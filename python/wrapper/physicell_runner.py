from __future__ import annotations

import concurrent.futures
import json
import logging
import os
import re
import shutil
import signal
import subprocess
import time
import uuid
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

try:
    import resource  # POSIX-only; unavailable on Windows.
except ImportError:  # pragma: no cover
    resource = None  # type: ignore[assignment]


LOGGER = logging.getLogger(__name__)


@dataclass
class SimulationResult:
    run_dir: Path
    exit_code: int
    wall_time: float
    success: bool
    output_files: Dict[int, Path] = field(default_factory=dict)


class PhysiCellRunner:
    def __init__(self, binary_path, config_path, output_dir, timeout_seconds=3600):
        self.binary_path = Path(binary_path).expanduser().resolve()
        self.config_path = Path(config_path).expanduser().resolve()
        self.output_dir = Path(output_dir).expanduser().resolve()
        self.timeout_seconds = float(timeout_seconds)

        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.cleanup_failed_runs = True
        self.poll_interval_seconds = 1.0

        # Optional memory limit in MB. Set env PHYSICELL_MEMORY_LIMIT_MB=...
        self.memory_limit_mb = self._load_memory_limit_from_env()

        if not self.binary_path.exists():
            raise FileNotFoundError(f"PhysiCell binary not found: {self.binary_path}")
        if not os.access(self.binary_path, os.X_OK):
            raise PermissionError(f"PhysiCell binary is not executable: {self.binary_path}")
        if not self.config_path.exists():
            raise FileNotFoundError(f"PhysiCell config not found: {self.config_path}")

    def run(self, intervention_json_path) -> SimulationResult:
        start_time = time.perf_counter()
        run_dir = self._new_run_dir()
        output_subdir = run_dir / "output"
        output_subdir.mkdir(parents=True, exist_ok=True)

        config_copy = run_dir / "config.xml"
        intervention_copy = run_dir / "intervention.json"
        stdout_log = run_dir / "stdout.log"
        stderr_log = run_dir / "stderr.log"

        try:
            self._copy_and_patch_config(self.config_path, config_copy, output_subdir)
            self._materialize_intervention(intervention_json_path, intervention_copy)

            primary_cmd = self._primary_command(config_copy, intervention_copy, output_subdir)
            exit_code, timed_out, mem_exceeded = self._execute_command(
                cmd=primary_cmd,
                stdout_log=stdout_log,
                stderr_log=stderr_log,
            )
            wall_time = time.perf_counter() - start_time

            output_files = self._collect_output_files(output_subdir)
            if exit_code != 0 and not timed_out and not mem_exceeded and not output_files:
                LOGGER.warning(
                    "Primary CLI failed for %s (exit=%s). Retrying with legacy positional CLI.",
                    run_dir,
                    exit_code,
                )
                # Keep first attempt logs and append legacy attempt logs.
                legacy_stdout = run_dir / "stdout_legacy.log"
                legacy_stderr = run_dir / "stderr_legacy.log"
                legacy_cmd = self._legacy_command(config_copy, intervention_copy)
                exit_code, timed_out, mem_exceeded = self._execute_command(
                    cmd=legacy_cmd,
                    stdout_log=legacy_stdout,
                    stderr_log=legacy_stderr,
                )
                wall_time = time.perf_counter() - start_time
                output_files = self._collect_output_files(output_subdir)

            success = (exit_code == 0) and self._output_has_files(output_subdir)
            result = SimulationResult(
                run_dir=run_dir,
                exit_code=exit_code,
                wall_time=wall_time,
                success=success,
                output_files=output_files,
            )

            if not success and self.cleanup_failed_runs:
                reason = "timeout" if timed_out else ("memory_limit" if mem_exceeded else "failed_exit")
                self._cleanup_failed_run(result, reason=reason)
            return result

        except Exception as exc:
            wall_time = time.perf_counter() - start_time
            LOGGER.exception("Simulation run failed before completion in %s: %s", run_dir, exc)
            result = SimulationResult(
                run_dir=run_dir,
                exit_code=-1,
                wall_time=wall_time,
                success=False,
                output_files={},
            )
            if self.cleanup_failed_runs:
                self._cleanup_failed_run(result, reason=f"exception:{type(exc).__name__}")
            return result

    def run_batch(self, intervention_jsons: list, max_parallel: int = 4) -> List[SimulationResult]:
        if max_parallel < 1:
            raise ValueError("max_parallel must be >= 1")

        results: List[Optional[SimulationResult]] = [None] * len(intervention_jsons)

        with concurrent.futures.ThreadPoolExecutor(max_workers=max_parallel) as pool:
            future_map = {
                pool.submit(self.run, intervention_json): idx
                for idx, intervention_json in enumerate(intervention_jsons)
            }
            for future in concurrent.futures.as_completed(future_map):
                idx = future_map[future]
                try:
                    results[idx] = future.result()
                except Exception as exc:  # pragma: no cover
                    LOGGER.exception("Unhandled batch worker failure at index %d: %s", idx, exc)
                    failed_dir = self._new_run_dir(prefix="batch_failure")
                    results[idx] = SimulationResult(
                        run_dir=failed_dir,
                        exit_code=-1,
                        wall_time=0.0,
                        success=False,
                        output_files={},
                    )

        return [r for r in results if r is not None]

    def run_batch_slurm(
        self,
        intervention_jsons: list,
        max_parallel: int = 4,
        poll_interval_seconds: int = 15,
        sbatch_args: Optional[List[str]] = None,
    ) -> List[SimulationResult]:
        if max_parallel < 1:
            raise ValueError("max_parallel must be >= 1")
        sbatch_args = sbatch_args or []

        prepared: List[Dict[str, Any]] = []
        for idx, intervention_json in enumerate(intervention_jsons):
            item = self._prepare_slurm_run(intervention_json)
            item["index"] = idx
            prepared.append(item)
        results: List[Optional[SimulationResult]] = [None] * len(prepared)

        pending = list(prepared)
        active: Dict[str, Dict[str, Any]] = {}

        while pending or active:
            while pending and len(active) < max_parallel:
                item = pending.pop(0)
                idx = item["index"]
                run_dir: Path = item["run_dir"]
                script_path = self._write_slurm_script(
                    run_dir=run_dir,
                    config_path=item["config_path"],
                    intervention_path=item["intervention_path"],
                    output_path=item["output_path"],
                )
                cmd = ["sbatch", "--parsable", *sbatch_args, str(script_path)]
                submit = subprocess.run(cmd, capture_output=True, text=True, check=False)
                if submit.returncode != 0:
                    LOGGER.error("sbatch failed for %s: %s", run_dir, submit.stderr.strip())
                    result = SimulationResult(run_dir=run_dir, exit_code=-1, wall_time=0.0, success=False, output_files={})
                    if self.cleanup_failed_runs:
                        self._cleanup_failed_run(result, reason="sbatch_submit_failed")
                    results[idx] = result
                    continue

                job_id = submit.stdout.strip().split(";")[0]
                active[job_id] = {
                    "index": idx,
                    "run_dir": run_dir,
                    "output_path": item["output_path"],
                    "start_time": time.perf_counter(),
                }
                LOGGER.info("Submitted SLURM job %s for run %s", job_id, run_dir)

            if not active:
                continue

            time.sleep(max(1, poll_interval_seconds))

            for job_id in list(active.keys()):
                ctx = active[job_id]
                elapsed = time.perf_counter() - ctx["start_time"]

                # Local timeout guard even when the scheduler allows longer.
                if elapsed > self.timeout_seconds:
                    subprocess.run(["scancel", job_id], capture_output=True, text=True, check=False)
                    result = SimulationResult(
                        run_dir=ctx["run_dir"],
                        exit_code=-1,
                        wall_time=elapsed,
                        success=False,
                        output_files={},
                    )
                    if self.cleanup_failed_runs:
                        self._cleanup_failed_run(result, reason="slurm_timeout")
                    results[ctx["index"]] = result
                    del active[job_id]
                    continue

                state, exit_code, terminal = self._query_slurm_state(job_id)
                if not terminal:
                    continue

                output_files = self._collect_output_files(ctx["output_path"])
                success = (state == "COMPLETED") and (exit_code == 0) and self._output_has_files(ctx["output_path"])
                result = SimulationResult(
                    run_dir=ctx["run_dir"],
                    exit_code=exit_code,
                    wall_time=elapsed,
                    success=success,
                    output_files=output_files,
                )
                if not success and self.cleanup_failed_runs:
                    self._cleanup_failed_run(result, reason=f"slurm_{state.lower()}")
                results[ctx["index"]] = result
                del active[job_id]

        return [r for r in results if r is not None]

    def _load_memory_limit_from_env(self) -> Optional[int]:
        raw = os.environ.get("PHYSICELL_MEMORY_LIMIT_MB")
        if not raw:
            return None
        try:
            value = int(raw)
            return value if value > 0 else None
        except ValueError:
            LOGGER.warning("Ignoring invalid PHYSICELL_MEMORY_LIMIT_MB value: %s", raw)
            return None

    def _new_run_dir(self, prefix: str = "run") -> Path:
        ts = datetime.utcnow().strftime("%Y%m%d_%H%M%S")
        uid = uuid.uuid4().hex[:12]
        run_dir = self.output_dir / f"{prefix}_{ts}_{uid}"
        run_dir.mkdir(parents=True, exist_ok=False)
        return run_dir

    def _copy_and_patch_config(self, src: Path, dst: Path, output_subdir: Path) -> None:
        tree = ET.parse(src)
        root = tree.getroot()

        # Prefer the canonical save/folder setting when present.
        folder_node = root.find("./save/folder")
        if folder_node is None:
            folder_node = root.find(".//folder")
        if folder_node is not None:
            folder_node.text = str(output_subdir)

        tree.write(dst, encoding="utf-8", xml_declaration=True)

    def _materialize_intervention(self, intervention_source: Any, dst: Path) -> None:
        payload: Any

        if isinstance(intervention_source, (str, os.PathLike, Path)):
            src_path = Path(intervention_source).expanduser()
            if src_path.exists():
                with src_path.open("r", encoding="utf-8") as f:
                    payload = json.load(f)
            else:
                # Allow inline JSON string as input for convenience.
                payload = json.loads(str(intervention_source))
        elif isinstance(intervention_source, (dict, list)):
            payload = intervention_source
        else:
            raise TypeError(f"Unsupported intervention source type: {type(intervention_source)}")

        with dst.open("w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2, sort_keys=True)
            f.write("\n")

    def _primary_command(self, config_path: Path, intervention_path: Path, output_path: Path) -> List[str]:
        return [
            str(self.binary_path),
            "--config",
            str(config_path),
            "--intervention",
            str(intervention_path),
            "--output",
            str(output_path),
        ]

    def _legacy_command(self, config_path: Path, intervention_path: Path) -> List[str]:
        return [str(self.binary_path), str(config_path), str(intervention_path)]

    def _execute_command(self, cmd: List[str], stdout_log: Path, stderr_log: Path) -> Tuple[int, bool, bool]:
        timed_out = False
        mem_exceeded = False

        with stdout_log.open("a", encoding="utf-8") as out_f, stderr_log.open("a", encoding="utf-8") as err_f:
            preexec = self._build_preexec_fn()
            process = subprocess.Popen(
                cmd,
                stdout=out_f,
                stderr=err_f,
                preexec_fn=preexec,
                text=True,
            )
            start = time.perf_counter()

            while True:
                ret = process.poll()
                if ret is not None:
                    return ret, timed_out, mem_exceeded

                elapsed = time.perf_counter() - start
                if elapsed > self.timeout_seconds:
                    timed_out = True
                    self._terminate_process_tree(process)
                    return -1, timed_out, mem_exceeded

                if self.memory_limit_mb is not None:
                    rss_mb = self._read_rss_mb(process.pid)
                    if rss_mb is not None and rss_mb > self.memory_limit_mb:
                        mem_exceeded = True
                        self._terminate_process_tree(process)
                        return -1, timed_out, mem_exceeded

                time.sleep(self.poll_interval_seconds)

    def _build_preexec_fn(self):
        if os.name != "posix":
            return None

        memory_limit_mb = self.memory_limit_mb

        def _preexec():
            os.setsid()
            if memory_limit_mb is not None and resource is not None:
                limit_bytes = int(memory_limit_mb) * 1024 * 1024
                try:
                    resource.setrlimit(resource.RLIMIT_AS, (limit_bytes, limit_bytes))
                except (ValueError, OSError):
                    pass

        return _preexec

    def _terminate_process_tree(self, process: subprocess.Popen) -> None:
        if process.poll() is not None:
            return
        try:
            if os.name == "posix":
                os.killpg(process.pid, signal.SIGTERM)
            else:  # pragma: no cover
                process.terminate()
            process.wait(timeout=10)
        except Exception:
            try:
                if os.name == "posix":
                    os.killpg(process.pid, signal.SIGKILL)
                else:  # pragma: no cover
                    process.kill()
            except Exception:
                pass

    def _read_rss_mb(self, pid: int) -> Optional[float]:
        status_file = Path(f"/proc/{pid}/status")
        if not status_file.exists():
            return None
        try:
            for line in status_file.read_text(encoding="utf-8", errors="ignore").splitlines():
                if line.startswith("VmRSS:"):
                    parts = line.split()
                    kb = float(parts[1])
                    return kb / 1024.0
        except Exception:
            return None
        return None

    def _collect_output_files(self, output_dir: Path) -> Dict[int, Path]:
        files: Dict[int, Path] = {}
        if not output_dir.exists():
            return files

        xml_pattern = re.compile(r"^output(\d+)\.xml$")
        for path in sorted(output_dir.glob("output*.xml")):
            match = xml_pattern.match(path.name)
            if match:
                files[int(match.group(1))] = path

        # Fallback mapping if no periodic outputs were produced.
        if not files:
            initial = output_dir / "initial.xml"
            final = output_dir / "final.xml"
            if initial.exists():
                files[0] = initial
            if final.exists():
                files[max(files.keys(), default=0) + 1] = final
        return files

    def _output_has_files(self, output_dir: Path) -> bool:
        if not output_dir.exists():
            return False
        return any(output_dir.iterdir())

    def _cleanup_failed_run(self, result: SimulationResult, reason: str) -> None:
        run_dir = result.run_dir
        output_dir = run_dir / "output"
        if output_dir.exists():
            shutil.rmtree(output_dir, ignore_errors=True)
        marker = run_dir / "FAILED.txt"
        marker.write_text(
            f"reason={reason}\nexit_code={result.exit_code}\nwall_time={result.wall_time:.6f}\n",
            encoding="utf-8",
        )

    def _prepare_slurm_run(self, intervention_source: Any) -> Dict[str, Any]:
        run_dir = self._new_run_dir(prefix="slurm_run")
        output_subdir = run_dir / "output"
        output_subdir.mkdir(parents=True, exist_ok=True)
        config_copy = run_dir / "config.xml"
        intervention_copy = run_dir / "intervention.json"

        self._copy_and_patch_config(self.config_path, config_copy, output_subdir)
        self._materialize_intervention(intervention_source, intervention_copy)

        return {
            "run_dir": run_dir,
            "config_path": config_copy,
            "intervention_path": intervention_copy,
            "output_path": output_subdir,
            "index": None,  # filled by caller
        }

    def _write_slurm_script(
        self,
        run_dir: Path,
        config_path: Path,
        intervention_path: Path,
        output_path: Path,
    ) -> Path:
        script_path = run_dir / "run.slurm.sh"
        cmd_primary = " ".join(
            [
                self._shq(self.binary_path),
                "--config",
                self._shq(config_path),
                "--intervention",
                self._shq(intervention_path),
                "--output",
                self._shq(output_path),
            ]
        )
        cmd_legacy = " ".join(
            [self._shq(self.binary_path), self._shq(config_path), self._shq(intervention_path)]
        )

        script = "\n".join(
            [
                "#!/usr/bin/env bash",
                "set -euo pipefail",
                f"cd {self._shq(Path.cwd())}",
                f"{cmd_primary} || {cmd_legacy}",
                "",
            ]
        )
        script_path.write_text(script, encoding="utf-8")
        script_path.chmod(0o755)
        return script_path

    def _query_slurm_state(self, job_id: str) -> Tuple[str, int, bool]:
        cmd = ["sacct", "-j", job_id, "--format=JobIDRaw,State,ExitCode", "-n", "-P"]
        probe = subprocess.run(cmd, capture_output=True, text=True, check=False)
        if probe.returncode == 0 and probe.stdout.strip():
            state, exit_code, terminal = self._parse_sacct_output(job_id, probe.stdout)
            if state:
                return state, exit_code, terminal

        # Fallback: if sacct is slow/unavailable, use squeue for liveness.
        sq = subprocess.run(["squeue", "-j", job_id, "-h", "-o", "%T"], capture_output=True, text=True, check=False)
        if sq.returncode == 0 and sq.stdout.strip():
            return sq.stdout.strip().splitlines()[0], -1, False
        return "UNKNOWN", -1, True

    def _parse_sacct_output(self, job_id: str, output: str) -> Tuple[str, int, bool]:
        terminal_states = {
            "COMPLETED",
            "FAILED",
            "TIMEOUT",
            "CANCELLED",
            "OUT_OF_MEMORY",
            "NODE_FAIL",
            "PREEMPTED",
        }
        for line in output.strip().splitlines():
            fields = line.split("|")
            if len(fields) < 3:
                continue
            jid, state_raw, exit_raw = fields[0].strip(), fields[1].strip(), fields[2].strip()
            if jid != job_id:
                continue
            state = state_raw.split()[0]
            exit_code = self._parse_exit_code(exit_raw)
            return state, exit_code, state in terminal_states
        return "", -1, False

    def _parse_exit_code(self, raw: str) -> int:
        # SLURM format is usually "<code>:<signal>".
        if ":" in raw:
            code_part = raw.split(":", 1)[0]
            try:
                return int(code_part)
            except ValueError:
                return -1
        try:
            return int(raw)
        except ValueError:
            return -1

    def _shq(self, value: Path) -> str:
        text = str(value)
        return "'" + text.replace("'", "'\"'\"'") + "'"
