#!/usr/bin/env python3
"""
Reality Check 3 (RC3) — SHH inhibition paradox test
====================================================

3 arms x 5 seeds (paired by seed):
  Arm A: Control
  Arm B: SHH inhibition only (vismodegib proxy)
  Arm C: SHH inhibition + drug

Timeline (single continuous run):
  0      -> T_PRE=20160   : barrier establishment
  T_PRE  -> T_END=40320   : intervention phase

Required (quorum >= 4/5 seeds):
  Final tumor rank: C < A < B
  Explicit checks:  B > A, C < A
"""
from __future__ import annotations

import argparse
import math
import os
import shutil
import subprocess
import sys
import textwrap
import time
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(PROJECT_ROOT))

from python.wrapper.output_parser import OutputParser

BINARY = PROJECT_ROOT / "stroma_world"
BASE_CONFIG = PROJECT_ROOT / "config" / "PhysiCell_settings.xml"

SEEDS = [42, 99, 137, 256, 1001]
NUM_REPLICATES = len(SEEDS)
QUORUM = 4

T_PRE = 20160.0
T_END = 40320.0
MAX_TIME = T_END
SAVE_INTERVAL = 360

DRUG_CONCENTRATION = 1.0

SLURM_PARTITION = os.environ.get("RC3_SLURM_PARTITION", os.environ.get("SLURM_PARTITION", "cpu384g"))
SLURM_CPUS = int(os.environ.get("RC3_SLURM_CPUS", os.environ.get("SLURM_CPUS", "32")))
SLURM_MEM = os.environ.get("RC3_SLURM_MEM", os.environ.get("SLURM_MEM", "0"))
SLURM_TIME = os.environ.get("RC3_SLURM_TIME", "06:00:00")
SLURM_POLL_INTERVAL = 20

TERMINAL_STATES = {
    "COMPLETED", "FAILED", "CANCELLED", "TIMEOUT",
    "OUT_OF_MEMORY", "NODE_FAIL", "PREEMPTED",
    "BOOT_FAIL", "DEADLINE",
}


@dataclass(frozen=True)
class ArmSpec:
    key: str
    name: str
    shh_start_time: float
    shh_strength: float
    drug_start_time: float
    drug_end_time: float
    drug_concentration: float


ARMS: List[ArmSpec] = [
    ArmSpec(
        key="A",
        name="Control",
        shh_start_time=1e18,
        shh_strength=0.0,
        drug_start_time=1e18,
        drug_end_time=1e18,
        drug_concentration=0.0,
    ),
    ArmSpec(
        key="B",
        name="SHH inhibition only",
        shh_start_time=T_PRE,
        shh_strength=1.0,
        drug_start_time=1e18,
        drug_end_time=1e18,
        drug_concentration=0.0,
    ),
    ArmSpec(
        key="C",
        name="SHH inhibition + drug",
        shh_start_time=T_PRE,
        shh_strength=1.0,
        drug_start_time=T_PRE,
        drug_end_time=T_END,
        drug_concentration=DRUG_CONCENTRATION,
    ),
]


@dataclass
class Snapshot:
    time: float
    n_tumor: int
    peri_ecm: float
    drug_at_tumor: float


@dataclass
class RunSpec:
    arm: ArmSpec
    replicate_index: int
    seed: int
    run_dir: Path
    config_path: Path


@dataclass
class RunResult:
    spec: RunSpec
    success: bool
    state: str = "UNKNOWN"
    exit_code: int = -1
    snap_pre: Optional[Snapshot] = None
    snap_final: Optional[Snapshot] = None


def _set_node(root: ET.Element, xpath: str, value: str) -> None:
    node = root.find(xpath)
    if node is not None:
        node.text = value


def _enforce_dirichlet(root: ET.Element, var_name: str, value: str) -> None:
    var = root.find(f".//microenvironment_setup/variable[@name='{var_name}']")
    if var is None:
        return
    dbc = var.find("./Dirichlet_boundary_condition")
    if dbc is not None:
        dbc.text = value
        dbc.set("enabled", "true")
    for bv in var.findall("./Dirichlet_options/boundary_value"):
        bv.text = value
        bv.set("enabled", "true")


def _patch_config(src: Path, dst: Path, output_dir: Path, seed: int, arm: ArmSpec) -> None:
    tree = ET.parse(src)
    root = tree.getroot()

    _set_node(root, "./save/folder", str(output_dir))
    _set_node(root, ".//overall/max_time", str(MAX_TIME))
    _set_node(root, ".//options/random_seed", str(seed))
    _set_node(root, ".//parallel/omp_num_threads", str(SLURM_CPUS))

    for n in root.findall(".//save//interval"):
        n.text = str(SAVE_INTERVAL)

    # Drug scheduling (reuses RC2 mechanism in custom.cpp)
    _set_node(root, ".//user_parameters/drug_start_time", str(arm.drug_start_time))
    _set_node(root, ".//user_parameters/drug_end_time", str(arm.drug_end_time))
    _set_node(root, ".//user_parameters/drug_concentration", str(arm.drug_concentration))

    # SHH inhibition intervention (new RC3 params)
    _set_node(root, ".//user_parameters/shh_inhibition_start_time", str(arm.shh_start_time))
    _set_node(root, ".//user_parameters/shh_inhibition_strength", str(arm.shh_strength))

    _enforce_dirichlet(root, "oxygen", "38")
    _enforce_dirichlet(root, "tgfb", "0")
    _enforce_dirichlet(root, "shh", "0")
    _enforce_dirichlet(root, "drug", "0")
    _enforce_dirichlet(root, "ecm_density", "0")

    tree.write(dst, encoding="utf-8", xml_declaration=True)


def _write_slurm_script(spec: RunSpec) -> Path:
    rep = spec.replicate_index + 1
    script = spec.run_dir / "run.slurm.sh"
    script.write_text(textwrap.dedent(f"""\
        #!/bin/bash
        #SBATCH --job-name=rc3_{spec.arm.key}_r{rep}_s{spec.seed}
        #SBATCH --partition={SLURM_PARTITION}
        #SBATCH --nodes=1
        #SBATCH --ntasks-per-node=1
        #SBATCH --cpus-per-task={SLURM_CPUS}
        #SBATCH --mem={SLURM_MEM}
        #SBATCH --time={SLURM_TIME}
        #SBATCH --output={spec.run_dir}/slurm_%j.out
        #SBATCH --error={spec.run_dir}/slurm_%j.err

        set -euo pipefail
        export OMP_NUM_THREADS={SLURM_CPUS}
        export OMP_PROC_BIND=spread
        export OMP_PLACES=cores

        cd {PROJECT_ROOT}
        echo "=== RC3 Arm {spec.arm.key} ({spec.arm.name}) Rep {rep} seed={spec.seed} ==="
        echo "Start: $(date)"
        {BINARY} {spec.config_path}
        echo "End:   $(date)"
    """), encoding="utf-8")
    script.chmod(0o755)
    return script


def _query_job(job_id: str) -> Tuple[str, int, bool]:
    try:
        r = subprocess.run(
            ["sacct", "-j", job_id, "--format=State,ExitCode", "--noheader", "--parsable2"],
            capture_output=True, text=True, check=False, timeout=30,
        )
        for line in r.stdout.strip().splitlines():
            parts = line.strip().split("|")
            if len(parts) >= 2:
                state = parts[0].strip().split()[0] if parts[0].strip() else "UNKNOWN"
                try:
                    exit_code = int(parts[1].split(":")[0])
                except (ValueError, IndexError):
                    exit_code = -1
                return state, exit_code, state in TERMINAL_STATES
    except Exception:
        pass

    try:
        r = subprocess.run(
            ["squeue", "-j", job_id, "--format=%T", "--noheader"],
            capture_output=True, text=True, check=False, timeout=15,
        )
        state = r.stdout.strip()
        if not state:
            return "COMPLETED", 0, True
        return state, 0, False
    except Exception:
        return "UNKNOWN", -1, False


def _sample_nearest(points: np.ndarray, voxel_coords: np.ndarray, values: np.ndarray) -> np.ndarray:
    if points.size == 0 or voxel_coords.size == 0:
        return np.array([], dtype=float)
    out = np.empty(points.shape[0], dtype=float)
    chunk = 500
    for i in range(0, points.shape[0], chunk):
        p = points[i : i + chunk]
        d2 = np.sum((p[:, None, :] - voxel_coords[None, :, :]) ** 2, axis=2)
        out[i : i + chunk] = values[np.argmin(d2, axis=1)]
    return out


def _parse_snapshot(output_dir: Path, target_time: float) -> Optional[Snapshot]:
    parser = OutputParser(output_dir)
    xmls = sorted(output_dir.glob("output*.xml"))
    if not xmls:
        return None

    snaps: List[Tuple[float, dict]] = []
    for xml in xmls:
        s = parser._read_physicell_xml(xml)
        snaps.append((float(s["time"]), s))
    if not snaps:
        return None

    snap_time, snap = min(snaps, key=lambda x: abs(x[0] - target_time))

    matrix = snap["cell_matrix"]
    labels = snap["label_name_map"]
    micro_coords = snap["micro_coords"]
    micro_values = snap["micro_values"]

    def _row(name: str) -> Optional[np.ndarray]:
        e = labels.get(name)
        if e is None:
            return None
        idx = int(e["index"])
        if idx < 0 or idx >= matrix.shape[0]:
            return None
        return matrix[idx, :]

    cell_type = _row("cell_type")
    dead = _row("dead")
    death_model = _row("current_death_model")

    n_cells = matrix.shape[1]
    live_mask = np.ones(n_cells, dtype=bool)
    if dead is not None:
        live_mask &= dead <= 0.5
    if death_model is not None:
        live_mask &= np.rint(death_model).astype(int) != 100

    ctype_int = np.rint(cell_type).astype(int) if cell_type is not None else np.full(n_cells, -1, dtype=int)
    live_tumor = live_mask & (ctype_int == 0)
    n_tumor = int(np.sum(live_tumor))

    pos = parser._get_positions(matrix, labels)
    tumor_pos = pos[live_tumor] if (pos.size and live_tumor.any()) else np.empty((0, 3))
    centroid = np.nanmean(tumor_pos, axis=0) if tumor_pos.shape[0] > 0 else np.zeros(3)
    tumor_radius = float(np.max(np.linalg.norm(tumor_pos - centroid, axis=1))) if tumor_pos.shape[0] > 1 else 0.0

    peri_ecm = math.nan
    ecm = micro_values.get("ecm_density")
    if ecm is not None and micro_coords.size > 0 and tumor_pos.shape[0] > 0:
        vd = np.linalg.norm(micro_coords - centroid, axis=1)
        shell = (vd >= tumor_radius) & (vd <= tumor_radius + 100.0)
        if shell.any():
            peri_ecm = float(np.nanmean(ecm[shell]))

    drug_at_tumor = math.nan
    drug = micro_values.get("drug")
    if drug is not None and micro_coords.size > 0 and tumor_pos.shape[0] > 0:
        drug_vals = _sample_nearest(tumor_pos, micro_coords, drug)
        if drug_vals.size > 0:
            drug_at_tumor = float(np.nanmean(drug_vals))

    return Snapshot(
        time=snap_time,
        n_tumor=n_tumor,
        peri_ecm=peri_ecm,
        drug_at_tumor=drug_at_tumor,
    )


def _make_run_specs(work_dir: Path) -> List[RunSpec]:
    specs: List[RunSpec] = []
    for arm in ARMS:
        for i, seed in enumerate(SEEDS):
            run_dir = work_dir / f"arm_{arm.key}" / f"replicate_{i+1:02d}_seed{seed}"
            cfg = run_dir / "config.xml"
            specs.append(RunSpec(arm=arm, replicate_index=i, seed=seed, run_dir=run_dir, config_path=cfg))
    return specs


def _prepare_runs(work_dir: Path) -> List[RunSpec]:
    if work_dir.exists():
        shutil.rmtree(work_dir)
    work_dir.mkdir(parents=True)

    specs = _make_run_specs(work_dir)
    for spec in specs:
        spec.run_dir.mkdir(parents=True, exist_ok=True)
        out_dir = spec.run_dir / "output"
        out_dir.mkdir(parents=True, exist_ok=True)

        _patch_config(BASE_CONFIG, spec.config_path, out_dir, spec.seed, spec.arm)

        local_cfg = spec.run_dir / "config"
        local_cfg.mkdir(exist_ok=True)
        shutil.copy(PROJECT_ROOT / "config" / "tumor_calibration_knobs.json", local_cfg / "tumor_calibration_knobs.json")
        for extra in ["gene_params_default.json"]:
            src = PROJECT_ROOT / "config" / extra
            if src.exists():
                shutil.copy(src, local_cfg / extra)

    return specs


def _dry_run_report(specs: List[RunSpec]) -> None:
    print("=" * 78)
    print("RC3 DRY RUN PLAN (no jobs submitted)")
    print("=" * 78)
    print(f"Total runs: {len(specs)}  ({len(ARMS)} arms x {NUM_REPLICATES} replicates)")
    print(f"Seeds: {SEEDS}")
    print(f"Timeline: barrier 0-{T_PRE:.0f} min, intervention {T_PRE:.0f}-{T_END:.0f} min")
    print()
    for arm in ARMS:
        print(f"Arm {arm.key}: {arm.name}")
        print(f"  shh_start={arm.shh_start_time:.0f}, shh_strength={arm.shh_strength:.2f}, "
              f"drug_start={arm.drug_start_time:.0f}, drug_end={arm.drug_end_time:.0f}, "
              f"drug_conc={arm.drug_concentration:.2f}")
        for s in [x for x in specs if x.arm.key == arm.key]:
            print(f"    - {s.run_dir}")
        print()


def _evaluate_run(spec: RunSpec, state: str, exit_code: int) -> RunResult:
    ok = state == "COMPLETED" and exit_code == 0
    r = RunResult(spec=spec, success=ok, state=state, exit_code=exit_code)
    if not ok:
        return r

    output_dir = spec.run_dir / "output"
    if not list(output_dir.glob("output*.xml")):
        r.success = False
        return r

    r.snap_pre = _parse_snapshot(output_dir, T_PRE)
    r.snap_final = _parse_snapshot(output_dir, T_END)
    if r.snap_pre is None or r.snap_final is None:
        r.success = False
    return r


def _median(xs: List[float]) -> float:
    return float(np.median(np.array(xs, dtype=float))) if xs else float("nan")


def _report(results: List[RunResult]) -> int:
    print("\n" + "=" * 78)
    print("RC3 RESULTS")
    print("=" * 78)

    by_arm: Dict[str, List[RunResult]] = {a.key: [] for a in ARMS}
    by_seed_arm: Dict[Tuple[int, str], RunResult] = {}
    for r in results:
        by_arm[r.spec.arm.key].append(r)
        by_seed_arm[(r.spec.seed, r.spec.arm.key)] = r

    for arm in ARMS:
        arm_results = [r for r in by_arm[arm.key] if r.success and r.snap_final is not None]
        finals = [float(r.snap_final.n_tumor) for r in arm_results]
        peris = [float(r.snap_final.peri_ecm) for r in arm_results if math.isfinite(r.snap_final.peri_ecm)]
        print(f"\nArm {arm.key} ({arm.name})")
        print(f"  successful runs: {len(arm_results)}/{NUM_REPLICATES}")
        print(f"  final tumor counts: {finals if finals else 'n/a'}")
        print(f"  median final tumor count: {_median(finals):.2f}")
        if peris:
            print(f"  median final peritumoral ECM: {_median(peris):.4f}")
        if arm.key == "C":
            drug_vals = [float(r.snap_final.drug_at_tumor) for r in arm_results if math.isfinite(r.snap_final.drug_at_tumor)]
            if drug_vals:
                print(f"  median final drug near tumor: {_median(drug_vals):.4f}")

    rank_count = 0
    b_gt_a_count = 0
    c_lt_a_count = 0
    c_lt_b_count = 0
    compared = 0

    print("\nPer-seed paired final tumor counts:")
    for seed in SEEDS:
        ra = by_seed_arm.get((seed, "A"))
        rb = by_seed_arm.get((seed, "B"))
        rc = by_seed_arm.get((seed, "C"))
        if not (ra and rb and rc and ra.success and rb.success and rc.success and ra.snap_final and rb.snap_final and rc.snap_final):
            print(f"  seed={seed}: missing/failed run(s)")
            continue

        a = ra.snap_final.n_tumor
        b = rb.snap_final.n_tumor
        c = rc.snap_final.n_tumor
        compared += 1

        b_gt_a = b > a
        c_lt_a = c < a
        c_lt_b = c < b
        full_rank = c < a < b

        b_gt_a_count += int(b_gt_a)
        c_lt_a_count += int(c_lt_a)
        c_lt_b_count += int(c_lt_b)
        rank_count += int(full_rank)

        print(f"  seed={seed}: A={a}, B={b}, C={c} | "
              f"B>A={'Y' if b_gt_a else 'N'}  C<A={'Y' if c_lt_a else 'N'}  C<A<B={'Y' if full_rank else 'N'}")

    print("\nRank checks (quorum >= 4/5):")
    print(f"  C < A < B : {rank_count}/{NUM_REPLICATES}")
    print(f"  B > A     : {b_gt_a_count}/{NUM_REPLICATES}")
    print(f"  C < A     : {c_lt_a_count}/{NUM_REPLICATES}")

    pass_full = rank_count >= QUORUM
    pass_ba = b_gt_a_count >= QUORUM
    pass_ca = c_lt_a_count >= QUORUM
    overall_pass = pass_full and pass_ba and pass_ca

    print("\nOutcome:")
    print(f"  RC3 {'PASS' if overall_pass else 'FAIL'}")

    if not pass_ba:
        print("  Hint: If B ≤ A, mechanical confinement may be too weak or SHH->ECM coupling is broken.")
    if not pass_ca:
        print("  Hint: If C ≥ A, drug may not be benefiting from thinner stroma (check ECM->drug diffusion coupling).")
    if c_lt_b_count < QUORUM:
        print("  Hint: If C ≥ B, drug may be too weak or schedule is wrong.")

    if compared < NUM_REPLICATES:
        print(f"  Note: only {compared}/{NUM_REPLICATES} seeds had complete paired data.")

    return 0 if overall_pass else 1


def main() -> int:
    parser = argparse.ArgumentParser(description="Run Reality Check 3")
    parser.add_argument("--dry-run", action="store_true", help="Plan all runs and print paths without submitting")
    parser.add_argument("--evaluate-only", action="store_true",
                        help="Evaluate existing output (no submission, no cleanup)")
    args = parser.parse_args()

    if not BINARY.exists():
        print(f"ERROR: binary not found at {BINARY}")
        return 1
    if not BASE_CONFIG.exists():
        print(f"ERROR: config not found at {BASE_CONFIG}")
        return 1

    work_dir = PROJECT_ROOT / "build" / "reality_check_3"

    print("=" * 78)
    print("REALITY CHECK 3 — SHH inhibition paradox")
    print(f"Arms: {len(ARMS)} | Replicates/arm: {NUM_REPLICATES} | Total runs: {len(ARMS) * NUM_REPLICATES}")
    print(f"Timeline: barrier 0-{T_PRE:.0f} min | intervention {T_PRE:.0f}-{T_END:.0f} min")
    print(f"SLURM: partition={SLURM_PARTITION} cpus={SLURM_CPUS} mem={SLURM_MEM} time={SLURM_TIME}")
    print(f"Quorum: >= {QUORUM}/{NUM_REPLICATES}")
    print("=" * 78)

    if args.evaluate_only:
        if not work_dir.exists():
            print(f"ERROR: work_dir {work_dir} does not exist — nothing to evaluate")
            return 1
        specs = _make_run_specs(work_dir)
        print("Evaluate-only mode — skipping submission, reading existing output ...")
        results: List[RunResult] = []
        for spec in specs:
            print(f"Evaluating arm {spec.arm.key} rep {spec.replicate_index+1} seed={spec.seed} ...")
            results.append(_evaluate_run(spec, "COMPLETED", 0))
        return _report(results)

    specs = _prepare_runs(work_dir)

    if args.dry_run:
        _dry_run_report(specs)
        return 0

    job_to_spec: Dict[str, RunSpec] = {}
    for spec in specs:
        script = _write_slurm_script(spec)
        r = subprocess.run(["sbatch", "--parsable", str(script)], capture_output=True, text=True, check=False)
        if r.returncode != 0:
            print(f"ERROR: sbatch failed for arm {spec.arm.key} rep {spec.replicate_index+1} seed={spec.seed}: {r.stderr.strip()}")
            return 1
        job_id = r.stdout.strip().split(";")[0]
        job_to_spec[job_id] = spec
        print(f"Submitted arm {spec.arm.key} rep {spec.replicate_index+1} seed={spec.seed} -> job {job_id}")

    print(f"\nSubmitted {len(job_to_spec)} jobs. Polling ...\n")

    done: Dict[str, Tuple[str, int]] = {}
    t0 = time.perf_counter()
    while len(done) < len(job_to_spec):
        time.sleep(SLURM_POLL_INTERVAL)
        elapsed = time.perf_counter() - t0
        for job_id, spec in job_to_spec.items():
            if job_id in done:
                continue
            state, exit_code, terminal = _query_job(job_id)
            if terminal:
                done[job_id] = (state, exit_code)
                tag = "OK" if (state == "COMPLETED" and exit_code == 0) else "XX"
                print(f"{tag} arm {spec.arm.key} rep {spec.replicate_index+1} seed={spec.seed}: {state} exit={exit_code} [{elapsed:.0f}s]")

    total_wall = time.perf_counter() - t0
    print(f"\nAll jobs finished in {total_wall:.1f}s ({total_wall/60:.1f} min)")

    results: List[RunResult] = []
    for _, spec in sorted(job_to_spec.items(), key=lambda kv: (kv[1].arm.key, kv[1].replicate_index)):
        state, exit_code = done.get(_, ("UNKNOWN", -1))
        print(f"Evaluating arm {spec.arm.key} rep {spec.replicate_index+1} seed={spec.seed} ...")
        results.append(_evaluate_run(spec, state, exit_code))

    return _report(results)


if __name__ == "__main__":
    raise SystemExit(main())
