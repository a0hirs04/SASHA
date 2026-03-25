#!/usr/bin/env python3
"""
Launch 7 targeted RC2 contingency jobs (C1-C7) to hedge against
NO_EFFECTIVE_TREATMENT / WEAK_DRUG_EXPOSURE failures.

Default behavior:
  - Check Phoenix job state (default job id: 10674)
  - If still running: do nothing
  - If completed and verdict indicates PASS: do nothing
  - If completed and verdict indicates Ghost-Drug style failure: launch hedge set

Use --force-launch to submit hedge jobs immediately.
"""

from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import List

PROJECT_ROOT = Path(__file__).resolve().parent
RUNNER = PROJECT_ROOT / "run_reality_check_2.py"
DEFAULT_VERDICT = PROJECT_ROOT / "logs" / "rc2_regrowth_candidate_verdict_snippet.txt"
DEFAULT_WORK_ROOT = PROJECT_ROOT / "build" / "rc2_hedge_set"


@dataclass(frozen=True)
class HedgeJob:
    cid: str
    km: float
    ap: float
    thresh: float


HEDGE_JOBS: List[HedgeJob] = [
    HedgeJob("C1", 2.1, 0.00080, 0.025),
    HedgeJob("C2", 2.2, 0.00075, 0.025),
    HedgeJob("C3", 2.0, 0.00075, 0.030),
    HedgeJob("C4", 2.3, 0.00080, 0.025),
    HedgeJob("C5", 1.9, 0.00070, 0.030),
    HedgeJob("C6", 2.4, 0.00075, 0.020),
    HedgeJob("C7", 2.2, 0.00065, 0.030),
]


def _job_state(job_id: str) -> str:
    r = subprocess.run(
        ["sacct", "-j", job_id, "--format=State", "--noheader", "-P"],
        capture_output=True,
        text=True,
        check=False,
    )
    lines = [ln.strip() for ln in r.stdout.splitlines() if ln.strip()]
    if not lines:
        return "UNKNOWN"
    # first line is parent state
    return lines[0].split("|")[0].strip().upper()


def _verdict_text(path: Path) -> str:
    if not path.exists():
        return ""
    return path.read_text(encoding="utf-8", errors="ignore")


def _is_pass(verdict_text: str) -> bool:
    t = verdict_text.upper()
    return "TRUE BIOLOGICAL PASS" in t or re.search(r"\bPASS\b", t) is not None and "FALSE PASS" not in t


def _is_ghost_fail(verdict_text: str) -> bool:
    t = verdict_text.upper()
    trigger_tokens = [
        "NO_EFFECTIVE_TREATMENT",
        "WEAK_DRUG_EXPOSURE",
        "GHOST DRUG",
        "FALSE PASS",
    ]
    return any(tok in t for tok in trigger_tokens)


def _launch_job(job: HedgeJob, work_root: Path, seeds: str, quorum: int, detach: bool) -> None:
    job_prefix = f"rc2_hedge_{job.cid}_km{job.km}_ap{job.ap}_th{job.thresh}"
    work_dir = work_root / job_prefix
    log_dir = PROJECT_ROOT / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / f"{job_prefix}.log"

    cmd = [
        sys.executable,
        str(RUNNER),
        "--seeds",
        *seeds.split(),
        "--quorum",
        str(quorum),
        "--work-dir",
        str(work_dir),
        "--job-name-prefix",
        job_prefix,
        "--drug-kill-multiplier",
        str(job.km),
        "--abcb1-production-rate",
        str(job.ap),
        "--drug-stress-threshold",
        str(job.thresh),
    ]

    if detach:
        with log_path.open("w", encoding="utf-8") as fp:
            subprocess.Popen(
                cmd,
                cwd=str(PROJECT_ROOT),
                stdout=fp,
                stderr=subprocess.STDOUT,
            )
        print(f"[LAUNCHED] {job_prefix} (detached) -> {log_path}")
    else:
        print(f"[RUNNING] {job_prefix}")
        subprocess.run(cmd, cwd=str(PROJECT_ROOT), check=False)


def main() -> int:
    parser = argparse.ArgumentParser(description="Launch RC2 hedge set C1-C7")
    parser.add_argument("--phoenix-job-id", default="10674")
    parser.add_argument("--verdict-file", type=Path, default=DEFAULT_VERDICT)
    parser.add_argument("--work-root", type=Path, default=DEFAULT_WORK_ROOT)
    parser.add_argument("--seeds", default="42", help="Space-separated seeds, e.g. '42' or '42 99 137'")
    parser.add_argument("--quorum", type=int, default=1)
    parser.add_argument("--force-launch", action="store_true")
    parser.add_argument("--detach", action="store_true", default=True)
    parser.add_argument("--no-detach", action="store_true")
    args = parser.parse_args()

    detach = False if args.no_detach else args.detach

    if not RUNNER.exists():
        print(f"ERROR: missing runner: {RUNNER}")
        return 1

    if not args.force_launch:
        state = _job_state(args.phoenix_job_id)
        print(f"Phoenix job {args.phoenix_job_id} state: {state}")

        running_states = {"PENDING", "RUNNING", "CONFIGURING", "COMPLETING", "RESIZING"}
        if state in running_states:
            print("Phoenix still active. Hedge set held.")
            return 0

        verdict = _verdict_text(args.verdict_file)
        if not verdict.strip():
            print(f"Phoenix finished but verdict file not ready: {args.verdict_file}")
            print("No hedge launch yet.")
            return 0

        if _is_pass(verdict):
            print("Phoenix verdict indicates PASS. Hedge set not launched.")
            return 0

        if not _is_ghost_fail(verdict):
            print("Phoenix failed, but not with Ghost-Drug signature; hedge set not auto-launched.")
            print("Use --force-launch to override.")
            return 0

        print("Ghost-Drug signature detected; launching C1-C7 now.")

    args.work_root.mkdir(parents=True, exist_ok=True)
    for job in HEDGE_JOBS:
        _launch_job(job, args.work_root, args.seeds, args.quorum, detach)

    print("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
