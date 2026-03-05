#!/usr/bin/env python3
"""
watch_rc3.py — Live monitor for RC3 SLURM jobs
================================================
Refreshes every 30 seconds and shows:
  • SLURM job states  (RUNNING / PENDING / COMPLETED / FAILED)
  • Per-replicate progress bars  (snapshots out of 113)
  • Drug / SHH scheduling events as they fire
  • Completion message when all replicates finish

Usage:
    python watch_rc3.py          # runs until all 15 jobs finish
    python watch_rc3.py --once   # single snapshot then exit
"""

from __future__ import annotations

import os
import sys
import time
import subprocess
import xml.etree.ElementTree as ET
from datetime import datetime, timedelta
from pathlib import Path
from typing import Dict, List, Optional

# ── Constants (must match run_reality_check_3.py) ────────────────────────────
PROJECT_ROOT  = Path(__file__).resolve().parent
WORK_DIR      = PROJECT_ROOT / "build" / "reality_check_3"
REFRESH       = 30           # seconds between redraws

T_PRE         = 20160.0      # barrier → intervention boundary (min)
T_END         = 40320.0      # total simulation time (min)
SAVE_INTERVAL = 360          # minutes between saved snapshots
# snapshots: output at t=0, 360, 720, … 40320  →  40320/360 + 1 = 113
TOTAL_SNAPS   = int(T_END / SAVE_INTERVAL) + 1   # 113

# snapshot index where interventions start (t = T_PRE)
INTERV_SNAP   = int(T_PRE / SAVE_INTERVAL)        # 56

ARMS      = ["A", "B", "C"]
ARM_NAMES = {
    "A": "Control",
    "B": "SHH inhibition only",
    "C": "SHH inhibition + drug",
}
SEEDS  = [42, 99, 137, 256, 1001]
N_REPS = len(SEEDS)

TERMINAL_STATES = {
    "COMPLETED", "FAILED", "CANCELLED", "TIMEOUT",
    "OUT_OF_MEMORY", "NODE_FAIL", "PREEMPTED", "BOOT_FAIL", "DEADLINE",
}

# ── ANSI helpers ─────────────────────────────────────────────────────────────
RST    = "\033[0m"
BOLD   = "\033[1m"
DIM    = "\033[2m"
GREEN  = "\033[92m"
YELLOW = "\033[93m"
BLUE   = "\033[94m"
CYAN   = "\033[96m"
RED    = "\033[91m"
MAGENTA= "\033[95m"

def C(text: str, *codes: str) -> str:
    return "".join(codes) + text + RST

def state_fmt(s: str) -> str:
    PAD = 14
    if s == "RUNNING":    return C(f"{s:<{PAD}}", GREEN)
    if s == "PENDING":    return C(f"{s:<{PAD}}", YELLOW)
    if s == "COMPLETED":  return C(f"{s:<{PAD}}", BLUE)
    if s in TERMINAL_STATES: return C(f"{s:<{PAD}}", RED)
    return f"{s:<{PAD}}"

def progress_bar(n: int, total: int, width: int = 25) -> str:
    filled  = int(width * n / total) if total > 0 else 0
    pct     = 100.0 * n / total if total > 0 else 0.0
    filled_str  = C("█" * filled, GREEN if n >= total else CYAN)
    empty_str   = DIM + "░" * (width - filled) + RST
    return f"[{filled_str}{empty_str}] {n:3d}/{total} ({pct:5.1f}%)"

# ── Filesystem helpers ────────────────────────────────────────────────────────
def run_dir(arm: str, rep_idx: int, seed: int) -> Path:
    return WORK_DIR / f"arm_{arm}" / f"replicate_{rep_idx+1:02d}_seed{seed}"

def find_job_id(rd: Path) -> Optional[str]:
    """Recover job ID from the slurm_<JOBID>.out file SLURM creates."""
    files = list(rd.glob("slurm_*.out"))
    if not files:
        return None
    try:
        return files[0].stem.split("_")[1]
    except Exception:
        return None

def count_snaps(output_dir: Path) -> int:
    if not output_dir.exists():
        return 0
    return len(list(output_dir.glob("output*.xml")))

def latest_sim_time(output_dir: Path) -> Optional[float]:
    """Parse <current_time> from the most recent PhysiCell output XML."""
    xmls = sorted(output_dir.glob("output*.xml"))
    if not xmls:
        return None
    try:
        root = ET.parse(xmls[-1]).getroot()
        ct = root.find(".//metadata/current_time")
        return float(ct.text) if ct is not None else None
    except Exception:
        return None

# ── SLURM helpers ─────────────────────────────────────────────────────────────
def query_states(job_ids: List[str]) -> Dict[str, str]:
    if not job_ids:
        return {}
    states: Dict[str, str] = {}

    # squeue covers running / pending
    try:
        r = subprocess.run(
            ["squeue", "-j", ",".join(job_ids), "--format=%i|%T", "--noheader"],
            capture_output=True, text=True, timeout=15, check=False,
        )
        for line in r.stdout.strip().splitlines():
            parts = line.strip().split("|")
            if len(parts) >= 2:
                states[parts[0].strip()] = parts[1].strip()
    except Exception:
        pass

    # sacct covers already-finished jobs
    missing = [j for j in job_ids if j not in states]
    if missing:
        try:
            r = subprocess.run(
                ["sacct", "-j", ",".join(missing),
                 "--format=JobID,State", "--noheader", "--parsable2"],
                capture_output=True, text=True, timeout=15, check=False,
            )
            for line in r.stdout.strip().splitlines():
                parts = line.strip().split("|")
                if len(parts) >= 2 and "." not in parts[0]:
                    jid = parts[0].strip()
                    st  = parts[1].strip().split()[0]
                    if jid not in states:
                        states[jid] = st
        except Exception:
            pass

    return states

# ── Main ─────────────────────────────────────────────────────────────────────
def main() -> None:
    once = "--once" in sys.argv

    # hide cursor for clean redraws
    print("\033[?25l", end="", flush=True)

    # event tracking (keyed by "arm_repN")
    seen_interv_on:  Dict[str, bool] = {}
    seen_drug_off:   Dict[str, bool] = {}
    event_log:       List[str]       = []
    wall_start = time.time()

    try:
        _loop(once, wall_start, seen_interv_on, seen_drug_off, event_log)
    finally:
        print("\033[?25h", end="", flush=True)  # restore cursor

def _loop(
    once: bool,
    wall_start: float,
    seen_interv_on: Dict[str, bool],
    seen_drug_off: Dict[str, bool],
    event_log: List[str],
) -> None:

    while True:
        now = datetime.now()

        # ── gather data ───────────────────────────────────────────────────────
        runs = []
        for arm in ARMS:
            for i, seed in enumerate(SEEDS):
                rd      = run_dir(arm, i, seed)
                jid     = find_job_id(rd) if rd.exists() else None
                out_dir = rd / "output"
                n       = count_snaps(out_dir)
                sim_t   = latest_sim_time(out_dir)
                runs.append(dict(arm=arm, rep=i+1, seed=seed,
                                 rd=rd, jid=jid, out_dir=out_dir,
                                 n=n, sim_t=sim_t))

        job_ids = [r["jid"] for r in runs if r["jid"]]
        states  = query_states(job_ids)

        # ── detect scheduling events ──────────────────────────────────────────
        for r in runs:
            key    = f"{r['arm']}_r{r['rep']}"
            ts_str = now.strftime("%H:%M:%S")
            sim_str = f" (sim≈{r['sim_t']:.0f} min)" if r["sim_t"] is not None else ""

            # SHH inhibition / drug ON  (n crosses INTERV_SNAP)
            if r["n"] > INTERV_SNAP and not seen_interv_on.get(key):
                seen_interv_on[key] = True
                if r["arm"] == "B":
                    event_log.append(
                        f"{ts_str}  {C('🚫 SHH inhibition ON', MAGENTA)}"
                        f"  Arm B rep {r['rep']} seed={r['seed']}{sim_str}"
                    )
                elif r["arm"] == "C":
                    event_log.append(
                        f"{ts_str}  {C('💊 Drug ON', GREEN)}"
                        f"  Arm C rep {r['rep']} seed={r['seed']}{sim_str}"
                    )

            # Drug WITHDRAWN  (run finishes = drug ends at T_END)
            st = states.get(r["jid"] or "", "")
            if r["arm"] == "C" and st in TERMINAL_STATES and not seen_drug_off.get(key):
                seen_drug_off[key] = True
                event_log.append(
                    f"{ts_str}  {C('🛑 Drug WITHDRAWN', YELLOW)}"
                    f"  Arm C rep {r['rep']} seed={r['seed']}"
                )

        event_log[:] = event_log[-12:]  # rolling window

        # ── is_done predicate ─────────────────────────────────────────────────
        def is_done(r: dict) -> bool:
            st = states.get(r["jid"] or "", "")
            return st in TERMINAL_STATES or r["n"] >= TOTAL_SNAPS

        n_done   = sum(1 for r in runs if is_done(r))
        all_done = n_done == len(runs)

        # ── render ────────────────────────────────────────────────────────────
        # move cursor to top-left without clearing (avoids flicker)
        print("\033[H\033[J", end="")

        elapsed = timedelta(seconds=int(time.time() - wall_start))
        print(
            f"{C('RC3 LIVE MONITOR', BOLD)}  "
            f"{now.strftime('%Y-%m-%d %H:%M:%S')}  "
            f"elapsed {elapsed}  "
            f"{C(f'(refresh {REFRESH}s)', DIM)}"
        )
        print("─" * 78)

        for arm in ARMS:
            arm_runs  = [r for r in runs if r["arm"] == arm]
            n_arm_done = sum(1 for r in arm_runs if is_done(r))
            print(
                f"\n  {C(f'Arm {arm}', BOLD)}  {ARM_NAMES[arm]}"
                f"  {C(f'{n_arm_done}/{N_REPS} complete', DIM)}"
            )
            for r in arm_runs:
                jid = r["jid"] or "?"
                st  = states.get(jid, "PENDING" if jid != "?" else "WAITING")
                pb  = progress_bar(r["n"], TOTAL_SNAPS)
                sf  = state_fmt(st)

                # inline tag for active intervention
                key = f"{r['arm']}_r{r['rep']}"
                tag = ""
                if r["arm"] in ("B", "C") and seen_interv_on.get(key):
                    if r["arm"] == "B":
                        tag = f"  {C('▶ SHH-inhib', MAGENTA)}"
                    elif not seen_drug_off.get(key):
                        tag = f"  {C('▶ Drug active', GREEN)}"
                    else:
                        tag = f"  {C('■ Withdrawn', DIM)}"

                sim_t_val = r["sim_t"]
                sim_str = (
                    f"  {C(f't={sim_t_val:.0f}min', DIM)}"
                    if sim_t_val is not None else ""
                )
                print(f"    rep{r['rep']} seed={r['seed']:4d}  {sf}  {pb}{sim_str}{tag}")

        print()
        print("─" * 78)

        # overall progress
        print(f"\n  Overall  {progress_bar(n_done, len(runs), width=30)}\n")

        # event log
        if event_log:
            print(f"  {C('Scheduling events', BOLD)}")
            for ev in event_log:
                print(f"    {ev}")
        else:
            print(f"  {C('No scheduling events yet — interventions fire at T=20160 min', DIM)}")

        print()

        # ── completion ────────────────────────────────────────────────────────
        if all_done or once:
            if all_done:
                # count successes
                n_ok  = sum(1 for r in runs if states.get(r["jid"] or "", "") == "COMPLETED")
                n_bad = len(runs) - n_ok
                print("=" * 78)
                print(f"  {C('✔  All 15 replicates finished!', BOLD + GREEN)}")
                print(f"  {n_ok}/15 COMPLETED  |  {n_bad}/15 non-zero exit")
                print(f"  Results are printed at the end of the run_reality_check_3.py output.")
                print("=" * 78)
            break

        time.sleep(REFRESH)


if __name__ == "__main__":
    main()
