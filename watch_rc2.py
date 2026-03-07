#!/usr/bin/env python3
"""
watch_rc2.py — Live monitor for RC2 single-replicate SLURM job
===============================================================
Refreshes every 15 seconds and shows:
  - SLURM job state (RUNNING / PENDING / COMPLETED / FAILED)
  - Progress bar (snapshots out of 169)
  - Current simulation time & phase (barrier / drug / recovery)
  - Live tumor cell count from latest snapshot
  - Completion message when done

Usage:
    python watch_rc2.py          # runs until job finishes
    python watch_rc2.py --once   # single snapshot then exit
"""
from __future__ import annotations

import os
import sys
import time
import subprocess
import xml.etree.ElementTree as ET
from datetime import datetime, timedelta
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# ── Constants (must match run_reality_check_2.py) ─────────────────────────────
PROJECT_ROOT = Path(__file__).resolve().parent
WORK_DIR     = PROJECT_ROOT / "build" / "rc2_full_seed42" / "replicate_01_seed42"
OUTPUT_DIR   = WORK_DIR / "output"
REFRESH      = 15  # seconds between redraws

T_PRE       = 20160.0   # day 14: drug starts
T_TREAT_END = 40320.0   # day 28: drug withdrawn
T_POST      = 60480.0   # day 42: simulation ends
SAVE_INTERVAL = 360     # minutes between snapshots
TOTAL_SNAPS   = int(T_POST / SAVE_INTERVAL) + 1  # 169 (output at t=0..60480)

SEED = 42

TERMINAL_STATES = {
    "COMPLETED", "FAILED", "CANCELLED", "TIMEOUT",
    "OUT_OF_MEMORY", "NODE_FAIL", "PREEMPTED", "BOOT_FAIL", "DEADLINE",
}

# ── ANSI helpers ──────────────────────────────────────────────────────────────
RST     = "\033[0m"
BOLD    = "\033[1m"
DIM     = "\033[2m"
GREEN   = "\033[92m"
YELLOW  = "\033[93m"
BLUE    = "\033[94m"
CYAN    = "\033[96m"
RED     = "\033[91m"
MAGENTA = "\033[95m"
WHITE   = "\033[97m"

def C(text: str, *codes: str) -> str:
    return "".join(codes) + text + RST

def state_fmt(s: str) -> str:
    PAD = 14
    if s == "RUNNING":   return C(f"{s:<{PAD}}", GREEN)
    if s == "PENDING":   return C(f"{s:<{PAD}}", YELLOW)
    if s == "COMPLETED": return C(f"{s:<{PAD}}", BLUE)
    if s in TERMINAL_STATES: return C(f"{s:<{PAD}}", RED)
    return f"{s:<{PAD}}"

def progress_bar(n: int, total: int, width: int = 35) -> str:
    filled = int(width * n / total) if total > 0 else 0
    pct    = 100.0 * n / total if total > 0 else 0.0
    filled_str = C("█" * filled, GREEN if n >= total else CYAN)
    empty_str  = DIM + "░" * (width - filled) + RST
    return f"[{filled_str}{empty_str}] {n:3d}/{total} ({pct:5.1f}%)"

def phase_label(sim_t: float) -> str:
    if sim_t < T_PRE:
        day = sim_t / 1440.0
        return C(f"Phase 0: Barrier forming (day {day:.1f}/14)", DIM)
    elif sim_t < T_TREAT_END:
        day = sim_t / 1440.0
        return C(f"Phase 1: DRUG ACTIVE (day {day:.1f}/28)", GREEN + BOLD)
    elif sim_t <= T_POST:
        day = sim_t / 1440.0
        return C(f"Phase 2: Recovery (day {day:.1f}/42)", YELLOW)
    return ""

# ── Filesystem helpers ────────────────────────────────────────────────────────
def find_job_id() -> Optional[str]:
    files = sorted(WORK_DIR.glob("slurm_*.out"))
    if not files:
        return None
    try:
        return files[-1].stem.split("_")[1]
    except Exception:
        return None

def count_snaps() -> int:
    if not OUTPUT_DIR.exists():
        return 0
    return len(list(OUTPUT_DIR.glob("output*.xml")))

def latest_snap_info() -> Tuple[Optional[float], Optional[int]]:
    """Return (sim_time, n_agents) from most recent snapshot XML."""
    xmls = sorted(OUTPUT_DIR.glob("output*.xml"))
    if not xmls:
        return None, None
    try:
        root = ET.parse(xmls[-1]).getroot()
        ct = root.find(".//metadata/current_time")
        sim_t = float(ct.text) if ct is not None else None
        # try to get cell count from the filename index
        return sim_t, None
    except Exception:
        return None, None

def read_slurm_tail(n: int = 6) -> List[str]:
    """Last n lines of the SLURM stdout log."""
    files = sorted(WORK_DIR.glob("slurm_*.out"))
    if not files:
        return []
    try:
        lines = files[-1].read_text().strip().splitlines()
        return lines[-n:]
    except Exception:
        return []

# ── SLURM helpers ─────────────────────────────────────────────────────────────
def query_state(job_id: str) -> str:
    # squeue first
    try:
        r = subprocess.run(
            ["squeue", "-j", job_id, "--format=%T", "--noheader"],
            capture_output=True, text=True, timeout=10, check=False,
        )
        st = r.stdout.strip()
        if st:
            return st
    except Exception:
        pass
    # sacct fallback
    try:
        r = subprocess.run(
            ["sacct", "-j", job_id, "--format=State", "--noheader", "--parsable2"],
            capture_output=True, text=True, timeout=10, check=False,
        )
        for line in r.stdout.strip().splitlines():
            st = line.strip().split()[0]
            if st:
                return st
    except Exception:
        pass
    return "UNKNOWN"

def parse_agents_from_log() -> Optional[int]:
    """Parse 'total agents: NNN' from SLURM stdout."""
    lines = read_slurm_tail(10)
    for line in reversed(lines):
        if "total agents:" in line:
            try:
                return int(line.split("total agents:")[1].strip())
            except Exception:
                pass
    return None

# ── Event tracking ────────────────────────────────────────────────────────────
def main() -> None:
    once = "--once" in sys.argv
    print("\033[?25l", end="", flush=True)  # hide cursor

    event_log: List[str] = []
    seen_drug_on = False
    seen_drug_off = False
    wall_start = time.time()

    try:
        while True:
            now = datetime.now()
            elapsed = timedelta(seconds=int(time.time() - wall_start))

            jid   = find_job_id()
            state = query_state(jid) if jid else "WAITING"
            n     = count_snaps()
            sim_t, _ = latest_snap_info()
            agents = parse_agents_from_log()

            # detect phase transitions
            ts_str = now.strftime("%H:%M:%S")
            if sim_t is not None:
                if sim_t >= T_PRE and not seen_drug_on:
                    seen_drug_on = True
                    event_log.append(
                        f"{ts_str}  {C('💊 Drug ON', GREEN + BOLD)}"
                        f"  t={sim_t:.0f} min (day 14)"
                    )
                if sim_t >= T_TREAT_END and not seen_drug_off:
                    seen_drug_off = True
                    event_log.append(
                        f"{ts_str}  {C('🛑 Drug WITHDRAWN', YELLOW + BOLD)}"
                        f"  t={sim_t:.0f} min (day 28)"
                    )
            event_log[:] = event_log[-8:]

            # ── render ─────────────────────────────────────────────────────
            print("\033[H\033[J", end="")

            print(
                f"{C('RC2 LIVE MONITOR', BOLD + WHITE)}  "
                f"drug_kill_coefficient = {C('0.0005', CYAN + BOLD)}  "
                f"seed = {C(str(SEED), CYAN)}"
            )
            print(
                f"{now.strftime('%Y-%m-%d %H:%M:%S')}  "
                f"elapsed {elapsed}  "
                f"{C(f'(refresh {REFRESH}s)', DIM)}"
            )
            print("─" * 78)

            # Job info
            print(f"\n  SLURM job:  {jid or '?'}  {state_fmt(state)}")
            print(f"  Node:       cpu384g (32 CPUs, OMP_NUM_THREADS=32)")
            print(f"  Config:     drug_kill_coefficient = 0.0005")
            print()

            # Progress
            print(f"  Snapshots:  {progress_bar(n, TOTAL_SNAPS)}")
            if sim_t is not None:
                print(f"  Sim time:   {sim_t:.0f} / {T_POST:.0f} min  "
                      f"(day {sim_t/1440:.1f} / 42.0)")
                print(f"  Phase:      {phase_label(sim_t)}")
            if agents is not None:
                print(f"  Agents:     {C(str(agents), BOLD)}")

            # Phase timeline
            print()
            print("─" * 78)
            print(f"  {C('PHASE TIMELINE', BOLD)}")
            bar_w = 60
            if sim_t is not None:
                p0_w = int(bar_w * T_PRE / T_POST)
                p1_w = int(bar_w * (T_TREAT_END - T_PRE) / T_POST)
                p2_w = bar_w - p0_w - p1_w
                cursor = int(bar_w * sim_t / T_POST)
                cursor = min(cursor, bar_w - 1)

                bar = ""
                for i in range(bar_w):
                    if i == cursor:
                        bar += C("▶", BOLD + WHITE)
                    elif i < p0_w:
                        bar += C("─", DIM)
                    elif i < p0_w + p1_w:
                        bar += C("█", GREEN)
                    else:
                        bar += C("░", YELLOW)
                print(f"  [{bar}]")
                print(f"  {C('Barrier', DIM):>18s}  │  "
                      f"{C('Drug ON', GREEN):>18s}  │  "
                      f"{C('Recovery', YELLOW)}")
                print(f"  {'d0─────d14':>18s}  │  {'d14────d28':>18s}  │  d28────d42")
            else:
                print(f"  {C('Waiting for simulation data...', DIM)}")

            # Event log
            print()
            print("─" * 78)
            if event_log:
                print(f"  {C('Events', BOLD)}")
                for ev in event_log:
                    print(f"    {ev}")
            else:
                print(f"  {C('No phase transitions yet — drug fires at t=20160 min (day 14)', DIM)}")

            # SLURM log tail
            log_lines = read_slurm_tail(4)
            if log_lines:
                print()
                print(f"  {C('SLURM log (tail)', DIM)}")
                for line in log_lines:
                    print(f"    {C(line.strip(), DIM)}")

            print()

            # ── completion ─────────────────────────────────────────────────
            is_done = state in TERMINAL_STATES
            if is_done or once:
                if is_done:
                    print("=" * 78)
                    if state == "COMPLETED":
                        print(f"  {C('SIMULATION COMPLETE', BOLD + GREEN)}  "
                              f"({n} snapshots, wall {elapsed})")
                        print(f"  Run:  python evaluate_rc2.py")
                    else:
                        print(f"  {C(f'JOB {state}', BOLD + RED)}")
                        print(f"  Check:  cat {WORK_DIR}/slurm_{jid}.err")
                    print("=" * 78)
                break

            time.sleep(REFRESH)

    finally:
        print("\033[?25h", end="", flush=True)  # restore cursor


if __name__ == "__main__":
    main()
