#!/usr/bin/env python3
"""
RC2 Smoke Test — fast sanity check before committing to 6-hour HPC runs.
=========================================================================

Runs 1 replicate with compressed timescale:
  Phase 0:  3 days establish  (0     → 4 320 min)
  Phase 1:  3 days treatment  (4 320 → 8 640 min)   drug BC = 1.0
  Phase 2:  3 days washout    (8 640 → 12 960 min)   drug BC → 0

Checks (all must pass):
  S1  Drug field activates:    max(drug) at t=6480 > 0.1
  S2  Drug field deactivates:  max(drug_bc) at Dirichlet → 0 after drug_end_time
                               AND interior drug at t=12960 < drug at t=8640
  S3  Tumor changes:           |n_tumor(t=8640) - n_tumor(t=4320)| > 0
                               (any directional response, even weak)
  S4  No NaN / negative:       all cell coords finite, no negative cell counts

Expected wall time: 2–3 minutes on 128-core compute node via SLURM.
"""
from __future__ import annotations

import math
import os
import shutil
import subprocess
import sys
import textwrap
import time
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from python.wrapper.output_parser import OutputParser

# ── Timing ────────────────────────────────────────────────────────────────
T_ESTABLISH    = 4320.0    # 3 days
T_TREAT_START  = T_ESTABLISH
T_TREAT_END    = 8640.0    # 3 days drug
T_FINAL        = 12960.0   # 3 days washout
SAVE_INTERVAL  = 360       # every 6 hours → 36 snapshots

SEED           = 42

# SLURM settings
SLURM_PARTITION = "compute"
SLURM_CPUS      = 128
SLURM_MEM       = "128G"
SLURM_TIME      = "00:30:00"     # 30 min wall (generous for a 9-day sim)
SLURM_POLL      = 10             # seconds between polls

BINARY         = PROJECT_ROOT / "stroma_world"
BASE_CONFIG    = PROJECT_ROOT / "config" / "PhysiCell_settings.xml"
WORK_DIR       = PROJECT_ROOT / "build" / "rc2_smoke"

TERMINAL_STATES = {
    "COMPLETED", "FAILED", "CANCELLED", "TIMEOUT",
    "OUT_OF_MEMORY", "NODE_FAIL", "PREEMPTED",
    "BOOT_FAIL", "DEADLINE",
}


def _patch_config(src: Path, dst: Path, output_dir: Path) -> None:
    tree = ET.parse(src)
    root = tree.getroot()

    def _set(xpath: str, value: str) -> None:
        node = root.find(xpath)
        if node is not None:
            node.text = value

    _set("./save/folder", str(output_dir))
    _set(".//overall/max_time", str(T_FINAL))
    _set(".//options/random_seed", str(SEED))
    _set(".//parallel/omp_num_threads", str(SLURM_CPUS))

    for n in root.findall(".//save//interval"):
        n.text = str(SAVE_INTERVAL)

    _set(".//user_parameters/drug_start_time",   str(T_TREAT_START))
    _set(".//user_parameters/drug_end_time",     str(T_TREAT_END))
    _set(".//user_parameters/drug_concentration", "1.0")

    # Dirichlet BCs
    def _enforce_dirichlet(var_name: str, value: str) -> None:
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

    _enforce_dirichlet("oxygen",      "38")
    _enforce_dirichlet("tgfb",        "0")
    _enforce_dirichlet("shh",         "0")
    _enforce_dirichlet("drug",        "0")

    tree.write(dst, encoding="utf-8", xml_declaration=True)


def _find_nearest_snapshot(output_dir: Path, target_time: float) -> Optional[Path]:
    xmls = sorted(output_dir.glob("output*.xml"))
    if not xmls:
        return None
    best, best_dt = None, float("inf")
    for xml in xmls:
        try:
            tree = ET.parse(xml)
            t = float(tree.getroot().find(".//current_time").text)
            dt = abs(t - target_time)
            if dt < best_dt:
                best_dt = dt
                best = xml
        except Exception:
            continue
    return best


def _read_snapshot(output_dir: Path, target_time: float):
    """Return (time, n_tumor, max_drug_field, positions_ok, micro_drug_array)."""
    parser = OutputParser(output_dir)
    xmls = sorted(output_dir.glob("output*.xml"))
    if not xmls:
        return None

    # Find nearest
    best_xml, best_dt = None, float("inf")
    for xml in xmls:
        try:
            snap = parser._read_physicell_xml(xml)
            t = float(snap["time"])
            dt = abs(t - target_time)
            if dt < best_dt:
                best_dt = dt
                best_xml = xml
                best_snap = snap
        except Exception:
            continue
    if best_xml is None:
        return None

    snap = best_snap
    matrix = snap["cell_matrix"]
    labels = snap["label_name_map"]
    micro_values = snap["micro_values"]

    def _row(name):
        e = labels.get(name)
        if e is None:
            return None
        idx = int(e["index"])
        if idx < 0 or idx >= matrix.shape[0]:
            return None
        return matrix[idx, :]

    # Cell type + live mask
    cell_type = _row("cell_type")
    dead = _row("dead")
    n_cells = matrix.shape[1]

    live_mask = np.ones(n_cells, dtype=bool)
    if dead is not None:
        live_mask &= dead <= 0.5

    ctype_int = np.rint(cell_type).astype(int) if cell_type is not None else np.full(n_cells, -1, dtype=int)
    n_tumor = int(np.sum(live_mask & (ctype_int == 0)))

    # Drug field
    drug = micro_values.get("drug")
    max_drug = float(np.nanmax(drug)) if drug is not None and drug.size > 0 else 0.0
    mean_drug = float(np.nanmean(drug)) if drug is not None and drug.size > 0 else 0.0

    # Position sanity
    pos = parser._get_positions(matrix, labels)
    positions_ok = True
    if pos.size > 0:
        if np.any(~np.isfinite(pos)):
            positions_ok = False

    return {
        "time": float(snap["time"]),
        "n_tumor": n_tumor,
        "max_drug": max_drug,
        "mean_drug": mean_drug,
        "positions_ok": positions_ok,
        "n_cells": n_cells,
    }


def _query_job(job_id: str):
    """Return (state, exit_code, is_terminal)."""
    try:
        r = subprocess.run(
            ["sacct", "-j", job_id, "--format=State,ExitCode",
             "--noheader", "--parsable2"],
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


def main() -> int:
    if not BINARY.exists():
        print(f"ERROR: binary not found: {BINARY}")
        return 1

    # Clean workspace
    if WORK_DIR.exists():
        shutil.rmtree(WORK_DIR)
    WORK_DIR.mkdir(parents=True)
    output_dir = WORK_DIR / "output"
    output_dir.mkdir()

    config_path = WORK_DIR / "config.xml"
    _patch_config(BASE_CONFIG, config_path, output_dir)

    # Copy supporting config files
    local_cfg = WORK_DIR / "config"
    local_cfg.mkdir(exist_ok=True)
    for f in ["tumor_calibration_knobs.json", "gene_params_default.json"]:
        src = PROJECT_ROOT / "config" / f
        if src.exists():
            shutil.copy(src, local_cfg / f)

    print("=" * 60)
    print("  RC2 SMOKE TEST  (SLURM)")
    print(f"  Phases: establish 0–{T_ESTABLISH:.0f} | drug {T_TREAT_START:.0f}–{T_TREAT_END:.0f} | washout {T_TREAT_END:.0f}–{T_FINAL:.0f}")
    print(f"  Seed={SEED}  CPUs={SLURM_CPUS}  save_interval={SAVE_INTERVAL}")
    print("=" * 60)

    # Write SLURM script
    slurm_script = WORK_DIR / "smoke.slurm.sh"
    slurm_script.write_text(textwrap.dedent(f"""\
        #!/bin/bash
        #SBATCH --job-name=rc2_smoke
        #SBATCH --partition={SLURM_PARTITION}
        #SBATCH --nodes=1
        #SBATCH --ntasks-per-node=1
        #SBATCH --cpus-per-task={SLURM_CPUS}
        #SBATCH --mem={SLURM_MEM}
        #SBATCH --time={SLURM_TIME}
        #SBATCH --output={WORK_DIR}/slurm_%j.out
        #SBATCH --error={WORK_DIR}/slurm_%j.err

        set -euo pipefail
        export OMP_NUM_THREADS={SLURM_CPUS}
        export OMP_PROC_BIND=spread
        export OMP_PLACES=cores

        cd {PROJECT_ROOT}
        echo "=== RC2 SMOKE TEST  seed={SEED} ==="
        echo "Start: $(date)"
        {BINARY} {config_path}
        echo "End:   $(date)"
    """), encoding="utf-8")
    slurm_script.chmod(0o755)

    # Submit
    result = subprocess.run(
        ["sbatch", "--parsable", str(slurm_script)],
        capture_output=True, text=True, check=False,
    )
    if result.returncode != 0:
        print(f"  ERROR: sbatch failed: {result.stderr.strip()}")
        return 1

    job_id = result.stdout.strip().split(";")[0]
    print(f"  Submitted SLURM job {job_id}")

    # Poll until done
    t0 = time.perf_counter()
    while True:
        time.sleep(SLURM_POLL)
        state, exit_code, terminal = _query_job(job_id)
        elapsed = time.perf_counter() - t0
        n_snaps = len(list(output_dir.glob("output*.xml")))
        print(f"    [{elapsed:5.0f}s] {state}  snapshots={n_snaps}/~36")
        if terminal:
            break

    wall = time.perf_counter() - t0
    print(f"\n  Job {job_id} finished: {state}  exit={exit_code}  wall={wall:.1f}s")

    if state != "COMPLETED" or exit_code != 0:
        print(f"\n  *** SMOKE FAIL: SLURM job {state} exit={exit_code} ***")
        # Dump stderr file
        for errfile in sorted(WORK_DIR.glob("slurm_*.err")):
            print(f"\n  {errfile.name} (last 30 lines):")
            lines = errfile.read_text().strip().splitlines()
            for l in lines[-30:]:
                print(f"    {l}")
        return 1

    # Read DRUG_SCHEDULE messages from stderr
    drug_on_msg = []
    drug_off_msg = []
    for errfile in sorted(WORK_DIR.glob("slurm_*.err")):
        for l in errfile.read_text().splitlines():
            if "Drug ON" in l:
                drug_on_msg.append(l)
            if "Drug WITHDRAWN" in l:
                drug_off_msg.append(l)

    n_snaps = len(list(output_dir.glob("output*.xml")))
    print(f"  Snapshots: {n_snaps}")
    if drug_on_msg:
        print(f"  {drug_on_msg[0].strip()}")
    if drug_off_msg:
        print(f"  {drug_off_msg[0].strip()}")

    # Read checkpoints
    snap_mid_treat = _read_snapshot(output_dir, (T_TREAT_START + T_TREAT_END) / 2.0)
    snap_pre       = _read_snapshot(output_dir, T_TREAT_START)
    snap_treat_end = _read_snapshot(output_dir, T_TREAT_END)
    snap_final     = _read_snapshot(output_dir, T_FINAL)

    for label, s in [("pre", snap_pre), ("mid-treat", snap_mid_treat),
                     ("treat_end", snap_treat_end), ("final", snap_final)]:
        if s is None:
            print(f"  WARNING: could not read {label} snapshot")
        else:
            print(f"  {label:10s}  t={s['time']:8.0f}  tumor={s['n_tumor']:5d}  "
                  f"max_drug={s['max_drug']:.4f}  mean_drug={s['mean_drug']:.4f}  "
                  f"positions_ok={s['positions_ok']}")

    # ── S1: Drug field activates ─────────────────────────────────────────
    s1 = False
    s1_detail = ""
    if snap_mid_treat is not None:
        s1 = snap_mid_treat["max_drug"] > 0.1
        s1_detail = f"max_drug(mid-treat)={snap_mid_treat['max_drug']:.4f}"
    else:
        s1_detail = "no mid-treatment snapshot"

    # ── S2: Drug field deactivates (natural dissipation) ─────────────────
    s2 = False
    s2_detail = ""
    if snap_treat_end is not None and snap_final is not None:
        # Interior drug should be lower at t_final than at t_treat_end
        # (natural diffusion + decay/uptake clears it over 3 days)
        s2 = snap_final["mean_drug"] < snap_treat_end["mean_drug"]
        s2_detail = (f"mean_drug(treat_end)={snap_treat_end['mean_drug']:.4f} → "
                     f"mean_drug(final)={snap_final['mean_drug']:.4f}")
        # Also confirm drug_off message
        if not drug_off_msg:
            s2 = False
            s2_detail += "  (NO Drug WITHDRAWN log message!)"
    else:
        s2_detail = "missing snapshots"

    # ── S3: Tumor responds ───────────────────────────────────────────────
    s3 = False
    s3_detail = ""
    if snap_pre is not None and snap_treat_end is not None:
        delta = abs(snap_treat_end["n_tumor"] - snap_pre["n_tumor"])
        s3 = delta > 0
        s3_detail = (f"tumor(pre={snap_pre['n_tumor']}) → "
                     f"tumor(treat_end={snap_treat_end['n_tumor']})  |Δ|={delta}")
    else:
        s3_detail = "missing snapshots"

    # ── S4: No NaN / negative ────────────────────────────────────────────
    s4 = True
    s4_detail = "OK"
    for label, s in [("pre", snap_pre), ("treat_end", snap_treat_end), ("final", snap_final)]:
        if s is not None:
            if not s["positions_ok"]:
                s4 = False
                s4_detail = f"NaN positions at {label}"
            if s["n_tumor"] < 0:
                s4 = False
                s4_detail = f"negative tumor count at {label}"

    # ── Report ────────────────────────────────────────────────────────────
    checks = [
        ("S1  Drug field activates",    s1, s1_detail),
        ("S2  Drug field deactivates",  s2, s2_detail),
        ("S3  Tumor responds",          s3, s3_detail),
        ("S4  No NaN / negative",       s4, s4_detail),
    ]

    print()
    print("=" * 60)
    print("  SMOKE TEST RESULTS")
    print("=" * 60)
    all_ok = True
    for label, ok, detail in checks:
        mark = "PASS" if ok else "FAIL"
        if not ok:
            all_ok = False
        print(f"  [{mark}] {label}")
        print(f"         {detail}")

    print()
    print(f"  Wall time: {wall:.1f}s  ({wall/60:.1f} min)")

    if all_ok:
        print("\n  *** RC2 SMOKE TEST: PASS ***")
        return 0
    else:
        print("\n  *** RC2 SMOKE TEST: FAIL ***")
        # Dump SLURM stderr for debug
        for errfile in sorted(WORK_DIR.glob("slurm_*.err")):
            lines = errfile.read_text().strip().splitlines()
            if lines:
                print(f"\n  {errfile.name} (last 30 lines):")
                for l in lines[-30:]:
                    print(f"    {l}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
