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

Expected wall time: 2–5 minutes on 16 cores.
"""
from __future__ import annotations

import math
import os
import shutil
import subprocess
import sys
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
OMP_THREADS    = min(os.cpu_count() or 4, 16)

BINARY         = PROJECT_ROOT / "stroma_world"
BASE_CONFIG    = PROJECT_ROOT / "config" / "PhysiCell_settings.xml"
WORK_DIR       = PROJECT_ROOT / "build" / "rc2_smoke"


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
    _set(".//parallel/omp_num_threads", str(OMP_THREADS))

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
    _enforce_dirichlet("ecm_density", "0")

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
    print("  RC2 SMOKE TEST")
    print(f"  Phases: establish 0–{T_ESTABLISH:.0f} | drug {T_TREAT_START:.0f}–{T_TREAT_END:.0f} | washout {T_TREAT_END:.0f}–{T_FINAL:.0f}")
    print(f"  Seed={SEED}  OMP={OMP_THREADS}  save_interval={SAVE_INTERVAL}")
    print("=" * 60)

    # Run simulation
    t0 = time.perf_counter()
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(OMP_THREADS)
    env["OMP_PROC_BIND"] = "spread"
    env["OMP_PLACES"] = "cores"

    proc = subprocess.run(
        [str(BINARY), str(config_path)],
        cwd=str(PROJECT_ROOT),
        env=env,
        capture_output=True,
        text=True,
        timeout=1200,   # 20 min hard cap
    )
    wall = time.perf_counter() - t0

    # Check for stderr DRUG_SCHEDULE messages
    stderr_lines = proc.stderr.strip().split("\n") if proc.stderr else []
    drug_on_msg  = [l for l in stderr_lines if "Drug ON" in l]
    drug_off_msg = [l for l in stderr_lines if "Drug WITHDRAWN" in l]

    print(f"\n  Simulation finished in {wall:.1f}s  exit_code={proc.returncode}")
    print(f"  Snapshots: {len(list(output_dir.glob('output*.xml')))}")
    if drug_on_msg:
        print(f"  {drug_on_msg[0].strip()}")
    if drug_off_msg:
        print(f"  {drug_off_msg[0].strip()}")

    if proc.returncode != 0:
        print(f"\n  *** SMOKE FAIL: simulation crashed ***")
        print(f"  Last 20 stderr lines:")
        for l in stderr_lines[-20:]:
            print(f"    {l}")
        return 1

    # Read checkpoints
    snap_mid_treat = _read_snapshot(output_dir, (T_TREAT_START + T_TREAT_END) / 2.0)   # mid-treatment
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
        # Dump last stderr for debug
        if stderr_lines:
            print("\n  Last 30 stderr lines:")
            for l in stderr_lines[-30:]:
                print(f"    {l}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
