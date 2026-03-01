#!/usr/bin/env python3
"""
Evaluate Reality Check 1 criteria on completed simulation output.

Reads the snapshots directly (avoiding the slow O(n²) _max_pairwise_distance
in OutputParser.parse_final_state) and evaluates all 8 biological criteria.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(PROJECT_ROOT))

from python.wrapper.output_parser import OutputParser

SEEDS = [42, 99, 137, 256, 1001]
QUORUM = 4
BOUNDARY_O2 = 38.0

CRITERIA_NAMES = [
    ("barrier_self_assembles",    "1. Barrier self-assembles"),
    ("spatial_ecm_gradient",      "2. Spatial ECM gradient"),
    ("central_hypoxia",           "3. Central hypoxia develops"),
    ("emt_peripheral_not_global", "4. EMT peripheral, not global"),
    ("pdac_like_ratio",           "5. Tumor-stroma ratio PDAC-like"),
    ("no_numerical_instability",  "6. No numerical instability"),
    ("caf_activation",            "7. CAF count increases"),
    ("tgfb_gradient",             "8. TGF-b gradient near tumor"),
]

def _distance_to_boundary(coords, all_coords):
    mins = np.min(all_coords, axis=0)
    maxs = np.max(all_coords, axis=0)
    d_lo = coords - mins[None, :]
    d_hi = maxs[None, :] - coords
    d = np.minimum(d_lo, d_hi)
    return np.min(d, axis=1)

def evaluate_replicate(rep_dir: Path, seed: int):
    """Evaluate 8 criteria on a single replicate. Returns (criteria, details)."""
    output_dir = rep_dir / "output"
    parser = OutputParser(output_dir)

    print(f"    Parsing final snapshot for seed={seed}...", flush=True)
    final_xml = parser._find_final_snapshot_xml()
    snapshot = parser._read_physicell_xml(final_xml)

    matrix = snapshot["cell_matrix"]
    labels = snapshot["label_name_map"]
    micro_coords = snapshot["micro_coords"]
    micro_values = snapshot["micro_values"]
    sim_time = snapshot["time"]

    def _row(name):
        entry = labels.get(name)
        if entry is None:
            return None
        idx = int(entry["index"])
        if idx < 0 or idx >= matrix.shape[0]:
            return None
        return matrix[idx, :]

    cell_type = _row("cell_type")
    dead = _row("dead")
    death_model = _row("current_death_model")
    is_mes = _row("is_mesenchymal")
    is_activated = _row("is_activated")
    if is_activated is None:
        is_activated = _row("ACTA2")
    positions = parser._get_positions(matrix, labels)

    cell_type_int = np.rint(cell_type).astype(int) if cell_type is not None else np.zeros(matrix.shape[1], dtype=int)
    tumor_mask = cell_type_int == 0
    stroma_mask = cell_type_int == 1

    dead_mask = (dead > 0.5) if dead is not None else np.zeros_like(tumor_mask, dtype=bool)
    apoptotic_mask = (np.rint(death_model).astype(int) == 100) if death_model is not None else np.zeros_like(tumor_mask, dtype=bool)
    live_mask = ~(dead_mask | apoptotic_mask)

    n_tumor = int(np.sum(tumor_mask))
    n_stroma = int(np.sum(stroma_mask))
    n_activated = int(np.sum(stroma_mask & (is_activated > 0.5))) if is_activated is not None else 0

    tumor_pos = positions[tumor_mask] if positions.size else np.empty((0, 3))
    centroid = np.nanmean(tumor_pos, axis=0) if tumor_pos.size else np.zeros(3)
    tumor_radius = float(np.max(np.linalg.norm(tumor_pos - centroid, axis=1))) if tumor_pos.size else 0.0

    o2 = micro_values.get("oxygen")
    ecm = micro_values.get("ecm_density")
    tgfb = micro_values.get("tgfb")

    criteria = {}
    details = {}

    print(f"    t={sim_time:.0f}min  tumor={n_tumor}  stroma={n_stroma}  "
          f"activated_cafs={n_activated}  tumor_radius={tumor_radius:.0f}um", flush=True)

    # --- C1: Barrier self-assembles ---
    peri_ecm = math.nan
    if ecm is not None and micro_coords.size and tumor_pos.size:
        vd = np.linalg.norm(micro_coords - centroid, axis=1)
        shell = (vd >= tumor_radius) & (vd <= tumor_radius + 100.0)
        if shell.any():
            peri_ecm = float(np.nanmean(ecm[shell]))
    c1 = math.isfinite(peri_ecm) and peri_ecm > 0.5
    criteria["barrier_self_assembles"] = c1
    details["barrier_self_assembles"] = f"peritumoral_ecm={peri_ecm:.4f} {'>' if c1 else '<='} 0.5"

    # --- C2: ECM gradient ---
    ecm_near = ecm_edge = math.nan
    if ecm is not None and micro_coords.size and tumor_pos.size:
        vd = np.linalg.norm(micro_coords - centroid, axis=1)
        near = vd <= tumor_radius + 150.0
        bd = _distance_to_boundary(micro_coords, micro_coords)
        edge = bd <= 50.0
        if near.any(): ecm_near = float(np.nanmean(ecm[near]))
        if edge.any(): ecm_edge = float(np.nanmean(ecm[edge]))
    c2 = math.isfinite(ecm_near) and math.isfinite(ecm_edge) and ecm_near > ecm_edge
    criteria["spatial_ecm_gradient"] = c2
    details["spatial_ecm_gradient"] = f"ECM_near={ecm_near:.4f} {'>' if c2 else '<='} ECM_edge={ecm_edge:.4f}"

    # --- C3: Central hypoxia ---
    o2_centroid = math.nan
    o2_boundary = BOUNDARY_O2
    if o2 is not None and micro_coords.size:
        dists = np.linalg.norm(micro_coords - centroid, axis=1)
        nearest_voxel = np.argmin(dists)
        o2_centroid = float(o2[nearest_voxel])
        bd = _distance_to_boundary(micro_coords, micro_coords)
        edge = bd <= 30.0
        if edge.any():
            o2_boundary = float(np.nanmean(o2[edge]))
    threshold = 0.5 * o2_boundary
    c3 = math.isfinite(o2_centroid) and o2_centroid < threshold
    criteria["central_hypoxia"] = c3
    details["central_hypoxia"] = f"O2_centroid={o2_centroid:.2f} {'<' if c3 else '>='} 0.5*O2_boundary={threshold:.2f}"

    # --- C4: EMT peripheral > core ---
    emt_core = emt_periph = math.nan
    if is_mes is not None and tumor_pos.size >= 10:
        tm = is_mes[tumor_mask] > 0.5
        radial = np.linalg.norm(tumor_pos - centroid, axis=1)
        core_cut = float(np.quantile(radial, 0.40))
        periph_cut = float(np.quantile(radial, 0.60))
        core_m = radial <= core_cut
        periph_m = radial >= periph_cut
        if core_m.any(): emt_core = float(np.mean(tm[core_m]))
        if periph_m.any(): emt_periph = float(np.mean(tm[periph_m]))
    c4 = math.isfinite(emt_core) and math.isfinite(emt_periph) and emt_periph > emt_core
    criteria["emt_peripheral_not_global"] = c4
    details["emt_peripheral_not_global"] = f"EMT_periph={emt_periph:.4f} {'>' if c4 else '<='} EMT_core={emt_core:.4f}"

    # --- C5: Tumor fraction < 50% ---
    total = n_tumor + n_stroma
    frac = n_tumor / max(total, 1)
    c5 = frac < 0.50
    criteria["pdac_like_ratio"] = c5
    details["pdac_like_ratio"] = f"tumor_frac={frac:.4f} ({'<' if c5 else '>='} 0.50) tumor={n_tumor} stroma={n_stroma}"

    # --- C6: No numerical instability ---
    issues = []
    if ecm is not None:
        if np.any(~np.isfinite(ecm)): issues.append("ECM has NaN/Inf")
        if np.any(ecm > 1.0 + 1e-9): issues.append(f"ECM max={np.max(ecm):.6f}>1.0")
    for fname, fvals in micro_values.items():
        if fvals is not None:
            if np.any(~np.isfinite(fvals)): issues.append(f"{fname} has NaN/Inf")
            if np.any(fvals < -1e-12): issues.append(f"{fname} min={np.min(fvals):.6e}<0")
    if o2 is not None and np.any(o2 > BOUNDARY_O2 + 1e-6):
        issues.append(f"O2 max={np.max(o2):.4f}>{BOUNDARY_O2}")
    c6 = len(issues) == 0
    criteria["no_numerical_instability"] = c6
    details["no_numerical_instability"] = "clean" if c6 else f"ISSUES: {'; '.join(issues)}"

    # --- C7: CAF activation ---
    c7 = n_activated > 0
    criteria["caf_activation"] = c7
    details["caf_activation"] = f"activated_cafs={n_activated}"

    # --- C8: TGF-b gradient ---
    tgfb_near = tgfb_edge = math.nan
    if tgfb is not None and micro_coords.size and tumor_pos.size:
        vd = np.linalg.norm(micro_coords - centroid, axis=1)
        near = vd <= tumor_radius + 150.0
        bd = _distance_to_boundary(micro_coords, micro_coords)
        edge = bd <= 50.0
        if near.any(): tgfb_near = float(np.nanmean(tgfb[near]))
        if edge.any(): tgfb_edge = float(np.nanmean(tgfb[edge]))
    c8 = math.isfinite(tgfb_near) and math.isfinite(tgfb_edge) and tgfb_near > tgfb_edge
    criteria["tgfb_gradient"] = c8
    details["tgfb_gradient"] = f"TGFb_near={tgfb_near:.6f} {'>' if c8 else '<='} TGFb_edge={tgfb_edge:.6f}"

    return criteria, details


def main():
    work_dir = PROJECT_ROOT / "build" / "reality_check_1"

    print("=" * 72)
    print("  REALITY CHECK 1 — Evaluation of completed HPC runs")
    print("=" * 72)

    results = []
    for i, seed in enumerate(SEEDS):
        rep_dir = work_dir / f"replicate_{i+1:02d}_seed{seed}"
        output_dir = rep_dir / "output"
        xmls = sorted(output_dir.glob("output*.xml"))
        if not xmls:
            print(f"\n  Replicate {i+1} (seed={seed}) — NO OUTPUT FILES")
            results.append(None)
            continue
        print(f"\n  Evaluating replicate {i+1} (seed={seed})  [{len(xmls)} snapshots]...")
        try:
            criteria, details = evaluate_replicate(rep_dir, seed)
            results.append((criteria, details))
        except Exception as e:
            print(f"    ERROR: {e}")
            results.append(None)

    # Print per-replicate results
    print("\n" + "=" * 72)
    print("  PER-REPLICATE RESULTS")
    print("=" * 72)
    for i, seed in enumerate(SEEDS):
        hdr = f"Replicate {i+1}  (seed={seed})"
        if results[i] is None:
            print(f"\n  {hdr}  — FAILED")
            continue
        criteria, details = results[i]
        all_pass = all(criteria.get(k, False) for k, _ in CRITERIA_NAMES)
        tag = "ALL PASS" if all_pass else "SOME FAIL"
        print(f"\n  {hdr}  — {tag}")
        for key, label in CRITERIA_NAMES:
            ok = criteria.get(key, False)
            detail = details.get(key, "")
            mark = "PASS" if ok else "FAIL"
            print(f"    [{mark}] {label}: {detail}")

    # Aggregate
    print("\n" + "=" * 72)
    print(f"  AGGREGATE  (>={QUORUM}/5 replicates must pass each criterion)")
    print("=" * 72)

    any_failure = False
    for key, label in CRITERIA_NAMES:
        passing = sum(
            1 for r in results
            if r is not None and r[0].get(key, False)
        )
        ok = passing >= QUORUM
        mark = "PASS" if ok else "FAIL"
        if not ok:
            any_failure = True
        print(f"  [{mark}] {label}  ({passing}/5 replicates)")

    print()
    if any_failure:
        print("  *** REALITY CHECK 1: FAIL ***")
        return 1
    else:
        print("  *** REALITY CHECK 1: PASS ***")
        return 0


if __name__ == "__main__":
    sys.exit(main())
