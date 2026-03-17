#!/usr/bin/env python3
"""
analyze_withdrawal.py — Detailed cell-state analysis during RC2 withdrawal phase.
Prints per-cell pressure decomposition (simple_pressure, mechanical_pressure,
neighbor count, local ECM) at key timepoints.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Dict, Any, Optional

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(PROJECT_ROOT))
from python.wrapper.output_parser import OutputParser

# ── Configuration ─────────────────────────────────────────────────────────────
OUT_DIR = PROJECT_ROOT / "build" / "rc2_fixG_seed42" / "replicate_01_seed42" / "output"
SAVE_INTERVAL = 360  # minutes between snapshots
HIF1A_BOOST = 0.14   # current hif1a_emt_boost
ECM_CAP = 0.10        # current ecm_emt_cap
EMT_OFF = 0.25        # current emt_off_threshold

# Days to sample
SUMMARY_DAYS = [14, 21, 24, 28, 31, 35, 38, 42]
DETAIL_DAYS  = [28, 31, 35, 42]


# ── Helper functions ──────────────────────────────────────────────────────────
def _row(matrix, labels, name):
    """Extract a single row from the cell matrix by label name."""
    entry = labels.get(name)
    if entry is None:
        return None
    idx = int(entry["index"])
    if idx < 0 or idx >= matrix.shape[0]:
        return None
    return matrix[idx, :]


def _sample_field_at_positions(positions, micro_coords, field_values):
    """Sample a microenvironment field at given cell positions using nearest voxel."""
    if positions.shape[0] == 0 or micro_coords.shape[0] == 0:
        return np.array([])
    # NumPy-only nearest-neighbor lookup (avoids SciPy dependency)
    pos2 = positions[:, :2]
    vox2 = micro_coords[:, :2]
    sampled = np.empty(pos2.shape[0], dtype=float)
    for i, p in enumerate(pos2):
        d2 = np.sum((vox2 - p) ** 2, axis=1)
        sampled[i] = field_values[int(np.argmin(d2))]
    return sampled


def _count_neighbors(positions, radius=30.0):
    """Count neighbors within `radius` µm for each cell."""
    if positions.shape[0] == 0:
        return np.array([], dtype=int)
    pos2 = positions[:, :2]
    r2 = radius * radius
    counts = np.zeros(pos2.shape[0], dtype=int)
    for i in range(pos2.shape[0]):
        d2 = np.sum((pos2 - pos2[i]) ** 2, axis=1)
        # include self then subtract 1
        counts[i] = int(np.count_nonzero(d2 <= r2) - 1)
    return counts


# ── Parse a single snapshot ───────────────────────────────────────────────────
def parse_snapshot(parser: OutputParser, idx: int) -> Dict[str, Any]:
    """Parse one output snapshot and return structured summary."""
    snap = parser._read_physicell_xml(OUT_DIR / f"output{idx:08d}.xml")
    matrix = snap["cell_matrix"]
    labels = snap["label_name_map"]
    micro_coords = snap["micro_coords"]
    micro_values = snap["micro_values"]
    n_cells = matrix.shape[1]

    # Live mask
    cell_type = _row(matrix, labels, "cell_type")
    dead = _row(matrix, labels, "dead")
    death_model = _row(matrix, labels, "current_death_model")

    live_mask = np.ones(n_cells, dtype=bool)
    if dead is not None:
        live_mask &= dead <= 0.5
    if death_model is not None:
        live_mask &= np.rint(death_model).astype(int) != 100

    ctype = np.rint(cell_type).astype(int) if cell_type is not None else np.full(n_cells, -1)
    live_tumor = live_mask & (ctype == 0)
    n_tumor = int(np.sum(live_tumor))

    result: Dict[str, Any] = {"n_tumor": n_tumor}

    if n_tumor == 0:
        result.update({
            "zeb1p_abcb1p": 0, "zeb1n_abcb1p": 0,
            "zeb1p_abcb1n": 0, "zeb1n_abcb1n": 0,
            "mean_tgfb": math.nan, "mean_ecm": math.nan,
            "mean_o2": math.nan, "mean_hif1a": math.nan,
            "mean_pressure": math.nan, "mean_drug": math.nan,
            "mean_drug_field": math.nan, "mean_nrf2": math.nan,
            "mean_simple_pressure": math.nan,
        })
        return result

    # Custom data rows
    zeb1 = _row(matrix, labels, "zeb1_active")
    abcb1 = _row(matrix, labels, "abcb1_active")
    mech_press = _row(matrix, labels, "mechanical_pressure")
    intracellular_drug = _row(matrix, labels, "intracellular_drug")
    hif1a = _row(matrix, labels, "hif1a_active")
    nrf2 = _row(matrix, labels, "nrf2_active")
    simple_press = _row(matrix, labels, "pressure")  # row 20: PhysiCell simple_pressure

    # ZEB1/ABCB1 quadrants
    z = (zeb1[live_tumor] > 0.5) if zeb1 is not None else np.zeros(n_tumor, dtype=bool)
    a = (abcb1[live_tumor] > 0.5) if abcb1 is not None else np.zeros(n_tumor, dtype=bool)
    result["zeb1p_abcb1p"] = int(np.sum(z & a))
    result["zeb1n_abcb1p"] = int(np.sum(~z & a))
    result["zeb1p_abcb1n"] = int(np.sum(z & ~a))
    result["zeb1n_abcb1n"] = int(np.sum(~z & ~a))

    # Positions
    positions = matrix[1:4, :].T  # rows 1-3 = x,y,z
    tumor_pos = positions[live_tumor]

    # Sample microenvironment
    tgfb = micro_values.get("tgfb")
    ecm = micro_values.get("ecm_density")
    drug_field = micro_values.get("drug")
    o2_field = micro_values.get("oxygen")

    if tgfb is not None and micro_coords.size > 0:
        result["mean_tgfb"] = float(np.mean(_sample_field_at_positions(tumor_pos, micro_coords, tgfb)))
    else:
        result["mean_tgfb"] = math.nan

    if ecm is not None and micro_coords.size > 0:
        result["mean_ecm"] = float(np.mean(_sample_field_at_positions(tumor_pos, micro_coords, ecm)))
    else:
        result["mean_ecm"] = math.nan

    if o2_field is not None and micro_coords.size > 0:
        result["mean_o2"] = float(np.mean(_sample_field_at_positions(tumor_pos, micro_coords, o2_field)))
    else:
        result["mean_o2"] = math.nan

    if drug_field is not None and micro_coords.size > 0:
        result["mean_drug_field"] = float(np.mean(_sample_field_at_positions(tumor_pos, micro_coords, drug_field)))
    else:
        result["mean_drug_field"] = math.nan

    # Cell-level custom data means
    result["mean_hif1a"] = float(np.mean(hif1a[live_tumor])) if hif1a is not None else math.nan
    result["mean_pressure"] = float(np.mean(mech_press[live_tumor])) if mech_press is not None else math.nan
    result["mean_simple_pressure"] = float(np.mean(simple_press[live_tumor])) if simple_press is not None else math.nan
    result["mean_drug"] = float(np.mean(intracellular_drug[live_tumor])) if intracellular_drug is not None else math.nan
    result["mean_nrf2"] = float(np.mean(nrf2[live_tumor])) if nrf2 is not None else math.nan

    # ── Per-cell detail arrays (for detailed breakdown) ───────────────────
    per_cell = []
    for i in range(n_tumor):
        cell_idx = np.where(live_tumor)[0][i]
        pos = positions[cell_idx]
        cell_info = {
            "pos": pos,
            "zeb1": bool(z[i]),
            "abcb1": bool(a[i]),
            "simple_pressure": float(simple_press[cell_idx]) if simple_press is not None else math.nan,
            "mechanical_pressure": float(mech_press[cell_idx]) if mech_press is not None else math.nan,
            "hif1a": float(hif1a[cell_idx]) if hif1a is not None else math.nan,
            "drug_i": float(intracellular_drug[cell_idx]) if intracellular_drug is not None else math.nan,
        }
        # Local ECM at this cell's position
        if ecm is not None and micro_coords.size > 0:
            cell_info["local_ecm"] = float(_sample_field_at_positions(
                pos.reshape(1, -1), micro_coords, ecm)[0])
        else:
            cell_info["local_ecm"] = math.nan
        per_cell.append(cell_info)

    # Neighbor counts
    neighbor_counts = _count_neighbors(tumor_pos, radius=30.0)
    for i, nc in enumerate(neighbor_counts):
        per_cell[i]["neighbors_30um"] = int(nc)

    result["per_cell"] = per_cell
    return result


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    if not OUT_DIR.exists():
        print(f"ERROR: Output directory not found: {OUT_DIR}")
        return 1

    parser = OutputParser(OUT_DIR)

    print("=" * 80)
    print("  RC2 WITHDRAWAL PHASE — DETAILED CELL-STATE ANALYSIS")
    print(f"  Run: hif1a_boost={HIF1A_BOOST}, ecm_cap={ECM_CAP}, off_thresh={EMT_OFF}")
    print(f"  Fix-G: crowding_base=0.3, tumor_CI=15.0, stromal_CI=0.8, dkc=0.0002")
    print("=" * 80)

    # ── Summary table ─────────────────────────────────────────────────────
    results = {}
    for day in SUMMARY_DAYS:
        t_min = day * 1440.0
        idx = int(t_min / SAVE_INTERVAL)
        try:
            r = parse_snapshot(parser, idx)
            results[day] = r
        except Exception as e:
            print(f"  WARNING: Could not parse day {day} (snap {idx}): {e}")

    print()
    print(f"  {'Day':>4}  {'Total':>5}  {'Z+A+':>5}  {'Z-A+':>5}  {'Z+A-':>5}  "
          f"{'Z-A-':>5}  {'TGFβ':>7}  {'ECM':>6}  {'O₂':>6}  {'HIF1α':>5}  "
          f"{'s_press':>7}  {'m_press':>7}  Phase")
    print(f"  {'────':>4}  {'─────':>5}  {'────':>5}  {'────':>5}  {'────':>5}  "
          f"{'────':>5}  {'─────':>7}  {'────':>6}  {'──':>6}  {'─────':>5}  "
          f"{'───────':>7}  {'───────':>7}  ─────")

    for day in SUMMARY_DAYS:
        r = results.get(day)
        if r is None:
            continue
        phase = "barrier" if day < 14 else ("drug-ON" if day <= 28 else "WITHDRAWAL")
        tgfb_s = f"{r['mean_tgfb']:.4f}" if math.isfinite(r.get('mean_tgfb', math.nan)) else "N/A"
        ecm_s = f"{r['mean_ecm']:.3f}" if math.isfinite(r.get('mean_ecm', math.nan)) else "N/A"
        o2_s = f"{r['mean_o2']:.1f}" if math.isfinite(r.get('mean_o2', math.nan)) else "N/A"
        hif_s = f"{r['mean_hif1a']:.2f}" if math.isfinite(r.get('mean_hif1a', math.nan)) else "N/A"
        sp_s = f"{r['mean_simple_pressure']:.4f}" if math.isfinite(r.get('mean_simple_pressure', math.nan)) else "N/A"
        mp_s = f"{r['mean_pressure']:.4f}" if math.isfinite(r.get('mean_pressure', math.nan)) else "N/A"

        print(f"  {day:4d}  {r['n_tumor']:5d}  {r['zeb1p_abcb1p']:5d}  "
              f"{r['zeb1n_abcb1p']:5d}  {r['zeb1p_abcb1n']:5d}  {r['zeb1n_abcb1n']:5d}  "
              f"{tgfb_s:>7}  {ecm_s:>6}  {o2_s:>6}  {hif_s:>5}  {sp_s:>7}  {mp_s:>7}  {phase}")

    # ── Detailed per-cell breakdown ───────────────────────────────────────
    print()
    print("=" * 80)
    print("  PER-CELL PRESSURE DECOMPOSITION (withdrawal phase)")
    print("  Formula: mechanical_pressure = simple_pressure × (crowding_base + ecm × collagen_frac × compaction_strength)")
    print("           crowding_base=0.3, compaction_strength=0.6, collagen_frac=0.4")
    print("           tumor contact_inhibition_threshold=15.0, stromal=0.8")
    print("=" * 80)

    for day in DETAIL_DAYS:
        r = results.get(day)
        if r is None:
            continue
        phase = "drug-ON" if day <= 28 else "WITHDRAWAL"
        n = r["n_tumor"]
        per_cell = r.get("per_cell", [])

        print(f"\n  ┌─ DAY {day} ({phase}) ─ {n} tumor cells ─────────────────────")

        if n == 0:
            print(f"  │ No survivors.")
            print(f"  └──────────────────────────────────────────────────")
            continue

        # Quadrant summary
        print(f"  │ ZEB1+: {r['zeb1p_abcb1n'] + r['zeb1p_abcb1p']}   ZEB1-: {r['zeb1n_abcb1n'] + r['zeb1n_abcb1p']}   "
              f"ABCB1+: {r['zeb1p_abcb1p'] + r['zeb1n_abcb1p']}   ABCB1-: {r['zeb1p_abcb1n'] + r['zeb1n_abcb1n']}")

        # Signal estimate
        ecm_boost = min(0.4 * max(0, r.get('mean_ecm', 0) - 0.3) / 0.3, ECM_CAP)
        hif_boost = HIF1A_BOOST if r.get('mean_hif1a', 0) > 0.5 else 0.0
        tgfb_val = r.get('mean_tgfb', 0) if math.isfinite(r.get('mean_tgfb', math.nan)) else 0
        signal = tgfb_val + hif_boost + ecm_boost
        signal_status = f"> {EMT_OFF} → ZEB1 ON" if signal > EMT_OFF else f"< {EMT_OFF} → ZEB1 REVERTS"
        print(f"  │ Signal: tgfb({tgfb_val:.4f}) + hif({hif_boost:.2f}) + ecm({ecm_boost:.3f}) = {signal:.4f}  ({signal_status})")

        # Per-cell table header
        print(f"  │")
        print(f"  │  {'#':>3}  {'x':>7}  {'y':>7}  {'Z':>1}{'A':>2}  "
              f"{'s_press':>8}  {'m_press':>8}  {'nbrs':>4}  {'ECM':>7}  "
              f"{'ECM_amp':>7}  {'drug_i':>8}  {'HIF':>4}")
        print(f"  │  {'───':>3}  {'───────':>7}  {'───────':>7}  {'──':>2}{' ':>1}  "
              f"{'────────':>8}  {'────────':>8}  {'────':>4}  {'───────':>7}  "
              f"{'───────':>7}  {'────────':>8}  {'────':>4}")

        for i, c in enumerate(per_cell):
            z_flag = "+" if c["zeb1"] else "-"
            a_flag = "+" if c["abcb1"] else "-"
            sp = c["simple_pressure"]
            mp = c["mechanical_pressure"]
            ecm_local = c["local_ecm"]
            nbrs = c["neighbors_30um"]
            drug_i = c["drug_i"]
            hif = c["hif1a"]

            # Compute ECM amplification factor: crowding_base + ecm * collagen_frac * compaction_strength
            ecm_amp = 0.3 + ecm_local * 0.4 * 0.6 if math.isfinite(ecm_local) else math.nan

            sp_s = f"{sp:.5f}" if math.isfinite(sp) else "N/A"
            mp_s = f"{mp:.5f}" if math.isfinite(mp) else "N/A"
            ecm_s = f"{ecm_local:.4f}" if math.isfinite(ecm_local) else "N/A"
            amp_s = f"{ecm_amp:.4f}" if math.isfinite(ecm_amp) else "N/A"
            drug_s = f"{drug_i:.6f}" if math.isfinite(drug_i) else "N/A"
            hif_s = f"{hif:.1f}" if math.isfinite(hif) else "N/A"

            print(f"  │  {i+1:3d}  {c['pos'][0]:7.1f}  {c['pos'][1]:7.1f}  "
                  f"{z_flag}{a_flag:>2}  {sp_s:>8}  {mp_s:>8}  {nbrs:4d}  "
                  f"{ecm_s:>7}  {amp_s:>7}  {drug_s:>8}  {hif_s:>4}")

        # Summary stats
        sps = [c["simple_pressure"] for c in per_cell if math.isfinite(c["simple_pressure"])]
        mps = [c["mechanical_pressure"] for c in per_cell if math.isfinite(c["mechanical_pressure"])]
        ecms = [c["local_ecm"] for c in per_cell if math.isfinite(c["local_ecm"])]
        nbrs_all = [c["neighbors_30um"] for c in per_cell]

        print(f"  │")
        if sps:
            print(f"  │  simple_pressure:  min={min(sps):.4f}  mean={np.mean(sps):.4f}  max={max(sps):.4f}")
        if mps:
            print(f"  │  mech_pressure:    min={min(mps):.4f}  mean={np.mean(mps):.4f}  max={max(mps):.4f}  "
                  f"(tumor CI @ 15.0, stroma CI @ 0.8)")
        if ecms:
            print(f"  │  local ECM:        min={min(ecms):.4f}  mean={np.mean(ecms):.4f}  max={max(ecms):.4f}")
        if nbrs_all:
            print(f"  │  neighbors(30µm):  min={min(nbrs_all)}  mean={np.mean(nbrs_all):.1f}  max={max(nbrs_all)}")

        # Diagnosis: is pressure from crowding or ECM amplification?
        if sps and ecms:
            mean_sp = np.mean(sps)
            mean_ecm_v = np.mean(ecms)
            ecm_amp_factor = 0.3 + mean_ecm_v * 0.4 * 0.6
            pure_crowd = mean_sp * 0.3  # pressure from crowding_base alone
            ecm_contrib = mean_sp * mean_ecm_v * 0.4 * 0.6  # pressure from ECM amplification
            print(f"  │")
            print(f"  │  DECOMPOSITION (mean):")
            print(f"  │    crowding_base contribution:  {pure_crowd:.5f}  ({100*pure_crowd/(pure_crowd+ecm_contrib) if (pure_crowd+ecm_contrib) > 0 else 0:.0f}%)")
            print(f"  │    ECM amplification:           {ecm_contrib:.5f}  ({100*ecm_contrib/(pure_crowd+ecm_contrib) if (pure_crowd+ecm_contrib) > 0 else 0:.0f}%)")
            print(f"  │    total mech_pressure ≈ sp × ({0.3:.2f} + {mean_ecm_v:.3f}×0.4×0.6) = {mean_sp:.4f} × {ecm_amp_factor:.4f} = {mean_sp*ecm_amp_factor:.5f}")
            if np.mean(mps) > 5.0 and np.mean(nbrs_all) < 3:
                print(f"  │")
                print(f"  │  ⚠  HIGH PRESSURE ({np.mean(mps):.3f}) with FEW NEIGHBORS ({np.mean(nbrs_all):.1f})")
                print(f"  │     → Pressure is from ECM stiffness, NOT cell crowding")
                print(f"  │     → Consider lowering crowding_base_pressure or mechanical_compaction_strength")

        print(f"  └──────────────────────────────────────────────────")

    # ── Final diagnosis ───────────────────────────────────────────────────
    print()
    print("=" * 80)
    print("  DIAGNOSIS")
    print("=" * 80)

    d28 = results.get(28, {})
    d42 = results.get(42, {})
    total_28 = d28.get("n_tumor", 0)
    total_42 = d42.get("n_tumor", 0)
    zeb1_28 = d28.get("zeb1p_abcb1n", 0) + d28.get("zeb1p_abcb1p", 0)

    print(f"\n  At d28: {total_28} survivors, {zeb1_28} are ZEB1+ ({100*zeb1_28/total_28:.0f}%)" if total_28 > 0 else "\n  At d28: 0 survivors")
    print(f"  At d42: {total_42} survivors")

    if total_42 > total_28:
        print(f"\n  ✓ REGROWTH DETECTED: {total_28} → {total_42}")
    elif total_42 > 0:
        print(f"\n  ✗ Still shrinking: {total_28} → {total_42} (no regrowth)")
    else:
        print(f"\n  ✗ EXTINCTION: {total_28} → 0")

    # Signal at d28
    if total_28 > 0:
        mean_tgfb_28 = d28.get("mean_tgfb", math.nan)
        mean_ecm_28 = d28.get("mean_ecm", math.nan)
        mean_hif_28 = d28.get("mean_hif1a", math.nan)
        mean_o2_28 = d28.get("mean_o2", math.nan)
        if math.isfinite(mean_tgfb_28) and math.isfinite(mean_ecm_28):
            ecm_boost = min(0.4 * max(0, mean_ecm_28 - 0.3) / 0.3, ECM_CAP)
            hif_boost = HIF1A_BOOST if (math.isfinite(mean_hif_28) and mean_hif_28 > 0.5) else 0.0
            total_sig = mean_tgfb_28 + hif_boost + ecm_boost
            print(f"\n  At d28 survivor positions:")
            print(f"    TGF-β    = {mean_tgfb_28:.5f}")
            print(f"    ECM      = {mean_ecm_28:.5f} → ECM boost = {ecm_boost:.3f}")
            print(f"    HIF1α    = {mean_hif_28:.2f} → HIF1α boost = {hif_boost:.3f}")
            if math.isfinite(mean_o2_28):
                print(f"    O₂       = {mean_o2_28:.2f} mmHg {'(HYPOXIC)' if mean_o2_28 < 15 else '(normoxic)'}")
            print(f"    Signal   = {mean_tgfb_28:.4f} + {hif_boost:.3f} + {ecm_boost:.3f} = {total_sig:.4f}")
            print(f"    emt_off_threshold = {EMT_OFF}")
            if total_sig > EMT_OFF:
                print(f"    Signal ({total_sig:.3f}) > emt_off_threshold ({EMT_OFF}) → ZEB1 LOCKED ON")
            else:
                print(f"    Signal ({total_sig:.3f}) < emt_off_threshold ({EMT_OFF}) → ZEB1 CAN REVERT ✓")

    print()
    print("=" * 80)
    return 0


if __name__ == "__main__":
    sys.exit(main())
