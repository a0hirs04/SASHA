#!/usr/bin/env python3
"""
Diagnostic script for RC2 seed 42 — reads existing output, prints drug/apoptosis metrics.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(PROJECT_ROOT))

from python.wrapper.output_parser import OutputParser
from python.wrapper.workdir_utils import default_reality_check_dir

OUTPUT_DIR = default_reality_check_dir(PROJECT_ROOT, "reality_check_2") / "replicate_01_seed42" / "output"
TIMEPOINTS = {
    "t_pre (day 14)": 20160.0,
    "t_end (day 28)": 40320.0,
}


def _sample_nearest(points, voxel_coords, values):
    if points.size == 0 or voxel_coords.size == 0:
        return np.array([], dtype=float)
    out = np.empty(points.shape[0], dtype=float)
    chunk = 500
    for i in range(0, points.shape[0], chunk):
        p = points[i : i + chunk]
        d2 = np.sum((p[:, None, :] - voxel_coords[None, :, :]) ** 2, axis=2)
        out[i : i + chunk] = values[np.argmin(d2, axis=1)]
    return out


def _row(matrix, labels, name):
    e = labels.get(name)
    if e is None:
        return None
    idx = int(e["index"])
    if idx < 0 or idx >= matrix.shape[0]:
        return None
    return matrix[idx, :]


def diagnose_snapshot(parser, xml_path):
    snap = parser._read_physicell_xml(xml_path)
    t = snap["time"]
    matrix = snap["cell_matrix"]
    labels = snap["label_name_map"]
    micro_coords = snap["micro_coords"]
    micro_values = snap["micro_values"]

    cell_type = _row(matrix, labels, "cell_type")
    dead = _row(matrix, labels, "dead")
    death_model = _row(matrix, labels, "current_death_model")
    death_rates = _row(matrix, labels, "death_rates")  # apoptosis rate (row 28)

    n_cells = matrix.shape[1]
    ctype_int = np.rint(cell_type).astype(int) if cell_type is not None else np.zeros(n_cells, int)

    tumor_mask = ctype_int == 0
    dead_mask = (dead > 0.5) if dead is not None else np.zeros(n_cells, bool)
    apoptotic_mask = (np.rint(death_model).astype(int) == 100) if death_model is not None else np.zeros(n_cells, bool)
    live_tumor = tumor_mask & ~dead_mask & ~apoptotic_mask
    dead_tumor = tumor_mask & dead_mask

    n_live = int(np.sum(live_tumor))
    n_dead = int(np.sum(dead_tumor))
    n_total_tumor = int(np.sum(tumor_mask))

    # Peritumoral ECM
    pos_row = labels.get("position")
    if pos_row:
        pidx = int(pos_row["index"])
        positions = matrix[pidx : pidx + 3, :].T
    else:
        positions = np.empty((0, 3))

    tumor_pos = positions[live_tumor] if live_tumor.any() else np.empty((0, 3))
    centroid = np.nanmean(tumor_pos, axis=0) if tumor_pos.shape[0] > 0 else np.zeros(3)
    tumor_radius = float(np.max(np.linalg.norm(tumor_pos - centroid, axis=1))) if tumor_pos.shape[0] > 1 else 0.0

    ecm = micro_values.get("ecm_density")
    peri_ecm = math.nan
    if ecm is not None and micro_coords.size > 0 and tumor_pos.shape[0] > 0:
        vd = np.linalg.norm(micro_coords - centroid, axis=1)
        shell = (vd >= tumor_radius) & (vd <= tumor_radius + 100.0)
        if shell.any():
            peri_ecm = float(np.nanmean(ecm[shell]))

    # Extracellular drug near tumor
    drug_field = micro_values.get("drug")
    drug_at_tumor = math.nan
    if drug_field is not None and micro_coords.size > 0 and tumor_pos.shape[0] > 0:
        sampled = _sample_nearest(tumor_pos, micro_coords, drug_field)
        drug_at_tumor = float(np.nanmean(sampled))

    # Intracellular drug
    ic_drug = _row(matrix, labels, "intracellular_drug")
    ic_mean = ic_max = math.nan
    ic_count_above = 0
    if ic_drug is not None and n_live > 0:
        vals = ic_drug[live_tumor]
        ic_mean = float(np.nanmean(vals))
        ic_max = float(np.nanmax(vals))
        ic_count_above = int(np.sum(vals > 0.01))

    # Apoptosis rate (death_rates label, index 28 = apoptosis rate row)
    apop_rate_mean = apop_rate_max = math.nan
    apop_frac = 0.0
    if death_rates is not None and n_live > 0:
        vals = death_rates[live_tumor]
        apop_rate_mean = float(np.nanmean(vals))
        apop_rate_max = float(np.nanmax(vals))

    if n_total_tumor > 0:
        apop_frac = n_dead / n_total_tumor

    # ABCB1 / NRF2
    abcb1 = _row(matrix, labels, "abcb1_active")
    nrf2 = _row(matrix, labels, "nrf2_active")
    frac_abcb1 = float(np.mean(abcb1[live_tumor] > 0.5)) if (abcb1 is not None and n_live > 0) else math.nan
    frac_nrf2 = float(np.mean(nrf2[live_tumor] > 0.5)) if (nrf2 is not None and n_live > 0) else math.nan

    print(f"  Live tumor cells:         {n_live}")
    print(f"  Dead tumor cells:         {n_dead}")
    print(f"  Total tumor cells:        {n_total_tumor}")
    print(f"  Peritumoral ECM (mean):   {peri_ecm:.4f}")
    print(f"  Tumor radius:             {tumor_radius:.1f} um")
    print()
    print(f"  Extracellular drug at tumor (mean): {drug_at_tumor:.6f}")
    print()
    print(f"  Intracellular drug (live tumor):")
    print(f"    mean:   {ic_mean:.6f}")
    print(f"    max:    {ic_max:.6f}")
    print(f"    >0.01:  {ic_count_above} / {n_live} ({100*ic_count_above/max(n_live,1):.1f}%)")
    print()
    print(f"  Apoptosis rate (live tumor):")
    print(f"    mean:   {apop_rate_mean:.6f} 1/min")
    print(f"    max:    {apop_rate_max:.6f} 1/min")
    print(f"  Dead fraction (all tumor): {apop_frac:.4f} ({n_dead}/{n_total_tumor})")
    print()
    print(f"  ABCB1 active fraction:     {frac_abcb1:.4f}")
    print(f"  NRF2 active fraction:      {frac_nrf2:.4f}")

    return n_live


def main():
    if not OUTPUT_DIR.exists():
        print(f"ERROR: Output directory not found: {OUTPUT_DIR}")
        return 1

    parser = OutputParser(OUTPUT_DIR)
    xmls = sorted(OUTPUT_DIR.glob("output*.xml"))
    if not xmls:
        print("ERROR: No snapshot XMLs found")
        return 1

    # Map time -> xml path
    snap_times = []
    for xml in xmls:
        snap = parser._read_physicell_xml(xml)
        snap_times.append((float(snap["time"]), xml))

    results = {}
    for label, target_time in TIMEPOINTS.items():
        best_time, best_xml = min(snap_times, key=lambda s: abs(s[0] - target_time))
        print("=" * 60)
        print(f"  {label}  (actual t={best_time:.0f}, file={best_xml.name})")
        print("=" * 60)
        results[label] = diagnose_snapshot(parser, best_xml)
        print()

    print("=" * 60)
    print("  DIAGNOSIS")
    print("=" * 60)
    labels = list(TIMEPOINTS.keys())
    n_pre = results[labels[0]]
    n_end = results[labels[1]]
    if n_end == 0:
        print(f"  Tumor ERADICATED: {n_pre} -> {n_end} (100% kill)")
        print("  Root cause: intracellular drug never decays (drug_natural_decay_rate=0.0)")
    elif n_end < 0.9 * n_pre:
        print(f"  Partial response: {n_pre} -> {n_end} ({100*(1-n_end/max(n_pre,1)):.1f}% reduction)")
    else:
        print(f"  Minimal response: {n_pre} -> {n_end}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
