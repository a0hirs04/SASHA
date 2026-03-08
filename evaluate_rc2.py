#!/usr/bin/env python3
"""
Evaluate RC2 criteria on the completed rc2_full_seed42 run.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(PROJECT_ROOT))
from python.wrapper.output_parser import OutputParser

# ---------------------------------------------------------------------------
OUT_DIR = PROJECT_ROOT / "build" / "rc2_full_seed42" / "replicate_01_seed42" / "output"
T_PRE       = 20160.0   # day 14
T_TREAT_END = 40320.0   # day 28
T_POST      = 60480.0   # day 42
SAVE_INTERVAL = 360


def _row(matrix, labels, name):
    e = labels.get(name)
    if e is None:
        return None
    idx = int(e["index"])
    if idx < 0 or idx >= matrix.shape[0]:
        return None
    return matrix[idx, :]


def _sample_ecm_at_positions(positions, micro_coords, ecm):
    """Sample ECM density at cell positions via nearest voxel."""
    vals = []
    mc = micro_coords[:, :2] if micro_coords.shape[1] > 2 else micro_coords
    for i in range(0, positions.shape[0], 500):
        p = positions[i:i+500, :2]
        d2 = np.sum((p[:, None, :] - mc[None, :, :]) ** 2, axis=2)
        vals.append(ecm[np.argmin(d2, axis=1)])
    return np.concatenate(vals) if vals else np.array([])


def parse_snap(parser, snap):
    matrix = snap["cell_matrix"]
    labels = snap["label_name_map"]
    micro_coords = snap["micro_coords"]
    micro_values = snap["micro_values"]

    cell_type = _row(matrix, labels, "cell_type")
    dead = _row(matrix, labels, "dead")
    death_model = _row(matrix, labels, "current_death_model")
    n_cells = matrix.shape[1]

    live_mask = np.ones(n_cells, dtype=bool)
    if dead is not None:
        live_mask &= dead <= 0.5
    if death_model is not None:
        live_mask &= np.rint(death_model).astype(int) != 100

    ctype = np.rint(cell_type).astype(int) if cell_type is not None else np.full(n_cells, -1)
    live_tumor = live_mask & (ctype == 0)
    live_stroma = live_mask & (ctype == 1)

    acta2 = _row(matrix, labels, "ZEB1")  # stromal offset map
    live_caf = live_stroma & (acta2 > 0.5) if acta2 is not None else np.zeros_like(live_stroma)

    pos = parser._get_positions(matrix, labels)
    tumor_pos = pos[live_tumor] if pos.size and live_tumor.any() else np.empty((0, 3))
    centroid = np.nanmean(tumor_pos, axis=0) if tumor_pos.shape[0] > 0 else np.zeros(3)
    tumor_radius = (
        float(np.max(np.linalg.norm(tumor_pos - centroid, axis=1)))
        if tumor_pos.shape[0] > 1 else 0.0
    )

    ecm = micro_values.get("ecm_density")
    peri_ecm = math.nan
    ecm_at_tumor = math.nan
    if ecm is not None and micro_coords.size > 0 and tumor_pos.shape[0] > 0:
        vd = np.linalg.norm(micro_coords - centroid, axis=1)
        shell = (vd >= tumor_radius) & (vd <= tumor_radius + 100.0)
        if shell.any():
            peri_ecm = float(np.nanmean(ecm[shell]))
        ecm_vals = _sample_ecm_at_positions(tumor_pos, micro_coords, ecm)
        if ecm_vals.size > 0:
            ecm_at_tumor = float(np.mean(ecm_vals))

    abcb1 = _row(matrix, labels, "abcb1_active")
    if abcb1 is None:
        abcb1 = _row(matrix, labels, "ABCB1")
    frac_abcb1 = (
        float(np.mean(abcb1[live_tumor] > 0.5))
        if (abcb1 is not None and live_tumor.any()) else math.nan
    )

    # ZEB1 fraction among live tumor cells
    zeb1 = _row(matrix, labels, "zeb1_active")
    if zeb1 is None:
        zeb1 = _row(matrix, labels, "ZEB1")
    frac_zeb1 = (
        float(np.mean(zeb1[live_tumor] > 0.5))
        if (zeb1 is not None and live_tumor.any()) else math.nan
    )

    return {
        "time": float(snap["time"]),
        "n_tumor": int(np.sum(live_tumor)),
        "n_stroma": int(np.sum(live_stroma)),
        "n_caf": int(np.sum(live_caf)),
        "peri_ecm": peri_ecm,
        "ecm_at_tumor": ecm_at_tumor,
        "frac_abcb1": frac_abcb1,
        "frac_zeb1": frac_zeb1,
        "tumor_radius": tumor_radius,
    }


def main():
    parser = OutputParser(OUT_DIR)
    xmls = sorted(OUT_DIR.glob("output*.xml"))
    print(f"Total snapshots: {len(xmls)}")

    snap_idx_pre   = int(T_PRE / SAVE_INTERVAL)        # 56
    snap_idx_treat = int(T_TREAT_END / SAVE_INTERVAL)   # 112
    snap_idx_post  = int(T_POST / SAVE_INTERVAL)        # 168

    def get_snap(idx):
        fname = OUT_DIR / f"output{idx:08d}.xml"
        if fname.exists():
            return parser._read_physicell_xml(fname)
        raise FileNotFoundError(f"Snapshot not found: {fname}")

    print()
    print("=" * 72)
    print("  RC2 EVALUATION — seed 42 (single replicate)")
    print("=" * 72)

    print(f"\n  Parsing snap {snap_idx_pre} (day 14, pre-treatment)...")
    s_pre = parse_snap(parser, get_snap(snap_idx_pre))
    print(f"    tumor={s_pre['n_tumor']}  stroma={s_pre['n_stroma']}  "
          f"caf={s_pre['n_caf']}  peri_ecm={s_pre['peri_ecm']:.4f}")

    print(f"  Parsing snap {snap_idx_treat} (day 28, treatment end)...")
    s_treat = parse_snap(parser, get_snap(snap_idx_treat))
    print(f"    tumor={s_treat['n_tumor']}  stroma={s_treat['n_stroma']}  "
          f"caf={s_treat['n_caf']}  peri_ecm={s_treat['peri_ecm']:.4f}")

    print(f"  Parsing snap {snap_idx_post} (day 42, post-withdrawal)...")
    s_post = parse_snap(parser, get_snap(snap_idx_post))
    print(f"    tumor={s_post['n_tumor']}  stroma={s_post['n_stroma']}  "
          f"caf={s_post['n_caf']}  peri_ecm={s_post['peri_ecm']:.4f}")

    # ---- Extended tumor timeline ----
    print()
    print("=" * 72)
    print("  TUMOR TIMELINE")
    print("=" * 72)
    timeline_days = [14, 17, 21, 24, 28, 31, 35, 42]
    print(f"  {'Day':>5}  {'Tumor':>6}  {'ABCB1%':>7}  {'ZEB1%':>6}  Phase")
    print(f"  {'---':>5}  {'-----':>6}  {'------':>7}  {'-----':>6}  -----")
    for day in timeline_days:
        t_min = day * 1440.0
        idx = int(t_min / SAVE_INTERVAL)
        fname = OUT_DIR / f"output{idx:08d}.xml"
        if not fname.exists():
            print(f"  {day:5d}  {'N/A':>6}")
            continue
        s = parse_snap(parser, get_snap(idx))
        phase = "barrier" if day < 14 else ("drug-ON" if day <= 28 else "regrowth")
        abcb1_str = f"{s['frac_abcb1']:.1%}" if math.isfinite(s['frac_abcb1']) else "N/A"
        zeb1_str = f"{s['frac_zeb1']:.1%}" if math.isfinite(s['frac_zeb1']) else "N/A"
        print(f"  {day:5d}  {s['n_tumor']:6d}  {abcb1_str:>7}  {zeb1_str:>6}  {phase}")

    print()
    print("=" * 72)
    print("  CRITERIA EVALUATION")
    print("=" * 72)

    # RC2-1: Partial response (>=10% tumor reduction)
    reduction = 1.0 - s_treat["n_tumor"] / max(s_pre["n_tumor"], 1)
    c1 = reduction >= 0.10
    mark = "PASS" if c1 else "FAIL"
    print(f"  [{mark}] RC2-1  Partial response")
    print(f"         tumor {s_pre['n_tumor']} → {s_treat['n_tumor']}  "
          f"reduction={reduction:.1%}  (need ≥10%)")

    # RC2-2: Not eradicated
    c2 = s_treat["n_tumor"] > 0
    mark = "PASS" if c2 else "FAIL"
    print(f"  [{mark}] RC2-2  Not eradicated")
    print(f"         survivors at treatment end = {s_treat['n_tumor']}")

    # RC2-3: Barrier persists (ECM >= 80% of pre)
    threshold = 0.80 * s_pre["peri_ecm"] if math.isfinite(s_pre["peri_ecm"]) else math.nan
    c3 = (math.isfinite(s_treat["peri_ecm"])
          and math.isfinite(threshold)
          and s_treat["peri_ecm"] >= threshold)
    mark = "PASS" if c3 else "FAIL"
    print(f"  [{mark}] RC2-3  Barrier persists")
    print(f"         peri_ecm {s_pre['peri_ecm']:.4f} → {s_treat['peri_ecm']:.4f}  "
          f"(need ≥{threshold:.4f})")

    # RC2-4: Spatial sanctuary
    c4 = (math.isfinite(s_treat["ecm_at_tumor"])
          and math.isfinite(s_pre["ecm_at_tumor"])
          and s_treat["ecm_at_tumor"] > s_pre["ecm_at_tumor"])
    mark = "PASS" if c4 else "FAIL"
    print(f"  [{mark}] RC2-4  Spatial sanctuary")
    print(f"         ecm@survivors={s_treat['ecm_at_tumor']:.4f} vs "
          f"ecm@tumor(pre)={s_pre['ecm_at_tumor']:.4f}")

    # RC2-5: ABCB1 emerges (SOFT)
    c5 = (math.isfinite(s_treat["frac_abcb1"])
          and math.isfinite(s_pre["frac_abcb1"])
          and s_treat["frac_abcb1"] > s_pre["frac_abcb1"])
    mark = "PASS" if c5 else "INFO"
    print(f"  [{mark}] RC2-5  ABCB1 resistance emerges [SOFT]")
    print(f"         frac_abcb1 {s_pre['frac_abcb1']:.4f} → {s_treat['frac_abcb1']:.4f}")

    # RC2-6: Regrowth
    c6 = s_post["n_tumor"] > s_treat["n_tumor"]
    mark = "PASS" if c6 else "FAIL"
    print(f"  [{mark}] RC2-6  Regrowth after withdrawal")
    print(f"         tumor {s_treat['n_tumor']} → {s_post['n_tumor']}")

    # Summary
    print()
    print("=" * 72)
    hard = [c1, c2, c3, c4, c6]
    n_pass = sum(hard)
    overall = all(hard)
    if overall:
        print("  *** RC2 (seed 42): ALL HARD CRITERIA PASS ✓ ***")
    else:
        print(f"  *** RC2 (seed 42): {n_pass}/5 HARD CRITERIA PASS ***")
        if not c1:
            print("    ✗ Drug is not effective enough (or too effective)")
        if not c2:
            print("    ✗ Tumor was eradicated — drug too strong")
        if not c3:
            print("    ✗ Barrier collapsed during treatment")
        if not c4:
            print("    ✗ No spatial sanctuary effect detected")
        if not c6:
            print("    ✗ No regrowth after drug withdrawal")
    if c5:
        print("  ABCB1 resistance: DETECTED ✓")
    else:
        print("  ABCB1 resistance: not detected (soft criterion, informational)")
    print("=" * 72)

    return 0 if overall else 1


if __name__ == "__main__":
    sys.exit(main())
