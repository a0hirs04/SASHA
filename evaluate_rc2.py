#!/usr/bin/env python3
"""Mechanistic RC2 evaluator.

Produces a detailed, interpretation-first report for each RC2 run:
  1) Tumor timeline with growth rates
  2) Explained RC2 criteria breakdown
  3) Resistance dynamics (ABCB1 / NRF2)
  4) Mechanism check (selection vs loophole)
  5) Drug dynamics (intra + extra)
  6) EMT / ZEB1 state
  7) Spatial / ECM analysis with p95 radius
  8) Final verdict (true biological pass vs false pass)
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(PROJECT_ROOT))
from python.wrapper.output_parser import OutputParser
from python.wrapper.workdir_utils import default_reality_check_dir


DEFAULT_OUT_DIR = (
    default_reality_check_dir(PROJECT_ROOT, "reality_check_2")
    / "replicate_01_seed42"
    / "output"
)
SAVE_INTERVAL = 360.0


def _row(matrix, labels, name):
    entry = labels.get(name)
    if entry is None:
        return None
    idx = int(entry["index"])
    if idx < 0 or idx >= matrix.shape[0]:
        return None
    return matrix[idx, :]


def _sample_field_at_positions(positions, micro_coords, field_vals):
    if positions.size == 0 or micro_coords.size == 0:
        return np.array([], dtype=float)
    out = []
    vox2 = micro_coords[:, :2] if micro_coords.shape[1] > 2 else micro_coords
    for i in range(0, positions.shape[0], 500):
        p = positions[i : i + 500, :2]
        d2 = np.sum((p[:, None, :] - vox2[None, :, :]) ** 2, axis=2)
        out.append(field_vals[np.argmin(d2, axis=1)])
    return np.concatenate(out) if out else np.array([], dtype=float)


def _fmt(v, nd=4):
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    if v is None or (isinstance(v, float) and (math.isnan(v) or math.isinf(v))):
        return "nan"
    return f"{float(v):.{nd}f}"


def _infer_seed_label(out_dir: Path) -> str:
    m = re.search(r"seed(\d+)", out_dir.as_posix())
    return m.group(1) if m else "?"


def _load_snapshot(parser: OutputParser, out_dir: Path, target_time_min: float):
    idx = int(round(target_time_min / SAVE_INTERVAL))
    for i in range(idx, -1, -1):
        f = out_dir / f"output{i:08d}.xml"
        if f.exists() and f.stat().st_size > 0:
            try:
                return parser._read_physicell_xml(f), i
            except Exception:
                pass
    raise FileNotFoundError(f"No readable snapshot at or before t={target_time_min} (idx={idx})")


def _parse_snapshot(parser: OutputParser, snap: dict):
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
    tumor_mask = ctype == 0
    live_tumor = live_mask & tumor_mask
    dead_tumor = tumor_mask & ~live_mask

    pos = parser._get_positions(matrix, labels)
    tumor_pos = pos[live_tumor] if pos.size and np.any(live_tumor) else np.empty((0, 3))
    centroid = np.nanmean(tumor_pos, axis=0) if tumor_pos.shape[0] > 0 else np.zeros(3)
    if tumor_pos.shape[0] > 1:
        d = np.linalg.norm(tumor_pos - centroid, axis=1)
        radius_p95 = float(np.percentile(d, 95))
        radius_max = float(np.max(d))
    else:
        radius_p95 = 0.0
        radius_max = 0.0

    # Signals / states
    abcb1 = _row(matrix, labels, "abcb1_active")
    if abcb1 is None:
        abcb1 = _row(matrix, labels, "ABCB1")
    nrf2 = _row(matrix, labels, "nrf2_active")
    if nrf2 is None:
        nrf2 = _row(matrix, labels, "NRF2")
    zeb1 = _row(matrix, labels, "zeb1_active")
    if zeb1 is None:
        zeb1 = _row(matrix, labels, "ZEB1")
    ic_drug = _row(matrix, labels, "intracellular_drug")

    frac_abcb1 = float(np.mean(abcb1[live_tumor] > 0.5)) if (abcb1 is not None and np.any(live_tumor)) else math.nan
    mean_abcb1 = float(np.mean(abcb1[live_tumor])) if (abcb1 is not None and np.any(live_tumor)) else math.nan
    mean_nrf2 = float(np.mean(nrf2[live_tumor])) if (nrf2 is not None and np.any(live_tumor)) else math.nan
    frac_zeb1 = float(np.mean(zeb1[live_tumor] > 0.5)) if (zeb1 is not None and np.any(live_tumor)) else math.nan
    mean_ic_drug = float(np.mean(ic_drug[live_tumor])) if (ic_drug is not None and np.any(live_tumor)) else math.nan

    # Extracellular drug near live tumor
    ec_drug_field = micro_values.get("drug")
    mean_ec_drug = math.nan
    if ec_drug_field is not None and micro_coords.size > 0 and tumor_pos.shape[0] > 0:
        sampled = _sample_field_at_positions(tumor_pos, micro_coords, ec_drug_field)
        if sampled.size > 0:
            mean_ec_drug = float(np.mean(sampled))

    # ECM metrics
    ecm_field = micro_values.get("ecm_density")
    peri_ecm = math.nan
    ecm_at_survivors = math.nan
    if ecm_field is not None and micro_coords.size > 0 and tumor_pos.shape[0] > 0:
        vd = np.linalg.norm(micro_coords - centroid, axis=1)
        shell = (vd >= radius_p95) & (vd <= radius_p95 + 100.0)
        if np.any(shell):
            peri_ecm = float(np.nanmean(ecm_field[shell]))
        sampled_ecm = _sample_field_at_positions(tumor_pos, micro_coords, ecm_field)
        if sampled_ecm.size > 0:
            ecm_at_survivors = float(np.mean(sampled_ecm))

    return {
        "time_min": float(snap["time"]),
        "n_live_tumor": int(np.sum(live_tumor)),
        "n_dead_tumor": int(np.sum(dead_tumor)),
        "frac_abcb1": frac_abcb1,
        "mean_abcb1": mean_abcb1,
        "mean_nrf2": mean_nrf2,
        "frac_zeb1": frac_zeb1,
        "mean_ic_drug": mean_ic_drug,
        "mean_ec_drug": mean_ec_drug,
        "peri_ecm": peri_ecm,
        "ecm_at_survivors": ecm_at_survivors,
        "tumor_radius_p95": radius_p95,
        "tumor_radius_max": radius_max,
    }


def _criteria(pre, treat, post):
    reduction = 1.0 - treat["n_live_tumor"] / max(pre["n_live_tumor"], 1)
    c1 = reduction >= 0.10
    c2 = treat["n_live_tumor"] > 0

    c3 = (
        math.isfinite(pre["peri_ecm"]) and math.isfinite(treat["peri_ecm"]) and
        treat["peri_ecm"] >= 0.80 * pre["peri_ecm"]
    )
    c4 = (
        math.isfinite(pre["ecm_at_survivors"]) and math.isfinite(treat["ecm_at_survivors"]) and
        treat["ecm_at_survivors"] > pre["ecm_at_survivors"]
    )
    c5 = (
        math.isfinite(pre["frac_abcb1"]) and math.isfinite(treat["frac_abcb1"]) and
        treat["frac_abcb1"] > pre["frac_abcb1"]
    )
    c6 = post["n_live_tumor"] > treat["n_live_tumor"]

    return {
        "RC2-1": (c1, reduction),
        "RC2-2": (c2, treat["n_live_tumor"]),
        "RC2-3": (c3, (pre["peri_ecm"], treat["peri_ecm"], 0.80 * pre["peri_ecm"] if math.isfinite(pre["peri_ecm"]) else math.nan)),
        "RC2-4": (c4, (pre["ecm_at_survivors"], treat["ecm_at_survivors"])),
        "RC2-5": (c5, (pre["frac_abcb1"], treat["frac_abcb1"])),
        "RC2-6": (c6, (treat["n_live_tumor"], post["n_live_tumor"])),
    }


def _report_one(run_name: str, out_dir: Path, timing: dict):
    parser = OutputParser(out_dir)

    day_to_min = {
        14: float(timing["drug_start_time"]),
        21: 21.0 * 1440.0,
        28: float(timing["drug_end_time"]),
        31: 31.0 * 1440.0,
        35: 35.0 * 1440.0,
        42: float(timing["max_time"]),
    }
    days = [14, 21, 28, 31, 35, 42]

    snaps = {}
    for d in days:
        snap, _ = _load_snapshot(parser, out_dir, day_to_min[d])
        snaps[d] = _parse_snapshot(parser, snap)

    pre = snaps[14]
    treat = snaps[28]
    post = snaps[42]
    crit = _criteria(pre, treat, post)

    print("\n" + "=" * 100)
    print(f"RC2 MECHANISTIC REPORT :: {run_name}")
    print(f"Output: {out_dir}")
    print("=" * 100)

    # 1) Tumor timeline
    print("\n1) Tumor Timeline (core signal)")
    print("  Day   LiveTumor  DeadTumor   Δcells/day")
    prev_day = None
    prev_live = None
    for d in days:
        live = snaps[d]["n_live_tumor"]
        dead = snaps[d]["n_dead_tumor"]
        if prev_day is None:
            rate = math.nan
        else:
            rate = (live - prev_live) / (d - prev_day)
        print(f"  {d:>3}   {live:>9}  {dead:>9}   {_fmt(rate, 2):>9}")
        prev_day = d
        prev_live = live

    shrink = treat["n_live_tumor"] < pre["n_live_tumor"]
    regrow = post["n_live_tumor"] > treat["n_live_tumor"]
    post_rate_1 = (snaps[31]["n_live_tumor"] - snaps[28]["n_live_tumor"]) / 3.0
    post_rate_2 = (snaps[42]["n_live_tumor"] - snaps[31]["n_live_tumor"]) / 11.0
    explosive = regrow and post_rate_2 > max(1.0, 1.5 * post_rate_1)
    print(f"\n  Interpretation: treatment shrink={'YES' if shrink else 'NO'}; "
          f"withdrawal regrowth={'YES' if regrow else 'NO'}; "
          f"regrowth profile={'explosive' if explosive else 'gradual/mixed'}.")

    # 2) Criteria breakdown
    print("\n2) RC2 Criteria Breakdown (explained)")
    c1, red = crit["RC2-1"]
    print(f"  RC2-1 Partial response: {'PASS' if c1 else 'FAIL'} | reduction={red:.2%} "
          f"(d14 {pre['n_live_tumor']} -> d28 {treat['n_live_tumor']}) | "
          f"{'Tumor reduced by >=10% during treatment.' if c1 else 'Tumor did not reduce enough during treatment.'}")

    c2, n_treat = crit["RC2-2"]
    print(f"  RC2-2 Not eradicated: {'PASS' if c2 else 'FAIL'} | survivors@d28={n_treat} | "
          f"{'Survivors remained after treatment.' if c2 else 'Treatment eradicated tumor.'}")

    c3, (peri_pre, peri_tr, peri_thr) = crit["RC2-3"]
    print(f"  RC2-3 Barrier persists: {'PASS' if c3 else 'FAIL'} | peri_ecm d14={_fmt(peri_pre)} d28={_fmt(peri_tr)} need>={_fmt(peri_thr)} | "
          f"{'Barrier remained above persistence threshold.' if c3 else 'Barrier collapsed below persistence threshold.'}")

    c4, (ecm_pre_surv, ecm_tr_surv) = crit["RC2-4"]
    print(f"  RC2-4 Spatial sanctuary: {'PASS' if c4 else 'FAIL'} | ecm@tumor d14={_fmt(ecm_pre_surv)} d28={_fmt(ecm_tr_surv)} | "
          f"{'Survivors occupy higher-ECM niches at treatment end.' if c4 else 'No ECM-enriched sanctuary signal at treatment end.'}")

    c5, (ab_pre, ab_tr) = crit["RC2-5"]
    print(f"  RC2-5 ABCB1 emergence (soft): {'PASS' if c5 else 'FAIL'} | frac_abcb1 d14={_fmt(ab_pre)} d28={_fmt(ab_tr)} | "
          f"{'Resistance fraction increased during treatment.' if c5 else 'Resistance fraction did not increase during treatment.'}")

    c6, (n_tr2, n_po2) = crit["RC2-6"]
    print(f"  RC2-6 Regrowth: {'PASS' if c6 else 'FAIL'} | d28={n_tr2} d42={n_po2} | "
          f"{'Tumor regrew after withdrawal.' if c6 else 'No regrowth after withdrawal.'}")

    # 3) Resistance dynamics
    print("\n3) Resistance Dynamics (critical)")
    print("  Day   ABCB1+%   mean_ABCB1   mean_NRF2")
    for d in [28, 29, 30, 31]:
        snap, _ = _load_snapshot(parser, out_dir, d * 1440.0)
        s = _parse_snapshot(parser, snap)
        print(f"  {d:>3}   {_fmt(100.0 * s['frac_abcb1'], 2):>7}%   {_fmt(s['mean_abcb1'], 4):>10}   {_fmt(s['mean_nrf2'], 4):>9}")
    ab28 = snaps[28]["mean_abcb1"]
    ab31 = _parse_snapshot(parser, _load_snapshot(parser, out_dir, 31 * 1440.0)[0])["mean_abcb1"]
    print(f"\n  Interpretation: ABCB1 persistence post-withdrawal is "
          f"{'maintained' if (math.isfinite(ab28) and math.isfinite(ab31) and ab31 >= 0.7 * max(ab28, 1e-9)) else 'weak/collapsing'}.")

    # 4) Mechanism check
    print("\n4) Mechanism Check")
    pref_surv = (math.isfinite(ab_pre) and math.isfinite(ab_tr) and (ab_tr - ab_pre) > 0.10)
    if regrow:
        ab42 = snaps[42]["mean_abcb1"]
        if math.isfinite(ab42) and ab42 > 0.5:
            regrowth_driver = "resistant lineage-dominant"
        elif math.isfinite(ab42) and ab42 < 0.2:
            regrowth_driver = "likely mixed/random survivors"
        else:
            regrowth_driver = "mixed lineage contribution"
    else:
        regrowth_driver = "no regrowth"

    stable_res = (
        math.isfinite(snaps[28]["mean_abcb1"]) and
        math.isfinite(snaps[42]["mean_abcb1"]) and
        snaps[42]["mean_abcb1"] >= 0.7 * max(snaps[28]["mean_abcb1"], 1e-9)
    )

    print(f"  Preferential resistant survival: {'YES' if pref_surv else 'NO/WEAK'}")
    print(f"  Regrowth driver: {regrowth_driver}")
    print(f"  Resistance stability: {'stable' if stable_res else 'transient/declining'}")

    # 5) Drug dynamics
    print("\n5) Drug Dynamics")
    print("  Day   mean_intra_drug   mean_extra_drug")
    intra = []
    extra = []
    for d in [28, 29, 30, 31]:
        snap, _ = _load_snapshot(parser, out_dir, d * 1440.0)
        s = _parse_snapshot(parser, snap)
        intra.append(s["mean_ic_drug"])
        extra.append(s["mean_ec_drug"])
        print(f"  {d:>3}      {_fmt(s['mean_ic_drug'], 6):>12}      {_fmt(s['mean_ec_drug'], 6):>12}")
    intra_clear = all(math.isfinite(intra[i]) and math.isfinite(intra[i + 1]) and intra[i + 1] <= intra[i] + 1e-12 for i in range(3))
    extra_clear = all(math.isfinite(extra[i]) and math.isfinite(extra[i + 1]) and extra[i + 1] <= extra[i] + 1e-12 for i in range(3))
    print(f"\n  Interpretation: intracellular clearing={'YES' if intra_clear else 'NO/MIXED'}; "
          f"extracellular clearing={'YES' if extra_clear else 'NO/MIXED'}.")

    # 6) EMT/ZEB1
    print("\n6) EMT / ZEB1 State")
    for d in [28, 31, 42]:
        print(f"  Day {d}: ZEB1+ fraction = {_fmt(100.0 * snaps[d]['frac_zeb1'], 2)}%")
    z28 = snaps[28]["frac_zeb1"]
    z42 = snaps[42]["frac_zeb1"]
    if math.isfinite(z28) and math.isfinite(z42):
        if z42 > 0.70:
            zmsg = "cells remain predominantly mesenchymal"
        elif z42 + 0.10 < z28:
            zmsg = "epithelial reversion occurred"
        else:
            zmsg = "partial/limited reversion"
    else:
        zmsg = "insufficient ZEB1 data"
    print(f"  Interpretation: {zmsg}.")

    # 7) Spatial / ECM
    print("\n7) Spatial / ECM Analysis")
    print(f"  ECM at survivors d28: {_fmt(treat['ecm_at_survivors'])}")
    print(f"  ECM at tumor d14:     {_fmt(pre['ecm_at_survivors'])}")
    print(f"  Tumor radius p95 d14/d28/d42: {_fmt(pre['tumor_radius_p95'],2)} / {_fmt(treat['tumor_radius_p95'],2)} / {_fmt(post['tumor_radius_p95'],2)}")
    print(f"  Sanctuary zones: {'YES' if c4 else 'NO'}")
    print(f"  Barrier intact: {'YES' if c3 else 'NO'}")

    # 8) Final verdict
    hard = [crit["RC2-1"][0], crit["RC2-2"][0], crit["RC2-3"][0], crit["RC2-4"][0], crit["RC2-6"][0]]
    hard_score = sum(hard)

    suspicious = []
    if crit["RC2-3"][0] and math.isfinite(pre["peri_ecm"]) and abs(pre["peri_ecm"]) < 1e-9:
        suspicious.append("Barrier pass may be loophole: peri-ECM pre-treatment is ~0.")
    if crit["RC2-4"][0] and math.isfinite(pre["ecm_at_survivors"]) and math.isfinite(treat["ecm_at_survivors"]) and abs(treat["ecm_at_survivors"] - pre["ecm_at_survivors"]) < 1e-6:
        suspicious.append("Sanctuary pass is numerically tiny; verify true spatial separation.")
    if all(math.isfinite(v) and abs(v) < 1e-8 for v in intra):
        suspicious.append("Intracellular drug ~0 throughout d28-31 survivors (possible survivorship filtering).")

    true_bio_pass = (hard_score == 5) and (len(suspicious) == 0) and pref_surv

    print("\n8) Final Verdict")
    print(f"  Overall RC2 hard score: {hard_score}/5")
    print(f"  Biological validity: {'TRUE BIOLOGICAL PASS' if true_bio_pass else 'CONDITIONAL / POSSIBLE FALSE PASS'}")
    if hard_score == 5:
        print(f"  Success mechanism: treatment debulking + sanctuary-preserved survivors + post-withdrawal regrowth.")
    else:
        print(f"  Failure mechanism: one or more hard RC2 gates not satisfied.")
    if suspicious:
        print("  Risks / flags:")
        for s in suspicious:
            print(f"    - {s}")
    else:
        print("  Risks / flags: none material.")

    return hard_score, true_bio_pass


def _jobs_from_manifest(manifest_path: Path):
    m = json.loads(manifest_path.read_text())
    jobs = []
    for name, v in m.get("variants", {}).items():
        out = Path(v["work_dir"]) / "output"
        timing = v.get("timing", {
            "drug_start_time": 20160.0,
            "drug_end_time": 40320.0,
            "max_time": 60480.0,
        })
        jobs.append((name, out, timing))
    return jobs


def main():
    ap = argparse.ArgumentParser(description="Mechanistic RC2 evaluator")
    ap.add_argument("--out-dir", type=Path, default=None,
                    help="Single RC2 output dir (contains output*.xml)")
    ap.add_argument("--manifest", type=Path, default=None,
                    help="Manifest JSON from a sweep launcher")
    ap.add_argument("--all-known", action="store_true",
                    help="Evaluate baseline + wave1 + wave2 + wave3 manifests")
    args = ap.parse_args()

    jobs = []

    if args.all_known:
        jobs.append((
            "baseline_10620",
            (PROJECT_ROOT / "build" / "rc2_full_seed42" / "replicate_01_seed42" / "output").resolve(),
            {"drug_start_time": 20160.0, "drug_end_time": 40320.0, "max_time": 60480.0},
        ))
        for mp in [
            Path("/work/a0hirs04/PROJECT-NORTHSTAR/build/rc2_countermeasures/manifest.json"),
            Path("/work/a0hirs04/PROJECT-NORTHSTAR/build/rc2_wave2/manifest.json"),
            Path("/work/a0hirs04/PROJECT-NORTHSTAR/build/rc2_wave3/manifest.json"),
        ]:
            if mp.exists():
                jobs.extend(_jobs_from_manifest(mp))
    elif args.manifest is not None:
        jobs = _jobs_from_manifest(args.manifest.resolve())
    else:
        out_dir = args.out_dir.resolve() if args.out_dir is not None else DEFAULT_OUT_DIR.resolve()
        jobs = [(
            f"seed_{_infer_seed_label(out_dir)}",
            out_dir,
            {"drug_start_time": 20160.0, "drug_end_time": 40320.0, "max_time": 60480.0},
        )]

    if not jobs:
        print("No jobs found to evaluate.")
        return 1

    summary = []
    for name, out_dir, timing in jobs:
        if not out_dir.exists():
            print(f"\n[SKIP] {name}: output dir missing: {out_dir}")
            continue
        try:
            hard_score, true_bio = _report_one(name, out_dir, timing)
            summary.append((name, hard_score, true_bio))
        except Exception as e:
            print(f"\n[ERROR] {name}: {e}")

    print("\n" + "=" * 100)
    print("SUMMARY")
    print("=" * 100)
    for name, hard, tb in summary:
        print(f"  {name:<28} hard={hard}/5   biological={'TRUE' if tb else 'CONDITIONAL/FALSE'}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
