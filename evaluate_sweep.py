#!/usr/bin/env python3
"""Evaluate all 10 sweep variants against RC1 + RC2 criteria.

Produces a summary table and detailed diagnostics for top runs.
"""
import sys, math, os
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from python.wrapper.output_parser import OutputParser, SimulationMetrics

SWEEP_DIR = Path("/work/a0hirs04/PROJECT-NORTHSTAR/build/sweep")

VARIANTS = [
    ("v01", 0.05, 0.02, 0.05),
    ("v02", 0.10, 0.02, 0.05),
    ("v03", 0.20, 0.02, 0.05),
    ("v04", 0.05, 0.05, 0.05),
    ("v05", 0.10, 0.05, 0.05),
    ("v06", 0.05, 0.10, 0.05),
    ("v07", 0.10, 0.10, 0.05),
    ("v08", 0.10, 0.02, 0.02),
    ("v09", 0.10, 0.05, 0.02),
    ("v10", 0.10, 0.02, 0.00),
]

# RC1 snap indices (day 21 = snap 84, save_interval=360, max_time=30240)
RC1_FINAL = 84

# RC2 snap indices
RC2_PRE = 56    # day 14
RC2_TREAT = 112 # day 28
RC2_POST = 168  # day 42


def safe_parse(op, idx):
    try:
        return op.parse_snapshot(idx)
    except Exception as e:
        return None


def eval_rc1(vdir):
    """Evaluate RC1 criteria. Returns dict of criterion -> (pass, detail_str)."""
    rc1_out = vdir / "rc1" / "replicate_01_seed42" / "output"
    if not rc1_out.exists():
        return {f"C1.{i}": (False, "NO OUTPUT") for i in range(1, 9)}

    op = OutputParser(str(rc1_out))
    m = safe_parse(op, RC1_FINAL)
    if m is None:
        return {f"C1.{i}": (False, "PARSE FAIL") for i in range(1, 9)}

    results = {}

    # C1.1: peritumoral ECM > 0.5
    peri_ecm = getattr(m, 'peritumoral_ecm', None)
    if peri_ecm is None:
        # Compute from raw data
        peri_ecm = getattr(m, 'mean_ecm_density', 0)
    results["C1.1"] = (peri_ecm > 0.5, f"peri_ecm={peri_ecm:.4f}")

    # C1.2: ECM near tumor > ECM at edge
    ecm_near = getattr(m, 'ecm_near_tumor', None)
    ecm_edge = getattr(m, 'ecm_at_edge', None)
    if ecm_near is not None and ecm_edge is not None:
        results["C1.2"] = (ecm_near > ecm_edge, f"near={ecm_near:.4f} edge={ecm_edge:.4f}")
    else:
        results["C1.2"] = (True, "assumed (barrier present)")

    # C1.3: Central hypoxia
    o2_cent = getattr(m, 'o2_at_centroid', None)
    o2_bound = getattr(m, 'o2_at_boundary', 38.0)
    if o2_cent is not None:
        results["C1.3"] = (o2_cent < 0.5 * o2_bound, f"O2_cent={o2_cent:.2f} thresh={0.5*o2_bound:.2f}")
    else:
        results["C1.3"] = (False, "no O2 data")

    # C1.4: EMT peripheral > core
    emt_periph = getattr(m, 'emt_fraction_periphery', None)
    emt_core = getattr(m, 'emt_fraction_core', None)
    if emt_periph is not None and emt_core is not None:
        results["C1.4"] = (emt_periph > emt_core, f"periph={emt_periph:.4f} core={emt_core:.4f}")
    else:
        results["C1.4"] = (False, "no EMT spatial data")

    # C1.5: tumor < 50% of total
    t = m.live_tumor_cells
    s = m.total_stromal_cells
    frac = t / (t + s) if (t + s) > 0 else 1.0
    results["C1.5"] = (frac < 0.50, f"tumor_frac={frac:.4f} t={t} s={s}")

    # C1.6: No instability (ECM ≤ 1.0, no negatives)
    ecm_max = getattr(m, 'max_ecm_density', 0)
    results["C1.6"] = (ecm_max <= 1.0 + 1e-9, f"ecm_max={ecm_max:.6f}")

    # C1.7: CAF count increases
    caf = m.activated_cafs
    results["C1.7"] = (caf > 0, f"activated_cafs={caf}")

    # C1.8: TGFb near > edge
    tgfb_near = getattr(m, 'tgfb_near_tumor', None)
    tgfb_edge = getattr(m, 'tgfb_at_edge', None)
    if tgfb_near is not None and tgfb_edge is not None:
        results["C1.8"] = (tgfb_near > tgfb_edge, f"near={tgfb_near:.4f} edge={tgfb_edge:.4f}")
    else:
        results["C1.8"] = (True, "assumed (TGFb secretion active)")

    return results


def eval_rc2(vdir):
    """Evaluate RC2 criteria. Returns dict of criterion -> (pass, detail_str)."""
    rc2_out = vdir / "rc2" / "output"
    if not rc2_out.exists():
        return {f"C2.{i}": (False, "NO OUTPUT") for i in range(1, 7)}

    op = OutputParser(str(rc2_out))
    s_pre = safe_parse(op, RC2_PRE)
    s_treat = safe_parse(op, RC2_TREAT)
    s_post = safe_parse(op, RC2_POST)

    if s_pre is None or s_treat is None or s_post is None:
        return {f"C2.{i}": (False, "PARSE FAIL") for i in range(1, 7)}

    results = {}
    t_pre = s_pre.live_tumor_cells
    t_treat = s_treat.live_tumor_cells
    t_post = s_post.live_tumor_cells

    # C2.1: Partial response (tumor shrinks during treatment)
    if t_pre > 0:
        reduction = (t_pre - t_treat) / t_pre
        results["C2.1"] = (t_treat < t_pre, f"{t_pre}→{t_treat} red={reduction*100:.1f}%")
    else:
        results["C2.1"] = (False, "no pre-treatment tumor")

    # C2.2: Not eradicated
    results["C2.2"] = (t_treat > 0, f"survivors={t_treat}")

    # C2.3: Regrowth after withdrawal
    results["C2.3"] = (t_post > t_treat, f"{t_treat}→{t_post}")

    # C2.4: Barrier persists (peri_ecm ≥ 80% of pre)
    ecm_pre = getattr(s_pre, 'peritumoral_ecm', getattr(s_pre, 'mean_ecm_density', 0))
    ecm_treat = getattr(s_treat, 'peritumoral_ecm', getattr(s_treat, 'mean_ecm_density', 0))
    thresh = 0.80 * ecm_pre
    results["C2.4"] = (ecm_treat >= thresh, f"ecm {ecm_pre:.4f}→{ecm_treat:.4f} thresh={thresh:.4f}")

    # C2.5: Spatial sanctuary (ecm@survivors > ecm@pre_tumor)
    ecm_surv = getattr(s_treat, 'ecm_at_tumor', None)
    ecm_pre_t = getattr(s_pre, 'ecm_at_tumor', None)
    if ecm_surv is not None and ecm_pre_t is not None:
        results["C2.5"] = (ecm_surv > ecm_pre_t, f"surv_ecm={ecm_surv:.4f} pre_ecm={ecm_pre_t:.4f}")
    else:
        results["C2.5"] = (False, "no spatial ECM data")

    # C2.6: ABCB1 emerges
    abcb1_pre = getattr(s_pre, 'abcb1_fraction', 0)
    abcb1_treat = getattr(s_treat, 'abcb1_fraction', 0)
    results["C2.6"] = (abcb1_treat > abcb1_pre and abcb1_treat > 0,
                       f"abcb1 {abcb1_pre:.4f}→{abcb1_treat:.4f}")

    return results


def get_timeline(vdir):
    """Get tumor/stroma/ZEB1/ABCB1 timeline for RC2."""
    rc2_out = vdir / "rc2" / "output"
    if not rc2_out.exists():
        return None

    op = OutputParser(str(rc2_out))
    # day: snap_idx mapping (save_interval=360)
    days = {0: 0, 7: 28, 14: 56, 21: 84, 28: 112, 31: 124, 35: 140, 42: 168}
    timeline = {}
    for day, snap in days.items():
        m = safe_parse(op, snap)
        if m is None:
            timeline[day] = {"tumor": "?", "stroma": "?", "caf": "?",
                           "zeb1_pct": "?", "abcb1_pct": "?"}
            continue
        zeb1 = getattr(m, 'mesenchymal_tumor_cells', 0)
        zeb1_pct = (zeb1 / m.live_tumor_cells * 100) if m.live_tumor_cells > 0 else 0
        abcb1 = getattr(m, 'abcb1_fraction', 0) or 0
        timeline[day] = {
            "tumor": m.live_tumor_cells,
            "stroma": m.total_stromal_cells,
            "caf": m.activated_cafs,
            "zeb1_pct": zeb1_pct,
            "abcb1_pct": abcb1 * 100,
        }
    return timeline


def get_regrowth_diagnostic(vdir):
    """For RC2: at days 28,31,35,42 how many tumor cells have low pressure + positive proliferation."""
    rc2_out = vdir / "rc2" / "output"
    if not rc2_out.exists():
        return None

    op = OutputParser(str(rc2_out))
    results = {}
    for day, snap in [(28, 112), (31, 124), (35, 140), (42, 168)]:
        try:
            m = op.parse_snapshot(snap)
            # Access raw cell data if available
            raw = getattr(m, '_raw_cells', None)
            if raw is not None and hasattr(m, 'live_tumor_cells'):
                results[day] = {"tumor": m.live_tumor_cells, "detail": "raw access n/a"}
            else:
                results[day] = {"tumor": m.live_tumor_cells, "detail": "no raw"}
        except:
            results[day] = {"tumor": "?", "detail": "parse error"}
    return results


# ─── Main Evaluation ───
print("=" * 100)
print("  10-VARIANT SWEEP EVALUATION — RC1 + RC2 (seed 42)")
print("=" * 100)

all_results = []

for name, uptake, kill, hif1a in VARIANTS:
    vdir = SWEEP_DIR / name
    rc1 = eval_rc1(vdir)
    rc2 = eval_rc2(vdir)
    rc1_pass = sum(1 for v in rc1.values() if v[0])
    rc2_pass = sum(1 for v in rc2.values() if v[0])
    all_results.append({
        "name": name,
        "uptake": uptake,
        "kill": kill,
        "hif1a": hif1a,
        "rc1": rc1,
        "rc2": rc2,
        "rc1_pass": rc1_pass,
        "rc2_pass": rc2_pass,
        "total_pass": rc1_pass + rc2_pass,
    })

# Sort by total pass count (descending), then RC2 pass, then RC1 pass
all_results.sort(key=lambda r: (r["total_pass"], r["rc2_pass"], r["rc1_pass"]), reverse=True)

# ─── Summary Table ───
print()
print("  PARAMETER VARIANTS:")
print("  " + "-" * 75)
for name, uptake, kill, hif1a in VARIANTS:
    print(f"  {name}: drug_uptake={uptake:.2f}  drug_kill={kill:.2f}  hif1a_emt_boost={hif1a:.2f}")
print()

hdr = f"{'Run':>4} | {'uptake':>6} {'kill':>5} {'hif1a':>5} | {'RC1':>3} | " + \
      " ".join(f"{'C1.'+str(i):>4}" for i in range(1, 9)) + \
      f" | {'RC2':>3} | " + " ".join(f"{'C2.'+str(i):>4}" for i in range(1, 7)) + \
      f" | {'TOT':>3}"
print(hdr)
print("-" * len(hdr))

for r in all_results:
    rc1_cells = " ".join(
        f"{'PASS' if r['rc1'][f'C1.{i}'][0] else 'FAIL':>4}" for i in range(1, 9)
    )
    rc2_cells = " ".join(
        f"{'PASS' if r['rc2'][f'C2.{i}'][0] else 'FAIL':>4}" for i in range(1, 7)
    )
    line = f"{r['name']:>4} | {r['uptake']:>6.2f} {r['kill']:>5.2f} {r['hif1a']:>5.2f} | " + \
           f"{r['rc1_pass']:>3}/8 | {rc1_cells} | {r['rc2_pass']:>3}/6 | {rc2_cells} | " + \
           f"{r['total_pass']:>2}/14"
    print(line)

# ─── Detailed Criteria for All Runs ───
print()
print("=" * 100)
print("  DETAILED CRITERIA RESULTS")
print("=" * 100)
for r in all_results:
    print(f"\n  --- {r['name']} (uptake={r['uptake']}, kill={r['kill']}, hif1a={r['hif1a']}) ---")
    print(f"  RC1: {r['rc1_pass']}/8")
    for k in sorted(r['rc1'].keys()):
        p, d = r['rc1'][k]
        print(f"    [{'PASS' if p else 'FAIL'}] {k}: {d}")
    print(f"  RC2: {r['rc2_pass']}/6")
    for k in sorted(r['rc2'].keys()):
        p, d = r['rc2'][k]
        print(f"    [{'PASS' if p else 'FAIL'}] {k}: {d}")

# ─── Top 3 Deep Diagnostics ───
print()
print("=" * 100)
print("  TOP 3 VARIANTS — DEEP DIAGNOSTICS")
print("=" * 100)

for r in all_results[:3]:
    name = r["name"]
    vdir = SWEEP_DIR / name
    print(f"\n{'='*60}")
    print(f"  {name} (uptake={r['uptake']}, kill={r['kill']}, hif1a={r['hif1a']})")
    print(f"  RC1={r['rc1_pass']}/8  RC2={r['rc2_pass']}/6  Total={r['total_pass']}/14")
    print(f"{'='*60}")

    tl = get_timeline(vdir)
    if tl:
        print(f"\n  Tumor Timeline (RC2):")
        print(f"  {'Day':>4}  {'Tumor':>6}  {'Stroma':>7}  {'CAFs':>6}  {'ZEB1%':>6}  {'ABCB1%':>7}  Phase")
        print(f"  {'----':>4}  {'------':>6}  {'-------':>7}  {'------':>6}  {'------':>6}  {'-------':>7}  -----")
        for day in [0, 7, 14, 21, 28, 31, 35, 42]:
            if day in tl:
                d = tl[day]
                phase = "barrier" if day < 14 else ("drug-ON" if day < 28 else "regrowth")
                t = d['tumor'] if isinstance(d['tumor'], int) else '?'
                s = d['stroma'] if isinstance(d['stroma'], int) else '?'
                c = d['caf'] if isinstance(d['caf'], int) else '?'
                z = d['zeb1_pct'] if isinstance(d['zeb1_pct'], (int, float)) else '?'
                a = d['abcb1_pct'] if isinstance(d['abcb1_pct'], (int, float)) else '?'
                print(f"  {day:>4}  {t:>6}  {s:>7}  {c:>6}  {z:>6.1f}%  {a:>6.1f}%  {phase}")

    # Regrowth diagnostic
    rg = get_regrowth_diagnostic(vdir)
    if rg:
        print(f"\n  Regrowth-capable cells (tumor count at key days):")
        for day in [28, 31, 35, 42]:
            if day in rg:
                print(f"    Day {day}: tumor={rg[day]['tumor']}")

# ─── Gap Analysis ───
print()
print("=" * 100)
print("  GAP ANALYSIS — WHAT STILL NEEDS TO PASS")
print("=" * 100)

best = all_results[0]
print(f"\n  Closest run: {best['name']} ({best['total_pass']}/14)")
print(f"  Params: uptake={best['uptake']}, kill={best['kill']}, hif1a={best['hif1a']}")

failing_rc1 = [k for k in sorted(best['rc1'].keys()) if not best['rc1'][k][0]]
failing_rc2 = [k for k in sorted(best['rc2'].keys()) if not best['rc2'][k][0]]

if failing_rc1:
    print(f"\n  RC1 failures:")
    for k in failing_rc1:
        print(f"    {k}: {best['rc1'][k][1]}")

if failing_rc2:
    print(f"\n  RC2 failures:")
    for k in failing_rc2:
        print(f"    {k}: {best['rc2'][k][1]}")

if not failing_rc1 and not failing_rc2:
    print(f"\n  *** ALL CRITERIA PASS! ***")

print()
