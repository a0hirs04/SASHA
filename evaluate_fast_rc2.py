#!/usr/bin/env python3
"""Evaluate fast_rc2 screening results for fixL variants.

Reads the fast_screen.csv written by the C++ fast-screen logger and applies
pass/fail criteria as a proxy for RC2-6 (post-withdrawal regrowth).

Pass criteria (within 2000 min after drug_end_time):
  - ZEB1 fraction decreases by >=10%  (EMT reversion happening)
  - AND mean cycle rate increases      (proliferation recovering)
  - AND tumor count stabilizes or increases (not still declining)

Fail criteria:
  - ZEB1 fraction stays flat or increases
  - AND cycle rate remains suppressed

Usage:
  python3 evaluate_fast_rc2.py [--base-dir DIR]
"""
import argparse
import csv
import sys
from pathlib import Path

DEFAULT_BASE = Path("/work/a0hirs04/PROJECT-NORTHSTAR/build/fixL")
VARIANTS = ["a", "b", "c", "d"]
T_TREAT_END = 40320.0
EVAL_WINDOW = 2000.0  # minutes after drug_end to evaluate


def load_csv(csv_path):
    """Load fast_screen.csv and return list of row dicts."""
    rows = []
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append({
                "time": float(row["time"]),
                "tumor_count": int(row["tumor_count"]),
                "zeb1_count": int(row["zeb1_count"]),
                "zeb1_fraction": float(row["zeb1_fraction"]),
                "mean_cycle_rate": float(row["mean_cycle_rate"]),
                "emt_to_epi": int(row["emt_to_epi"]),
            })
    return rows


def evaluate_variant(base_dir, vname):
    """Evaluate a single variant's fast_rc2 results."""
    out_dir = base_dir / vname / "fast_rc2" / "output"
    csv_path = out_dir / "fast_screen.csv"
    fail_marker = out_dir / "FAST_SCREEN_FAIL.txt"

    result = {
        "variant": vname,
        "status": "UNKNOWN",
        "reason": "",
        "baseline": {},
        "final": {},
        "zeb1_change_pct": None,
        "cycle_change_pct": None,
        "count_change_pct": None,
    }

    # Check for early termination marker
    if fail_marker.exists():
        result["status"] = "FAIL"
        result["reason"] = "Early terminated by C++ fast-screen"
        with open(fail_marker) as f:
            result["reason"] += f" — {f.read().strip()}"
        # Still try to load CSV for reporting
        if not csv_path.exists():
            return result

    if not csv_path.exists():
        result["status"] = "NO_DATA"
        result["reason"] = f"No fast_screen.csv in {out_dir}"
        return result

    rows = load_csv(csv_path)
    if len(rows) < 5:
        result["status"] = "INCOMPLETE"
        result["reason"] = f"Only {len(rows)} data points (need >=5)"
        return result

    # Baseline: first row (just after drug withdrawal)
    baseline = rows[0]
    # Final: last row
    final = rows[-1]

    result["baseline"] = baseline
    result["final"] = final

    # Calculate changes
    if baseline["zeb1_fraction"] > 0:
        result["zeb1_change_pct"] = (
            (final["zeb1_fraction"] - baseline["zeb1_fraction"])
            / baseline["zeb1_fraction"] * 100
        )
    else:
        result["zeb1_change_pct"] = 0.0

    if baseline["mean_cycle_rate"] > 0:
        result["cycle_change_pct"] = (
            (final["mean_cycle_rate"] - baseline["mean_cycle_rate"])
            / baseline["mean_cycle_rate"] * 100
        )
    else:
        result["cycle_change_pct"] = (
            100.0 if final["mean_cycle_rate"] > 0 else 0.0
        )

    if baseline["tumor_count"] > 0:
        result["count_change_pct"] = (
            (final["tumor_count"] - baseline["tumor_count"])
            / baseline["tumor_count"] * 100
        )
    else:
        result["count_change_pct"] = 0.0

    # Already marked FAIL by early termination
    if result["status"] == "FAIL":
        return result

    # Apply pass/fail criteria
    zeb1_decreased = result["zeb1_change_pct"] <= -10.0
    cycle_increased = result["cycle_change_pct"] > 0.0
    count_stable = final["tumor_count"] >= baseline["tumor_count"]

    zeb1_stuck = result["zeb1_change_pct"] > -1.0
    cycle_suppressed = result["cycle_change_pct"] <= 5.0
    count_declining = final["tumor_count"] < baseline["tumor_count"]

    if zeb1_decreased and cycle_increased and count_stable:
        result["status"] = "PASS"
        result["reason"] = (
            f"ZEB1 {result['zeb1_change_pct']:+.1f}%, "
            f"cycle {result['cycle_change_pct']:+.1f}%, "
            f"count {baseline['tumor_count']}→{final['tumor_count']}"
        )
    elif zeb1_stuck and cycle_suppressed and count_declining:
        result["status"] = "FAIL"
        result["reason"] = (
            f"ZEB1 {result['zeb1_change_pct']:+.1f}% (stuck), "
            f"cycle {result['cycle_change_pct']:+.1f}% (suppressed), "
            f"count {baseline['tumor_count']}→{final['tumor_count']} (declining)"
        )
    elif zeb1_decreased or cycle_increased:
        result["status"] = "MARGINAL"
        result["reason"] = (
            f"ZEB1 {result['zeb1_change_pct']:+.1f}%, "
            f"cycle {result['cycle_change_pct']:+.1f}%, "
            f"count {baseline['tumor_count']}→{final['tumor_count']} "
            "(partial recovery signals)"
        )
    else:
        result["status"] = "FAIL"
        result["reason"] = (
            f"ZEB1 {result['zeb1_change_pct']:+.1f}%, "
            f"cycle {result['cycle_change_pct']:+.1f}%, "
            f"count {baseline['tumor_count']}→{final['tumor_count']}"
        )

    return result


def print_timeseries(base_dir, vname):
    """Print condensed timeseries for a variant."""
    csv_path = base_dir / vname / "fast_rc2" / "output" / "fast_screen.csv"
    if not csv_path.exists():
        return
    rows = load_csv(csv_path)
    if not rows:
        return

    print(f"\n  Timeseries for {vname}:")
    print(f"  {'time':>8} {'tumor':>7} {'ZEB1%':>7} {'cycle':>10} {'E→epi':>6}")
    print(f"  {'-'*42}")
    # Show every 5th row + first + last
    indices = [0] + list(range(4, len(rows), 5))
    if len(rows) - 1 not in indices:
        indices.append(len(rows) - 1)
    for i in sorted(set(indices)):
        r = rows[i]
        print(f"  {r['time']:>8.0f} {r['tumor_count']:>7d} "
              f"{r['zeb1_fraction']*100:>6.1f}% {r['mean_cycle_rate']:>10.6f} "
              f"{r['emt_to_epi']:>6d}")


def main():
    parser = argparse.ArgumentParser(
        description="Evaluate fixL fast_rc2 screening results")
    parser.add_argument("--base-dir", type=Path, default=DEFAULT_BASE)
    parser.add_argument("--verbose", "-v", action="store_true",
                        help="Show per-variant timeseries")
    args = parser.parse_args()

    print("=" * 68)
    print("fixL Fast RC2 Screen Evaluation")
    print(f"Base: {args.base_dir}")
    print("=" * 68)

    results = []
    for vname in VARIANTS:
        r = evaluate_variant(args.base_dir, vname)
        results.append(r)

    # Summary table
    print(f"\n{'Var':<4} {'Status':<10} {'ZEB1%':>8} {'Cycle%':>8} "
          f"{'Count':>12} {'Reason'}")
    print("-" * 72)
    for r in results:
        zc = f"{r['zeb1_change_pct']:+.1f}" if r["zeb1_change_pct"] is not None else "-"
        cc = f"{r['cycle_change_pct']:+.1f}" if r["cycle_change_pct"] is not None else "-"
        if r["baseline"] and r["final"]:
            cnt = f"{r['baseline']['tumor_count']}→{r['final']['tumor_count']}"
        else:
            cnt = "-"
        print(f"{r['variant']:<4} {r['status']:<10} {zc:>8} {cc:>8} "
              f"{cnt:>12} {r['reason']}")

    if args.verbose:
        for r in results:
            print_timeseries(args.base_dir, r["variant"])

    # Summary
    passed = [r for r in results if r["status"] == "PASS"]
    marginal = [r for r in results if r["status"] == "MARGINAL"]
    failed = [r for r in results if r["status"] == "FAIL"]

    print(f"\nPASS: {len(passed)}  MARGINAL: {len(marginal)}  FAIL: {len(failed)}")

    if passed:
        print(f"\nPromote to full RC2: {', '.join(r['variant'] for r in passed)}")
        print(f"  python3 launch_fixL.py --promote")
    elif marginal:
        print(f"\nMarginal variants (consider promoting): "
              f"{', '.join(r['variant'] for r in marginal)}")
    else:
        print("\nAll variants failed fast screening. Iterate on parameters.")

    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
