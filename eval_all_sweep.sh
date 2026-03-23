#!/usr/bin/env bash
# Evaluate all 10 sweep variants using official evaluators
set -uo pipefail

SWEEP="/work/a0hirs04/PROJECT-NORTHSTAR/build/sweep"
PROJECT="/home/a0hirs04/PROJECT-NORTHSTAR"

echo "================================================================"
echo "  FULL SWEEP EVALUATION — 10 VARIANTS × (RC1 + RC2)"
echo "================================================================"

for v in v01 v02 v03 v04 v05 v06 v07 v08 v09 v10; do
    echo ""
    echo "================================================================"
    echo "  VARIANT: $v"
    echo "================================================================"

    # Get parameters from config
    uptake=$(grep -oP 'drug_uptake_rate[^>]*>\K[^<]+' "$SWEEP/$v/rc2/config.xml" 2>/dev/null || echo "?")
    kill=$(grep -oP 'drug_kill_coefficient[^>]*>\K[^<]+' "$SWEEP/$v/rc2/config.xml" 2>/dev/null || echo "?")
    hif1a=$(grep -oP 'hif1a_emt_boost[^>]*>\K[^<]+' "$SWEEP/$v/rc1/replicate_01_seed42/config.xml" 2>/dev/null || echo "?")
    echo "  Params: uptake=$uptake  kill=$kill  hif1a=$hif1a"

    echo ""
    echo "  --- RC1 ---"
    python3 "$PROJECT/run_reality_check_1.py" \
        --work-dir "$SWEEP/$v/rc1" \
        --seeds 42 --quorum 1 --evaluate-only 2>&1 | \
        grep -E "^\s+\[(PASS|FAIL)\]|REALITY CHECK 1"

    echo ""
    echo "  --- RC2 ---"
    python3 "$PROJECT/evaluate_rc2.py" \
        --out-dir "$SWEEP/$v/rc2/output" 2>&1 | \
        grep -E "^\s+\[(PASS|FAIL|INFO)\]|RC2 \(seed|HARD CRITERIA|ABCB1 resistance"
done

echo ""
echo "================================================================"
echo "  DONE"
echo "================================================================"
