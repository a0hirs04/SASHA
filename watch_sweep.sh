#!/usr/bin/env bash
# watch_sweep.sh — Monitor 10-variant RC1+RC2 parameter sweep
# Usage: watch -n 30 bash watch_sweep.sh

set -uo pipefail

SWEEP_DIR="/work/a0hirs04/PROJECT-NORTHSTAR/build/sweep"
MANIFEST="$SWEEP_DIR/jobs.txt"
PROJECT="/home/a0hirs04/PROJECT-NORTHSTAR"

RC1_TOTAL=85   # 21 days @ 360 min intervals + initial
RC2_TOTAL=169  # 42 days @ 360 min intervals + initial

VARIANTS=(v01 v02 v03 v04 v05 v06 v07 v08 v09 v10)

# Parameter labels for display
declare -A PARAMS
PARAMS[v01]="uptake=0.05 kill=0.02 hif=0.05"
PARAMS[v02]="uptake=0.10 kill=0.02 hif=0.05"
PARAMS[v03]="uptake=0.20 kill=0.02 hif=0.05"
PARAMS[v04]="uptake=0.05 kill=0.05 hif=0.05"
PARAMS[v05]="uptake=0.10 kill=0.05 hif=0.05"
PARAMS[v06]="uptake=0.05 kill=0.10 hif=0.05"
PARAMS[v07]="uptake=0.10 kill=0.10 hif=0.05"
PARAMS[v08]="uptake=0.10 kill=0.02 hif=0.02"
PARAMS[v09]="uptake=0.10 kill=0.05 hif=0.02"
PARAMS[v10]="uptake=0.10 kill=0.02 hif=0.00"

count_snaps() {
    ls "$1"/output*.xml 2>/dev/null | wc -l
}

job_state() {
    local jid=$1
    local st
    st=$(squeue -j "$jid" -h -o "%T" 2>/dev/null)
    if [[ -z "$st" ]]; then
        st=$(sacct -j "$jid" -n -o State -X 2>/dev/null | head -1 | tr -d ' ')
    fi
    echo "${st:-UNKNOWN}"
}

pct_bar() {
    local n=$1 total=$2
    local pct=$(( n * 100 / (total > 0 ? total : 1) ))
    if (( pct >= 100 )); then echo "DONE"
    elif (( pct >= 50 )); then echo "${pct}%"
    else echo "${pct}%"
    fi
}

echo "=================================================================="
echo "  10-VARIANT SWEEP — $(date '+%H:%M:%S %Z')"
echo "=================================================================="
echo ""
printf "  %-4s  %-35s  %-12s  %-12s\n" "V#" "Parameters" "RC1" "RC2"
printf "  %-4s  %-35s  %-12s  %-12s\n" "----" "-----------------------------------" "------------" "------------"

all_done=true
for v in "${VARIANTS[@]}"; do
    rc1_dir="$SWEEP_DIR/$v/rc1/replicate_01_seed42"
    rc2_dir="$SWEEP_DIR/$v/rc2"

    # Count snapshots
    rc1_snaps=0; rc2_snaps=0
    [[ -d "$rc1_dir/output" ]] && rc1_snaps=$(count_snaps "$rc1_dir/output")
    [[ -d "$rc2_dir/output" ]] && rc2_snaps=$(count_snaps "$rc2_dir/output")

    # Job states from manifest
    rc1_jid=$(grep "^$v	rc1	" "$MANIFEST" 2>/dev/null | cut -f3)
    rc2_jid=$(grep "^$v	rc2	" "$MANIFEST" 2>/dev/null | cut -f3)

    rc1_st="?"; rc2_st="?"
    [[ -n "$rc1_jid" ]] && rc1_st=$(job_state "$rc1_jid")
    [[ -n "$rc2_jid" ]] && rc2_st=$(job_state "$rc2_jid")

    [[ "$rc1_st" != "COMPLETED" || "$rc2_st" != "COMPLETED" ]] && all_done=false

    rc1_info="${rc1_snaps}/${RC1_TOTAL} $(pct_bar $rc1_snaps $RC1_TOTAL)"
    rc2_info="${rc2_snaps}/${RC2_TOTAL} $(pct_bar $rc2_snaps $RC2_TOTAL)"

    printf "  %-4s  %-35s  %-12s  %-12s\n" "$v" "${PARAMS[$v]}" "$rc1_info" "$rc2_info"
done

# If all done, run quick evaluation
if $all_done; then
    echo ""
    echo "  ALL JOBS COMPLETE — EVALUATING..."
    echo ""
    printf "  %-4s  %-8s  %-8s  %-8s  %-8s  %-8s  %-8s  %-8s  %-8s  %-8s  %-8s  %-8s  %-8s  %-8s\n" \
        "V#" "RC1-C1" "RC1-C2" "RC1-C3" "RC1-C4" "RC1-C5" "RC1-C6" "RC1-C7" "RC1-C8" "RC2-C1" "RC2-C2" "RC2-C3" "RC2-C4" "RC2-C6"
    printf "  %-4s  %-8s  %-8s  %-8s  %-8s  %-8s  %-8s  %-8s  %-8s  %-8s  %-8s  %-8s  %-8s  %-8s\n" \
        "----" "--------" "--------" "--------" "--------" "--------" "--------" "--------" "--------" "--------" "--------" "--------" "--------" "--------"

    for v in "${VARIANTS[@]}"; do
        rc1_work="$SWEEP_DIR/$v/rc1"
        rc2_out="$SWEEP_DIR/$v/rc2/output"

        # Run evaluations and capture results
        rc1_result=$(cd "$PROJECT" && python3 evaluate_rc1.py --work-dir "$rc1_work" --seeds 42 --quorum 1 2>&1 || true)
        rc2_result=$(cd "$PROJECT" && python3 evaluate_rc2.py --out-dir "$rc2_out" 2>&1 || true)

        # Parse RC1 criteria (8 criteria)
        rc1_c=()
        for i in 1 2 3 4 5 6 7 8; do
            if echo "$rc1_result" | grep -q "\[PASS\] $i\."; then
                rc1_c+=("PASS")
            else
                rc1_c+=("FAIL")
            fi
        done

        # Parse RC2 criteria
        rc2_c=()
        for label in "RC2-1" "RC2-2" "RC2-3" "RC2-4" "RC2-6"; do
            if echo "$rc2_result" | grep -q "\[PASS\] $label"; then
                rc2_c+=("PASS")
            else
                rc2_c+=("FAIL")
            fi
        done

        # Color coding via markers
        line=$(printf "  %-4s" "$v")
        for c in "${rc1_c[@]}"; do
            if [[ "$c" == "PASS" ]]; then
                line+=$(printf "  %-8s" "PASS")
            else
                line+=$(printf "  *%-7s" "FAIL")
            fi
        done
        for c in "${rc2_c[@]}"; do
            if [[ "$c" == "PASS" ]]; then
                line+=$(printf "  %-8s" "PASS")
            else
                line+=$(printf "  *%-7s" "FAIL")
            fi
        done
        echo "$line"
    done
fi

echo ""
echo "  Evaluate manually:  python3 evaluate_rc1.py --work-dir $SWEEP_DIR/vNN/rc1 --seeds 42 --quorum 1"
echo "                      python3 evaluate_rc2.py --out-dir $SWEEP_DIR/vNN/rc2/output"
echo "=================================================================="
