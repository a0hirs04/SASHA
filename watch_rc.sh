#!/usr/bin/env bash
# watch_rc.sh — Combined RC1 + RC2 live monitor
# Usage:  watch -n 30 bash watch_rc.sh

set -uo pipefail

RC1_DIR="/work/a0hirs04/PROJECT-NORTHSTAR/build/rc1_final"
RC2_DIR="/work/a0hirs04/PROJECT-NORTHSTAR/build/rc2_final"
SEEDS=(42 99 137 256 1001)

# RC1: 21 days = 30240 min, interval 360 -> 84 snapshots + initial = 85
RC1_TOTAL=85
# RC2: 42 days = 60480 min, interval 360 -> 168 snapshots + initial = 169
RC2_TOTAL=169

bar() {
    local n=$1 total=$2 w=25
    local filled=$(( n * w / total ))
    local empty=$(( w - filled ))
    local pct=$(( n * 100 / total ))
    printf "[%s%s] %3d%%" "$(printf '#%.0s' $(seq 1 $filled 2>/dev/null) )" \
           "$(printf -- '-%.0s' $(seq 1 $empty 2>/dev/null) )" "$pct"
}

job_state() {
    local jid=$1
    squeue -j "$jid" -h -o "%T" 2>/dev/null || sacct -j "$jid" -n -o State -X 2>/dev/null | head -1 | tr -d ' '
}

count_snaps() {
    local dir="$1"
    ls "$dir"/output*.xml 2>/dev/null | wc -l
}

# ── Cell count helper via Python ──
cell_counts_at_snap() {
    local outdir="$1" snap_idx="$2"
    python3 -c "
import sys
sys.path.insert(0, '/home/a0hirs04/PROJECT-NORTHSTAR')
from python.wrapper.output_parser import OutputParser
op = OutputParser('$outdir')
try:
    m = op.parse_snapshot($snap_idx)
    print(f'{m.live_tumor_cells} {m.total_stromal_cells} {m.activated_cafs}')
except:
    print('-- -- --')
" 2>/dev/null || echo "-- -- --"
}

echo "============================================================"
echo "  RC1 + RC2 — Live Monitor   $(date '+%H:%M:%S %Z')"
echo "============================================================"

# ════════════════════ RC1 ════════════════════
echo ""
echo "  ┌─ RC1: Barrier Self-Assembly (21 days, 5 seeds) ─────────"

rc1_all_done=true
for i in "${!SEEDS[@]}"; do
    seed=${SEEDS[$i]}
    rep_idx=$((i + 1))
    rep_dir="$RC1_DIR/replicate_$(printf '%02d' $rep_idx)_seed${seed}"
    out_dir="$rep_dir/output"

    # Find SLURM job ID
    jid=$(ls "$rep_dir"/slurm_*.out 2>/dev/null | head -1 | grep -oP '\d+(?=\.out)' || echo "?")

    snaps=0
    if [[ -d "$out_dir" ]]; then
        snaps=$(count_snaps "$out_dir")
    fi

    state="?"
    if [[ "$jid" != "?" ]]; then
        state=$(job_state "$jid")
    fi
    [[ "$state" != "COMPLETED" ]] && rc1_all_done=false

    printf "  │ Seed %-5d  %s  %3d/%d snaps  job %s (%s)\n" \
        "$seed" "$(bar $snaps $RC1_TOTAL)" "$snaps" "$RC1_TOTAL" "$jid" "$state"
done

# RC1 cell counts at final snapshot (day 21 = snap 84)
if $rc1_all_done; then
    echo "  │"
    echo "  │ FINAL CELL COUNTS (Day 21):"
    printf "  │   %-8s %6s %6s %6s\n" "Seed" "Tumor" "Stroma" "CAFs"
    for i in "${!SEEDS[@]}"; do
        seed=${SEEDS[$i]}
        rep_idx=$((i + 1))
        out="$RC1_DIR/replicate_$(printf '%02d' $rep_idx)_seed${seed}/output"
        read -r t s c <<< "$(cell_counts_at_snap "$out" 84)"
        printf "  │   %-8d %6s %6s %6s\n" "$seed" "$t" "$s" "$c"
    done
fi
echo "  └─────────────────────────────────────────────────────────"

# ════════════════════ RC2 ════════════════════
echo ""
echo "  ┌─ RC2: Drug Response (42 days, 5 seeds) ────────────────"

rc2_all_done=true
for i in "${!SEEDS[@]}"; do
    seed=${SEEDS[$i]}
    rep_idx=$((i + 1))
    rep_dir="$RC2_DIR/replicate_$(printf '%02d' $rep_idx)_seed${seed}"
    out_dir="$rep_dir/output"

    jid=$(ls "$rep_dir"/slurm_*.out 2>/dev/null | head -1 | grep -oP '\d+(?=\.out)' || echo "?")

    snaps=0
    if [[ -d "$out_dir" ]]; then
        snaps=$(count_snaps "$out_dir")
    fi

    state="?"
    if [[ "$jid" != "?" ]]; then
        state=$(job_state "$jid")
    fi
    [[ "$state" != "COMPLETED" ]] && rc2_all_done=false

    # Phase indicator
    phase="barrier"
    [[ $snaps -gt 56 ]] && phase="drug-ON"
    [[ $snaps -gt 112 ]] && phase="regrowth"

    printf "  │ Seed %-5d  %s  %3d/%d snaps  %-9s  job %s (%s)\n" \
        "$seed" "$(bar $snaps $RC2_TOTAL)" "$snaps" "$RC2_TOTAL" "$phase" "$jid" "$state"
done

# RC2 cell counts at key timepoints (day14=56, day28=112, day42=168)
echo "  │"
echo "  │ CELL COUNTS AT KEY TIMEPOINTS (seed 42):"
out42="$RC2_DIR/replicate_01_seed42/output"
snaps42=$(count_snaps "$out42" 2>/dev/null || echo 0)
printf "  │   %-10s %6s %6s %6s  %s\n" "Timepoint" "Tumor" "Stroma" "CAFs" "Status"
for label_snap_status in "Day_14:56:pre-treatment" "Day_28:112:treatment-end" "Day_42:168:post-withdrawal"; do
    IFS=: read -r label snap status <<< "$label_snap_status"
    if [[ "$snaps42" -gt "$snap" ]]; then
        read -r t s c <<< "$(cell_counts_at_snap "$out42" "$snap")"
    else
        t="--"; s="--"; c="--"
    fi
    printf "  │   %-10s %6s %6s %6s  %s\n" "${label/_/ }" "$t" "$s" "$c" "$status"
done

# Quick criteria check for seed 42
if [[ "$snaps42" -gt 112 ]]; then
    read -r t14 _ _ <<< "$(cell_counts_at_snap "$out42" 56)"
    read -r t28 _ _ <<< "$(cell_counts_at_snap "$out42" 112)"
    if [[ "$t14" != "--" && "$t28" != "--" && "$t14" -gt 0 ]]; then
        reduction=$(python3 -c "print(f'{(1 - $t28/$t14)*100:.1f}')" 2>/dev/null || echo "?")
        echo "  │"
        echo "  │ EARLY RC2 CRITERIA (seed 42):"
        if (( t28 < t14 )); then
            echo "  │   RC2-1 Partial response: $t14 -> $t28 (${reduction}% reduction) [LIKELY PASS]"
        else
            echo "  │   RC2-1 Partial response: $t14 -> $t28 (${reduction}% reduction) [FAIL]"
        fi
        [[ "$t28" -gt 0 ]] && echo "  │   RC2-2 Not eradicated: $t28 survivors [PASS]"
    fi
fi

if [[ "$snaps42" -gt 168 ]]; then
    read -r t28_2 _ _ <<< "$(cell_counts_at_snap "$out42" 112)"
    read -r t42 _ _ <<< "$(cell_counts_at_snap "$out42" 168)"
    if [[ "$t28_2" != "--" && "$t42" != "--" ]]; then
        if (( t42 > t28_2 )); then
            echo "  │   RC2-6 Regrowth: $t28_2 -> $t42 [PASS]"
        else
            echo "  │   RC2-6 Regrowth: $t28_2 -> $t42 [FAIL]"
        fi
    fi
fi

echo "  └─────────────────────────────────────────────────────────"

# ════════════════════ SUMMARY ════════════════════
echo ""
if $rc1_all_done && $rc2_all_done; then
    echo "  >>> ALL JOBS COMPLETE! Run full evaluation:"
    echo "      tail -100 /tmp/rc1_run.log"
    echo "      tail -100 /tmp/rc2_run.log"
elif $rc1_all_done; then
    echo "  >>> RC1 COMPLETE. RC2 still running."
    echo "      RC1 results: tail -100 /tmp/rc1_run.log"
elif $rc2_all_done; then
    echo "  >>> RC2 COMPLETE. RC1 still running."
    echo "      RC2 results: tail -100 /tmp/rc2_run.log"
else
    echo "  >>> Both still running. RC1 ~30-40 min, RC2 ~1.5 hrs."
fi
echo "============================================================"
