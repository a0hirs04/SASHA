#!/usr/bin/env bash
# watch_mega_sweep.sh — Monitor 100-variant mega sweep
# Usage:  watch -n 60 bash watch_mega_sweep.sh
set -uo pipefail

SWEEP="/work/a0hirs04/PROJECT-NORTHSTAR/build/mega_sweep"
MANIFEST="$SWEEP/jobs.txt"

# Snapshot counts for progress tracking
RC1_TOTAL=85    # 21 days / 360 min = 84 + initial
RC2_TOTAL=169   # 42 days / 360 min = 168 + initial
RC3_TOTAL=225   # 56 days / 360 min = 224 + initial

count_snaps() {
    ls "$1"/output*.xml 2>/dev/null | wc -l
}

echo "════════════════════════════════════════════════════════════"
echo "  MEGA SWEEP — 100 Variants   $(date '+%H:%M:%S %Z')"
echo "════════════════════════════════════════════════════════════"

if [[ ! -f "$MANIFEST" ]]; then
    echo "  No manifest found at $MANIFEST"
    echo "  Run launch_mega_sweep.py first."
    exit 1
fi

# Count job states
total=0; running=0; pending=0; completed=0; failed=0; other=0
while IFS=$'\t' read -r name rc jid; do
    total=$((total + 1))
    state=$(squeue -j "$jid" -h -o "%T" 2>/dev/null)
    if [[ -z "$state" ]]; then
        state=$(sacct -j "$jid" -n -o State -X 2>/dev/null | head -1 | tr -d ' ')
    fi
    case "$state" in
        RUNNING)   running=$((running + 1)) ;;
        PENDING)   pending=$((pending + 1)) ;;
        COMPLETED) completed=$((completed + 1)) ;;
        FAILED|TIMEOUT|CANCELLED*) failed=$((failed + 1)) ;;
        *)         other=$((other + 1)) ;;
    esac
done < "$MANIFEST"

echo ""
echo "  Jobs: $total total | $completed done | $running running | $pending queued | $failed failed"
pct=0
[[ $total -gt 0 ]] && pct=$(( completed * 100 / total ))
echo "  Progress: ${pct}%"
echo ""

# Show per-variant RC1/RC2 completion status
echo "  ┌─ Variant Status (RC1 / RC2) ─────────────────────────"
printf "  │ %-6s  %4s  %4s  %-20s\n" "Name" "RC1" "RC2" "Parameters"

for v in $(seq -w 1 100); do
    vname="v${v}"
    vdir="$SWEEP/$vname"
    [[ -d "$vdir" ]] || continue

    # RC1 status
    rc1_out="$vdir/rc1/replicate_01_seed42/output"
    rc1_snaps=0
    [[ -d "$rc1_out" ]] && rc1_snaps=$(count_snaps "$rc1_out")
    if [[ $rc1_snaps -ge $RC1_TOTAL ]]; then
        rc1_st="DONE"
    elif [[ $rc1_snaps -gt 0 ]]; then
        rc1_st="${rc1_snaps}/${RC1_TOTAL}"
    else
        rc1_st="wait"
    fi

    # RC2 status
    rc2_out="$vdir/rc2/output"
    rc2_snaps=0
    [[ -d "$rc2_out" ]] && rc2_snaps=$(count_snaps "$rc2_out")
    if [[ $rc2_snaps -ge $RC2_TOTAL ]]; then
        rc2_st="DONE"
    elif [[ $rc2_snaps -gt 0 ]]; then
        rc2_st="${rc2_snaps}/${RC2_TOTAL}"
    else
        rc2_st="wait"
    fi

    # Get params from variants.tsv
    params=""
    if [[ -f "$SWEEP/variants.tsv" ]]; then
        params=$(grep "^${vname}\b" "$SWEEP/variants.tsv" 2>/dev/null | \
                 awk -F'\t' '{printf "u=%.2f k=%.2f h=%.2f c=%.2f o=%.2f t=%.1f", $2,$3,$4,$5,$6,$7}')
    fi

    printf "  │ %-6s  %4s  %4s  %s\n" "$vname" "$rc1_st" "$rc2_st" "$params"
done | head -60

nvars=$(ls -d "$SWEEP"/v??? 2>/dev/null | wc -l)
if [[ $nvars -gt 60 ]]; then
    echo "  │ ... ($((nvars - 60)) more variants)"
fi

echo "  └──────────────────────────────────────────────────────"

# Estimate remaining time
if [[ $total -gt 0 && $completed -lt $total ]]; then
    remaining=$((total - completed))
    # Rough estimate: ~1 hr per 8 jobs
    est_hrs=$(python3 -c "print(f'{$remaining / 8 * 1.0:.1f}')" 2>/dev/null || echo "?")
    echo ""
    echo "  Remaining: ~$remaining jobs, est. ~${est_hrs} hours"
fi

if [[ $completed -eq $total && $total -gt 0 ]]; then
    echo ""
    echo "  >>> ALL JOBS COMPLETE! Run evaluation:"
    echo "      bash /home/a0hirs04/PROJECT-NORTHSTAR/eval_mega_sweep.sh"
fi

echo "════════════════════════════════════════════════════════════"
