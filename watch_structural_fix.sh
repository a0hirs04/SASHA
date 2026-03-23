#!/bin/bash
# Watch all structural fix runs: main + backup1 + backup2
WORK="/work/a0hirs04/PROJECT-NORTHSTAR/build"
RUNS=("structural_fix" "backup1" "backup2")
LABELS=("Main (off=0.10, tgfb=0.8)" "Backup1 (off=0.12, tgfb=0.8)" "Backup2 (off=0.10, tgfb=1.2)")

while true; do
    clear
    echo "=== Structural Fix Monitor  $(date) ==="
    echo ""

    # SLURM status
    echo "── SLURM Jobs ──"
    squeue -u "$USER" -o "%.10i %.20j %.2t %.10M %.6D %R" 2>/dev/null | grep -E "JOBID|sf_|bk_" || echo "  (no jobs in queue)"
    echo ""

    all_done=true
    for i in "${!RUNS[@]}"; do
        run="${RUNS[$i]}"
        label="${LABELS[$i]}"
        base="$WORK/$run"

        rc1_out="$base/rc1/replicate_01_seed42/output"
        rc2_out="$base/rc2/output"
        rc1_snaps=$(ls "$rc1_out"/output*.xml 2>/dev/null | wc -l)
        rc2_snaps=$(ls "$rc2_out"/output*.xml 2>/dev/null | wc -l)
        rc1_pct=$(( rc1_snaps * 100 / 85 ))
        rc2_pct=$(( rc2_snaps * 100 / 169 ))

        [[ $rc1_snaps -ge 85 ]] && rc1_st="DONE" || rc1_st="running"
        [[ $rc2_snaps -ge 169 ]] && rc2_st="DONE" || rc2_st="running"

        if [[ $rc1_snaps -lt 85 ]] || [[ $rc2_snaps -lt 169 ]]; then
            all_done=false
        fi

        echo "── $label ──"
        printf "  RC1: %3d/85  (%3d%%)  [%s]\n" "$rc1_snaps" "$rc1_pct" "$rc1_st"
        printf "  RC2: %3d/169 (%3d%%)  [%s]\n" "$rc2_snaps" "$rc2_pct" "$rc2_st"
        echo ""
    done

    if $all_done; then
        echo "══════════════════════════════════"
        echo "  ALL 3 RUNS COMPLETE"
        echo "══════════════════════════════════"
        for run in "${RUNS[@]}"; do
            base="$WORK/$run"
            echo "  python3 run_reality_check_1.py --work-dir $base/rc1 --seeds 42 --quorum 1 --evaluate-only"
            echo "  python3 evaluate_rc2.py --out-dir $base/rc2/output"
        done
        break
    fi

    sleep 60
done
