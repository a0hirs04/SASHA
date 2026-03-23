#!/bin/bash
# Watch ecm_reversion_weight sweep (8 variants x RC1+RC2)
WORK="/work/a0hirs04/PROJECT-NORTHSTAR/build/reversion_weight"

while true; do
    clear
    echo "=== ECM Reversion Weight Sweep  $(date) ==="
    echo ""

    echo "── SLURM Jobs ──"
    squeue -u "$USER" -o "%.10i %.22j %.2t %.10M" 2>/dev/null | grep -E "JOBID|rw_" || echo "  (no jobs in queue)"
    echo ""

    printf "%-4s  %10s %8s %6s %6s  %7s %8s\n" "Var" "ecm_rev_wt" "emt_off" "tgfb" "hif1a" "RC1" "RC2"
    echo "--------------------------------------------------------------"

    all_done=true
    for v in v1 v2 v3 v4 v5 v6 v7 v8; do
        rc1_snaps=$(ls "$WORK/$v/rc1/replicate_01_seed42/output"/output*.xml 2>/dev/null | wc -l)
        rc2_snaps=$(ls "$WORK/$v/rc2/output"/output*.xml 2>/dev/null | wc -l)
        [[ $rc1_snaps -ge 85 ]] && rc1_st="DONE" || { rc1_st="${rc1_snaps}/85"; all_done=false; }
        [[ $rc2_snaps -ge 169 ]] && rc2_st="DONE" || { rc2_st="${rc2_snaps}/169"; all_done=false; }

        # Read params from config
        cfg="$WORK/$v/rc1/replicate_01_seed42/PhysiCell_settings.xml"
        rw=$(grep "ecm_reversion_weight" "$cfg" 2>/dev/null | sed 's/.*>//' | sed 's/<.*//')
        eo=$(grep "emt_off_threshold" "$cfg" 2>/dev/null | sed 's/.*>//' | sed 's/<.*//')
        tg=$(grep "tgfb_secretion_rate" "$cfg" 2>/dev/null | grep -v caf | head -1 | sed 's/.*>//' | sed 's/<.*//')
        hi=$(grep "hif1a_emt_boost" "$cfg" 2>/dev/null | sed 's/.*>//' | sed 's/<.*//')

        printf "%-4s  %10s %8s %6s %6s  %7s %8s\n" "$v" "$rw" "$eo" "$tg" "$hi" "$rc1_st" "$rc2_st"
    done

    echo ""
    if $all_done; then
        echo "══ ALL 8 VARIANTS COMPLETE ══"
        break
    fi

    sleep 60
done
