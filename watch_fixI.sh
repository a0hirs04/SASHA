#!/bin/bash
# Watch fixI sweep: 2 persistence timer (RC1+RC2) + 4 induction threshold (RC1 only)
WORK="/work/a0hirs04/PROJECT-NORTHSTAR/build/fixI"

while true; do
    clear
    echo "=== fixI Monitor  $(date) ==="
    echo ""

    echo "── SLURM Jobs ──"
    squeue -u "$USER" -o "%.10i %.22j %.2t %.10M" 2>/dev/null | grep -E "JOBID|fI_" || echo "  (no jobs in queue)"
    echo ""

    echo "══ SET 1: Persistence Timer (reverted structural fix) ══"
    printf "%-5s  %11s %10s  %7s %8s\n" "Var" "persist_min" "drug_kill" "RC1" "RC2"
    echo "------------------------------------------------------"
    all_done=true
    for v in pt1 pt2; do
        rc1_snaps=$(ls "$WORK/$v/rc1/replicate_01_seed42/output"/output*.xml 2>/dev/null | wc -l)
        rc2_snaps=$(ls "$WORK/$v/rc2/output"/output*.xml 2>/dev/null | wc -l)
        [[ $rc1_snaps -ge 85 ]] && rc1_st="DONE" || { rc1_st="${rc1_snaps}/85"; all_done=false; }
        [[ $rc2_snaps -ge 169 ]] && rc2_st="DONE" || { rc2_st="${rc2_snaps}/169"; all_done=false; }

        cfg="$WORK/$v/rc1/replicate_01_seed42/PhysiCell_settings.xml"
        pt=$(grep "emt_persistence_time" "$cfg" 2>/dev/null | sed 's/.*>//' | sed 's/<.*//')
        dk=$(grep "drug_kill_coefficient" "$cfg" 2>/dev/null | sed 's/.*>//' | sed 's/<.*//')

        printf "%-5s  %11s %10s  %7s %8s\n" "$v" "$pt" "$dk" "$rc1_st" "$rc2_st"
    done

    echo ""
    echo "══ SET 2: Induction Threshold Fix (structural fix kept) ══"
    printf "%-5s  %10s  %7s\n" "Var" "threshold" "RC1"
    echo "-------------------------------"
    for v in th18 th19 th20 th22; do
        rc1_snaps=$(ls "$WORK/$v/rc1/replicate_01_seed42/output"/output*.xml 2>/dev/null | wc -l)
        [[ $rc1_snaps -ge 85 ]] && rc1_st="DONE" || { rc1_st="${rc1_snaps}/85"; all_done=false; }

        cfg="$WORK/$v/rc1/replicate_01_seed42/PhysiCell_settings.xml"
        th=$(grep "emt_induction_threshold" "$cfg" 2>/dev/null | sed 's/.*>//' | sed 's/<.*//')

        printf "%-5s  %10s  %7s\n" "$v" "$th" "$rc1_st"
    done

    echo ""
    if $all_done; then
        echo "══ ALL FIXL VARIANTS COMPLETE ══"
        break
    fi

    sleep 60
done
