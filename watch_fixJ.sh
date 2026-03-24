#!/bin/bash
# Watch fixJ: 20-variant combined sweep (persistence timer + structural fix)
WORK="/work/a0hirs04/PROJECT-NORTHSTAR/build/fixJ"

while true; do
    clear
    echo "=== fixJ Monitor (20 variants)  $(date) ==="
    echo ""

    echo "── SLURM Queue ──"
    squeue -u "$USER" -o "%.10i %.22j %.2t %.10M" 2>/dev/null | grep -E "JOBID|fJ_" || echo "  (no fJ jobs)"
    echo ""

    echo "══ RC1 Progress (85 snapshots = DONE) ══"
    printf "%-5s  %7s %8s %7s  %7s %s\n" "Var" "persist" "ecm_rev" "drug_k" "RC1" "RC2"
    echo "------------------------------------------------"

    all_done=true
    for v in v01 v02 v03 v04 v05 v06 v07 v08 v09 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20; do
        rc1_snaps=$(ls "$WORK/$v/rc1/replicate_01_seed42/output"/output*.xml 2>/dev/null | wc -l)
        rc2_snaps=$(ls "$WORK/$v/rc2/output"/output*.xml 2>/dev/null | wc -l)

        [[ $rc1_snaps -ge 85 ]] && rc1_st="DONE" || { rc1_st="${rc1_snaps}/85"; all_done=false; }

        # Read params from config if available
        cfg="$WORK/$v/rc1/replicate_01_seed42/PhysiCell_settings.xml"
        if [[ -f "$cfg" ]]; then
            pt=$(grep "emt_persistence_time" "$cfg" 2>/dev/null | sed 's/.*>//' | sed 's/<.*//')
            rv=$(grep "ecm_reversion_weight" "$cfg" 2>/dev/null | sed 's/.*>//' | sed 's/<.*//')
            dk=$(grep "drug_kill_coefficient" "$cfg" 2>/dev/null | sed 's/.*>//' | sed 's/<.*//')
        else
            pt="-"; rv="-"; dk="-"
        fi

        # Show RC2 if it exists
        rc2_info="-"
        if [[ $rc2_snaps -gt 0 ]]; then
            [[ $rc2_snaps -ge 169 ]] && rc2_info="DONE" || { rc2_info="${rc2_snaps}/169"; all_done=false; }
        fi

        printf "%-5s  %7s %8s %7s  %7s %s\n" "$v" "$pt" "$rv" "$dk" "$rc1_st" "$rc2_info"
    done

    echo ""
    if $all_done; then
        echo "══ ALL VARIANTS COMPLETE ══"
        break
    fi

    sleep 60
done
