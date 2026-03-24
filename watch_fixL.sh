#!/bin/bash
# watch_fixL.sh: Monitor fixL sweep (drug toxicity trap fix)
# Shows fast_rc2 screening + full RC1/RC2 status
WORK="/work/a0hirs04/PROJECT-NORTHSTAR/build/fixL"

while true; do
    clear
    echo "=== fixL Monitor $(date) ==="
    echo ""
    echo "── SLURM Queue ──"
    squeue -u "$USER" -o "%.10i %.22j %.2t %.10M" 2>/dev/null |
    grep -E "JOBID|fL" || echo "  (no fL jobs)"
    echo ""

    # ── Fast RC2 Screening ──
    has_fast=false
    for v in a b c d; do
        if [[ -d "$WORK/$v/fast_rc2" ]]; then has_fast=true; break; fi
    done

    if $has_fast; then
        echo "=== Fast RC2 Screening ==="
        printf "%-5s %10s %11s  %8s %s\n" "Var" "drug_kill" "decay_rate" "Status" "Detail"
        echo "-------------------------------------------------------------"
        for v in a b c d; do
            fast_dir="$WORK/$v/fast_rc2/output"
            cfg="$WORK/$v/fast_rc2/PhysiCell_settings.xml"
            if [[ -f "$cfg" ]]; then
                dk=$(grep "drug_kill_coefficient" "$cfg" 2>/dev/null | sed 's/.*>//' | sed 's/<.*//')
                dr=$(grep "drug_natural_decay_rate" "$cfg" 2>/dev/null | sed 's/.*>//' | sed 's/<.*//')
            else
                dk="-"; dr="-"
            fi

            if [[ -f "$fast_dir/FAST_SCREEN_FAIL.txt" ]]; then
                fast_st="FAIL"
                detail=$(head -1 "$fast_dir/FAST_SCREEN_FAIL.txt" 2>/dev/null)
            elif [[ -f "$fast_dir/fast_screen.csv" ]]; then
                lines=$(wc -l < "$fast_dir/fast_screen.csv" 2>/dev/null)
                lines=$((lines - 1))  # subtract header
                fast_snaps=$(ls "$fast_dir"/output*.xml 2>/dev/null | wc -l)
                # Check if simulation finished (max_time reached)
                if [[ -f "$fast_dir/final.xml" ]] || [[ $fast_snaps -ge 42 ]]; then
                    fast_st="DONE"
                    detail="${lines} metrics logged"
                else
                    fast_st="RUN"
                    detail="${lines} metrics, ${fast_snaps} snaps"
                fi
            elif [[ -d "$WORK/$v/fast_rc2" ]]; then
                fast_st="QUEUE"
                detail="waiting"
            else
                fast_st="-"
                detail=""
            fi
            printf "%-5s %10s %11s  %8s %s\n" "$v" "$dk" "$dr" "$fast_st" "$detail"
        done
        echo ""
    fi

    # ── Full RC1/RC2 Progress ──
    echo "=== Full RC1/RC2 Progress ==="
    printf "%-5s %10s %11s  %7s %7s\n" "Var" "drug_kill" "decay_rate" "RC1" "RC2"
    echo "----------------------------------------------"
    all_done=true
    any_full=false
    for v in a b c d; do
        rc1_snaps=$(ls "$WORK/$v/rc1/replicate_01_seed42/output"/output*.xml 2>/dev/null | wc -l)
        rc2_snaps=$(ls "$WORK/$v/rc2/output"/output*.xml 2>/dev/null | wc -l)
        full_rc2_snaps=$(ls "$WORK/$v/full_rc2/output"/output*.xml 2>/dev/null | wc -l)

        # Use full_rc2 if it exists (from --promote), else rc2
        if [[ $full_rc2_snaps -gt 0 ]]; then
            rc2_snaps=$full_rc2_snaps
        fi

        if [[ $rc1_snaps -gt 0 ]] || [[ $rc2_snaps -gt 0 ]]; then
            any_full=true
        fi

        # RC1 status
        if [[ $rc1_snaps -ge 85 ]]; then
            rc1_st="DONE"
        elif [[ $rc1_snaps -gt 0 ]]; then
            rc1_st="${rc1_snaps}/85"
            all_done=false
        else
            rc1_st="-"
        fi

        # RC2 status
        if [[ $rc2_snaps -ge 169 ]]; then
            rc2_st="DONE"
        elif [[ $rc2_snaps -gt 0 ]]; then
            rc2_st="${rc2_snaps}/169"
            all_done=false
        else
            rc2_st="-"
        fi

        # Read params from config if available
        cfg="$WORK/$v/rc2/PhysiCell_settings.xml"
        [[ ! -f "$cfg" ]] && cfg="$WORK/$v/full_rc2/PhysiCell_settings.xml"
        if [[ -f "$cfg" ]]; then
            dk=$(grep "drug_kill_coefficient" "$cfg" 2>/dev/null | sed 's/.*>//' | sed 's/<.*//')
            dr=$(grep "drug_natural_decay_rate" "$cfg" 2>/dev/null | sed 's/.*>//' | sed 's/<.*//')
        else
            dk="-"; dr="-"
        fi

        printf "%-5s %10s %11s  %7s %7s\n" "$v" "$dk" "$dr" "$rc1_st" "$rc2_st"
    done
    echo ""

    if $any_full && $all_done; then
        echo "══ ALL FULL RUNS COMPLETE ══"
        break
    fi

    # If only fast runs exist and they're all done, also break
    if $has_fast && ! $any_full; then
        fast_all_done=true
        for v in a b c d; do
            fast_dir="$WORK/$v/fast_rc2/output"
            if [[ -d "$WORK/$v/fast_rc2" ]] && \
               [[ ! -f "$fast_dir/FAST_SCREEN_FAIL.txt" ]] && \
               [[ ! -f "$fast_dir/final.xml" ]]; then
                fast_snaps=$(ls "$fast_dir"/output*.xml 2>/dev/null | wc -l)
                [[ $fast_snaps -lt 42 ]] && fast_all_done=false
            fi
        done
        if $fast_all_done; then
            echo "══ ALL FAST SCREENS COMPLETE ══"
            echo "Run: python3 evaluate_fast_rc2.py -v"
            break
        fi
    fi

    sleep 60
done
