#!/bin/bash
#watch_fixK.sh: 8-jobs combined sweep (persistence timer + structural fix)^U

WORK="/work/a0hirs04/PROJECT-NORTHSTAR/build/fixK"

while true; do
    clear
    echo "===fixK Monitor A-D variants $(date) ==="
    echo ""
    echo "── SLURM Queue ──"
    squeue -u "$USER" -o "%.10i %.22j %.2t %.10M" 2>/dev/null | 
    grep -E "JOBID|fK_" || echo "no fK jobs"
    echo ""
    echo "=== Progress ==="
    printf "%-5s %-7s %7s %s\n" "Var" "RC1" "RC2"
    echo "-----------------------------"
    all_done=true
    for v in a b c d; do
        rc1_snaps=$(ls "$WORK/$v/rc1/output"/output*.xml 2>/dev/null | wc -l)
        rc2_snaps=$(ls "$WORK/$v/rc2/output"/output*.xml 2>/dev/null | wc -l)

        #RC1 status
        if [[ $rc1_snaps -ge 85 ]]; then
            rc1_st="DONE"
        else
            rc1_st="${rc1_snaps}/85"
            all_done=false
        fi
        
        #RC2 status
        if [[ $rc2_snaps -ge 169 ]]; then
            rc2_st="DONE"
        elif [[ $rc2_snaps -gt 0 ]]; then
            rc2_st="${rc2_snaps}/169"
            all_done=false
        

        else
            rc2_st="-"
        fi
        printf "%-5s %-7s %7s %s\n" "$v" "$rc1_st" "$rc2_st"
        done
    echo ""
    if $all_done; then
        echo "══ ALL VARIANTS COMPLETE ══"
        break
    fi
    sleep 60
done
