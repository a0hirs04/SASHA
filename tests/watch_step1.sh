#!/bin/bash
# Step-1 validation live monitor — run with: bash tests/watch_step1.sh
# Watches SLURM jobs 17696 (rc1_baseline_seed42) and 17697 (rc2_probe_seed42)
# Press Ctrl+C to stop.

JOB_RC1=17700
JOB_RC2=17701

DIR_RC1="/home/a0hirs04/PROJECT-NORTHSTAR/build/step1_validation/rc1_baseline_seed42"
DIR_RC2="/home/a0hirs04/PROJECT-NORTHSTAR/build/step1_validation/rc2_probe_seed42"

SNAPS_RC1=56   # 20160 / 360
SNAPS_RC2=64   # 23040 / 360

REFRESH=20

bar() {
    local n=$1 total=$2 width=25
    local filled=$(( n * width / total ))
    local b=""
    for ((i=0; i<filled; i++));  do b+="█"; done
    for ((i=filled; i<width; i++)); do b+="░"; done
    printf "%s" "$b"
}

job_state() {
    local jid=$1
    # Check squeue first (fast for running jobs)
    local st
    st=$(squeue -j "$jid" --format="%T" --noheader 2>/dev/null | head -1 | tr -d ' ')
    if [ -n "$st" ]; then
        printf "%s" "$st"
        return
    fi
    # Fall back to sacct for terminal states
    st=$(sacct -j "$jid" --format=State --noheader --parsable2 2>/dev/null \
         | head -1 | cut -d'|' -f1 | tr -d ' ')
    [ -n "$st" ] && printf "%s" "$st" || printf "UNKNOWN"
}

state_color() {
    case "$1" in
        RUNNING)   printf "\033[32m" ;;   # green
        PENDING)   printf "\033[33m" ;;   # yellow
        COMPLETED) printf "\033[36m" ;;   # cyan
        FAILED|CANCELLED|TIMEOUT) printf "\033[31m" ;;  # red
        *)         printf "\033[0m"  ;;
    esac
}

RESET="\033[0m"
BOLD="\033[1m"
DIM="\033[2m"

while true; do
    clear

    ST_RC1=$(job_state $JOB_RC1)
    ST_RC2=$(job_state $JOB_RC2)

    N_RC1=$(ls "$DIR_RC1"/output/output*.xml 2>/dev/null | wc -l)
    N_RC2=$(ls "$DIR_RC2"/output/output*.xml 2>/dev/null | wc -l)

    PCT_RC1=$(( N_RC1 * 100 / SNAPS_RC1 ))
    PCT_RC2=$(( N_RC2 * 100 / SNAPS_RC2 ))

    BAR_RC1=$(bar $N_RC1 $SNAPS_RC1)
    BAR_RC2=$(bar $N_RC2 $SNAPS_RC2)

    # Elapsed from sacct
    ELAPSED_RC1=$(sacct -j $JOB_RC1 --format=Elapsed --noheader --parsable2 2>/dev/null \
                  | head -1 | cut -d'|' -f1 | tr -d ' ')
    ELAPSED_RC2=$(sacct -j $JOB_RC2 --format=Elapsed --noheader --parsable2 2>/dev/null \
                  | head -1 | cut -d'|' -f1 | tr -d ' ')

    printf "${BOLD}═══════════════════════════════════════════════════════════════${RESET}\n"
    printf "${BOLD}  STEP-1 VALIDATION MONITOR   %s${RESET}\n" "$(date '+%H:%M:%S  %a %b %d')"
    printf "${BOLD}═══════════════════════════════════════════════════════════════${RESET}\n"
    echo ""

    # ── Job table ────────────────────────────────────────────────────────
    printf "  ${BOLD}%-8s  %-24s  %-12s  %-10s  %s${RESET}\n" \
           "JOB ID" "NAME" "STATE" "ELAPSED" "NODE"
    printf "  %s\n" "───────────────────────────────────────────────────────────"

    for row in \
        "$JOB_RC1|step1_rc1_s42|$ST_RC1|${ELAPSED_RC1:----}" \
        "$JOB_RC2|step1_rc2probe_s42|$ST_RC2|${ELAPSED_RC2:----}"
    do
        IFS='|' read -r jid jname jstate jelapsed <<< "$row"
        node=$(squeue -j "$jid" --format="%N" --noheader 2>/dev/null | head -1 | tr -d ' ')
        [ -z "$node" ] && node="-"
        COL=$(state_color "$jstate")
        printf "  %-8s  %-24s  ${COL}%-12s${RESET}  %-10s  %s\n" \
               "$jid" "$jname" "$jstate" "$jelapsed" "$node"
    done
    echo ""

    # ── Progress bars ─────────────────────────────────────────────────────
    printf "  ${BOLD}SNAPSHOT PROGRESS${RESET}  (interval=360 min)\n"
    printf "  %s\n" "───────────────────────────────────────────────────────────"

    # RC1 baseline (day 0–14, 56 snaps)
    COL1=$(state_color "$ST_RC1")
    printf "  ${BOLD}17696${RESET}  rc1_baseline_seed42\n"
    printf "         day 0–14  (max 20160 min)   target: %d snaps\n" $SNAPS_RC1
    printf "         [%s]  %3d / %d  (%d%%)  ${COL1}%s${RESET}\n" \
           "$BAR_RC1" $N_RC1 $SNAPS_RC1 $PCT_RC1 "$ST_RC1"
    echo ""

    # RC2 probe (day 0–16, 64 snaps)
    COL2=$(state_color "$ST_RC2")
    printf "  ${BOLD}17697${RESET}  rc2_probe_seed42\n"
    printf "         day 0–16  (max 23040 min)   target: %d snaps\n" $SNAPS_RC2
    printf "         [%s]  %3d / %d  (%d%%)  ${COL2}%s${RESET}\n" \
           "$BAR_RC2" $N_RC2 $SNAPS_RC2 $PCT_RC2 "$ST_RC2"
    echo ""

    # ── Drug schedule events ──────────────────────────────────────────────
    printf "  ${BOLD}DRUG SCHEDULE EVENTS${RESET}\n"
    printf "  %s\n" "───────────────────────────────────────────────────────────"
    for entry in \
        "rc1_baseline_seed42|$DIR_RC1" \
        "rc2_probe_seed42|$DIR_RC2"
    do
        IFS='|' read -r ename edir <<< "$entry"
        errfile=$(ls "$edir"/slurm_*.err 2>/dev/null | head -1)
        if [ -n "$errfile" ]; then
            msgs=$(grep "DRUG_SCHEDULE" "$errfile" 2>/dev/null)
            if [ -n "$msgs" ]; then
                printf "  ${DIM}%s:${RESET}\n" "$ename"
                echo "$msgs" | sed 's/^/    /'
            else
                printf "  ${DIM}%s: (none yet)${RESET}\n" "$ename"
            fi
        else
            printf "  ${DIM}%s: (no log file yet)${RESET}\n" "$ename"
        fi
    done
    echo ""

    # ── Recent simulation output lines ────────────────────────────────────
    printf "  ${BOLD}RECENT SIM OUTPUT${RESET}\n"
    printf "  %s\n" "───────────────────────────────────────────────────────────"
    for entry in \
        "17696 rc1_baseline|$DIR_RC1" \
        "17697 rc2_probe  |$DIR_RC2"
    do
        IFS='|' read -r elabel edir <<< "$entry"
        outfile=$(ls "$edir"/slurm_*.out 2>/dev/null | head -1)
        if [ -n "$outfile" ]; then
            last=$(tail -3 "$outfile" 2>/dev/null | grep -v "^$" | tail -1)
            printf "  ${DIM}%s:${RESET}  %s\n" "$elabel" "${last:-(no output yet)}"
        else
            printf "  ${DIM}%s:${RESET}  (no log yet)\n" "$elabel"
        fi
    done
    echo ""

    # ── Completion check ──────────────────────────────────────────────────
    rc1_done=false
    rc2_done=false
    [[ "$ST_RC1" =~ ^(COMPLETED|FAILED|CANCELLED|TIMEOUT)$ ]] && rc1_done=true
    [[ "$ST_RC2" =~ ^(COMPLETED|FAILED|CANCELLED|TIMEOUT)$ ]] && rc2_done=true

    if $rc1_done && $rc2_done; then
        printf "  ${BOLD}\033[36m*** BOTH JOBS DONE — run evaluation ***${RESET}\n"
    elif $rc1_done; then
        printf "  ${DIM}17696 finished.  Waiting for 17697...${RESET}\n"
    elif $rc2_done; then
        printf "  ${DIM}17697 finished.  Waiting for 17696...${RESET}\n"
    else
        printf "  ${DIM}Both jobs running.  Next refresh in ${REFRESH}s  (Ctrl+C to stop)${RESET}\n"
    fi

    echo ""
    sleep $REFRESH
done
