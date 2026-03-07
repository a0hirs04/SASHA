#!/usr/bin/env bash
# ═══════════════════════════════════════════════════════════════════════
#  PROJECT NORTHSTAR — Live Job Monitor
#  Usage:  bash watch_job.sh [SLURM_JOB_ID]
#          (auto-detects your running jobs if no ID given)
# ═══════════════════════════════════════════════════════════════════════
set -uo pipefail

# ──────────────────────── Configuration ────────────────────────────────
PROJECT_ROOT="$(cd "$(dirname "$0")" && pwd)"
DEFAULT_OUT_DIR="$PROJECT_ROOT/build/rc2_full_seed42/replicate_01_seed42/output"

SAVE_INTERVAL=360       # minutes per snapshot
TOTAL_TIME=60480        # total sim minutes (42 days)
TOTAL_SNAPS=168         # TOTAL_TIME / SAVE_INTERVAL

SNAP_PRE=56             # end of BARRIER  (day 14)
SNAP_TREAT=112          # end of DRUG-ON  (day 28)

POLL=10                 # seconds between refreshes
BAR_WIDTH=62            # characters for progress bar

# ──────────────────────── Colors / Symbols ─────────────────────────────
RST="\033[0m"
BOLD="\033[1m"
DIM="\033[2m"
RED="\033[1;31m"
GRN="\033[1;32m"
YEL="\033[1;33m"
CYN="\033[1;36m"
BLU="\033[1;34m"
MAG="\033[1;35m"
WHT="\033[1;37m"

SPIN_CHARS=('⠋' '⠙' '⠹' '⠸' '⠼' '⠴' '⠦' '⠧' '⠇' '⠏')
SPARK_CHARS=('▁' '▂' '▃' '▄' '▅' '▆' '▇' '█')

# ──────────────────────── Auto-detect Job ──────────────────────────────
if [ -n "${1:-}" ]; then
    JOB_ID="$1"
else
    JOB_ID=$(squeue -u "$USER" --format="%i %j" --noheader 2>/dev/null \
             | grep rc2_full | head -1 | awk '{print $1}' | tr -d ' ')
    if [ -z "$JOB_ID" ]; then
        JOB_ID=$(squeue -u "$USER" --format="%i" --noheader 2>/dev/null | head -1 | tr -d ' ')
    fi
    if [ -z "$JOB_ID" ]; then
        echo -e "${RED}ERROR:${RST} No running jobs found. Pass a job ID: bash watch_job.sh <JOB_ID>"
        exit 1
    fi
    echo -e "${DIM}Auto-detected job: ${JOB_ID}${RST}"
fi

OUT_DIR="$DEFAULT_OUT_DIR"

# ──────────────────────── Helper Functions ─────────────────────────────

job_info() {
    squeue -j "$JOB_ID" --format="%T|%M|%l|%P|%D|%N|%j|%C|%m" --noheader 2>/dev/null | head -1
}

job_state_sacct() {
    sacct -j "$JOB_ID" --format=State --noheader 2>/dev/null | head -1 | tr -d ' '
}

phase_of() {
    local s=$1
    if   (( s <= SNAP_PRE ));   then echo "BARRIER"
    elif (( s <= SNAP_TREAT )); then echo "DRUG-ON"
    else                              echo "REGROWTH"
    fi
}

phase_color() {
    case "$1" in
        BARRIER)  echo -ne "$BLU" ;;
        DRUG-ON)  echo -ne "$MAG" ;;
        REGROWTH) echo -ne "$CYN" ;;
        *)        echo -ne "$WHT" ;;
    esac
}

state_color() {
    case "$1" in
        RUNNING)    echo -ne "$GRN" ;;
        COMPLETED)  echo -ne "$CYN" ;;
        PENDING)    echo -ne "$YEL" ;;
        *)          echo -ne "$RED" ;;
    esac
}

tumor_count() {
    python3 -c "
import sys, numpy as np
sys.path.insert(0, '$PROJECT_ROOT')
from python.wrapper.output_parser import OutputParser
p = OutputParser('$OUT_DIR')
s = p._read_physicell_xml('$1')
m = s['cell_matrix']; l = s['label_name_map']
ct = m[int(l['cell_type']['index']), :]
dead = m[int(l['dead']['index']), :] if 'dead' in l else np.zeros(m.shape[1])
dm = m[int(l['current_death_model']['index']), :] if 'current_death_model' in l else np.zeros(m.shape[1])
tumor = (np.rint(ct).astype(int)==0) & (dead<0.5) & (np.rint(dm).astype(int)!=100)
stroma = (np.rint(ct).astype(int)==1) & (dead<0.5) & (np.rint(dm).astype(int)!=100)
print(f'{int(np.sum(tumor))}|{int(np.sum(stroma))}|{int(m.shape[1])}')
" 2>/dev/null || echo "?|?|?"
}

sparkline() {
    local arr=("$@")
    local n=${#arr[@]}
    if (( n < 2 )); then return; fi
    local min=${arr[0]} max=${arr[0]}
    for v in "${arr[@]}"; do
        (( v < min )) && min=$v
        (( v > max )) && max=$v
    done
    local range=$(( max - min ))
    (( range == 0 )) && range=1
    local start=0
    (( n > 30 )) && start=$(( n - 30 ))
    for (( i=start; i<n; i++ )); do
        local idx=$(( (arr[i] - min) * 7 / range ))
        (( idx > 7 )) && idx=7
        (( idx < 0 )) && idx=0
        printf "%s" "${SPARK_CHARS[$idx]}"
    done
}

build_progress_bar() {
    local n_snaps=$1
    local seg1=$(( SNAP_PRE * BAR_WIDTH / TOTAL_SNAPS ))
    local seg2=$(( (SNAP_TREAT - SNAP_PRE) * BAR_WIDTH / TOTAL_SNAPS ))
    local filled=$(( n_snaps * BAR_WIDTH / TOTAL_SNAPS ))
    for (( i=0; i<BAR_WIDTH; i++ )); do
        if (( i == seg1 )) || (( i == seg1 + seg2 )); then
            printf "${WHT}│${RST}"
        elif (( i < filled )); then
            if   (( i < seg1 ));          then printf "${BLU}█${RST}"
            elif (( i < seg1 + seg2 ));   then printf "${MAG}█${RST}"
            else                               printf "${CYN}█${RST}"
            fi
        else
            printf "${DIM}░${RST}"
        fi
    done
}

fmt_time() {
    local s=$1
    printf "%02d:%02d:%02d" $((s/3600)) $(((s%3600)/60)) $((s%60))
}

# ──────────────────────── State Variables ──────────────────────────────
LAST_TUMOR_XML=""
TUMOR_LIVE="?"
STROMA_LIVE="?"
TOTAL_CELLS="?"
declare -a TUMOR_HISTORY=()
START_TIME=$(date +%s)
SPIN_IDX=0
TICK=0
DISK="..."

# ──────────────────────── Banner ───────────────────────────────────────
clear
echo ""
echo -e "${BLU}${BOLD}"
cat << 'BANNER'
    ██╗   ██╗ ██████╗ ██████╗ ████████╗██╗  ██╗███████╗████████╗ █████╗ ██████╗
    ███╗  ██║██╔═══██╗██╔══██╗╚══██╔══╝██║  ██║██╔════╝╚══██╔══╝██╔══██╗██╔══██╗
    ████╗ ██║██║   ██║██████╔╝   ██║   ███████║███████╗   ██║   ███████║██████╔╝
    ██╔██╗██║██║   ██║██╔══██╗   ██║   ██╔══██║╚════██║   ██║   ██╔══██║██╔══██╗
    ██║╚████║╚██████╔╝██║  ██║   ██║   ██║  ██║███████║   ██║   ██║  ██║██║  ██║
    ╚═╝ ╚═══╝ ╚═════╝ ╚═╝  ╚═╝   ╚═╝   ╚═╝  ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═╝
BANNER
echo -e "${RST}"

# Print placeholder lines for in-place update
DISPLAY_LINES=17
for (( i=0; i<DISPLAY_LINES; i++ )); do echo ""; done

# ──────────────────────── Main Loop ────────────────────────────────────
while true; do
    NOW=$(date +%s)
    WALL_ELAPSED=$(( NOW - START_TIME ))
    SPIN_IDX=$(( (SPIN_IDX + 1) % ${#SPIN_CHARS[@]} ))
    SPINNER="${SPIN_CHARS[$SPIN_IDX]}"
    TICK=$(( TICK + 1 ))

    # ── Job info ───────────────────────────────────────────────────
    INFO=$(job_info)
    if [ -n "$INFO" ]; then
        IFS='|' read -r STATE JOB_ELAPSED TIMELIMIT PARTITION NODES NODELIST JOBNAME NCPUS MEM <<< "$INFO"
        STATE=$(echo "$STATE" | tr -d ' ')
    else
        STATE=$(job_state_sacct)
        JOB_ELAPSED="N/A"; TIMELIMIT="N/A"; PARTITION="?"; NODELIST="?"; JOBNAME="?"; NCPUS="?"; MEM="?"
    fi

    # ── Snapshots ──────────────────────────────────────────────────
    N_SNAPS=$(find "$OUT_DIR" -maxdepth 1 -name "output*.xml" 2>/dev/null | wc -l)
    SIM_DAY=$(echo "scale=1; $N_SNAPS * $SAVE_INTERVAL / 1440" | bc 2>/dev/null || echo "?")
    PHASE=$(phase_of "$N_SNAPS")
    PCT=$(( N_SNAPS * 100 / TOTAL_SNAPS ))

    # ── ETA ────────────────────────────────────────────────────────
    if (( N_SNAPS > 3 && WALL_ELAPSED > 0 )); then
        REMAINING=$(( TOTAL_SNAPS - N_SNAPS ))
        ETA_SECS=$(( WALL_ELAPSED * REMAINING / N_SNAPS ))
        ETA_STR=$(fmt_time "$ETA_SECS")
        SPEED=$(echo "scale=2; $N_SNAPS * 60 / $WALL_ELAPSED" | bc 2>/dev/null || echo "?")
    else
        ETA_STR="calculating..."
        SPEED="?"
    fi

    # ── Tumor count (on new snapshot) ──────────────────────────────
    LATEST_XML=$(find "$OUT_DIR" -maxdepth 1 -name "output*.xml" 2>/dev/null | sort | tail -1)
    if [ -n "${LATEST_XML:-}" ] && [ "${LATEST_XML:-}" != "${LAST_TUMOR_XML:-}" ]; then
        COUNTS=$(tumor_count "$LATEST_XML")
        IFS='|' read -r TUMOR_LIVE STROMA_LIVE TOTAL_CELLS <<< "$COUNTS"
        LAST_TUMOR_XML="$LATEST_XML"
        if [ "$TUMOR_LIVE" != "?" ]; then
            TUMOR_HISTORY+=("$TUMOR_LIVE")
        fi
    fi

    # ── Disk usage (every 6 ticks) ─────────────────────────────────
    if (( TICK % 6 == 1 )); then
        DISK=$(du -sh "$OUT_DIR" 2>/dev/null | cut -f1)
    fi

    # ── Sparkline ──────────────────────────────────────────────────
    SPARK=""
    if (( ${#TUMOR_HISTORY[@]} > 1 )); then
        SPARK=$(sparkline "${TUMOR_HISTORY[@]}")
    fi

    # ── Wall time bar ──────────────────────────────────────────────
    TLIMIT_BAR=""
    if [[ "${TIMELIMIT:-}" =~ ^([0-9]+):([0-9]+):([0-9]+)$ ]]; then
        TL_SECS=$(( ${BASH_REMATCH[1]} * 3600 + ${BASH_REMATCH[2]} * 60 + ${BASH_REMATCH[3]} ))
        JE_SECS=0
        if [[ "${JOB_ELAPSED:-}" =~ ^([0-9]+):([0-9]+):([0-9]+)$ ]]; then
            JE_SECS=$(( ${BASH_REMATCH[1]} * 3600 + ${BASH_REMATCH[2]} * 60 + ${BASH_REMATCH[3]} ))
        elif [[ "${JOB_ELAPSED:-}" =~ ^([0-9]+):([0-9]+)$ ]]; then
            JE_SECS=$(( ${BASH_REMATCH[1]} * 60 + ${BASH_REMATCH[2]} ))
        fi
        if (( TL_SECS > 0 )); then
            TIME_PCT=$(( JE_SECS * 100 / TL_SECS ))
            TIME_FILLED=$(( JE_SECS * 20 / TL_SECS ))
            TLIMIT_BAR="["
            for (( i=0; i<20; i++ )); do
                if (( i < TIME_FILLED )); then
                    if (( TIME_PCT > 80 )); then TLIMIT_BAR+=$(printf "${RED}█${RST}")
                    elif (( TIME_PCT > 60 )); then TLIMIT_BAR+=$(printf "${YEL}█${RST}")
                    else TLIMIT_BAR+=$(printf "${GRN}█${RST}")
                    fi
                else
                    TLIMIT_BAR+=$(printf "${DIM}░${RST}")
                fi
            done
            TLIMIT_BAR+="] ${TIME_PCT}%"
        fi
    fi

    # ── SVG count ────────────────────────────────────────────────
    N_SVGS=$(find "$OUT_DIR" -maxdepth 1 -name "snapshot*.svg" 2>/dev/null | wc -l)

    # ── Draw ───────────────────────────────────────────────────────
    printf "\033[${DISPLAY_LINES}A"

    SC=$(state_color "$STATE")
    PC=$(phase_color "$PHASE")

    printf "  ${DIM}━━━━━━━━━━━━━━━━━━━━━━━━ Job Info ━━━━━━━━━━━━━━━━━━━━━━━━━${RST}\n"
    printf "  ${BOLD}Job ID:${RST}  %-10s  ${BOLD}Name:${RST}  %-16s  ${BOLD}State:${RST} %b● ${BOLD}%s${RST} %s\n" \
        "$JOB_ID" "$JOBNAME" "$SC" "$STATE" "$SPINNER"
    printf "  ${BOLD}Partition:${RST} %-8s  ${BOLD}Node:${RST}  %-16s  ${BOLD}CPUs:${RST}  %-4s ${BOLD}Mem:${RST} %s\n" \
        "$PARTITION" "$NODELIST" "$NCPUS" "$MEM"
    printf "  ${BOLD}Elapsed:${RST}  %-10s  ${BOLD}Limit:${RST} %-16s  %b\n" \
        "$JOB_ELAPSED" "$TIMELIMIT" "$TLIMIT_BAR"

    printf "  ${DIM}━━━━━━━━━━━━━━━━━━━━ Simulation Progress ━━━━━━━━━━━━━━━━━━${RST}\n"
    printf "  ["
    build_progress_bar "$N_SNAPS"
    printf "]\n"
    printf "  ${BLU}■${RST} BARRIER (d0-14)  ${WHT}│${RST}  ${MAG}■${RST} DRUG-ON (d14-28)  ${WHT}│${RST}  ${CYN}■${RST} REGROWTH (d28-42)\n"
    echo ""
    printf "  ${BOLD}Day:${RST} %-6s  ${BOLD}Snap:${RST} %3d / %-3d [%b%3d%%${RST}]  ${BOLD}Phase:${RST} %b%-10s${RST}  ${BOLD}ETA:${RST} %s\n" \
        "$SIM_DAY" "$N_SNAPS" "$TOTAL_SNAPS" "$PC" "$PCT" "$PC" "$PHASE" "$ETA_STR"
    printf "  ${BOLD}Speed:${RST} %s snap/min %28s ${BOLD}SVGs:${RST} %-5s ${BOLD}Disk:${RST} %s\n" \
        "$SPEED" "" "$N_SVGS" "${DISK:-...}"

    printf "  ${DIM}━━━━━━━━━━━━━━━━━━━━━━ Cell Counts ━━━━━━━━━━━━━━━━━━━━━━━━${RST}\n"
    printf "  🧬 ${BOLD}Tumor:${RST} %-8s  🔬 ${BOLD}Stroma:${RST} %-8s  📊 ${BOLD}Total:${RST} %-8s\n" \
        "$TUMOR_LIVE" "$STROMA_LIVE" "$TOTAL_CELLS"
    if [ -n "$SPARK" ]; then
        printf "  ${BOLD}Tumor trend:${RST} ${GRN}%s${RST}\n" "$SPARK"
    else
        printf "  ${DIM}Tumor trend: (collecting data...)${RST}\n"
    fi
    echo ""
    printf "  ${DIM}Monitor: %-8s │ Refresh: ${POLL}s │ %s${RST}\n" \
        "$(fmt_time $WALL_ELAPSED)" "$(date '+%H:%M:%S')"
    printf "  ${DIM}Ctrl+C to detach (job continues running)${RST}\n"

    # ── Exit on terminal state ─────────────────────────────────────
    if [ "$STATE" = "COMPLETED" ] || [ "$STATE" = "FAILED" ] || \
       [[ "$STATE" == CANCELLED* ]] || [ "$STATE" = "TIMEOUT" ]; then
        echo ""
        if [ "$STATE" = "COMPLETED" ]; then
            echo -e "  ${GRN}✓ Job completed!${RST}  Run: python3 evaluate_rc2.py"
        else
            echo -e "  ${RED}✗ Job ended: ${STATE}${RST}"
            ERR=$(find "$(dirname "$OUT_DIR")" -name "slurm_${JOB_ID}.err" 2>/dev/null | head -1)
            [ -n "$ERR" ] && [ -f "$ERR" ] && tail -5 "$ERR" 2>/dev/null | sed 's/^/    /'
        fi
        break
    fi

    sleep "$POLL"
done
