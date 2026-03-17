#!/bin/bash
# watch_fixJ.sh — Monitor RC1 + RC2 Fix-J jobs
# Usage:
#   bash watch_fixJ.sh [poll_seconds] [rc2_jobid] [rc1_jobid]
set -u

POLL=${1:-30}
RC2_JOB=${2:-9860}
RC1_JOB=${3:-9859}

RC2_DIR="/home/a0hirs04/PROJECT-NORTHSTAR/build/rc2_fixJ_seed42/replicate_01_seed42/output"
RC1_DIR="/home/a0hirs04/PROJECT-NORTHSTAR/build/rc1_fixJ_seed42/replicate_01_seed42/output"
RC2_EXPECTED=169
RC1_EXPECTED=85

RED='\033[0;31m'
GRN='\033[0;32m'
YLW='\033[0;33m'
CYN='\033[0;36m'
RST='\033[0m'

bar() {
    local cur=$1 total=$2 width=30
    local denom=$(( total > 0 ? total : 1 ))
    local pct=$(( cur * 100 / denom ))
    local filled=$(( cur * width / denom ))
    [[ $filled -gt $width ]] && filled=$width
    local empty=$(( width - filled ))
    printf "[${GRN}%s${RST}%s] %3d%%" \
        "$(printf '#%.0s' $(seq 1 $filled 2>/dev/null))" \
        "$(printf '.%.0s' $(seq 1 $empty 2>/dev/null))" \
        "$pct"
}

get_state() {
    local jid=${1:-}
    [[ -z "$jid" ]] && { echo "UNKNOWN"; return; }
    sacct -j "$jid" --format=State --noheader -P 2>/dev/null | head -1 | tr -d ' '
}

get_node() {
    local jid=${1:-}
    [[ -z "$jid" ]] && { echo "N/A"; return; }
    sacct -j "$jid" --format=NodeList --noheader -P 2>/dev/null | head -1 | tr -d ' '
}

get_elapsed() {
    local jid=${1:-}
    [[ -z "$jid" ]] && { echo "N/A"; return; }
    sacct -j "$jid" --format=Elapsed --noheader -P 2>/dev/null | head -1 | tr -d ' '
}

snap_count() {
    local dir=$1
    ls "$dir"/output*.xml 2>/dev/null | wc -l
}

latest_day() {
    local dir=$1
    local last
    last=$(ls "$dir"/output*.xml 2>/dev/null | sort | tail -1)
    if [[ -z "$last" ]]; then
        echo "N/A"
        return
    fi
    local idx
    idx=$(basename "$last" .xml | sed 's/output//' | sed 's/^0*//')
    [[ -z "$idx" ]] && idx=0
    local minutes=$(( 10#$idx * 360 ))
    local day=$(( minutes / 1440 ))
    echo "d${day}"
}

echo -e "\n${CYN}═══════════════════════════════════════════════════════════${RST}"
echo -e "${CYN}  Fix-J Job Monitor — RC1 (${RC1_JOB}) + RC2 (${RC2_JOB})${RST}"
echo -e "${CYN}  caf_base_death=0.0  caf_crowded_death=0.00005${RST}"
echo -e "${CYN}  CI=25  MCS=0.6  stromal_CI=0.8${RST}"
echo -e "${CYN}═══════════════════════════════════════════════════════════${RST}"

rc1_done=0
rc2_done=0

while [[ $rc1_done -eq 0 || $rc2_done -eq 0 ]]; do
    clear 2>/dev/null || true
    echo -e "${CYN}═══════════════════════════════════════════════════════════${RST}"
    echo -e "${CYN}  Fix-J Job Monitor          $(date '+%H:%M:%S')${RST}"
    echo -e "${CYN}  crowded CAF death=0.00005 | frontier CAF death=0.0${RST}"
    echo -e "${CYN}═══════════════════════════════════════════════════════════${RST}"

    rc1_state=$(get_state "$RC1_JOB")
    rc1_elapsed=$(get_elapsed "$RC1_JOB")
    rc1_node=$(get_node "$RC1_JOB")
    rc1_snaps=$(snap_count "$RC1_DIR")
    rc1_day=$(latest_day "$RC1_DIR")

    echo -e "\n  ${YLW}RC1${RST} (job ${RC1_JOB}) — 21-day untreated"
    echo -e "  State: ${GRN}${rc1_state}${RST}  Node: ${rc1_node}  Elapsed: ${rc1_elapsed}"
    printf "  Snapshots: %d / %d  (%s)  " "$rc1_snaps" "$RC1_EXPECTED" "$rc1_day"
    bar "$rc1_snaps" "$RC1_EXPECTED"
    echo ""

    if [[ "$rc1_state" == "COMPLETED" || "$rc1_state" == "FAILED" || "$rc1_state" == "CANCELLED" || "$rc1_state" == "TIMEOUT" ]]; then
        rc1_done=1
    fi

    rc2_state=$(get_state "$RC2_JOB")
    rc2_elapsed=$(get_elapsed "$RC2_JOB")
    rc2_node=$(get_node "$RC2_JOB")
    rc2_snaps=$(snap_count "$RC2_DIR")
    rc2_day=$(latest_day "$RC2_DIR")

    echo -e "\n  ${YLW}RC2${RST} (job ${RC2_JOB}) — 42-day drug+withdrawal"
    echo -e "  State: ${GRN}${rc2_state}${RST}  Node: ${rc2_node}  Elapsed: ${rc2_elapsed}"
    printf "  Snapshots: %d / %d  (%s)  " "$rc2_snaps" "$RC2_EXPECTED" "$rc2_day"
    bar "$rc2_snaps" "$RC2_EXPECTED"
    echo ""

    if [[ "$rc2_state" == "COMPLETED" || "$rc2_state" == "FAILED" || "$rc2_state" == "CANCELLED" || "$rc2_state" == "TIMEOUT" ]]; then
        rc2_done=1
    fi

    if [[ -z "$RC2_JOB" && $rc2_snaps -ge $RC2_EXPECTED ]]; then rc2_done=1; fi
    if [[ -z "$RC1_JOB" && $rc1_snaps -ge $RC1_EXPECTED ]]; then rc1_done=1; fi

    if [[ $rc1_done -eq 1 && $rc2_done -eq 1 ]]; then
        break
    fi

    echo -e "\n  Polling every ${POLL}s... (Ctrl+C to stop)"
    sleep "$POLL"
done

echo -e "\n${CYN}═══════════════════════════════════════════════════════════${RST}"
echo -e "${CYN}  BOTH JOBS FINISHED${RST}"
echo -e "${CYN}═══════════════════════════════════════════════════════════${RST}"
echo -e "  RC1 (${RC1_JOB}): $(get_state "$RC1_JOB")  —  $rc1_snaps snapshots"
echo -e "  RC2 (${RC2_JOB}): $(get_state "$RC2_JOB")  —  $rc2_snaps snapshots"
echo -e "\n  Run evaluations:"
echo -e "    python3 evaluate_rc1.py"
echo -e "    python3 evaluate_rc2.py"
echo ""
