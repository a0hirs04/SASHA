#!/bin/bash
# ─────────────────────────────────────────────────────────────────────
#  watch_fixA.sh — Monitor RC1 + RC2 Fix-A jobs in parallel
#  Usage:  bash watch_fixA.sh [poll_interval_seconds]
# ─────────────────────────────────────────────────────────────────────
set -u

POLL=${1:-30}
RC2_JOB=9530
RC1_JOB=9531
RC2_DIR="/home/a0hirs04/PROJECT-NORTHSTAR/build/rc2_fixA_seed42/replicate_01_seed42/output"
RC1_DIR="/home/a0hirs04/PROJECT-NORTHSTAR/build/rc1_fixA_seed42/replicate_01_seed42/output"

# Expected snapshots: RC2 = 42 days = 60480 min / 360 = 168 + initial = 169
#                     RC1 = 21 days = 30240 min / 360 =  84 + initial =  85
RC2_EXPECTED=169
RC1_EXPECTED=85

RED='\033[0;31m'
GRN='\033[0;32m'
YLW='\033[0;33m'
CYN='\033[0;36m'
RST='\033[0m'

bar() {
    local cur=$1 total=$2 width=30
    local pct=$(( cur * 100 / (total > 0 ? total : 1) ))
    local filled=$(( cur * width / (total > 0 ? total : 1) ))
    [[ $filled -gt $width ]] && filled=$width
    local empty=$(( width - filled ))
    printf "[${GRN}%s${RST}%s] %3d%%" \
        "$(printf '#%.0s' $(seq 1 $filled 2>/dev/null))" \
        "$(printf '.%.0s' $(seq 1 $empty 2>/dev/null))" \
        "$pct"
}

get_state() {
    local jid=$1
    sacct -j "$jid" --format=State --noheader -P 2>/dev/null | head -1 | tr -d ' '
}

get_node() {
    local jid=$1
    sacct -j "$jid" --format=NodeList --noheader -P 2>/dev/null | head -1 | tr -d ' '
}

get_elapsed() {
    local jid=$1
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
    # Extract index from filename, compute sim day
    local idx
    idx=$(basename "$last" .xml | sed 's/output//' | sed 's/^0*//' )
    [[ -z "$idx" ]] && idx=0
    local minutes=$(( 10#$idx * 360 ))
    local day=$(( minutes / 1440 ))
    echo "d${day}"
}

# ── Initial node check ────────────────────────────────────────────
echo -e "\n${CYN}═══════════════════════════════════════════════════════════${RST}"
echo -e "${CYN}  Fix-A Job Monitor — RC1 ($RC1_JOB) + RC2 ($RC2_JOB)${RST}"
echo -e "${CYN}═══════════════════════════════════════════════════════════${RST}"

for JOB in $RC2_JOB $RC1_JOB; do
    NODE=$(get_node "$JOB")
    STATE=$(get_state "$JOB")
    if [[ "$NODE" == *"login"* || "$NODE" == *"head"* ]]; then
        echo -e "  ${RED}⚠  Job $JOB on LOGIN node ($NODE) — DANGER!${RST}"
    elif [[ -z "$NODE" || "$NODE" == "None" ]]; then
        echo -e "  ${YLW}⚠  Job $JOB — node unknown (state=$STATE)${RST}"
    else
        echo -e "  ${GRN}✓  Job $JOB on compute node: $NODE  (state=$STATE)${RST}"
    fi
done
echo ""

# ── Main monitoring loop ──────────────────────────────────────────
rc1_done=0
rc2_done=0

while [[ $rc1_done -eq 0 || $rc2_done -eq 0 ]]; do
    clear 2>/dev/null || true
    echo -e "${CYN}═══════════════════════════════════════════════════════════${RST}"
    echo -e "${CYN}  Fix-A Job Monitor          $(date '+%H:%M:%S')${RST}"
    echo -e "${CYN}═══════════════════════════════════════════════════════════${RST}"

    # ── RC2 ───────────────────────────────────────────────────────
    rc2_state=$(get_state $RC2_JOB)
    rc2_elapsed=$(get_elapsed $RC2_JOB)
    rc2_node=$(get_node $RC2_JOB)
    rc2_snaps=$(snap_count "$RC2_DIR")
    rc2_day=$(latest_day "$RC2_DIR")

    echo -e "\n  ${YLW}RC2${RST} (job $RC2_JOB) — 42-day drug+withdrawal"
    echo -e "  State: ${GRN}${rc2_state}${RST}  Node: ${rc2_node}  Elapsed: ${rc2_elapsed}"
    printf "  Snapshots: %d / %d  (%s)  " "$rc2_snaps" "$RC2_EXPECTED" "$rc2_day"
    bar "$rc2_snaps" "$RC2_EXPECTED"
    echo ""

    if [[ "$rc2_state" == "COMPLETED" ]]; then
        echo -e "  ${GRN}✓ RC2 COMPLETED${RST}"
        rc2_done=1
    elif [[ "$rc2_state" == "FAILED" || "$rc2_state" == "CANCELLED" || "$rc2_state" == "TIMEOUT" ]]; then
        echo -e "  ${RED}✗ RC2 ${rc2_state}${RST}"
        rc2_done=1
    fi

    # ── RC1 ───────────────────────────────────────────────────────
    rc1_state=$(get_state $RC1_JOB)
    rc1_elapsed=$(get_elapsed $RC1_JOB)
    rc1_node=$(get_node $RC1_JOB)
    rc1_snaps=$(snap_count "$RC1_DIR")
    rc1_day=$(latest_day "$RC1_DIR")

    echo -e "\n  ${YLW}RC1${RST} (job $RC1_JOB) — 21-day barrier-only"
    echo -e "  State: ${GRN}${rc1_state}${RST}  Node: ${rc1_node}  Elapsed: ${rc1_elapsed}"
    printf "  Snapshots: %d / %d  (%s)  " "$rc1_snaps" "$RC1_EXPECTED" "$rc1_day"
    bar "$rc1_snaps" "$RC1_EXPECTED"
    echo ""

    if [[ "$rc1_state" == "COMPLETED" ]]; then
        echo -e "  ${GRN}✓ RC1 COMPLETED${RST}"
        rc1_done=1
    elif [[ "$rc1_state" == "FAILED" || "$rc1_state" == "CANCELLED" || "$rc1_state" == "TIMEOUT" ]]; then
        echo -e "  ${RED}✗ RC1 ${rc1_state}${RST}"
        rc1_done=1
    fi

    # ── Key parameters reminder ──────────────────────────────────
    echo -e "\n  ${CYN}Parameters: ecm_emt_cap=0.10, dkc=0.0002, ggp=0.95, emt_death=0.00005${RST}"

    if [[ $rc1_done -eq 1 && $rc2_done -eq 1 ]]; then
        break
    fi

    echo -e "\n  Polling every ${POLL}s... (Ctrl+C to stop)"
    sleep "$POLL"
done

# ── Final summary ────────────────────────────────────────────────
echo -e "\n${CYN}═══════════════════════════════════════════════════════════${RST}"
echo -e "${CYN}  BOTH JOBS FINISHED${RST}"
echo -e "${CYN}═══════════════════════════════════════════════════════════${RST}"
echo -e "  RC2 ($RC2_JOB): $(get_state $RC2_JOB)  —  $rc2_snaps snapshots"
echo -e "  RC1 ($RC1_JOB): $(get_state $RC1_JOB)  —  $rc1_snaps snapshots"
echo ""
echo "  To evaluate:"
echo "    python3 evaluate_rc2.py   (update OUT_DIR to rc2_fixA_seed42)"
echo "    python3 evaluate_rc1.py   (update work_dir to rc1_fixA_seed42)"
echo ""
