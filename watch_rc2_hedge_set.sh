#!/usr/bin/env bash
# Dashboard watcher for RC2 hedge set C1-C7.
# Usage:
#   bash watch_rc2_hedge_set.sh
#   bash watch_rc2_hedge_set.sh --once

set -uo pipefail

PROJECT_ROOT="/home/a0hirs04/PROJECT-NORTHSTAR"
LOG_DIR="$PROJECT_ROOT/logs"
BUILD_ROOT="$PROJECT_ROOT/build/rc2_hedge_set"
if [[ ! -d "$BUILD_ROOT" && -d "/work/a0hirs04/PROJECT-NORTHSTAR/build/rc2_hedge_set" ]]; then
  BUILD_ROOT="/work/a0hirs04/PROJECT-NORTHSTAR/build/rc2_hedge_set"
fi
POLL=20
ONCE="false"

if [[ "${1:-}" == "--once" ]]; then
  ONCE="true"
fi

JOBS=(
  "C1|rc2_hedge_C1_km2.1_ap0.0008_th0.025"
  "C2|rc2_hedge_C2_km2.2_ap0.00075_th0.025"
  "C3|rc2_hedge_C3_km2.0_ap0.00075_th0.03"
  "C4|rc2_hedge_C4_km2.3_ap0.0008_th0.025"
  "C5|rc2_hedge_C5_km1.9_ap0.0007_th0.03"
  "C6|rc2_hedge_C6_km2.4_ap0.00075_th0.02"
  "C7|rc2_hedge_C7_km2.2_ap0.00065_th0.03"
)

job_line_state() {
  local prefix="$1"
  local state="NOT_LAUNCHED"
  local reps="-"
  local snaps="-"

  local log_file="$LOG_DIR/${prefix}.log"
  local rep1_out="$BUILD_ROOT/${prefix}/replicate_01_seed42/output"

  # Active SLURM replicas by job name prefix
  local sq
  sq=$(squeue -h -o "%j|%T|%M|%N" 2>/dev/null | awk -F'|' -v p="$prefix" '$1 ~ ("^"p"_r")')

  if [[ -n "$sq" ]]; then
    state="RUNNING"
    reps=$(echo "$sq" | wc -l)
  fi

  if [[ -d "$rep1_out" ]]; then
    snaps=$(find "$rep1_out" -maxdepth 1 -name 'output*.xml' 2>/dev/null | wc -l)
  fi

  if [[ -f "$log_file" ]]; then
    if grep -q "\*\*\* REALITY CHECK 2: PASS \*\*\*" "$log_file"; then
      state="PASS"
    elif grep -q "\*\*\* REALITY CHECK 2: FAIL \*\*\*" "$log_file"; then
      state="FAIL"
    elif [[ "$state" != "RUNNING" ]]; then
      if grep -q "Submitted rep" "$log_file"; then
        state="SUBMITTED"
      else
        state="LOG_PRESENT"
      fi
    fi
  fi

  printf "%-3s %-46s %-12s %-4s %-5s\n" "$cid" "$prefix" "$state" "$reps" "$snaps"
}

all_terminal() {
  local n_term=0
  local n_all=0
  while IFS='|' read -r cid prefix; do
    [[ -z "$cid" ]] && continue
    n_all=$((n_all+1))
    local log_file="$LOG_DIR/${prefix}.log"
    if [[ -f "$log_file" ]] && (grep -q "\*\*\* REALITY CHECK 2: PASS \*\*\*" "$log_file" || grep -q "\*\*\* REALITY CHECK 2: FAIL \*\*\*" "$log_file"); then
      n_term=$((n_term+1))
    fi
  done < <(printf "%s\n" "${JOBS[@]}")

  [[ "$n_all" -gt 0 && "$n_term" -eq "$n_all" ]]
}

while true; do
  clear
  echo "RC2 HEDGE SET WATCH (C1-C7)"
  echo "--------------------------------------------------------------------------"
  echo "Time: $(date '+%F %T')"
  echo "Columns: ID | job_prefix | state | active_reps | rep1_snaps"
  echo "--------------------------------------------------------------------------"

  pass_n=0
  fail_n=0
  run_n=0
  launch_n=0

  while IFS='|' read -r cid prefix; do
    [[ -z "$cid" ]] && continue
    line=$(cid="$cid" job_line_state "$prefix")
    echo "$line"

    st=$(echo "$line" | awk '{print $3}')
    case "$st" in
      PASS) pass_n=$((pass_n+1));;
      FAIL) fail_n=$((fail_n+1));;
      RUNNING|SUBMITTED) run_n=$((run_n+1));;
      NOT_LAUNCHED|LOG_PRESENT) launch_n=$((launch_n+1));;
    esac
  done < <(printf "%s\n" "${JOBS[@]}")

  echo "--------------------------------------------------------------------------"
  echo "Summary: PASS=$pass_n FAIL=$fail_n ACTIVE=$run_n PENDING_LAUNCH=$launch_n"
  echo "Logs: $LOG_DIR/rc2_hedge_*.log"
  echo "Build root: $BUILD_ROOT"

  if all_terminal; then
    echo "All hedge jobs reached terminal verdicts (PASS/FAIL)."
    break
  fi

  if [[ "$ONCE" == "true" ]]; then
    break
  fi

  echo "Refresh: ${POLL}s (Ctrl+C to stop)"
  sleep "$POLL"
done
