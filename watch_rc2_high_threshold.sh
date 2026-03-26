#!/usr/bin/env bash
# Focused polling protocol for high-threshold hedge jobs: C3, C5, C7.
# - Monitors terminal completion
# - Runs strict evaluator immediately on completion
# - Writes strict report + compact verdict snippet
#
# Usage:
#   bash watch_rc2_high_threshold.sh
#   bash watch_rc2_high_threshold.sh --once

set -uo pipefail

PROJECT_ROOT="/home/a0hirs04/PROJECT-NORTHSTAR"
WORK_ROOT="$PROJECT_ROOT/build/rc2_hedge_set"
if [[ ! -d "$WORK_ROOT" && -d "/work/a0hirs04/PROJECT-NORTHSTAR/build/rc2_hedge_set" ]]; then
  WORK_ROOT="/work/a0hirs04/PROJECT-NORTHSTAR/build/rc2_hedge_set"
fi
LOG_ROOT="$PROJECT_ROOT/logs"
POLL=20
ONCE="false"

if [[ "${1:-}" == "--once" ]]; then
  ONCE="true"
fi

JOBS=(
  "C3|rc2_hedge_C3_km2.0_ap0.00075_th0.03"
  "C5|rc2_hedge_C5_km1.9_ap0.0007_th0.03"
  "C7|rc2_hedge_C7_km2.2_ap0.00065_th0.03"
)

state_for_prefix() {
  local prefix="$1"

  local sq
  sq=$(squeue -h -o "%j|%T" 2>/dev/null | awk -F'|' -v p="$prefix" '$1 ~ ("^"p"_r") {print $2; exit}')
  if [[ -n "$sq" ]]; then
    echo "$sq"
    return
  fi

  local sa
  sa=$(sacct -X -n -P --format=JobName,State 2>/dev/null | awk -F'|' -v p="$prefix" '$1 ~ ("^"p"_r") {print $2; exit}')
  if [[ -n "$sa" ]]; then
    echo "$sa"
  else
    echo "UNKNOWN"
  fi
}

is_terminal() {
  local s="$1"
  [[ "$s" == "COMPLETED" || "$s" == "FAILED" || "$s" == TIMEOUT* || "$s" == CANCELLED* || "$s" == OUT_OF_MEMORY* || "$s" == NODE_FAIL* ]]
}

evaluate_if_ready() {
  local cid="$1"
  local prefix="$2"

  local out_dir="$WORK_ROOT/$prefix/replicate_01_seed42/output"
  local strict_file="$LOG_ROOT/${prefix}_strict_eval.txt"
  local snippet_file="$LOG_ROOT/${prefix}_verdict_snippet.txt"
  local done_flag="$LOG_ROOT/.${prefix}.strict.done"

  [[ -f "$done_flag" ]] && return 0
  [[ -d "$out_dir" ]] || return 0

  local n_xml
  n_xml=$(find "$out_dir" -maxdepth 1 -name 'output*.xml' 2>/dev/null | wc -l)
  if [[ "$n_xml" -lt 5 ]]; then
    return 0
  fi

  python3 "$PROJECT_ROOT/evaluate_rc2.py" --out-dir "$out_dir" > "$strict_file" 2>&1 || true

  {
    echo "RC2 strict verdict :: $cid"
    echo "job_prefix=$prefix"
    echo "=================================================="
    grep -E 'RC2-6 Regrowth:|Hard Score:|Classification:|Reason:|Risk Flags:|NO_EFFECTIVE_TREATMENT|seed_42\s+hard=' "$strict_file" || true
  } > "$snippet_file"

  touch "$done_flag"

  echo "[STRICT READY] $cid -> $snippet_file"

  if grep -q 'RC2-6 Regrowth: PASS' "$strict_file"; then
    if grep -Eq 'Risk Flags:\s+NONE|Risk Flags:\s+LOW_PRE_BARRIER' "$strict_file"; then
      echo "[SUCCESS SIGNAL] $cid meets requested signal (RC2-6 PASS + acceptable risk flags)"
    fi
  fi
}

all_done() {
  local done=0
  local total=0
  while IFS='|' read -r cid prefix; do
    [[ -z "$cid" ]] && continue
    total=$((total+1))
    [[ -f "$LOG_ROOT/.${prefix}.strict.done" ]] && done=$((done+1))
  done < <(printf "%s\n" "${JOBS[@]}")
  [[ "$total" -gt 0 && "$done" -eq "$total" ]]
}

while true; do
  clear
  echo "RC2 HIGH-THRESHOLD POLLING (C3/C5/C7)"
  echo "--------------------------------------------------------------------------"
  echo "Time: $(date '+%F %T')"
  echo "Columns: ID | prefix | state | rep1_snaps | strict_status"
  echo "--------------------------------------------------------------------------"

  while IFS='|' read -r cid prefix; do
    [[ -z "$cid" ]] && continue
    state=$(state_for_prefix "$prefix")

    out_dir="$WORK_ROOT/$prefix/replicate_01_seed42/output"
    snaps="-"
    if [[ -d "$out_dir" ]]; then
      snaps=$(find "$out_dir" -maxdepth 1 -name 'output*.xml' 2>/dev/null | wc -l)
    fi

    strict_status="PENDING"
    if [[ -f "$LOG_ROOT/.${prefix}.strict.done" ]]; then
      strict_status="READY"
    fi

    printf "%-3s %-42s %-12s %-9s %-8s\n" "$cid" "$prefix" "$state" "$snaps" "$strict_status"

    if is_terminal "$state"; then
      evaluate_if_ready "$cid" "$prefix"
    fi
  done < <(printf "%s\n" "${JOBS[@]}")

  echo "--------------------------------------------------------------------------"
  echo "Snippets: $LOG_ROOT/*_verdict_snippet.txt"

  if all_done; then
    echo "All C3/C5/C7 strict verdicts are ready."
    break
  fi

  if [[ "$ONCE" == "true" ]]; then
    break
  fi

  echo "Refresh: ${POLL}s (Ctrl+C to stop)"
  sleep "$POLL"
done
