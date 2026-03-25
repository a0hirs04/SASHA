#!/bin/bash
# Usage: bash watch_stage1_micro.sh <SLURM_JOB_ID>

set -euo pipefail

JOB_ID="${1:?Usage: bash watch_stage1_micro.sh <JOB_ID>}"
ROOT="/home/a0hirs04/PROJECT-NORTHSTAR"
WORK_DIR="$ROOT/build/stage1_micro_${JOB_ID}"
OUT_LOG="$ROOT/logs/stage1_micro_${JOB_ID}.out"
ERR_LOG="$ROOT/logs/stage1_micro_${JOB_ID}.err"
SUMMARY="$WORK_DIR/stage1_summary.json"

job_state() {
  local jid="$1" st
  st=$(squeue -j "$jid" --format="%T" --noheader 2>/dev/null | head -1 | tr -d ' ')
  if [ -n "$st" ]; then printf "%s" "$st"; return; fi
  st=$(sacct -j "$jid" --format=State --noheader 2>/dev/null | head -1 | tr -d ' ')
  printf "%s" "${st:-UNKNOWN}"
}

while true; do
  clear 2>/dev/null || true
  STATE=$(job_state "$JOB_ID")
  echo "=== Stage 1 Micro Watch $(date) ==="
  echo "Job:      $JOB_ID"
  echo "State:    $STATE"
  echo "Work dir: $WORK_DIR"
  echo

  if [ -f "$SUMMARY" ]; then
    echo "--- stage1_summary.json ---"
    cat "$SUMMARY"
    echo
  else
    echo "Summary not ready yet."
    if [ -f "$OUT_LOG" ]; then
      echo
      echo "--- recent stdout ---"
      tail -n 30 "$OUT_LOG" || true
    fi
    if [ -f "$ERR_LOG" ]; then
      echo
      echo "--- recent stderr ---"
      tail -n 20 "$ERR_LOG" || true
    fi
  fi

  case "$STATE" in
    COMPLETED|FAILED|CANCELLED|TIMEOUT)
      echo
      echo "Final state: $STATE"
      echo "Stdout: $OUT_LOG"
      echo "Stderr: $ERR_LOG"
      break
      ;;
  esac

  sleep 20
done
