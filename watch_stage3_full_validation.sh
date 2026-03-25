#!/bin/bash
# Usage: bash watch_stage3_full_validation.sh <job_id1> [job_id2 ...]

set -euo pipefail
if [ "$#" -lt 1 ]; then
  echo "Usage: bash watch_stage3_full_validation.sh <job_id1> [job_id2 ...]"
  exit 1
fi

job_state() {
  local jid="$1" st
  st=$(squeue -j "$jid" --format="%T" --noheader 2>/dev/null | head -1 | tr -d ' ')
  if [ -n "$st" ]; then printf "%s" "$st"; return; fi
  st=$(sacct -j "$jid" --format=State --noheader 2>/dev/null | head -1 | tr -d ' ')
  printf "%s" "${st:-UNKNOWN}"
}

while true; do
  clear
  echo "=== Stage 3 Full Validation Watch $(date) ==="
  echo
  all_done=true
  for jid in "$@"; do
    st=$(job_state "$jid")
    printf "  %-10s %s\n" "$jid" "$st"
    case "$st" in
      COMPLETED|FAILED|CANCELLED|TIMEOUT|OUT_OF_MEMORY|NODE_FAIL|PREEMPTED|BOOT_FAIL|DEADLINE) ;;
      *) all_done=false ;;
    esac
  done
  echo
  if $all_done; then
    echo "All Stage 3 jobs reached terminal states."
    break
  fi
  sleep 20
done
