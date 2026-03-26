#!/usr/bin/env bash
# Watch Phoenix RC2 run and optionally trigger hedge set at completion.
# Usage:
#   bash watch_rc2_phoenix.sh [JOB_ID] [OUT_DIR]
#   bash watch_rc2_phoenix.sh 10674 /work/.../output --auto-hedge

set -uo pipefail

JOB_ID="${1:-10674}"
OUT_DIR="${2:-/work/a0hirs04/PROJECT-NORTHSTAR/build/rc2_regrowth_candidate/replicate_01_seed42/output}"
AUTO_HEDGE="false"
if [[ "${3:-}" == "--auto-hedge" || "${1:-}" == "--auto-hedge" || "${2:-}" == "--auto-hedge" ]]; then
  AUTO_HEDGE="true"
fi

PROJECT_ROOT="/home/a0hirs04/PROJECT-NORTHSTAR"
VERDICT_FILE="$PROJECT_ROOT/logs/rc2_regrowth_candidate_verdict_snippet.txt"
STRICT_FILE="$PROJECT_ROOT/logs/rc2_regrowth_candidate_strict_eval.txt"
POLL=15
SAVE_INTERVAL=360
TOTAL_SNAPS=168
SNAP_PRE=56
SNAP_TREAT=112

phase_of() {
  local s=$1
  if   (( s <= SNAP_PRE ));   then echo "BARRIER"
  elif (( s <= SNAP_TREAT )); then echo "DRUG-ON"
  else                             echo "REGROWTH"
  fi
}

job_info() {
  squeue -j "$JOB_ID" --format="%T|%M|%l|%P|%N|%j|%C|%m" --noheader 2>/dev/null | head -1
}

job_state_sacct() {
  sacct -j "$JOB_ID" --format=State --noheader -P 2>/dev/null | head -1 | tr -d ' '
}

clear
while true; do
  INFO="$(job_info)"
  if [ -n "$INFO" ]; then
    IFS='|' read -r STATE JOB_ELAPSED TIMELIMIT PARTITION NODELIST JOBNAME NCPUS MEM <<< "$INFO"
    STATE=$(echo "$STATE" | tr -d ' ')
  else
    STATE="$(job_state_sacct)"
    JOB_ELAPSED="N/A"; TIMELIMIT="N/A"; PARTITION="?"; NODELIST="?"; JOBNAME="?"; NCPUS="?"; MEM="?"
  fi

  if [ -d "$OUT_DIR" ]; then
    N_SNAPS=$(find "$OUT_DIR" -maxdepth 1 -name 'output*.xml' 2>/dev/null | wc -l)
    DAY=$(echo "scale=1; $N_SNAPS * $SAVE_INTERVAL / 1440" | bc 2>/dev/null || echo "?")
    PHASE=$(phase_of "$N_SNAPS")
    PCT=$(( N_SNAPS * 100 / TOTAL_SNAPS ))
    DISK=$(du -sh "$OUT_DIR" 2>/dev/null | cut -f1)
  else
    N_SNAPS=0
    DAY="?"
    PHASE="?"
    PCT=0
    DISK="?"
  fi

  clear
  echo "RC2 PHOENIX WATCH"
  echo "--------------------------------------------------------------"
  echo "Job: $JOB_ID  Name: $JOBNAME  State: $STATE"
  echo "Partition: $PARTITION  Node: $NODELIST  CPUs: $NCPUS  Mem: $MEM"
  echo "Elapsed: $JOB_ELAPSED  Limit: $TIMELIMIT"
  echo ""
  echo "Progress: day=$DAY  snap=$N_SNAPS/$TOTAL_SNAPS  phase=$PHASE  done=${PCT}%"
  echo "Output dir: $OUT_DIR"
  echo "Disk: $DISK"
  echo ""
  echo "Verdict file: $VERDICT_FILE"
  if [ -s "$VERDICT_FILE" ]; then
    tail -n 6 "$VERDICT_FILE"
  else
    echo "(not ready)"
  fi
  echo ""
  echo "Refresh: ${POLL}s (Ctrl+C to stop)"

  if [ "$STATE" = "COMPLETED" ] || [ "$STATE" = "FAILED" ] || [[ "$STATE" == CANCELLED* ]] || [ "$STATE" = "TIMEOUT" ]; then
    echo ""
    echo "Final job state: $STATE"
    sacct -j "$JOB_ID" --format=JobID,State,Elapsed,NodeList -P 2>/dev/null || true
    echo ""
    if [ -f "$STRICT_FILE" ]; then
      echo "Strict eval tail:"
      tail -n 30 "$STRICT_FILE" || true
    fi

    if [ "$AUTO_HEDGE" = "true" ]; then
      echo ""
      echo "Auto-hedge enabled: invoking launch_rc2_hedge_set.py"
      python3 "$PROJECT_ROOT/launch_rc2_hedge_set.py" --phoenix-job-id "$JOB_ID" | tee -a "$PROJECT_ROOT/logs/rc2_phoenix_autohedge.log"
    fi
    break
  fi

  sleep "$POLL"
done
