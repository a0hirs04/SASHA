#!/usr/bin/env bash
# Rich watch dashboard for promoted RC2 winner run
# Usage: bash watch_rc2_transition_winner.sh [JOB_ID] [OUT_DIR]

set -uo pipefail

JOB_ID="${1:-10666}"
OUT_DIR="${2:-/work/a0hirs04/PROJECT-NORTHSTAR/build/rc2_full_transition_winner/replicate_01_seed42/output}"
PROJECT_ROOT="/home/a0hirs04/PROJECT-NORTHSTAR"

SAVE_INTERVAL=360
TOTAL_SNAPS=168
SNAP_PRE=56
SNAP_TREAT=112
POLL=10

if [ ! -d "$OUT_DIR" ]; then
  echo "ERROR: output dir not found: $OUT_DIR"
  exit 1
fi

phase_of() {
  local s=$1
  if   (( s <= SNAP_PRE ));   then echo "BARRIER"
  elif (( s <= SNAP_TREAT )); then echo "DRUG-ON"
  else                              echo "REGROWTH"
  fi
}

job_info() {
  squeue -j "$JOB_ID" --format="%T|%M|%l|%P|%N|%j|%C|%m" --noheader 2>/dev/null | head -1
}

job_state_sacct() {
  sacct -j "$JOB_ID" --format=State --noheader 2>/dev/null | head -1 | tr -d ' '
}

fmt_hms() {
  local s=$1
  printf "%02d:%02d:%02d" $((s/3600)) $(((s%3600)/60)) $((s%60))
}

snapshot_metrics() {
  python3 - <<PY 2>/dev/null || echo "?|?|?|?|?|?|?"
import sys, numpy as np
sys.path.insert(0, "$PROJECT_ROOT")
from python.wrapper.output_parser import OutputParser

p = OutputParser("$OUT_DIR")
s = p._read_physicell_xml("$1")
m = s["cell_matrix"]
l = s["label_name_map"]

def row(name):
    e = l.get(name)
    if e is None:
        return None
    i = int(e["index"])
    if i < 0 or i >= m.shape[0]:
        return None
    return m[i,:]

ct = row("cell_type")
dead = row("dead")
dm = row("current_death_model")
ab = row("abcb1_active")
if ab is None:
    ab = row("ABCB1")
ic = row("intracellular_drug")

n = m.shape[1]
ctype = np.rint(ct).astype(int) if ct is not None else np.full(n, -1)
tumor = ctype == 0
stroma = ctype == 1
live = np.ones(n, dtype=bool)
if dead is not None:
    live &= dead <= 0.5
if dm is not None:
    live &= np.rint(dm).astype(int) != 100

lt = live & tumor
ls = live & stroma
n_live_t = int(np.sum(lt))
n_dead_t = int(np.sum(tumor & ~live))
n_live_s = int(np.sum(ls))
n_total = int(n)

if n_live_t > 0:
    ab_pct = float(np.mean(ab[lt] > 0.5) * 100.0) if ab is not None else float('nan')
    ab_mean = float(np.mean(ab[lt])) if ab is not None else float('nan')
    ic_mean = float(np.mean(ic[lt])) if ic is not None else float('nan')
else:
    ab_pct, ab_mean, ic_mean = float('nan'), float('nan'), float('nan')

print(f"{n_live_t}|{n_dead_t}|{n_live_s}|{n_total}|{ab_pct:.2f}|{ab_mean:.4f}|{ic_mean:.6f}")
PY
}

START_TIME=$(date +%s)
LAST_XML=""
TUMOR_LIVE="?"
TUMOR_DEAD="?"
STROMA_LIVE="?"
TOTAL_CELLS="?"
ABCB1_PCT="?"
ABCB1_MEAN="?"
IC_MEAN="?"

clear
echo "RC2 TRANSITION WINNER WATCH"
echo "job: $JOB_ID"
echo "out: $OUT_DIR"
echo ""

while true; do
  INFO=$(job_info)
  if [ -n "$INFO" ]; then
    IFS='|' read -r STATE JOB_ELAPSED TIMELIMIT PARTITION NODELIST JOBNAME NCPUS MEM <<< "$INFO"
    STATE=$(echo "$STATE" | tr -d ' ')
  else
    STATE=$(job_state_sacct)
    JOB_ELAPSED="N/A"; TIMELIMIT="N/A"; PARTITION="?"; NODELIST="?"; JOBNAME="?"; NCPUS="?"; MEM="?"
  fi

  N_SNAPS=$(find "$OUT_DIR" -maxdepth 1 -name 'output*.xml' 2>/dev/null | wc -l)
  DAY=$(echo "scale=1; $N_SNAPS * $SAVE_INTERVAL / 1440" | bc 2>/dev/null || echo "?")
  PHASE=$(phase_of "$N_SNAPS")
  PCT=$(( N_SNAPS * 100 / TOTAL_SNAPS ))

  NOW=$(date +%s)
  ELAPSED=$(( NOW - START_TIME ))
  if (( N_SNAPS > 3 && ELAPSED > 0 )); then
    ETA=$(( ELAPSED * (TOTAL_SNAPS - N_SNAPS) / N_SNAPS ))
    SPEED=$(echo "scale=2; $N_SNAPS * 60 / $ELAPSED" | bc 2>/dev/null || echo "?")
  else
    ETA=0
    SPEED="?"
  fi

  LATEST_XML=$(find "$OUT_DIR" -maxdepth 1 -name 'output*.xml' 2>/dev/null | sort | tail -1)
  if [ -n "${LATEST_XML:-}" ] && [ "$LATEST_XML" != "$LAST_XML" ]; then
    MET=$(snapshot_metrics "$LATEST_XML")
    IFS='|' read -r TUMOR_LIVE TUMOR_DEAD STROMA_LIVE TOTAL_CELLS ABCB1_PCT ABCB1_MEAN IC_MEAN <<< "$MET"
    LAST_XML="$LATEST_XML"
  fi

  DISK=$(du -sh "$OUT_DIR" 2>/dev/null | cut -f1)

  clear
  echo "RC2 TRANSITION WINNER WATCH"
  echo "--------------------------------------------------------------"
  echo "Job: $JOB_ID  Name: $JOBNAME  State: $STATE"
  echo "Partition: $PARTITION  Node: $NODELIST  CPUs: $NCPUS  Mem: $MEM"
  echo "Elapsed: $JOB_ELAPSED  Limit: $TIMELIMIT"
  echo ""
  echo "Progress: day=$DAY  snap=$N_SNAPS/$TOTAL_SNAPS  phase=$PHASE  done=${PCT}%"
  echo "Rate: $SPEED snap/min  ETA: $(fmt_hms "$ETA")  Disk: ${DISK:-?}"
  echo ""
  echo "Tumor live/dead: $TUMOR_LIVE / $TUMOR_DEAD"
  echo "Stroma live: $STROMA_LIVE  Total cells: $TOTAL_CELLS"
  echo "ABCB1+%: $ABCB1_PCT  mean_ABCB1: $ABCB1_MEAN  ic_mean: $IC_MEAN"
  echo "Latest XML: ${LATEST_XML:-none}"
  echo ""
  echo "Refresh: ${POLL}s   (Ctrl+C to detach)"

  if [ "$STATE" = "COMPLETED" ] || [ "$STATE" = "FAILED" ] || [[ "$STATE" == CANCELLED* ]] || [ "$STATE" = "TIMEOUT" ]; then
    echo ""
    echo "Final job state: $STATE"
    sacct -j "$JOB_ID" --format=JobID,State,Elapsed,NodeList -P 2>/dev/null || true
    echo ""
    echo "Run strict evaluator:"
    echo "  python3 /home/a0hirs04/PROJECT-NORTHSTAR/evaluate_rc2.py --out-dir $OUT_DIR"
    break
  fi

  sleep "$POLL"
done
