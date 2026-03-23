#!/bin/bash
# Watch script for Frankenstein RC2 run (job 10230, seed 42)
# Usage: bash watch_frankenstein.sh
#   or:  watch -n 30 bash watch_frankenstein.sh

OUT_DIR="/work/a0hirs04/PROJECT-NORTHSTAR/build/rc2_frankenstein/replicate_01_seed42/output"
SLURM_LOG="/work/a0hirs04/PROJECT-NORTHSTAR/build/rc2_frankenstein/replicate_01_seed42/slurm_10230.out"
PROJECT="/home/a0hirs04/PROJECT-NORTHSTAR"
TOTAL_SNAPS=168  # 60480 min / 360 min per save

echo "============================================================"
echo "  FRANKENSTEIN RC2 — Live Monitor"
echo "============================================================"

# Job status
JOB_STATE=$(squeue -j 10230 --noheader --format="%T" 2>/dev/null)
if [ -z "$JOB_STATE" ]; then
    JOB_STATE="COMPLETED (or not found)"
fi
echo "  Job 10230:  $JOB_STATE"

# Snapshot count
N_SNAPS=$(ls "$OUT_DIR"/output*.xml 2>/dev/null | wc -l)
PCT=$(( 100 * N_SNAPS / TOTAL_SNAPS ))
BAR_LEN=30
FILLED=$(( BAR_LEN * N_SNAPS / TOTAL_SNAPS ))
EMPTY=$(( BAR_LEN - FILLED ))
BAR=$(printf '%0.s#' $(seq 1 $FILLED 2>/dev/null))$(printf '%0.s-' $(seq 1 $EMPTY 2>/dev/null))
echo "  Snapshots:  $N_SNAPS / $TOTAL_SNAPS  [$BAR] ${PCT}%"

# Estimate phase
if [ "$N_SNAPS" -gt 0 ]; then
    if [ "$N_SNAPS" -le 56 ]; then
        PHASE="Phase 0: Barrier forming (day 0-14)"
    elif [ "$N_SNAPS" -le 112 ]; then
        PHASE="Phase 1: DRUG ON (day 14-28)"
    else
        PHASE="Phase 2: Drug withdrawn / regrowth (day 28-42)"
    fi
    echo "  Phase:      $PHASE"
fi

# Latest snapshot timestamp
if [ "$N_SNAPS" -gt 0 ]; then
    LATEST=$(ls -t "$OUT_DIR"/output*.xml 2>/dev/null | head -1)
    LATEST_AGE=$(( ($(date +%s) - $(stat -c %Y "$LATEST")) ))
    echo "  Latest:     $(basename $LATEST)  (${LATEST_AGE}s ago)"
fi

# Wall time from SLURM log
if [ -f "$SLURM_LOG" ]; then
    START_LINE=$(head -2 "$SLURM_LOG" | grep "Start:")
    if [ -n "$START_LINE" ]; then
        echo "  Started:    ${START_LINE#Start: }"
    fi
fi

echo ""

# --- Cell counts at key timepoints (day 14, 28, 42) ---
# Snap indices: day14=56, day28=112, day42=168
echo "  CELL COUNTS AT KEY TIMEPOINTS"
echo "  -----------------------------------------------------------"
printf "  %-10s  %6s  %6s  %6s  %s\n" "Timepoint" "Tumor" "Stroma" "CAFs" "Status"
printf "  %-10s  %6s  %6s  %6s  %s\n" "----------" "------" "------" "------" "------"

python3 - "$OUT_DIR" "$PROJECT" <<'PYEOF'
import sys, os
out_dir = sys.argv[1]
project = sys.argv[2]
sys.path.insert(0, project)

from python.wrapper.output_parser import OutputParser
import numpy as np

checkpoints = [
    ("Day 14", 56,  "pre-treatment"),
    ("Day 28", 112, "treatment end"),
    ("Day 42", 168, "post-withdrawal"),
]

parser = OutputParser(out_dir)

def _row(matrix, labels, name):
    e = labels.get(name)
    if e is None:
        return None
    idx = int(e["index"])
    if idx < 0 or idx >= matrix.shape[0]:
        return None
    return matrix[idx, :]

for label, snap_idx, desc in checkpoints:
    # Find the snapshot file (or closest earlier one)
    found = False
    for try_idx in range(snap_idx, max(snap_idx - 3, -1), -1):
        fname = os.path.join(out_dir, f"output{try_idx:08d}.xml")
        if os.path.exists(fname) and os.path.getsize(fname) > 0:
            try:
                snap = parser._read_physicell_xml(fname)
                matrix = snap["cell_matrix"]
                labels = snap["label_name_map"]

                cell_type = _row(matrix, labels, "cell_type")
                dead = _row(matrix, labels, "dead")
                death_model = _row(matrix, labels, "current_death_model")
                n_cells = matrix.shape[1]

                live_mask = np.ones(n_cells, dtype=bool)
                if dead is not None:
                    live_mask &= dead <= 0.5
                if death_model is not None:
                    live_mask &= np.rint(death_model).astype(int) != 100

                ctype = np.rint(cell_type).astype(int) if cell_type is not None else np.full(n_cells, -1)
                live_tumor = live_mask & (ctype == 0)
                live_stroma = live_mask & (ctype == 1)

                # CAF count (acta2_active stored at ZEB1 label offset for stromal cells)
                acta2 = _row(matrix, labels, "ZEB1")
                live_caf = live_stroma & (acta2 > 0.5) if acta2 is not None else np.zeros_like(live_stroma)

                n_tumor = int(np.sum(live_tumor))
                n_stroma = int(np.sum(live_stroma))
                n_caf = int(np.sum(live_caf))

                print(f"  {label:<10s}  {n_tumor:6d}  {n_stroma:6d}  {n_caf:6d}  {desc}")
                found = True
                break
            except Exception:
                pass
    if not found:
        print(f"  {label:<10s}  {'--':>6s}  {'--':>6s}  {'--':>6s}  (not yet available)")

# Quick RC2 criteria preview if day 28 and day 14 are both available
try:
    snap14 = None
    snap28 = None
    snap42 = None
    for try_idx in range(56, 53, -1):
        f = os.path.join(out_dir, f"output{try_idx:08d}.xml")
        if os.path.exists(f) and os.path.getsize(f) > 0:
            snap14 = parser._read_physicell_xml(f)
            break
    for try_idx in range(112, 109, -1):
        f = os.path.join(out_dir, f"output{try_idx:08d}.xml")
        if os.path.exists(f) and os.path.getsize(f) > 0:
            snap28 = parser._read_physicell_xml(f)
            break
    for try_idx in range(168, 165, -1):
        f = os.path.join(out_dir, f"output{try_idx:08d}.xml")
        if os.path.exists(f) and os.path.getsize(f) > 0:
            snap42 = parser._read_physicell_xml(f)
            break

    if snap14 and snap28:
        mat14 = snap14["cell_matrix"]
        lbl14 = snap14["label_name_map"]
        ct14 = _row(mat14, lbl14, "cell_type")
        d14 = _row(mat14, lbl14, "dead")
        dm14 = _row(mat14, lbl14, "current_death_model")
        n14 = mat14.shape[1]
        lm14 = np.ones(n14, dtype=bool)
        if d14 is not None: lm14 &= d14 <= 0.5
        if dm14 is not None: lm14 &= np.rint(dm14).astype(int) != 100
        c14 = np.rint(ct14).astype(int) if ct14 is not None else np.full(n14, -1)
        n_pre = int(np.sum(lm14 & (c14 == 0)))

        mat28 = snap28["cell_matrix"]
        lbl28 = snap28["label_name_map"]
        ct28 = _row(mat28, lbl28, "cell_type")
        d28 = _row(mat28, lbl28, "dead")
        dm28 = _row(mat28, lbl28, "current_death_model")
        n28 = mat28.shape[1]
        lm28 = np.ones(n28, dtype=bool)
        if d28 is not None: lm28 &= d28 <= 0.5
        if dm28 is not None: lm28 &= np.rint(dm28).astype(int) != 100
        c28 = np.rint(ct28).astype(int) if ct28 is not None else np.full(n28, -1)
        n_end = int(np.sum(lm28 & (c28 == 0)))

        print("")
        print("  EARLY CRITERIA CHECK")
        print("  -----------------------------------------------------------")
        reduction = 1.0 - n_end / max(n_pre, 1)
        tag = "PASS" if reduction >= 0.10 else "FAIL"
        print(f"  RC2-1 Partial response: {n_pre} -> {n_end} ({reduction:.1%} reduction) [{tag}]")
        tag = "PASS" if n_end > 0 else "FAIL"
        print(f"  RC2-2 Not eradicated:   {n_end} survivors [{tag}]")

        if snap42:
            mat42 = snap42["cell_matrix"]
            lbl42 = snap42["label_name_map"]
            ct42 = _row(mat42, lbl42, "cell_type")
            d42 = _row(mat42, lbl42, "dead")
            dm42 = _row(mat42, lbl42, "current_death_model")
            n42c = mat42.shape[1]
            lm42 = np.ones(n42c, dtype=bool)
            if d42 is not None: lm42 &= d42 <= 0.5
            if dm42 is not None: lm42 &= np.rint(dm42).astype(int) != 100
            c42 = np.rint(ct42).astype(int) if ct42 is not None else np.full(n42c, -1)
            n_post = int(np.sum(lm42 & (c42 == 0)))
            tag = "PASS" if n_post > n_end else "FAIL"
            print(f"  RC2-6 Regrowth:         {n_end} -> {n_post} [{tag}]")
except Exception:
    pass
PYEOF

echo ""

# If job is done, prompt for full evaluation
if [ "$JOB_STATE" = "COMPLETED (or not found)" ] && [ "$N_SNAPS" -ge 160 ]; then
    echo "  >>> Job finished! Run full evaluation with:"
    echo "  python3 $PROJECT/evaluate_rc2.py \\"
    echo "      --out-dir $OUT_DIR"
    echo ""

    if [ -f "$SLURM_LOG" ]; then
        echo "  SLURM log tail:"
        tail -5 "$SLURM_LOG"
    fi
fi

echo "============================================================"
