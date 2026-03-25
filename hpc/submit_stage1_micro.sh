#!/bin/bash
#SBATCH --job-name=ns_s1_micro
#SBATCH --partition=cpu384g
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=logs/stage1_micro_%j.out
#SBATCH --error=logs/stage1_micro_%j.err

set -euo pipefail

if command -v module >/dev/null 2>&1; then
  module purge || true
  module load gcc/12 2>/dev/null || true
  module load python/3.10 2>/dev/null || true
fi

PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$PROJECT_ROOT"
mkdir -p logs build

VENV_ACTIVATE="$PROJECT_ROOT/venv/bin/activate"
if [ -f "$VENV_ACTIVATE" ]; then
  # shellcheck disable=SC1090
  source "$VENV_ACTIVATE"
fi

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
export OMP_PROC_BIND="${OMP_PROC_BIND:-spread}"
export OMP_PLACES="${OMP_PLACES:-cores}"
export PYTHONPATH="$PROJECT_ROOT:${PYTHONPATH:-}"

PHYSI_DIR="$PROJECT_ROOT/Stroma_world/PhysiCell"
BUILD_BIN="$PROJECT_ROOT/build/bin/release/stroma_world"

echo "[stage1] build start: $(date)"
make BUILD=release PHYSICELL_DIR="$PHYSI_DIR" release
echo "[stage1] build done:  $(date)"

WORK_DIR="$PROJECT_ROOT/build/stage1_micro_${SLURM_JOB_ID}"

set +e
python3 run_stage1_micro_sim.py \
  --binary "$BUILD_BIN" \
  --base-config "$PROJECT_ROOT/config/PhysiCell_settings.xml" \
  --work-dir "$WORK_DIR" \
  --seed 42 \
  --threads "${SLURM_CPUS_PER_TASK:-16}" \
  --save-interval 360
RC=$?
set -e

echo "[stage1] summary: $WORK_DIR/stage1_summary.json"
echo "[stage1] done rc=$RC: $(date)"
exit $RC
