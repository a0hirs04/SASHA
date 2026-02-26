#!/bin/bash
#SBATCH --job-name=stroma_world_ea
#SBATCH --partition=cpu384g
# On Zurada, switch to another partition if needed (e.g., batch/compute).
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=logs/ea_%j.out
#SBATCH --error=logs/ea_%j.err

set -euo pipefail

# University of Louisville Zurada environment defaults.
if command -v module >/dev/null 2>&1; then
  module purge || true
  module load gcc/12 2>/dev/null || true
fi

PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
export ZURADA_PROJECT_ROOT="${ZURADA_PROJECT_ROOT:-$PROJECT_ROOT}"

# Activate project venv (Python 3.9 — no cluster python module available).
VENV_ACTIVATE="$ZURADA_PROJECT_ROOT/venv/bin/activate"
if [ -f "$VENV_ACTIVATE" ]; then
  # shellcheck source=/dev/null
  source "$VENV_ACTIVATE"
fi
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-32}"
export OMP_PROC_BIND="${OMP_PROC_BIND:-spread}"
export OMP_PLACES="${OMP_PLACES:-cores}"
export PYTHONPATH="$ZURADA_PROJECT_ROOT:${PYTHONPATH:-}"

cd "$ZURADA_PROJECT_ROOT"
mkdir -p logs output

python python/run_ea.py \
  --config config/ea_config_hpc.json \
  --physicell-binary ./stroma_world \
  --physicell-config config/PhysiCell_settings.xml \
  --output-dir output/ea_"$SLURM_JOB_ID" \
  --parallel "${SLURM_CPUS_PER_TASK:-32}"
