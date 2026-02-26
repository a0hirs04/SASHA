#!/bin/bash
#SBATCH --job-name=stroma_world_ea_array
#SBATCH --partition=cpu384g
# On Zurada, switch to another partition if needed (e.g., batch/compute).
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --array=1-50%50
#SBATCH --output=logs/ea_array_%A_%a.out
#SBATCH --error=logs/ea_array_%A_%a.err

set -euo pipefail

# University of Louisville Zurada environment defaults.
if command -v module >/dev/null 2>&1; then
  module purge || true
  module load gcc/12
  module load python/3.10
fi

PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
export ZURADA_PROJECT_ROOT="${ZURADA_PROJECT_ROOT:-$PROJECT_ROOT}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-32}"
export OMP_PROC_BIND="${OMP_PROC_BIND:-spread}"
export OMP_PLACES="${OMP_PLACES:-cores}"
export PYTHONPATH="$ZURADA_PROJECT_ROOT:${PYTHONPATH:-}"

cd "$ZURADA_PROJECT_ROOT"
mkdir -p logs output/ea_array_configs

SEED="${SLURM_ARRAY_TASK_ID}"
BASE_CONFIG="config/ea_config_hpc.json"
ARRAY_CONFIG="output/ea_array_configs/ea_config_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.json"
RUN_OUTPUT="output/ea_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"

python - "$BASE_CONFIG" "$ARRAY_CONFIG" "$SEED" <<'PY'
import json
import sys
from pathlib import Path

base_cfg = Path(sys.argv[1])
out_cfg = Path(sys.argv[2])
seed = int(sys.argv[3])

payload = {}
if base_cfg.exists():
    with base_cfg.open("r", encoding="utf-8") as f:
        payload = json.load(f)
if not isinstance(payload, dict):
    payload = {}

payload["random_seed"] = seed
out_cfg.parent.mkdir(parents=True, exist_ok=True)
with out_cfg.open("w", encoding="utf-8") as f:
    json.dump(payload, f, indent=2)
    f.write("\n")
PY

python python/run_ea.py \
  --config "$ARRAY_CONFIG" \
  --physicell-binary ./stroma_world \
  --physicell-config config/PhysiCell_settings.xml \
  --output-dir "$RUN_OUTPUT" \
  --parallel "${SLURM_CPUS_PER_TASK:-32}" \
  --seed "$SEED"
