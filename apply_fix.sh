#!/bin/bash
set -euo pipefail

ROOT="/home/a0hirs04/PROJECT-NORTHSTAR"
cd "$ROOT"

mkdir -p logs build

echo "[apply_fix] Submitting Stage 1 micro pipeline on HPC compute nodes..."
JOB_ID=$(sbatch --parsable hpc/submit_stage1_micro.sh)
JOB_ID=${JOB_ID%%;*}

echo "[apply_fix] Submitted Stage 1 job: $JOB_ID"
echo "[apply_fix] Watch: bash watch_stage1_micro.sh $JOB_ID"
