#!/bin/bash
# Usage: bash watch_stage2_rc2_single.sh <SLURM_JOB_ID>
set -euo pipefail
JOB_ID="${1:?Usage: bash watch_stage2_rc2_single.sh <JOB_ID>}"
exec bash watch_rc2.sh "$JOB_ID"
