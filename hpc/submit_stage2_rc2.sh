#!/bin/bash
#SBATCH --job-name=ns_s2_rc2
#SBATCH --partition=cpu384g
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=06:00:00
#SBATCH --output=logs/stage2_rc2_%j.out
#SBATCH --error=logs/stage2_rc2_%j.err

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

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-32}"
export OMP_PROC_BIND="${OMP_PROC_BIND:-spread}"
export OMP_PLACES="${OMP_PLACES:-cores}"
export PYTHONPATH="$PROJECT_ROOT:${PYTHONPATH:-}"

PHYSI_DIR="$PROJECT_ROOT/Stroma_world/PhysiCell"
BUILD_BIN="$PROJECT_ROOT/build/bin/release/stroma_world"

echo "[stage2] build start: $(date)"
make BUILD=release PHYSICELL_DIR="$PHYSI_DIR" release
echo "[stage2] build done:  $(date)"

# Keep root-level copy in sync for any legacy tooling.
cp "$BUILD_BIN" "$PROJECT_ROOT/stroma_world"

# ── Simulation layout ─────────────────────────────────────────────────────────
SEED=42
WORK_DIR="$PROJECT_ROOT/build/rc2_full_seed42"
RUN_DIR="$WORK_DIR/replicate_01_seed42"
OUT_DIR="$RUN_DIR/output"
CONFIG_SRC="$PROJECT_ROOT/config/PhysiCell_settings.xml"
CONFIG_DST="$RUN_DIR/config.xml"

rm -rf "$WORK_DIR"
mkdir -p "$OUT_DIR"

# Copy auxiliary config files expected by the binary.
LOCAL_CFG="$RUN_DIR/config"
mkdir -p "$LOCAL_CFG"
for f in tumor_calibration_knobs.json gene_params_default.json; do
  [ -f "$PROJECT_ROOT/config/$f" ] && cp "$PROJECT_ROOT/config/$f" "$LOCAL_CFG/$f" || true
done

# Patch XML config: seed, timing, domain, threads, drug schedule.
python3 - "$CONFIG_SRC" "$CONFIG_DST" "$OUT_DIR" "$SEED" \
          "${SLURM_CPUS_PER_TASK:-32}" <<'PYEOF'
import sys, xml.etree.ElementTree as ET

src, dst, out_dir, seed, cpus = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5]

T_PRE      = 20160.0   # day 14  — drug ON
T_TREAT_END= 40320.0   # day 28  — drug OFF / withdrawal
T_POST     = 60480.0   # day 42  — end of regrowth window
DRUG_CONC  = 0.10      # matched to sweep-validated concentration
SAVE_INT   = 360       # 6 h between snapshots

tree = ET.parse(src)
root = tree.getroot()

def _set(xpath, val):
    n = root.find(xpath)
    if n is not None:
        n.text = str(val)

_set("./save/folder",                       out_dir)
_set(".//overall/max_time",                 str(T_POST))
_set(".//options/random_seed",              seed)
_set(".//parallel/omp_num_threads",         cpus)
_set(".//user_parameters/drug_start_time",  str(T_PRE))
_set(".//user_parameters/drug_end_time",    str(T_TREAT_END))
_set(".//user_parameters/drug_concentration", str(DRUG_CONC))

for n in root.findall(".//save//interval"):
    n.text = str(SAVE_INT)

for var_name, val in [("oxygen","38"),("tgfb","0"),("shh","0"),
                      ("drug","0")]:
    var = root.find(f".//microenvironment_setup/variable[@name='{var_name}']")
    if var is None:
        continue
    dbc = var.find("./Dirichlet_boundary_condition")
    if dbc is not None:
        dbc.text = val
        dbc.set("enabled", "true")
    for bv in var.findall("./Dirichlet_options/boundary_value"):
        bv.text = val
        bv.set("enabled", "true")

tree.write(dst, encoding="utf-8", xml_declaration=True)
print(f"[stage2] config written → {dst}")
PYEOF

# ── Run simulation ─────────────────────────────────────────────────────────────
echo "[stage2] simulation start: $(date)"
"$BUILD_BIN" "$CONFIG_DST"
echo "[stage2] simulation done:  $(date)"
echo "[stage2] output: $OUT_DIR"
echo "[stage2] job: ${SLURM_JOB_ID}"
