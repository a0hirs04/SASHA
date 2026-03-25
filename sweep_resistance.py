#!/usr/bin/env python3
"""
sweep_resistance.py — Parallel parameter sweep for drug resistance tuning.

Launches 4 SLURM Stage 1 micro-sim jobs (4 kill multipliers x 1 production rate).
All jobs run simultaneously on HPC compute nodes.

Usage:
    python3 sweep_resistance.py
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent
BINARY = PROJECT_ROOT / "stroma_world"
BASE_CONFIG = PROJECT_ROOT / "config" / "PhysiCell_settings.xml"

# Output to /work for large HPC scratch
SWEEP_DIR = Path("/work/a0hirs04/PROJECT-NORTHSTAR/build/resistance_sweep")

SEED = 42
PARTITION = "cpu384g"
CPUS = 16
MEM = "64G"
WALL_TIME = "02:00:00"

# Stage 1 micro-sim timing
DAY_MIN = 1440.0
T_DRUG_START = 14.0 * DAY_MIN  # 20160
T_WITHDRAW = 28.0 * DAY_MIN    # 40320
T_END = 31.0 * DAY_MIN         # 44640
SAVE_INTERVAL = 360

# Sweep grid — Round 3: transition zone
# Round 1 (km=2-10, ap=0.0005): total kill — drug too fast for resistance
# Round 2 (km=0.1-1.0, ap=0.001-0.01): total survival — 100% ABCB1, ic=0
# Round 3: narrow the transition zone between kill and survival
KILL_MULTIPLIERS = [2.0, 2.2, 2.5, 2.8]
ABCB1_PRODUCTION_RATES = [0.0008]
FIXED_STRESS_THRESHOLD = 0.02


def _set(root: ET.Element, xpath: str, value: str) -> None:
    node = root.find(xpath)
    if node is not None:
        node.text = value


def _set_domain(root: ET.Element) -> None:
    _set(root, ".//domain/x_min", "-200")
    _set(root, ".//domain/x_max", "200")
    _set(root, ".//domain/y_min", "-200")
    _set(root, ".//domain/y_max", "200")
    _set(root, ".//domain/z_min", "-10")
    _set(root, ".//domain/z_max", "10")
    _set(root, ".//domain/dx", "20")
    _set(root, ".//domain/dy", "20")
    _set(root, ".//domain/dz", "20")
    _set(root, ".//domain/use_2D", "true")


def patch_config(
    src: Path,
    dst: Path,
    out_dir: Path,
    kill_multiplier: float,
    abcb1_production_rate: float,
) -> None:
    tree = ET.parse(src)
    root = tree.getroot()

    # Stage 1 micro-sim geometry/runtime (from run_stage1_micro_sim.py)
    _set_domain(root)
    _set(root, "./save/folder", str(out_dir))
    _set(root, ".//overall/max_time", str(T_END))
    _set(root, ".//overall/dt_diffusion", "0.1")
    _set(root, ".//overall/dt_mechanics", "0.5")
    _set(root, ".//overall/dt_phenotype", "6")
    _set(root, ".//options/random_seed", str(SEED))
    _set(root, ".//parallel/omp_num_threads", str(CPUS))

    # Disable SVG, keep full snapshots
    _set(root, ".//save/SVG/enable", "false")
    _set(root, ".//save/legacy_data/enable", "false")
    for n in root.findall(".//save//interval"):
        n.text = str(SAVE_INTERVAL)

    # Drug schedule
    _set(root, ".//user_parameters/drug_start_time", str(T_DRUG_START))
    _set(root, ".//user_parameters/drug_end_time", str(T_WITHDRAW))
    _set(root, ".//user_parameters/drug_concentration", "0.10")

    # Small cell counts
    _set(root, ".//user_parameters/number_of_tumor_cells", "60")
    _set(root, ".//user_parameters/number_of_stromal_cells", "20")
    _set(root, ".//user_parameters/tumor_cluster_radius", "60.0")
    _set(root, ".//user_parameters/stroma_inner_radius", "80.0")
    _set(root, ".//user_parameters/stroma_outer_radius", "120.0")

    # Deterministic boundaries + high-ECM
    for var_name, value in [
        ("oxygen", "38"),
        ("tgfb", "0"),
        ("shh", "0"),
        ("drug", "0"),
        ("ecm_density", "0.8"),
    ]:
        var = root.find(f".//microenvironment_setup/variable[@name='{var_name}']")
        if var is None:
            continue
        init = var.find("./initial_condition")
        if init is not None:
            init.text = value
        dbc = var.find("./Dirichlet_boundary_condition")
        if dbc is not None:
            dbc.text = value
            dbc.set("enabled", "true")
        for bv in var.findall("./Dirichlet_options/boundary_value"):
            bv.text = value
            bv.set("enabled", "true")

    # === SWEEP PARAMETERS ===
    _set(root, ".//user_parameters/drug_stress_threshold", str(FIXED_STRESS_THRESHOLD))
    _set(root, ".//user_parameters/drug_kill_multiplier", str(kill_multiplier))
    _set(root, ".//user_parameters/abcb1_production_rate", str(abcb1_production_rate))

    tree.write(dst, encoding="utf-8", xml_declaration=True)


def make_slurm_script(variant_name: str, variant_dir: Path, config_path: Path) -> str:
    return f"""#!/bin/bash
#SBATCH --job-name=rs_{variant_name}
#SBATCH --partition={PARTITION}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task={CPUS}
#SBATCH --mem={MEM}
#SBATCH --time={WALL_TIME}
#SBATCH --output={variant_dir}/slurm_%j.out
#SBATCH --error={variant_dir}/slurm_%j.err

set -euo pipefail
export OMP_NUM_THREADS={CPUS}
export OMP_PROC_BIND=spread
export OMP_PLACES=cores

cd {PROJECT_ROOT}
echo "=== rs_{variant_name} seed={SEED} ==="
echo "drug_stress_threshold=$(grep -oP 'drug_stress_threshold[^>]*>\\K[^<]+' {config_path})"
echo "drug_kill_multiplier=$(grep -oP 'drug_kill_multiplier[^>]*>\\K[^<]+' {config_path})"
echo "Start: $(date)"
{BINARY} {config_path}
echo "End:   $(date)"
"""


def main() -> int:
    if not BINARY.exists():
        print(f"ERROR: binary not found: {BINARY}")
        print("Run 'make BUILD=release' first.")
        return 1
    if not BASE_CONFIG.exists():
        print(f"ERROR: base config not found: {BASE_CONFIG}")
        return 1

    # Clean and create sweep directory
    if SWEEP_DIR.exists():
        shutil.rmtree(SWEEP_DIR)
    SWEEP_DIR.mkdir(parents=True, exist_ok=True)

    manifest = {"sweep_params": ["drug_kill_multiplier", "abcb1_production_rate"], "variants": {}}
    job_ids = []

    for km in KILL_MULTIPLIERS:
        for ap in ABCB1_PRODUCTION_RATES:
            variant_name = f"km{km}_ap{ap}"
            variant_dir = SWEEP_DIR / variant_name
            out_dir = variant_dir / "output"
            out_dir.mkdir(parents=True, exist_ok=True)

            # Copy auxiliary config files
            local_cfg = variant_dir / "config"
            local_cfg.mkdir(exist_ok=True)
            for extra in ["tumor_calibration_knobs.json", "gene_params_default.json"]:
                src = PROJECT_ROOT / "config" / extra
                if src.exists():
                    shutil.copy(src, local_cfg / extra)

            # Patch config
            config_path = variant_dir / "config.xml"
            patch_config(BASE_CONFIG, config_path, out_dir, km, ap)

            # Write SLURM script
            slurm_path = variant_dir / "run.slurm.sh"
            slurm_path.write_text(make_slurm_script(variant_name, variant_dir, config_path))
            slurm_path.chmod(0o755)

            # Submit
            result = subprocess.run(
                ["sbatch", "--parsable", str(slurm_path)],
                capture_output=True, text=True,
                cwd=str(PROJECT_ROOT),
            )
            if result.returncode != 0:
                print(f"ERROR submitting {variant_name}: {result.stderr.strip()}")
                return 1

            job_id = result.stdout.strip().split(";")[0]
            job_ids.append((variant_name, job_id))

            manifest["variants"][variant_name] = {
                "drug_kill_multiplier": km,
                "abcb1_production_rate": ap,
                "drug_stress_threshold": FIXED_STRESS_THRESHOLD,
                "job_id": job_id,
                "work_dir": str(variant_dir),
            }

            print(f"  {variant_name:20s}  job={job_id}  km={km}  ap={ap}")

    # Write manifest
    manifest_path = SWEEP_DIR / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2))

    print(f"\n{'=' * 60}")
    print(f"  {len(job_ids)} jobs submitted to {PARTITION}")
    print(f"  Manifest: {manifest_path}")
    print(f"  Monitor:  squeue -u $USER")
    print(f"  Evaluate: python3 evaluate_resistance_sweep.py")
    print(f"{'=' * 60}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
