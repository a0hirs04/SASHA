#!/usr/bin/env python3
"""
Full RC2 run — seed 42 only, 42-day schedule.
Phase 0: day 0-14  (barrier)
Phase 1: day 14-28 (drug ON)
Phase 2: day 28-42 (drug OFF / regrowth)
"""
from __future__ import annotations

import shutil
import subprocess
import sys
import textwrap
import xml.etree.ElementTree as ET
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent
BINARY = PROJECT_ROOT / "stroma_world"
BASE_CONFIG = PROJECT_ROOT / "config" / "PhysiCell_settings.xml"

SEED = 42
SLURM_PARTITION = "cpu384g"
SLURM_CPUS = 32
SLURM_MEM = "128G"
SLURM_TIME = "06:00:00"
SAVE_INTERVAL = 360  # minutes

# Phase boundaries (minutes)
T_PRE = 20160.0       # day 14
T_TREAT_END = 40320.0 # day 28
T_POST = 60480.0      # day 42


def _patch_config(src, dst, output_dir, seed, max_time, drug_start, drug_end, drug_conc):
    tree = ET.parse(src)
    root = tree.getroot()

    def _set(xpath, value):
        node = root.find(xpath)
        if node is not None:
            node.text = str(value)

    _set("./save/folder", str(output_dir))
    _set(".//overall/max_time", str(max_time))
    _set(".//options/random_seed", str(seed))
    _set(".//parallel/omp_num_threads", str(SLURM_CPUS))

    for n in root.findall(".//save//interval"):
        n.text = str(SAVE_INTERVAL)

    _set(".//user_parameters/drug_start_time", str(drug_start))
    _set(".//user_parameters/drug_end_time", str(drug_end))
    _set(".//user_parameters/drug_concentration", str(drug_conc))

    for var_name, val in [("oxygen", "38"), ("tgfb", "0"), ("shh", "0"),
                          ("drug", "0")]:
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


def main():
    if not BINARY.exists():
        print(f"ERROR: binary not found: {BINARY}")
        return 1

    work_dir = PROJECT_ROOT / "build" / "rc2_full_seed42"
    if work_dir.exists():
        shutil.rmtree(work_dir)
    work_dir.mkdir(parents=True)

    run_dir = work_dir / "replicate_01_seed42"
    run_dir.mkdir(parents=True)
    out_dir = run_dir / "output"
    out_dir.mkdir()

    config_path = run_dir / "config.xml"
    _patch_config(BASE_CONFIG, config_path, out_dir, SEED,
                  max_time=T_POST,
                  drug_start=T_PRE,
                  drug_end=T_TREAT_END,
                  drug_conc=1.0)

    local_cfg = run_dir / "config"
    local_cfg.mkdir(exist_ok=True)
    for f in ["tumor_calibration_knobs.json", "gene_params_default.json"]:
        src = PROJECT_ROOT / "config" / f
        if src.exists():
            shutil.copy(src, local_cfg / f)

    script = run_dir / "run.slurm.sh"
    script.write_text(textwrap.dedent(f"""\
        #!/bin/bash
        #SBATCH --job-name=rc2_full_s42
        #SBATCH --partition={SLURM_PARTITION}
        #SBATCH --nodes=1
        #SBATCH --ntasks-per-node=1
        #SBATCH --cpus-per-task={SLURM_CPUS}
        #SBATCH --mem={SLURM_MEM}
        #SBATCH --time={SLURM_TIME}
        #SBATCH --output={run_dir}/slurm_%j.out
        #SBATCH --error={run_dir}/slurm_%j.err

        set -euo pipefail
        export OMP_NUM_THREADS={SLURM_CPUS}
        export OMP_PROC_BIND=spread
        export OMP_PLACES=cores

        cd {PROJECT_ROOT}
        echo "=== RC2 Full — seed 42 ==="
        echo "Start: $(date)"
        {BINARY} {config_path}
        echo "End:   $(date)"
    """))
    script.chmod(0o755)

    result = subprocess.run(
        ["sbatch", "--parsable", str(script)],
        capture_output=True, text=True, check=False,
    )
    if result.returncode != 0:
        print(f"ERROR: sbatch failed: {result.stderr.strip()}")
        return 1

    job_id = result.stdout.strip().split(";")[0]
    print(f"  Submitted RC2 full (seed 42) -> SLURM job {job_id}")
    print(f"  Output: {out_dir}")
    print(f"  Monitor: bash watch_rc2.sh {job_id}")
    print(f"  Evaluate: python diagnose_rc2_full.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
