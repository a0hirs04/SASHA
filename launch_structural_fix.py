#!/usr/bin/env python3
"""Launch RC1 + RC2 with the structural EMT reversion fix.

Binary already rebuilt with reversion_signal decoupling.
XML already has v09+restored EMT params:
  ecm_emt_cap=0.30, emt_off_threshold=0.15, hif1a_emt_boost=0.02,
  drug_uptake_rate=0.10, drug_kill_coefficient=0.05
"""
import os
import shutil
import subprocess
import textwrap
import xml.etree.ElementTree as ET
from pathlib import Path

PROJECT_ROOT = Path("/home/a0hirs04/PROJECT-NORTHSTAR")
BINARY = PROJECT_ROOT / "stroma_world"
BASE_CONFIG = PROJECT_ROOT / "config" / "PhysiCell_settings.xml"
OUT_BASE = Path("/work/a0hirs04/PROJECT-NORTHSTAR/build/structural_fix")

SEED = 42
PARTITION = os.environ.get("RC_SLURM_PARTITION", "cpu384g")
CPUS = 32
RC1_MAX = 30240     # 21 days
RC2_MAX = 60480     # 42 days
T_PRE = 20160       # day 14
T_TREAT_END = 40320 # day 28
SAVE_INTERVAL = 360


def patch_config(src, dst, output_dir, max_time, is_rc2=False):
    tree = ET.parse(src)
    root = tree.getroot()

    def _set(xpath, value):
        node = root.find(xpath)
        if node is not None:
            node.text = str(value)

    _set("./save/folder", str(output_dir))
    _set(".//overall/max_time", str(max_time))
    _set(".//options/random_seed", str(SEED))
    _set(".//parallel/omp_num_threads", str(CPUS))

    for n in root.findall(".//save//interval"):
        n.text = str(SAVE_INTERVAL)

    if is_rc2:
        _set(".//user_parameters/drug_start_time", str(T_PRE))
        _set(".//user_parameters/drug_end_time", str(T_TREAT_END))
        _set(".//user_parameters/drug_concentration", "1.0")
    else:
        _set(".//user_parameters/drug_start_time", str(max_time + 10000))
        _set(".//user_parameters/drug_concentration", "0")

    for var_name, value in [("oxygen", "38"), ("tgfb", "0"), ("shh", "0"),
                            ("drug", "0"), ("ecm_density", "0")]:
        var = root.find(f".//microenvironment_setup/variable[@name='{var_name}']")
        if var is None:
            continue
        dbc = var.find("./Dirichlet_boundary_condition")
        if dbc is not None:
            dbc.text = value
            dbc.set("enabled", "true")
        for bv in var.findall("./Dirichlet_options/boundary_value"):
            bv.text = value
            bv.set("enabled", "true")

    tree.write(dst, encoding="utf-8", xml_declaration=True)


def setup_and_launch(name, max_time, is_rc2, slurm_time):
    run_dir = OUT_BASE / name
    out_dir = run_dir / "output"
    run_dir.mkdir(parents=True, exist_ok=True)
    out_dir.mkdir(parents=True, exist_ok=True)

    cfg = run_dir / "PhysiCell_settings.xml"
    patch_config(BASE_CONFIG, cfg, out_dir, max_time, is_rc2)

    script = run_dir / "run.slurm.sh"
    script.write_text(textwrap.dedent(f"""\
        #!/bin/bash
        #SBATCH --job-name=sf_{name}
        #SBATCH --partition={PARTITION}
        #SBATCH --nodes=1
        #SBATCH --ntasks-per-node=1
        #SBATCH --cpus-per-task={CPUS}
        #SBATCH --mem=0
        #SBATCH --time={slurm_time}
        #SBATCH --output={run_dir}/slurm_%j.out
        #SBATCH --error={run_dir}/slurm_%j.err

        set -euo pipefail
        export OMP_NUM_THREADS={CPUS}
        export OMP_PROC_BIND=spread
        export OMP_PLACES=cores

        cd {PROJECT_ROOT}
        echo "=== structural_fix {name} seed={SEED} ==="
        ./stroma_world {cfg}
        echo "=== DONE ==="
    """))

    result = subprocess.run(["sbatch", str(script)], capture_output=True, text=True)
    job_id = result.stdout.strip().split()[-1] if result.stdout.strip() else "FAILED"
    print(f"  {name}: job {job_id}")
    return job_id


def main():
    OUT_BASE.mkdir(parents=True, exist_ok=True)
    print("Launching structural fix runs...")

    # RC1: 21 days, no drug
    rc1_dir = OUT_BASE / "rc1" / "replicate_01_seed42"
    rc1_dir.mkdir(parents=True, exist_ok=True)
    rc1_job = setup_and_launch("rc1/replicate_01_seed42", RC1_MAX, is_rc2=False, slurm_time="04:00:00")

    # RC2: 42 days, drug day 14-28
    rc2_job = setup_and_launch("rc2", RC2_MAX, is_rc2=True, slurm_time="08:00:00")

    print(f"\nOutputs: {OUT_BASE}")
    print(f"Monitor: squeue -u $USER")
    print(f"\nEvaluate when done:")
    print(f"  python3 run_reality_check_1.py --work-dir {OUT_BASE}/rc1 --seeds 42 --quorum 1 --evaluate-only")
    print(f"  python3 evaluate_rc2.py --out-dir {OUT_BASE}/rc2/output")


if __name__ == "__main__":
    main()
