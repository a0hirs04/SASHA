#!/usr/bin/env python3
"""Launch 8 ecm_reversion_weight variants — RC1 + RC2 each (16 jobs total).

All share v011 base: drug_uptake=0.10, drug_kill=0.05, hif1a_emt_boost=0.02,
ecm_emt_cap=0.30, ecm_emt_require_caf=1.0, emt_activation_delay=360, tgfb=0.8.
"""
import os
import subprocess
import textwrap
import xml.etree.ElementTree as ET
from pathlib import Path

PROJECT_ROOT = Path("/home/a0hirs04/PROJECT-NORTHSTAR")
BASE_CONFIG = PROJECT_ROOT / "config" / "PhysiCell_settings.xml"
WORK_BASE = Path("/work/a0hirs04/PROJECT-NORTHSTAR/build/structural_fix/rc1")

SEED = 42
PARTITION = os.environ.get("RC_SLURM_PARTITION", "cpu384g")
CPUS = 32
RC1_MAX = 30240
RC2_MAX = 60480
T_PRE = 20160
T_TREAT_END = 40320
SAVE_INTERVAL = 360

VARIANTS = {
    "v1": {"ecm_reversion_weight": 0.20, "emt_off_threshold": 0.10},
    "v2": {"ecm_reversion_weight": 0.25, "emt_off_threshold": 0.11},
    "v3": {"ecm_reversion_weight": 0.30, "emt_off_threshold": 0.12},
    "v4": {"ecm_reversion_weight": 0.35, "emt_off_threshold": 0.14},
    "v5": {"ecm_reversion_weight": 0.40, "emt_off_threshold": 0.15},
    "v6": {"ecm_reversion_weight": 0.50, "emt_off_threshold": 0.18},
    "v7": {"ecm_reversion_weight": 0.30, "emt_off_threshold": 0.13, "tgfb_secretion_rate": 1.2},
    "v8": {"ecm_reversion_weight": 0.25, "emt_off_threshold": 0.13, "hif1a_emt_boost": 0.04},
}

OUT_BASE = Path("/work/a0hirs04/PROJECT-NORTHSTAR/build/reversion_weight")


def patch_config(src, dst, output_dir, max_time, is_rc2, overrides):
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

    for param, value in overrides.items():
        _set(f".//user_parameters/{param}", str(value))

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


def submit(vname, arm, max_time, is_rc2, slurm_time, overrides):
    run_dir = OUT_BASE / vname / arm
    out_dir = run_dir / "output"
    run_dir.mkdir(parents=True, exist_ok=True)
    out_dir.mkdir(parents=True, exist_ok=True)

    cfg = run_dir / "PhysiCell_settings.xml"
    patch_config(BASE_CONFIG, cfg, out_dir, max_time, is_rc2, overrides)

    script = run_dir / "run.slurm.sh"
    script.write_text(textwrap.dedent(f"""\
        #!/bin/bash
        #SBATCH --job-name=rw_{vname}_{arm}
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
        echo "=== reversion_weight {vname} {arm} seed={SEED} ==="
        ./stroma_world {cfg}
        echo "=== DONE ==="
    """))

    result = subprocess.run(["sbatch", str(script)], capture_output=True, text=True)
    job_id = result.stdout.strip().split()[-1] if result.stdout.strip() else "FAILED"
    return job_id


def main():
    OUT_BASE.mkdir(parents=True, exist_ok=True)
    print("Launching ecm_reversion_weight sweep (8 variants x 2 arms = 16 jobs)\n")
    print(f"{'Variant':<6} {'ecm_rev_wt':>10} {'emt_off':>8} {'tgfb':>6} {'hif1a':>6}  {'RC1 job':>8} {'RC2 job':>8}")
    print("-" * 70)

    for vname, overrides in VARIANTS.items():
        rc1_job = submit(vname, "rc1/replicate_01_seed42", RC1_MAX, False, "04:00:00", overrides)
        rc2_job = submit(vname, "rc2", RC2_MAX, True, "08:00:00", overrides)
        print(f"{vname:<6} {overrides.get('ecm_reversion_weight', 0.3):>10.2f} "
              f"{overrides.get('emt_off_threshold', 0.15):>8.2f} "
              f"{overrides.get('tgfb_secretion_rate', 0.8):>6.1f} "
              f"{overrides.get('hif1a_emt_boost', 0.02):>6.2f}  "
              f"{rc1_job:>8} {rc2_job:>8}")

    print(f"\nOutputs: {OUT_BASE}")
    print("Monitor: squeue -u $USER | grep rw_")


if __name__ == "__main__":
    main()
