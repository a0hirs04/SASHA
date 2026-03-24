#!/usr/bin/env python3
"""fixK: Final sweep — v03 base (RC1 8/8 proven) + raised emt_off_threshold.

v03 base: persist=720, ecm_rev=0.35, drug_kill=0.05, threshold=0.30
The key tweak: raise emt_off_threshold so cells can revert after treatment.

Math: reversion_signal ≈ 0.135 with ECM maxed.
  - emt_off=0.10: 0.135 > 0.10 → NO reversion → RC2-6 FAIL
  - emt_off=0.14: 0.135 < 0.14 → reversion → RC2-6 should PASS

4 variants × (RC1 + RC2) = 8 SLURM jobs (fits QOS limit).
"""
import os
import subprocess
import textwrap
import xml.etree.ElementTree as ET
from collections import OrderedDict
from pathlib import Path

PROJECT_ROOT = Path("/home/a0hirs04/PROJECT-NORTHSTAR")
BASE_CONFIG = PROJECT_ROOT / "config" / "PhysiCell_settings.xml"
OUT_BASE = Path("/work/a0hirs04/PROJECT-NORTHSTAR/build/fixK")
BINARY = Path("/work/a0hirs04/PROJECT-NORTHSTAR/build/stroma_world_induct")

SEED = 42
PARTITION = os.environ.get("RC_SLURM_PARTITION", "cpu384g")
CPUS = 32
RC1_MAX = 30240
RC2_MAX = 60480
T_PRE = 20160
T_TREAT_END = 40320
SAVE_INTERVAL = 360

# v03 base (proven RC1 8/8) + fixed params
FIXED = {
    "emt_persistence_time": 720.0,
    "ecm_reversion_weight": 0.35,
    "emt_induction_threshold": 0.30,
    "ecm_emt_cap": 0.30,
    "drug_uptake_rate": 0.10,
    "hif1a_emt_boost": 0.02,
    "ecm_emt_require_caf_contact": 1.0,
    "emt_activation_delay": 360.0,
}

VARIANTS = OrderedDict([
    ("a", {"emt_off_threshold": 0.12, "drug_kill_coefficient": 0.05}),
    ("b", {"emt_off_threshold": 0.14, "drug_kill_coefficient": 0.05}),
    ("c", {"emt_off_threshold": 0.15, "drug_kill_coefficient": 0.05}),
    ("d", {"emt_off_threshold": 0.14, "drug_kill_coefficient": 0.04}),
])


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

    all_params = {**FIXED, **overrides}
    for param, value in all_params.items():
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
        #SBATCH --job-name=fK_{vname}_{arm[:3]}
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

        echo "=== fixK {vname} {arm} seed={SEED} ==="
        {BINARY} {cfg}
        echo "=== DONE ==="
    """))

    result = subprocess.run(["sbatch", str(script)], capture_output=True, text=True)
    job_id = result.stdout.strip().split()[-1] if result.stdout.strip() else "FAILED"
    return job_id


def main():
    OUT_BASE.mkdir(parents=True, exist_ok=True)
    print("=" * 68)
    print("fixK: Final sweep — v03 base + raised emt_off_threshold")
    print(f"Binary: {BINARY}")
    print("=" * 68)
    print(f"{'Var':<4} {'off_thresh':>10} {'drug_kill':>10}  {'RC1 job':>8} {'RC2 job':>8}")
    print("-" * 50)

    for vname, overrides in VARIANTS.items():
        rc1_job = submit(vname, "rc1/replicate_01_seed42", RC1_MAX, False,
                         "04:00:00", overrides)
        rc2_job = submit(vname, "rc2", RC2_MAX, True, "08:00:00", overrides)
        ot = overrides["emt_off_threshold"]
        dk = overrides["drug_kill_coefficient"]
        print(f"{vname:<4} {ot:>10.2f} {dk:>10.2f}  {rc1_job:>8} {rc2_job:>8}")

    print(f"\nOutput: {OUT_BASE}")


if __name__ == "__main__":
    main()
