#!/usr/bin/env python3
"""fixL: Fix post-withdrawal drug toxicity trap.

Two C++ bug fixes:
  1. Intracellular drug decay now always active (was gated on local_drug<=0)
  2. Efflux pump works post-withdrawal (was gated on local_drug>0)

Plus XML parameter: drug_natural_decay_rate > 0 for metabolic clearance.

Uses fixJ v03 base (proven RC1 8/8): persist=720, ecm_rev=0.35, drug_kill=0.05.
4 variants sweeping drug_natural_decay_rate + 1 stronger drug variant.

Modes:
  --fast    : Run fast_rc2 screening only (drug_end + 5000 min, 16 CPUs)
  --full    : Run full RC1 + RC2 (default, original behavior)
  --promote : Run full RC2 only for variants that passed fast screening
"""
import argparse
import os
import subprocess
import textwrap
import xml.etree.ElementTree as ET
from collections import OrderedDict
from pathlib import Path

PROJECT_ROOT = Path("/home/a0hirs04/PROJECT-NORTHSTAR")
BASE_CONFIG = PROJECT_ROOT / "config" / "PhysiCell_settings.xml"
OUT_BASE = Path("/work/a0hirs04/PROJECT-NORTHSTAR/build/fixL")
BINARY = Path("/work/a0hirs04/PROJECT-NORTHSTAR/build/stroma_world_fixL")
BINARY_FAST = Path("/work/a0hirs04/PROJECT-NORTHSTAR/build/stroma_world_fixL2")

SEED = 42
PARTITION = os.environ.get("RC_SLURM_PARTITION", "cpu384g")
CPUS_FULL = 32
CPUS_FAST = 16
RC1_MAX = 30240
RC2_MAX = 60480
T_PRE = 20160
T_TREAT_END = 40320
SAVE_INTERVAL = 360
FAST_SAVE_INTERVAL = 120
FAST_MAX = T_TREAT_END + 5000  # drug_end + 5000 = 45320

# v03 base (proven RC1 8/8)
FIXED = {
    "emt_persistence_time": 720.0,
    "ecm_reversion_weight": 0.35,
    "emt_induction_threshold": 0.30,
    "ecm_emt_cap": 0.30,
    "emt_off_threshold": 0.10,
    "drug_uptake_rate": 0.10,
    "hif1a_emt_boost": 0.02,
    "ecm_emt_require_caf_contact": 1.0,
    "emt_activation_delay": 360.0,
}

VARIANTS = OrderedDict([
    ("a", {"drug_kill_coefficient": 0.05, "drug_natural_decay_rate": 0.02}),
    ("b", {"drug_kill_coefficient": 0.05, "drug_natural_decay_rate": 0.03}),
    ("c", {"drug_kill_coefficient": 0.05, "drug_natural_decay_rate": 0.05}),
    ("d", {"drug_kill_coefficient": 0.08, "drug_natural_decay_rate": 0.03}),
])


def patch_config(src, dst, output_dir, max_time, is_rc2, overrides,
                 save_interval=SAVE_INTERVAL, cpus=CPUS_FULL,
                 fast_screen=False):
    tree = ET.parse(src)
    root = tree.getroot()

    def _set(xpath, value):
        node = root.find(xpath)
        if node is not None:
            node.text = str(value)

    _set("./save/folder", str(output_dir))
    _set(".//overall/max_time", str(max_time))
    _set(".//options/random_seed", str(SEED))
    _set(".//parallel/omp_num_threads", str(cpus))

    for n in root.findall(".//save//interval"):
        n.text = str(save_interval)

    all_params = {**FIXED, **overrides}
    if fast_screen:
        all_params["fast_screen_mode"] = 1.0
    for param, value in all_params.items():
        # Ensure the parameter node exists; create if needed
        node = root.find(f".//user_parameters/{param}")
        if node is None:
            up = root.find(".//user_parameters")
            if up is not None:
                node = ET.SubElement(up, param, type="double", units="dimensionless")
        if node is not None:
            node.text = str(value)

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


def submit(vname, arm, max_time, is_rc2, slurm_time, overrides,
           save_interval=SAVE_INTERVAL, cpus=CPUS_FULL, fast_screen=False,
           job_prefix="fL", binary=None):
    if binary is None:
        binary = BINARY
    run_dir = OUT_BASE / vname / arm
    out_dir = run_dir / "output"
    run_dir.mkdir(parents=True, exist_ok=True)
    out_dir.mkdir(parents=True, exist_ok=True)

    cfg = run_dir / "PhysiCell_settings.xml"
    patch_config(BASE_CONFIG, cfg, out_dir, max_time, is_rc2, overrides,
                 save_interval=save_interval, cpus=cpus,
                 fast_screen=fast_screen)

    script = run_dir / "run.slurm.sh"
    script.write_text(textwrap.dedent(f"""\
        #!/bin/bash
        #SBATCH --job-name={job_prefix}_{vname}_{arm[:3]}
        #SBATCH --partition={PARTITION}
        #SBATCH --nodes=1
        #SBATCH --ntasks-per-node=1
        #SBATCH --cpus-per-task={cpus}
        #SBATCH --mem=0
        #SBATCH --time={slurm_time}
        #SBATCH --output={run_dir}/slurm_%j.out
        #SBATCH --error={run_dir}/slurm_%j.err

        set -euo pipefail
        export OMP_NUM_THREADS={cpus}
        export OMP_PROC_BIND=spread
        export OMP_PLACES=cores

        echo "=== fixL {vname} {arm} seed={SEED} ==="
        {binary} {cfg}
        echo "=== DONE ==="
    """))

    result = subprocess.run(["sbatch", str(script)], capture_output=True, text=True)
    job_id = result.stdout.strip().split()[-1] if result.stdout.strip() else "FAILED"
    return job_id


def mode_fast():
    """Submit fast_rc2 screening runs for all variants."""
    print("=" * 68)
    print("fixL FAST SCREEN: post-withdrawal dynamics screening")
    print(f"Binary: {BINARY}")
    print(f"max_time: {FAST_MAX} (drug_end + 5000)")
    print(f"CPUs: {CPUS_FAST}, save_interval: {FAST_SAVE_INTERVAL}")
    print("=" * 68)
    print(f"{'Var':<4} {'drug_kill':>10} {'decay_rate':>11}  {'Job':>8}")
    print("-" * 40)

    for vname, overrides in VARIANTS.items():
        job = submit(vname, "fast_rc2", FAST_MAX, True, "02:00:00", overrides,
                     save_interval=FAST_SAVE_INTERVAL, cpus=CPUS_FAST,
                     fast_screen=True, job_prefix="fLf", binary=BINARY_FAST)
        dk = overrides["drug_kill_coefficient"]
        dr = overrides["drug_natural_decay_rate"]
        print(f"{vname:<4} {dk:>10.2f} {dr:>11.2f}  {job:>8}")

    print(f"\nOutput: {OUT_BASE}/<var>/fast_rc2/output/")
    print("Evaluate: python3 evaluate_fast_rc2.py")


def mode_full():
    """Submit full RC1 + RC2 for all variants (original behavior)."""
    print("=" * 68)
    print("fixL FULL: RC1 + RC2 for all variants")
    print(f"Binary: {BINARY}")
    print("=" * 68)
    print(f"{'Var':<4} {'drug_kill':>10} {'decay_rate':>11}  {'RC1 job':>8} {'RC2 job':>8}")
    print("-" * 50)

    for vname, overrides in VARIANTS.items():
        rc1_job = submit(vname, "rc1/replicate_01_seed42", RC1_MAX, False,
                         "04:00:00", overrides)
        rc2_job = submit(vname, "rc2", RC2_MAX, True, "08:00:00", overrides)
        dk = overrides["drug_kill_coefficient"]
        dr = overrides["drug_natural_decay_rate"]
        print(f"{vname:<4} {dk:>10.2f} {dr:>11.2f}  {rc1_job:>8} {rc2_job:>8}")

    print(f"\nOutput: {OUT_BASE}")
    print("Monitor: bash watch_fixL.sh")


def mode_promote():
    """Submit full RC1 + RC2 only for variants that passed fast screening."""
    print("=" * 68)
    print("fixL PROMOTE: full runs for fast-screen PASS variants")
    print("=" * 68)

    passed = []
    for vname in VARIANTS:
        fail_marker = OUT_BASE / vname / "fast_rc2" / "output" / "FAST_SCREEN_FAIL.txt"
        csv_file = OUT_BASE / vname / "fast_rc2" / "output" / "fast_screen.csv"
        if fail_marker.exists():
            print(f"  {vname}: FAIL (early terminated)")
        elif csv_file.exists():
            passed.append(vname)
            print(f"  {vname}: PASS (promoting to full RC1+RC2)")
        else:
            print(f"  {vname}: no fast_rc2 results found")

    if not passed:
        print("\nNo variants passed fast screening. Nothing to promote.")
        return

    print(f"\n{'Var':<4} {'drug_kill':>10} {'decay_rate':>11}  {'RC1 job':>8} {'RC2 job':>8}")
    print("-" * 50)

    for vname in passed:
        overrides = VARIANTS[vname]
        rc1_job = submit(vname, "full_rc2/../rc1/replicate_01_seed42", RC1_MAX,
                         False, "04:00:00", overrides)
        rc2_job = submit(vname, "full_rc2", RC2_MAX, True, "08:00:00",
                         overrides)
        dk = overrides["drug_kill_coefficient"]
        dr = overrides["drug_natural_decay_rate"]
        print(f"{vname:<4} {dk:>10.2f} {dr:>11.2f}  {rc1_job:>8} {rc2_job:>8}")

    print(f"\nOutput: {OUT_BASE}")


def main():
    parser = argparse.ArgumentParser(description="fixL launcher")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--fast", action="store_true",
                       help="Run fast_rc2 screening only")
    group.add_argument("--full", action="store_true",
                       help="Run full RC1 + RC2 (default)")
    group.add_argument("--promote", action="store_true",
                       help="Run full RC1+RC2 for fast-screen PASS variants")
    args = parser.parse_args()

    OUT_BASE.mkdir(parents=True, exist_ok=True)

    if args.fast:
        mode_fast()
    elif args.promote:
        mode_promote()
    else:
        mode_full()


if __name__ == "__main__":
    main()
