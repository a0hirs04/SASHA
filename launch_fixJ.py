#!/usr/bin/env python3
"""fixJ: 20-variant combined sweep (persistence timer + structural fix).

Phase 1 (default): submit 20 RC1-only jobs.
Phase 2 (--phase2 v01 v02 ...): submit RC2 for specified winning variants.

Binary: stroma_world_induct (persistence timer + ecm_reversion_weight).
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
OUT_BASE = Path("/work/a0hirs04/PROJECT-NORTHSTAR/build/fixJ")
BINARY = Path("/work/a0hirs04/PROJECT-NORTHSTAR/build/stroma_world_induct")

SEED = 42
PARTITION = os.environ.get("RC_SLURM_PARTITION", "cpu384g")
CPUS = 32
RC1_MAX = 30240   # 21 days
RC2_MAX = 60480   # 42 days
T_PRE = 20160     # drug start (day 14)
T_TREAT_END = 40320  # drug end (day 28)
SAVE_INTERVAL = 360

# ═══════════════════════════════════════════════════════════════════
# Fixed parameters (same for all 20 variants)
# ═══════════════════════════════════════════════════════════════════
FIXED = {
    "emt_induction_threshold": 0.30,
    "ecm_emt_cap": 0.30,
    "emt_off_threshold": 0.10,
    "drug_uptake_rate": 0.10,
    "hif1a_emt_boost": 0.02,
    "ecm_emt_require_caf_contact": 1.0,
    "emt_activation_delay": 360.0,
}

# ═══════════════════════════════════════════════════════════════════
# 20 Variants: persist_time × ecm_reversion_weight × drug_kill
# ═══════════════════════════════════════════════════════════════════
VARIANTS = OrderedDict()

# Grid v01-v16: persistence × ecm_reversion_weight at drug_kill=0.05
grid = [
    ("v01",  720, 0.15, 0.05),
    ("v02",  720, 0.25, 0.05),
    ("v03",  720, 0.35, 0.05),
    ("v04", 1080, 0.15, 0.05),
    ("v05", 1080, 0.25, 0.05),
    ("v06", 1080, 0.35, 0.05),
    ("v07", 1440, 0.15, 0.05),
    ("v08", 1440, 0.25, 0.05),
    ("v09", 1440, 0.35, 0.05),
    ("v10", 2160, 0.15, 0.05),
    ("v11", 2160, 0.25, 0.05),
    ("v12", 2160, 0.35, 0.05),
    ("v13", 2880, 0.20, 0.05),
    ("v14", 2880, 0.30, 0.05),
    ("v15", 4320, 0.20, 0.05),
    ("v16", 4320, 0.30, 0.05),
]

# Drug sensitivity v17-v20
drug_sens = [
    ("v17", 1440, 0.25, 0.04),
    ("v18", 1440, 0.25, 0.06),
    ("v19", 2160, 0.25, 0.04),
    ("v20", 2160, 0.25, 0.06),
]

for name, persist, rev, drug in grid + drug_sens:
    VARIANTS[name] = {
        "emt_persistence_time": float(persist),
        "ecm_reversion_weight": rev,
        "drug_kill_coefficient": drug,
    }


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

    # Apply fixed + variant overrides
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

    # Dirichlet BCs
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


def submit_job(vname, arm, max_time, is_rc2, slurm_time, overrides):
    run_dir = OUT_BASE / vname / arm
    out_dir = run_dir / "output"
    run_dir.mkdir(parents=True, exist_ok=True)
    out_dir.mkdir(parents=True, exist_ok=True)

    cfg = run_dir / "PhysiCell_settings.xml"
    patch_config(BASE_CONFIG, cfg, out_dir, max_time, is_rc2, overrides)

    script = run_dir / "run.slurm.sh"
    script.write_text(textwrap.dedent(f"""\
        #!/bin/bash
        #SBATCH --job-name=fJ_{vname}_{arm[:3]}
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

        echo "=== fixJ {vname} {arm} seed={SEED} ==="
        {BINARY} {cfg}
        echo "=== DONE ==="
    """))

    result = subprocess.run(["sbatch", str(script)], capture_output=True, text=True)
    job_id = result.stdout.strip().split()[-1] if result.stdout.strip() else "FAILED"
    return job_id


def phase1():
    """Submit 20 RC1-only jobs."""
    OUT_BASE.mkdir(parents=True, exist_ok=True)
    print("=" * 72)
    print("fixJ Phase 1: 20 RC1-only jobs (persistence timer + structural fix)")
    print(f"Binary: {BINARY}")
    print("=" * 72)
    print(f"{'V#':<5} {'persist':>7} {'ecm_rev':>8} {'drug_k':>7}  {'RC1 job':>8}")
    print("-" * 42)

    for vname, overrides in VARIANTS.items():
        job = submit_job(vname, "rc1/replicate_01_seed42", RC1_MAX, False,
                         "04:00:00", overrides)
        pt = overrides["emt_persistence_time"]
        rv = overrides["ecm_reversion_weight"]
        dk = overrides["drug_kill_coefficient"]
        print(f"{vname:<5} {pt:>7.0f} {rv:>8.2f} {dk:>7.2f}  {job:>8}")

    print(f"\nOutput: {OUT_BASE}")
    print("Monitor: bash watch_fixJ.sh")


def phase2(winners):
    """Submit RC2 for specified winning variants."""
    print("=" * 72)
    print(f"fixJ Phase 2: RC2 for winners {winners}")
    print("=" * 72)
    print(f"{'V#':<5} {'persist':>7} {'ecm_rev':>8} {'drug_k':>7}  {'RC2 job':>8}")
    print("-" * 42)

    for vname in winners:
        if vname not in VARIANTS:
            print(f"  WARNING: {vname} not in variant list, skipping")
            continue
        overrides = VARIANTS[vname]
        job = submit_job(vname, "rc2", RC2_MAX, True, "08:00:00", overrides)
        pt = overrides["emt_persistence_time"]
        rv = overrides["ecm_reversion_weight"]
        dk = overrides["drug_kill_coefficient"]
        print(f"{vname:<5} {pt:>7.0f} {rv:>8.2f} {dk:>7.2f}  {job:>8}")


def main():
    parser = argparse.ArgumentParser(description="fixJ launcher")
    parser.add_argument("--phase2", nargs="+", metavar="VARIANT",
                        help="Submit RC2 for these winning variants (e.g., --phase2 v08 v11)")
    args = parser.parse_args()

    if args.phase2:
        phase2(args.phase2)
    else:
        phase1()


if __name__ == "__main__":
    main()
