#!/usr/bin/env python3
"""Launch fixI sweep: 2 approaches x multiple variants = 8 SLURM jobs.

Set 1 — Persistence timer (4 jobs):
  Reverted structural fix (reversion = full induction signal).
  Added emt_persistence_time so EMT can't revert for N minutes.
  v011 base params. Two persistence times, each RC1+RC2.
    pt1: emt_persistence_time=7200 (5 days)
    pt2: emt_persistence_time=4320 (3 days)

Set 2 — Induction threshold fix (4 jobs):
  Kept structural fix (ecm_reversion_weight=0.3).
  Lower emt_induction_threshold so peripheral TGFb can trigger EMT.
  RC1 only (test C1.4).
    th18: threshold=0.18
    th19: threshold=0.19
    th20: threshold=0.20
    th22: threshold=0.22
"""
import os
import subprocess
import textwrap
import xml.etree.ElementTree as ET
from pathlib import Path

PROJECT_ROOT = Path("/home/a0hirs04/PROJECT-NORTHSTAR")
BASE_CONFIG = PROJECT_ROOT / "config" / "PhysiCell_settings.xml"
OUT_BASE = Path("/work/a0hirs04/PROJECT-NORTHSTAR/build/fixI")

# Binaries
BIN_PERSIST = Path("/work/a0hirs04/PROJECT-NORTHSTAR/build/stroma_world_persist")
BIN_INDUCT = Path("/work/a0hirs04/PROJECT-NORTHSTAR/build/stroma_world_induct")

SEED = 42
PARTITION = os.environ.get("RC_SLURM_PARTITION", "cpu384g")
CPUS = 32
RC1_MAX = 30240
RC2_MAX = 60480
T_PRE = 20160
T_TREAT_END = 40320
SAVE_INTERVAL = 360

# ── Set 1: Persistence timer variants ──────────────────────────────
# v011 base: ecm_emt_cap=0.30, emt_off_threshold=0.10, drug_kill=0.05,
#            drug_uptake=0.10, hif1a_emt_boost=0.02, ecm_emt_require_caf=1.0,
#            ecm_reversion_weight irrelevant (code uses full signal)
PERSIST_VARIANTS = {
    "pt1": {
        "emt_persistence_time": 7200.0,
        "emt_off_threshold": 0.10,
        "drug_kill_coefficient": 0.05,
        "ecm_emt_cap": 0.30,
        "ecm_emt_require_caf_contact": 1.0,
    },
    "pt2": {
        "emt_persistence_time": 4320.0,
        "emt_off_threshold": 0.10,
        "drug_kill_coefficient": 0.05,
        "ecm_emt_cap": 0.30,
        "ecm_emt_require_caf_contact": 1.0,
    },
}

# ── Set 2: Induction threshold variants ────────────────────────────
# Structural fix active (ecm_reversion_weight=0.3). Lower threshold
# so peripheral cells (TGFb-only ~0.07-0.15) can cross.
INDUCT_VARIANTS = {
    "th18": {
        "emt_induction_threshold": 0.18,
        "ecm_reversion_weight": 0.3,
        "emt_off_threshold": 0.10,
        "ecm_emt_cap": 0.30,
        "ecm_emt_require_caf_contact": 1.0,
        "emt_persistence_time": 0.0,  # disabled for this approach
    },
    "th19": {
        "emt_induction_threshold": 0.19,
        "ecm_reversion_weight": 0.3,
        "emt_off_threshold": 0.10,
        "ecm_emt_cap": 0.30,
        "ecm_emt_require_caf_contact": 1.0,
        "emt_persistence_time": 0.0,
    },
    "th20": {
        "emt_induction_threshold": 0.20,
        "ecm_reversion_weight": 0.3,
        "emt_off_threshold": 0.10,
        "ecm_emt_cap": 0.30,
        "ecm_emt_require_caf_contact": 1.0,
        "emt_persistence_time": 0.0,
    },
    "th22": {
        "emt_induction_threshold": 0.22,
        "ecm_reversion_weight": 0.3,
        "emt_off_threshold": 0.10,
        "ecm_emt_cap": 0.30,
        "ecm_emt_require_caf_contact": 1.0,
        "emt_persistence_time": 0.0,
    },
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


def submit(vname, arm, max_time, is_rc2, slurm_time, overrides, binary):
    run_dir = OUT_BASE / vname / arm
    out_dir = run_dir / "output"
    run_dir.mkdir(parents=True, exist_ok=True)
    out_dir.mkdir(parents=True, exist_ok=True)

    cfg = run_dir / "PhysiCell_settings.xml"
    patch_config(BASE_CONFIG, cfg, out_dir, max_time, is_rc2, overrides)

    script = run_dir / "run.slurm.sh"
    script.write_text(textwrap.dedent(f"""\
        #!/bin/bash
        #SBATCH --job-name=fI_{vname}_{arm[:3]}
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

        echo "=== fixI {vname} {arm} seed={SEED} ==="
        {binary} {cfg}
        echo "=== DONE ==="
    """))

    result = subprocess.run(["sbatch", str(script)], capture_output=True, text=True)
    job_id = result.stdout.strip().split()[-1] if result.stdout.strip() else "FAILED"
    return job_id


def main():
    OUT_BASE.mkdir(parents=True, exist_ok=True)

    # ── Set 1: Persistence timer ──
    print("=" * 65)
    print("SET 1: Persistence Timer (reverted structural fix)")
    print(f"Binary: {BIN_PERSIST}")
    print("=" * 65)
    print(f"{'Var':<5} {'persist_min':>11} {'drug_kill':>10}  {'RC1 job':>8} {'RC2 job':>8}")
    print("-" * 55)
    for vname, overrides in PERSIST_VARIANTS.items():
        rc1_job = submit(vname, "rc1/replicate_01_seed42", RC1_MAX, False, "04:00:00",
                         overrides, BIN_PERSIST)
        rc2_job = submit(vname, "rc2", RC2_MAX, True, "08:00:00",
                         overrides, BIN_PERSIST)
        pt = overrides.get("emt_persistence_time", 7200)
        dk = overrides.get("drug_kill_coefficient", 0.05)
        print(f"{vname:<5} {pt:>11.0f} {dk:>10.2f}  {rc1_job:>8} {rc2_job:>8}")

    # ── Set 2: Induction threshold ──
    print()
    print("=" * 65)
    print("SET 2: Induction Threshold Fix (structural fix kept)")
    print(f"Binary: {BIN_INDUCT}")
    print("=" * 65)
    print(f"{'Var':<5} {'threshold':>10}  {'RC1 job':>8}")
    print("-" * 30)
    for vname, overrides in INDUCT_VARIANTS.items():
        rc1_job = submit(vname, "rc1/replicate_01_seed42", RC1_MAX, False, "04:00:00",
                         overrides, BIN_INDUCT)
        th = overrides.get("emt_induction_threshold", 0.30)
        print(f"{vname:<5} {th:>10.2f}  {rc1_job:>8}")

    print(f"\nOutputs: {OUT_BASE}")
    print("Monitor: squeue -u $USER | grep fI_")


if __name__ == "__main__":
    main()
