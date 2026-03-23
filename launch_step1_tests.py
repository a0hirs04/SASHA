#!/usr/bin/env python3
"""STEP 1: Two diagnostic runs to test hypotheses about C1.4 fix.

Run A: v09 params + EMT fix (ecm_emt_cap=0.30, emt_off_threshold=0.15)
       RC2 ONLY, seed 42.  Hypothesis: RC2-6 regrowth may fail.

Run B: v09 params + tgfb_secretion_rate=2.0 (keep current EMT params)
       RC1 + RC2, seed 42.  Hypothesis: stronger TGF-beta drives
       peripheral EMT without needing ecm_emt_cap changes.
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
TEST_DIR = Path("/work/a0hirs04/PROJECT-NORTHSTAR/build/step1_tests")

SEED = 42
PARTITION = os.environ.get("RC_SLURM_PARTITION", "cpu384g")
CPUS = 32
RC1_MAX = 30240     # 21 days
RC2_MAX = 60480     # 42 days
T_PRE = 20160       # day 14
T_TREAT_END = 40320 # day 28
SAVE_INTERVAL = 360


def patch_config(src, dst, output_dir, seed, max_time, params, is_rc2=False):
    tree = ET.parse(src)
    root = tree.getroot()

    def _set(xpath, value):
        node = root.find(xpath)
        if node is not None:
            node.text = str(value)

    _set("./save/folder", str(output_dir))
    _set(".//overall/max_time", str(max_time))
    _set(".//options/random_seed", str(seed))
    _set(".//parallel/omp_num_threads", str(CPUS))

    for n in root.findall(".//save//interval"):
        n.text = str(SAVE_INTERVAL)

    # Patch user parameters
    for key, value in params.items():
        _set(f".//user_parameters/{key}", str(value))

    # Drug timing
    if is_rc2:
        _set(".//user_parameters/drug_start_time", str(T_PRE))
        _set(".//user_parameters/drug_end_time", str(T_TREAT_END))
        _set(".//user_parameters/drug_concentration", "1.0")
    else:
        _set(".//user_parameters/drug_start_time", str(max_time + 10000))
        _set(".//user_parameters/drug_concentration", "0")

    # Enforce Dirichlet BCs
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


def write_slurm(rep_dir, config_path, job_name, slurm_time):
    script = rep_dir / "run.slurm.sh"
    script.write_text(textwrap.dedent(f"""\
        #!/bin/bash
        #SBATCH --job-name={job_name}
        #SBATCH --partition={PARTITION}
        #SBATCH --nodes=1
        #SBATCH --ntasks-per-node=1
        #SBATCH --cpus-per-task={CPUS}
        #SBATCH --mem=0
        #SBATCH --time={slurm_time}
        #SBATCH --output={rep_dir}/slurm_%j.out
        #SBATCH --error={rep_dir}/slurm_%j.err

        set -euo pipefail
        export OMP_NUM_THREADS={CPUS}
        export OMP_PROC_BIND=spread
        export OMP_PLACES=cores

        cd {PROJECT_ROOT}
        echo "=== {job_name} seed={SEED} ==="
        echo "Start: $(date)"
        {BINARY} {config_path}
        echo "End:   $(date)"
    """), encoding="utf-8")
    script.chmod(0o755)
    return script


def submit(script):
    r = subprocess.run(["sbatch", "--parsable", str(script)],
                       capture_output=True, text=True)
    return r.stdout.strip()


def main():
    if TEST_DIR.exists():
        shutil.rmtree(TEST_DIR)
    TEST_DIR.mkdir(parents=True)

    jobs = []

    # ═══════════════════════════════════════════════════════════════
    # Run A: v09 + EMT fix → RC2 ONLY
    # Hypothesis: Restoring Phase-1 EMT params may break RC2-6 regrowth
    # ═══════════════════════════════════════════════════════════════
    run_a_params = {
        'drug_uptake_rate': 0.10,
        'drug_kill_coefficient': 0.05,
        'hif1a_emt_boost': 0.02,
        'ecm_emt_cap': 0.30,         # restored from Phase-1
        'emt_off_threshold': 0.15,   # restored from Phase-1
    }
    a_dir = TEST_DIR / "runA_rc2"
    a_out = a_dir / "output"
    a_out.mkdir(parents=True)
    a_cfg = a_dir / "config.xml"
    patch_config(BASE_CONFIG, a_cfg, a_out, SEED, RC2_MAX, run_a_params, is_rc2=True)
    a_script = write_slurm(a_dir, a_cfg, "runA_rc2", "04:00:00")
    jid = submit(a_script)
    print(f"  Run A (RC2): job {jid}")
    print(f"    Params: v09 drug + EMT fix (cap=0.30, off=0.15)")
    print(f"    Test: Does RC2-6 regrowth still work?")
    jobs.append(("runA", "rc2", jid))

    # ═══════════════════════════════════════════════════════════════
    # Run B: v09 + tgfb=2.0 → RC1 + RC2
    # Hypothesis: Higher TGF-beta at periphery drives EMT without
    # needing ecm_emt_cap/emt_off_threshold changes
    # ═══════════════════════════════════════════════════════════════
    run_b_params = {
        'drug_uptake_rate': 0.10,
        'drug_kill_coefficient': 0.05,
        'hif1a_emt_boost': 0.02,
        'ecm_emt_cap': 0.20,          # keep current
        'emt_off_threshold': 0.25,    # keep current
        'tgfb_secretion_rate': 2.0,   # boosted from 0.8
    }

    # RC1
    b_rc1_work = TEST_DIR / "runB_rc1"
    b_rc1_rep = b_rc1_work / f"replicate_01_seed{SEED}"
    b_rc1_out = b_rc1_rep / "output"
    b_rc1_out.mkdir(parents=True)
    b_rc1_cfg = b_rc1_rep / "config.xml"
    patch_config(BASE_CONFIG, b_rc1_cfg, b_rc1_out, SEED, RC1_MAX,
                 run_b_params, is_rc2=False)
    b_rc1_script = write_slurm(b_rc1_rep, b_rc1_cfg, "runB_rc1", "02:00:00")
    jid1 = submit(b_rc1_script)
    print(f"\n  Run B (RC1): job {jid1}")
    print(f"    Params: v09 drug + tgfb=2.0 (current EMT params)")
    print(f"    Test: Does higher TGF-beta fix C1.4?")
    jobs.append(("runB", "rc1", jid1))

    # RC2
    b_rc2_dir = TEST_DIR / "runB_rc2"
    b_rc2_out = b_rc2_dir / "output"
    b_rc2_out.mkdir(parents=True)
    b_rc2_cfg = b_rc2_dir / "config.xml"
    patch_config(BASE_CONFIG, b_rc2_cfg, b_rc2_out, SEED, RC2_MAX,
                 run_b_params, is_rc2=True)
    b_rc2_script = write_slurm(b_rc2_dir, b_rc2_cfg, "runB_rc2", "04:00:00")
    jid2 = submit(b_rc2_script)
    print(f"  Run B (RC2): job {jid2}")
    print(f"    Params: same as RC1")
    print(f"    Test: Does tgfb=2.0 break RC2?")
    jobs.append(("runB", "rc2", jid2))

    # Write manifest
    manifest = TEST_DIR / "jobs.txt"
    with open(manifest, "w") as f:
        for name, rc, jid in jobs:
            f.write(f"{name}\t{rc}\t{jid}\n")

    print(f"\n  === {len(jobs)} jobs submitted to {TEST_DIR} ===")
    print(f"  Estimated completion:")
    print(f"    Run A (RC2):  ~1.5-2 hours")
    print(f"    Run B (RC1):  ~30-40 minutes")
    print(f"    Run B (RC2):  ~1.5-2 hours")
    print(f"\n  Evaluate with:")
    print(f"    # Run A (RC2):")
    print(f"    python3 {PROJECT_ROOT}/evaluate_rc2.py --out-dir {TEST_DIR}/runA_rc2/output")
    print(f"    # Run B (RC1):")
    print(f"    python3 {PROJECT_ROOT}/run_reality_check_1.py \\")
    print(f"        --work-dir {TEST_DIR}/runB_rc1 --seeds 42 --quorum 1 --evaluate-only")
    print(f"    # Run B (RC2):")
    print(f"    python3 {PROJECT_ROOT}/evaluate_rc2.py --out-dir {TEST_DIR}/runB_rc2/output")


if __name__ == "__main__":
    main()
