#!/usr/bin/env python3
"""Launch 10-variant parameter sweep for RC1 + RC2.

Each variant patches 3 parameters: drug_uptake_rate, drug_kill_coefficient, hif1a_emt_boost.
Submits 20 SLURM jobs (10 RC1 + 10 RC2), all seed=42.
"""
import os
import shutil
import subprocess
import sys
import textwrap
import xml.etree.ElementTree as ET
from pathlib import Path

PROJECT_ROOT = Path("/home/a0hirs04/PROJECT-NORTHSTAR")
BINARY = PROJECT_ROOT / "stroma_world"
BASE_CONFIG = PROJECT_ROOT / "config" / "PhysiCell_settings.xml"
SWEEP_DIR = Path("/work/a0hirs04/PROJECT-NORTHSTAR/build/sweep")

SEED = 42
SLURM_PARTITION = os.environ.get("RC_SLURM_PARTITION", "cpu384g")
SLURM_CPUS = 32
SLURM_MEM = 0
SLURM_TIME_RC1 = "02:00:00"
SLURM_TIME_RC2 = "04:00:00"

RC1_MAX_TIME = 30240    # 21 days
RC2_MAX_TIME = 60480    # 42 days
T_PRE = 20160           # day 14
T_TREAT_END = 40320     # day 28
SAVE_INTERVAL = 360

VARIANTS = [
    # (name, drug_uptake_rate, drug_kill_coefficient, hif1a_emt_boost)
    ("v01", 0.05, 0.02, 0.05),
    ("v02", 0.10, 0.02, 0.05),
    ("v03", 0.20, 0.02, 0.05),
    ("v04", 0.05, 0.05, 0.05),
    ("v05", 0.10, 0.05, 0.05),
    ("v06", 0.05, 0.10, 0.05),
    ("v07", 0.10, 0.10, 0.05),
    ("v08", 0.10, 0.02, 0.02),
    ("v09", 0.10, 0.05, 0.02),
    ("v10", 0.10, 0.02, 0.00),
]


def patch_config(src: Path, dst: Path, output_dir: Path,
                 seed: int, max_time: float,
                 drug_uptake_rate: float,
                 drug_kill_coefficient: float,
                 hif1a_emt_boost: float,
                 is_rc2: bool = False) -> None:
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

    # Variant parameters
    _set(".//user_parameters/drug_uptake_rate", str(drug_uptake_rate))
    _set(".//user_parameters/drug_kill_coefficient", str(drug_kill_coefficient))
    _set(".//user_parameters/hif1a_emt_boost", str(hif1a_emt_boost))

    # Drug timing
    if is_rc2:
        _set(".//user_parameters/drug_start_time", str(T_PRE))
        _set(".//user_parameters/drug_end_time", str(T_TREAT_END))
        _set(".//user_parameters/drug_concentration", "1.0")
    else:
        # RC1: no drug
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


def write_slurm(rep_dir: Path, config_path: Path, job_name: str,
                slurm_time: str) -> Path:
    script = rep_dir / "run.slurm.sh"
    script.write_text(textwrap.dedent(f"""\
        #!/bin/bash
        #SBATCH --job-name={job_name}
        #SBATCH --partition={SLURM_PARTITION}
        #SBATCH --nodes=1
        #SBATCH --ntasks-per-node=1
        #SBATCH --cpus-per-task={SLURM_CPUS}
        #SBATCH --mem={SLURM_MEM}
        #SBATCH --time={slurm_time}
        #SBATCH --output={rep_dir}/slurm_%j.out
        #SBATCH --error={rep_dir}/slurm_%j.err

        set -euo pipefail
        export OMP_NUM_THREADS={SLURM_CPUS}
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


def submit(script: Path) -> str:
    r = subprocess.run(["sbatch", "--parsable", str(script)],
                       capture_output=True, text=True)
    return r.stdout.strip()


def main():
    # Clean and recreate sweep dir
    if SWEEP_DIR.exists():
        shutil.rmtree(SWEEP_DIR)
    SWEEP_DIR.mkdir(parents=True)

    jobs = []

    for name, uptake, kill, hif1a in VARIANTS:
        vdir = SWEEP_DIR / name

        # RC1 — evaluator expects work-dir/replicate_01_seed42/output/
        rc1_work = vdir / "rc1"
        rc1_rep = rc1_work / f"replicate_01_seed{SEED}"
        rc1_out = rc1_rep / "output"
        rc1_out.mkdir(parents=True, exist_ok=True)
        rc1_cfg = rc1_rep / "config.xml"
        patch_config(BASE_CONFIG, rc1_cfg, rc1_out, SEED, RC1_MAX_TIME,
                     uptake, kill, hif1a, is_rc2=False)
        rc1_script = write_slurm(rc1_rep, rc1_cfg,
                                 f"{name}_rc1", SLURM_TIME_RC1)
        jid1 = submit(rc1_script)
        print(f"  {name} RC1 -> job {jid1}")
        jobs.append((name, "rc1", jid1))

        # RC2 — evaluator takes --out-dir pointing to output/
        rc2_rep = vdir / "rc2"
        rc2_out = rc2_rep / "output"
        rc2_out.mkdir(parents=True, exist_ok=True)
        rc2_cfg = rc2_rep / "config.xml"
        patch_config(BASE_CONFIG, rc2_cfg, rc2_out, SEED, RC2_MAX_TIME,
                     uptake, kill, hif1a, is_rc2=True)
        rc2_script = write_slurm(rc2_rep, rc2_cfg,
                                 f"{name}_rc2", SLURM_TIME_RC2)
        jid2 = submit(rc2_script)
        print(f"  {name} RC2 -> job {jid2}")
        jobs.append((name, "rc2", jid2))

    # Write job manifest for watch script
    manifest = SWEEP_DIR / "jobs.txt"
    with open(manifest, "w") as f:
        for name, rc, jid in jobs:
            f.write(f"{name}\t{rc}\t{jid}\n")

    print(f"\n  All {len(jobs)} jobs submitted. Monitor with:")
    print(f"    watch -n 30 bash {PROJECT_ROOT}/watch_sweep.sh")


if __name__ == "__main__":
    main()
