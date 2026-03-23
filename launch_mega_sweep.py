#!/usr/bin/env python3
"""Launch 100-variant parameter sweep for RC1-RC3.

Each variant tests a unique combination of parameters targeting:
  - RC1 C1.4 fix (EMT spatial pattern: peripheral > core)
  - RC2 C2.1-C2.6 (drug response calibration)
  - RC3 C3.1-C3.6 (vismodegib paradox: SHH inhibition arms)

Per variant: 5 SLURM jobs (RC1, RC2, RC3-ArmA, RC3-ArmB, RC3-ArmC)
Total: 500 SLURM jobs.

RC4/RC5 require ECM-degradation intervention (hyaluronidase) which needs
code-level support not yet parameterizable via XML. These checks are
deferred until a winning RC1-RC3 config is found and ECM degradation
is implemented.

Groups:
  A (27): EMT-cap restoration + drug grid
  B (18): TGF-beta boost approach
  C (20): Drug-strength grid with Phase-1 EMT fix
  D (18): Combined TGF-beta + EMT fine-tuning
  E (10): HIF1-alpha exploration
  F  (7): Edge cases + golden config
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
SWEEP_DIR = Path("/work/a0hirs04/PROJECT-NORTHSTAR/build/mega_sweep")

SEED = 42
PARTITION = os.environ.get("RC_SLURM_PARTITION", "cpu384g")
CPUS = 32

RC1_MAX_TIME = 30240     # 21 days
RC2_MAX_TIME = 60480     # 42 days
RC3_MAX_TIME = 80640     # 56 days (8 weeks)
T_PRE = 20160            # day 14
T_TREAT_END = 40320      # day 28 (also: RC3 intervention start)
RC3_DRUG_END = 60480     # day 42 (RC3 drug: weeks 4-6 = days 28-42)
SAVE_INTERVAL = 360

SLURM_TIME_RC1 = "02:00:00"
SLURM_TIME_RC2 = "04:00:00"
SLURM_TIME_RC3 = "06:00:00"


# ═══════════════════════════════════════════════════════════════════════
#  100 VARIANTS
#  Each: (name, drug_uptake, drug_kill, hif1a_emt_boost,
#         ecm_emt_cap, emt_off_threshold, tgfb_secretion_rate)
# ═══════════════════════════════════════════════════════════════════════
def generate_variants():
    """Generate 100 unique parameter configurations."""
    seen = set()
    variants = []

    def _add(uptake, kill, hif, cap, off, tgfb):
        key = (uptake, kill, hif, cap, off, tgfb)
        if key not in seen:
            seen.add(key)
            name = f"v{len(variants)+1:03d}"
            variants.append((name, uptake, kill, hif, cap, off, tgfb))

    # ── GROUP A (27): EMT-cap restoration + drug calibration ──
    # Strategy: Restore ecm_emt_cap / emt_off_threshold from Phase-1-RC1
    # Drug: vary uptake around v09 sweet spot
    for cap in [0.25, 0.30, 0.35]:
        for off in [0.10, 0.15, 0.20]:
            for uptake in [0.08, 0.10, 0.15]:
                _add(uptake, 0.05, 0.02, cap, off, 0.8)

    # ── GROUP B (18): TGF-beta boost ──
    # Strategy: Increase tgfb to drive peripheral EMT via higher
    # concentration at tumor-stroma interface
    for tgfb in [1.2, 1.6, 2.0]:
        for cap in [0.20, 0.25, 0.30]:
            for off in [0.15, 0.25]:
                _add(0.10, 0.05, 0.02, cap, off, tgfb)

    # ── GROUP C (20): Drug-strength grid with Phase-1 EMT fix ──
    # Strategy: Fix EMT (Phase-1), then sweep drug space broadly
    for uptake in [0.05, 0.08, 0.10, 0.15, 0.20]:
        for kill in [0.03, 0.04, 0.05, 0.07]:
            _add(uptake, kill, 0.02, 0.30, 0.15, 0.8)

    # ── GROUP D (18): Combined TGF-beta + EMT fine-tuning ──
    # Strategy: Both approaches together for maximum coverage
    for tgfb in [1.2, 1.5, 2.0]:
        for cap in [0.28, 0.30, 0.33]:
            for off in [0.12, 0.15]:
                _add(0.10, 0.05, 0.02, cap, off, tgfb)

    # ── GROUP E (10): HIF1-alpha exploration ──
    # Strategy: Explore hif1a_emt_boost range with fixed EMT params
    for hif in [0.00, 0.03, 0.06, 0.08, 0.10]:
        for cap in [0.25, 0.30]:
            _add(0.10, 0.05, hif, cap, 0.15, 0.8)

    # ── GROUP F: Edge cases + golden config ──
    _add(0.25, 0.03, 0.02, 0.30, 0.15, 0.8)   # very high uptake, low kill
    _add(0.10, 0.05, 0.02, 0.40, 0.08, 0.8)   # very high cap, very low off
    _add(0.10, 0.05, 0.00, 0.30, 0.15, 2.5)   # zero hif, high tgfb
    _add(0.12, 0.04, 0.02, 0.30, 0.15, 1.0)   # moderate everything
    _add(0.08, 0.06, 0.02, 0.30, 0.15, 0.8)   # low uptake, high kill
    _add(0.15, 0.05, 0.02, 0.30, 0.15, 1.2)   # moderate uplift
    _add(0.10, 0.05, 0.02, 0.30, 0.10, 1.5)   # low off + tgfb boost
    _add(0.10, 0.03, 0.02, 0.30, 0.15, 0.8)   # gentle kill
    _add(0.10, 0.08, 0.02, 0.30, 0.15, 0.8)   # strong kill
    _add(0.20, 0.05, 0.02, 0.30, 0.15, 0.8)   # high uptake
    _add(0.10, 0.05, 0.02, 0.35, 0.12, 0.8)   # high cap, low off
    _add(0.10, 0.05, 0.04, 0.30, 0.15, 1.2)   # moderate hif, tgfb boost
    _add(0.12, 0.05, 0.02, 0.28, 0.15, 1.6)   # near-sweet-spot combo
    _add(0.10, 0.06, 0.02, 0.30, 0.15, 0.8)   # slightly higher kill
    _add(0.10, 0.05, 0.02, 0.30, 0.18, 0.8)   # slightly higher off

    # Pad to 100 if needed with additional explorations
    extras = [
        (0.10, 0.05, 0.01, 0.30, 0.15, 0.8),
        (0.10, 0.05, 0.05, 0.30, 0.15, 0.8),
        (0.10, 0.05, 0.02, 0.22, 0.15, 1.8),
        (0.10, 0.05, 0.02, 0.27, 0.13, 0.8),
        (0.13, 0.05, 0.02, 0.30, 0.15, 0.8),
        (0.10, 0.045, 0.02, 0.30, 0.15, 0.8),
        (0.10, 0.055, 0.02, 0.30, 0.15, 0.8),
        (0.10, 0.05, 0.02, 0.32, 0.14, 0.8),
        (0.10, 0.05, 0.02, 0.30, 0.15, 1.4),
        (0.10, 0.05, 0.02, 0.30, 0.15, 1.8),
        (0.10, 0.05, 0.02, 0.30, 0.15, 0.6),
        (0.11, 0.05, 0.02, 0.30, 0.15, 0.8),
        (0.09, 0.05, 0.02, 0.30, 0.15, 0.8),
        (0.10, 0.05, 0.02, 0.30, 0.15, 0.9),
        (0.10, 0.05, 0.02, 0.30, 0.15, 1.1),
    ]
    for e in extras:
        if len(variants) >= 100:
            break
        _add(*e)

    return variants[:100]


# ═══════════════════════════════════════════════════════════════════════
#  CONFIG PATCHING
# ═══════════════════════════════════════════════════════════════════════
def patch_config(src, dst, output_dir, seed, max_time,
                 drug_uptake_rate, drug_kill_coefficient,
                 hif1a_emt_boost, ecm_emt_cap, emt_off_threshold,
                 tgfb_secretion_rate,
                 is_rc2=False,
                 rc3_mode=None):
    """Patch XML config for a specific run.

    rc3_mode: None (RC1/RC2), 'control', 'shh_only', 'shh_drug'
    """
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

    # ── Variant parameters ──
    _set(".//user_parameters/drug_uptake_rate", str(drug_uptake_rate))
    _set(".//user_parameters/drug_kill_coefficient", str(drug_kill_coefficient))
    _set(".//user_parameters/hif1a_emt_boost", str(hif1a_emt_boost))
    _set(".//user_parameters/ecm_emt_cap", str(ecm_emt_cap))
    _set(".//user_parameters/emt_off_threshold", str(emt_off_threshold))
    _set(".//user_parameters/tgfb_secretion_rate", str(tgfb_secretion_rate))

    # ── Drug timing ──
    if is_rc2:
        # RC2: drug days 14-28
        _set(".//user_parameters/drug_start_time", str(T_PRE))
        _set(".//user_parameters/drug_end_time", str(T_TREAT_END))
        _set(".//user_parameters/drug_concentration", "1.0")
    elif rc3_mode == 'shh_drug':
        # RC3 Arm C: SHH off + drug, drug days 28-42
        _set(".//user_parameters/drug_start_time", str(T_TREAT_END))
        _set(".//user_parameters/drug_end_time", str(RC3_DRUG_END))
        _set(".//user_parameters/drug_concentration", "1.0")
    else:
        # RC1, RC3-ArmA, RC3-ArmB: no drug
        _set(".//user_parameters/drug_start_time", str(max_time + 10000))
        _set(".//user_parameters/drug_concentration", "0")

    # ── SHH inhibition (RC3) ──
    if rc3_mode in ('shh_only', 'shh_drug'):
        # SHH inhibition starts at day 28 (week 4)
        _set(".//user_parameters/shh_inhibition_start_time", str(T_TREAT_END))
        _set(".//user_parameters/shh_inhibition_strength", "1.0")
    else:
        # No SHH inhibition
        _set(".//user_parameters/shh_inhibition_start_time", "1e18")
        _set(".//user_parameters/shh_inhibition_strength", "0.0")

    # ── Dirichlet BCs ──
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
    # Parse args
    include_rc3 = "--include-rc3" in sys.argv
    dry_run = "--dry-run" in sys.argv

    variants = generate_variants()
    print(f"  Generated {len(variants)} unique parameter configurations")

    if SWEEP_DIR.exists():
        shutil.rmtree(SWEEP_DIR)
    SWEEP_DIR.mkdir(parents=True)

    jobs = []
    total_rc1 = 0
    total_rc2 = 0
    total_rc3 = 0

    for name, uptake, kill, hif, cap, off, tgfb in variants:
        vdir = SWEEP_DIR / name
        params = dict(drug_uptake_rate=uptake, drug_kill_coefficient=kill,
                      hif1a_emt_boost=hif, ecm_emt_cap=cap,
                      emt_off_threshold=off, tgfb_secretion_rate=tgfb)

        # ── RC1: Natural history (21 days, no drug) ──
        rc1_work = vdir / "rc1"
        rc1_rep = rc1_work / f"replicate_01_seed{SEED}"
        rc1_out = rc1_rep / "output"
        rc1_out.mkdir(parents=True, exist_ok=True)
        rc1_cfg = rc1_rep / "config.xml"
        patch_config(BASE_CONFIG, rc1_cfg, rc1_out, SEED, RC1_MAX_TIME,
                     uptake, kill, hif, cap, off, tgfb, is_rc2=False)
        rc1_script = write_slurm(rc1_rep, rc1_cfg, f"{name}_rc1", SLURM_TIME_RC1)
        if not dry_run:
            jid = submit(rc1_script)
            jobs.append((name, "rc1", jid))
        total_rc1 += 1

        # ── RC2: Drug response (42 days) ──
        rc2_dir = vdir / "rc2"
        rc2_out = rc2_dir / "output"
        rc2_out.mkdir(parents=True, exist_ok=True)
        rc2_cfg = rc2_dir / "config.xml"
        patch_config(BASE_CONFIG, rc2_cfg, rc2_out, SEED, RC2_MAX_TIME,
                     uptake, kill, hif, cap, off, tgfb, is_rc2=True)
        rc2_script = write_slurm(rc2_dir, rc2_cfg, f"{name}_rc2", SLURM_TIME_RC2)
        if not dry_run:
            jid = submit(rc2_script)
            jobs.append((name, "rc2", jid))
        total_rc2 += 1

        if include_rc3:
            # ── RC3-ArmA: Control (56 days, no intervention) ──
            rc3a_dir = vdir / "rc3_armA"
            rc3a_out = rc3a_dir / "output"
            rc3a_out.mkdir(parents=True, exist_ok=True)
            rc3a_cfg = rc3a_dir / "config.xml"
            patch_config(BASE_CONFIG, rc3a_cfg, rc3a_out, SEED, RC3_MAX_TIME,
                         uptake, kill, hif, cap, off, tgfb,
                         is_rc2=False, rc3_mode='control')
            rc3a_script = write_slurm(rc3a_dir, rc3a_cfg,
                                      f"{name}_r3a", SLURM_TIME_RC3)
            if not dry_run:
                jid = submit(rc3a_script)
                jobs.append((name, "rc3a", jid))
            total_rc3 += 1

            # ── RC3-ArmB: SHH inhibition only (56 days) ──
            rc3b_dir = vdir / "rc3_armB"
            rc3b_out = rc3b_dir / "output"
            rc3b_out.mkdir(parents=True, exist_ok=True)
            rc3b_cfg = rc3b_dir / "config.xml"
            patch_config(BASE_CONFIG, rc3b_cfg, rc3b_out, SEED, RC3_MAX_TIME,
                         uptake, kill, hif, cap, off, tgfb,
                         is_rc2=False, rc3_mode='shh_only')
            rc3b_script = write_slurm(rc3b_dir, rc3b_cfg,
                                      f"{name}_r3b", SLURM_TIME_RC3)
            if not dry_run:
                jid = submit(rc3b_script)
                jobs.append((name, "rc3b", jid))
            total_rc3 += 1

            # ── RC3-ArmC: SHH inhibition + drug (56 days) ──
            rc3c_dir = vdir / "rc3_armC"
            rc3c_out = rc3c_dir / "output"
            rc3c_out.mkdir(parents=True, exist_ok=True)
            rc3c_cfg = rc3c_dir / "config.xml"
            patch_config(BASE_CONFIG, rc3c_cfg, rc3c_out, SEED, RC3_MAX_TIME,
                         uptake, kill, hif, cap, off, tgfb,
                         is_rc2=False, rc3_mode='shh_drug')
            rc3c_script = write_slurm(rc3c_dir, rc3c_cfg,
                                      f"{name}_r3c", SLURM_TIME_RC3)
            if not dry_run:
                jid = submit(rc3c_script)
                jobs.append((name, "rc3c", jid))
            total_rc3 += 1

    # Write job manifest
    manifest = SWEEP_DIR / "jobs.txt"
    with open(manifest, "w") as f:
        for name, rc, jid in jobs:
            f.write(f"{name}\t{rc}\t{jid}\n")

    # Write variant table
    vtable = SWEEP_DIR / "variants.tsv"
    with open(vtable, "w") as f:
        f.write("name\tdrug_uptake\tdrug_kill\thif1a_emt\tecm_emt_cap\t"
                "emt_off_thresh\ttgfb_secr\n")
        for name, uptake, kill, hif, cap, off, tgfb in variants:
            f.write(f"{name}\t{uptake}\t{kill}\t{hif}\t{cap}\t{off}\t{tgfb}\n")

    total = total_rc1 + total_rc2 + total_rc3
    print(f"\n  ════════════════════════════════════════════════════════")
    print(f"  MEGA SWEEP: {len(variants)} variants submitted")
    print(f"  ════════════════════════════════════════════════════════")
    print(f"    RC1 jobs:  {total_rc1}")
    print(f"    RC2 jobs:  {total_rc2}")
    if include_rc3:
        print(f"    RC3 jobs:  {total_rc3} (3 arms × {len(variants)})")
    print(f"    TOTAL:     {total} SLURM jobs")
    print(f"  ────────────────────────────────────────────────────────")
    print(f"  Output: {SWEEP_DIR}")
    print(f"  Variants table: {vtable}")
    print(f"  Job manifest:   {manifest}")
    print(f"  ────────────────────────────────────────────────────────")

    # Estimate completion time
    # Assumptions: 8 concurrent slots, RC1=30min, RC2=90min, RC3=150min
    rc1_hrs = (total_rc1 / 8) * 0.5
    rc2_hrs = (total_rc2 / 8) * 1.5
    rc3_hrs = (total_rc3 / 8) * 2.5 if include_rc3 else 0
    # Jobs interleave, so total is roughly max of the two phases
    total_hrs = rc1_hrs + rc2_hrs + rc3_hrs
    # In practice ~60% of sequential due to interleaving
    est_hrs = total_hrs * 0.6

    print(f"\n  ESTIMATED COMPLETION TIME:")
    print(f"    Sequential estimate: {total_hrs:.0f} hours")
    print(f"    With interleaving:   ~{est_hrs:.0f} hours ({est_hrs/24:.1f} days)")
    print(f"    (assumes 8 concurrent SLURM slots)")

    print(f"\n  Monitor with:")
    print(f"    watch -n 60 bash {PROJECT_ROOT}/watch_mega_sweep.sh")
    print(f"\n  Evaluate with:")
    print(f"    bash {PROJECT_ROOT}/eval_mega_sweep.sh")

    if dry_run:
        print(f"\n  *** DRY RUN — no jobs were submitted ***")


if __name__ == "__main__":
    main()
