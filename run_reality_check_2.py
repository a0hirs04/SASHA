#!/usr/bin/env python3
"""
Reality Check 2 — Drug Response, Sanctuary, Resistance & Regrowth
==================================================================

Runs a 42-day single-replicate simulation per seed:
  Phase 0 (day  0–14, t=      0–20160): No drug.  Barrier establishes.
  Phase 1 (day 14–28, t=20160–40320):  Drug at max dose (Dirichlet BC = 1.0).
  Phase 2 (day 28–42, t=40320–60480):  Drug withdrawn.  Observe regrowth.

Drug is applied/withdrawn by the drug scheduling code in custom_function
(drug_start_time, drug_end_time XML params control the Dirichlet BCs).

CRITERIA (must pass in >=4 of 5 replicates unless noted):
  RC2-1  Partial response:   tumor(t_treat_end) <= 0.90 * tumor(t_pre)
                              (meaningful ≥10% reduction, not mere fluctuation)
  RC2-2  Not eradicated:     tumor(t_treat_end) > 0
  RC2-3  Barrier persists:   peri_ecm(t_treat_end) >= 0.80 * peri_ecm(t_pre)
  RC2-4  Spatial sanctuary:  mean_ecm@survivors(t_treat_end) >
                              mean_ecm@all_tumor(t_pre)
  RC2-5  ABCB1 emerges:      frac_abcb1(t_treat_end) > frac_abcb1(t_pre)
                              [SOFT — informational, does not block PASS]
  RC2-6  Regrowth:           tumor(t_post) > tumor(t_treat_end)

FAIL DEBUGGING GUIDE (per spec):
  Drug too effective (eradication)  → ECM→drug diffusion too weak.
                                       Increase ECM impedance weights.
  Drug ineffective (no reduction)   → M7→M4 kill signal broken.  Check E24.
  No regrowth                       → Survivors lost secretion capacity.
                                       Check Rule 5 daughter cells.
  Barrier collapsed                 → CAFs died.  Check Rule 31 (drug-indifferent)
                                       and Rule 32 (zero death rate).
"""
from __future__ import annotations

import argparse
import math
import os
import shutil
import subprocess
import sys
import textwrap
import time
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(PROJECT_ROOT))

from python.wrapper.output_parser import OutputParser
from python.wrapper.workdir_utils import default_reality_check_dir

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
BINARY       = PROJECT_ROOT / "stroma_world"
BASE_CONFIG  = PROJECT_ROOT / "config" / "PhysiCell_settings.xml"

SEEDS            = [42, 99, 137, 256, 1001]
NUM_REPLICATES   = len(SEEDS)
QUORUM           = 4

# Phase boundaries (minutes)
T_PRE           = 20160.0   # end of barrier-forming phase  (day 14)
T_TREAT_END     = 40320.0   # end of treatment phase        (day 28)
T_POST          = 60480.0   # end of regrowth observation   (day 42)
MAX_TIME        = T_POST

DRUG_CONCENTRATION = 1.0    # Dirichlet BC value when drug is on
SAVE_INTERVAL      = 360    # minutes between snapshots (~168 total)

PARTIAL_RESPONSE_MIN_REDUCTION = 0.10   # RC2-1: ≥10% tumor reduction required
BOUNDARY_BARRIER_FRACTION = 0.80        # RC2-3: ECM must stay ≥ 80% of pre-treatment

# HPC SLURM
SLURM_PARTITION    = os.environ.get("RC2_SLURM_PARTITION", os.environ.get("SLURM_PARTITION", "cpu384g"))
SLURM_CPUS         = int(os.environ.get("RC2_SLURM_CPUS", os.environ.get("SLURM_CPUS", "32")))
SLURM_MEM          = os.environ.get("RC2_SLURM_MEM", os.environ.get("SLURM_MEM", "0"))
SLURM_TIME         = os.environ.get("RC2_SLURM_TIME", "06:00:00")
SLURM_POLL_INTERVAL = 20   # seconds

TERMINAL_STATES = {
    "COMPLETED", "FAILED", "CANCELLED", "TIMEOUT",
    "OUT_OF_MEMORY", "NODE_FAIL", "PREEMPTED",
    "BOOT_FAIL", "DEADLINE",
}

# Stromal custom-data column mapping (same offset trick as RC1)
STROMAL_TO_TUMOR_LABEL: Dict[str, str] = {
    "acta2_active":   "ZEB1",
    "gli1_active":    "CDH1",
}


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------
@dataclass
class Snapshot:
    time: float
    n_tumor: int
    n_caf: int
    peri_ecm: float            # mean ECM in peritumoral shell (0–100 µm beyond radius)
    frac_abcb1: float          # fraction of live tumor with abcb1_active > 0.5
    ecm_at_tumor: float        # mean ECM density at live tumor cell positions
    tumor_radius: float


@dataclass
class ReplicateResult:
    seed: int
    run_dir: Path
    success: bool
    wall_time_s: float = 0.0
    snap_pre: Optional[Snapshot] = None
    snap_treat_end: Optional[Snapshot] = None
    snap_post: Optional[Snapshot] = None
    criteria: Dict[str, bool] = field(default_factory=dict)
    details: Dict[str, str] = field(default_factory=dict)


CRITERIA_NAMES = [
    ("partial_response",   "RC2-1  Partial response (tumor reduced during treatment)"),
    ("not_eradicated",     "RC2-2  Not eradicated (survivors at treatment end)"),
    ("barrier_persists",   "RC2-3  Barrier persists (ECM ≥80% of pre-treatment)"),
    ("spatial_sanctuary",  "RC2-4  Spatial sanctuary (survivors in high-ECM zones)"),
    ("abcb1_emerges",      "RC2-5  ABCB1 resistance emerges [SOFT]"),
    ("regrowth",           "RC2-6  Regrowth after withdrawal"),
]


# ---------------------------------------------------------------------------
# Config patching
# ---------------------------------------------------------------------------
def _patch_config(
    src: Path,
    dst: Path,
    output_dir: Path,
    seed: int,
    drug_kill_multiplier: Optional[float] = None,
    abcb1_production_rate: Optional[float] = None,
    drug_stress_threshold: Optional[float] = None,
) -> None:
    tree = ET.parse(src)
    root = tree.getroot()

    def _set(xpath: str, value: str) -> None:
        node = root.find(xpath)
        if node is not None:
            node.text = value

    _set("./save/folder", str(output_dir))
    _set(".//overall/max_time", str(MAX_TIME))
    _set(".//options/random_seed", str(seed))
    _set(".//parallel/omp_num_threads", str(SLURM_CPUS))

    # Save interval
    for n in root.findall(".//save//interval"):
        n.text = str(SAVE_INTERVAL)

    # Drug timing
    _set(".//user_parameters/drug_start_time",   str(T_PRE))
    _set(".//user_parameters/drug_end_time",     str(T_TREAT_END))
    _set(".//user_parameters/drug_concentration", str(DRUG_CONCENTRATION))

    # Optional sweep/transition overrides
    if drug_kill_multiplier is not None:
        _set(".//user_parameters/drug_kill_multiplier", str(drug_kill_multiplier))
    if abcb1_production_rate is not None:
        _set(".//user_parameters/abcb1_production_rate", str(abcb1_production_rate))
    if drug_stress_threshold is not None:
        _set(".//user_parameters/drug_stress_threshold", str(drug_stress_threshold))

    # Enable Dirichlet BCs for all substrates (drug starts at 0; C++ will raise it)
    def _enforce_dirichlet(var_name: str, value: str) -> None:
        var = root.find(f".//microenvironment_setup/variable[@name='{var_name}']")
        if var is None:
            return
        dbc = var.find("./Dirichlet_boundary_condition")
        if dbc is not None:
            dbc.text = value
            dbc.set("enabled", "true")
        for bv in var.findall("./Dirichlet_options/boundary_value"):
            bv.text = value
            bv.set("enabled", "true")

    _enforce_dirichlet("oxygen",      "38")
    _enforce_dirichlet("tgfb",        "0")
    _enforce_dirichlet("shh",         "0")
    _enforce_dirichlet("drug",        "0")   # C++ will raise to drug_conc at T_PRE

    tree.write(dst, encoding="utf-8", xml_declaration=True)


# ---------------------------------------------------------------------------
# SLURM helpers
# ---------------------------------------------------------------------------
def _write_slurm_script(
    rep_dir: Path,
    config_path: Path,
    rep_idx: int,
    seed: int,
    job_name_prefix: str = "rc2",
    drug_pipeline_debug: bool = False,
) -> Path:
    script = rep_dir / "run.slurm.sh"
    debug_env_line = "export RC2_DRUG_DEBUG=1\n" if drug_pipeline_debug else ""
    script.write_text(textwrap.dedent(f"""\
        #!/bin/bash
        #SBATCH --job-name={job_name_prefix}_r{rep_idx+1}_s{seed}
        #SBATCH --partition={SLURM_PARTITION}
        #SBATCH --nodes=1
        #SBATCH --ntasks-per-node=1
        #SBATCH --cpus-per-task={SLURM_CPUS}
        #SBATCH --mem={SLURM_MEM}
        #SBATCH --time={SLURM_TIME}
        #SBATCH --output={rep_dir}/slurm_%j.out
        #SBATCH --error={rep_dir}/slurm_%j.err

        set -euo pipefail
        export OMP_NUM_THREADS={SLURM_CPUS}
        export OMP_PROC_BIND=spread
        export OMP_PLACES=cores
        {debug_env_line.rstrip()}

        cd {PROJECT_ROOT}
        echo "=== Reality Check 2 -- Replicate {rep_idx+1} -- seed={seed} ==="
        echo "Start: $(date)"
        {BINARY} {config_path}
        echo "End:   $(date)"
    """), encoding="utf-8")
    script.chmod(0o755)
    return script


def _query_job(job_id: str) -> Tuple[str, int, bool]:
    try:
        r = subprocess.run(
            ["sacct", "-j", job_id, "--format=State,ExitCode",
             "--noheader", "--parsable2"],
            capture_output=True, text=True, check=False, timeout=30,
        )
        for line in r.stdout.strip().splitlines():
            parts = line.strip().split("|")
            if len(parts) >= 2:
                state = parts[0].strip().split()[0] if parts[0].strip() else "UNKNOWN"
                try:
                    exit_code = int(parts[1].split(":")[0])
                except (ValueError, IndexError):
                    exit_code = -1
                return state, exit_code, state in TERMINAL_STATES
    except Exception:
        pass
    try:
        r = subprocess.run(
            ["squeue", "-j", job_id, "--format=%T", "--noheader"],
            capture_output=True, text=True, check=False, timeout=15,
        )
        state = r.stdout.strip()
        if not state:
            return "COMPLETED", 0, True
        return state, 0, False
    except Exception:
        return "UNKNOWN", -1, False


# ---------------------------------------------------------------------------
# Spatial helpers
# ---------------------------------------------------------------------------
def _sample_nearest(
    points: np.ndarray,
    voxel_coords: np.ndarray,
    values: np.ndarray,
) -> np.ndarray:
    if points.size == 0 or voxel_coords.size == 0:
        return np.array([], dtype=float)
    out = np.empty(points.shape[0], dtype=float)
    chunk = 500
    for i in range(0, points.shape[0], chunk):
        p = points[i : i + chunk]
        d2 = np.sum((p[:, None, :] - voxel_coords[None, :, :]) ** 2, axis=2)
        out[i : i + chunk] = values[np.argmin(d2, axis=1)]
    return out


def _distance_to_boundary(coords: np.ndarray, all_coords: np.ndarray) -> np.ndarray:
    mins = np.min(all_coords, axis=0)
    maxs = np.max(all_coords, axis=0)
    d = np.minimum(coords - mins[None, :], maxs[None, :] - coords)
    return np.min(d, axis=1)


# ---------------------------------------------------------------------------
# Snapshot parser — extract one Snapshot from a PhysiCell XML + mat
# ---------------------------------------------------------------------------
def _parse_snapshot(output_dir: Path, target_time: float) -> Optional[Snapshot]:
    parser = OutputParser(output_dir)
    xmls = sorted(output_dir.glob("output*.xml"))
    if not xmls:
        return None

    snaps = []
    for xml in xmls:
        snap = parser._read_physicell_xml(xml)
        snaps.append((float(snap["time"]), snap))
    if not snaps:
        return None

    snap_time, snap = min(snaps, key=lambda s: abs(s[0] - target_time))

    matrix      = snap["cell_matrix"]
    labels      = snap["label_name_map"]
    micro_coords = snap["micro_coords"]
    micro_values = snap["micro_values"]

    def _row(name: str) -> Optional[np.ndarray]:
        e = labels.get(name)
        if e is None:
            return None
        idx = int(e["index"])
        if idx < 0 or idx >= matrix.shape[0]:
            return None
        return matrix[idx, :]

    # Live masks
    cell_type  = _row("cell_type")
    dead       = _row("dead")
    death_model = _row("current_death_model")

    n_cells = matrix.shape[1]
    live_mask = np.ones(n_cells, dtype=bool)
    if dead is not None:
        live_mask &= dead <= 0.5
    if death_model is not None:
        live_mask &= np.rint(death_model).astype(int) != 100

    tumor_type  = 0
    stroma_type = 1
    ctype_int   = np.rint(cell_type).astype(int) if cell_type is not None else np.full(n_cells, -1, dtype=int)

    live_tumor  = live_mask & (ctype_int == tumor_type)
    live_stroma = live_mask & (ctype_int == stroma_type)

    # CAF count (acta2_active via offset map)
    acta2_row = _row("ZEB1")   # acta2_active is at offset 9 = ZEB1 label
    if acta2_row is not None:
        live_caf = live_stroma & (acta2_row > 0.5)
    else:
        live_caf = np.zeros_like(live_stroma)

    n_tumor = int(np.sum(live_tumor))
    n_caf   = int(np.sum(live_caf))

    # Positions
    pos = parser._get_positions(matrix, labels)
    tumor_pos = pos[live_tumor] if (pos.size and live_tumor.any()) else np.empty((0, 3))
    centroid  = np.nanmean(tumor_pos, axis=0) if tumor_pos.shape[0] > 0 else np.zeros(3)
    tumor_radius = (
        float(np.max(np.linalg.norm(tumor_pos - centroid, axis=1)))
        if tumor_pos.shape[0] > 1 else 0.0
    )

    # ECM density
    ecm = micro_values.get("ecm_density")
    peri_ecm = math.nan
    ecm_at_tumor = math.nan
    if ecm is not None and micro_coords.size > 0:
        if tumor_pos.shape[0] > 0:
            # Peritumoral shell: tumor_radius to tumor_radius+100 µm
            vd = np.linalg.norm(micro_coords - centroid, axis=1)
            shell = (vd >= tumor_radius) & (vd <= tumor_radius + 100.0)
            if shell.any():
                peri_ecm = float(np.nanmean(ecm[shell]))
            # Mean ECM at tumor cell locations
            ecm_at_pos = _sample_nearest(tumor_pos, micro_coords, ecm)
            if ecm_at_pos.size > 0:
                ecm_at_tumor = float(np.nanmean(ecm_at_pos))

    # ABCB1 fraction among live tumor cells
    abcb1 = _row("abcb1_active")
    frac_abcb1 = math.nan
    if abcb1 is not None and live_tumor.any():
        frac_abcb1 = float(np.mean(abcb1[live_tumor] > 0.5))

    return Snapshot(
        time=snap_time,
        n_tumor=n_tumor,
        n_caf=n_caf,
        peri_ecm=peri_ecm,
        frac_abcb1=frac_abcb1,
        ecm_at_tumor=ecm_at_tumor,
        tumor_radius=tumor_radius,
    )


# ---------------------------------------------------------------------------
# Evaluate 6 RC2 criteria
# ---------------------------------------------------------------------------
def _evaluate(r: ReplicateResult) -> ReplicateResult:
    output_dir = r.run_dir / "output"
    print(f"    Parsing pre-treatment snapshot (t≈{T_PRE:.0f})...")
    r.snap_pre       = _parse_snapshot(output_dir, T_PRE)
    print(f"    Parsing treatment-end snapshot  (t≈{T_TREAT_END:.0f})...")
    r.snap_treat_end = _parse_snapshot(output_dir, T_TREAT_END)
    print(f"    Parsing post-withdrawal snapshot (t≈{T_POST:.0f})...")
    r.snap_post      = _parse_snapshot(output_dir, T_POST)

    sp = r.snap_pre
    st = r.snap_treat_end
    sw = r.snap_post

    if sp is None or st is None or sw is None:
        r.success = False
        return r

    # ── RC2-1: Partial response (≥10% reduction) ────────────────────────
    reduction = 1.0 - st.n_tumor / max(sp.n_tumor, 1)
    c1 = reduction >= PARTIAL_RESPONSE_MIN_REDUCTION
    r.criteria["partial_response"] = c1
    r.details["partial_response"] = (
        f"tumor(pre={sp.n_tumor}) → tumor(treat_end={st.n_tumor})  "
        f"reduction={reduction:.1%}  threshold≥{PARTIAL_RESPONSE_MIN_REDUCTION:.0%}  "
        f"{'✓' if c1 else '✗'}"
    )

    # ── RC2-2: Not eradicated ────────────────────────────────────────────
    c2 = st.n_tumor > 0
    r.criteria["not_eradicated"] = c2
    r.details["not_eradicated"] = (
        f"survivors at treatment end = {st.n_tumor}  "
        f"({'>0 ✓' if c2 else '=0 — ERADICATED ✗'})"
    )

    # ── RC2-3: Barrier persists ──────────────────────────────────────────
    barrier_threshold = (
        BOUNDARY_BARRIER_FRACTION * sp.peri_ecm
        if math.isfinite(sp.peri_ecm) else math.nan
    )
    c3 = (
        math.isfinite(st.peri_ecm)
        and math.isfinite(barrier_threshold)
        and st.peri_ecm >= barrier_threshold
    )
    r.criteria["barrier_persists"] = c3
    r.details["barrier_persists"] = (
        f"peri_ecm(pre={sp.peri_ecm:.4f}) → peri_ecm(treat_end={st.peri_ecm:.4f})  "
        f"threshold={barrier_threshold:.4f}  "
        f"({'≥ ✓' if c3 else '< ✗'})"
    )

    # ── RC2-4: Spatial sanctuary ─────────────────────────────────────────
    # Surviving cells (treat_end) live in denser ECM than the pre-treatment
    # tumor population average → they were sheltered behind the barrier.
    c4 = (
        math.isfinite(st.ecm_at_tumor)
        and math.isfinite(sp.ecm_at_tumor)
        and st.ecm_at_tumor > sp.ecm_at_tumor
    )
    r.criteria["spatial_sanctuary"] = c4
    r.details["spatial_sanctuary"] = (
        f"ecm@survivors(treat_end={st.ecm_at_tumor:.4f}) vs "
        f"ecm@tumor(pre={sp.ecm_at_tumor:.4f})  "
        f"({'> ✓' if c4 else '≤ ✗'})"
    )

    # ── RC2-5: ABCB1 emerges ─────────────────────────────────────────────
    c5 = (
        math.isfinite(st.frac_abcb1)
        and math.isfinite(sp.frac_abcb1)
        and st.frac_abcb1 > sp.frac_abcb1
    )
    r.criteria["abcb1_emerges"] = c5
    r.details["abcb1_emerges"] = (
        f"frac_abcb1(pre={sp.frac_abcb1:.4f}) → "
        f"frac_abcb1(treat_end={st.frac_abcb1:.4f})  "
        f"({'increased ✓' if c5 else 'NOT increased ✗'})"
    )

    # ── RC2-6: Regrowth ──────────────────────────────────────────────────
    c6 = sw.n_tumor > st.n_tumor
    r.criteria["regrowth"] = c6
    r.details["regrowth"] = (
        f"tumor(treat_end={st.n_tumor}) → tumor(post={sw.n_tumor})  "
        f"({'regrew ✓' if c6 else 'did NOT regrow ✗'})"
    )

    return r


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> int:
    default_work_dir = default_reality_check_dir(PROJECT_ROOT, "reality_check_2")
    parser = argparse.ArgumentParser(description="Run Reality Check 2")
    parser.add_argument(
        "--seeds",
        nargs="+",
        type=int,
        default=SEEDS,
        help="Seed list to run/evaluate (default: full 5-seed gate)",
    )
    parser.add_argument(
        "--work-dir",
        type=Path,
        default=default_work_dir,
        help="Root directory for RC2 outputs",
    )
    parser.add_argument(
        "--quorum",
        type=int,
        default=None,
        help="Pass threshold; defaults to min(default quorum, number of seeds)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print planned replicate directories without submitting jobs",
    )
    parser.add_argument(
        "--evaluate-only",
        action="store_true",
        help="Evaluate existing outputs only; do not delete or submit jobs",
    )
    parser.add_argument(
        "--job-name-prefix",
        type=str,
        default="rc2",
        help="SLURM job-name prefix for submitted replicates",
    )
    parser.add_argument(
        "--drug-kill-multiplier",
        type=float,
        default=None,
        help="Optional override for user_parameters.drug_kill_multiplier",
    )
    parser.add_argument(
        "--abcb1-production-rate",
        type=float,
        default=None,
        help="Optional override for user_parameters.abcb1_production_rate",
    )
    parser.add_argument(
        "--drug-stress-threshold",
        type=float,
        default=None,
        help="Optional override for user_parameters.drug_stress_threshold",
    )
    parser.add_argument(
        "--drug-pipeline-debug",
        action="store_true",
        help="Enable focused RC2 drug pipeline diagnostics (sets RC2_DRUG_DEBUG=1 in SLURM job)",
    )
    args = parser.parse_args()

    if args.dry_run and args.evaluate_only:
        print("ERROR: --dry-run and --evaluate-only are mutually exclusive")
        return 1

    seeds = list(dict.fromkeys(args.seeds))
    if not seeds:
        print("ERROR: at least one seed is required")
        return 1

    num_replicates = len(seeds)
    quorum = args.quorum if args.quorum is not None else min(QUORUM, num_replicates)
    if quorum < 1 or quorum > num_replicates:
        print(f"ERROR: quorum must be in [1, {num_replicates}], got {quorum}")
        return 1

    work_dir = args.work_dir.resolve()
    rep_dirs = [
        work_dir / f"replicate_{i+1:02d}_seed{seed}"
        for i, seed in enumerate(seeds)
    ]

    print("=" * 72)
    print("  REALITY CHECK 2 — Drug Response, Sanctuary, Resistance & Regrowth")
    print(f"  Replicates: {num_replicates}   Seeds: {seeds}")
    print(f"  Phases (min):  barrier 0–{T_PRE:.0f}  |  drug {T_PRE:.0f}–{T_TREAT_END:.0f}  |  recovery {T_TREAT_END:.0f}–{T_POST:.0f}")
    print(f"  HPC: {SLURM_PARTITION} partition, {SLURM_CPUS} CPUs/node")
    print(f"  Quorum: >={quorum}/{num_replicates}")
    print(f"  Work dir: {work_dir}")
    if args.drug_kill_multiplier is not None:
        print(f"  Override: drug_kill_multiplier={args.drug_kill_multiplier}")
    if args.abcb1_production_rate is not None:
        print(f"  Override: abcb1_production_rate={args.abcb1_production_rate}")
    if args.drug_stress_threshold is not None:
        print(f"  Override: drug_stress_threshold={args.drug_stress_threshold}")
    print("=" * 72)

    if args.dry_run:
        print("  Dry run only; no jobs submitted.")
        for i, seed in enumerate(seeds):
            print(f"    rep {i+1} seed={seed} -> {rep_dirs[i]}")
        return 0

    if not args.evaluate_only:
        if not BINARY.exists():
            print(f"ERROR: binary not found at {BINARY}")
            return 1
        if not BASE_CONFIG.exists():
            print(f"ERROR: config not found at {BASE_CONFIG}")
            return 1
        if work_dir.exists():
            shutil.rmtree(work_dir)
        work_dir.mkdir(parents=True)
    elif not work_dir.exists():
        print(f"ERROR: work dir does not exist for evaluate-only mode: {work_dir}")
        return 1

    # ---- Prepare & submit ----
    job_ids: Dict[str, int] = {}
    completed: Dict[int, Tuple[str, int]] = {}
    total_wall = 0.0

    if args.evaluate_only:
        print("\n  Evaluate-only mode — skipping submission and reading existing output.\n")
        completed = {i: ("COMPLETED", 0) for i in range(num_replicates)}
    else:
        for i, seed in enumerate(seeds):
            rep_dir = rep_dirs[i]
            rep_dir.mkdir(parents=True, exist_ok=True)
            output_dir = rep_dir / "output"
            output_dir.mkdir(parents=True, exist_ok=True)

            config_path = rep_dir / "config.xml"
            _patch_config(
                BASE_CONFIG,
                config_path,
                output_dir,
                seed,
                drug_kill_multiplier=args.drug_kill_multiplier,
                abcb1_production_rate=args.abcb1_production_rate,
                drug_stress_threshold=args.drug_stress_threshold,
            )

            # Copy tumour calibration knobs into run dir so binary can find them
            local_cfg = rep_dir / "config"
            local_cfg.mkdir(exist_ok=True)
            shutil.copy(
                PROJECT_ROOT / "config" / "tumor_calibration_knobs.json",
                local_cfg / "tumor_calibration_knobs.json",
            )
            for extra in ["gene_params_default.json"]:
                src = PROJECT_ROOT / "config" / extra
                if src.exists():
                    shutil.copy(src, local_cfg / extra)

            script = _write_slurm_script(
                rep_dir,
                config_path,
                i,
                seed,
                job_name_prefix=args.job_name_prefix,
                drug_pipeline_debug=args.drug_pipeline_debug,
            )

            result = subprocess.run(
                ["sbatch", "--parsable", str(script)],
                capture_output=True, text=True, check=False,
            )
            if result.returncode != 0:
                print(f"  ERROR: sbatch failed for rep {i+1}: {result.stderr.strip()}")
                return 1

            job_id = result.stdout.strip().split(";")[0]
            job_ids[job_id] = i
            print(f"  Submitted rep {i+1} (seed={seed}) → SLURM job {job_id}")

        print(f"\n  All {num_replicates} jobs submitted.  Polling ...\n")

        # ---- Poll ----
        t0 = time.perf_counter()

        while len(completed) < len(job_ids):
            time.sleep(SLURM_POLL_INTERVAL)
            elapsed = time.perf_counter() - t0

            for job_id, rep_idx in list(job_ids.items()):
                if rep_idx in completed:
                    continue
                state, exit_code, terminal = _query_job(job_id)
                if terminal:
                    completed[rep_idx] = (state, exit_code)
                    seed = seeds[rep_idx]
                    tag = "OK" if (state == "COMPLETED" and exit_code == 0) else "XX"
                    print(f"  {tag} Rep {rep_idx+1} (seed={seed}): {state}  exit={exit_code}  [{elapsed:.0f}s]")

            still_running = num_replicates - len(completed)
            if still_running > 0:
                for job_id, rep_idx in job_ids.items():
                    if rep_idx in completed:
                        continue
                    n_snaps = len(list(rep_dirs[rep_idx].glob("output/output*.xml")))
                    if n_snaps > 0:
                        print(f"    rep {rep_idx+1}: {n_snaps} snapshots  [{elapsed:.0f}s]")

        total_wall = time.perf_counter() - t0
        print(f"\n  All jobs done in {total_wall:.1f}s ({total_wall/60:.1f} min)\n")

    # ---- Evaluate ----
    results: List[ReplicateResult] = []

    for i, seed in enumerate(seeds):
        state, exit_code = completed.get(i, ("UNKNOWN", -1))
        success = state == "COMPLETED" and exit_code == 0
        output_dir = rep_dirs[i] / "output"
        if success and not list(output_dir.glob("output*.xml")):
            success = False

        r = ReplicateResult(seed=seed, run_dir=rep_dirs[i], success=success,
                            wall_time_s=total_wall / max(num_replicates, 1))
        if success:
            print(f"\n  Evaluating rep {i+1} (seed={seed}) ...")
            try:
                r = _evaluate(r)
            except Exception as exc:
                import traceback
                print(f"    ERROR: {exc}")
                traceback.print_exc()
                r.success = False
        results.append(r)

    # ---- Per-replicate summary ----
    print()
    print("=" * 72)
    print("  PER-REPLICATE RESULTS")
    print("=" * 72)

    for i, r in enumerate(results):
        hdr = f"Replicate {i+1}  (seed={r.seed})"
        if not r.success:
            print(f"\n  {hdr}  — SIMULATION FAILED")
            continue
        all_pass = all(r.criteria.get(k, False) for k, _ in CRITERIA_NAMES)
        tag = "ALL PASS" if all_pass else "SOME FAIL"
        print(f"\n  {hdr}  — {tag}")
        for key, label in CRITERIA_NAMES:
            ok = r.criteria.get(key, False)
            mark = "PASS" if ok else "FAIL"
            print(f"    [{mark}] {label}")
            print(f"           {r.details.get(key, '')}")
        if r.snap_pre and r.snap_treat_end and r.snap_post:
            sp, st, sw = r.snap_pre, r.snap_treat_end, r.snap_post
            print(f"    Timecourse:")
            print(f"      t={sp.time:.0f}  tumor={sp.n_tumor:4d}  caf={sp.n_caf:4d}  "
                  f"peri_ecm={sp.peri_ecm:.4f}  abcb1={sp.frac_abcb1:.4f}  "
                  f"ecm@tumor={sp.ecm_at_tumor:.4f}")
            print(f"      t={st.time:.0f}  tumor={st.n_tumor:4d}  caf={st.n_caf:4d}  "
                  f"peri_ecm={st.peri_ecm:.4f}  abcb1={st.frac_abcb1:.4f}  "
                  f"ecm@tumor={st.ecm_at_tumor:.4f}")
            print(f"      t={sw.time:.0f}  tumor={sw.n_tumor:4d}  caf={sw.n_caf:4d}  "
                  f"peri_ecm={sw.peri_ecm:.4f}  abcb1={sw.frac_abcb1:.4f}  "
                  f"ecm@tumor={sw.ecm_at_tumor:.4f}")

    # ---- Aggregate ----
    print()
    print("=" * 72)
    print(f"  AGGREGATE  (>={quorum}/{num_replicates} replicates must pass)")
    print("=" * 72)

    SOFT_CRITERIA = {"abcb1_emerges"}   # informational, doesn't block PASS
    any_failure = False
    for key, label in CRITERIA_NAMES:
        passing = sum(1 for r in results if r.success and r.criteria.get(key, False))
        is_soft = key in SOFT_CRITERIA
        ok = passing >= quorum
        if is_soft:
            mark = "INFO" if not ok else "PASS"
        else:
            mark = "PASS" if ok else "FAIL"
            if not ok:
                any_failure = True
        print(f"  [{mark}] {label}  ({passing}/{num_replicates})")

    print()
    print(f"  Total wall time: {total_wall:.1f}s  ({total_wall/60:.1f} min)")
    if any_failure:
        print("\n  *** REALITY CHECK 2: FAIL ***")
        return 1
    else:
        print("\n  *** REALITY CHECK 2: PASS ***")
        return 0


if __name__ == "__main__":
    sys.exit(main())
