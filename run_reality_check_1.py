#!/usr/bin/env python3
"""
Reality Check 1 — Emergent Stroma Barrier Self-Assembly
=======================================================

HPC-accelerated version.  Submits all 5 replicates as parallel SLURM
jobs (128 CPUs each on the compute partition), then collects outputs
and evaluates 8 biological criteria.

Setup
-----
- 50 tumor cells at domain centre (canonical PDAC genotype)
- 200 quiescent PSCs in surrounding ring
- No drug.  No intervention.
- Simulated time: 7 days (10 080 min)
- 5 replicates with different random seeds

Criteria (must pass in ≥4 of 5 replicates)
-------------------------------------------
1. Barrier self-assembles           – mean peritumoral ECM > 0.5
2. Spatial ECM gradient             – ECM near tumor > ECM at edge
3. Central hypoxia develops         – O₂ centroid < 0.5 × O₂ boundary
4. EMT is peripheral, not global    – EMT fraction periphery > core
5. PDAC-like tumor-stroma ratio     – tumor fraction < 50% by count
6. No numerical instability         – no NaN, no negative conc, ECM ≤ 1.0
7. CAF count increases from 0       – activated_cafs_end > 0
8. TGF-β highest near tumor         – mean TGF-β near tumor > at edge
"""
from __future__ import annotations

import json
import math
import os
import re
import shutil
import subprocess
import sys
import textwrap
import time
import xml.etree.ElementTree as ET
import random
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

# ---------------------------------------------------------------------------
# Project root & imports
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(PROJECT_ROOT))

from python.wrapper.output_parser import OutputParser, SimulationMetrics

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
BINARY = PROJECT_ROOT / "stroma_world"
BASE_CONFIG = PROJECT_ROOT / "config" / "PhysiCell_settings.xml"

NUM_REPLICATES = 5
SEEDS = [42, 99, 137, 256, 1001]
SIM_MAX_TIME_MIN = 10_080.0          # 7 days (28 snapshot intervals of 360 min)
QUORUM = 4                            # >=4 of 5 must pass
BOUNDARY_O2 = 38.0                    # mmHg

# Time-series checkpoints (minutes)
TARGET_TIMES = [1440, 2880, 4320, 7200, 10080]

# HPC SLURM settings
SLURM_PARTITION = "compute"
SLURM_CPUS = 128                      # full node
SLURM_MEM = "128G"
SLURM_TIME = "04:00:00"               # 4h wall time limit per replicate
SLURM_POLL_INTERVAL = 15              # seconds between sacct polls


@dataclass
class ReplicateResult:
    seed: int
    run_dir: Path
    success: bool
    wall_time_s: float
    metrics: Optional[SimulationMetrics] = None
    criteria: Dict[str, bool] = field(default_factory=dict)
    details: Dict[str, str] = field(default_factory=dict)
    spatial: Dict[str, float] = field(default_factory=dict)
    timeseries: List[Dict[str, Any]] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Config patching
# ---------------------------------------------------------------------------

def _patch_config(
    src: Path,
    dst: Path,
    output_dir: Path,
    seed: int,
    max_time: float,
    omp_threads: int = SLURM_CPUS,
) -> None:
    """Copy base XML and patch seed/max_time/output folder for RC1."""
    tree = ET.parse(src)
    root = tree.getroot()

    # Set output folder
    folder = root.find("./save/folder")
    if folder is not None:
        folder.text = str(output_dir)

    # Set max time
    mt = root.find(".//overall/max_time")
    if mt is not None:
        mt.text = str(max_time)

    # Set random seed
    rs = root.find(".//options/random_seed")
    if rs is not None:
        rs.text = str(seed)

    # Set OMP threads to match SLURM allocation
    omp = root.find(".//parallel/omp_num_threads")
    if omp is not None:
        omp.text = str(omp_threads)

    # Reduce output frequency (every 360 min instead of 60)
    for interval_node in root.findall(".//save//interval"):
        interval_node.text = "360"

    # Enforce Dirichlet boundaries for all 5 substrates.
    # RC1 keeps drug concentration at 0, but BC stays enabled.
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

    _enforce_dirichlet("oxygen", "38")
    _enforce_dirichlet("tgfb", "0")
    _enforce_dirichlet("shh", "0")
    _enforce_dirichlet("drug", "0")
    _enforce_dirichlet("ecm_density", "0")

    # Push drug_start_time beyond sim end
    ds = root.find(".//user_parameters/drug_start_time")
    if ds is not None:
        ds.text = str(max_time + 10000)

    dc = root.find(".//user_parameters/drug_concentration")
    if dc is not None:
        dc.text = "0"

    tree.write(dst, encoding="utf-8", xml_declaration=True)


def _compute_live_population_snapshot(snapshot: Dict[str, Any]) -> Dict[str, Any]:
    """Count live populations directly from the snapshot cell table."""
    matrix = snapshot["cell_matrix"]
    labels = snapshot["label_name_map"]

    def _row(name: str):
        entry = labels.get(name)
        if entry is None:
            return None
        idx = int(entry["index"])
        if idx < 0 or idx >= matrix.shape[0]:
            return None
        return matrix[idx, :]

    cell_type = _row("cell_type")
    if cell_type is None:
        n = matrix.shape[1]
        cell_type_int = np.full(n, -1, dtype=int)
    else:
        cell_type_int = np.rint(cell_type).astype(int)

    dead = _row("dead")
    death_model = _row("current_death_model")
    live_mask = np.ones(cell_type_int.shape[0], dtype=bool)
    if dead is not None:
        live_mask &= dead <= 0.5
    if death_model is not None:
        live_mask &= np.rint(death_model).astype(int) != 100

    tumor_mask = cell_type_int == 0
    stroma_mask = cell_type_int == 1
    live_tumor_mask = tumor_mask & live_mask
    live_stroma_mask = stroma_mask & live_mask

    acta2 = _row("acta2_active")
    if acta2 is None:
        acta2 = _row("is_activated")
    if acta2 is None:
        acta2 = _row("ACTA2")
    live_caf_mask = live_stroma_mask & (acta2 > 0.5) if acta2 is not None else np.zeros_like(live_stroma_mask)

    is_mes = _row("is_mesenchymal")
    if is_mes is None:
        is_mes = _row("zeb1_active")
    if is_mes is not None and np.any(live_tumor_mask):
        emt_fraction_live = float(np.mean(is_mes[live_tumor_mask] > 0.5))
    else:
        emt_fraction_live = math.nan

    return {
        "live_mask": live_mask,
        "tumor_mask": tumor_mask,
        "stroma_mask": stroma_mask,
        "live_tumor_mask": live_tumor_mask,
        "live_stroma_mask": live_stroma_mask,
        "live_caf_mask": live_caf_mask,
        "n_live_tumor": int(np.sum(live_tumor_mask)),
        "n_live_stroma": int(np.sum(live_stroma_mask)),
        "n_live_caf": int(np.sum(live_caf_mask)),
        "emt_fraction_live_tumor": emt_fraction_live,
    }


def _print_t360_proliferation_diagnostic(
    output_dir: Path,
    config_path: Path,
    seed: int,
) -> None:
    """One-time diagnostic at first snapshot (~t=360 min) for tumor proliferation rates."""
    parser = OutputParser(output_dir)
    xml_files = sorted(output_dir.glob("output*.xml"))
    if not xml_files:
        print(f"[DIAG_T360_PROLIF] No snapshots in {output_dir}")
        return

    snapshots: List[Tuple[float, Path, Dict[str, Any]]] = []
    for xml in xml_files:
        snap = parser._read_physicell_xml(xml)
        snapshots.append((float(snap["time"]), xml, snap))
    snap_time, xml_path, snap = min(snapshots, key=lambda s: abs(s[0] - 360.0))

    matrix = snap["cell_matrix"]
    labels = snap["label_name_map"]

    def _row(name: str):
        entry = labels.get(name)
        if entry is None:
            return None
        idx = int(entry["index"])
        if idx < 0 or idx >= matrix.shape[0]:
            return None
        return matrix[idx, :]

    population = _compute_live_population_snapshot(snap)
    live_tumor_mask = population["live_tumor_mask"]

    rate = _row("current_cycle_phase_exit_rate")
    cell_id = _row("cell_id")
    parent_id = _row("parent_id")
    if cell_id is None:
        cell_id = _row("ID")

    tumor_indices = np.where(live_tumor_mask)[0].tolist()
    rng = random.Random(seed + 360)
    rng.shuffle(tumor_indices)
    sampled_tumors = tumor_indices[: min(10, len(tumor_indices))]

    daughter_indices: List[int] = []
    if parent_id is not None:
        daughter_indices = [i for i in tumor_indices if parent_id[i] >= 0.0]
        rng.shuffle(daughter_indices)
        daughter_indices = daughter_indices[: min(3, len(daughter_indices))]

    # Read configured rates to identify source explicitly.
    base_rate = math.nan
    generic_rate = math.nan
    try:
        cfg_root = ET.parse(config_path).getroot()
        base_node = cfg_root.find(".//user_parameters/base_proliferation_rate")
        generic_node = cfg_root.find(".//user_parameters/proliferation_rate")
        if base_node is not None and base_node.text is not None:
            base_rate = float(base_node.text)
        if generic_node is not None and generic_node.text is not None:
            generic_rate = float(generic_node.text)
    except Exception:
        pass

    print(f"[DIAG_T360_PROLIF] snapshot={xml_path.name} time={snap_time:.1f} min")
    print(
        "[DIAG_T360_PROLIF] module4 parameter source: "
        "base_proliferation_rate (read_xml_double_or_default), "
        f"base_proliferation_rate={base_rate}, proliferation_rate={generic_rate}"
    )
    print("[DIAG_T360_PROLIF] 10 random live tumor cells (idx, cell_id, parent_id, transition_rate):")
    for idx in sampled_tumors:
        cid = float(cell_id[idx]) if cell_id is not None else float(idx)
        pid = float(parent_id[idx]) if parent_id is not None else math.nan
        tr = float(rate[idx]) if rate is not None else math.nan
        print(f"[DIAG_T360_PROLIF] tumor idx={idx} cell_id={cid:.0f} parent_id={pid:.0f} transition_rate={tr:.6f}")

    print("[DIAG_T360_PROLIF] 3 random live daughter tumor cells (parent_id != -1):")
    if not daughter_indices:
        print("[DIAG_T360_PROLIF] daughters: none in snapshot")
    for idx in daughter_indices:
        cid = float(cell_id[idx]) if cell_id is not None else float(idx)
        pid = float(parent_id[idx]) if parent_id is not None else math.nan
        tr = float(rate[idx]) if rate is not None else math.nan
        print(f"[DIAG_T360_PROLIF] daughter idx={idx} cell_id={cid:.0f} parent_id={pid:.0f} transition_rate={tr:.6f}")


# ---------------------------------------------------------------------------
# SLURM helpers
# ---------------------------------------------------------------------------

def _write_slurm_script(
    rep_dir: Path,
    config_path: Path,
    intervention_path: Path,
    rep_idx: int,
    seed: int,
) -> Path:
    """Write a SLURM batch script for one replicate."""
    script = rep_dir / "run.slurm.sh"
    script.write_text(textwrap.dedent(f"""\
        #!/bin/bash
        #SBATCH --job-name=rc1_rep{rep_idx+1}_s{seed}
        #SBATCH --partition={SLURM_PARTITION}
        #SBATCH --nodes=1
        #SBATCH --ntasks-per-node=1
        #SBATCH --cpus-per-task={SLURM_CPUS}
        #SBATCH --mem={SLURM_MEM}
        #SBATCH --time={SLURM_TIME}
        #SBATCH --output={rep_dir}/slurm_%j.out
        #SBATCH --error={rep_dir}/slurm_%j.err

        set -euo pipefail

        # Environment
        if command -v module >/dev/null 2>&1; then
            module purge || true
            module load gcc/12 2>/dev/null || true
        fi

        export OMP_NUM_THREADS={SLURM_CPUS}
        export OMP_PROC_BIND=spread
        export OMP_PLACES=cores

        cd {PROJECT_ROOT}

        echo "=== Reality Check 1 -- Replicate {rep_idx+1} -- seed={seed} ==="
        echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"
        echo "Binary: {BINARY}"
        echo "Config: {config_path}"
        echo "Start:  $(date)"

        # Run simulation
        {BINARY} {config_path} {intervention_path}
        EXIT_CODE=$?

        echo "Exit:   $EXIT_CODE"
        echo "End:    $(date)"
        exit $EXIT_CODE
    """), encoding="utf-8")
    script.chmod(0o755)
    return script


TERMINAL_STATES = {
    "COMPLETED", "FAILED", "CANCELLED", "TIMEOUT",
    "OUT_OF_MEMORY", "NODE_FAIL", "PREEMPTED",
    "BOOT_FAIL", "DEADLINE",
}


def _query_job(job_id: str) -> Tuple[str, int, bool]:
    """Query sacct for job state.  Returns (state, exit_code, is_terminal)."""
    try:
        result = subprocess.run(
            ["sacct", "-j", job_id, "--format=State,ExitCode",
             "--noheader", "--parsable2"],
            capture_output=True, text=True, check=False, timeout=30,
        )
        for line in result.stdout.strip().splitlines():
            parts = line.strip().split("|")
            if len(parts) >= 2:
                state = parts[0].strip().split()[0] if parts[0].strip() else "UNKNOWN"
                try:
                    exit_code = int(parts[1].split(":")[0])
                except (ValueError, IndexError):
                    exit_code = -1
                terminal = state in TERMINAL_STATES
                return state, exit_code, terminal
    except Exception:
        pass

    # Fallback: squeue
    try:
        result = subprocess.run(
            ["squeue", "-j", job_id, "--format=%T", "--noheader"],
            capture_output=True, text=True, check=False, timeout=15,
        )
        state = result.stdout.strip()
        if not state:
            return "COMPLETED", 0, True
        return state, 0, False
    except Exception:
        return "UNKNOWN", -1, False


# ---------------------------------------------------------------------------
# Spatial analysis helpers
# ---------------------------------------------------------------------------

def _sample_nearest(
    points: np.ndarray,
    voxel_coords: np.ndarray,
    values: np.ndarray,
) -> np.ndarray:
    """For each point, find the nearest voxel and return its value."""
    if points.size == 0 or voxel_coords.size == 0:
        return np.array([], dtype=float)
    out = np.empty(points.shape[0], dtype=float)
    chunk = 500
    for i in range(0, points.shape[0], chunk):
        p = points[i:i+chunk]
        d2 = np.sum((p[:, None, :] - voxel_coords[None, :, :]) ** 2, axis=2)
        idx = np.argmin(d2, axis=1)
        out[i:i+chunk] = values[idx]
    return out


def _distance_to_boundary(coords: np.ndarray, all_coords: np.ndarray) -> np.ndarray:
    """Min distance from each coord to the domain boundary."""
    mins = np.min(all_coords, axis=0)
    maxs = np.max(all_coords, axis=0)
    d_lo = coords - mins[None, :]
    d_hi = maxs[None, :] - coords
    d = np.minimum(d_lo, d_hi)
    return np.min(d, axis=1)


# ---------------------------------------------------------------------------
# Evaluate 8 criteria on a single replicate
# ---------------------------------------------------------------------------

def _evaluate(result: ReplicateResult) -> ReplicateResult:
    output_dir = result.run_dir / "output"
    parser = OutputParser(output_dir)

    final_xml = parser._find_final_snapshot_xml()
    snapshot = parser._read_physicell_xml(final_xml)
    matrix = snapshot["cell_matrix"]
    labels = snapshot["label_name_map"]
    micro_coords = snapshot["micro_coords"]
    micro_values = snapshot["micro_values"]
    result.metrics = parser.parse_final_state()
    m = result.metrics

    # ---------- cell masks ----------
    def _row(name):
        entry = labels.get(name)
        if entry is None:
            return None
        idx = int(entry["index"])
        if idx < 0 or idx >= matrix.shape[0]:
            return None
        return matrix[idx, :]

    cell_type = _row("cell_type")
    dead = _row("dead")
    death_model = _row("current_death_model")
    is_mes = _row("is_mesenchymal")
    positions = parser._get_positions(matrix, labels)

    pop = _compute_live_population_snapshot(snapshot)
    tumor_mask = pop["live_tumor_mask"]
    stroma_mask = pop["live_stroma_mask"]
    live_mask = pop["live_mask"]
    live_tumor_mask = pop["live_tumor_mask"]

    tumor_pos = positions[tumor_mask] if positions.size else np.empty((0, 3))
    centroid = np.nanmean(tumor_pos, axis=0) if tumor_pos.size else np.zeros(3)
    tumor_radius = (
        float(np.max(np.linalg.norm(tumor_pos - centroid, axis=1)))
        if tumor_pos.size else 0.0
    )

    o2 = micro_values.get("oxygen")
    ecm = micro_values.get("ecm_density")
    tgfb = micro_values.get("tgfb")

    # ===== CRITERION 1: Barrier self-assembles =====
    peri_ecm = math.nan
    if ecm is not None and micro_coords.size and tumor_pos.size:
        voxel_dist = np.linalg.norm(micro_coords - centroid, axis=1)
        shell_mask = (
            (voxel_dist >= tumor_radius)
            & (voxel_dist <= tumor_radius + 100.0)
        )
        if shell_mask.any():
            peri_ecm = float(np.nanmean(ecm[shell_mask]))
    c1 = math.isfinite(peri_ecm) and peri_ecm > 0.5
    result.spatial["peritumoral_ecm"] = peri_ecm
    result.criteria["barrier_self_assembles"] = c1
    result.details["barrier_self_assembles"] = (
        f"peritumoral_ecm={peri_ecm:.4f} {'>' if c1 else '<='} 0.5"
    )

    # ===== CRITERION 2: ECM near tumor > ECM at edge =====
    ecm_near = ecm_edge = math.nan
    if ecm is not None and micro_coords.size and tumor_pos.size:
        voxel_dist = np.linalg.norm(micro_coords - centroid, axis=1)
        near_mask = voxel_dist <= tumor_radius + 150.0
        bd = _distance_to_boundary(micro_coords, micro_coords)
        edge_mask = bd <= 50.0
        if near_mask.any():
            ecm_near = float(np.nanmean(ecm[near_mask]))
        if edge_mask.any():
            ecm_edge = float(np.nanmean(ecm[edge_mask]))
    c2 = (
        math.isfinite(ecm_near) and math.isfinite(ecm_edge)
        and ecm_near > ecm_edge
    )
    result.spatial["ecm_near_tumor"] = ecm_near
    result.spatial["ecm_at_edge"] = ecm_edge
    result.criteria["spatial_ecm_gradient"] = c2
    result.details["spatial_ecm_gradient"] = (
        f"ECM_near={ecm_near:.4f} {'>' if c2 else '<='} ECM_edge={ecm_edge:.4f}"
    )

    # ===== CRITERION 3: Central hypoxia =====
    o2_centroid = math.nan
    o2_boundary_actual = BOUNDARY_O2
    if o2 is not None and micro_coords.size:
        o2_centroid = float(
            _sample_nearest(centroid[None, :], micro_coords, o2)[0]
        )
        bd = _distance_to_boundary(micro_coords, micro_coords)
        edge_mask = bd <= 30.0
        if edge_mask.any():
            o2_boundary_actual = float(np.nanmean(o2[edge_mask]))
    threshold_o2 = 0.5 * o2_boundary_actual
    c3 = math.isfinite(o2_centroid) and o2_centroid < threshold_o2
    result.spatial["o2_centroid"] = o2_centroid
    result.spatial["o2_boundary"] = o2_boundary_actual
    result.criteria["central_hypoxia"] = c3
    result.details["central_hypoxia"] = (
        f"O2_centroid={o2_centroid:.2f} {'<' if c3 else '>='} "
        f"0.5*O2_boundary={threshold_o2:.2f}"
    )

    # ===== CRITERION 4: EMT peripheral > core =====
    emt_core = emt_periph = math.nan
    if is_mes is not None and tumor_pos.size >= 10:
        tumor_mesenchymal = is_mes[tumor_mask] > 0.5
        radial = np.linalg.norm(tumor_pos - centroid, axis=1)
        core_cut = float(np.quantile(radial, 0.40))
        periph_cut = float(np.quantile(radial, 0.60))
        core_m = radial <= core_cut
        periph_m = radial >= periph_cut
        if core_m.any():
            emt_core = float(np.mean(tumor_mesenchymal[core_m]))
        if periph_m.any():
            emt_periph = float(np.mean(tumor_mesenchymal[periph_m]))
    c4 = (
        math.isfinite(emt_core) and math.isfinite(emt_periph)
        and emt_periph > emt_core
    )
    result.spatial["emt_core_fraction"] = emt_core
    result.spatial["emt_periphery_fraction"] = emt_periph
    result.criteria["emt_peripheral_not_global"] = c4
    result.details["emt_peripheral_not_global"] = (
        f"EMT_periph={emt_periph:.4f} {'>' if c4 else '<='} "
        f"EMT_core={emt_core:.4f}"
    )

    # ===== CRITERION 5: Tumor fraction < 50% =====
    total_cells = pop["n_live_tumor"] + pop["n_live_stroma"]
    tumor_frac = pop["n_live_tumor"] / max(total_cells, 1)
    c5 = tumor_frac < 0.50
    result.spatial["tumor_fraction"] = tumor_frac
    result.criteria["pdac_like_ratio"] = c5
    result.details["pdac_like_ratio"] = (
        f"tumor_frac={tumor_frac:.4f} ({'<' if c5 else '>='} 0.50)  "
        f"tumor={pop['n_live_tumor']}  stroma={pop['n_live_stroma']}"
    )

    # ===== CRITERION 6: No numerical instability =====
    issues: List[str] = []
    for attr in (
        "mean_ecm_density", "max_ecm_density", "tumor_extent",
        "drug_penetration", "hypoxic_fraction",
    ):
        v = getattr(m, attr, None)
        if v is not None and not math.isfinite(v):
            issues.append(f"{attr}=NaN")
    if ecm is not None and np.any(ecm > 1.0 + 1e-9):
        issues.append(f"ECM max={np.max(ecm):.6f} > 1.0")
    for fname, fvals in micro_values.items():
        if fvals is not None and np.any(fvals < -1e-12):
            issues.append(
                f"{fname} has negative values (min={np.min(fvals):.6e})"
            )
    if o2 is not None and np.any(o2 > BOUNDARY_O2 + 1e-6):
        issues.append(f"O2 max={np.max(o2):.4f} > boundary {BOUNDARY_O2}")
    c6 = len(issues) == 0
    result.criteria["no_numerical_instability"] = c6
    result.details["no_numerical_instability"] = (
        "clean" if c6 else f"ISSUES: {'; '.join(issues)}"
    )

    # ===== CRITERION 7: CAF activation =====
    c7 = pop["n_live_caf"] > 0
    result.criteria["caf_activation"] = c7
    result.details["caf_activation"] = f"activated_cafs={pop['n_live_caf']}"

    # ===== CRITERION 8: TGF-beta gradient =====
    tgfb_near = tgfb_edge = math.nan
    if tgfb is not None and micro_coords.size and tumor_pos.size:
        voxel_dist = np.linalg.norm(micro_coords - centroid, axis=1)
        near_mask = voxel_dist <= tumor_radius + 150.0
        bd = _distance_to_boundary(micro_coords, micro_coords)
        edge_mask = bd <= 50.0
        if near_mask.any():
            tgfb_near = float(np.nanmean(tgfb[near_mask]))
        if edge_mask.any():
            tgfb_edge = float(np.nanmean(tgfb[edge_mask]))
    c8 = (
        math.isfinite(tgfb_near) and math.isfinite(tgfb_edge)
        and tgfb_near > tgfb_edge
    )
    result.spatial["tgfb_near_tumor"] = tgfb_near
    result.spatial["tgfb_at_edge"] = tgfb_edge
    result.criteria["tgfb_gradient"] = c8
    result.details["tgfb_gradient"] = (
        f"TGFb_near={tgfb_near:.4f} {'>' if c8 else '<='} "
        f"TGFb_edge={tgfb_edge:.4f}"
    )

    return result


def _collect_time_series(output_dir: Path) -> List[Dict[str, Any]]:
    """Collect coarse time-series at target minutes using nearest snapshots."""
    parser = OutputParser(output_dir)
    xml_files = sorted(output_dir.glob("output*.xml"))
    snapshots: List[Tuple[float, Path, Dict[str, Any]]] = []

    for xml in xml_files:
        snap = parser._read_physicell_xml(xml)
        snapshots.append((float(snap["time"]), xml, snap))

    if not snapshots:
        return []

    timeseries: List[Dict[str, Any]] = []
    for target in TARGET_TIMES:
        closest = min(snapshots, key=lambda s: abs(s[0] - target))
        snap_time, xml_path, snap = closest
        pop = _compute_live_population_snapshot(snap)

        matrix = snap["cell_matrix"]
        labels = snap["label_name_map"]
        micro_coords = snap["micro_coords"]
        micro_values = snap["micro_values"]

        def _row(name):
            entry = labels.get(name)
            if entry is None:
                return None
            idx = int(entry["index"])
            if idx < 0 or idx >= matrix.shape[0]:
                return None
            return matrix[idx, :]

        cell_type = _row("cell_type")
        is_mes = _row("is_mesenchymal")
        positions = parser._get_positions(matrix, labels)

        tumor_mask = pop["live_tumor_mask"]
        stroma_mask = pop["live_stroma_mask"]
        tumor_pos = positions[tumor_mask] if positions.size else np.empty((0, 3))
        centroid = np.nanmean(tumor_pos, axis=0) if tumor_pos.size else np.zeros(3)
        tumor_radius = (
            float(np.max(np.linalg.norm(tumor_pos - centroid, axis=1)))
            if tumor_pos.size else 0.0
        )

        ecm = micro_values.get("ecm_density")
        tgfb = micro_values.get("tgfb")

        peri_ecm = math.nan
        if ecm is not None and micro_coords.size and tumor_pos.size:
            voxel_dist = np.linalg.norm(micro_coords - centroid, axis=1)
            shell_mask = (
                (voxel_dist >= tumor_radius)
                & (voxel_dist <= tumor_radius + 100.0)
            )
            if shell_mask.any():
                peri_ecm = float(np.nanmean(ecm[shell_mask]))

        max_tgfb_at_psc = math.nan
        if tgfb is not None and micro_coords.size and stroma_mask.any():
            stroma_pos = positions[stroma_mask]
            nearest = _sample_nearest(stroma_pos, micro_coords, tgfb)
            if nearest.size:
                max_tgfb_at_psc = float(np.nanmax(nearest))

        emt_fraction = pop["emt_fraction_live_tumor"]

        timeseries.append({
            "target_time": float(target),
            "snapshot_time": float(snap_time),
            "total_tumor_cells": pop["n_live_tumor"],
            "total_caf_count": pop["n_live_caf"],
            "max_tgfb_at_psc": max_tgfb_at_psc,
            "mean_peritumoral_ecm": peri_ecm,
            "emt_fraction": emt_fraction,
        })

    return timeseries


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

CRITERIA_NAMES = [
    ("barrier_self_assembles",    "1. Barrier self-assembles"),
    ("spatial_ecm_gradient",      "2. Spatial ECM gradient"),
    ("central_hypoxia",           "3. Central hypoxia develops"),
    ("emt_peripheral_not_global", "4. EMT peripheral, not global"),
    ("pdac_like_ratio",           "5. Tumor-stroma ratio PDAC-like"),
    ("no_numerical_instability",  "6. No numerical instability"),
    ("caf_activation",            "7. CAF count increases"),
    ("tgfb_gradient",             "8. TGF-b gradient near tumor"),
]


def main() -> int:
    if not BINARY.exists():
        print(f"ERROR: binary not found at {BINARY}")
        return 1
    if not BASE_CONFIG.exists():
        print(f"ERROR: config not found at {BASE_CONFIG}")
        return 1

    work_dir = PROJECT_ROOT / "build" / "reality_check_1"
    if work_dir.exists():
        shutil.rmtree(work_dir)
    work_dir.mkdir(parents=True)

    print("=" * 72)
    print("  REALITY CHECK 1 -- Emergent Stroma Barrier Self-Assembly (HPC)")
    print(f"  Replicates: {NUM_REPLICATES}   Seeds: {SEEDS}")
    print(f"  Sim time: {SIM_MAX_TIME_MIN:.0f} min "
          f"({SIM_MAX_TIME_MIN/60/24:.1f} days)")
    print(f"  HPC: {SLURM_PARTITION} partition, {SLURM_CPUS} CPUs/node")
    print(f"  Quorum: >={QUORUM}/{NUM_REPLICATES}")
    print("=" * 72)

    # ---- Prepare and submit all replicates in parallel via SLURM ----
    job_ids: Dict[str, int] = {}   # job_id -> rep_idx
    rep_dirs: List[Path] = []

    for i, seed in enumerate(SEEDS[:NUM_REPLICATES]):
        rep_dir = work_dir / f"replicate_{i + 1:02d}_seed{seed}"
        rep_dir.mkdir(parents=True, exist_ok=True)
        output_dir = rep_dir / "output"
        output_dir.mkdir(parents=True, exist_ok=True)
        rep_dirs.append(rep_dir)

        config_path = rep_dir / "config.xml"
        _patch_config(
            BASE_CONFIG, config_path, output_dir, seed,
            SIM_MAX_TIME_MIN, SLURM_CPUS,
        )

        intervention = rep_dir / "intervention.json"
        intervention.write_text(json.dumps({
            "knob_interventions": [],
            "drug_delivery": {"enabled": False},
        }, indent=2), encoding="utf-8")

        script = _write_slurm_script(
            rep_dir, config_path, intervention, i, seed,
        )

        # Submit via sbatch
        result = subprocess.run(
            ["sbatch", "--parsable", str(script)],
            capture_output=True, text=True, check=False,
        )
        if result.returncode != 0:
            print(f"  ERROR: sbatch failed for rep {i+1}: "
                  f"{result.stderr.strip()}")
            return 1

        job_id = result.stdout.strip().split(";")[0]
        job_ids[job_id] = i
        print(f"  Submitted rep {i+1} (seed={seed}) -> SLURM job {job_id}")

    print(f"\n  All {NUM_REPLICATES} jobs submitted.  "
          f"Polling for completion ...\n")

    # ---- Poll until all jobs finish ----
    t0 = time.perf_counter()
    completed: Dict[int, Tuple[str, int]] = {}  # rep_idx -> (state, exit_code)

    while len(completed) < len(job_ids):
        time.sleep(SLURM_POLL_INTERVAL)
        elapsed = time.perf_counter() - t0

        for job_id, rep_idx in list(job_ids.items()):
            if rep_idx in completed:
                continue
            state, exit_code, terminal = _query_job(job_id)
            if terminal:
                completed[rep_idx] = (state, exit_code)
                seed = SEEDS[rep_idx]
                status = "OK" if (
                    state == "COMPLETED" and exit_code == 0
                ) else "XX"
                print(f"  {status} Rep {rep_idx+1} (seed={seed}): "
                      f"{state} exit={exit_code}  [{elapsed:.0f}s]")

        # Progress indicator for still-running jobs
        n_done = len(completed)
        still_running = NUM_REPLICATES - n_done
        if still_running > 0:
            for job_id, rep_idx in job_ids.items():
                if rep_idx in completed:
                    continue
                snap_count = len(
                    list(rep_dirs[rep_idx].glob("output/output*.xml"))
                )
                if snap_count > 0:
                    print(f"    rep {rep_idx+1}: {snap_count} snapshots "
                          f"so far ...  [{elapsed:.0f}s elapsed]")

    total_wall = time.perf_counter() - t0
    print(f"\n  All jobs completed in {total_wall:.1f}s "
          f"({total_wall/60:.1f} min)\n")

    # ---- Evaluate criteria ----
    results: List[ReplicateResult] = []
    t360_diag_printed = False
    for i, seed in enumerate(SEEDS[:NUM_REPLICATES]):
        rep_dir = rep_dirs[i]
        state, exit_code = completed.get(i, ("UNKNOWN", -1))
        success = state == "COMPLETED" and exit_code == 0

        # Check output files exist
        output_dir = rep_dir / "output"
        xmls = sorted(output_dir.glob("output*.xml"))
        if success and not xmls:
            success = False

        r = ReplicateResult(
            seed=seed,
            run_dir=rep_dir,
            success=success,
            wall_time_s=total_wall / NUM_REPLICATES,
        )
        if success:
            try:
                r = _evaluate(r)
                r.timeseries = _collect_time_series(output_dir)
                if not t360_diag_printed:
                    _print_t360_proliferation_diagnostic(
                        output_dir=output_dir,
                        config_path=rep_dir / "config.xml",
                        seed=seed,
                    )
                    t360_diag_printed = True
            except Exception as exc:
                print(f"  WARNING: evaluation failed for rep {i+1}: {exc}")
                r.success = False
        results.append(r)

    # ---- Print per-replicate results ----
    print()
    for i, r in enumerate(results):
        hdr = f"Replicate {i+1}  (seed={r.seed})"
        if not r.success:
            print(f"  {hdr}  -- SIMULATION FAILED")
            continue
        all_pass = all(r.criteria.get(k, False) for k, _ in CRITERIA_NAMES)
        tag = "ALL PASS" if all_pass else "SOME FAIL"
        print(f"  {hdr}  -- {tag}")
        for key, label in CRITERIA_NAMES:
            ok = r.criteria.get(key, False)
            detail = r.details.get(key, "")
            mark = "PASS" if ok else "FAIL"
            print(f"    [{mark}] {label}: {detail}")
        if r.timeseries:
            print("    Time series (t_min → snapshot):")
            for ts in r.timeseries:
                def _fmt(v):
                    return "nan" if v is None or (isinstance(v, float) and not math.isfinite(v)) else f"{v:.4f}" if isinstance(v, float) else str(v)
                print(
                    f"      t={ts['target_time']:.0f} (snap {ts['snapshot_time']:.0f}) | "
                    f"tumor={ts['total_tumor_cells']:5d} caf={ts['total_caf_count']:4d} | "
                    f"TGFb_psc_max={_fmt(ts['max_tgfb_at_psc'])} "
                    f"peri_ecm={_fmt(ts['mean_peritumoral_ecm'])} "
                    f"emt_frac={_fmt(ts['emt_fraction'])}"
                )
        print()

    # ---- Aggregate ----
    print("=" * 72)
    print(f"  AGGREGATE  (>={QUORUM}/{NUM_REPLICATES} replicates "
          f"must pass each criterion)")
    print("=" * 72)

    any_failure = False
    for key, label in CRITERIA_NAMES:
        passing = sum(
            1 for r in results
            if r.success and r.criteria.get(key, False)
        )
        ok = passing >= QUORUM
        mark = "PASS" if ok else "FAIL"
        if not ok:
            any_failure = True
        print(f"  [{mark}] {label}  ({passing}/{NUM_REPLICATES} replicates)")

    print()
    print(f"  Total wall time: {total_wall:.1f}s  ({total_wall/60:.1f} min)")
    if any_failure:
        print("\n  *** REALITY CHECK 1: FAIL ***")
        return 1
    else:
        print("\n  *** REALITY CHECK 1: PASS ***")
        return 0


if __name__ == "__main__":
    sys.exit(main())
