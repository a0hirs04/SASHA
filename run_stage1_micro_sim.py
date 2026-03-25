#!/usr/bin/env python3
"""
Stage 1 micro-sim gate for resistance persistence (fast, compute-node friendly).

Protocol:
- 50-100 tumor cells (default 80), no stroma
- High-ECM microenvironment
- Drug ON from day 14 to day 28
- Drug OFF day 28 to day 31

Outputs at day 28/29/30/31:
- ABCB1+ fraction (%) among live tumor cells
- Mean intracellular drug
- Live tumor count (survival)

Exit code:
- 0 => Stage 1 PASS (no collapse)
- 2 => Stage 1 FAIL (collapse detected)
- 1 => runtime/config error
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(PROJECT_ROOT))
from python.wrapper.output_parser import OutputParser

DAY_MIN = 1440.0
T_WITHDRAW = 28.0 * DAY_MIN
T_END = 31.0 * DAY_MIN
T_DRUG_START = 14.0 * DAY_MIN
CHECK_DAYS = [28, 29, 30, 31]


def _set(root: ET.Element, xpath: str, value: str) -> None:
    node = root.find(xpath)
    if node is not None:
        node.text = value


def _set_domain(root: ET.Element, xmin: float, xmax: float, ymin: float, ymax: float, dx: float) -> None:
    _set(root, ".//domain/x_min", str(xmin))
    _set(root, ".//domain/x_max", str(xmax))
    _set(root, ".//domain/y_min", str(ymin))
    _set(root, ".//domain/y_max", str(ymax))
    _set(root, ".//domain/z_min", "-10")
    _set(root, ".//domain/z_max", "10")
    _set(root, ".//domain/dx", str(dx))
    _set(root, ".//domain/dy", str(dx))
    _set(root, ".//domain/dz", "20")
    _set(root, ".//domain/use_2D", "true")


def _patch_config(
    src: Path,
    dst: Path,
    out_dir: Path,
    seed: int,
    threads: int,
    save_interval: int,
    drug_concentration: float,
) -> None:
    tree = ET.parse(src)
    root = tree.getroot()

    # Fast micro-sim geometry/runtime.
    _set_domain(root, -200.0, 200.0, -200.0, 200.0, 20.0)
    _set(root, "./save/folder", str(out_dir))
    _set(root, ".//overall/max_time", str(T_END))
    _set(root, ".//overall/dt_diffusion", "0.1")
    _set(root, ".//overall/dt_mechanics", "0.5")
    _set(root, ".//overall/dt_phenotype", "6")
    _set(root, ".//options/random_seed", str(seed))
    _set(root, ".//parallel/omp_num_threads", str(threads))

    # Speed: disable SVG, keep full snapshots for parser metrics.
    _set(root, ".//save/SVG/enable", "false")
    _set(root, ".//save/legacy_data/enable", "false")
    for n in root.findall(".//save//interval"):
        n.text = str(save_interval)

    # Stage 1 protocol: short pre-phase, then treatment to day 28, then withdrawal.
    _set(root, ".//user_parameters/drug_start_time", str(T_DRUG_START))
    _set(root, ".//user_parameters/drug_end_time", str(T_WITHDRAW))
    _set(root, ".//user_parameters/drug_concentration", str(drug_concentration))

    # 50-100 cells only.
    _set(root, ".//user_parameters/number_of_tumor_cells", "60")
    _set(root, ".//user_parameters/number_of_stromal_cells", "20")
    _set(root, ".//user_parameters/tumor_cluster_radius", "60.0")
    _set(root, ".//user_parameters/stroma_inner_radius", "80.0")
    _set(root, ".//user_parameters/stroma_outer_radius", "120.0")

    # Deterministic boundaries + high-ECM region.
    for var_name, value in [
        ("oxygen", "38"),
        ("tgfb", "0"),
        ("shh", "0"),
        ("drug", "0"),
        ("ecm_density", "0.8"),
    ]:
        var = root.find(f".//microenvironment_setup/variable[@name='{var_name}']")
        if var is None:
            continue
        init = var.find("./initial_condition")
        if init is not None:
            init.text = value
        dbc = var.find("./Dirichlet_boundary_condition")
        if dbc is not None:
            dbc.text = value
            dbc.set("enabled", "true")
        for bv in var.findall("./Dirichlet_options/boundary_value"):
            bv.text = value
            bv.set("enabled", "true")

    tree.write(dst, encoding="utf-8", xml_declaration=True)


def _row(matrix: np.ndarray, labels: dict, name: str):
    entry = labels.get(name)
    if entry is None:
        return None
    idx = int(entry["index"])
    if idx < 0 or idx >= matrix.shape[0]:
        return None
    return matrix[idx, :]


def _nearest_snapshot(parser: OutputParser, xmls: list[Path], target_time: float) -> tuple[Path, float]:
    best_xml = xmls[0]
    best_t = 0.0
    best_dt = float("inf")
    for xml in xmls:
        snap = parser._read_physicell_xml(xml)
        t = float(snap["time"])
        dt = abs(t - target_time)
        if dt < best_dt:
            best_dt = dt
            best_xml = xml
            best_t = t
    return best_xml, best_t


def _metrics_for_snapshot(parser: OutputParser, xml_path: Path) -> dict[str, float]:
    snap = parser._read_physicell_xml(xml_path)
    matrix = snap["cell_matrix"]
    labels = snap["label_name_map"]

    cell_type = _row(matrix, labels, "cell_type")
    dead = _row(matrix, labels, "dead")
    death_model = _row(matrix, labels, "current_death_model")
    abcb1 = _row(matrix, labels, "abcb1_active")
    ic_drug = _row(matrix, labels, "intracellular_drug")

    n_cells = matrix.shape[1]
    ctype = np.rint(cell_type).astype(int) if cell_type is not None else np.full(n_cells, -1)
    tumor_mask = ctype == 0
    dead_mask = (dead > 0.5) if dead is not None else np.zeros(n_cells, dtype=bool)
    apoptotic_mask = (np.rint(death_model).astype(int) == 100) if death_model is not None else np.zeros(n_cells, dtype=bool)
    live_tumor = tumor_mask & ~dead_mask & ~apoptotic_mask

    n_live = int(np.sum(live_tumor))
    if n_live == 0:
        return {
            "n_live": 0.0,
            "abcb1_frac": float("nan"),
            "ic_mean": float("nan"),
        }

    return {
        "n_live": float(n_live),
        "abcb1_frac": float(np.mean(abcb1[live_tumor] > 0.5)) if abcb1 is not None else float("nan"),
        "ic_mean": float(np.mean(ic_drug[live_tumor])) if ic_drug is not None else float("nan"),
    }


def _fmt(x: float) -> str:
    if np.isnan(x):
        return "nan"
    return f"{x:.4f}"


def main() -> int:
    ap = argparse.ArgumentParser(description="Stage 1 micro-sim resistance persistence gate")
    ap.add_argument("--binary", type=Path, default=PROJECT_ROOT / "stroma_world")
    ap.add_argument("--base-config", type=Path, default=PROJECT_ROOT / "config" / "PhysiCell_settings.xml")
    ap.add_argument("--work-dir", type=Path, default=PROJECT_ROOT / "build" / "stage1_micro")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--threads", type=int, default=int(os.environ.get("OMP_NUM_THREADS", "16")))
    ap.add_argument("--save-interval", type=int, default=360)
    ap.add_argument("--drug-concentration", type=float, default=0.10)
    args = ap.parse_args()

    binary = args.binary.resolve()
    base_config = args.base_config.resolve()
    work_dir = args.work_dir.resolve()
    out_dir = work_dir / "output"

    if not binary.exists():
        print(f"ERROR: binary not found: {binary}")
        return 1
    if not base_config.exists():
        print(f"ERROR: base config not found: {base_config}")
        return 1

    if work_dir.exists():
        shutil.rmtree(work_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    config_path = work_dir / "config.xml"
    _patch_config(
        base_config,
        config_path,
        out_dir,
        args.seed,
        args.threads,
        args.save_interval,
        args.drug_concentration,
    )

    local_cfg = work_dir / "config"
    local_cfg.mkdir(exist_ok=True)
    for extra in ["tumor_calibration_knobs.json", "gene_params_default.json"]:
        src = PROJECT_ROOT / "config" / extra
        if src.exists():
            shutil.copy(src, local_cfg / extra)

    print("=" * 72)
    print("Stage 1 micro-sim: resistance persistence gate")
    print(f"  binary: {binary}")
    print(f"  output: {out_dir}")
    print(f"  seed:   {args.seed}")
    print(f"  drug:   {args.drug_concentration}")
    print("=" * 72)

    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = str(args.threads)
    result = subprocess.run([str(binary), str(config_path)], cwd=PROJECT_ROOT, env=env)
    if result.returncode != 0:
        print(f"ERROR: simulation exited with code {result.returncode}")
        return result.returncode

    parser = OutputParser(out_dir)
    xmls = sorted(out_dir.glob("output*.xml"))
    if not xmls:
        print("ERROR: no snapshots found")
        return 1

    rows: dict[int, dict[str, float]] = {}
    print("\nStage 1 metrics (live tumor cells):")
    print("  day   t(min)   n_live   ABCB1+%   ic_mean")
    print("  ----  -------  ------   -------   -------")
    for day in CHECK_DAYS:
        t_target = day * DAY_MIN
        xml, t_actual = _nearest_snapshot(parser, xmls, t_target)
        m = _metrics_for_snapshot(parser, xml)
        rows[day] = m
        abcb1_pct = 100.0 * m["abcb1_frac"] if not np.isnan(m["abcb1_frac"]) else float("nan")
        print(
            f"  {day:>3d}  {t_actual:>7.0f}  {int(m['n_live']):>6d}   "
            f"{_fmt(abcb1_pct):>7}   {_fmt(m['ic_mean']):>7}"
        )

    abcb1_28 = rows[28]["abcb1_frac"]
    abcb1_29 = rows[29]["abcb1_frac"]
    abcb1_30 = rows[30]["abcb1_frac"]
    abcb1_31 = rows[31]["abcb1_frac"]

    n28 = rows[28]["n_live"]
    n31 = rows[31]["n_live"]

    # Collapse detector: if post-withdrawal day31 drops below half of day28,
    # resistance memory failed.
    pass_abcb1_memory = (
        not np.isnan(abcb1_28)
        and not np.isnan(abcb1_31)
        and abcb1_31 >= 0.5 * abcb1_28
        and abcb1_29 >= 0.5 * abcb1_28
        and abcb1_30 >= 0.5 * abcb1_28
    )
    pass_survival = (not np.isnan(n28)) and (not np.isnan(n31)) and (n28 > 0.0) and (n31 >= 0.7 * n28)

    print("\nStage 1 gates:")
    print(f"  ABCB1 memory (d29-31 >= 50% of d28): {'PASS' if pass_abcb1_memory else 'FAIL'}")
    print(f"  Survival retention (d31 >= 70% of d28): {'PASS' if pass_survival else 'FAIL'}")

    overall = pass_abcb1_memory and pass_survival
    print(f"\nStage 1 overall: {'PASS' if overall else 'FAIL'}")

    summary = {
        "overall_pass": bool(overall),
        "rows": {
            str(day): {
                "n_live": float(rows[day]["n_live"]),
                "abcb1_frac": float(rows[day]["abcb1_frac"]),
                "ic_mean": float(rows[day]["ic_mean"]),
            }
            for day in CHECK_DAYS
        },
        "gates": {
            "abcb1_memory": bool(pass_abcb1_memory),
            "survival_retention": bool(pass_survival),
        },
    }
    (work_dir / "stage1_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    return 0 if overall else 2


if __name__ == "__main__":
    sys.exit(main())
