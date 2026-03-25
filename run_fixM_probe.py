#!/usr/bin/env python3
"""
run_fixM_probe.py

Runs a single-seed RC2 probe through day 31 and reports the withdrawal window
metrics required by fixM:
  - ABCB1+ fraction at days 28, 29, 30, 31
  - Intracellular drug at same timepoints
  - PASS/FAIL gate for post-withdrawal resistance persistence
"""
from __future__ import annotations

import argparse
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


T_PRE = 20160.0
T_TREAT_END = 40320.0
T_DAY_31 = 44640.0
CHECK_DAYS = [28, 29, 30, 31]


def _set(root, xpath: str, value: str) -> None:
    node = root.find(xpath)
    if node is not None:
        node.text = value


def _patch_config(src: Path, dst: Path, out_dir: Path, seed: int, threads: int, save_interval: int) -> None:
    tree = ET.parse(src)
    root = tree.getroot()

    _set(root, "./save/folder", str(out_dir))
    _set(root, ".//overall/max_time", str(T_DAY_31))
    _set(root, ".//options/random_seed", str(seed))
    _set(root, ".//parallel/omp_num_threads", str(threads))

    for n in root.findall(".//save//interval"):
        n.text = str(save_interval)

    _set(root, ".//user_parameters/drug_start_time", str(T_PRE))
    _set(root, ".//user_parameters/drug_end_time", str(T_TREAT_END))
    _set(root, ".//user_parameters/drug_concentration", "1.0")

    # Keep boundary conditions deterministic.
    for var_name, value in [
        ("oxygen", "38"),
        ("tgfb", "0"),
        ("shh", "0"),
        ("drug", "0"),
        ("ecm_density", "0"),
    ]:
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


def _row(matrix, labels, name: str):
    entry = labels.get(name)
    if entry is None:
        return None
    idx = int(entry["index"])
    if idx < 0 or idx >= matrix.shape[0]:
        return None
    return matrix[idx, :]


def _nearest_snapshot(parser: OutputParser, xmls: list[Path], target_time: float):
    best_xml = None
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
            "n_live": 0,
            "abcb1_frac": float("nan"),
            "ic_mean": float("nan"),
            "ic_max": float("nan"),
        }

    abcb1_frac = float(np.mean(abcb1[live_tumor] > 0.5)) if abcb1 is not None else float("nan")
    ic_mean = float(np.mean(ic_drug[live_tumor])) if ic_drug is not None else float("nan")
    ic_max = float(np.max(ic_drug[live_tumor])) if ic_drug is not None else float("nan")

    return {
        "n_live": float(n_live),
        "abcb1_frac": abcb1_frac,
        "ic_mean": ic_mean,
        "ic_max": ic_max,
    }


def _fmt(x: float) -> str:
    if np.isnan(x):
        return "nan"
    return f"{x:.4f}"


def main() -> int:
    ap = argparse.ArgumentParser(description="Run fixM day28-31 probe")
    ap.add_argument("--binary", type=Path, default=PROJECT_ROOT / "stroma_world")
    ap.add_argument("--base-config", type=Path, default=PROJECT_ROOT / "config" / "PhysiCell_settings.xml")
    ap.add_argument("--work-dir", type=Path, default=PROJECT_ROOT / "build" / "fixM_probe_seed42")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--save-interval", type=int, default=360)
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
    _patch_config(base_config, config_path, out_dir, args.seed, args.threads, args.save_interval)

    local_cfg = work_dir / "config"
    local_cfg.mkdir(exist_ok=True)
    for extra in ["tumor_calibration_knobs.json", "gene_params_default.json"]:
        src = PROJECT_ROOT / "config" / extra
        if src.exists():
            shutil.copy(src, local_cfg / extra)

    print("=" * 72)
    print("fixM probe: running single-seed day28-31 withdrawal check")
    print(f"  binary: {binary}")
    print(f"  output: {out_dir}")
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
    print("\nDay-window metrics (live tumor cells):")
    print("  day   t(min)   n_live   ABCB1+%   ic_mean   ic_max")
    print("  ----  -------  ------   -------   -------   -------")
    for day in CHECK_DAYS:
        t_target = day * 1440.0
        xml, t_actual = _nearest_snapshot(parser, xmls, t_target)
        m = _metrics_for_snapshot(parser, xml)
        rows[day] = m
        abcb1_pct = 100.0 * m["abcb1_frac"] if not np.isnan(m["abcb1_frac"]) else float("nan")
        print(
            f"  {day:>3d}  {t_actual:>7.0f}  {int(m['n_live']):>6d}   "
            f"{_fmt(abcb1_pct):>7}   {_fmt(m['ic_mean']):>7}   {_fmt(m['ic_max']):>7}"
        )

    # PASS criteria from directive
    abcb1_28 = rows[28]["abcb1_frac"]
    abcb1_29 = rows[29]["abcb1_frac"]
    abcb1_30 = rows[30]["abcb1_frac"]
    abcb1_31 = rows[31]["abcb1_frac"]

    ic_29 = rows[29]["ic_mean"]
    ic_30 = rows[30]["ic_mean"]
    ic_31 = rows[31]["ic_mean"]

    pass_abcb1_level = (abcb1_29 >= 0.50) and (abcb1_30 >= 0.50)
    pass_abcb1_gradual = (abcb1_29 >= 0.5 * abcb1_28) and (abcb1_30 >= 0.5 * abcb1_29) and (abcb1_31 >= 0.5 * abcb1_30)
    pass_ic_decay = (ic_30 <= ic_29) and (ic_31 <= ic_30)

    print("\nPASS checks:")
    print(f"  ABCB1 day29/day30 >= 50%:      {'PASS' if pass_abcb1_level else 'FAIL'}")
    print(f"  ABCB1 gradual decay (<=50/day): {'PASS' if pass_abcb1_gradual else 'FAIL'}")
    print(f"  Intracellular drug decreases:   {'PASS' if pass_ic_decay else 'FAIL'}")

    overall = pass_abcb1_level and pass_abcb1_gradual and pass_ic_decay
    print(f"\nfixM withdrawal probe overall: {'PASS' if overall else 'FAIL'}")
    return 0 if overall else 2


if __name__ == "__main__":
    sys.exit(main())
