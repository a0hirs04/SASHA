from __future__ import annotations

import logging
import math
import re
from collections import Counter
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np
import pandas as pd

try:
    import matplotlib.pyplot as plt
except ImportError:  # pragma: no cover
    plt = None  # type: ignore[assignment]

try:
    from ..ea.evolutionary_algorithm import EAResult
    from ..wrapper.output_parser import OutputParser
except Exception:  # pragma: no cover
    from python.ea.evolutionary_algorithm import EAResult  # type: ignore
    from python.wrapper.output_parser import OutputParser  # type: ignore


LOGGER = logging.getLogger(__name__)

GENE_NAME_PATTERN = re.compile(r"^[A-Z][A-Z0-9_]{1,24}$")


def _require_matplotlib() -> None:
    if plt is None:
        raise ImportError("matplotlib is required for plotting. Install with `pip install matplotlib`.")


def _apply_publication_style() -> None:
    _require_matplotlib()
    try:
        plt.style.use("seaborn-v0_8-whitegrid")
    except Exception:  # pragma: no cover
        plt.style.use("default")
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.titlesize": 12,
            "axes.labelsize": 10,
            "legend.fontsize": 9,
            "figure.dpi": 130,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "grid.alpha": 0.25,
            "lines.linewidth": 2.0,
        }
    )


def _output_base(output_path: Path) -> Path:
    path = Path(output_path).expanduser().resolve()
    if path.suffix.lower() in {".png", ".pdf"}:
        path = path.with_suffix("")
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def _save_dual_format(fig, output_path: Path) -> Tuple[Path, Path]:
    base = _output_base(output_path)
    png_path = base.with_suffix(".png")
    pdf_path = base.with_suffix(".pdf")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    return png_path, pdf_path


def _safe_float(value: Any, default: float = math.nan) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _get_fitness_std(df: pd.DataFrame) -> pd.Series:
    for col in ("std_fitness", "fitness_std", "stdev_fitness", "fitness_stdev"):
        if col in df.columns:
            return pd.to_numeric(df[col], errors="coerce").fillna(0.0)

    # Fallback estimate when explicit std is unavailable.
    if {"best_fitness", "mean_fitness", "worst_fitness"}.issubset(df.columns):
        best = pd.to_numeric(df["best_fitness"], errors="coerce")
        mean = pd.to_numeric(df["mean_fitness"], errors="coerce")
        worst = pd.to_numeric(df["worst_fitness"], errors="coerce")
        return pd.concat([(best - mean).abs(), (mean - worst).abs()], axis=1).max(axis=1).fillna(0.0)
    return pd.Series(np.zeros(len(df)))


def _snapshot_path_for_timestep(parser: OutputParser, timestep: int) -> Path:
    numbered = parser._numbered_output_xmls()  # noqa: SLF001 (intentional reuse)
    if numbered:
        by_step: Dict[int, Path] = {}
        for path in numbered:
            match = re.match(r"^output(\d+)\.xml$", path.name)
            if match:
                by_step[int(match.group(1))] = path

        if timestep in by_step:
            return by_step[timestep]
        nearest_step = min(by_step.keys(), key=lambda s: abs(s - timestep))
        LOGGER.warning("Requested timestep %d not found; using nearest available timestep %d.", timestep, nearest_step)
        return by_step[nearest_step]

    xmls = parser._list_timeseries_xmls()  # noqa: SLF001
    if not xmls:
        raise FileNotFoundError(f"No timestep XML files found in {parser.output_dir}")
    idx = int(np.clip(timestep, 0, len(xmls) - 1))
    return xmls[idx]


def _get_row(snapshot: Dict[str, Any], label_name: str) -> np.ndarray | None:
    label_map = snapshot.get("label_name_map", {})
    cell_matrix = snapshot.get("cell_matrix")
    if cell_matrix is None or not isinstance(label_map, dict):
        return None

    label = label_map.get(label_name)
    if label is None:
        return None

    idx = int(label.get("index", -1))
    if idx < 0 or idx >= cell_matrix.shape[0]:
        return None
    return cell_matrix[idx, :]


def _get_positions(snapshot: Dict[str, Any]) -> np.ndarray:
    label_map = snapshot.get("label_name_map", {})
    cell_matrix = snapshot.get("cell_matrix")
    if cell_matrix is None or not isinstance(label_map, dict):
        return np.empty((0, 3), dtype=float)

    pos_label = label_map.get("position")
    if pos_label is None:
        return np.empty((0, 3), dtype=float)

    start = int(pos_label.get("index", -1))
    size = int(pos_label.get("size", 0))
    if start < 0 or size < 2:
        return np.empty((0, 3), dtype=float)

    stop = min(start + min(size, 3), cell_matrix.shape[0])
    if stop - start < 2:
        return np.empty((0, 3), dtype=float)

    block = cell_matrix[start:stop, :].T
    if block.shape[1] == 2:
        block = np.column_stack([block, np.zeros(block.shape[0])])
    return block


def _extract_gene_matrix(snapshot: Dict[str, Any]) -> Tuple[List[str], np.ndarray]:
    label_map = snapshot.get("label_name_map", {})
    cell_matrix = snapshot.get("cell_matrix")
    if cell_matrix is None or not isinstance(label_map, dict):
        return [], np.empty((0, 0), dtype=float)

    rows: List[Tuple[int, str, np.ndarray]] = []
    for name, entry in label_map.items():
        if not GENE_NAME_PATTERN.match(str(name)):
            continue

        idx = int(entry.get("index", -1))
        size = int(entry.get("size", 1))
        if size != 1 or idx < 0 or idx >= cell_matrix.shape[0]:
            continue

        values = cell_matrix[idx, :]
        finite = values[np.isfinite(values)]
        if finite.size == 0:
            continue

        frac_unit_interval = np.mean((finite >= -1e-6) & (finite <= 1.0 + 1e-6))
        if frac_unit_interval < 0.95:
            continue
        rows.append((idx, str(name), values))

    rows.sort(key=lambda x: x[0])
    if not rows:
        return [], np.empty((0, 0), dtype=float)

    names = [r[1] for r in rows]
    data = np.vstack([r[2] for r in rows])
    return names, np.clip(data, 0.0, 1.0)


def plot_fitness_history(fitness_csv: Path, output_path: Path) -> None:
    _apply_publication_style()

    fitness_csv = Path(fitness_csv).expanduser().resolve()
    if not fitness_csv.exists():
        raise FileNotFoundError(f"Fitness CSV not found: {fitness_csv}")

    df = pd.read_csv(fitness_csv)
    required = {"best_fitness", "mean_fitness", "worst_fitness"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns in {fitness_csv}: {sorted(missing)}")

    gen = pd.to_numeric(df["generation"], errors="coerce") if "generation" in df.columns else pd.Series(np.arange(len(df)))
    best = pd.to_numeric(df["best_fitness"], errors="coerce")
    mean = pd.to_numeric(df["mean_fitness"], errors="coerce")
    worst = pd.to_numeric(df["worst_fitness"], errors="coerce")
    std = _get_fitness_std(df)

    fig, ax = plt.subplots(figsize=(9.5, 5.5))
    ax.plot(gen, best, color="#2166ac", label="Best")
    ax.plot(gen, mean, color="#1b9e77", label="Mean")
    ax.plot(gen, worst, color="#b2182b", label="Worst")
    lower = np.clip((mean - std).to_numpy(), 0.0, 1.0)
    upper = np.clip((mean + std).to_numpy(), 0.0, 1.0)
    ax.fill_between(gen.to_numpy(), lower, upper, color="#1b9e77", alpha=0.2, label="Std. dev.")
    ax.set_xlabel("Generation")
    ax.set_ylabel("Fitness")
    ax.set_title("EA Fitness History")
    ax.set_ylim(0.0, 1.0)
    ax.legend(loc="best")
    fig.tight_layout()

    _save_dual_format(fig, output_path)
    plt.close(fig)


def plot_intervention_frequency(population_history: list, output_path: Path) -> None:
    _apply_publication_style()

    if not population_history:
        raise ValueError("population_history is empty.")

    final_population = population_history[-1]
    if not isinstance(final_population, list):
        raise ValueError("Final population is not a list.")

    counts: Counter[str] = Counter()
    for individual in final_population:
        if not isinstance(individual, (list, tuple)):
            continue
        seen = set()
        for intervention in individual:
            if not isinstance(intervention, dict):
                continue
            knob = str(intervention.get("knob", "")).strip()
            if not knob:
                knob = str(intervention.get("gene", "")).strip()
            if not knob or knob in seen:
                continue
            seen.add(knob)
            counts[knob] += 1

    if not counts:
        raise ValueError("No valid interventions found in final population.")

    knobs, values = zip(*sorted(counts.items(), key=lambda kv: (-kv[1], kv[0])))
    x = np.arange(len(knobs))

    fig, ax = plt.subplots(figsize=(10, 5.5))
    bars = ax.bar(x, values, color="#4c78a8")
    ax.set_xticks(x)
    ax.set_xticklabels(knobs, rotation=30, ha="right")
    ax.set_ylabel("Individuals targeting knob")
    ax.set_title("Final Population Knob Frequency")
    ax.grid(axis="y", alpha=0.25)

    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.05, str(val), ha="center", va="bottom", fontsize=8)

    fig.tight_layout()
    _save_dual_format(fig, output_path)
    plt.close(fig)


def plot_gene_network_state(simulation_output: Path, timestep: int, output_path: Path) -> None:
    _apply_publication_style()

    parser = OutputParser(Path(simulation_output))
    xml_path = _snapshot_path_for_timestep(parser, timestep)
    snapshot = parser._read_physicell_xml(xml_path)  # noqa: SLF001

    gene_names, matrix = _extract_gene_matrix(snapshot)
    if matrix.size == 0:
        raise ValueError(f"No gene-like rows in [0,1] found for {xml_path}")

    n_genes, n_cells = matrix.shape
    max_cells = 500
    if n_cells > max_cells:
        sample_idx = np.linspace(0, n_cells - 1, max_cells, dtype=int)
        matrix = matrix[:, sample_idx]

    fig_height = max(5.5, 0.25 * n_genes + 2.5)
    fig, ax = plt.subplots(figsize=(12, fig_height))
    im = ax.imshow(matrix, aspect="auto", interpolation="nearest", cmap="magma", vmin=0.0, vmax=1.0)
    ax.set_yticks(np.arange(len(gene_names)))
    ax.set_yticklabels(gene_names)
    ax.set_xlabel(f"Cells ({matrix.shape[1]} shown)")
    ax.set_ylabel("Genes")
    ax.set_title(f"Gene Network State Heatmap (timestep={timestep})")
    cbar = fig.colorbar(im, ax=ax, fraction=0.015, pad=0.02)
    cbar.set_label("Expression level")
    fig.tight_layout()

    _save_dual_format(fig, output_path)
    plt.close(fig)


def plot_simulation_timeseries(timeseries_df: pd.DataFrame, output_path: Path) -> None:
    _apply_publication_style()

    if timeseries_df is None or timeseries_df.empty:
        raise ValueError("timeseries_df is empty.")

    df = timeseries_df.copy()
    for col in ["time", "tumor_count", "stroma_count", "mean_ecm", "drug_penetration", "activated_cafs"]:
        if col not in df.columns:
            raise ValueError(f"Missing required timeseries column: {col}")

    time = pd.to_numeric(df["time"], errors="coerce")
    tumor_count = pd.to_numeric(df["tumor_count"], errors="coerce")
    mean_ecm = pd.to_numeric(df["mean_ecm"], errors="coerce")
    drug_penetration = pd.to_numeric(df["drug_penetration"], errors="coerce")
    stroma_count = pd.to_numeric(df["stroma_count"], errors="coerce")
    activated_cafs = pd.to_numeric(df["activated_cafs"], errors="coerce")
    activated_fraction = (activated_cafs / stroma_count.replace(0, np.nan)).clip(0.0, 1.0)

    fig, axes = plt.subplots(4, 1, figsize=(10, 11), sharex=True)

    axes[0].plot(time, tumor_count, color="#2166ac")
    axes[0].set_ylabel("Tumor cells")
    axes[0].set_title("Tumor Cell Count")

    axes[1].plot(time, mean_ecm, color="#1b9e77")
    axes[1].set_ylabel("Mean ECM")
    axes[1].set_title("Mean ECM Density")

    axes[2].plot(time, drug_penetration, color="#f46d43")
    axes[2].set_ylabel("Drug")
    axes[2].set_title("Drug Penetration")

    axes[3].plot(time, activated_fraction, color="#762a83")
    axes[3].set_ylabel("CAF fraction")
    axes[3].set_xlabel("Time")
    axes[3].set_title("Activated CAF Fraction")
    axes[3].set_ylim(0.0, 1.0)

    for ax in axes:
        ax.grid(alpha=0.25)

    fig.tight_layout()
    _save_dual_format(fig, output_path)
    plt.close(fig)


def plot_stroma_barrier(simulation_output: Path, output_path: Path) -> None:
    _apply_publication_style()

    parser = OutputParser(Path(simulation_output))
    final_xml = parser._find_final_snapshot_xml()  # noqa: SLF001
    snapshot = parser._read_physicell_xml(final_xml)  # noqa: SLF001

    coords = np.asarray(snapshot.get("micro_coords", []), dtype=float)
    micro_values = snapshot.get("micro_values", {})
    ecm = None
    if isinstance(micro_values, dict):
        ecm = micro_values.get("ecm_density")
        if ecm is None:
            for key, values in micro_values.items():
                if "ecm" in str(key).lower():
                    ecm = values
                    break

    if coords.size == 0 or ecm is None or len(ecm) == 0:
        raise ValueError(f"ECM microenvironment data missing in {final_xml}")

    ecm = np.asarray(ecm, dtype=float)
    x = np.asarray(coords[:, 0], dtype=float)
    y = np.asarray(coords[:, 1], dtype=float)

    fig, ax = plt.subplots(figsize=(8.5, 7))

    xr = np.round(x, 6)
    yr = np.round(y, 6)
    xu = np.unique(xr)
    yu = np.unique(yr)

    if xu.size * yu.size == ecm.size:
        x_idx = {val: i for i, val in enumerate(xu)}
        y_idx = {val: i for i, val in enumerate(yu)}
        grid = np.full((yu.size, xu.size), np.nan, dtype=float)
        for xv, yv, ev in zip(xr, yr, ecm):
            grid[y_idx[yv], x_idx[xv]] = ev
        im = ax.imshow(
            grid,
            origin="lower",
            extent=[float(xu.min()), float(xu.max()), float(yu.min()), float(yu.max())],
            aspect="equal",
            cmap="viridis",
            vmin=0.0,
            vmax=1.0,
        )
    else:
        im = ax.tricontourf(x, y, ecm, levels=30, cmap="viridis", vmin=0.0, vmax=1.0)

    cell_types = _get_row(snapshot, "cell_type")
    positions = _get_positions(snapshot)
    if cell_types is not None and positions.size:
        mask = np.rint(cell_types).astype(int) == 0
        tumor_xy = positions[mask, :2]
        if tumor_xy.size:
            ax.scatter(
                tumor_xy[:, 0],
                tumor_xy[:, 1],
                s=8,
                c="white",
                edgecolors="black",
                linewidths=0.25,
                alpha=0.75,
                label="Tumor cells",
            )
            ax.legend(loc="upper right")

    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("ECM density")
    ax.set_xlabel("x (micron)")
    ax.set_ylabel("y (micron)")
    ax.set_title("Stroma Barrier: ECM Density with Tumor Overlay")
    fig.tight_layout()

    _save_dual_format(fig, output_path)
    plt.close(fig)


def compare_interventions(results: list[EAResult], labels: list[str], output_path: Path) -> None:
    _apply_publication_style()

    if not results:
        raise ValueError("results is empty.")
    if len(labels) != len(results):
        raise ValueError("labels length must match results length.")

    def as_mapping(result: Any) -> Dict[str, Any]:
        if isinstance(result, dict):
            return result
        return {
            "best_fitness": getattr(result, "best_fitness", math.nan),
            "best_individual": getattr(result, "best_individual", []),
            "fitness_history": getattr(result, "fitness_history", []),
            "runtime_seconds": getattr(result, "runtime_seconds", math.nan),
        }

    normalized = [as_mapping(r) for r in results]

    fig, axes = plt.subplots(2, 1, figsize=(11.5, 9), gridspec_kw={"height_ratios": [2.2, 1.3]})

    # Panel 1: best-fitness trajectory per run.
    ax0 = axes[0]
    plotted = False
    for label, result in zip(labels, normalized):
        history = result.get("fitness_history", []) or []
        if not history:
            continue
        gens = [_safe_float(h.get("generation"), i) for i, h in enumerate(history)]
        best = [_safe_float(h.get("best_fitness"), math.nan) for h in history]
        ax0.plot(gens, best, label=label)
        plotted = True
    ax0.set_title("Best Fitness Trajectories")
    ax0.set_ylabel("Best fitness")
    ax0.set_ylim(0.0, 1.0)
    ax0.grid(alpha=0.25)
    if plotted:
        ax0.legend(loc="best")

    # Panel 2: final scores and intervention size.
    ax1 = axes[1]
    x = np.arange(len(normalized))
    best_scores = np.array([_safe_float(r.get("best_fitness"), math.nan) for r in normalized], dtype=float)
    final_means = np.array(
        [
            _safe_float((r.get("fitness_history") or [{}])[-1].get("mean_fitness"), math.nan)
            if (r.get("fitness_history") or [])
            else math.nan
            for r in normalized
        ],
        dtype=float,
    )
    intervention_counts = np.array([len(r.get("best_individual", []) or []) for r in normalized], dtype=float)

    width = 0.38
    ax1.bar(x - width / 2, best_scores, width=width, color="#2166ac", label="Best fitness")
    ax1.bar(x + width / 2, final_means, width=width, color="#1b9e77", label="Final mean fitness")
    ax1.set_xticks(x)
    ax1.set_xticklabels(labels, rotation=20, ha="right")
    ax1.set_ylabel("Fitness")
    ax1.set_ylim(0.0, 1.0)
    ax1.grid(axis="y", alpha=0.25)

    ax1b = ax1.twinx()
    ax1b.plot(x, intervention_counts, "o--", color="#333333", label="#targets")
    ax1b.set_ylabel("Interventions in best strategy")

    handles_a, labels_a = ax1.get_legend_handles_labels()
    handles_b, labels_b = ax1b.get_legend_handles_labels()
    ax1.legend(handles_a + handles_b, labels_a + labels_b, loc="upper right")
    ax1.set_title("Final Run Comparison")

    fig.tight_layout()
    _save_dual_format(fig, output_path)
    plt.close(fig)
