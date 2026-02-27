#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import logging
import random
import sys
from dataclasses import asdict
from pathlib import Path
from typing import Any, Dict, List

# Ensure project-root imports work when running as a script.
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from python.ea.evolutionary_algorithm import EAConfig, EAResult, StromaWorldEA
from python.wrapper.output_parser import OutputParser
from python.wrapper.physicell_runner import PhysiCellRunner

try:
    import matplotlib.pyplot as plt
except ImportError:  # pragma: no cover
    plt = None


LOGGER = logging.getLogger("run_ea")

KNOB_ANALOG_LABELS = {
    "tgfb_secretion_rate": "tgfb-secretion modulator-like",
    "shh_secretion_rate": "shh-secretion modulator-like",
    "efflux_induction_delay": "efflux-kinetics modulator-like",
    "efflux_strength": "efflux-blocker-like",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the Stroma World evolutionary algorithm.")
    parser.add_argument("--config", required=True, help="Path to EA JSON config.")
    parser.add_argument("--physicell-binary", default=None, help="Path to PhysiCell executable.")
    parser.add_argument("--physicell-config", default=None, help="Path to PhysiCell XML config.")
    parser.add_argument("--output-dir", required=True, help="Output directory for EA artifacts.")
    parser.add_argument(
        "--parallel",
        type=int,
        default=None,
        help="Override number of parallel simulations to evaluate at once.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Override EA random seed (useful for array jobs).",
    )
    parser.add_argument("--resume", default=None, help="Checkpoint file to resume from.")
    parser.add_argument("--slurm", action="store_true", help="Use SLURM batch execution backend.")
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Generate and print 3 random intervention JSONs without running simulations.",
    )
    parser.add_argument(
        "--checkpoint-interval",
        type=int,
        default=1,
        help="Checkpoint interval in generations (default: 1).",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Logging verbosity.",
    )
    return parser.parse_args()


def configure_logging(level: str) -> None:
    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )


def load_ea_config(config_path: Path) -> EAConfig:
    with config_path.open("r", encoding="utf-8") as f:
        payload = json.load(f)
    if not isinstance(payload, dict):
        raise ValueError(f"EA config JSON must be an object: {config_path}")

    fields = set(EAConfig.__dataclass_fields__.keys())  # type: ignore[attr-defined]
    kwargs = {k: v for k, v in payload.items() if k in fields}
    unknown = sorted(set(payload.keys()) - fields)
    if unknown:
        LOGGER.warning("Ignoring unknown EA config keys: %s", ", ".join(unknown))
    return EAConfig(**kwargs)


def apply_cli_overrides(cfg: EAConfig, args: argparse.Namespace, output_dir: Path) -> EAConfig:
    cfg = EAConfig(**asdict(cfg))
    if args.physicell_binary:
        cfg.binary_path = str(Path(args.physicell_binary).expanduser().resolve())
    if args.physicell_config:
        cfg.config_path = str(Path(args.physicell_config).expanduser().resolve())

    cfg.simulation_output_dir = str((output_dir / "sim_runs").resolve())
    cfg.stats_csv_path = str((output_dir / "fitness_history.csv").resolve())
    if args.parallel is not None:
        if args.parallel < 1:
            raise ValueError("--parallel must be >= 1")
        cfg.parallel_sims = int(args.parallel)
    if args.seed is not None:
        cfg.random_seed = int(args.seed)
    if args.slurm:
        cfg.use_slurm = True
    return cfg


def generate_random_individual_payloads(cfg: EAConfig, n: int = 3) -> List[Dict[str, Any]]:
    rng = random.Random(cfg.random_seed)
    payloads: List[Dict[str, Any]] = []
    max_allowed = max(1, min(cfg.max_interventions, len(cfg.targetable_knobs)))

    for j in range(n):
        k = rng.randint(1, max_allowed)
        knobs = rng.sample(cfg.targetable_knobs, k=k)
        interventions = []
        for i, knob in enumerate(knobs):
            effect = rng.choices(["INHIBIT", "ACTIVATE"], weights=[0.8, 0.2], k=1)[0]
            strength = round(rng.uniform(0.1, 1.0), 6)
            interventions.append(
                {
                    "knob": knob,
                    "effect": effect,
                    "strength": strength,
                    "name": f"{effect}_{knob}_{i}",
                }
            )
        payloads.append(
            {
                "calibration_profile": cfg.calibration_profile,
                "knob_interventions": interventions,
                "candidate_index": j,
            }
        )
    return payloads


def save_best_intervention(result: EAResult, outpath: Path, calibration_profile: str) -> None:
    payload = {
        "best_fitness": result.best_fitness,
        "calibration_profile": calibration_profile,
        "knob_interventions": result.best_individual,
    }
    outpath.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def save_fitness_history_csv(result: EAResult, outpath: Path) -> None:
    fieldnames = ["generation", "best_fitness", "mean_fitness", "worst_fitness", "best_individual"]
    with outpath.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in result.fitness_history:
            writer.writerow({k: row.get(k) for k in fieldnames})


def save_ea_result_json(result: EAResult, outpath: Path) -> None:
    payload = {
        "best_individual": result.best_individual,
        "best_fitness": result.best_fitness,
        "fitness_history": result.fitness_history,
        "population_history": result.population_history,
        "runtime_seconds": result.runtime_seconds,
    }
    outpath.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def save_fitness_plot(result: EAResult, plots_dir: Path) -> None:
    plots_dir.mkdir(parents=True, exist_ok=True)
    if plt is None:
        LOGGER.warning("matplotlib is not installed; skipping fitness curve plot.")
        return
    if not result.fitness_history:
        LOGGER.warning("No fitness history available; skipping fitness curve plot.")
        return

    generations = [row.get("generation", 0) for row in result.fitness_history]
    best = [row.get("best_fitness", 0.0) for row in result.fitness_history]
    mean = [row.get("mean_fitness", 0.0) for row in result.fitness_history]
    worst = [row.get("worst_fitness", 0.0) for row in result.fitness_history]

    fig = plt.figure(figsize=(9, 5))
    ax = fig.add_subplot(111)
    ax.plot(generations, best, label="best", linewidth=2.0)
    ax.plot(generations, mean, label="mean", linewidth=1.8)
    ax.plot(generations, worst, label="worst", linewidth=1.5)
    ax.set_xlabel("Generation")
    ax.set_ylabel("Fitness")
    ax.set_title("EA Fitness Progress")
    ax.set_ylim(0.0, 1.0)
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(plots_dir / "fitness_curve.png", dpi=150)
    plt.close(fig)


def build_interpretation(best_individual: List[Dict[str, Any]]) -> str:
    knobs = {entry.get("knob", "") for entry in best_individual}
    fragments: List[str] = []

    if knobs & {"tgfb_secretion_rate", "shh_secretion_rate"}:
        fragments.append("reshape stromal signaling pressure")
    if "efflux_induction_delay" in knobs:
        fragments.append("slow resistance induction timing")
    if "efflux_strength" in knobs:
        fragments.append("limit efflux-mediated drug resistance")

    if not fragments:
        return "combine pathway perturbations to reduce barrier-supported tumor fitness."
    if len(fragments) == 1:
        return f"{fragments[0]}."
    if len(fragments) == 2:
        return f"{fragments[0]}, and {fragments[1]}."
    return ", ".join(fragments[:-1]) + f", and {fragments[-1]}."


def format_human_summary(result: EAResult) -> str:
    lines = [f"Best strategy (fitness={result.best_fitness:.2f}):"]
    for i, intervention in enumerate(result.best_individual, start=1):
        knob = intervention.get("knob", "UNKNOWN")
        effect = str(intervention.get("effect", "INHIBIT")).upper()
        strength = float(intervention.get("strength", 0.0))
        pct = int(round(strength * 100.0))
        verb = "Inhibit" if effect == "INHIBIT" else "Activate"
        alias = KNOB_ANALOG_LABELS.get(knob, "targeted-like")
        lines.append(f"{i}. {verb} {knob} at {pct}% strength ({alias})")

    interpretation = build_interpretation(result.best_individual)
    lines.append(f"Interpretation: {interpretation}")
    return "\n".join(lines)


def main() -> int:
    args = parse_args()
    configure_logging(args.log_level)

    config_path = Path(args.config).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    cfg = load_ea_config(config_path)
    cfg = apply_cli_overrides(cfg, args, output_dir)

    if args.dry_run:
        payloads = generate_random_individual_payloads(cfg, n=3)
        print("Dry-run: generated random intervention candidates")
        for payload in payloads:
            idx = payload.pop("candidate_index")
            print(f"\nCandidate {idx + 1}:")
            print(json.dumps(payload, indent=2))
        return 0

    # Explicitly initialize shared components requested by pipeline contract.
    runner = PhysiCellRunner(
        binary_path=cfg.binary_path,
        config_path=cfg.config_path,
        output_dir=cfg.simulation_output_dir,
        timeout_seconds=cfg.timeout_seconds,
    )
    _parser_probe = OutputParser(output_dir)
    _ = _parser_probe  # explicit init requested by pipeline contract.

    ea = StromaWorldEA(cfg)
    ea.runner = runner
    if args.slurm:
        ea.use_slurm = True
        LOGGER.info("SLURM backend enabled: evaluations use run_batch_slurm().")

    checkpoint_path = output_dir / "checkpoint_latest.pkl"
    start_generation = 0
    if args.resume:
        resume_path = Path(args.resume).expanduser().resolve()
        if not resume_path.exists():
            raise FileNotFoundError(f"Checkpoint not found: {resume_path}")
        start_generation = ea.load_checkpoint(resume_path)
        LOGGER.info("Resumed from checkpoint %s at generation %d", resume_path, start_generation)

    result = ea.run(
        start_generation=start_generation,
        checkpoint_path=checkpoint_path,
        checkpoint_interval=max(1, args.checkpoint_interval),
    )

    best_intervention_path = output_dir / "best_intervention.json"
    fitness_history_path = output_dir / "fitness_history.csv"
    ea_result_path = output_dir / "ea_result.json"
    plots_dir = output_dir / "plots"
    summary_path = output_dir / "summary.txt"

    save_best_intervention(result, best_intervention_path, cfg.calibration_profile)
    save_fitness_history_csv(result, fitness_history_path)
    save_ea_result_json(result, ea_result_path)
    save_fitness_plot(result, plots_dir)

    summary = format_human_summary(result)
    summary_path.write_text(summary + "\n", encoding="utf-8")
    print(summary)

    LOGGER.info("Saved best intervention to %s", best_intervention_path)
    LOGGER.info("Saved fitness history CSV to %s", fitness_history_path)
    LOGGER.info("Saved full EA result JSON to %s", ea_result_path)
    LOGGER.info("Saved summary to %s", summary_path)
    LOGGER.info("Latest checkpoint: %s", checkpoint_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
