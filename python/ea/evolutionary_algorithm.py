from __future__ import annotations

import concurrent.futures
import copy
import csv
import json
import logging
import math
import pickle
import random
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

try:
    from deap import base, creator, tools  # type: ignore[import-not-found]
except ImportError:  # pragma: no cover
    base = None  # type: ignore[assignment]
    creator = None  # type: ignore[assignment]
    tools = None  # type: ignore[assignment]

try:
    from ..wrapper.output_parser import OutputParser, SimulationMetrics
    from ..wrapper.physicell_runner import PhysiCellRunner
    from .fitness import (
        FitnessConfig,
        compute_fitness,
        compute_fitness_detailed,
        set_fitness_config,
        validate_metrics,
    )
except Exception:  # pragma: no cover
    from python.wrapper.output_parser import OutputParser, SimulationMetrics  # type: ignore
    from python.wrapper.physicell_runner import PhysiCellRunner  # type: ignore
    from python.ea.fitness import (  # type: ignore
        FitnessConfig,
        compute_fitness,
        compute_fitness_detailed,
        set_fitness_config,
        validate_metrics,
    )
try:
    from .knob_schema import TARGETABLE_KNOBS as SCHEMA_TARGETABLE_KNOBS, validate_partition as validate_knob_partition
except Exception:  # pragma: no cover
    from python.ea.knob_schema import TARGETABLE_KNOBS as SCHEMA_TARGETABLE_KNOBS, validate_partition as validate_knob_partition  # type: ignore


LOGGER = logging.getLogger(__name__)

DEFAULT_TARGETABLE_KNOBS = list(SCHEMA_TARGETABLE_KNOBS)

ALLOWED_EFFECTS = ("INHIBIT", "ACTIVATE")


@dataclass
class EAConfig:
    population_size: int = 50
    generations: int = 100
    mutation_rate: float = 0.2
    crossover_rate: float = 0.7
    tournament_size: int = 3
    max_interventions: int = 4
    parallel_sims: int = 4
    targetable_knobs: List[str] = field(default_factory=lambda: list(DEFAULT_TARGETABLE_KNOBS))
    calibration_profile: str = "AsPC-1"

    # Runtime / integration settings.
    binary_path: str = "/work/a0hirs04/PhysiCell/stroma_world"
    config_path: str = "/work/a0hirs04/PhysiCell/config/PhysiCell_settings.xml"
    simulation_output_dir: str = "/work/a0hirs04/PhysiCell/ea_runs"
    timeout_seconds: int = 3600
    random_seed: Optional[int] = None
    keep_population_history: bool = True
    use_slurm: bool = False
    slurm_poll_interval_seconds: int = 15
    slurm_sbatch_args: List[str] = field(default_factory=list)

    # EA bookkeeping.
    elitism_fraction: float = 0.05
    stats_csv_path: str = "/home/a0hirs04/PROJECT-NORTHSTAR/python/ea/ea_generation_stats.csv"

    # Fitness normalization constants (mirrors FitnessConfig).
    initial_tumor_cells: int = 50
    max_expected_drug: float = 1.0
    domain_size: float = 2000.0


@dataclass
class EAResult:
    best_individual: List[Dict[str, Any]]
    best_fitness: float
    fitness_history: List[Dict[str, Any]]
    population_history: List[List[List[Dict[str, Any]]]]
    runtime_seconds: float


class StromaWorldEA:
    def __init__(self, config: EAConfig):
        if base is None or creator is None or tools is None:
            raise ImportError("DEAP is required. Install with `pip install deap`.")

        self.config = config
        validate_knob_partition()
        self._rng = random.Random(config.random_seed)
        if config.random_seed is not None:
            random.seed(config.random_seed)

        self.runner = PhysiCellRunner(
            binary_path=config.binary_path,
            config_path=config.config_path,
            output_dir=config.simulation_output_dir,
            timeout_seconds=config.timeout_seconds,
        )
        self.use_slurm = bool(config.use_slurm)
        self.slurm_poll_interval_seconds = int(config.slurm_poll_interval_seconds)
        self.slurm_sbatch_args = list(config.slurm_sbatch_args)

        self._fitness_config = FitnessConfig(
            initial_tumor_cells=config.initial_tumor_cells,
            max_expected_drug=config.max_expected_drug,
            domain_size=config.domain_size,
        )
        set_fitness_config(self._fitness_config)

        self.toolbox = base.Toolbox()
        self._setup_deap_toolbox()

        self.population: List[Any] = []
        self.fitness_history: List[Dict[str, Any]] = []
        self.population_history: List[List[List[Dict[str, Any]]]] = []

        self.stats_csv_path = Path(config.stats_csv_path).expanduser().resolve()
        self.stats_csv_path.parent.mkdir(parents=True, exist_ok=True)
        self._init_stats_csv()

    def initialize_population(self) -> list:
        return self.toolbox.population(n=self.config.population_size)

    def evaluate(self, individual) -> float:
        """
        Evaluate one intervention strategy by full in-silico simulation:
        1) serialize intervention policy
        2) run PhysiCell
        3) parse final state
        4) compute Cycle 1 fitness
        """
        intervention_payload = self._individual_to_json(individual)
        sim_result = self.runner.run(intervention_payload)
        return self._score_sim_result(sim_result)

    def mutate(self, individual):
        """
        Mutation strategy:
        - For each intervention, with probability mutation_rate:
          40% mutate strength
          20% mutate target knob
          10% mutate effect type
          15% add intervention
          15% remove intervention
        """
        if len(individual) == 0:
            individual.append(self._random_intervention())

        pending_add = 0
        pending_remove = 0

        # Per-intervention mutations.
        for entry in list(individual):
            if self._rng.random() >= self.config.mutation_rate:
                continue
            roll = self._rng.random()
            if roll < 0.40:
                delta = self._rng.gauss(0.0, 0.1)
                entry["strength"] = self._clamp(entry["strength"] + delta, 0.1, 1.0)
            elif roll < 0.60:
                current_knobs = {e["knob"] for e in individual if e is not entry}
                entry["knob"] = self._random_knob(exclude=current_knobs)
            elif roll < 0.70:
                entry["effect"] = self._random_effect(exclude=entry.get("effect"))
            elif roll < 0.85:
                pending_add += 1
            else:
                pending_remove += 1

        for _ in range(pending_add):
            if len(individual) >= self.config.max_interventions:
                break
            current_knobs = {e["knob"] for e in individual}
            individual.append(self._random_intervention(exclude_knobs=current_knobs))

        for _ in range(pending_remove):
            if len(individual) <= 1:
                break
            drop_idx = self._rng.randrange(len(individual))
            individual.pop(drop_idx)

        normalized = self._normalize_interventions(individual)
        individual[:] = normalized
        return (individual,)

    def crossover(self, ind1, ind2):
        """
        Uniform crossover at intervention level:
        - pool interventions from both parents
        - each child samples uniformly from pooled interventions
        - no duplicate knob targets in resulting children
        """
        pool = [copy.deepcopy(x) for x in ind1] + [copy.deepcopy(x) for x in ind2]
        if not pool:
            return ind1, ind2

        child1 = self._sample_child_from_pool(pool)
        child2 = self._sample_child_from_pool(pool)

        ind1[:] = child1
        ind2[:] = child2
        return ind1, ind2

    def run(
        self,
        start_generation: int = 0,
        checkpoint_path: Optional[Path | str] = None,
        checkpoint_interval: int = 1,
    ) -> EAResult:
        start = time.perf_counter()

        if not self.population:
            self.population = self.initialize_population()

        self._evaluate_population(self.population)
        if start_generation <= 0:
            self._record_generation_stats(generation=0, population=self.population)

        elite_count = max(1, int(round(self.config.population_size * self.config.elitism_fraction)))
        elite_count = min(elite_count, self.config.population_size)

        first_generation = max(1, int(start_generation) + 1)
        for generation in range(first_generation, self.config.generations + 1):
            parents = self.toolbox.select(self.population, len(self.population))
            offspring = list(map(self.toolbox.clone, parents))

            for i in range(1, len(offspring), 2):
                if self._rng.random() < self.config.crossover_rate:
                    self.toolbox.mate(offspring[i - 1], offspring[i])
                    if hasattr(offspring[i - 1].fitness, "values"):
                        del offspring[i - 1].fitness.values
                    if hasattr(offspring[i].fitness, "values"):
                        del offspring[i].fitness.values

            for individual in offspring:
                self.toolbox.mutate(individual)
                if hasattr(individual.fitness, "values"):
                    del individual.fitness.values

            self._evaluate_population(offspring)

            elites = list(map(self.toolbox.clone, tools.selBest(self.population, elite_count)))
            non_elite_count = self.config.population_size - elite_count
            survivors = tools.selBest(offspring, non_elite_count)
            self.population = elites + survivors

            self._record_generation_stats(generation=generation, population=self.population)

            if checkpoint_path is not None and checkpoint_interval > 0:
                if generation % checkpoint_interval == 0:
                    self.save_checkpoint(generation=generation, filepath=checkpoint_path)

        best = tools.selBest(self.population, 1)[0]
        runtime = time.perf_counter() - start
        if checkpoint_path is not None:
            self.save_checkpoint(generation=self.config.generations, filepath=checkpoint_path)
        return EAResult(
            best_individual=self._plain_interventions(best),
            best_fitness=float(best.fitness.values[0]),
            fitness_history=list(self.fitness_history),
            population_history=list(self.population_history),
            runtime_seconds=runtime,
        )

    def save_checkpoint(self, generation, filepath):
        path = Path(filepath).expanduser().resolve()
        path.parent.mkdir(parents=True, exist_ok=True)

        payload = {
            "generation": int(generation),
            "config": asdict(self.config),
            "population": [self._plain_interventions(ind) for ind in self.population],
            "population_fitness": [
                float(ind.fitness.values[0]) if ind.fitness.valid else None for ind in self.population
            ],
            "fitness_history": self.fitness_history,
            "population_history": self.population_history,
            "random_state": random.getstate(),
        }
        with path.open("wb") as f:
            pickle.dump(payload, f)

    def load_checkpoint(self, filepath):
        path = Path(filepath).expanduser().resolve()
        with path.open("rb") as f:
            payload = pickle.load(f)

        saved_population = payload.get("population", [])
        saved_fitness = payload.get("population_fitness", [])

        self.population = []
        for idx, interventions in enumerate(saved_population):
            ind = self._individual_from_plain(interventions)
            fit = saved_fitness[idx] if idx < len(saved_fitness) else None
            if fit is not None:
                ind.fitness.values = (float(fit),)
            self.population.append(ind)

        self.fitness_history = list(payload.get("fitness_history", []))
        self.population_history = list(payload.get("population_history", []))

        random_state = payload.get("random_state")
        if random_state is not None:
            random.setstate(random_state)

        return int(payload.get("generation", 0))

    def generate_summary_report(self, result: EAResult, output_path: Optional[Path | str] = None) -> str:
        lines = [
            "StromaWorld EA Summary",
            "======================",
            f"Best fitness: {result.best_fitness:.6f}",
            f"Runtime (s): {result.runtime_seconds:.2f}",
            f"Intervention count: {len(result.best_individual)}",
            "",
            "Best intervention strategy:",
        ]
        for i, intervention in enumerate(result.best_individual, start=1):
            lines.append(
                f"{i}. knob={intervention['knob']} effect={intervention['effect']} strength={intervention['strength']:.3f}"
            )

        if result.fitness_history:
            last = result.fitness_history[-1]
            lines.extend(
                [
                    "",
                    "Final generation stats:",
                    f"- generation: {last.get('generation')}",
                    f"- best_fitness: {last.get('best_fitness'):.6f}",
                    f"- mean_fitness: {last.get('mean_fitness'):.6f}",
                    f"- worst_fitness: {last.get('worst_fitness'):.6f}",
                ]
            )

        report = "\n".join(lines) + "\n"
        if output_path is not None:
            path = Path(output_path).expanduser().resolve()
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(report, encoding="utf-8")
        return report

    # ---- Internal helpers -------------------------------------------------

    def _setup_deap_toolbox(self) -> None:
        if not hasattr(creator, "StromaFitnessMax"):
            creator.create("StromaFitnessMax", base.Fitness, weights=(1.0,))
        if not hasattr(creator, "StromaInterventionIndividual"):
            creator.create("StromaInterventionIndividual", list, fitness=creator.StromaFitnessMax)

        self.toolbox.register("individual", self._make_random_individual)
        self.toolbox.register("population", tools.initRepeat, list, self.toolbox.individual)
        self.toolbox.register("evaluate", self._evaluate_for_deap)
        self.toolbox.register("mate", self.crossover)
        self.toolbox.register("mutate", self.mutate)
        self.toolbox.register("select", tools.selTournament, tournsize=self.config.tournament_size)
        self.toolbox.register("clone", copy.deepcopy)

    def _make_random_individual(self):
        max_allowed = max(1, min(self.config.max_interventions, len(self.config.targetable_knobs)))
        n_interventions = self._rng.randint(1, max_allowed)
        knobs = self._rng.sample(self.config.targetable_knobs, k=n_interventions)
        interventions = [self._random_intervention(force_knob=k) for k in knobs]
        interventions = self._normalize_interventions(interventions)
        return creator.StromaInterventionIndividual(interventions)

    def _evaluate_for_deap(self, individual):
        return (self.evaluate(individual),)

    def _evaluate_population(self, population: Sequence[Any]) -> None:
        invalid = [ind for ind in population if not ind.fitness.valid]
        if not invalid:
            return

        if self.use_slurm:
            self._evaluate_population_slurm(invalid)
            return

        with concurrent.futures.ThreadPoolExecutor(max_workers=self.config.parallel_sims) as pool:
            futures = {pool.submit(self.evaluate, ind): ind for ind in invalid}
            for future in concurrent.futures.as_completed(futures):
                ind = futures[future]
                try:
                    fit = float(future.result())
                except Exception:
                    LOGGER.exception("EA evaluation failed for individual; assigning fitness=0.")
                    fit = 0.0
                ind.fitness.values = (fit,)

    def _evaluate_population_slurm(self, invalid: Sequence[Any]) -> None:
        payloads = [self._individual_to_json(ind) for ind in invalid]
        try:
            sim_results = self.runner.run_batch_slurm(
                payloads,
                max_parallel=self.config.parallel_sims,
                poll_interval_seconds=self.slurm_poll_interval_seconds,
                sbatch_args=self.slurm_sbatch_args,
            )
        except Exception:
            LOGGER.exception("SLURM batch evaluation failed; assigning fitness=0 for batch.")
            for ind in invalid:
                ind.fitness.values = (0.0,)
            return

        for ind, sim_result in zip(invalid, sim_results):
            fit = self._score_sim_result(sim_result)
            ind.fitness.values = (fit,)

        if len(sim_results) < len(invalid):
            for ind in invalid[len(sim_results) :]:
                ind.fitness.values = (0.0,)

    def _score_sim_result(self, sim_result) -> float:
        if not sim_result.success:
            return 0.0

        try:
            parser = OutputParser(sim_result.run_dir / "output")
            metrics = parser.parse_final_state()
        except Exception:
            LOGGER.exception("Failed to parse run output from %s", sim_result.run_dir)
            return 0.0

        if not validate_metrics(metrics):
            return 0.0
        return float(compute_fitness(metrics))

    def _record_generation_stats(self, generation: int, population: Sequence[Any]) -> None:
        fitness_values = [float(ind.fitness.values[0]) for ind in population]
        best = max(fitness_values) if fitness_values else 0.0
        mean = sum(fitness_values) / len(fitness_values) if fitness_values else 0.0
        worst = min(fitness_values) if fitness_values else 0.0
        best_ind = tools.selBest(population, 1)[0] if population else None

        row = {
            "generation": generation,
            "best_fitness": best,
            "mean_fitness": mean,
            "worst_fitness": worst,
            "best_individual": json.dumps(self._plain_interventions(best_ind) if best_ind is not None else []),
        }
        self.fitness_history.append(row)

        if self.config.keep_population_history:
            snapshot = [self._plain_interventions(ind) for ind in population]
            self.population_history.append(snapshot)

        self._append_stats_csv(row)

    def _init_stats_csv(self) -> None:
        with self.stats_csv_path.open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(
                f,
                fieldnames=["generation", "best_fitness", "mean_fitness", "worst_fitness", "best_individual"],
            )
            writer.writeheader()

    def _append_stats_csv(self, row: Dict[str, Any]) -> None:
        with self.stats_csv_path.open("a", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(
                f,
                fieldnames=["generation", "best_fitness", "mean_fitness", "worst_fitness", "best_individual"],
            )
            writer.writerow(row)

    def _sample_child_from_pool(self, pool: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        self._rng.shuffle(pool)
        child: List[Dict[str, Any]] = []
        used_knobs = set()
        for intervention in pool:
            if self._rng.random() < 0.5 and intervention["knob"] not in used_knobs:
                child.append(copy.deepcopy(intervention))
                used_knobs.add(intervention["knob"])
            if len(child) >= self.config.max_interventions:
                break

        if not child:
            pick = copy.deepcopy(self._rng.choice(pool))
            child = [pick]

        return self._normalize_interventions(child)

    def _random_intervention(
        self,
        exclude_knobs: Optional[set] = None,
        force_knob: Optional[str] = None,
    ) -> Dict[str, Any]:
        knob = force_knob if force_knob is not None else self._random_knob(exclude=exclude_knobs)
        effect = self._random_effect()
        strength = self._rng.uniform(0.1, 1.0)
        return {
            "knob": knob,
            "effect": effect,
            "strength": round(float(strength), 6),
        }

    def _random_knob(self, exclude: Optional[set] = None) -> str:
        exclude = exclude or set()
        candidates = [k for k in self.config.targetable_knobs if k not in exclude]
        if not candidates:
            candidates = list(self.config.targetable_knobs)
        return self._rng.choice(candidates)

    def _random_effect(self, exclude: Optional[str] = None) -> str:
        effects = list(ALLOWED_EFFECTS)
        weights = [0.8, 0.2]  # INHIBIT / ACTIVATE
        if exclude in effects and len(effects) > 1:
            idx = effects.index(exclude)
            effects.pop(idx)
            weights.pop(idx)
        total = sum(weights)
        probs = [w / total for w in weights]
        return self._rng.choices(effects, weights=probs, k=1)[0]

    def _normalize_interventions(self, interventions: Sequence[Dict[str, Any]]) -> List[Dict[str, Any]]:
        seen = set()
        normalized: List[Dict[str, Any]] = []
        for entry in interventions:
            knob = str(entry.get("knob", "")).strip()
            effect = str(entry.get("effect", "INHIBIT")).strip().upper()
            strength = float(entry.get("strength", 0.5))

            if knob not in self.config.targetable_knobs:
                raise ValueError(
                    f"Partition violation: knob '{knob}' is not targetable. "
                    "EA may only use tgfb_secretion_rate, shh_secretion_rate, "
                    "efflux_induction_delay, efflux_strength."
                )
            if effect not in ALLOWED_EFFECTS:
                effect = "INHIBIT"
            if knob in seen:
                continue

            seen.add(knob)
            normalized.append(
                {
                    "knob": knob,
                    "effect": effect,
                    "strength": round(self._clamp(strength, 0.1, 1.0), 6),
                }
            )
            if len(normalized) >= self.config.max_interventions:
                break

        if not normalized:
            normalized.append(self._random_intervention())
        return normalized

    def _individual_to_json(self, individual: Sequence[Dict[str, Any]]) -> Dict[str, Any]:
        interventions = []
        for idx, entry in enumerate(self._normalize_interventions(individual)):
            knob = entry["knob"]
            effect = entry["effect"]
            strength = float(entry["strength"])
            interventions.append(
                {
                    "knob": knob,
                    "effect": effect,
                    "strength": self._clamp(strength, 0.0, 1.0),
                    "name": f"{effect}_{knob}_{idx}",
                }
            )
        return {
            "calibration_profile": self.config.calibration_profile,
            "knob_interventions": interventions,
        }

    def _individual_from_plain(self, interventions: Sequence[Dict[str, Any]]):
        normalized = self._normalize_interventions(interventions)
        return creator.StromaInterventionIndividual(normalized)

    def _plain_interventions(self, individual: Optional[Sequence[Dict[str, Any]]]) -> List[Dict[str, Any]]:
        if individual is None:
            return []
        return [dict(entry) for entry in individual]

    @staticmethod
    def _clamp(x: float, lo: float, hi: float) -> float:
        if math.isnan(x):
            return lo
        if x < lo:
            return lo
        if x > hi:
            return hi
        return x
