"""Genetic operators for knob-based EA individuals."""

from __future__ import annotations

import copy
import random
from typing import Dict, List, Optional, Tuple

from .population import DEFAULT_TARGETABLE_KNOBS


def mutate_individual(
    individual: List[Dict],
    mutation_rate: float = 0.2,
    max_interventions: int = 4,
    targetable_knobs: Optional[List[str]] = None,
) -> List[Dict]:
    knobs = targetable_knobs if targetable_knobs is not None else DEFAULT_TARGETABLE_KNOBS
    ind = copy.deepcopy(individual)

    for i in range(len(ind)):
        if random.random() >= mutation_rate:
            continue

        r = random.random()
        if r < 0.40:
            delta = random.gauss(0.0, 0.1)
            new_strength = float(ind[i]["strength"]) + delta
            ind[i]["strength"] = round(max(0.1, min(1.0, new_strength)), 4)
        elif r < 0.60:
            current_knobs = {iv["knob"] for iv in ind}
            candidates = [k for k in knobs if k != ind[i]["knob"] and k not in current_knobs]
            if candidates:
                ind[i]["knob"] = random.choice(candidates)
        elif r < 0.70:
            ind[i]["effect"] = "ACTIVATE" if ind[i]["effect"] == "INHIBIT" else "INHIBIT"
        elif r < 0.85:
            if len(ind) < max_interventions:
                existing_knobs = {iv["knob"] for iv in ind}
                candidates = [k for k in knobs if k not in existing_knobs]
                if candidates:
                    ind.append(
                        {
                            "knob": random.choice(candidates),
                            "effect": "INHIBIT" if random.random() < 0.8 else "ACTIVATE",
                            "strength": round(random.uniform(0.1, 1.0), 4),
                        }
                    )
        else:
            if len(ind) > 1:
                del ind[random.randrange(len(ind))]
                break

    return ind


def crossover_individuals(
    ind1: List[Dict],
    ind2: List[Dict],
    max_interventions: int = 4,
) -> Tuple[List[Dict], List[Dict]]:
    pool: Dict[str, Dict] = {}
    for iv in ind1:
        pool[iv["knob"]] = copy.deepcopy(iv)
    for iv in ind2:
        pool[iv["knob"]] = copy.deepcopy(iv)

    pool_list = list(pool.values())
    random.shuffle(pool_list)

    def _sample_child() -> List[Dict]:
        seen: set = set()
        child: List[Dict] = []
        for iv in pool_list:
            if iv["knob"] not in seen and len(child) < max_interventions:
                child.append(copy.deepcopy(iv))
                seen.add(iv["knob"])
        if not child and pool_list:
            child = [copy.deepcopy(random.choice(pool_list))]
        return child

    child1 = _sample_child()
    random.shuffle(pool_list)
    child2 = _sample_child()
    return child1, child2


def tournament_select(
    population: List,
    fitnesses: List[float],
    tournament_size: int = 3,
) -> object:
    contestants = random.sample(range(len(population)), min(tournament_size, len(population)))
    winner_idx = max(contestants, key=lambda i: fitnesses[i])
    return population[winner_idx]


def elitist_survivors(
    population: List,
    fitnesses: List[float],
    elite_fraction: float = 0.05,
) -> Tuple[List, List[float]]:
    n_elite = max(1, int(len(population) * elite_fraction))
    ranked = sorted(zip(population, fitnesses), key=lambda x: x[1], reverse=True)
    elites = ranked[:n_elite]
    elite_inds = [copy.deepcopy(e[0]) for e in elites]
    elite_fits = [e[1] for e in elites]
    return elite_inds, elite_fits
