"""
Population utilities for knob-based EA individuals.

Each individual is a list of interventions:
    [{"knob": "...", "effect": "INHIBIT|ACTIVATE", "strength": 0..1}, ...]
"""

from __future__ import annotations

import random
from typing import Dict, List, Optional

from .knob_schema import TARGETABLE_KNOBS as DEFAULT_TARGETABLE_KNOBS, validate_partition as validate_knob_partition

ALLOWED_EFFECTS: tuple = ("INHIBIT", "ACTIVATE")
INHIBIT_WEIGHT: float = 0.8

validate_knob_partition()


def random_intervention(
    targetable_knobs: Optional[List[str]] = None,
    inhibit_weight: float = INHIBIT_WEIGHT,
) -> Dict:
    knobs = targetable_knobs if targetable_knobs is not None else DEFAULT_TARGETABLE_KNOBS
    knob = random.choice(knobs)
    effect = "INHIBIT" if random.random() < inhibit_weight else "ACTIVATE"
    strength = round(random.uniform(0.1, 1.0), 4)
    return {"knob": knob, "effect": effect, "strength": strength}


def random_individual(
    max_interventions: int = 4,
    min_interventions: int = 1,
    targetable_knobs: Optional[List[str]] = None,
    inhibit_weight: float = INHIBIT_WEIGHT,
) -> List[Dict]:
    knobs = targetable_knobs if targetable_knobs is not None else DEFAULT_TARGETABLE_KNOBS
    n = random.randint(min_interventions, min(max_interventions, len(knobs)))
    selected = random.sample(knobs, n)
    individual = []
    for knob in selected:
        effect = "INHIBIT" if random.random() < inhibit_weight else "ACTIVATE"
        strength = round(random.uniform(0.1, 1.0), 4)
        individual.append({"knob": knob, "effect": effect, "strength": strength})
    return individual


def initialize_population(
    population_size: int,
    max_interventions: int = 4,
    min_interventions: int = 1,
    targetable_knobs: Optional[List[str]] = None,
    inhibit_weight: float = INHIBIT_WEIGHT,
) -> List[List[Dict]]:
    return [
        random_individual(
            max_interventions=max_interventions,
            min_interventions=min_interventions,
            targetable_knobs=targetable_knobs,
            inhibit_weight=inhibit_weight,
        )
        for _ in range(population_size)
    ]


def validate_individual(
    individual: List[Dict],
    targetable_knobs: Optional[List[str]] = None,
    max_interventions: int = 4,
) -> bool:
    knobs = set(targetable_knobs) if targetable_knobs is not None else set(DEFAULT_TARGETABLE_KNOBS)
    if not individual or len(individual) > max_interventions:
        return False

    seen_knobs: set = set()
    for entry in individual:
        if not isinstance(entry, dict):
            return False
        if not {"knob", "effect", "strength"}.issubset(entry):
            return False
        if entry["knob"] not in knobs:
            return False
        if entry["effect"] not in ALLOWED_EFFECTS:
            return False
        s = entry["strength"]
        if not (isinstance(s, (int, float)) and 0.0 <= float(s) <= 1.0):
            return False
        if entry["knob"] in seen_knobs:
            return False
        seen_knobs.add(entry["knob"])
    return True


def individual_to_json_payload(
    individual: List[Dict],
    calibration_profile: str = "AsPC-1",
) -> dict:
    entries = []
    for i, iv in enumerate(individual):
        entries.append(
            {
                "knob": iv["knob"],
                "effect": iv["effect"],
                "strength": float(iv["strength"]),
                "name": f"EA_{iv['knob']}_{i}",
            }
        )
    return {"calibration_profile": calibration_profile, "knob_interventions": entries}
