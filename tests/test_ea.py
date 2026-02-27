from __future__ import annotations

import pytest

pytest.importorskip("deap")

from python.ea.evolutionary_algorithm import ALLOWED_EFFECTS, EAConfig, StromaWorldEA


def _assert_individual_valid(individual, cfg: EAConfig) -> None:
    assert 1 <= len(individual) <= cfg.max_interventions
    knobs = [entry["knob"] for entry in individual]
    assert len(knobs) == len(set(knobs))

    for entry in individual:
        assert entry["knob"] in cfg.targetable_knobs
        assert entry["effect"] in ALLOWED_EFFECTS
        assert 0.1 <= float(entry["strength"]) <= 1.0


@pytest.fixture
def ea_instance(ea_test_paths, tmp_path):
    cfg = EAConfig(
        population_size=8,
        generations=1,
        mutation_rate=1.0,
        crossover_rate=1.0,
        tournament_size=2,
        max_interventions=4,
        parallel_sims=1,
        binary_path=str(ea_test_paths["binary"]),
        config_path=str(ea_test_paths["config"]),
        simulation_output_dir=str(ea_test_paths["sim_output_dir"]),
        stats_csv_path=str(ea_test_paths["stats_csv"]),
        random_seed=7,
        keep_population_history=True,
    )
    return StromaWorldEA(cfg)


def test_individual_generation_valid_interventions_only(ea_instance):
    population = ea_instance.initialize_population()
    assert len(population) == ea_instance.config.population_size
    for individual in population:
        _assert_individual_valid(individual, ea_instance.config)


def test_mutation_stays_within_bounds(ea_instance):
    individual = ea_instance._individual_from_plain(  # noqa: SLF001
        [
            {"knob": "tgfb_secretion_rate", "effect": "INHIBIT", "strength": 0.5},
            {"knob": "efflux_strength", "effect": "INHIBIT", "strength": 0.7},
        ]
    )

    for _ in range(25):
        ea_instance.mutate(individual)
        _assert_individual_valid(individual, ea_instance.config)


def test_crossover_no_duplicate_gene_targets(ea_instance):
    ind1 = ea_instance._individual_from_plain(  # noqa: SLF001
        [
            {"knob": "tgfb_secretion_rate", "effect": "INHIBIT", "strength": 0.4},
            {"knob": "efflux_strength", "effect": "INHIBIT", "strength": 0.8},
            {"knob": "efflux_induction_delay", "effect": "INHIBIT", "strength": 0.6},
        ]
    )
    ind2 = ea_instance._individual_from_plain(  # noqa: SLF001
        [
            {"knob": "tgfb_secretion_rate", "effect": "ACTIVATE", "strength": 0.3},
            {"knob": "shh_secretion_rate", "effect": "INHIBIT", "strength": 0.5},
            {"knob": "efflux_strength", "effect": "INHIBIT", "strength": 0.9},
        ]
    )

    child1, child2 = ea_instance.crossover(ind1, ind2)
    _assert_individual_valid(child1, ea_instance.config)
    _assert_individual_valid(child2, ea_instance.config)


def test_ea_runs_one_generation_with_mock_fitness(ea_instance, monkeypatch):
    def fake_evaluate(individual):
        avg_strength = sum(float(entry["strength"]) for entry in individual) / len(individual)
        return max(0.0, min(1.0, avg_strength))

    monkeypatch.setattr(ea_instance, "evaluate", fake_evaluate)

    result = ea_instance.run()
    assert 0.0 <= result.best_fitness <= 1.0
    assert len(result.best_individual) >= 1
    assert len(result.fitness_history) >= 2  # generation 0 + generation 1


def test_payload_emits_knob_interventions_only(ea_instance):
    individual = ea_instance._individual_from_plain(  # noqa: SLF001
        [
            {"knob": "tgfb_secretion_rate", "effect": "INHIBIT", "strength": 0.5},
            {"knob": "efflux_strength", "effect": "ACTIVATE", "strength": 0.3},
        ]
    )
    payload = ea_instance._individual_to_json(individual)  # noqa: SLF001
    assert "knob_interventions" in payload
    assert "interventions" not in payload
    knobs = {entry["knob"] for entry in payload["knob_interventions"]}
    assert knobs.issubset(set(ea_instance.config.targetable_knobs))
