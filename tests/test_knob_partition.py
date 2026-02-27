from __future__ import annotations

import pytest

pytest.importorskip("deap")

from python.ea.evolutionary_algorithm import EAConfig, StromaWorldEA
from python.ea.knob_schema import (
    FIXED_KNOBS,
    LEGACY_GENE_BRIDGE,
    OBSERVABLE_KNOBS,
    TARGETABLE_KNOBS,
    map_legacy_gene_to_knob,
    validate_partition,
)


def test_partition_membership_and_counts():
    validate_partition()
    assert len(FIXED_KNOBS) == 5
    assert len(OBSERVABLE_KNOBS) == 4
    assert len(TARGETABLE_KNOBS) == 4
    assert set(TARGETABLE_KNOBS) == {
        "tgfb_secretion_rate",
        "shh_secretion_rate",
        "efflux_induction_delay",
        "efflux_strength",
    }


def test_reject_fixed_and_observable_overrides_from_ea_path(ea_test_paths):
    cfg = EAConfig(
        population_size=4,
        generations=1,
        binary_path=str(ea_test_paths["binary"]),
        config_path=str(ea_test_paths["config"]),
        simulation_output_dir=str(ea_test_paths["sim_output_dir"]),
        stats_csv_path=str(ea_test_paths["stats_csv"]),
        random_seed=3,
    )
    ea = StromaWorldEA(cfg)

    with pytest.raises(ValueError, match="Partition violation"):
        ea._individual_from_plain(  # noqa: SLF001
            [{"knob": "tgfb_brake_sensitivity", "effect": "INHIBIT", "strength": 1.0}]
        )

    with pytest.raises(ValueError, match="Partition violation"):
        ea._individual_from_plain(  # noqa: SLF001
            [{"knob": "emt_induction_threshold", "effect": "INHIBIT", "strength": 1.0}]
        )


def test_legacy_bridge_behavior():
    assert map_legacy_gene_to_knob("TGFB1") == LEGACY_GENE_BRIDGE["TGFB1"]
    assert map_legacy_gene_to_knob("SHH") == LEGACY_GENE_BRIDGE["SHH"]
    assert map_legacy_gene_to_knob("ABCB1") == LEGACY_GENE_BRIDGE["ABCB1"]

    for non_bridgeable in ["EGFR", "HAS2", "MMP2", "GLI1", "BCL_XL"]:
        with pytest.raises(ValueError, match="Partition violation"):
            map_legacy_gene_to_knob(non_bridgeable)
