from __future__ import annotations

from typing import Dict, List

TARGETABLE_KNOBS: List[str] = [
    "tgfb_secretion_rate",
    "shh_secretion_rate",
    "efflux_induction_delay",
    "efflux_strength",
]

FIXED_KNOBS: List[str] = [
    "tgfb_brake_sensitivity",
    "proliferation_rate",
    "checkpoint_integrity",
    "apoptosis_resistance",
    "mechanical_compaction_strength",
]

OBSERVABLE_KNOBS: List[str] = [
    "emt_induction_threshold",
    "emt_phenotype_extent",
    "hypoxia_response_threshold",
    "hypoxia_phenotype_shift",
]

LEGACY_GENE_BRIDGE: Dict[str, str] = {
    "TGFB1": "tgfb_secretion_rate",
    "SHH": "shh_secretion_rate",
    "ABCB1": "efflux_strength",
}


def all_knobs() -> List[str]:
    return list(FIXED_KNOBS + OBSERVABLE_KNOBS + TARGETABLE_KNOBS)


def validate_partition() -> None:
    knobs = all_knobs()
    if len(knobs) != 13:
        raise ValueError("Expected exactly 13 knobs in schema.")
    if len(FIXED_KNOBS) != 5 or len(OBSERVABLE_KNOBS) != 4 or len(TARGETABLE_KNOBS) != 4:
        raise ValueError("Expected partition counts Fixed=5, Observable=4, Targetable=4.")
    expected_targetable = {"tgfb_secretion_rate", "shh_secretion_rate", "efflux_induction_delay", "efflux_strength"}
    if set(TARGETABLE_KNOBS) != expected_targetable:
        raise ValueError("Targetable knobs must be exactly 1a/1b/7a/7b.")


def map_legacy_gene_to_knob(gene: str) -> str:
    try:
        return LEGACY_GENE_BRIDGE[gene]
    except KeyError as exc:
        raise ValueError(
            f"Partition violation: legacy gene '{gene}' is not bridgeable. "
            "Allowed: TGFB1, SHH, ABCB1."
        ) from exc
