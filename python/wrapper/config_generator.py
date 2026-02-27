from __future__ import annotations

import json
import logging
import shutil
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any, Dict, List, Tuple

try:
    from ..ea.knob_schema import map_legacy_gene_to_knob
except Exception:  # pragma: no cover
    from python.ea.knob_schema import map_legacy_gene_to_knob  # type: ignore

LOGGER = logging.getLogger(__name__)


def _clamp(x: float, lo: float, hi: float) -> float:
    if x < lo:
        return lo
    if x > hi:
        return hi
    return x


class ConfigGenerator:
    """
    Generate per-run PhysiCell inputs from an EA individual.

    The generator copies a base XML config and base gene-parameter JSON into a
    run directory, emits an intervention JSON payload, and optionally patches
    drug-related XML settings based on intervention content.
    """

    # Knobs treated as affecting drug-delivery context.
    DEFAULT_DRUG_RELATED_KNOBS = {
        "tgfb_secretion_rate",
        "shh_secretion_rate",
        "efflux_induction_delay",
        "efflux_strength",
    }

    def __init__(self, base_config_xml: Path, base_gene_params: Path):
        self.base_config_xml = Path(base_config_xml).expanduser().resolve()
        self.base_gene_params = Path(base_gene_params).expanduser().resolve()

        if not self.base_config_xml.exists():
            raise FileNotFoundError(f"Base PhysiCell XML config not found: {self.base_config_xml}")
        if not self.base_gene_params.exists():
            raise FileNotFoundError(f"Base gene params JSON not found: {self.base_gene_params}")

        with self.base_gene_params.open("r", encoding="utf-8") as f:
            self.gene_params_template = json.load(f)
        if not isinstance(self.gene_params_template, dict):
            raise ValueError(f"Base gene params must be a JSON object: {self.base_gene_params}")

        self.drug_related_knobs = set(self.DEFAULT_DRUG_RELATED_KNOBS)

    def generate(self, individual: List[Dict[str, Any]], run_dir: Path) -> Tuple[Path, Path]:
        """
        Build per-run files from a single EA individual.

        Returns:
            (config_xml_path, intervention_json_path)
        """
        run_dir = Path(run_dir).expanduser().resolve()
        run_dir.mkdir(parents=True, exist_ok=True)

        config_xml_path = run_dir / "config.xml"
        intervention_json_path = run_dir / "intervention.json"
        run_gene_params_path = run_dir / "gene_params.json"
        run_output_dir = run_dir / "output"
        run_output_dir.mkdir(parents=True, exist_ok=True)

        shutil.copy2(self.base_config_xml, config_xml_path)
        shutil.copy2(self.base_gene_params, run_gene_params_path)

        normalized_interventions = self._normalize_individual(individual)
        drug_protocol = self.generate_drug_protocol(normalized_interventions)

        self._patch_xml_config(
            config_xml_path=config_xml_path,
            output_dir=run_output_dir,
            drug_protocol=drug_protocol,
        )

        intervention_payload = {
            "knob_interventions": normalized_interventions,
            "drug_delivery": drug_protocol,
            "gene_params_file": str(run_gene_params_path),
        }
        with intervention_json_path.open("w", encoding="utf-8") as f:
            json.dump(intervention_payload, f, indent=2, sort_keys=True)
            f.write("\n")

        return config_xml_path, intervention_json_path

    def generate_drug_protocol(self, individual: List[Dict[str, Any]]) -> Dict[str, Any]:
        """
        Build drug-delivery protocol from an intervention list.

        Rule:
        - INHIBIT intervention on a drug-related knob contributes dose
          (strength * max_drug_concentration)
        - interventions marked as non-drug delivery are skipped
        """
        max_drug_concentration = 1.0
        start_time = 2880.0
        mode = "constant_boundary"

        targets: List[Dict[str, Any]] = []
        boundary_concentration = 0.0

        for entry in self._normalize_individual(individual):
            knob = entry["knob"]
            effect = entry["effect"]
            strength = float(entry["strength"])

            if effect != "INHIBIT":
                continue
            if knob not in self.drug_related_knobs:
                continue

            # Allow callers to explicitly mark non-drug interventions.
            delivery_mode = str(entry.get("delivery", "")).strip().lower()
            if delivery_mode in {"genetic", "none"}:
                continue
            if entry.get("requires_drug") is False:
                continue

            dose = _clamp(strength, 0.0, 1.0) * max_drug_concentration
            boundary_concentration = max(boundary_concentration, dose)
            targets.append(
                {
                    "knob": knob,
                    "effect": effect,
                    "strength": strength,
                    "dose": dose,
                }
            )

        enabled = boundary_concentration > 0.0
        return {
            "enabled": enabled,
            "mode": mode,
            "start_time": start_time,
            "max_drug_concentration": max_drug_concentration,
            "boundary_concentration": boundary_concentration,
            "targets": targets,
        }

    def _normalize_individual(self, individual: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        if not isinstance(individual, list):
            raise TypeError("individual must be a list of intervention dicts")

        normalized: List[Dict[str, Any]] = []
        for idx, raw in enumerate(individual):
            if not isinstance(raw, dict):
                LOGGER.warning("Skipping intervention[%d]: expected dict, got %s", idx, type(raw).__name__)
                continue

            knob = str(raw.get("knob", "")).strip()
            if not knob:
                legacy_gene = str(raw.get("gene", "")).strip()
                if legacy_gene:
                    try:
                        knob = map_legacy_gene_to_knob(legacy_gene)
                    except ValueError:
                        LOGGER.warning("Skipping intervention[%d]: unsupported legacy gene '%s'", idx, legacy_gene)
                        continue
            if not knob:
                LOGGER.warning("Skipping intervention[%d]: missing knob", idx)
                continue

            effect = str(raw.get("effect", "INHIBIT")).strip().upper()
            if effect not in {"INHIBIT", "ACTIVATE"}:
                effect = "INHIBIT"

            try:
                strength = float(raw.get("strength", 0.0))
            except (TypeError, ValueError):
                strength = 0.0
            strength = _clamp(strength, 0.0, 1.0)

            item = {
                "knob": knob,
                "effect": effect,
                "strength": strength,
            }

            # Preserve optional metadata fields if present.
            if "name" in raw:
                item["name"] = str(raw["name"])
            if "delivery" in raw:
                item["delivery"] = str(raw["delivery"])
            if "requires_drug" in raw:
                item["requires_drug"] = bool(raw["requires_drug"])

            normalized.append(item)

        return normalized

    def _patch_xml_config(self, config_xml_path: Path, output_dir: Path, drug_protocol: Dict[str, Any]) -> None:
        tree = ET.parse(config_xml_path)
        root = tree.getroot()

        # Keep each run self-contained.
        folder_node = root.find("./save/folder")
        if folder_node is None:
            folder_node = root.find(".//folder")
        if folder_node is not None:
            folder_node.text = str(output_dir)

        self._apply_drug_protocol_to_xml(root, drug_protocol)

        tree.write(config_xml_path, encoding="utf-8", xml_declaration=True)

    def _apply_drug_protocol_to_xml(self, root: ET.Element, drug_protocol: Dict[str, Any]) -> None:
        enabled = bool(drug_protocol.get("enabled", False))
        boundary_conc = float(drug_protocol.get("boundary_concentration", 0.0))
        start_time = float(drug_protocol.get("start_time", 2880.0))

        # Update user parameters where available.
        self._set_user_parameter(root, "drug_start_time", start_time)
        self._set_user_parameter(root, "drug_concentration", boundary_conc)

        # Locate "drug" variable and set boundary condition concentration.
        drug_var = root.find(".//microenvironment_setup/variable[@name='drug']")
        if drug_var is None:
            LOGGER.warning("No microenvironment variable named 'drug' found; skipping XML drug patch.")
            return

        dbc = drug_var.find("./Dirichlet_boundary_condition")
        if dbc is not None:
            dbc.text = f"{boundary_conc:.6g}"
            dbc.set("enabled", "true" if enabled else "false")

        # Update explicit boundary entries for the drug substrate.
        for bv in drug_var.findall("./Dirichlet_options/boundary_value"):
            bv.text = f"{boundary_conc:.6g}"
            bv.set("enabled", "true" if enabled else "false")

        ic = drug_var.find("./initial_condition")
        if ic is not None and not enabled:
            # Keep drug absent at t=0 when no protocol is active.
            ic.text = "0"

    @staticmethod
    def _set_user_parameter(root: ET.Element, name: str, value: float) -> None:
        node = root.find(f".//user_parameters/{name}")
        if node is not None:
            node.text = f"{float(value):.6g}"
