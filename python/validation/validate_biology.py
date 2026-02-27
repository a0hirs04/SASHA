from __future__ import annotations

import argparse
import csv
import json
import logging
import math
import shutil
import sys
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple
import xml.etree.ElementTree as ET

# Ensure project-root imports work when running as a script.
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

try:
    from ..ea.fitness import compute_fitness
    from ..wrapper.output_parser import OutputParser, SimulationMetrics
    from ..wrapper.physicell_runner import PhysiCellRunner
except Exception:  # pragma: no cover
    from python.ea.fitness import compute_fitness  # type: ignore
    from python.wrapper.output_parser import OutputParser, SimulationMetrics  # type: ignore
    from python.wrapper.physicell_runner import PhysiCellRunner  # type: ignore


LOGGER = logging.getLogger(__name__)


@dataclass
class ScenarioResult:
    scenario: str
    success: bool
    criterion: str
    details: str
    run_dir: str
    exit_code: int
    wall_time: float
    fitness: float
    metrics: Dict[str, Any]
    extras: Dict[str, Any]


@dataclass
class ValidationSummary:
    total: int
    passed: int
    failed: int
    scenario_results: List[ScenarioResult]


@dataclass
class ScenarioSpec:
    name: str
    description: str
    criterion: str
    intervention_payload: Dict[str, Any]
    user_parameter_overrides: Dict[str, Any]
    variable_overrides: Dict[str, Dict[str, Any]]
    evaluator: Callable[[SimulationMetrics, Dict[str, Any], Dict[str, Any]], Tuple[bool, str]]


class BiologyValidator:
    def __init__(
        self,
        binary_path: Path | str,
        config_path: Path | str,
        output_dir: Path | str,
        timeout_seconds: int = 7200,
        sim_max_time: Optional[float] = None,
    ):
        self.binary_path = Path(binary_path).expanduser().resolve()
        self.base_config_path = Path(config_path).expanduser().resolve()
        self.output_dir = Path(output_dir).expanduser().resolve()
        self.timeout_seconds = int(timeout_seconds)
        self.sim_max_time: Optional[float] = sim_max_time

        if not self.binary_path.exists():
            raise FileNotFoundError(f"PhysiCell binary not found: {self.binary_path}")
        if not self.base_config_path.exists():
            raise FileNotFoundError(f"PhysiCell config not found: {self.base_config_path}")

        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.initial_tumor_cells = self._read_user_param_int(self.base_config_path, "number_of_tumor_cells", default=50)
        self._baseline_cache: Dict[str, Any] = {}

    def run_all(self) -> ValidationSummary:
        specs = self._scenario_specs()
        scenario_results: List[ScenarioResult] = []
        for spec in specs:
            LOGGER.info("Running validation scenario: %s", spec.name)
            result = self._run_scenario(spec)
            scenario_results.append(result)
            self._baseline_cache[spec.name] = {
                "fitness": result.fitness,
                "metrics": result.metrics,
                "extras": result.extras,
            }

        passed = sum(1 for r in scenario_results if r.success)
        failed = len(scenario_results) - passed
        return ValidationSummary(
            total=len(scenario_results),
            passed=passed,
            failed=failed,
            scenario_results=scenario_results,
        )

    def save_summary(self, summary: ValidationSummary, output_prefix: Path | str) -> Tuple[Path, Path]:
        output_prefix = Path(output_prefix).expanduser().resolve()
        output_prefix.parent.mkdir(parents=True, exist_ok=True)
        json_path = output_prefix.with_suffix(".json")
        csv_path = output_prefix.with_suffix(".csv")

        payload = {
            "total": summary.total,
            "passed": summary.passed,
            "failed": summary.failed,
            "scenario_results": [asdict(r) for r in summary.scenario_results],
        }
        json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

        fieldnames = [
            "scenario",
            "success",
            "criterion",
            "details",
            "run_dir",
            "exit_code",
            "wall_time",
            "fitness",
        ]
        with csv_path.open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for result in summary.scenario_results:
                writer.writerow(
                    {
                        "scenario": result.scenario,
                        "success": result.success,
                        "criterion": result.criterion,
                        "details": result.details,
                        "run_dir": result.run_dir,
                        "exit_code": result.exit_code,
                        "wall_time": f"{result.wall_time:.3f}",
                        "fitness": f"{result.fitness:.6f}",
                    }
                )
        return json_path, csv_path

    def _run_scenario(self, spec: ScenarioSpec) -> ScenarioResult:
        scenario_dir = self.output_dir / f"{datetime.utcnow().strftime('%Y%m%d_%H%M%S')}_{spec.name.lower()}"
        scenario_dir.mkdir(parents=True, exist_ok=True)

        config_copy = scenario_dir / "scenario_config.xml"
        shutil.copy2(self.base_config_path, config_copy)
        self._patch_config(
            config_copy,
            user_parameter_overrides=spec.user_parameter_overrides,
            variable_overrides=spec.variable_overrides,
        )

        runner = PhysiCellRunner(
            binary_path=self.binary_path,
            config_path=config_copy,
            output_dir=scenario_dir / "sim_runs",
            timeout_seconds=self.timeout_seconds,
        )

        sim_result = runner.run(spec.intervention_payload)
        if not sim_result.success:
            return ScenarioResult(
                scenario=spec.name,
                success=False,
                criterion=spec.criterion,
                details=f"simulation failed (exit_code={sim_result.exit_code})",
                run_dir=str(sim_result.run_dir),
                exit_code=sim_result.exit_code,
                wall_time=sim_result.wall_time,
                fitness=0.0,
                metrics={},
                extras={},
            )

        parser = OutputParser(sim_result.run_dir / "output")
        metrics = parser.parse_final_state()
        fitness = float(compute_fitness(metrics))
        extras = self._compute_extras(parser, metrics)

        ok, details = spec.evaluator(metrics, extras, self._baseline_cache)
        return ScenarioResult(
            scenario=spec.name,
            success=ok,
            criterion=spec.criterion,
            details=details,
            run_dir=str(sim_result.run_dir),
            exit_code=sim_result.exit_code,
            wall_time=sim_result.wall_time,
            fitness=fitness,
            metrics=asdict(metrics),
            extras=extras,
        )

    def _compute_extras(self, parser: OutputParser, metrics: SimulationMetrics) -> Dict[str, Any]:
        extras: Dict[str, Any] = {}

        try:
            df = parser.parse_timeseries()
            extras["timeseries_rows"] = int(len(df))
            if len(df) >= 2:
                extras["ecm_delta"] = float(df["mean_ecm"].iloc[-1] - df["mean_ecm"].iloc[0])
                extras["activated_cafs_delta"] = float(df["activated_cafs"].iloc[-1] - df["activated_cafs"].iloc[0])
                extras["tumor_count_delta"] = float(df["tumor_count"].iloc[-1] - df["tumor_count"].iloc[0])
            else:
                extras["ecm_delta"] = math.nan
                extras["activated_cafs_delta"] = math.nan
                extras["tumor_count_delta"] = math.nan
        except Exception as exc:  # pragma: no cover
            LOGGER.warning("Timeseries parsing failed: %s", exc)
            extras["timeseries_rows"] = 0

        # Gene-level summary for resistance scenario (NRF2 / ABCB1 in surviving tumor cells).
        extras["mean_nrf2_surviving_tumor"] = self._mean_gene_in_surviving_tumor(parser, "NRF2")
        extras["mean_abcb1_surviving_tumor"] = self._mean_gene_in_surviving_tumor(parser, "ABCB1")
        extras["activated_caf_fraction"] = (
            float(metrics.activated_cafs) / float(metrics.total_stromal_cells)
            if metrics.total_stromal_cells > 0
            else math.nan
        )
        return extras

    def _mean_gene_in_surviving_tumor(self, parser: OutputParser, gene_name: str) -> float:
        try:
            final_xml = parser._find_final_snapshot_xml()  # noqa: SLF001
            snapshot = parser._read_physicell_xml(final_xml)  # noqa: SLF001
            matrix = snapshot["cell_matrix"]
            label_name_map = snapshot["label_name_map"]
        except Exception:  # pragma: no cover
            return math.nan

        gene_row = self._row_by_label(matrix, label_name_map, gene_name)
        cell_type = self._row_by_label(matrix, label_name_map, "cell_type")
        dead = self._row_by_label(matrix, label_name_map, "dead")
        death_model = self._row_by_label(matrix, label_name_map, "current_death_model")

        if gene_row is None or cell_type is None:
            return math.nan

        tumor_mask = np_rint(cell_type) == 0
        if dead is not None:
            tumor_mask &= dead <= 0.5
        if death_model is not None:
            tumor_mask &= np_rint(death_model) != 100

        if not tumor_mask.any():
            return math.nan
        return float(gene_row[tumor_mask].mean())

    @staticmethod
    def _row_by_label(matrix, label_name_map: Dict[str, Dict[str, Any]], label: str):
        entry = label_name_map.get(label)
        if entry is None:
            return None
        idx = int(entry.get("index", -1))
        if idx < 0 or idx >= matrix.shape[0]:
            return None
        return matrix[idx, :]

    def _patch_config(
        self,
        config_path: Path,
        user_parameter_overrides: Dict[str, Any],
        variable_overrides: Dict[str, Dict[str, Any]],
    ) -> None:
        tree = ET.parse(config_path)
        root = tree.getroot()

        if self.sim_max_time is not None:
            for tag in ("max_time", "overall_time"):
                node = root.find(f".//{tag}")
                if node is not None:
                    node.text = str(self.sim_max_time)

        for key, value in user_parameter_overrides.items():
            node = root.find(f".//user_parameters/{key}")
            if node is None:
                continue
            node.text = str(value)

        for variable_name, overrides in variable_overrides.items():
            var_node = root.find(f".//microenvironment_setup/variable[@name='{variable_name}']")
            if var_node is None:
                continue

            if "Dirichlet_boundary_condition" in overrides:
                dbc = var_node.find("./Dirichlet_boundary_condition")
                if dbc is not None:
                    cfg = overrides["Dirichlet_boundary_condition"]
                    if "value" in cfg:
                        dbc.text = str(cfg["value"])
                    if "enabled" in cfg:
                        dbc.set("enabled", "true" if bool(cfg["enabled"]) else "false")

            if "Dirichlet_options" in overrides:
                bcfg = overrides["Dirichlet_options"]
                for bnode in var_node.findall("./Dirichlet_options/boundary_value"):
                    if "value" in bcfg:
                        bnode.text = str(bcfg["value"])
                    if "enabled" in bcfg:
                        bnode.set("enabled", "true" if bool(bcfg["enabled"]) else "false")

            if "initial_condition" in overrides:
                inode = var_node.find("./initial_condition")
                if inode is not None:
                    inode.text = str(overrides["initial_condition"])

        tree.write(config_path, encoding="utf-8", xml_declaration=True)

    @staticmethod
    def _read_user_param_int(config_path: Path, name: str, default: int) -> int:
        try:
            tree = ET.parse(config_path)
            root = tree.getroot()
            node = root.find(f".//user_parameters/{name}")
            if node is not None and node.text is not None:
                return int(float(node.text.strip()))
        except Exception:  # pragma: no cover
            pass
        return int(default)

    def _scenario_specs(self) -> List[ScenarioSpec]:
        no_cytotoxic_user = {"drug_start_time": 1.0e9, "drug_concentration": 0.0}
        no_cytotoxic_var = {
            "drug": {
                "Dirichlet_boundary_condition": {"enabled": False, "value": 0.0},
                "Dirichlet_options": {"enabled": False, "value": 0.0},
                "initial_condition": 0.0,
            }
        }
        cytotoxic_user = {"drug_start_time": 0.0, "drug_concentration": 1.0}
        cytotoxic_var = {
            "drug": {
                "Dirichlet_boundary_condition": {"enabled": True, "value": 1.0},
                "Dirichlet_options": {"enabled": True, "value": 1.0},
                "initial_condition": 1.0,
            }
        }
        high_tgfb_var = {
            "tgfb": {
                "Dirichlet_boundary_condition": {"enabled": True, "value": 1.0},
                "Dirichlet_options": {"enabled": True, "value": 1.0},
                "initial_condition": 1.0,
            }
        }

        specs: List[ScenarioSpec] = [
            ScenarioSpec(
                name="BASELINE_NO_CYTOTOXIC",
                description="AsPC-1 baseline with no cytotoxic exposure.",
                criterion="reference baseline is produced for SHH paradox comparison",
                intervention_payload={"calibration_profile": "AsPC-1", "knob_interventions": []},
                user_parameter_overrides=no_cytotoxic_user,
                variable_overrides=no_cytotoxic_var,
                evaluator=self._eval_baseline_no_cytotoxic,
            ),
            ScenarioSpec(
                name="SHH_OFF_NO_CYTOTOXIC",
                description="Knob 1b inhibition without cytotoxic context.",
                criterion="worse tumor outcome than BASELINE_NO_CYTOTOXIC",
                intervention_payload={
                    "calibration_profile": "AsPC-1",
                    "knob_interventions": [
                        {"knob": "shh_secretion_rate", "effect": "INHIBIT", "strength": 1.0, "name": "SHH_OFF"}
                    ],
                },
                user_parameter_overrides=no_cytotoxic_user,
                variable_overrides=no_cytotoxic_var,
                evaluator=self._eval_shh_off_no_cytotoxic,
            ),
            ScenarioSpec(
                name="SHH_OFF_WITH_CYTOTOXIC",
                description="Knob 1b inhibition with cytotoxic enabled.",
                criterion="improved tumor kill vs SHH_OFF_NO_CYTOTOXIC",
                intervention_payload={
                    "calibration_profile": "AsPC-1",
                    "knob_interventions": [
                        {"knob": "shh_secretion_rate", "effect": "INHIBIT", "strength": 1.0, "name": "SHH_OFF"}
                    ],
                },
                user_parameter_overrides=cytotoxic_user,
                variable_overrides=cytotoxic_var,
                evaluator=self._eval_shh_off_with_cytotoxic,
            ),
            ScenarioSpec(
                name="KNOB2_ASPC1_CONTROL",
                description="AsPC-1 under high TGF-beta exposure for Knob 2 comparator.",
                criterion="reference control for KNOB2_PANC1_CONTROL",
                intervention_payload={"calibration_profile": "AsPC-1", "knob_interventions": []},
                user_parameter_overrides=no_cytotoxic_user,
                variable_overrides={**no_cytotoxic_var, **high_tgfb_var},
                evaluator=self._eval_knob2_aspc1_control,
            ),
            ScenarioSpec(
                name="KNOB2_PANC1_CONTROL",
                description="PANC-1 high Knob 2 sensitivity should show stronger TGF-beta brake.",
                criterion="tumor burden lower than KNOB2_ASPC1_CONTROL",
                intervention_payload={"calibration_profile": "PANC-1", "knob_interventions": []},
                user_parameter_overrides=no_cytotoxic_user,
                variable_overrides={**no_cytotoxic_var, **high_tgfb_var},
                evaluator=self._eval_knob2_panc1_control,
            ),
        ]
        return specs

    @staticmethod
    def _eval_baseline_no_cytotoxic(
        metrics: SimulationMetrics, extras: Dict[str, Any], baseline: Dict[str, Any]
    ) -> Tuple[bool, str]:
        ok = metrics.total_tumor_cells > 0
        return ok, (
            f"tumor={metrics.total_tumor_cells}, "
            f"live_tumor={metrics.live_tumor_cells}, "
            f"fitness={compute_fitness(metrics):.3f}"
        )

    @staticmethod
    def _eval_shh_off_no_cytotoxic(
        metrics: SimulationMetrics, extras: Dict[str, Any], baseline: Dict[str, Any]
    ) -> Tuple[bool, str]:
        ref = baseline.get("BASELINE_NO_CYTOTOXIC", {}).get("metrics", {})
        ref_fit = baseline.get("BASELINE_NO_CYTOTOXIC", {}).get("fitness", math.nan)
        if not ref:
            return False, "missing BASELINE_NO_CYTOTOXIC reference"

        ref_tumor = float(ref.get("total_tumor_cells", math.nan))
        cur_fit = float(compute_fitness(metrics))
        worse_burden = metrics.total_tumor_cells > ref_tumor
        worse_fitness = math.isfinite(float(ref_fit)) and (cur_fit < float(ref_fit))
        ok = worse_burden or worse_fitness
        details = (
            f"tumor={metrics.total_tumor_cells} vs baseline={ref_tumor:.1f}, "
            f"fitness={cur_fit:.3f} vs baseline={float(ref_fit):.3f}"
        )
        return ok, details

    @staticmethod
    def _eval_shh_off_with_cytotoxic(
        metrics: SimulationMetrics, extras: Dict[str, Any], baseline: Dict[str, Any]
    ) -> Tuple[bool, str]:
        ref = baseline.get("SHH_OFF_NO_CYTOTOXIC", {}).get("metrics", {})
        ref_fit = baseline.get("SHH_OFF_NO_CYTOTOXIC", {}).get("fitness", math.nan)
        if not ref:
            return False, "missing SHH_OFF_NO_CYTOTOXIC reference"

        ref_live = float(ref.get("live_tumor_cells", math.nan))
        ref_tumor = float(ref.get("total_tumor_cells", math.nan))
        cur_fit = float(compute_fitness(metrics))
        better_kill = metrics.live_tumor_cells < ref_live or metrics.total_tumor_cells < ref_tumor
        better_fitness = math.isfinite(float(ref_fit)) and (cur_fit > float(ref_fit))
        ok = better_kill and better_fitness
        details = (
            f"live_tumor={metrics.live_tumor_cells} vs shh_off_no_cyto={ref_live:.1f}, "
            f"tumor={metrics.total_tumor_cells} vs shh_off_no_cyto={ref_tumor:.1f}, "
            f"fitness={cur_fit:.3f} vs shh_off_no_cyto={float(ref_fit):.3f}"
        )
        return ok, details

    @staticmethod
    def _eval_knob2_aspc1_control(
        metrics: SimulationMetrics, extras: Dict[str, Any], baseline: Dict[str, Any]
    ) -> Tuple[bool, str]:
        ok = metrics.total_tumor_cells > 0
        return ok, f"aspc1_high_tgfb_total_tumor={metrics.total_tumor_cells}"

    @staticmethod
    def _eval_knob2_panc1_control(
        metrics: SimulationMetrics, extras: Dict[str, Any], baseline: Dict[str, Any]
    ) -> Tuple[bool, str]:
        aspc1 = baseline.get("KNOB2_ASPC1_CONTROL", {}).get("metrics", {})
        if not aspc1:
            return False, "missing KNOB2_ASPC1_CONTROL reference"
        aspc1_tumor = float(aspc1.get("total_tumor_cells", math.nan))
        ok = metrics.total_tumor_cells < aspc1_tumor
        details = f"panc1_total_tumor={metrics.total_tumor_cells} vs aspc1_total_tumor={aspc1_tumor:.1f}"
        return ok, details


def np_rint(values):
    # Local helper to avoid importing numpy at module import-time in minimal environments.
    import numpy as np

    return np.rint(values).astype(int)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run biological validation scenarios for Stroma World.")
    parser.add_argument("--physicell-binary", required=True, help="Path to stroma_world executable.")
    parser.add_argument("--physicell-config", "--base-config", dest="physicell_config", required=True, help="Path to PhysiCell_settings.xml.")
    parser.add_argument("--output-dir", required=True, help="Directory for validation outputs.")
    parser.add_argument("--timeout-seconds", type=int, default=7200, help="Per-simulation timeout in seconds.")
    parser.add_argument(
        "--sim-max-time",
        type=float,
        default=None,
        help="Override max_time (minutes) in the PhysiCell config for each validation run. "
             "Use a short value (e.g. 1440) to limit output size.",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Logging level.",
    )
    return parser.parse_args()


def configure_logging(level: str) -> None:
    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )


def main() -> int:
    args = parse_args()
    configure_logging(args.log_level)

    validator = BiologyValidator(
        binary_path=Path(args.physicell_binary),
        config_path=Path(args.physicell_config),
        output_dir=Path(args.output_dir),
        timeout_seconds=args.timeout_seconds,
        sim_max_time=args.sim_max_time,
    )
    summary = validator.run_all()

    prefix = Path(args.output_dir).expanduser().resolve() / "biology_validation_summary"
    json_path, csv_path = validator.save_summary(summary, prefix)

    print(f"Biology validation complete: {summary.passed}/{summary.total} passed")
    print(f"JSON summary: {json_path}")
    print(f"CSV summary: {csv_path}")
    for scenario in summary.scenario_results:
        status = "PASS" if scenario.success else "FAIL"
        print(f"- [{status}] {scenario.scenario}: {scenario.details}")
    return 0 if summary.failed == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
