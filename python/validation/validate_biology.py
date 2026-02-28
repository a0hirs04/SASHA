from __future__ import annotations

import argparse
import csv
import json
import logging
import math
import statistics
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
    anchor_id: int
    scenario: str
    success: bool
    criterion: str
    details: str
    dependencies: List[str]
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
    tests_total: int
    tests_passed: int
    tests_failed: int
    checks_total: int
    checks_passed: int
    checks_failed: int
    median_wall_time_per_run_seconds: float
    estimated_full_validation_wall_time_seconds: float
    scenario_results: List[ScenarioResult]


@dataclass
class ScenarioSpec:
    anchor_id: int
    name: str
    description: str
    criterion: str
    intervention_payload: Dict[str, Any]
    user_parameter_overrides: Dict[str, Any]
    variable_overrides: Dict[str, Dict[str, Any]]
    dependencies: List[str]
    evaluator: Callable[[SimulationMetrics, Dict[str, Any], Dict[str, Any]], Tuple[bool, str]]
    virtual: bool = False  # True = cache-only; no PhysiCell simulation is launched


class BiologyValidator:
    def __init__(
        self,
        binary_path: Path | str,
        config_path: Path | str,
        output_dir: Path | str,
        timeout_seconds: int = 7200,
        sim_max_time: Optional[float] = None,
        replicates_per_arm: int = 3,
        random_seed_base: int = 1001,
    ):
        self.binary_path = Path(binary_path).expanduser().resolve()
        self.base_config_path = Path(config_path).expanduser().resolve()
        self.output_dir = Path(output_dir).expanduser().resolve()
        self.timeout_seconds = int(timeout_seconds)
        self.sim_max_time: Optional[float] = sim_max_time
        self.replicates_per_arm = max(1, int(replicates_per_arm))
        self.random_seed_base = int(random_seed_base)

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
        scenario_by_name: Dict[str, ScenarioResult] = {}
        for spec in specs:
            failed_dependencies = [
                dep
                for dep in spec.dependencies
                if dep not in scenario_by_name or not scenario_by_name[dep].success
            ]
            if failed_dependencies:
                dep_msg = ", ".join(failed_dependencies)
                LOGGER.warning(
                    "Skipping validation scenario %s (anchor %d): dependency failure: %s",
                    spec.name,
                    spec.anchor_id,
                    dep_msg,
                )
                skipped = ScenarioResult(
                    anchor_id=spec.anchor_id,
                    scenario=spec.name,
                    success=False,
                    criterion=spec.criterion,
                    details=f"skipped due to failed dependencies: {dep_msg}",
                    dependencies=list(spec.dependencies),
                    run_dir="",
                    exit_code=-2,
                    wall_time=0.0,
                    fitness=0.0,
                    metrics={},
                    extras={"skipped": True, "failed_dependencies": failed_dependencies},
                )
                scenario_results.append(skipped)
                scenario_by_name[spec.name] = skipped
                self._baseline_cache[spec.name] = {
                    "fitness": skipped.fitness,
                    "metrics": skipped.metrics,
                    "extras": skipped.extras,
                }
                continue

            LOGGER.info("Running validation scenario: %s", spec.name)
            if spec.virtual:
                result = self._run_virtual_check(spec)
            else:
                result = self._run_scenario(spec)
            scenario_results.append(result)
            scenario_by_name[spec.name] = result
            self._baseline_cache[spec.name] = {
                "fitness": result.fitness,
                "metrics": result.metrics,
                "extras": result.extras,
            }

        passed = sum(1 for r in scenario_results if r.success)
        failed = len(scenario_results) - passed
        tests_total = 0
        tests_passed = 0
        checks_total = 0
        checks_passed = 0
        per_replica_wall: List[float] = []
        for result in scenario_results:
            anchor_tests = result.extras.get("anchor_tests", [])
            if isinstance(anchor_tests, list):
                tests_total += len(anchor_tests)
                tests_passed += sum(1 for t in anchor_tests if bool(t.get("passed", False)))
            check_tests = result.extras.get("check_tests", [])
            if isinstance(check_tests, list):
                checks_total += len(check_tests)
                checks_passed += sum(1 for t in check_tests if bool(t.get("passed", False)))
            wall_values = result.extras.get("replicate_wall_times", [])
            if isinstance(wall_values, list):
                for w in wall_values:
                    if isinstance(w, (int, float)) and math.isfinite(float(w)) and float(w) >= 0.0:
                        per_replica_wall.append(float(w))
        tests_failed = tests_total - tests_passed
        checks_failed = checks_total - checks_passed
        median_wall = float(statistics.median(per_replica_wall)) if per_replica_wall else math.nan
        estimated_full_wall = (median_wall * 60.0) if math.isfinite(median_wall) else math.nan
        skipped_count = sum(1 for r in scenario_results if r.exit_code == -2)
        # Anchor directional test count is stable at 34; Reality Check tests are counted separately.
        if skipped_count == 0 and tests_total != 34:
            LOGGER.warning(
                "Directional test matrix size is %d (expected 34). Check anchor test wiring.",
                tests_total,
            )
        return ValidationSummary(
            total=len(scenario_results),
            passed=passed,
            failed=failed,
            tests_total=tests_total,
            tests_passed=tests_passed,
            tests_failed=tests_failed,
            checks_total=checks_total,
            checks_passed=checks_passed,
            checks_failed=checks_failed,
            median_wall_time_per_run_seconds=median_wall,
            estimated_full_validation_wall_time_seconds=estimated_full_wall,
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
            "tests_total": summary.tests_total,
            "tests_passed": summary.tests_passed,
            "tests_failed": summary.tests_failed,
            "checks_total": summary.checks_total,
            "checks_passed": summary.checks_passed,
            "checks_failed": summary.checks_failed,
            "median_wall_time_per_run_seconds": summary.median_wall_time_per_run_seconds,
            "estimated_full_validation_wall_time_seconds": summary.estimated_full_validation_wall_time_seconds,
            "replicates_per_arm": self.replicates_per_arm,
            "scenario_results": [asdict(r) for r in summary.scenario_results],
        }
        json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

        fieldnames = [
            "anchor_id",
            "scenario",
            "success",
            "criterion",
            "details",
            "dependencies",
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
                        "anchor_id": result.anchor_id,
                        "scenario": result.scenario,
                        "success": result.success,
                        "criterion": result.criterion,
                        "details": result.details,
                        "dependencies": ";".join(result.dependencies),
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

        metrics_replicates: List[SimulationMetrics] = []
        extras_replicates: List[Dict[str, Any]] = []
        fitness_replicates: List[float] = []
        run_dirs: List[str] = []
        wall_times: List[float] = []

        for rep_idx in range(self.replicates_per_arm):
            rep_dir = scenario_dir / f"replicate_{rep_idx + 1:02d}"
            rep_dir.mkdir(parents=True, exist_ok=True)
            config_copy = rep_dir / "scenario_config.xml"
            shutil.copy2(self.base_config_path, config_copy)

            user_overrides = dict(spec.user_parameter_overrides)
            user_overrides["random_seed"] = self.random_seed_base + spec.anchor_id * 100 + rep_idx
            self._patch_config(
                config_copy,
                user_parameter_overrides=user_overrides,
                variable_overrides=spec.variable_overrides,
            )

            runner = PhysiCellRunner(
                binary_path=self.binary_path,
                config_path=config_copy,
                output_dir=rep_dir / "sim_runs",
                timeout_seconds=self.timeout_seconds,
            )

            sim_result = runner.run(spec.intervention_payload)
            if not sim_result.success:
                return ScenarioResult(
                    anchor_id=spec.anchor_id,
                    scenario=spec.name,
                    success=False,
                    criterion=spec.criterion,
                    details=(
                        f"simulation failed on replicate {rep_idx + 1}/{self.replicates_per_arm} "
                        f"(exit_code={sim_result.exit_code})"
                    ),
                    dependencies=list(spec.dependencies),
                    run_dir=str(sim_result.run_dir),
                    exit_code=sim_result.exit_code,
                    wall_time=sim_result.wall_time,
                    fitness=0.0,
                    metrics={},
                    extras={"replicate_index": rep_idx + 1},
                )

            parser = OutputParser(sim_result.run_dir / "output")
            metrics = parser.parse_final_state()
            fitness = float(compute_fitness(metrics))
            extras = self._compute_extras(parser, metrics)

            metrics_replicates.append(metrics)
            fitness_replicates.append(fitness)
            extras_replicates.append(extras)
            run_dirs.append(str(sim_result.run_dir))
            wall_times.append(float(sim_result.wall_time))

        agg_metrics = self._median_metrics(metrics_replicates)
        agg_fitness = self._median_scalar(fitness_replicates)
        agg_extras = self._median_extras(extras_replicates)
        agg_extras["replicate_count"] = int(self.replicates_per_arm)
        agg_extras["replicate_run_dirs"] = run_dirs
        agg_extras["replicate_wall_times"] = wall_times
        agg_extras["median_wall_time"] = self._median_scalar(wall_times)

        ok, details = spec.evaluator(agg_metrics, agg_extras, self._baseline_cache)
        return ScenarioResult(
            anchor_id=spec.anchor_id,
            scenario=spec.name,
            success=ok,
            criterion=spec.criterion,
            details=details,
            dependencies=list(spec.dependencies),
            run_dir=str(scenario_dir),
            exit_code=0,
            wall_time=float(sum(wall_times)),
            fitness=agg_fitness,
            metrics=asdict(agg_metrics),
            extras=agg_extras,
        )

    def _run_virtual_check(self, spec: ScenarioSpec) -> ScenarioResult:
        """Run a Reality Check without launching a PhysiCell simulation.

        The evaluator receives empty metrics/extras and reads purely from
        self._baseline_cache (populated by earlier anchor scenarios).
        """
        extras: Dict[str, Any] = {}
        try:
            ok, details = spec.evaluator(None, extras, self._baseline_cache)
        except Exception as exc:  # pragma: no cover
            ok = False
            details = f"check evaluator raised {type(exc).__name__}: {exc}"
        return ScenarioResult(
            anchor_id=spec.anchor_id,
            scenario=spec.name,
            success=ok,
            criterion=spec.criterion,
            details=details,
            dependencies=list(spec.dependencies),
            run_dir="",
            exit_code=0 if ok else 1,
            wall_time=0.0,
            fitness=0.0,
            metrics={},
            extras=extras,
        )

    @staticmethod
    def _median_scalar(values: List[float]) -> float:
        finite = [float(v) for v in values if isinstance(v, (int, float)) and math.isfinite(float(v))]
        if not finite:
            return math.nan
        return float(statistics.median(finite))

    def _median_metrics(self, metrics_list: List[SimulationMetrics]) -> SimulationMetrics:
        if not metrics_list:
            raise ValueError("cannot aggregate empty metrics list")

        fields = SimulationMetrics.__dataclass_fields__
        aggregated: Dict[str, Any] = {}
        for name in fields:
            values = [getattr(m, name) for m in metrics_list]
            if name == "source_file":
                aggregated[name] = values[-1]
                continue
            if name in {"total_tumor_cells", "total_stromal_cells", "live_tumor_cells", "activated_cafs", "mesenchymal_tumor_cells"}:
                aggregated[name] = int(round(self._median_scalar([float(v) for v in values])))
                continue
            aggregated[name] = self._median_scalar([float(v) for v in values])
        return SimulationMetrics(**aggregated)

    def _median_extras(self, extras_list: List[Dict[str, Any]]) -> Dict[str, Any]:
        if not extras_list:
            return {}

        keys: set[str] = set()
        for extra in extras_list:
            keys.update(extra.keys())

        out: Dict[str, Any] = {}
        for key in keys:
            values = [extra.get(key) for extra in extras_list if key in extra]
            if not values:
                continue
            if all(isinstance(v, bool) for v in values):
                out[key] = bool(round(self._median_scalar([1.0 if v else 0.0 for v in values])))
                continue
            if all(isinstance(v, (int, float)) for v in values):
                out[key] = self._median_scalar([float(v) for v in values])
                continue
            # Keep latest non-scalar payloads (lists/maps) for debug context.
            out[key] = values[-1]
        return out

    def _compute_extras(self, parser: OutputParser, metrics: SimulationMetrics) -> Dict[str, Any]:
        extras: Dict[str, Any] = {}

        try:
            df = parser.parse_timeseries()
            extras["timeseries_rows"] = int(len(df))
            if len(df) >= 2:
                extras["ecm_delta"] = float(df["mean_ecm"].iloc[-1] - df["mean_ecm"].iloc[0])
                extras["activated_cafs_delta"] = float(df["activated_cafs"].iloc[-1] - df["activated_cafs"].iloc[0])
                extras["tumor_count_delta"] = float(df["tumor_count"].iloc[-1] - df["tumor_count"].iloc[0])
                extras["tumor_count_min"] = float(df["tumor_count"].min())
                extras["tumor_count_end"] = float(df["tumor_count"].iloc[-1])
                extras["tumor_regrowth_after_nadir"] = bool(
                    extras["tumor_count_end"] > extras["tumor_count_min"]
                )
                nadir_idx = int(df["tumor_count"].idxmin())
                extras["mean_ecm_at_nadir"] = float(df["mean_ecm"].iloc[nadir_idx]) if "mean_ecm" in df.columns else math.nan
                extras["activated_cafs_at_nadir"] = (
                    float(df["activated_cafs"].iloc[nadir_idx]) if "activated_cafs" in df.columns else math.nan
                )
                extras.update(self._timeseries_snapshots(df))
                if "activated_cafs" in df.columns and "time" in df.columns:
                    first_on = df[df["activated_cafs"] > 0]
                    if len(first_on) > 0:
                        extras["time_to_first_caf_activation"] = float(first_on["time"].iloc[0] - df["time"].iloc[0])
                        extras["total_sim_time"] = float(df["time"].iloc[-1] - df["time"].iloc[0])
                    else:
                        extras["time_to_first_caf_activation"] = math.nan
                        extras["total_sim_time"] = float(df["time"].iloc[-1] - df["time"].iloc[0])
            else:
                extras["ecm_delta"] = math.nan
                extras["activated_cafs_delta"] = math.nan
                extras["tumor_count_delta"] = math.nan
                extras["tumor_count_min"] = math.nan
                extras["tumor_count_end"] = math.nan
                extras["tumor_regrowth_after_nadir"] = False
                extras["time_to_first_caf_activation"] = math.nan
                extras["total_sim_time"] = math.nan
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
        extras.update(self._compute_spatial_emt_metrics(parser))
        extras.update(self._compute_sanctuary_metrics(parser))
        extras.update(self._compute_snapshot_profiles(parser))
        return extras

    @staticmethod
    def _timeseries_snapshots(df) -> Dict[str, float]:
        out: Dict[str, float] = {}
        if len(df) == 0:
            return out

        idx_map = {
            "early": 0,
            "mid": len(df) // 2,
            "late": len(df) - 1,
        }
        columns = [
            "mean_ecm",
            "drug_penetration",
            "hypoxic_fraction",
            "tumor_count",
            "tumor_extent",
            "mesenchymal_tumor_cells",
            "activated_cafs",
        ]
        for col in columns:
            if col not in df.columns:
                continue
            for stage, idx in idx_map.items():
                out[f"{col}_{stage}"] = float(df[col].iloc[idx])
        return out

    def _compute_spatial_emt_metrics(self, parser: OutputParser) -> Dict[str, Any]:
        result: Dict[str, Any] = {
            "core_mesenchymal_fraction": math.nan,
            "peripheral_mesenchymal_fraction": math.nan,
            "emt_peripheral_minus_core": math.nan,
            "core_cdh1_mean": math.nan,
            "peripheral_cdh1_mean": math.nan,
        }
        try:
            final_xml = parser._find_final_snapshot_xml()  # noqa: SLF001
            snapshot = parser._read_physicell_xml(final_xml)  # noqa: SLF001
            matrix = snapshot["cell_matrix"]
            label_name_map = snapshot["label_name_map"]
        except Exception:  # pragma: no cover
            return result

        cell_type = self._row_by_label(matrix, label_name_map, "cell_type")
        dead = self._row_by_label(matrix, label_name_map, "dead")
        death_model = self._row_by_label(matrix, label_name_map, "current_death_model")
        is_mesenchymal = self._row_by_label(matrix, label_name_map, "is_mesenchymal")
        cdh1 = self._row_by_label(matrix, label_name_map, "CDH1")
        positions = self._positions_by_label(matrix, label_name_map, "position")

        if cell_type is None or is_mesenchymal is None or positions.size == 0:
            return result

        import numpy as np

        tumor_mask = np_rint(cell_type) == 0
        if dead is not None:
            tumor_mask &= dead <= 0.5
        if death_model is not None:
            tumor_mask &= np_rint(death_model) != 100

        if int(np.sum(tumor_mask)) < 10:
            return result

        tumor_positions = positions[tumor_mask, :]
        tumor_mesenchymal = is_mesenchymal[tumor_mask] > 0.5
        centroid = np.nanmean(tumor_positions, axis=0)
        radial = np.linalg.norm(tumor_positions - centroid, axis=1)
        core_cutoff = float(np.quantile(radial, 0.40))
        peripheral_cutoff = float(np.quantile(radial, 0.60))
        core_mask = radial <= core_cutoff
        peripheral_mask = radial >= peripheral_cutoff

        if not core_mask.any() or not peripheral_mask.any():
            return result

        core_frac = float(np.mean(tumor_mesenchymal[core_mask]))
        peripheral_frac = float(np.mean(tumor_mesenchymal[peripheral_mask]))
        result["core_mesenchymal_fraction"] = core_frac
        result["peripheral_mesenchymal_fraction"] = peripheral_frac
        result["emt_peripheral_minus_core"] = float(peripheral_frac - core_frac)
        if cdh1 is not None:
            tumor_cdh1 = cdh1[tumor_mask]
            result["core_cdh1_mean"] = float(np.mean(tumor_cdh1[core_mask]))
            result["peripheral_cdh1_mean"] = float(np.mean(tumor_cdh1[peripheral_mask]))
        return result

    def _compute_sanctuary_metrics(self, parser: OutputParser) -> Dict[str, Any]:
        result: Dict[str, Any] = {
            "live_mean_local_ecm": math.nan,
            "all_mean_local_ecm": math.nan,
            "live_mean_local_drug": math.nan,
            "all_mean_local_drug": math.nan,
            "live_fraction_in_top_ecm_quartile": math.nan,
            "live_tumor_count_snapshot": 0,
            "dead_mean_local_ecm": math.nan,
            "dead_mean_local_drug": math.nan,
            "mean_distance_survivor_from_boundary": math.nan,
            "mean_distance_dead_from_boundary": math.nan,
        }
        try:
            final_xml = parser._find_final_snapshot_xml()  # noqa: SLF001
            snapshot = parser._read_physicell_xml(final_xml)  # noqa: SLF001
            matrix = snapshot["cell_matrix"]
            label_name_map = snapshot["label_name_map"]
            micro_coords = snapshot["micro_coords"]
            micro_values = snapshot["micro_values"]
        except Exception:  # pragma: no cover
            return result

        cell_type = self._row_by_label(matrix, label_name_map, "cell_type")
        dead = self._row_by_label(matrix, label_name_map, "dead")
        death_model = self._row_by_label(matrix, label_name_map, "current_death_model")
        positions = self._positions_by_label(matrix, label_name_map, "position")
        ecm_vals = micro_values.get("ecm_density")
        drug_vals = micro_values.get("drug")

        if cell_type is None or positions.size == 0:
            return result
        if ecm_vals is None or drug_vals is None or micro_coords.size == 0:
            return result

        import numpy as np

        tumor_mask = np_rint(cell_type) == 0
        live_tumor_mask = tumor_mask.copy()
        if dead is not None:
            live_tumor_mask &= dead <= 0.5
        if death_model is not None:
            live_tumor_mask &= np_rint(death_model) != 100
        dead_tumor_mask = tumor_mask & (~live_tumor_mask)

        if not tumor_mask.any():
            return result

        tumor_positions = positions[tumor_mask, :]
        live_positions = positions[live_tumor_mask, :]
        result["live_tumor_count_snapshot"] = int(np.sum(live_tumor_mask))

        tumor_ecm = self._sample_nearest_voxel_values(tumor_positions, micro_coords, ecm_vals)
        tumor_drug = self._sample_nearest_voxel_values(tumor_positions, micro_coords, drug_vals)
        result["all_mean_local_ecm"] = float(np.nanmean(tumor_ecm))
        result["all_mean_local_drug"] = float(np.nanmean(tumor_drug))

        if live_positions.size == 0:
            return result

        live_ecm = self._sample_nearest_voxel_values(live_positions, micro_coords, ecm_vals)
        live_drug = self._sample_nearest_voxel_values(live_positions, micro_coords, drug_vals)
        result["live_mean_local_ecm"] = float(np.nanmean(live_ecm))
        result["live_mean_local_drug"] = float(np.nanmean(live_drug))
        top_ecm_cutoff = float(np.quantile(tumor_ecm, 0.75))
        result["live_fraction_in_top_ecm_quartile"] = float(np.mean(live_ecm >= top_ecm_cutoff))

        if dead_tumor_mask.any():
            dead_positions = positions[dead_tumor_mask, :]
            dead_ecm = self._sample_nearest_voxel_values(dead_positions, micro_coords, ecm_vals)
            dead_drug = self._sample_nearest_voxel_values(dead_positions, micro_coords, drug_vals)
            result["dead_mean_local_ecm"] = float(np.nanmean(dead_ecm))
            result["dead_mean_local_drug"] = float(np.nanmean(dead_drug))

            d_survivor = self._distance_to_domain_boundary(live_positions, micro_coords)
            d_dead = self._distance_to_domain_boundary(dead_positions, micro_coords)
            result["mean_distance_survivor_from_boundary"] = float(np.nanmean(d_survivor))
            result["mean_distance_dead_from_boundary"] = float(np.nanmean(d_dead))
        return result

    def _compute_snapshot_profiles(self, parser: OutputParser) -> Dict[str, Any]:
        out: Dict[str, Any] = {}
        try:
            snapshots = parser._list_timeseries_xmls()  # noqa: SLF001
            if not snapshots:
                return out
            idx_map = {
                "early": 0,
                "mid": len(snapshots) // 2,
                "late": len(snapshots) - 1,
            }
            for stage, idx in idx_map.items():
                snap = parser._read_physicell_xml(snapshots[idx])  # noqa: SLF001
                stage_metrics = self._snapshot_metrics(snap)
                for key, value in stage_metrics.items():
                    out[f"{key}_{stage}"] = value
            return out
        except Exception:  # pragma: no cover
            return out

    def _snapshot_metrics(self, snapshot: Dict[str, Any]) -> Dict[str, float]:
        import numpy as np

        matrix = snapshot["cell_matrix"]
        labels = snapshot["label_name_map"]
        micro_coords = snapshot["micro_coords"]
        micro_vals = snapshot["micro_values"]

        cell_type = self._row_by_label(matrix, labels, "cell_type")
        dead = self._row_by_label(matrix, labels, "dead")
        death_model = self._row_by_label(matrix, labels, "current_death_model")
        hif1a = self._row_by_label(matrix, labels, "HIF1A")
        cdh1 = self._row_by_label(matrix, labels, "CDH1")
        is_mes = self._row_by_label(matrix, labels, "is_mesenchymal")
        pos = self._positions_by_label(matrix, labels, "position")
        oxygen = micro_vals.get("oxygen")
        ecm = micro_vals.get("ecm_density")
        drug = micro_vals.get("drug")

        out = {
            "oxygen_centroid": math.nan,
            "oxygen_periphery": math.nan,
            "hif1a_core_fraction": math.nan,
            "hif1a_periphery_fraction": math.nan,
            "emt_core_fraction": math.nan,
            "emt_periphery_fraction": math.nan,
            "cdh1_core_mean": math.nan,
            "cdh1_periphery_mean": math.nan,
            "ecm_peri_tumor_shell": math.nan,
            "ecm_domain_boundary": math.nan,
            "barrier_maturity_index": math.nan,
            "drug_centroid": math.nan,
            "drug_penetration_half_depth": math.nan,
        }
        if cell_type is None or pos.size == 0:
            return out

        tumor_mask = np_rint(cell_type) == 0
        if dead is not None:
            tumor_mask &= dead <= 0.5
        if death_model is not None:
            tumor_mask &= np_rint(death_model) != 100
        if not tumor_mask.any():
            return out

        tumor_pos = pos[tumor_mask, :]
        centroid = np.nanmean(tumor_pos, axis=0)
        radial = np.linalg.norm(tumor_pos - centroid, axis=1)
        radius = float(np.max(radial)) if radial.size else 0.0
        core_mask = radial <= (0.5 * radius)
        periph_mask = radial >= (0.5 * radius)

        if oxygen is not None and micro_coords.size:
            out["oxygen_centroid"] = float(self._sample_nearest_voxel_values(centroid[None, :], micro_coords, oxygen)[0])
            if periph_mask.any():
                periph_o2 = self._sample_nearest_voxel_values(tumor_pos[periph_mask, :], micro_coords, oxygen)
                out["oxygen_periphery"] = float(np.nanmean(periph_o2))

        if hif1a is not None:
            tumor_hif = hif1a[tumor_mask] > 0.5
            if core_mask.any():
                out["hif1a_core_fraction"] = float(np.mean(tumor_hif[core_mask]))
            if periph_mask.any():
                out["hif1a_periphery_fraction"] = float(np.mean(tumor_hif[periph_mask]))

        if is_mes is not None:
            tumor_mes = is_mes[tumor_mask] > 0.5
            if core_mask.any():
                out["emt_core_fraction"] = float(np.mean(tumor_mes[core_mask]))
            if periph_mask.any():
                out["emt_periphery_fraction"] = float(np.mean(tumor_mes[periph_mask]))

        if cdh1 is not None:
            tumor_cdh1 = cdh1[tumor_mask]
            if core_mask.any():
                out["cdh1_core_mean"] = float(np.mean(tumor_cdh1[core_mask]))
            if periph_mask.any():
                out["cdh1_periphery_mean"] = float(np.mean(tumor_cdh1[periph_mask]))

        if ecm is not None and micro_coords.size:
            voxel_dist = np.linalg.norm(micro_coords - centroid, axis=1)
            shell_mask = (voxel_dist >= radius) & (voxel_dist <= radius + 50.0)
            if shell_mask.any():
                out["ecm_peri_tumor_shell"] = float(np.nanmean(ecm[shell_mask]))
                out["barrier_maturity_index"] = out["ecm_peri_tumor_shell"]

            boundary_dist = self._distance_to_domain_boundary(micro_coords, micro_coords)
            if boundary_dist.size:
                edge_cut = float(np.quantile(boundary_dist, 0.2))
                edge_mask = boundary_dist <= edge_cut
                if edge_mask.any():
                    out["ecm_domain_boundary"] = float(np.nanmean(ecm[edge_mask]))

        if drug is not None and micro_coords.size:
            out["drug_centroid"] = float(self._sample_nearest_voxel_values(centroid[None, :], micro_coords, drug)[0])
            out["drug_penetration_half_depth"] = self._drug_half_depth(micro_coords, drug)

        return out

    @staticmethod
    def _distance_to_domain_boundary(points, reference_voxels):
        import numpy as np

        pts = np.asarray(points, dtype=float)
        vox = np.asarray(reference_voxels, dtype=float)
        if pts.size == 0 or vox.size == 0:
            return np.asarray([], dtype=float)

        mins = np.min(vox, axis=0)
        maxs = np.max(vox, axis=0)
        d_to_min = pts - mins[None, :]
        d_to_max = maxs[None, :] - pts
        d = np.minimum(d_to_min, d_to_max)
        return np.min(d, axis=1)

    def _drug_half_depth(self, voxel_coords, drug_values) -> float:
        import numpy as np

        dist = self._distance_to_domain_boundary(voxel_coords, voxel_coords)
        if dist.size == 0:
            return math.nan
        edge_cut = float(np.quantile(dist, 0.2))
        boundary_mask = dist <= edge_cut
        if not boundary_mask.any():
            return math.nan
        boundary_drug = float(np.nanmean(drug_values[boundary_mask]))
        if not math.isfinite(boundary_drug) or boundary_drug <= 0.0:
            return math.nan
        half_mask = drug_values >= (0.5 * boundary_drug)
        if not half_mask.any():
            return 0.0
        return float(np.nanmax(dist[half_mask]))

    @staticmethod
    def _positions_by_label(matrix, label_name_map: Dict[str, Dict[str, Any]], label: str):
        import numpy as np

        entry = label_name_map.get(label)
        if entry is None:
            return np.empty((0, 3), dtype=float)
        idx = int(entry.get("index", -1))
        size = int(entry.get("size", 0))
        if idx < 0 or size < 3:
            return np.empty((0, 3), dtype=float)
        if idx + 3 > matrix.shape[0]:
            return np.empty((0, 3), dtype=float)
        return matrix[idx : idx + 3, :].T

    @staticmethod
    def _sample_nearest_voxel_values(points, voxels, values):
        import numpy as np

        if points.size == 0:
            return np.asarray([], dtype=float)
        d2 = np.sum((points[:, None, :] - voxels[None, :, :]) ** 2, axis=2)
        nearest_idx = np.argmin(d2, axis=1)
        return np.asarray(values[nearest_idx], dtype=float)

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

    @staticmethod
    def _make_test(test_id: str, passed: bool, criterion: str, fail_if: str) -> Dict[str, Any]:
        return {
            "id": test_id,
            "passed": bool(passed),
            "criterion": criterion,
            "fail_if": fail_if,
        }

    def _finalize_anchor_tests(
        self,
        extras: Dict[str, Any],
        tests: List[Dict[str, Any]],
        derived_notes: Optional[List[str]] = None,
    ) -> Tuple[bool, str]:
        extras["anchor_tests"] = tests
        if derived_notes:
            extras["anchor_derived_notes"] = derived_notes
        failed = [t for t in tests if not bool(t["passed"])]
        if not failed:
            return True, f"{len(tests)}/{len(tests)} directional tests passed"
        fail_ids = ", ".join(t["id"] for t in failed)
        fail_hints = " | ".join(t["fail_if"] for t in failed)
        details = (
            f"{len(tests)-len(failed)}/{len(tests)} directional tests passed; "
            f"failed: {fail_ids}; debug: {fail_hints}"
        )
        return False, details

    @staticmethod
    def _cached_metrics(baseline: Dict[str, Any], name: str) -> Optional[Dict[str, Any]]:
        m = baseline.get(name, {}).get("metrics", {})
        return m if isinstance(m, dict) and m else None

    @staticmethod
    def _cached_extras(baseline: Dict[str, Any], name: str) -> Dict[str, Any]:
        e = baseline.get(name, {}).get("extras", {})
        return e if isinstance(e, dict) else {}

    @staticmethod
    def _viability(metrics: Dict[str, Any]) -> float:
        total = float(metrics.get("total_tumor_cells", math.nan))
        live = float(metrics.get("live_tumor_cells", math.nan))
        if not math.isfinite(total) or total <= 0.0 or not math.isfinite(live):
            return math.nan
        return live / total

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
        high_stroma_user = {"number_of_stromal_cells": 400}
        shh_off = [
            {"knob": "shh_secretion_rate", "effect": "INHIBIT", "strength": 1.0, "name": "SHH_OFF"}
        ]
        ha_degrade_user  = {**cytotoxic_user, "ha_degrade_strength": 0.9}
        col_degrade_user = {**cytotoxic_user, "col_degrade_strength": 0.9}
        both_degrade_user = {**cytotoxic_user, "ha_degrade_strength": 0.9, "col_degrade_strength": 0.9}
        ecm_degrade_only_user = {**no_cytotoxic_user, "ha_degrade_strength": 0.9, "col_degrade_strength": 0.9}

        specs: List[ScenarioSpec] = [
            ScenarioSpec(
                anchor_id=1,
                name="ANCHOR_1_SELF_ASSEMBLY",
                description="Barrier self-assembles from tumor->stroma paracrine signaling.",
                criterion="CAF activation rises and ECM accumulates without external instruction",
                intervention_payload={"calibration_profile": "AsPC-1", "knob_interventions": []},
                user_parameter_overrides=no_cytotoxic_user,
                variable_overrides=no_cytotoxic_var,
                dependencies=[],
                evaluator=self._eval_anchor1_self_assembly,
            ),
            ScenarioSpec(
                anchor_id=3,
                name="ANCHOR_3_SHH_PARADOX",
                description="SHH inhibition alone should worsen tumor outcome vs baseline.",
                criterion="tumor burden and extent increase vs Anchor 1 without cytotoxic",
                intervention_payload={
                    "calibration_profile": "AsPC-1",
                    "knob_interventions": shh_off,
                },
                user_parameter_overrides=no_cytotoxic_user,
                variable_overrides=no_cytotoxic_var,
                dependencies=["ANCHOR_1_SELF_ASSEMBLY"],
                evaluator=self._eval_anchor3_shh_paradox,
            ),
            ScenarioSpec(
                anchor_id=2,
                name="ANCHOR_2_DRUG_PENETRATION_MATURITY",
                description="Drug penetration declines as ECM barrier matures over time.",
                criterion="ECM rises early->mid->late while tumor drug penetration drops monotonically",
                intervention_payload={"calibration_profile": "AsPC-1", "knob_interventions": []},
                user_parameter_overrides=cytotoxic_user,
                variable_overrides=cytotoxic_var,
                dependencies=["ANCHOR_1_SELF_ASSEMBLY"],
                evaluator=self._eval_anchor2_penetration_maturity,
            ),
            ScenarioSpec(
                anchor_id=4,
                name="ANCHOR_4_SMAD4_ASYMMETRY",
                description="SMAD4-intact profile should show stronger TGF-beta growth braking.",
                criterion="PANC-1 (high knob2) grows slower than AsPC-1 baseline under high TGF-beta",
                intervention_payload={
                    "calibration_profile": "PANC-1",
                    "knob_interventions": [],
                },
                user_parameter_overrides=no_cytotoxic_user,
                variable_overrides={**no_cytotoxic_var, **high_tgfb_var},
                dependencies=["ANCHOR_1_SELF_ASSEMBLY"],
                evaluator=self._eval_anchor4_smad4_asymmetry,
            ),
            ScenarioSpec(
                anchor_id=5,
                name="ANCHOR_5_CENTRAL_HYPOXIA",
                description="Hypoxia should emerge in the interior as barrier density increases.",
                criterion="hypoxic fraction rises with ECM maturation over time",
                intervention_payload={"calibration_profile": "AsPC-1", "knob_interventions": []},
                user_parameter_overrides=no_cytotoxic_user,
                variable_overrides=no_cytotoxic_var,
                dependencies=["ANCHOR_1_SELF_ASSEMBLY"],
                evaluator=self._eval_anchor5_central_hypoxia,
            ),
            ScenarioSpec(
                anchor_id=6,
                name="ANCHOR_6_PERIPHERAL_EMT",
                description="EMT should be spatially enriched at tumor periphery.",
                criterion="peripheral mesenchymal fraction is higher than core fraction",
                intervention_payload={"calibration_profile": "AsPC-1", "knob_interventions": []},
                user_parameter_overrides=no_cytotoxic_user,
                variable_overrides=no_cytotoxic_var,
                dependencies=["ANCHOR_5_CENTRAL_HYPOXIA"],
                evaluator=self._eval_anchor6_peripheral_emt,
            ),
            ScenarioSpec(
                anchor_id=7,
                name="ANCHOR_7_ECM_DEGRADATION_LIMITS",
                description="Barrier opening without cytotoxic is insufficient; with cytotoxic it improves kill.",
                criterion="SHH-off + cytotoxic improves kill vs SHH-off alone while residual tumor persists",
                intervention_payload={
                    "calibration_profile": "AsPC-1",
                    "knob_interventions": shh_off,
                },
                user_parameter_overrides=cytotoxic_user,
                variable_overrides=cytotoxic_var,
                dependencies=["ANCHOR_3_SHH_PARADOX", "ANCHOR_2_DRUG_PENETRATION_MATURITY"],
                evaluator=self._eval_anchor7_ecm_degradation_limits,
            ),
            ScenarioSpec(
                anchor_id=7,
                name="ANCHOR_7B_ECM_DEGRADE_ONLY",
                description="ECM depletion without cytotoxic: barrier removal alone is insufficient for tumor control.",
                criterion="data collection for C4 treatment ranking (0 directional tests)",
                intervention_payload={"calibration_profile": "AsPC-1", "knob_interventions": []},
                user_parameter_overrides=ecm_degrade_only_user,
                variable_overrides=no_cytotoxic_var,
                dependencies=["ANCHOR_7_ECM_DEGRADATION_LIMITS"],
                evaluator=self._eval_anchor7b_ecm_degrade_only,
            ),
            ScenarioSpec(
                anchor_id=8,
                name="ANCHOR_8_HA_DEPLETED_DRUG",
                description="HA-depleted arm: hyaluronidase reduces HA diffusion barrier without altering collagen.",
                criterion="data collection for three-arm independence test (0 directional tests in this arm)",
                intervention_payload={"calibration_profile": "AsPC-1", "knob_interventions": []},
                user_parameter_overrides=ha_degrade_user,
                variable_overrides=cytotoxic_var,
                dependencies=["ANCHOR_7_ECM_DEGRADATION_LIMITS"],
                evaluator=self._eval_anchor8_arm_data,
            ),
            ScenarioSpec(
                anchor_id=8,
                name="ANCHOR_8_COL_DEPLETED_DRUG",
                description="Collagen-depleted arm: collagenase reduces mechanical stiffness without altering HA.",
                criterion="data collection for three-arm independence test (0 directional tests in this arm)",
                intervention_payload={"calibration_profile": "AsPC-1", "knob_interventions": []},
                user_parameter_overrides=col_degrade_user,
                variable_overrides=cytotoxic_var,
                dependencies=["ANCHOR_7_ECM_DEGRADATION_LIMITS"],
                evaluator=self._eval_anchor8_arm_data,
            ),
            ScenarioSpec(
                anchor_id=8,
                name="ANCHOR_8_BOTH_DEPLETED_DRUG",
                description="HA + collagen depletion: independence test shows HA drives penetration and collagen drives confinement.",
                criterion="pen(HA)>pen(COL) and extent(COL)>extent(HA) and combined additivity",
                intervention_payload={"calibration_profile": "AsPC-1", "knob_interventions": []},
                user_parameter_overrides=both_degrade_user,
                variable_overrides=cytotoxic_var,
                dependencies=["ANCHOR_8_HA_DEPLETED_DRUG", "ANCHOR_8_COL_DEPLETED_DRUG"],
                evaluator=self._eval_anchor8_independence,
            ),
            ScenarioSpec(
                anchor_id=9,
                name="ANCHOR_9_BARRIER_DENSITY_PROGNOSTIC",
                description="Higher PSC density should worsen response for same drug regimen.",
                criterion="higher stromal seeding yields lower penetration and higher residual tumor",
                intervention_payload={"calibration_profile": "AsPC-1", "knob_interventions": []},
                user_parameter_overrides={**cytotoxic_user, **high_stroma_user},
                variable_overrides=cytotoxic_var,
                dependencies=["ANCHOR_2_DRUG_PENETRATION_MATURITY"],
                evaluator=self._eval_anchor9_density_prognostic,
            ),
            ScenarioSpec(
                anchor_id=10,
                name="ANCHOR_10_SPATIAL_SANCTUARY",
                description="Residual survivors should cluster in dense ECM sanctuary under strong drug pressure.",
                criterion="survivors localize to high-ECM/low-drug niches and show regrowth tendency",
                intervention_payload={"calibration_profile": "AsPC-1", "knob_interventions": []},
                user_parameter_overrides={**cytotoxic_user, **high_stroma_user},
                variable_overrides=cytotoxic_var,
                dependencies=[
                    "ANCHOR_1_SELF_ASSEMBLY",
                    "ANCHOR_2_DRUG_PENETRATION_MATURITY",
                    "ANCHOR_3_SHH_PARADOX",
                    "ANCHOR_4_SMAD4_ASYMMETRY",
                    "ANCHOR_5_CENTRAL_HYPOXIA",
                    "ANCHOR_6_PERIPHERAL_EMT",
                    "ANCHOR_7_ECM_DEGRADATION_LIMITS",
                    "ANCHOR_8_BOTH_DEPLETED_DRUG",
                    "ANCHOR_9_BARRIER_DENSITY_PROGNOSTIC",
                ],
                evaluator=self._eval_anchor10_spatial_sanctuary,
            ),
            # ----------------------------------------------------------------
            # §4 Reality Check Suite — integration gates.
            # All are virtual (no simulation run); they read from the baseline
            # cache populated by earlier anchor runs.
            # The EA does not start until all 5 Checks pass.
            # ----------------------------------------------------------------
            ScenarioSpec(
                anchor_id=11,
                name="CHECK_1_NATURAL_HISTORY",
                description="Whole-system natural history: tumor grows, stroma assembles, hypoxia and EMT emerge.",
                criterion="All four spontaneous PDAC hallmarks present simultaneously in one coherent run",
                intervention_payload={},
                user_parameter_overrides={},
                variable_overrides={},
                dependencies=["ANCHOR_1_SELF_ASSEMBLY", "ANCHOR_5_CENTRAL_HYPOXIA", "ANCHOR_6_PERIPHERAL_EMT"],
                evaluator=self._eval_check1_natural_history,
                virtual=True,
            ),
            ScenarioSpec(
                anchor_id=12,
                name="CHECK_2_CHEMO_FAILURE",
                description="Single-agent chemotherapy partially responds but fails to eradicate; barrier persists, survivors regrow, resistance emerges.",
                criterion="partial response + barrier persistence ≥80% + regrowth after withdrawal + ABCB1 emergence",
                intervention_payload={},
                user_parameter_overrides={},
                variable_overrides={},
                dependencies=[
                    "ANCHOR_1_SELF_ASSEMBLY",
                    "ANCHOR_2_DRUG_PENETRATION_MATURITY",
                    "ANCHOR_7_ECM_DEGRADATION_LIMITS",
                    "ANCHOR_10_SPATIAL_SANCTUARY",
                ],
                evaluator=self._eval_check2_chemo_failure,
                virtual=True,
            ),
            ScenarioSpec(
                anchor_id=13,
                name="CHECK_3_VISMODEGIB_PARADOX",
                description="Three-way rank: SHH+drug best, control middle, SHH-only worst.",
                criterion="tumor_count(SHH+drug) < tumor_count(control) < tumor_count(SHH-only)",
                intervention_payload={},
                user_parameter_overrides={},
                variable_overrides={},
                dependencies=["ANCHOR_1_SELF_ASSEMBLY", "ANCHOR_3_SHH_PARADOX", "ANCHOR_7_ECM_DEGRADATION_LIMITS"],
                evaluator=self._eval_check3_vismodegib_paradox,
                virtual=True,
            ),
            ScenarioSpec(
                anchor_id=14,
                name="CHECK_4_FITNESS_RANKING",
                description="Treatment hierarchy: ECM-degrade+drug > drug-only > ECM-degrade-only > SHH-only.",
                criterion="tc(A8_BOTH) < tc(A2) < tc(ECM_DEGRADE_ONLY) < tc(A3) — four pairwise comparisons",
                intervention_payload={},
                user_parameter_overrides={},
                variable_overrides={},
                dependencies=[
                    "ANCHOR_2_DRUG_PENETRATION_MATURITY",
                    "ANCHOR_3_SHH_PARADOX",
                    "ANCHOR_7B_ECM_DEGRADE_ONLY",
                    "ANCHOR_8_BOTH_DEPLETED_DRUG",
                ],
                evaluator=self._eval_check4_fitness_ranking,
                virtual=True,
            ),
            ScenarioSpec(
                anchor_id=15,
                name="CHECK_5_SANCTUARY_REGROWTH",
                description="Sanctuary survivors regrow after drug withdrawal; barrier re-assembles.",
                criterion="All four sanctuary-regrowth criteria from Anchor 10 satisfied concurrently",
                intervention_payload={},
                user_parameter_overrides={},
                variable_overrides={},
                dependencies=[
                    "CHECK_1_NATURAL_HISTORY",
                    "CHECK_2_CHEMO_FAILURE",
                    "CHECK_3_VISMODEGIB_PARADOX",
                    "CHECK_4_FITNESS_RANKING",
                    "ANCHOR_10_SPATIAL_SANCTUARY",
                ],
                evaluator=self._eval_check5_sanctuary_regrowth,
                virtual=True,
            ),
        ]
        return specs

    def _eval_anchor1_self_assembly(
        self, metrics: SimulationMetrics, extras: Dict[str, Any], baseline: Dict[str, Any]
    ) -> Tuple[bool, str]:
        tests = [
            self._make_test(
                "P1a",
                float(extras.get("activated_caf_fraction", math.nan)) > 0.0,
                "CAF activation fraction at endpoint > CAF activation fraction at t=0",
                "CAF activation loop broken: check tumor secretion -> PSC receptor activation -> ACTA2 transition.",
            ),
            self._make_test(
                "P1b",
                float(extras.get("ecm_peri_tumor_shell_late", math.nan))
                > float(extras.get("ecm_domain_boundary_late", math.nan)),
                "ECM near tumor > ECM at domain boundary",
                "ECM gradient missing/inverted: check paracrine diffusion and ACTA2 -> HAS2/COL1A1 linkage.",
            ),
            self._make_test(
                "P1c",
                float(extras.get("time_to_first_caf_activation", math.nan))
                < float(extras.get("total_sim_time", math.nan)),
                "Time to first CAF activation < total simulation duration",
                "Activation only at simulation end: thresholds/marginal triggering likely miswired.",
            ),
        ]
        return self._finalize_anchor_tests(extras, tests)

    def _eval_anchor3_shh_paradox(
        self, metrics: SimulationMetrics, extras: Dict[str, Any], baseline: Dict[str, Any]
    ) -> Tuple[bool, str]:
        ref = self._cached_metrics(baseline, "ANCHOR_1_SELF_ASSEMBLY")
        if ref is None:
            return False, "missing ANCHOR_1_SELF_ASSEMBLY reference"

        tests = [
            self._make_test(
                "P3a",
                float(metrics.total_tumor_cells) > float(ref.get("total_tumor_cells", math.nan)),
                "Tumor cell count (SHH-blocked) > Tumor cell count (Control)",
                "SHH paradox absent: check stroma-as-restraint feedback and SHH wiring.",
            ),
            self._make_test(
                "P3b",
                float(metrics.tumor_extent) > float(ref.get("tumor_extent", math.nan)),
                "Tumor spatial extent (SHH-blocked) > Tumor spatial extent (Control)",
                "Mechanical confinement too weak: inspect compaction/contact-inhibition coupling.",
            ),
            self._make_test(
                "P3c",
                float(metrics.mean_ecm_density) < float(ref.get("mean_ecm_density", math.nan)),
                "Barrier density (SHH-blocked) < Barrier density (Control)",
                "SHH->stroma arm ineffective: verify GLI1/HAS2/COL1A1 response and diffusion.",
            ),
        ]
        return self._finalize_anchor_tests(extras, tests)

    def _eval_anchor2_penetration_maturity(
        self, metrics: SimulationMetrics, extras: Dict[str, Any], baseline: Dict[str, Any]
    ) -> Tuple[bool, str]:
        tests = [
            self._make_test(
                "P2a",
                float(extras.get("drug_centroid_late", math.nan))
                < float(extras.get("drug_centroid_mid", math.nan))
                < float(extras.get("drug_centroid_early", math.nan)),
                "Drug at centroid late < mid < early",
                "Centroid drug not decreasing with maturity: check ECM->diffusion coupling.",
            ),
            self._make_test(
                "P2b",
                float(extras.get("drug_penetration_half_depth_late", math.nan))
                < float(extras.get("drug_penetration_half_depth_early", math.nan)),
                "Drug penetration half-depth late < early",
                "Drug front not retreating: diffusion barrier physics likely disconnected.",
            ),
            self._make_test(
                "P2c",
                float(extras.get("barrier_maturity_index_early", math.nan))
                < float(extras.get("barrier_maturity_index_mid", math.nan))
                < float(extras.get("barrier_maturity_index_late", math.nan)),
                "Barrier maturity index early < mid < late",
                "Barrier not maturing across snapshots: anchor preconditions invalid.",
            ),
        ]
        return self._finalize_anchor_tests(extras, tests)

    def _eval_anchor4_smad4_asymmetry(
        self, metrics: SimulationMetrics, extras: Dict[str, Any], baseline: Dict[str, Any]
    ) -> Tuple[bool, str]:
        aspc1 = self._cached_metrics(baseline, "ANCHOR_1_SELF_ASSEMBLY")
        if aspc1 is None:
            return False, "missing ANCHOR_1_SELF_ASSEMBLY reference"

        tests = [
            self._make_test(
                "P4a",
                float(metrics.total_tumor_cells) < float(aspc1.get("total_tumor_cells", math.nan)),
                "Tumor count (SMAD4-intact) < Tumor count (SMAD4-null)",
                "SMAD4 asymmetry broken: TGF-beta brake not engaging in intact context.",
            ),
            self._make_test(
                "P4b",
                float(extras.get("tumor_count_delta", math.nan))
                < float(baseline.get("ANCHOR_1_SELF_ASSEMBLY", {}).get("extras", {}).get("tumor_count_delta", math.nan)),
                "Tumor proliferation proxy (delta count) SMAD4-intact < SMAD4-null",
                "Growth brake not mechanistically visible: inspect Knob 2 -> RB1/TGF arm wiring.",
            ),
            self._make_test(
                "P4c",
                float(metrics.mean_ecm_density) <= float(aspc1.get("mean_ecm_density", math.nan)),
                "Barrier density (SMAD4-intact) <= Barrier density (SMAD4-null)",
                "Unexpected barrier increase in SMAD4-intact branch: inspect feedback coupling.",
            ),
        ]
        return self._finalize_anchor_tests(extras, tests)

    def _eval_anchor5_central_hypoxia(
        self, metrics: SimulationMetrics, extras: Dict[str, Any], baseline: Dict[str, Any]
    ) -> Tuple[bool, str]:
        tests = [
            self._make_test(
                "P5a",
                float(extras.get("oxygen_centroid_late", math.nan))
                < float(extras.get("oxygen_centroid_early", math.nan)),
                "O2 at centroid late < O2 at centroid early",
                "Central oxygen not declining: inspect ECM->oxygen diffusion coupling.",
            ),
            self._make_test(
                "P5b",
                float(extras.get("oxygen_centroid_late", math.nan))
                < float(extras.get("oxygen_periphery_late", math.nan)),
                "O2 at centroid late < O2 at periphery late",
                "No center-periphery oxygen gradient: check domain geometry/boundary conditions.",
            ),
            self._make_test(
                "P5c",
                float(extras.get("hif1a_core_fraction_late", math.nan))
                > float(extras.get("hif1a_periphery_fraction_late", math.nan)),
                "HIF1A fraction core > periphery at late stage",
                "HIF1A gradient inverted/flat: inspect oxygen sensing thresholds and update chain.",
            ),
        ]
        return self._finalize_anchor_tests(extras, tests)

    def _eval_anchor6_peripheral_emt(
        self, metrics: SimulationMetrics, extras: Dict[str, Any], baseline: Dict[str, Any]
    ) -> Tuple[bool, str]:
        tests = [
            self._make_test(
                "P6a",
                float(extras.get("emt_periphery_fraction_late", math.nan))
                > float(extras.get("emt_core_fraction_late", math.nan)),
                "EMT fraction periphery > core",
                "Spatial EMT restriction failed: check emt induction threshold and TGF-beta gradient.",
            ),
            self._make_test(
                "P6b",
                float(extras.get("cdh1_core_mean_late", math.nan))
                > float(extras.get("emt_core_fraction_late", math.nan)),
                "Core epithelial tendency > core mesenchymal tendency",
                "Whole-tumor EMT drift: core no longer predominantly epithelial.",
            ),
            self._make_test(
                "P6c",
                float(extras.get("cdh1_core_mean_late", math.nan))
                > float(extras.get("cdh1_periphery_mean_late", math.nan)),
                "CDH1 center > CDH1 periphery gradient",
                "CDH1 gradient absent/inverted: check ZEB1/CDH1 wiring and spatial signals.",
            ),
        ]
        return self._finalize_anchor_tests(extras, tests)

    def _eval_anchor7_ecm_degradation_limits(
        self, metrics: SimulationMetrics, extras: Dict[str, Any], baseline: Dict[str, Any]
    ) -> Tuple[bool, str]:
        shh_no_cyto = self._cached_metrics(baseline, "ANCHOR_3_SHH_PARADOX")
        control_nodrug = self._cached_metrics(baseline, "ANCHOR_1_SELF_ASSEMBLY")
        drug_only = self._cached_metrics(baseline, "ANCHOR_2_DRUG_PENETRATION_MATURITY")
        if shh_no_cyto is None or drug_only is None or control_nodrug is None:
            return False, "missing ANCHOR_3_SHH_PARADOX or ANCHOR_2_DRUG_PENETRATION_MATURITY reference"

        viab_degrade_only = self._viability(shh_no_cyto)
        viab_control = self._viability(control_nodrug)
        viab_degrade_drug = self._viability(asdict(metrics))
        viab_control_drug = self._viability(drug_only)
        ecm_recovery = float(extras.get("ecm_delta", math.nan))
        tests = [
            self._make_test(
                "P7a",
                viab_degrade_only >= viab_control,
                "Viability (Degrade-only) >= Viability (Control)",
                "Degrade-only unexpectedly cytotoxic: inspect unintended death channels.",
            ),
            self._make_test(
                "P7b",
                float(metrics.drug_penetration) > float(drug_only.get("drug_penetration", math.nan)),
                "Drug at centroid (Degrade+Drug) > Drug at centroid (Control+Drug)",
                "No penetration gain from ECM degradation: inspect ECM->diffusion coupling.",
            ),
            self._make_test(
                "P7c",
                viab_degrade_drug < viab_control_drug,
                "Viability (Degrade+Drug) < Viability (Control+Drug)",
                "Combination not outperforming drug-alone: barrier opening not functionally connected.",
            ),
            self._make_test(
                "P7d",
                ecm_recovery > 0.0,
                "ECM recovery rate after degradation > 0",
                "Barrier not rebuilding: check CAF persistence and ECM production loop.",
            ),
            self._make_test(
                "P7e",
                viab_degrade_drug > 0.0,
                "Viability (Degrade+Drug) > 0 (residual survivors remain)",
                "Complete eradication indicates resistance/sanctuary machinery too weak.",
            ),
        ]
        return self._finalize_anchor_tests(extras, tests)

    def _eval_anchor7b_ecm_degrade_only(
        self, metrics: SimulationMetrics, extras: Dict[str, Any], baseline: Dict[str, Any]
    ) -> Tuple[bool, str]:
        # Data collection arm only — directional tests run in CHECK_4_FITNESS_RANKING.
        extras["anchor_tests"] = []
        return True, "ECM-degrade-only data collected (0 directional tests; ranking evaluated in Check 4)"

    def _eval_anchor8_arm_data(
        self, metrics: SimulationMetrics, extras: Dict[str, Any], baseline: Dict[str, Any]
    ) -> Tuple[bool, str]:
        # Data collection arm only — directional tests run in ANCHOR_8_BOTH_DEPLETED_DRUG.
        extras["anchor_tests"] = []
        return True, "arm data collected (0 directional tests; independence evaluated in BOTH_DEPLETED_DRUG arm)"

    def _eval_anchor8_independence(
        self, metrics: SimulationMetrics, extras: Dict[str, Any], baseline: Dict[str, Any]
    ) -> Tuple[bool, str]:
        ha_arm  = self._cached_metrics(baseline, "ANCHOR_8_HA_DEPLETED_DRUG")
        col_arm = self._cached_metrics(baseline, "ANCHOR_8_COL_DEPLETED_DRUG")
        if ha_arm is None or col_arm is None:
            return False, "missing HA-depleted or COL-depleted arm cache for Anchor 8 independence test"
        both = asdict(metrics)

        pen_ha   = float(ha_arm.get("drug_penetration", math.nan))
        pen_col  = float(col_arm.get("drug_penetration", math.nan))
        pen_both = float(both.get("drug_penetration", math.nan))
        ext_ha   = float(ha_arm.get("tumor_extent", math.nan))
        ext_col  = float(col_arm.get("tumor_extent", math.nan))

        tests = [
            self._make_test(
                "P8a",
                pen_ha > pen_col,
                "Drug penetration (HA-depleted) > Drug penetration (Collagen-depleted)",
                "HA not the dominant diffusion barrier: check ecm_ha_fraction weighting in barrier formula.",
            ),
            self._make_test(
                "P8b",
                ext_col > ext_ha,
                "Tumor spatial extent (Collagen-depleted) > Tumor spatial extent (HA-depleted)",
                "Collagen not dominant mechanical confinement: check mech_col_frac weighting.",
            ),
            self._make_test(
                "P8c",
                pen_both > pen_ha,
                "Drug penetration (Both-depleted) > Drug penetration (HA-depleted alone)",
                "Dual depletion no better than HA alone: collagen contributes to diffusion barrier.",
            ),
            self._make_test(
                "P8d",
                pen_both > pen_col,
                "Drug penetration (Both-depleted) > Drug penetration (Collagen-depleted alone)",
                "Dual depletion no better than collagen alone: HA contribution to barrier not additive.",
            ),
        ]
        return self._finalize_anchor_tests(extras, tests)

    def _eval_anchor9_density_prognostic(
        self, metrics: SimulationMetrics, extras: Dict[str, Any], baseline: Dict[str, Any]
    ) -> Tuple[bool, str]:
        ref = self._cached_metrics(baseline, "ANCHOR_2_DRUG_PENETRATION_MATURITY")
        if ref is None:
            return False, "missing ANCHOR_2_DRUG_PENETRATION_MATURITY reference"

        tests = [
            self._make_test(
                "P9a",
                float(metrics.live_tumor_cells) > float(ref.get("live_tumor_cells", math.nan)),
                "Viability (high PSC) > Viability (medium PSC control)",
                "Drug efficacy insensitive to barrier density: inspect ECM->drug coupling.",
            ),
            self._make_test(
                "P9b",
                float(metrics.stroma_barrier_score) > float(ref.get("stroma_barrier_score", math.nan)),
                "Pre-treatment barrier density (high PSC) > medium PSC control",
                "PSC sweep not changing barrier thickness; anchor setup invalid.",
            ),
            self._make_test(
                "P9c",
                (
                    (float(metrics.stroma_barrier_score) - float(ref.get("stroma_barrier_score", math.nan)))
                    * (float(metrics.live_tumor_cells) - float(ref.get("live_tumor_cells", math.nan)))
                ) > 0.0,
                "Barrier density and post-treatment viability move in same direction",
                "Expected positive barrier-response relationship absent.",
            ),
        ]
        return self._finalize_anchor_tests(extras, tests)

    def _eval_anchor10_spatial_sanctuary(
        self, metrics: SimulationMetrics, extras: Dict[str, Any], baseline: Dict[str, Any]
    ) -> Tuple[bool, str]:
        tests = [
            self._make_test(
                "P10a",
                float(extras.get("live_mean_local_ecm", math.nan))
                > float(extras.get("dead_mean_local_ecm", math.nan)),
                "Mean ECM at survivor locations > mean ECM at dead-cell locations",
                "Survivors not enriched in dense ECM: sanctuary effect missing.",
            ),
            self._make_test(
                "P10b",
                float(extras.get("mean_distance_survivor_from_boundary", math.nan))
                > float(extras.get("mean_distance_dead_from_boundary", math.nan)),
                "Survivors are farther from boundary than dead cells",
                "Distance shielding absent: drug penetration likely too uniform.",
            ),
            self._make_test(
                "P10d",
                float(extras.get("tumor_count_end", math.nan))
                > float(extras.get("tumor_count_min", math.nan)),
                "Tumor count after withdrawal window > count at nadir",
                "Residual population not regrowing: check genotype preservation/proliferation recovery.",
            ),
            self._make_test(
                "P10e",
                float(extras.get("mean_ecm_late", math.nan))
                > float(extras.get("mean_ecm_at_nadir", math.nan)),
                "ECM density after regrowth window > ECM density at nadir",
                "Barrier not re-emerging during regrowth: CAF->ECM loop likely broken.",
            ),
        ]
        derived = [
            "P10c is implied by concurrent P10a and P10b.",
            "P10f tracked as directional note: activated_cafs_late > activated_cafs_at_nadir.",
        ]
        if "activated_cafs_late" in extras and "activated_cafs_at_nadir" in extras:
            extras["p10f_directional"] = (
                float(extras.get("activated_cafs_late", math.nan))
                > float(extras.get("activated_cafs_at_nadir", math.nan))
            )
        return self._finalize_anchor_tests(extras, tests, derived_notes=derived)

    # ── §4 Reality Check evaluators (virtual — no simulation launched) ─────────

    def _eval_check1_natural_history(
        self, metrics: Optional[SimulationMetrics], extras: Dict[str, Any], baseline: Dict[str, Any]
    ) -> Tuple[bool, str]:
        a1 = self._cached_extras(baseline, "ANCHOR_1_SELF_ASSEMBLY")
        a5 = self._cached_extras(baseline, "ANCHOR_5_CENTRAL_HYPOXIA")
        a6 = self._cached_extras(baseline, "ANCHOR_6_PERIPHERAL_EMT")
        tests = [
            self._make_test(
                "C1a",
                float(a1.get("activated_caf_fraction", math.nan)) > 0.0,
                "Stroma self-assembles: CAF activation fraction > 0",
                "CAF activation absent from Anchor 1; check tumor->PSC paracrine loop.",
            ),
            self._make_test(
                "C1b",
                float(a1.get("ecm_peri_tumor_shell_late", math.nan))
                > float(a1.get("ecm_domain_boundary_late", math.nan)),
                "ECM gradient present: peri-tumor ECM > domain-boundary ECM",
                "No ECM gradient from Anchor 1; spatial ECM accumulation missing.",
            ),
            self._make_test(
                "C1c",
                float(a5.get("oxygen_centroid_late", math.nan))
                < float(a5.get("oxygen_periphery_late", math.nan)),
                "Central hypoxia: core oxygen < peripheral oxygen",
                "No oxygen gradient from Anchor 5; ECM->diffusion coupling missing.",
            ),
            self._make_test(
                "C1d",
                float(a6.get("emt_periphery_fraction_late", math.nan))
                > float(a6.get("emt_core_fraction_late", math.nan)),
                "Peripheral EMT: peripheral mesenchymal fraction > core mesenchymal fraction",
                "EMT not spatially enriched from Anchor 6; check ZEB1 vs oxygen gradient.",
            ),
        ]
        extras["check_tests"] = tests
        failed = [t for t in tests if not bool(t["passed"])]
        if not failed:
            return True, f"CHECK 1 PASS — {len(tests)}/{len(tests)} natural history hallmarks confirmed"
        fail_ids = ", ".join(t["id"] for t in failed)
        return False, f"CHECK 1 FAIL — {len(tests) - len(failed)}/{len(tests)} hallmarks; failed: {fail_ids}"

    def _eval_check2_chemo_failure(
        self, metrics: Optional[SimulationMetrics], extras: Dict[str, Any], baseline: Dict[str, Any]
    ) -> Tuple[bool, str]:
        a1_m = self._cached_metrics(baseline, "ANCHOR_1_SELF_ASSEMBLY")
        a2_m = self._cached_metrics(baseline, "ANCHOR_2_DRUG_PENETRATION_MATURITY")
        a7_m = self._cached_metrics(baseline, "ANCHOR_7_ECM_DEGRADATION_LIMITS")
        a10_x = self._cached_extras(baseline, "ANCHOR_10_SPATIAL_SANCTUARY")
        a1_x = self._cached_extras(baseline, "ANCHOR_1_SELF_ASSEMBLY")
        a7_x = self._cached_extras(baseline, "ANCHOR_7_ECM_DEGRADATION_LIMITS")
        if a1_m is None or a2_m is None or a7_m is None:
            return False, "missing ANCHOR_1 / ANCHOR_2 / ANCHOR_7 metrics for Check 2"

        # C2a: Partial response — drug reduces tumor count vs untreated control
        tc_a2 = float(a2_m.get("total_tumor_cells", math.nan))
        tc_a1 = float(a1_m.get("total_tumor_cells", math.nan))
        # C2b: Barrier persistence — ECM survives drug treatment at ≥80% of pre-treatment density
        barrier_a2 = float(a2_m.get("stroma_barrier_score", math.nan))
        barrier_a1 = float(a1_m.get("stroma_barrier_score", math.nan))
        barrier_ratio = barrier_a2 / barrier_a1 if (math.isfinite(barrier_a1) and barrier_a1 > 0) else math.nan
        # C2c: Regrowth after withdrawal — survivors rebuild (from A10's time series)
        tc_end = float(a10_x.get("tumor_count_end", math.nan))
        tc_min = float(a10_x.get("tumor_count_min", math.nan))
        # C2d: Resistance emergence — ABCB1 rises under drug pressure
        abcb1_a7 = float(a7_x.get("mean_abcb1_surviving_tumor", math.nan))
        abcb1_a1 = float(a1_x.get("mean_abcb1_surviving_tumor", math.nan))

        tests = [
            self._make_test(
                "C2a",
                tc_a2 < tc_a1,
                "Drug reduces tumor count: tc(A2: drug) < tc(A1: control)",
                "Drug ineffective: tumor count with drug not lower than untreated control.",
            ),
            self._make_test(
                "C2b",
                math.isfinite(barrier_ratio) and barrier_ratio >= 0.8,
                f"Barrier persists through treatment: stroma(A2)/stroma(A1) = {barrier_ratio:.2f} >= 0.80",
                "Barrier demolished by drug: ECM did not persist at ≥80% of untreated density.",
            ),
            self._make_test(
                "C2c",
                math.isfinite(tc_end) and math.isfinite(tc_min) and tc_end > tc_min,
                "Sanctuary survivors regrow after drug withdrawal: tc(endpoint) > tc(nadir) in A10",
                "No regrowth after withdrawal: surviving population static; check genotype preservation.",
            ),
            self._make_test(
                "C2d",
                abcb1_a7 > abcb1_a1,
                "ABCB1 resistance emerges: mean ABCB1 in A7 survivors > A1 baseline",
                "No efflux pump induction: Rule 11 (efflux_induction_delay) temporal dynamics not working.",
            ),
        ]
        extras["check_tests"] = tests
        failed = [t for t in tests if not bool(t["passed"])]
        if not failed:
            return True, f"CHECK 2 PASS — {len(tests)}/{len(tests)} chemo-failure criteria confirmed"
        fail_ids = ", ".join(t["id"] for t in failed)
        return False, f"CHECK 2 FAIL — {len(tests) - len(failed)}/{len(tests)} criteria; failed: {fail_ids}"

    def _eval_check3_vismodegib_paradox(
        self, metrics: Optional[SimulationMetrics], extras: Dict[str, Any], baseline: Dict[str, Any]
    ) -> Tuple[bool, str]:
        shh_drug_m = self._cached_metrics(baseline, "ANCHOR_7_ECM_DEGRADATION_LIMITS")  # SHH-off + drug
        control_m  = self._cached_metrics(baseline, "ANCHOR_1_SELF_ASSEMBLY")            # no drug, no intervention
        shh_only_m = self._cached_metrics(baseline, "ANCHOR_3_SHH_PARADOX")             # SHH-off, no drug
        if shh_drug_m is None or control_m is None or shh_only_m is None:
            return False, "missing A1/A3/A7 metrics for Check 3 three-way rank"

        tc_shh_drug = float(shh_drug_m.get("total_tumor_cells", math.nan))
        tc_control  = float(control_m.get("total_tumor_cells", math.nan))
        tc_shh_only = float(shh_only_m.get("total_tumor_cells", math.nan))
        tests = [
            self._make_test(
                "C3a",
                tc_shh_drug < tc_control,
                "Tumor count (SHH-off+drug) < Tumor count (control)",
                "SHH inhibition + drug not better than untreated: check vismodegib->barrier-opening->drug chain.",
            ),
            self._make_test(
                "C3b",
                tc_control < tc_shh_only,
                "Tumor count (control) < Tumor count (SHH-off alone)",
                "SHH inhibition without drug not worse than untreated: SHH paradox missing.",
            ),
        ]
        extras["check_tests"] = tests
        failed = [t for t in tests if not bool(t["passed"])]
        if not failed:
            return True, (
                f"CHECK 3 PASS — three-way rank confirmed: "
                f"SHH+drug({tc_shh_drug:.0f}) < control({tc_control:.0f}) < SHH-only({tc_shh_only:.0f})"
            )
        fail_ids = ", ".join(t["id"] for t in failed)
        return False, f"CHECK 3 FAIL — three-way rank broken; failed: {fail_ids}"

    def _eval_check4_fitness_ranking(
        self, metrics: Optional[SimulationMetrics], extras: Dict[str, Any], baseline: Dict[str, Any]
    ) -> Tuple[bool, str]:
        a8_m = self._cached_metrics(baseline, "ANCHOR_8_BOTH_DEPLETED_DRUG")
        a2_m = self._cached_metrics(baseline, "ANCHOR_2_DRUG_PENETRATION_MATURITY")
        a7b_m = self._cached_metrics(baseline, "ANCHOR_7B_ECM_DEGRADE_ONLY")
        a3_m = self._cached_metrics(baseline, "ANCHOR_3_SHH_PARADOX")
        if any(m is None for m in (a8_m, a2_m, a7b_m, a3_m)):
            return False, "missing A2/A3/A7B/A8_BOTH metrics for Check 4 fitness ranking"

        # Four-strategy ranking by total tumor count (ascending = better outcome):
        # ECM-degrade+drug < drug-only < ECM-degrade-only < SHH-only
        tc_a8 = float(a8_m.get("total_tumor_cells", math.nan))
        tc_a2 = float(a2_m.get("total_tumor_cells", math.nan))
        tc_7b = float(a7b_m.get("total_tumor_cells", math.nan))
        tc_a3 = float(a3_m.get("total_tumor_cells", math.nan))

        tests = [
            self._make_test(
                "C4a",
                tc_a8 < tc_a2,
                f"tc(ECM-degrade+drug)={tc_a8:.0f} < tc(drug-only)={tc_a2:.0f}",
                "Combination not better than standard of care: check barrier decomposition + drug synergy.",
            ),
            self._make_test(
                "C4b",
                tc_a2 < tc_7b,
                f"tc(drug-only)={tc_a2:.0f} < tc(ECM-degrade-only)={tc_7b:.0f}",
                "Drug alone not better than barrier removal alone: drug mechanism non-functional.",
            ),
            self._make_test(
                "C4c",
                tc_7b < tc_a3,
                f"tc(ECM-degrade-only)={tc_7b:.0f} < tc(SHH-only)={tc_a3:.0f}",
                "Barrier removal not better than actively harmful SHH-only: ECM degrade ineffective.",
            ),
            self._make_test(
                "C4d",
                tc_a2 < tc_a3,
                f"tc(drug-only)={tc_a2:.0f} < tc(SHH-only)={tc_a3:.0f}",
                "Standard of care not better than SHH-only: drug or SHH paradox mechanism broken.",
            ),
        ]
        extras["check_tests"] = tests
        failed = [t for t in tests if not bool(t["passed"])]
        if not failed:
            return True, (
                f"CHECK 4 PASS — full ranking: ECM+drug({tc_a8:.0f}) < drug({tc_a2:.0f}) "
                f"< ECM-only({tc_7b:.0f}) < SHH-only({tc_a3:.0f})"
            )
        fail_ids = ", ".join(t["id"] for t in failed)
        return False, f"CHECK 4 FAIL — fitness hierarchy broken; failed: {fail_ids}"

    def _eval_check5_sanctuary_regrowth(
        self, metrics: Optional[SimulationMetrics], extras: Dict[str, Any], baseline: Dict[str, Any]
    ) -> Tuple[bool, str]:
        # Gate: all four prerequisite checks must have passed.
        for gate in (
            "CHECK_1_NATURAL_HISTORY",
            "CHECK_2_CHEMO_FAILURE",
            "CHECK_3_VISMODEGIB_PARADOX",
            "CHECK_4_FITNESS_RANKING",
        ):
            gate_tests = self._cached_extras(baseline, gate).get("check_tests", [])
            if not gate_tests or not all(bool(t.get("passed", False)) for t in gate_tests):
                extras["check_tests"] = []
                return False, f"CHECK 5 skipped: prerequisite {gate} did not pass"

        a10_extras = self._cached_extras(baseline, "ANCHOR_10_SPATIAL_SANCTUARY")
        a10_tests  = a10_extras.get("anchor_tests", [])
        a10_by_id  = {t["id"]: t for t in a10_tests if isinstance(t, dict)}
        tests = [
            self._make_test(
                "C5a",
                bool(a10_by_id.get("P10a", {}).get("passed", False)),
                "Survivors localize to dense ECM sanctuary (P10a passed in Anchor 10)",
                "ECM sanctuary enrichment absent: P10a failed; check ECM->drug coupling.",
            ),
            self._make_test(
                "C5b",
                bool(a10_by_id.get("P10b", {}).get("passed", False)),
                "Survivors farther from boundary than dead cells (P10b passed in Anchor 10)",
                "Distance shielding absent: P10b failed; drug penetration too uniform.",
            ),
            self._make_test(
                "C5c",
                bool(a10_by_id.get("P10d", {}).get("passed", False)),
                "Tumor regrowth after drug withdrawal (P10d passed in Anchor 10)",
                "No post-withdrawal regrowth: P10d failed; check genotype preservation/proliferation.",
            ),
            self._make_test(
                "C5d",
                bool(a10_by_id.get("P10e", {}).get("passed", False)),
                "ECM barrier re-assembles during regrowth window (P10e passed in Anchor 10)",
                "Barrier not re-emerging: P10e failed; check CAF->ECM loop.",
            ),
        ]
        extras["check_tests"] = tests
        failed = [t for t in tests if not bool(t["passed"])]
        if not failed:
            return True, f"CHECK 5 PASS — {len(tests)}/{len(tests)} sanctuary-regrowth criteria confirmed"
        fail_ids = ", ".join(t["id"] for t in failed)
        return False, f"CHECK 5 FAIL — {len(tests) - len(failed)}/{len(tests)} criteria; failed: {fail_ids}"


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
    parser.add_argument(
        "--replicates-per-arm",
        type=int,
        default=3,
        help="Replicate runs per scenario arm (decision protocol recommends >=3).",
    )
    parser.add_argument(
        "--random-seed-base",
        type=int,
        default=1001,
        help="Base random seed used to deterministically fan out replicate seeds.",
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
        replicates_per_arm=args.replicates_per_arm,
        random_seed_base=args.random_seed_base,
    )
    summary = validator.run_all()

    prefix = Path(args.output_dir).expanduser().resolve() / "biology_validation_summary"
    json_path, csv_path = validator.save_summary(summary, prefix)

    print(
        "Biology validation complete: "
        f"{summary.passed}/{summary.total} anchors passed; "
        f"{summary.tests_passed}/{summary.tests_total} directional tests passed"
    )
    print(f"JSON summary: {json_path}")
    print(f"CSV summary: {csv_path}")
    if math.isfinite(summary.median_wall_time_per_run_seconds):
        est_minutes = summary.estimated_full_validation_wall_time_seconds / 60.0
        print(
            "Estimated pre-EA validation budget (~60 runs): "
            f"{est_minutes:.2f} min at median {summary.median_wall_time_per_run_seconds:.2f}s/run"
        )
    for scenario in summary.scenario_results:
        status = "PASS" if scenario.success else "FAIL"
        print(f"- [{status}] Anchor {scenario.anchor_id} {scenario.scenario}: {scenario.details}")
    return 0 if summary.failed == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
