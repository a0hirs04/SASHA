from __future__ import annotations

from pathlib import Path

from python.validation.validate_biology import BiologyValidator, ScenarioResult, ScenarioSpec


def _make_result(spec: ScenarioSpec, success: bool, details: str = "ok") -> ScenarioResult:
    return ScenarioResult(
        anchor_id=spec.anchor_id,
        scenario=spec.name,
        success=success,
        criterion=spec.criterion,
        details=details,
        dependencies=list(spec.dependencies),
        run_dir="/tmp/fake",
        exit_code=0 if success else 1,
        wall_time=0.01,
        fitness=0.5 if success else 0.0,
        metrics={"total_tumor_cells": 50},
        extras={},
    )


def _make_validator(ea_test_paths: dict, tmp_path: Path) -> BiologyValidator:
    return BiologyValidator(
        binary_path=ea_test_paths["binary"],
        config_path=ea_test_paths["config"],
        output_dir=tmp_path / "validation_runs",
        timeout_seconds=5,
        sim_max_time=10.0,
    )


def test_anchor_suite_has_expected_order_and_dependencies(ea_test_paths, tmp_path):
    validator = _make_validator(ea_test_paths, tmp_path)
    specs = validator._scenario_specs()  # noqa: SLF001

    # 13 simulation anchors (1-10, Anchor-7B, three Anchor-8 arms) + 5 virtual Reality Checks = 18
    assert len(specs) == 18
    assert set(spec.anchor_id for spec in specs) == set(range(1, 16))
    assert specs[0].name == "ANCHOR_1_SELF_ASSEMBLY"
    assert specs[1].name == "ANCHOR_3_SHH_PARADOX"
    assert specs[-1].name == "CHECK_5_SANCTUARY_REGROWTH"

    by_name = {spec.name: spec for spec in specs}

    # Anchor 6 preserves its single-parent dependency
    assert by_name["ANCHOR_6_PERIPHERAL_EMT"].dependencies == ["ANCHOR_5_CENTRAL_HYPOXIA"]

    # Anchor 8 three-arm independence structure
    assert by_name["ANCHOR_8_HA_DEPLETED_DRUG"].dependencies == ["ANCHOR_7_ECM_DEGRADATION_LIMITS"]
    assert by_name["ANCHOR_8_COL_DEPLETED_DRUG"].dependencies == ["ANCHOR_7_ECM_DEGRADATION_LIMITS"]
    assert by_name["ANCHOR_8_BOTH_DEPLETED_DRUG"].dependencies == [
        "ANCHOR_8_HA_DEPLETED_DRUG",
        "ANCHOR_8_COL_DEPLETED_DRUG",
    ]

    # Anchor 10 includes the new Anchor 8 capstone arm
    assert "ANCHOR_9_BARRIER_DENSITY_PROGNOSTIC" in by_name["ANCHOR_10_SPATIAL_SANCTUARY"].dependencies
    assert "ANCHOR_8_BOTH_DEPLETED_DRUG" in by_name["ANCHOR_10_SPATIAL_SANCTUARY"].dependencies

    # Reality Checks are virtual and have the right gate dependencies
    assert by_name["CHECK_1_NATURAL_HISTORY"].virtual is True
    assert by_name["CHECK_5_SANCTUARY_REGROWTH"].virtual is True

    # CHECK_2 needs A1 (baseline), A2 (drug), A7 (treatment), A10 (regrowth)
    assert by_name["CHECK_2_CHEMO_FAILURE"].dependencies == [
        "ANCHOR_1_SELF_ASSEMBLY",
        "ANCHOR_2_DRUG_PENETRATION_MATURITY",
        "ANCHOR_7_ECM_DEGRADATION_LIMITS",
        "ANCHOR_10_SPATIAL_SANCTUARY",
    ]

    # CHECK_4 ranks four treatment strategies
    assert by_name["CHECK_4_FITNESS_RANKING"].dependencies == [
        "ANCHOR_2_DRUG_PENETRATION_MATURITY",
        "ANCHOR_3_SHH_PARADOX",
        "ANCHOR_7B_ECM_DEGRADE_ONLY",
        "ANCHOR_8_BOTH_DEPLETED_DRUG",
    ]

    assert by_name["CHECK_5_SANCTUARY_REGROWTH"].dependencies == [
        "CHECK_1_NATURAL_HISTORY",
        "CHECK_2_CHEMO_FAILURE",
        "CHECK_3_VISMODEGIB_PARADOX",
        "CHECK_4_FITNESS_RANKING",
        "ANCHOR_10_SPATIAL_SANCTUARY",
    ]


def test_dependency_gate_skips_all_downstream_when_anchor1_fails(ea_test_paths, tmp_path, monkeypatch):
    validator = _make_validator(ea_test_paths, tmp_path)
    calls: list[str] = []

    def fake_run(spec: ScenarioSpec) -> ScenarioResult:
        calls.append(spec.name)
        if spec.name == "ANCHOR_1_SELF_ASSEMBLY":
            return _make_result(spec, success=False, details="forced fail")
        return _make_result(spec, success=True)

    monkeypatch.setattr(validator, "_run_scenario", fake_run)
    monkeypatch.setattr(validator, "_run_virtual_check", fake_run)

    summary = validator.run_all()

    assert calls == ["ANCHOR_1_SELF_ASSEMBLY"]
    assert summary.total == 18
    assert summary.passed == 0
    assert summary.failed == 18
    for result in summary.scenario_results[1:]:
        assert result.exit_code == -2
        assert "failed dependencies" in result.details


def test_dependency_gate_skips_only_affected_branches(ea_test_paths, tmp_path, monkeypatch):
    validator = _make_validator(ea_test_paths, tmp_path)
    calls: list[str] = []

    def fake_run(spec: ScenarioSpec) -> ScenarioResult:
        calls.append(spec.name)
        if spec.name == "ANCHOR_3_SHH_PARADOX":
            return _make_result(spec, success=False, details="forced fail")
        return _make_result(spec, success=True)

    monkeypatch.setattr(validator, "_run_scenario", fake_run)
    monkeypatch.setattr(validator, "_run_virtual_check", fake_run)

    summary = validator.run_all()
    by_name = {r.scenario: r for r in summary.scenario_results}

    # A3 failure propagates to A7, A7B, and all three A8 arms
    assert "ANCHOR_7_ECM_DEGRADATION_LIMITS" not in calls
    assert "ANCHOR_7B_ECM_DEGRADE_ONLY" not in calls
    assert "ANCHOR_8_HA_DEPLETED_DRUG" not in calls
    assert "ANCHOR_8_COL_DEPLETED_DRUG" not in calls
    assert "ANCHOR_8_BOTH_DEPLETED_DRUG" not in calls

    # Branches independent of A3 still run
    assert "ANCHOR_2_DRUG_PENETRATION_MATURITY" in calls
    assert "ANCHOR_9_BARRIER_DENSITY_PROGNOSTIC" in calls

    # A10 and checks that need A3/A7/A8 are skipped
    assert "ANCHOR_10_SPATIAL_SANCTUARY" not in calls

    assert by_name["ANCHOR_7_ECM_DEGRADATION_LIMITS"].exit_code == -2
    assert by_name["ANCHOR_7B_ECM_DEGRADE_ONLY"].exit_code == -2
    assert by_name["ANCHOR_8_HA_DEPLETED_DRUG"].exit_code == -2
    assert by_name["ANCHOR_8_COL_DEPLETED_DRUG"].exit_code == -2
    assert by_name["ANCHOR_8_BOTH_DEPLETED_DRUG"].exit_code == -2
    assert by_name["ANCHOR_10_SPATIAL_SANCTUARY"].exit_code == -2
