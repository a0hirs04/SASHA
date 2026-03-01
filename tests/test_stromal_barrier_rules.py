from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[1]
BUILD_TEST_DIR = PROJECT_ROOT / "build" / "tests"
BUILD_LOG_DIR = PROJECT_ROOT / "build" / "test_logs"
OBJ_DIR = PROJECT_ROOT / "build" / "obj" / "release"

CPP_TEST_SOURCES = {
    "module3_stromal_activation_test": PROJECT_ROOT / "tests" / "module3_stromal_activation_test.cpp",
    "module6_ecm_production_test": PROJECT_ROOT / "tests" / "module6_ecm_production_test.cpp",
    "e12_caf_tgfb_feedback_closing_test": PROJECT_ROOT / "tests" / "e12_caf_tgfb_feedback_closing_test.cpp",
    "microenvironment_five_fields_config_test": PROJECT_ROOT / "tests" / "microenvironment_five_fields_config_test.cpp",
    "e11_e13_gli1_boost_edges_test": PROJECT_ROOT / "tests" / "e11_e13_gli1_boost_edges_test.cpp",
    "rule30_caf_division_boost_test": PROJECT_ROOT / "tests" / "rule30_caf_division_boost_test.cpp",
    "stromal_cell_definition_actor2_test": PROJECT_ROOT / "tests" / "stromal_cell_definition_actor2_test.cpp",
    "ecm_diffusion_coupling_test": PROJECT_ROOT / "tests" / "ecm_diffusion_coupling_test.cpp",
    "rule36_ecm_traps_paracrine_test": PROJECT_ROOT / "tests" / "rule36_ecm_traps_paracrine_test.cpp",
    "e23_ecm_tgfb_trapping_edge_test": PROJECT_ROOT / "tests" / "e23_ecm_tgfb_trapping_edge_test.cpp",
    "e26_mechanical_to_contact_inhibition_feedback_closing_test": PROJECT_ROOT / "tests" / "e26_mechanical_to_contact_inhibition_feedback_closing_test.cpp",
    "e27_e28_compaction_secondary_edges_test": PROJECT_ROOT / "tests" / "e27_e28_compaction_secondary_edges_test.cpp",
    "e39_stroma_removal_releases_confinement_test": PROJECT_ROOT / "tests" / "e39_stroma_removal_releases_confinement_test.cpp",
}

RULE_MATRIX = [
    {
        "rules": "22,23",
        "binary": "module3_stromal_activation_test",
        "tokens": [
            "PASS Rule22_A_TGFB_alone",
            "PASS Rule22_B_SHH_alone",
            "PASS Rule22_C_combined_above_threshold",
            "PASS Rule22_D_combined_below_threshold",
            "PASS Rule23_irreversible_ACTA2",
            "PASS module3_stromal_activation_test",
        ],
    },
    {
        "rules": "25,26",
        "binary": "module6_ecm_production_test",
        "tokens": [
            "PASS Rule25_caf_produces_HA_and_collagen",
            "PASS Rule26_gli1_boosts_ecm_production",
            "PASS module6_ecm_production_test",
        ],
    },
    {
        "rules": "27",
        "binary": "e12_caf_tgfb_feedback_closing_test",
        "tokens": ["PASS E12 scenario"],
    },
    {
        "rules": "28,29",
        "binary": "microenvironment_five_fields_config_test",
        "tokens": [
            "PASS Rule28_ecm_non_diffusive",
            "PASS Rule29_ecm_near_permanent",
            "PASS microenvironment_five_fields_config_test",
        ],
    },
    {
        "rules": "30",
        "binary": "rule30_caf_division_boost_test",
        "tokens": [
            "PASS Rule30_caf_division_gli1_boost",
            "PASS rule30_caf_division_boost_test",
        ],
    },
    {
        "rules": "31,32",
        "binary": "stromal_cell_definition_actor2_test",
        "tokens": [
            "PASS Rule31_caf_drug_indifferent",
            "PASS Rule32_caf_zero_death",
            "PASS stromal_cell_definition_actor2_test",
        ],
    },
    {
        "rules": "33,34,35",
        "binary": "ecm_diffusion_coupling_test",
        "tokens": [
            "PASS Rule33_monotonic_diffusion_impedance",
            "PASS Rule34_ha_impedes_more_than_collagen",
            "PASS Rule35_collagen_stiffer_than_HA",
            "PASS ecm_diffusion_coupling_test",
        ],
    },
    {
        "rules": "36",
        "binary": "rule36_ecm_traps_paracrine_test",
        "tokens": [
            "PASS Rule36_ecm_traps_paracrine_signals",
            "PASS rule36_ecm_traps_paracrine_test",
        ],
    },
    {
        "rules": "37",
        "binary": "e26_mechanical_to_contact_inhibition_feedback_closing_test",
        "tokens": ["PASS E26 integration"],
    },
    {
        "rules": "38",
        "binary": "e27_e28_compaction_secondary_edges_test",
        "tokens": ["PASS E27/E28 scenario"],
    },
    {
        "rules": "39",
        "binary": "e39_stroma_removal_releases_confinement_test",
        "tokens": [
            "PASS Rule39_stroma_removal_worsens_tumor_growth",
            "PASS e39_stroma_removal_releases_confinement_test",
        ],
    },
]


def _run(cmd: list[str], cwd: Path) -> None:
    subprocess.run(cmd, cwd=cwd, check=True)


def _ensure_release_objects() -> None:
    _run(
        [
            "make",
            "-j",
            str(max(1, os.cpu_count() or 1)),
            "BUILD=release",
            "PHYSICELL_DIR=Stroma_world/PhysiCell",
            "build-mode",
        ],
        PROJECT_ROOT,
    )


def _compile_cpp_test(binary_name: str) -> Path:
    source = CPP_TEST_SOURCES[binary_name]
    assert source.exists(), f"Missing test source: {source}"

    BUILD_TEST_DIR.mkdir(parents=True, exist_ok=True)
    binary = BUILD_TEST_DIR / binary_name

    objects = sorted(OBJ_DIR.rglob("*.o"))
    objects = [o for o in objects if not str(o).replace("\\", "/").endswith("build/generated/main.o")]
    assert objects, "No release objects found in build/obj/release; build step failed"

    cmd = [
        "g++",
        "-I.",
        "-IStroma_world/PhysiCell",
        "-IStroma_world/PhysiCell/BioFVM",
        "-std=c++17",
        "-march=native",
        "-m64",
        "-fomit-frame-pointer",
        "-mfpmath=both",
        "-fopenmp",
        "-O2",
        str(source),
        *[str(o) for o in objects],
        "-o",
        str(binary),
        "-fopenmp",
    ]
    _run(cmd, PROJECT_ROOT)
    return binary


def _run_binary(binary: Path) -> str:
    BUILD_LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_path = BUILD_LOG_DIR / f"{binary.name}.log"

    proc = subprocess.run(
        [str(binary)],
        cwd=PROJECT_ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    output = proc.stdout
    log_path.write_text(output, encoding="utf-8")

    assert proc.returncode == 0, (
        f"{binary.name} failed with exit code {proc.returncode}. "
        f"See {log_path}"
    )
    return output


def test_stromal_barrier_rules() -> None:
    _ensure_release_objects()

    missing_tokens: list[str] = []
    executed: list[str] = []

    for item in RULE_MATRIX:
        binary_name = item["binary"]
        binary = _compile_cpp_test(binary_name)
        output = _run_binary(binary)
        executed.append(f"rules={item['rules']}:{binary_name}")

        for token in item["tokens"]:
            if token not in output:
                missing_tokens.append(
                    f"{binary_name}: missing token '{token}'"
                )

    if missing_tokens:
        joined = "\n".join(missing_tokens)
        pytest.fail(f"Stromal barrier rules missing expected PASS markers:\n{joined}")

    # Smoke assertion on execution breadth.
    assert len(executed) == len(RULE_MATRIX)
