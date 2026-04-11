from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from dataclasses import asdict
from pathlib import Path

import pytest

from quarkConstraints.paper_0710_1869.artifacts import (
    write_default_paper_0710_1869_kaon_artifact_exports,
)
from quarkConstraints.paper_0710_1869.benchmarks import (
    PAPER_0710_1869_BENCHMARK_SCHEMA_ID,
    default_paper_0710_1869_pr1_benchmark,
    paper_0710_1869_pr1_benchmarks,
)
from quarkConstraints.paper_0710_1869.eft_deltaf2.matching_kkgluon import (
    PAPER_0710_1869_DELTAF2_MATCHING_OBSERVABLE_SUPPORT_STATUS_ID,
)
from quarkConstraints.paper_0710_1869.eft_deltaf2.operators import (
    PAPER_0710_1869_DELTAF2_OPERATOR_BASIS_ID,
    PAPER_0710_1869_DELTAF2_OPERATOR_DEFINITION_Q4_LR_ID,
    PAPER_0710_1869_DELTAF2_OPERATOR_DEFINITION_Q5_LR_ID,
    PAPER_0710_1869_DELTAF2_PROJECTOR_NORMALIZATION_ID,
    PAPER_0710_1869_DELTAF2_PROJECTOR_NORMALIZATION_NOTE,
)
from quarkConstraints.paper_0710_1869.eft_deltaf2.rg_inputs import (
    PAPER_0710_1869_DELTAF2_RG_LR_BASIS_CONTRACT_ID,
    PAPER_0710_1869_DELTAF2_RG_LR_BASIS_STATUS_ID,
    PAPER_0710_1869_DELTAF2_RG_LR_BMU_BASIS_ID,
    PAPER_0710_1869_DELTAF2_RG_SUPPORTED_OPERATOR_SUBSET_ID,
)
from quarkConstraints.paper_0710_1869.inputs import (
    PAPER_0710_1869_EQ3_SCHEMA_ID,
    PAPER_0710_1869_SOURCE_SCHEMA_ID,
    PAPER_0710_1869_TABLE_I_SCHEMA_ID,
    Paper07101869Eq3Example,
    Paper07101869SourceRef,
    Paper07101869TableIInputs,
    default_paper_0710_1869_eq3_example,
    default_paper_0710_1869_table_i_inputs,
)
from quarkConstraints.paper_0710_1869.validation import module_has_forbidden_import
from quarkConstraints.paper_0710_1869.verifier import TolerancePolicy

REPO_ROOT = Path(__file__).resolve().parents[1]
PACKAGE_ROOT = REPO_ROOT / "quarkConstraints" / "paper_0710_1869"
BENCHMARK_SCRIPT = REPO_ROOT / "scripts" / "benchmark_quark_0710_1869.py"
GOLDEN_ARTIFACT_DIR = (
    REPO_ROOT / "tests" / "golden" / "paper_0710_1869" / "default_kaon_np_only"
)
EXPECTED_ARTIFACT_FILENAMES = {
    "wilsons": "wilsons.json",
    "hadronic": "hadronic.json",
    "observables": "observables.json",
    "provenance": "provenance.json",
}
EXPECTED_ARTIFACT_ROWS = {"M12_K_NP.re", "M12_K_NP.im", "Delta_m_K_NP"}
EXPECTED_SUPPORTED_OPERATORS = ["Q1_VLL", "Q1_VRR"]
EXPECTED_FULL_RG_SUPPORTED_OPERATORS = ["Q1_VLL", "Q1_VRR", "Q4_LR", "Q5_LR"]
EXPECTED_GUARDED_LR_OPERATORS = ["Q4_LR", "Q5_LR"]
EXPECTED_GUARDED_LR_DEFINITION_IDS = [
    PAPER_0710_1869_DELTAF2_OPERATOR_DEFINITION_Q4_LR_ID,
    PAPER_0710_1869_DELTAF2_OPERATOR_DEFINITION_Q5_LR_ID,
]
EXPECTED_LR_PAPER_OPERATOR_ORDER = ["Q4_LR", "Q5_LR"]
EXPECTED_LR_BMU_OPERATOR_ORDER = ["Q1_LR_BMU", "Q2_LR_BMU"]
EXPECTED_PAPER_TO_BMU_OPERATOR_MAP = [[0.0, 2.0], [1.0, 0.0]]
EXPECTED_BMU_TO_PAPER_OPERATOR_MAP = [[0.0, 1.0], [0.5, 0.0]]
EXPECTED_PAPER_TO_BMU_WILSON_MAP = [[0.0, 0.5], [1.0, 0.0]]
EXPECTED_BMU_TO_PAPER_WILSON_MAP = [[0.0, 1.0], [2.0, 0.0]]
EXPECTED_MATCHING_Q1_VLL = {"real": -1.1004498898491606e-12, "imag": -7.601544537269597e-13}
EXPECTED_RG_Q1_VLL = {"real": -1.509786848232995e-12, "imag": -1.0429109107548845e-12}
EXPECTED_M12_K_NP = {"real": -1.3495753042583394e-14, "imag": -9.322420653906486e-15}
EXPECTED_DELTA_M_K_NP_GEV = -2.6991506085166787e-14
EXPECTED_SYNTHETIC_LR_INPUT = [
    {"real": 1.25, "imag": -0.5},
    {"real": -0.75, "imag": 0.25},
]
EXPECTED_SYNTHETIC_BMU_INPUT = [
    {"real": -0.375, "imag": 0.125},
    {"real": 1.25, "imag": -0.5},
]
EXPECTED_CUSTOM_TOTAL_SCOPE_ID = "kaon.np_only.custom_total.q1_plus_lr.v1"
EXPECTED_CUSTOM_TOTAL_INTERPRETATION_ID = "kaon.np_only.custom_total.q1_plus_lr.v1"
EXPECTED_CUSTOM_TOTAL_M12_OBSERVABLE_ID = "M12_K_NP_CUSTOM_TOTAL"
EXPECTED_CUSTOM_TOTAL_DELTA_M_OBSERVABLE_ID = "Delta_m_K_NP_CUSTOM_TOTAL"
EXPECTED_CUSTOM_B_Q1_PROBE_CONFIG = {
    "B_d": {
        "scope_id": "bd.np_only.custom_q1.m12.v1",
        "interpretation_id": "bd.np_only.custom_q1.v1",
        "m12_id": "M12_Bd_NP_CUSTOM_Q1",
        "delta_m_id": "Delta_m_Bd_NP_CUSTOM_Q1",
    },
    "B_s": {
        "scope_id": "bs.np_only.custom_q1.m12.v1",
        "interpretation_id": "bs.np_only.custom_q1.v1",
        "m12_id": "M12_Bs_NP_CUSTOM_Q1",
        "delta_m_id": "Delta_m_Bs_NP_CUSTOM_Q1",
    },
}


def _assert_lr_status_semantics(status_id: object) -> None:
    assert isinstance(status_id, str)
    assert status_id == PAPER_0710_1869_DELTAF2_RG_LR_BASIS_STATUS_ID
    assert ".lr_rg_active." in status_id
    assert ".custom_lr_hadronic_active." in status_id
    assert ".custom_lr_only_observable_active." in status_id
    assert ".custom_combined_observable_active." in status_id
    assert status_id.endswith(".default_export_q1_only.v6")


def _canonical_tolerance_policy() -> dict[str, object]:
    return asdict(TolerancePolicy())


def _run_acceptance_benchmark(*args: str) -> dict[str, object]:
    completed = subprocess.run(
        [sys.executable, str(BENCHMARK_SCRIPT), "--emit-json", *args],
        check=False,
        capture_output=True,
        text=True,
        cwd=REPO_ROOT,
    )
    assert completed.stdout, completed.stderr
    return json.loads(completed.stdout)


def _canonical_payload_sha256(value: object) -> str:
    canonical_json = json.dumps(value, indent=2, sort_keys=True) + "\n"
    return hashlib.sha256(canonical_json.encode("utf-8")).hexdigest()


def _assert_complex_payload(value: object, expected: dict[str, float]) -> None:
    assert isinstance(value, dict)
    assert float(value["real"]) == pytest.approx(expected["real"], rel=0.0, abs=1.0e-24)
    assert float(value["imag"]) == pytest.approx(expected["imag"], rel=0.0, abs=1.0e-24)


def _assert_complex_vector_payload(
    value: object,
    expected: list[dict[str, float]],
) -> None:
    assert isinstance(value, list)
    assert len(value) == len(expected)
    for observed, expected_item in zip(value, expected, strict=True):
        _assert_complex_payload(observed, expected_item)


def _assert_custom_b_observable_probe(summary: dict[str, object], system_id: str) -> None:
    system_key = system_id.lower().replace("_", "")
    probe = summary[f"custom_{system_key}_observable_probe"]
    checks = probe["checks"]
    expected = EXPECTED_CUSTOM_B_Q1_PROBE_CONFIG[system_id]

    assert probe["status"] == "ok"
    assert probe["deterministic"] is True
    assert probe["system_id"] == system_id
    assert probe["declared_scheme_id"] == probe["scheme_id"]
    assert float(probe["declared_mu_had_GeV"]) == pytest.approx(
        float(probe["mu_had_GeV"]),
        rel=0.0,
        abs=1.0e-12,
    )
    assert float(probe["matching_scale_GeV"]) == pytest.approx(
        float(probe["declared_mu_had_GeV"]),
        rel=0.0,
        abs=1.0e-12,
    )
    expected_m12 = (
        (
            complex(
                float(probe["probe_q1_vll"]["real"]),
                float(probe["probe_q1_vll"]["imag"]),
            )
            + complex(
                float(probe["probe_q1_vrr"]["real"]),
                float(probe["probe_q1_vrr"]["imag"]),
            )
        )
        * float(probe["q1_matrix_element_GeV4"])
    ) / (2.0 * float(probe["meson_mass_GeV"]))
    expected_delta_m = 2.0 * float(expected_m12.real)

    _assert_complex_payload(
        probe["M12_NP_GeV"],
        {"real": expected_m12.real, "imag": expected_m12.imag},
    )
    assert float(probe["delta_m_NP_GeV"]) == pytest.approx(
        expected_delta_m,
        rel=0.0,
        abs=1.0e-24,
    )
    assert checks["system_id_matches_probe"] is True
    assert checks["scope_id_is_custom_b_surface"] is True
    assert checks["interpretation_is_custom_b_surface"] is True
    assert checks["observable_ids_are_custom_b_surface"] is True
    assert checks["declared_scheme_matches_payload"] is True
    assert checks["declared_mu_had_matches_payload"] is True
    assert checks["matching_scale_matches_payload"] is True
    assert checks["m12_matches_hand_calculation"] is True
    assert checks["delta_m_matches_hand_calculation"] is True
    assert checks["probe_payload_is_deterministic"] is True
    assert checks["probe_payload_is_cross_process_deterministic"] is True
    assert checks["probe_payload_sha256_matches_repeat"] is True
    assert checks["probe_payload_sha256_matches_cross_process"] is True
    assert checks["q1_matrix_element_present"] is True
    assert checks["lr_coefficients_remain_zero"] is True
    assert checks["default_observables_unchanged"] is True
    assert checks["default_artifacts_unchanged"] is True
    assert probe["default_observables_unchanged"] is True
    assert probe["default_artifacts_unchanged"] is True
    assert probe["cross_process_deterministic"] is True

    payload = probe["payload"]
    assert probe["payload_sha256"] == _canonical_payload_sha256(payload)
    assert probe["same_process_repeat_payload_sha256"] == probe["payload_sha256"]
    assert probe["cross_process_payload_sha256"] == probe["payload_sha256"]
    assert payload["system_id"] == system_id
    assert payload["observable_scope_id"] == expected["scope_id"]
    assert payload["interpretation"] == expected["interpretation_id"]
    assert payload["m12_observable_id"] == expected["m12_id"]
    assert payload["delta_m_observable_id"] == expected["delta_m_id"]


def test_table_i_inputs_freeze_extracted_values_and_metadata() -> None:
    table_inputs = default_paper_0710_1869_table_i_inputs()

    assert table_inputs.schema_id == PAPER_0710_1869_TABLE_I_SCHEMA_ID
    assert table_inputs.mode_id == "paper_0710_1869"
    assert [sector.sector_id for sector in table_inputs.sectors] == ["Q", "u", "d"]

    expected_pairs = {
        "Q": [(0.64, 0.002), (0.59, 0.01), (0.46, 0.2)],
        "u": [(0.68, 0.0007), (0.53, 0.06), (-0.06, 0.8)],
        "d": [(0.65, 0.002), (0.60, 0.008), (0.58, 0.02)],
    }
    expected_details = {
        "Q": "c_Q/f_Q eigenvalues",
        "u": "c_u/f_u eigenvalues",
        "d": "c_d/f_d eigenvalues",
    }

    for sector in table_inputs.sectors:
        assert sector.source.schema_id == PAPER_0710_1869_SOURCE_SCHEMA_ID
        assert sector.source.paper_id == "arXiv:0710.1869"
        assert sector.source.source_kind == "table"
        assert sector.source.locator_label == "Table I"
        assert sector.source.detail == expected_details[sector.sector_id]
        assert sector.source.citation == "arXiv:0710.1869, Table I"
        assert [
            (entry.c_value, entry.f_value)
            for entry in sector.entries
        ] == expected_pairs[sector.sector_id]


def test_eq3_example_freezes_extracted_values_and_notes() -> None:
    example = default_paper_0710_1869_eq3_example()

    assert example.schema_id == PAPER_0710_1869_EQ3_SCHEMA_ID
    assert example.mode_id == "paper_0710_1869"
    assert example.source.schema_id == PAPER_0710_1869_SOURCE_SCHEMA_ID
    assert example.source.paper_id == "arXiv:0710.1869"
    assert example.source.source_kind == "equation"
    assert example.source.locator_label == "Eq. (3)"
    assert example.source.detail == "example parameters (a, r, theta12, theta23, theta13, delta)"
    assert example.source.citation == "arXiv:0710.1869, Eq. (3)"
    assert "delta is carried as the quoted paper parameter" in (example.source.notes or "")
    assert example.a == 0.8
    assert example.r == 0.3
    assert example.theta12_deg == 115.0
    assert example.theta23_deg == 65.0
    assert example.theta13_deg == 70.0
    assert example.delta == 0.6


def test_pr1_benchmark_wraps_sourced_inputs_without_extra_logic() -> None:
    benchmark = default_paper_0710_1869_pr1_benchmark()
    payload = benchmark.as_dict()

    assert benchmark.schema_id == PAPER_0710_1869_BENCHMARK_SCHEMA_ID
    assert benchmark.mode_id == "paper_0710_1869"
    assert benchmark.paper_id == "arXiv:0710.1869"
    assert benchmark.benchmark_id == "pr1.table_i_eq3_example.v1"
    assert benchmark.label == "table_i_eq3_example"
    assert benchmark.extraction_version == 1
    assert benchmark.status == "sourced_structural_only"
    assert "No fit, EFT matching, RG running, or observable logic" in benchmark.notes
    assert payload["table_i_inputs"]["schema_id"] == PAPER_0710_1869_TABLE_I_SCHEMA_ID
    assert payload["eq3_example"]["schema_id"] == PAPER_0710_1869_EQ3_SCHEMA_ID
    assert payload["table_i_inputs"]["sectors"][1]["entries"][2] == {
        "generation": 3,
        "c_value": -0.06,
        "f_value": 0.8,
    }


def test_pr1_benchmark_catalog_is_closed_and_singleton() -> None:
    benchmarks = paper_0710_1869_pr1_benchmarks()

    assert len(benchmarks) == 1
    assert benchmarks[0] == default_paper_0710_1869_pr1_benchmark()


def test_inputs_and_benchmarks_modules_do_not_use_absolute_or_repo_v1_imports() -> None:
    forbidden_modules = {
        "quarkConstraints",
        "deltaf2",
        "model",
        "couplings",
        "fit",
        "proxies",
        "scan",
        "quarkConstraints.deltaf2",
        "quarkConstraints.model",
        "quarkConstraints.couplings",
        "quarkConstraints.benchmarks",
        "quarkConstraints.fit",
        "quarkConstraints.proxies",
        "quarkConstraints.scan",
        "quarkConstraints.validation",
    }

    assert not module_has_forbidden_import(PACKAGE_ROOT / "inputs.py", forbidden_modules)
    assert not module_has_forbidden_import(PACKAGE_ROOT / "benchmarks.py", forbidden_modules)


def test_importing_benchmark_modules_does_not_load_repo_v1_runtime_modules() -> None:
    script = """
import importlib
import sys

importlib.import_module("quarkConstraints.paper_0710_1869.inputs")
importlib.import_module("quarkConstraints.paper_0710_1869.benchmarks")
forbidden = sorted(
        name
        for name in sys.modules
        if name in {
            "quarkConstraints.benchmarks",
            "quarkConstraints.deltaf2",
            "quarkConstraints.model",
            "quarkConstraints.couplings",
            "quarkConstraints.fit",
            "quarkConstraints.proxies",
            "quarkConstraints.scan",
            "quarkConstraints.validation",
        }
)
if forbidden:
    raise SystemExit(",".join(forbidden))
"""
    completed = subprocess.run(
        [sys.executable, "-c", script],
        check=False,
        capture_output=True,
        text=True,
        cwd=REPO_ROOT,
    )

    assert completed.returncode == 0, completed.stderr or completed.stdout


def test_closed_contract_rejects_wrong_source_paper_id_and_nested_types() -> None:
    with pytest.raises(ValueError, match="paper_id"):
        Paper07101869SourceRef(paper_id="arXiv:0000.0000")

    with pytest.raises(ValueError, match="source must be"):
        Paper07101869Eq3Example(source="bad")  # type: ignore[arg-type]

    with pytest.raises(ValueError, match="sectors must contain only"):
        Paper07101869TableIInputs(sectors=("bad", "bad", "bad"))  # type: ignore[arg-type]

    with pytest.raises(ValueError, match="table_i_inputs must be"):
        default_paper_0710_1869_pr1_benchmark().__class__(
            table_i_inputs="bad"  # type: ignore[arg-type]
        )


def test_acceptance_benchmark_reports_artifact_verifier_status() -> None:
    summary = _run_acceptance_benchmark()
    artifact_summary = summary["artifacts"]
    checks = artifact_summary["checks"]
    lr_summary = summary["lr_contract_freeze"]
    lr_checks = lr_summary["checks"]

    assert artifact_summary["status"] == "ok"
    assert artifact_summary["import_isolation_ok"] is True
    assert artifact_summary["schema_ok"] is True
    assert artifact_summary["numeric_match_ok"] is True
    assert artifact_summary["scope_ok"] is True
    assert set(artifact_summary["bundle_ids"]) == {
        "wilsons",
        "hadronic",
        "observables",
        "provenance",
    }
    assert artifact_summary["point_id"]
    assert artifact_summary["tolerance_policy"] == _canonical_tolerance_policy()
    assert artifact_summary["tolerances"] == _canonical_tolerance_policy()
    assert set(artifact_summary["observable_rows"]) == EXPECTED_ARTIFACT_ROWS
    assert checks["observable_source_wilson_bundle_linked"] is True
    assert checks["observable_source_hadronic_bundle_linked"] is True
    assert checks["hadronic_mu_had_matches_wilson_scale"] is True
    assert checks["hadronic_operator_normalization_matches_wilson"] is True
    assert checks["observable_operator_normalization_matches_wilson"] is True
    assert checks["writer_outputs_are_deterministic"] is True

    assert lr_summary["status"] == "ok"
    assert lr_summary["deterministic"] is True
    assert lr_summary["paper_operator_basis_id"] == PAPER_0710_1869_DELTAF2_OPERATOR_BASIS_ID
    assert lr_summary["paper_lr_operator_names"] == EXPECTED_GUARDED_LR_OPERATORS
    assert lr_summary["paper_lr_definition_ids"] == EXPECTED_GUARDED_LR_DEFINITION_IDS
    assert (
        lr_summary["projector_normalization_id"]
        == PAPER_0710_1869_DELTAF2_PROJECTOR_NORMALIZATION_ID
    )
    assert (
        lr_summary["projector_normalization_note"]
        == PAPER_0710_1869_DELTAF2_PROJECTOR_NORMALIZATION_NOTE
    )
    assert lr_summary["lr_basis_contract_id"] == PAPER_0710_1869_DELTAF2_RG_LR_BASIS_CONTRACT_ID
    _assert_lr_status_semantics(lr_summary["lr_basis_status_id"])
    assert lr_summary["bmu_lr_basis_id"] == PAPER_0710_1869_DELTAF2_RG_LR_BMU_BASIS_ID
    assert lr_summary["mapping_matrix_frozen"] is True
    assert lr_summary["paper_operator_order"] == EXPECTED_LR_PAPER_OPERATOR_ORDER
    assert lr_summary["bmu_lr_operator_order"] == EXPECTED_LR_BMU_OPERATOR_ORDER
    assert lr_summary["paper_to_bmu_operator_map_matrix"] == EXPECTED_PAPER_TO_BMU_OPERATOR_MAP
    assert lr_summary["bmu_to_paper_operator_map_matrix"] == EXPECTED_BMU_TO_PAPER_OPERATOR_MAP
    assert lr_summary["paper_to_bmu_wilson_map_matrix"] == EXPECTED_PAPER_TO_BMU_WILSON_MAP
    assert lr_summary["bmu_to_paper_wilson_map_matrix"] == EXPECTED_BMU_TO_PAPER_WILSON_MAP
    assert (
        lr_summary["matching_lr_support_contract_id"]
        == PAPER_0710_1869_DELTAF2_RG_LR_BASIS_CONTRACT_ID
    )
    _assert_lr_status_semantics(lr_summary["matching_lr_basis_status_id"])
    assert (
        lr_summary["matching_lr_support_status_id"]
        == PAPER_0710_1869_DELTAF2_MATCHING_OBSERVABLE_SUPPORT_STATUS_ID
    )
    assert lr_summary["matching_lr_support_status_id"] != lr_summary["lr_basis_status_id"]
    assert (
        lr_summary["rg_supported_operator_subset_id"]
        == PAPER_0710_1869_DELTAF2_RG_SUPPORTED_OPERATOR_SUBSET_ID
    )
    assert (
        lr_summary["matching_supported_observable_operator_names"]
        == EXPECTED_SUPPORTED_OPERATORS
    )
    assert (
        lr_summary["matching_unsupported_observable_operator_names"]
        == EXPECTED_GUARDED_LR_OPERATORS
    )
    assert lr_summary["rg_supported_operator_names"] == EXPECTED_FULL_RG_SUPPORTED_OPERATORS
    assert lr_summary["rg_unsupported_operator_names"] == []
    assert lr_checks["source_payloads_are_deterministic"] is True
    assert lr_checks["paper_lr_operator_order_is_frozen"] is True
    assert lr_checks["paper_lr_definition_ids_are_explicit"] is True
    assert lr_checks["projector_normalization_id_present"] is True
    assert lr_checks["projector_normalization_note_mentions_projector_convention"] is True
    assert lr_checks["lr_contract_references_paper_basis"] is True
    assert lr_checks["lr_contract_references_paper_lr_definitions"] is True
    assert lr_checks["lr_definitions_are_frozen"] is True
    assert lr_checks["mapping_matrix_is_frozen"] is True
    assert lr_checks["paper_to_bmu_operator_map_is_exact"] is True
    assert lr_checks["bmu_to_paper_operator_map_is_exact"] is True
    assert lr_checks["paper_to_bmu_wilson_map_is_exact"] is True
    assert lr_checks["bmu_to_paper_wilson_map_is_exact"] is True
    assert lr_checks["operator_map_round_trip_is_identity"] is True
    assert lr_checks["wilson_map_round_trip_is_identity"] is True
    assert lr_checks["operator_basis_vectors_map_as_frozen"] is True
    assert lr_checks["wilson_basis_vectors_map_as_frozen"] is True
    assert lr_summary["lr_running_activated"] is True
    assert lr_summary["rg_lr_basis_map_supported"] is True
    assert lr_checks["matching_references_lr_contract"] is True
    assert lr_checks["matching_references_lr_contract_status"] is True
    assert lr_checks["matching_supported_subset_stays_q1_only"] is True
    assert lr_checks["matching_unsupported_subset_is_lr_only"] is True
    assert lr_checks["rg_contract_references_lr_contract"] is True
    assert lr_checks["rg_contract_references_lr_contract_status"] is True
    assert lr_summary["rg_supported_operator_names"] == EXPECTED_FULL_RG_SUPPORTED_OPERATORS
    assert lr_summary["rg_unsupported_operator_names"] == []


def test_acceptance_benchmark_reports_frozen_lr_contract_state() -> None:
    summary = _run_acceptance_benchmark()
    rg_summary = summary["eft_rg_lo"]
    checks = rg_summary["checks"]
    lr_contract_freeze = rg_summary["lr_contract_freeze"]
    lr_running_state = rg_summary["lr_running_state"]

    assert checks["lr_contract_id_present"] is True
    assert checks["lr_contract_status_present"] is True
    assert checks["lr_contract_ids_match_freeze_contract"] is True
    assert checks["lr_contract_status_matches_freeze_contract"] is True
    assert checks["lr_contract_definitions_frozen"] is True
    assert checks["lr_contract_mapping_matrix_is_frozen"] is True
    assert lr_contract_freeze["contract_present"] is True
    assert lr_contract_freeze["status_present"] is True
    assert lr_contract_freeze["definitions_frozen"] is True
    assert lr_contract_freeze["mapping_matrix_frozen"] is True
    assert lr_contract_freeze["running_active"] is True
    assert lr_contract_freeze["contract_matches_rg_export"] is True
    assert lr_contract_freeze["status_matches_rg_export"] is True
    assert lr_contract_freeze["paper_operator_order"] == EXPECTED_LR_PAPER_OPERATOR_ORDER
    assert lr_contract_freeze["bmu_lr_operator_order"] == EXPECTED_LR_BMU_OPERATOR_ORDER
    assert lr_running_state["contract_id"] == PAPER_0710_1869_DELTAF2_RG_LR_BASIS_CONTRACT_ID
    _assert_lr_status_semantics(lr_running_state["status_id"])
    assert (
        lr_running_state["supported_operator_subset_id"]
        == PAPER_0710_1869_DELTAF2_RG_SUPPORTED_OPERATOR_SUBSET_ID
    )
    assert lr_running_state["lr_running_activated"] is True
    assert lr_running_state["lr_basis_map_supported"] is True
    assert lr_running_state["supported_operator_names"] == EXPECTED_FULL_RG_SUPPORTED_OPERATORS
    assert lr_running_state["unsupported_operator_names"] == []
    assert (
        lr_contract_freeze["paper_to_bmu_operator_map_matrix"]
        == EXPECTED_PAPER_TO_BMU_OPERATOR_MAP
    )
    assert (
        lr_contract_freeze["bmu_to_paper_operator_map_matrix"]
        == EXPECTED_BMU_TO_PAPER_OPERATOR_MAP
    )
    assert (
        lr_contract_freeze["paper_to_bmu_wilson_map_matrix"]
        == EXPECTED_PAPER_TO_BMU_WILSON_MAP
    )
    assert (
        lr_contract_freeze["bmu_to_paper_wilson_map_matrix"]
        == EXPECTED_BMU_TO_PAPER_WILSON_MAP
    )
    assert lr_contract_freeze["operator_map_round_trip_is_identity"] is True
    assert lr_contract_freeze["wilson_map_round_trip_is_identity"] is True
    assert lr_contract_freeze["operator_basis_vectors_map_as_frozen"] is True
    assert lr_contract_freeze["wilson_basis_vectors_map_as_frozen"] is True
    assert len(lr_contract_freeze["paper_operator_definition_ids"]) == 2
    assert all(
        isinstance(item, str) and item
        for item in lr_contract_freeze["paper_operator_definition_ids"]
    )


def test_acceptance_benchmark_reports_synthetic_lr_probe() -> None:
    summary = _run_acceptance_benchmark()
    probe = summary["eft_rg_lr_probe"]
    checks = probe["checks"]

    assert probe["status"] == "ok"
    assert probe["lr_running_activated"] is True
    assert probe["lr_basis_map_supported"] is True
    assert probe["public_run_error"] is None
    _assert_complex_vector_payload(
        probe["probe_input_paper_lr_wilsons"],
        EXPECTED_SYNTHETIC_LR_INPUT,
    )
    _assert_complex_vector_payload(
        probe["mapped_input_bmu_lr_wilsons"],
        EXPECTED_SYNTHETIC_BMU_INPUT,
    )
    assert len(probe["evolved_bmu_lr_wilsons"]) == 2
    assert len(probe["mapped_back_paper_lr_wilsons"]) == 2
    assert len(probe["public_evolved_paper_lr_wilsons"]) == 2
    assert checks["segment_product_matches_bmu_total"] is True
    assert checks["paper_block_matches_winv_u_bmu_w"] is True
    assert checks["mapped_back_vector_matches_conjugated_block"] is True
    assert checks["public_rg_accepts_nonzero_lr"] is True
    assert checks["public_rg_matches_winv_u_bmu_w"] is True


def test_acceptance_benchmark_reports_custom_lr_hadronic_probe() -> None:
    summary = _run_acceptance_benchmark()
    hadronic_summary = summary["hadronic_inputs"]
    probe = summary["custom_lr_hadronic_probe"]
    checks = probe["checks"]

    assert hadronic_summary["status"] == "ok"
    assert hadronic_summary["checks"]["supported_operator_subset_is_pr5a"] is True
    assert hadronic_summary["checks"]["unsupported_operator_subset_is_lr_only"] is True

    assert probe["status"] == "ok"
    assert probe["default_bundle_unchanged"] is True
    assert probe["default_lr_fields_absent"] is True
    assert "custom" in str(probe["input_provenance_mode_id"]).lower()
    assert probe["contract_id"]
    assert probe["q4_formula_id"]
    assert probe["q5_formula_id"]
    assert probe["declared_scheme_id"] == probe["scheme_id"]
    assert probe["declared_scheme_id"] != hadronic_summary["scheme_id"]
    assert float(probe["declared_mu_had_GeV"]) == pytest.approx(
        float(probe["mu_had_GeV"]),
        rel=0.0,
        abs=1.0e-12,
    )
    assert float(probe["declared_mu_had_GeV"]) != pytest.approx(
        float(hadronic_summary["mu_had_GeV"]),
        rel=0.0,
        abs=1.0e-12,
    )
    assert float(probe["b4_source_scale_GeV"]) == pytest.approx(
        float(probe["declared_mu_had_GeV"]),
        rel=0.0,
        abs=1.0e-12,
    )
    assert float(probe["b5_source_scale_GeV"]) == pytest.approx(
        float(probe["declared_mu_had_GeV"]),
        rel=0.0,
        abs=1.0e-12,
    )
    assert float(probe["r_chi_source_scale_GeV"]) == pytest.approx(
        float(probe["declared_mu_had_GeV"]),
        rel=0.0,
        abs=1.0e-12,
    )
    assert probe["b4_source_scheme_id"] == probe["declared_scheme_id"]
    assert probe["b5_source_scheme_id"] == probe["declared_scheme_id"]
    assert probe["r_chi_source_scheme_id"] == probe["declared_scheme_id"]

    expected_q4 = (
        2.0
        * float(probe["R_chi_mu_had"])
        * (float(probe["m_K0_GeV"]) ** 2)
        * (float(probe["f_K_GeV"]) ** 2)
        * float(probe["B4_mu_had"])
    )
    expected_q5 = (
        (2.0 / 3.0)
        * float(probe["R_chi_mu_had"])
        * (float(probe["m_K0_GeV"]) ** 2)
        * (float(probe["f_K_GeV"]) ** 2)
        * float(probe["B5_mu_had"])
    )
    assert float(probe["q4_matrix_element_GeV4"]) == pytest.approx(
        expected_q4,
        rel=0.0,
        abs=1.0e-15,
    )
    assert float(probe["q5_matrix_element_GeV4"]) == pytest.approx(
        expected_q5,
        rel=0.0,
        abs=1.0e-15,
    )
    assert checks["default_bundle_unchanged"] is True
    assert checks["default_lr_fields_absent"] is True
    assert checks["custom_mode_only"] is True
    assert checks["contract_id_present"] is True
    assert checks["q4_formula_id_present"] is True
    assert checks["q5_formula_id_present"] is True
    assert checks["declared_scheme_matches_payload"] is True
    assert checks["declared_mu_had_matches_payload"] is True
    assert checks["source_scheme_matches_declared_scheme"] is True
    assert checks["source_scale_matches_declared_mu_had"] is True
    assert checks["q4_matches_bv2004_eq5"] is True
    assert checks["q5_matches_bv2004_eq5"] is True
    assert checks["observables_remain_blocked"] is True


def test_acceptance_benchmark_reports_custom_lr_observable_probe() -> None:
    summary = _run_acceptance_benchmark()
    probe = summary["custom_lr_observable_probe"]
    checks = probe["checks"]

    assert probe["status"] == "ok"
    assert probe["declared_scheme_id"] == probe["scheme_id"]
    assert float(probe["declared_mu_had_GeV"]) == pytest.approx(
        float(probe["mu_had_GeV"]),
        rel=0.0,
        abs=1.0e-12,
    )
    expected_m12 = (
        complex(
            float(probe["probe_q4_lr"]["real"]),
            float(probe["probe_q4_lr"]["imag"]),
        )
        * float(probe["q4_matrix_element_GeV4"])
        + complex(
            float(probe["probe_q5_lr"]["real"]),
            float(probe["probe_q5_lr"]["imag"]),
        )
        * float(probe["q5_matrix_element_GeV4"])
    ) / (2.0 * 0.497611)
    expected_delta_m = 2.0 * float(expected_m12.real)
    _assert_complex_payload(
        probe["M12_K_LR_NP_GeV"],
        {"real": expected_m12.real, "imag": expected_m12.imag},
    )
    assert float(probe["delta_m_K_LR_NP_GeV"]) == pytest.approx(
        expected_delta_m,
        rel=0.0,
        abs=1.0e-24,
    )
    assert checks["declared_scheme_matches_payload"] is True
    assert checks["declared_mu_had_matches_payload"] is True
    assert checks["m12_matches_hand_calculation"] is True
    assert checks["delta_m_matches_hand_calculation"] is True
    assert checks["q4_matrix_element_present"] is True
    assert checks["q5_matrix_element_present"] is True
    assert checks["default_observables_unchanged"] is True
    assert checks["default_artifacts_unchanged"] is True
    assert probe["default_observables_unchanged"] is True
    assert probe["default_artifacts_unchanged"] is True


def test_acceptance_benchmark_reports_custom_combined_observable_probe() -> None:
    summary = _run_acceptance_benchmark()
    probe = summary["custom_combined_observable_probe"]
    checks = probe["checks"]

    assert probe["status"] == "ok"
    assert probe["declared_scheme_id"] == probe["scheme_id"]
    assert float(probe["declared_mu_had_GeV"]) == pytest.approx(
        float(probe["mu_had_GeV"]),
        rel=0.0,
        abs=1.0e-12,
    )
    assert checks["scope_id_is_custom_total"] is True
    assert checks["interpretation_is_custom_total"] is True
    assert checks["observable_ids_are_custom_total"] is True
    assert checks["declared_scheme_matches_payload"] is True
    assert checks["declared_mu_had_matches_payload"] is True
    assert checks["q1_component_matches_hand_calculation"] is True
    assert checks["lr_component_matches_hand_calculation"] is True
    assert checks["total_matches_hand_calculation"] is True
    assert checks["delta_m_matches_hand_calculation"] is True
    assert checks["total_equals_component_sum"] is True
    assert checks["total_matches_existing_q1_plus_lr_surfaces"] is True
    assert checks["custom_q1_bundle_active"] is True
    assert checks["epsilon_k_blocked"] is True
    assert checks["default_observables_unchanged"] is True
    assert checks["default_artifacts_unchanged"] is True
    assert probe["default_observables_unchanged"] is True
    assert probe["default_artifacts_unchanged"] is True

    payload = probe["payload"]
    assert payload["observable_scope_id"] == EXPECTED_CUSTOM_TOTAL_SCOPE_ID
    assert payload["interpretation"] == EXPECTED_CUSTOM_TOTAL_INTERPRETATION_ID
    assert payload["m12_observable_id"] == EXPECTED_CUSTOM_TOTAL_M12_OBSERVABLE_ID
    assert payload["delta_m_observable_id"] == EXPECTED_CUSTOM_TOTAL_DELTA_M_OBSERVABLE_ID

    q1_piece = complex(
        float(probe["M12_K_NP_Q1_GeV"]["real"]),
        float(probe["M12_K_NP_Q1_GeV"]["imag"]),
    )
    lr_piece = complex(
        float(probe["M12_K_LR_NP_GeV"]["real"]),
        float(probe["M12_K_LR_NP_GeV"]["imag"]),
    )
    total_piece = complex(
        float(probe["M12_K_NP_TOTAL_GeV"]["real"]),
        float(probe["M12_K_NP_TOTAL_GeV"]["imag"]),
    )
    assert total_piece == pytest.approx(q1_piece + lr_piece, rel=0.0, abs=1.0e-24)
    assert float(probe["delta_m_K_NP_TOTAL_GeV"]) == pytest.approx(
        2.0 * total_piece.real,
        rel=0.0,
        abs=1.0e-24,
    )


def test_acceptance_benchmark_reports_custom_bd_observable_probe() -> None:
    summary = _run_acceptance_benchmark()
    _assert_custom_b_observable_probe(summary, "B_d")


def test_acceptance_benchmark_reports_custom_bs_observable_probe() -> None:
    summary = _run_acceptance_benchmark()
    _assert_custom_b_observable_probe(summary, "B_s")


def test_acceptance_benchmark_keeps_default_kaon_reference_numbers_fixed() -> None:
    summary = _run_acceptance_benchmark()

    matching_payload = summary["eft_matching"]["payload"]
    rg_payload = summary["eft_rg_lo"]["payload"]
    observable_payload = summary["observables"]["payload"]

    _assert_complex_payload(matching_payload["coefficients"]["Q1_VLL"], EXPECTED_MATCHING_Q1_VLL)
    _assert_complex_payload(rg_payload["coefficients"]["Q1_VLL"], EXPECTED_RG_Q1_VLL)
    _assert_complex_payload(observable_payload["M12_K_NP_GeV"], EXPECTED_M12_K_NP)
    assert "Delta_m_K_NP" in observable_payload["observables"]
    assert float(observable_payload["delta_m_K_NP_GeV"]) == pytest.approx(
        EXPECTED_DELTA_M_K_NP_GEV,
        rel=0.0,
        abs=1.0e-24,
    )


def test_acceptance_benchmark_exports_canonical_artifacts_deterministically(tmp_path: Path) -> None:
    benchmark_dir = tmp_path / "benchmark"
    repeat_dir = tmp_path / "repeat"
    writer_dir = tmp_path / "writer"

    summary = _run_acceptance_benchmark("--export-artifacts-dir", str(benchmark_dir))
    repeat_summary = _run_acceptance_benchmark("--export-artifacts-dir", str(repeat_dir))
    writer_paths = write_default_paper_0710_1869_kaon_artifact_exports(writer_dir)
    artifact_summary = summary["artifacts"]

    assert artifact_summary["status"] == "ok"
    writer_files = {
        "wilsons": writer_paths.wilson_path,
        "hadronic": writer_paths.hadronic_path,
        "observables": writer_paths.observable_path,
        "provenance": writer_paths.provenance_path,
    }
    for key, filename in EXPECTED_ARTIFACT_FILENAMES.items():
        exported = benchmark_dir / filename
        repeated = repeat_dir / filename
        writer_export = writer_files[key]
        golden = GOLDEN_ARTIFACT_DIR / filename
        assert exported.exists(), filename
        assert repeated.exists(), filename
        assert writer_export.exists(), filename
        assert golden.exists(), filename
        assert exported.read_text(encoding="utf-8") == repeated.read_text(encoding="utf-8")
        assert exported.read_text(encoding="utf-8") == writer_export.read_text(encoding="utf-8")
        assert exported.read_text(encoding="utf-8") == golden.read_text(encoding="utf-8")

    assert artifact_summary["file_sha256"] == repeat_summary["artifacts"]["file_sha256"]
    assert artifact_summary["file_sha256"] == artifact_summary["writer_file_sha256"]
