from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
import shutil
from contextlib import contextmanager
from dataclasses import asdict
from pathlib import Path

import pytest

from quarkConstraints.paper_0710_1869.artifacts import (
    write_default_paper_0710_1869_kaon_artifact_exports,
    write_strict_paper_0710_1869_kaon_artifact_exports,
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
STRICT_PAPER_ARTIFACT_ROOT_ENV = (
    "QUARKCONSTRAINTS_PAPER_0710_1869_STRICT_ARTIFACT_ROOT"
)
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
EXPECTED_SUPPORTED_OPERATORS = ["Q1_VLL", "Q1_VRR", "Q4_LR", "Q5_LR"]
EXPECTED_FULL_RG_SUPPORTED_OPERATORS = ["Q1_VLL", "Q1_VRR", "Q4_LR", "Q5_LR"]
EXPECTED_GUARDED_LR_OPERATORS = ["Q4_LR", "Q5_LR"]
EXPECTED_UNSUPPORTED_OPERATORS: list[str] = []
EXPECTED_GUARDED_LR_DEFINITION_IDS = [
    PAPER_0710_1869_DELTAF2_OPERATOR_DEFINITION_Q4_LR_ID,
    PAPER_0710_1869_DELTAF2_OPERATOR_DEFINITION_Q5_LR_ID,
]
EXPECTED_LR_PAPER_OPERATOR_ORDER = ["Q4_LR", "Q5_LR"]
EXPECTED_LR_BMU_OPERATOR_ORDER = ["Q1_LR_BMU", "Q2_LR_BMU"]
EXPECTED_PAPER_TO_BMU_OPERATOR_MAP = [[0.0, -2.0], [1.0, 0.0]]
EXPECTED_BMU_TO_PAPER_OPERATOR_MAP = [[0.0, 1.0], [-0.5, 0.0]]
EXPECTED_PAPER_TO_BMU_WILSON_MAP = [[0.0, -0.5], [1.0, 0.0]]
EXPECTED_BMU_TO_PAPER_WILSON_MAP = [[0.0, 1.0], [-2.0, 0.0]]
EXPECTED_MATCHING_Q1_VLL = {"real": -1.3193697123462574e-11, "imag": -9.113770397032225e-12}
EXPECTED_RG_Q1_VLL = {"real": -9.616540103444602e-12, "imag": -6.642788423632198e-12}
EXPECTED_M12_K_NP = {"real": -2.1490194213852247e-14, "imag": -1.4844716687059954e-14}
EXPECTED_DELTA_M_K_NP_GEV = -4.2980388427704495e-14
EXPECTED_SYNTHETIC_LR_INPUT = [
    {"real": 1.25, "imag": -0.5},
    {"real": -0.75, "imag": 0.25},
]
EXPECTED_SYNTHETIC_BMU_INPUT = [
    {"real": 0.375, "imag": -0.125},
    {"real": 1.25, "imag": -0.5},
]
EXPECTED_CUSTOM_TOTAL_SCOPE_ID = "kaon.np_only.custom_total.q1_plus_lr.v1"
EXPECTED_CUSTOM_TOTAL_INTERPRETATION_ID = "kaon.np_only.custom_total.q1_plus_lr.v1"
EXPECTED_CUSTOM_TOTAL_M12_OBSERVABLE_ID = "M12_K_NP_CUSTOM_TOTAL"
EXPECTED_CUSTOM_TOTAL_DELTA_M_OBSERVABLE_ID = "Delta_m_K_NP_CUSTOM_TOTAL"
EXPECTED_KAON_LR_R_CHI_MASS_SCHEME_ID = "pdg.2024.msbar.running_masses.at_2gev.v1"
EXPECTED_KAON_LR_R_CHI_ACTIVE_FLAVOR_POLICY_ID = (
    "pdg.2024.quark_masses.n_l_4.at_2gev.explicit.v1"
)
EXPECTED_KAON_LR_R_CHI_NO_HIDDEN_CONVERSION_POLICY_ID = (
    "none.freeze_pdg2024_msbar_nl4_inputs_at_2gev.v1"
)
EXPECTED_KAON_LR_R_CHI_M_S_2GEV_GEV = 0.09274
EXPECTED_KAON_LR_R_CHI_M_D_2GEV_GEV = 0.00469
EXPECTED_KAON_LR_R_CHI_M_K0_GEV = 0.497611
EXPECTED_KAON_LR_R_CHI_EXACT_VALUE = 26.085222120747908
EXPECTED_KAON_LR_R_CHI_FREEZE_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.eft_deltaf2.kaon_lr_r_chi_freeze.v1"
)
EXPECTED_KAON_LR_R_CHI_SUMMARY_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.eft_deltaf2.kaon_lr_r_chi_summary.v1"
)
EXPECTED_DEFAULT_LR_HADRONIC_BUNDLE_ID = "hadronic.kaon.lr.default.etm2013_ms_2gev.v1"
EXPECTED_DEFAULT_LR_HADRONIC_SOURCE_ID = (
    "hadronic.kaon.lr.default.etm2013_ms_2gev.aggregate.v1"
)
EXPECTED_DEFAULT_LR_HADRONIC_INPUT_POLICY_ID = (
    "default_source.etm2013.table1.ms_2gev.no_hidden_conversion.v1"
)
EXPECTED_CUSTOM_LR_HADRONIC_INPUT_POLICY_ID = (
    "custom_input_only.default_lr_bundle_frozen_separately.no_auto_consumption.v1"
)
EXPECTED_DEFAULT_LR_HADRONIC_B4_SOURCE_ID = (
    "hadronic.kaon.lr.b4.etm2013.table1.ms_2gev.v1"
)
EXPECTED_DEFAULT_LR_HADRONIC_B5_SOURCE_ID = (
    "hadronic.kaon.lr.b5.etm2013.table1.ms_2gev.v1"
)
EXPECTED_DEFAULT_LR_HADRONIC_R_CHI_SOURCE_ID = "hadronic.kaon.lr.r_chi.pdg2024_msbar_nl4.v1"
EXPECTED_DEFAULT_LR_HADRONIC_MASS_SOURCE_ID = "pdg.2024.k0.mass.v1"
EXPECTED_DEFAULT_LR_HADRONIC_DECAY_CONSTANT_SOURCE_ID = "pdg.2024.fkplus.eq72.14.v1"
EXPECTED_DEFAULT_LR_HADRONIC_ETM_CITATION = (
    "ETM Collaboration, JHEP 03 (2013) 089, arXiv:1207.1287"
)
EXPECTED_DEFAULT_LR_HADRONIC_B4_VALUE = 0.78
EXPECTED_DEFAULT_LR_HADRONIC_B5_VALUE = 0.57
EXPECTED_DEFAULT_LR_HADRONIC_PROVENANCE_IDS = [
    EXPECTED_DEFAULT_LR_HADRONIC_SOURCE_ID,
    EXPECTED_DEFAULT_LR_HADRONIC_B4_SOURCE_ID,
    EXPECTED_DEFAULT_LR_HADRONIC_B5_SOURCE_ID,
    EXPECTED_DEFAULT_LR_HADRONIC_R_CHI_SOURCE_ID,
    EXPECTED_DEFAULT_LR_HADRONIC_MASS_SOURCE_ID,
    EXPECTED_DEFAULT_LR_HADRONIC_DECAY_CONSTANT_SOURCE_ID,
]
EXPECTED_DEFAULT_LR_HADRONIC_SUPPORTED_OPERATORS = ["Q4_LR", "Q5_LR"]
EXPECTED_DEFAULT_LR_HADRONIC_UNSUPPORTED_OPERATORS = ["Q1_VLL", "Q1_VRR"]
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
    "D0": {
        "scope_id": "d0.np_only.custom_q1.m12.v1",
        "interpretation_id": "d0.np_only.custom_q1.v1",
        "m12_id": "M12_D0_NP",
        "delta_m_id": "Delta_m_D0_NP",
        "require_nonzero_probe": True,
    },
}


def _assert_lr_status_semantics(status_id: object) -> None:
    assert isinstance(status_id, str)
    assert status_id == PAPER_0710_1869_DELTAF2_RG_LR_BASIS_STATUS_ID
    assert ".lr_rg_active." in status_id
    assert ".custom_lr_hadronic_active." in status_id
    assert ".custom_lr_only_observable_active." in status_id
    assert ".custom_combined_observable_active." in status_id
    assert status_id.endswith(".default_export_lr_capable.v7")


def _canonical_tolerance_policy() -> dict[str, object]:
    return asdict(TolerancePolicy())


def _run_acceptance_benchmark(
    *args: str,
    strict_canonical_root: Path | None = None,
) -> dict[str, object]:
    env = os.environ.copy()
    if strict_canonical_root is not None:
        env[STRICT_PAPER_ARTIFACT_ROOT_ENV] = str(strict_canonical_root)
    completed = subprocess.run(
        [sys.executable, str(BENCHMARK_SCRIPT), "--emit-json", *args],
        check=False,
        capture_output=True,
        text=True,
        cwd=REPO_ROOT,
        env=env,
    )
    assert completed.stdout, completed.stderr
    return json.loads(completed.stdout)


@contextmanager
def _temporarily_hide_strict_canonical_exports(tmp_path: Path):
    canonical_dir = tmp_path / "strict_paper_kaon"
    if canonical_dir.exists():
        shutil.rmtree(canonical_dir)
    yield canonical_dir


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


def _assert_etm_source_metadata_text(
    metadata_text: str,
    *,
    operator_label: str,
    bag_value: float,
) -> None:
    assert EXPECTED_DEFAULT_LR_HADRONIC_ETM_CITATION in metadata_text
    assert "Table 1" in metadata_text
    assert "Buras" in metadata_text
    assert "2 GeV" in metadata_text
    assert operator_label in metadata_text
    assert f"{bag_value:.2f}" in metadata_text


def _normalized_metadata_text(value: object) -> str:
    lowered = str(value).strip().lower().replace("-", " ").replace("_", " ")
    return " ".join(lowered.split())


def _is_current_custom_lr_input_policy(policy_id: object) -> bool:
    lowered = str(policy_id).strip().lower()
    return bool(lowered) and lowered == EXPECTED_CUSTOM_LR_HADRONIC_INPUT_POLICY_ID.lower()


def _has_current_custom_lr_note_core(notes: object) -> bool:
    lowered = _normalized_metadata_text(notes)
    return bool(lowered) and "defaults not frozen" not in lowered


def _is_current_custom_lr_contract_notes(notes: object) -> bool:
    lowered = _normalized_metadata_text(notes)
    return (
        _has_current_custom_lr_note_core(notes)
        and "frozen default lr hadronic bundle remains separate" in lowered
        and "not auto consumed on this custom surface" in lowered
    )


def _is_current_custom_lr_bundle_notes(notes: object) -> bool:
    lowered = _normalized_metadata_text(notes)
    return (
        _has_current_custom_lr_note_core(notes)
        and "frozen default lr hadronic bundle remains separate" in lowered
        and (
            "lr only and custom combined observable surfaces still "
            "require explicit custom lr inputs"
            in lowered
        )
    )


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
    if bool(expected.get("require_nonzero_probe", False)):
        assert checks["probe_payload_is_nonzero"] is True
        assert abs(expected_m12) > 0.0
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
        == EXPECTED_UNSUPPORTED_OPERATORS
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
    # C-2: LR coefficients are exported and flow through matching instead of staying guarded.
    assert lr_checks["matching_supported_subset_is_q1_lr"] is True
    assert lr_checks["matching_unsupported_subset_is_empty"] is True
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
    assert hadronic_summary["checks"]["default_export_matches_builder_default"] is True
    assert hadronic_summary["checks"]["supported_operator_subset_is_pr5a"] is True
    assert hadronic_summary["checks"]["unsupported_operator_subset_is_lr_only"] is True

    assert probe["status"] == "ok"
    assert probe["default_bundle_unchanged"] is True
    assert probe["default_lr_fields_absent"] is True
    assert "custom" in str(probe["input_provenance_mode_id"]).lower()
    assert _is_current_custom_lr_input_policy(probe["input_policy_id"])
    assert _is_current_custom_lr_input_policy(probe["contract_input_policy_id"])
    assert probe["input_policy_id"] == probe["contract_input_policy_id"]
    assert _is_current_custom_lr_bundle_notes(probe["bundle_notes"])
    assert _is_current_custom_lr_contract_notes(probe["contract_notes"])
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

    # M-27: PL/PR projectors add the audited /4 relative to the old Eq. (5) pin.
    expected_q4 = (
        0.5
        * float(probe["R_chi_mu_had"])
        * (float(probe["m_K0_GeV"]) ** 2)
        * (float(probe["f_K_GeV"]) ** 2)
        * float(probe["B4_mu_had"])
    )
    expected_q5 = (
        (1.0 / 6.0)
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
    assert checks["input_policy_id_present"] is True
    assert checks["contract_input_policy_id_present"] is True
    assert checks["input_policy_id_matches_contract"] is True
    assert (
        checks["input_policy_id_is_current_custom_frozen_default_bundle_no_auto_consumption"]
        is True
    )
    assert checks["bundle_notes_match_current_custom_lr_contract"] is True
    assert checks["contract_notes_match_current_custom_lr_contract"] is True
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
    assert checks["default_q1_only_observable_surface_rejects_custom_lr_bundle"] is True


def test_acceptance_benchmark_reports_default_lr_hadronic_probe() -> None:
    summary = _run_acceptance_benchmark()
    hadronic_summary = summary["hadronic_inputs"]
    probe = summary["default_lr_hadronic_probe"]
    checks = probe["checks"]

    assert probe["status"] == "ok"
    assert probe["deterministic"] is True
    assert probe["cross_process_deterministic"] is True
    assert probe["payload_sha256"] == _canonical_payload_sha256(probe["payload"])
    assert probe["same_process_repeat_payload_sha256"] == probe["payload_sha256"]
    assert probe["cross_process_payload_sha256"] == probe["payload_sha256"]
    assert probe["system_id"] == "kaon"
    assert float(probe["mu_had_GeV"]) == pytest.approx(2.0, rel=0.0, abs=1.0e-12)
    assert probe["bundle_id"] == EXPECTED_DEFAULT_LR_HADRONIC_BUNDLE_ID
    assert probe["source_id"] == EXPECTED_DEFAULT_LR_HADRONIC_SOURCE_ID
    assert probe["input_policy_id"] == EXPECTED_DEFAULT_LR_HADRONIC_INPUT_POLICY_ID
    assert probe["provenance_ids"] == EXPECTED_DEFAULT_LR_HADRONIC_PROVENANCE_IDS
    assert (
        probe["payload"]["supported_operator_names"]
        == EXPECTED_DEFAULT_LR_HADRONIC_SUPPORTED_OPERATORS
    )
    assert (
        probe["payload"]["unsupported_operator_names"]
        == EXPECTED_DEFAULT_LR_HADRONIC_UNSUPPORTED_OPERATORS
    )
    assert float(probe["B4_mu_had"]) == pytest.approx(
        EXPECTED_DEFAULT_LR_HADRONIC_B4_VALUE,
        rel=0.0,
        abs=1.0e-15,
    )
    assert float(probe["B5_mu_had"]) == pytest.approx(
        EXPECTED_DEFAULT_LR_HADRONIC_B5_VALUE,
        rel=0.0,
        abs=1.0e-15,
    )
    assert float(probe["R_chi_mu_had"]) == pytest.approx(
        EXPECTED_KAON_LR_R_CHI_EXACT_VALUE,
        rel=0.0,
        abs=1.0e-15,
    )
    assert probe["b4_source_id"] == EXPECTED_DEFAULT_LR_HADRONIC_B4_SOURCE_ID
    assert probe["b5_source_id"] == EXPECTED_DEFAULT_LR_HADRONIC_B5_SOURCE_ID
    assert probe["r_chi_source_id"] == EXPECTED_DEFAULT_LR_HADRONIC_R_CHI_SOURCE_ID
    assert probe["mass_source_id"] == EXPECTED_DEFAULT_LR_HADRONIC_MASS_SOURCE_ID
    assert (
        probe["decay_constant_source_id"]
        == EXPECTED_DEFAULT_LR_HADRONIC_DECAY_CONSTANT_SOURCE_ID
    )
    assert probe["b4_source_scheme_id"] == probe["scheme_id"]
    assert probe["b5_source_scheme_id"] == probe["scheme_id"]
    assert probe["r_chi_source_scheme_id"] == EXPECTED_KAON_LR_R_CHI_MASS_SCHEME_ID
    assert probe["r_chi_source_scheme_id"] != probe["scheme_id"]
    assert float(probe["b4_source_scale_GeV"]) == pytest.approx(
        float(probe["mu_had_GeV"]),
        rel=0.0,
        abs=1.0e-12,
    )
    assert float(probe["b5_source_scale_GeV"]) == pytest.approx(
        float(probe["mu_had_GeV"]),
        rel=0.0,
        abs=1.0e-12,
    )
    assert float(probe["r_chi_source_scale_GeV"]) == pytest.approx(
        float(probe["mu_had_GeV"]),
        rel=0.0,
        abs=1.0e-12,
    )
    _assert_etm_source_metadata_text(
        str(probe["b4_source_metadata_text"]),
        operator_label="B4",
        bag_value=EXPECTED_DEFAULT_LR_HADRONIC_B4_VALUE,
    )
    _assert_etm_source_metadata_text(
        str(probe["b5_source_metadata_text"]),
        operator_label="B5",
        bag_value=EXPECTED_DEFAULT_LR_HADRONIC_B5_VALUE,
    )

    # M-27: PL/PR projectors add the audited /4 relative to the old Eq. (5) pin.
    expected_q4 = (
        0.5
        * float(probe["R_chi_mu_had"])
        * (float(probe["m_K0_GeV"]) ** 2)
        * (float(probe["f_K_GeV"]) ** 2)
        * float(probe["B4_mu_had"])
    )
    expected_q5 = (
        (1.0 / 6.0)
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
    assert hadronic_summary["status"] == "ok"
    assert hadronic_summary["checks"]["default_export_matches_builder_default"] is True
    assert hadronic_summary["checks"]["supported_operator_subset_is_pr5a"] is True
    assert hadronic_summary["checks"]["unsupported_operator_subset_is_lr_only"] is True
    for key in (
        "default_export_matches_builder_default",
        "operator_basis_id_is_exact",
        "operator_normalization_id_is_exact",
        "hamiltonian_convention_id_is_exact",
        "q4_formula_id_is_exact",
        "q5_formula_id_is_exact",
        "b4_source_transformation_id_is_none",
        "b5_source_transformation_id_is_none",
        "contract_operator_basis_id_is_exact",
        "contract_operator_normalization_id_is_exact",
        "contract_hamiltonian_convention_id_is_exact",
    ):
        assert checks[key] is True
    assert checks["probe_payload_is_deterministic"] is True
    assert checks["probe_payload_is_cross_process_deterministic"] is True
    assert checks["probe_payload_sha256_matches_repeat"] is True
    assert checks["probe_payload_sha256_matches_cross_process"] is True
    assert checks["system_is_kaon"] is True
    assert checks["mu_had_is_2gev"] is True
    assert checks["bundle_id_is_exact"] is True
    assert checks["source_id_is_exact"] is True
    assert checks["input_policy_id_is_exact"] is True
    assert checks["provenance_ids_match_frozen_source_package"] is True
    assert checks["supported_operator_subset_is_lr_only"] is True
    assert checks["unsupported_operator_subset_is_q1_only"] is True
    assert checks["b4_value_matches_etm2013_table1"] is True
    assert checks["b5_value_matches_etm2013_table1"] is True
    assert checks["r_chi_matches_frozen_default"] is True
    assert checks["b4_source_id_is_exact"] is True
    assert checks["b5_source_id_is_exact"] is True
    assert checks["r_chi_source_id_is_exact"] is True
    assert checks["mass_source_id_is_exact"] is True
    assert checks["decay_constant_source_id_is_exact"] is True
    assert checks["b4_source_scheme_matches_bundle"] is True
    assert checks["b5_source_scheme_matches_bundle"] is True
    assert checks["r_chi_source_scheme_matches_frozen_mass_scheme"] is True
    assert checks["r_chi_source_scheme_stays_distinct_from_bundle_scheme"] is True
    assert checks["b4_source_scale_matches_bundle"] is True
    assert checks["b5_source_scale_matches_bundle"] is True
    assert checks["r_chi_source_scale_matches_bundle"] is True
    assert checks["b4_source_mentions_etm2013_table1_ms_2gev"] is True
    assert checks["b5_source_mentions_etm2013_table1_ms_2gev"] is True
    assert checks["q4_matches_bv2004_eq5"] is True
    assert checks["q5_matches_bv2004_eq5"] is True
    assert checks["drifted_b4_value_rejected"] is True
    assert checks["drifted_b5_value_rejected"] is True
    assert checks["drifted_source_id_rejected"] is True
    assert checks["drifted_b4_source_id_rejected"] is True
    assert checks["drifted_b5_source_id_rejected"] is True
    assert checks["drifted_r_chi_source_id_rejected"] is True
    assert checks["drifted_bundle_scheme_rejected"] is True
    assert checks["forced_r_chi_source_bundle_scheme_rejected"] is True
    assert checks["drifted_mu_had_rejected"] is True
    assert checks["hidden_conversion_policy_drift_rejected"] is True
    assert checks["hidden_r_chi_conversion_drift_rejected"] is True
    assert checks["custom_lr_only_rejects_default_lr_bundle"] is True
    assert checks["custom_combined_rejects_default_lr_bundle"] is True
    assert checks["custom_lr_surfaces_still_require_custom_inputs"] is True
    assert checks["default_q1_hadronic_bundle_stays_q1_only"] is True
    assert checks["default_observables_unchanged"] is True
    assert checks["default_artifacts_unchanged"] is True
    assert "input_provenance_mode_id" in str(probe["custom_lr_only_rejection_error"])
    assert "input_provenance_mode_id" in str(probe["custom_combined_rejection_error"])
    assert probe["default_q1_hadronic_bundle_stays_q1_only"] is True
    assert probe["default_observables_unchanged"] is True
    assert probe["default_artifacts_unchanged"] is True


def test_acceptance_benchmark_reports_lr_r_chi_freeze_probe() -> None:
    summary = _run_acceptance_benchmark()
    hadronic_summary = summary["hadronic_inputs"]
    probe = summary["lr_r_chi_freeze_probe"]
    checks = probe["checks"]

    assert probe["status"] == "ok"
    assert probe["deterministic"] is True
    assert probe["cross_process_deterministic"] is True
    assert probe["payload_sha256"] == _canonical_payload_sha256(probe["payload"])
    assert probe["same_process_repeat_payload_sha256"] == probe["payload_sha256"]
    assert probe["cross_process_payload_sha256"] == probe["payload_sha256"]
    assert probe["system_id"] == "kaon"
    assert float(probe["mu_had_GeV"]) == pytest.approx(2.0, rel=0.0, abs=1.0e-12)
    assert probe["freeze_id"]
    assert probe["source_id"]
    assert "derived" in str(probe["input_provenance_mode_id"]).lower()
    assert "freeze_only" in str(probe["input_policy_id"]).lower()
    assert probe["object_schema_id"] == EXPECTED_KAON_LR_R_CHI_FREEZE_SCHEMA_ID
    assert probe["summary_schema_id"] == EXPECTED_KAON_LR_R_CHI_SUMMARY_SCHEMA_ID
    assert probe["operator_scheme_id"] == probe["operator_renormalization_scheme_id"]
    assert probe["mass_scheme_id"] == probe["mass_renormalization_scheme_id"]
    assert probe["mass_scheme_id"] == EXPECTED_KAON_LR_R_CHI_MASS_SCHEME_ID
    assert probe["operator_scheme_id"] != probe["mass_scheme_id"]
    assert (
        probe["mass_active_flavor_policy_id"]
        == EXPECTED_KAON_LR_R_CHI_ACTIVE_FLAVOR_POLICY_ID
    )
    assert probe["derivation_formula_id"]
    assert probe["derivation_formula_source_id"]
    assert (
        probe["no_hidden_conversion_policy_id"]
        == EXPECTED_KAON_LR_R_CHI_NO_HIDDEN_CONVERSION_POLICY_ID
    )
    assert float(probe["m_K0_GeV"]) == pytest.approx(
        EXPECTED_KAON_LR_R_CHI_M_K0_GEV,
        rel=0.0,
        abs=1.0e-15,
    )
    assert float(probe["m_s_mu_had_GeV"]) == pytest.approx(
        EXPECTED_KAON_LR_R_CHI_M_S_2GEV_GEV,
        rel=0.0,
        abs=1.0e-15,
    )
    assert float(probe["m_d_mu_had_GeV"]) == pytest.approx(
        EXPECTED_KAON_LR_R_CHI_M_D_2GEV_GEV,
        rel=0.0,
        abs=1.0e-15,
    )
    assert probe["provenance_ids"] == [
        probe["source_id"],
        probe["kaon_mass_source_id"],
        probe["strange_mass_source_id"],
        probe["down_mass_source_id"],
    ]
    assert probe["default_lr_hadronic_available"] is True
    assert any(probe["default_lr_hadronic_exports_present"].values()) is True
    assert hadronic_summary["status"] == "ok"
    assert hadronic_summary["checks"]["default_export_matches_builder_default"] is True
    assert hadronic_summary["checks"]["supported_operator_subset_is_pr5a"] is True
    assert hadronic_summary["checks"]["unsupported_operator_subset_is_lr_only"] is True
    assert float(probe["strange_mass_source_scale_GeV"]) == pytest.approx(
        float(probe["mu_had_GeV"]),
        rel=0.0,
        abs=1.0e-12,
    )
    assert float(probe["down_mass_source_scale_GeV"]) == pytest.approx(
        float(probe["mu_had_GeV"]),
        rel=0.0,
        abs=1.0e-12,
    )
    expected_r_chi = (
        float(probe["m_K0_GeV"])
        / (float(probe["m_s_mu_had_GeV"]) + float(probe["m_d_mu_had_GeV"]))
    ) ** 2
    assert expected_r_chi == pytest.approx(
        EXPECTED_KAON_LR_R_CHI_EXACT_VALUE,
        rel=0.0,
        abs=1.0e-15,
    )
    assert float(probe["R_chi_mu_had"]) == pytest.approx(
        EXPECTED_KAON_LR_R_CHI_EXACT_VALUE,
        rel=0.0,
        abs=1.0e-15,
    )
    assert checks["probe_payload_is_deterministic"] is True
    assert checks["probe_payload_is_cross_process_deterministic"] is True
    assert checks["probe_payload_sha256_matches_repeat"] is True
    assert checks["probe_payload_sha256_matches_cross_process"] is True
    assert checks["object_schema_is_exact"] is True
    assert checks["summary_schema_is_exact"] is True
    assert checks["object_and_summary_schema_are_distinct"] is True
    assert checks["system_is_kaon"] is True
    assert checks["mass_scale_is_2gev"] is True
    assert checks["formula_matches_bv2004"] is True
    assert checks["mass_sources_frozen"] is True
    assert checks["mass_sources_at_declared_scale"] is True
    assert checks["mass_scheme_is_explicit"] is True
    assert checks["mass_scheme_id_is_exact"] is True
    assert checks["active_flavor_policy_is_frozen"] is True
    assert checks["active_flavor_policy_is_explicit_n_l_4"] is True
    assert checks["derivation_formula_frozen"] is True
    assert checks["provenance_mode_is_derived_default"] is True
    assert checks["input_policy_is_r_chi_freeze_only"] is True
    assert checks["no_hidden_conversion_policy"] is True
    assert checks["no_hidden_conversion_policy_is_explicit"] is True
    assert checks["m_s_mu_had_matches_exact_freeze"] is True
    assert checks["m_d_mu_had_matches_exact_freeze"] is True
    assert checks["m_k0_gev_matches_exact_freeze"] is True
    assert checks["r_chi_mu_had_matches_exact_freeze"] is True
    assert checks["operator_vs_mass_scheme_not_aliased"] is True
    assert checks["custom_lr_only_rejects_r_chi_freeze"] is True
    assert checks["custom_combined_rejects_r_chi_freeze"] is True
    assert checks["custom_lr_only_rejection_message_is_explicit_type_guard"] is True
    assert checks["custom_combined_rejection_message_is_explicit_type_guard"] is True
    assert checks["custom_lr_surfaces_still_require_custom_inputs"] is True
    for key in (
        "default_hadronic_bundle_stays_q1_only",
        "default_observables_unchanged",
        "default_artifacts_unchanged",
        "default_observable_values_match_frozen_reference",
        "default_artifact_exports_match_tracked_files",
    ):
        assert checks[key] is True
    assert (
        probe["custom_lr_only_rejection_error"]
        == "ValueError: hadronic_inputs must be a Paper07101869KaonLRHadronicInputs"
    )
    assert (
        probe["custom_combined_rejection_error"]
        == "ValueError: lr_hadronic_inputs must be a Paper07101869KaonLRHadronicInputs"
    )
    assert probe["default_hadronic_bundle_stays_q1_only"] is True
    assert probe["default_observables_unchanged"] is True
    assert probe["default_artifacts_unchanged"] is True


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


def test_acceptance_benchmark_reports_custom_d0_observable_probe() -> None:
    summary = _run_acceptance_benchmark()
    _assert_custom_b_observable_probe(summary, "D0")


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


def test_acceptance_benchmark_reports_strict_paper_artifact_summary_and_exports(
    tmp_path: Path,
) -> None:
    benchmark_dir = tmp_path / "benchmark"
    repeat_dir = tmp_path / "repeat"
    writer_dir = tmp_path / "writer"

    with _temporarily_hide_strict_canonical_exports(tmp_path) as canonical_dir:
        assert not canonical_dir.exists()

        summary = _run_acceptance_benchmark(
            "--export-artifacts-dir",
            str(benchmark_dir),
            strict_canonical_root=canonical_dir,
        )
        strict_summary = summary.get("strict_paper_artifacts")
        if not strict_summary:
            pytest.fail("strict-paper artifact summary path not exposed yet")

        assert summary["artifacts"]["status"] == "ok"
        assert strict_summary["status"] == "ok"
        assert strict_summary["canonical_root_source"] == "env_override"
        assert strict_summary["canonical_root_path"] == str(canonical_dir)
        checks = strict_summary["checks"]
        assert checks["point_id_matches_strict_paper_contract"] is True
        assert checks["bundle_ids_match_strict_paper_contract"] is True
        assert checks["strict_numerics_match_default_export"] is True
        assert checks["strict_provenance_ids_are_disjoint_from_default"] is True
        assert checks["strict_source_bundle_ids_are_distinct_from_default"] is True
        assert checks["canonical_export_files_present"] is False
        assert checks["canonical_export_tree_is_complete_or_absent"] is True
        assert checks["effective_strict_exports_match_current_export"] is False
        assert strict_summary["canonical_file_sha256"] == {}
        assert (
            strict_summary["matching_coefficients"]
            == summary["artifacts"]["matching_coefficients"]
        )

        write_strict_paper_0710_1869_kaon_artifact_exports(canonical_dir)

        corrupted_provenance_file = canonical_dir / EXPECTED_ARTIFACT_FILENAMES["provenance"]
        corrupted_provenance_file.unlink()

        partial_summary = _run_acceptance_benchmark(
            "--export-artifacts-dir",
            str(repeat_dir),
            strict_canonical_root=canonical_dir,
        )
        partial_strict_summary = partial_summary.get("strict_paper_artifacts")
        if not partial_strict_summary:
            pytest.fail(
                "strict-paper artifact summary path not exposed after partial canonical write"
            )

        partial_checks = partial_strict_summary["checks"]
        assert partial_strict_summary["canonical_root_source"] == "env_override"
        assert partial_strict_summary["canonical_root_path"] == str(canonical_dir)
        assert partial_summary["strict_paper_artifacts"]["status"] == "failed"
        assert partial_checks["canonical_export_files_present"] is False
        assert partial_checks["canonical_export_tree_is_complete_or_absent"] is False
        assert partial_checks["effective_strict_exports_match_current_export"] is False

        write_strict_paper_0710_1869_kaon_artifact_exports(canonical_dir)
        repeat_summary = _run_acceptance_benchmark(
            "--export-artifacts-dir",
            str(repeat_dir),
            strict_canonical_root=canonical_dir,
        )
        repeat_strict_summary = repeat_summary.get("strict_paper_artifacts")
        if not repeat_strict_summary:
            pytest.fail("strict-paper artifact summary path not exposed after canonical write")

        repeat_checks = repeat_strict_summary["checks"]
        assert repeat_strict_summary["canonical_root_source"] == "env_override"
        assert repeat_strict_summary["canonical_root_path"] == str(canonical_dir)
        assert repeat_strict_summary["repo_tracked_canonical_root_path"] == str(
            REPO_ROOT / "results" / "paper_0710_1869" / "strict_paper_kaon"
        )
        assert repeat_checks["canonical_export_files_present"] is True
        assert repeat_checks["effective_strict_exports_match_current_export"] is True
        assert repeat_checks["canonical_export_tree_is_complete_or_absent"] is True
        assert repeat_strict_summary["canonical_file_sha256"] == repeat_strict_summary["file_sha256"]
        assert (
            repeat_strict_summary["matching_coefficients"]
            == repeat_summary["artifacts"]["matching_coefficients"]
        )

    writer_paths = write_strict_paper_0710_1869_kaon_artifact_exports(writer_dir)
    strict_export_dir = benchmark_dir / "strict_paper"
    repeat_strict_export_dir = repeat_dir / "strict_paper"
    writer_files = {
        "wilsons": writer_paths.wilson_path,
        "hadronic": writer_paths.hadronic_path,
        "observables": writer_paths.observable_path,
        "provenance": writer_paths.provenance_path,
    }
    for key, filename in EXPECTED_ARTIFACT_FILENAMES.items():
        exported = strict_export_dir / filename
        repeated = repeat_strict_export_dir / filename
        writer_export = writer_files[key]
        assert exported.exists(), filename
        assert repeated.exists(), filename
        assert writer_export.exists(), filename
        assert exported.read_text(encoding="utf-8") == repeated.read_text(encoding="utf-8")
        assert exported.read_text(encoding="utf-8") == writer_export.read_text(encoding="utf-8")

    assert strict_summary["file_sha256"] == repeat_strict_summary["file_sha256"]
    assert strict_summary["file_sha256"] == strict_summary["writer_file_sha256"]
