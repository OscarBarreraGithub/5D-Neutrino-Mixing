"""Contract-first tests for the paper-mode coupling slice.

These tests are written to be forward-compatible while the paper-mode coupling
API is still landing. They are strict once the coupling module exists.
"""

from __future__ import annotations

import importlib
import inspect
import json
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]

FORBIDDEN_REPO_V1_MODULES = {
    "quarkConstraints.benchmarks",
    "quarkConstraints.couplings",
    "quarkConstraints.deltaf2",
    "quarkConstraints.fit",
    "quarkConstraints.model",
    "quarkConstraints.proxies",
    "quarkConstraints.scan",
    "quarkConstraints.scales",
    "quarkConstraints.validation",
}


def _has_paper_couplings_module() -> bool:
    return (REPO_ROOT / "quarkConstraints" / "paper_0710_1869" / "couplings.py").exists()


def _load_couplings_module():
    if not _has_paper_couplings_module():
        pytest.skip("paper_0710_1869 coupling slice not implemented yet")
    return __import__("quarkConstraints.paper_0710_1869.couplings", fromlist=["*"])


def _get_contract_summary_callable(module):
    for name in (
        "coupling_contract_summary",
        "build_coupling_contract_summary",
        "paper_0710_1869_coupling_contract_summary",
    ):
        candidate = getattr(module, name, None)
        if callable(candidate):
            return candidate
    return None


def _load_kkgluon_module():
    if not _has_paper_couplings_module():
        pytest.skip("paper_0710_1869 coupling slice not implemented yet")
    return __import__("quarkConstraints.paper_0710_1869.kkgluon", fromlist=["*"])


def _physical_bulk_state() -> object:
    benchmarks_module = importlib.import_module("quarkConstraints.paper_0710_1869.benchmarks")
    model_module = importlib.import_module("quarkConstraints.paper_0710_1869.model")

    seed = benchmarks_module.Paper07101869BenchmarkSpurionSeed(
        up_singular_values=(0.45, 1.2, 3.4),
        down_singular_values=(0.15, 0.55, 1.7),
        overall_scale=0.8,
        up_left=model_module.Paper07101869RotationParameters.from_degrees(
            theta12_deg=19.0,
            theta13_deg=4.0,
            theta23_deg=11.0,
            delta=0.3,
        ),
        up_right=model_module.Paper07101869RotationParameters.from_degrees(
            theta12_deg=13.0,
            theta13_deg=7.0,
            theta23_deg=5.0,
            delta=0.1,
        ),
        down_left=model_module.Paper07101869RotationParameters.from_degrees(
            theta12_deg=7.0,
            theta13_deg=2.0,
            theta23_deg=9.0,
            delta=0.5,
        ),
        down_right=model_module.Paper07101869RotationParameters.from_degrees(
            theta12_deg=17.0,
            theta13_deg=3.0,
            theta23_deg=6.0,
            delta=0.2,
        ),
        notes="couplings-test-seed",
    )
    benchmark = benchmarks_module.default_paper_0710_1869_pr1_benchmark()
    physical_point = benchmarks_module.build_paper_0710_1869_seeded_physical_point(
        benchmark,
        seed,
        metadata={"test_case": "couplings_acceptance"},
    )
    return model_module.derive_paper_0710_1869_physical_bulk_state(physical_point)


def _physical_kk_gluon_adapter_callable(module):
    for name in (
        "build_paper_0710_1869_physical_kk_gluon_couplings",
        "build_paper_0710_1869_kk_gluon_couplings_from_physical_bulk_state",
        "build_paper_0710_1869_kk_gluon_couplings_from_physical_state",
    ):
        candidate = getattr(module, name, None)
        if callable(candidate):
            return candidate

    candidate = getattr(module, "build_paper_0710_1869_kk_gluon_couplings", None)
    if callable(candidate):
        try:
            parameters = inspect.signature(candidate).parameters
        except (TypeError, ValueError):
            parameters = {}
        if any(name in parameters for name in ("physical_bulk_state", "bulk_state")):
            return candidate
    return None


def _invoke_physical_kk_gluon_adapter(callable_obj, physical_bulk_state):
    try:
        parameters = inspect.signature(callable_obj).parameters
    except (TypeError, ValueError):
        parameters = {}
    for name in ("physical_bulk_state", "bulk_state", "physical_state"):
        if name in parameters:
            return callable_obj(**{name: physical_bulk_state})
    raise AssertionError(
        "physical bulk-state adapter exists but exposes no supported keyword parameter"
    )


def test_importing_paper_couplings_does_not_load_repo_v1_modules() -> None:
    if not _has_paper_couplings_module():
        pytest.skip("paper_0710_1869 coupling slice not implemented yet")

    script = f"""
import importlib
import json
import sys

importlib.import_module("quarkConstraints.paper_0710_1869.couplings")

forbidden = sorted(
    name for name in sys.modules
    if name in {sorted(FORBIDDEN_REPO_V1_MODULES)!r}
)
print(json.dumps(forbidden))
"""
    completed = subprocess.run(
        [sys.executable, "-c", script],
        check=True,
        capture_output=True,
        text=True,
        cwd=str(REPO_ROOT),
    )
    assert json.loads(completed.stdout) == []


def test_importing_paper_kkgluon_does_not_load_repo_v1_modules() -> None:
    if not _has_paper_couplings_module():
        pytest.skip("paper_0710_1869 coupling slice not implemented yet")

    script = f"""
import importlib
import json
import sys

importlib.import_module("quarkConstraints.paper_0710_1869.kkgluon")

forbidden = sorted(
    name for name in sys.modules
    if name in {sorted(FORBIDDEN_REPO_V1_MODULES)!r}
)
print(json.dumps(forbidden))
"""
    completed = subprocess.run(
        [sys.executable, "-c", script],
        check=True,
        capture_output=True,
        text=True,
        cwd=str(REPO_ROOT),
    )
    assert json.loads(completed.stdout) == []


def test_couplings_contract_summary_exists_and_is_deterministic() -> None:
    module = _load_couplings_module()
    summary_fn = _get_contract_summary_callable(module)
    assert summary_fn is not None, (
        "paper coupling slice exists but exposes no contract summary callable; "
        "expected one of: coupling_contract_summary / build_coupling_contract_summary / "
        "paper_0710_1869_coupling_contract_summary"
    )

    first = summary_fn()
    second = summary_fn()
    assert json.dumps(first, sort_keys=True) == json.dumps(second, sort_keys=True)


def test_mu_gs_semantics_and_basis_no_fcnc_under_universal_profiles() -> None:
    """Run behavioral coupling checks if the coupling slice exposes the expected API.

    This test is intentionally conditional: the coupling slice can land with only a
    contract summary first, and add the heavier coupling machinery next.
    """

    module = _load_couplings_module()

    normalize_fn = getattr(module, "kk_gluon_normalization", None)
    coupling_fn = getattr(module, "kk_gluon_coupling_matrix", None)
    if not callable(normalize_fn) or not callable(coupling_fn):
        pytest.skip("behavioral coupling API not exposed yet (expected normalization + matrix)")

    scales_mod = __import__("quarkConstraints.paper_0710_1869.scales", fromlist=["*"])
    Paper07101869ScalePoint = scales_mod.Paper07101869ScalePoint

    # mu_gs semantics: normalization must depend on mu_gs, not mu_match, when mu_gs differs.
    base = Paper07101869ScalePoint(
        Lambda_IR_GeV=3000.0,
        m_g1_GeV=4200.0,
        mu_match_GeV=3000.0,
        mu_gs_GeV=4200.0,
    )
    variant_mu_match = Paper07101869ScalePoint(
        Lambda_IR_GeV=3000.0,
        m_g1_GeV=4200.0,
        mu_match_GeV=2500.0,
        mu_gs_GeV=4200.0,
    )
    variant_mu_gs = Paper07101869ScalePoint(
        Lambda_IR_GeV=3000.0,
        m_g1_GeV=4200.0,
        mu_match_GeV=3000.0,
        mu_gs_GeV=2500.0,
    )

    norm_base = normalize_fn(base)
    norm_mu_match = normalize_fn(variant_mu_match)
    norm_mu_gs = normalize_fn(variant_mu_gs)

    assert norm_base["alpha_s_mu_gs"] == pytest.approx(norm_mu_match["alpha_s_mu_gs"])
    assert norm_base["g_s_mu_gs"] == pytest.approx(norm_mu_match["g_s_mu_gs"])
    assert json.dumps(norm_base, sort_keys=True) == json.dumps(norm_mu_match, sort_keys=True), (
        "kk_gluon_normalization should be controlled by mu_gs, not mu_match"
    )
    assert norm_base["alpha_s_mu_gs"] != pytest.approx(norm_mu_gs["alpha_s_mu_gs"])
    assert norm_base["g_s_mu_gs"] != pytest.approx(norm_mu_gs["g_s_mu_gs"])
    assert json.dumps(norm_base, sort_keys=True) != json.dumps(norm_mu_gs, sort_keys=True), (
        "kk_gluon_normalization should change when mu_gs changes"
    )

    # Universal profiles + identity rotations: coupling matrix must be diagonal (no FCNC).
    profiles = [1.0, 1.0, 1.0]
    identity = [
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ]
    matrix = coupling_fn(profiles=profiles, u_left=identity, u_right=identity, scales=base)
    # Require strict zeros on off-diagonals within floating tolerance.
    off_diag = []
    for i in range(3):
        for j in range(3):
            if i != j:
                off_diag.append(float(matrix[i][j]))
    assert max(abs(value) for value in off_diag) <= 1e-14


def test_default_kkgluon_benchmark_summary_is_present_and_deterministic() -> None:
    module = _load_kkgluon_module()
    summary_fn = getattr(module, "default_paper_0710_1869_kk_gluon_benchmark_summary", None)
    assert callable(summary_fn)

    first = summary_fn().as_dict()
    second = summary_fn().as_dict()
    assert json.dumps(first, sort_keys=True) == json.dumps(second, sort_keys=True)
    assert first["benchmark_status"] == "sourced_structural_only"
    assert first["kk_gluon_normalization_id"] == "explicit_mu_gs.g_s_sqrt_4pi_alpha_s.v1"
    assert first["mu_gs_semantics_id"] == "alpha_s.high_precision.msbar.at_mu_gs_GeV.v1"
    assert first["propagator_mass_rule_id"] == (
        "propagator_mass.scale_point_propagator_mass_GeV.v1"
    )
    assert len(first["matrix_summaries"]) == 4
    matrix_summaries = {item["label"]: item for item in first["matrix_summaries"]}
    assert matrix_summaries["left_up_aligned"]["raw_offdiag_fro_norm"] == pytest.approx(0.0)
    assert matrix_summaries["left_down_aligned"]["raw_offdiag_fro_norm"] > 0.0
    assert matrix_summaries["left_down_aligned"]["max_abs_imag_raw"] > 0.0
    assert matrix_summaries["right_up_diagonal"]["raw_offdiag_fro_norm"] == pytest.approx(0.0)
    assert matrix_summaries["right_down_diagonal"]["raw_offdiag_fro_norm"] == pytest.approx(0.0)
    assert (
        matrix_summaries["right_up_diagonal"]["universal_subtracted_offdiag_fro_norm"]
        == pytest.approx(0.0)
    )
    assert (
        matrix_summaries["right_down_diagonal"]["universal_subtracted_offdiag_fro_norm"]
        == pytest.approx(0.0)
    )


def test_physical_bulk_state_adapter_is_deterministic_and_distinct_from_reference_path() -> None:
    module = _load_kkgluon_module()
    adapter = _physical_kk_gluon_adapter_callable(module)
    if not callable(adapter):
        pytest.skip("physical bulk-state KK-gluon adapter not exposed yet")

    physical_bulk_state = _physical_bulk_state()
    first = _invoke_physical_kk_gluon_adapter(adapter, physical_bulk_state)
    second = _invoke_physical_kk_gluon_adapter(adapter, physical_bulk_state)
    benchmark_default = module.default_paper_0710_1869_kk_gluon_couplings()

    assert isinstance(first, module.Paper07101869PhysicalKKGluonCouplings)
    assert isinstance(second, module.Paper07101869PhysicalKKGluonCouplings)
    assert json.dumps(first.as_dict(), sort_keys=True) == json.dumps(second.as_dict(), sort_keys=True)
    assert first.physical_profile_status_id == module.PAPER_0710_1869_PHYSICAL_KK_GLUON_STATUS_ID
    assert second.physical_profile_status_id == module.PAPER_0710_1869_PHYSICAL_KK_GLUON_STATUS_ID
    assert first.left_up_aligned.basis_id != benchmark_default.left_up_aligned.basis_id
    assert first.left_down_aligned.basis_id != benchmark_default.left_down_aligned.basis_id
    assert first.right_up.basis_id != benchmark_default.right_up.basis_id
    assert first.right_down.basis_id != benchmark_default.right_down.basis_id
    assert second.left_up_aligned.basis_id != benchmark_default.left_up_aligned.basis_id
    assert second.left_down_aligned.basis_id != benchmark_default.left_down_aligned.basis_id
    assert second.right_up.basis_id != benchmark_default.right_up.basis_id
    assert second.right_down.basis_id != benchmark_default.right_down.basis_id
    assert (
        json.dumps(first.summary().as_dict(), sort_keys=True)
        != json.dumps(benchmark_default.summary().as_dict(), sort_keys=True)
    )
