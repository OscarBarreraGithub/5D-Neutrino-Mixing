"""Contract-first tests for the paper-mode coupling slice.

These tests are written to be forward-compatible while the paper-mode coupling
API is still landing. They are strict once the coupling module exists.
"""

from __future__ import annotations

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
