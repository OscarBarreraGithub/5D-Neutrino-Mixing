"""Production tests for K002 (Delta m_K)."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest
from scipy import constants
import yaml

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import point_builder
from flavor_catalog_constraints.base import ConstraintProtocol, Severity
from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.deltaf2 import (
    DEFAULT_DELTA_F2_INPUTS_V1,
    DELTA_M_K,
    compute_delta_f2_wilsons,
    evaluate_delta_mk,
    evaluate_delta_mk_with_running,
)

_PID = "K002"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "kaon" / "K002.yaml"
_HBAR_GEV_SECONDS = float(constants.hbar / constants.electron_volt * 1.0e-9)
_UNIT_SCALE_HBAR_PER_SECOND = 1.0e10


def _yaml_pdg_block():
    with open(_SIDECAR) as handle:
        return yaml.safe_load(handle)["pdg_or_equivalent"]


def _sd_couplings(
    left: complex,
    right: complex,
    M_KK: float = 3000.0,
) -> QuarkMassBasisCouplings:
    """Minimal valid mass-basis couplings with only the s-d slot populated."""
    zeros = np.zeros((3, 3), dtype=np.complex128)
    left_down = zeros.copy()
    right_down = zeros.copy()
    left_down[0, 1] = left
    left_down[1, 0] = np.conj(left)
    right_down[0, 1] = right
    right_down[1, 0] = np.conj(right)
    return QuarkMassBasisCouplings(
        M_KK=M_KK,
        xi_KK=1.0,
        alpha_s=0.09,
        g_s=1.0,
        left_overlap=zeros,
        right_up_overlap=zeros,
        right_down_overlap=zeros,
        left_up=zeros,
        left_down=left_down,
        right_up=zeros,
        right_down=right_down,
    )


def _kaon_wilsons(couplings: QuarkMassBasisCouplings):
    wilsons = compute_delta_f2_wilsons(
        couplings,
        inputs=DEFAULT_DELTA_F2_INPUTS_V1,
    )
    return next(w for w in wilsons if w.input.key == "epsilon_k")


def _audited_delta_mk_with_running(couplings: QuarkMassBasisCouplings):
    return evaluate_delta_mk_with_running(_kaon_wilsons(couplings), mu_had=2.0)


def _gev_from_yaml_rate(value: float) -> float:
    return float(value) * _UNIT_SCALE_HBAR_PER_SECOND * _HBAR_GEV_SECONDS


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.HARD
    assert constraint.family == "kaon"
    assert constraint.observable == "Delta m_K"


def test_anchor_matches_yaml_and_core_budget():
    constraint = fcc.get(_PID)
    pdg = _yaml_pdg_block()
    exp = pdg["pdg_fit_assuming_cpt"]
    expected_value_gev = _gev_from_yaml_rate(exp["value"])
    expected_uncertainty_gev = _gev_from_yaml_rate(exp["uncertainty"])

    assert constraint.anchor.experimental.value == pytest.approx(exp["value"])
    assert constraint.anchor.experimental.uncertainty == pytest.approx(exp["uncertainty"])
    assert constraint.anchor.experimental.units == exp["units"]
    assert constraint.anchor.experimental.source_url == exp["source_url"]
    assert constraint.anchor.value == pytest.approx(expected_value_gev)
    assert constraint.anchor.uncertainty == pytest.approx(expected_uncertainty_gev)
    assert constraint.anchor.budget == pytest.approx(expected_value_gev / 2.0)

    core_budget = DELTA_M_K / 2.0
    assert core_budget == pytest.approx(1.742e-15)
    assert abs(DELTA_M_K - constraint.anchor.value) < constraint.anchor.uncertainty


def test_evaluate_without_input_degrades_gracefully():
    result = fcc.get(_PID).evaluate(point_builder.empty_point())

    assert result.process_id == _PID
    assert result.passes is True
    assert result.predicted is None
    assert result.ratio is None
    assert result.experimental == pytest.approx(fcc.get(_PID).anchor.value)
    assert result.sm_prediction is None
    assert result.budget == pytest.approx(fcc.get(_PID).anchor.budget)
    assert result.diagnostics["missing_extra"] == "quark_mass_basis_couplings"


def test_evaluate_runs_end_to_end_with_real_couplings_and_real_finite_fields():
    couplings = _sd_couplings(left=1.0e-4 + 0.5e-4j, right=1.0e-4 + 0.2e-4j)
    point = point_builder.build_from_quark_couplings(couplings)
    result = fcc.get(_PID).evaluate(point)

    assert result.process_id == _PID
    for value in (
        result.predicted,
        result.ratio,
        result.budget,
        result.experimental,
    ):
        assert isinstance(value, float)
        assert math.isfinite(value)
    assert result.sm_prediction is None
    assert isinstance(result.diagnostics["left_sd_coupling"], complex)
    assert isinstance(result.diagnostics["wilson_coefficients"]["C4_LR"], complex)
    for key in (
        "abs_m12_np_gev",
        "hadronic_scale_gev",
        "matching_scale_gev",
        "m_kk_gev",
        "delta_m_k_exp_gev",
        "delta_m_k_uncertainty_gev",
        "core_delta_m_k_exp_gev",
        "core_m12_budget_gev",
        "core_ratio_to_exp_without_catalog_override",
        "core_catalog_budget_pull_sigma",
        "catalog_core_budget_relative_difference",
    ):
        assert isinstance(result.diagnostics[key], float)
        assert math.isfinite(result.diagnostics[key])
    assert result.diagnostics["qcd_running_applied"] is True
    assert result.diagnostics["hadronic_scale_gev"] == pytest.approx(2.0)
    assert result.diagnostics["sm_subtracted"] is False
    assert result.diagnostics["long_distance_dominated"] is True


def test_reference_couplings_show_qcd_running_enhancement():
    couplings = _sd_couplings(left=1.0e-4 + 0.5e-4j, right=1.0e-4 + 0.2e-4j)
    wilsons = _kaon_wilsons(couplings)

    unrun = evaluate_delta_mk(wilsons)
    run = evaluate_delta_mk_with_running(wilsons, mu_had=2.0)

    assert unrun.ratio_to_exp == pytest.approx(0.008753373215325102)
    assert run.ratio_to_exp == pytest.approx(0.09343760390731619)
    assert run.ratio_to_exp > 10.0 * unrun.ratio_to_exp


@pytest.mark.parametrize(
    ("couplings", "expected_pass"),
    [
        (_sd_couplings(left=1.0e-4 + 0.5e-4j, right=1.0e-4 + 0.2e-4j), True),
        (_sd_couplings(left=1.0e-3 + 0.5e-3j, right=1.0e-3 + 0.2e-3j), False),
    ],
)
def test_pass_fail_and_numbers_match_audited_evaluate_delta_mk(
    couplings: QuarkMassBasisCouplings,
    expected_pass: bool,
):
    constraint = fcc.get(_PID)
    point = point_builder.build_from_quark_couplings(couplings)
    result = constraint.evaluate(point)
    audited = _audited_delta_mk_with_running(couplings)
    audited_budget = audited.abs_m12_np / audited.ratio_to_exp

    assert result.passes is expected_pass
    assert result.predicted == pytest.approx(audited.abs_m12_np)
    assert result.budget == pytest.approx(constraint.anchor.budget)
    assert result.budget == pytest.approx(audited_budget, rel=3.0e-5)
    assert result.ratio == pytest.approx(result.predicted / result.budget)
    assert result.ratio == pytest.approx(audited.ratio_to_exp, rel=3.0e-5)
    assert result.diagnostics["core_ratio_to_exp_without_catalog_override"] == (
        pytest.approx(audited.ratio_to_exp)
    )
    if expected_pass:
        assert result.ratio < 0.1
    else:
        assert result.ratio > 9.0


def test_evaluate_is_pure_and_deterministic():
    couplings = _sd_couplings(left=1.0e-4 + 0.5e-4j, right=1.0e-4 + 0.2e-4j)
    before_left_down = couplings.left_down.copy()
    point = point_builder.build_from_quark_couplings(couplings)
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    np.testing.assert_array_equal(couplings.left_down, before_left_down)
