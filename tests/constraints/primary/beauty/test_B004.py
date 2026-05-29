"""Production tests for B004 (phi_s in B_s -> J/psi phi)."""

from __future__ import annotations

import cmath
import math
from pathlib import Path

import numpy as np
import pytest
import yaml

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import point_builder
from flavor_catalog_constraints.anchors import Anchor, AnchorError
from flavor_catalog_constraints.base import ConstraintProtocol, Severity
from flavor_catalog_constraints.primary.beauty import B004 as b004_module
from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.deltaf2 import (
    B_1_BS,
    B_4_BS,
    B_5_BS,
    DEFAULT_DELTA_F2_INPUTS_V1,
    DELTA_M_BS_SM,
    F_BS,
    M_BS,
    M_B_QUARK,
    M_S_QUARK_BS,
    compute_delta_f2_wilsons,
    compute_m12_np,
    evaluate_bs_mixing_with_running,
    _evolve_wilsons,
)

_PID = "B004"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "beauty" / "B004.yaml"


def _yaml():
    with open(_SIDECAR) as handle:
        return yaml.safe_load(handle)


def _yaml_pdg_block():
    return _yaml()["pdg_or_equivalent"]


def _wrap_phase(value: float) -> float:
    return float(math.atan2(math.sin(value), math.cos(value)))


def _bs_couplings(
    left: complex,
    right: complex,
    M_KK: float = 3000.0,
) -> QuarkMassBasisCouplings:
    """Minimal valid mass-basis couplings with only the s-b slot populated."""
    zeros = np.zeros((3, 3), dtype=np.complex128)
    left_down = zeros.copy()
    right_down = zeros.copy()
    left_down[1, 2] = left
    left_down[2, 1] = np.conj(left)
    right_down[1, 2] = right
    right_down[2, 1] = np.conj(right)
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


def _bs_wilsons(couplings: QuarkMassBasisCouplings):
    wilsons = compute_delta_f2_wilsons(
        couplings,
        inputs=DEFAULT_DELTA_F2_INPUTS_V1,
    )
    return next(w for w in wilsons if w.input.key == "b_s")


def _budget_from_yaml_and_core():
    pdg = _yaml_pdg_block()
    exp = pdg["canonical_hflav_average"]
    sm = pdg["standard_model_reference"]
    exp_sigma = float(exp["uncertainty"])
    sm_upper = abs(float(sm["upper_uncertainty"]))
    sm_lower = abs(float(sm["lower_uncertainty"]))
    combined_upper = math.sqrt(exp_sigma * exp_sigma + sm_upper * sm_upper)
    combined_lower = math.sqrt(exp_sigma * exp_sigma + sm_lower * sm_lower)
    return {
        "central_residual": abs(_wrap_phase(float(sm["value"]) - float(exp["value"]))),
        "combined_upper": combined_upper,
        "combined_lower": combined_lower,
        "hard_budget": max(combined_upper, combined_lower),
        "m12_sm": DELTA_M_BS_SM / 2.0,
    }


def _direct_core_phi_s(
    couplings: QuarkMassBasisCouplings,
    constraint,
):
    wilsons = _bs_wilsons(couplings)
    evolved = _evolve_wilsons(wilsons, mu_had=2.0)
    m12_np = compute_m12_np(
        evolved,
        F_BS,
        M_BS,
        M_B_QUARK,
        M_S_QUARK_BS,
        B_1_BS,
        B_4_BS,
        B_5_BS,
    )
    magnitude = evaluate_bs_mixing_with_running(wilsons, mu_had=2.0)
    m12_ratio = m12_np / constraint.anchor.budget_band.m12_sm_gev
    phi_s_np = cmath.phase(1.0 + m12_ratio)
    predicted = _wrap_phase(constraint.anchor.sm_value + phi_s_np)
    signed_residual = _wrap_phase(predicted - constraint.anchor.value)
    budget = (
        constraint.anchor.budget_band.combined_sigma_upper
        if signed_residual >= 0.0
        else constraint.anchor.budget_band.combined_sigma_lower
    )
    residual = abs(signed_residual)
    ratio = residual / budget
    return m12_np, magnitude, phi_s_np, predicted, residual, budget, ratio


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.HARD
    assert constraint.family == "beauty"
    assert constraint.observable == "phi_s(B_s -> J/psi phi)"


def test_anchor_matches_yaml_and_budget_band():
    constraint = fcc.get(_PID)
    pdg = _yaml_pdg_block()
    exp = pdg["canonical_hflav_average"]
    mode = pdg["mode_specific_jpsi_kk_average"]
    latest = pdg["latest_lhcb_jpsi_kk_input"]
    sm = pdg["standard_model_reference"]
    budget = _budget_from_yaml_and_core()

    assert constraint.anchor.experimental.value == pytest.approx(exp["value"])
    assert constraint.anchor.experimental.uncertainty == pytest.approx(
        exp["uncertainty"]
    )
    assert constraint.anchor.experimental.units == exp["units"]
    assert constraint.anchor.experimental.source_url == exp["source_url"]
    assert constraint.anchor.mode_specific_jpsi_kk.value == pytest.approx(
        mode["value"]
    )
    assert constraint.anchor.mode_specific_jpsi_kk.uncertainty == pytest.approx(
        mode["uncertainty"]
    )
    assert constraint.anchor.latest_lhcb_jpsi_kk.value == pytest.approx(
        latest["value"]
    )
    assert constraint.anchor.standard_model.value == pytest.approx(sm["value"])
    assert constraint.anchor.standard_model.upper_uncertainty == pytest.approx(
        sm["upper_uncertainty"]
    )
    assert constraint.anchor.standard_model.lower_uncertainty == pytest.approx(
        sm["lower_uncertainty"]
    )
    assert constraint.anchor.standard_model.source_url == sm["source_url"]

    assert constraint.anchor.budget_band.central_residual == pytest.approx(
        budget["central_residual"]
    )
    assert constraint.anchor.budget_band.combined_sigma_upper == pytest.approx(
        budget["combined_upper"]
    )
    assert constraint.anchor.budget_band.combined_sigma_lower == pytest.approx(
        budget["combined_lower"]
    )
    assert constraint.anchor.budget == pytest.approx(budget["hard_budget"])
    assert constraint.anchor.budget == pytest.approx(0.016025292508658584)
    assert constraint.anchor.budget_band.m12_sm_gev == pytest.approx(
        budget["m12_sm"]
    )


def test_phase_anchor_metadata_is_loudly_validated(monkeypatch: pytest.MonkeyPatch):
    bad_anchor = Anchor(
        process_id=_PID,
        block_key="bad_phi_s",
        value=-0.041,
        uncertainty=0.016,
        observable="phi_s^ccs",
        units="degrees",
    )
    monkeypatch.setattr(b004_module, "load_anchor", lambda *args, **kwargs: bad_anchor)

    with pytest.raises(AnchorError, match="must use units"):
        b004_module._load_phi_s_anchor(_PID)


def test_evaluate_without_input_degrades_gracefully():
    constraint = fcc.get(_PID)
    result = constraint.evaluate(point_builder.empty_point())

    assert result.process_id == _PID
    assert result.passes is True
    assert result.predicted is None
    assert result.ratio is None
    assert result.experimental == pytest.approx(constraint.anchor.value)
    assert result.sm_prediction == pytest.approx(constraint.anchor.sm_value)
    assert result.budget == pytest.approx(constraint.anchor.budget)
    assert result.diagnostics["missing_extra"] == "quark_mass_basis_couplings"
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"][0]


def test_sm_limit_uses_yaml_phi_s_reference():
    constraint = fcc.get(_PID)
    couplings = _bs_couplings(left=0.0j, right=0.0j)
    result = constraint.evaluate(point_builder.build_from_quark_couplings(couplings))

    assert result.predicted == pytest.approx(-0.0368)
    assert result.predicted == pytest.approx(result.sm_prediction)
    assert result.ratio == pytest.approx(
        abs(result.predicted - constraint.anchor.value)
        / constraint.anchor.budget_band.combined_sigma_upper
    )
    assert result.passes is True
    assert result.diagnostics["phi_s_np_rad"] == pytest.approx(0.0)
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"][0]


def test_evaluate_runs_end_to_end_with_real_finite_fields_and_complex_diagnostics():
    couplings = _bs_couplings(left=1.0e-3 + 1.0e-3j, right=0.0j)
    point = point_builder.build_from_quark_couplings(couplings)
    result = fcc.get(_PID).evaluate(point)

    assert result.process_id == _PID
    for value in (
        result.predicted,
        result.ratio,
        result.budget,
        result.sm_prediction,
        result.experimental,
    ):
        assert isinstance(value, float)
        assert math.isfinite(value)
    for key in (
        "m12_np_gev",
        "m12_np_over_m12_sm",
        "total_mixing_ratio",
        "left_sb_coupling",
        "right_sb_coupling",
    ):
        assert isinstance(result.diagnostics[key], complex)
    assert isinstance(result.diagnostics["wilson_coefficients"]["C1_VLL"], complex)
    for key in (
        "abs_m12_np_gev",
        "im_m12_np_gev",
        "core_abs_m12_np_gev",
        "m12_sm_gev",
        "re_m12_np_over_m12_sm",
        "im_m12_np_over_m12_sm",
        "phi_s_np_rad",
        "phi_s_np_deg",
        "phi_s_total_rad",
        "phi_s_residual_rad",
        "hadronic_scale_gev",
        "matching_scale_gev",
        "m_kk_gev",
        "central_sm_exp_residual",
        "hard_veto_budget",
        "budget_combined_sigma_upper",
        "budget_combined_sigma_lower",
        "budget_experimental_sigma",
        "budget_sm_phase_sigma_upper",
        "budget_sm_phase_sigma_lower",
        "delta_m_bs_sm_gev",
        "core_delta_m_bs_exp_gev",
        "core_legacy_m12_budget",
        "core_legacy_ratio_to_budget",
    ):
        assert isinstance(result.diagnostics[key], float)
        assert math.isfinite(result.diagnostics[key])
    assert result.diagnostics["qcd_running_applied"] is True
    assert result.diagnostics["hadronic_scale_gev"] == pytest.approx(2.0)
    assert result.diagnostics["core_input_key"] == "b_s"
    assert result.diagnostics["down_sector_indices"] == (1, 2)
    assert result.diagnostics["phase_uses_complex_m12_not_abs"] is True
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"][0]


def test_numbers_match_direct_running_complex_m12_phase_evaluator():
    couplings = _bs_couplings(left=1.0e-3 + 1.0e-3j, right=0.0j)
    constraint = fcc.get(_PID)
    point = point_builder.build_from_quark_couplings(couplings)
    result = constraint.evaluate(point)
    m12_np, magnitude, phi_s_np, predicted, residual, budget, ratio = (
        _direct_core_phi_s(couplings, constraint)
    )

    assert result.predicted == pytest.approx(predicted)
    assert result.ratio == pytest.approx(ratio)
    assert result.budget == pytest.approx(budget)
    assert result.diagnostics["phi_s_residual_rad"] == pytest.approx(residual)
    assert result.diagnostics["phi_s_np_rad"] == pytest.approx(phi_s_np)
    assert result.diagnostics["m12_np_gev"].real == pytest.approx(m12_np.real)
    assert result.diagnostics["m12_np_gev"].imag == pytest.approx(m12_np.imag)
    assert result.diagnostics["abs_m12_np_gev"] == pytest.approx(abs(m12_np))
    assert result.diagnostics["core_abs_m12_np_gev"] == pytest.approx(
        magnitude.abs_m12_np
    )
    assert abs(m12_np) == pytest.approx(magnitude.abs_m12_np)


@pytest.mark.parametrize(
    ("couplings", "expected_pass"),
    [
        (_bs_couplings(left=1.0e-3 + 1.0e-3j, right=0.0j), True),
        (_bs_couplings(left=1.0e-2 + 1.0e-2j, right=0.0j), False),
    ],
)
def test_safe_point_passes_and_large_np_point_fails(
    couplings: QuarkMassBasisCouplings,
    expected_pass: bool,
):
    constraint = fcc.get(_PID)
    result = constraint.evaluate(point_builder.build_from_quark_couplings(couplings))

    assert result.passes is expected_pass
    if expected_pass:
        assert result.ratio <= 1.0
    else:
        assert result.ratio > 1.0


def test_evaluate_is_pure_and_deterministic():
    couplings = _bs_couplings(left=1.0e-3 + 1.0e-3j, right=0.0j)
    before_left_down = couplings.left_down.copy()
    before_right_down = couplings.right_down.copy()
    point = point_builder.build_from_quark_couplings(couplings)
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    np.testing.assert_array_equal(couplings.left_down, before_left_down)
    np.testing.assert_array_equal(couplings.right_down, before_right_down)
