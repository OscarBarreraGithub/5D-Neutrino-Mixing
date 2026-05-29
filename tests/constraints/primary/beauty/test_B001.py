"""Production tests for B001 (Delta m_d in B_d mixing)."""

from __future__ import annotations

from dataclasses import replace
import math
from pathlib import Path

import numpy as np
import pytest
import yaml

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import point_builder
from flavor_catalog_constraints.anchors import AnchorError
from flavor_catalog_constraints.base import ConstraintProtocol, Severity
from flavor_catalog_constraints.primary.beauty import B001
from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.deltaf2 import (
    DEFAULT_DELTA_F2_INPUTS_V1,
    DELTA_M_BD_SM,
    compute_delta_f2_wilsons,
    evaluate_bd_mixing_with_running,
)

_PID = "B001"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "beauty" / "B001.yaml"


def _yaml():
    with open(_SIDECAR) as handle:
        return yaml.safe_load(handle)


def _yaml_pdg_block():
    return _yaml()["pdg_or_equivalent"]


def _yaml_sm_prediction():
    return _yaml_pdg_block()["standard_model_prediction"]


def _yaml_code_inputs():
    return _yaml()["auxiliary_code_inputs"]["deltaf2_bd_constants"]


def _bd_couplings(
    left: complex,
    right: complex,
    M_KK: float = 3000.0,
) -> QuarkMassBasisCouplings:
    """Minimal valid mass-basis couplings with only the d-b slot populated."""
    zeros = np.zeros((3, 3), dtype=np.complex128)
    left_down = zeros.copy()
    right_down = zeros.copy()
    left_down[0, 2] = left
    left_down[2, 0] = np.conj(left)
    right_down[0, 2] = right
    right_down[2, 0] = np.conj(right)
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


def _bd_wilsons(couplings: QuarkMassBasisCouplings):
    wilsons = compute_delta_f2_wilsons(
        couplings,
        inputs=DEFAULT_DELTA_F2_INPUTS_V1,
    )
    return next(w for w in wilsons if w.input.key == "b_d")


def _budget_from_yaml_and_core():
    exp = _yaml_pdg_block()["canonical_experimental_average"]
    sm = _yaml_sm_prediction()
    code_values = _yaml_code_inputs()["values"]
    exp_gev = float(code_values["DELTA_M_BD_EXP_GeV"])
    gev_per_ps_inverse = exp_gev / float(exp["value"])
    exp_sigma_gev = float(exp["uncertainty"]) * gev_per_ps_inverse
    sm_sigma_gev = float(sm["uncertainty"]) * gev_per_ps_inverse
    combined_sigma = math.sqrt(exp_sigma_gev * exp_sigma_gev + sm_sigma_gev * sm_sigma_gev)
    central_budget = abs(exp_gev - DELTA_M_BD_SM) / 2.0
    loose_budget = (abs(exp_gev - DELTA_M_BD_SM) + combined_sigma) / 2.0
    tight_budget = max(abs(exp_gev - DELTA_M_BD_SM) - combined_sigma, exp_sigma_gev) / 2.0
    return {
        "exp_gev": exp_gev,
        "gev_per_ps_inverse": gev_per_ps_inverse,
        "combined_sigma": combined_sigma,
        "exp_sigma_gev": exp_sigma_gev,
        "sm_sigma_gev": sm_sigma_gev,
        "sm_sigma_ps_inv": float(sm["uncertainty"]),
        "central_budget": central_budget,
        "tight_budget": tight_budget,
        "loose_budget": loose_budget,
        "sm_ps_inv": DELTA_M_BD_SM / gev_per_ps_inverse,
    }


def _expected_sm_sigma_from_yaml():
    sm = _yaml_sm_prediction()
    asymmetric_side = max(
        abs(float(sm["asymmetric_uncertainty_plus"])),
        abs(float(sm["asymmetric_uncertainty_minus"])),
    )
    return math.hypot(asymmetric_side, abs(float(sm["second_uncertainty"])))


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.HARD
    assert constraint.family == "beauty"
    assert constraint.observable == "Delta m_d"


def test_anchor_matches_yaml():
    constraint = fcc.get(_PID)
    exp = _yaml_pdg_block()["canonical_experimental_average"]
    sm = _yaml_sm_prediction()
    code = _yaml_code_inputs()
    budget = _budget_from_yaml_and_core()

    assert constraint.anchor.experimental.value == pytest.approx(exp["value"])
    assert constraint.anchor.experimental.uncertainty == pytest.approx(exp["uncertainty"])
    assert constraint.anchor.experimental.source_url == exp["source_url"]
    assert constraint.anchor.code_inputs.source == code["source"]
    assert constraint.anchor.code_inputs.lines == tuple(code["lines"])
    assert constraint.anchor.code_inputs.delta_m_bd_exp_gev == pytest.approx(
        code["values"]["DELTA_M_BD_EXP_GeV"]
    )
    assert constraint.anchor.code_inputs.f_bd_gev == pytest.approx(
        code["values"]["F_BD_GeV"]
    )
    assert constraint.anchor.standard_model.value == pytest.approx(sm["value"])
    assert constraint.anchor.standard_model.uncertainty == pytest.approx(
        sm["uncertainty"]
    )
    assert constraint.anchor.standard_model.source_url == sm["source_url"]
    assert constraint.anchor.sm_value == pytest.approx(DELTA_M_BD_SM)
    assert constraint.anchor.budget_band.gev_per_ps_inverse == pytest.approx(
        budget["gev_per_ps_inverse"]
    )
    assert constraint.anchor.budget_band.sm_delta_m_ps_inv == pytest.approx(
        budget["sm_ps_inv"]
    )
    assert constraint.anchor.budget_band.sm_delta_m_sigma_ps_inv == pytest.approx(
        budget["sm_sigma_ps_inv"]
    )
    assert constraint.anchor.budget_band.sm_delta_m_sigma_policy == (
        sm["uncertainty_policy"]
    )
    assert float(sm["uncertainty"]) == pytest.approx(
        _expected_sm_sigma_from_yaml(),
        rel=5.0e-3,
        abs=5.0e-4,
    )
    assert constraint.anchor.central_budget == pytest.approx(budget["central_budget"])
    assert constraint.anchor.budget_band.tight_budget == pytest.approx(
        budget["tight_budget"]
    )
    assert constraint.anchor.budget_band.combined_delta_m_sigma_gev == pytest.approx(
        budget["combined_sigma"]
    )
    assert constraint.anchor.budget == pytest.approx(budget["loose_budget"])
    assert constraint.anchor.budget == pytest.approx(3.369899779459366e-14)


def test_missing_or_mismatched_sm_uncertainty_anchor_fails_loudly():
    constraint = fcc.get(_PID)
    sm_sub = _yaml_sm_prediction()
    kwargs = {
        "process_id": _PID,
        "experimental": constraint.anchor.experimental,
        "standard_model_sub": sm_sub,
        "code_inputs": constraint.anchor.code_inputs,
    }

    missing = replace(constraint.anchor.standard_model, uncertainty=None)
    with pytest.raises(AnchorError, match="SM Delta m_d uncertainty anchor"):
        B001._build_budget_band(standard_model=missing, **kwargs)

    mismatched = replace(
        constraint.anchor.standard_model,
        uncertainty=constraint.anchor.standard_model.uncertainty * 2.0,
    )
    with pytest.raises(AnchorError, match="does not match component quadrature"):
        B001._build_budget_band(standard_model=mismatched, **kwargs)


def test_evaluate_without_input_degrades_gracefully():
    result = fcc.get(_PID).evaluate(point_builder.empty_point())

    assert result.process_id == _PID
    assert result.passes is True
    assert result.predicted is None
    assert result.ratio is None
    assert result.experimental == pytest.approx(fcc.get(_PID).anchor.delta_m_experimental_gev)
    assert result.sm_prediction == pytest.approx(fcc.get(_PID).anchor.sm_value)
    assert result.budget == pytest.approx(fcc.get(_PID).anchor.budget)
    assert result.diagnostics["missing_extra"] == "quark_mass_basis_couplings"


def test_evaluate_runs_end_to_end_with_real_couplings_and_real_finite_fields():
    couplings = _bd_couplings(left=1.0e-3 + 0.4e-3j, right=0.7e-3 + 0.2e-3j)
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
    assert isinstance(result.diagnostics["left_db_coupling"], complex)
    assert isinstance(result.diagnostics["wilson_coefficients"]["C4_LR"], complex)
    for key in (
        "abs_m12_np_gev",
        "hadronic_scale_gev",
        "matching_scale_gev",
        "m_kk_gev",
        "central_np_budget",
        "tight_band_np_budget",
        "loose_band_np_budget",
        "hard_veto_np_budget",
        "budget_combined_delta_m_sigma_gev",
        "budget_experimental_delta_m_sigma_gev",
        "budget_sm_delta_m_sigma_gev",
        "experimental_delta_m_ps_inv",
        "sm_delta_m_ps_inv",
        "gev_per_ps_inverse",
        "core_default_budget",
        "core_default_ratio",
    ):
        assert isinstance(result.diagnostics[key], float)
        assert math.isfinite(result.diagnostics[key])
    assert result.diagnostics["qcd_running_applied"] is True
    assert result.diagnostics["hadronic_scale_gev"] == pytest.approx(2.0)
    assert result.diagnostics["core_input_key"] == "b_d"
    assert result.diagnostics["down_sector_indices"] == (0, 2)


@pytest.mark.parametrize(
    ("couplings", "expected_pass"),
    [
        (_bd_couplings(left=1.0e-3 + 0.4e-3j, right=0.7e-3 + 0.2e-3j), True),
        (_bd_couplings(left=5.0e-3 + 2.0e-3j, right=3.5e-3 + 1.0e-3j), False),
    ],
)
def test_pass_fail_and_numbers_match_running_evaluator_amplitude(
    couplings: QuarkMassBasisCouplings,
    expected_pass: bool,
):
    constraint = fcc.get(_PID)
    point = point_builder.build_from_quark_couplings(couplings)
    result = constraint.evaluate(point)
    audited = evaluate_bd_mixing_with_running(_bd_wilsons(couplings), mu_had=2.0)

    assert result.passes is expected_pass
    assert result.predicted == pytest.approx(audited.abs_m12_np)
    assert result.budget == pytest.approx(constraint.anchor.budget)
    assert result.ratio == pytest.approx(audited.abs_m12_np / constraint.anchor.budget)
    assert result.diagnostics["core_default_budget"] == pytest.approx(audited.budget)
    assert result.diagnostics["core_default_ratio"] == pytest.approx(
        audited.ratio_to_budget
    )
    if expected_pass:
        assert result.predicted == pytest.approx(2.1778918692437868e-14)
        assert result.ratio == pytest.approx(0.6462779345898491)
        assert result.ratio <= 1.0
    else:
        assert result.ratio > 10.0


def test_evaluate_is_pure_and_deterministic():
    couplings = _bd_couplings(left=1.0e-3 + 0.4e-3j, right=0.7e-3 + 0.2e-3j)
    before_left_down = couplings.left_down.copy()
    before_right_down = couplings.right_down.copy()
    point = point_builder.build_from_quark_couplings(couplings)
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    np.testing.assert_array_equal(couplings.left_down, before_left_down)
    np.testing.assert_array_equal(couplings.right_down, before_right_down)
