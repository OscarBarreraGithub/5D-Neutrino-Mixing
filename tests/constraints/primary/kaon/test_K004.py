"""Production tests for K004 (K+ -> pi+ nu nubar)."""

from __future__ import annotations

import math
from pathlib import Path
import re

import numpy as np
import pytest
import yaml

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import point_builder
from flavor_catalog_constraints.base import ConstraintProtocol, Severity
from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.rare_kaon_snd import evaluate_kplus_to_piplus_nunu

_PID = "K004"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "kaon" / "K004.yaml"
_ASYM_RE = re.compile(r"^\s*\+?(?P<upper>[0-9.eE+-]+)\s*/\s*-(?P<lower>[0-9.eE+-]+)\s*$")


def _yaml_pdg_block():
    with open(_SIDECAR) as handle:
        return yaml.safe_load(handle)["pdg_or_equivalent"]


def _scaled(block):
    return float(block["value"]) * float(block.get("scale", 1.0))


def _scaled_uncertainty_pair(block):
    scale = float(block.get("scale", 1.0))
    uncertainty = block["uncertainty"]
    if isinstance(uncertainty, str):
        match = _ASYM_RE.match(uncertainty)
        if match is not None:
            return (
                float(match.group("upper")) * scale,
                float(match.group("lower")) * scale,
            )
    sigma = float(uncertainty) * scale
    return sigma, sigma


def _sd_couplings(
    left: complex,
    right: complex = 0.0j,
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


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.HARD
    assert constraint.family == "kaon"
    assert constraint.observable == "BR(K+ -> pi+ nu nubar)"


def test_anchor_matches_yaml_and_budget_band():
    constraint = fcc.get(_PID)
    pdg = _yaml_pdg_block()
    exp = pdg["latest_experimental_value"]
    sm = pdg["sm_prediction_buras_venturini_2022"]
    bgs = pdg["sm_prediction_brod_gorbahn_stamou_2021"]
    exp_up, exp_down = _scaled_uncertainty_pair(exp)
    sm_up, sm_down = _scaled_uncertainty_pair(sm)
    combined_up = math.sqrt(exp_up * exp_up + sm_up * sm_up)
    combined_down = math.sqrt(exp_down * exp_down + sm_down * sm_down)

    assert constraint.anchor.experimental.value == pytest.approx(_scaled(exp))
    assert constraint.anchor.experimental.uncertainty_upper == pytest.approx(exp_up)
    assert constraint.anchor.experimental.uncertainty_lower == pytest.approx(exp_down)
    assert constraint.anchor.experimental.source_url == exp["source_url"]
    assert constraint.anchor.standard_model.value == pytest.approx(_scaled(sm))
    assert constraint.anchor.standard_model.uncertainty == pytest.approx(sm_up)
    assert constraint.anchor.standard_model.source_url == sm["source_url"]
    assert constraint.anchor.validation_standard_model.value == pytest.approx(
        _scaled(bgs)
    )
    assert constraint.anchor.validation_standard_model.source_url == bgs["source_url"]
    assert constraint.anchor.budget_band.central_residual == pytest.approx(
        _scaled(exp) - _scaled(sm)
    )
    assert constraint.anchor.budget_band.combined_sigma_upper == pytest.approx(
        combined_up
    )
    assert constraint.anchor.budget_band.combined_sigma_lower == pytest.approx(
        combined_down
    )
    assert constraint.anchor.budget == pytest.approx(max(combined_up, combined_down))


def test_evaluate_without_input_degrades_gracefully():
    result = fcc.get(_PID).evaluate(point_builder.empty_point())

    assert result.process_id == _PID
    assert result.passes is True
    assert result.predicted is None
    assert result.ratio is None
    assert result.experimental == pytest.approx(fcc.get(_PID).anchor.value)
    assert result.sm_prediction == pytest.approx(
        fcc.get(_PID).sm_result.branching_fraction
    )
    assert result.budget == pytest.approx(fcc.get(_PID).anchor.budget)
    assert result.diagnostics["missing_extra"] == "quark_mass_basis_couplings"


def test_sm_limit_branching_fraction_matches_short_distance_reference():
    constraint = fcc.get(_PID)
    couplings = _sd_couplings(left=0.0j, right=0.0j)
    point = point_builder.build_from_quark_couplings(couplings)
    result = constraint.evaluate(point)
    direct = evaluate_kplus_to_piplus_nunu(couplings)
    sm_anchor = constraint.anchor.sm_value

    assert result.predicted == pytest.approx(8.472598449611133e-11)
    assert result.predicted == pytest.approx(result.sm_prediction)
    assert result.predicted == pytest.approx(direct.sm_branching_fraction)
    assert result.diagnostics["sm_anchor_branching_fraction"] == pytest.approx(
        sm_anchor
    )
    assert result.diagnostics["sm_formula_minus_anchor"] == pytest.approx(
        result.sm_prediction - sm_anchor
    )
    assert result.diagnostics["delta_em_correction"] == pytest.approx(0.997)
    assert result.passes is True


def test_evaluate_runs_end_to_end_with_real_finite_fields_and_complex_diagnostics():
    couplings = _sd_couplings(left=1.0e-5 + 0.2e-5j, right=0.5e-5j)
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
        "left_sd_coupling",
        "right_sd_coupling",
        "lambda_c",
        "lambda_t",
        "x_eff_top",
        "x_np_total",
    ):
        assert isinstance(result.diagnostics[key], complex)
    for key in (
        "m_kk_gev",
        "matching_scale_gev",
        "kappa_plus",
        "p_c",
        "x_t",
        "lambda_wolfenstein",
        "np_shift_branching_fraction",
        "budget_combined_sigma_upper",
        "budget_combined_sigma_lower",
    ):
        assert isinstance(result.diagnostics[key], float)
        assert math.isfinite(result.diagnostics[key])
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]


@pytest.mark.parametrize(
    ("couplings", "expected_pass"),
    [
        (_sd_couplings(left=1.0e-5), True),
        (_sd_couplings(left=5.0e-2), False),
    ],
)
def test_safe_point_passes_and_large_np_point_fails(
    couplings: QuarkMassBasisCouplings,
    expected_pass: bool,
):
    constraint = fcc.get(_PID)
    point = point_builder.build_from_quark_couplings(couplings)
    result = constraint.evaluate(point)
    direct = evaluate_kplus_to_piplus_nunu(couplings)

    assert result.passes is expected_pass
    assert result.predicted == pytest.approx(direct.branching_fraction)
    assert result.sm_prediction == pytest.approx(direct.sm_branching_fraction)
    if expected_pass:
        assert result.ratio < 1.0
    else:
        assert result.ratio > 10.0


def test_optional_kk_ew_mass_extra_changes_matching_scale():
    couplings = _sd_couplings(left=1.0e-3)
    default_point = point_builder.build_from_quark_couplings(couplings)
    ew_point = point_builder.make_point(
        quark_mass_basis_couplings=couplings,
        kk_ew_mass_gev=6000.0,
    )
    default_result = fcc.get(_PID).evaluate(default_point)
    ew_result = fcc.get(_PID).evaluate(ew_point)

    assert default_result.diagnostics["m_kk_gev"] == pytest.approx(3000.0)
    assert ew_result.diagnostics["m_kk_gev"] == pytest.approx(6000.0)
    assert ew_result.diagnostics["kk_ew_mass_extra_used"] is True
    assert abs(ew_result.diagnostics["x_np_total"]) == pytest.approx(
        abs(default_result.diagnostics["x_np_total"]) / 4.0
    )


def test_evaluate_is_pure_and_deterministic():
    couplings = _sd_couplings(left=1.0e-5 + 0.2e-5j, right=0.5e-5j)
    before_left_down = couplings.left_down.copy()
    before_right_down = couplings.right_down.copy()
    point = point_builder.build_from_quark_couplings(couplings)
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    np.testing.assert_array_equal(couplings.left_down, before_left_down)
    np.testing.assert_array_equal(couplings.right_down, before_right_down)
