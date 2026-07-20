"""Production tests for B005 (B_s -> mu+ mu-)."""

from __future__ import annotations

import math
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest
import yaml

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import point_builder
from flavor_catalog_constraints.anchors import AnchorError
from flavor_catalog_constraints.base import ConstraintProtocol, Severity
from flavor_catalog_constraints.primary.beauty import B005 as b005_module
from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.model import RotationParameters, ckm_like_unitary
from quarkConstraints.rare_b_dilepton import (
    RareBDileptonWilsonCoefficients,
    ckm_factors,
    evaluate_bq_to_mumu,
)
from tests.rare_b_phase3d_helpers import (
    core_wilsons_from_rs_coeff,
    rs_coeff,
    sample_rare_b_point,
    scaled_rare_b_point,
    sm_limit_rare_b_point,
)

_PID = "B005"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "beauty" / "B005.yaml"


def _yaml_pdg_block():
    with open(_SIDECAR) as handle:
        return yaml.safe_load(handle)["pdg_or_equivalent"]


def _bs_couplings(
    left: complex,
    right: complex = 0.0j,
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


def _manual_sm_bs_mumu(inputs) -> float:
    matrix = ckm_like_unitary(
        RotationParameters(
            theta12=inputs.theta12,
            theta13=inputs.theta13,
            theta23=inputs.theta23,
            delta=inputs.delta,
        )
    )
    meson = inputs.bs
    lambda_t = complex(matrix[2, 2] * np.conjugate(matrix[2, 1]))
    tau_gev_inv = meson.lifetime_ps * 1.0e-12 / inputs.hbar_gev_s
    beta = math.sqrt(
        1.0
        - 4.0 * inputs.muon_mass_gev**2 / meson.meson_mass_gev**2
    )
    time_factor = (1.0 + meson.width_difference_y) / (
        1.0 - meson.width_difference_y**2
    )
    return float(
        tau_gev_inv
        * inputs.gf_gev_minus2**2
        * inputs.alpha_em_mz**2
        / (16.0 * math.pi**3)
        * meson.decay_constant_gev**2
        * meson.meson_mass_gev
        * inputs.muon_mass_gev**2
        * abs(lambda_t) ** 2
        * beta
        * abs(inputs.c10_sm) ** 2
        * time_factor
    )


def _synthetic_wilsons(inputs, **wilson_overrides) -> RareBDileptonWilsonCoefficients:
    factors = ckm_factors("b_s", inputs)
    values = {
        "c9_np": 0.0j,
        "c10_np": 0.0j,
        "c9p_np": 0.0j,
        "c10p_np": 0.0j,
        "cs_np": 0.0j,
        "csp_np": 0.0j,
        "cp_np": 0.0j,
        "cpp_np": 0.0j,
    }
    values.update(wilson_overrides)
    return RareBDileptonWilsonCoefficients(
        model_label="test_rare_b_dilepton_wilsons",
        operator_convention="unit-test",
        matching_assumption="unit-test direct Wilson coefficients",
        transition_key="b_s",
        M_KK=3000.0,
        matching_scale=3000.0,
        lambda_t=factors.lambda_t,
        left_qb_coupling=0.0j,
        right_qb_coupling=0.0j,
        left_qb_overlap=0.0j,
        right_qb_overlap=0.0j,
        left_quark_delta=0.0j,
        right_quark_delta=0.0j,
        muon_left_delta=0.0,
        muon_right_delta=0.0,
        muon_vector_delta=0.0,
        muon_axial_delta=0.0,
        **values,
    )


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.HARD
    assert constraint.family == "beauty"
    assert constraint.observable == "BR(B_s -> mu+ mu-)"


def test_anchor_matches_yaml_and_budget_band():
    constraint = fcc.get(_PID)
    pdg = _yaml_pdg_block()
    exp = pdg["canonical_experimental_average"]
    hist = pdg["hflav_historical_average"]
    sm = pdg["standard_model_prediction"]
    combined = math.sqrt(exp["uncertainty"] ** 2 + sm["uncertainty"] ** 2)
    central = abs(exp["value"] - sm["value"])

    assert constraint.anchor.experimental.value == pytest.approx(exp["value"])
    assert constraint.anchor.experimental.uncertainty == pytest.approx(
        exp["uncertainty"]
    )
    assert constraint.anchor.experimental.units == exp["units"]
    assert constraint.anchor.experimental.snapshot_path == exp["snapshot_path"]
    assert constraint.anchor.historical_average.value == pytest.approx(hist["value"])
    assert constraint.anchor.historical_average.uncertainty == pytest.approx(
        hist["uncertainty"]
    )
    assert constraint.anchor.standard_model.value == pytest.approx(sm["value"])
    assert constraint.anchor.standard_model.uncertainty == pytest.approx(
        sm["uncertainty"]
    )
    assert constraint.anchor.standard_model.source_url == sm["source_url"]
    assert constraint.anchor.budget_band.central_residual == pytest.approx(central)
    assert constraint.anchor.budget_band.combined_sigma == pytest.approx(combined)
    assert constraint.anchor.budget == pytest.approx(central + combined)
    assert constraint.anchor.budget == pytest.approx(5.95465734053883e-10)


def test_missing_sm_uncertainty_anchor_fails_loudly():
    constraint = fcc.get(_PID)
    bad_sm = replace(constraint.anchor.standard_model, uncertainty=None)

    with pytest.raises(AnchorError, match="standard_model uncertainty is required"):
        b005_module._build_budget_band(
            process_id=_PID,
            experimental=constraint.anchor.experimental,
            standard_model=bad_sm,
        )


def test_evaluate_without_input_degrades_gracefully():
    constraint = fcc.get(_PID)
    result = constraint.evaluate(point_builder.empty_point())
    old_style_result = constraint.evaluate(
        point_builder.build_from_quark_couplings(_bs_couplings(left=1.0e-2))
    )

    assert result.process_id == _PID
    assert result.passes is True
    assert result.predicted is None
    assert result.ratio is None
    assert result.experimental == pytest.approx(constraint.anchor.value)
    assert result.sm_prediction == pytest.approx(constraint.sm_result.branching_fraction)
    assert result.budget == pytest.approx(constraint.anchor.budget)
    assert result.diagnostics["missing_extra"] == "rs_semileptonic_wilsons"
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]
    assert old_style_result.passes is True
    assert old_style_result.predicted is None
    assert old_style_result.diagnostics["evaluated"] is False
    assert old_style_result.diagnostics["missing_extra"] == "rs_semileptonic_wilsons"
    assert old_style_result.diagnostics["legacy_quark_mass_basis_couplings_present"] is True


def test_sm_limit_branching_fraction_matches_independent_formula_and_anchor():
    constraint = fcc.get(_PID)
    result = constraint.evaluate(sm_limit_rare_b_point())
    manual = _manual_sm_bs_mumu(constraint.sm_inputs)

    assert result.predicted == pytest.approx(manual)
    assert result.predicted == pytest.approx(3.6483515801441554e-09)
    assert result.predicted == pytest.approx(result.sm_prediction)
    assert abs(result.predicted - constraint.anchor.sm_value) < (
        constraint.anchor.standard_model.uncertainty
    )
    assert result.diagnostics["sm_anchor_branching_fraction"] == pytest.approx(
        constraint.anchor.sm_value
    )
    assert result.diagnostics["time_integrated_factor"] == pytest.approx(
        1.0 / (1.0 - constraint.sm_inputs.bs.width_difference_y)
    )
    assert result.diagnostics["a_delta_gamma"] == pytest.approx(1.0)
    assert result.diagnostics["c9_np"] == pytest.approx(0.0j)
    assert result.diagnostics["c10_np"] == pytest.approx(0.0j)
    assert result.diagnostics["wilson_prefactor_reused"] is False
    assert result.diagnostics["second_mkk_suppression_applied"] is False
    assert result.passes is True


def test_np_time_integration_uses_amplitude_dependent_a_delta_gamma():
    inputs = fcc.get(_PID).sm_inputs
    y_s = inputs.bs.width_difference_y
    wilsons = _synthetic_wilsons(
        inputs,
        c10_np=-inputs.c10_sm + 1.0 + 1.0j,
        cs_np=0.01j,
    )

    result = evaluate_bq_to_mumu(wilsons, transition="b_s", inputs=inputs)
    pseudoscalar = result.pseudoscalar_amplitude
    scalar = result.scalar_amplitude
    expected_a_delta_gamma = (
        (pseudoscalar * pseudoscalar).real - (scalar * scalar).real
    ) / (abs(pseudoscalar) ** 2 + abs(scalar) ** 2)
    expected_time_factor = (1.0 + expected_a_delta_gamma * y_s) / (1.0 - y_s**2)

    assert result.a_delta_gamma == pytest.approx(expected_a_delta_gamma)
    assert result.diagnostics["a_delta_gamma"] == pytest.approx(
        expected_a_delta_gamma
    )
    assert abs(result.a_delta_gamma - 1.0) > 0.1
    assert result.time_integrated_factor == pytest.approx(expected_time_factor)
    assert result.branching_fraction == pytest.approx(
        result.prompt_branching_fraction * result.time_integrated_factor
    )
    assert result.time_integrated_factor != pytest.approx(1.0 / (1.0 - y_s))
    assert (1.0 / (1.0 + y_s)) <= result.time_integrated_factor <= (
        1.0 / (1.0 - y_s)
    )

    lower_edge = evaluate_bq_to_mumu(
        _synthetic_wilsons(inputs, c10_np=-inputs.c10_sm + 1.0j),
        transition="b_s",
        inputs=inputs,
    )
    assert lower_edge.a_delta_gamma == pytest.approx(-1.0)
    assert lower_edge.time_integrated_factor == pytest.approx(1.0 / (1.0 + y_s))


def test_evaluate_runs_end_to_end_with_real_finite_fields_and_complex_diagnostics():
    point = sample_rare_b_point()
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
        "left_qb_coupling",
        "right_qb_coupling",
        "lambda_t",
        "c10_total",
        "c10_leptonic_np",
        "c9_effective_np",
        "c9_np",
        "c10_np",
    ):
        assert isinstance(result.diagnostics[key], complex)
    for key in (
        "m_kk_gev",
        "matching_scale_gev",
        "muon_vector_delta",
        "muon_axial_delta",
        "prompt_branching_fraction",
        "time_integrated_factor",
        "budget_combined_sigma",
        "budget_central_residual",
        "total_minus_experiment_ratio_to_combined_sigma",
    ):
        assert isinstance(result.diagnostics[key], float)
        assert math.isfinite(result.diagnostics[key])
    assert result.diagnostics["c9_does_not_enter_leptonic_rate"] is True
    assert result.diagnostics["rs_semileptonic_wilsons_present"] is True
    assert result.diagnostics["rs_semileptonic_matching_status"].endswith(
        "no_second_1_over_M_KK_squared"
    )
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]


@pytest.mark.parametrize(
    ("scale", "expected_pass"),
    [
        (1.0, True),
        (1.0e5, False),
    ],
)
def test_safe_point_passes_and_large_np_point_fails(
    scale: float,
    expected_pass: bool,
):
    constraint = fcc.get(_PID)
    point = scaled_rare_b_point(scale=scale)
    result = constraint.evaluate(point)
    direct = evaluate_bq_to_mumu(
        core_wilsons_from_rs_coeff(
            rs_coeff(point, transition="b_s", lepton="mu"),
            matching_scale_gev=point.extras["kk_ew_mass_gev"],
        ),
        transition="b_s",
        inputs=constraint.sm_inputs,
    )

    assert result.passes is expected_pass
    assert result.predicted == pytest.approx(direct.branching_fraction)
    assert result.sm_prediction == pytest.approx(direct.sm_branching_fraction)
    assert result.ratio == pytest.approx(
        abs(direct.np_shift_branching_fraction) / constraint.anchor.budget
    )
    if expected_pass:
        assert result.ratio < 1.0
    else:
        assert result.ratio > 10.0


def test_optional_kk_ew_mass_extra_changes_matching_scale():
    default_point = scaled_rare_b_point(scale=20.0)
    ew_point = scaled_rare_b_point(scale=20.0, kk_ew_mass_gev=6000.0)
    default_result = fcc.get(_PID).evaluate(default_point)
    ew_result = fcc.get(_PID).evaluate(ew_point)

    assert default_result.diagnostics["m_kk_gev"] == pytest.approx(3000.0)
    assert ew_result.diagnostics["m_kk_gev"] == pytest.approx(6000.0)
    assert ew_result.diagnostics["kk_ew_mass_extra_used"] is True
    assert ew_result.predicted == pytest.approx(default_result.predicted)
    assert ew_result.diagnostics["c10_leptonic_np"] == pytest.approx(
        default_result.diagnostics["c10_leptonic_np"]
    )
    assert ew_result.diagnostics["second_mkk_suppression_applied"] is False


def test_evaluate_is_pure_and_deterministic():
    point = sample_rare_b_point()
    before_bundle = point.extras["rs_semileptonic_wilsons"]
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    assert point.extras["rs_semileptonic_wilsons"] is before_bundle
