"""Production tests for B006 (B0 -> mu+ mu-)."""

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
from flavor_catalog_constraints.primary.beauty import B006 as b006_module
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

_PID = "B006"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "beauty" / "B006.yaml"


def _yaml_pdg_block():
    with open(_SIDECAR) as handle:
        return yaml.safe_load(handle)["pdg_or_equivalent"]


def _input_measurement(source: str):
    for entry in _yaml_pdg_block()["input_measurements"]:
        if entry["source"] == source:
            return entry
    raise AssertionError(f"no B006 input_measurements entry with source={source!r}")


def _bd_couplings(
    left: complex,
    right: complex = 0.0j,
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


def _manual_sm_bd_mumu(inputs) -> float:
    matrix = ckm_like_unitary(
        RotationParameters(
            theta12=inputs.theta12,
            theta13=inputs.theta13,
            theta23=inputs.theta23,
            delta=inputs.delta,
        )
    )
    meson = inputs.bd
    lambda_t = complex(matrix[2, 2] * np.conjugate(matrix[2, 0]))
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
    factors = ckm_factors("b_d", inputs)
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
        transition_key="b_d",
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
    assert constraint.observable == "BR(B0 -> mu+ mu-)"


def test_anchor_matches_yaml_and_budget_band():
    constraint = fcc.get(_PID)
    pdg = _yaml_pdg_block()
    exp = pdg["canonical_experimental_limit"]
    hflav = pdg["hflav_equivalent_ratio"]
    sm = pdg["standard_model_prediction"]
    cms = _input_measurement("CMS2023:BdMuMu")
    lhcb = _input_measurement("LHCb2022:BdMuMu")
    active_limit = min(
        exp["upper_limit"],
        cms["upper_limit_90cl"],
        lhcb["upper_limit_90cl"],
    )

    assert constraint.anchor.experimental_limit.value == pytest.approx(
        exp["upper_limit"]
    )
    assert constraint.anchor.experimental_limit.units == exp["units"]
    assert constraint.anchor.experimental_limit.snapshot_path == exp["snapshot_path"]
    assert constraint.anchor.experimental_limit.confidence_level == pytest.approx(
        exp["confidence_level"]
    )
    assert constraint.anchor.hflav_ratio_limit.value == pytest.approx(
        hflav["upper_limit"]
    )
    assert constraint.anchor.hflav_ratio_limit.units == hflav["units"]
    assert constraint.anchor.standard_model.value == pytest.approx(sm["value"])
    assert constraint.anchor.standard_model.uncertainty == pytest.approx(
        sm["uncertainty"]
    )
    assert constraint.anchor.standard_model.source_url == sm["source_url"]
    assert constraint.anchor.cms_limit.value == pytest.approx(cms["upper_limit_90cl"])
    assert constraint.anchor.cms_limit.source_url == cms["source_url"]
    assert constraint.anchor.lhcb_limit.value == pytest.approx(
        lhcb["upper_limit_90cl"]
    )
    assert constraint.anchor.lhcb_limit.source_url == lhcb["source_url"]
    assert constraint.anchor.budget_band.active_limit == pytest.approx(active_limit)
    assert constraint.anchor.budget == pytest.approx(active_limit - sm["value"])
    assert constraint.anchor.budget == pytest.approx(4.4e-11)
    assert constraint.anchor.budget_band.limit_minus_formula_sm == pytest.approx(
        active_limit - constraint.sm_result.branching_fraction
    )


def test_b006_anchors_route_through_scaffold_load_anchor(monkeypatch):
    calls: list[tuple[str, ...]] = []
    original_load_anchor = b006_module.load_anchor

    def spy_load_anchor(*args, **kwargs):
        calls.append(tuple(kwargs["candidates"]))
        return original_load_anchor(*args, **kwargs)

    monkeypatch.setattr(b006_module, "load_anchor", spy_load_anchor)
    anchor = b006_module._load_bd_mumu_anchor(
        _PID,
        formula_sm=fcc.get(_PID).sm_result.branching_fraction,
    )

    assert calls == [
        ("canonical_experimental_limit",),
        ("hflav_equivalent_ratio",),
        ("standard_model_prediction",),
        ("input_measurements[0]",),
        ("input_measurements[1]",),
    ]
    assert anchor.value == pytest.approx(fcc.get(_PID).anchor.value)


def test_limit_below_sm_anchor_fails_loudly():
    constraint = fcc.get(_PID)
    bad_sm = replace(
        constraint.anchor.standard_model,
        value=constraint.anchor.budget_band.active_limit * 1.01,
    )

    with pytest.raises(AnchorError, match="limit-SM budget"):
        b006_module._build_budget_band(
            process_id=_PID,
            experimental_limit=constraint.anchor.experimental_limit,
            cms_limit=constraint.anchor.cms_limit,
            lhcb_limit=constraint.anchor.lhcb_limit,
            standard_model=bad_sm,
            formula_sm=constraint.sm_result.branching_fraction,
        )


def test_evaluate_without_input_degrades_gracefully():
    constraint = fcc.get(_PID)
    result = constraint.evaluate(point_builder.empty_point())
    old_style_result = constraint.evaluate(
        point_builder.build_from_quark_couplings(_bd_couplings(left=1.0e-2))
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
    manual = _manual_sm_bd_mumu(constraint.sm_inputs)

    assert result.predicted == pytest.approx(manual)
    assert result.predicted == pytest.approx(1.0466846288738303e-10)
    assert result.predicted == pytest.approx(result.sm_prediction)
    assert abs(result.predicted - constraint.anchor.sm_value) < (
        constraint.anchor.standard_model.uncertainty
    )
    assert result.diagnostics["sm_anchor_branching_fraction"] == pytest.approx(
        constraint.anchor.sm_value
    )
    assert result.diagnostics["lambda_t"] == pytest.approx(
        ckm_factors("b_d", constraint.sm_inputs).lambda_t
    )
    assert result.diagnostics["time_integrated_factor"] == pytest.approx(1.0)
    assert result.diagnostics["a_delta_gamma"] == pytest.approx(1.0)
    assert result.diagnostics["c9_np"] == pytest.approx(0.0j)
    assert result.diagnostics["c10_np"] == pytest.approx(0.0j)
    assert result.diagnostics["wilson_prefactor_reused"] is False
    assert result.diagnostics["second_mkk_suppression_applied"] is False
    assert result.ratio == pytest.approx(0.0)
    assert result.passes is True


def test_bd_time_integration_keeps_amplitude_a_delta_gamma_but_y_d_is_zero():
    inputs = fcc.get(_PID).sm_inputs
    wilsons = _synthetic_wilsons(
        inputs,
        c10_np=-inputs.c10_sm + 1.0 + 1.0j,
        cs_np=0.01j,
    )

    result = evaluate_bq_to_mumu(wilsons, transition="b_d", inputs=inputs)
    pseudoscalar = result.pseudoscalar_amplitude
    scalar = result.scalar_amplitude
    expected_a_delta_gamma = (
        (pseudoscalar * pseudoscalar).real - (scalar * scalar).real
    ) / (abs(pseudoscalar) ** 2 + abs(scalar) ** 2)

    assert inputs.bd.width_difference_y == pytest.approx(0.0)
    assert result.a_delta_gamma == pytest.approx(expected_a_delta_gamma)
    assert abs(result.a_delta_gamma - 1.0) > 0.1
    assert result.time_integrated_factor == pytest.approx(1.0)
    assert result.branching_fraction == pytest.approx(result.prompt_branching_fraction)


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
        "hard_veto_np_shift_budget",
        "total_limit_ratio",
        "upward_excess_over_sm_anchor",
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
        (1.0e-3, True),
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
            rs_coeff(point, transition="b_d", lepton="mu"),
            matching_scale_gev=point.extras["kk_ew_mass_gev"],
        ),
        transition="b_d",
        inputs=constraint.sm_inputs,
    )
    expected_ratio = (
        max(0.0, direct.branching_fraction - constraint.anchor.sm_value)
        / constraint.anchor.budget
    )

    assert result.passes is expected_pass
    assert result.predicted == pytest.approx(direct.branching_fraction)
    assert result.sm_prediction == pytest.approx(direct.sm_branching_fraction)
    assert result.ratio == pytest.approx(expected_ratio)
    if expected_pass:
        assert result.ratio < 1.0
        assert result.predicted < constraint.anchor.budget_band.active_limit
    else:
        assert result.ratio > 1.0
        assert result.predicted > constraint.anchor.budget_band.active_limit


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
