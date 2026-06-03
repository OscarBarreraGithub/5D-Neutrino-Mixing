"""Production tests for C007 (D+ -> pi+ mu+ mu-)."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest
import yaml

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import point_builder
from flavor_catalog_constraints.base import ConstraintProtocol, Severity
from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.rs_semileptonic_wilsons import RSSemileptonicWilsonCoefficients
from tests.rare_charm_phase3d_helpers import (
    rare_charm_point_with_wilsons,
    rs_coeff,
    sample_rare_charm_point,
    sm_limit_rare_charm_point,
)

_PID = "C007"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "charm" / "C007.yaml"


def _yaml_pdg_block():
    with open(_SIDECAR) as handle:
        return yaml.safe_load(handle)["pdg_or_equivalent"]


def _uc_couplings(
    left: complex,
    right: complex = 0.0j,
    M_KK: float = 3000.0,
) -> QuarkMassBasisCouplings:
    """Minimal valid mass-basis couplings with only the u-c slot populated."""
    zeros = np.zeros((3, 3), dtype=np.complex128)
    left_up = zeros.copy()
    right_up = zeros.copy()
    left_up[0, 1] = left
    left_up[1, 0] = np.conj(left)
    right_up[0, 1] = right
    right_up[1, 0] = np.conj(right)
    return QuarkMassBasisCouplings(
        M_KK=M_KK,
        xi_KK=1.0,
        alpha_s=0.09,
        g_s=1.0,
        left_overlap=zeros,
        right_up_overlap=zeros,
        right_down_overlap=zeros,
        left_up=left_up,
        left_down=zeros.copy(),
        right_up=right_up,
        right_down=zeros.copy(),
    )


def _manual_fplus(q2: float, inputs) -> float:
    ff = inputs.form_factor
    x = q2 / ff.pole_mass_gev**2
    return ff.fplus_0 / ((1.0 - x) * (1.0 - ff.fplus_shape_a * x))


def _manual_fzero(q2: float, inputs) -> float:
    ff = inputs.form_factor
    x = q2 / ff.pole_mass_gev**2
    return ff.fplus_0 / (1.0 - x / ff.fzero_shape_b)


def _manual_kallen(a: float, b: float, c: float) -> float:
    return a * a + b * b + c * c - 2.0 * (a * b + a * c + b * c)


def _manual_smooth_dtopi_mumu_proxy(
    coeff: RSSemileptonicWilsonCoefficients,
    inputs,
) -> float:
    c9_semileptonic = complex(coeff.c9_np + coeff.c9p_np)
    c10_semileptonic = complex(
        inputs.rare_charm.c10_sm + coeff.c10_np + coeff.c10p_np
    )

    m_d = inputs.dplus_mass_gev
    m_pi = inputs.piplus_mass_gev
    m_mu = inputs.rare_charm.muon.mass_gev
    m_d2 = m_d * m_d
    m_pi2 = m_pi * m_pi
    q2_min = 4.0 * m_mu * m_mu
    q2_max = (m_d - m_pi) ** 2
    tau = inputs.dplus_lifetime_ps * 1.0e-12 / inputs.rare_charm.hbar_gev_s
    lambda_b = complex(coeff.lambda_ckm)
    gamma0 = (
        inputs.rare_charm.gf_gev_minus2**2
        * inputs.rare_charm.alpha_em_mz**2
        * abs(lambda_b) ** 2
        / ((4.0 * math.pi) ** 5 * m_d**3)
    )

    nodes, weights = np.polynomial.legendre.leggauss(inputs.quadrature_points)
    mid = 0.5 * (q2_max + q2_min)
    half_width = 0.5 * (q2_max - q2_min)
    total = 0.0
    for node, weight in zip(nodes, weights):
        q2 = mid + half_width * float(node)
        lam = _manual_kallen(m_d2, m_pi2, q2)
        beta2 = 1.0 - 4.0 * m_mu * m_mu / q2
        sqrt_lam = math.sqrt(max(0.0, lam))
        beta = math.sqrt(max(0.0, beta2))
        fplus = _manual_fplus(q2, inputs)
        fzero = _manual_fzero(q2, inputs)
        vector_amp = fplus * c9_semileptonic
        axial_amp = fplus * c10_semileptonic
        pseudoscalar_amp = (
            -m_mu
            * (fplus - (m_d2 - m_pi2) / q2 * (fzero - fplus))
            * c10_semileptonic
        )
        a_term = gamma0 * sqrt_lam * beta * (
            2.0 * q2 * abs(pseudoscalar_amp) ** 2
            + 0.5 * lam * (abs(axial_amp) ** 2 + abs(vector_amp) ** 2)
            + 4.0
            * m_mu
            * (m_d2 - m_pi2 + q2)
            * (axial_amp * np.conjugate(pseudoscalar_amp)).real
            + 8.0 * m_mu * m_mu * m_d2 * abs(axial_amp) ** 2
        )
        c_term = gamma0 * sqrt_lam * beta * (
            -0.5 * lam * beta2 * (abs(vector_amp) ** 2 + abs(axial_amp) ** 2)
        )
        total += float(weight) * 2.0 * tau * (a_term + c_term / 3.0)
    return float(half_width * total)


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.HARD
    assert constraint.family == "charm"
    assert constraint.observable == "BR(D+ -> pi+ mu+ mu-) smooth full-q2 proxy"


def test_anchor_matches_yaml_and_budget():
    constraint = fcc.get(_PID)
    pdg = _yaml_pdg_block()
    current = pdg["canonical_current_limit"]
    lhcb = pdg["lhcb_current_result"]
    previous = pdg["lhcb_previous_result"]
    scope = pdg["lhcb_2021_search_scope"]
    sm_scale = pdg["lhcb_2021_short_distance_sm_scale"]
    fplus0 = pdg["dtopi_form_factor_fplus0"]
    pole = pdg["dtopi_form_factor_pole_mass"]
    shape_a = pdg["dtopi_form_factor_fplus_shape_a"]
    shape_b = pdg["dtopi_form_factor_fzero_shape_b"]
    theory = pdg["theory_context"]
    rs = pdg["rs_baseline"]

    assert constraint.anchor.current_limit.value == pytest.approx(current["value"])
    assert constraint.anchor.current_limit.confidence_level == pytest.approx(
        current["confidence_level"]
    )
    assert constraint.anchor.current_limit.source_url == current["source_url"]
    assert constraint.anchor.lhcb_current_limit.value == pytest.approx(lhcb["value"])
    assert constraint.anchor.lhcb_current_limit.companion_95cl_limit == pytest.approx(
        lhcb["companion_95cl_limit"]
    )
    assert constraint.anchor.previous_nonresonant_limit.value == pytest.approx(
        previous["value"]
    )
    assert constraint.anchor.previous_nonresonant_limit.companion_95cl_limit == (
        pytest.approx(previous["companion_95cl_limit"])
    )
    assert constraint.anchor.search_scope.value == pytest.approx(scope["value"])
    assert constraint.anchor.short_distance_sm_scale.value == pytest.approx(
        sm_scale["value"]
    )
    assert constraint.anchor.short_distance_sm_scale.qualifier == sm_scale["qualifier"]
    assert constraint.anchor.form_factor.fplus0.value == pytest.approx(fplus0["value"])
    assert constraint.anchor.form_factor.fplus0.uncertainty == pytest.approx(
        fplus0["uncertainty"]
    )
    assert constraint.anchor.form_factor.fplus0.snapshot_path == fplus0["snapshot_path"]
    assert constraint.anchor.form_factor.pole_mass.value == pytest.approx(pole["value"])
    assert constraint.anchor.form_factor.pole_mass.uncertainty == pytest.approx(
        pole["uncertainty"]
    )
    assert constraint.anchor.form_factor.fplus_shape_a.value == pytest.approx(
        shape_a["value"]
    )
    assert constraint.anchor.form_factor.fplus_shape_a.uncertainty == pytest.approx(
        shape_a["uncertainty"]
    )
    assert constraint.anchor.form_factor.fzero_shape_b.value == pytest.approx(
        shape_b["value"]
    )
    assert constraint.anchor.form_factor.fzero_shape_b.uncertainty == pytest.approx(
        shape_b["uncertainty"]
    )
    assert constraint.anchor.form_factor.normalization_uncertainty_fraction == (
        pytest.approx(2.0 * fplus0["uncertainty"] / fplus0["value"])
    )
    assert constraint.anchor.theory_context.long_distance_note == (
        theory["long_distance_note"]
    )
    assert constraint.anchor.rs_baseline.generic_rs_kk_gluon_scale_tev == pytest.approx(
        rs["generic_rs_kk_gluon_scale_TeV"]
    )
    assert constraint.anchor.budget == pytest.approx(current["value"])
    assert constraint.anchor.budget == pytest.approx(6.7e-8)
    assert constraint.anchor.sm_value == pytest.approx(1.0e-12)


def test_anchor_loading_fails_loudly_for_missing_candidate_or_value():
    with pytest.raises(fcc.AnchorError, match="none of the expected anchor keys"):
        fcc.load_anchor(_PID, family="charm", candidates=("missing_c007_anchor",))
    with pytest.raises(fcc.AnchorError, match="has no 'missing_limit' field"):
        fcc.load_anchor(
            _PID,
            family="charm",
            candidates=("canonical_current_limit",),
            value_key="missing_limit",
        )


def test_evaluate_without_input_degrades_gracefully():
    constraint = fcc.get(_PID)
    old_style_point = point_builder.build_from_quark_couplings(_uc_couplings(left=1.0))
    result = constraint.evaluate(old_style_point)

    assert result.process_id == _PID
    assert result.passes is True
    assert result.predicted is None
    assert result.ratio is None
    assert result.experimental == pytest.approx(constraint.anchor.value)
    assert result.sm_prediction == pytest.approx(constraint.anchor.sm_value)
    assert result.budget == pytest.approx(constraint.anchor.budget)
    assert result.diagnostics["evaluated"] is False
    assert result.diagnostics["missing_extra"] == "rs_semileptonic_wilsons"
    assert result.diagnostics["legacy_quark_mass_basis_couplings_present"] is True
    assert result.diagnostics["full_q2_proxy_constrained"] is True
    assert result.diagnostics["lhcb_nonresonant_dimuon_window_applied"] is False
    assert "full_kinematic_q2" in result.diagnostics["q2_treatment"]
    assert result.diagnostics["long_distance_resonance_dominated"] is True
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]


def test_zero_np_prediction_uses_yaml_short_distance_sm_scale():
    constraint = fcc.get(_PID)
    point = sm_limit_rare_charm_point()
    result = constraint.evaluate(point)

    assert constraint.smooth_sm_result.branching_fraction == pytest.approx(0.0)
    assert result.sm_prediction == pytest.approx(1.0e-12)
    assert result.predicted == pytest.approx(1.0e-12)
    assert result.ratio == pytest.approx(1.0e-12 / constraint.anchor.budget)
    assert result.diagnostics["evaluated"] is True
    assert result.diagnostics["rs_semileptonic_wilsons_present"] is True
    assert abs(result.diagnostics["c9_semileptonic_np"]) < 1.0e-12
    assert abs(result.diagnostics["c10_semileptonic_np"]) < 1.0e-12
    assert result.diagnostics["catalog_sm_short_distance_scale_added_incoherently"] is True
    assert result.diagnostics["long_distance_not_subtracted"] is True
    assert result.passes is True


def test_rigorous_wilson_path_matches_independent_manual_integration():
    constraint = fcc.get(_PID)
    point = sample_rare_charm_point()
    result = constraint.evaluate(point)
    manual = _manual_smooth_dtopi_mumu_proxy(
        rs_coeff(point, lepton="mu"),
        constraint.sd_inputs,
    )

    assert result.diagnostics["smooth_full_q2_np_proxy_branching_fraction"] == (
        pytest.approx(manual, rel=1e-12, abs=0.0)
    )
    assert manual > 0.0
    assert result.predicted == pytest.approx(manual + constraint.anchor.sm_value)
    assert result.ratio == pytest.approx(result.predicted / constraint.anchor.budget)
    assert result.diagnostics["fplus_0"] == pytest.approx(0.67)
    assert result.diagnostics["form_factor_pole_mass_gev"] == pytest.approx(1.90)
    assert result.diagnostics["form_factor_fplus0_anchor_block"] == (
        "dtopi_form_factor_fplus0"
    )
    assert result.diagnostics["form_factor_normalization_uncertainty_fraction"] == (
        pytest.approx(2.0 * 0.03 / 0.67)
    )
    assert result.diagnostics["budget_form_factor_uncertainty_policy"].startswith(
        "f_+(0) normalization uncertainty"
    )


def test_evaluate_runs_end_to_end_with_real_finite_fields_and_complex_diagnostics():
    point = sample_rare_charm_point()
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
        "left_uc_coupling",
        "right_uc_coupling",
        "lambda_b",
        "c9_semileptonic_np",
        "c10_semileptonic_np",
        "c9_np",
        "c10_np",
    ):
        assert isinstance(result.diagnostics[key], complex)
    for key in (
        "m_kk_gev",
        "matching_scale_gev",
        "q2_min_gev2",
        "q2_max_gev2",
        "smooth_full_q2_np_proxy_branching_fraction",
        "predicted_full_q2_proxy_total",
    ):
        assert isinstance(result.diagnostics[key], float)
        assert math.isfinite(result.diagnostics[key])
    assert result.diagnostics["evaluated"] is True
    assert result.diagnostics["rs_semileptonic_wilsons_present"] is True
    assert result.diagnostics["wilson_prefactor_reused"] is False
    assert result.diagnostics["second_mkk_suppression_applied"] is False
    assert result.diagnostics["semileptonic_primed_combination_is_plus"] is True
    assert result.diagnostics["resonance_amplitudes_included"] is False
    assert result.diagnostics["full_q2_proxy_constrained"] is True
    assert result.diagnostics["lhcb_nonresonant_dimuon_window_applied"] is False
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]


@pytest.mark.parametrize(
    ("point", "expected_pass"),
    [
        (rare_charm_point_with_wilsons(lepton="mu", c9_np=10.0), True),
        (rare_charm_point_with_wilsons(lepton="mu", c9_np=1.0e4), False),
    ],
)
def test_safe_point_passes_and_large_np_point_fails(
    point,
    expected_pass: bool,
):
    constraint = fcc.get(_PID)
    result = constraint.evaluate(point)
    manual = _manual_smooth_dtopi_mumu_proxy(
        rs_coeff(point, lepton="mu"),
        constraint.sd_inputs,
    )

    assert result.passes is expected_pass
    assert result.predicted == pytest.approx(manual + constraint.anchor.sm_value)
    assert result.ratio == pytest.approx(result.predicted / constraint.anchor.budget)
    if expected_pass:
        assert result.ratio < 1.0
    else:
        assert result.ratio > 1.0


def test_optional_kk_ew_mass_extra_is_diagnostic_only_no_second_suppression():
    default_point = rare_charm_point_with_wilsons(lepton="mu", c9_np=100.0)
    extras = dict(default_point.extras)
    extras["kk_ew_mass_gev"] = 6000.0
    ew_point = point_builder.make_point(**extras)
    default_result = fcc.get(_PID).evaluate(default_point)
    ew_result = fcc.get(_PID).evaluate(ew_point)

    assert default_result.diagnostics["m_kk_gev"] == pytest.approx(3000.0)
    assert ew_result.diagnostics["m_kk_gev"] == pytest.approx(6000.0)
    assert ew_result.diagnostics["kk_ew_mass_extra_used"] is True
    assert ew_result.predicted == pytest.approx(default_result.predicted)
    assert ew_result.diagnostics["c9_semileptonic_np"] == pytest.approx(
        default_result.diagnostics["c9_semileptonic_np"]
    )
    assert ew_result.diagnostics["second_mkk_suppression_applied"] is False


def test_evaluate_is_pure_and_deterministic():
    constraint = fcc.get(_PID)
    point = sample_rare_charm_point()
    before_bundle = point.extras["rs_semileptonic_wilsons"]

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    assert point.extras["rs_semileptonic_wilsons"] is before_bundle
