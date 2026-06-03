"""Production tests for C005 (D0 -> e+ e-)."""

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
from quarkConstraints.model import RotationParameters, ckm_like_unitary
from quarkConstraints.rare_charm_dilepton import evaluate_d0_to_ll
from tests.rare_charm_phase3d_helpers import (
    core_dilepton_wilsons_from_rs_coeff,
    rare_charm_point_with_wilsons,
    rs_coeff,
    sample_rare_charm_point,
    sm_limit_rare_charm_point,
)

_PID = "C005"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "charm" / "C005.yaml"


def _yaml():
    with open(_SIDECAR) as handle:
        return yaml.safe_load(handle)


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


def _manual_sm_sd_d0_ee(inputs) -> float:
    matrix = ckm_like_unitary(
        RotationParameters(
            theta12=inputs.theta12,
            theta13=inputs.theta13,
            theta23=inputs.theta23,
            delta=inputs.delta,
        )
    )
    meson = inputs.d0
    electron = inputs.electron
    lambda_b = complex(np.conjugate(matrix[1, 2]) * matrix[0, 2])
    tau_gev_inv = meson.lifetime_ps * 1.0e-12 / inputs.hbar_gev_s
    beta = math.sqrt(
        1.0
        - 4.0 * electron.mass_gev**2 / meson.meson_mass_gev**2
    )
    return float(
        tau_gev_inv
        * inputs.gf_gev_minus2**2
        * inputs.alpha_em_mz**2
        / (16.0 * math.pi**3)
        * meson.decay_constant_gev**2
        * meson.meson_mass_gev
        * electron.mass_gev**2
        * abs(lambda_b) ** 2
        * beta
        * abs(inputs.c10_sm) ** 2
    )


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.HARD
    assert constraint.family == "charm"
    assert constraint.observable == "BR(D0 -> e+ e-) short-distance"


def test_anchor_matches_yaml_and_budget():
    constraint = fcc.get(_PID)
    data = _yaml()
    pdg = data["pdg_or_equivalent"]
    current = pdg["canonical_current_limit"]
    belle = pdg["belle_current_result"]
    babar = pdg["babar_independent_result"]
    theory = data["theory_context"]["rare_charm_model_dependence"]
    rs = data["paper_era_reference"]["rs_baseline"]

    assert constraint.anchor.current_limit.value == pytest.approx(current["value"])
    assert constraint.anchor.current_limit.confidence_level == pytest.approx(
        current["confidence_level"]
    )
    assert constraint.anchor.current_limit.source_url == current["source_url"]
    assert constraint.anchor.belle_current_limit.value == pytest.approx(belle["value"])
    assert constraint.anchor.belle_current_limit.confidence_level == pytest.approx(
        belle["confidence_level"]
    )
    assert constraint.anchor.babar_independent_limit.value == pytest.approx(
        babar["value"]
    )
    assert constraint.anchor.babar_independent_limit.confidence_level == pytest.approx(
        babar["confidence_level"]
    )
    assert constraint.anchor.theory_context.source == theory["source"]
    assert constraint.anchor.theory_context.topic == theory["topic"]
    assert constraint.anchor.rs_baseline.generic_rs_kk_gluon_scale_tev == pytest.approx(
        rs["generic_rs_kk_gluon_scale_TeV"]
    )
    assert (
        constraint.anchor.rs_baseline.composite_pseudo_goldstone_kk_gluon_scale_tev
        == pytest.approx(rs["composite_pseudo_goldstone_kk_gluon_scale_TeV"])
    )
    assert constraint.anchor.budget == pytest.approx(current["value"])
    assert constraint.anchor.budget == pytest.approx(7.9e-8)


def test_anchor_loading_fails_loudly_for_missing_candidate_or_value():
    with pytest.raises(fcc.AnchorError, match="none of the expected anchor keys"):
        fcc.load_anchor(_PID, family="charm", candidates=("missing_c005_anchor",))
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
    assert result.sm_prediction == pytest.approx(constraint.sm_result.branching_fraction)
    assert result.budget == pytest.approx(constraint.anchor.budget)
    assert result.diagnostics["evaluated"] is False
    assert result.diagnostics["missing_extra"] == "rs_semileptonic_wilsons"
    assert result.diagnostics["legacy_quark_mass_basis_couplings_present"] is True
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]


def test_sm_sd_limit_matches_manual_formula_and_zero_short_distance_policy():
    constraint = fcc.get(_PID)
    point = sm_limit_rare_charm_point()
    result = constraint.evaluate(point)
    manual = _manual_sm_sd_d0_ee(constraint.sm_inputs)

    assert result.predicted == pytest.approx(manual)
    assert result.predicted == pytest.approx(0.0)
    assert result.predicted == pytest.approx(result.sm_prediction)
    assert result.diagnostics["evaluated"] is True
    assert result.diagnostics["rs_semileptonic_wilsons_present"] is True
    assert abs(result.diagnostics["c10_leptonic_np"]) < 1.0e-12
    assert result.diagnostics["sm_long_distance_numeric_anchor_available"] is False
    assert result.diagnostics["long_distance_not_subtracted"] is True
    assert result.passes is True


def test_electron_mode_prediction_matches_independent_core_wilson_consumption():
    constraint = fcc.get(_PID)
    point = sample_rare_charm_point()
    result = constraint.evaluate(point)
    direct = evaluate_d0_to_ll(
        core_dilepton_wilsons_from_rs_coeff(rs_coeff(point, lepton="e")),
        lepton="e",
    )

    assert result.predicted == pytest.approx(
        direct.branching_fraction,
        rel=1e-12,
        abs=0.0,
    )
    assert abs(result.diagnostics["c10_leptonic_np"]) > 0.0
    assert result.ratio == pytest.approx(result.predicted / constraint.anchor.budget)


def test_electron_prediction_has_expected_helicity_suppression_vs_muon():
    constraint = fcc.get(_PID)
    point = rare_charm_point_with_wilsons(lepton="e", c10_np=100.0)
    electron_result = constraint.evaluate(point)
    muon_direct = evaluate_d0_to_ll(
        core_dilepton_wilsons_from_rs_coeff(rs_coeff(point, lepton="e")),
        lepton="mu",
    )
    inputs = constraint.sm_inputs
    beta_e = math.sqrt(
        1.0
        - 4.0 * inputs.electron.mass_gev**2 / inputs.d0.meson_mass_gev**2
    )
    beta_mu = math.sqrt(
        1.0
        - 4.0 * inputs.muon.mass_gev**2 / inputs.d0.meson_mass_gev**2
    )
    expected_ratio = (
        inputs.electron.mass_gev**2
        * beta_e
        / (inputs.muon.mass_gev**2 * beta_mu)
    )

    assert electron_result.predicted / muon_direct.branching_fraction == pytest.approx(
        expected_ratio
    )
    assert expected_ratio == pytest.approx(
        electron_result.diagnostics["electron_to_muon_helicity_suppression_m2"]
        / beta_mu
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
        "lepton_vector_delta",
        "lepton_axial_delta",
        "np_shift_branching_fraction",
        "electron_to_muon_helicity_suppression_m2",
    ):
        assert isinstance(result.diagnostics[key], float)
        assert math.isfinite(result.diagnostics[key])
    assert result.diagnostics["evaluated"] is True
    assert result.diagnostics["rs_semileptonic_wilsons_present"] is True
    assert result.diagnostics["wilson_prefactor_reused"] is False
    assert result.diagnostics["second_mkk_suppression_applied"] is False
    assert result.diagnostics["c9_does_not_enter_leptonic_rate"] is True
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]


@pytest.mark.parametrize(
    ("point", "expected_pass"),
    [
        (rare_charm_point_with_wilsons(lepton="e", c10_np=100.0), True),
        (rare_charm_point_with_wilsons(lepton="e", c10_np=2.0e7), False),
    ],
)
def test_safe_point_passes_and_large_np_point_fails(
    point,
    expected_pass: bool,
):
    constraint = fcc.get(_PID)
    result = constraint.evaluate(point)
    direct = evaluate_d0_to_ll(
        core_dilepton_wilsons_from_rs_coeff(rs_coeff(point, lepton="e")),
        lepton="e",
    )

    assert result.passes is expected_pass
    assert result.predicted == pytest.approx(
        direct.branching_fraction,
        rel=1e-12,
        abs=0.0,
    )
    assert result.sm_prediction == pytest.approx(direct.sm_branching_fraction)
    assert result.ratio == pytest.approx(
        direct.branching_fraction / constraint.anchor.budget
    )
    if expected_pass:
        assert result.ratio < 1.0
    else:
        assert result.ratio > 10.0


def test_optional_kk_ew_mass_extra_is_diagnostic_only_no_second_suppression():
    default_point = rare_charm_point_with_wilsons(lepton="e", c10_np=1000.0)
    extras = dict(default_point.extras)
    extras["kk_ew_mass_gev"] = 6000.0
    ew_point = point_builder.make_point(**extras)
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
    point = sample_rare_charm_point()
    before_bundle = point.extras["rs_semileptonic_wilsons"]
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    assert point.extras["rs_semileptonic_wilsons"] is before_bundle
