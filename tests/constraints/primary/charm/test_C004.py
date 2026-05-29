"""Production tests for C004 (D0 -> mu+ mu-)."""

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

_PID = "C004"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "charm" / "C004.yaml"


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


def _manual_sm_sd_d0_mumu(inputs) -> float:
    matrix = ckm_like_unitary(
        RotationParameters(
            theta12=inputs.theta12,
            theta13=inputs.theta13,
            theta23=inputs.theta23,
            delta=inputs.delta,
        )
    )
    meson = inputs.d0
    muon = inputs.muon
    lambda_b = complex(np.conjugate(matrix[1, 2]) * matrix[0, 2])
    tau_gev_inv = meson.lifetime_ps * 1.0e-12 / inputs.hbar_gev_s
    beta = math.sqrt(
        1.0
        - 4.0 * muon.mass_gev**2 / meson.meson_mass_gev**2
    )
    return float(
        tau_gev_inv
        * inputs.gf_gev_minus2**2
        * inputs.alpha_em_mz**2
        / (16.0 * math.pi**3)
        * meson.decay_constant_gev**2
        * meson.meson_mass_gev
        * muon.mass_gev**2
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
    assert constraint.observable == "BR(D0 -> mu+ mu-) short-distance"


def test_anchor_matches_yaml_and_budget():
    constraint = fcc.get(_PID)
    pdg = _yaml_pdg_block()
    current = pdg["canonical_current_limit"]
    cms = pdg["cms_current_95cl_limit"]
    lhcb = pdg["lhcb_previous_limit"]
    ld = pdg["standard_model_long_distance_context"]
    expected_vmd = (
        ld["relation_to_d0_gammagamma_factor"]
        * ld["d0_gammagamma_vmd_estimate"]
    )

    assert constraint.anchor.current_limit.value == pytest.approx(current["value"])
    assert constraint.anchor.current_limit.confidence_level == pytest.approx(
        current["confidence_level"]
    )
    assert constraint.anchor.current_limit.source_url == current["source_url"]
    assert constraint.anchor.cms_95_limit.value == pytest.approx(cms["value"])
    assert constraint.anchor.cms_95_limit.confidence_level == pytest.approx(
        cms["confidence_level"]
    )
    assert constraint.anchor.lhcb_previous_limit.value == pytest.approx(lhcb["value"])
    assert constraint.anchor.lhcb_previous_limit.confidence_level == pytest.approx(
        lhcb["confidence_level"]
    )
    assert constraint.anchor.long_distance_context.relation_to_d0_gammagamma_factor == (
        pytest.approx(ld["relation_to_d0_gammagamma_factor"])
    )
    assert constraint.anchor.long_distance_context.quoted_minimum_branching_fraction == (
        pytest.approx(ld["quoted_minimum_branching_fraction"])
    )
    assert constraint.anchor.long_distance_context.d0_gammagamma_vmd_estimate == (
        pytest.approx(ld["d0_gammagamma_vmd_estimate"])
    )
    assert constraint.anchor.long_distance_context.vmd_implied_branching_fraction == (
        pytest.approx(expected_vmd)
    )
    assert constraint.anchor.budget == pytest.approx(current["value"])
    assert constraint.anchor.budget == pytest.approx(2.1e-9)


def test_anchor_loading_fails_loudly_for_missing_candidate_or_value():
    with pytest.raises(fcc.AnchorError, match="none of the expected anchor keys"):
        fcc.load_anchor(_PID, family="charm", candidates=("missing_c004_anchor",))
    with pytest.raises(fcc.AnchorError, match="has no 'missing_limit' field"):
        fcc.load_anchor(
            _PID,
            family="charm",
            candidates=("canonical_current_limit",),
            value_key="missing_limit",
        )


def test_evaluate_without_input_degrades_gracefully():
    constraint = fcc.get(_PID)
    result = constraint.evaluate(point_builder.empty_point())

    assert result.process_id == _PID
    assert result.passes is True
    assert result.predicted is None
    assert result.ratio is None
    assert result.experimental == pytest.approx(constraint.anchor.value)
    assert result.sm_prediction == pytest.approx(constraint.sm_result.branching_fraction)
    assert result.budget == pytest.approx(constraint.anchor.budget)
    assert result.diagnostics["missing_extra"] == "quark_mass_basis_couplings"
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]


def test_sm_sd_limit_matches_manual_formula_and_yaml_ld_context():
    constraint = fcc.get(_PID)
    couplings = _uc_couplings(left=0.0j, right=0.0j)
    point = point_builder.build_from_quark_couplings(couplings)
    result = constraint.evaluate(point)
    manual = _manual_sm_sd_d0_mumu(constraint.sm_inputs)
    pdg = _yaml_pdg_block()
    ld = pdg["standard_model_long_distance_context"]
    expected_vmd = (
        ld["relation_to_d0_gammagamma_factor"]
        * ld["d0_gammagamma_vmd_estimate"]
    )

    assert result.predicted == pytest.approx(manual)
    assert result.predicted == pytest.approx(0.0)
    assert result.predicted == pytest.approx(result.sm_prediction)
    assert result.diagnostics["sm_long_distance_quoted_minimum"] == pytest.approx(
        ld["quoted_minimum_branching_fraction"]
    )
    assert result.diagnostics["sm_long_distance_vmd_implied_branching_fraction"] == (
        pytest.approx(expected_vmd)
    )
    assert expected_vmd == pytest.approx(9.45e-13)
    assert result.diagnostics["long_distance_not_subtracted"] is True
    assert result.passes is True


def test_evaluate_runs_end_to_end_with_real_finite_fields_and_complex_diagnostics():
    couplings = _uc_couplings(left=1.0e-2 + 0.2e-2j, right=0.5e-2j)
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
        "short_distance_plus_vmd_long_distance",
        "short_distance_plus_vmd_ratio_to_limit",
    ):
        assert isinstance(result.diagnostics[key], float)
        assert math.isfinite(result.diagnostics[key])
    assert result.diagnostics["c9_does_not_enter_leptonic_rate"] is True
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]


@pytest.mark.parametrize(
    ("couplings", "expected_pass"),
    [
        (_uc_couplings(left=1.0e-2), True),
        (_uc_couplings(left=5.0), False),
    ],
)
def test_safe_point_passes_and_large_np_point_fails(
    couplings: QuarkMassBasisCouplings,
    expected_pass: bool,
):
    constraint = fcc.get(_PID)
    point = point_builder.build_from_quark_couplings(couplings)
    result = constraint.evaluate(point)
    direct = evaluate_d0_to_ll(couplings, lepton="mu")

    assert result.passes is expected_pass
    assert result.predicted == pytest.approx(direct.branching_fraction)
    assert result.sm_prediction == pytest.approx(direct.sm_branching_fraction)
    assert result.ratio == pytest.approx(
        direct.branching_fraction / constraint.anchor.budget
    )
    if expected_pass:
        assert result.ratio < 1.0
    else:
        assert result.ratio > 10.0


def test_optional_kk_ew_mass_extra_changes_matching_scale():
    couplings = _uc_couplings(left=1.0)
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
    assert abs(ew_result.diagnostics["c10_leptonic_np"]) == pytest.approx(
        abs(default_result.diagnostics["c10_leptonic_np"]) / 4.0
    )


def test_evaluate_is_pure_and_deterministic():
    couplings = _uc_couplings(left=1.0e-2 + 0.2e-2j, right=0.5e-2j)
    before_left_up = couplings.left_up.copy()
    before_right_up = couplings.right_up.copy()
    point = point_builder.build_from_quark_couplings(couplings)
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    np.testing.assert_array_equal(couplings.left_up, before_left_up)
    np.testing.assert_array_equal(couplings.right_up, before_right_up)
