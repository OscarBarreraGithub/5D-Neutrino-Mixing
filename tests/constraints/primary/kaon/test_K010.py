"""Production tests for K010 (K_S -> pi0 e+ e-)."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest
import yaml

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import point_builder
from flavor_catalog_constraints.base import ConstraintProtocol, Severity
from flavor_catalog_constraints.physics_adapters.rare_kaon_dilepton_ks import (
    KSHORT_PI0EE_A_S_BRANCHING_COEFFICIENT,
)
from flavor_catalog_constraints.primary.kaon import K010 as k010_module
from quarkConstraints.couplings import QuarkMassBasisCouplings
from tests.rare_kaon_phase3d_helpers import (
    core_y7_wilsons_from_rs_coeff,
    rs_coeff,
    sample_rare_kaon_point,
    scaled_rare_kaon_point,
    sm_limit_rare_kaon_point,
)

_PID = "K010"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "kaon" / "K010.yaml"


def _yaml():
    with open(_SIDECAR) as handle:
        return yaml.safe_load(handle)


def _value_entry(observable: str):
    for entry in _yaml()["pdg_or_equivalent"]["values"]:
        if entry["observable"] == observable:
            return entry
    raise AssertionError(f"no K010 value entry for {observable!r}")


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


def _independent_ks_pi0ee_full_rate(
    chpt,
    lambda_y7v_np_proxy: complex = 0.0j,
) -> dict[str, float | complex]:
    a_s_np_proxy_complex = lambda_y7v_np_proxy / 1.0e-4
    a_s_np_proxy = a_s_np_proxy_complex.real
    positive_a_s_effective = chpt.a_s_abs + a_s_np_proxy
    negative_a_s_effective = -chpt.a_s_abs + a_s_np_proxy
    positive_branching_fraction = (
        chpt.rate_coefficient * positive_a_s_effective**2
    )
    negative_branching_fraction = (
        chpt.rate_coefficient * negative_a_s_effective**2
    )
    return {
        "a_s_np_proxy_complex": a_s_np_proxy_complex,
        "a_s_np_proxy": a_s_np_proxy,
        "positive_a_s_effective": positive_a_s_effective,
        "negative_a_s_effective": negative_a_s_effective,
        "positive_a_s_branching_fraction": positive_branching_fraction,
        "negative_a_s_branching_fraction": negative_branching_fraction,
        "branching_fraction": positive_branching_fraction,
        "sm_branching_fraction": chpt.rate_coefficient * chpt.a_s_abs**2,
    }


def _independent_budget_ratio(predicted: float, constraint) -> tuple[float, float]:
    budget = (
        constraint.anchor.full_rate.uncertainty_upper
        if predicted >= constraint.anchor.value
        else constraint.anchor.full_rate.uncertainty_lower
    )
    ratio = abs(predicted - constraint.anchor.value) / budget
    return budget, ratio


def _independent_conservative_branch(
    independent: dict[str, float | complex],
    constraint,
) -> dict[str, float | str]:
    branches = []
    for label in ("positive", "negative"):
        predicted = float(independent[f"{label}_a_s_branching_fraction"])
        budget, ratio = _independent_budget_ratio(predicted, constraint)
        branches.append(
            {
                "branch": label,
                "predicted": predicted,
                "budget": budget,
                "ratio": ratio,
                "a_s_effective": float(independent[f"{label}_a_s_effective"]),
            }
        )
    return min(branches, key=lambda item: float(item["ratio"]))


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.HARD
    assert constraint.family == "kaon"
    assert constraint.observable == "BR(K_S -> pi0 e+ e-)"


def test_anchor_matches_yaml_and_budget_band():
    constraint = fcc.get(_PID)
    pdg = _yaml()["pdg_or_equivalent"]
    partial = pdg["canonical_measurement"]
    full = pdg["extrapolated_total_rate"]
    observed = _value_entry(k010_module._OBSERVED_CANDIDATES_OBS)
    background = _value_entry(k010_module._EXPECTED_BACKGROUND_OBS)
    partial_up = math.sqrt(
        partial["uncertainty_stat_positive"] ** 2 + partial["uncertainty_syst"] ** 2
    )
    partial_down = math.sqrt(
        partial["uncertainty_stat_negative"] ** 2 + partial["uncertainty_syst"] ** 2
    )
    expected_a_s = math.sqrt(
        float(full["value"]) / KSHORT_PI0EE_A_S_BRANCHING_COEFFICIENT
    )

    assert constraint.anchor.partial_rate.value == pytest.approx(float(partial["value"]))
    assert constraint.anchor.partial_rate.uncertainty_upper == pytest.approx(partial_up)
    assert constraint.anchor.partial_rate.uncertainty_lower == pytest.approx(partial_down)
    assert constraint.anchor.partial_rate.source_url == partial["source_url"]
    assert constraint.anchor.full_rate.value == pytest.approx(float(full["value"]))
    assert constraint.anchor.full_rate.uncertainty_upper == pytest.approx(
        float(full["uncertainty_positive"])
    )
    assert constraint.anchor.full_rate.uncertainty_lower == pytest.approx(
        float(full["uncertainty_negative"])
    )
    assert constraint.anchor.full_rate.assumptions == full["assumptions"]
    assert constraint.anchor.observed_candidates.value == pytest.approx(
        float(observed["value"])
    )
    assert constraint.anchor.expected_background.value == pytest.approx(
        float(background["value"])
    )
    assert constraint.anchor.a_s_inputs.rate_coefficient == pytest.approx(
        KSHORT_PI0EE_A_S_BRANCHING_COEFFICIENT
    )
    assert constraint.anchor.a_s_inputs.a_s_abs == pytest.approx(expected_a_s)
    assert constraint.anchor.budget == pytest.approx(
        max(full["uncertainty_positive"], full["uncertainty_negative"])
    )
    assert constraint.anchor.budget_band.lower_edge == pytest.approx(
        full["value"] - full["uncertainty_negative"]
    )
    assert constraint.anchor.budget_band.upper_edge == pytest.approx(
        full["value"] + full["uncertainty_positive"]
    )


def test_k010_anchor_loud_fail_probe():
    with pytest.raises(k010_module.AnchorError):
        k010_module._load_value_anchor(_PID, "not a K010 observable")
    with pytest.raises(k010_module.AnchorError):
        k010_module._load_branching_subblock_anchor(_PID, "not_a_k010_block")


def test_evaluate_without_input_degrades_gracefully():
    constraint = fcc.get(_PID)
    result = constraint.evaluate(point_builder.empty_point())
    old_style_result = constraint.evaluate(
        point_builder.build_from_quark_couplings(_sd_couplings(left=1.0e-2))
    )

    assert result.process_id == _PID
    assert result.passes is True
    assert result.predicted is None
    assert result.ratio is None
    assert result.experimental == pytest.approx(fcc.get(_PID).anchor.value)
    assert result.sm_prediction == pytest.approx(
        fcc.get(_PID).sm_result.sm_branching_fraction
    )
    assert result.budget == pytest.approx(fcc.get(_PID).anchor.budget)
    assert result.diagnostics["evaluated"] is False
    assert result.diagnostics["missing_extra"] == "rs_semileptonic_wilsons"
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]
    assert old_style_result.passes is True
    assert old_style_result.predicted is None
    assert old_style_result.diagnostics["evaluated"] is False
    assert old_style_result.diagnostics["missing_extra"] == "rs_semileptonic_wilsons"
    assert old_style_result.diagnostics["legacy_quark_mass_basis_couplings_present"] is True


def test_sm_limit_branching_fraction_matches_independent_recomputation():
    constraint = fcc.get(_PID)
    point = sm_limit_rare_kaon_point()
    result = constraint.evaluate(point)
    independent = _independent_ks_pi0ee_full_rate(constraint.chpt_inputs)
    selected = _independent_conservative_branch(independent, constraint)

    assert result.predicted == pytest.approx(selected["predicted"])
    assert result.predicted == pytest.approx(result.sm_prediction)
    assert result.predicted == pytest.approx(constraint.anchor.full_rate.value)
    assert result.diagnostics["a_s_abs"] == pytest.approx(
        constraint.chpt_inputs.a_s_abs
    )
    assert result.diagnostics["a_s_effective"] == pytest.approx(
        selected["a_s_effective"]
    )
    assert result.diagnostics["positive_a_s_effective"] == pytest.approx(
        independent["positive_a_s_effective"]
    )
    assert result.diagnostics["negative_a_s_effective"] == pytest.approx(
        independent["negative_a_s_effective"]
    )
    assert result.diagnostics["positive_a_s_branching_fraction"] == pytest.approx(
        independent["positive_a_s_branching_fraction"]
    )
    assert result.diagnostics["negative_a_s_branching_fraction"] == pytest.approx(
        independent["negative_a_s_branching_fraction"]
    )
    assert result.diagnostics["a_s_np_proxy"] == pytest.approx(0.0)
    assert result.diagnostics["lambda_y7v_np_proxy"] == pytest.approx(0.0j)
    assert result.diagnostics["wilson_prefactor_reused"] is False
    assert result.diagnostics["second_mkk_suppression_applied"] is False
    assert result.diagnostics["a_s_sign_envelope_used_for_hard_verdict"] is True
    assert result.diagnostics["selected_a_s_branch"] == "positive"
    assert result.ratio == pytest.approx(0.0)
    assert result.passes is True


def test_evaluate_runs_end_to_end_with_real_finite_fields_and_complex_diagnostics():
    point = sample_rare_kaon_point()
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
        "lambda_y7v_np_proxy",
        "lambda_y7a_np_proxy",
        "a_s_np_proxy_complex",
    ):
        assert isinstance(result.diagnostics[key], complex)
    for key in (
        "m_kk_gev",
        "matching_scale_gev",
        "rate_coefficient_per_a_s_squared",
        "a_s_abs",
        "a_s_effective",
        "positive_a_s_effective",
        "negative_a_s_effective",
        "positive_a_s_branching_fraction",
        "negative_a_s_branching_fraction",
        "positive_a_s_ratio",
        "negative_a_s_ratio",
        "a_s_np_proxy",
        "np_shift_branching_fraction",
        "budget_lower_edge",
        "budget_upper_edge",
    ):
        assert isinstance(result.diagnostics[key], float)
        assert math.isfinite(result.diagnostics[key])
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]
    assert result.diagnostics["ks_a_s_driven_full_rate"] is True
    assert result.diagnostics["cp_conserving_ks_mode"] is True
    assert result.diagnostics["direct_cp_kl_formula_used"] is False
    assert result.diagnostics["uses_k008_y7v_y7a_rs_proxy"] is False
    assert result.diagnostics["uses_k008_y7v_y7a_core_inputs"] is True
    assert result.diagnostics["uses_vector_y7v_proxy_only_for_a_s"] is False
    assert result.diagnostics["uses_vector_y7v_only_for_a_s"] is True
    assert result.diagnostics["a_s_sign_ambiguity_evaluated"] is True
    assert result.diagnostics["a_s_sign_envelope_used_for_hard_verdict"] is True
    assert result.diagnostics["selected_a_s_branch"] in ("positive", "negative")
    assert result.diagnostics["semileptonic_qcd_running_applied"] is False
    assert result.diagnostics["semileptonic_qcd_running_effect_fraction"] == 0.0
    assert result.diagnostics["rs_semileptonic_wilsons_present"] is True
    assert result.diagnostics["rs_semileptonic_matching_status"].endswith(
        "no_second_1_over_M_KK_squared"
    )
    assert result.diagnostics["rare_kaon_proxy_reused"] is False


@pytest.mark.parametrize(
    ("scale", "expected_pass"),
    [
        (0.1, True),
        (1.0e5, False),
    ],
)
def test_safe_point_passes_and_large_np_point_fails(
    scale: float,
    expected_pass: bool,
):
    constraint = fcc.get(_PID)
    point = scaled_rare_kaon_point(scale=scale)
    result = constraint.evaluate(point)
    wilsons = core_y7_wilsons_from_rs_coeff(
        rs_coeff(point, lepton="e"),
        inputs=constraint.sm_inputs,
        matching_scale_gev=point.extras["kk_ew_mass_gev"],
    )
    independent = _independent_ks_pi0ee_full_rate(
        constraint.chpt_inputs,
        lambda_y7v_np_proxy=wilsons.lambda_y7v_np_proxy,
    )
    selected = _independent_conservative_branch(independent, constraint)

    assert result.passes is expected_pass
    assert result.predicted == pytest.approx(selected["predicted"])
    assert result.sm_prediction == pytest.approx(
        independent["sm_branching_fraction"]
    )
    assert result.ratio == pytest.approx(selected["ratio"])
    assert result.budget == pytest.approx(selected["budget"])
    assert result.diagnostics["positive_a_s_branching_fraction"] == pytest.approx(
        independent["positive_a_s_branching_fraction"]
    )
    assert result.diagnostics["negative_a_s_branching_fraction"] == pytest.approx(
        independent["negative_a_s_branching_fraction"]
    )
    assert result.diagnostics["selected_a_s_branch"] == selected["branch"]
    if expected_pass:
        assert result.ratio < 1.0
    else:
        assert result.ratio > 10.0


def test_optional_kk_ew_mass_extra_is_diagnostic_only_for_rigorous_wilsons():
    default_point = scaled_rare_kaon_point(scale=20.0)
    ew_point = scaled_rare_kaon_point(scale=20.0, kk_ew_mass_gev=6000.0)
    default_result = fcc.get(_PID).evaluate(default_point)
    ew_result = fcc.get(_PID).evaluate(ew_point)

    assert default_result.diagnostics["m_kk_gev"] == pytest.approx(3000.0)
    assert ew_result.diagnostics["m_kk_gev"] == pytest.approx(6000.0)
    assert ew_result.diagnostics["kk_ew_mass_extra_used"] is True
    assert ew_result.predicted == pytest.approx(default_result.predicted)
    assert ew_result.diagnostics["lambda_y7v_np_proxy"] == pytest.approx(
        default_result.diagnostics["lambda_y7v_np_proxy"]
    )
    assert ew_result.diagnostics["a_s_np_proxy"] == pytest.approx(
        default_result.diagnostics["a_s_np_proxy"]
    )
    assert ew_result.diagnostics["second_mkk_suppression_applied"] is False


def test_evaluate_is_pure_and_deterministic():
    point = sample_rare_kaon_point()
    before_bundle = point.extras["rs_semileptonic_wilsons"]
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    assert point.extras["rs_semileptonic_wilsons"] is before_bundle
