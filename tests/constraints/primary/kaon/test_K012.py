"""Production tests for K012 (K_S -> mu+ mu- short-distance proxy)."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest
import yaml

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import point_builder
from flavor_catalog_constraints.anchors import Anchor, AnchorError
from flavor_catalog_constraints.base import ConstraintProtocol, Severity
from flavor_catalog_constraints.physics_adapters.rare_kaon_dilepton_kshort_mumu import (
    kshort_mumu_lifetime_inputs_default,
)
from flavor_catalog_constraints.primary.kaon import K012 as k012_module
from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.model import RotationParameters, ckm_like_unitary
from quarkConstraints.rare_kaon_dilepton import (
    evaluate_klong_mumu_short_distance,
)
from tests.rare_kaon_phase3d_helpers import (
    core_y_wilsons_from_rs_coeff,
    rare_kaon_point_with_muon_y_np_total,
    rs_coeff,
    sample_rare_kaon_point,
    scaled_rare_kaon_point,
    sm_limit_rare_kaon_point,
)

_PID = "K012"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "kaon" / "K012.yaml"


def _yaml():
    with open(_SIDECAR) as handle:
        return yaml.safe_load(handle)


def _entry(value_id: str):
    for entry in _yaml()["pdg_or_equivalent"]:
        if entry["value_id"] == value_id:
            return entry
    raise AssertionError(f"no K012 value_id {value_id!r}")


def _limit(value: str) -> float:
    return float(value.split()[0].replace("<", ""))


def _approximate(value: str) -> float:
    return float(value.replace("approximately", "").strip())


def _percent_limit(value: str) -> float:
    return float(value.replace("<", "").replace("%", "").strip()) / 100.0


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


def _independent_kshort_mumu_sd(inputs, lifetime_inputs, y_np_total=0.0j) -> float:
    matrix = ckm_like_unitary(
        RotationParameters(
            theta12=inputs.theta12,
            theta13=inputs.theta13,
            theta23=inputs.theta23,
            delta=inputs.delta,
        )
    )
    lam = float(abs(matrix[0, 1]))
    lambda_c = complex(np.conjugate(matrix[1, 1]) * matrix[1, 0])
    lambda_t = complex(np.conjugate(matrix[2, 1]) * matrix[2, 0])
    kappa_mu = inputs.kappa_mu_ref * (lam / inputs.kappa_lambda_ref) ** 8
    amplitude = (
        (lambda_t * inputs.y_t + complex(y_np_total)).imag / lam**5
        - lambda_c.imag / lam * inputs.p_c_y
    )
    return float(lifetime_inputs.lifetime_ratio * kappa_mu * amplitude * amplitude)


def _independent_sm_kshort_mumu_sd(inputs, lifetime_inputs) -> float:
    return _independent_kshort_mumu_sd(inputs, lifetime_inputs)


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.HARD
    assert constraint.family == "kaon"
    assert constraint.observable == "BR(K_S -> mu+ mu-)_SD"


def test_anchor_matches_yaml_and_budget_band():
    constraint = fcc.get(_PID)
    headline = _entry("PDG2025:K012:headline_limit")
    standalone = _entry("LHCb2020:K012:standalone_2016_2018_limit")
    combined = _entry("LHCb2020:K012:combined_limit")
    run1 = _entry("LHCb2017:K012:run1_limit")
    sm_estimate = _entry("Chobanova2018:K012:sm_estimate")
    hadronic = _entry(
        "DeryGhoshGrossmanSchacht2021:K012:ell0_hadronic_uncertainty"
    )

    assert constraint.anchor.headline_limit.value == pytest.approx(
        _limit(headline["value"])
    )
    assert constraint.anchor.headline_limit.is_upper_limit is True
    assert constraint.anchor.headline_limit.source_url == headline["source_url"]
    assert constraint.anchor.lhcb_standalone_limit.value == pytest.approx(
        _limit(standalone["value"])
    )
    assert constraint.anchor.lhcb_combined_limit.value == pytest.approx(
        _limit(combined["value"])
    )
    assert constraint.anchor.run1_limit.value == pytest.approx(_limit(run1["value"]))
    assert constraint.anchor.standard_model_total_estimate.value == pytest.approx(
        _approximate(sm_estimate["value"])
    )
    assert constraint.anchor.standard_model_total_estimate.is_approximate is True
    assert constraint.anchor.hadronic_uncertainty_target.value == pytest.approx(
        _percent_limit(hadronic["value"])
    )
    assert constraint.anchor.budget == pytest.approx(_limit(headline["value"]))
    assert constraint.anchor.budget_band.sm_sd_formula_value == pytest.approx(
        constraint.sm_result.branching_fraction
    )
    assert constraint.anchor.budget_band.sm_total_estimate_value == pytest.approx(
        _approximate(sm_estimate["value"])
    )
    assert constraint.anchor.budget_band.sm_total_estimate_is_short_distance is False
    assert constraint.anchor.budget_band.sm_subtracted is False
    assert constraint.anchor.budget_band.long_distance_dominated is True


def test_k012_value_anchors_route_through_scaffold_load_anchor(monkeypatch):
    calls: list[tuple[str, ...]] = []
    original_load_anchor = k012_module.load_anchor

    def spy_load_anchor(*args, **kwargs):
        calls.append(tuple(kwargs["candidates"]))
        return original_load_anchor(*args, **kwargs)

    monkeypatch.setattr(k012_module, "load_anchor", spy_load_anchor)
    anchor = k012_module._load_k012_anchor(
        _PID,
        sm_formula_value=fcc.get(_PID).sm_result.branching_fraction,
    )

    assert calls == [
        ("pdg_or_equivalent[0]",),
        ("pdg_or_equivalent[1]",),
        ("pdg_or_equivalent[2]",),
        ("pdg_or_equivalent[3]",),
        ("pdg_or_equivalent[4]",),
        ("pdg_or_equivalent[5]",),
    ]
    assert anchor.value == pytest.approx(fcc.get(_PID).anchor.value)

    def mismatched_load_anchor(*args, **kwargs):
        return Anchor(
            process_id=_PID,
            block_key="pdg_or_equivalent[99]",
            value=1.0,
            uncertainty=None,
        )

    monkeypatch.setattr(k012_module, "load_anchor", mismatched_load_anchor)
    with pytest.raises(AnchorError, match="load_anchor selected"):
        k012_module._load_value_anchor(
            _PID,
            "PDG2025:K012:headline_limit",
            expect_upper_limit=True,
            expected_units="branching fraction",
        )


def test_k012_anchor_loud_fail_probe():
    with pytest.raises(k012_module.AnchorError):
        k012_module._load_value_anchor(
            _PID,
            "not a K012 value_id",
            expect_upper_limit=True,
        )
    with pytest.raises(k012_module.AnchorError):
        k012_module._load_value_anchor(
            _PID,
            "Chobanova2018:K012:sm_estimate",
            expect_upper_limit=True,
            expected_units="branching fraction",
        )


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
        fcc.get(_PID).sm_result.branching_fraction
    )
    assert result.budget == pytest.approx(fcc.get(_PID).anchor.budget)
    assert result.diagnostics["evaluated"] is False
    assert result.diagnostics["missing_extra"] == "rs_semileptonic_wilsons"
    assert result.diagnostics["long_distance_dominated_total_rate"] is True
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
    independent = _independent_sm_kshort_mumu_sd(
        constraint.sm_inputs,
        constraint.lifetime_inputs,
    )
    direct_klong = evaluate_klong_mumu_short_distance(
        core_y_wilsons_from_rs_coeff(
            rs_coeff(point, lepton="mu"),
            inputs=constraint.sm_inputs,
            matching_scale_gev=point.extras["kk_ew_mass_gev"],
        ),
        inputs=constraint.sm_inputs,
    )

    assert result.predicted == pytest.approx(independent)
    assert result.predicted == pytest.approx(result.sm_prediction)
    assert result.predicted == pytest.approx(1.864194959306828e-13)
    assert result.predicted < (
        direct_klong.branching_fraction * constraint.lifetime_inputs.lifetime_ratio
    )
    assert result.diagnostics["klong_short_distance_branching_fraction"] == (
        pytest.approx(direct_klong.branching_fraction)
    )
    assert result.diagnostics["uses_imaginary_ks_ell0_projection"] is True
    assert result.diagnostics["kshort_ell0_projection_source"] == (
        "Im[-lambda_c Y_c + lambda_t C10]^2"
    )
    assert result.diagnostics["sm_context_total_not_sd_anchor"] is True
    assert result.diagnostics["kshort_to_klong_lifetime_ratio"] == pytest.approx(
        kshort_mumu_lifetime_inputs_default().lifetime_ratio
    )
    assert result.diagnostics["sm_context_total_branching_fraction"] == pytest.approx(
        constraint.anchor.standard_model_total_estimate.value
    )
    assert result.diagnostics["sm_context_is_short_distance"] is False
    assert result.diagnostics["y_np_total"] == pytest.approx(0.0j)
    assert result.diagnostics["wilson_prefactor_reused"] is False
    assert result.diagnostics["second_mkk_suppression_applied"] is False
    assert result.predicted < constraint.anchor.value
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
        "lambda_c",
        "lambda_t",
        "y_eff",
        "y_np_total",
    ):
        assert isinstance(result.diagnostics[key], complex)
    for key in (
        "m_kk_gev",
        "matching_scale_gev",
        "kappa_mu",
        "p_c_y",
        "y_t",
        "short_distance_amplitude",
        "np_shift_branching_fraction",
        "klong_short_distance_branching_fraction",
        "kshort_to_klong_lifetime_ratio",
    ):
        assert isinstance(result.diagnostics[key], float)
        assert math.isfinite(result.diagnostics[key])
    assert result.diagnostics["uses_k006_muon_sd_evaluator"] is False
    assert result.diagnostics["uses_k006_muon_sd_wilson_proxy"] is False
    assert result.diagnostics["uses_k006_muon_sd_y_core_inputs"] is True
    assert result.diagnostics["uses_lifetime_rescaled_short_distance_proxy"] is False
    assert result.diagnostics["uses_imaginary_ks_ell0_projection"] is True
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]
    assert result.diagnostics["rs_semileptonic_wilsons_present"] is True
    assert result.diagnostics["rs_semileptonic_matching_status"].endswith(
        "no_second_1_over_M_KK_squared"
    )
    assert result.diagnostics["rare_kaon_proxy_reused"] is False


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
    point = scaled_rare_kaon_point(scale=scale)
    result = constraint.evaluate(point)
    wilsons = core_y_wilsons_from_rs_coeff(
        rs_coeff(point, lepton="mu"),
        inputs=constraint.sm_inputs,
        matching_scale_gev=point.extras["kk_ew_mass_gev"],
    )
    direct_klong = evaluate_klong_mumu_short_distance(
        wilsons,
        inputs=constraint.sm_inputs,
    )
    expected = _independent_kshort_mumu_sd(
        constraint.sm_inputs,
        constraint.lifetime_inputs,
        y_np_total=wilsons.y_np_total,
    )

    assert result.passes is expected_pass
    assert result.predicted == pytest.approx(expected)
    assert result.sm_prediction == pytest.approx(
        constraint.sm_result.branching_fraction
    )
    assert result.diagnostics["klong_short_distance_branching_fraction"] == (
        pytest.approx(direct_klong.branching_fraction)
    )
    assert result.ratio == pytest.approx(expected / constraint.anchor.budget)
    if expected_pass:
        assert result.ratio < 1.0
    else:
        assert result.ratio > 10.0


def test_real_only_np_coupling_leaves_ks_imaginary_sd_near_sm():
    constraint = fcc.get(_PID)
    point = rare_kaon_point_with_muon_y_np_total(
        2.0e-3,
        inputs=constraint.sm_inputs,
    )
    result = constraint.evaluate(point)

    assert result.passes is True
    assert result.predicted == pytest.approx(constraint.sm_result.branching_fraction)
    assert result.ratio < 1.0
    assert result.diagnostics["y_np_total"].imag == pytest.approx(0.0)
    assert result.diagnostics["uses_imaginary_ks_ell0_projection"] is True


def test_optional_kk_ew_mass_extra_is_diagnostic_only_for_rigorous_wilsons():
    default_point = scaled_rare_kaon_point(scale=20.0)
    ew_point = scaled_rare_kaon_point(scale=20.0, kk_ew_mass_gev=6000.0)
    default_result = fcc.get(_PID).evaluate(default_point)
    ew_result = fcc.get(_PID).evaluate(ew_point)

    assert default_result.diagnostics["m_kk_gev"] == pytest.approx(3000.0)
    assert ew_result.diagnostics["m_kk_gev"] == pytest.approx(6000.0)
    assert ew_result.diagnostics["kk_ew_mass_extra_used"] is True
    assert ew_result.predicted == pytest.approx(default_result.predicted)
    assert ew_result.diagnostics["y_np_total"] == pytest.approx(
        default_result.diagnostics["y_np_total"]
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
