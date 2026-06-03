"""Production tests for K004 (K+ -> pi+ nu nubar)."""

from __future__ import annotations

import math
from pathlib import Path
import re

import pytest
import yaml

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import point_builder
from flavor_catalog_constraints.base import ConstraintProtocol, Severity
from flavor_catalog_constraints.physics_adapters.rare_kaon import (
    kplus_piplus_nunu_from_rs_semileptonic_wilsons,
)
from tests.constraints.primary.nunu_phase4d_helpers import (
    direct_contact_x_np,
    nunu_block,
    rigorous_point,
    scalar_x_np,
    sm_limit_point,
)

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


def test_absent_rs_semileptonic_wilsons_degrades_gracefully():
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
    assert result.diagnostics["evaluated"] is False
    assert result.diagnostics["missing_extra"] == "rs_semileptonic_wilsons"
    assert "needs_human_physics" not in result.diagnostics


def test_legacy_quark_proxy_only_is_unevaluated_not_real_pass():
    result = fcc.get(_PID).evaluate(
        point_builder.make_point(quark_mass_basis_couplings=object())
    )

    assert result.passes is True
    assert result.predicted is None
    assert result.ratio is None
    assert result.diagnostics["evaluated"] is False
    assert result.diagnostics["missing_extra"] == "rs_semileptonic_wilsons"
    assert result.diagnostics["legacy_quark_mass_basis_couplings_present"] is True


def test_invalid_rs_semileptonic_wilsons_is_unevaluated_not_real_pass():
    result = fcc.get(_PID).evaluate(
        point_builder.make_point(rs_semileptonic_wilsons=object())
    )

    assert result.passes is True
    assert result.predicted is None
    assert result.ratio is None
    assert result.diagnostics["evaluated"] is False
    assert result.diagnostics["invalid_extra"] == "rs_semileptonic_wilsons"


def test_sm_limit_universal_quark_c_recovers_committed_sm_branching_fraction():
    constraint = fcc.get(_PID)
    point = sm_limit_point()
    result = constraint.evaluate(point)
    block = nunu_block(point, "s_to_d")
    x_left, x_right, x_total = scalar_x_np(block)

    assert x_left == pytest.approx(0.0j, abs=1.0e-18)
    assert x_right == pytest.approx(0.0j, abs=1.0e-18)
    assert x_total == pytest.approx(0.0j, abs=1.0e-18)
    assert result.diagnostics["x_np_total"] == pytest.approx(0.0j, abs=1.0e-18)
    assert result.predicted == pytest.approx(8.472598449611133e-11)
    assert result.predicted == pytest.approx(result.sm_prediction)
    assert result.diagnostics["sm_anchor_branching_fraction"] == pytest.approx(
        constraint.anchor.sm_value
    )
    assert result.diagnostics["sm_formula_minus_anchor"] == pytest.approx(
        result.sm_prediction - constraint.anchor.sm_value
    )
    assert result.diagnostics["delta_em_correction"] == pytest.approx(0.997)
    assert result.diagnostics["p_c"] == pytest.approx(0.404)
    assert result.passes is True


def test_rigorous_nonzero_xnp_shifts_br_and_matches_direct_contact():
    constraint = fcc.get(_PID)
    point = rigorous_point()
    result = constraint.evaluate(point)
    direct = kplus_piplus_nunu_from_rs_semileptonic_wilsons(
        point.extras["rs_semileptonic_wilsons"],
        inputs=constraint.sm_inputs,
    )
    x_left, x_right, x_total = direct_contact_x_np(point, "s_to_d")

    assert abs(x_total) > 1.0e-4
    assert result.diagnostics["x_np_left"] == pytest.approx(x_left)
    assert result.diagnostics["x_np_right"] == pytest.approx(x_right)
    assert result.diagnostics["x_np_total"] == pytest.approx(x_total)
    assert result.predicted == pytest.approx(direct.branching_fraction)
    assert result.predicted != pytest.approx(result.sm_prediction)
    assert result.passes is False
    assert result.ratio > 1.0
    assert result.diagnostics["nunu_mapping"] == "X_NP=C/g_SM^2"
    assert result.diagnostics["wilson_prefactor_reused"] is False
    assert result.diagnostics["second_mkk_suppression_applied"] is False
    assert result.diagnostics["legacy_one_z_proxy_reused"] is False
    assert "needs_human_physics" not in result.diagnostics


def test_evaluate_runs_end_to_end_with_real_finite_fields_and_complex_diagnostics():
    result = fcc.get(_PID).evaluate(rigorous_point())

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
    assert result.diagnostics["evaluated"] is True


def test_kk_ew_mass_extra_is_diagnostic_only_no_second_mkk_suppression():
    base_point = rigorous_point()
    ew_point = point_builder.make_point(
        rs_semileptonic_wilsons=base_point.extras["rs_semileptonic_wilsons"],
        kk_ew_mass_gev=6000.0,
    )
    default_result = fcc.get(_PID).evaluate(base_point)
    ew_result = fcc.get(_PID).evaluate(ew_point)

    assert default_result.diagnostics["m_kk_gev"] == pytest.approx(3000.0)
    assert ew_result.diagnostics["m_kk_gev"] == pytest.approx(6000.0)
    assert ew_result.diagnostics["kk_ew_mass_extra_used"] is True
    assert ew_result.predicted == pytest.approx(default_result.predicted)
    assert ew_result.diagnostics["x_np_total"] == pytest.approx(
        default_result.diagnostics["x_np_total"]
    )


def test_majorana_and_dirac_active_nu_rates_match():
    dirac = fcc.get(_PID).evaluate(rigorous_point(alpha=0.0, beta=0.0))
    majorana = fcc.get(_PID).evaluate(rigorous_point(alpha=1.1, beta=-0.7))

    assert majorana.predicted == pytest.approx(dirac.predicted, rel=0.0, abs=1.0e-22)
    assert majorana.diagnostics["x_np_total"] == pytest.approx(
        dirac.diagnostics["x_np_total"],
        abs=1.0e-18,
    )
    assert majorana.diagnostics["majorana_dirac_rate_factor"] == pytest.approx(1.0)


def test_evaluate_is_pure_and_deterministic():
    point = rigorous_point()
    before_bundle = point.extras["rs_semileptonic_wilsons"]
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    assert point.extras["rs_semileptonic_wilsons"] is before_bundle
