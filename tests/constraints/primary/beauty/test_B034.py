"""Production tests for B034 (B_s -> phi phi penguin CP phase stub)."""

from __future__ import annotations

from dataclasses import replace
import math
from pathlib import Path

import pytest
import yaml

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import point_builder
from flavor_catalog_constraints.base import ConstraintProtocol, Severity
from flavor_catalog_constraints.physics_adapters.charmless_b_phiphi import (
    compare_bs_phiphi_phase_to_room,
)
from flavor_catalog_constraints.primary.beauty import B034 as b034_module

_PID = "B034"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "beauty" / "B034.yaml"

_PDG_BLOCK = "pdg_beta_s_b_to_ssbars"
_LHCB_BLOCK = "lhcb2023_precision_and_combination"
_BETA_S = "beta_s(b -> s sbar s)"
_PHI_S = "phi_s^{s sbar s}, LHCb combined note"
_LAMBDA = "|lambda|, LHCb combined note"
_LHCB_RUN12_PHI_S = "phi_s^{s sbar s}, LHCb Run 1 + Run 2 combination"
_LHCB_RUN12_LAMBDA = "|lambda|, LHCb Run 1 + Run 2 combination"


def _yaml_pdg_block():
    with open(_SIDECAR) as handle:
        return yaml.safe_load(handle)["pdg_or_equivalent"]


def _observable(block_key: str, name: str):
    entries = _yaml_pdg_block()[block_key]["observables"]
    return next(entry for entry in entries if entry["name"] == name)


def _expected_phi_from_beta() -> tuple[float, float]:
    beta = _observable(_PDG_BLOCK, _BETA_S)
    return -2.0e-2 * float(beta["value"]), 2.0e-2 * float(beta["uncertainty"])


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.INFO
    assert constraint.family == "beauty"
    assert constraint.observable == "phi_s^sss(Bs -> phi phi)"


def test_anchor_matches_yaml_and_non_vetoing_uncertainty_room():
    constraint = fcc.get(_PID)
    pdg = _yaml_pdg_block()
    phi_s = _observable(_PDG_BLOCK, _PHI_S)
    beta_s = _observable(_PDG_BLOCK, _BETA_S)
    lam = _observable(_PDG_BLOCK, _LAMBDA)
    lhcb_phi_s = _observable(_LHCB_BLOCK, _LHCB_RUN12_PHI_S)
    lhcb_lam = _observable(_LHCB_BLOCK, _LHCB_RUN12_LAMBDA)
    phi_from_beta, phi_sigma_from_beta = _expected_phi_from_beta()

    assert constraint.anchor.phi_s_sss.value == pytest.approx(phi_s["value"])
    assert constraint.anchor.phi_s_sss.uncertainty == pytest.approx(
        phi_s["uncertainty"]
    )
    assert constraint.anchor.phi_s_sss.source_url == pdg[_PDG_BLOCK]["source_url"]
    assert constraint.anchor.phi_s_sss.convention == pdg[_PDG_BLOCK]["convention"]
    assert constraint.anchor.beta_s.value == pytest.approx(beta_s["value"])
    assert constraint.anchor.beta_s.uncertainty == pytest.approx(
        beta_s["uncertainty"]
    )
    assert constraint.anchor.lambda_abs.value == pytest.approx(lam["value"])
    assert constraint.anchor.lambda_abs.uncertainty == pytest.approx(
        lam["uncertainty"]
    )
    assert constraint.anchor.lhcb_run12_phi_s_sss.value == pytest.approx(
        lhcb_phi_s["value"]
    )
    assert constraint.anchor.lhcb_run12_phi_s_sss.uncertainty == pytest.approx(
        lhcb_phi_s["uncertainty"]
    )
    assert constraint.anchor.lhcb_run12_lambda_abs.value == pytest.approx(
        lhcb_lam["value"]
    )
    assert constraint.anchor.lhcb_run12_lambda_abs.uncertainty == pytest.approx(
        lhcb_lam["uncertainty"]
    )
    assert constraint.anchor.phi_s_from_beta_convention == pytest.approx(
        phi_from_beta
    )
    assert constraint.anchor.phi_s_uncertainty_from_beta_convention == pytest.approx(
        phi_sigma_from_beta
    )
    assert constraint.anchor.value == pytest.approx(phi_s["value"])
    assert constraint.anchor.budget == pytest.approx(float(phi_s["uncertainty"]))
    assert constraint.anchor.budget == pytest.approx(0.069)


def test_anchor_loading_routes_list_entries_through_load_anchor(monkeypatch):
    calls: list[tuple[str, ...]] = []
    original_load_anchor = b034_module.load_anchor

    def spy_load_anchor(*args, **kwargs):
        calls.append(tuple(kwargs["candidates"]))
        return original_load_anchor(*args, **kwargs)

    monkeypatch.setattr(b034_module, "load_anchor", spy_load_anchor)
    anchor = b034_module._load_b034_anchor(_PID)

    assert calls == [
        ("pdg_beta_s_b_to_ssbars.observables[1]",),
        ("pdg_beta_s_b_to_ssbars.observables[0]",),
        ("pdg_beta_s_b_to_ssbars.observables[2]",),
        ("lhcb2023_precision_and_combination.observables[2]",),
        ("lhcb2023_precision_and_combination.observables[3]",),
    ]
    assert anchor.value == pytest.approx(fcc.get(_PID).anchor.value)

    with pytest.raises(fcc.AnchorError, match="has no 'value' field"):
        fcc.load_anchor(_PID, family="beauty", candidates=(_PDG_BLOCK,))

    with pytest.raises(fcc.AnchorError, match="has no observable"):
        b034_module._load_b034_observable(
            _yaml_pdg_block(),
            _PDG_BLOCK,
            "missing B034 observable",
            process_id=_PID,
            expected_units="rad",
        )


def test_anchor_loading_rejects_mismatched_load_anchor_block_key(monkeypatch):
    original_load_anchor = b034_module.load_anchor

    def mismatched_load_anchor(*args, **kwargs):
        anchor = original_load_anchor(*args, **kwargs)
        return replace(anchor, block_key="wrong_block")

    monkeypatch.setattr(b034_module, "load_anchor", mismatched_load_anchor)
    with pytest.raises(fcc.AnchorError, match="load_anchor selected 'wrong_block'"):
        b034_module._load_b034_anchor(_PID)


def test_stub_numeric_validation_matches_independent_yaml_recompute():
    constraint = fcc.get(_PID)
    result = constraint.evaluate(point_builder.empty_point())
    phi_s = _observable(_PDG_BLOCK, _PHI_S)
    phi_from_beta, phi_sigma_from_beta = _expected_phi_from_beta()
    expected_budget = float(phi_s["uncertainty"])
    expected_ratio = abs(float(phi_s["value"])) / expected_budget

    assert result.predicted is None
    assert result.sm_prediction is None
    assert result.experimental == pytest.approx(phi_s["value"])
    assert result.budget == pytest.approx(expected_budget)
    assert result.ratio == pytest.approx(expected_ratio)
    assert result.ratio == pytest.approx(0.074 / 0.069)
    assert result.diagnostics["phi_s_sss_value"] == pytest.approx(phi_s["value"])
    assert result.diagnostics["phi_s_sss_uncertainty"] == pytest.approx(
        phi_s["uncertainty"]
    )
    assert result.diagnostics["phi_s_from_beta_convention_rad"] == pytest.approx(
        phi_from_beta
    )
    assert result.diagnostics[
        "phi_s_uncertainty_from_beta_convention_rad"
    ] == pytest.approx(phi_sigma_from_beta)
    assert result.diagnostics["no_hadronic_penguin_calculation"] is True
    assert result.diagnostics["no_delta_b1_penguin_matching"] is True
    assert result.diagnostics["bs_mixing_phase_reused"] is False
    assert result.diagnostics["budget_is_measurement_uncertainty"] is True


def test_adapter_room_comparison_has_pass_fail_behavior_without_amplitudes():
    safe = compare_bs_phiphi_phase_to_room(
        measured_phi_s=-0.02,
        experimental_uncertainty=0.07,
        documented_np_room_abs=0.05,
    )
    excluded_advisory = compare_bs_phiphi_phase_to_room(
        measured_phi_s=0.20,
        experimental_uncertainty=0.07,
        documented_np_room_abs=0.05,
    )

    assert safe.passes is True
    assert safe.ratio_to_room == pytest.approx(0.4)
    assert excluded_advisory.passes is False
    assert excluded_advisory.ratio_to_room == pytest.approx(4.0)


def test_evaluate_is_info_non_vetoing_and_uses_no_point_inputs():
    constraint = fcc.get(_PID)
    empty_result = constraint.evaluate(point_builder.empty_point())
    irrelevant_point = point_builder.make_point(
        kk_gluon_mass_gev=100.0,
        kk_ew_mass_gev=200.0,
    )
    irrelevant_result = constraint.evaluate(irrelevant_point)

    assert empty_result == irrelevant_result
    assert empty_result.process_id == _PID
    assert empty_result.severity is Severity.INFO
    assert empty_result.passes is False
    assert empty_result.predicted is None
    assert empty_result.sm_prediction is None
    assert empty_result.experimental == pytest.approx(constraint.anchor.value)
    assert empty_result.budget == pytest.approx(constraint.anchor.budget)
    assert empty_result.diagnostics["non_vetoing"] is True
    assert empty_result.diagnostics["severity_policy"] == (
        "INFO/non-vetoing even when advisory passes=False"
    )
    assert empty_result.diagnostics["parameter_point_inputs_used"] == ()
    assert empty_result.diagnostics["full_sm_phi_s_sss_prediction_available"] is False
    assert empty_result.diagnostics["rs_delta_b1_penguin_amplitude_available"] is False
    assert "NEEDS-HUMAN-PHYSICS" in empty_result.diagnostics["needs_human_physics_sm"]
    assert "NEEDS-HUMAN-PHYSICS" in empty_result.diagnostics["needs_human_physics_np"]


def test_evaluate_has_real_finite_numeric_fields():
    result = fcc.get(_PID).evaluate(point_builder.empty_point())

    assert result.predicted is None
    assert result.sm_prediction is None
    for value in (result.experimental, result.ratio, result.budget):
        assert isinstance(value, float)
        assert math.isfinite(value)
    for key in (
        "phi_s_sss_value",
        "phi_s_sss_uncertainty",
        "phi_s_sss_abs",
        "phi_s_room_abs",
        "beta_s_value_10minus2_rad",
        "beta_s_uncertainty_10minus2_rad",
        "phi_s_from_beta_convention_rad",
        "phi_s_uncertainty_from_beta_convention_rad",
        "lambda_abs_value",
        "lambda_abs_uncertainty",
        "lhcb_run12_phi_s_value",
        "lhcb_run12_phi_s_uncertainty",
        "lhcb_run12_lambda_abs_value",
        "lhcb_run12_lambda_abs_uncertainty",
    ):
        assert isinstance(result.diagnostics[key], float)
        assert math.isfinite(result.diagnostics[key])


def test_evaluate_is_pure_and_deterministic():
    point = point_builder.make_point(kk_gluon_mass_gev=3000.0)
    before_extras = dict(point.extras)
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    assert dict(point.extras) == before_extras
