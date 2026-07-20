"""Production tests for B033 (B0 -> phi K_S penguin CP asymmetry stub)."""

from __future__ import annotations

import math
from dataclasses import replace
from pathlib import Path

import pytest
import yaml

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import point_builder
from flavor_catalog_constraints.base import ConstraintProtocol, Severity
from flavor_catalog_constraints.physics_adapters.charmless_b_phiks import (
    compare_sphiks_to_sin2beta_reference,
)
from flavor_catalog_constraints.primary.beauty import B033 as b033_module

_PID = "B033"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "beauty" / "B033.yaml"

_MEASUREMENT_BLOCK = "canonical_hflav_phiK0"
_REFERENCE_BLOCK = "hflav_comparison_sin2beta"
_S_PHI_KS = "sin(2 beta_eff)(phi K0)"
_C_PHI_K0 = "C_CP(phi K0)"
_SIN2BETA = "sin(2 beta), all charmonium"
_SIN2BETA_JPSI = "sin(2 beta), J/psi K_S mode"
_DELTA_S = "difference sin(2 beta_eff)(phi K0) - sin(2 beta all charmonium)"


def _yaml_pdg_block():
    with open(_SIDECAR) as handle:
        return yaml.safe_load(handle)["pdg_or_equivalent"]


def _observable(block_key: str, name: str):
    entries = _yaml_pdg_block()[block_key]["observables"]
    return next(entry for entry in entries if entry["name"] == name)


def _expected_delta_s() -> tuple[float, float, float]:
    s_phi = _observable(_MEASUREMENT_BLOCK, _S_PHI_KS)
    sin2beta = _observable(_REFERENCE_BLOCK, _SIN2BETA)
    delta = float(s_phi["value"]) - float(sin2beta["value"])
    sigma = math.sqrt(
        float(s_phi["uncertainty"]) ** 2 + float(sin2beta["uncertainty"]) ** 2
    )
    ratio = abs(delta) / sigma
    return float(delta), float(sigma), float(ratio)


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.INFO
    assert constraint.family == "beauty"
    assert constraint.observable == "S_phiK_S"


def test_anchor_matches_yaml_and_reference_budget():
    constraint = fcc.get(_PID)
    pdg = _yaml_pdg_block()
    s_phi = _observable(_MEASUREMENT_BLOCK, _S_PHI_KS)
    c_phi = _observable(_MEASUREMENT_BLOCK, _C_PHI_K0)
    sin2beta = _observable(_REFERENCE_BLOCK, _SIN2BETA)
    sin2beta_jpsi = _observable(_REFERENCE_BLOCK, _SIN2BETA_JPSI)
    delta_yaml = _observable(_REFERENCE_BLOCK, _DELTA_S)
    delta, sigma, ratio = _expected_delta_s()

    assert constraint.anchor.s_phi_ks.value == pytest.approx(s_phi["value"])
    assert constraint.anchor.s_phi_ks.uncertainty == pytest.approx(
        s_phi["uncertainty"]
    )
    assert constraint.anchor.s_phi_ks.source_url == pdg[_MEASUREMENT_BLOCK][
        "source_url"
    ]
    assert constraint.anchor.s_phi_ks.convention == pdg[_MEASUREMENT_BLOCK][
        "convention"
    ]
    assert constraint.anchor.c_phi_k0.value == pytest.approx(c_phi["value"])
    assert constraint.anchor.c_phi_k0.uncertainty == pytest.approx(
        c_phi["uncertainty"]
    )
    assert constraint.anchor.sin2beta_all_charmonium.value == pytest.approx(
        sin2beta["value"]
    )
    assert constraint.anchor.sin2beta_all_charmonium.uncertainty == pytest.approx(
        sin2beta["uncertainty"]
    )
    assert constraint.anchor.sin2beta_jpsi_ks.value == pytest.approx(
        sin2beta_jpsi["value"]
    )
    assert constraint.anchor.sin2beta_jpsi_ks.uncertainty == pytest.approx(
        sin2beta_jpsi["uncertainty"]
    )
    assert constraint.anchor.delta_s_yaml.value == pytest.approx(delta_yaml["value"])
    assert constraint.anchor.value == pytest.approx(s_phi["value"])
    assert constraint.anchor.sm_value == pytest.approx(sin2beta["value"])
    assert constraint.anchor.comparison.delta_s == pytest.approx(delta)
    assert constraint.anchor.budget == pytest.approx(sigma)
    assert constraint.anchor.comparison.ratio_to_reference_uncertainty == pytest.approx(
        ratio
    )
    assert constraint.anchor.budget == pytest.approx(0.12050311199301038)


def test_anchor_loading_routes_list_entries_through_load_anchor(monkeypatch):
    calls: list[tuple[str, ...]] = []
    original_load_anchor = b033_module.load_anchor

    def spy_load_anchor(*args, **kwargs):
        calls.append(tuple(kwargs["candidates"]))
        return original_load_anchor(*args, **kwargs)

    monkeypatch.setattr(b033_module, "load_anchor", spy_load_anchor)
    anchor = b033_module._load_b033_anchor(_PID)

    assert calls == [
        ("canonical_hflav_phiK0.observables[0]",),
        ("canonical_hflav_phiK0.observables[1]",),
        ("hflav_comparison_sin2beta.observables[0]",),
        ("hflav_comparison_sin2beta.observables[1]",),
        ("hflav_comparison_sin2beta.observables[2]",),
    ]
    assert anchor.value == pytest.approx(fcc.get(_PID).anchor.value)

    with pytest.raises(fcc.AnchorError, match="has no 'value' field"):
        fcc.load_anchor(
            _PID,
            family="beauty",
            candidates=("canonical_hflav_phiK0",),
        )

    with pytest.raises(fcc.AnchorError, match="has no observable"):
        b033_module._load_hflav_observable(
            _yaml_pdg_block(),
            _MEASUREMENT_BLOCK,
            "missing B033 observable",
            process_id=_PID,
        )


def test_anchor_loading_rejects_mismatched_load_anchor_block_key(monkeypatch):
    original_load_anchor = b033_module.load_anchor

    def mismatched_load_anchor(*args, **kwargs):
        anchor = original_load_anchor(*args, **kwargs)
        return replace(anchor, block_key="wrong_block")

    monkeypatch.setattr(b033_module, "load_anchor", mismatched_load_anchor)
    with pytest.raises(fcc.AnchorError, match="load_anchor selected 'wrong_block'"):
        b033_module._load_b033_anchor(_PID)


def test_stub_numeric_validation_matches_independent_yaml_recompute():
    constraint = fcc.get(_PID)
    result = constraint.evaluate(point_builder.empty_point())
    s_phi = _observable(_MEASUREMENT_BLOCK, _S_PHI_KS)
    sin2beta = _observable(_REFERENCE_BLOCK, _SIN2BETA)
    delta, sigma, ratio = _expected_delta_s()

    assert result.predicted is None
    assert result.experimental == pytest.approx(s_phi["value"])
    assert result.sm_prediction == pytest.approx(sin2beta["value"])
    assert result.budget == pytest.approx(sigma)
    assert result.ratio == pytest.approx(ratio)
    assert result.ratio == pytest.approx(0.24895622877980005)
    assert result.diagnostics["delta_s_value"] == pytest.approx(delta)
    assert result.diagnostics["delta_s_uncertainty"] == pytest.approx(sigma)
    assert result.diagnostics["delta_s_reference_ratio"] == pytest.approx(ratio)
    assert result.diagnostics["no_hadronic_penguin_calculation"] is True
    assert result.diagnostics["no_delta_b1_penguin_matching"] is True
    assert result.diagnostics["b002_clean_mixing_phase_reused"] is False


def test_adapter_reference_comparison_has_pass_fail_behavior_without_amplitudes():
    delta, sigma, ratio = _expected_delta_s()
    safe = compare_sphiks_to_sin2beta_reference(
        measured_s_phi_ks=0.74,
        measured_uncertainty=0.12,
        sin2beta_reference=0.710,
        sin2beta_uncertainty=0.011,
    )
    excluded_advisory = compare_sphiks_to_sin2beta_reference(
        measured_s_phi_ks=1.0,
        measured_uncertainty=0.01,
        sin2beta_reference=0.0,
        sin2beta_uncertainty=0.01,
    )

    assert safe.passes is True
    assert safe.delta_s == pytest.approx(delta)
    assert safe.delta_s_uncertainty == pytest.approx(sigma)
    assert safe.ratio_to_reference_uncertainty == pytest.approx(ratio)
    assert excluded_advisory.passes is False
    assert excluded_advisory.ratio_to_reference_uncertainty > 10.0


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
    assert empty_result.passes is True
    assert empty_result.predicted is None
    assert empty_result.experimental == pytest.approx(constraint.anchor.value)
    assert empty_result.sm_prediction == pytest.approx(constraint.anchor.sm_value)
    assert empty_result.budget == pytest.approx(constraint.anchor.budget)
    assert empty_result.diagnostics["non_vetoing"] is True
    assert empty_result.diagnostics["parameter_point_inputs_used"] == ()
    assert empty_result.diagnostics["full_sm_sphi_ks_prediction_available"] is False
    assert empty_result.diagnostics["rs_delta_b1_penguin_amplitude_available"] is False
    assert "NEEDS-HUMAN-PHYSICS" in empty_result.diagnostics["needs_human_physics_sm"]
    assert "NEEDS-HUMAN-PHYSICS" in empty_result.diagnostics["needs_human_physics_np"]


def test_evaluate_has_real_finite_numeric_fields():
    result = fcc.get(_PID).evaluate(point_builder.empty_point())

    assert result.predicted is None
    for value in (result.sm_prediction, result.experimental, result.ratio, result.budget):
        assert isinstance(value, float)
        assert math.isfinite(value)
    for key in (
        "s_phi_ks_value",
        "s_phi_ks_uncertainty",
        "sin2beta_reference_value",
        "sin2beta_reference_uncertainty",
        "delta_s_value",
        "delta_s_yaml_central",
        "delta_s_uncertainty",
        "delta_s_reference_ratio",
        "c_phi_k0_value",
        "c_phi_k0_uncertainty",
        "sin2beta_jpsi_ks_value",
        "sin2beta_jpsi_ks_uncertainty",
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
