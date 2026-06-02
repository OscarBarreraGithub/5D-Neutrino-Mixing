"""Production tests for CR012 (spin-1 diboson high-mass resonance)."""

from __future__ import annotations

import math
from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import point_builder
from flavor_catalog_constraints.anchors import AnchorError
from flavor_catalog_constraints.base import ConstraintProtocol, Severity
from flavor_catalog_constraints.primary.collider_rs import CR012 as cr012_module
from quarkConstraints import collider_resonance as core

_PID = "CR012"
_ACTIVE_VALUE_ID = "PDG2025:CR012:HVTB_Wprime_WZ_mass_lower"
_VPRIME_VV_VALUE_ID = "CMS2023:CR012:HVTB_Vprime_VV_mass_lower"
_VV_VH_VALUE_ID = "CMS2023:CR012:HVTB_Vprime_VV_VH_mass_lower"
_ZPRIME_INTERVAL_VALUE_ID = "CMS2023:CR012:HVTB_Zprime_WW_excluded_intervals"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "collider_rs" / "CR012.yaml"


def _yaml_pdg_block():
    with open(_SIDECAR) as handle:
        return yaml.safe_load(handle)["pdg_or_equivalent"]


def _entry(value_id: str):
    matches = [
        entry
        for entry in _yaml_pdg_block()["values"]
        if entry.get("value_id") == value_id
    ]
    assert len(matches) == 1
    return matches[0]


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.HARD
    assert constraint.family == "collider_rs"
    assert constraint.observable == "m(V_KK^(1) spin-1 -> WW/WZ/ZZ)"


def test_anchor_matches_yaml_and_loud_fail_probe():
    constraint = fcc.get(_PID)
    active = _entry(_ACTIVE_VALUE_ID)
    vprime_vv = _entry(_VPRIME_VV_VALUE_ID)
    vv_vh = _entry(_VV_VH_VALUE_ID)
    zprime = _entry(_ZPRIME_INTERVAL_VALUE_ID)

    assert constraint.anchor.active_limit.value_tev == pytest.approx(active["value"])
    assert constraint.anchor.active_limit.value_gev == pytest.approx(
        active["value"] * 1000.0
    )
    assert constraint.anchor.active_limit.units == active["units"]
    assert constraint.anchor.active_limit.source_url == active["source_url"]
    assert constraint.anchor.active_limit.display == active["display"]
    assert constraint.anchor.active_limit.experiment == active["experiment"]
    assert constraint.anchor.active_limit.model == active["model"]
    assert constraint.anchor.active_limit.decay_mode == active["decay_mode"]
    assert constraint.anchor.active_limit.confidence_level == pytest.approx(
        active["confidence_level"]
    )
    assert constraint.anchor.budget == pytest.approx(active["value"])
    assert constraint.anchor.mass_degenerate_vv.value_tev == pytest.approx(
        vprime_vv["value"]
    )
    assert constraint.anchor.mass_degenerate_vv.decay_mode == vprime_vv["decay_mode"]
    assert constraint.anchor.vv_vh_combined.value_tev == pytest.approx(vv_vh["value"])
    assert constraint.anchor.vv_vh_combined.decay_mode == vv_vh["decay_mode"]
    assert constraint.anchor.strongest_pure_diboson_limit.value_id == (
        _VPRIME_VV_VALUE_ID
    )
    assert constraint.anchor.zprime_ww_intervals.excluded_ranges_tev == tuple(
        tuple(float(value) for value in interval)
        for interval in zprime["excluded_ranges_TeV"]
    )
    assert constraint.anchor.zprime_ww_intervals.expected_limit_tev == pytest.approx(
        zprime["expected_limit"]
    )
    assert constraint.anchor.zprime_ww_intervals.limit_type == zprime["limit_type"]

    with pytest.raises(AnchorError):
        cr012_module._load_mass_limit_anchor("not-a-real-value-id", process_id=_PID)
    with pytest.raises(AnchorError):
        cr012_module._load_interval_anchor("not-a-real-value-id", process_id=_PID)


def test_evaluate_without_input_degrades_gracefully():
    result = fcc.get(_PID).evaluate(point_builder.empty_point())

    assert result.process_id == _PID
    assert result.passes is True
    assert result.predicted is None
    assert result.ratio is None
    assert result.experimental == pytest.approx(fcc.get(_PID).anchor.value)
    assert result.budget == pytest.approx(fcc.get(_PID).anchor.budget)
    assert result.diagnostics["missing_extras"] == (
        "kk_ew_mass_gev",
        "quark_mass_basis_couplings",
    )
    assert result.diagnostics["active_value_id"] == _ACTIVE_VALUE_ID
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]
    assert result.diagnostics["sigma_times_br_recast_status"] == (
        "NEEDS-HUMAN-PHYSICS"
    )
    assert result.diagnostics["branching_surface_recast_status"] == (
        "NEEDS-HUMAN-PHYSICS"
    )


def test_mass_limit_numerics_validate_against_core_not_adapter():
    constraint = fcc.get(_PID)
    point = point_builder.make_point(kk_ew_mass_gev=5000.0)
    result = constraint.evaluate(point)
    active = _entry(_ACTIVE_VALUE_ID)

    limit = core.ColliderResonanceLimit(
        process_id=_PID,
        resonance="V_KK^(1) spin-1",
        final_state="WW/WZ/ZZ",
        limit_kind=core.MASS_LOWER_BOUND,
        value=float(active["value"]),
        units=active["units"],
    )
    prediction = core.ColliderResonancePrediction(
        resonance="V_KK^(1) spin-1",
        final_state="WW/WZ/ZZ",
        mass_tev=5.0,
    )
    direct = core.evaluate_resonance_limit(prediction, limit)

    assert result.predicted == pytest.approx(5.0)
    assert result.predicted == pytest.approx(direct.predicted_mass_tev)
    assert result.experimental == pytest.approx(float(active["value"]))
    assert result.ratio == pytest.approx(float(active["value"]) / 5.0)
    assert result.ratio == pytest.approx(direct.ratio_to_budget)
    assert result.budget == pytest.approx(direct.budget)
    assert result.passes is True
    assert result.diagnostics["active_value_id"] == _ACTIVE_VALUE_ID
    assert result.diagnostics["active_model"] == active["model"]


def test_mass_falls_back_to_couplings_when_explicit_ew_mass_absent():
    point = point_builder.make_point(
        quark_mass_basis_couplings=SimpleNamespace(M_KK=5200.0)
    )
    result = fcc.get(_PID).evaluate(point)

    assert result.predicted == pytest.approx(5.2)
    assert result.diagnostics["m_kk_gev"] == pytest.approx(5200.0)
    assert result.diagnostics["mass_source"] == "quark_mass_basis_couplings.M_KK"
    assert result.passes is True


def test_safe_point_passes_and_excluded_point_fails():
    constraint = fcc.get(_PID)
    limit_gev = constraint.anchor.budget * 1000.0

    safe = constraint.evaluate(
        point_builder.make_point(kk_ew_mass_gev=limit_gev + 100.0)
    )
    excluded = constraint.evaluate(
        point_builder.make_point(kk_ew_mass_gev=limit_gev - 100.0)
    )

    assert safe.passes is True
    assert safe.predicted == pytest.approx((limit_gev + 100.0) / 1000.0)
    assert safe.ratio < 1.0
    assert excluded.passes is False
    assert excluded.predicted == pytest.approx((limit_gev - 100.0) / 1000.0)
    assert excluded.ratio > 1.0
    assert excluded.experimental == pytest.approx(constraint.anchor.budget)
    assert excluded.diagnostics["mass_proxy"] == (
        "m_spin1_diboson = kk_ew_mass_gev or M_KK"
    )
    assert excluded.diagnostics["mass_source"] == "kk_ew_mass_gev"
    assert "NEEDS-HUMAN-PHYSICS" in excluded.diagnostics["needs_human_physics"]
    assert excluded.diagnostics["sigma_times_br_recast_status"] == (
        "NEEDS-HUMAN-PHYSICS"
    )
    assert excluded.diagnostics["pure_diboson_mass_limits_tev"][
        _VPRIME_VV_VALUE_ID
    ] == pytest.approx(_entry(_VPRIME_VV_VALUE_ID)["value"])
    assert excluded.diagnostics["nonactive_vv_vh_combined_limit_tev"] == pytest.approx(
        _entry(_VV_VH_VALUE_ID)["value"]
    )
    assert excluded.diagnostics["zprime_ww_excluded_ranges_tev"] == tuple(
        tuple(float(value) for value in interval)
        for interval in _entry(_ZPRIME_INTERVAL_VALUE_ID)["excluded_ranges_TeV"]
    )


def test_invalid_ew_mass_is_non_crashing_failure():
    point = point_builder.make_point(kk_ew_mass_gev=-1.0)
    result = fcc.get(_PID).evaluate(point)

    assert result.passes is False
    assert result.predicted is None
    assert result.ratio is None
    assert result.experimental == pytest.approx(fcc.get(_PID).anchor.value)
    assert result.diagnostics["exception_type"] == "ValueError"
    assert result.diagnostics["sigma_times_br_recast_status"] == (
        "NEEDS-HUMAN-PHYSICS"
    )


def test_evaluate_is_pure_and_deterministic():
    point = point_builder.make_point(kk_ew_mass_gev=4600.0)
    before = dict(point.extras)
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    assert point.extras == before
    for value in (first.predicted, first.experimental, first.ratio, first.budget):
        assert isinstance(value, float)
        assert math.isfinite(value)
