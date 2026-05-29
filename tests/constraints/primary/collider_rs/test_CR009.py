"""Production tests for CR009 (Drell-Yan contact-operator scale bound)."""

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
from flavor_catalog_constraints.primary.collider_rs import CR009 as cr009_module
from quarkConstraints import collider_resonance as core

_PID = "CR009"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "collider_rs" / "CR009.yaml"


def _yaml_pdg_block():
    with open(_SIDECAR) as handle:
        return yaml.safe_load(handle)["pdg_or_equivalent"]


def _values():
    return list(_yaml_pdg_block()["values"])


def _entry(value_id: str):
    matches = [entry for entry in _values() if entry.get("value_id") == value_id]
    assert len(matches) == 1
    return matches[0]


def _active_entry():
    return min(_values(), key=lambda entry: float(entry["value"]))


def _strongest_entry():
    return max(_values(), key=lambda entry: float(entry["value"]))


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.HARD
    assert constraint.family == "collider_rs"
    assert constraint.observable == "Lambda(llqq contact operator)"


def test_anchor_matches_yaml_and_loud_fail_probe():
    constraint = fcc.get(_PID)
    active = _active_entry()

    assert constraint.anchor.active_limit.value_id == active["value_id"]
    assert constraint.anchor.active_limit.value_tev == pytest.approx(active["value"])
    assert constraint.anchor.active_limit.units == active["units"]
    assert constraint.anchor.active_limit.source_url == active["source_url"]
    assert constraint.anchor.active_limit.display == active["display"]
    assert constraint.anchor.active_limit.experiment == active["experiment"]
    assert constraint.anchor.active_limit.model_assumptions == (
        active["model_assumptions"]
    )
    assert constraint.anchor.active_limit.value_id == (
        "PDG2025:CR009:CMS_LL_destructive_range_endpoint"
    )
    assert constraint.anchor.active_limit.value_tev == pytest.approx(23.9)
    assert constraint.anchor.budget == pytest.approx(active["value"])
    assert "conservative DY contact policy" in constraint.anchor.budget_policy
    assert "weakest PDG2025" in constraint.anchor.budget_policy
    assert "helicity/interference matching" in constraint.anchor.budget_policy
    assert constraint.anchor.parent_source == _yaml_pdg_block()["canonical_source"]
    assert constraint.anchor.parent_source_url == _yaml_pdg_block()["source_url"]
    assert constraint.anchor.contact_matching_relation == (
        "4*pi/Lambda^2 ~ |g_q g_l|/M_V^2"
    )

    for entry in _values():
        anchor = next(
            limit
            for limit in constraint.anchor.all_limits
            if limit.value_id == entry["value_id"]
        )
        assert anchor.value_tev == pytest.approx(entry["value"])
        assert anchor.model_assumptions == entry["model_assumptions"]

    with pytest.raises(AnchorError):
        cr009_module._load_contact_limit_anchor(
            "not-a-real-value-id",
            process_id=_PID,
        )


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
        "kk_gluon_mass_gev",
        "quark_mass_basis_couplings",
    )
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]
    assert result.diagnostics["eft_recast_status"] == "NEEDS-HUMAN-PHYSICS"
    assert result.diagnostics["helicity_interference_matching_status"] == (
        "NEEDS-HUMAN-PHYSICS"
    )


def test_contact_scale_numerics_validate_against_core_not_adapter():
    constraint = fcc.get(_PID)
    point = point_builder.make_point(kk_ew_mass_gev=40000.0)
    result = constraint.evaluate(point)
    active = _active_entry()

    limit = core.ColliderResonanceLimit(
        process_id=_PID,
        resonance="llqq contact operator",
        final_state="ee + mumu high-mass tail",
        limit_kind=core.MASS_LOWER_BOUND,
        value=float(active["value"]),
        units=active["units"],
    )
    prediction = core.ColliderResonancePrediction(
        resonance="llqq contact operator",
        final_state="ee + mumu high-mass tail",
        mass_tev=40.0,
    )
    direct = core.evaluate_resonance_limit(prediction, limit)

    assert result.predicted == pytest.approx(40.0)
    assert result.predicted == pytest.approx(direct.predicted_mass_tev)
    assert result.experimental == pytest.approx(float(active["value"]))
    assert result.ratio == pytest.approx(float(active["value"]) / 40.0)
    assert result.ratio == pytest.approx(direct.ratio_to_budget)
    assert result.budget == pytest.approx(direct.budget)
    assert result.passes is True
    assert result.diagnostics["active_value_id"] == active["value_id"]
    assert result.diagnostics["active_value_id"] == (
        "PDG2025:CR009:CMS_LL_destructive_range_endpoint"
    )
    assert result.diagnostics["helicity_interference_matching_status"] == (
        "NEEDS-HUMAN-PHYSICS"
    )
    strongest = _strongest_entry()
    assert result.diagnostics["stronger_contact_limits_tev"][strongest["value_id"]] == (
        pytest.approx(strongest["value"])
    )
    assert result.diagnostics["contact_limits_by_benchmark"][strongest["value_id"]] == {
        "value_tev": pytest.approx(strongest["value"]),
        "experiment": strongest["experiment"],
        "helicity": strongest["model_assumptions"]["helicity"],
        "interference": strongest["model_assumptions"]["interference"],
    }


def test_contact_scale_falls_back_to_gluon_mass_then_couplings():
    from_gluon = fcc.get(_PID).evaluate(
        point_builder.make_point(kk_gluon_mass_gev=41000.0)
    )
    from_couplings = fcc.get(_PID).evaluate(
        point_builder.make_point(
            quark_mass_basis_couplings=SimpleNamespace(M_KK=42000.0)
        )
    )

    assert from_gluon.predicted == pytest.approx(41.0)
    assert from_gluon.diagnostics["scale_source"] == "kk_gluon_mass_gev"
    assert from_gluon.passes is True
    assert from_couplings.predicted == pytest.approx(42.0)
    assert from_couplings.diagnostics["scale_source"] == (
        "quark_mass_basis_couplings.M_KK"
    )
    assert from_couplings.passes is True


def test_safe_point_passes_and_excluded_point_fails():
    constraint = fcc.get(_PID)
    limit_gev = constraint.anchor.budget * 1000.0

    safe = constraint.evaluate(point_builder.make_point(kk_ew_mass_gev=limit_gev + 1000.0))
    excluded = constraint.evaluate(
        point_builder.make_point(kk_ew_mass_gev=limit_gev - 1000.0)
    )

    assert safe.passes is True
    assert safe.predicted == pytest.approx((limit_gev + 1000.0) / 1000.0)
    assert safe.ratio < 1.0
    assert excluded.passes is False
    assert excluded.predicted == pytest.approx((limit_gev - 1000.0) / 1000.0)
    assert excluded.ratio > 1.0
    assert excluded.experimental == pytest.approx(constraint.anchor.budget)
    assert excluded.diagnostics["contact_scale_proxy"] == (
        "Lambda_RS = kk_ew_mass_gev, kk_gluon_mass_gev, or M_KK"
    )
    assert "NEEDS-HUMAN-PHYSICS" in excluded.diagnostics["needs_human_physics"]
    assert excluded.diagnostics["eft_recast_status"] == "NEEDS-HUMAN-PHYSICS"
    assert excluded.diagnostics["helicity_interference_matching_status"] == (
        "NEEDS-HUMAN-PHYSICS"
    )
    assert excluded.diagnostics["sigma_or_likelihood_recast_status"] == (
        "NEEDS-HUMAN-PHYSICS"
    )


def test_invalid_contact_scale_is_non_crashing_failure():
    result = fcc.get(_PID).evaluate(point_builder.make_point(kk_ew_mass_gev=-1.0))

    assert result.passes is False
    assert result.predicted is None
    assert result.ratio is None
    assert result.experimental == pytest.approx(fcc.get(_PID).anchor.value)
    assert result.diagnostics["exception_type"] == "ValueError"
    assert result.diagnostics["eft_recast_status"] == "NEEDS-HUMAN-PHYSICS"
    assert result.diagnostics["helicity_interference_matching_status"] == (
        "NEEDS-HUMAN-PHYSICS"
    )


def test_evaluate_is_pure_and_deterministic():
    point = point_builder.make_point(kk_ew_mass_gev=39000.0)
    before = dict(point.extras)
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    assert point.extras == before
    for value in (first.predicted, first.experimental, first.ratio, first.budget):
        assert isinstance(value, float)
        assert math.isfinite(value)
