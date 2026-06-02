"""Production tests for CR014 (top-philic vector four-top limit)."""

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
from flavor_catalog_constraints.primary.collider_rs import CR014 as cr014_module
from quarkConstraints import collider_resonance as core

_PID = "CR014"
_ACTIVE_VALUE_ID = (
    "CMSB2G25005:CR014:top_philic_vector_50pct_width_excluded_up_to"
)
_ATLAS_4TOP_XSEC_ID = "PDG2025:CR014:ATLAS_4top_cross_section"
_CMS_4TOP_XSEC_ID = "PDG2025:CR014:CMS_4top_cross_section"
_ATLAS_TOP_PHILIC_RANGE_ID = (
    "ATLAS2024:CR014:top_philic_Zprime_cross_section_limit_range"
)
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "collider_rs" / "CR014.yaml"


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


def _point_with_mkk(m_kk_gev: float):
    return point_builder.make_point(
        quark_mass_basis_couplings=SimpleNamespace(M_KK=m_kk_gev)
    )


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.HARD
    assert constraint.family == "collider_rs"
    assert constraint.observable == "m(top-philic vector Z' -> t tbar t tbar)"


def test_anchor_matches_yaml_and_loud_fail_probe():
    constraint = fcc.get(_PID)
    active = _entry(_ACTIVE_VALUE_ID)
    atlas_xsec = _entry(_ATLAS_4TOP_XSEC_ID)
    cms_xsec = _entry(_CMS_4TOP_XSEC_ID)
    atlas_range = _entry(_ATLAS_TOP_PHILIC_RANGE_ID)

    assert constraint.anchor.active_limit.value_gev == pytest.approx(active["value"])
    assert constraint.anchor.active_limit.value_tev == pytest.approx(
        active["value"] / 1000.0
    )
    assert constraint.anchor.active_limit.units == active["units"]
    assert constraint.anchor.active_limit.source_url == active["source_url"]
    assert constraint.anchor.active_limit.display == active["display"]
    assert constraint.anchor.active_limit.experiment == active["experiment"]
    assert constraint.anchor.active_limit.source_key == active["source_key"]
    assert constraint.anchor.active_limit.width_over_mass == pytest.approx(
        active["assumptions"]["width_over_mass"]
    )
    assert constraint.anchor.active_limit.resonance_type == (
        active["assumptions"]["resonance_type"]
    )
    assert constraint.anchor.active_limit.coupling_pattern == (
        active["assumptions"]["coupling_pattern"]
    )
    assert constraint.anchor.budget == pytest.approx(active["value"] / 1000.0)
    assert constraint.anchor.width_over_mass == pytest.approx(0.50)
    assert constraint.anchor.all_applicable_mass_limits[0].value_id == _ACTIVE_VALUE_ID

    assert constraint.anchor.atlas_four_top_cross_section.value_fb == pytest.approx(
        atlas_xsec["value"]
    )
    assert constraint.anchor.atlas_four_top_cross_section.uncertainty_plus == (
        pytest.approx(atlas_xsec["uncertainty_plus"])
    )
    assert constraint.anchor.cms_four_top_cross_section.value_fb == pytest.approx(
        cms_xsec["value"]
    )
    assert constraint.anchor.cms_four_top_cross_section.uncertainty_minus == (
        pytest.approx(cms_xsec["uncertainty_minus"])
    )
    assert (
        constraint.anchor.atlas_top_philic_cross_section_range.value_min_fb
        == pytest.approx(atlas_range["value_min"])
    )
    assert (
        constraint.anchor.atlas_top_philic_cross_section_range.value_max_fb
        == pytest.approx(atlas_range["value_max"])
    )

    with pytest.raises(AnchorError):
        cr014_module._load_mass_limit_anchor("not-a-real-value-id", process_id=_PID)
    with pytest.raises(AnchorError):
        cr014_module._load_cross_section_anchor("not-a-real-value-id", process_id=_PID)


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
    assert result.diagnostics["active_limit_gev"] == pytest.approx(
        _entry(_ACTIVE_VALUE_ID)["value"]
    )
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]
    assert result.diagnostics["width_dependence_recast_status"] == (
        "NEEDS-HUMAN-PHYSICS"
    )
    assert result.diagnostics["sm_four_top_background_recast_status"] == (
        "NEEDS-HUMAN-PHYSICS"
    )


def test_mass_limit_numerics_validate_against_core_not_adapter():
    constraint = fcc.get(_PID)
    point = point_builder.make_point(kk_ew_mass_gev=1000.0)
    result = constraint.evaluate(point)
    active = _entry(_ACTIVE_VALUE_ID)
    budget_tev = float(active["value"]) / 1000.0

    limit = core.ColliderResonanceLimit(
        process_id=_PID,
        resonance="top-philic vector mediator Z'",
        final_state="t tbar t tbar two-lepton",
        limit_kind=core.MASS_LOWER_BOUND,
        value=budget_tev,
        units="TeV",
    )
    prediction = core.ColliderResonancePrediction(
        resonance="top-philic vector mediator Z'",
        final_state="t tbar t tbar two-lepton",
        mass_tev=1.0,
    )
    direct = core.evaluate_resonance_limit(prediction, limit)

    assert result.predicted == pytest.approx(1.0)
    assert result.predicted == pytest.approx(direct.predicted_mass_tev)
    assert result.experimental == pytest.approx(budget_tev)
    assert result.ratio == pytest.approx(budget_tev / 1.0)
    assert result.ratio == pytest.approx(direct.ratio_to_budget)
    assert result.budget == pytest.approx(direct.budget)
    assert result.passes is True
    assert result.diagnostics["active_value_id"] == _ACTIVE_VALUE_ID
    assert result.diagnostics["active_width_over_mass"] == pytest.approx(
        active["assumptions"]["width_over_mass"]
    )


def test_mass_falls_back_to_couplings_when_explicit_ew_mass_absent():
    result = fcc.get(_PID).evaluate(_point_with_mkk(900.0))

    assert result.predicted == pytest.approx(0.9)
    assert result.diagnostics["mass_source"] == "quark_mass_basis_couplings.M_KK"
    assert result.diagnostics["m_top_philic_vector_proxy_gev"] == pytest.approx(900.0)
    assert result.passes is True


def test_safe_point_passes_and_excluded_point_fails():
    constraint = fcc.get(_PID)
    limit_gev = constraint.anchor.value_gev

    safe = constraint.evaluate(point_builder.make_point(kk_ew_mass_gev=limit_gev + 50.0))
    excluded = constraint.evaluate(
        point_builder.make_point(kk_ew_mass_gev=limit_gev - 50.0)
    )

    assert safe.passes is True
    assert safe.predicted == pytest.approx((limit_gev + 50.0) / 1000.0)
    assert safe.ratio < 1.0
    assert excluded.passes is False
    assert excluded.predicted == pytest.approx((limit_gev - 50.0) / 1000.0)
    assert excluded.ratio > 1.0
    assert excluded.experimental == pytest.approx(constraint.anchor.budget)
    assert excluded.diagnostics["mass_proxy"] == (
        "m_Zprime_top_philic = kk_ew_mass_gev or M_KK"
    )
    assert excluded.diagnostics["mass_source"] == "kk_ew_mass_gev"
    assert "NEEDS-HUMAN-PHYSICS" in excluded.diagnostics["needs_human_physics"]
    assert excluded.diagnostics["sigma_times_br_recast_status"] == (
        "NEEDS-HUMAN-PHYSICS"
    )
    assert excluded.diagnostics["four_top_acceptance_recast_status"] == (
        "NEEDS-HUMAN-PHYSICS"
    )
    assert excluded.diagnostics["all_applicable_mass_limits_gev"][
        _ACTIVE_VALUE_ID
    ] == pytest.approx(_entry(_ACTIVE_VALUE_ID)["value"])
    assert excluded.diagnostics["cms_four_top_cross_section_fb"] == pytest.approx(
        _entry(_CMS_4TOP_XSEC_ID)["value"]
    )
    assert excluded.diagnostics["atlas_top_philic_cross_section_limit_range_fb"] == (
        pytest.approx(_entry(_ATLAS_TOP_PHILIC_RANGE_ID)["value_min"]),
        pytest.approx(_entry(_ATLAS_TOP_PHILIC_RANGE_ID)["value_max"]),
    )


def test_invalid_ew_mass_is_non_crashing_failure():
    result = fcc.get(_PID).evaluate(point_builder.make_point(kk_ew_mass_gev=-1.0))

    assert result.passes is False
    assert result.predicted is None
    assert result.ratio is None
    assert result.experimental == pytest.approx(fcc.get(_PID).anchor.value)
    assert result.diagnostics["exception_type"] == "ValueError"
    assert result.diagnostics["width_dependence_recast_status"] == (
        "NEEDS-HUMAN-PHYSICS"
    )


def test_evaluate_is_pure_and_deterministic():
    point = point_builder.make_point(kk_ew_mass_gev=1000.0)
    before = dict(point.extras)
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    assert point.extras == before
    for value in (first.predicted, first.experimental, first.ratio, first.budget):
        assert isinstance(value, float)
        assert math.isfinite(value)
