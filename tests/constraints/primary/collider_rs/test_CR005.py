"""Production tests for CR005 (neutral EW KK dilepton resonance)."""

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
from flavor_catalog_constraints.primary.collider_rs import CR005 as cr005_module
from quarkConstraints import collider_resonance as core

_PID = "CR005"
_ACTIVE_VALUE_ID = "CMS2021:CR005:ZSSM_ll_mass_lower_limit"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "collider_rs" / "CR005.yaml"


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
    assert constraint.observable == "m((gamma^(1), Z^(1))_KK -> l+l-)"


def test_anchor_matches_yaml_and_loud_fail_probe():
    constraint = fcc.get(_PID)
    active = _entry(_ACTIVE_VALUE_ID)
    cms_zpsi = _entry("CMS2021:CR005:Zpsi_ll_mass_lower_limit")
    atlas_zssm = _entry("ATLAS2019:CR005:ZSSM_ll_mass_lower_limit")
    atlas_xsec = _entry("ATLAS2019:CR005:fiducial_xsec_limit_6tev_zero_width")
    pdg_xsec = _entry("PDG2024:CR005:dilepton_zprime_xsec_summary_limit")

    assert constraint.anchor.canonical_value_id == _ACTIVE_VALUE_ID
    assert constraint.anchor.active_limit.value_tev == pytest.approx(active["value"])
    assert constraint.anchor.active_limit.units == active["units"]
    assert constraint.anchor.active_limit.source_url == active["source_url"]
    assert constraint.anchor.active_limit.display == active["display"]
    assert constraint.anchor.active_limit.model_assumption == active["model_assumption"]
    assert constraint.anchor.cms_zpsi.value_tev == pytest.approx(cms_zpsi["value"])
    assert constraint.anchor.atlas_zssm.value_tev == pytest.approx(atlas_zssm["value"])
    assert constraint.anchor.atlas_xsec_6tev.value_fb == pytest.approx(
        atlas_xsec["value"]
    )
    assert constraint.anchor.atlas_xsec_6tev.mass_point_tev == pytest.approx(
        atlas_xsec["mass_point"]["value"]
    )
    assert constraint.anchor.pdg_xsec_summary.value_fb == pytest.approx(
        pdg_xsec["value"]
    )
    assert constraint.anchor.budget == pytest.approx(active["value"])

    with pytest.raises(AnchorError):
        cr005_module._load_value_anchor("not-a-real-value-id", process_id=_PID)


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
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]


def test_mass_limit_numerics_validate_against_core_not_adapter():
    constraint = fcc.get(_PID)
    point = point_builder.make_point(kk_ew_mass_gev=6000.0)
    result = constraint.evaluate(point)
    active = _entry(_ACTIVE_VALUE_ID)

    limit = core.ColliderResonanceLimit(
        process_id=_PID,
        resonance="(gamma^(1), Z^(1))_KK",
        final_state="ee + mumu",
        limit_kind=core.MASS_LOWER_BOUND,
        value=float(active["value"]),
        units=active["units"],
    )
    prediction = core.ColliderResonancePrediction(
        resonance="(gamma^(1), Z^(1))_KK",
        final_state="ee + mumu",
        mass_tev=6.0,
    )
    direct = core.evaluate_resonance_limit(prediction, limit)

    assert result.predicted == pytest.approx(6.0)
    assert result.predicted == pytest.approx(direct.predicted_mass_tev)
    assert result.experimental == pytest.approx(float(active["value"]))
    assert result.ratio == pytest.approx(float(active["value"]) / 6.0)
    assert result.ratio == pytest.approx(direct.ratio_to_budget)
    assert result.budget == pytest.approx(direct.budget)
    assert result.passes is True


def test_mass_falls_back_to_couplings_when_explicit_ew_mass_absent():
    point = point_builder.make_point(
        quark_mass_basis_couplings=SimpleNamespace(M_KK=6200.0)
    )
    result = fcc.get(_PID).evaluate(point)

    assert result.predicted == pytest.approx(6.2)
    assert result.diagnostics["m_kk_gev"] == pytest.approx(6200.0)
    assert result.diagnostics["mass_source"] == "quark_mass_basis_couplings.M_KK"
    assert result.passes is True


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
    assert excluded.diagnostics["mass_proxy"] == (
        "m_(gamma/Z KK) = kk_ew_mass_gev or M_KK"
    )
    assert "NEEDS-HUMAN-PHYSICS" in excluded.diagnostics["needs_human_physics"]
    assert excluded.diagnostics["sigma_times_br_recast_status"] == (
        "NEEDS-HUMAN-PHYSICS"
    )


def test_invalid_ew_mass_is_non_crashing_failure():
    point = point_builder.make_point(kk_ew_mass_gev=-1.0)
    result = fcc.get(_PID).evaluate(point)

    assert result.passes is False
    assert result.predicted is None
    assert result.ratio is None
    assert result.experimental == pytest.approx(fcc.get(_PID).anchor.value)
    assert result.diagnostics["exception_type"] == "ValueError"


def test_evaluate_is_pure_and_deterministic():
    point = point_builder.make_point(kk_ew_mass_gev=6500.0)
    before = dict(point.extras)
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    assert point.extras == before
    for value in (first.predicted, first.experimental, first.ratio, first.budget):
        assert isinstance(value, float)
        assert math.isfinite(value)
