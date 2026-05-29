"""Production tests for CR006 (charged-current EW KK/W' resonance)."""

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
from flavor_catalog_constraints.primary.collider_rs import CR006 as cr006_module
from quarkConstraints import collider_resonance as core

_PID = "CR006"
_ACTIVE_VALUE_ID = "PDG2025:CR006:Wprime_SSM_enu_mass_lower_bound"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "collider_rs" / "CR006.yaml"


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
    assert constraint.observable == "m(W_KK^(1) -> ell nu, tb)"


def test_anchor_matches_yaml_and_loud_fail_probe():
    constraint = fcc.get(_PID)
    active = _entry(_ACTIVE_VALUE_ID)
    muon = _entry("PDG2025:CR006:Wprime_SSM_munu_mass_lower_bound")
    combined = _entry("CMS2022:CR006:Wprime_SSM_enu_munu_combined_mass_lower_bound")
    tb_r = _entry("CMS2024:CR006:Wprime_R_tb_mass_lower_bound")
    tb_l = _entry("CMS2024:CR006:Wprime_L_tb_mass_lower_bound")

    assert constraint.anchor.active_limit.value_tev == pytest.approx(active["value"])
    assert constraint.anchor.active_limit.units == active["units"]
    assert constraint.anchor.active_limit.source_url == active["source_url"]
    assert constraint.anchor.active_limit.display == active["display"]
    assert constraint.anchor.active_limit.experiment == active["experiment"]
    assert constraint.anchor.pdg_cms_munu.value_tev == pytest.approx(muon["value"])
    assert constraint.anchor.cms_combined_lnu.value_tev == pytest.approx(
        combined["value"]
    )
    assert constraint.anchor.cms_tb_r.value_tev == pytest.approx(tb_r["value"])
    assert constraint.anchor.cms_tb_r.model_assumption == tb_r["model_assumption"]
    assert constraint.anchor.cms_tb_l.value_tev == pytest.approx(tb_l["value"])
    assert constraint.anchor.cms_tb_l.model_assumption == tb_l["model_assumption"]
    assert constraint.anchor.budget == pytest.approx(active["value"])

    with pytest.raises(AnchorError):
        cr006_module._load_value_anchor("not-a-real-value-id", process_id=_PID)


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
        resonance="W_KK^(1)",
        final_state="ell nu",
        limit_kind=core.MASS_LOWER_BOUND,
        value=float(active["value"]),
        units=active["units"],
    )
    prediction = core.ColliderResonancePrediction(
        resonance="W_KK^(1)",
        final_state="ell nu",
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
        quark_mass_basis_couplings=SimpleNamespace(M_KK=6500.0)
    )
    result = fcc.get(_PID).evaluate(point)

    assert result.predicted == pytest.approx(6.5)
    assert result.diagnostics["m_kk_gev"] == pytest.approx(6500.0)
    assert result.diagnostics["mass_source"] == "quark_mass_basis_couplings.M_KK"
    assert result.passes is True


def test_safe_point_passes_and_excluded_point_fails():
    constraint = fcc.get(_PID)
    limit_gev = constraint.anchor.budget * 1000.0

    safe = constraint.evaluate(
        point_builder.make_point(kk_ew_mass_gev=limit_gev + 1000.0)
    )
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
    assert excluded.diagnostics["mass_proxy"] == "m_WKK = kk_ew_mass_gev or M_KK"
    assert "NEEDS-HUMAN-PHYSICS" in excluded.diagnostics["needs_human_physics"]
    assert excluded.diagnostics["sigma_times_br_recast_status"] == (
        "NEEDS-HUMAN-PHYSICS"
    )
    assert excluded.diagnostics["tb_mass_limits_tev"][
        "CMS2024:CR006:Wprime_R_tb_mass_lower_bound"
    ] == pytest.approx(4.3)


def test_invalid_ew_mass_is_non_crashing_failure():
    point = point_builder.make_point(kk_ew_mass_gev=-1.0)
    result = fcc.get(_PID).evaluate(point)

    assert result.passes is False
    assert result.predicted is None
    assert result.ratio is None
    assert result.experimental == pytest.approx(fcc.get(_PID).anchor.value)
    assert result.diagnostics["exception_type"] == "ValueError"


def test_evaluate_is_pure_and_deterministic():
    point = point_builder.make_point(kk_ew_mass_gev=7000.0)
    before = dict(point.extras)
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    assert point.extras == before
    for value in (first.predicted, first.experimental, first.ratio, first.budget):
        assert isinstance(value, float)
        assert math.isfinite(value)
