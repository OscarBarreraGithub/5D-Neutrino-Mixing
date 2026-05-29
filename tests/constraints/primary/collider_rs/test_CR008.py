"""Production tests for CR008 (singlet VLQ T pair-production mass limit)."""

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
from flavor_catalog_constraints.primary.collider_rs import CR008 as cr008_module
from quarkConstraints import collider_resonance as core

_PID = "CR008"
_ACTIVE_VALUE_ID = "ATLAS2024:CR008:T_singlet_pair_mass_limit"
_CMS_ENVELOPE_VALUE_ID = "CMS2023:CR008:T_pair_all_third_generation_decays_envelope"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "collider_rs" / "CR008.yaml"


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
    assert constraint.observable == "m(T singlet pair -> Wb/Ht/Zt)"


def test_anchor_matches_yaml_and_loud_fail_probe():
    constraint = fcc.get(_PID)
    active = _entry(_ACTIVE_VALUE_ID)
    cms_envelope = _entry(_CMS_ENVELOPE_VALUE_ID)

    assert constraint.anchor.active_limit.value_tev == pytest.approx(active["value"])
    assert constraint.anchor.active_limit.value_gev == pytest.approx(
        active["value"] * 1000.0
    )
    assert constraint.anchor.active_limit.units == active["units"]
    assert constraint.anchor.active_limit.source_url == active["source_url"]
    assert constraint.anchor.active_limit.display == active["display"]
    assert constraint.anchor.active_limit.experiment == active["experiment"]
    assert constraint.anchor.active_limit.sqrt_s == active["sqrt_s"]
    assert constraint.anchor.active_limit.model_assumptions == (
        active["model_assumptions"]
    )
    assert constraint.anchor.active_limit.pdg_crosscheck_source_url == (
        active["pdg_crosscheck"]["source_url"]
    )
    assert constraint.anchor.cms_all_third_generation.value_tev == pytest.approx(
        cms_envelope["value"]
    )
    assert constraint.anchor.cms_all_third_generation.model_assumptions == (
        cms_envelope["model_assumptions"]
    )
    assert constraint.anchor.budget == pytest.approx(active["value"])

    with pytest.raises(AnchorError):
        cr008_module._load_mass_limit_anchor("not-a-real-value-id", process_id=_PID)


def test_evaluate_without_input_degrades_gracefully():
    result = fcc.get(_PID).evaluate(point_builder.empty_point())

    assert result.process_id == _PID
    assert result.passes is True
    assert result.predicted is None
    assert result.ratio is None
    assert result.experimental == pytest.approx(fcc.get(_PID).anchor.value)
    assert result.budget == pytest.approx(fcc.get(_PID).anchor.budget)
    assert result.diagnostics["missing_extra"] == "quark_mass_basis_couplings"
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]
    assert result.diagnostics["sigma_times_br_recast_status"] == (
        "NEEDS-HUMAN-PHYSICS"
    )


def test_mass_limit_numerics_validate_against_core_not_adapter():
    constraint = fcc.get(_PID)
    point = _point_with_mkk(1600.0)
    result = constraint.evaluate(point)
    active = _entry(_ACTIVE_VALUE_ID)

    limit = core.ColliderResonanceLimit(
        process_id=_PID,
        resonance="T pair",
        final_state="Wb/Ht/Zt singlet",
        limit_kind=core.MASS_LOWER_BOUND,
        value=float(active["value"]),
        units=active["units"],
    )
    prediction = core.ColliderResonancePrediction(
        resonance="T pair",
        final_state="Wb/Ht/Zt singlet",
        mass_tev=1.6,
    )
    direct = core.evaluate_resonance_limit(prediction, limit)

    assert result.predicted == pytest.approx(1.6)
    assert result.predicted == pytest.approx(direct.predicted_mass_tev)
    assert result.experimental == pytest.approx(float(active["value"]))
    assert result.ratio == pytest.approx(float(active["value"]) / 1.6)
    assert result.ratio == pytest.approx(direct.ratio_to_budget)
    assert result.budget == pytest.approx(direct.budget)
    assert result.passes is True
    assert result.diagnostics["active_value_id"] == _ACTIVE_VALUE_ID
    assert result.diagnostics["active_branching_fractions"] == (
        active["model_assumptions"]["branching_fractions"]
    )


@pytest.mark.parametrize(
    ("m_kk_gev", "expected_pass"),
    [
        (1500.0, True),
        (1200.0, False),
    ],
)
def test_safe_point_passes_and_excluded_point_fails(
    m_kk_gev: float,
    expected_pass: bool,
):
    result = fcc.get(_PID).evaluate(_point_with_mkk(m_kk_gev))

    assert result.passes is expected_pass
    assert result.predicted == pytest.approx(m_kk_gev / 1000.0)
    assert result.experimental == pytest.approx(_entry(_ACTIVE_VALUE_ID)["value"])
    if expected_pass:
        assert result.ratio < 1.0
    else:
        assert result.ratio > 1.0
    assert result.diagnostics["mass_proxy"] == "m_T = M_KK"
    assert result.diagnostics["mass_source"] == "quark_mass_basis_couplings.M_KK"
    assert result.diagnostics["sigma_times_br_recast_status"] == (
        "NEEDS-HUMAN-PHYSICS"
    )
    assert result.diagnostics["all_mass_limits_tev"][
        _CMS_ENVELOPE_VALUE_ID
    ] == pytest.approx(_entry(_CMS_ENVELOPE_VALUE_ID)["value"])


def test_invalid_coupling_mass_is_non_crashing_failure():
    result = fcc.get(_PID).evaluate(_point_with_mkk(-1.0))

    assert result.passes is False
    assert result.predicted is None
    assert result.ratio is None
    assert result.experimental == pytest.approx(_entry(_ACTIVE_VALUE_ID)["value"])
    assert result.diagnostics["exception_type"] == "ValueError"
    assert result.diagnostics["sigma_times_br_recast_status"] == (
        "NEEDS-HUMAN-PHYSICS"
    )


def test_evaluate_is_pure_and_deterministic():
    point = _point_with_mkk(1500.0)
    before = dict(point.extras)
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    assert point.extras == before
    for value in (first.predicted, first.experimental, first.ratio, first.budget):
        assert isinstance(value, float)
        assert math.isfinite(value)
