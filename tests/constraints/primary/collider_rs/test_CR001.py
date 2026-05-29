"""Production tests for CR001 (KK-gluon -> ttbar resonance)."""

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
from flavor_catalog_constraints.primary.collider_rs import CR001 as cr001_module
from quarkConstraints import collider_resonance as core

_PID = "CR001"
_ACTIVE_VALUE_ID = "CMS2026:CR001:gkk_ttbar_mass_exclusion"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "collider_rs" / "CR001.yaml"


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
    assert constraint.observable == "m(g_KK -> t tbar)"


def test_anchor_matches_yaml_and_loud_fail_probe():
    constraint = fcc.get(_PID)
    active = _entry(_ACTIVE_VALUE_ID)
    cms_2019 = _entry("CMS2019:CR001:gkk_ttbar_mass_exclusion")
    pdg = _entry("PDGLive2026:CR001:s071kkg_aaboud2018bi")

    assert constraint.anchor.active_limit.value == pytest.approx(active["value"])
    assert constraint.anchor.active_limit.units == active["units"]
    assert constraint.anchor.active_limit.source_url == active["source_url"]
    assert constraint.anchor.active_limit.display == active["display"]
    assert constraint.anchor.active_limit.mass_interval_low == pytest.approx(
        active["benchmark_assumptions"]["mass_interval_low"]
    )
    assert constraint.anchor.cms_2019.value == pytest.approx(cms_2019["value"])
    assert constraint.anchor.pdg_2026.value == pytest.approx(pdg["value"])
    assert constraint.anchor.budget == pytest.approx(5.5)

    with pytest.raises(AnchorError):
        cr001_module._load_mass_limit_anchor("not-a-real-value-id", process_id=_PID)


def test_evaluate_without_input_degrades_gracefully():
    result = fcc.get(_PID).evaluate(point_builder.empty_point())

    assert result.process_id == _PID
    assert result.passes is True
    assert result.predicted is None
    assert result.ratio is None
    assert result.experimental == pytest.approx(fcc.get(_PID).anchor.value)
    assert result.budget == pytest.approx(fcc.get(_PID).anchor.budget)
    assert result.diagnostics["missing_extras"] == (
        "kk_gluon_mass_gev",
        "quark_mass_basis_couplings",
    )


def test_mass_limit_numerics_validate_against_core_not_adapter():
    constraint = fcc.get(_PID)
    point = point_builder.make_point(kk_gluon_mass_gev=6000.0)
    result = constraint.evaluate(point)
    active = _entry(_ACTIVE_VALUE_ID)

    limit = core.ColliderResonanceLimit(
        process_id=_PID,
        resonance="g_KK^(1)",
        final_state="t tbar",
        limit_kind=core.MASS_LOWER_BOUND,
        value=float(active["value"]),
        units=active["units"],
    )
    prediction = core.kk_gluon_prediction_from_m_kk_gev(6000.0)
    direct = core.evaluate_resonance_limit(prediction, limit)

    assert result.predicted == pytest.approx(6.0)
    assert result.predicted == pytest.approx(direct.predicted_mass_tev)
    assert result.ratio == pytest.approx(5.5 / 6.0)
    assert result.ratio == pytest.approx(direct.ratio_to_budget)
    assert result.budget == pytest.approx(direct.budget)
    assert result.passes is True


def test_mass_can_be_derived_from_quark_couplings_m_kk():
    point = point_builder.build_from_quark_couplings(SimpleNamespace(M_KK=6200.0))
    result = fcc.get(_PID).evaluate(point)

    assert result.predicted == pytest.approx(6.2)
    assert result.diagnostics["m_kk_gev"] == pytest.approx(6200.0)
    assert result.diagnostics["mass_source"] == "kk_gluon_mass_gev"
    assert result.passes is True


def test_mass_falls_back_to_couplings_when_explicit_extra_absent():
    point = point_builder.make_point(
        quark_mass_basis_couplings=SimpleNamespace(M_KK=6300.0)
    )
    result = fcc.get(_PID).evaluate(point)

    assert result.predicted == pytest.approx(6.3)
    assert result.diagnostics["m_kk_gev"] == pytest.approx(6300.0)
    assert result.diagnostics["mass_source"] == "quark_mass_basis_couplings.M_KK"
    assert result.passes is True


@pytest.mark.parametrize(
    ("m_kk_gev", "expected_pass"),
    [
        (6500.0, True),
        (3000.0, False),
    ],
)
def test_safe_point_passes_and_excluded_point_fails(
    m_kk_gev: float,
    expected_pass: bool,
):
    result = fcc.get(_PID).evaluate(
        point_builder.make_point(kk_gluon_mass_gev=m_kk_gev)
    )

    assert result.passes is expected_pass
    assert result.predicted == pytest.approx(m_kk_gev / 1000.0)
    assert result.experimental == pytest.approx(5.5)
    if expected_pass:
        assert result.ratio < 1.0
    else:
        assert result.ratio > 1.0
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]


def test_evaluate_is_pure_and_deterministic():
    point = point_builder.make_point(kk_gluon_mass_gev=6500.0)
    before = dict(point.extras)
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    assert point.extras == before
    for value in (first.predicted, first.experimental, first.ratio, first.budget):
        assert isinstance(value, float)
        assert math.isfinite(value)
