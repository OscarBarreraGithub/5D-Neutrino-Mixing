"""Production tests for CR002 (VLQ T_5/3 pair-production mass limit)."""

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
from flavor_catalog_constraints.primary.collider_rs import CR002 as cr002_module
from quarkConstraints import collider_resonance as core

_PID = "CR002"
_ACTIVE_VALUE_ID = "PDG2026:CR002:ATLAS2023:X53_pair_Wt"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "collider_rs" / "CR002.yaml"


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
    assert constraint.observable == "m(T_5/3 pair -> tW tW)"


def test_anchor_matches_yaml_and_loud_fail_probe():
    constraint = fcc.get(_PID)
    active = _entry(_ACTIVE_VALUE_ID)
    degenerate = _entry("ATLAS2023:CR002:mass_degenerate_doublet")
    cms_rh = _entry("CMS2019:CR002:X53_RH")
    atlas_pair_only = _entry("ATLAS2018:CR002:T53_pair_only")
    atlas_pair_plus_single = _entry(
        "ATLAS2018:CR002:T53_pair_plus_single_unit_coupling"
    )

    assert constraint.anchor.active_limit.value_gev == pytest.approx(active["value"])
    assert constraint.anchor.active_limit.units == active["units"]
    assert constraint.anchor.active_limit.source_url == active["source_url"]
    assert constraint.anchor.active_limit.display == active["display"]
    assert constraint.anchor.active_limit.model_assumptions == tuple(
        active["model_assumptions"]
    )
    assert constraint.anchor.atlas_2023_degenerate.value_gev == pytest.approx(
        degenerate["value"]
    )
    assert constraint.anchor.cms_2019_rh.value_gev == pytest.approx(cms_rh["value"])
    assert constraint.anchor.atlas_2018_pair_only.value_gev == pytest.approx(
        atlas_pair_only["value"]
    )
    assert constraint.anchor.atlas_2018_pair_plus_single.value_gev == pytest.approx(
        atlas_pair_plus_single["value"]
    )
    assert constraint.anchor.budget == pytest.approx(active["value"] / 1000.0)

    with pytest.raises(AnchorError):
        cr002_module._load_mass_limit_anchor("not-a-real-value-id", process_id=_PID)


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


def test_mass_limit_numerics_validate_against_core_not_adapter():
    constraint = fcc.get(_PID)
    point = _point_with_mkk(1600.0)
    result = constraint.evaluate(point)
    active = _entry(_ACTIVE_VALUE_ID)

    limit = core.ColliderResonanceLimit(
        process_id=_PID,
        resonance="T_5/3 pair",
        final_state="tW tW",
        limit_kind=core.MASS_LOWER_BOUND,
        value=float(active["value"]) / 1000.0,
        units="TeV",
    )
    prediction = core.ColliderResonancePrediction(
        resonance="T_5/3 pair",
        final_state="tW tW",
        mass_tev=1.6,
    )
    direct = core.evaluate_resonance_limit(prediction, limit)

    assert result.predicted == pytest.approx(1.6)
    assert result.predicted == pytest.approx(direct.predicted_mass_tev)
    assert result.experimental == pytest.approx(1.46)
    assert result.ratio == pytest.approx(1.46 / 1.6)
    assert result.ratio == pytest.approx(direct.ratio_to_budget)
    assert result.budget == pytest.approx(direct.budget)
    assert result.passes is True


@pytest.mark.parametrize(
    ("m_kk_gev", "expected_pass"),
    [
        (1800.0, True),
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
    assert result.experimental == pytest.approx(1.46)
    if expected_pass:
        assert result.ratio < 1.0
    else:
        assert result.ratio > 1.0
    assert result.diagnostics["mass_proxy"] == "m_VLQ = M_KK"
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]


def test_invalid_coupling_mass_is_non_crashing_failure():
    point = point_builder.make_point(
        quark_mass_basis_couplings=SimpleNamespace(M_KK=-1.0)
    )
    result = fcc.get(_PID).evaluate(point)

    assert result.passes is False
    assert result.predicted is None
    assert result.ratio is None
    assert result.experimental == pytest.approx(1.46)
    assert result.diagnostics["exception_type"] == "ValueError"


def test_evaluate_is_pure_and_deterministic():
    point = _point_with_mkk(1800.0)
    before = dict(point.extras)
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    assert point.extras == before
    for value in (first.predicted, first.experimental, first.ratio, first.budget):
        assert isinstance(value, float)
        assert math.isfinite(value)
