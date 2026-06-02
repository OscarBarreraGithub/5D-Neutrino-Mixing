"""Production tests for CR010 ((T,B)-doublet VLQ pair-production limit)."""

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
from flavor_catalog_constraints.primary.collider_rs import CR010 as cr010_module
from quarkConstraints import collider_resonance as core

_PID = "CR010"
_ACTIVE_T_VALUE_ID = "PDG2025:CR010:TB_doublet_T_mass_lower_limit"
_ACTIVE_B_VALUE_ID = "PDG2025:CR010:TB_doublet_B_mass_lower_limit"
_T_WB_VALUE_ID = "PDG2025:CR010:T_to_Wb_100pct_mass_lower_limit"
_B_HB_VALUE_ID = "PDG2025:CR010:B_to_Hb_100pct_mass_lower_limit"
_B_WT_VALUE_ID = "PDG2025:CR010:B_to_Wt_100pct_mass_lower_limit"
_B_ZB_VALUE_ID = "PDG2025:CR010:B_to_Zb_100pct_mass_lower_limit"
_CMS_T_UNIFORM_VALUE_ID = (
    "CMS2023:CR010:T_all_third_generation_decays_uniform_lower_limit"
)
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "collider_rs" / "CR010.yaml"


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
    assert constraint.observable == (
        "m((T,B) doublet pair -> mixed third-generation final states)"
    )


def test_anchor_matches_yaml_and_loud_fail_probe():
    constraint = fcc.get(_PID)
    active_t = _entry(_ACTIVE_T_VALUE_ID)
    active_b = _entry(_ACTIVE_B_VALUE_ID)
    t_wb = _entry(_T_WB_VALUE_ID)
    b_hb = _entry(_B_HB_VALUE_ID)
    b_wt = _entry(_B_WT_VALUE_ID)
    b_zb = _entry(_B_ZB_VALUE_ID)
    cms_uniform = _entry(_CMS_T_UNIFORM_VALUE_ID)

    assert constraint.anchor.t_doublet_limit.value_tev == pytest.approx(
        active_t["value"]
    )
    assert constraint.anchor.t_doublet_limit.value_gev == pytest.approx(
        active_t["value"] * 1000.0
    )
    assert constraint.anchor.t_doublet_limit.units == active_t["units"]
    assert constraint.anchor.t_doublet_limit.source_url == active_t["source_url"]
    assert constraint.anchor.t_doublet_limit.display == active_t["display"]
    assert constraint.anchor.t_doublet_limit.experiment == active_t["experiment"]
    assert constraint.anchor.t_doublet_limit.source_key == active_t["source_key"]
    assert constraint.anchor.t_doublet_limit.original_result_url == (
        active_t["original_result_url"]
    )
    assert constraint.anchor.t_doublet_limit.model_assumption == (
        active_t["model_assumption"]
    )
    assert constraint.anchor.b_doublet_limit.value_tev == pytest.approx(
        active_b["value"]
    )
    assert constraint.anchor.b_doublet_limit.source_url == active_b["source_url"]
    assert constraint.anchor.b_doublet_limit.model_assumption == (
        active_b["model_assumption"]
    )
    assert constraint.anchor.budget == pytest.approx(
        max(active_t["value"], active_b["value"])
    )
    assert constraint.anchor.active_value_ids == (
        _ACTIVE_T_VALUE_ID,
        _ACTIVE_B_VALUE_ID,
    )

    assert constraint.anchor.t_wb_100pct.value_tev == pytest.approx(t_wb["value"])
    assert constraint.anchor.b_hb_100pct.value_tev == pytest.approx(b_hb["value"])
    assert constraint.anchor.b_wt_100pct.value_tev == pytest.approx(b_wt["value"])
    assert constraint.anchor.b_zb_100pct.value_tev == pytest.approx(b_zb["value"])
    assert constraint.anchor.cms_t_uniform.value_tev == pytest.approx(
        cms_uniform["value"]
    )
    assert "not a simultaneous" in constraint.anchor.t_wb_100pct.model_assumption
    assert "not a simultaneous" in constraint.anchor.b_wt_100pct.model_assumption

    with pytest.raises(AnchorError):
        cr010_module._load_mass_limit_anchor("not-a-real-value-id", process_id=_PID)


def test_evaluate_without_input_degrades_gracefully():
    result = fcc.get(_PID).evaluate(point_builder.empty_point())

    assert result.process_id == _PID
    assert result.passes is True
    assert result.predicted is None
    assert result.ratio is None
    assert result.experimental == pytest.approx(fcc.get(_PID).anchor.value)
    assert result.budget == pytest.approx(fcc.get(_PID).anchor.budget)
    assert result.diagnostics["missing_extra"] == "quark_mass_basis_couplings"
    assert result.diagnostics["active_value_ids"] == (
        _ACTIVE_T_VALUE_ID,
        _ACTIVE_B_VALUE_ID,
    )
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]
    assert result.diagnostics["sigma_times_br_recast_status"] == (
        "NEEDS-HUMAN-PHYSICS"
    )


def test_mass_limit_numerics_validate_against_core_not_adapter():
    constraint = fcc.get(_PID)
    point = _point_with_mkk(1600.0)
    result = constraint.evaluate(point)
    active_t = _entry(_ACTIVE_T_VALUE_ID)
    active_b = _entry(_ACTIVE_B_VALUE_ID)
    active_budget = max(float(active_t["value"]), float(active_b["value"]))

    limit = core.ColliderResonanceLimit(
        process_id=_PID,
        resonance="(T,B) doublet pair",
        final_state="mixed W/Z/H third-generation final states",
        limit_kind=core.MASS_LOWER_BOUND,
        value=active_budget,
        units="TeV",
    )
    prediction = core.ColliderResonancePrediction(
        resonance="(T,B) doublet pair",
        final_state="mixed W/Z/H third-generation final states",
        mass_tev=1.6,
    )
    direct = core.evaluate_resonance_limit(prediction, limit)

    assert result.predicted == pytest.approx(1.6)
    assert result.predicted == pytest.approx(direct.predicted_mass_tev)
    assert result.experimental == pytest.approx(active_budget)
    assert result.ratio == pytest.approx(active_budget / 1.6)
    assert result.ratio == pytest.approx(direct.ratio_to_budget)
    assert result.budget == pytest.approx(direct.budget)
    assert result.passes is True
    assert result.diagnostics["active_value_ids"] == (
        _ACTIVE_T_VALUE_ID,
        _ACTIVE_B_VALUE_ID,
    )
    assert result.diagnostics["active_t_limit_tev"] == pytest.approx(
        active_t["value"]
    )
    assert result.diagnostics["active_b_limit_tev"] == pytest.approx(
        active_b["value"]
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
    assert result.experimental == pytest.approx(_entry(_ACTIVE_T_VALUE_ID)["value"])
    if expected_pass:
        assert result.ratio < 1.0
    else:
        assert result.ratio > 1.0
    assert result.diagnostics["mass_proxy"] == "m_T = m_B = M_KK"
    assert result.diagnostics["mass_source"] == "quark_mass_basis_couplings.M_KK"
    assert result.diagnostics["sigma_times_br_recast_status"] == (
        "NEEDS-HUMAN-PHYSICS"
    )
    assert result.diagnostics["all_mass_limits_tev"][
        _CMS_T_UNIFORM_VALUE_ID
    ] == pytest.approx(_entry(_CMS_T_UNIFORM_VALUE_ID)["value"])
    assert result.diagnostics["individual_endpoint_limits_tev"][
        _T_WB_VALUE_ID
    ] == pytest.approx(_entry(_T_WB_VALUE_ID)["value"])


def test_invalid_coupling_mass_is_non_crashing_failure():
    result = fcc.get(_PID).evaluate(_point_with_mkk(-1.0))

    assert result.passes is False
    assert result.predicted is None
    assert result.ratio is None
    assert result.experimental == pytest.approx(_entry(_ACTIVE_T_VALUE_ID)["value"])
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
