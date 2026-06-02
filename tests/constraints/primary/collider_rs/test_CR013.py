"""Production tests for CR013 (diphoton spin-2 KK-graviton limit)."""

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
from flavor_catalog_constraints.primary.collider_rs import CR013 as cr013_module
from quarkConstraints.scales import GAUGE_KK_ROOT_NN, SPIN2_GRAVITON_KK_ROOT
from quarkConstraints import collider_resonance as core

_PID = "CR013"
_ACTIVE_VALUE_ID = "PDG2025:CR013:CMS2024_RSG_diphoton_kMPl_0p1"
_CMS_2024_RANGE_ID = (
    "CMS2024:CR013:RSG_diphoton_exclusion_range_ktilde_0p01_0p2"
)
_ATLAS_2021_KTILDE_0P1_ID = "ATLAS2021:CR013:RSG_diphoton_kMPl_0p1"
_ATLAS_2021_KTILDE_0P05_ID = "ATLAS2021:CR013:RSG_diphoton_kMPl_0p05"
_ATLAS_2021_KTILDE_0P01_ID = "ATLAS2021:CR013:RSG_diphoton_kMPl_0p01"
_CMS_2018_RANGE_ID = (
    "CMS2018:CR013:RSG_diphoton_exclusion_range_ktilde_0p01_0p2"
)
_ATLAS_2015_KTILDE_0P1_ID = "ATLAS2015:CR013:RSG_diphoton_kMPl_0p1"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "collider_rs" / "CR013.yaml"
_GRAVITON_TO_GAUGE_ROOT_RATIO = SPIN2_GRAVITON_KK_ROOT / GAUGE_KK_ROOT_NN


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
    assert constraint.observable == "m(G_KK^(1) -> gamma gamma)"


def test_anchor_matches_yaml_and_loud_fail_probe():
    constraint = fcc.get(_PID)
    active = _entry(_ACTIVE_VALUE_ID)
    cms_range = _entry(_CMS_2024_RANGE_ID)
    atlas_2021_0p1 = _entry(_ATLAS_2021_KTILDE_0P1_ID)
    atlas_2021_0p05 = _entry(_ATLAS_2021_KTILDE_0P05_ID)
    atlas_2021_0p01 = _entry(_ATLAS_2021_KTILDE_0P01_ID)
    cms_2018_range = _entry(_CMS_2018_RANGE_ID)
    atlas_2015_0p1 = _entry(_ATLAS_2015_KTILDE_0P1_ID)

    assert constraint.anchor.active_limit.value_tev == pytest.approx(active["value"])
    assert constraint.anchor.active_limit.units == active["units"]
    assert constraint.anchor.active_limit.source_url == active["source_url"]
    assert constraint.anchor.active_limit.display == active["display"]
    assert constraint.anchor.active_limit.experiment == active["experiment"]
    assert constraint.anchor.active_limit.analysis_year == active["analysis_year"]
    assert constraint.anchor.active_limit.model == active["model"]
    assert constraint.anchor.benchmark_ktilde == pytest.approx(0.1)
    assert constraint.anchor.budget == pytest.approx(active["value"])

    assert constraint.anchor.cms_2024_range.value_tev == pytest.approx(
        cms_range["value"]
    )
    assert constraint.anchor.cms_2024_range.range_min_tev == pytest.approx(
        cms_range["range_min"]
    )
    assert constraint.anchor.cms_2024_range.range_max_tev == pytest.approx(
        cms_range["range_max"]
    )
    assert constraint.anchor.atlas_2021_ktilde_0p1.value_tev == pytest.approx(
        atlas_2021_0p1["value"]
    )
    assert constraint.anchor.atlas_2021_ktilde_0p05.value_tev == pytest.approx(
        atlas_2021_0p05["value"]
    )
    assert constraint.anchor.atlas_2021_ktilde_0p01.value_tev == pytest.approx(
        atlas_2021_0p01["value"]
    )
    assert constraint.anchor.cms_2018_range.value_tev == pytest.approx(
        cms_2018_range["value"]
    )
    assert constraint.anchor.atlas_2015_ktilde_0p1.value_tev == pytest.approx(
        atlas_2015_0p1["value"]
    )
    assert constraint.anchor.strongest_reported_endpoint.value_tev == pytest.approx(
        cms_range["value"]
    )

    with pytest.raises(AnchorError):
        cr013_module._load_mass_limit_anchor("not-a-real-value-id", process_id=_PID)


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
    assert result.diagnostics["active_value_id"] == _ACTIVE_VALUE_ID
    assert result.diagnostics["benchmark_ktilde"] == pytest.approx(0.1)
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]
    assert result.diagnostics["sigma_times_br_recast_status"] == (
        "NEEDS-HUMAN-PHYSICS"
    )
    assert result.diagnostics["coupling_recast_status"] == "NEEDS-HUMAN-PHYSICS"


def test_mass_limit_numerics_validate_against_core_not_adapter():
    constraint = fcc.get(_PID)
    point = point_builder.make_point(kk_ew_mass_gev=3100.0)
    result = constraint.evaluate(point)
    active = _entry(_ACTIVE_VALUE_ID)
    expected_mass_tev = 3.1 * _GRAVITON_TO_GAUGE_ROOT_RATIO

    limit = core.ColliderResonanceLimit(
        process_id=_PID,
        resonance="G_KK^(1)",
        final_state="gamma gamma",
        limit_kind=core.MASS_LOWER_BOUND,
        value=float(active["value"]),
        units=active["units"],
    )
    prediction = core.ColliderResonancePrediction(
        resonance="G_KK^(1)",
        final_state="gamma gamma",
        mass_tev=expected_mass_tev,
    )
    direct = core.evaluate_resonance_limit(prediction, limit)

    assert result.predicted == pytest.approx(expected_mass_tev)
    assert result.predicted == pytest.approx(direct.predicted_mass_tev)
    assert result.experimental == pytest.approx(float(active["value"]))
    assert result.ratio == pytest.approx(float(active["value"]) / expected_mass_tev)
    assert result.ratio == pytest.approx(direct.ratio_to_budget)
    assert result.budget == pytest.approx(direct.budget)
    assert result.passes is True
    assert result.diagnostics["benchmark_ktilde"] == pytest.approx(0.1)


def test_graviton_spectrum_mapping_sources_use_existing_declared_extras():
    ew_point = point_builder.make_point(kk_ew_mass_gev=3000.0)
    gluon_point = point_builder.make_point(kk_gluon_mass_gev=3200.0)
    couplings_point = point_builder.make_point(
        quark_mass_basis_couplings=SimpleNamespace(M_KK=1400.0)
    )
    constraint = fcc.get(_PID)

    ew_result = constraint.evaluate(ew_point)
    gluon_result = constraint.evaluate(gluon_point)
    couplings_result = constraint.evaluate(couplings_point)

    assert ew_result.predicted == pytest.approx(
        3.0 * _GRAVITON_TO_GAUGE_ROOT_RATIO
    )
    assert ew_result.diagnostics["mass_source"] == "kk_ew_mass_gev"
    assert ew_result.diagnostics["lambda_ir_gev"] == pytest.approx(
        3000.0 / GAUGE_KK_ROOT_NN
    )
    assert gluon_result.predicted == pytest.approx(
        3.2 * _GRAVITON_TO_GAUGE_ROOT_RATIO
    )
    assert gluon_result.diagnostics["mass_source"] == "kk_gluon_mass_gev"
    assert couplings_result.predicted == pytest.approx(
        1.4 * SPIN2_GRAVITON_KK_ROOT
    )
    assert couplings_result.diagnostics["mass_source"] == (
        "quark_mass_basis_couplings.M_KK"
    )
    assert couplings_result.diagnostics["lambda_ir_gev"] == pytest.approx(1400.0)
    assert couplings_result.diagnostics["mass_proxy"] == (
        "m_GKK = SPIN2_GRAVITON_KK_ROOT * Lambda_IR"
    )
    assert couplings_result.diagnostics["spectrum_mapping_status"] == (
        "NEEDS-HUMAN-PHYSICS"
    )
    assert couplings_result.diagnostics["spin0_resonance_recast_status"].startswith(
        "NEEDS-HUMAN-PHYSICS"
    )


def test_safe_point_passes_and_excluded_point_fails():
    constraint = fcc.get(_PID)
    limit_lambda_ir_gev = constraint.anchor.budget * 1000.0 / SPIN2_GRAVITON_KK_ROOT

    safe = constraint.evaluate(_point_with_mkk(limit_lambda_ir_gev + 500.0))
    excluded = constraint.evaluate(_point_with_mkk(limit_lambda_ir_gev - 100.0))

    assert safe.passes is True
    assert safe.predicted == pytest.approx(
        (limit_lambda_ir_gev + 500.0) * SPIN2_GRAVITON_KK_ROOT / 1000.0
    )
    assert safe.ratio < 1.0
    assert excluded.passes is False
    assert excluded.predicted == pytest.approx(
        (limit_lambda_ir_gev - 100.0) * SPIN2_GRAVITON_KK_ROOT / 1000.0
    )
    assert excluded.ratio > 1.0
    assert excluded.experimental == pytest.approx(constraint.anchor.budget)
    assert excluded.diagnostics["active_value_id"] == _ACTIVE_VALUE_ID
    assert excluded.diagnostics["mass_proxy"] == (
        "m_GKK = SPIN2_GRAVITON_KK_ROOT * Lambda_IR"
    )
    assert excluded.diagnostics["lambda_ir_source"] == (
        "quark_mass_basis_couplings.M_KK"
    )
    assert "NEEDS-HUMAN-PHYSICS" in excluded.diagnostics["needs_human_physics"]
    assert excluded.diagnostics["sigma_times_br_recast_status"] == (
        "NEEDS-HUMAN-PHYSICS"
    )
    assert excluded.diagnostics["branching_surface_recast_status"] == (
        "NEEDS-HUMAN-PHYSICS"
    )
    assert excluded.diagnostics["all_diphoton_mass_limits_tev"][
        _CMS_2024_RANGE_ID
    ] == pytest.approx(_entry(_CMS_2024_RANGE_ID)["value"])


def test_invalid_graviton_mass_proxy_is_non_crashing_failure():
    point = point_builder.make_point(kk_ew_mass_gev=-1.0)
    result = fcc.get(_PID).evaluate(point)

    assert result.passes is False
    assert result.predicted is None
    assert result.ratio is None
    assert result.experimental == pytest.approx(fcc.get(_PID).anchor.value)
    assert result.diagnostics["exception_type"] == "ValueError"
    assert result.diagnostics["benchmark_ktilde"] == pytest.approx(0.1)


def test_evaluate_is_pure_and_deterministic():
    point = point_builder.make_point(kk_ew_mass_gev=3600.0)
    before = dict(point.extras)
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    assert point.extras == before
    for value in (first.predicted, first.experimental, first.ratio, first.budget):
        assert isinstance(value, float)
        assert math.isfinite(value)
