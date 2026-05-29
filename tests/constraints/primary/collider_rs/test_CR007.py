"""Production tests for CR007 (spin-2 KK-graviton resonance)."""

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
from flavor_catalog_constraints.primary.collider_rs import CR007 as cr007_module
from quarkConstraints.scales import GAUGE_KK_ROOT_NN, SPIN2_GRAVITON_KK_ROOT
from quarkConstraints import collider_resonance as core

_PID = "CR007"
_ACTIVE_VALUE_ID = "PDG2025:CR007:bulk_graviton_diboson_mass_kMpl_0p5"
_DIPHOTON_VALUE_ID = "PDG2025:CR007:RS_graviton_diphoton_mass_kMpl_0p1"
_DILEPTON_VALUE_ID = "PDG2025:CR007:RS_graviton_dilepton_mass_kMpl_0p1"
_ATLAS_BULK_VALUE_ID = "ATLAS2018:CR007:bulk_graviton_combined_mass_kMpl_1"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "collider_rs" / "CR007.yaml"
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


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.HARD
    assert constraint.family == "collider_rs"
    assert constraint.observable == "m(G_KK^(1) -> WW, ZZ)"


def test_anchor_matches_yaml_and_loud_fail_probe():
    constraint = fcc.get(_PID)
    active = _entry(_ACTIVE_VALUE_ID)
    diphoton = _entry(_DIPHOTON_VALUE_ID)
    dilepton = _entry(_DILEPTON_VALUE_ID)
    atlas_bulk = _entry(_ATLAS_BULK_VALUE_ID)

    assert constraint.anchor.active_limit.value_tev == pytest.approx(active["value"])
    assert constraint.anchor.active_limit.units == active["units"]
    assert constraint.anchor.active_limit.source_url == active["source_url"]
    assert constraint.anchor.active_limit.display == active["display"]
    assert constraint.anchor.active_limit.experiment == active["experiment"]
    assert constraint.anchor.active_limit.channel == active["channel"]
    assert constraint.anchor.active_limit.model_assumptions == active["model_assumptions"]
    assert constraint.anchor.generic_diphoton.value_tev == pytest.approx(
        diphoton["value"]
    )
    assert constraint.anchor.generic_diphoton.channel == diphoton["channel"]
    assert constraint.anchor.generic_dilepton.value_tev == pytest.approx(
        dilepton["value"]
    )
    assert constraint.anchor.atlas_bulk_combined.value_tev == pytest.approx(
        atlas_bulk["value"]
    )
    assert constraint.anchor.budget == pytest.approx(active["value"])

    with pytest.raises(AnchorError):
        cr007_module._load_value_anchor("not-a-real-value-id", process_id=_PID)


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


def test_mass_limit_numerics_validate_against_core_not_adapter():
    constraint = fcc.get(_PID)
    point = point_builder.make_point(kk_ew_mass_gev=2000.0)
    result = constraint.evaluate(point)
    active = _entry(_ACTIVE_VALUE_ID)
    expected_mass_tev = 2.0 * _GRAVITON_TO_GAUGE_ROOT_RATIO

    limit = core.ColliderResonanceLimit(
        process_id=_PID,
        resonance="G_KK^(1)",
        final_state="WW, ZZ",
        limit_kind=core.MASS_LOWER_BOUND,
        value=float(active["value"]),
        units=active["units"],
    )
    prediction = core.ColliderResonancePrediction(
        resonance="G_KK^(1)",
        final_state="WW, ZZ",
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


def test_graviton_spectrum_mapping_sources_use_existing_declared_extras():
    ew_point = point_builder.make_point(kk_ew_mass_gev=2300.0)
    gluon_point = point_builder.make_point(kk_gluon_mass_gev=2400.0)
    couplings_point = point_builder.make_point(
        quark_mass_basis_couplings=SimpleNamespace(M_KK=2500.0)
    )
    constraint = fcc.get(_PID)

    ew_result = constraint.evaluate(ew_point)
    gluon_result = constraint.evaluate(gluon_point)
    couplings_result = constraint.evaluate(couplings_point)

    assert ew_result.predicted == pytest.approx(
        2.3 * _GRAVITON_TO_GAUGE_ROOT_RATIO
    )
    assert ew_result.diagnostics["mass_source"] == "kk_ew_mass_gev"
    assert ew_result.diagnostics["lambda_ir_gev"] == pytest.approx(
        2300.0 / GAUGE_KK_ROOT_NN
    )
    assert gluon_result.predicted == pytest.approx(
        2.4 * _GRAVITON_TO_GAUGE_ROOT_RATIO
    )
    assert gluon_result.diagnostics["mass_source"] == "kk_gluon_mass_gev"
    assert couplings_result.predicted == pytest.approx(
        2.5 * SPIN2_GRAVITON_KK_ROOT
    )
    assert couplings_result.diagnostics["mass_source"] == (
        "quark_mass_basis_couplings.M_KK"
    )
    assert couplings_result.diagnostics["lambda_ir_gev"] == pytest.approx(2500.0)
    assert couplings_result.diagnostics["graviton_spin2_root"] == pytest.approx(
        SPIN2_GRAVITON_KK_ROOT
    )
    assert couplings_result.diagnostics["mass_proxy"] == (
        "m_GKK = SPIN2_GRAVITON_KK_ROOT * Lambda_IR"
    )
    assert couplings_result.diagnostics["spectrum_mapping_status"] == (
        "NEEDS-HUMAN-PHYSICS"
    )


def test_point_builder_couplings_mkk_is_treated_as_lambda_ir():
    point = point_builder.build_from_quark_couplings(SimpleNamespace(M_KK=3000.0))
    result = fcc.get(_PID).evaluate(point)

    assert result.predicted == pytest.approx(3.0 * SPIN2_GRAVITON_KK_ROOT)
    assert result.diagnostics["mass_source"] == "quark_mass_basis_couplings.M_KK"
    assert result.diagnostics["lambda_ir_source"] == "quark_mass_basis_couplings.M_KK"
    assert result.diagnostics["lambda_ir_gev"] == pytest.approx(3000.0)


def test_safe_point_passes_and_excluded_point_fails():
    constraint = fcc.get(_PID)
    limit_lambda_ir_gev = constraint.anchor.budget * 1000.0 / SPIN2_GRAVITON_KK_ROOT

    safe = constraint.evaluate(
        point_builder.make_point(
            quark_mass_basis_couplings=SimpleNamespace(
                M_KK=limit_lambda_ir_gev + 1000.0
            )
        )
    )
    excluded = constraint.evaluate(
        point_builder.make_point(
            quark_mass_basis_couplings=SimpleNamespace(
                M_KK=limit_lambda_ir_gev - 100.0
            )
        )
    )

    assert safe.passes is True
    assert safe.predicted == pytest.approx(
        (limit_lambda_ir_gev + 1000.0) * SPIN2_GRAVITON_KK_ROOT / 1000.0
    )
    assert safe.ratio < 1.0
    assert excluded.passes is False
    assert excluded.predicted == pytest.approx(
        (limit_lambda_ir_gev - 100.0) * SPIN2_GRAVITON_KK_ROOT / 1000.0
    )
    assert excluded.ratio > 1.0
    assert excluded.experimental == pytest.approx(constraint.anchor.budget)
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
    assert excluded.diagnostics["generic_rs_mass_limits_tev"][
        _DIPHOTON_VALUE_ID
    ] == pytest.approx(_entry(_DIPHOTON_VALUE_ID)["value"])


def test_invalid_graviton_mass_proxy_is_non_crashing_failure():
    point = point_builder.make_point(kk_ew_mass_gev=-1.0)
    result = fcc.get(_PID).evaluate(point)

    assert result.passes is False
    assert result.predicted is None
    assert result.ratio is None
    assert result.experimental == pytest.approx(fcc.get(_PID).anchor.value)
    assert result.diagnostics["exception_type"] == "ValueError"


def test_evaluate_is_pure_and_deterministic():
    point = point_builder.make_point(kk_ew_mass_gev=2600.0)
    before = dict(point.extras)
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    assert point.extras == before
    for value in (first.predicted, first.experimental, first.ratio, first.budget):
        assert isinstance(value, float)
        assert math.isfinite(value)
