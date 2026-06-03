"""Production tests for K020 (K+ -> pi+ mu+ e- LFV)."""

from __future__ import annotations

import math
from pathlib import Path

import pytest
import yaml

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import point_builder
from flavor_catalog_constraints.base import ConstraintLevel, ConstraintProtocol, Severity
from flavor_catalog_constraints.secondary.kaon import K020 as k020_module
from tests.constraints.lfv_rare_phase4c_helpers import (
    diagonal_lfv_rare_point,
    lfv_coeff,
    lfv_live_rare_point,
    manual_kaon_ktopi_rate,
    scaled_lfv_rare_point,
)

_PID = "K020"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = (
    _REPO_ROOT
    / "flavor_catalog"
    / "processes"
    / "secondary"
    / "kaon"
    / "K020.yaml"
)


def _yaml():
    with open(_SIDECAR) as handle:
        return yaml.safe_load(handle)


def _value_entry(value_id: str):
    for entry in _yaml()["pdg_or_equivalent"]["values"]:
        if entry["value_id"] == value_id:
            return entry
    raise AssertionError(f"no K020 value entry for {value_id!r}")


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.HARD
    assert constraint.level is ConstraintLevel.SECONDARY
    assert constraint.family == "kaon"
    assert constraint.observable == "BR(K+ -> pi+ mu+ e-) LFV"


def test_anchor_matches_yaml_and_routes_through_scaffold_load_anchor(monkeypatch):
    calls: list[tuple[str, ...]] = []
    original_load_anchor = k020_module.load_anchor

    def spy_load_anchor(*args, **kwargs):
        calls.append(tuple(kwargs["candidates"]))
        return original_load_anchor(*args, **kwargs)

    monkeypatch.setattr(k020_module, "load_anchor", spy_load_anchor)
    anchor = k020_module._load_k020_anchor(_PID)
    current = _value_entry(k020_module._CURRENT_LIMIT_VALUE_ID)
    sher = _value_entry(k020_module._SHER_E865_VALUE_ID)
    opposite = _value_entry(k020_module._OPPOSITE_CHARGE_VALUE_ID)

    assert calls == [
        ("pdg_or_equivalent.values[0]",),
        ("pdg_or_equivalent.values[1]",),
        ("pdg_or_equivalent.values[2]",),
    ]
    assert anchor.current_limit.value == pytest.approx(float(current["value"]))
    assert anchor.current_limit.observable == current["observable"]
    assert anchor.current_limit.limit_type == current["limit_type"]
    assert anchor.current_limit.confidence_level == current["cl"]
    assert anchor.current_limit.units == current["units"]
    assert anchor.current_limit.source_url == current["source_url"]
    assert anchor.sher_e865_limit.value == pytest.approx(float(sher["value"]))
    assert anchor.sher_e865_limit.source_url == sher["source_url"]
    assert anchor.opposite_charge_limit.value == pytest.approx(float(opposite["value"]))
    assert anchor.opposite_charge_limit.source_url == opposite["source_url"]
    assert anchor.value == pytest.approx(1.3e-11)
    assert anchor.budget == pytest.approx(anchor.value)


def test_k020_anchor_loud_fail_probe(monkeypatch):
    broken = _yaml()
    first = dict(broken["pdg_or_equivalent"]["values"][0])
    first["cl"] = "95% CL"
    broken["pdg_or_equivalent"]["values"][0] = first

    def fake_load_full_yaml(*_args, **_kwargs):
        return broken

    monkeypatch.setattr(k020_module, "load_full_yaml", fake_load_full_yaml)
    with pytest.raises(k020_module.AnchorError):
        k020_module._load_k020_anchor(_PID)


def test_absent_rs_semileptonic_wilsons_degrades_gracefully():
    constraint = fcc.get(_PID)
    result = constraint.evaluate(point_builder.empty_point())

    assert result.passes is True
    assert result.predicted is None
    assert result.ratio is None
    assert result.sm_prediction == pytest.approx(0.0)
    assert result.experimental == pytest.approx(constraint.anchor.value)
    assert result.budget == pytest.approx(constraint.anchor.budget)
    assert result.diagnostics["evaluated"] is False
    assert result.diagnostics["missing_extra"] == "rs_semileptonic_wilsons"


def test_invalid_rs_semileptonic_wilsons_is_unevaluated_not_real_pass():
    result = fcc.get(_PID).evaluate(
        point_builder.make_point(rs_semileptonic_wilsons=object())
    )

    assert result.passes is True
    assert result.predicted is None
    assert result.ratio is None
    assert result.diagnostics["evaluated"] is False
    assert result.diagnostics["invalid_extra"] == "rs_semileptonic_wilsons"
    assert result.diagnostics["exception_type"] == "ValueError"


def test_diagonal_rs_ew_builder_gives_zero_tree_lfv_and_resolves_proxy_flag():
    result = fcc.get(_PID).evaluate(diagonal_lfv_rare_point())
    coeff = lfv_coeff(diagonal_lfv_rare_point(), "s_to_d")

    assert coeff.c9_lfv_np == pytest.approx(0.0j, abs=1.0e-24)
    assert result.predicted == pytest.approx(0.0)
    assert result.ratio == pytest.approx(0.0)
    assert result.passes is True
    assert result.diagnostics["evaluated"] is True
    assert result.diagnostics["tree_level_matching_status"] == (
        "rigorous_tree_light_z_lfv_llqq_from_rs_semileptonic_wilsons"
    )
    assert "diagonal charged-lepton fit" in result.diagnostics["lfv_tree_level_note"]
    assert "NEEDS-HUMAN-PHYSICS" not in str(result.diagnostics)
    assert "proxy" not in result.notes.lower()


def test_lfv_live_toy_matches_independent_rate_and_mkk_scaling():
    constraint = fcc.get(_PID)
    low_point = lfv_live_rare_point(3000.0)
    high_point = lfv_live_rare_point(6000.0)
    low = constraint.evaluate(low_point)
    high = constraint.evaluate(high_point)
    coeff = lfv_coeff(low_point, "s_to_d")
    manual = manual_kaon_ktopi_rate(
        coeff,
        constraint.sm_inputs,
        charge_mode="muplus_eminus",
    )

    assert abs(coeff.c9_lfv_np) > 0.0
    assert low.predicted == pytest.approx(manual, rel=2.0e-10, abs=1.0e-30)
    assert low.ratio == pytest.approx(low.predicted / constraint.anchor.budget)
    assert high.predicted / low.predicted == pytest.approx(
        (3000.0 / 6000.0) ** 8,
        rel=1.0e-8,
    )
    assert low.diagnostics["charge_mode"] == "muplus_eminus"
    assert low.diagnostics["quark_vector_current_uses_left_plus_right"] is True
    assert low.diagnostics["wilson_prefactor_reused"] is False
    assert low.diagnostics["second_mkk_suppression_applied"] is False


def test_amplified_lfv_live_toy_bites_limit():
    result = fcc.get(_PID).evaluate(
        scaled_lfv_rare_point(lfv_live_rare_point(), "s_to_d", 2.0e5)
    )

    assert result.passes is False
    assert result.ratio > 1.0


def test_evaluate_runs_end_to_end_with_real_fields_and_complex_diagnostics():
    result = fcc.get(_PID).evaluate(lfv_live_rare_point())

    for value in (
        result.predicted,
        result.sm_prediction,
        result.experimental,
        result.ratio,
        result.budget,
    ):
        assert isinstance(value, float)
        assert math.isfinite(value)
    for key in (
        "rs_semileptonic_lambda_ckm",
        "y_vector_lfv",
        "y_axial_lfv",
    ):
        assert isinstance(result.diagnostics[key], complex)
    assert result.diagnostics["rs_semileptonic_wilsons_present"] is True
    assert result.diagnostics["c9_c10_to_rare_kaon_effective_inputs"] is True
    assert result.diagnostics["wilson_coefficients"]


def test_evaluate_is_pure_and_deterministic():
    point = lfv_live_rare_point()
    before_bundle = point.extras["rs_semileptonic_wilsons"]
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    assert point.extras["rs_semileptonic_wilsons"] is before_bundle
