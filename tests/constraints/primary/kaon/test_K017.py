"""Production tests for K017 (charged-kaon leptonic LFU ratio R_K)."""

from __future__ import annotations

from dataclasses import replace
import math
from pathlib import Path

import pytest
import yaml

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import anchors, point_builder
from flavor_catalog_constraints.base import ConstraintProtocol, Severity
from flavor_catalog_constraints.primary.kaon import K017 as k017_module
from quarkConstraints.leptonic_tree import (
    evaluate_leptonic_lfu_ratio,
    kplus_enu_over_munu_inputs_from_sm_ratio_anchor,
)

_PID = "K017"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "kaon" / "K017.yaml"


def _yaml_pdg_block():
    with open(_SIDECAR) as handle:
        return yaml.safe_load(handle)["pdg_or_equivalent"]


def _scaled(block) -> float:
    return float(block["value"]) * float(block.get("scale", 1.0))


def _scaled_uncertainty(block) -> float:
    return float(block["uncertainty"]) * float(block.get("scale", 1.0))


def _core_inputs_from_yaml_sm_anchor():
    sm = _yaml_pdg_block()["sm_prediction"]
    return kplus_enu_over_munu_inputs_from_sm_ratio_anchor(
        sm_ratio=_scaled(sm),
        constants_citation=f"K017.yaml SM anchor: {sm['source']}",
    )


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.HARD
    assert constraint.family == "kaon"
    assert constraint.observable == "R_K = Gamma(K -> e nu) / Gamma(K -> mu nu)"


def test_anchor_matches_yaml_and_loud_fail_probe():
    constraint = fcc.get(_PID)
    pdg = _yaml_pdg_block()
    exp = pdg["canonical_ratio"]
    sm = pdg["sm_prediction"]
    na62 = pdg["dominant_experimental_input"]
    kloe = pdg["supporting_experimental_input_kloe"]
    combined = math.sqrt(_scaled_uncertainty(exp) ** 2 + _scaled_uncertainty(sm) ** 2)
    central = abs(_scaled(exp) - _scaled(sm))
    na62_scale = float(na62["scale"])
    kloe_scale = float(kloe["scale"])

    assert constraint.anchor.experimental.value == pytest.approx(_scaled(exp))
    assert constraint.anchor.experimental.uncertainty == pytest.approx(
        _scaled_uncertainty(exp)
    )
    assert constraint.anchor.experimental.source_url == exp["source_url"]
    assert constraint.anchor.standard_model.value == pytest.approx(_scaled(sm))
    assert constraint.anchor.standard_model.uncertainty == pytest.approx(
        _scaled_uncertainty(sm)
    )
    assert constraint.anchor.standard_model.source_url == sm["source_url"]
    assert constraint.anchor.dominant_experimental_input.value == pytest.approx(
        float(na62["value"]) * na62_scale
    )
    assert constraint.anchor.dominant_experimental_input.combined_uncertainty == pytest.approx(
        math.sqrt(float(na62["stat_uncertainty"]) ** 2 + float(na62["syst_uncertainty"]) ** 2)
        * na62_scale
    )
    assert constraint.anchor.supporting_experimental_input.value == pytest.approx(
        float(kloe["value"]) * kloe_scale
    )
    assert constraint.anchor.supporting_experimental_input.combined_uncertainty == pytest.approx(
        math.sqrt(float(kloe["stat_uncertainty"]) ** 2 + float(kloe["syst_uncertainty"]) ** 2)
        * kloe_scale
    )
    assert constraint.anchor.budget_band.central_residual == pytest.approx(central)
    assert constraint.anchor.budget_band.combined_sigma == pytest.approx(combined)
    assert constraint.anchor.budget == pytest.approx(central + combined)
    assert constraint.anchor.budget == pytest.approx(2.005538513813748e-07)

    with pytest.raises(anchors.AnchorError):
        anchors.load_anchor(
            _PID,
            family="kaon",
            candidates=("no_such_block",),
        )


def test_ratio_anchors_route_through_scaffold_load_anchor(monkeypatch):
    calls: list[tuple[str, ...]] = []
    original_load_anchor = k017_module.load_anchor

    def spy_load_anchor(*args, **kwargs):
        calls.append(tuple(kwargs["candidates"]))
        return original_load_anchor(*args, **kwargs)

    monkeypatch.setattr(k017_module, "load_anchor", spy_load_anchor)
    anchor = k017_module._load_k017_anchor(_PID)

    assert ("canonical_ratio",) in calls
    assert ("sm_prediction",) in calls
    assert ("dominant_experimental_input",) in calls
    assert ("supporting_experimental_input_kloe",) in calls
    assert anchor.experimental.value == pytest.approx(2.488e-5)
    assert anchor.standard_model.value == pytest.approx(2.477e-5)

    def mismatched_load_anchor(*args, **kwargs):
        anchor = original_load_anchor(*args, **kwargs)
        return replace(anchor, block_key="wrong_block")

    monkeypatch.setattr(k017_module, "load_anchor", mismatched_load_anchor)
    with pytest.raises(anchors.AnchorError, match="load_anchor selected 'wrong_block'"):
        k017_module._load_k017_anchor(_PID)


def test_evaluate_without_kk_mass_uses_sm_and_flags_missing_proxy_input():
    constraint = fcc.get(_PID)
    result = constraint.evaluate(point_builder.empty_point())

    assert result.process_id == _PID
    assert result.passes is True
    assert result.predicted == pytest.approx(constraint.sm_result.ratio)
    assert result.sm_prediction == pytest.approx(constraint.sm_result.sm_ratio)
    assert result.experimental == pytest.approx(constraint.anchor.value)
    assert result.ratio == pytest.approx(0.0)
    assert result.budget == pytest.approx(constraint.anchor.budget)
    assert result.diagnostics["lepton_nonuniversal_proxy_evaluated"] is False
    assert result.diagnostics["missing_extra"] == "kk_ew_mass_gev"
    assert result.diagnostics["sm_anchor_ratio"] == pytest.approx(
        constraint.anchor.sm_value
    )
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]


def test_sm_limit_ratio_matches_core_evaluator_and_yaml_anchor():
    constraint = fcc.get(_PID)
    result = constraint.evaluate(point_builder.make_point(kk_ew_mass_gev=1.0e9))
    core_inputs = _core_inputs_from_yaml_sm_anchor()
    core_sm = evaluate_leptonic_lfu_ratio(inputs=core_inputs)

    assert core_sm.diagnostics["tree_ratio_without_radiation"] == pytest.approx(
        2.5689629165499624e-05
    )
    assert core_inputs.radiative_correction_multiplier == pytest.approx(
        0.9642023183917867
    )
    assert constraint.sm_result.ratio == pytest.approx(core_sm.ratio)
    assert constraint.sm_result.ratio == pytest.approx(2.477e-05)
    assert constraint.sm_result.ratio == pytest.approx(constraint.anchor.sm_value)
    assert result.sm_prediction == pytest.approx(core_sm.sm_ratio)
    assert result.diagnostics["sm_anchor_ratio"] == pytest.approx(
        constraint.anchor.sm_value
    )
    assert result.passes is True


def test_lepton_nonuniversal_proxy_matches_core_evaluator_from_yaml_anchor():
    constraint = fcc.get(_PID)
    m_kk = 20.0
    result = constraint.evaluate(point_builder.make_point(kk_ew_mass_gev=m_kk))
    core_inputs = _core_inputs_from_yaml_sm_anchor()
    expected = evaluate_leptonic_lfu_ratio(m_kk_gev=m_kk, inputs=core_inputs)

    assert result.predicted == pytest.approx(expected.ratio)
    assert result.sm_prediction == pytest.approx(expected.sm_ratio)
    assert result.diagnostics["numerator_np_amplitude_ratio"] == pytest.approx(
        expected.numerator_np_amplitude_ratio
    )
    assert result.diagnostics["denominator_np_amplitude_ratio"] == pytest.approx(
        expected.denominator_np_amplitude_ratio
    )
    assert result.diagnostics["numerator_amplitude_multiplier"] == pytest.approx(
        expected.numerator_amplitude_multiplier
    )
    assert result.diagnostics["denominator_amplitude_multiplier"] == pytest.approx(
        expected.denominator_amplitude_multiplier
    )
    assert result.ratio == pytest.approx(
        abs(expected.np_shift_ratio) / constraint.anchor.budget
    )
    assert result.diagnostics["kk_ew_mass_extra_used"] is True
    assert result.diagnostics["uses_lepton_nonuniversal_proxy"] is True
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]


def test_evaluate_runs_end_to_end_with_real_fields_and_complex_diagnostics():
    result = fcc.get(_PID).evaluate(point_builder.make_point(kk_ew_mass_gev=20.0))

    for value in (
        result.predicted,
        result.ratio,
        result.budget,
        result.sm_prediction,
        result.experimental,
    ):
        assert isinstance(value, float)
        assert math.isfinite(value)
    for key in ("numerator_np_amplitude_ratio", "denominator_np_amplitude_ratio"):
        assert isinstance(result.diagnostics[key], complex)
    for key in (
        "m_kk_gev",
        "matching_scale_gev",
        "numerator_amplitude_multiplier",
        "denominator_amplitude_multiplier",
        "np_shift_ratio",
        "budget_combined_sigma",
        "hard_veto_np_shift_budget",
        "meson_mass_gev",
        "tree_ratio_without_radiation",
        "radiative_correction_multiplier",
    ):
        assert isinstance(result.diagnostics[key], float)
        assert math.isfinite(result.diagnostics[key])


@pytest.mark.parametrize(
    ("m_kk_gev", "expected_pass"),
    [
        (3000.0, True),
        (3.0, False),
    ],
)
def test_safe_point_passes_and_large_proxy_point_fails(
    m_kk_gev: float,
    expected_pass: bool,
):
    result = fcc.get(_PID).evaluate(point_builder.make_point(kk_ew_mass_gev=m_kk_gev))

    assert result.passes is expected_pass
    if expected_pass:
        assert result.ratio < 1.0
    else:
        assert result.ratio > 1.0


def test_evaluate_is_pure_and_deterministic():
    point = point_builder.make_point(kk_ew_mass_gev=20.0)
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    assert point.get_extra("kk_ew_mass_gev") == 20.0
