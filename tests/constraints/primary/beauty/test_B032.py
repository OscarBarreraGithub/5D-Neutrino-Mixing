"""Production tests for B032 (charmless B -> pi K nonleptonic stub)."""

from __future__ import annotations

import math
from dataclasses import replace
from pathlib import Path

import pytest
import yaml

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import point_builder
from flavor_catalog_constraints.base import ConstraintProtocol, Severity
from flavor_catalog_constraints.physics_adapters.charmless_b_pik import (
    compare_b_to_pi_k_np_room_to_measurement,
)
from flavor_catalog_constraints.primary.beauty import B032 as b032_module

_PID = "B032"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "beauty" / "B032.yaml"

_CHARGED_PI0 = "A_CP(B+ -> K+ pi0)"
_NEUTRAL_KPLUS_PIMINUS = "A_CP(B0 -> K+ pi-)"
_BELLE_II_SUM_RULE = "Belle II Kpi sum-rule test"


def _yaml_pdg_block():
    with open(_SIDECAR) as handle:
        return yaml.safe_load(handle)["pdg_or_equivalent"]


def _find_entry(entries, observable):
    return next(entry for entry in entries if entry["observable"] == observable)


def _expected_delta_from_yaml() -> tuple[float, float, float, float]:
    direct = _yaml_pdg_block()["direct_cp_asymmetries"]
    charged = _find_entry(direct, _CHARGED_PI0)
    neutral = _find_entry(direct, _NEUTRAL_KPLUS_PIMINUS)
    delta_percent = float(charged["value_percent"]) - float(neutral["value_percent"])
    sigma_percent = math.sqrt(
        float(charged["uncertainty_percent"]) ** 2
        + float(neutral["uncertainty_percent"]) ** 2
    )
    return (
        float(delta_percent),
        float(sigma_percent),
        float(delta_percent / 100.0),
        float(sigma_percent / 100.0),
    )


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.INFO
    assert constraint.family == "beauty"
    assert constraint.observable == "Delta A_CP(B -> K pi)"


def test_anchor_matches_yaml_and_budget():
    constraint = fcc.get(_PID)
    pdg = _yaml_pdg_block()
    charged = _find_entry(pdg["direct_cp_asymmetries"], _CHARGED_PI0)
    neutral = _find_entry(pdg["direct_cp_asymmetries"], _NEUTRAL_KPLUS_PIMINUS)
    sum_rule = _find_entry(pdg["post_2008_measurements"], _BELLE_II_SUM_RULE)
    delta_percent, sigma_percent, delta, sigma = _expected_delta_from_yaml()

    assert len(constraint.anchor.branching_fractions) == 4
    assert len(constraint.anchor.direct_cp_asymmetries) == 4
    assert len(constraint.anchor.time_dependent_cp) == 2
    assert len(constraint.anchor.post_2008_measurements) == 2
    assert constraint.anchor.delta_acp_kpi.charged_pi0.value_percent == pytest.approx(
        charged["value_percent"]
    )
    assert (
        constraint.anchor.delta_acp_kpi.neutral_kplus_piminus.value_percent
        == pytest.approx(neutral["value_percent"])
    )
    assert constraint.anchor.delta_acp_kpi.value_percent == pytest.approx(
        delta_percent
    )
    assert constraint.anchor.delta_acp_kpi.uncertainty_percent == pytest.approx(
        sigma_percent
    )
    assert constraint.anchor.value == pytest.approx(delta)
    assert constraint.anchor.uncertainty == pytest.approx(sigma)
    assert constraint.anchor.budget == pytest.approx(abs(delta))
    assert constraint.anchor.budget == pytest.approx(0.1101)
    assert (
        constraint.anchor.delta_acp_kpi.charged_pi0.source_url
        == charged["source_url"]
    )
    assert (
        b032_module._find_observable_anchor(
            constraint.anchor.post_2008_measurements,
            _BELLE_II_SUM_RULE,
            process_id=_PID,
            block_name="post_2008_measurements",
        ).value
        == pytest.approx(sum_rule["value"])
    )


def test_anchor_loading_routes_list_entries_through_load_anchor(monkeypatch):
    calls: list[tuple[str, ...]] = []
    original_load_anchor = b032_module.load_anchor

    def spy_load_anchor(*args, **kwargs):
        calls.append(tuple(kwargs["candidates"]))
        return original_load_anchor(*args, **kwargs)

    monkeypatch.setattr(b032_module, "load_anchor", spy_load_anchor)
    anchor = b032_module._load_b032_anchor(_PID)

    assert calls == [
        ("direct_cp_asymmetries[0]",),
        ("direct_cp_asymmetries[1]",),
        ("direct_cp_asymmetries[2]",),
        ("direct_cp_asymmetries[3]",),
        ("branching_fractions[0]",),
        ("branching_fractions[1]",),
        ("branching_fractions[2]",),
        ("branching_fractions[3]",),
        ("time_dependent_cp[0]",),
        ("time_dependent_cp[1]",),
        ("post_2008_measurements[0]",),
        ("post_2008_measurements[1]",),
    ]
    assert anchor.value == pytest.approx(fcc.get(_PID).anchor.value)

    with pytest.raises(fcc.AnchorError, match="not a mapping"):
        fcc.load_anchor(_PID, family="beauty", candidates=("direct_cp_asymmetries",))

    with pytest.raises(fcc.AnchorError, match="has no observable"):
        b032_module._find_observable_anchor(
            anchor.direct_cp_asymmetries,
            "missing B032 observable",
            process_id=_PID,
            block_name="direct_cp_asymmetries",
        )


def test_anchor_loading_rejects_mismatched_load_anchor_block_key(monkeypatch):
    original_load_anchor = b032_module.load_anchor

    def mismatched_load_anchor(*args, **kwargs):
        anchor = original_load_anchor(*args, **kwargs)
        return replace(anchor, block_key="wrong_block")

    monkeypatch.setattr(b032_module, "load_anchor", mismatched_load_anchor)
    with pytest.raises(fcc.AnchorError, match="load_anchor selected 'wrong_block'"):
        b032_module._load_b032_anchor(_PID)


def test_stub_numeric_validation_matches_independent_yaml_recompute():
    constraint = fcc.get(_PID)
    delta_percent, sigma_percent, delta, sigma = _expected_delta_from_yaml()
    expected_budget = abs(delta)
    expected_ratio = abs(delta) / expected_budget

    result = constraint.evaluate(point_builder.empty_point())

    assert result.predicted is None
    assert result.sm_prediction is None
    assert result.experimental == pytest.approx(delta)
    assert result.budget == pytest.approx(expected_budget)
    assert result.ratio == pytest.approx(expected_ratio)
    assert result.ratio == pytest.approx(1.0)
    assert result.diagnostics["delta_acp_value_percent"] == pytest.approx(
        delta_percent
    )
    assert result.diagnostics["delta_acp_uncertainty_percent"] == pytest.approx(
        sigma_percent
    )
    assert result.diagnostics["delta_acp_uncertainty"] == pytest.approx(sigma)
    assert result.diagnostics["no_hadronic_amplitude"] is True
    assert result.diagnostics["no_penguin_calculation"] is True


def test_adapter_room_comparison_has_pass_fail_behavior_without_amplitudes():
    measurement = 0.1101
    uncertainty = 0.012394756956825766

    safe = compare_b_to_pi_k_np_room_to_measurement(
        measured_observable=measurement,
        experimental_uncertainty=uncertainty,
        documented_np_room_abs=measurement,
    )
    overfilled = compare_b_to_pi_k_np_room_to_measurement(
        measured_observable=measurement,
        experimental_uncertainty=uncertainty,
        documented_np_room_abs=measurement / 2.0,
    )

    assert safe.passes is True
    assert safe.ratio_to_room == pytest.approx(1.0)
    assert overfilled.passes is False
    assert overfilled.ratio_to_room == pytest.approx(2.0)


def test_evaluate_is_info_non_vetoing_and_uses_no_point_inputs():
    constraint = fcc.get(_PID)
    empty_result = constraint.evaluate(point_builder.empty_point())
    irrelevant_point = point_builder.make_point(
        kk_gluon_mass_gev=100.0,
        kk_ew_mass_gev=200.0,
    )
    irrelevant_result = constraint.evaluate(irrelevant_point)

    assert empty_result == irrelevant_result
    assert empty_result.process_id == _PID
    assert empty_result.severity is Severity.INFO
    assert empty_result.passes is True
    assert empty_result.predicted is None
    assert empty_result.sm_prediction is None
    assert empty_result.experimental == pytest.approx(constraint.anchor.value)
    assert empty_result.budget == pytest.approx(constraint.anchor.budget)
    assert empty_result.ratio == pytest.approx(1.0)
    assert empty_result.diagnostics["non_vetoing"] is True
    assert empty_result.diagnostics["parameter_point_inputs_used"] == ()
    assert "NEEDS-HUMAN-PHYSICS" in empty_result.diagnostics["needs_human_physics_sm"]
    assert "NEEDS-HUMAN-PHYSICS" in empty_result.diagnostics["needs_human_physics_np"]


def test_evaluate_has_real_finite_numeric_fields():
    result = fcc.get(_PID).evaluate(point_builder.empty_point())

    assert result.predicted is None
    assert result.sm_prediction is None
    for value in (result.experimental, result.ratio, result.budget):
        assert isinstance(value, float)
        assert math.isfinite(value)
    for key in (
        "delta_acp_value_percent",
        "delta_acp_uncertainty_percent",
        "delta_acp_value",
        "delta_acp_uncertainty",
        "measurement_abs",
        "documented_np_room_abs",
        "measurement_to_np_room_ratio",
        "charged_pi0_acp_percent",
        "charged_pi0_acp_uncertainty_percent",
        "neutral_kplus_piminus_acp_percent",
        "neutral_kplus_piminus_acp_uncertainty_percent",
        "belle_ii_kpi_sum_rule_value",
    ):
        assert isinstance(result.diagnostics[key], float)
        assert math.isfinite(result.diagnostics[key])


def test_evaluate_is_pure_and_deterministic():
    point = point_builder.make_point(kk_gluon_mass_gev=3000.0)
    before_extras = dict(point.extras)
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    assert dict(point.extras) == before_extras
