"""Production tests for EW003 (semileptonic CKM tensions)."""

from __future__ import annotations

import math
from pathlib import Path

import pytest
import yaml

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import point_builder
from flavor_catalog_constraints.anchors import AnchorError
from flavor_catalog_constraints.base import ConstraintProtocol, Severity
from flavor_catalog_constraints.primary.top_higgs_ew import EW003 as ew003_module
from tests.constraints.charged_current_phase5b_helpers import (
    charged_with_epsilon,
    universal_charged_point,
)

_PID = "EW003"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "top_higgs_ew" / "EW003.yaml"


def _yaml():
    with open(_SIDECAR) as handle:
        return yaml.safe_load(handle)


def _entry(observable: str):
    for item in _yaml()["pdg_or_equivalent"]:
        if item["observable"] == observable:
            return item
    raise AssertionError(f"missing observable {observable}")


def _uncertainty_pair(entry) -> tuple[float, float]:
    uncertainty = entry.get("uncertainty")
    if uncertainty is not None:
        sigma = float(uncertainty)
        return sigma, sigma
    components = entry["uncertainty_components"]
    upper_terms = []
    lower_terms = []
    for key, value in components.items():
        label = str(key).lower()
        number = float(value)
        if "plus" in label:
            upper_terms.append(number)
        elif "minus" in label:
            lower_terms.append(number)
        else:
            upper_terms.append(number)
            lower_terms.append(number)
    upper = math.sqrt(sum(value * value for value in upper_terms))
    lower = math.sqrt(sum(value * value for value in lower_terms))
    return upper, lower


def _manual_pull(inclusive_observable: str, exclusive_observable: str) -> float:
    inclusive = _entry(inclusive_observable)
    exclusive = _entry(exclusive_observable)
    inc_upper, inc_lower = _uncertainty_pair(inclusive)
    exc_upper, exc_lower = _uncertainty_pair(exclusive)
    difference = float(inclusive["value"]) - float(exclusive["value"])
    if difference >= 0.0:
        inc_sigma = inc_lower
        exc_sigma = exc_upper
    else:
        inc_sigma = inc_upper
        exc_sigma = exc_lower
    return abs(difference) / math.sqrt(inc_sigma * inc_sigma + exc_sigma * exc_sigma)


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.SOFT
    assert constraint.family == "top_higgs_ew"
    assert constraint.observable == "|V_cb| and |V_ub| inclusive-exclusive tension"


def test_anchor_matches_yaml_and_uncertainty_components():
    constraint = fcc.get(_PID)
    vcb_inc = _entry("PDG 2024 |V_cb| inclusive")
    vcb_exc = _entry("PDG 2024 |V_cb| exclusive")
    vub_inc = _entry("PDG 2024 |V_ub| inclusive")
    vub_exc = _entry("PDG 2024 |V_ub| exclusive")
    budget = _entry("PDG 2024 |V_cb| inclusive-exclusive marginal consistency")
    vub_inc_up, vub_inc_down = _uncertainty_pair(vub_inc)
    vub_exc_up, vub_exc_down = _uncertainty_pair(vub_exc)

    assert constraint.anchor.pdg_vcb_inclusive.value == pytest.approx(vcb_inc["value"])
    assert constraint.anchor.pdg_vcb_inclusive.uncertainty == pytest.approx(
        vcb_inc["uncertainty"]
    )
    assert constraint.anchor.pdg_vcb_inclusive.source_url == vcb_inc["source_url"]
    assert constraint.anchor.pdg_vcb_exclusive.value == pytest.approx(vcb_exc["value"])
    assert constraint.anchor.pdg_vcb_exclusive.uncertainty == pytest.approx(
        vcb_exc["uncertainty"]
    )
    assert constraint.anchor.pdg_vub_inclusive.value == pytest.approx(vub_inc["value"])
    assert constraint.anchor.pdg_vub_inclusive.uncertainty_upper == pytest.approx(
        vub_inc_up
    )
    assert constraint.anchor.pdg_vub_inclusive.uncertainty_lower == pytest.approx(
        vub_inc_down
    )
    assert constraint.anchor.pdg_vub_exclusive.uncertainty_upper == pytest.approx(
        vub_exc_up
    )
    assert constraint.anchor.pdg_vub_exclusive.uncertainty_lower == pytest.approx(
        vub_exc_down
    )
    assert constraint.anchor.pdg_vcb_consistency_sigma.value == pytest.approx(
        budget["value"]
    )
    assert constraint.anchor.pdg_vcb_consistency_sigma.source_url == budget["source_url"]
    assert constraint.anchor.budget == pytest.approx(3.0)


def test_ew003_list_anchors_route_through_scaffold_load_anchor(monkeypatch):
    calls: list[tuple[str, ...]] = []
    original_load_anchor = ew003_module.load_anchor

    def spy_load_anchor(*args, **kwargs):
        calls.append(tuple(kwargs["candidates"]))
        return original_load_anchor(*args, **kwargs)

    monkeypatch.setattr(ew003_module, "load_anchor", spy_load_anchor)
    anchor = ew003_module._load_ew003_anchor(_PID)

    assert calls == [
        ("pdg_or_equivalent[0]",),
        ("pdg_or_equivalent[1]",),
        ("pdg_or_equivalent[2]",),
        ("pdg_or_equivalent[3]",),
        ("pdg_or_equivalent[4]",),
        ("pdg_or_equivalent[5]",),
        ("pdg_or_equivalent[11]",),
        ("pdg_or_equivalent[6]",),
        ("pdg_or_equivalent[7]",),
        ("pdg_or_equivalent[8]",),
        ("pdg_or_equivalent[9]",),
        ("pdg_or_equivalent[10]",),
    ]
    assert anchor.pdg_vcb_inclusive.value == pytest.approx(
        fcc.get(_PID).anchor.pdg_vcb_inclusive.value
    )


def test_anchor_loader_fails_loudly_on_missing_observable():
    with pytest.raises(AnchorError):
        ew003_module._load_value_anchor(
            "not an EW003 observable",
            process_id=_PID,
            expected_units="10^-3",
        )


def test_numerical_pulls_match_independent_yaml_recomputation():
    constraint = fcc.get(_PID)
    result = constraint.evaluate(point_builder.empty_point())
    vcb_pull = _manual_pull(
        "PDG 2024 |V_cb| inclusive",
        "PDG 2024 |V_cb| exclusive",
    )
    vub_pull = _manual_pull(
        "PDG 2024 |V_ub| inclusive",
        "PDG 2024 |V_ub| exclusive",
    )
    max_pull = max(vcb_pull, vub_pull)

    assert vcb_pull == pytest.approx(3.0728851183895105)
    assert vub_pull == pytest.approx(1.4270051199399387)
    assert result.predicted == pytest.approx(max_pull)
    assert result.experimental == pytest.approx(3.0)
    assert result.budget == pytest.approx(3.0)
    assert result.ratio == pytest.approx(max_pull / 3.0)
    assert result.sm_prediction == pytest.approx(0.0)


def test_absent_charged_current_path_keeps_data_level_pull_non_crashing():
    constraint = fcc.get(_PID)
    result = constraint.evaluate(point_builder.empty_point())

    assert result.predicted == pytest.approx(3.0728851183895105)
    assert result.experimental == pytest.approx(3.0)
    assert result.ratio == pytest.approx(3.0728851183895105 / 3.0)
    assert result.diagnostics["charged_current_diagnostics_available"] is False
    assert (
        result.diagnostics["charged_current_missing_optional_extra"]
        == "rs_charged_current"
    )
    assert result.diagnostics["parameter_point_used_for_scalar_pull"] is False


def test_charged_current_epsilon_diagnostics_preserve_scalar_pull():
    constraint = fcc.get(_PID)
    baseline = constraint.evaluate(point_builder.empty_point())
    result = constraint.evaluate(universal_charged_point())
    diagnostics = result.diagnostics["charged_current_diagnostics"]

    assert result.predicted == pytest.approx(baseline.predicted)
    assert result.experimental == pytest.approx(baseline.experimental)
    assert result.ratio == pytest.approx(baseline.ratio)
    assert result.passes is baseline.passes
    assert result.diagnostics["charged_current_diagnostics_available"] is True
    assert isinstance(diagnostics["epsilon_cb_light_average_e_mu"], complex)
    assert isinstance(diagnostics["epsilon_ub_light_average_e_mu"], complex)
    assert diagnostics["epsilon_cb"]["up_flavor"] == "c"
    assert diagnostics["epsilon_cb"]["down_flavor"] == "b"
    assert diagnostics["epsilon_ub"]["up_flavor"] == "u"
    assert diagnostics["epsilon_ub"]["down_flavor"] == "b"
    assert "m_wprime_gev" in diagnostics


def test_universal_charged_current_shift_cancels_in_inclusive_exclusive_ratio():
    constraint = fcc.get(_PID)
    baseline = constraint.evaluate(point_builder.empty_point())
    updates = {}
    for lepton_index in range(3):
        updates[(1, 2, lepton_index)] = 0.2 + 0.0j
        updates[(0, 2, lepton_index)] = 0.2 + 0.0j
    charged = charged_with_epsilon(
        universal_charged_point().extras["rs_charged_current"],
        updates,
    )
    result = constraint.evaluate(point_builder.make_point(rs_charged_current=charged))
    diagnostics = result.diagnostics["charged_current_diagnostics"]

    assert diagnostics["epsilon_cb_light_average_e_mu"] == pytest.approx(0.2 + 0.0j)
    assert diagnostics["epsilon_ub_light_average_e_mu"] == pytest.approx(0.2 + 0.0j)
    assert diagnostics["epsilon_cb_abs_1_plus"] == pytest.approx(1.2)
    assert diagnostics["epsilon_ub_abs_1_plus"] == pytest.approx(1.2)
    assert diagnostics["universal_cc_ratio_multiplier_vcb"] == pytest.approx(1.0)
    assert diagnostics["universal_cc_ratio_multiplier_vub"] == pytest.approx(1.0)
    assert diagnostics["universal_cc_pull_multiplier_vcb"] == pytest.approx(1.0)
    assert diagnostics["universal_cc_pull_multiplier_vub"] == pytest.approx(1.0)
    assert diagnostics["universal_cc_pull_shift_sigma_vcb"] == pytest.approx(0.0)
    assert diagnostics["universal_cc_pull_shift_sigma_vub"] == pytest.approx(0.0)
    assert diagnostics["universal_cc_cancellation_documented"] is True
    assert result.predicted == pytest.approx(baseline.predicted)
    assert result.ratio == pytest.approx(baseline.ratio)


def test_covariance_scheme_needs_human_note_is_present():
    result = fcc.get(_PID).evaluate(point_builder.empty_point())

    assert result.diagnostics["matching_coverage"] == "PARTIAL"
    assert result.diagnostics["ew003_covariance_scheme_input_supplied"] is False
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]
    assert "covariance/scheme" in result.diagnostics["needs_human_physics"]
    assert "diagnostic only" in result.notes


def test_safe_subobservable_passes_and_tension_subobservable_fails_soft_budget():
    result = fcc.get(_PID).evaluate(point_builder.empty_point())
    pulls = result.diagnostics["primary_pulls"]
    vcb = pulls["PDG 2024 |V_cb|"]
    vub = pulls["PDG 2024 |V_ub|"]

    assert vub["passes_budget"] is True
    assert vub["ratio_to_budget"] < 1.0
    assert vcb["passes_budget"] is False
    assert vcb["ratio_to_budget"] > 1.0
    assert result.passes is False
    assert result.severity is Severity.SOFT
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]


def test_evaluate_has_real_finite_fields_and_no_parameter_point_dependence():
    point = point_builder.make_point(kk_ew_mass_gev=6000.0)
    result = fcc.get(_PID).evaluate(point)

    for value in (
        result.predicted,
        result.ratio,
        result.budget,
        result.sm_prediction,
        result.experimental,
    ):
        assert isinstance(value, float)
        assert math.isfinite(value)
    assert result.diagnostics["parameter_point_used"] is False
    assert result.diagnostics["qcd_running_applied"] is False


def test_evaluate_is_pure_and_deterministic():
    constraint = fcc.get(_PID)
    empty = point_builder.empty_point()
    with_extra = point_builder.make_point(kk_ew_mass_gev=6000.0)

    first = constraint.evaluate(empty)
    second = constraint.evaluate(empty)
    third = constraint.evaluate(with_extra)

    assert first == second
    assert first == third
