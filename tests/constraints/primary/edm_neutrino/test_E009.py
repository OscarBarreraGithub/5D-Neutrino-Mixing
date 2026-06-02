"""Production tests for E009 (Weinberg three-gluon operator stub)."""

from __future__ import annotations

from dataclasses import replace
import math
from pathlib import Path

import pytest
import yaml

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import anchors, point_builder
from flavor_catalog_constraints.base import ConstraintProtocol, Severity
from flavor_catalog_constraints.physics_adapters.weinberg_operator import (
    compare_weinberg_coefficient_to_reference_bound,
)
from flavor_catalog_constraints.primary.edm_neutrino import E009 as e009_module

_PID = "E009"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "edm_neutrino" / "E009.yaml"
_X1E_MINUS_26_E_CM = 1.0e-26
_MEV_TO_GEV = 1.0e-3
_PR_RESPONSE_VALUE_ID = "PospelovRitz2005:E009:weinberg_normalization"
_PR_BOUND_VALUE_ID = "PDG2026-PospelovRitz2005:E009:w_bound"
_HH_RESPONSE_VALUE_ID = "HaischHala2019:E009:o6_normalization"
_HH_BOUND_VALUE_ID = "PDG2026-HaischHala2019:E009:c6_bound"


def _yaml_document():
    with open(_SIDECAR) as handle:
        return yaml.safe_load(handle)


def _yaml_pdg_block():
    return _yaml_document()["pdg_or_equivalent"]


def _value_entry(pdg, value_id: str):
    for item in pdg["values"]:
        if item["value_id"] == value_id:
            return item
    raise AssertionError(f"missing test value_id {value_id}")


def _primary_measurement_value_e_cm(block) -> float:
    return float(block["value"]) * _X1E_MINUS_26_E_CM


def _primary_measurement_total_uncertainty_e_cm(block) -> float:
    stat = float(block["statistical_uncertainty"])
    syst = float(block["systematic_uncertainty"])
    return math.sqrt(stat * stat + syst * syst) * _X1E_MINUS_26_E_CM


def _unrounded_bound_gev_minus2(document, response_value_id: str) -> float:
    pdg = document["pdg_or_equivalent"]
    conversion = document["auxiliary_theory_inputs"]["current_neutron_limit_conversion"]
    response = _value_entry(pdg, response_value_id)
    coefficient_gev = float(response["coefficient_value"]) * _MEV_TO_GEV
    return float(conversion["neutron_limit_gev_inverse"]) / coefficient_gev


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.INFO
    assert constraint.family == "edm_neutrino"
    assert constraint.observable == "Weinberg three-gluon operator reference bounds"


def test_anchor_matches_yaml_and_loud_fail_probe():
    constraint = fcc.get(_PID)
    pdg = _yaml_pdg_block()
    measured = pdg["measured_experimental_anchor"]
    neutron = measured["neutron_edm"]
    primary = measured["primary_experiment"]
    pr_response = _value_entry(pdg, _PR_RESPONSE_VALUE_ID)
    pr_bound = _value_entry(pdg, _PR_BOUND_VALUE_ID)
    hh_response = _value_entry(pdg, _HH_RESPONSE_VALUE_ID)
    hh_bound = _value_entry(pdg, _HH_BOUND_VALUE_ID)

    assert constraint.anchor.neutron_edm_limit.value == pytest.approx(neutron["value"])
    assert constraint.anchor.neutron_edm_limit.source_url == neutron["source_url"]
    assert constraint.anchor.neutron_edm_limit.limit_operator == neutron["limit_operator"]
    assert constraint.anchor.neutron_edm_limit.confidence_level == neutron["confidence_level"]
    assert constraint.anchor.neutron_edm_limit.table_value == pytest.approx(
        neutron["table_value"]
    )
    assert constraint.anchor.primary_experiment.central_value_e_cm == pytest.approx(
        _primary_measurement_value_e_cm(primary)
    )
    assert constraint.anchor.primary_experiment.total_uncertainty_e_cm == pytest.approx(
        _primary_measurement_total_uncertainty_e_cm(primary)
    )
    assert constraint.anchor.pospelov_ritz_response.coefficient_gev == pytest.approx(
        pr_response["coefficient_value"] * _MEV_TO_GEV
    )
    assert constraint.anchor.pospelov_ritz_w_bound.value == pytest.approx(
        pr_bound["derived_bound_value"]
    )
    assert constraint.anchor.haisch_hala_response.coefficient_gev == pytest.approx(
        hh_response["coefficient_value"] * _MEV_TO_GEV
    )
    assert constraint.anchor.haisch_hala_c6_bound.value == pytest.approx(
        hh_bound["derived_bound_value"]
    )
    assert constraint.anchor.haisch_hala_response.coefficient_relative_uncertainty == (
        pytest.approx(hh_response["coefficient_relative_uncertainty"])
    )
    assert constraint.anchor.reference_bound.value_id == _HH_BOUND_VALUE_ID
    assert constraint.anchor.budget == pytest.approx(hh_bound["derived_bound_value"])

    with pytest.raises(anchors.AnchorError, match="missing E009 value_id"):
        e009_module._load_value_anchor(
            _PID,
            value_id="no_such_e009_value_id",
            value_key="derived_bound_value",
        )
    with pytest.raises(anchors.AnchorError, match="none of the expected anchor keys"):
        e009_module._load_nested_anchor(
            _PID,
            parent_key="measured_experimental_anchor",
            nested_key="no_such_e009_nested_block",
        )


def test_e009_anchors_route_through_scaffold_load_anchor(monkeypatch):
    calls: list[tuple[str, ...]] = []
    original_load_anchor = e009_module.load_anchor

    def spy_load_anchor(*args, **kwargs):
        calls.append(tuple(kwargs["candidates"]))
        return original_load_anchor(*args, **kwargs)

    monkeypatch.setattr(e009_module, "load_anchor", spy_load_anchor)
    anchor = e009_module._load_e009_anchor(_PID)

    assert ("measured_experimental_anchor.neutron_edm",) in calls
    assert ("measured_experimental_anchor.primary_experiment",) in calls
    assert (f"values.{_PR_BOUND_VALUE_ID}",) in calls
    assert (f"values.{_HH_BOUND_VALUE_ID}",) in calls
    assert anchor.value == pytest.approx(fcc.get(_PID).anchor.value)


def test_e009_anchor_rejects_mismatched_load_anchor_block_key(monkeypatch):
    original_load_anchor = e009_module.load_anchor

    def mismatched_load_anchor(*args, **kwargs):
        anchor = original_load_anchor(*args, **kwargs)
        if kwargs["candidates"] == (f"values.{_HH_BOUND_VALUE_ID}",):
            return replace(anchor, block_key="wrong_block")
        return anchor

    monkeypatch.setattr(e009_module, "load_anchor", mismatched_load_anchor)
    with pytest.raises(anchors.AnchorError, match="load_anchor selected 'wrong_block'"):
        e009_module._load_e009_anchor(_PID)


def test_stub_numeric_validation_matches_independent_yaml_recompute():
    constraint = fcc.get(_PID)
    document = _yaml_document()
    pdg = document["pdg_or_equivalent"]
    pr_bound = _value_entry(pdg, _PR_BOUND_VALUE_ID)
    hh_bound = _value_entry(pdg, _HH_BOUND_VALUE_ID)
    pr_unrounded = _unrounded_bound_gev_minus2(document, _PR_RESPONSE_VALUE_ID)
    hh_unrounded = _unrounded_bound_gev_minus2(document, _HH_RESPONSE_VALUE_ID)

    result = constraint.evaluate(point_builder.empty_point())

    assert constraint.anchor.pospelov_ritz_w_bound.value == pytest.approx(
        pr_unrounded,
        rel=0.02,
    )
    assert constraint.anchor.haisch_hala_c6_bound.value == pytest.approx(
        hh_unrounded,
        rel=0.04,
    )
    assert result.predicted is None
    assert result.sm_prediction is None
    assert result.ratio is None
    assert result.experimental == pytest.approx(hh_bound["derived_bound_value"])
    assert result.budget == pytest.approx(hh_bound["derived_bound_value"])
    assert result.diagnostics["pospelov_ritz_w_bound_gev_minus2"] == pytest.approx(
        pr_bound["derived_bound_value"]
    )
    assert result.diagnostics["haisch_hala_c6_bound_gev_minus2"] == pytest.approx(
        hh_bound["derived_bound_value"]
    )
    assert result.diagnostics["pospelov_ritz_response_coefficient_gev"] == pytest.approx(
        _value_entry(pdg, _PR_RESPONSE_VALUE_ID)["coefficient_value"] * _MEV_TO_GEV
    )
    assert result.diagnostics["haisch_hala_response_coefficient_gev"] == pytest.approx(
        _value_entry(pdg, _HH_RESPONSE_VALUE_ID)["coefficient_value"] * _MEV_TO_GEV
    )
    assert result.diagnostics["no_hadronic_calculation"] is True
    assert result.diagnostics["no_rs_cp_odd_gluonic_matching"] is True


def test_adapter_reference_comparison_has_pass_fail_behavior_without_matching():
    bound = fcc.get(_PID).anchor.budget

    safe = compare_weinberg_coefficient_to_reference_bound(
        coefficient_gev_minus2=0.5 * bound,
        reference_bound_gev_minus2=bound,
    )
    excluded = compare_weinberg_coefficient_to_reference_bound(
        coefficient_gev_minus2=2.0 * bound,
        reference_bound_gev_minus2=bound,
    )

    assert safe.passes is True
    assert safe.ratio_to_bound == pytest.approx(0.5)
    assert excluded.passes is False
    assert excluded.ratio_to_bound == pytest.approx(2.0)


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
    assert empty_result.ratio is None
    assert empty_result.experimental == pytest.approx(constraint.anchor.value)
    assert empty_result.budget == pytest.approx(constraint.anchor.budget)
    assert empty_result.diagnostics["non_vetoing"] is True
    assert empty_result.diagnostics["parameter_point_inputs_used"] == ()
    assert empty_result.diagnostics["hadronic_matrix_elements_available"] is False
    assert empty_result.diagnostics["rs_cp_odd_gluonic_matching_available"] is False
    assert "NEEDS-HUMAN-PHYSICS" in (
        empty_result.diagnostics["needs_human_physics_hadronic_matrix_elements"]
    )
    assert "NEEDS-HUMAN-PHYSICS" in (
        empty_result.diagnostics["needs_human_physics_rs_cp_odd_gluonic_matching"]
    )


def test_evaluate_has_real_finite_numeric_fields():
    result = fcc.get(_PID).evaluate(point_builder.empty_point())

    assert result.predicted is None
    assert result.sm_prediction is None
    assert result.ratio is None
    for value in (result.experimental, result.budget):
        assert isinstance(value, float)
        assert math.isfinite(value)
    for key in (
        "reference_bound_gev_minus2",
        "pospelov_ritz_w_bound_gev_minus2",
        "pospelov_ritz_response_coefficient_gev",
        "haisch_hala_c6_bound_gev_minus2",
        "haisch_hala_response_coefficient_gev",
        "haisch_hala_response_relative_uncertainty",
        "neutron_edm_limit_e_cm",
        "neutron_edm_limit_gev_inverse",
        "primary_experiment_central_e_cm",
        "primary_experiment_total_uncertainty_e_cm",
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
