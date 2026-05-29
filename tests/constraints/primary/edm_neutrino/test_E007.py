"""Production tests for E007 (Ra-225/Xe-129 atomic EDM stub)."""

from __future__ import annotations

import math
from pathlib import Path

import pytest
import yaml

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import anchors, point_builder
from flavor_catalog_constraints.base import ConstraintProtocol, Severity
from flavor_catalog_constraints.physics_adapters.atomic_edm import (
    compare_atomic_edm_to_limit,
)

_PID = "E007"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "edm_neutrino" / "E007.yaml"


def _yaml_pdg_block():
    with open(_SIDECAR) as handle:
        return yaml.safe_load(handle)["pdg_or_equivalent"]


def _yaml_post_context():
    with open(_SIDECAR) as handle:
        return yaml.safe_load(handle)["post_2008_context"]


def _xe_best_total_uncertainty(block) -> float:
    stat = float(block["statistical_uncertainty"])
    syst = float(block["systematic_uncertainty"])
    return math.sqrt(stat * stat + syst * syst)


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.INFO
    assert constraint.family == "edm_neutrino"
    assert constraint.observable == "|d_Ra|, |d_Xe|"


def test_anchor_matches_yaml_and_loud_fail_probe():
    constraint = fcc.get(_PID)
    pdg = _yaml_pdg_block()
    post = _yaml_post_context()
    ra = pdg["ra225_current_direct_limit"]
    ra_first = pdg["ra225_first_measurement"]
    xe = pdg["xe129_best_direct_limit"]
    xe_independent = pdg["xe129_independent_heidelberg_limit"]

    assert constraint.anchor.ra225_current_direct_limit.value == pytest.approx(
        ra["value"]
    )
    assert constraint.anchor.ra225_current_direct_limit.block_key == (
        "ra225_current_direct_limit"
    )
    assert constraint.anchor.ra225_current_direct_limit.source_url == ra["source_url"]
    assert constraint.anchor.ra225_current_direct_limit.limit_operator == "<"
    assert constraint.anchor.ra225_current_direct_limit.confidence_level == "95%"
    assert constraint.anchor.ra225_current_direct_limit.improvement_factor == pytest.approx(
        ra["improvement_factor"]
    )
    assert constraint.anchor.ra225_first_measurement.value == pytest.approx(
        ra_first["value"]
    )
    assert constraint.anchor.xe129_best_direct_limit.value == pytest.approx(
        xe["limit_value"]
    )
    assert constraint.anchor.xe129_best_direct_limit.central_value_e_cm == pytest.approx(
        xe["central_value"]
    )
    assert (
        constraint.anchor.xe129_best_direct_limit.statistical_uncertainty_e_cm
        == pytest.approx(xe["statistical_uncertainty"])
    )
    assert (
        constraint.anchor.xe129_best_direct_limit.systematic_uncertainty_e_cm
        == pytest.approx(xe["systematic_uncertainty"])
    )
    assert constraint.anchor.xe129_best_direct_limit.total_uncertainty_e_cm == pytest.approx(
        _xe_best_total_uncertainty(xe)
    )
    assert constraint.anchor.xe129_independent_heidelberg_limit.value == pytest.approx(
        xe_independent["limit_value"]
    )
    assert (
        constraint.anchor.xe129_independent_heidelberg_limit.central_value_e_cm
        == pytest.approx(xe_independent["central_value"])
    )
    assert constraint.anchor.budget == pytest.approx(xe["limit_value"])
    assert constraint.anchor.reference_limit.block_key == "xe129_best_direct_limit"
    assert constraint.anchor.argonne_ra225_program_context.numeric_values[
        "half_life_days"
    ] == pytest.approx(post["argonne_ra225_program_context"]["half_life_days"])

    with pytest.raises(anchors.AnchorError, match="none of the expected anchor keys"):
        anchors.load_anchor(
            _PID,
            family="edm_neutrino",
            candidates=("no_such_e007_block",),
        )
    with pytest.raises(anchors.AnchorError, match="has no 'value' field"):
        anchors.load_anchor(
            _PID,
            family="edm_neutrino",
            candidates=("xe129_best_direct_limit",),
        )


def test_stub_numeric_validation_matches_independent_yaml_recompute():
    constraint = fcc.get(_PID)
    pdg = _yaml_pdg_block()
    ra = pdg["ra225_current_direct_limit"]
    xe = pdg["xe129_best_direct_limit"]
    selected_limit = min(float(ra["value"]), float(xe["limit_value"]))
    expected_ratio = abs(float(xe["central_value"])) / float(xe["limit_value"])

    result = constraint.evaluate(point_builder.empty_point())

    assert result.predicted is None
    assert result.sm_prediction is None
    assert result.experimental == pytest.approx(selected_limit)
    assert result.budget == pytest.approx(selected_limit)
    assert result.ratio == pytest.approx(expected_ratio)
    assert result.ratio == pytest.approx(0.1)
    assert result.diagnostics["selected_limit_block"] == "xe129_best_direct_limit"
    assert result.diagnostics["ra225_current_limit_e_cm"] == pytest.approx(ra["value"])
    assert result.diagnostics["xe129_best_limit_e_cm"] == pytest.approx(
        xe["limit_value"]
    )
    assert result.diagnostics["xe129_best_central_value_e_cm"] == pytest.approx(
        xe["central_value"]
    )
    assert result.diagnostics["xe129_best_total_uncertainty_e_cm"] == pytest.approx(
        _xe_best_total_uncertainty(xe)
    )
    assert result.diagnostics["no_schiff_moment_calculation"] is True
    assert result.diagnostics["no_atomic_structure_calculation"] is True
    assert result.diagnostics["no_rs_cp_odd_matching"] is True


def test_adapter_limit_comparison_has_pass_fail_behavior_without_matching():
    limit = fcc.get(_PID).anchor.budget

    safe = compare_atomic_edm_to_limit(
        atomic_edm_e_cm=0.5 * limit,
        experimental_limit_e_cm=limit,
        isotope="129Xe",
    )
    excluded = compare_atomic_edm_to_limit(
        atomic_edm_e_cm=2.0 * limit,
        experimental_limit_e_cm=limit,
        isotope="129Xe",
    )

    assert safe.passes is True
    assert safe.ratio_to_limit == pytest.approx(0.5)
    assert excluded.passes is False
    assert excluded.ratio_to_limit == pytest.approx(2.0)


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
    assert empty_result.ratio == pytest.approx(0.1)
    assert empty_result.diagnostics["non_vetoing"] is True
    assert empty_result.diagnostics["parameter_point_inputs_used"] == ()
    assert (
        empty_result.diagnostics["schiff_moment_atomic_structure_available"]
        is False
    )
    assert empty_result.diagnostics["rs_cp_odd_matching_available"] is False
    assert "NEEDS-HUMAN-PHYSICS" in (
        empty_result.diagnostics["needs_human_physics_schiff_atomic_structure"]
    )
    assert "NEEDS-HUMAN-PHYSICS" in (
        empty_result.diagnostics["needs_human_physics_rs_cp_odd_matching"]
    )


def test_evaluate_has_real_finite_numeric_fields():
    result = fcc.get(_PID).evaluate(point_builder.empty_point())

    assert result.predicted is None
    assert result.sm_prediction is None
    for value in (result.experimental, result.ratio, result.budget):
        assert isinstance(value, float)
        assert math.isfinite(value)
    for key in (
        "selected_limit_e_cm",
        "selected_measurement_to_limit_ratio",
        "ra225_current_limit_e_cm",
        "ra225_current_improvement_factor",
        "ra225_first_measurement_limit_e_cm",
        "xe129_best_limit_e_cm",
        "xe129_best_central_value_e_cm",
        "xe129_best_statistical_uncertainty_e_cm",
        "xe129_best_systematic_uncertainty_e_cm",
        "xe129_best_total_uncertainty_e_cm",
        "xe129_best_measurement_to_limit_ratio",
        "xe129_best_improvement_factor",
        "xe129_independent_limit_e_cm",
        "xe129_independent_central_value_e_cm",
        "xe129_independent_total_uncertainty_e_cm",
        "argonne_ra225_half_life_days",
        "ra225_projected_statistical_sensitivity_e_cm",
        "ra225_projected_systematic_uncertainty_e_cm",
        "heidelberg_xe_upgrade_starting_limit_e_cm",
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
