"""Production tests for E006 (Mercury-199 EDM stub)."""

from __future__ import annotations

import math
from pathlib import Path

import pytest
import yaml

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import anchors, point_builder
from flavor_catalog_constraints.base import ConstraintProtocol, Severity
from flavor_catalog_constraints.physics_adapters.mercury_edm import (
    compare_mercury_edm_measurement_to_limit,
)
from flavor_catalog_constraints.primary.edm_neutrino import E006 as e006_module

_PID = "E006"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "edm_neutrino" / "E006.yaml"
_X1E_MINUS_30_E_CM = 1.0e-30


def _yaml_pdg_block():
    with open(_SIDECAR) as handle:
        return yaml.safe_load(handle)["pdg_or_equivalent"]


def _measurement_value_e_cm(block) -> float:
    return float(block["central_value_from_arxiv_v4"]) * _X1E_MINUS_30_E_CM


def _measurement_stat_e_cm(block) -> float:
    return float(block["statistical_uncertainty"]) * _X1E_MINUS_30_E_CM


def _measurement_syst_e_cm(block) -> float:
    return float(block["systematic_uncertainty"]) * _X1E_MINUS_30_E_CM


def _measurement_total_uncertainty_e_cm(block) -> float:
    stat = float(block["statistical_uncertainty"])
    syst = float(block["systematic_uncertainty"])
    return math.sqrt(stat * stat + syst * syst) * _X1E_MINUS_30_E_CM


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.INFO
    assert constraint.family == "edm_neutrino"
    assert constraint.observable == "|d_Hg|"


def test_anchor_matches_yaml_and_loud_fail_probe():
    constraint = fcc.get(_PID)
    pdg = _yaml_pdg_block()
    limit = pdg["canonical_direct_limit"]
    improvement = pdg["post_2008_improvement_factor"]
    cross_ref = pdg["pdg_cross_reference_neutron"]
    theory = pdg["theory_translation_context"]

    assert constraint.anchor.limit.experimental.value == pytest.approx(limit["value"])
    assert constraint.anchor.limit.experimental.block_key == "canonical_direct_limit"
    assert constraint.anchor.limit.experimental.source_url == limit["source_url"]
    assert constraint.anchor.limit.experimental.units == "e cm"
    assert constraint.anchor.budget == pytest.approx(7.4e-30)
    assert constraint.anchor.limit.limit_operator == limit["limit_operator"]
    assert constraint.anchor.limit.confidence_level == limit["confidence_level"]
    assert constraint.anchor.limit.central_value_e_cm == pytest.approx(
        _measurement_value_e_cm(limit)
    )
    assert constraint.anchor.limit.statistical_uncertainty_e_cm == pytest.approx(
        _measurement_stat_e_cm(limit)
    )
    assert constraint.anchor.limit.systematic_uncertainty_e_cm == pytest.approx(
        _measurement_syst_e_cm(limit)
    )
    assert constraint.anchor.limit.total_uncertainty_e_cm == pytest.approx(
        _measurement_total_uncertainty_e_cm(limit)
    )
    assert constraint.anchor.improvement_factor.value == pytest.approx(
        improvement["value"]
    )
    assert constraint.anchor.pdg_cross_reference.graner_row_value == pytest.approx(
        cross_ref["graner_row_value"]
    )
    assert constraint.anchor.pdg_cross_reference.sahoo_row_value == pytest.approx(
        cross_ref["sahoo_row_value"]
    )
    assert constraint.anchor.theory_translation_context.theta_bar_limit == pytest.approx(
        theory["theta_bar_limit"]
    )

    with pytest.raises(anchors.AnchorError, match="none of the expected anchor keys"):
        anchors.load_anchor(
            _PID,
            family="edm_neutrino",
            candidates=("no_such_e006_block",),
        )
    with pytest.raises(anchors.AnchorError, match="has no 'missing_hg_limit' field"):
        anchors.load_anchor(
            _PID,
            family="edm_neutrino",
            candidates=("canonical_direct_limit",),
            value_key="missing_hg_limit",
        )


def test_anchor_rejects_improvement_factor_as_limit(monkeypatch):
    improvement = _yaml_pdg_block()["post_2008_improvement_factor"]

    def improvement_load_anchor(*args, **kwargs):
        return anchors.Anchor(
            process_id=_PID,
            block_key="post_2008_improvement_factor",
            value=improvement["value"],
            uncertainty=None,
            observable=improvement["observable"],
            units=improvement["units"],
            source=improvement["source"],
            source_url=improvement["source_url"],
            year=improvement["year"],
            snapshot_path=improvement["snapshot_path"],
        )

    monkeypatch.setattr(e006_module, "load_anchor", improvement_load_anchor)
    with pytest.raises(anchors.AnchorError, match="canonical_direct_limit"):
        e006_module._load_e006_anchor(_PID)


def test_stub_numeric_validation_matches_independent_yaml_recompute():
    constraint = fcc.get(_PID)
    limit = _yaml_pdg_block()["canonical_direct_limit"]
    expected_limit = float(limit["value"])
    measurement = _measurement_value_e_cm(limit)
    total_uncertainty = _measurement_total_uncertainty_e_cm(limit)
    expected_ratio = abs(measurement) / expected_limit

    result = constraint.evaluate(point_builder.empty_point())

    assert result.predicted is None
    assert result.sm_prediction is None
    assert result.experimental == pytest.approx(expected_limit)
    assert result.budget == pytest.approx(expected_limit)
    assert result.ratio == pytest.approx(expected_ratio)
    assert result.ratio == pytest.approx(2.20 / 7.4)
    assert result.diagnostics["measurement_central_e_cm"] == pytest.approx(measurement)
    assert result.diagnostics["measurement_abs_e_cm"] == pytest.approx(abs(measurement))
    assert result.diagnostics["total_uncertainty_e_cm"] == pytest.approx(
        total_uncertainty
    )
    assert result.diagnostics["no_schiff_moment_calculation"] is True
    assert result.diagnostics["no_atomic_structure_calculation"] is True
    assert result.diagnostics["no_rs_cp_odd_matching"] is True


def test_adapter_limit_comparison_has_pass_fail_behavior_without_hg_matching():
    limit = 7.4e-30

    safe = compare_mercury_edm_measurement_to_limit(
        measured_mercury_edm_e_cm=0.5 * limit,
        experimental_limit_e_cm=limit,
    )
    over_limit = compare_mercury_edm_measurement_to_limit(
        measured_mercury_edm_e_cm=2.0 * limit,
        experimental_limit_e_cm=limit,
    )

    assert safe.passes is True
    assert safe.ratio_to_limit == pytest.approx(0.5)
    assert over_limit.passes is False
    assert over_limit.ratio_to_limit == pytest.approx(2.0)


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
    assert empty_result.ratio == pytest.approx(2.20 / 7.4)
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
        "experimental_limit_e_cm",
        "measurement_central_e_cm",
        "measurement_abs_e_cm",
        "measurement_to_limit_ratio",
        "statistical_uncertainty_e_cm",
        "systematic_uncertainty_e_cm",
        "total_uncertainty_e_cm",
        "improvement_factor",
        "pdg_cross_reference_graner_row_value",
        "pdg_cross_reference_sahoo_row_value",
        "theory_translation_neutron_edm_limit",
        "theory_translation_proton_edm_limit",
        "theory_translation_theta_bar_limit",
        "theory_translation_up_down_cedm_difference_limit",
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
