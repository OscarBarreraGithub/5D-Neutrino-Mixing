"""Production tests for E008 (quark cEDM stub)."""

from __future__ import annotations

from dataclasses import replace
import math
from pathlib import Path

import pytest
import yaml

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import anchors, point_builder
from flavor_catalog_constraints.base import ConstraintProtocol, Severity
from flavor_catalog_constraints.physics_adapters.quark_cedm import (
    compare_quark_cedm_to_reference_bound,
)
from flavor_catalog_constraints.primary.edm_neutrino import E008 as e008_module

_PID = "E008"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "edm_neutrino" / "E008.yaml"
_X1E_MINUS_30_E_CM = 1.0e-30


def _yaml_pdg_block():
    with open(_SIDECAR) as handle:
        return yaml.safe_load(handle)["pdg_or_equivalent"]


def _neutron_unrounded_bound_cm(pdg) -> float:
    neutron_limit = float(pdg["experimental_anchors"]["neutron_edm"]["value"])
    coefficient = float(
        pdg["qcedm_translations"]["neutron_combination"][
            "coefficient_for_qcedm_combination"
        ]
    )
    return neutron_limit / coefficient


def _mercury_unrounded_bound_cm(pdg) -> float:
    mercury_limit = (
        float(pdg["experimental_anchors"]["mercury_199_edm"]["limit_value"])
        * _X1E_MINUS_30_E_CM
    )
    coefficient = float(
        pdg["qcedm_translations"]["mercury_isovector_combination"][
            "coefficient_for_isovector_qcedm"
        ]
    )
    return mercury_limit / coefficient


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.INFO
    assert constraint.family == "edm_neutrino"
    assert constraint.observable == "quark chromo-EDM reference bounds"


def test_anchor_matches_yaml_and_loud_fail_probe():
    constraint = fcc.get(_PID)
    pdg = _yaml_pdg_block()
    neutron_exp = pdg["experimental_anchors"]["neutron_edm"]
    mercury_exp = pdg["experimental_anchors"]["mercury_199_edm"]
    neutron = pdg["qcedm_translations"]["neutron_combination"]
    mercury = pdg["qcedm_translations"]["mercury_isovector_combination"]

    assert constraint.anchor.neutron_edm_limit.limit_value_e_cm == pytest.approx(
        neutron_exp["value"]
    )
    assert constraint.anchor.neutron_edm_limit.source_url == neutron_exp["source_url"]
    assert constraint.anchor.mercury_edm_limit.limit_value_e_cm == pytest.approx(
        mercury_exp["limit_value"] * _X1E_MINUS_30_E_CM
    )
    assert constraint.anchor.mercury_edm_limit.source_url == mercury_exp["source_url"]
    assert constraint.anchor.neutron_combination.value == pytest.approx(
        neutron["derived_bound_value"]
    )
    assert constraint.anchor.neutron_combination.coefficient_value == pytest.approx(
        neutron["coefficient_for_qcedm_combination"]
    )
    assert constraint.anchor.neutron_combination.derived_bound_observable == (
        neutron["derived_bound_observable"]
    )
    assert constraint.anchor.mercury_isovector_combination.value == pytest.approx(
        mercury["derived_bound_value"]
    )
    assert (
        constraint.anchor.mercury_isovector_combination.coefficient_value
        == pytest.approx(mercury["coefficient_for_isovector_qcedm"])
    )
    assert constraint.anchor.reference_bound.block_key == (
        "qcedm_translations.mercury_isovector_combination"
    )
    assert constraint.anchor.budget == pytest.approx(mercury["derived_bound_value"])

    with pytest.raises(anchors.AnchorError, match="none of the expected anchor keys"):
        e008_module._load_qcedm_translation_bound(
            _PID,
            nested_key="no_such_e008_translation",
            coefficient_key="coefficient_for_qcedm_combination",
        )


def test_qcedm_translation_anchors_route_through_scaffold_load_anchor(monkeypatch):
    calls: list[tuple[str, ...]] = []
    original_load_anchor = e008_module.load_anchor

    def spy_load_anchor(*args, **kwargs):
        calls.append(tuple(kwargs["candidates"]))
        return original_load_anchor(*args, **kwargs)

    monkeypatch.setattr(e008_module, "load_anchor", spy_load_anchor)
    anchor = e008_module._load_e008_anchor(_PID)

    assert (
        "qcedm_translations.neutron_combination",
    ) in calls
    assert (
        "qcedm_translations.mercury_isovector_combination",
    ) in calls
    assert anchor.value == pytest.approx(fcc.get(_PID).anchor.value)


def test_qcedm_translation_anchor_rejects_mismatched_load_anchor_block_key(monkeypatch):
    original_load_anchor = e008_module.load_anchor

    def mismatched_load_anchor(*args, **kwargs):
        anchor = original_load_anchor(*args, **kwargs)
        if kwargs["candidates"] == (
            "qcedm_translations.mercury_isovector_combination",
        ):
            return replace(anchor, block_key="wrong_block")
        return anchor

    monkeypatch.setattr(e008_module, "load_anchor", mismatched_load_anchor)
    with pytest.raises(anchors.AnchorError, match="load_anchor selected 'wrong_block'"):
        e008_module._load_e008_anchor(_PID)


def test_stub_numeric_validation_matches_independent_yaml_recompute():
    constraint = fcc.get(_PID)
    pdg = _yaml_pdg_block()
    neutron = pdg["qcedm_translations"]["neutron_combination"]
    mercury = pdg["qcedm_translations"]["mercury_isovector_combination"]
    neutron_unrounded = _neutron_unrounded_bound_cm(pdg)
    mercury_unrounded = _mercury_unrounded_bound_cm(pdg)

    result = constraint.evaluate(point_builder.empty_point())

    assert constraint.anchor.neutron_unrounded_bound_cm == pytest.approx(
        neutron_unrounded
    )
    assert constraint.anchor.mercury_unrounded_bound_cm == pytest.approx(
        mercury_unrounded
    )
    assert constraint.anchor.neutron_combination.value == pytest.approx(
        neutron_unrounded,
        rel=0.03,
    )
    assert constraint.anchor.mercury_isovector_combination.value == pytest.approx(
        mercury_unrounded,
        rel=0.05,
    )
    assert result.predicted is None
    assert result.sm_prediction is None
    assert result.ratio is None
    assert result.experimental == pytest.approx(mercury["derived_bound_value"])
    assert result.budget == pytest.approx(mercury["derived_bound_value"])
    assert result.diagnostics["neutron_qcedm_bound_cm"] == pytest.approx(
        neutron["derived_bound_value"]
    )
    assert result.diagnostics["mercury_qcedm_bound_cm"] == pytest.approx(
        mercury["derived_bound_value"]
    )
    assert result.diagnostics["neutron_qcedm_unrounded_bound_cm"] == pytest.approx(
        neutron_unrounded
    )
    assert result.diagnostics["mercury_qcedm_unrounded_bound_cm"] == pytest.approx(
        mercury_unrounded
    )
    assert result.diagnostics["no_hadronic_calculation"] is True
    assert result.diagnostics["no_rs_cp_odd_cedm_matching"] is True


def test_adapter_reference_comparison_has_pass_fail_behavior_without_matching():
    bound = fcc.get(_PID).anchor.budget

    safe = compare_quark_cedm_to_reference_bound(
        qcedm_combination_cm=0.5 * bound,
        reference_bound_cm=bound,
    )
    excluded = compare_quark_cedm_to_reference_bound(
        qcedm_combination_cm=2.0 * bound,
        reference_bound_cm=bound,
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
    assert (
        empty_result.diagnostics["rs_cp_odd_quark_cedm_matching_available"]
        is False
    )
    assert "NEEDS-HUMAN-PHYSICS" in (
        empty_result.diagnostics["needs_human_physics_hadronic_matrix_elements"]
    )
    assert "NEEDS-HUMAN-PHYSICS" in (
        empty_result.diagnostics["needs_human_physics_rs_cp_odd_matching"]
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
        "reference_bound_cm",
        "neutron_qcedm_bound_cm",
        "neutron_qcedm_unrounded_bound_cm",
        "neutron_qcedm_coefficient",
        "neutron_edm_limit_e_cm",
        "mercury_qcedm_bound_cm",
        "mercury_qcedm_unrounded_bound_cm",
        "mercury_qcedm_coefficient",
        "mercury_edm_limit_e_cm",
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
