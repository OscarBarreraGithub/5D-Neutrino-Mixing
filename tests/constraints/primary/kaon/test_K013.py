"""Production tests for K013 (K_L -> pi0 gamma gamma stub)."""

from __future__ import annotations

import math
from pathlib import Path

import pytest
import yaml

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import point_builder
from flavor_catalog_constraints.anchors import Anchor, AnchorError
from flavor_catalog_constraints.base import ConstraintProtocol, Severity
from flavor_catalog_constraints.physics_adapters.radiative_kaon import (
    compare_kl_pi0gammagamma_np_room_to_measurement,
)
from flavor_catalog_constraints.primary.kaon import K013 as k013_module

_PID = "K013"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "kaon" / "K013.yaml"


def _yaml_pdg_block():
    with open(_SIDECAR) as handle:
        return yaml.safe_load(handle)["pdg_or_equivalent"]


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.INFO
    assert constraint.family == "kaon"
    assert constraint.observable == "BR(K_L -> pi0 gamma gamma)"


def test_anchor_matches_yaml_and_budget():
    constraint = fcc.get(_PID)
    pdg = _yaml_pdg_block()
    inputs = pdg["average_inputs"]

    assert constraint.anchor.block_key == "pdg_or_equivalent"
    assert constraint.anchor.value == pytest.approx(pdg["value"])
    assert constraint.anchor.uncertainty == pytest.approx(pdg["uncertainty"])
    assert constraint.anchor.units == pdg["units"]
    assert constraint.anchor.source == pdg["source"]
    assert constraint.anchor.source_url == pdg["source_url"]
    assert constraint.anchor.snapshot_path == pdg["snapshot_path"]
    assert constraint.anchor.display == pdg["display"]
    assert constraint.anchor.sha256 == pdg["sha256"]
    assert constraint.anchor.budget == pytest.approx(abs(float(pdg["value"])))
    assert len(constraint.anchor.average_inputs) == len(inputs)

    for parsed, raw in zip(constraint.anchor.average_inputs, inputs):
        assert parsed.label == raw["label"]
        assert parsed.experiment == raw["experiment"]
        assert parsed.value == pytest.approx(raw["value"])
        assert parsed.uncertainties == pytest.approx(
            tuple(float(value) for value in raw["uncertainties"])
        )
        assert (
            parsed.original_reported_branching_fraction
            == raw["original_reported_branching_fraction"]
        )
        assert parsed.vector_exchange_parameter == raw["vector_exchange_parameter"]
        assert parsed.snapshot_path == raw["snapshot_path"]
        assert parsed.sha256 == raw["sha256"]


def test_k013_flat_anchor_routes_through_scaffold_load_anchor(monkeypatch):
    calls: list[tuple[str, ...]] = []
    original_load_anchor = k013_module.load_anchor

    def spy_load_anchor(*args, **kwargs):
        calls.append(tuple(kwargs["candidates"]))
        return original_load_anchor(*args, **kwargs)

    monkeypatch.setattr(k013_module, "load_anchor", spy_load_anchor)
    anchor = k013_module._load_k013_anchor(_PID)

    assert calls == [("pdg_or_equivalent",)]
    assert anchor.value == pytest.approx(fcc.get(_PID).anchor.value)

    def mismatched_load_anchor(*args, **kwargs):
        return Anchor(
            process_id=_PID,
            block_key="not_the_flat_pdg_block",
            value=1.0,
            uncertainty=0.1,
        )

    monkeypatch.setattr(k013_module, "load_anchor", mismatched_load_anchor)
    with pytest.raises(AnchorError, match="load_anchor selected"):
        k013_module._load_flat_pdg_anchor(_PID)


def test_anchor_loading_fails_loudly_for_missing_candidate_or_value():
    with pytest.raises(fcc.AnchorError, match="none of the expected anchor keys"):
        fcc.load_anchor(_PID, family="kaon", candidates=("missing_k013_anchor",))
    with pytest.raises(k013_module.AnchorError, match="has no 'missing_br' field"):
        k013_module._load_flat_pdg_anchor(_PID, value_key="missing_br")


def test_stub_numeric_validation_matches_independent_yaml_recompute():
    constraint = fcc.get(_PID)
    pdg = _yaml_pdg_block()
    expected_value = float(pdg["value"])
    expected_uncertainty = float(pdg["uncertainty"])
    expected_budget = abs(expected_value)
    expected_ratio = abs(expected_value) / expected_budget
    expected_naive_significance = expected_value / expected_uncertainty

    result = constraint.evaluate(point_builder.empty_point())

    assert result.predicted is None
    assert result.sm_prediction is None
    assert result.experimental == pytest.approx(expected_value)
    assert result.budget == pytest.approx(expected_budget)
    assert result.ratio == pytest.approx(expected_ratio)
    assert result.ratio == pytest.approx(1.0)
    assert result.diagnostics["experimental_uncertainty"] == pytest.approx(
        expected_uncertainty
    )
    assert expected_naive_significance == pytest.approx(
        abs(result.experimental) / result.diagnostics["experimental_uncertainty"]
    )
    assert result.diagnostics["no_chpt_calculation"] is True
    assert result.diagnostics["no_np_matching"] is True


def test_adapter_room_comparison_has_pass_fail_behavior_without_chpt():
    measurement = 1.273e-6
    uncertainty = 0.033e-6

    safe = compare_kl_pi0gammagamma_np_room_to_measurement(
        measured_branching_fraction=measurement,
        experimental_uncertainty=uncertainty,
        documented_np_room_abs=measurement,
    )
    overfilled = compare_kl_pi0gammagamma_np_room_to_measurement(
        measured_branching_fraction=measurement,
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
    assert "non-vetoing only" in empty_result.diagnostics["passes_semantics"]
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
        "experimental_uncertainty",
        "measurement_abs",
        "documented_np_room_abs",
        "measurement_to_np_room_ratio",
        "pdg_average_input_count",
    ):
        assert isinstance(result.diagnostics[key], float)
        assert math.isfinite(result.diagnostics[key])
    for entry in result.diagnostics["pdg_average_inputs"]:
        assert isinstance(entry["value"], float)
        assert math.isfinite(entry["value"])
        for uncertainty in entry["uncertainties"]:
            assert isinstance(uncertainty, float)
            assert math.isfinite(uncertainty)


def test_evaluate_is_pure_and_deterministic():
    point = point_builder.make_point(kk_gluon_mass_gev=3000.0)
    before_extras = dict(point.extras)
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    assert dict(point.extras) == before_extras
