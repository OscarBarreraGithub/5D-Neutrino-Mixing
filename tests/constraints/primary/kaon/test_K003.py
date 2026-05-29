"""Production tests for K003 (Re(epsilon'/epsilon) stub)."""

from __future__ import annotations

import math
from pathlib import Path
import re

import pytest
import yaml

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import point_builder
from flavor_catalog_constraints.anchors import Anchor, AnchorError
from flavor_catalog_constraints.base import ConstraintProtocol, Severity
from flavor_catalog_constraints.physics_adapters.kaon_direct_cp import (
    compare_epsilon_prime_np_room_to_measurement,
)
from flavor_catalog_constraints.primary.kaon import K003 as k003_module

_PID = "K003"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "kaon" / "K003.yaml"
_NUMBER_RE = re.compile(r"[0-9]+(?:\.[0-9]*)?(?:[eE][+-]?[0-9]+)?")


def _yaml():
    with open(_SIDECAR) as handle:
        return yaml.safe_load(handle)


def _yaml_pdg_block():
    return _yaml()["pdg_or_equivalent"]


def _yaml_supporting():
    return {entry["key"]: entry for entry in _yaml()["supporting_values"]}


def _scale(units: str) -> float:
    if units == "dimensionless":
        return 1.0
    if units == "10^-3":
        return 1.0e-3
    if units == "10^-4":
        return 1.0e-4
    raise AssertionError(f"unexpected units {units!r}")


def _scaled_value(block) -> float:
    return float(block["value"]) * _scale(block["units"])


def _scaled_uncertainty(block) -> float:
    raw = block["uncertainty"]
    scale = _scale(block["units"])
    if isinstance(raw, str):
        components = [float(match.group(0)) * scale for match in _NUMBER_RE.finditer(raw)]
        return math.sqrt(sum(component * component for component in components))
    return float(raw) * scale


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.INFO
    assert constraint.family == "kaon"
    assert constraint.observable == "Re(epsilon'/epsilon)"


def test_anchor_matches_yaml_and_budget():
    constraint = fcc.get(_PID)
    pdg = _yaml_pdg_block()
    supporting = _yaml_supporting()

    assert constraint.anchor.block_key == "pdg_or_equivalent"
    assert constraint.anchor.value == pytest.approx(pdg["value"])
    assert constraint.anchor.uncertainty == pytest.approx(pdg["uncertainty"])
    assert constraint.anchor.units == pdg["units"]
    assert constraint.anchor.source == pdg["source"]
    assert constraint.anchor.source_url == pdg["source_url"]
    assert constraint.anchor.snapshot_path == pdg["snapshot_path"]
    assert constraint.anchor.display_value == pdg["display_value"]
    assert constraint.anchor.sha256 == pdg["sha256_of_text_snapshot"]
    assert constraint.anchor.budget == pytest.approx(abs(float(pdg["value"])))

    for key in ("ktev_2011", "na48_2002", "rbc_ukqcd_2020"):
        parsed = constraint.anchor.supporting(key)
        raw = supporting[key]
        assert parsed.value == pytest.approx(_scaled_value(raw))
        assert parsed.uncertainty == pytest.approx(_scaled_uncertainty(raw))
        assert parsed.source_url == raw["source_url"]
        assert parsed.snapshot_path == raw["snapshot_path"]

    assert constraint.anchor.sm_context.key == "rbc_ukqcd_2020"
    assert constraint.anchor.sm_context.value == pytest.approx(21.7e-4)


def test_k003_flat_anchor_routes_through_scaffold_load_anchor(monkeypatch):
    calls: list[tuple[str, ...]] = []
    original_load_anchor = k003_module.load_anchor

    def spy_load_anchor(*args, **kwargs):
        calls.append(tuple(kwargs["candidates"]))
        return original_load_anchor(*args, **kwargs)

    monkeypatch.setattr(k003_module, "load_anchor", spy_load_anchor)
    anchor = k003_module._load_k003_anchor(_PID)

    assert calls == [("pdg_or_equivalent",)]
    assert anchor.value == pytest.approx(fcc.get(_PID).anchor.value)

    def mismatched_load_anchor(*args, **kwargs):
        return Anchor(
            process_id=_PID,
            block_key="not_the_flat_pdg_block",
            value=1.0,
            uncertainty=0.1,
        )

    monkeypatch.setattr(k003_module, "load_anchor", mismatched_load_anchor)
    with pytest.raises(AnchorError, match="load_anchor selected"):
        k003_module._load_flat_pdg_anchor(_PID)


def test_anchor_loading_fails_loudly_for_missing_candidate_or_value():
    with pytest.raises(fcc.AnchorError, match="none of the expected anchor keys"):
        fcc.load_anchor(_PID, family="kaon", candidates=("missing_k003_anchor",))
    with pytest.raises(k003_module.AnchorError, match="has no 'missing_value' field"):
        k003_module._load_flat_pdg_anchor(_PID, value_key="missing_value")


def test_stub_numeric_validation_matches_independent_yaml_recompute():
    constraint = fcc.get(_PID)
    pdg = _yaml_pdg_block()
    supporting = _yaml_supporting()
    rbc = supporting["rbc_ukqcd_2020"]
    expected_value = float(pdg["value"])
    expected_uncertainty = float(pdg["uncertainty"])
    expected_budget = abs(expected_value)
    expected_ratio = abs(expected_value) / expected_budget
    expected_sm_context = _scaled_value(rbc)
    expected_sm_context_uncertainty = _scaled_uncertainty(rbc)

    result = constraint.evaluate(point_builder.empty_point())

    assert result.predicted is None
    assert result.sm_prediction == pytest.approx(expected_sm_context)
    assert result.experimental == pytest.approx(expected_value)
    assert result.budget == pytest.approx(expected_budget)
    assert result.ratio == pytest.approx(expected_ratio)
    assert result.ratio == pytest.approx(1.0)
    assert result.diagnostics["experimental_uncertainty"] == pytest.approx(
        expected_uncertainty
    )
    assert result.diagnostics["sm_context_value"] == pytest.approx(expected_sm_context)
    assert result.diagnostics["sm_context_uncertainty"] == pytest.approx(
        expected_sm_context_uncertainty
    )
    assert result.diagnostics["no_penguin_calculation"] is True
    assert result.diagnostics["sm_qcd_electroweak_penguin_cancellation"] is True
    assert result.diagnostics["rs_delta_s1_penguin_matching_available"] is False


def test_adapter_room_comparison_has_pass_fail_behavior_without_penguins():
    measurement = 1.66e-3
    uncertainty = 0.23e-3

    safe = compare_epsilon_prime_np_room_to_measurement(
        measured_re_epsilon_prime_over_epsilon=measurement,
        experimental_uncertainty=uncertainty,
        documented_np_room_abs=measurement,
    )
    overfilled = compare_epsilon_prime_np_room_to_measurement(
        measured_re_epsilon_prime_over_epsilon=measurement,
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
    assert empty_result.sm_prediction == pytest.approx(constraint.anchor.sm_context.value)
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
    for value in (result.sm_prediction, result.experimental, result.ratio, result.budget):
        assert isinstance(value, float)
        assert math.isfinite(value)
    for key in (
        "experimental_uncertainty",
        "measurement_abs",
        "documented_np_room_abs",
        "measurement_to_np_room_ratio",
        "sm_context_value",
        "sm_context_uncertainty",
        "ktev_2011_value",
        "na48_2002_value",
        "aebischer_buras_2020_octet_value",
        "aebischer_buras_2020_nonet_value",
    ):
        assert isinstance(result.diagnostics[key], float)
        assert math.isfinite(result.diagnostics[key])
    for entry in result.diagnostics["supporting_values"]:
        assert isinstance(entry["value"], float)
        assert math.isfinite(entry["value"])
        assert isinstance(entry["scale"], float)
        assert math.isfinite(entry["scale"])
        for uncertainty in entry["uncertainty_components"]:
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
