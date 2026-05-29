"""Production tests for C003 (direct charm Delta A_CP stub)."""

from __future__ import annotations

import math
from pathlib import Path

import pytest
import yaml

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import point_builder
from flavor_catalog_constraints.base import ConstraintProtocol, Severity
from flavor_catalog_constraints.physics_adapters.charm_direct_cp import (
    compare_delta_acp_np_room_to_measurement,
)

_PID = "C003"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "charm" / "C003.yaml"


def _yaml_pdg_block():
    with open(_SIDECAR) as handle:
        return yaml.safe_load(handle)["pdg_or_equivalent"]


def _scaled_lhcb_value(block) -> float:
    return float(block["value"]) * float(block["scale"])


def _scaled_lhcb_uncertainty(block) -> float:
    return float(block["uncertainty"]) * float(block["scale"])


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.INFO
    assert constraint.family == "charm"
    assert constraint.observable == "Delta A_CP(D0 -> K+K-, pi+pi-)"


def test_anchor_matches_yaml_and_budget():
    constraint = fcc.get(_PID)
    pdg = _yaml_pdg_block()
    lhcb = pdg["lhcb2019_discovery_average"]
    direct = pdg["hflav_direct_cpv_average"]
    indirect = pdg["hflav_indirect_cpv_companion"]
    no_cpv = pdg["hflav_no_cpv_test"]
    cfw = pdg["cfw2008_context"]
    expected_value = _scaled_lhcb_value(lhcb)
    expected_uncertainty = _scaled_lhcb_uncertainty(lhcb)

    assert constraint.anchor.lhcb_discovery.raw_value == pytest.approx(lhcb["value"])
    assert constraint.anchor.lhcb_discovery.raw_uncertainty == pytest.approx(
        lhcb["uncertainty"]
    )
    assert constraint.anchor.lhcb_discovery.scale == pytest.approx(float(lhcb["scale"]))
    assert constraint.anchor.value == pytest.approx(expected_value)
    assert constraint.anchor.uncertainty == pytest.approx(expected_uncertainty)
    assert constraint.anchor.lhcb_discovery.value_percent == pytest.approx(
        lhcb["value_percent"]
    )
    assert constraint.anchor.lhcb_discovery.uncertainty_percent == pytest.approx(
        lhcb["uncertainty_percent"]
    )
    assert constraint.anchor.source_url == lhcb["source_url"]
    assert constraint.anchor.hflav_direct_cpv.value_percent == pytest.approx(
        direct["value"]
    )
    assert constraint.anchor.hflav_direct_cpv.uncertainty_percent == pytest.approx(
        direct["uncertainty"]
    )
    assert constraint.anchor.hflav_indirect_cpv.value_percent == pytest.approx(
        indirect["value"]
    )
    assert constraint.anchor.hflav_no_cpv.delta_chi2 == pytest.approx(
        no_cpv["delta_chi2"]
    )
    assert constraint.anchor.hflav_no_cpv.significance_sigma == pytest.approx(
        no_cpv["significance_sigma"]
    )
    assert constraint.anchor.cfw_context.value_summary == cfw["value_summary"]
    assert constraint.anchor.budget == pytest.approx(abs(expected_value))
    assert constraint.anchor.budget == pytest.approx(1.54e-3)


def test_anchor_loading_fails_loudly_for_missing_candidate_or_value():
    with pytest.raises(fcc.AnchorError, match="none of the expected anchor keys"):
        fcc.load_anchor(_PID, family="charm", candidates=("missing_c003_anchor",))
    with pytest.raises(fcc.AnchorError, match="has no 'missing_delta_acp' field"):
        fcc.load_anchor(
            _PID,
            family="charm",
            candidates=("lhcb2019_discovery_average",),
            value_key="missing_delta_acp",
        )


def test_stub_numeric_validation_matches_independent_yaml_recompute():
    constraint = fcc.get(_PID)
    lhcb = _yaml_pdg_block()["lhcb2019_discovery_average"]
    expected_value = _scaled_lhcb_value(lhcb)
    expected_uncertainty = _scaled_lhcb_uncertainty(lhcb)
    expected_budget = abs(expected_value)
    expected_ratio = abs(expected_value) / expected_budget
    expected_naive_significance = abs(float(lhcb["value"])) / float(
        lhcb["uncertainty"]
    )

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
    assert result.diagnostics["lhcb_discovery_significance_sigma_naive"] == (
        pytest.approx(expected_naive_significance)
    )
    assert result.diagnostics["hflav_no_cpv_significance_sigma"] == pytest.approx(5.3)
    assert result.diagnostics["no_penguin_calculation"] is True


def test_adapter_room_comparison_has_pass_fail_behavior_without_penguins():
    measurement = 1.54e-3
    uncertainty = 0.29e-3

    safe = compare_delta_acp_np_room_to_measurement(
        measured_delta_acp=measurement,
        experimental_uncertainty=uncertainty,
        documented_np_room_abs=measurement,
    )
    overfilled = compare_delta_acp_np_room_to_measurement(
        measured_delta_acp=measurement,
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
        "experimental_value_raw_x1e_minus4",
        "experimental_uncertainty_raw_x1e_minus4",
        "experimental_scale",
        "experimental_value_percent",
        "experimental_uncertainty_percent",
        "experimental_uncertainty",
        "measurement_abs",
        "documented_np_room_abs",
        "measurement_to_np_room_ratio",
        "lhcb_discovery_significance_sigma_naive",
        "hflav_direct_cpv_value_percent",
        "hflav_direct_cpv_uncertainty_percent",
        "hflav_indirect_cpv_value_percent",
        "hflav_indirect_cpv_uncertainty_percent",
        "hflav_no_cpv_delta_chi2",
        "hflav_no_cpv_significance_sigma",
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
