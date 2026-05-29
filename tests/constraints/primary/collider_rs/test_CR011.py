"""Production tests for CR011 (longitudinal VBS)."""

from __future__ import annotations

import math
from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import point_builder
from flavor_catalog_constraints.anchors import AnchorError
from flavor_catalog_constraints.base import ConstraintProtocol, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters import vbs_longitudinal as adapter
from flavor_catalog_constraints.primary.collider_rs import CR011 as cr011_module

_PID = "CR011"
_ACTIVE_VALUE_ID = "ATLAS2025:CR011:fiducial_sigma_WLWL_same_sign_upper_limit"
_CMS_VALUE_ID = "CMS2020:CR011:fiducial_sigma_WLWL_same_sign_upper_limit"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "collider_rs" / "CR011.yaml"


def _yaml_pdg_block():
    with open(_SIDECAR) as handle:
        return yaml.safe_load(handle)["pdg_or_equivalent"]


def _entry(value_id: str):
    matches = [
        entry
        for entry in _yaml_pdg_block()["values"]
        if entry.get("value_id") == value_id
    ]
    assert len(matches) == 1
    return matches[0]


def _raw_prediction(sigma_fb: float, source: str = "unit-test human recast"):
    return {
        adapter.HUMAN_SUPPLIED_SIGMA_RAW_KEY: sigma_fb,
        adapter.HUMAN_SUPPLIED_SIGMA_SOURCE_RAW_KEY: source,
    }


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.INFO
    assert constraint.family == "collider_rs"
    assert constraint.observable == "sigma_fid(pp -> jj W_L^+/- W_L^+/-)"


def test_anchor_matches_yaml_and_loud_fail_probe():
    constraint = fcc.get(_PID)
    active = _entry(_ACTIVE_VALUE_ID)
    cms = _entry(_CMS_VALUE_ID)

    assert constraint.anchor.active_limit.value_fb == pytest.approx(active["value"])
    assert constraint.anchor.active_limit.units == active["units"]
    assert constraint.anchor.active_limit.source_url == active["source_url"]
    assert constraint.anchor.active_limit.display == active["display"]
    assert constraint.anchor.active_limit.experiment == active["experiment"]
    assert constraint.anchor.active_limit.cl == active["cl"]
    assert constraint.anchor.active_limit.expected_limit_value == pytest.approx(
        active["expected_limit"]["value"]
    )
    assert constraint.anchor.cms_2020.value_fb == pytest.approx(cms["value"])
    assert constraint.anchor.cms_2020.source_url == cms["source_url"]
    assert constraint.anchor.cms_2020.expected_limit_value == pytest.approx(
        cms["expected_limit"]["value"]
    )
    assert constraint.anchor.budget == pytest.approx(active["value"])

    with pytest.raises(AnchorError):
        cr011_module._load_value_anchor("not-a-real-value-id", process_id=_PID)


def test_evaluate_without_human_prediction_records_bound_non_vetoing():
    result = fcc.get(_PID).evaluate(point_builder.empty_point())

    assert result.process_id == _PID
    assert result.severity is Severity.INFO
    assert result.passes is True
    assert result.predicted is None
    assert result.ratio is None
    assert result.experimental == pytest.approx(fcc.get(_PID).anchor.value)
    assert result.budget == pytest.approx(fcc.get(_PID).anchor.budget)
    assert result.diagnostics["prediction_status"] == "missing-human-recast"
    assert result.diagnostics["human_prediction_raw_key"] == (
        adapter.HUMAN_SUPPLIED_SIGMA_RAW_KEY
    )
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics[
        "needs_human_physics_sm_eft"
    ]
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics[
        "needs_human_physics_rs_matching"
    ]


def test_human_prediction_numeric_comparison_is_independent_of_adapter():
    constraint = fcc.get(_PID)
    active = _entry(_ACTIVE_VALUE_ID)
    predicted = 0.225
    point = ParameterPoint(raw=_raw_prediction(predicted), extras={})
    result = constraint.evaluate(point)

    expected_ratio = predicted / float(active["value"])

    assert result.predicted == pytest.approx(predicted)
    assert result.experimental == pytest.approx(float(active["value"]))
    assert result.budget == pytest.approx(float(active["value"]))
    assert result.ratio == pytest.approx(expected_ratio)
    assert result.ratio == pytest.approx(0.5)
    assert result.passes is True
    assert result.diagnostics["prediction_status"] == "human-supplied"
    assert result.diagnostics["human_supplied_prediction_used"] is True


@pytest.mark.parametrize(
    ("predicted_sigma_fb", "expected_pass"),
    [
        (0.10, True),
        (0.90, False),
    ],
)
def test_safe_point_passes_and_excluded_point_fails_advisory_info(
    predicted_sigma_fb: float,
    expected_pass: bool,
):
    constraint = fcc.get(_PID)
    point = ParameterPoint(raw=_raw_prediction(predicted_sigma_fb), extras={})
    result = constraint.evaluate(point)

    assert result.passes is expected_pass
    assert result.severity is Severity.INFO
    assert result.predicted == pytest.approx(predicted_sigma_fb)
    assert result.ratio == pytest.approx(predicted_sigma_fb / constraint.anchor.budget)
    if expected_pass:
        assert result.ratio <= 1.0
    else:
        assert result.ratio > 1.0
        assert result.diagnostics["severity_policy"] == (
            "INFO/non-vetoing even when advisory passes=False"
        )


def test_object_raw_prediction_path_and_invalid_input_are_non_crashing():
    constraint = fcc.get(_PID)

    object_result = constraint.evaluate(
        ParameterPoint(
            raw=SimpleNamespace(
                cr011_human_fiducial_sigma_fb=0.20,
                cr011_human_fiducial_sigma_source="object-style human recast",
            ),
            extras={},
        )
    )
    invalid_result = constraint.evaluate(
        ParameterPoint(raw=_raw_prediction(-0.1), extras={})
    )

    assert object_result.predicted == pytest.approx(0.20)
    assert object_result.ratio == pytest.approx(0.20 / constraint.anchor.budget)
    assert object_result.diagnostics["human_prediction_source"] == (
        "object-style human recast"
    )
    assert invalid_result.passes is True
    assert invalid_result.predicted is None
    assert invalid_result.ratio is None
    assert invalid_result.experimental == pytest.approx(constraint.anchor.value)
    assert invalid_result.budget == pytest.approx(constraint.anchor.budget)
    assert invalid_result.diagnostics["exception_type"] == "ValueError"
    assert "NEEDS-HUMAN-PHYSICS" in invalid_result.diagnostics[
        "needs_human_physics_sm_eft"
    ]
    assert "NEEDS-HUMAN-PHYSICS" in invalid_result.diagnostics[
        "needs_human_physics_rs_matching"
    ]


def test_evaluate_is_pure_and_deterministic():
    raw = _raw_prediction(0.225)
    point = ParameterPoint(raw=raw, extras={})
    before = dict(raw)
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    assert raw == before
    for value in (first.predicted, first.experimental, first.ratio, first.budget):
        assert isinstance(value, float)
        assert math.isfinite(value)
    assert first.sm_prediction is None
