"""Production tests for B025 (R_D in B -> D tau nu)."""

from __future__ import annotations

import math
from pathlib import Path

import pytest
import yaml

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import point_builder
from flavor_catalog_constraints.anchors import AnchorError
from flavor_catalog_constraints.base import ConstraintProtocol, Severity
from flavor_catalog_constraints.primary.beauty import B025 as b025_module
from tests.constraints.charged_current_phase5b_helpers import (
    charged_with_epsilon,
    sample_charged_point,
    universal_charged_point,
)

_PID = "B025"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "beauty" / "B025.yaml"


def _yaml_pdg_block():
    with open(_SIDECAR) as handle:
        return yaml.safe_load(handle)["pdg_or_equivalent"]


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.HARD
    assert constraint.family == "beauty"
    assert constraint.observable == "R_D = Gamma(B -> D tau nu) / Gamma(B -> D l nu)"


def test_anchor_matches_yaml_and_budget_band():
    constraint = fcc.get(_PID)
    pdg = _yaml_pdg_block()
    exp = pdg["canonical_average"]
    sm = pdg["sm_reference_used_by_hflav"]
    joint = pdg["joint_fit_context"]
    central = abs(float(exp["value"]) - float(sm["rd_value"]))
    combined = math.sqrt(float(exp["uncertainty"]) ** 2 + float(sm["rd_uncertainty"]) ** 2)

    assert constraint.anchor.experimental.value == pytest.approx(exp["value"])
    assert constraint.anchor.experimental.uncertainty == pytest.approx(exp["uncertainty"])
    assert constraint.anchor.experimental.source_url == exp["source_url"]
    assert constraint.anchor.standard_model.value == pytest.approx(sm["rd_value"])
    assert constraint.anchor.standard_model.uncertainty == pytest.approx(
        sm["rd_uncertainty"]
    )
    assert constraint.anchor.standard_model.source_url == sm["source_url"]
    assert constraint.anchor.joint_fit_context.rdstar_value == pytest.approx(
        joint["rdstar_value"]
    )
    assert constraint.anchor.joint_fit_context.correlation_rd_rdstar == pytest.approx(
        joint["correlation_rd_rdstar"]
    )
    assert constraint.anchor.budget_band.central_residual == pytest.approx(central)
    assert constraint.anchor.budget_band.combined_sigma == pytest.approx(combined)
    assert constraint.anchor.budget == pytest.approx(central + combined)
    assert constraint.anchor.budget == pytest.approx(0.08633105012119291)
    assert "B025.yaml" in constraint.anchor.budget_band.source


def test_anchor_loads_through_scaffold_and_fails_loudly(monkeypatch):
    calls: list[tuple[str, ...]] = []
    original_load_anchor = b025_module.load_anchor

    def spy_load_anchor(*args, **kwargs):
        calls.append(tuple(kwargs["candidates"]))
        return original_load_anchor(*args, **kwargs)

    monkeypatch.setattr(b025_module, "load_anchor", spy_load_anchor)
    anchor = b025_module._load_b025_anchor(_PID)

    assert calls == [("canonical_average",), ("sm_reference_used_by_hflav",)]
    assert anchor.value == pytest.approx(fcc.get(_PID).anchor.value)
    assert anchor.sm_value == pytest.approx(0.296)

    def missing_load_anchor(*args, **kwargs):
        raise AnchorError("forced missing B025 anchor")

    monkeypatch.setattr(b025_module, "load_anchor", missing_load_anchor)
    with pytest.raises(AnchorError, match="forced missing B025 anchor"):
        b025_module._load_b025_anchor(_PID)

    with pytest.raises(AnchorError):
        original_load_anchor(
            _PID,
            family="beauty",
            candidates=("not_a_b025_block",),
        )


def test_evaluate_without_input_degrades_gracefully():
    constraint = fcc.get(_PID)
    result = constraint.evaluate(point_builder.empty_point())

    assert result.process_id == _PID
    assert result.passes is True
    assert result.predicted is None
    assert result.ratio is None
    assert result.experimental == pytest.approx(constraint.anchor.value)
    assert result.sm_prediction == pytest.approx(constraint.anchor.sm_value)
    assert result.budget == pytest.approx(constraint.anchor.budget)
    assert result.diagnostics["evaluated"] is False
    assert result.diagnostics["missing_extra"] == "rs_charged_current"
    assert result.diagnostics["matching_coverage"] == "PARTIAL"
    assert "PARTIAL" in result.diagnostics["b025_partial_matching_status"]
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]


def test_sm_limit_rd_matches_yaml_sm_reference_and_manual_budget():
    constraint = fcc.get(_PID)
    result = constraint.evaluate(universal_charged_point())
    expected_ratio = abs(constraint.anchor.sm_value - constraint.anchor.value) / (
        abs(constraint.anchor.value - constraint.anchor.sm_value)
        + math.sqrt(
            constraint.anchor.experimental.uncertainty**2
            + constraint.anchor.standard_model.uncertainty**2
        )
    )

    assert result.predicted == pytest.approx(0.296)
    assert result.predicted == pytest.approx(result.sm_prediction)
    assert result.predicted == pytest.approx(constraint.anchor.sm_value)
    assert result.ratio == pytest.approx(expected_ratio)
    assert result.ratio == pytest.approx(0.7181657110965687)
    assert result.passes is True
    assert result.diagnostics["sm_formula_minus_anchor"] == pytest.approx(0.0)
    assert result.diagnostics["epsilon_cb_tau"] == 0.0j
    assert result.diagnostics["epsilon_cb_light_average_e_mu"] == 0.0j
    assert result.diagnostics["np_shift_rd"] == pytest.approx(0.0)


def test_vector_lfu_prediction_matches_independent_formula():
    constraint = fcc.get(_PID)
    charged = charged_with_epsilon(
        universal_charged_point().extras["rs_charged_current"],
        {(1, 2, 2): 0.5 + 0.0j, (1, 2, 0): 0.0j, (1, 2, 1): 0.0j},
    )
    point = point_builder.make_point(rs_charged_current=charged)
    result = constraint.evaluate(point)
    manual = constraint.anchor.sm_value * abs(1.0 + 0.5) ** 2

    assert result.predicted == pytest.approx(manual)
    assert result.diagnostics["epsilon_cb_tau"] == pytest.approx(0.5 + 0.0j)
    assert result.diagnostics["epsilon_cb_light_average_e_mu"] == pytest.approx(0.0j)
    assert result.diagnostics["tau_rate_multiplier"] == pytest.approx(2.25)
    assert result.diagnostics["light_rate_multiplier"] == pytest.approx(1.0)
    assert result.diagnostics["np_shift_rd"] == pytest.approx(
        result.predicted - constraint.anchor.sm_value
    )
    assert result.passes is False
    assert "mass proxy" not in result.notes


def test_evaluate_runs_end_to_end_with_real_finite_fields_and_complex_diagnostics():
    result = fcc.get(_PID).evaluate(sample_charged_point())

    assert result.process_id == _PID
    for value in (
        result.predicted,
        result.ratio,
        result.budget,
        result.sm_prediction,
        result.experimental,
    ):
        assert isinstance(value, float)
        assert math.isfinite(value)
    for key in ("epsilon_cb_tau", "epsilon_cb_light_average_e_mu"):
        assert isinstance(result.diagnostics[key], complex)
    for key in (
        "kk_ew_mass_gev",
        "m_w_gev",
        "m_wprime_gev",
        "tau_rate_multiplier",
        "light_rate_multiplier",
        "np_shift_rd",
        "budget_combined_sigma",
    ):
        assert isinstance(result.diagnostics[key], float)
        assert math.isfinite(result.diagnostics[key])
    assert result.diagnostics["matching_coverage"] == "PARTIAL"
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]


@pytest.mark.parametrize(
    ("epsilon_cb_tau", "expected_pass"),
    [
        (0.0, True),
        (0.5, False),
    ],
)
def test_safe_point_passes_and_large_vector_shift_fails(
    epsilon_cb_tau: float,
    expected_pass: bool,
):
    constraint = fcc.get(_PID)
    charged = charged_with_epsilon(
        universal_charged_point().extras["rs_charged_current"],
        {(1, 2, 2): complex(epsilon_cb_tau), (1, 2, 0): 0.0j, (1, 2, 1): 0.0j},
    )
    result = constraint.evaluate(point_builder.make_point(rs_charged_current=charged))
    manual = constraint.anchor.sm_value * abs(1.0 + epsilon_cb_tau) ** 2

    assert result.passes is expected_pass
    assert result.predicted == pytest.approx(manual)
    assert result.ratio == pytest.approx(
        abs(result.predicted - constraint.anchor.value) / constraint.anchor.budget
    )
    assert (result.ratio < 1.0) is expected_pass


def test_universal_vector_shift_cancels_in_lfu_ratio_and_partial_flag_stays():
    constraint = fcc.get(_PID)
    charged = charged_with_epsilon(
        universal_charged_point().extras["rs_charged_current"],
        {
            (1, 2, 0): 0.2 + 0.0j,
            (1, 2, 1): 0.2 + 0.0j,
            (1, 2, 2): 0.2 + 0.0j,
        },
    )
    result = constraint.evaluate(point_builder.make_point(rs_charged_current=charged))

    assert result.predicted == pytest.approx(constraint.anchor.sm_value)
    assert result.diagnostics["tau_rate_multiplier"] == pytest.approx(1.44)
    assert result.diagnostics["light_rate_multiplier"] == pytest.approx(1.44)
    assert result.diagnostics["matching_coverage"] == "PARTIAL"
    assert "PARTIAL" in result.diagnostics["b025_partial_matching_status"]


def test_evaluate_is_pure_and_deterministic():
    point = sample_charged_point()
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    assert point.get_extra("rs_charged_current") is not None
