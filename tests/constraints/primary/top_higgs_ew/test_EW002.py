"""Production tests for EW002 (first-row CKM unitarity)."""

from __future__ import annotations

import math
from pathlib import Path

import pytest
import yaml

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import point_builder
from flavor_catalog_constraints.anchors import Anchor, AnchorError
from flavor_catalog_constraints.base import ConstraintProtocol, Severity
from flavor_catalog_constraints.primary.top_higgs_ew import EW002 as ew002_module
from tests.constraints.charged_current_phase5b_helpers import (
    charged_with_epsilon,
    sample_charged_point,
    universal_charged_point,
)

_PID = "EW002"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "top_higgs_ew" / "EW002.yaml"


def _yaml_entries_by_id():
    with open(_SIDECAR) as handle:
        entries = yaml.safe_load(handle)["pdg_or_equivalent"]
    return {entry["value_id"]: entry for entry in entries if "value_id" in entry}


def _controlled_constraint(
    *,
    vud: float,
    vus: float,
    vub: float,
    budget: float,
) -> tuple[ew002_module.Constraint, float]:
    first_row_sum = float(vud * vud + vus * vus + vub * vub)
    constraint = object.__new__(ew002_module.Constraint)
    constraint.anchor = ew002_module.EW002Anchor(
        vud=Anchor(
            process_id=_PID,
            block_key="test:vud",
            value=float(vud),
            uncertainty=1.0e-6,
        ),
        vus=Anchor(
            process_id=_PID,
            block_key="test:vus",
            value=float(vus),
            uncertainty=1.0e-6,
        ),
        first_row_sum=ew002_module.FirstRowSumAnchor(
            anchor=Anchor(
                process_id=_PID,
                block_key="test:first_row_sum",
                value=first_row_sum,
                uncertainty=None,
            ),
            uncertainty_components={
                "Vud_squared": 1.0e-6,
                "Vus_squared": 1.0e-6,
                "combined_quoted_elsewhere": float(budget),
            },
            budget=float(budget),
            budget_source="test controlled CKM first-row budget",
        ),
    )
    return constraint, first_row_sum


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.SOFT
    assert constraint.family == "top_higgs_ew"
    assert constraint.observable == "Delta_CKM first-row CKM unitarity"


def test_anchor_matches_yaml_and_loud_fail_probe():
    constraint = fcc.get(_PID)
    entries = _yaml_entries_by_id()
    vud = entries["PDG2025:EW002:Vud_superallowed"]
    vus = entries["PDG2025:EW002:Vus_kaon_average"]
    first_row = entries["PDG2025:EW002:first_row_sum"]

    assert constraint.anchor.vud.value == pytest.approx(vud["value"])
    assert constraint.anchor.vud.uncertainty == pytest.approx(vud["uncertainty"])
    assert constraint.anchor.vud.source_url == vud["source_url"]
    assert constraint.anchor.vus.value == pytest.approx(vus["value"])
    assert constraint.anchor.vus.uncertainty == pytest.approx(vus["uncertainty"])
    assert constraint.anchor.vus.source_url == vus["source_url"]
    assert constraint.anchor.first_row_sum.value == pytest.approx(first_row["value"])
    assert constraint.anchor.first_row_sum.source_url == first_row["source_url"]
    assert constraint.anchor.first_row_sum.budget == pytest.approx(
        first_row["uncertainties"]["combined_quoted_elsewhere"]
    )
    for key in ("Vud_squared", "Vus_squared", "combined_quoted_elsewhere"):
        assert constraint.anchor.first_row_sum.uncertainty_components[key] == (
            pytest.approx(first_row["uncertainties"][key])
        )

    with pytest.raises(AnchorError):
        ew002_module._load_value_anchor(_PID, "PDG2025:EW002:not_present")


def test_list_anchors_route_through_scaffold_load_anchor_and_fail_loudly(monkeypatch):
    calls: list[tuple[str, ...]] = []
    original_load_anchor = ew002_module.load_anchor

    def spy_load_anchor(*args, **kwargs):
        calls.append(tuple(kwargs["candidates"]))
        return original_load_anchor(*args, **kwargs)

    monkeypatch.setattr(ew002_module, "load_anchor", spy_load_anchor)
    anchor = ew002_module._load_ew002_anchor(_PID)

    assert calls == [
        ("pdg_or_equivalent[0]",),
        ("pdg_or_equivalent[3]",),
        ("pdg_or_equivalent[4]",),
    ]
    assert anchor.value == pytest.approx(fcc.get(_PID).anchor.value)

    with pytest.raises(AnchorError):
        ew002_module._load_value_anchor(_PID, "PDG2025:EW002:not_present")

    def mismatched_load_anchor(*args, **kwargs):
        return Anchor(
            process_id=_PID,
            block_key="pdg_or_equivalent[99]",
            value=1.0,
            uncertainty=0.1,
        )

    monkeypatch.setattr(ew002_module, "load_anchor", mismatched_load_anchor)
    with pytest.raises(AnchorError, match="load_anchor selected"):
        ew002_module._load_value_anchor(_PID, "PDG2025:EW002:Vud_superallowed")


def test_numerical_delta_ckm_matches_independent_yaml_recomputation():
    constraint = fcc.get(_PID)
    entries = _yaml_entries_by_id()
    first_row = entries["PDG2025:EW002:first_row_sum"]
    budget = first_row["uncertainties"]["combined_quoted_elsewhere"]
    expected_delta = float(first_row["value"] - 1.0)
    expected_pull = abs(expected_delta) / budget

    result = constraint.evaluate(universal_charged_point())

    assert result.predicted == pytest.approx(1.0)
    assert result.sm_prediction == pytest.approx(1.0)
    assert result.experimental == pytest.approx(first_row["value"])
    assert result.budget == pytest.approx(0.0007)
    assert result.diagnostics["sm_vs_data_delta_ckm"] == pytest.approx(expected_delta)
    assert result.diagnostics["np_shift_delta_ckm"] == pytest.approx(0.0, abs=1.0e-18)
    assert result.ratio == pytest.approx(expected_pull)
    assert result.ratio == pytest.approx(2.4285714285714284)


def test_vud_vus_diagnostic_arithmetic_is_independent():
    constraint = fcc.get(_PID)
    entries = _yaml_entries_by_id()
    vud = entries["PDG2025:EW002:Vud_superallowed"]
    vus = entries["PDG2025:EW002:Vus_kaon_average"]
    vud_squared = float(vud["value"] * vud["value"])
    vus_squared = float(vus["value"] * vus["value"])
    sum_without_vub = float(vud_squared + vus_squared)
    propagated_without_vub = math.sqrt(
        (2.0 * vud["value"] * vud["uncertainty"]) ** 2
        + (2.0 * vus["value"] * vus["uncertainty"]) ** 2
    )
    result = constraint.evaluate(point_builder.empty_point())

    assert result.diagnostics["vud_squared_from_yaml_value"] == pytest.approx(
        vud_squared
    )
    assert result.diagnostics["vus_squared_from_yaml_value"] == pytest.approx(
        vus_squared
    )
    assert result.diagnostics["sum_vud_vus_squared_without_vub"] == pytest.approx(
        sum_without_vub
    )
    assert result.diagnostics["sum_vud_vus_squared_without_vub"] == pytest.approx(
        0.998348245
    )
    assert result.diagnostics["pdg_first_row_sum_minus_vud_vus_squared"] == (
        pytest.approx(entries["PDG2025:EW002:first_row_sum"]["value"] - sum_without_vub)
    )
    manual_propagated_without_vub = math.sqrt(
        (2.0 * vud["value"] * vud["uncertainty"]) ** 2
        + (2.0 * vus["value"] * vus["uncertainty"]) ** 2
        + (2.0 * 1.0e-12 * 1.0e-12) ** 2
    )
    assert manual_propagated_without_vub == pytest.approx(propagated_without_vub)


def test_soft_tension_reports_np_gap_and_real_finite_fields():
    result = fcc.get(_PID).evaluate(sample_charged_point())

    assert result.process_id == _PID
    assert result.passes is False
    for value in (
        result.predicted,
        result.ratio,
        result.budget,
        result.sm_prediction,
        result.experimental,
    ):
        assert isinstance(value, float)
        assert math.isfinite(value)
    charged = sample_charged_point().extras["rs_charged_current"]
    expected_delta_np = (
        2.0
        * fcc.get(_PID).anchor.vud.value**2
        * charged.epsilon[0, 0, 0].real
        + 2.0
        * fcc.get(_PID).anchor.vus.value**2
        * (0.5 * (charged.epsilon[0, 1, 0] + charged.epsilon[0, 1, 1])).real
    )
    assert result.diagnostics["evaluated"] is True
    assert result.diagnostics["np_shift_delta_ckm"] == pytest.approx(expected_delta_np)
    assert result.diagnostics["np_shift_delta_ckm"] == pytest.approx(
        5.254734244800866e-05
    )
    assert result.predicted == pytest.approx(1.000052547342448)
    assert result.diagnostics["epsilon_us_light_average_e_mu"] == pytest.approx(
        0.5 * (charged.epsilon[0, 1, 0] + charged.epsilon[0, 1, 1])
    )
    assert result.diagnostics["vub_value_block_in_ew002_yaml"] is False
    assert "mass proxy" not in result.notes


def test_absent_charged_current_degrades_non_vetoing():
    result = fcc.get(_PID).evaluate(point_builder.empty_point())

    assert result.passes is True
    assert result.predicted is None
    assert result.ratio is None
    assert result.diagnostics["evaluated"] is False
    assert result.diagnostics["missing_extra"] == "rs_charged_current"
    assert result.diagnostics["sm_vs_data_pull_sigma"] == pytest.approx(
        2.4285714285714284
    )


def test_safe_unitarity_sum_passes_and_observed_soft_tension_fails():
    safe_constraint, safe_sum = _controlled_constraint(
        vud=0.8,
        vus=0.6,
        vub=0.0,
        budget=0.01,
    )
    excluded_constraint, excluded_sum = _controlled_constraint(
        vud=0.8,
        vus=0.5,
        vub=0.0,
        budget=0.02,
    )

    safe = safe_constraint.evaluate(universal_charged_point())
    excluded = excluded_constraint.evaluate(universal_charged_point())
    expected_safe_delta = float(0.8 * 0.8 + 0.6 * 0.6 + 0.0 * 0.0 - 1.0)
    expected_excluded_delta = float(1.0 - excluded_sum)
    expected_excluded_pull = abs(expected_excluded_delta) / 0.02

    assert safe_sum == pytest.approx(1.0)
    assert safe.passes is True
    assert safe.predicted == pytest.approx(safe_sum)
    assert safe.diagnostics["delta_ckm"] == pytest.approx(expected_safe_delta)
    assert safe.ratio == pytest.approx(0.0)

    assert excluded_sum == pytest.approx(0.89)
    assert excluded.passes is False
    assert excluded.predicted == pytest.approx(1.0)
    assert excluded.diagnostics["delta_ckm"] == pytest.approx(expected_excluded_delta)
    assert excluded.ratio == pytest.approx(expected_excluded_pull)

    shifted = charged_with_epsilon(
        universal_charged_point().extras["rs_charged_current"],
        {(0, 0, 0): 1.0e-3 + 0.0j},
    )
    shifted_result = safe_constraint.evaluate(
        point_builder.make_point(rs_charged_current=shifted)
    )
    assert shifted_result.diagnostics["epsilon_ud_e"] == pytest.approx(1.0e-3 + 0.0j)
    assert shifted_result.diagnostics["np_shift_delta_ckm"] == pytest.approx(
        2.0 * 0.8 * 0.8 * 1.0e-3
    )


def test_evaluate_is_pure_and_deterministic():
    point = point_builder.empty_point()
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    assert point.extras == {}
