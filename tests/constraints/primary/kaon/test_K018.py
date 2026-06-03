"""Production tests for K018 (K_l3 |V_us| f_+(0) extraction)."""

from __future__ import annotations

import math
from pathlib import Path

import pytest
import yaml

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import point_builder
from flavor_catalog_constraints.anchors import Anchor, AnchorError
from flavor_catalog_constraints.base import ConstraintProtocol, Severity
from flavor_catalog_constraints.primary.kaon import K018 as k018_module
from tests.constraints.charged_current_phase5b_helpers import (
    charged_with_epsilon,
    sample_charged_point,
    universal_charged_point,
)
from quarkConstraints.ckm_extraction import (
    extract_vus_from_kl3,
    vus_consistency_pull,
)

_PID = "K018"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "kaon" / "K018.yaml"
_FPLUS_VUS_VALUE_ID = "PDG2025:K018:average_fplus_vus"
_FPLUS_VALUE_ID = "FLAG2024:K018:fplus_Nf211"
_PDG_DERIVED_VUS_VALUE_ID = "PDG2025:K018:Vus_from_Kl3"
_FLAG_DERIVED_VUS_VALUE_ID = "FLAG2024:K018:Vus_from_fplus_Nf211"


def _entries_by_id():
    with open(_SIDECAR) as handle:
        entries = yaml.safe_load(handle)["pdg_or_equivalent"]
    return {entry["value_id"]: entry for entry in entries if "value_id" in entry}


def _value_anchor(
    *,
    value_id: str,
    block_key: str,
    value: float,
    uncertainty: float,
    uncertainty_total: float | None = None,
    uncertainty_breakdown: dict[str, float] | None = None,
) -> k018_module.K018ValueAnchor:
    return k018_module.K018ValueAnchor(
        anchor=Anchor(
            process_id=_PID,
            block_key=block_key,
            value=float(value),
            uncertainty=float(uncertainty),
            units="dimensionless",
        ),
        value_id=value_id,
        entry_index=0,
        uncertainty_total=uncertainty_total,
        uncertainty_breakdown={} if uncertainty_breakdown is None else uncertainty_breakdown,
        nf=None,
        value_summary=None,
    )


def _controlled_constraint(
    *,
    fplus_vus: float,
    fplus: float,
    reference_vus: float,
    budget: float,
) -> k018_module.Constraint:
    constraint = object.__new__(k018_module.Constraint)
    constraint.anchor = k018_module.K018Anchor(
        fplus_vus=_value_anchor(
            value_id=_FPLUS_VUS_VALUE_ID,
            block_key="test:fplus_vus",
            value=fplus_vus,
            uncertainty=1.0e-5,
        ),
        fplus_zero=_value_anchor(
            value_id=_FPLUS_VALUE_ID,
            block_key="test:fplus",
            value=fplus,
            uncertainty=1.0e-5,
        ),
        pdg_derived_vus=_value_anchor(
            value_id=_PDG_DERIVED_VUS_VALUE_ID,
            block_key="test:pdg_vus",
            value=reference_vus,
            uncertainty=budget,
            uncertainty_total=budget,
            uncertainty_breakdown={"radiative_isospin": 0.0},
        ),
        flag_derived_vus=_value_anchor(
            value_id=_FLAG_DERIVED_VUS_VALUE_ID,
            block_key="test:flag_vus",
            value=reference_vus,
            uncertainty=budget,
        ),
        budget=float(budget),
        budget_source="test controlled K018 budget",
    )
    return constraint


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.SOFT
    assert constraint.family == "kaon"
    assert constraint.observable == "K_l3 |V_us| f_+(0) extraction"


def test_anchor_matches_yaml_and_loud_fail_probe():
    constraint = fcc.get(_PID)
    entries = _entries_by_id()
    fplus_vus = entries[_FPLUS_VUS_VALUE_ID]
    fplus = entries[_FPLUS_VALUE_ID]
    pdg_vus = entries[_PDG_DERIVED_VUS_VALUE_ID]
    flag_vus = entries[_FLAG_DERIVED_VUS_VALUE_ID]

    assert constraint.anchor.fplus_vus.value == pytest.approx(fplus_vus["value"])
    assert constraint.anchor.fplus_vus.uncertainty == pytest.approx(
        fplus_vus["uncertainty"]
    )
    assert constraint.anchor.fplus_vus.source_url == fplus_vus["source_url"]
    assert constraint.anchor.fplus_zero.value == pytest.approx(fplus["value"])
    assert constraint.anchor.fplus_zero.uncertainty == pytest.approx(
        fplus["uncertainty"]
    )
    assert constraint.anchor.fplus_zero.source_url == fplus["source_url"]
    assert constraint.anchor.pdg_derived_vus.value == pytest.approx(pdg_vus["value"])
    assert constraint.anchor.pdg_derived_vus.uncertainty_total == pytest.approx(
        pdg_vus["uncertainty_total"]
    )
    assert constraint.anchor.pdg_derived_vus.uncertainty_breakdown == {
        key: pytest.approx(value)
        for key, value in pdg_vus["uncertainty_breakdown"].items()
    }
    assert constraint.anchor.flag_derived_vus.value == pytest.approx(flag_vus["value"])
    assert constraint.anchor.flag_derived_vus.source_url == flag_vus["source_url"]
    assert constraint.anchor.budget == pytest.approx(pdg_vus["uncertainty_total"])

    with pytest.raises(AnchorError):
        k018_module._load_value_anchor("PDG2025:K018:not_present", process_id=_PID)


def test_list_anchors_route_through_scaffold_load_anchor(monkeypatch):
    calls: list[tuple[str, ...]] = []
    original_load_anchor = k018_module.load_anchor

    def spy_load_anchor(*args, **kwargs):
        calls.append(tuple(kwargs["candidates"]))
        return original_load_anchor(*args, **kwargs)

    monkeypatch.setattr(k018_module, "load_anchor", spy_load_anchor)
    anchor = k018_module._load_k018_anchor(_PID)

    assert calls == [
        ("pdg_or_equivalent[0]",),
        ("pdg_or_equivalent[10]",),
        ("pdg_or_equivalent[9]",),
        ("pdg_or_equivalent[11]",),
    ]
    assert anchor.value == pytest.approx(fcc.get(_PID).anchor.value)


def test_numerical_extraction_matches_core_ckm_recomputation():
    constraint = fcc.get(_PID)
    radiative_isospin = constraint.anchor.pdg_derived_vus.uncertainty_breakdown[
        "radiative_isospin"
    ]
    extraction = extract_vus_from_kl3(
        fplus_vus_product=constraint.anchor.fplus_vus.value,
        fplus_vus_uncertainty=constraint.anchor.fplus_vus.uncertainty,
        fplus_zero=constraint.anchor.fplus_zero.value,
        fplus_zero_uncertainty=constraint.anchor.fplus_zero.uncertainty,
        radiative_isospin_uncertainty=radiative_isospin,
    )
    consistency = vus_consistency_pull(
        extracted_vus=extraction.vus,
        reference_vus=constraint.anchor.pdg_derived_vus.value,
        budget=constraint.anchor.budget,
    )

    result = constraint.evaluate(universal_charged_point())

    assert result.predicted == pytest.approx(extraction.vus)
    assert result.predicted == pytest.approx(0.22330377397401527)
    assert result.experimental == pytest.approx(consistency.reference_vus)
    assert result.budget == pytest.approx(consistency.budget)
    assert result.ratio == pytest.approx(consistency.pull_sigma)
    assert result.ratio == pytest.approx(0.0071207056891883)
    assert result.diagnostics["extracted_vus_uncertainty"] == pytest.approx(
        extraction.total_uncertainty
    )
    assert result.diagnostics["extracted_vus_uncertainty_components"][
        "experimental"
    ] == pytest.approx(extraction.experimental_uncertainty)
    assert result.diagnostics["extracted_vus_uncertainty_components"][
        "lattice"
    ] == pytest.approx(extraction.lattice_uncertainty)
    assert result.diagnostics["extracted_vus_uncertainty_components"][
        "radiative_isospin"
    ] == pytest.approx(extraction.radiative_isospin_uncertainty)
    assert result.diagnostics["delta_vus_vs_pdg_derived"] == pytest.approx(
        consistency.delta_vus
    )
    assert result.diagnostics["epsilon_us_light_average_e_mu"] == 0.0j
    assert result.diagnostics["np_shift_vus"] == pytest.approx(0.0)
    assert result.diagnostics["flag_derived_vus"] == pytest.approx(
        constraint.anchor.flag_derived_vus.value
    )
    assert result.passes is consistency.passes


def test_rigorous_epsilon_path_shifts_vus_and_keeps_real_finite_fields():
    constraint = fcc.get(_PID)
    point = sample_charged_point()
    charged = point.extras["rs_charged_current"]
    epsilon_us = 0.5 * (charged.epsilon[0, 1, 0] + charged.epsilon[0, 1, 1])
    sm_vus = 0.22330377397401527
    expected = sm_vus * abs(1.0 + epsilon_us)

    result = constraint.evaluate(point)

    for value in (result.predicted, result.ratio, result.budget, result.experimental):
        assert isinstance(value, float)
        assert math.isfinite(value)
    assert result.sm_prediction is None
    assert result.diagnostics["evaluated"] is True
    assert result.diagnostics["parameter_point_used"] is True
    assert result.predicted == pytest.approx(expected)
    assert result.predicted == pytest.approx(0.2233096753995752)
    assert result.diagnostics["epsilon_us_light_average_e_mu"] == pytest.approx(
        epsilon_us
    )
    assert result.diagnostics["np_shift_vus"] == pytest.approx(expected - sm_vus)
    assert "PARTIAL" in result.diagnostics["radiative_isospin_mode_weight_status"]
    assert "mass proxy" not in result.notes


def test_absent_charged_current_degrades_non_vetoing():
    result = fcc.get(_PID).evaluate(point_builder.empty_point())

    assert result.passes is True
    assert result.predicted is None
    assert result.ratio is None
    assert result.diagnostics["evaluated"] is False
    assert result.diagnostics["missing_extra"] == "rs_charged_current"
    assert result.diagnostics["np_shift_vus"] == pytest.approx(0.0)


def test_safe_consistency_point_passes_and_controlled_tension_fails():
    safe = _controlled_constraint(
        fplus_vus=0.2,
        fplus=1.0,
        reference_vus=0.2,
        budget=0.01,
    ).evaluate(universal_charged_point())
    excluded = _controlled_constraint(
        fplus_vus=0.2,
        fplus=1.0,
        reference_vus=0.16,
        budget=0.01,
    ).evaluate(universal_charged_point())

    assert safe.passes is True
    assert safe.predicted == pytest.approx(0.2)
    assert safe.ratio == pytest.approx(0.0)
    assert excluded.passes is False
    assert excluded.predicted == pytest.approx(0.2)
    assert excluded.ratio == pytest.approx(4.0)

    shifted = charged_with_epsilon(
        universal_charged_point().extras["rs_charged_current"],
        {(0, 1, 0): 0.1 + 0.0j, (0, 1, 1): 0.1 + 0.0j},
    )
    shifted_result = _controlled_constraint(
        fplus_vus=0.2,
        fplus=1.0,
        reference_vus=0.2,
        budget=0.01,
    ).evaluate(point_builder.make_point(rs_charged_current=shifted))
    assert shifted_result.predicted == pytest.approx(0.22)
    assert shifted_result.passes is False


def test_evaluate_is_pure_and_deterministic():
    point = sample_charged_point()
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)
    empty = constraint.evaluate(point_builder.empty_point())

    assert first == second
    assert empty.diagnostics["evaluated"] is False
    assert point.get_extra("rs_charged_current") is not None
