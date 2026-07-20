"""Production tests for L023 (muon-neutrino trident production)."""

from __future__ import annotations

import math
from dataclasses import replace
from pathlib import Path
from statistics import NormalDist

import pytest
import yaml

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import point_builder
from flavor_catalog_constraints.base import ConstraintProtocol, Severity
from flavor_catalog_constraints.primary.charged_lepton import L023 as l023_module

_PID = "L023"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = (
    _REPO_ROOT / "flavor_catalog" / "processes" / "charged_lepton" / "L023.yaml"
)
_CCFR_CL_VALUE_ID = "Altmannshofer2014:L023:ccfr_exclusion_cl"


def _yaml_pdg_block():
    with open(_SIDECAR) as handle:
        return yaml.safe_load(handle)["pdg_or_equivalent"]


def _find_measurement(experiment: str):
    return next(
        entry
        for entry in _yaml_pdg_block()["measured_observables"]
        if entry["experiment"] == experiment
    )


def _find_value(value_id: str):
    return next(
        entry for entry in _yaml_pdg_block()["values"] if entry["value_id"] == value_id
    )


def _point_from_mapping(mapping):
    return point_builder.make_point(lepton_mass_basis_couplings=dict(mapping))


def _manual_trident_ratio(delta_c_vector: float, delta_c_axial: float, sin2: float):
    c_vector_sm = 1.0 + 4.0 * sin2
    c_axial_sm = 1.0
    return (
        (c_vector_sm + delta_c_vector) ** 2
        + (c_axial_sm + delta_c_axial) ** 2
    ) / (c_vector_sm**2 + c_axial_sm**2)


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.HARD
    assert constraint.family == "charged_lepton"
    assert constraint.observable == "sigma(nu_mu N -> nu_mu N mu+ mu-) / sigma_SM"


def test_anchor_matches_yaml_and_budget_band():
    constraint = fcc.get(_PID)
    ccfr = _find_measurement("CCFR")
    charmii = _find_measurement("CHARM-II")
    nutev = _find_measurement("NuTeV")
    cl_entry = _find_value(_CCFR_CL_VALUE_ID)
    cl_fraction = float(cl_entry["value"]) / 100.0
    gaussian_z = NormalDist().inv_cdf(0.5 + 0.5 * cl_fraction)
    expected_budget = gaussian_z * float(ccfr["uncertainty"])

    assert len(constraint.anchor.measurements) == 3
    assert constraint.anchor.active_measurement.experiment == "CCFR"
    assert constraint.anchor.active_measurement.value == pytest.approx(ccfr["value"])
    assert constraint.anchor.active_measurement.uncertainty == pytest.approx(
        ccfr["uncertainty"]
    )
    assert constraint.anchor.active_measurement.source_url == ccfr["source_url"]
    assert constraint.anchor.ccfr_exclusion_cl.value_percent == pytest.approx(
        cl_entry["value"]
    )
    assert constraint.anchor.ccfr_exclusion_cl.source_url == cl_entry["source_url"]
    assert constraint.anchor.budget_band.gaussian_z == pytest.approx(gaussian_z)
    assert constraint.anchor.budget == pytest.approx(expected_budget)
    assert constraint.anchor.budget == pytest.approx(0.548789915671215)
    assert constraint.anchor.budget_band.upper_edge == pytest.approx(
        ccfr["value"] + expected_budget
    )
    assert constraint.anchor.budget_band.lower_edge == pytest.approx(
        ccfr["value"] - expected_budget
    )
    assert l023_module._find_measurement(
        constraint.anchor.measurements,
        "CHARM-II",
        process_id=_PID,
    ).value == pytest.approx(charmii["value"])
    nutev_anchor = l023_module._find_measurement(
        constraint.anchor.measurements,
        "NuTeV",
        process_id=_PID,
    )
    assert nutev_anchor.uncertainty_upper == pytest.approx(nutev["uncertainty_plus"])
    assert nutev_anchor.uncertainty_lower == pytest.approx(nutev["uncertainty_minus"])


def test_anchor_loading_routes_list_entries_through_load_anchor(monkeypatch):
    calls: list[tuple[str, ...]] = []
    original_load_anchor = l023_module.load_anchor

    def spy_load_anchor(*args, **kwargs):
        calls.append(tuple(kwargs["candidates"]))
        return original_load_anchor(*args, **kwargs)

    monkeypatch.setattr(l023_module, "load_anchor", spy_load_anchor)
    anchor = l023_module._load_l023_anchor(_PID)

    assert calls == [
        ("measured_observables[0]",),
        ("measured_observables[1]",),
        ("measured_observables[2]",),
        ("values[0]",),
    ]
    assert anchor.value == pytest.approx(fcc.get(_PID).anchor.value)

    with pytest.raises(fcc.AnchorError, match="not a mapping"):
        fcc.load_anchor(
            _PID,
            family="charged_lepton",
            candidates=("measured_observables",),
        )


def test_anchor_loading_rejects_mismatched_load_anchor_block_key(monkeypatch):
    original_load_anchor = l023_module.load_anchor

    def mismatched_load_anchor(*args, **kwargs):
        anchor = original_load_anchor(*args, **kwargs)
        return replace(anchor, block_key="wrong_block")

    monkeypatch.setattr(l023_module, "load_anchor", mismatched_load_anchor)
    with pytest.raises(fcc.AnchorError, match="load_anchor selected 'wrong_block'"):
        l023_module._load_l023_anchor(_PID)


def test_evaluate_without_input_degrades_gracefully():
    result = fcc.get(_PID).evaluate(point_builder.empty_point())

    assert result.process_id == _PID
    assert result.passes is True
    assert result.predicted is None
    assert result.ratio is None
    assert result.sm_prediction == pytest.approx(1.0)
    assert result.experimental == pytest.approx(fcc.get(_PID).anchor.value)
    assert result.budget == pytest.approx(fcc.get(_PID).anchor.budget)
    assert result.notes.startswith("NOT EVALUATED -")
    assert result.diagnostics["evaluated"] is False
    assert result.diagnostics["missing_extra"] == "lepton_mass_basis_couplings"
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]


def test_invalid_lepton_input_is_unevaluated_not_real_pass():
    result = fcc.get(_PID).evaluate(_point_from_mapping({"source": "empty"}))

    assert result.passes is False
    assert result.predicted is None
    assert result.ratio is None
    assert result.notes.startswith("NOT EVALUATED -")
    assert result.diagnostics["evaluated"] is True
    assert result.diagnostics["invalid_extra"] == "lepton_mass_basis_couplings"
    assert result.diagnostics["exception_type"] in {"KeyError", "ValueError"}


def test_sm_limit_and_proxy_numerics_match_independent_formula():
    constraint = fcc.get(_PID)
    sm_result = constraint.evaluate(_point_from_mapping({"trident_delta_c_vector": 0.0}))
    shifted = constraint.evaluate(
        _point_from_mapping(
            {
                "trident_delta_c_vector": 0.2,
                "trident_delta_c_axial": -0.1,
                "source": "L023 test direct effective-shift proxy",
            }
        )
    )

    assert sm_result.predicted == pytest.approx(1.0)
    assert sm_result.predicted == pytest.approx(sm_result.sm_prediction)
    assert sm_result.passes is True
    sin2 = shifted.diagnostics["sin2_theta_w"]
    expected = _manual_trident_ratio(0.2, -0.1, sin2)
    expected_pull = expected - constraint.anchor.value
    expected_ratio = abs(expected_pull) / constraint.anchor.budget

    assert shifted.predicted == pytest.approx(expected)
    assert shifted.predicted == pytest.approx(1.1317562936528025)
    assert shifted.ratio == pytest.approx(expected_ratio)
    assert shifted.diagnostics["delta_c_vector"] == pytest.approx(0.2)
    assert shifted.diagnostics["delta_c_axial"] == pytest.approx(-0.1)
    assert "NEEDS-HUMAN-PHYSICS" in shifted.diagnostics["needs_human_physics"]


@pytest.mark.parametrize(
    ("lepton", "expected_pass"),
    [
        ({"trident_delta_c_vector": 0.1}, True),
        ({"trident_delta_c_vector": 1.0}, False),
    ],
)
def test_safe_point_passes_and_large_np_point_fails(lepton, expected_pass):
    constraint = fcc.get(_PID)
    result = constraint.evaluate(_point_from_mapping(lepton))
    expected = _manual_trident_ratio(
        lepton["trident_delta_c_vector"],
        0.0,
        result.diagnostics["sin2_theta_w"],
    )

    assert result.passes is expected_pass
    assert result.predicted == pytest.approx(expected)
    assert result.sm_prediction == pytest.approx(1.0)
    if expected_pass:
        assert result.ratio < 1.0
    else:
        assert result.ratio > 1.0


def test_heavy_zprime_proxy_accepts_kk_ew_mass_override():
    constraint = fcc.get(_PID)
    default_result = constraint.evaluate(
        _point_from_mapping(
            {
                "g_nu_mu": 1.0,
                "g_mu_vector": 1.0,
                "m_kk_gev": 3000.0,
                "source": "L023 zprime proxy",
            }
        )
    )
    override_result = constraint.evaluate(
        point_builder.make_point(
            lepton_mass_basis_couplings={
                "g_nu_mu": 1.0,
                "g_mu_vector": 1.0,
                "m_kk_gev": 3000.0,
                "source": "L023 zprime proxy with override",
            },
            kk_ew_mass_gev=6000.0,
        )
    )

    assert default_result.diagnostics["delta_c_vector"] == pytest.approx(
        2.0 * 246.0**2 / 3000.0**2
    )
    assert override_result.diagnostics["delta_c_vector"] == pytest.approx(
        default_result.diagnostics["delta_c_vector"] / 4.0
    )
    assert override_result.diagnostics["kk_ew_mass_extra_used"] is True
    assert override_result.diagnostics["kk_ew_mass_override_used"] is True


def test_evaluate_runs_end_to_end_with_real_finite_fields():
    result = fcc.get(_PID).evaluate(
        _point_from_mapping(
            {
                "trident_delta_c_vector": 0.2,
                "trident_delta_c_axial": 0.05,
            }
        )
    )

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
    for key in (
        "sin2_theta_w",
        "c_vector_sm",
        "c_axial_sm",
        "c_vector_total",
        "c_axial_total",
        "delta_c_vector",
        "delta_c_axial",
        "budget_gaussian_z",
        "budget_hard_veto_upper",
        "budget_hard_veto_lower",
        "pull_to_active_measurement",
    ):
        assert isinstance(result.diagnostics[key], float)
        assert math.isfinite(result.diagnostics[key])


def test_evaluate_is_pure_and_deterministic():
    lepton = {
        "trident_delta_c_vector": 0.2,
        "trident_delta_c_axial": -0.1,
        "source": "L023 purity test",
    }
    point = _point_from_mapping(lepton)
    before_extras = dict(point.extras)
    before_lepton = dict(point.extras["lepton_mass_basis_couplings"])
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    assert point.extras == before_extras
    assert point.extras["lepton_mass_basis_couplings"] == before_lepton
