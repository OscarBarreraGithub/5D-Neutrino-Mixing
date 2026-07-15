"""Production tests for L006 (muonium-antimuonium conversion)."""

from __future__ import annotations

import math
from pathlib import Path

import pytest
import yaml

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import anchors, point_builder
from flavor_catalog_constraints.base import ConstraintProtocol, Severity

_PID = "L006"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = (
    _REPO_ROOT / "flavor_catalog" / "processes" / "charged_lepton" / "L006.yaml"
)


def _yaml_pdg_block():
    with open(_SIDECAR) as handle:
        return yaml.safe_load(handle)["pdg_or_equivalent"]


def _point_from_mapping(mapping):
    return point_builder.make_point(lepton_mass_basis_couplings=dict(mapping))


def _manual_probability_from_coupling(
    g_mmbar_over_gf: complex,
    *,
    probability_limit: float,
    coupling_limit_over_gf: float,
) -> float:
    return float(
        probability_limit * abs(g_mmbar_over_gf) ** 2 / coupling_limit_over_gf**2
    )


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.HARD
    assert constraint.family == "charged_lepton"
    assert constraint.observable == "P(M -> Mbar)"


def test_anchor_matches_yaml_and_loud_fail_probe():
    constraint = fcc.get(_PID)
    pdg = _yaml_pdg_block()
    probability = pdg["primary_current_probability_limit"]
    coupling = pdg["pdg_effective_coupling_limit"]
    original = pdg["original_experiment"]

    assert constraint.anchor.value == pytest.approx(probability["value"])
    assert constraint.anchor.budget == pytest.approx(probability["value"])
    assert constraint.anchor.budget == pytest.approx(8.3e-11)
    assert constraint.anchor.probability_limit.source_url == probability["source_url"]
    assert constraint.anchor.probability_metadata.relation == probability["relation"]
    assert constraint.anchor.probability_metadata.confidence_level == (
        probability["confidence_level"]
    )
    assert constraint.anchor.probability_metadata.conditions == probability["conditions"]
    assert constraint.anchor.coupling_limit_over_gf == pytest.approx(
        coupling["value"]
    )
    assert constraint.anchor.effective_coupling_limit.source_url == (
        coupling["source_url"]
    )
    assert constraint.anchor.coupling_metadata.relation == coupling["relation"]
    assert constraint.anchor.original_experiment.value == pytest.approx(
        original["value"]
    )
    assert constraint.anchor.original_experiment.source_url == original["source_url"]
    assert constraint.anchor.original_metadata.confidence_level == (
        original["confidence_level"]
    )

    with pytest.raises(anchors.AnchorError):
        anchors.load_anchor(
            _PID,
            family="charged_lepton",
            candidates=("no_such_block",),
        )


def test_evaluate_without_input_degrades_gracefully():
    result = fcc.get(_PID).evaluate(point_builder.empty_point())

    assert result.process_id == _PID
    assert result.passes is True
    assert result.predicted is None
    assert result.ratio is None
    assert result.sm_prediction == pytest.approx(0.0)
    assert result.experimental == pytest.approx(fcc.get(_PID).anchor.value)
    assert result.budget == pytest.approx(fcc.get(_PID).anchor.budget)
    assert result.notes.startswith("NOT EVALUATED --")
    assert result.diagnostics["evaluated"] is False
    assert result.diagnostics["missing_extra"] == "lepton_mass_basis_couplings"
    assert result.diagnostics["budget_limit_status"] == "observed_experimental_bound"
    assert "needs_human_physics" not in result.diagnostics


def test_invalid_lepton_input_is_unevaluated_not_real_pass():
    result = fcc.get(_PID).evaluate(
        point_builder.make_point(lepton_mass_basis_couplings={"source": "empty"})
    )

    assert result.passes is False
    assert result.predicted is None
    assert result.ratio is None
    assert result.notes.startswith("NOT EVALUATED --")
    assert result.diagnostics["evaluated"] is True
    assert result.diagnostics["invalid_extra"] == "lepton_mass_basis_couplings"
    assert result.diagnostics["exception_type"] == "KeyError"


def test_proxy_numerics_match_independent_anchor_calibration():
    constraint = fcc.get(_PID)
    g_proxy = 1.5e-4 + 2.0e-4j
    point = _point_from_mapping(
        {
            "g_mmbar_over_gf": g_proxy,
            "source": "L006 test effective-coupling proxy",
        }
    )
    result = constraint.evaluate(point)
    expected_probability = _manual_probability_from_coupling(
        g_proxy,
        probability_limit=constraint.anchor.budget,
        coupling_limit_over_gf=constraint.anchor.coupling_limit_over_gf,
    )
    expected_ratio = expected_probability / constraint.anchor.budget

    assert result.predicted == pytest.approx(expected_probability)
    assert result.predicted == pytest.approx(5.763888888888889e-13)
    assert result.ratio == pytest.approx(expected_ratio)
    assert result.sm_prediction == pytest.approx(0.0)
    assert result.budget == pytest.approx(constraint.anchor.budget)
    assert result.experimental == pytest.approx(constraint.anchor.value)
    assert result.diagnostics["effective_coupling_over_gf"] == pytest.approx(
        g_proxy
    )
    assert result.diagnostics["effective_coupling_abs_over_gf"] == pytest.approx(
        abs(g_proxy)
    )
    assert result.diagnostics["probability_per_coupling_ratio_squared"] == (
        pytest.approx(
            constraint.anchor.budget / constraint.anchor.coupling_limit_over_gf**2
        )
    )
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]
    assert result.diagnostics["input_kind"] == "g_mmbar_over_gf"
    assert result.diagnostics["direct_probability_proxy_used"] is False


def test_direct_probability_proxy_is_compared_to_macs_probability_limit():
    constraint = fcc.get(_PID)
    result = constraint.evaluate(
        _point_from_mapping(
            {
                "conversion_probability": 4.15e-11,
                "source": "L006 test direct probability proxy",
            }
        )
    )

    assert result.predicted == pytest.approx(4.15e-11)
    assert result.ratio == pytest.approx(0.5)
    assert result.passes is True
    assert result.diagnostics["effective_coupling_over_gf"] is None
    assert result.diagnostics["direct_probability_proxy_used"] is True
    assert result.diagnostics["input_kind"] == "conversion_probability"


def test_evaluate_runs_end_to_end_with_real_finite_fields_and_complex_diagnostics():
    result = fcc.get(_PID).evaluate(
        _point_from_mapping(
            {
                "G_C_over_G_F": 2.0e-4 - 1.0e-4j,
                "source": "L006 finite field probe",
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
    assert isinstance(result.diagnostics["effective_coupling_over_gf"], complex)
    for key in (
        "probability_limit",
        "effective_coupling_limit_over_gf",
        "probability_per_coupling_ratio_squared",
        "effective_coupling_abs_over_gf",
    ):
        value = result.diagnostics[key]
        assert isinstance(value, float)
        assert math.isfinite(value)
    assert result.diagnostics["evaluated"] is True
    assert result.diagnostics["budget_limit_status"] == "observed_experimental_bound"
    assert result.diagnostics["budget_verdict_role"] == "HARD observed upper-limit veto"


@pytest.mark.parametrize(
    ("lepton", "expected_pass"),
    [
        ({"g_mmbar_over_gf": 1.0e-4, "source": "safe L006 proxy"}, True),
        ({"g_mmbar_over_gf": 1.0e-2, "source": "excluded L006 proxy"}, False),
    ],
)
def test_safe_point_passes_and_large_np_point_fails(lepton, expected_pass: bool):
    result = fcc.get(_PID).evaluate(_point_from_mapping(lepton))

    assert result.passes is expected_pass
    if expected_pass:
        assert result.ratio < 1.0
    else:
        assert result.ratio > 1.0


def test_evaluate_is_pure_and_deterministic():
    lepton = {
        "g_mmbar_over_gf": 2.0e-4,
        "source": "L006 purity proxy",
    }
    point = _point_from_mapping(lepton)
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    assert point.get_extra("lepton_mass_basis_couplings") == lepton
