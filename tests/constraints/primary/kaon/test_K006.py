"""Production tests for K006 (K_L -> mu+ mu- short-distance)."""

from __future__ import annotations

import math
from pathlib import Path
import re

import numpy as np
import pytest
import yaml

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import point_builder
from flavor_catalog_constraints.base import ConstraintProtocol, Severity
from flavor_catalog_constraints.primary.kaon import K006 as k006_module
from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.model import RotationParameters, ckm_like_unitary
from quarkConstraints.rare_kaon_dilepton import evaluate_klong_mumu_short_distance

_PID = "K006"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "kaon" / "K006.yaml"
_LIMIT_RE = re.compile(r"^\s*<\s*(?P<value>[0-9.eE+-]+)\s*$")


def _yaml():
    with open(_SIDECAR) as handle:
        return yaml.safe_load(handle)


def _supporting_entry(observable: str):
    for entry in _yaml()["supporting_numeric_values"]:
        if entry["observable"] == observable:
            return entry
    raise AssertionError(f"no supporting entry for {observable!r}")


def _limit_value(entry) -> float:
    match = _LIMIT_RE.match(entry["value"])
    assert match is not None
    return float(match.group("value"))


def _sd_couplings(
    left: complex,
    right: complex = 0.0j,
    M_KK: float = 3000.0,
) -> QuarkMassBasisCouplings:
    """Minimal valid mass-basis couplings with only the s-d slot populated."""
    zeros = np.zeros((3, 3), dtype=np.complex128)
    left_down = zeros.copy()
    right_down = zeros.copy()
    left_down[0, 1] = left
    left_down[1, 0] = np.conj(left)
    right_down[0, 1] = right
    right_down[1, 0] = np.conj(right)
    return QuarkMassBasisCouplings(
        M_KK=M_KK,
        xi_KK=1.0,
        alpha_s=0.09,
        g_s=1.0,
        left_overlap=zeros,
        right_up_overlap=zeros,
        right_down_overlap=zeros,
        left_up=zeros,
        left_down=left_down,
        right_up=zeros,
        right_down=right_down,
    )


def _independent_sm_short_distance(inputs) -> float:
    matrix = ckm_like_unitary(
        RotationParameters(
            theta12=inputs.theta12,
            theta13=inputs.theta13,
            theta23=inputs.theta23,
            delta=inputs.delta,
        )
    )
    lam = float(abs(matrix[0, 1]))
    lambda_c = complex(np.conjugate(matrix[1, 1]) * matrix[1, 0])
    lambda_t = complex(np.conjugate(matrix[2, 1]) * matrix[2, 0])
    kappa_mu = inputs.kappa_mu_ref * (lam / inputs.kappa_lambda_ref) ** 8
    amplitude = (
        (lambda_t * inputs.y_t).real / lam**5
        + lambda_c.real / lam * inputs.p_c_y
    )
    return float(kappa_mu * amplitude * amplitude)


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.HARD
    assert constraint.family == "kaon"
    assert constraint.observable == "BR(K_L -> mu+ mu-)_SD"


def test_anchor_matches_yaml_and_budget_band():
    constraint = fcc.get(_PID)
    pdg = _yaml()["pdg_or_equivalent"]
    sm_sd = pdg["standard_model_short_distance_prediction"]
    e871 = _supporting_entry("BNL E871 BR(K_L -> mu+ mu-) measurement")
    sd_bound = _supporting_entry(
        "Conservative upper bound on B(K_L -> mu+ mu-)_short"
    )
    limit = _limit_value(sd_bound)

    assert constraint.anchor.total_rate.value == pytest.approx(float(pdg["value"]))
    assert constraint.anchor.total_rate.uncertainty == pytest.approx(
        float(pdg["uncertainty"])
    )
    assert constraint.anchor.total_rate.source_url == pdg["source_url"]
    assert constraint.anchor.e871_measurement.value == pytest.approx(
        float(e871["value"])
    )
    assert constraint.anchor.e871_measurement.uncertainty == pytest.approx(
        float(e871["uncertainty"])
    )
    assert constraint.anchor.e871_measurement.source_url == e871["source_url"]
    assert constraint.anchor.standard_model_short_distance.value == pytest.approx(
        float(sm_sd["value"])
    )
    assert constraint.anchor.standard_model_short_distance.uncertainty == (
        pytest.approx(float(sm_sd["uncertainty"]))
    )
    assert constraint.anchor.standard_model_short_distance.source_url == (
        sm_sd["source_url"]
    )
    assert abs(
        constraint.sm_result.branching_fraction
        - constraint.anchor.standard_model_short_distance.value
    ) < constraint.anchor.standard_model_short_distance.uncertainty
    assert constraint.anchor.short_distance_bound.value == pytest.approx(limit)
    assert constraint.anchor.short_distance_bound.is_upper_limit is True
    assert constraint.anchor.short_distance_bound.source_url == sd_bound["source_url"]
    assert constraint.anchor.budget == pytest.approx(limit)
    assert constraint.anchor.budget_band.sm_formula_value == pytest.approx(
        constraint.sm_result.branching_fraction
    )
    assert constraint.anchor.budget_band.sm_anchor_value == pytest.approx(
        float(sm_sd["value"])
    )
    assert constraint.anchor.budget_band.sm_anchor_uncertainty == pytest.approx(
        float(sm_sd["uncertainty"])
    )
    assert constraint.anchor.budget_band.sm_formula_minus_anchor == pytest.approx(
        constraint.sm_result.branching_fraction - float(sm_sd["value"])
    )
    assert constraint.anchor.budget_band.bound_minus_sm_formula == pytest.approx(
        limit - constraint.sm_result.branching_fraction
    )
    assert constraint.anchor.budget_band.long_distance_dominated is True
    assert constraint.anchor.budget_band.sm_subtracted is False


def test_k006_anchors_route_through_scaffold_load_anchor(monkeypatch):
    calls: list[tuple[str, ...]] = []
    original_load_anchor = k006_module.load_anchor

    def spy_load_anchor(*args, **kwargs):
        calls.append(tuple(kwargs["candidates"]))
        return original_load_anchor(*args, **kwargs)

    monkeypatch.setattr(k006_module, "load_anchor", spy_load_anchor)
    anchor = k006_module._load_k006_anchor(
        _PID,
        sm_formula_value=fcc.get(_PID).sm_result.branching_fraction,
    )

    assert calls == [
        ("pdg_or_equivalent",),
        ("supporting_numeric_values[0]",),
        ("standard_model_short_distance_prediction",),
        ("supporting_numeric_values[2]",),
    ]
    assert anchor.short_distance_bound.value == pytest.approx(
        fcc.get(_PID).anchor.value
    )


def test_k006_anchor_loud_fail_probe():
    with pytest.raises(k006_module.AnchorError):
        k006_module._load_supporting_anchor(
            _PID,
            block_key="supporting_numeric_values[2]",
            expect_upper_limit=True,
            observable="not the K006 short-distance bound",
        )


def test_evaluate_without_input_degrades_gracefully():
    result = fcc.get(_PID).evaluate(point_builder.empty_point())

    assert result.process_id == _PID
    assert result.passes is True
    assert result.predicted is None
    assert result.ratio is None
    assert result.experimental == pytest.approx(fcc.get(_PID).anchor.value)
    assert result.sm_prediction == pytest.approx(
        fcc.get(_PID).sm_result.branching_fraction
    )
    assert result.budget == pytest.approx(fcc.get(_PID).anchor.budget)
    assert result.diagnostics["missing_extra"] == "quark_mass_basis_couplings"
    assert result.diagnostics["long_distance_dominated_total_rate"] is True


def test_sm_limit_branching_fraction_matches_short_distance_reference():
    constraint = fcc.get(_PID)
    couplings = _sd_couplings(left=0.0j, right=0.0j)
    point = point_builder.build_from_quark_couplings(couplings)
    result = constraint.evaluate(point)
    direct = evaluate_klong_mumu_short_distance(couplings)
    independent = _independent_sm_short_distance(constraint.sm_inputs)
    sm_anchor = constraint.anchor.standard_model_short_distance

    assert result.predicted == pytest.approx(constraint.sm_result.branching_fraction)
    assert result.predicted == pytest.approx(independent)
    assert result.predicted == pytest.approx(result.sm_prediction)
    assert result.predicted == pytest.approx(direct.sm_branching_fraction)
    assert abs(result.predicted - sm_anchor.value) < sm_anchor.uncertainty
    assert result.diagnostics["sm_anchor_branching_fraction"] == pytest.approx(
        sm_anchor.value
    )
    assert result.diagnostics["sm_anchor_uncertainty"] == pytest.approx(
        sm_anchor.uncertainty
    )
    assert result.predicted < constraint.anchor.short_distance_bound.value
    assert result.diagnostics["sm_short_distance_amplitude"] == pytest.approx(
        direct.diagnostics["sm_short_distance_amplitude"]
    )
    assert result.diagnostics["kappa_mu"] == pytest.approx(direct.kappa_mu)
    assert result.diagnostics["p_c_y"] == pytest.approx(direct.p_c_y)
    assert result.diagnostics["y_t"] == pytest.approx(direct.y_t)
    assert result.passes is True


def test_evaluate_runs_end_to_end_with_real_finite_fields_and_complex_diagnostics():
    couplings = _sd_couplings(left=1.0e-5 + 0.2e-5j, right=0.5e-5)
    point = point_builder.build_from_quark_couplings(couplings)
    result = fcc.get(_PID).evaluate(point)

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
        "left_sd_coupling",
        "right_sd_coupling",
        "lambda_c",
        "lambda_t",
        "y_eff",
        "y_np_total",
    ):
        assert isinstance(result.diagnostics[key], complex)
    for key in (
        "m_kk_gev",
        "matching_scale_gev",
        "kappa_mu",
        "p_c_y",
        "y_t",
        "lambda_wolfenstein",
        "short_distance_amplitude",
        "np_shift_branching_fraction",
        "budget_bound_minus_sm_formula",
    ):
        assert isinstance(result.diagnostics[key], float)
        assert math.isfinite(result.diagnostics[key])
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]
    assert result.diagnostics["short_distance_only"] is True
    assert result.diagnostics["right_handed_enters_with_minus_sign"] is True


@pytest.mark.parametrize(
    ("couplings", "expected_pass"),
    [
        (_sd_couplings(left=1.0e-5), True),
        (_sd_couplings(left=5.0e-2), False),
    ],
)
def test_safe_point_passes_and_large_np_point_fails(
    couplings: QuarkMassBasisCouplings,
    expected_pass: bool,
):
    constraint = fcc.get(_PID)
    point = point_builder.build_from_quark_couplings(couplings)
    result = constraint.evaluate(point)
    direct = evaluate_klong_mumu_short_distance(couplings)

    assert result.passes is expected_pass
    assert result.predicted == pytest.approx(direct.branching_fraction)
    assert result.sm_prediction == pytest.approx(direct.sm_branching_fraction)
    assert result.ratio == pytest.approx(
        direct.branching_fraction / constraint.anchor.budget
    )
    if expected_pass:
        assert result.ratio < 1.0
    else:
        assert result.ratio > 10.0


def test_optional_kk_ew_mass_extra_changes_matching_scale():
    couplings = _sd_couplings(left=1.0e-3)
    default_point = point_builder.build_from_quark_couplings(couplings)
    ew_point = point_builder.make_point(
        quark_mass_basis_couplings=couplings,
        kk_ew_mass_gev=6000.0,
    )
    default_result = fcc.get(_PID).evaluate(default_point)
    ew_result = fcc.get(_PID).evaluate(ew_point)

    assert default_result.diagnostics["m_kk_gev"] == pytest.approx(3000.0)
    assert ew_result.diagnostics["m_kk_gev"] == pytest.approx(6000.0)
    assert ew_result.diagnostics["kk_ew_mass_extra_used"] is True
    assert abs(ew_result.diagnostics["y_np_total"]) == pytest.approx(
        abs(default_result.diagnostics["y_np_total"]) / 4.0
    )


def test_evaluate_is_pure_and_deterministic():
    couplings = _sd_couplings(left=1.0e-5 + 0.2e-5j, right=0.5e-5)
    before_left_down = couplings.left_down.copy()
    before_right_down = couplings.right_down.copy()
    point = point_builder.build_from_quark_couplings(couplings)
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    np.testing.assert_array_equal(couplings.left_down, before_left_down)
    np.testing.assert_array_equal(couplings.right_down, before_right_down)
