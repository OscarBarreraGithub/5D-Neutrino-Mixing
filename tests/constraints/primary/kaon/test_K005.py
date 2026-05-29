"""Production tests for K005 (K_L -> pi0 nu nubar)."""

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
from flavor_catalog_constraints.primary.kaon import K005 as k005_module
from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.rare_kaon_snd import evaluate_klong_to_pi0_nunu

_PID = "K005"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "kaon" / "K005.yaml"
_LIMIT_RE = re.compile(r"^\s*<\s*(?P<value>[0-9.eE+-]+)\s*$")


def _yaml():
    with open(_SIDECAR) as handle:
        return yaml.safe_load(handle)


def _source_entry(section: str, source_needle: str):
    for entry in _yaml()[section]:
        if source_needle in entry["source"]:
            return entry
    raise AssertionError(f"no {section} entry with source containing {source_needle}")


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


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.HARD
    assert constraint.family == "kaon"
    assert constraint.observable == "BR(K_L -> pi0 nu nubar)"


def test_anchor_matches_yaml_and_budget_band():
    constraint = fcc.get(_PID)
    koto = _source_entry("experimental_inputs", "KOTO")
    pdg = _source_entry("pdg_or_equivalent", "PDG")
    bv = _source_entry("theory_inputs", "Buras and Venturini")
    bgs = _source_entry("theory_inputs", "Brod, Gorbahn, Stamou")
    koto_limit = _limit_value(koto)

    assert constraint.anchor.experimental.value == pytest.approx(koto_limit)
    assert constraint.anchor.experimental.source_url == koto["source_url"]
    assert constraint.anchor.experimental.confidence_level == koto["confidence_level"]
    assert constraint.anchor.pdg_limit.value == pytest.approx(_limit_value(pdg))
    assert constraint.anchor.pdg_limit.source_url == pdg["source_url"]
    assert constraint.anchor.standard_model.value == pytest.approx(float(bv["value"]))
    assert constraint.anchor.standard_model.uncertainty == pytest.approx(
        float(bv["uncertainty"])
    )
    assert constraint.anchor.standard_model.source_url == bv["source_url"]
    assert constraint.anchor.validation_standard_model.value == pytest.approx(
        float(bgs["value"])
    )
    assert constraint.anchor.validation_standard_model.uncertainty == pytest.approx(
        float(bgs["uncertainty"])
    )
    assert constraint.anchor.budget == pytest.approx(koto_limit)
    assert constraint.anchor.budget_band.limit_minus_sm_anchor == pytest.approx(
        koto_limit - float(bv["value"])
    )
    assert constraint.anchor.budget_band.sm_subtracted is False
    assert constraint.anchor.budget_band.confidence_level == "90%"


def test_k005_list_anchors_route_through_scaffold_load_anchor(monkeypatch):
    calls: list[tuple[str, ...]] = []
    original_load_anchor = k005_module.load_anchor

    def spy_load_anchor(*args, **kwargs):
        calls.append(tuple(kwargs["candidates"]))
        return original_load_anchor(*args, **kwargs)

    monkeypatch.setattr(k005_module, "load_anchor", spy_load_anchor)
    anchor = k005_module._load_k005_anchor(_PID)

    assert calls == [
        ("experimental_inputs[0]",),
        ("pdg_or_equivalent[0]",),
        ("theory_inputs[1]",),
        ("theory_inputs[0]",),
    ]
    assert anchor.experimental.value == pytest.approx(fcc.get(_PID).anchor.value)


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


def test_sm_limit_branching_fraction_matches_short_distance_reference():
    constraint = fcc.get(_PID)
    couplings = _sd_couplings(left=0.0j, right=0.0j)
    point = point_builder.build_from_quark_couplings(couplings)
    result = constraint.evaluate(point)
    direct = evaluate_klong_to_pi0_nunu(couplings)
    sm_anchor = constraint.anchor.sm_value

    assert result.predicted == pytest.approx(2.95375989343059e-11)
    assert result.predicted == pytest.approx(result.sm_prediction)
    assert result.predicted == pytest.approx(direct.sm_branching_fraction)
    assert result.predicted == pytest.approx(direct.branching_fraction)
    assert abs(result.predicted - sm_anchor) < constraint.anchor.standard_model.uncertainty
    assert result.diagnostics["sm_anchor_branching_fraction"] == pytest.approx(
        sm_anchor
    )
    assert result.diagnostics["sm_formula_minus_anchor"] == pytest.approx(
        result.sm_prediction - sm_anchor
    )
    assert result.diagnostics["kappa_l"] == pytest.approx(2.266439872701711e-10)
    assert result.diagnostics["imaginary_projection_only"] is True
    assert result.passes is True


def test_evaluate_runs_end_to_end_with_real_finite_fields_and_complex_diagnostics():
    couplings = _sd_couplings(left=1.0e-5 + 0.2e-5j, right=0.5e-5j)
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
        "lambda_t",
        "x_eff_top",
        "x_np_total",
    ):
        assert isinstance(result.diagnostics[key], complex)
    for key in (
        "m_kk_gev",
        "matching_scale_gev",
        "kappa_l",
        "x_t",
        "lambda_wolfenstein",
        "im_x_eff_top",
        "np_shift_branching_fraction",
        "budget_limit_minus_sm_anchor",
    ):
        assert isinstance(result.diagnostics[key], float)
        assert math.isfinite(result.diagnostics[key])
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]
    assert result.diagnostics["qcd_running_applied"] is False
    assert result.diagnostics["budget_sm_subtracted"] is False


@pytest.mark.parametrize(
    ("couplings", "expected_pass"),
    [
        (_sd_couplings(left=1.0e-5j), True),
        (_sd_couplings(left=1.0e-1j), False),
    ],
)
def test_safe_point_passes_and_large_np_point_fails(
    couplings: QuarkMassBasisCouplings,
    expected_pass: bool,
):
    constraint = fcc.get(_PID)
    point = point_builder.build_from_quark_couplings(couplings)
    result = constraint.evaluate(point)
    direct = evaluate_klong_to_pi0_nunu(couplings)

    assert result.passes is expected_pass
    assert result.predicted == pytest.approx(direct.branching_fraction)
    assert result.sm_prediction == pytest.approx(direct.sm_branching_fraction)
    assert result.ratio == pytest.approx(direct.branching_fraction / constraint.anchor.budget)
    if expected_pass:
        assert result.ratio < 1.0
    else:
        assert result.ratio > 10.0


def test_optional_kk_ew_mass_extra_changes_matching_scale():
    couplings = _sd_couplings(left=1.0e-3j)
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
    assert abs(ew_result.diagnostics["x_np_total"]) == pytest.approx(
        abs(default_result.diagnostics["x_np_total"]) / 4.0
    )


def test_evaluate_is_pure_and_deterministic():
    couplings = _sd_couplings(left=1.0e-5 + 0.2e-5j, right=0.5e-5j)
    before_left_down = couplings.left_down.copy()
    before_right_down = couplings.right_down.copy()
    point = point_builder.build_from_quark_couplings(couplings)
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    np.testing.assert_array_equal(couplings.left_down, before_left_down)
    np.testing.assert_array_equal(couplings.right_down, before_right_down)
