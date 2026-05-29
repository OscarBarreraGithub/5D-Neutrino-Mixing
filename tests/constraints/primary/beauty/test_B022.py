"""Production tests for B022 (B+ -> K+ nu nubar)."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest
import yaml

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import point_builder
from flavor_catalog_constraints.anchors import Anchor, AnchorError
from flavor_catalog_constraints.base import ConstraintProtocol, Severity
from flavor_catalog_constraints.primary.beauty import B022 as b022_module
from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.rare_b_nunu import (
    default_sm_inputs,
    evaluate_bplus_to_kplus_nunu,
    sm_inputs_with_bplus_kplus_normalization,
)

_PID = "B022"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "beauty" / "B022.yaml"


def _yaml_pdg_block():
    with open(_SIDECAR) as handle:
        return yaml.safe_load(handle)["pdg_or_equivalent"]


def _find_entry(entries, key: str, expected: str):
    for entry in entries:
        if str(entry.get(key)) == expected:
            return entry
    raise AssertionError(f"missing {key}={expected!r}")


def _sb_couplings(
    left: complex,
    right: complex = 0.0j,
    M_KK: float = 3000.0,
) -> QuarkMassBasisCouplings:
    """Minimal valid mass-basis couplings with only the s-b slot populated."""
    zeros = np.zeros((3, 3), dtype=np.complex128)
    left_down = zeros.copy()
    right_down = zeros.copy()
    left_down[1, 2] = left
    left_down[2, 1] = np.conj(left)
    right_down[1, 2] = right
    right_down[2, 1] = np.conj(right)
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


def _core_inputs_for_constraint():
    constraint = fcc.get(_PID)
    return sm_inputs_with_bplus_kplus_normalization(
        constraint.anchor.sm_value,
        inputs=default_sm_inputs(),
    )


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.HARD
    assert constraint.family == "beauty"
    assert constraint.observable == "BR(B+ -> K+ nu nubar)"


def test_anchor_matches_yaml_and_budget_band():
    constraint = fcc.get(_PID)
    pdg = _yaml_pdg_block()
    exp = _find_entry(
        pdg["observables"],
        "name",
        "Belle II evidence measurement",
    )
    sm = _find_entry(
        pdg["values"],
        "value_id",
        "HPQCD2023:B022:sm_prediction",
    )
    exp_up = math.hypot(
        float(exp["statistical_uncertainty"]),
        float(exp["systematic_uncertainty_positive"]),
    )
    exp_down = math.hypot(
        float(exp["statistical_uncertainty"]),
        float(exp["systematic_uncertainty_negative"]),
    )
    sm_sigma = float(sm["uncertainty"])
    combined_up = math.sqrt(exp_up * exp_up + sm_sigma * sm_sigma)
    combined_down = math.sqrt(exp_down * exp_down + sm_sigma * sm_sigma)

    assert constraint.anchor.experimental.value == pytest.approx(exp["value"])
    assert constraint.anchor.experimental.uncertainty_upper == pytest.approx(exp_up)
    assert constraint.anchor.experimental.uncertainty_lower == pytest.approx(exp_down)
    assert constraint.anchor.experimental.source_url == exp["source_url"]
    assert constraint.anchor.experimental.snapshot_path == exp["snapshot_path"]
    assert constraint.anchor.standard_model.value == pytest.approx(float(sm["value"]))
    assert constraint.anchor.standard_model.uncertainty == pytest.approx(sm_sigma)
    assert constraint.anchor.standard_model.source_url == sm["source_url"]
    assert constraint.anchor.standard_model.snapshot_path == sm["snapshot_path"]

    assert constraint.anchor.budget_band.central_residual == pytest.approx(
        float(exp["value"]) - float(sm["value"])
    )
    assert constraint.anchor.budget_band.combined_sigma_upper == pytest.approx(
        combined_up
    )
    assert constraint.anchor.budget_band.combined_sigma_lower == pytest.approx(
        combined_down
    )
    assert constraint.anchor.budget == pytest.approx(max(combined_up, combined_down))
    assert constraint.anchor.budget == pytest.approx(7.080741486595878e-6)


def test_value_anchors_route_through_scaffold_load_anchor_and_fail_loudly(
    monkeypatch,
):
    calls: list[tuple[str, ...]] = []
    original_load_anchor = b022_module.anchor_scaffold.load_anchor

    def spy_load_anchor(*args, **kwargs):
        calls.append(tuple(kwargs["candidates"]))
        return original_load_anchor(*args, **kwargs)

    monkeypatch.setattr(b022_module.anchor_scaffold, "load_anchor", spy_load_anchor)
    anchor = b022_module._load_b022_anchor(_PID)

    assert calls == [
        ("pdg_or_equivalent.observables[2]",),
        ("pdg_or_equivalent.values[0]",),
    ]
    assert anchor.value == pytest.approx(fcc.get(_PID).anchor.value)
    assert anchor.sm_value == pytest.approx(5.58e-6)

    def missing_load_anchor(*args, **kwargs):
        raise AnchorError("forced missing load_anchor")

    monkeypatch.setattr(
        b022_module.anchor_scaffold,
        "load_anchor",
        missing_load_anchor,
    )
    with pytest.raises(AnchorError, match="forced missing load_anchor"):
        b022_module._load_b022_anchor(_PID)

    def mismatched_load_anchor(*args, **kwargs):
        return Anchor(
            process_id=_PID,
            block_key="pdg_or_equivalent.values[99]",
            value=1.0,
            uncertainty=0.1,
        )

    monkeypatch.setattr(
        b022_module.anchor_scaffold,
        "load_anchor",
        mismatched_load_anchor,
    )
    with pytest.raises(AnchorError, match="load_anchor selected"):
        b022_module._load_b022_anchor(_PID)


def test_evaluate_without_input_degrades_gracefully():
    constraint = fcc.get(_PID)
    result = constraint.evaluate(point_builder.empty_point())

    assert result.process_id == _PID
    assert result.passes is True
    assert result.predicted is None
    assert result.ratio is None
    assert result.experimental == pytest.approx(constraint.anchor.value)
    assert result.sm_prediction == pytest.approx(
        constraint.sm_result.branching_fraction
    )
    assert result.budget == pytest.approx(constraint.anchor.budget)
    assert result.diagnostics["missing_extra"] == "quark_mass_basis_couplings"
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]


def test_sm_limit_branching_fraction_and_top_function_validation():
    constraint = fcc.get(_PID)
    couplings = _sb_couplings(left=0.0j, right=0.0j)
    point = point_builder.build_from_quark_couplings(couplings)
    result = constraint.evaluate(point)
    direct = evaluate_bplus_to_kplus_nunu(couplings, inputs=_core_inputs_for_constraint())

    x = (163.5 / 80.379) ** 2
    x0 = x / 8.0 * (
        (x + 2.0) / (x - 1.0)
        + (3.0 * x - 6.0) / ((x - 1.0) ** 2) * math.log(x)
    )
    expected_xt = 0.994 * x0

    assert result.predicted == pytest.approx(5.58e-6)
    assert result.predicted == pytest.approx(result.sm_prediction)
    assert result.predicted == pytest.approx(direct.branching_fraction)
    assert result.diagnostics["x_t"] == pytest.approx(expected_xt)
    assert result.diagnostics["c_l_sm"].real == pytest.approx(
        -expected_xt / 0.23122
    )
    assert result.ratio == pytest.approx(
        abs(result.predicted - constraint.anchor.value)
        / constraint.anchor.budget_band.combined_sigma_lower
    )
    assert result.passes is False
    assert result.ratio > 2.7


def test_evaluate_runs_end_to_end_with_real_finite_fields_and_complex_diagnostics():
    couplings = _sb_couplings(left=-0.7 + 0.05j, right=0.02j)
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
        "left_sb_coupling",
        "right_sb_coupling",
        "lambda_t_bs",
        "c_l_sm",
        "c_l_total",
        "c_r_total",
        "x_eff_left",
        "x_eff_right",
        "x_np_total",
    ):
        assert isinstance(result.diagnostics[key], complex)
    for key in (
        "m_kk_gev",
        "matching_scale_gev",
        "x_t",
        "epsilon",
        "eta",
        "r_k",
        "r_kstar",
        "np_shift_branching_fraction",
        "budget_combined_sigma_upper",
        "budget_combined_sigma_lower",
    ):
        assert isinstance(result.diagnostics[key], float)
        assert math.isfinite(result.diagnostics[key])
    assert result.diagnostics["down_sector_indices"] == (1, 2)
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]


@pytest.mark.parametrize(
    ("couplings", "expected_pass"),
    [
        (_sb_couplings(left=-0.7), True),
        (_sb_couplings(left=-2.0), False),
    ],
)
def test_belle_like_point_passes_and_large_np_point_fails(
    couplings: QuarkMassBasisCouplings,
    expected_pass: bool,
):
    constraint = fcc.get(_PID)
    point = point_builder.build_from_quark_couplings(couplings)
    result = constraint.evaluate(point)
    direct = evaluate_bplus_to_kplus_nunu(couplings, inputs=_core_inputs_for_constraint())

    assert result.passes is expected_pass
    assert result.predicted == pytest.approx(direct.branching_fraction)
    assert result.sm_prediction == pytest.approx(direct.sm_branching_fraction)
    if expected_pass:
        assert result.ratio < 1.0
    else:
        assert result.ratio > 1.0


def test_charged_long_distance_piece_is_not_rescaled_by_short_distance_np():
    couplings = _sb_couplings(left=-0.7)
    direct = evaluate_bplus_to_kplus_nunu(couplings, inputs=_core_inputs_for_constraint())
    br_ld = 6.09e-7
    br_sm = 5.58e-6

    assert direct.r_k == pytest.approx(4.00211528673656)
    assert direct.branching_fraction == pytest.approx(
        br_ld + (br_sm - br_ld) * direct.r_k
    )
    assert direct.branching_fraction == pytest.approx(2.0503514974408267e-5)
    assert direct.branching_fraction != pytest.approx(br_sm * direct.r_k)


def test_optional_kk_ew_mass_extra_changes_matching_scale():
    couplings = _sb_couplings(left=-0.7)
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
    couplings = _sb_couplings(left=-0.7 + 0.05j, right=0.02j)
    before_left_down = couplings.left_down.copy()
    before_right_down = couplings.right_down.copy()
    point = point_builder.build_from_quark_couplings(couplings)
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    np.testing.assert_array_equal(couplings.left_down, before_left_down)
    np.testing.assert_array_equal(couplings.right_down, before_right_down)
