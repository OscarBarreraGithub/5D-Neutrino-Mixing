"""Production tests for B023 (B -> K* nu nubar)."""

from __future__ import annotations

import math
from pathlib import Path
import re

import numpy as np
import pytest
import yaml

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import point_builder
from flavor_catalog_constraints.anchors import Anchor, AnchorError
from flavor_catalog_constraints.base import ConstraintProtocol, Severity
from flavor_catalog_constraints.primary.beauty import B023 as b023_module
from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.rare_b_nunu import (
    RARE_B_NUNU_KSTAR_ETA_COEFFICIENT,
    compute_rare_b_nunu_wilsons,
    short_distance_response,
)

_PID = "B023"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "beauty" / "B023.yaml"
_LIMIT_RE = re.compile(r"^\s*<\s*(?P<value>[0-9.eE+-]+)\s*$")


def _yaml_pdg_block():
    with open(_SIDECAR) as handle:
        return yaml.safe_load(handle)["pdg_or_equivalent"]


def _value_entry(value_id: str):
    for entry in _yaml_pdg_block()["values"]:
        if str(entry.get("value_id")) == value_id:
            return entry
    raise AssertionError(f"missing value_id={value_id!r}")


def _limit_value(entry) -> float:
    match = _LIMIT_RE.match(entry["value"])
    assert match is not None
    return float(match.group("value"))


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


def _direct_kstar_branching(
    couplings: QuarkMassBasisCouplings,
    *,
    sm_branching_fraction: float,
    m_kk_gev: float | None = None,
) -> tuple[float, object, object]:
    wilsons = compute_rare_b_nunu_wilsons(couplings, m_kk_gev=m_kk_gev)
    response = short_distance_response(wilsons.x_np_left, wilsons.x_np_right)
    return float(sm_branching_fraction * response.r_kstar), response, wilsons


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.HARD
    assert constraint.family == "beauty"
    assert constraint.observable == "BR(B -> K* nu nubar)"


def test_anchor_matches_yaml_and_budget_band():
    constraint = fcc.get(_PID)
    exp = _value_entry("Belle2017:B023:combined_vector_limit")
    sm = _value_entry("Buras2015:B023:sm_b0_kstar0_nunu")
    limit = _limit_value(exp)
    sm_value = float(sm["value"])

    assert constraint.anchor.experimental.value == pytest.approx(limit)
    assert constraint.anchor.experimental.value_raw == exp["value"]
    assert constraint.anchor.experimental.is_upper_limit is True
    assert constraint.anchor.experimental.confidence_level == exp["confidence_level"]
    assert constraint.anchor.experimental.source_url == exp["source_url"]
    assert constraint.anchor.experimental.snapshot_path == exp["snapshot_path"]
    assert constraint.anchor.standard_model.value == pytest.approx(sm_value)
    assert constraint.anchor.standard_model.uncertainty == pytest.approx(
        float(sm["uncertainty"])
    )
    assert constraint.anchor.standard_model.source_url == sm["source_url"]
    assert constraint.anchor.standard_model.snapshot_path == sm["snapshot_path"]
    assert constraint.anchor.budget == pytest.approx(limit)
    assert constraint.anchor.budget == pytest.approx(2.7e-5)
    assert constraint.anchor.budget_band.limit_minus_sm_anchor == pytest.approx(
        limit - sm_value
    )
    assert constraint.anchor.budget_band.limit_minus_sm_anchor == pytest.approx(1.78e-5)
    assert constraint.anchor.budget_band.sm_subtracted is False
    assert constraint.anchor.budget_band.confidence_level == "90% CL"


def test_value_anchors_route_through_scaffold_load_anchor_and_fail_loudly(
    monkeypatch,
):
    calls: list[tuple[str, ...]] = []
    original_load_anchor = b023_module.anchor_scaffold.load_anchor

    def spy_load_anchor(*args, **kwargs):
        calls.append(tuple(kwargs["candidates"]))
        return original_load_anchor(*args, **kwargs)

    monkeypatch.setattr(b023_module.anchor_scaffold, "load_anchor", spy_load_anchor)
    anchor = b023_module._load_b023_anchor(_PID)

    assert calls == [
        ("pdg_or_equivalent.values[0]",),
        ("pdg_or_equivalent.values[1]",),
    ]
    assert anchor.value == pytest.approx(fcc.get(_PID).anchor.value)
    assert anchor.sm_value == pytest.approx(9.2e-6)

    def missing_load_anchor(*args, **kwargs):
        raise AnchorError("forced missing load_anchor")

    monkeypatch.setattr(
        b023_module.anchor_scaffold,
        "load_anchor",
        missing_load_anchor,
    )
    with pytest.raises(AnchorError, match="forced missing load_anchor"):
        b023_module._load_b023_anchor(_PID)

    def mismatched_load_anchor(*args, **kwargs):
        return Anchor(
            process_id=_PID,
            block_key="pdg_or_equivalent.values[99]",
            value=1.0,
            uncertainty=None,
        )

    monkeypatch.setattr(
        b023_module.anchor_scaffold,
        "load_anchor",
        mismatched_load_anchor,
    )
    with pytest.raises(AnchorError, match="load_anchor selected"):
        b023_module._load_b023_anchor(_PID)


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


def test_sm_limit_branching_fraction_matches_independent_core_response():
    constraint = fcc.get(_PID)
    couplings = _sb_couplings(left=0.0j, right=0.0j)
    point = point_builder.build_from_quark_couplings(couplings)
    result = constraint.evaluate(point)
    direct, response, _ = _direct_kstar_branching(
        couplings,
        sm_branching_fraction=constraint.anchor.sm_value,
    )

    assert result.predicted == pytest.approx(9.2e-6)
    assert result.predicted == pytest.approx(result.sm_prediction)
    assert result.predicted == pytest.approx(direct)
    assert result.predicted == pytest.approx(9.0e-6, rel=0.03)
    assert result.diagnostics["epsilon"] == pytest.approx(1.0)
    assert result.diagnostics["eta"] == pytest.approx(0.0)
    assert result.diagnostics["r_kstar"] == pytest.approx(1.0)
    assert response.r_kstar == pytest.approx(1.0)
    assert result.diagnostics["kstar_eta_coefficient"] == pytest.approx(
        RARE_B_NUNU_KSTAR_ETA_COEFFICIENT
    )
    assert result.ratio == pytest.approx(
        result.predicted / constraint.anchor.budget
    )
    assert result.passes is True


def test_evaluate_runs_end_to_end_with_real_finite_fields_and_complex_diagnostics():
    couplings = _sb_couplings(left=-0.4 + 0.05j, right=0.05 - 0.03j)
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
        "budget_limit_minus_sm_anchor",
    ):
        assert isinstance(result.diagnostics[key], float)
        assert math.isfinite(result.diagnostics[key])
    assert result.diagnostics["down_sector_indices"] == (1, 2)
    assert result.diagnostics["budget_sm_subtracted"] is False
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]


@pytest.mark.parametrize(
    ("couplings", "expected_pass"),
    [
        (_sb_couplings(left=-0.4), True),
        (_sb_couplings(left=-0.7), False),
    ],
)
def test_safe_point_passes_and_large_np_point_fails(
    couplings: QuarkMassBasisCouplings,
    expected_pass: bool,
):
    constraint = fcc.get(_PID)
    point = point_builder.build_from_quark_couplings(couplings)
    result = constraint.evaluate(point)
    direct, response, _ = _direct_kstar_branching(
        couplings,
        sm_branching_fraction=constraint.anchor.sm_value,
    )

    assert result.passes is expected_pass
    assert result.predicted == pytest.approx(direct)
    assert result.sm_prediction == pytest.approx(constraint.anchor.sm_value)
    assert result.ratio == pytest.approx(direct / constraint.anchor.budget)
    assert result.diagnostics["r_kstar"] == pytest.approx(response.r_kstar)
    if expected_pass:
        assert result.ratio < 1.0
    else:
        assert result.ratio > 1.0


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
    couplings = _sb_couplings(left=-0.4 + 0.05j, right=0.05 - 0.03j)
    before_left_down = couplings.left_down.copy()
    before_right_down = couplings.right_down.copy()
    point = point_builder.build_from_quark_couplings(couplings)
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    np.testing.assert_array_equal(couplings.left_down, before_left_down)
    np.testing.assert_array_equal(couplings.right_down, before_right_down)
