"""Production tests for B008 (B_s/B0 -> tau+ tau-)."""

from __future__ import annotations

from dataclasses import replace
import math
from pathlib import Path

import numpy as np
import pytest
import yaml

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import point_builder
from flavor_catalog_constraints.anchors import AnchorError
from flavor_catalog_constraints.base import ConstraintLevel, ConstraintProtocol, Severity
from flavor_catalog_constraints.secondary.beauty import B008 as b008_module
from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.model import RotationParameters, ckm_like_unitary
from quarkConstraints.rare_b_dilepton import evaluate_bq_to_mumu
from tests.rare_b_phase3d_helpers import (
    core_wilsons_from_rs_coeff,
    rs_coeff,
    sample_rare_b_point,
    scaled_rare_b_point,
    sm_limit_rare_b_point,
)

_PID = "B008"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = (
    _REPO_ROOT
    / "flavor_catalog"
    / "processes"
    / "secondary"
    / "beauty"
    / "B008.yaml"
)


def _yaml_doc():
    with open(_SIDECAR) as handle:
        return yaml.safe_load(handle)


def _value_row(row_id: str):
    for row in _yaml_doc()["pdg_or_equivalent"]["values"]:
        if row["id"] == row_id:
            return row
    raise AssertionError(f"no B008 pdg_or_equivalent.values row with id={row_id!r}")


def _sm_row(source_key: str, observable: str):
    for row in _yaml_doc()["auxiliary_theory_inputs"]:
        if row.get("source_key") == source_key and row.get("observable") == observable:
            return row
    raise AssertionError(
        f"no B008 SM row with source_key={source_key!r}, observable={observable!r}"
    )


def _bq_couplings(
    *,
    bs_left: complex = 0.0j,
    bd_left: complex = 0.0j,
    bs_right: complex = 0.0j,
    bd_right: complex = 0.0j,
    M_KK: float = 3000.0,
) -> QuarkMassBasisCouplings:
    """Minimal valid mass-basis couplings with b-s and b-d slots populated."""
    zeros = np.zeros((3, 3), dtype=np.complex128)
    left_down = zeros.copy()
    right_down = zeros.copy()
    left_down[1, 2] = bs_left
    left_down[2, 1] = np.conj(bs_left)
    left_down[0, 2] = bd_left
    left_down[2, 0] = np.conj(bd_left)
    right_down[1, 2] = bs_right
    right_down[2, 1] = np.conj(bs_right)
    right_down[0, 2] = bd_right
    right_down[2, 0] = np.conj(bd_right)
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


def _manual_sm_tautau(inputs, transition: str) -> float:
    matrix = ckm_like_unitary(
        RotationParameters(
            theta12=inputs.theta12,
            theta13=inputs.theta13,
            theta23=inputs.theta23,
            delta=inputs.delta,
        )
    )
    meson = inputs.meson(transition)
    lambda_t = complex(
        matrix[2, 2] * np.conjugate(matrix[2, meson.light_down_index])
    )
    tau_gev_inv = meson.lifetime_ps * 1.0e-12 / inputs.hbar_gev_s
    beta = math.sqrt(
        1.0 - 4.0 * inputs.muon_mass_gev**2 / meson.meson_mass_gev**2
    )
    time_factor = (1.0 + meson.width_difference_y) / (
        1.0 - meson.width_difference_y**2
    )
    return float(
        tau_gev_inv
        * inputs.gf_gev_minus2**2
        * inputs.alpha_em_mz**2
        / (16.0 * math.pi**3)
        * meson.decay_constant_gev**2
        * meson.meson_mass_gev
        * inputs.muon_mass_gev**2
        * abs(lambda_t) ** 2
        * beta
        * abs(inputs.c10_sm) ** 2
        * time_factor
    )


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.HARD
    assert constraint.level is ConstraintLevel.SECONDARY
    assert constraint.family == "beauty"
    assert constraint.observable == "BR(B_s -> tau+ tau-), BR(B0 -> tau+ tau-)"


def test_anchor_matches_yaml_and_budget_rows():
    constraint = fcc.get(_PID)
    bs_limit = _value_row("pdg2026_bs_tautau_95cl")
    bd_limit = _value_row("pdg2026_bd_tautau_95cl")
    bs_sm = _sm_row("BobethEtAl2014_BqTauTauSM", "SM BR(B_s -> tau+ tau-)")
    bd_sm = _sm_row("BobethEtAl2014_BqTauTauSM", "SM BR(B_d -> tau+ tau-)")

    assert constraint.anchor.bs_limit.value == pytest.approx(bs_limit["upper_limit"])
    assert constraint.anchor.bs_limit.units == bs_limit["units"]
    assert constraint.anchor.bs_limit.confidence_level == pytest.approx(
        bs_limit["confidence_level"]
    )
    assert constraint.anchor.bs_limit.snapshot_path == bs_limit["snapshot_path"]
    assert constraint.anchor.bd_limit.value == pytest.approx(bd_limit["upper_limit"])
    assert constraint.anchor.bd_limit.units == bd_limit["units"]
    assert constraint.anchor.bd_limit.confidence_level == pytest.approx(
        bd_limit["confidence_level"]
    )
    assert constraint.anchor.bd_limit.source_url == bd_limit["source_url"]
    assert constraint.anchor.bs_standard_model.value == pytest.approx(bs_sm["value"])
    assert constraint.anchor.bs_standard_model.uncertainty == pytest.approx(
        bs_sm["uncertainty"]
    )
    assert constraint.anchor.bs_standard_model.source_url == bs_sm["source_url"]
    assert constraint.anchor.bd_standard_model.value == pytest.approx(bd_sm["value"])
    assert constraint.anchor.bd_standard_model.uncertainty == pytest.approx(
        bd_sm["uncertainty"]
    )
    assert constraint.anchor.bs_budget.hard_veto_budget == pytest.approx(
        bs_limit["upper_limit"]
    )
    assert constraint.anchor.bd_budget.hard_veto_budget == pytest.approx(
        bd_limit["upper_limit"]
    )
    assert constraint.anchor.bs_budget.sm_anchor_to_limit_ratio == pytest.approx(
        bs_sm["value"] / bs_limit["upper_limit"]
    )
    assert constraint.anchor.bd_budget.sm_anchor_to_limit_ratio == pytest.approx(
        bd_sm["value"] / bd_limit["upper_limit"]
    )


def test_b008_anchors_route_through_scaffold_load_anchor(monkeypatch):
    calls: list[tuple[str, ...]] = []
    original_load_anchor = b008_module.load_anchor

    def spy_load_anchor(*args, **kwargs):
        calls.append(tuple(kwargs["candidates"]))
        return original_load_anchor(*args, **kwargs)

    constraint = fcc.get(_PID)
    monkeypatch.setattr(b008_module, "load_anchor", spy_load_anchor)
    anchor = b008_module._load_b008_anchor(
        _PID,
        formula_sm_bs=constraint.bs_sm_result.branching_fraction,
        formula_sm_bd=constraint.bd_sm_result.branching_fraction,
    )

    assert calls == [
        ("values[0]",),
        ("values[1]",),
        ("auxiliary_theory_inputs[0]",),
        ("auxiliary_theory_inputs[1]",),
    ]
    assert anchor.bs_limit.value == pytest.approx(constraint.anchor.bs_limit.value)
    assert anchor.bd_limit.value == pytest.approx(constraint.anchor.bd_limit.value)


def test_loud_fail_for_mismatched_limit_anchor_block(monkeypatch):
    original_load_anchor = b008_module.load_anchor

    def wrong_block_load_anchor(*args, **kwargs):
        anchor = original_load_anchor(*args, **kwargs)
        return replace(anchor, block_key="values[999]")

    monkeypatch.setattr(b008_module, "load_anchor", wrong_block_load_anchor)

    with pytest.raises(AnchorError, match=r"expected 'values\[0\]' for B008 anchor"):
        b008_module._limit_anchor_from_value_row(
            _PID,
            row_id="pdg2026_bs_tautau_95cl",
        )


def test_loud_fail_for_missing_limit_row_and_limit_below_sm():
    constraint = fcc.get(_PID)

    with pytest.raises(AnchorError, match="expected exactly one"):
        b008_module._value_row_by_id(_PID, "not_a_real_row")

    bad_sm = replace(
        constraint.anchor.bs_standard_model,
        anchor=replace(
            constraint.anchor.bs_standard_model.anchor,
            value=constraint.anchor.bs_limit.value * 1.01,
        ),
    )
    with pytest.raises(AnchorError, match="limit must exceed the SM anchor"):
        b008_module._build_channel_budget(
            process_id=_PID,
            channel_key="b_s",
            limit=constraint.anchor.bs_limit,
            standard_model=bad_sm,
            formula_sm=constraint.bs_sm_result.branching_fraction,
        )


def test_evaluate_without_input_degrades_gracefully():
    constraint = fcc.get(_PID)
    result = constraint.evaluate(point_builder.empty_point())
    old_style_result = constraint.evaluate(
        point_builder.build_from_quark_couplings(_bq_couplings(bs_left=1.0))
    )

    assert result.process_id == _PID
    assert result.passes is True
    assert result.predicted is None
    assert result.ratio is None
    assert result.experimental == pytest.approx(constraint.anchor.value)
    assert result.sm_prediction == pytest.approx(constraint.anchor.sm_value)
    assert result.budget == pytest.approx(constraint.anchor.budget)
    assert result.diagnostics["missing_extra"] == "rs_semileptonic_wilsons"
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]
    assert old_style_result.passes is True
    assert old_style_result.predicted is None
    assert old_style_result.diagnostics["evaluated"] is False
    assert old_style_result.diagnostics["missing_extra"] == "rs_semileptonic_wilsons"
    assert old_style_result.diagnostics["legacy_quark_mass_basis_couplings_present"] is True


def test_sm_limit_branching_fractions_match_independent_formula_and_anchors():
    constraint = fcc.get(_PID)
    point = sm_limit_rare_b_point()
    result = constraint.evaluate(point)
    channels = result.diagnostics["channels"]
    manual_bs = _manual_sm_tautau(constraint.sm_inputs, "b_s")
    manual_bd = _manual_sm_tautau(constraint.sm_inputs, "b_d")

    assert channels["b_s"]["predicted_branching_fraction"] == pytest.approx(manual_bs)
    assert channels["b_d"]["predicted_branching_fraction"] == pytest.approx(manual_bd)
    assert manual_bs == pytest.approx(7.737985313074278e-07)
    assert manual_bd == pytest.approx(2.1909728194059746e-08)
    assert abs(manual_bs - constraint.anchor.bs_standard_model.value) < (
        constraint.anchor.bs_standard_model.uncertainty
    )
    assert abs(manual_bd - constraint.anchor.bd_standard_model.value) < (
        constraint.anchor.bd_standard_model.uncertainty
    )
    assert result.predicted == pytest.approx(manual_bs)
    assert result.sm_prediction == pytest.approx(manual_bs)
    assert result.diagnostics["active_channel"] == "b_s"
    assert channels["b_s"]["c9_np"] == pytest.approx(0.0j)
    assert channels["b_d"]["c10_np"] == pytest.approx(0.0j)
    assert channels["b_s"]["wilson_prefactor_reused"] is False
    assert channels["b_d"]["second_mkk_suppression_applied"] is False
    assert result.passes is True


def test_evaluate_runs_end_to_end_with_real_fields_and_complex_diagnostics():
    result = fcc.get(_PID).evaluate(sample_rare_b_point())

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
    for channel in ("b_s", "b_d"):
        diagnostics = result.diagnostics["channels"][channel]
        for key in (
            "left_qb_coupling",
            "right_qb_coupling",
            "lambda_t",
            "c10_total",
            "c10_leptonic_np",
            "c9_effective_np",
            "c9_np",
            "c10_np",
        ):
            assert isinstance(diagnostics[key], complex)
        for key in (
            "m_kk_gev",
            "matching_scale_gev",
            "muon_vector_delta",
            "muon_axial_delta",
            "prompt_branching_fraction",
            "time_integrated_factor",
            "hard_veto_np_budget",
            "total_limit_ratio",
            "upward_excess_over_sm_anchor",
        ):
            assert isinstance(diagnostics[key], float)
            assert math.isfinite(diagnostics[key])
        assert diagnostics["c9_does_not_enter_leptonic_rate"] is True
        assert diagnostics["rs_semileptonic_wilsons_present"] is True
        assert diagnostics["rs_semileptonic_matching_status"].endswith(
            "no_second_1_over_M_KK_squared"
        )
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]


@pytest.mark.parametrize(
    ("scale", "expected_pass"),
    [
        (1.0, True),
        (1.0e4, False),
    ],
)
def test_safe_point_passes_and_large_np_point_fails(
    scale: float,
    expected_pass: bool,
):
    constraint = fcc.get(_PID)
    point = scaled_rare_b_point(scale=scale)
    result = constraint.evaluate(point)
    direct_bs = evaluate_bq_to_mumu(
        core_wilsons_from_rs_coeff(
            rs_coeff(point, transition="b_s", lepton="tau"),
            matching_scale_gev=point.extras["kk_ew_mass_gev"],
        ),
        transition="b_s",
        inputs=constraint.sm_inputs,
    )
    direct_bd = evaluate_bq_to_mumu(
        core_wilsons_from_rs_coeff(
            rs_coeff(point, transition="b_d", lepton="tau"),
            matching_scale_gev=point.extras["kk_ew_mass_gev"],
        ),
        transition="b_d",
        inputs=constraint.sm_inputs,
    )
    expected_bs_ratio = (
        max(0.0, direct_bs.branching_fraction - constraint.anchor.bs_standard_model.value)
        / constraint.anchor.bs_budget.hard_veto_budget
    )
    expected_bd_ratio = (
        max(0.0, direct_bd.branching_fraction - constraint.anchor.bd_standard_model.value)
        / constraint.anchor.bd_budget.hard_veto_budget
    )

    assert result.passes is expected_pass
    assert result.ratio == pytest.approx(max(expected_bs_ratio, expected_bd_ratio))
    assert result.diagnostics["channels"]["b_s"]["predicted_branching_fraction"] == (
        pytest.approx(direct_bs.branching_fraction)
    )
    assert result.diagnostics["channels"]["b_d"]["predicted_branching_fraction"] == (
        pytest.approx(direct_bd.branching_fraction)
    )
    if expected_pass:
        assert result.ratio < 1.0
    else:
        assert result.ratio > 1.0
        assert result.diagnostics["active_channel"] == "b_d"


def test_optional_kk_ew_mass_extra_changes_matching_scale():
    default_point = scaled_rare_b_point(scale=20.0)
    ew_point = scaled_rare_b_point(scale=20.0, kk_ew_mass_gev=6000.0)
    default_result = fcc.get(_PID).evaluate(default_point)
    ew_result = fcc.get(_PID).evaluate(ew_point)

    for channel in ("b_s", "b_d"):
        default_channel = default_result.diagnostics["channels"][channel]
        ew_channel = ew_result.diagnostics["channels"][channel]
        assert default_channel["m_kk_gev"] == pytest.approx(3000.0)
        assert ew_channel["m_kk_gev"] == pytest.approx(6000.0)
        assert ew_channel["c10_leptonic_np"] == pytest.approx(
            default_channel["c10_leptonic_np"]
        )
        assert ew_channel["second_mkk_suppression_applied"] is False
    assert ew_result.diagnostics["kk_ew_mass_extra_used"] is True
    assert ew_result.predicted == pytest.approx(default_result.predicted)


def test_evaluate_is_pure_and_deterministic():
    point = sample_rare_b_point()
    before_bundle = point.extras["rs_semileptonic_wilsons"]
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    assert point.extras["rs_semileptonic_wilsons"] is before_bundle
