"""Production tests for B022 (B+ -> K+ nu nubar)."""

from __future__ import annotations

import math
from pathlib import Path

import pytest
import yaml

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import point_builder
from flavor_catalog_constraints.anchors import Anchor, AnchorError
from flavor_catalog_constraints.base import ConstraintProtocol, Severity
from flavor_catalog_constraints.primary.beauty import B022 as b022_module
from flavor_catalog_constraints.physics_adapters.rare_b_nunu import (
    bplus_kplus_nunu_from_rs_semileptonic_wilsons,
)
from tests.constraints.primary.nunu_phase4d_helpers import (
    direct_contact_x_np,
    nunu_block,
    rigorous_point,
    scalar_x_np,
    sm_limit_point,
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
    budget_raises = abs(float(exp["value"]) - float(sm["value"])) + combined_up
    budget_lowers = abs(float(exp["value"]) - float(sm["value"])) + combined_down
    assert constraint.anchor.budget_band.budget_raises_branching_fraction == (
        pytest.approx(budget_raises)
    )
    assert constraint.anchor.budget_band.budget_lowers_branching_fraction == (
        pytest.approx(budget_lowers)
    )
    assert constraint.anchor.budget == pytest.approx(max(budget_raises, budget_lowers))
    assert constraint.anchor.budget == pytest.approx(2.4500741486595877e-5)


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


def test_absent_rs_semileptonic_wilsons_degrades_gracefully():
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
    assert result.diagnostics["evaluated"] is False
    assert result.diagnostics["missing_extra"] == "rs_semileptonic_wilsons"
    assert result.diagnostics["budget_policy_id"] == (
        "b022_belleii2023_hpqcd2023_np_shift_one_sigma_v1"
    )
    assert result.diagnostics["confidence_level"] == "68.27% one_sigma_sensitivity"
    assert "needs_human_physics" not in result.diagnostics


def test_legacy_quark_proxy_only_is_unevaluated_not_real_pass():
    result = fcc.get(_PID).evaluate(
        point_builder.make_point(quark_mass_basis_couplings=object())
    )

    assert result.passes is True
    assert result.predicted is None
    assert result.ratio is None
    assert result.diagnostics["evaluated"] is False
    assert result.diagnostics["missing_extra"] == "rs_semileptonic_wilsons"
    assert result.diagnostics["legacy_quark_mass_basis_couplings_present"] is True


def test_invalid_rs_semileptonic_wilsons_is_unevaluated_not_real_pass():
    result = fcc.get(_PID).evaluate(
        point_builder.make_point(rs_semileptonic_wilsons=object())
    )

    assert result.passes is True
    assert result.predicted is None
    assert result.ratio is None
    assert result.diagnostics["evaluated"] is False
    assert result.diagnostics["invalid_extra"] == "rs_semileptonic_wilsons"


def test_sm_limit_universal_quark_c_recovers_committed_sm_branching_fraction():
    constraint = fcc.get(_PID)
    point = sm_limit_point()
    result = constraint.evaluate(point)
    block = nunu_block(point, "b_to_s")
    x_left, x_right, x_total = scalar_x_np(block)

    assert x_left == pytest.approx(0.0j, abs=1.0e-18)
    assert x_right == pytest.approx(0.0j, abs=1.0e-18)
    assert x_total == pytest.approx(0.0j, abs=1.0e-18)
    assert result.diagnostics["x_np_total"] == pytest.approx(0.0j, abs=1.0e-18)
    assert result.predicted == pytest.approx(5.58e-6)
    assert result.predicted == pytest.approx(result.sm_prediction)
    assert result.ratio == pytest.approx(0.0, abs=1.0e-12)
    assert result.diagnostics["np_shift_from_sm_branching_fraction"] == pytest.approx(
        0.0,
        abs=1.0e-18,
    )
    assert result.diagnostics["budget_policy_id"] == (
        "b022_belleii2023_hpqcd2023_np_shift_one_sigma_v1"
    )
    assert result.diagnostics["confidence_level"] == "68.27% one_sigma_sensitivity"
    assert result.passes is True


def test_rigorous_nonzero_xnp_shifts_br_and_matches_direct_contact():
    constraint = fcc.get(_PID)
    point = rigorous_point()
    result = constraint.evaluate(point)
    direct = bplus_kplus_nunu_from_rs_semileptonic_wilsons(
        point.extras["rs_semileptonic_wilsons"],
        inputs=constraint.sm_inputs,
    )
    x_left, x_right, x_total = direct_contact_x_np(point, "b_to_s")

    assert abs(x_total) > 1.0e-2
    assert result.diagnostics["x_np_left"] == pytest.approx(x_left)
    assert result.diagnostics["x_np_right"] == pytest.approx(x_right)
    assert result.diagnostics["x_np_total"] == pytest.approx(x_total)
    assert result.predicted == pytest.approx(direct.branching_fraction)
    assert result.predicted != pytest.approx(result.sm_prediction)
    assert result.diagnostics["nunu_mapping"] == "X_NP=C/g_SM^2"
    assert result.diagnostics["wilson_prefactor_reused"] is False
    assert result.diagnostics["second_mkk_suppression_applied"] is False
    assert result.diagnostics["legacy_one_z_proxy_reused"] is False
    assert "needs_human_physics" not in result.diagnostics


def test_evaluate_runs_end_to_end_with_real_finite_fields_and_complex_diagnostics():
    result = fcc.get(_PID).evaluate(rigorous_point())

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
        "budget_raises_branching_fraction",
        "budget_lowers_branching_fraction",
        "np_shift_from_sm_branching_fraction",
    ):
        assert isinstance(result.diagnostics[key], float)
        assert math.isfinite(result.diagnostics[key])
    assert result.diagnostics["down_sector_indices"] == (1, 2)
    assert result.diagnostics["evaluated"] is True


def test_kk_ew_mass_extra_is_diagnostic_only_no_second_mkk_suppression():
    base_point = rigorous_point()
    ew_point = point_builder.make_point(
        rs_semileptonic_wilsons=base_point.extras["rs_semileptonic_wilsons"],
        kk_ew_mass_gev=6000.0,
    )
    default_result = fcc.get(_PID).evaluate(base_point)
    ew_result = fcc.get(_PID).evaluate(ew_point)

    assert default_result.diagnostics["m_kk_gev"] == pytest.approx(3000.0)
    assert ew_result.diagnostics["m_kk_gev"] == pytest.approx(6000.0)
    assert ew_result.diagnostics["kk_ew_mass_extra_used"] is True
    assert ew_result.predicted == pytest.approx(default_result.predicted)
    assert ew_result.diagnostics["x_np_total"] == pytest.approx(
        default_result.diagnostics["x_np_total"]
    )


def test_majorana_and_dirac_active_nu_rates_match():
    dirac = fcc.get(_PID).evaluate(rigorous_point(alpha=0.0, beta=0.0))
    majorana = fcc.get(_PID).evaluate(rigorous_point(alpha=1.1, beta=-0.7))

    assert majorana.predicted == pytest.approx(dirac.predicted, rel=0.0, abs=1.0e-18)
    assert majorana.diagnostics["x_np_total"] == pytest.approx(
        dirac.diagnostics["x_np_total"],
        abs=1.0e-18,
    )
    assert majorana.diagnostics["majorana_dirac_rate_factor"] == pytest.approx(1.0)


def test_evaluate_is_pure_and_deterministic():
    point = rigorous_point()
    before_bundle = point.extras["rs_semileptonic_wilsons"]
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    assert point.extras["rs_semileptonic_wilsons"] is before_bundle
