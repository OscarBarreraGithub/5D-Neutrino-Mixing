"""Production tests for C006 (D0 -> e+- mu-+ LFV)."""

from __future__ import annotations

import math
from pathlib import Path

import pytest
import yaml

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import point_builder
from flavor_catalog_constraints.base import ConstraintProtocol, Severity
from tests.constraints.lfv_rare_phase4c_helpers import (
    diagonal_lfv_rare_point,
    lfv_coeff,
    lfv_live_rare_point,
    manual_d0_emu_rate,
    scaled_lfv_rare_point,
)

_PID = "C006"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "charm" / "C006.yaml"


def _yaml():
    with open(_SIDECAR) as handle:
        return yaml.safe_load(handle)


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.HARD
    assert constraint.family == "charm"
    assert constraint.observable == "BR(D0 -> e+- mu-+) LFV"


def test_anchor_matches_yaml_and_budget():
    constraint = fcc.get(_PID)
    data = _yaml()
    pdg = data["pdg_or_equivalent"]
    current = pdg["canonical_current_limit"]
    lhcb = pdg["lhcb_primary_result"]
    belle = pdg["previous_limits"]["belle_2010"]
    babar = pdg["previous_limits"]["babar_2012"]
    rs = data["paper_era_reference"]["rs_baseline"]

    assert constraint.anchor.current_limit.value == pytest.approx(current["value"])
    assert constraint.anchor.current_limit.confidence_level == pytest.approx(
        current["confidence_level"]
    )
    assert constraint.anchor.current_limit.source_url == current["source_url"]
    assert constraint.anchor.lhcb_primary_limit.value == pytest.approx(lhcb["value"])
    assert constraint.anchor.lhcb_primary_limit.integrated_luminosity_fb_inv == (
        pytest.approx(lhcb["integrated_luminosity_fb_inv"])
    )
    assert constraint.anchor.lhcb_primary_limit.collision_energy_tev == tuple(
        float(item) for item in lhcb["collision_energy_TeV"]
    )
    assert constraint.anchor.belle_previous_limit.value == pytest.approx(belle["value"])
    assert constraint.anchor.babar_previous_limit.value == pytest.approx(babar["value"])
    assert constraint.anchor.rs_baseline.generic_rs_kk_gluon_scale_tev == pytest.approx(
        rs["generic_rs_kk_gluon_scale_TeV"]
    )
    assert (
        constraint.anchor.rs_baseline.composite_pseudo_goldstone_kk_gluon_scale_tev
        == pytest.approx(rs["composite_pseudo_goldstone_kk_gluon_scale_TeV"])
    )
    assert constraint.anchor.budget == pytest.approx(current["value"])
    assert constraint.anchor.budget == pytest.approx(1.3e-8)


def test_anchor_loading_fails_loudly_for_missing_candidate_or_value():
    with pytest.raises(fcc.AnchorError, match="none of the expected anchor keys"):
        fcc.load_anchor(_PID, family="charm", candidates=("missing_c006_anchor",))
    with pytest.raises(fcc.AnchorError, match="has no 'missing_limit' field"):
        fcc.load_anchor(
            _PID,
            family="charm",
            candidates=("canonical_current_limit",),
            value_key="missing_limit",
        )


def test_absent_rs_semileptonic_wilsons_degrades_gracefully():
    constraint = fcc.get(_PID)
    result = constraint.evaluate(point_builder.empty_point())

    assert result.process_id == _PID
    assert result.passes is True
    assert result.predicted is None
    assert result.ratio is None
    assert result.sm_prediction == pytest.approx(0.0)
    assert result.experimental == pytest.approx(constraint.anchor.value)
    assert result.budget == pytest.approx(constraint.anchor.budget)
    assert result.notes.startswith("NOT EVALUATED -")
    assert result.diagnostics["evaluated"] is False
    assert result.diagnostics["missing_extra"] == "rs_semileptonic_wilsons"


def test_invalid_rs_semileptonic_wilsons_is_unevaluated_not_real_pass():
    result = fcc.get(_PID).evaluate(
        point_builder.make_point(rs_semileptonic_wilsons=object())
    )

    assert result.passes is True
    assert result.predicted is None
    assert result.ratio is None
    assert result.notes.startswith("NOT EVALUATED -")
    assert result.diagnostics["evaluated"] is False
    assert result.diagnostics["invalid_extra"] == "rs_semileptonic_wilsons"
    assert result.diagnostics["exception_type"] == "ValueError"


def test_diagonal_rs_ew_builder_gives_zero_tree_lfv_and_resolves_proxy_flag():
    result = fcc.get(_PID).evaluate(diagonal_lfv_rare_point())
    coeff = lfv_coeff(diagonal_lfv_rare_point(), "c_to_u")

    assert coeff.c9_lfv_np == pytest.approx(0.0j, abs=1.0e-24)
    assert result.predicted == pytest.approx(0.0)
    assert result.ratio == pytest.approx(0.0)
    assert result.passes is True
    assert result.diagnostics["evaluated"] is True
    assert result.diagnostics["tree_level_matching_status"] == (
        "rigorous_tree_light_z_lfv_llqq_from_rs_semileptonic_wilsons"
    )
    assert "diagonal charged-lepton fit" in result.diagnostics["lfv_tree_level_note"]
    assert "NEEDS-HUMAN-PHYSICS" not in str(result.diagnostics)
    assert "proxy" not in result.notes.lower()


def test_lfv_live_toy_matches_independent_rate_and_mkk_scaling():
    constraint = fcc.get(_PID)
    low_point = lfv_live_rare_point(3000.0)
    high_point = lfv_live_rare_point(6000.0)
    low = constraint.evaluate(low_point)
    high = constraint.evaluate(high_point)
    coeff = lfv_coeff(low_point, "c_to_u")
    manual = manual_d0_emu_rate(coeff, constraint.sm_inputs)

    assert abs(coeff.c9_lfv_np) > 0.0
    assert low.predicted == pytest.approx(manual, rel=1.0e-12, abs=1.0e-30)
    assert low.ratio == pytest.approx(low.predicted / constraint.anchor.budget)
    assert high.predicted / low.predicted == pytest.approx(
        (3000.0 / 6000.0) ** 8,
        rel=1.0e-8,
    )
    assert low.diagnostics["charge_conjugate_modes_included"] is True
    assert low.diagnostics["wilson_prefactor_reused"] is False
    assert low.diagnostics["second_mkk_suppression_applied"] is False


def test_amplified_lfv_live_toy_bites_limit():
    result = fcc.get(_PID).evaluate(
        scaled_lfv_rare_point(lfv_live_rare_point(), "c_to_u", 1.0e5)
    )

    assert result.passes is False
    assert result.ratio > 1.0


def test_evaluate_runs_end_to_end_with_real_finite_fields_and_complex_diagnostics():
    result = fcc.get(_PID).evaluate(lfv_live_rare_point())

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
        "rs_semileptonic_lambda_ckm",
        "c9_lfv_combination",
        "c10_lfv_combination",
        "c9_lfv_np",
        "c10_lfv_np",
    ):
        assert isinstance(result.diagnostics[key], complex)
    assert result.diagnostics["evaluated"] is True
    assert result.diagnostics["rs_semileptonic_wilsons_present"] is True
    assert result.diagnostics["wilson_coefficients"]


def test_evaluate_is_pure_and_deterministic():
    point = lfv_live_rare_point()
    before_bundle = point.extras["rs_semileptonic_wilsons"]
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    assert point.extras["rs_semileptonic_wilsons"] is before_bundle
