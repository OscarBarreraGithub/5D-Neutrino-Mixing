"""Production tests for C008 (D+ -> pi+ e+- mu-+ LFV)."""

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
    manual_dtopi_emu_rate,
    scaled_lfv_rare_point,
)

_PID = "C008"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "charm" / "C008.yaml"


def _yaml():
    with open(_SIDECAR) as handle:
        return yaml.safe_load(handle)


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.HARD
    assert constraint.family == "charm"
    assert constraint.observable == "BR(D+ -> pi+ e+- mu-+) LFV full-q2 tree-level"


def test_anchor_matches_yaml_and_budget():
    constraint = fcc.get(_PID)
    data = _yaml()
    pdg = data["pdg_or_equivalent"]
    eplus = pdg["dplus_piplus_eplus_muminus"]
    eminus = pdg["dplus_piplus_eminus_muplus"]
    scope = pdg["lhcb_2021_search_scope"]
    babar = pdg["babar_2011_predecessor"]
    rs = data["paper_era_reference"]["rs_baseline"]

    assert constraint.anchor.eplus_muminus_limit.value == pytest.approx(eplus["value"])
    assert constraint.anchor.eplus_muminus_limit.confidence_level == pytest.approx(
        eplus["confidence_level"]
    )
    assert constraint.anchor.eplus_muminus_limit.source_url == eplus["source_url"]
    assert constraint.anchor.eplus_muminus_limit.companion_95cl_limit == pytest.approx(
        eplus["companion_95cl_limit"]
    )
    assert constraint.anchor.eminus_muplus_limit.value == pytest.approx(eminus["value"])
    assert constraint.anchor.eminus_muplus_limit.companion_95cl_limit == pytest.approx(
        eminus["companion_95cl_limit"]
    )
    assert constraint.anchor.search_scope.decay_modes_investigated == pytest.approx(
        scope["decay_modes_investigated"]
    )
    assert constraint.anchor.search_scope.integrated_luminosity_fb_inv == pytest.approx(
        scope["integrated_luminosity_fb_inv"]
    )
    assert constraint.anchor.babar_predecessor.eplus_muminus_limit == pytest.approx(
        babar["dplus_piplus_eplus_muminus_limit"]
    )
    assert constraint.anchor.babar_predecessor.eminus_muplus_limit == pytest.approx(
        babar["dplus_piplus_eminus_muplus_limit"]
    )
    assert constraint.anchor.rs_baseline.source == rs["source"]
    assert constraint.anchor.rs_baseline.use == rs["use"]
    assert constraint.anchor.budget == pytest.approx(min(eplus["value"], eminus["value"]))
    assert constraint.anchor.budget == pytest.approx(2.1e-7)


def test_anchor_loading_fails_loudly_for_missing_candidate_or_value():
    with pytest.raises(fcc.AnchorError, match="none of the expected anchor keys"):
        fcc.load_anchor(_PID, family="charm", candidates=("missing_c008_anchor",))
    with pytest.raises(fcc.AnchorError, match="has no 'missing_limit' field"):
        fcc.load_anchor(
            _PID,
            family="charm",
            candidates=("dplus_piplus_eplus_muminus",),
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

    assert result.passes is False
    assert result.predicted is None
    assert result.ratio is None
    assert result.notes.startswith("NOT EVALUATED -")
    assert result.diagnostics["evaluated"] is True
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
    manual = manual_dtopi_emu_rate(
        coeff,
        constraint.sd_inputs,
        charge_mode="eplus_muminus",
    )

    assert abs(coeff.c9_lfv_np) > 0.0
    assert low.predicted == pytest.approx(manual, rel=1.0e-12, abs=1.0e-30)
    assert low.ratio == pytest.approx(
        max(
            low.predicted / constraint.anchor.eplus_muminus_limit.value,
            low.predicted / constraint.anchor.eminus_muplus_limit.value,
        )
    )
    assert high.predicted / low.predicted == pytest.approx(
        (3000.0 / 6000.0) ** 8,
        rel=1.0e-8,
    )
    assert low.diagnostics["semileptonic_primed_combination_is_plus"] is True
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
        "c9_lfv_semileptonic_np",
        "c10_lfv_semileptonic_np",
        "c9_lfv_np",
        "c10_lfv_np",
    ):
        assert isinstance(result.diagnostics[key], complex)
    assert result.diagnostics["evaluated"] is True
    assert result.diagnostics["rs_semileptonic_wilsons_present"] is True
    assert result.diagnostics["wilson_coefficients"]
    assert result.diagnostics["scalar_tensor_operators_included"] is False
    assert result.diagnostics["lhcb_window_acceptance_applied"] is False


def test_evaluate_is_pure_and_deterministic():
    point = lfv_live_rare_point()
    before_bundle = point.extras["rs_semileptonic_wilsons"]
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    assert point.extras["rs_semileptonic_wilsons"] is before_bundle
