"""Production tests for K001 (epsilon_K)."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest
import yaml

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import point_builder
from flavor_catalog_constraints.base import ConstraintProtocol, Severity
from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.deltaf2 import (
    DEFAULT_DELTA_F2_INPUTS_V1,
    compute_delta_f2_wilsons,
    evaluate_epsilon_k,
    evaluate_epsilon_k_with_running,
)

_PID = "K001"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "kaon" / "K001.yaml"


def _yaml_pdg_block():
    with open(_SIDECAR) as handle:
        return yaml.safe_load(handle)["pdg_or_equivalent"]


def _sd_couplings(
    left: complex,
    right: complex,
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


def _kaon_wilsons(couplings: QuarkMassBasisCouplings):
    wilsons = compute_delta_f2_wilsons(
        couplings,
        inputs=DEFAULT_DELTA_F2_INPUTS_V1,
    )
    return next(w for w in wilsons if w.input.key == "epsilon_k")


def _audited_epsilon_k_with_running(
    couplings: QuarkMassBasisCouplings,
    budget: float,
):
    return evaluate_epsilon_k_with_running(
        _kaon_wilsons(couplings),
        mu_had=2.0,
        epsilon_k_np_budget_override=budget,
    )


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.HARD
    assert constraint.family == "kaon"
    assert constraint.observable == "epsilon_K"


def test_anchor_matches_yaml():
    constraint = fcc.get(_PID)
    pdg = _yaml_pdg_block()
    exp = pdg["canonical_experimental_value"]
    sm = pdg["standard_model_reference"]
    flag = pdg["flag_bag_parameters"]

    assert constraint.anchor.experimental.value == pytest.approx(exp["value"])
    assert constraint.anchor.experimental.uncertainty == pytest.approx(exp["uncertainty"])
    assert constraint.anchor.experimental.source_url == exp["source_url"]
    assert constraint.anchor.standard_model.value == pytest.approx(sm["value"])
    assert constraint.anchor.standard_model.source_url == sm["source_url"]
    assert constraint.anchor.flag_bag_parameters.source_url == flag["source_url"]
    assert constraint.anchor.flag_bag_parameters.b_k_msbar_2gev_nf21 == (
        flag["B_K_MSbar_2GeV_Nf21"]
    )
    central_budget = abs(exp["value"] - sm["value"])
    bgs_sigma = math.sqrt(
        sum(value * value for value in sm["grouped_uncertainties"].values())
    )
    combined_sigma = math.sqrt(
        bgs_sigma * bgs_sigma
        + exp["uncertainty"] * exp["uncertainty"]
        + (0.15e-3) * (0.15e-3)
    )
    loose_budget = abs(exp["value"] - (sm["value"] - combined_sigma))

    assert constraint.anchor.central_budget == pytest.approx(central_budget)
    assert constraint.anchor.budget_band.sm_theory_sigma == pytest.approx(bgs_sigma)
    assert constraint.anchor.budget_band.combined_sigma == pytest.approx(
        combined_sigma
    )
    assert constraint.anchor.budget_band.tight_budget == pytest.approx(
        exp["uncertainty"]
    )
    assert constraint.anchor.budget_band.loose_budget == pytest.approx(loose_budget)
    assert constraint.anchor.budget == pytest.approx(loose_budget)
    assert constraint.anchor.budget == pytest.approx(3.0370868171657726e-4)


def test_evaluate_without_input_degrades_gracefully():
    result = fcc.get(_PID).evaluate(point_builder.empty_point())

    assert result.process_id == _PID
    assert result.passes is True
    assert result.predicted is None
    assert result.ratio is None
    assert result.experimental == pytest.approx(fcc.get(_PID).anchor.value)
    assert result.sm_prediction == pytest.approx(fcc.get(_PID).anchor.sm_value)
    assert result.budget == pytest.approx(fcc.get(_PID).anchor.budget)
    assert result.diagnostics["missing_extra"] == "quark_mass_basis_couplings"


def test_evaluate_runs_end_to_end_with_real_couplings_and_real_finite_fields():
    couplings = _sd_couplings(left=1.0e-5 + 0.5e-5j, right=1.0e-5 + 0.2e-5j)
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
    assert isinstance(result.diagnostics["left_sd_coupling"], complex)
    assert isinstance(result.diagnostics["wilson_coefficients"]["C4_LR"], complex)
    for key in (
        "im_m12_np_gev",
        "hadronic_scale_gev",
        "matching_scale_gev",
        "m_kk_gev",
        "central_np_budget",
        "tight_band_np_budget",
        "loose_band_np_budget",
        "hard_veto_np_budget",
        "budget_combined_sigma",
        "budget_sm_theory_sigma",
        "budget_experimental_sigma",
        "budget_sm_choice_sigma",
        "budget_sm_at_loose_edge",
        "budget_sm_at_tight_edge",
    ):
        assert isinstance(result.diagnostics[key], float)
        assert math.isfinite(result.diagnostics[key])
    assert result.diagnostics["qcd_running_applied"] is True
    assert result.diagnostics["hadronic_scale_gev"] == pytest.approx(2.0)


def test_reference_couplings_show_qcd_running_enhancement():
    couplings = _sd_couplings(left=1.0e-5 + 0.5e-5j, right=1.0e-5 + 0.2e-5j)
    wilsons = _kaon_wilsons(couplings)
    central_budget = fcc.get(_PID).anchor.central_budget

    unrun = evaluate_epsilon_k(
        wilsons,
        epsilon_k_np_budget_override=central_budget,
    )
    run = evaluate_epsilon_k_with_running(
        wilsons,
        mu_had=2.0,
        epsilon_k_np_budget_override=central_budget,
    )

    # Re-pinned after B3 (GGMS Eq. 8 O4/O5 un-swap + 1/(2 m_M)); fixed _sd_couplings
    # (no SVD pipeline) so B2 does not enter.  The run (post-RG, physical) ratio
    # moves ~x1.73 as expected; the unrun ratio moves more (the swap annihilated
    # the leading chiral term pre-running, PLAN §0.2), both legitimate snapshots
    # downstream of the test_epsilon_k_physics literature-anchored O4/O5 pins.
    assert unrun.ratio_to_budget == pytest.approx(1.4067965606347435)
    assert run.ratio_to_budget == pytest.approx(4.932284791072436)


@pytest.mark.parametrize(
    ("couplings", "expected_pass"),
    [
        # expected_pass re-pinned True->False after B3 (GGMS O4/O5 un-swap +
        # 1/(2 m_M)): the corrected kaon LR contribution is ~x1.7 larger, so this
        # point now crosses the epsilon_K veto (ratio ~1.09 > 1).
        (_sd_couplings(left=1.0e-5 + 0.5e-5j, right=1.0e-5 + 0.2e-5j), False),
        (_sd_couplings(left=1.0e-4 + 0.5e-4j, right=1.0e-4 + 0.2e-4j), False),
    ],
)
def test_pass_fail_and_numbers_match_audited_evaluate_epsilon_k(
    couplings: QuarkMassBasisCouplings,
    expected_pass: bool,
):
    constraint = fcc.get(_PID)
    point = point_builder.build_from_quark_couplings(couplings)
    result = constraint.evaluate(point)
    audited = _audited_epsilon_k_with_running(couplings, constraint.anchor.budget)

    assert result.passes is expected_pass
    assert result.predicted == pytest.approx(audited.epsilon_k_np)
    assert result.ratio == pytest.approx(audited.ratio_to_budget)
    assert result.budget == pytest.approx(audited.epsilon_k_np_budget)
    assert result.diagnostics["im_m12_np_gev"] == pytest.approx(audited.im_m12_np)
    if expected_pass:
        assert result.ratio <= 1.0
    else:
        # bound relaxed 10.0 -> 1.0 after B3 (GGMS O4/O5 un-swap + 1/(2 m_M)):
        # the smaller couplings0 point now crosses the veto at ratio ~1.09
        # (couplings1 remains ~108), so both points exceed the budget.
        assert result.ratio > 1.0


def test_evaluate_is_pure_and_deterministic():
    couplings = _sd_couplings(left=1.0e-5 + 0.5e-5j, right=1.0e-5 + 0.2e-5j)
    before_left_down = couplings.left_down.copy()
    point = point_builder.build_from_quark_couplings(couplings)
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    np.testing.assert_array_equal(couplings.left_down, before_left_down)
