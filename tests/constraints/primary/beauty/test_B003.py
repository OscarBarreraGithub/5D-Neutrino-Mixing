"""Production tests for B003 (Delta m_s / B_s mixing)."""

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
from flavor_catalog_constraints.primary.beauty import B003 as b003_module
from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.deltaf2 import (
    DEFAULT_DELTA_F2_INPUTS_V1,
    DELTA_M_BS_EXP,
    DELTA_M_BS_SM,
    compute_delta_f2_wilsons,
    evaluate_bs_mixing_with_running,
)

_PID = "B003"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "beauty" / "B003.yaml"
_HBAR_GEV_PER_PS = 6.582119569e-13


def _yaml():
    with open(_SIDECAR) as handle:
        return yaml.safe_load(handle)


def _yaml_pdg_block():
    return _yaml()["pdg_or_equivalent"]


def _yaml_auxiliary_theory_block():
    return _yaml()["auxiliary_theory_inputs"]


def _bs_couplings(
    left: complex,
    right: complex,
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


def _bs_wilsons(couplings: QuarkMassBasisCouplings):
    wilsons = compute_delta_f2_wilsons(
        couplings,
        inputs=DEFAULT_DELTA_F2_INPUTS_V1,
    )
    return next(w for w in wilsons if w.input.key == "b_s")


def _audited_bs_with_running(couplings: QuarkMassBasisCouplings):
    return evaluate_bs_mixing_with_running(_bs_wilsons(couplings), mu_had=4.18)


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.HARD
    assert constraint.family == "beauty"
    assert constraint.observable == "Delta m_s"


def test_anchor_matches_yaml():
    constraint = fcc.get(_PID)
    pdg = _yaml_pdg_block()
    aux = _yaml_auxiliary_theory_block()
    exp = pdg["canonical_hflav_recommended"]
    flag = aux["flag_2024_bmixing"]
    fbs = flag["f_Bs_sqrt_Bhat_Bs"]
    xi = flag["xi"]

    assert constraint.anchor.experimental.value == pytest.approx(exp["value"])
    assert constraint.anchor.experimental.uncertainty == pytest.approx(
        exp["uncertainty"]
    )
    assert constraint.anchor.experimental.units == exp["units"]
    assert constraint.anchor.experimental.source_url == exp["source_url"]
    assert constraint.anchor.flag_bmixing.source_url == flag["source_url"]
    assert constraint.anchor.flag_bmixing.f_bs_sqrt_bhat_bs_value_mev == (
        pytest.approx(fbs["value"])
    )
    assert constraint.anchor.flag_bmixing.f_bs_sqrt_bhat_bs_uncertainty_mev == (
        pytest.approx(fbs["uncertainty"])
    )
    assert constraint.anchor.flag_bmixing.xi_value == pytest.approx(xi["value"])
    assert constraint.anchor.flag_bmixing.xi_uncertainty == pytest.approx(
        xi["uncertainty"]
    )

    exp_delta_m_gev = exp["value"] * _HBAR_GEV_PER_PS
    exp_sigma_delta_m_gev = exp["uncertainty"] * _HBAR_GEV_PER_PS
    sm_theory_sigma = DELTA_M_BS_SM * 2.0 * fbs["uncertainty"] / fbs["value"]
    combined_sigma = math.sqrt(
        exp_sigma_delta_m_gev * exp_sigma_delta_m_gev
        + sm_theory_sigma * sm_theory_sigma
    )
    central_budget = abs(exp_delta_m_gev - DELTA_M_BS_SM) / 2.0
    loose_budget = (abs(exp_delta_m_gev - DELTA_M_BS_SM) + combined_sigma) / 2.0
    core_default_budget = loose_budget
    core_legacy_budget = DELTA_M_BS_EXP / 2.0

    assert constraint.anchor.budget_band.experimental_delta_m_gev == pytest.approx(
        exp_delta_m_gev
    )
    assert constraint.anchor.budget_band.sm_delta_m_gev == pytest.approx(
        DELTA_M_BS_SM
    )
    assert constraint.anchor.central_budget == pytest.approx(central_budget)
    assert constraint.anchor.budget_band.sm_theory_sigma_delta_m_gev == (
        pytest.approx(sm_theory_sigma)
    )
    assert constraint.anchor.budget_band.combined_sigma_delta_m_gev == pytest.approx(
        combined_sigma
    )
    assert constraint.anchor.budget_band.loose_budget == pytest.approx(loose_budget)
    assert constraint.anchor.budget == pytest.approx(loose_budget)
    assert constraint.anchor.budget == pytest.approx(2.635167648629676e-13)
    assert constraint.anchor.budget_band.core_default_m12_budget == pytest.approx(
        core_default_budget
    )
    assert constraint.anchor.budget_band.core_legacy_m12_budget == pytest.approx(
        core_legacy_budget
    )


@pytest.mark.parametrize(
    ("metadata", "bad_value", "error_match"),
    [
        ("units", "GeV", "must use units"),
        ("observable", "DeltaGamma_s", "must have observable"),
    ],
)
def test_delta_ms_anchor_metadata_is_validated_before_budget_construction(
    monkeypatch: pytest.MonkeyPatch,
    metadata: str,
    bad_value: str,
    error_match: str,
):
    anchor_kwargs = {
        "process_id": _PID,
        "block_key": "bad_anchor",
        "value": 17.766,
        "uncertainty": 0.006,
        "observable": "Delta m_s",
        "units": "ps^-1",
    }
    anchor_kwargs[metadata] = bad_value
    bad_anchor = Anchor(**anchor_kwargs)

    monkeypatch.setattr(b003_module, "load_anchor", lambda *args, **kwargs: bad_anchor)

    def fail_if_reached(*args, **kwargs):
        pytest.fail("budget inputs loaded before experimental anchor validation")

    monkeypatch.setattr(b003_module, "_load_flag_bmixing_inputs", fail_if_reached)

    with pytest.raises(AnchorError, match=error_match):
        b003_module._load_bs_mixing_anchor(_PID)


def test_missing_auxiliary_anchor_error_names_auxiliary_context(
    monkeypatch: pytest.MonkeyPatch,
):
    data = _yaml()
    data["auxiliary_theory_inputs"] = {"unexpected_auxiliary_block": {}}
    monkeypatch.setattr(b003_module, "load_full_yaml", lambda *args, **kwargs: data)

    with pytest.raises(AnchorError) as exc_info:
        b003_module._load_flag_bmixing_inputs(_PID)

    message = str(exc_info.value)
    assert "auxiliary_theory_inputs" in message
    assert "pdg_or_equivalent" not in message


def test_evaluate_without_input_degrades_gracefully():
    constraint = fcc.get(_PID)
    result = constraint.evaluate(point_builder.empty_point())

    assert result.process_id == _PID
    assert result.passes is True
    assert result.predicted is None
    assert result.ratio is None
    assert result.experimental == pytest.approx(
        constraint.anchor.budget_band.experimental_delta_m_gev
    )
    assert result.sm_prediction == pytest.approx(
        constraint.anchor.budget_band.sm_delta_m_gev
    )
    assert result.budget == pytest.approx(constraint.anchor.budget)
    assert "Delta m_s in GeV" in result.notes
    assert "|M12^NP| in GeV" in result.notes
    assert result.diagnostics["missing_extra"] == "quark_mass_basis_couplings"


def test_evaluate_runs_end_to_end_with_real_couplings_and_real_finite_fields():
    couplings = _bs_couplings(left=1.0e-3 + 0.5e-3j, right=1.0e-3 + 0.2e-3j)
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
    assert isinstance(result.diagnostics["left_sb_coupling"], complex)
    assert isinstance(result.diagnostics["wilson_coefficients"]["C4_LR"], complex)
    for key in (
        "abs_m12_np_gev",
        "delta_m_np_gev",
        "hadronic_scale_gev",
        "matching_scale_gev",
        "m_kk_gev",
        "central_np_budget",
        "loose_band_np_budget",
        "hard_veto_np_budget",
        "experimental_delta_m_ps_inv",
        "experimental_delta_m_gev",
        "experimental_m12_gev",
        "sm_delta_m_ps_inv",
        "sm_delta_m_gev",
        "sm_m12_gev",
        "budget_combined_sigma_delta_m_gev",
        "budget_sm_theory_sigma_delta_m_gev",
        "budget_experimental_sigma_delta_m_gev",
        "flag_f_bs_sqrt_bhat_relative_sigma",
        "flag_f_bs_sqrt_bhat_bs_mev",
        "flag_f_bs_sqrt_bhat_bs_uncertainty_mev",
        "core_delta_m_exp_gev",
        "core_delta_m_sm_gev",
        "core_default_m12_budget",
        "core_legacy_m12_budget",
        "core_default_ratio_to_budget",
        "core_legacy_ratio_to_budget",
    ):
        assert isinstance(result.diagnostics[key], float)
        assert math.isfinite(result.diagnostics[key])
    assert result.diagnostics["qcd_running_applied"] is True
    assert result.diagnostics["hadronic_scale_gev"] == pytest.approx(4.18)
    assert result.diagnostics["budget_policy_id"] == (
        "b_s_delta_m_flag2024_hflav2024_one_sigma_v1"
    )
    assert result.diagnostics["confidence_level"] == "68.27% one_sigma_sensitivity"
    assert result.diagnostics["system"] == "B_s"
    assert result.sm_prediction == pytest.approx(
        result.diagnostics["sm_delta_m_gev"]
    )
    assert result.experimental == pytest.approx(
        result.diagnostics["experimental_delta_m_gev"]
    )
    assert result.diagnostics["sm_m12_gev"] == pytest.approx(
        result.sm_prediction / 2.0
    )
    assert result.diagnostics["experimental_m12_gev"] == pytest.approx(
        result.experimental / 2.0
    )
    assert "|M12^NP| in GeV" in result.notes
    assert "Delta m_s in GeV" in result.notes


def test_reference_couplings_use_uncertainty_budget_not_core_legacy_budget():
    couplings = _bs_couplings(left=1.0e-3 + 0.5e-3j, right=1.0e-3 + 0.2e-3j)
    constraint = fcc.get(_PID)
    result = constraint.evaluate(point_builder.build_from_quark_couplings(couplings))
    audited = _audited_bs_with_running(couplings)

    assert audited.budget == pytest.approx(2.6351676486296766e-13)
    # The promoted core default now matches the B003 catalog one-sigma budget.
    # The former full-Delta-m envelope is retained only in diagnostics.
    assert result.budget == pytest.approx(2.635167648629676e-13)
    # Re-pinned after M-6: B_s Wilsons run to m_b=4.18 GeV.
    assert result.predicted == pytest.approx(3.6767250485257285e-14)
    assert result.ratio == pytest.approx(0.13952528031518893)
    assert result.predicted == pytest.approx(audited.abs_m12_np)
    assert result.diagnostics["core_default_m12_budget"] == pytest.approx(
        audited.budget
    )
    assert result.diagnostics["core_default_ratio_to_budget"] == pytest.approx(
        audited.ratio_to_budget
    )
    assert result.diagnostics["core_legacy_m12_budget"] == pytest.approx(
        DELTA_M_BS_EXP / 2.0
    )
    assert result.diagnostics["core_legacy_ratio_to_budget"] == pytest.approx(
        result.predicted / (DELTA_M_BS_EXP / 2.0)
    )
    assert result.ratio == pytest.approx(result.predicted / constraint.anchor.budget)


@pytest.mark.parametrize(
    ("couplings", "expected_pass"),
    [
        (_bs_couplings(left=1.0e-3 + 0.5e-3j, right=1.0e-3 + 0.2e-3j), True),
        (_bs_couplings(left=1.0e-2 + 0.5e-2j, right=1.0e-2 + 0.2e-2j), False),
    ],
)
def test_pass_fail_and_numbers_match_running_evaluator_prediction(
    couplings: QuarkMassBasisCouplings,
    expected_pass: bool,
):
    constraint = fcc.get(_PID)
    point = point_builder.build_from_quark_couplings(couplings)
    result = constraint.evaluate(point)
    audited = _audited_bs_with_running(couplings)

    assert result.passes is expected_pass
    assert result.predicted == pytest.approx(audited.abs_m12_np)
    assert result.budget == pytest.approx(constraint.anchor.budget)
    assert result.ratio == pytest.approx(audited.abs_m12_np / constraint.anchor.budget)
    assert result.diagnostics["core_default_m12_budget"] == pytest.approx(
        audited.budget
    )
    assert result.diagnostics["core_default_ratio_to_budget"] == pytest.approx(
        audited.ratio_to_budget
    )
    assert result.diagnostics["core_legacy_m12_budget"] == pytest.approx(
        DELTA_M_BS_EXP / 2.0
    )
    assert result.diagnostics["core_legacy_ratio_to_budget"] == pytest.approx(
        result.predicted / (DELTA_M_BS_EXP / 2.0)
    )
    for value in (audited.abs_m12_np, audited.ratio_to_budget, audited.budget):
        assert isinstance(value, float)
        assert math.isfinite(value)
    if expected_pass:
        assert result.ratio <= 1.0
    else:
        assert result.ratio > 10.0


def test_evaluate_is_pure_and_deterministic():
    couplings = _bs_couplings(left=1.0e-3 + 0.5e-3j, right=1.0e-3 + 0.2e-3j)
    before_left_down = couplings.left_down.copy()
    before_right_down = couplings.right_down.copy()
    point = point_builder.build_from_quark_couplings(couplings)
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    np.testing.assert_array_equal(couplings.left_down, before_left_down)
    np.testing.assert_array_equal(couplings.right_down, before_right_down)
