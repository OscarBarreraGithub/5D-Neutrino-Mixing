"""Production tests for C001 (neutral D0 mixing)."""

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
    evaluate_d0_mixing_with_running,
)

_PID = "C001"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "charm" / "C001.yaml"
_MU_HAD_GEV = 3.0


def _yaml_pdg_block():
    with open(_SIDECAR) as handle:
        return yaml.safe_load(handle)["pdg_or_equivalent"]


def _uc_couplings(
    left: complex,
    right: complex,
    M_KK: float = 3000.0,
) -> QuarkMassBasisCouplings:
    """Minimal valid mass-basis couplings with only the u-c slot populated."""
    zeros = np.zeros((3, 3), dtype=np.complex128)
    left_up = zeros.copy()
    right_up = zeros.copy()
    left_up[0, 1] = left
    left_up[1, 0] = np.conj(left)
    right_up[0, 1] = right
    right_up[1, 0] = np.conj(right)
    return QuarkMassBasisCouplings(
        M_KK=M_KK,
        xi_KK=1.0,
        alpha_s=0.09,
        g_s=1.0,
        left_overlap=zeros,
        right_up_overlap=zeros,
        right_down_overlap=zeros,
        left_up=left_up,
        left_down=zeros.copy(),
        right_up=right_up,
        right_down=zeros.copy(),
    )


def _d0_wilsons(couplings: QuarkMassBasisCouplings):
    wilsons = compute_delta_f2_wilsons(
        couplings,
        inputs=DEFAULT_DELTA_F2_INPUTS_V1,
    )
    return next(w for w in wilsons if w.input.key == "d")


def _audited_d0_with_running(couplings: QuarkMassBasisCouplings):
    return evaluate_d0_mixing_with_running(
        _d0_wilsons(couplings),
        mu_had=_MU_HAD_GEV,
    )


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.HARD
    assert constraint.family == "charm"
    assert constraint.observable == "D0-D0bar mixing"


def test_anchor_matches_yaml():
    constraint = fcc.get(_PID)
    pdg = _yaml_pdg_block()
    x_d = pdg["x_D"]
    y_d = pdg["y_D"]
    delta_m_d = pdg["Delta_m_D"]

    assert constraint.anchor.x_d.value == pytest.approx(x_d["value"])
    assert constraint.anchor.x_d.uncertainty == pytest.approx(x_d["uncertainty"])
    assert constraint.anchor.x_d.source_url == x_d["source_url"]
    assert constraint.anchor.y_d.value == pytest.approx(y_d["value"])
    assert constraint.anchor.y_d.uncertainty == pytest.approx(y_d["uncertainty"])
    assert constraint.anchor.y_d.source_url == y_d["source_url"]
    assert constraint.anchor.delta_m_d_native.value == pytest.approx(
        delta_m_d["value"]
    )
    assert constraint.anchor.delta_m_d_native.uncertainty == pytest.approx(
        delta_m_d["uncertainty"]
    )
    assert constraint.anchor.delta_m_d_gev.value == pytest.approx(
        delta_m_d["value_GeV"]
    )
    assert constraint.anchor.delta_m_d_gev.uncertainty == pytest.approx(
        delta_m_d["uncertainty_GeV"]
    )
    assert constraint.anchor.delta_m_d_gev.source_url == delta_m_d["source_url"]
    assert constraint.anchor.budget == pytest.approx(delta_m_d["value_GeV"] / 2.0)
    assert constraint.anchor.budget == pytest.approx(3.281e-15)


def test_anchor_loading_fails_loudly_for_missing_candidate_or_value():
    with pytest.raises(fcc.AnchorError, match="none of the expected anchor keys"):
        fcc.load_anchor(_PID, family="charm", candidates=("missing_c001_anchor",))
    with pytest.raises(fcc.AnchorError, match="has no 'missing_value_GeV' field"):
        fcc.load_anchor(
            _PID,
            family="charm",
            candidates=("Delta_m_D",),
            value_key="missing_value_GeV",
        )


def test_evaluate_without_input_degrades_gracefully():
    result = fcc.get(_PID).evaluate(point_builder.empty_point())

    assert result.process_id == _PID
    assert result.passes is True
    assert result.predicted is None
    assert result.ratio is None
    assert result.experimental == pytest.approx(fcc.get(_PID).anchor.value)
    assert result.sm_prediction is None
    assert result.budget == pytest.approx(fcc.get(_PID).anchor.budget)
    assert result.diagnostics["missing_extra"] == "quark_mass_basis_couplings"


def test_evaluate_runs_end_to_end_with_real_couplings_and_real_finite_fields():
    couplings = _uc_couplings(
        left=1.0e-4 + 0.5e-4j,
        right=0.7e-4 + 0.2e-4j,
    )
    point = point_builder.build_from_quark_couplings(couplings)
    result = fcc.get(_PID).evaluate(point)

    assert result.process_id == _PID
    for value in (
        result.predicted,
        result.ratio,
        result.budget,
        result.experimental,
    ):
        assert isinstance(value, float)
        assert math.isfinite(value)
    assert result.sm_prediction is None
    assert isinstance(result.diagnostics["left_uc_coupling"], complex)
    assert isinstance(result.diagnostics["wilson_coefficients"]["C4_LR"], complex)
    for key in (
        "abs_m12_np_gev",
        "hadronic_scale_gev",
        "matching_scale_gev",
        "m_kk_gev",
        "delta_m_d_budget_gev",
        "delta_m_d_experimental_gev",
        "delta_m_d_uncertainty_gev",
        "delta_m_d_native_value",
        "delta_m_d_native_uncertainty",
        "x_d_percent",
        "x_d_uncertainty_percent",
        "y_d_percent",
        "y_d_uncertainty_percent",
    ):
        assert isinstance(result.diagnostics[key], float)
        assert math.isfinite(result.diagnostics[key])
    assert result.diagnostics["qcd_running_applied"] is True
    assert result.diagnostics["hadronic_scale_gev"] == pytest.approx(_MU_HAD_GEV)


@pytest.mark.parametrize(
    ("couplings", "expected_pass"),
    [
        (
            _uc_couplings(left=1.0e-4 + 0.5e-4j, right=0.7e-4 + 0.2e-4j),
            True,
        ),
        (
            _uc_couplings(left=2.0e-3 + 1.0e-3j, right=1.4e-3 + 0.4e-3j),
            False,
        ),
    ],
)
def test_pass_fail_and_numbers_match_audited_evaluate_d0_mixing(
    couplings: QuarkMassBasisCouplings,
    expected_pass: bool,
):
    constraint = fcc.get(_PID)
    point = point_builder.build_from_quark_couplings(couplings)
    result = constraint.evaluate(point)
    audited = _audited_d0_with_running(couplings)
    yaml_delta_m_d = constraint.anchor.delta_m_d_gev.value
    expected_budget = yaml_delta_m_d / 2.0
    expected_ratio = audited.abs_m12_np / expected_budget

    assert isinstance(audited.abs_m12_np, float)
    assert isinstance(audited.ratio_to_budget, float)
    assert isinstance(audited.budget, float)
    assert math.isfinite(audited.abs_m12_np)
    assert math.isfinite(audited.ratio_to_budget)
    assert math.isfinite(audited.budget)
    assert result.passes is expected_pass
    assert result.predicted == pytest.approx(audited.abs_m12_np)
    # The direct core object carries quarkConstraints' older built-in D0
    # budget.  C001 intentionally overrides that with the YAML anchor, so the
    # apples-to-apples ratio oracle is direct |M12^NP| divided by the same YAML
    # budget used by the constraint.
    assert result.ratio == pytest.approx(expected_ratio)
    assert result.budget == pytest.approx(expected_budget)
    assert result.budget == pytest.approx(3.281e-15)
    if expected_pass:
        assert result.ratio <= 1.0
    else:
        assert result.ratio > 10.0


def test_evaluate_is_pure_and_deterministic():
    couplings = _uc_couplings(
        left=1.0e-4 + 0.5e-4j,
        right=0.7e-4 + 0.2e-4j,
    )
    before_left_up = couplings.left_up.copy()
    point = point_builder.build_from_quark_couplings(couplings)
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    np.testing.assert_array_equal(couplings.left_up, before_left_up)
