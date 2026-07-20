"""Tests for quark mass-basis KK-gluon couplings."""

import sys
from pathlib import Path
from math import log, sqrt

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from quarkConstraints.benchmarks import default_quark_targets, default_spurion_seed
from quarkConstraints.couplings import (
    COUPLING_POLICY_FIXED_GSSTAR_3,
    COUPLING_POLICY_PERTURBATIVE_4D_LEGACY,
    COUPLING_POLICY_RS_VOLUME_SQRT2L_PHYSICAL,
    OPERATOR_CONVENTION_FIXED_GSSTAR_3,
    OPERATOR_CONVENTION_PERTURBATIVE_4D_LEGACY,
    OPERATOR_CONVENTION_RS_VOLUME_SQRT2L_PHYSICAL,
    compute_quark_kk_gluon_couplings,
)
from quarkConstraints.fit import fit_quark_sector
from quarkConstraints.scales import KKGluonMassConventionError


def _fit_solution(r_value: float):
    seed = default_spurion_seed()
    return fit_quark_sector(
        default_quark_targets(),
        r=r_value,
        seed=seed,
        overall_scale=seed.overall_scale,
        max_nfev=120,
    )


def _offdiag_norm(matrix: np.ndarray) -> float:
    arr = np.asarray(matrix, dtype=np.complex128).copy()
    np.fill_diagonal(arr, 0.0)
    return float(np.linalg.norm(arr, ord="fro"))


def test_mass_basis_couplings_have_expected_shapes_and_are_finite():
    result = _fit_solution(0.25).result
    # perturbative g_s (legacy repo_v1 behavior)
    couplings = compute_quark_kk_gluon_couplings(result, g_s_star=None)

    for matrix in (
        couplings.left_up,
        couplings.left_down,
        couplings.right_up,
        couplings.right_down,
        couplings.left_overlap,
        couplings.right_up_overlap,
        couplings.right_down_overlap,
    ):
        assert matrix.shape == (3, 3)
        assert np.all(np.isfinite(matrix.real))
        assert np.all(np.isfinite(matrix.imag))
        assert np.allclose(matrix, matrix.conjugate().T, atol=1e-12)

    np.testing.assert_allclose(couplings.ckm_matrix, result.ckm_matrix)
    assert couplings.ckm_source == "QuarkFitResult.ckm_matrix (U_L_u^dagger U_L_d)"
    assert couplings.ckm_matrix.flags.writeable is False

    predicted_left_up = result.ckm_matrix @ couplings.left_overlap @ result.ckm_matrix.conjugate().T
    assert np.allclose(couplings.left_up / couplings.g_s, predicted_left_up, atol=1e-10)


def test_explicit_mkk_override_changes_scale_and_running():
    result = _fit_solution(0.25).result
    # perturbative g_s (legacy repo_v1 behavior)
    default_couplings = compute_quark_kk_gluon_couplings(result, g_s_star=None)
    overridden = compute_quark_kk_gluon_couplings(
        result,
        m_kk_physical_gev=6000.0,
        lambda_ir_gev=result.point.Lambda_IR,
        xi_KK=2.0,
        g_s_star=None,
    )

    assert np.isclose(default_couplings.M_KK, result.point.Lambda_IR)
    assert np.isclose(overridden.M_KK, 6000.0)
    assert np.isclose(overridden.xi_KK, 2.0)
    assert np.isclose(overridden.m_kk_physical_gev, 6000.0)
    assert np.isclose(overridden.lambda_ir_gev, result.point.Lambda_IR)
    assert overridden.alpha_s < default_couplings.alpha_s
    assert overridden.g_s < default_couplings.g_s


def test_kk_gluon_coupling_policies_have_distinct_operator_conventions():
    result = _fit_solution(0.25).result
    legacy = compute_quark_kk_gluon_couplings(result, M_KK=6000.0, g_s_star=None)
    fixed = compute_quark_kk_gluon_couplings(result, M_KK=6000.0, g_s_star=3.0)
    physical = compute_quark_kk_gluon_couplings(
        result,
        M_KK=6000.0,
        coupling_policy_id=COUPLING_POLICY_RS_VOLUME_SQRT2L_PHYSICAL,
    )

    assert legacy.coupling_policy_id == COUPLING_POLICY_PERTURBATIVE_4D_LEGACY
    assert fixed.coupling_policy_id == COUPLING_POLICY_FIXED_GSSTAR_3
    assert physical.coupling_policy_id == COUPLING_POLICY_RS_VOLUME_SQRT2L_PHYSICAL
    assert legacy.operator_convention_id == OPERATOR_CONVENTION_PERTURBATIVE_4D_LEGACY
    assert fixed.operator_convention_id == OPERATOR_CONVENTION_FIXED_GSSTAR_3
    assert physical.operator_convention_id == OPERATOR_CONVENTION_RS_VOLUME_SQRT2L_PHYSICAL
    assert len(
        {
            legacy.operator_convention_id,
            fixed.operator_convention_id,
            physical.operator_convention_id,
        }
    ) == 3
    assert legacy.g_s_multiplier == pytest.approx(1.0)
    assert fixed.g_s_multiplier == pytest.approx(3.0 / fixed.g_s_4d)
    expected_volume = sqrt(2.0 * log(result.point.k / result.point.Lambda_IR))
    assert physical.g_s_multiplier == pytest.approx(expected_volume)
    assert legacy.g_eff == pytest.approx(legacy.g_s_4d)
    assert fixed.g_eff == pytest.approx(3.0)
    assert physical.g_eff == pytest.approx(physical.g_s_4d * expected_volume)


def test_explicit_perturbative_policy_rejects_gsstar_override():
    result = _fit_solution(0.25).result

    with pytest.raises(ValueError, match="perturbative_4d_legacy requires g_s_star=None"):
        compute_quark_kk_gluon_couplings(
            result,
            M_KK=6000.0,
            g_s_star=3.0,
            coupling_policy_id=COUPLING_POLICY_PERTURBATIVE_4D_LEGACY,
        )


def test_mismatched_physical_mass_and_lambda_ir_fails_loudly():
    result = _fit_solution(0.25).result

    with pytest.raises(KKGluonMassConventionError, match="mass convention mismatch"):
        compute_quark_kk_gluon_couplings(
            result,
            m_kk_physical_gev=6000.0,
            lambda_ir_gev=result.point.Lambda_IR,
            xi_KK=1.0,
            g_s_star=None,
        )


def test_small_r_reduces_down_sector_left_handed_offdiagonals():
    targets = default_quark_targets()
    seed = default_spurion_seed()
    small_r = fit_quark_sector(
        targets,
        r=0.05,
        seed=seed,
        overall_scale=seed.overall_scale,
        max_nfev=120,
    ).result
    large_r = fit_quark_sector(
        targets,
        r=1.0,
        seed=seed,
        overall_scale=seed.overall_scale,
        max_nfev=120,
    ).result

    # perturbative g_s (legacy repo_v1 behavior)
    small_couplings = compute_quark_kk_gluon_couplings(small_r, g_s_star=None)
    large_couplings = compute_quark_kk_gluon_couplings(large_r, g_s_star=None)

    assert _offdiag_norm(small_couplings.left_down) < _offdiag_norm(large_couplings.left_down)
    assert abs(small_couplings.left_down[1, 2]) < abs(large_couplings.left_down[1, 2])
