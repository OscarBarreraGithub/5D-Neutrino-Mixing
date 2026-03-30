"""Tests for quark mass-basis KK-gluon couplings."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from quarkConstraints.benchmarks import default_quark_targets, default_spurion_seed
from quarkConstraints.couplings import compute_quark_kk_gluon_couplings
from quarkConstraints.fit import fit_quark_sector


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
    couplings = compute_quark_kk_gluon_couplings(result)

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

    predicted_left_up = result.ckm_matrix @ couplings.left_overlap @ result.ckm_matrix.conjugate().T
    assert np.allclose(couplings.left_up / couplings.g_s, predicted_left_up, atol=1e-10)


def test_explicit_mkk_override_changes_scale_and_running():
    result = _fit_solution(0.25).result
    default_couplings = compute_quark_kk_gluon_couplings(result)
    overridden = compute_quark_kk_gluon_couplings(result, M_KK=6000.0, xi_KK=7.0)

    assert np.isclose(default_couplings.M_KK, result.point.Lambda_IR)
    assert np.isclose(overridden.M_KK, 6000.0)
    assert np.isclose(overridden.xi_KK, 7.0)
    assert overridden.alpha_s < default_couplings.alpha_s
    assert overridden.g_s < default_couplings.g_s


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

    small_couplings = compute_quark_kk_gluon_couplings(small_r)
    large_couplings = compute_quark_kk_gluon_couplings(large_r)

    assert _offdiag_norm(small_couplings.left_down) < _offdiag_norm(large_couplings.left_down)
    assert abs(small_couplings.left_down[1, 2]) < abs(large_couplings.left_down[1, 2])
