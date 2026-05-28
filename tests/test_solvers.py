import numpy as np
import pytest
from scipy.special import jn_zeros

from solvers.bessel import solve_kk
from warpConfig.baseParams import get_warp_params


def _geometry(lambda_ir=1.0, epsilon=1.0e-15):
    return get_warp_params(k=lambda_ir / epsilon, Lambda_IR=lambda_ir)


def test_fermion_c_half_pp_ironly_matches_j0_roots():
    geometry = _geometry(lambda_ir=1.0)

    masses, extras = solve_kk(
        species="fermion",
        bc="++",
        geometry=geometry,
        c=0.5,
        n_roots=3,
        exact=False,
    )

    expected_roots = jn_zeros(0, 3)
    assert extras["nu"] == pytest.approx(0.0)
    assert np.allclose(extras["x"], expected_roots, rtol=1e-12, atol=1e-12)
    assert np.allclose(masses, expected_roots, rtol=1e-12, atol=1e-12)


def test_fermion_boundary_conditions_select_different_roots():
    geometry = _geometry(lambda_ir=1.0)

    pp_masses, pp_extras = solve_kk(
        species="fermion",
        bc="++",
        geometry=geometry,
        c=0.5,
        n_roots=3,
        exact=False,
    )
    mm_masses, mm_extras = solve_kk(
        species="fermion",
        bc="--",
        geometry=geometry,
        c=0.5,
        n_roots=3,
        exact=False,
    )

    assert pp_extras["nu"] == pytest.approx(0.0)
    assert mm_extras["nu"] == pytest.approx(1.0)
    assert np.allclose(pp_masses, jn_zeros(0, 3), rtol=1e-12, atol=1e-12)
    assert np.allclose(mm_masses, jn_zeros(1, 3), rtol=1e-12, atol=1e-12)
    assert mm_masses[0] > pp_masses[0]
    assert not np.allclose(pp_masses, mm_masses)


def test_gauge_nn_exact_first_kk_mass_is_standard_rs_root():
    lambda_ir = 3000.0
    geometry = _geometry(lambda_ir=lambda_ir)

    masses, extras = solve_kk(
        species="gauge",
        bc="NN",
        geometry=geometry,
        n_roots=1,
        exact=True,
    )

    assert extras["nu"] == pytest.approx(0.0)
    assert extras["x"][0] == pytest.approx(2.45, rel=1e-3)
    assert masses[0] == pytest.approx(2.45 * lambda_ir, rel=1e-3)
