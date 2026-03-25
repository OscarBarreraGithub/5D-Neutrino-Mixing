import numpy as np
from scipy.optimize import brentq
from scipy.special import jv

from solvers.bessel import solve_kk
from warpConfig.baseParams import get_warp_params


def _jv_roots(nu: float, n_roots: int) -> np.ndarray:
    """Find the first positive zeros of J_nu(x) for test validation."""
    roots = []
    step = 0.05
    x_lo = 1.0e-6
    f_lo = jv(nu, x_lo)
    x_hi = x_lo

    while len(roots) < n_roots and x_hi < 50.0:
        x_hi += step
        f_hi = jv(nu, x_hi)
        if f_lo == 0.0 or f_lo * f_hi < 0.0:
            roots.append(brentq(lambda x: jv(nu, x), x_lo, x_hi))
        x_lo = x_hi
        f_lo = f_hi

    if len(roots) != n_roots:
        raise RuntimeError(f"Failed to bracket {n_roots} roots for J_{nu}")
    return np.array(roots)


def test_fermion_pp_negative_c_uses_alpha_plus_one_order():
    """For c < -1/2, the (++) branch should use nu = alpha + 1."""
    c = -0.75
    alpha = abs(c + 0.5)
    expected_nu = alpha + 1.0
    geometry = get_warp_params(Lambda_IR=1.0)

    masses, extras = solve_kk(
        species="fermion",
        bc="++",
        geometry=geometry,
        c=c,
        n_roots=3,
        exact=False,
    )

    expected_x = _jv_roots(expected_nu, 3)

    assert np.isclose(extras["nu"], expected_nu)
    assert np.allclose(extras["x"], expected_x, rtol=1e-10, atol=1e-12)
    assert np.allclose(masses, expected_x, rtol=1e-10, atol=1e-12)


def test_fermion_pp_positive_c_keeps_alpha_minus_one_order():
    """For c > -1/2, the historical nu = alpha - 1 branch remains unchanged."""
    c = 0.40
    alpha = abs(c + 0.5)
    expected_nu = alpha - 1.0
    geometry = get_warp_params(Lambda_IR=1.0)

    _, extras = solve_kk(
        species="fermion",
        bc="++",
        geometry=geometry,
        c=c,
        n_roots=3,
        exact=False,
    )

    expected_x = _jv_roots(expected_nu, 3)

    assert np.isclose(extras["nu"], expected_nu)
    assert np.allclose(extras["x"], expected_x, rtol=1e-10, atol=1e-12)
