import numpy as np
import pytest
from scipy.optimize import brentq
from scipy.special import jv

from solvers.bessel import _F_exact, _nu_for, solve_kk
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


def _exact_quantization_roots(nu: float, eps: float, n_roots: int) -> np.ndarray:
    """Brute-force sign scan of the solver's exact quantization function."""
    roots = []
    F = _F_exact(nu, eps)
    step = 0.02
    x_lo = 1.0e-5
    f_lo = F(x_lo)
    x_hi = x_lo

    while len(roots) < n_roots and x_hi < 80.0:
        x_hi += step
        f_hi = F(x_hi)
        if np.isfinite(f_lo) and np.isfinite(f_hi) and np.sign(f_lo) != np.sign(f_hi):
            roots.append(brentq(F, x_lo, x_hi, xtol=1e-13, rtol=1e-13))
        x_lo = x_hi
        f_lo = f_hi

    if len(roots) != n_roots:
        raise RuntimeError(f"Failed to bracket {n_roots} exact roots for nu={nu}")
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


def test_deep_exact_gauge_nn_roots_are_sorted_and_do_not_skip():
    """C-4: exact gauge NN tower should include the sixth root near 18.118."""
    geometry = get_warp_params(Lambda_IR=3000.0)
    n_roots = 8
    _, extras = solve_kk(
        species="gauge",
        bc="NN",
        geometry=geometry,
        n_roots=n_roots,
        exact=True,
    )
    expected = _exact_quantization_roots(extras["nu"], geometry["epsilon"], n_roots)

    assert np.all(np.diff(extras["x"]) > 0.0)
    assert np.allclose(extras["x"], expected, rtol=1e-10, atol=1e-10)
    assert extras["x"][5] == pytest.approx(18.118, abs=0.002)


def test_deep_exact_fermion_pp_roots_are_sorted_and_do_not_skip():
    """C-4: fermion ++ c=0.27 fifth root is ~15.29, not the old 52.98 outlier."""
    geometry = get_warp_params(Lambda_IR=3000.0)
    c = 0.27
    n_roots = 8
    _, extras = solve_kk(
        species="fermion",
        bc="++",
        geometry=geometry,
        c=c,
        n_roots=n_roots,
        exact=True,
    )
    nu, _ = _nu_for("fermion", "++", c)
    expected = _exact_quantization_roots(nu, geometry["epsilon"], n_roots)

    assert np.all(np.diff(extras["x"]) > 0.0)
    assert np.allclose(extras["x"], expected, rtol=1e-10, atol=1e-10)
    assert extras["x"][4] == pytest.approx(15.29, abs=0.002)
