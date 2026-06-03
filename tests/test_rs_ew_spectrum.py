import math

import numpy as np
import pytest
from numpy.polynomial.legendre import leggauss
from scipy.optimize import brentq
from scipy.special import j0, j1, y0, y1

from quarkConstraints.rs_ew_spectrum import (
    RSEWOverlapConvergenceError,
    RSEWSpectrum,
    g0_squared,
    kk_ew_mass,
)
from solvers.bessel import solve_kk
from warpConfig.baseParams import get_warp_params
from warpConfig.wavefuncs import f_IR, f_UV


LAMBDA_IR_GEV = 3000.0
EPSILON_RS = 1.0e-15
K_GEV = LAMBDA_IR_GEV / EPSILON_RS
QUADRATURE_ORDER = 1024


def _independent_gauge_nn_roots(epsilon, n_roots):
    def cross_product(x):
        return float(j0(x) * y0(epsilon * x) - j0(epsilon * x) * y0(x))

    roots = []
    x_prev = 1.0e-10
    f_prev = cross_product(x_prev)
    x_limit = (n_roots + 2.0) * math.pi
    step = math.pi / 32.0

    while len(roots) < n_roots and x_prev < x_limit:
        x_next = min(x_prev + step, x_limit)
        f_next = cross_product(x_next)
        if np.sign(f_prev) != np.sign(f_next):
            root = brentq(cross_product, x_prev, x_next, xtol=1.0e-13, rtol=1.0e-13)
            if not roots or root > roots[-1] + 1.0e-8:
                roots.append(float(root))
        x_prev, f_prev = x_next, f_next

    if len(roots) != n_roots:
        raise AssertionError(f"found {len(roots)} independent roots, expected {n_roots}")
    return np.asarray(roots, dtype=float)


def _independent_a_partials(c, roots, epsilon, quadrature_order):
    nodes, weights = leggauss(int(quadrature_order))
    y = 0.5 * (nodes + 1.0)
    y_weights = 0.5 * weights
    log_eps = math.log(float(epsilon))
    warp_log = -log_eps
    t = np.exp(log_eps * (1.0 - y))

    bvals = -j0(roots) / y0(roots)
    u = roots[:, None] * t[None, :]
    raw = t[None, :] * (j1(u) + bvals[:, None] * y1(u))
    norm_sq = warp_log * np.sum(y_weights[None, :] * raw * raw, axis=1)
    profile_norms = 1.0 / np.sqrt(norm_sq)
    chi = math.sqrt(warp_log) * profile_norms[:, None] * raw
    chi_ir = math.sqrt(warp_log) * profile_norms * (j1(roots) + bvals * y1(roots))

    density_y = warp_log * np.asarray(g0_squared(c, t, epsilon), dtype=float)
    omega = np.sum(y_weights[None, :] * density_y[None, :] * chi, axis=1)
    return np.cumsum((roots[0] / roots) ** 2 * chi_ir * omega)


@pytest.fixture(scope="module")
def spectrum():
    return RSEWSpectrum.build(
        lambda_ir_gev=LAMBDA_IR_GEV,
        k_gev=K_GEV,
        n_gauge_modes=512,
        quadrature_order=QUADRATURE_ORDER,
    )


def test_exact_gauge_nn_root_and_kk_ew_mass_use_physical_root(spectrum):
    geometry = get_warp_params(k=K_GEV, Lambda_IR=LAMBDA_IR_GEV)
    masses, extras = solve_kk("gauge", "NN", geometry, n_roots=1, exact=True)
    direct_root = float(extras["x"][0])

    assert direct_root == pytest.approx(2.450509663813736, rel=0.0, abs=1.0e-12)
    assert 2.449 < direct_root < 2.451
    assert abs(direct_root - 2.404825557695773) > 4.0e-2
    assert spectrum.gauge_roots_x[0] == pytest.approx(direct_root, rel=0.0, abs=1.0e-12)
    assert masses[0] == pytest.approx(direct_root * LAMBDA_IR_GEV, rel=0.0, abs=1.0e-9)
    assert spectrum.kk_ew_mass_gev == pytest.approx(7351.528991441208, rel=0.0, abs=1.0e-9)
    assert kk_ew_mass(lambda_ir_gev=LAMBDA_IR_GEV, k_gev=K_GEV) == pytest.approx(
        spectrum.kk_ew_mass_gev,
        rel=0.0,
        abs=1.0e-9,
    )


def test_gauge_nn_tower_is_strictly_increasing_unique_and_not_skipped(spectrum):
    expected_first_six = np.asarray(
        [
            2.450509663813736,
            5.567547286728064,
            8.70195714458206,
            11.840260538783498,
            14.980018189299924,
            18.120466878759633,
        ],
        dtype=float,
    )

    roots = spectrum.gauge_roots_x
    assert np.all(np.diff(roots) > 1.0e-8)
    assert np.unique(np.round(roots, decimals=12)).size == roots.size
    assert roots[:6] == pytest.approx(expected_first_six, rel=0.0, abs=5.0e-12)
    assert roots[5] == pytest.approx(18.120466878759633, rel=0.0, abs=5.0e-12)
    assert roots[5] < 19.0
    assert np.max(np.abs(np.diff(roots[:8]) - math.pi)) < 3.0e-2


def test_zero_mode_endpoint_identities_match_existing_overlap_factors(spectrum):
    residuals = []
    for c in [-0.4, 0.0, 0.5, 0.7, 1.2]:
        ir_residual = abs(
            g0_squared(c, 1.0, spectrum.epsilon)
            - 2.0 * float(f_IR(c, spectrum.epsilon)) ** 2
        )
        uv_residual = abs(
            g0_squared(c, spectrum.epsilon, spectrum.epsilon)
            - 2.0 * float(f_UV(c, spectrum.epsilon)) ** 2
        )
        residuals.extend([ir_residual, uv_residual])

    assert max(residuals) < 3.0e-15


def test_a_convergence_uses_doubled_truncation_and_forced_failure_reaches_512(spectrum):
    c_value = 0.2
    independent_roots = _independent_gauge_nn_roots(EPSILON_RS, 512)
    manual_partials = _independent_a_partials(
        c_value,
        independent_roots,
        EPSILON_RS,
        QUADRATURE_ORDER,
    )

    a_16 = float(manual_partials[15])
    a_32 = float(manual_partials[31])
    manual_rel = abs(a_32 - a_16) / max(abs(a_32), 1.0e-12)
    result = spectrum.a_with_diagnostics(c_value)

    assert spectrum.gauge_roots_x[:512] == pytest.approx(
        independent_roots,
        rel=0.0,
        abs=5.0e-11,
    )
    assert a_16 == pytest.approx(21.859749049606535, rel=1.0e-12)
    assert a_32 == pytest.approx(21.860547338064457, rel=1.0e-12)
    assert manual_rel == pytest.approx(3.65173133854731e-05, rel=1.0e-12)
    assert result.raw_value == pytest.approx(a_32, rel=1.0e-12)
    assert result.modes == 32
    assert result.previous_modes == 16
    assert result.relative_delta < 1.0e-3

    with pytest.raises(RSEWOverlapConvergenceError, match="N=512"):
        spectrum.a_with_diagnostics(c_value, rel_tol=0.0)


def test_universal_c_subtraction_precursor_vanishes_for_species_set(spectrum):
    universal_c = 0.65
    a_ref = spectrum.a(universal_c)
    species_bulk_masses = {
        "c_Q": [universal_c, universal_c, universal_c],
        "c_u": [universal_c, universal_c, universal_c],
        "c_d": [universal_c, universal_c, universal_c],
        "c_L": [universal_c, universal_c, universal_c],
        "c_E": [universal_c, universal_c, universal_c],
    }
    subtracted = [
        spectrum.a(c, a_ref=a_ref)
        for values in species_bulk_masses.values()
        for c in values
    ]

    assert a_ref == pytest.approx(-1.4571399330265025, rel=1.0e-12)
    assert np.allclose(subtracted, 0.0, rtol=0.0, atol=0.0)


def test_ir_localized_b_right_precursor_has_larger_a_than_light_singlet(spectrum):
    a_light = spectrum.a(0.65)
    a_b_right = spectrum.a(0.2)

    assert a_light == pytest.approx(-1.4571399330265025, rel=1.0e-12)
    assert a_b_right == pytest.approx(21.860547338064457, rel=1.0e-12)
    assert a_b_right > a_light


def test_overlap_outputs_are_finite_and_cache_deterministic(spectrum):
    uncached = spectrum.a(0.2)
    cached = spectrum.a(0.2)
    omega_first = spectrum.omega(0.2, max_modes=64)
    omega_cached = spectrum.omega(0.2, max_modes=64)
    values = np.concatenate(([uncached, cached], omega_first, omega_cached))

    assert math.isfinite(uncached)
    assert cached == pytest.approx(uncached, rel=0.0, abs=0.0)
    assert np.all(np.isfinite(values))
    assert np.array_equal(omega_first, omega_cached)
