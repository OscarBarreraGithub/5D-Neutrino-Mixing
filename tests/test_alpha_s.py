"""Tests for the qcd module.

Validates alpha_s running against PDG tabulated values, RunDec cross-checks,
and known analytic results.
"""

import numpy as np
import pytest

from qcd import alpha_s, alpha_s_array, match_alpha_s, ALPHA_S_MZ, M_Z, M_TOP_MS
from qcd.beta_function import beta_0, beta_1, beta_2, beta_coefficients


# ---------- Beta function coefficients ----------


def test_beta_0_nf5():
    """beta_0(5) = 11 - 10/3 = 23/3."""
    assert np.isclose(beta_0(5), 23.0 / 3.0)


def test_beta_0_nf6():
    """beta_0(6) = 11 - 4 = 7."""
    assert np.isclose(beta_0(6), 7.0)


def test_beta_0_pure_ym():
    """Pure Yang–Mills (n_f=0): beta_0 = 11."""
    assert np.isclose(beta_0(0), 11.0)


def test_beta_1_nf5():
    """beta_1(5) = 102 - 190/3 = 116/3."""
    assert np.isclose(beta_1(5), 116.0 / 3.0)


def test_beta_coefficients_shape():
    """beta_coefficients returns correct length array."""
    for nl in [1, 2, 3, 4]:
        c = beta_coefficients(5, nl)
        assert c.shape == (nl,)


def test_beta_coefficients_invalid():
    """Invalid n_f or n_loops raises ValueError."""
    with pytest.raises(ValueError):
        beta_coefficients(7, 3)
    with pytest.raises(ValueError):
        beta_coefficients(5, 5)


# ---------- alpha_s running ----------


def test_identity():
    """alpha_s at M_Z should return the input value."""
    assert np.isclose(alpha_s(M_Z), ALPHA_S_MZ, rtol=1e-8)


def test_asymptotic_freedom():
    """alpha_s decreases with increasing energy."""
    a_mz = alpha_s(M_Z)
    a_1tev = alpha_s(1000.0)
    a_10tev = alpha_s(10000.0)
    assert a_1tev < a_mz
    assert a_10tev < a_1tev


def test_mt_4loop():
    """alpha_s(m_t) ~ 0.1084 at 4-loop (MS-bar threshold)."""
    assert np.isclose(alpha_s(M_TOP_MS, n_loops=4), 0.1084, atol=0.002)


def test_1tev_4loop():
    """alpha_s(1 TeV) ~ 0.0885 at 4-loop (RunDec-style)."""
    assert np.isclose(alpha_s(1000.0, n_loops=4), 0.0885, atol=0.002)


def test_3tev_4loop():
    """alpha_s(3 TeV) ~ 0.0797 at 4-loop (RunDec-style)."""
    assert np.isclose(alpha_s(3000.0, n_loops=4), 0.0797, atol=0.002)


def test_10tev_4loop():
    """alpha_s(10 TeV) ~ 0.0718 at 4-loop with 3-loop matching."""
    assert np.isclose(alpha_s(10000.0, n_loops=4), 0.0718, atol=0.002)


def test_cms_896gev():
    """CMS/RunDec example: alpha_s(896 GeV) ≈ 0.0889."""
    assert np.isclose(alpha_s(896.0, n_loops=4), 0.0889, atol=0.0035)


def test_loop_order_convergence():
    """Higher loop orders give similar results (< 2% difference at 1 TeV)."""
    a1 = alpha_s(1000.0, n_loops=1)
    a2 = alpha_s(1000.0, n_loops=2)
    a3 = alpha_s(1000.0, n_loops=3)
    a4 = alpha_s(1000.0, n_loops=4)
    assert abs(a1 - a4) / a4 < 0.02
    assert abs(a2 - a4) / a4 < 0.005
    assert abs(a3 - a4) / a4 < 0.001


def test_1loop_analytic():
    """1-loop numerical result matches the exact analytic formula.

    Below the top threshold (n_f=5 throughout):
        alpha_s(mu) = alpha_s(M_Z) / (1 + beta_0*alpha_s(M_Z)/(4pi)*ln(mu^2/M_Z^2))
    """
    mu = 150.0  # below m_t, so n_f=5 the whole way
    b0 = beta_0(5)
    expected = ALPHA_S_MZ / (
        1.0 + b0 * ALPHA_S_MZ / (4.0 * np.pi) * np.log(mu**2 / M_Z**2)
    )
    assert np.isclose(alpha_s(mu, n_loops=1), expected, rtol=1e-6)


# ---------- Threshold handling ----------


def test_threshold_matching():
    """alpha_s satisfies the matching relation across the top threshold."""
    mt = M_TOP_MS
    a_below = alpha_s(mt - 0.01, n_loops=4)
    a_above = alpha_s(mt + 0.01, n_loops=4)
    a_below_to_6 = match_alpha_s(a_below, n_f_from=5, n_f_to=6, matching_loops=3)
    assert np.isclose(a_below_to_6, a_above, rtol=1e-4)


def test_no_thresholds():
    """Disabling thresholds gives a finite positive result."""
    result = alpha_s(150.0, n_loops=2, thresholds=[])
    assert 0.05 < result < 0.2


# ---------- Edge cases ----------


def test_negative_mu_raises():
    with pytest.raises(ValueError):
        alpha_s(-100.0)


def test_zero_mu_raises():
    with pytest.raises(ValueError):
        alpha_s(0.0)


def test_invalid_loop_order():
    with pytest.raises(ValueError):
        alpha_s(1000.0, n_loops=0)
    with pytest.raises(ValueError):
        alpha_s(1000.0, n_loops=5)


def test_alpha_s_array():
    """alpha_s_array returns correct shape and monotonic values."""
    mus = np.array([100.0, 1000.0, 10000.0])
    result = alpha_s_array(mus)
    assert result.shape == (3,)
    assert np.all(np.diff(result) < 0)


def test_matching_round_trip():
    """Matching down then up returns the original value (to matching order)."""
    a = alpha_s(1000.0, n_loops=4)
    a4 = match_alpha_s(a, n_f_from=6, n_f_to=5, matching_loops=3)
    a6 = match_alpha_s(a4, n_f_from=5, n_f_to=6, matching_loops=3)
    assert np.isclose(a, a6, rtol=1e-6)
