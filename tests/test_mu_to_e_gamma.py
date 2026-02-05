import numpy as np

from flavorConstraints.muToEGamma import (
    C_PAPER,
    check_mu_to_e_gamma,
    check_mu_to_e_gamma_raw,
    coefficient_from_br_limit,
)
from neutrinos.neutrinoValues import pmns_matrix
from yukawa import compute_all_yukawas


# PDG-style mixing values (as used in neutrinoValues.py)
SIN2_THETA12 = 0.307
SIN2_THETA23 = 0.534
SIN2_THETA13 = 0.0216
DELTA_CP = 1.21 * np.pi  # normal ordering value in neutrinoValues.py


def _pmns_from_sin2():
    theta12 = np.arcsin(np.sqrt(SIN2_THETA12))
    theta23 = np.arcsin(np.sqrt(SIN2_THETA23))
    theta13 = np.arcsin(np.sqrt(SIN2_THETA13))
    return pmns_matrix(theta12, theta23, theta13, DELTA_CP, 0.0, 0.0)


def test_paper_example_passes_paper_bound():
    """Paper example (Eq. 3.10) should satisfy the paper bound (Eq. 4.14)."""
    # Perez-Randall example: k * Y_N ~= (0.02, 0.03, 0.07)
    y_n_k = np.array([0.02, 0.03, 0.07])
    y_n_bar = 2.0 * y_n_k  # Ȳ_N = 2k Y_N
    U = _pmns_from_sin2()

    res = check_mu_to_e_gamma_raw(y_n_bar, U, M_KK=3000.0, C=C_PAPER)
    assert res["passes"] is True
    assert np.isclose(res["lhs"], 0.00140753, rtol=0.1)


def test_updated_megii_limit_is_stronger():
    """Published MEG II 2025 limit tightens C relative to paper value."""
    br_limit_megii_2025 = 1.5e-13
    c_updated = coefficient_from_br_limit(br_limit_megii_2025)

    assert c_updated < C_PAPER
    assert np.isclose(c_updated, 0.00193649, rtol=0.02)

    # Benchmark Yukawas from compute_all_yukawas() output (Perez-Randall point)
    y_n_bar = np.array([0.20416916, 0.43081761, 1.02432944])
    U = _pmns_from_sin2()

    res = check_mu_to_e_gamma_raw(y_n_bar, U, M_KK=3000.0, C=c_updated)
    assert res["passes"] is False
    assert np.isclose(res["lhs"], 0.07234407, rtol=0.02)


def test_m_kk_override_changes_lfv_rhs_scaling():
    """check_mu_to_e_gamma should respect explicit M_KK override."""
    result = compute_all_yukawas(
        Lambda_IR=3000.0,
        c_L=0.58,
        c_E=[0.75, 0.60, 0.50],
        c_N=0.27,
        M_N=1.22e18,
        lightest_nu_mass=0.002,
        ordering="normal",
    )

    no_override = check_mu_to_e_gamma(result, C=C_PAPER, reference_scale=3000.0)
    override = check_mu_to_e_gamma(
        result,
        C=C_PAPER,
        reference_scale=3000.0,
        M_KK_override=6000.0,
    )

    # RHS scales as (M_KK / reference_scale)^2
    assert np.isclose(override["rhs"], 4.0 * no_override["rhs"])
    assert np.isclose(override["lhs"], no_override["lhs"])
