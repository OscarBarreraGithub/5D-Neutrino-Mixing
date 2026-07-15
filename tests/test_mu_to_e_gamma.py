import numpy as np
import pytest

from flavorConstraints.muToEGamma import (
    C_PAPER,
    GAUGE_KK_ROOT_NN,
    PEREZ_RANDALL_LFV_M_KK_CONVENTION,
    PEREZ_RANDALL_LFV_XI_KK,
    assert_perez_randall_lfv_m_kk_convention,
    check_mu_to_e_gamma,
    check_mu_to_e_gamma_raw,
    coefficient_from_br_limit,
    default_m_kk_from_lambda_ir,
    perez_randall_lfv_m_kk_from_lambda_ir,
)
from neutrinos.neutrinoValues import pmns_matrix
from yukawa import compute_all_yukawas

# NuFIT 6.1 (2025), IC24 with SK atmospheric data, normal ordering.
SIN2_THETA12 = 0.3088
SIN2_THETA23 = 0.470
SIN2_THETA13 = 0.02248
DELTA_CP = np.deg2rad(212.0)


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
    assert np.isclose(res["lhs"], 0.00125309, rtol=0.1)


def test_updated_megii_limit_is_stronger():
    """Published MEG II 2025 limit tightens C relative to paper value."""
    br_limit_megii_2025 = 1.5e-13
    c_updated = coefficient_from_br_limit(br_limit_megii_2025)

    assert c_updated < C_PAPER
    assert np.isclose(c_updated, 0.00193649, rtol=0.02)

    # Benchmark Yukawas from compute_all_yukawas() output for the repo-local
    # paper-inspired validation point.
    y_n_bar = np.array([0.20416916, 0.43091265, 1.02237364])
    U = _pmns_from_sin2()

    res = check_mu_to_e_gamma_raw(y_n_bar, U, M_KK=3000.0, C=c_updated)
    assert res["passes"] is False
    assert np.isclose(res["lhs"], 0.06226618, rtol=0.02)


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


def test_default_lfv_check_uses_internal_lambda_ir_convention():
    """Default LFV check should not silently switch to the gauge KK root."""
    result = compute_all_yukawas(
        Lambda_IR=3000.0,
        c_L=0.58,
        c_E=[0.75, 0.60, 0.50],
        c_N=0.27,
        M_N=1.22e18,
        lightest_nu_mass=0.002,
        ordering="normal",
    )

    explicit_gauge_mkk = default_m_kk_from_lambda_ir(3000.0)
    assert np.isclose(explicit_gauge_mkk, GAUGE_KK_ROOT_NN * 3000.0)
    assert PEREZ_RANDALL_LFV_M_KK_CONVENTION == "perez_randall_geometric_lambda_ir_v1"
    assert np.isclose(PEREZ_RANDALL_LFV_XI_KK, 1.0)
    assert np.isclose(perez_randall_lfv_m_kk_from_lambda_ir(3000.0), 3000.0)
    assert np.isclose(
        assert_perez_randall_lfv_m_kk_convention(
            m_kk_gev=3000.0,
            Lambda_IR=3000.0,
        ),
        3000.0,
    )

    lfv = check_mu_to_e_gamma(result, C=C_PAPER, reference_scale=3000.0)
    assert np.isclose(lfv["rhs"], C_PAPER)


def test_perez_randall_lfv_mkk_assertion_rejects_physical_gauge_mass():
    """M-14 guard: do not feed the first-gauge Bessel mass to the calibrated LFV prefactor."""
    physical_gauge_mkk = default_m_kk_from_lambda_ir(3000.0)

    with pytest.raises(AssertionError, match=PEREZ_RANDALL_LFV_M_KK_CONVENTION):
        assert_perez_randall_lfv_m_kk_convention(
            m_kk_gev=physical_gauge_mkk,
            Lambda_IR=3000.0,
        )


def test_repo_benchmark_raw_core_values_and_br_oracle_are_distinct():
    result = compute_all_yukawas(
        Lambda_IR=3000.0,
        c_L=0.58,
        c_E=[0.75, 0.60, 0.50],
        c_N=0.27,
        M_N=1.22e18,
        lightest_nu_mass=0.002,
        ordering="normal",
    )
    U = _pmns_from_sin2()

    core = check_mu_to_e_gamma_raw(
        result.Y_N_bar,
        U,
        M_KK=3000.0,
        C=C_PAPER,
        reference_scale=3000.0,
    )
    c_consistent = coefficient_from_br_limit(1.5e-13)
    core_consistent = check_mu_to_e_gamma_raw(
        result.Y_N_bar,
        U,
        M_KK=3000.0,
        C=c_consistent,
        reference_scale=3000.0,
    )
    br_np = 4.0e-8 * float(core["lhs"]) ** 2

    np.testing.assert_allclose(
        result.Y_N_bar,
        [0.20416916, 0.43091265, 1.02237364],
        rtol=1.0e-8,
        atol=1.0e-10,
    )
    assert core["off_diagonal_12"] == pytest.approx(
        -0.034773700046830024 + 0.05165132075170267j
    )
    assert core["lhs"] == pytest.approx(0.062266115587389717)
    assert core["rhs"] == pytest.approx(0.02)
    assert core["ratio"] == pytest.approx(3.1133057793694858)
    assert core_consistent["ratio"] == pytest.approx(32.154083827064859)
    assert br_np == pytest.approx(1.5508276601368708e-10)
    assert br_np / 1.5e-13 == pytest.approx(1033.8851067579139)
    assert br_np > 1.5e-13
