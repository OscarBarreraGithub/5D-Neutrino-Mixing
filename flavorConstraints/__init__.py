"""flavorConstraints package — lepton flavor violation bounds."""

from .muToEGamma import (
    BR_LIMIT_PAPER,
    C_PAPER,
    GAUGE_KK_ROOT_NN,
    PEREZ_RANDALL_LFV_M_KK_CONVENTION,
    PEREZ_RANDALL_LFV_XI_KK,
    PREFAC_BR,
    assert_perez_randall_lfv_m_kk_convention,
    check_mu_to_e_gamma,
    check_mu_to_e_gamma_raw,
    coefficient_from_br_limit,
    default_m_kk_from_lambda_ir,
    perez_randall_lfv_m_kk_from_lambda_ir,
)

__all__ = [
    'check_mu_to_e_gamma',
    'check_mu_to_e_gamma_raw',
    'coefficient_from_br_limit',
    'default_m_kk_from_lambda_ir',
    'perez_randall_lfv_m_kk_from_lambda_ir',
    'assert_perez_randall_lfv_m_kk_convention',
    'PREFAC_BR',
    'BR_LIMIT_PAPER',
    'C_PAPER',
    'GAUGE_KK_ROOT_NN',
    'PEREZ_RANDALL_LFV_M_KK_CONVENTION',
    'PEREZ_RANDALL_LFV_XI_KK',
]
