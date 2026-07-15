"""μ→eγ lepton flavor violation constraint.

Checks whether the off-diagonal element (Ȳ_N Ȳ_N†)₁₂ respects the dipole-
operator bound from μ→eγ.

Physics:
    In the charged-lepton mass basis the neutrino Yukawa matrix is

        Y_N = U_PMNS · diag(Y_{N_1}, Y_{N_2}, Y_{N_3})

    with rescaled (dimensionless) Yukawas Ȳ_N = 2k Y_N.  The product

        Ȳ_N Ȳ_N† = (2k)² U diag(Y²_{N_i}) U†

    has off-diagonal entries that mediate lepton-flavor-violating processes.
    The μ→eγ rate is controlled by the (1,2) element, constrained by

        |(Ȳ_N Ȳ_N†)₁₂| ≤ C × (M_KK / 3 TeV)²

    with C = 0.02 from naive dimensional analysis of the brane-localized
    dipole operator (Perez & Randall, arXiv:0805.4652, Eq. 4.14).

Reference:
    Perez & Randall, arXiv:0805.4652
"""

from typing import Dict, Optional, Union

import numpy as np

PREFAC_BR = 4.0e-8
# Paper-era experimental limit used in Perez & Randall (MEGA):
BR_LIMIT_PAPER = 1.2e-11
# Paper bound quoted for IR-brane Higgs, Eq. (4.14).
# Note: sqrt(BR_LIMIT_PAPER / PREFAC_BR) ≈ 0.0173; the paper rounds up to 0.02.
C_PAPER = 0.02
PEREZ_RANDALL_LFV_M_KK_CONVENTION = "perez_randall_geometric_lambda_ir_v1"
PEREZ_RANDALL_LFV_XI_KK = 1.0
# Optional utility for converting the geometric IR scale into the first gauge KK
# mass, m_{g^(1)} = x_1 * Lambda_IR, with x_1 from the (NN) Bessel equation.
# This is not the repo's default LFV convention.
GAUGE_KK_ROOT_NN = 2.448687135269161


def perez_randall_lfv_m_kk_from_lambda_ir(
    Lambda_IR: float,
    *,
    xi_KK: float = PEREZ_RANDALL_LFV_XI_KK,
) -> float:
    """Return the KK scale paired with the repo's Perez-Randall LFV prefactor.

    The calibrated repo default uses the geometric IR scale directly,
    ``M_KK = Lambda_IR``.  Non-unit ``xi_KK`` values are explicit systematic
    variations and must be paired with consistent metadata/reference-scale choices.
    """
    if Lambda_IR <= 0:
        raise ValueError("Lambda_IR must be positive")
    if xi_KK <= 0:
        raise ValueError("xi_KK must be positive")
    return float(xi_KK * Lambda_IR)


def assert_perez_randall_lfv_m_kk_convention(
    *,
    m_kk_gev: float,
    Lambda_IR: float,
    xi_KK: float = PEREZ_RANDALL_LFV_XI_KK,
    rtol: float = 1.0e-12,
    atol: float = 1.0e-9,
) -> float:
    """Assert that ``m_kk_gev`` uses the named Perez-Randall LFV convention."""
    expected = perez_randall_lfv_m_kk_from_lambda_ir(Lambda_IR, xi_KK=xi_KK)
    if not np.isclose(float(m_kk_gev), expected, rtol=rtol, atol=atol):
        raise AssertionError(
            f"{PEREZ_RANDALL_LFV_M_KK_CONVENTION} requires M_KK={expected:.17g} GeV "
            f"for Lambda_IR={float(Lambda_IR):.17g} GeV and xi_KK={float(xi_KK):.17g}; "
            f"got {float(m_kk_gev):.17g} GeV"
        )
    return float(m_kk_gev)


def default_m_kk_from_lambda_ir(Lambda_IR: float, xi_KK: float = GAUGE_KK_ROOT_NN) -> float:
    """Map the geometric IR scale to an explicit physical KK-mass convention.

    Parameters
    ----------
    Lambda_IR : float
        Geometric IR scale, Lambda_IR = 1 / z_v.
    xi_KK : float, optional
        Chosen physical-mass convention. Default is the first gauge KK root for
        Neumann-Neumann boundary conditions.
    """
    if Lambda_IR <= 0:
        raise ValueError("Lambda_IR must be positive")
    if xi_KK <= 0:
        raise ValueError("xi_KK must be positive")
    return float(xi_KK * Lambda_IR)


def coefficient_from_br_limit(br_limit: float, prefactor: float = PREFAC_BR) -> float:
    """Compute C from a μ→eγ branching ratio limit.

    Uses the NDA estimate:
        BR(μ→eγ) ≈ prefactor × |(Ȳ_N Ȳ_N†)₁₂|² × (3 TeV / M_KK)⁴

    so C = sqrt(BR_limit / prefactor).
    """
    if br_limit <= 0:
        raise ValueError("br_limit must be positive")
    if prefactor <= 0:
        raise ValueError("prefactor must be positive")
    return float(np.sqrt(br_limit / prefactor))

def check_mu_to_e_gamma(
    yukawa_result,
    C: float = C_PAPER,
    reference_scale: float = 3000.0,
    M_KK_override: Optional[float] = None,
) -> Dict[str, Union[float, bool, complex, np.ndarray]]:
    """Check the μ→eγ dipole constraint on the neutrino Yukawa.

    Parameters
    ----------
    yukawa_result : yukawa.compute_yukawas.YukawaResult
        Output of ``compute_all_yukawas()``.
    C : float, optional
        Numerical coefficient in the bound.  Default 0.02 (Perez–Randall).
    reference_scale : float, optional
        Reference KK scale in GeV (denominator).  Default 3000 (= 3 TeV).
    M_KK_override : float, optional
        If provided, use this value (GeV) for M_KK instead of the value stored
        in ``yukawa_result.params['M_KK']`` or, as a fallback, the repo's
        internal LFV convention ``M_KK = Lambda_IR``.

    Returns
    -------
    dict with keys
        'lhs' : float
            |(Ȳ_N Ȳ_N†)₁₂|  (left-hand side of the bound).
        'rhs' : float
            C × (M_KK / reference_scale)²  (right-hand side of the bound).
        'passes' : bool
            True if lhs ≤ rhs.
        'ratio' : float
            lhs / rhs.  Values > 1 violate the bound.
        'off_diagonal_12' : complex
            The raw (1,2) matrix element before taking the absolute value.
        'product_matrix' : np.ndarray
            The full 3×3 Hermitian matrix Ȳ_N Ȳ_N†.
    """
    try:
        params = yukawa_result.params
        k = params['k']
        Lambda_IR = params['Lambda_IR']
    except Exception as exc:
        raise ValueError(
            "yukawa_result must include params with 'k' and 'Lambda_IR'."
        ) from exc
    # Internal LFV convention: unless a point stores an explicit M_KK or the
    # caller overrides it, use the geometric IR scale directly.
    M_KK = float(params.get('M_KK', Lambda_IR))
    if M_KK_override is not None:
        M_KK = float(M_KK_override)

    # Build the rescaled Yukawa matrix: Ȳ_N = 2k · Y_N_matrix
    Y_N_bar_matrix = 2.0 * k * yukawa_result.Y_N_matrix

    # Compute the Hermitian product Ȳ_N Ȳ_N†
    product = Y_N_bar_matrix @ Y_N_bar_matrix.conj().T

    # Extract the (e, μ) = (1, 2) element  →  0-indexed (0, 1)
    off_diagonal_12 = product[0, 1]
    lhs = float(np.abs(off_diagonal_12))

    # Bound: C × (M_KK / 3 TeV)²
    rhs = C * (M_KK / reference_scale) ** 2

    return {
        'lhs': lhs,
        'rhs': rhs,
        'passes': lhs <= rhs,
        'ratio': lhs / rhs if rhs > 0 else float('inf'),
        'off_diagonal_12': complex(off_diagonal_12),
        'product_matrix': product,
    }


def check_mu_to_e_gamma_raw(
    Y_N_bar: np.ndarray,
    pmns: np.ndarray,
    M_KK: float,
    C: float = C_PAPER,
    reference_scale: float = 3000.0,
) -> Dict[str, Union[float, bool, complex, np.ndarray]]:
    """Standalone μ→eγ check from raw arrays (no YukawaResult needed).

    Parameters
    ----------
    Y_N_bar : array-like of shape (3,)
        Rescaled neutrino Yukawa eigenvalues Ȳ_{N_i} = 2k Y_{N_i}.
    pmns : np.ndarray of shape (3, 3)
        PMNS mixing matrix.
    M_KK : float
        KK scale in GeV for the chosen LFV convention. The repo default uses
        ``M_KK = Lambda_IR``; physical first-KK masses should be passed
        explicitly and paired with a consistent reference scale.
    C : float, optional
        Numerical coefficient.  Default 0.02 (Perez–Randall).
    reference_scale : float, optional
        Reference KK scale in GeV.  Default 3000.

    Returns
    -------
    dict
        Same keys as ``check_mu_to_e_gamma``.
    """
    Y_N_bar = np.asarray(Y_N_bar, dtype=complex)
    pmns = np.asarray(pmns, dtype=complex)

    # Ȳ_N_matrix = U · diag(Ȳ_N)
    Y_N_bar_matrix = pmns @ np.diag(Y_N_bar)

    # Ȳ_N Ȳ_N† = U diag(Ȳ²_N) U†
    product = Y_N_bar_matrix @ Y_N_bar_matrix.conj().T

    off_diagonal_12 = product[0, 1]
    lhs = float(np.abs(off_diagonal_12))
    rhs = C * (M_KK / reference_scale) ** 2

    return {
        'lhs': lhs,
        'rhs': rhs,
        'passes': lhs <= rhs,
        'ratio': lhs / rhs if rhs > 0 else float('inf'),
        'off_diagonal_12': complex(off_diagonal_12),
        'product_matrix': product,
    }
