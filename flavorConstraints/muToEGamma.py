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

    with C = 0.028 from naive dimensional analysis of the dipole operator
    (Perez & Randall, arXiv:0805.4652).

Reference:
    Perez & Randall, arXiv:0805.4652
"""

from typing import Dict, Union

import numpy as np


def check_mu_to_e_gamma(
    yukawa_result,
    C: float = 0.028,
    reference_scale: float = 3000.0,
) -> Dict[str, Union[float, bool, complex, np.ndarray]]:
    """Check the μ→eγ dipole constraint on the neutrino Yukawa.

    Parameters
    ----------
    yukawa_result : yukawa.compute_yukawas.YukawaResult
        Output of ``compute_all_yukawas()``.
    C : float, optional
        Numerical coefficient in the bound.  Default 0.028.
    reference_scale : float, optional
        Reference KK scale in GeV (denominator).  Default 3000 (= 3 TeV).

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
        M_KK = params['Lambda_IR']
    except Exception as exc:
        raise ValueError(
            "yukawa_result must include params with 'k' and 'Lambda_IR'."
        ) from exc

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
    C: float = 0.028,
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
        KK mass scale in GeV (= Λ_IR).
    C : float, optional
        Numerical coefficient.  Default 0.028.
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
