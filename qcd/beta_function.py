"""QCD beta function coefficients in the MS-bar scheme.

The RG equation for the strong coupling is:

    mu^2 d(alpha_s)/d(mu^2) = beta(alpha_s)

    beta(alpha_s) = -(alpha_s^2 / (4*pi)) * [beta_0
                     + beta_1 * (alpha_s / (4*pi))
                     + beta_2 * (alpha_s / (4*pi))^2
                     + beta_3 * (alpha_s / (4*pi))^3 + ...]

References
----------
van Ritbergen, Vermaseren, Larin, Phys.Lett. B400 (1997) 379
    [hep-ph/9701390] — 4-loop coefficient.
PDG 2024, "Quantum Chromodynamics" review.
"""

import numpy as np

# Apery's constant  zeta(3)
_ZETA3 = 1.2020569031595942


def beta_0(n_f: int) -> float:
    """1-loop beta function coefficient.

    beta_0 = 11 - 2*n_f/3
    """
    return 11.0 - 2.0 * n_f / 3.0


def beta_1(n_f: int) -> float:
    """2-loop beta function coefficient.

    beta_1 = 102 - 38*n_f/3
    """
    return 102.0 - 38.0 * n_f / 3.0


def beta_2(n_f: int) -> float:
    """3-loop beta function coefficient.

    beta_2 = 2857/2 - 5033*n_f/18 + 325*n_f^2/54
    """
    return 2857.0 / 2.0 - 5033.0 * n_f / 18.0 + 325.0 * n_f**2 / 54.0


def beta_3(n_f: int) -> float:
    """4-loop beta function coefficient.

    From van Ritbergen, Vermaseren, Larin (1997).
    """
    return (
        149753.0 / 6.0 + 3564.0 * _ZETA3
        - (1078361.0 / 162.0 + 6508.0 * _ZETA3 / 27.0) * n_f
        + (50065.0 / 162.0 + 6472.0 * _ZETA3 / 81.0) * n_f**2
        + 1093.0 * n_f**3 / 729.0
    )


_BETA_FUNCS = [beta_0, beta_1, beta_2, beta_3]


def beta_coefficients(n_f: int, n_loops: int) -> np.ndarray:
    """Return the first *n_loops* beta function coefficients.

    Parameters
    ----------
    n_f : int
        Number of active quark flavors (0–6).
    n_loops : int
        Loop order (1–4).

    Returns
    -------
    np.ndarray
        Array [beta_0, beta_1, ...] of length *n_loops*.

    Raises
    ------
    ValueError
        If *n_f* or *n_loops* is out of range.
    """
    if not 0 <= n_f <= 6:
        raise ValueError(f"n_f must be 0..6, got {n_f}")
    if not 1 <= n_loops <= 4:
        raise ValueError(f"n_loops must be 1..4, got {n_loops}")
    return np.array([_BETA_FUNCS[i](n_f) for i in range(n_loops)])


def beta_rhs(alpha_s: float, n_f: int, n_loops: int) -> float:
    """Evaluate the beta function for the ODE integrator.

    Computes  mu^2 d(alpha_s)/d(mu^2) = -(alpha_s^2/(4pi)) sum_i beta_i (alpha_s/(4pi))^i.

    Parameters
    ----------
    alpha_s : float
        Current value of the strong coupling.
    n_f : int
        Number of active quark flavors.
    n_loops : int
        Loop order (1–4).

    Returns
    -------
    float
        Value of beta(alpha_s).
    """
    coeffs = beta_coefficients(n_f, n_loops)
    a_over_4pi = alpha_s / (4.0 * np.pi)

    series = sum(b * a_over_4pi**i for i, b in enumerate(coeffs))
    return -(alpha_s**2 / (4.0 * np.pi)) * series
