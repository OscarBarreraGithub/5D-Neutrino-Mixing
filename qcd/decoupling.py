"""Decoupling relations for alpha_s across heavy-quark thresholds.

Implements MS-bar decoupling constants through O(alpha_s^3) at mu = m_h.
The relation is (Chetyrkin–Kuhn–Steinhauser):

    alpha_s^(n_l)(mu) = alpha_s^(n_l+1)(mu) *
        [1 + c2 * (alpha_s/pi)^2 + c3 * (alpha_s/pi)^3 + ...]

with c2 = 11/72 and
    c3 = 564731/124416 - 82043/27648 * zeta(3) - (2633/31104) * n_l

We use mu = m_h by default, so log terms vanish. This matches the standard
RunDec choice and yields per-mille level corrections at the top threshold.
"""

from typing import Tuple

import numpy as np

# Apery's constant zeta(3)
_ZETA3 = 1.2020569031595942


def _coeffs_msbar(n_l: int) -> Tuple[float, float]:
    """Return (c2, c3) for MS-bar decoupling at mu = m_h.

    Parameters
    ----------
    n_l : int
        Number of light flavors (below threshold).
    """
    c2 = 11.0 / 72.0
    c3 = (
        564731.0 / 124416.0
        - 82043.0 * _ZETA3 / 27648.0
        - 2633.0 * n_l / 31104.0
    )
    return c2, c3


def match_alpha_s(
    alpha_s: float,
    n_f_from: int,
    n_f_to: int,
    matching_loops: int,
) -> float:
    """Match alpha_s across a heavy-quark threshold (MS-bar, mu = m_h).

    Parameters
    ----------
    alpha_s : float
        Coupling in the initial theory (n_f_from).
    n_f_from, n_f_to : int
        Flavor numbers before and after the threshold. Must differ by 1.
    matching_loops : int
        Matching order: 0/1 = no matching, 2 = include c2, 3 = include c3.

    Returns
    -------
    float
        Coupling in the target theory (n_f_to).
    """
    if n_f_from == n_f_to:
        return float(alpha_s)
    if abs(n_f_from - n_f_to) != 1:
        raise ValueError("Matching only supports adjacent thresholds (n_f +/- 1).")
    if matching_loops <= 1:
        return float(alpha_s)

    n_l = min(n_f_from, n_f_to)
    c2, c3 = _coeffs_msbar(n_l)
    if matching_loops < 3:
        c3 = 0.0

    a = alpha_s / np.pi

    if n_f_to < n_f_from:
        # Downward: alpha_s^(n_l) from alpha_s^(n_l+1)
        return float(alpha_s * (1.0 + c2 * a**2 + c3 * a**3))

    # Upward: invert the series to this order
    return float(alpha_s * (1.0 - c2 * a**2 - c3 * a**3))
