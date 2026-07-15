"""Decoupling relations for alpha_s and MS-bar masses across heavy-quark thresholds.

Implements MS-bar alpha_s decoupling constants through O(alpha_s^3) at
mu = m_h (Chetyrkin, Kniehl, Steinhauser, PRL 79 (1997) 2184):

    alpha_s^(n_l)(mu) = alpha_s^(n_l+1)(mu) *
        [1 + c2 * (alpha_s/pi)^2 + c3 * (alpha_s/pi)^3 + ...]

with c2 = 11/72 and
    c3 = 564731/124416 - 82043/27648 * zeta(3) - (2633/31104) * n_l.

Also implements MS-bar mass matching at mu = m_h to 2-loop and 3-loop
(Chetyrkin, Kniehl, Steinhauser 1997 / Bekavac et al. 2007 / RunDec):

    m^(n_l)(m_h) = m^(n_l+1)(m_h) * [1 + d2 * a^2 + d3 * a^3 + ...]

with d2 = 89/432 and d3(n_l) from the CKS closed form (see
_coeffs_msbar_mass).

We use mu = m_h by default, so log terms vanish. Corrections are per-mille
at the bottom threshold and well below per-mille at the charm threshold.
"""

from typing import Optional, Tuple

import numpy as np

# Apery's constant zeta(3) and zeta(4)
_ZETA3 = 1.2020569031595942
_ZETA4 = np.pi ** 4 / 90.0
# CKS Eq. (21): B4 = 16 Li4(1/2) - 13/2 zeta(4)
#               - 4 zeta(2) log^2(2) + 2/3 log^4(2)
_B4 = -1.7628000870737708


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


def _coeffs_msbar_mass(n_l: int) -> Tuple[float, float]:
    """Return (d2, d3) for MS-bar mass decoupling at mu = m_h.

    The series is

        m^(n_l)(m_h) = m^(n_l+1)(m_h) * [1 + d2 * a^2 + d3 * a^3 + ...]

    with a = alpha_s^(n_l)(m_h) / pi.

    The 2-loop coefficient is the well-known d2 = 89/432 (independent of
    n_l). The 3-loop coefficient is evaluated from the CKS closed form at
    mu = m_h, where the logarithms vanish.

    Parameters
    ----------
    n_l : int
        Number of light flavors (below the threshold).
    """
    d2 = 89.0 / 432.0
    # CKS Eq. (20), evaluated at mu = m_h:
    # d3 = 2951/2916 - 407/864 * zeta(3) + 5/4 * zeta(4) - B4/36
    #      + n_l * (1327/11664 - 2/27 * zeta(3))
    # Numerically, d3 = 1.8476267422 + 0.0247276094 * n_l.
    d3 = (
        2951.0 / 2916.0
        - 407.0 * _ZETA3 / 864.0
        + 5.0 * _ZETA4 / 4.0
        - _B4 / 36.0
        + n_l * (1327.0 / 11664.0 - 2.0 * _ZETA3 / 27.0)
    )
    return d2, d3


def match_msbar_mass(
    m: float,
    *,
    alpha_s: float,
    direction: str,
    n_f_high: int,
    matching_loops: int = 3,
    n_f_from: Optional[int] = None,
    n_f_to: Optional[int] = None,
) -> float:
    """Match an MS-bar quark mass across a heavy-quark threshold at mu = m_h.

    Convention: ``n_f_high`` is the higher of the two flavor counts (so the
    decoupled side has ``n_f_low = n_f_high - 1``); ``alpha_s`` is the
    coupling at mu = m_h **in the direction we are moving toward** (the
    convention used by :func:`match_alpha_s` after one applies the alpha_s
    matching first).

    Parameters
    ----------
    m : float
        Input MS-bar mass at mu = m_h (GeV).
    alpha_s : float
        alpha_s at mu = m_h on the target side of the threshold.
    direction : {'up', 'down'}
        ``'up'`` matches from ``n_l`` to ``n_l + 1`` (going to higher mu);
        ``'down'`` matches from ``n_l + 1`` to ``n_l`` (going to lower mu).
    n_f_high : int
        Higher of the two flavor counts spanning the threshold.
    matching_loops : int, optional
        Matching order: 0/1 = identity, 2 = include d2, 3 = include d3.
    n_f_from, n_f_to : int, optional
        Explicit flavor counts before and after matching. If both are supplied
        and are equal, the matching is an identity, matching
        :func:`match_alpha_s` behavior for equal-flavor crossings.
    """
    if (n_f_from is None) != (n_f_to is None):
        raise ValueError("n_f_from and n_f_to must be supplied together")
    if n_f_from is not None and n_f_to is not None:
        if n_f_from not in (3, 4, 5, 6) or n_f_to not in (3, 4, 5, 6):
            raise ValueError("n_f_from and n_f_to must be in {3,4,5,6}")
        if n_f_from == n_f_to:
            return float(m)
        if abs(n_f_from - n_f_to) != 1:
            raise ValueError("Matching only supports adjacent thresholds (n_f +/- 1).")
        expected_direction = "up" if n_f_to > n_f_from else "down"
        if direction != expected_direction:
            raise ValueError("direction is inconsistent with n_f_from/n_f_to")
        if n_f_high != max(n_f_from, n_f_to):
            raise ValueError("n_f_high is inconsistent with n_f_from/n_f_to")

    if matching_loops <= 1:
        return float(m)
    if n_f_high not in (4, 5, 6):
        raise ValueError("n_f_high must be in {4,5,6}")
    if direction not in ("up", "down"):
        raise ValueError("direction must be 'up' or 'down'")

    n_l = n_f_high - 1
    d2, d3 = _coeffs_msbar_mass(n_l)
    if matching_loops < 3:
        d3 = 0.0

    a = alpha_s / np.pi

    if direction == "down":
        # m^(n_l)(m_h) = m^(n_l+1)(m_h) * (1 + d2 a^2 + d3 a^3)
        return float(m * (1.0 + d2 * a ** 2 + d3 * a ** 3))
    # Upward: invert the series to this order.
    return float(m * (1.0 - d2 * a ** 2 - d3 * a ** 3))
