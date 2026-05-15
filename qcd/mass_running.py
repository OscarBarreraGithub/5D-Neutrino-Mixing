"""MS-bar quark mass RG running with threshold matching.

Implements the standard MS-bar running

    d ln m(mu) / d ln mu^2 = -gamma_m(alpha_s)

with the quark-mass anomalous dimension expanded as

    gamma_m(a) = sum_{n=0}^{N-1} gamma_m^(n) * a^{n+1},
    a = alpha_s(mu) / (4 pi),

through 4-loop order. Heavy-quark threshold matching is applied at
``m_b`` and (if needed) ``m_c`` using
:func:`qcd.decoupling.match_msbar_mass` (3-loop matching at mu = m_h).

References
----------
PDG 2024, "Quantum Chromodynamics" review.
Chetyrkin, Phys. Lett. B 404 (1997) 161 — 4-loop gamma_m.
Vermaseren, Larin, van Ritbergen, Phys. Lett. B 405 (1997) 327.
Chetyrkin, Kniehl, Steinhauser, PRL 79 (1997) 2184 — mass-matching at thresholds.
"""

from __future__ import annotations

from typing import Callable, List, Optional, Tuple

import numpy as np
from scipy.integrate import solve_ivp

from .beta_function import beta_rhs
from .constants import M_BOTTOM, M_CHARM, M_TOP_MS
from .decoupling import match_alpha_s, match_msbar_mass

_ZETA3 = 1.2020569031595942
_ZETA4 = np.pi ** 4 / 90.0
_ZETA5 = 1.0369277551433699


def _gamma_m_coeffs(n_f: int) -> Tuple[float, float, float, float]:
    """Return (gamma0, gamma1, gamma2, gamma3) at given n_f.

    Coefficients in the convention

        gamma_m(a) = sum_n gamma_n * a^{n+1},   a = alpha_s/(4 pi).

    Numerical values match Chetyrkin (1997) and the PDG QCD review.
    """
    nf = float(n_f)
    g0 = 4.0
    g1 = 202.0 / 3.0 - 20.0 / 9.0 * nf
    g2 = (
        1249.0
        - (2216.0 / 27.0 + 160.0 / 3.0 * _ZETA3) * nf
        - 140.0 / 81.0 * nf * nf
    )
    # 4-loop mass anomalous dimension (Chetyrkin 1997)
    g3 = (
        4603055.0 / 162.0
        + 135680.0 / 27.0 * _ZETA3
        - 8800.0 * _ZETA5
        + nf
        * (
            -91723.0 / 27.0
            - 34192.0 / 9.0 * _ZETA3
            + 880.0 * _ZETA4
            + 18400.0 / 9.0 * _ZETA5
        )
        + nf * nf
        * (
            5242.0 / 243.0
            + 800.0 / 9.0 * _ZETA3
            - 160.0 / 3.0 * _ZETA4
        )
        + nf * nf * nf
        * (-332.0 / 243.0 + 64.0 / 27.0 * _ZETA3)
    )
    return g0, g1, g2, g3


def _gamma_m(alpha_s: float, n_f: int, n_loops: int) -> float:
    """Return gamma_m(alpha_s) summed up to n_loops."""
    if n_loops < 1 or n_loops > 4:
        raise ValueError("n_loops must be 1..4")
    g0, g1, g2, g3 = _gamma_m_coeffs(n_f)
    a = alpha_s / (4.0 * np.pi)
    out = g0 * a
    if n_loops >= 2:
        out += g1 * a ** 2
    if n_loops >= 3:
        out += g2 * a ** 3
    if n_loops >= 4:
        out += g3 * a ** 4
    return float(out)


def _ordered_thresholds_between(
    mu_start: float, mu_end: float
) -> List[Tuple[float, int, int]]:
    """Return ordered list of mass thresholds strictly between mu_start and mu_end.

    Each tuple is (mu_thresh_GeV, n_f_below, n_f_above), ordered in the
    direction of running.
    """
    raw = [
        (M_CHARM, 3, 4),
        (M_BOTTOM, 4, 5),
        (M_TOP_MS, 5, 6),
    ]
    lo, hi = (mu_start, mu_end) if mu_end > mu_start else (mu_end, mu_start)
    crossings = [(m, nb, na) for (m, nb, na) in raw if lo < m < hi]
    crossings.sort(key=lambda x: x[0], reverse=(mu_end < mu_start))
    return crossings


def run_msbar_mass(
    m_ref: float,
    mu_ref: float,
    mu_target: float,
    n_f_ref: int,
    *,
    n_loops: int = 4,
    matching_loops: int = 3,
    alpha_s_callable: Optional[Callable[[float, int], float]] = None,
    rtol: float = 1e-10,
    atol: float = 1e-14,
) -> float:
    """Run an MS-bar quark mass from ``mu_ref`` to ``mu_target``.

    Threshold matching at ``m_b`` and ``m_c`` is performed using
    :func:`qcd.decoupling.match_msbar_mass` (3-loop matching at ``mu = m_h``).
    The running uses ``alpha_s`` from :mod:`qcd.running` with the same
    ``n_loops``-1 (capped at 3) matching prescription.

    Parameters
    ----------
    m_ref
        Mass at the reference scale (GeV).
    mu_ref
        Reference scale (GeV).
    mu_target
        Target scale (GeV). Must be strictly positive.
    n_f_ref
        Active MS-bar flavor count at ``mu_ref``.
    n_loops
        Loop order for both ``gamma_m`` and the underlying ``beta`` (1..4).
        Default 4.
    matching_loops
        Loop order for mass and alpha_s matching at thresholds (0..3).
        Default 3.
    alpha_s_callable
        Optional callable ``(mu, n_f) -> alpha_s`` for in-segment alpha_s
        evaluation. If None, alpha_s is integrated alongside the mass on
        each segment using :func:`qcd.beta_function.beta_rhs`.

    Returns
    -------
    float
        ``m(mu_target)`` in GeV.
    """
    if m_ref <= 0.0:
        raise ValueError("m_ref must be positive")
    if mu_ref <= 0.0 or mu_target <= 0.0:
        raise ValueError("mu_ref, mu_target must be positive")
    if n_f_ref not in (3, 4, 5, 6):
        raise ValueError("n_f_ref must be in {3,4,5,6}")

    if np.isclose(mu_ref, mu_target, rtol=1e-14):
        return float(m_ref)

    # Anchor alpha_s. We integrate alpha_s and the mass jointly on each
    # segment using beta_rhs, applying threshold matching at boundaries.
    from .running import alpha_s as _alpha_s

    alpha_at_ref = _alpha_s(
        mu_ref,
        n_loops=n_loops,
        matching_loops=min(matching_loops, max(0, n_loops - 1), 3),
    )

    crossings = _ordered_thresholds_between(mu_ref, mu_target)
    going_up = mu_target > mu_ref

    segments: List[Tuple[float, float, int, Optional[int]]] = []
    current_mu = mu_ref
    current_nf = n_f_ref
    for mu_thresh, nf_below, nf_above in crossings:
        next_nf = nf_above if going_up else nf_below
        segments.append((current_mu, mu_thresh, current_nf, next_nf))
        current_mu = mu_thresh
        current_nf = next_nf
    segments.append((current_mu, mu_target, current_nf, None))

    current_alpha = alpha_at_ref
    current_m = float(m_ref)
    eff_match_loops = min(matching_loops, max(0, n_loops - 1), 3)

    for mu_start, mu_end, nf, nf_next in segments:
        if np.isclose(mu_start, mu_end, rtol=1e-14):
            # Pure threshold step (zero-length segment) — just match.
            pass
        else:
            t_start = np.log(mu_start ** 2)
            t_end = np.log(mu_end ** 2)

            def rhs(t, y, _nf=nf, _nl=n_loops):
                a_s = y[0]
                m = y[1]
                d_alpha = beta_rhs(a_s, _nf, _nl)
                # d ln m / d t = - gamma_m(a_s)
                d_m = -_gamma_m(a_s, _nf, _nl) * m
                return [d_alpha, d_m]

            sol = solve_ivp(
                rhs,
                t_span=(t_start, t_end),
                y0=[current_alpha, current_m],
                method="RK45",
                rtol=rtol,
                atol=atol,
            )
            if not sol.success:
                raise RuntimeError(
                    f"mass-running ODE failed in segment "
                    f"[{mu_start}, {mu_end}] n_f={nf}: {sol.message}"
                )
            current_alpha = float(sol.y[0, -1])
            current_m = float(sol.y[1, -1])

        if nf_next is not None and eff_match_loops > 0:
            # Apply alpha_s and mass matching at mu = m_h (mu_end).
            current_alpha = match_alpha_s(
                current_alpha,
                n_f_from=nf,
                n_f_to=nf_next,
                matching_loops=eff_match_loops,
            )
            direction = "up" if (nf_next > nf) else "down"
            n_f_high = max(nf, nf_next)
            current_m = match_msbar_mass(
                current_m,
                alpha_s=current_alpha,
                direction=direction,
                n_f_high=n_f_high,
                matching_loops=eff_match_loops,
            )

    return float(current_m)
