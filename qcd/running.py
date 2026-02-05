"""QCD running coupling alpha_s(mu) via numerical RG evolution.

Evolves alpha_s from a reference point (M_Z by default) to an arbitrary
scale mu, crossing quark mass thresholds where n_f changes.

The integration uses t = ln(mu^2) as the evolution variable so the ODE is
    d(alpha_s)/dt = beta(alpha_s)
with no explicit mu dependence on the RHS.

At flavor thresholds alpha_s is matched with MS-bar decoupling constants
through the requested loop order (default: 3-loop matching for 4-loop
running). The default assumes matching at mu = m_h, so log terms vanish.

References
----------
PDG 2024, "Quantum Chromodynamics" review.
Chetyrkin, Kuhn, Sturm, Eur.Phys.J. C48 (2006) 107 — threshold matching.
"""

from typing import List, Optional, Tuple

import numpy as np
from scipy.integrate import solve_ivp

from .constants import ALPHA_S_MZ, M_Z, THRESHOLD_LIST
from .beta_function import beta_rhs
from .decoupling import match_alpha_s


def _n_f_at_scale(mu: float, thresholds: List[Tuple[float, int, int]]) -> int:
    """Determine the number of active flavors at scale *mu*."""
    if not thresholds:
        return 5  # default when thresholds disabled
    n_f = thresholds[0][1]  # n_f below the lowest threshold
    for mass, _, n_f_above in thresholds:
        if mu >= mass:
            n_f = n_f_above
        else:
            break
    return n_f


def alpha_s(
    mu: float,
    n_loops: int = 4,
    alpha_s_ref: float = ALPHA_S_MZ,
    mu_ref: float = M_Z,
    thresholds: Optional[List[Tuple[float, int, int]]] = None,
    matching_loops: Optional[int] = None,
    precision: Optional[str] = None,
    rtol: float = 1e-10,
    atol: float = 1e-12,
) -> float:
    """Compute alpha_s(mu) by integrating the MS-bar beta function.

    Parameters
    ----------
    mu : float
        Target energy scale (GeV).  Must be positive.
    n_loops : int, optional
        Loop order for the beta function (1–4).  Default 4 (N^3LO).
    alpha_s_ref : float, optional
        alpha_s at the reference scale.  Default 0.1180 (PDG 2024).
    mu_ref : float, optional
        Reference scale (GeV).  Default M_Z = 91.1876 GeV.
    thresholds : list of (float, int, int) or None, optional
        Flavor thresholds as (mass_GeV, n_f_below, n_f_above).
        Default uses PDG quark masses.  Pass ``[]`` to disable.
    matching_loops : int or None, optional
        Decoupling order at thresholds (0–3).  Default is n_loops-1,
        capped at 3 (i.e., 3-loop matching for 4-loop running).
    precision : str or None, optional
        Convenience preset that overrides *n_loops* and *matching_loops*:
        ``'low'`` selects 3-loop running with continuous threshold matching;
        ``'high'`` selects 4-loop running with 3-loop decoupling.
        Default ``None`` (use explicit parameters).
    rtol, atol : float, optional
        ODE integrator tolerances.

    Returns
    -------
    float
        The strong coupling alpha_s(mu).

    Raises
    ------
    ValueError
        If mu <= 0, n_loops not in {1,2,3,4}, or alpha_s_ref <= 0.
    RuntimeError
        If the ODE integration fails.

    Examples
    --------
    >>> from qcd import alpha_s
    >>> f"{alpha_s(1000.0):.4f}"   # 4-loop, 1 TeV
    '0.0884'
    >>> alpha_s(3000.0, precision='low')   # 3-loop, continuous matching
    >>> alpha_s(3000.0, precision='high')  # 4-loop, 3-loop decoupling
    """
    # Apply precision preset (overrides n_loops and matching_loops)
    if precision is not None:
        if precision == 'low':
            n_loops, matching_loops = 3, 0
        elif precision == 'high':
            n_loops, matching_loops = 4, 3
        else:
            raise ValueError(
                f"precision must be 'low' or 'high', got {precision!r}"
            )

    if mu <= 0:
        raise ValueError(f"mu must be positive, got {mu}")
    if alpha_s_ref <= 0:
        raise ValueError(f"alpha_s_ref must be positive, got {alpha_s_ref}")
    if not 1 <= n_loops <= 4:
        raise ValueError(f"n_loops must be 1..4, got {n_loops}")

    if thresholds is None:
        thresholds = list(THRESHOLD_LIST)
    thresholds = sorted(thresholds, key=lambda x: x[0])
    if matching_loops is None:
        matching_loops = max(0, min(n_loops - 1, 3))
    if matching_loops < 0 or matching_loops > 3:
        raise ValueError(f"matching_loops must be 0..3, got {matching_loops}")

    # Trivial case
    if np.isclose(mu, mu_ref, rtol=1e-12):
        return alpha_s_ref

    t_ref = np.log(mu_ref**2)
    t_target = np.log(mu**2)
    running_up = t_target > t_ref

    # Find thresholds between mu_ref and mu
    crossings = []
    for mass, nf_below, nf_above in thresholds:
        t_thresh = np.log(mass**2)
        if running_up and t_ref < t_thresh < t_target:
            crossings.append((t_thresh, nf_below, nf_above))
        elif not running_up and t_target < t_thresh < t_ref:
            crossings.append((t_thresh, nf_below, nf_above))

    crossings.sort(key=lambda x: x[0], reverse=(not running_up))

    # Build integration segments: (t_start, t_end, n_f, n_f_next)
    segments = []
    current_t = t_ref
    current_nf = _n_f_at_scale(mu_ref, thresholds)

    for t_thresh, nf_below, nf_above in crossings:
        next_nf = nf_above if running_up else nf_below
        segments.append((current_t, t_thresh, current_nf, next_nf))
        current_t = t_thresh
        current_nf = next_nf

    segments.append((current_t, t_target, current_nf, None))

    # Integrate each segment
    current_alpha = alpha_s_ref

    for t_start, t_end, nf, nf_next in segments:
        if np.isclose(t_start, t_end, rtol=1e-14):
            continue

        def rhs(t, y, _nf=nf, _nl=n_loops):
            return [beta_rhs(y[0], _nf, _nl)]

        sol = solve_ivp(
            rhs,
            t_span=(t_start, t_end),
            y0=[current_alpha],
            method='RK45',
            rtol=rtol,
            atol=atol,
        )

        if not sol.success:
            raise RuntimeError(
                f"ODE integration failed in segment "
                f"[mu={np.exp(t_start / 2):.2f}, mu={np.exp(t_end / 2):.2f}] "
                f"with n_f={nf}: {sol.message}"
            )

        current_alpha = float(sol.y[0, -1])
        if nf_next is not None and matching_loops > 0:
            current_alpha = match_alpha_s(
                current_alpha, n_f_from=nf, n_f_to=nf_next,
                matching_loops=matching_loops,
            )

    return current_alpha


def alpha_s_array(
    mu_values,
    n_loops: int = 4,
    alpha_s_ref: float = ALPHA_S_MZ,
    mu_ref: float = M_Z,
    thresholds: Optional[List[Tuple[float, int, int]]] = None,
    matching_loops: Optional[int] = None,
    precision: Optional[str] = None,
) -> np.ndarray:
    """Compute alpha_s at multiple scales.

    Parameters
    ----------
    mu_values : array-like
        Energy scales (GeV).
    n_loops : int, optional
        Loop order (1–4).  Default 4.
    alpha_s_ref : float, optional
        Reference value.  Default ALPHA_S_MZ.
    mu_ref : float, optional
        Reference scale (GeV).  Default M_Z.
    thresholds : list or None, optional
        Flavor thresholds.  Default from constants.py.
    matching_loops : int or None, optional
        Decoupling order at thresholds (0–3).
    precision : str or None, optional
        ``'low'`` or ``'high'`` preset.  See :func:`alpha_s`.

    Returns
    -------
    np.ndarray
        Array of alpha_s values, same shape as *mu_values*.
    """
    mu_arr = np.asarray(mu_values, dtype=float)
    result = np.empty_like(mu_arr)
    for i, mu_val in enumerate(mu_arr.flat):
        result.flat[i] = alpha_s(
            mu_val, n_loops=n_loops,
            alpha_s_ref=alpha_s_ref, mu_ref=mu_ref,
            thresholds=thresholds,
            matching_loops=matching_loops,
            precision=precision,
        )
    return result
