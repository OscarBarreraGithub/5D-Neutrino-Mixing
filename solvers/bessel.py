
"""
KK–Bessel Solver (RS models)
----------------------------
This module finds KK masses by solving Bessel-function boundary conditions
for gauge bosons and fermions in Randall–Sundrum backgrounds.

Design goals:
- Clear inputs (geometry + field/BC + numerics)
- Stable root finding (uses SciPy when available)
- Minimal, readable implementation with sensible defaults
- Exact ratio equation and fast IR-only approximation

Dependencies:
- numpy
- scipy.special (Bessel J/Y)
- scipy.optimize (brentq)

"""

from __future__ import annotations
import math
import warnings
from typing import Callable, Dict, List, Optional, Tuple, Any

import numpy as np
from scipy.special import jv, yv, jn_zeros  # jn_zeros only for integer order
from scipy.optimize import brentq

# =========================
# 0) GLOBAL / DEFAULT PARAMS
# =========================

# Numerics
DEFAULT_N_ROOTS  = 3
DEFAULT_TOL      = 1e-12
DEFAULT_EXACT    = True     # True = exact ratio condition; False = IR-only Jν(x)=0 approximation

# Small cutoffs
_TINY = 1e-14
_MIN_X = 1e-12


# =========================
# 1) GEOMETRY HANDLING
# =========================
def _validate_geometry(geo: Dict[str, Any]) -> Dict[str, float]:
    """
    Validate that the geometry dict has the necessary keys.
    It expects the output format of warpConfig.baseParams.get_warp_params.
    
    Required keys:
      - Lambda_IR
      - z_v
      - epsilon
    """
    # Copy to avoid mutating input
    g = geo.copy()
    
    # Ensure essential keys exist
    required = ["Lambda_IR", "z_v", "epsilon"]
    missing = [k for k in required if k not in g]
    if missing:
        raise ValueError(f"Geometry dict missing required keys: {missing}. ")
         
    return g



# =========================
# 2) BOUNDARY CONDITION EQUATIONS
# =========================
def _nu_for(species: str, bc: str, c: Optional[float]) -> Tuple[float, str]:
    """
    Determine the Bessel order ν required by the BC for the given species.
    Returns (nu, label) where label describes which ν was used.

    Strict keys required:
      - species: 'gauge' or 'fermion'
      - bc (gauge): 'NN'
      - bc (fermion): '++' or '--'
    """
    if species == "gauge":
        if bc != "NN":
            raise ValueError(f"Gauge species requires bc='NN'. Got '{bc}'.")
        return 0.0, "nu=0 (gauge, NN)"
    elif species == "fermion":
        if c is None:
            raise ValueError("Fermion requires bulk mass parameter c.")
        alpha = abs(c + 0.5)
        if bc == "++":
            # Do not clamp to 0. Negative orders are valid (e.g. c < 0.5)
            return alpha - 1.0, "nu=alpha-1 (fermion, ++)"
        elif bc == "--":
            return alpha, "nu=alpha (fermion, --)"
        else:
            raise ValueError(f"Fermion bc must be '++' or '--'. Got '{bc}'.")
    else:
        raise ValueError(f"Species must be 'gauge' or 'fermion'. Got '{species}'.")


def _F_exact(nu: float, eps: float) -> Callable[[float], float]:
    """
    Exact quantization equation in a numerically stable cross-product form:
        F(x) = J_ν(x) Y_ν(εx) - J_ν(εx) Y_ν(x) = 0
    This avoids explicit division by Y_ν and stays finite near zeros.
    """
    def F(x: float) -> float:
        if x <= _MIN_X:
            return np.sign(x) * 1.0  # avoid evaluating at 0
        return jv(nu, x) * yv(nu, eps * x) - jv(nu, eps * x) * yv(nu, x)
    return F


def _F_ironly(nu: float) -> Callable[[float], float]:
    """
    IR-only approximation: J_ν(x) = 0
    """
    def F(x: float) -> float:
        if x <= _MIN_X:
            return 1.0
        return jv(nu, x)
    return F


# =========================
# 3) ROOT BRACKETING UTILITIES
# =========================
def _approx_jv_zeros(nu: float, n_roots: int) -> np.ndarray:
    """
    Get rough initial guesses for zeros of J_ν(x).
    - If ν is (near) an integer, use scipy.special.jn_zeros directly.
    - Otherwise, use the large-n asymptotic spacing:
        x_{ν,n} ~ (n + ν/2 - 1/4) π,  n = 1,2,3,...
      which is decent even for small n as a seed.

    Returns an array of length n_roots with positive guesses.
    """
    nu_int = int(round(nu))
    if abs(nu - nu_int) < 1e-12:
        # exact integer (use abs because jn_zeros requires n >= 0)
        return jn_zeros(abs(nu_int), n_roots)
    # non-integer: asymptotic seeds
    n = np.arange(1, n_roots + 1, dtype=float)
    return (n + 0.5 * nu - 0.25) * math.pi


def _bracket_around(F: Callable[[float], float],
                    x0: float,
                    rel: float = 0.2,
                    max_expand: int = 25) -> Optional[Tuple[float, float]]:
    """
    Try to find a sign-change bracket around x0 by symmetric expansion.
    """
    a = max(x0 * (1.0 - rel), _MIN_X)
    b = x0 * (1.0 + rel)
    fa = F(a)
    fb = F(b)
    if np.sign(fa) == 0.0:
        return (max(a * 0.9, _MIN_X), a * 1.1)
    if np.sign(fb) == 0.0:
        return (max(b * 0.9, _MIN_X), b * 1.1)
    if np.sign(fa) != np.sign(fb):
        return (a, b)

    # Expand multiplicatively until sign change (or give up)
    for _ in range(max_expand):
        a = max(a * 0.8, _MIN_X)
        b = b * 1.25
        fa = F(a)
        fb = F(b)
        if np.sign(fa) != np.sign(fb):
            return (a, b)
    return None


def _scan_for_brackets(F: Callable[[float], float],
                       x_start: float,
                       step: float,
                       n_needed: int,
                       x_max: float = 200.0) -> List[Tuple[float, float]]:
    """
    Linear scan to collect brackets with sign changes of F.
    """
    brackets = []
    x_prev = max(x_start, _MIN_X)
    f_prev = F(x_prev)
    x = x_prev
    while len(brackets) < n_needed and x < x_max:
        x = x + step
        f = F(x)
        if np.sign(f_prev) == 0.0:
            # near exact zero, create a small bracket
            brackets.append((max(x_prev * 0.9, _MIN_X), x_prev * 1.1))
        elif np.sign(f_prev) != np.sign(f):
            brackets.append((x_prev, x))
        x_prev, f_prev = x, f
    return brackets


# =========================
# 4) MAIN SOLVER
# =========================
def solve_kk(species: str,
             bc: str,
             geometry: Dict[str, float],
             c: Optional[float] = None,
             # Numerics
             n_roots: int = DEFAULT_N_ROOTS,
             exact: bool = DEFAULT_EXACT,
             tol: float = DEFAULT_TOL,
             x_max: float = 200.0,
             return_extras: bool = True) -> Tuple[np.ndarray, Dict[str, Any]]:
    """
    Find the first `n_roots` KK masses for the given field and BC.

    Parameters
    ----------
    species : 'gauge' or 'fermion'
    bc      : 'NN' (for gauge), '++' or '--' (for fermion)
    geometry: dict from `warpConfig.baseParams.get_warp_params`. REQUIRED.
    c       : fermion bulk mass parameter (dimensionless). Required for fermions.
    n_roots : number of KK roots to find
    exact   : if True, solve exact ratio equation; if False, solve IR-only Jν(x)=0
    tol     : root-finder tolerance in x
    x_max   : upper scan limit in dimensionless x
    return_extras : if True, also return dict with 'x', 'b', and 'masses'

    Returns
    -------
    masses : np.ndarray of shape (n_roots,), in GeV (same units as Lambda_IR)
    extras : dict containing
             - 'x'      : the dimensionless roots
             - 'b'      : Bessel mixing coefficients b_n (see notes below)
             - 'nu'     : the Bessel order used
             - 'labels' : metadata strings
    Notes
    -----
    - Dimensionless variable x ≡ m z_v ⇒ m = x / z_v = x Lambda_IR.
    - Exact equation uses F(x) = J_ν(x) Y_ν(εx) - J_ν(εx) Y_ν(x) = 0 (stable).
    - IR-only approximation uses J_ν(x)=0 and is excellent for ε ≪ 1.
    - b_n determination:
        * exact=True : b_n = - J_ν(x) / Y_ν(x)   (IR ratio)
        * exact=False: b_n = - J_ν(εx) / Y_ν(εx) (UV ratio, since J_ν(x)=0 at IR)
    """
    # Validate and normalize the geometry dict
    geometry = _validate_geometry(geometry)
    
    z_v = geometry["z_v"]
    eps = geometry["epsilon"]
    Lam = geometry["Lambda_IR"]

    # Determine ν for this species/BC
    nu, nu_label = _nu_for(species, bc, c)

    # Build the function F(x)
    F = _F_exact(nu, eps) if exact else _F_ironly(nu)

    # Seed guesses from J_ν zeros (or asymptotic), then bracket for F
    seeds = _approx_jv_zeros(nu, max(n_roots, 5))
    brackets: List[Tuple[float, float]] = []

    # First try symmetric brackets around each seed
    for x0 in seeds:
        if len(brackets) >= n_roots:
            break
        br = _bracket_around(F, x0, rel=0.2, max_expand=30)
        if br is not None:
            brackets.append(br)

    # If we still need more brackets, do a linear scan
    if len(brackets) < n_roots:
        # Step ~ π is a decent spacing for consecutive roots
        extra = _scan_for_brackets(F, x_start=seeds[-1] if len(seeds) else 0.5, step=math.pi, n_needed=n_roots - len(brackets), x_max=x_max)
        brackets.extend(extra)

    if len(brackets) < n_roots:
        warnings.warn(f"Only found {len(brackets)} sign-change brackets up to x={x_max}. Returning fewer roots.")
        n_roots = len(brackets)

    # Solve each bracket with Brent
    xs = []
    for (a, b) in brackets[:n_roots]:
        try:
            root = brentq(F, a, b, xtol=tol, rtol=tol, maxiter=200)
        except ValueError:
            # Failsafe: nudge ends slightly and retry once
            aa = max(a * 0.99, _MIN_X)
            bb = b * 1.01
            root = brentq(F, aa, bb, xtol=tol, rtol=tol, maxiter=200)
        xs.append(root)
    xs = np.array(xs)

    # Map to masses: m = x * Λ
    masses = xs * Lam

    # Compute b_n per mode
    bvals = []
    for x in xs:
        if exact:
            # IR ratio
            denom = yv(nu, x)
            if abs(denom) < 1e-300:
                bvals.append(np.nan)
            else:
                bvals.append(- jv(nu, x) / denom)
        else:
            # UV ratio (since Jν(x)=0 at IR in the approximation)
            xp = eps * x
            denom = yv(nu, xp)
            if abs(denom) < 1e-300:
                bvals.append(np.nan)
            else:
                bvals.append(- jv(nu, xp) / denom)
    bvals = np.array(bvals)

    extras = dict(
        x=xs,
        b=bvals,
        nu=nu,
        labels=dict(nu_label=nu_label, species=species, bc=bc, exact=exact),
        geometry=geometry,
    )
    return masses, extras



