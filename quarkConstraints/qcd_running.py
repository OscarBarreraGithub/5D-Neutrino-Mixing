"""QCD renormalization group evolution for Delta F=2 Wilson coefficients.

This module implements leading-log QCD running of the ΔF=2 Wilson coefficients
from the KK matching scale (multi-TeV) down to the hadronic scale (μ = 2 GeV).
The anomalous dimension matrix follows Buras, Misiak, and Urban (NPB 2000).

Operator basis: {O1_VLL, O1_VRR, O4_LR, O5_LR}.  The VLL/VRR operators are
the BMU current-current operators.  The LR entries are the conventional scalar
O4/O5 operators used with the code's B4/B5 matrix elements.  Their LO running
is obtained from the BMU LR basis by the four-dimensional map
``Q1_LR^BMU = -2 O5_LR`` and ``Q2_LR^BMU = O4_LR``.

Key physics:
  - VLL and VRR operators run multiplicatively (no mixing with LR sector).
  - O4_LR and O5_LR mix under QCD running via a 2×2 anomalous dimension matrix;
    the large scalar-basis C5 -> C4 entry is the mapped BMU LR off-diagonal.
  - The LR operators receive significant enhancement from running, which
    strengthens the kaon mixing constraint.

Threshold matching at m_t, m_b, and m_c is implemented with continuous Wilson
coefficients (no finite matching corrections at leading order).
"""

from __future__ import annotations

import math

import numpy as np

from qcd.constants import M_TOP_MS

# ---------------------------------------------------------------------------
# Standard physical constants
# ---------------------------------------------------------------------------

_ALPHA_S_MZ_DEFAULT = 0.1179
_MZ_DEFAULT = 91.1876
_MT_DEFAULT = M_TOP_MS
_MB_DEFAULT = 4.18
_MC_DEFAULT = 1.27

# ---------------------------------------------------------------------------
# One-loop anomalous dimensions for ΔF=2 operators (BMU and mapped LR basis)
# ---------------------------------------------------------------------------

# VLL (same for VRR): single anomalous dimension.
# γ_VLL = 6(N_c-1)/N_c = 4 for N_c=3.
# Positive γ → coefficient suppressed at lower scales in the BMU Wilson
# evolution convention.
_GAMMA_VLL = 4.0

# BMU LR sector for coefficients ordered as [C1_LR^BMU, C2_LR^BMU] is
# [[2, 0], [12, -16]].  The code stores conventional scalar LR coefficients
# ordered as [C4_LR, C5_LR], with C_BMU = (-C5/2, C4).  Conjugating the BMU
# matrix gives the scalar-basis coefficient ADM below.
_GAMMA_LR = np.array([
    [-16.0, -6.0],
    [0.0, 2.0],
], dtype=float)


def _beta0(nf: int) -> float:
    """One-loop QCD beta function coefficient β₀ = (33 - 2n_f)/3."""
    return (33.0 - 2.0 * nf) / 3.0


def _nf_for_scale(mu: float, *, m_t: float, m_b: float, m_c: float) -> int:
    """Return the number of active quark flavors at scale μ."""
    if mu > m_t:
        return 6
    if mu > m_b:
        return 5
    if mu > m_c:
        return 4
    return 3


def run_alpha_s(
    mu: float,
    *,
    alpha_s_mz: float = _ALPHA_S_MZ_DEFAULT,
    m_z: float = _MZ_DEFAULT,
    m_t: float = _MT_DEFAULT,
    m_b: float = _MB_DEFAULT,
    m_c: float = _MC_DEFAULT,
) -> float:
    """One-loop running of α_s with threshold matching at m_t, m_b, and m_c.

    The running uses the one-loop formula:
        α_s(μ) = α_s(μ_ref) / (1 + (β₀ α_s(μ_ref)/(2π)) ln(μ/μ_ref))

    At each quark threshold (m_t, m_b, m_c), the coupling is matched continuously
    and the number of active flavors changes.

    Parameters
    ----------
    mu : float
        Target renormalization scale in GeV.
    alpha_s_mz : float
        α_s at the Z mass (default 0.1179).
    m_z : float
        Z boson mass in GeV (default 91.1876).
    m_t : float
        Top quark threshold in GeV (default: qcd.constants.M_TOP_MS).
    m_b : float
        Bottom quark mass threshold in GeV (default 4.18).
    m_c : float
        Charm quark mass threshold in GeV (default 1.27).

    Returns
    -------
    float
        α_s(μ) at one-loop accuracy.

    Raises
    ------
    ValueError
        If μ is non-positive or the running produces an unphysical result.
    """
    if mu <= 0.0:
        raise ValueError("mu must be positive")
    if alpha_s_mz <= 0.0:
        raise ValueError("alpha_s_mz must be positive")
    if not (m_c < m_b < m_z < m_t):
        raise ValueError("thresholds must satisfy m_c < m_b < m_z < m_t")

    # Build the ordered list of thresholds between mu and m_z.
    # We always start at m_z and run to mu, crossing thresholds as needed.

    # Step 1: run from m_z down to m_b (nf=5) if mu < m_z
    # or run from m_z up to mu if mu > m_z
    def _run_segment(alpha_ref: float, mu_ref: float, mu_target: float, nf: int) -> float:
        """Run α_s from mu_ref to mu_target with nf active flavors."""
        if mu_ref == mu_target:
            return alpha_ref
        b0 = _beta0(nf)
        log_ratio = math.log(mu_target / mu_ref)
        denominator = 1.0 + b0 * alpha_ref / (2.0 * math.pi) * log_ratio
        if denominator <= 0.0:
            raise ValueError(
                f"Landau pole encountered running alpha_s from {mu_ref:.1f} "
                f"to {mu_target:.1f} GeV with nf={nf}"
            )
        return alpha_ref / denominator

    if mu >= m_z:
        # Running upward from alpha_s^(5)(M_Z), crossing the top threshold if
        # needed.  Matching is continuous at LO.
        if mu <= m_t:
            return _run_segment(alpha_s_mz, m_z, mu, 5)
        alpha_at_top = _run_segment(alpha_s_mz, m_z, m_t, 5)
        return _run_segment(alpha_at_top, m_t, mu, 6)

    # Running downward from m_z
    alpha_current = alpha_s_mz
    mu_current = m_z

    # Thresholds to cross, ordered from high to low.  M_Z is below m_t, so the
    # downward path from the default reference point only crosses b and c.
    thresholds = [(m_b, 5, 4), (m_c, 4, 3)]

    for threshold, nf_above, nf_below in thresholds:
        if mu >= threshold:
            # Target is above this threshold; run with nf_above and stop
            return _run_segment(alpha_current, mu_current, mu, nf_above)
        # Run down to the threshold
        alpha_current = _run_segment(alpha_current, mu_current, threshold, nf_above)
        mu_current = threshold

    # Below all thresholds (mu < m_c): run with nf=3
    return _run_segment(alpha_current, mu_current, mu, 3)


def _vll_evolution_factor(
    alpha_s_low: float,
    alpha_s_high: float,
    nf: int,
) -> float:
    """Compute the VLL/VRR evolution factor for one nf segment.

    U_VLL = (α_s(μ_high)/α_s(μ_low))^(γ_VLL/(2β₀))
    """
    b0 = _beta0(nf)
    exponent = _GAMMA_VLL / (2.0 * b0)
    return (alpha_s_high / alpha_s_low) ** exponent


def _lr_evolution_matrix(
    alpha_s_low: float,
    alpha_s_high: float,
    nf: int,
) -> np.ndarray:
    """Compute the 2×2 LR evolution matrix for one nf segment.

    U_LR = M × diag((α_s_high/α_s_low)^(d_i/(2β₀))) × M⁻¹

    where d_i are the eigenvalues of γ_LR and M is the eigenvector matrix.
    """
    b0 = _beta0(nf)
    ratio = alpha_s_high / alpha_s_low

    # Diagonalize the anomalous dimension matrix
    eigenvalues, M = np.linalg.eig(_GAMMA_LR)

    # Build diagonal evolution matrix
    diag_factors = np.array([ratio ** (d / (2.0 * b0)) for d in eigenvalues])
    D = np.diag(diag_factors)

    # U = M D M^{-1}
    M_inv = np.linalg.inv(M)
    return M @ D @ M_inv


def evolve_deltaf2_wilsons(
    c_vll: complex,
    c_vrr: complex,
    c4_lr: complex,
    c5_lr: complex,
    mu_high: float,
    mu_low: float = 2.0,
    *,
    alpha_s_mz: float = _ALPHA_S_MZ_DEFAULT,
    m_z: float = _MZ_DEFAULT,
    m_t: float = _MT_DEFAULT,
    m_b: float = _MB_DEFAULT,
    m_c: float = _MC_DEFAULT,
) -> tuple[complex, complex, complex, complex]:
    """Evolve ΔF=2 Wilson coefficients from mu_high to mu_low using LO QCD RG.

    The evolution uses leading-log anomalous dimensions from Buras, Misiak,
    and Urban (NPB 2000). Threshold matching at m_t, m_b, and m_c is handled with
    continuous Wilson coefficients.

    Parameters
    ----------
    c_vll, c_vrr : complex
        VLL and VRR Wilson coefficients at mu_high.
    c4_lr, c5_lr : complex
        LR Wilson coefficients at mu_high.
    mu_high : float
        Matching scale (typically M_KK) in GeV.
    mu_low : float
        Hadronic scale in GeV (default 2.0).
    alpha_s_mz : float
        α_s(M_Z) (default 0.1179).
    m_z : float
        Z mass in GeV.
    m_t : float
        t quark threshold in GeV.
    m_b : float
        b quark threshold in GeV.
    m_c : float
        c quark threshold in GeV.

    Returns
    -------
    tuple[complex, complex, complex, complex]
        (c_vll_low, c_vrr_low, c4_lr_low, c5_lr_low) at scale mu_low.

    Raises
    ------
    ValueError
        If mu_high <= mu_low or other parameter constraints are violated.
    """
    if mu_high <= mu_low:
        raise ValueError(
            f"mu_high ({mu_high}) must be greater than mu_low ({mu_low})"
        )
    if mu_low <= 0.0:
        raise ValueError("mu_low must be positive")
    if mu_high <= 0.0:
        raise ValueError("mu_high must be positive")

    alpha_s_kwargs = dict(
        alpha_s_mz=alpha_s_mz, m_z=m_z, m_t=m_t, m_b=m_b, m_c=m_c,
    )

    # Determine which thresholds lie between mu_low and mu_high
    thresholds = sorted(
        [t for t in [m_t, m_b, m_c] if mu_low < t < mu_high],
        reverse=True,
    )

    # Build the sequence of segments: mu_high -> threshold_1 -> ... -> mu_low
    boundaries = [mu_high] + thresholds + [mu_low]

    # Initialize the cumulative evolution
    vll_factor = 1.0
    vrr_factor = 1.0
    lr_matrix = np.eye(2, dtype=float)

    for i in range(len(boundaries) - 1):
        mu_upper = boundaries[i]
        mu_lower = boundaries[i + 1]
        nf = _nf_for_scale(mu_upper, m_t=m_t, m_b=m_b, m_c=m_c)

        alpha_upper = run_alpha_s(mu_upper, **alpha_s_kwargs)
        alpha_lower = run_alpha_s(mu_lower, **alpha_s_kwargs)

        # VLL and VRR evolution (multiplicative)
        vll_factor *= _vll_evolution_factor(alpha_lower, alpha_upper, nf)
        vrr_factor *= _vll_evolution_factor(alpha_lower, alpha_upper, nf)

        # LR evolution (2×2 matrix)
        segment_lr = _lr_evolution_matrix(alpha_lower, alpha_upper, nf)
        lr_matrix = segment_lr @ lr_matrix

    # Apply evolution to Wilson coefficients
    c_vll_low = complex(c_vll * vll_factor)
    c_vrr_low = complex(c_vrr * vrr_factor)

    lr_vec = np.array([complex(c4_lr), complex(c5_lr)])
    lr_evolved = lr_matrix @ lr_vec
    c4_lr_low = complex(lr_evolved[0])
    c5_lr_low = complex(lr_evolved[1])

    return c_vll_low, c_vrr_low, c4_lr_low, c5_lr_low


__all__ = [
    "evolve_deltaf2_wilsons",
    "run_alpha_s",
]
