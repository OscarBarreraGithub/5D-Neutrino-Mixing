"""neutrinoValues.py

This module provides the oscillation inputs used to compute neutrino masses
and the PMNS matrix.

Default best-fit values are taken from NuFIT 6.1 (2025), using the
"IC24 with SK atmospheric data" variant posted on the NuFIT results page for
data available in November 2025.

Both normal and inverted mass orderings are supported.
"""

from __future__ import annotations

from typing import Tuple

import numpy as np

NUFIT_DATASET = "NuFIT 6.1 (2025), IC24 with SK atmospheric data"
NUFIT_REFERENCE = "https://www.nu-fit.org/?q=node/309"

# Constraint on the sum of neutrino masses from cosmology (in eV).
# This is kept as a separate phenomenological prior and is not part of the
# NuFIT oscillation-only dataset above.
sum_mass_constraint = 0.082  # eV

# NuFIT 6.1 best-fit mass-squared differences (in eV^2), IC24 with SK-atm.
# The published table quotes Δm^2_21 and Δm^2_3ℓ, with ℓ = 1 for NO and ℓ = 2
# for IO. We keep those direct quantities and provide legacy aliases below.
delta_m21_sq = 7.537e-5
delta_m3l_sq_normal = 2.511e-3
delta_m3l_sq_inverted = -2.483e-3

# Legacy aliases retained for compatibility with older imports.
# - For NO, this equals Δm^2_32 = Δm^2_31 - Δm^2_21.
# - For IO, the historical name is misleading; the value matches m1^2 - m3^2,
#   which was the quantity consumed by the old compute_masses() implementation.
delta_m32_sq_NH = delta_m3l_sq_normal - delta_m21_sq
delta_m32_sq_IH = abs(delta_m3l_sq_inverted) - delta_m21_sq

# NuFIT 6.1 best-fit mixing angles. θ12 is common to both orderings in the
# published table, while θ23, θ13, and δ_CP are ordering-dependent.
sin_sq_theta12 = 0.3088
sin_sq_theta23_normal = 0.470
sin_sq_theta23_inverted = 0.550
sin_sq_theta13_normal = 0.02248
sin_sq_theta13_inverted = 0.02262

# Backward-compatible aliases default to normal ordering.
sin_sq_theta23 = sin_sq_theta23_normal
sin_sq_theta13 = sin_sq_theta13_normal

# Calculate angles in radians.
theta12 = np.arcsin(np.sqrt(sin_sq_theta12))
theta23_normal = np.arcsin(np.sqrt(sin_sq_theta23_normal))
theta23_inverted = np.arcsin(np.sqrt(sin_sq_theta23_inverted))
theta13_normal = np.arcsin(np.sqrt(sin_sq_theta13_normal))
theta13_inverted = np.arcsin(np.sqrt(sin_sq_theta13_inverted))

# Backward-compatible aliases default to normal ordering.
theta23 = theta23_normal
theta13 = theta13_normal

# CP-violating phase best-fit values from NuFIT 6.1 (degrees converted to rad).
delta_CP_normal = np.deg2rad(212.0)
delta_CP_inverted = np.deg2rad(274.0)

__all__ = [
   "NUFIT_DATASET",
   "NUFIT_REFERENCE",
   "delta_m21_sq",
   "delta_m3l_sq_normal",
   "delta_m3l_sq_inverted",
   "delta_m32_sq_NH",
   "delta_m32_sq_IH",
   "compute_masses",
   "theta12",
   "theta23_normal",
   "theta23_inverted",
   "theta23",
   "theta13_normal",
   "theta13_inverted",
   "theta13",
   "delta_CP_normal",
   "delta_CP_inverted",
   "sin_sq_theta12",
   "sin_sq_theta23_normal",
   "sin_sq_theta23_inverted",
   "sin_sq_theta23",
   "sin_sq_theta13_normal",
   "sin_sq_theta13_inverted",
   "sin_sq_theta13",
   "pmns_matrix",
   "get_pmns"
]


def compute_masses(
   lightest_mass: float, ordering: str = "normal"
) -> Tuple[float, float, float, float]:
   """Compute neutrino masses and their sum.

   Args:
      lightest_mass: Lightest neutrino mass (in eV). Must be non-negative.
      ordering: 'normal' or 'inverted' (case-insensitive).

   Returns:
      Tuple of (m1, m2, m3, M_nu) where M_nu = m1 + m2 + m3.

   Raises:
      ValueError: if `lightest_mass` is negative or ordering is invalid.
   """
   if lightest_mass < 0:
      raise ValueError("lightest_mass must be non-negative")

   ordering_key = ordering.strip().lower()
   if ordering_key == "normal":
      m1 = float(lightest_mass)
      m2 = float(np.sqrt(m1 * m1 + delta_m21_sq))
      m3 = float(np.sqrt(m1 * m1 + delta_m3l_sq_normal))
   elif ordering_key == "inverted":
      m3 = float(lightest_mass)
      m2 = float(np.sqrt(m3 * m3 + abs(delta_m3l_sq_inverted)))
      m1 = float(np.sqrt(m2 * m2 - delta_m21_sq))
   else:
      raise ValueError("ordering must be 'normal' or 'inverted'")

   M_nu = m1 + m2 + m3
   return m1, m2, m3, M_nu

def get_pmns(ordering: str = "normal", alpha: float = 0.0, beta: float = 0.0):
   """
   Returns the PMNS matrix using the configured mixing angles and the
   CP phase appropriate for the mass ordering.

   Args:
      ordering: 'normal' or 'inverted' (case-insensitive). Determines which
              delta_CP value to use.
      alpha, beta: Majorana phases (radians).

   Returns:
      3x3 complex numpy array: the PMNS matrix.
   """
   ordering_key = ordering.strip().lower()
   if ordering_key == "normal":
      theta23_local = theta23_normal
      theta13_local = theta13_normal
      delta_cp = delta_CP_normal
   elif ordering_key == "inverted":
      theta23_local = theta23_inverted
      theta13_local = theta13_inverted
      delta_cp = delta_CP_inverted
   else:
      raise ValueError("ordering must be 'normal' or 'inverted'")

   return pmns_matrix(theta12, theta23_local, theta13_local, delta_cp, alpha, beta)

def pmns_matrix(theta12, theta23, theta13, delta_cp, alpha=0.0, beta=0.0):
    """
    Construct the PMNS matrix (PDG parameterization) with optional Majorana phases.

    Args:
        theta12, theta23, theta13: mixing angles (radians)
        delta_cp: Dirac CP phase δ (radians)
        alpha, beta: Majorana phases (radians); implemented as U * diag(1, e^{i α/2}, e^{i β/2})

    Returns:
        U: 3x3 numpy array (complex)
    """
    s12 = np.sin(theta12)
    c12 = np.cos(theta12)
    s23 = np.sin(theta23)
    c23 = np.cos(theta23)
    s13 = np.sin(theta13)
    c13 = np.cos(theta13)

    # CP phase term
    # Note: In PDG, the term is e^{-i \delta} in the (0,2) element.
    # And e^{i \delta} appears in the other terms.
    e_minus_i_delta = np.exp(-1j * delta_cp)
    e_plus_i_delta = np.exp(1j * delta_cp)

    # Construct the matrix elements
    # Row 1
    u00 = c12 * c13
    u01 = s12 * c13
    u02 = s13 * e_minus_i_delta

    # Row 2
    u10 = -s12 * c23 - c12 * s23 * s13 * e_plus_i_delta
    u11 = c12 * c23 - s12 * s23 * s13 * e_plus_i_delta
    u12 = s23 * c13

    # Row 3
    u20 = s12 * s23 - c12 * c23 * s13 * e_plus_i_delta
    u21 = -c12 * s23 - s12 * c23 * s13 * e_plus_i_delta
    u22 = c23 * c13

    U = np.array([
        [u00, u01, u02],
        [u10, u11, u12],
        [u20, u21, u22]
    ], dtype=complex)

    # Apply Majorana phases on the right: diag(1, e^{i α/2}, e^{i β/2})
    # These are often denoted as alpha21 and alpha31
    maj_phases = np.diag([1.0 + 0j, np.exp(1j * alpha / 2.0), np.exp(1j * beta / 2.0)])

    return np.dot(U, maj_phases)
