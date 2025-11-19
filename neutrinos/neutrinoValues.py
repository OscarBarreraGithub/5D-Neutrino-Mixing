"""neutrinoValues.py

This module provides the experimental constraints needed 
to compute neutrino masses.

Both normal and inverted mass orderings are supported.

"""

from __future__ import annotations

import numpy as np
from typing import Tuple
import cmath

# Constraint on the sum of neutrino masses from cosmology (in eV)
sum_mass_constraint = 0.082  # eV

# Central values of mass-squared differences (in eV^2)
delta_m21_sq = 7.53e-5   # Solar mass splitting
delta_m32_sq_NH = 2.455e-3   # Atmospheric splitting for Normal Hierarchy
delta_m32_sq_IH = 2.529e-3   # Atmospheric splitting for Inverted Hierarchy


# Mixing angles (from PDG 2024/2025 values)
# sin^2(theta) values
#Technically theta23 is for normal ordering, but the difference is negligible
sin_sq_theta12 = 0.307
sin_sq_theta23 = 0.534
sin_sq_theta13 = 0.0216

# Calculate angles in radians
theta12 = np.arcsin(np.sqrt(sin_sq_theta12))
theta23 = np.arcsin(np.sqrt(sin_sq_theta23))
theta13 = np.arcsin(np.sqrt(sin_sq_theta13))

# CP violating phase (Normal ordering value used here)
delta_CP_normal = 1.21 * np.pi  # radians
delta_CP_inverted = 1.58 * np.pi  # radians

__all__ = [
   "delta_m21_sq",
   "delta_m32_sq_NH",
   "delta_m32_sq_IH",
   "compute_masses",
   "theta12",
   "theta23",
   "theta13",
   "delta_CP_normal",
   "delta_CP_inverted",
   "sin_sq_theta12",
   "sin_sq_theta23",
   "sin_sq_theta13",
   "pmns_matrix",
   "get_pmns"
]


def compute_masses(lightest_mass: float, ordering: str = "normal") -> Tuple[float, float, float, float]:
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
      m3 = float(np.sqrt(m2 * m2 + delta_m32_sq_NH))
   elif ordering_key == "inverted":
      m3 = float(lightest_mass)
      m1 = float(np.sqrt(m3 * m3 + delta_m32_sq_IH))
      m2 = float(np.sqrt(m1 * m1 + delta_m21_sq))
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
      delta_cp = delta_CP_normal
   elif ordering_key == "inverted":
      delta_cp = delta_CP_inverted
   else:
      raise ValueError("ordering must be 'normal' or 'inverted'")

   return pmns_matrix(theta12, theta23, theta13, delta_cp, alpha, beta)

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

