"""massConfig.py

This module provides the experimental constraints needed 
to compute neutrino masses and their 1σ uncertainties.

Both normal and inverted mass orderings are supported.

Example:
   from massConfig import compute_masses
   m1, m2, m3, M_nu = compute_masses(1e-3, ordering='normal')
"""

from __future__ import annotations

import numpy as np
from typing import Tuple

# Constraint on the sum of neutrino masses from cosmology (in eV)
sum_mass_constraint = 0.082  # eV

# Central values of mass-squared differences (in eV^2)
delta_m21_sq = 7.53e-5   # Solar mass splitting
delta_m32_sq_NH = 2.455e-3   # Atmospheric splitting for Normal Hierarchy
delta_m32_sq_IH = 2.529e-3   # Atmospheric splitting for Inverted Hierarchy

# Experimental uncertainties (1σ) in eV^2
sigma_m21_sq = 0
sigma_m32_sq = 0

__all__ = [
   "delta_m21_sq",
   "delta_m32_sq_NH",
   "delta_m32_sq_IH",
   "sigma_m21_sq",
   "sigma_m32_sq",
   "compute_masses",
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
