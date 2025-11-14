"""massConstraints.py

Helpers to find allowed neutrino mass spectra given a sum constraint.

This module uses `compute_masses` from `massConfig.py` to compute
the three masses and their sum, then provides helpers to sweep the
lightest mass and return the allowed values under a maximum total
mass constraint (e.g. 0.082 eV).
"""

from __future__ import annotations

from typing import Tuple

import numpy as np

from massConfig import compute_masses, sum_mass_constraint


def find_allowed_lightest_masses(
    max_sum: float = sum_mass_constraint,
    ordering: str = "normal",
    step: float = 1e-4,
    max_lightest: float = 10.0,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Sweep the lightest neutrino mass and return allowed values.

    Args:
        max_sum: maximum allowed sum of neutrino masses (eV).
        ordering: 'normal' or 'inverted' (case-insensitive).
        step: increment for the lightest mass sweep (eV).
        max_lightest: upper bound for the sweep (eV). Sweep stops early
            once the total mass meets or exceeds `max_sum`.

        Returns:
            A tuple `(m1_values, m2_values, m3_values, Mtot_values)`
            where each element is a NumPy array.
            Values are included only while `Mtot < max_sum`.
    """
    if step <= 0:
        raise ValueError("step must be positive")
    if max_lightest <= 0:
        raise ValueError("max_lightest must be positive")
    
    m1_values = []
    m2_values = []
    m3_values = []
    Mtot_values = []

    # Use a simple loop to allow early stopping when constraint is hit.
    m = 0.0
    while m <= max_lightest:
        try:
            m1, m2, m3, Mtot = compute_masses(m, ordering=ordering)
        except ValueError:
            # propagate invalid ordering or negative mass errors
            raise

        if Mtot >= max_sum:
            break

        m1_values.append(m1)
        m2_values.append(m2)
        m3_values.append(m3)
        Mtot_values.append(Mtot)
        
        m += step

    return (
        np.array(m1_values),
        np.array(m2_values),
        np.array(m3_values),
        np.array(Mtot_values),
    )
