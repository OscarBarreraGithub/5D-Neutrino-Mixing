"""wavefuncs.py

This module provides functions to compute the IR and UV overlap factors
for 5D fermion zero modes in a warped geometry.
"""

from typing import Union
import numpy as np

def f_IR(c: Union[float, np.ndarray], epsilon: float) -> Union[float, np.ndarray]:
    """
    Computes the IR overlap factor f_IR(c).

    Formula:
        f_IR^2 = (0.5 - c) / (1 - epsilon**(1 - 2*c))

    Parameters
    ----------
    c : float or np.ndarray
        The bulk mass parameter (c = M/k).
    epsilon : float
        The warp factor epsilon = exp(-pi k rc).

    Returns
    -------
    float or np.ndarray
        The value of f_IR(c).
    """
    c_arr = np.asarray(c, dtype=float)
    res_sq = np.zeros_like(c_arr)
    
    # Handle the singularity at c = 0.5
    mask = ~np.isclose(c_arr, 0.5)
    
    # For c != 0.5
    if np.any(mask):
        c_safe = c_arr[mask]
        # Note: epsilon**(1-2c) can be very large or very small.
        # If c > 0.5, 1-2c < 0, exponent is negative, term is large (epsilon is small).
        # If c < 0.5, 1-2c > 0, exponent is positive, term is small.
        denom = 1.0 - epsilon**(1.0 - 2.0*c_safe)
        res_sq[mask] = (0.5 - c_safe) / denom

    # For c == 0.5, take the limit
    # Limit_{c->0.5} f^2 = 1 / (-2 * ln(epsilon))
    if np.any(~mask):
        res_sq[~mask] = 1.0 / (-2.0 * np.log(epsilon))
        
    # Ensure non-negative before sqrt (numerical noise might cause issues if result is effectively 0)
    res_sq = np.maximum(res_sq, 0.0)
    
    return np.sqrt(res_sq)


def f_UV(c: Union[float, np.ndarray], epsilon: float) -> Union[float, np.ndarray]:
    """
    Computes the UV overlap factor f_UV(c).

    Formula:
        f_UV^2 = (0.5 - c) / (epsilon**(2*c - 1) - 1)

    Parameters
    ----------
    c : float or np.ndarray
        The bulk mass parameter (c = M/k).
    epsilon : float
        The warp factor epsilon = exp(-pi k rc).

    Returns
    -------
    float or np.ndarray
        The value of f_UV(c).
    """
    c_arr = np.asarray(c, dtype=float)
    res_sq = np.zeros_like(c_arr)
    
    # Handle the singularity at c = 0.5
    mask = ~np.isclose(c_arr, 0.5)
    
    # For c != 0.5
    if np.any(mask):
        c_safe = c_arr[mask]
        denom = epsilon**(2.0*c_safe - 1.0) - 1.0
        res_sq[mask] = (0.5 - c_safe) / denom

    # For c == 0.5, take the limit
    # Limit_{c->0.5} f^2 = 1 / (-2 * ln(epsilon))
    if np.any(~mask):
        res_sq[~mask] = 1.0 / (-2.0 * np.log(epsilon))
        
    # Ensure non-negative
    res_sq = np.maximum(res_sq, 0.0)

    return np.sqrt(res_sq)

