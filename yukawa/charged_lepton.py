"""Charged lepton Yukawa computation.

Computes the diagonal Yukawa couplings for charged leptons (e, μ, τ) by inverting
the mass relation for an IR-localized Higgs in warped geometry.

Physics:
    The charged lepton mass in RS with IR-brane Higgs is:
        m_{E_i} = 2 v k · f_L · Y_{E_i} · f_{E_i}

    where:
        v = 174 GeV (electroweak VEV)
        k = AdS curvature (~M_Planck)
        f_L = IR overlap factor for lepton doublets
        f_{E_i} = IR overlap factor for RH charged lepton of generation i
        Y_{E_i} = 5D Yukawa coupling (dimension [mass]^{-1})

    Inverting:
        Y_{E_i} = m_{E_i} / (2 v k · f_L · f_{E_i})

    The rescaled (dimensionless) Yukawa is:
        Ȳ_{E_i} = 2k · Y_{E_i} = m_{E_i} / (v · f_L · f_{E_i})

    For naturalness, we want |Ȳ| ~ O(1), with perturbativity requiring |Ȳ| < 4.

Reference:
    Perez & Randall, arXiv:0805.4652, Eq. (5)
"""

from typing import Union, Dict, Tuple
import numpy as np

from .constants import LEPTON_MASSES


def compute_charged_lepton_yukawas(
    f_L: float,
    f_E: Union[np.ndarray, Tuple[float, float, float]],
    k: float,
    v: float = 174.0,
    lepton_masses: Tuple[float, float, float] = LEPTON_MASSES
) -> Dict[str, np.ndarray]:
    """
    Compute charged lepton Yukawa couplings from the IR-brane mass formula.

    Parameters
    ----------
    f_L : float
        IR overlap factor for lepton doublets (universal).
        Computed as f_IR(c_L, epsilon) from warpConfig.wavefuncs.
    f_E : array-like of shape (3,)
        IR overlap factors for RH charged leptons [f_e, f_μ, f_τ].
        Computed as f_IR(c_{E_i}, epsilon) for each generation.
    k : float
        AdS curvature scale (GeV). Typically ~M_Planck ≈ 1.22e19 GeV.
    v : float, optional
        Electroweak VEV (GeV). Default 174 GeV.
    lepton_masses : tuple of 3 floats, optional
        Physical charged lepton masses (m_e, m_μ, m_τ) in GeV.
        Default uses PDG 2024 values from constants.py.

    Returns
    -------
    dict with keys:
        'Y_E' : np.ndarray of shape (3,)
            5D Yukawa couplings Y_{E_i} (dimension [mass]^{-1}).
        'Y_E_bar' : np.ndarray of shape (3,)
            Rescaled (dimensionless) Yukawas Ȳ_{E_i} = 2k · Y_{E_i}.
        'masses_GeV' : np.ndarray of shape (3,)
            Input lepton masses used (for verification).

    Raises
    ------
    ValueError
        If any f-factor is zero or negative (would cause division by zero).

    Examples
    --------
    >>> from warpConfig.baseParams import get_warp_params, MPL
    >>> from warpConfig.wavefuncs import f_IR
    >>> params = get_warp_params(Lambda_IR=3000)
    >>> f_L = f_IR(0.58, params['epsilon'])
    >>> f_E = f_IR(np.array([0.75, 0.60, 0.50]), params['epsilon'])
    >>> result = compute_charged_lepton_yukawas(f_L, f_E, MPL)
    >>> print(result['Y_E_bar'])  # Should be O(1) to O(few)
    """
    # Convert to numpy arrays
    f_E_arr = np.asarray(f_E, dtype=float)
    m_E = np.asarray(lepton_masses, dtype=float)

    # Validate inputs
    if f_L <= 0:
        raise ValueError(f"f_L must be positive, got {f_L}")
    if np.any(f_E_arr <= 0):
        raise ValueError(f"All f_E values must be positive, got {f_E_arr}")
    if f_E_arr.shape != (3,):
        raise ValueError(f"f_E must have shape (3,), got {f_E_arr.shape}")

    # Compute 5D Yukawa couplings
    # Y_{E_i} = m_{E_i} / (2 v k · f_L · f_{E_i})
    Y_E = m_E / (2.0 * v * k * f_L * f_E_arr)

    # Compute rescaled (dimensionless) Yukawas
    # Ȳ_{E_i} = 2k · Y_{E_i} = m_{E_i} / (v · f_L · f_{E_i})
    Y_E_bar = 2.0 * k * Y_E  # Equivalently: m_E / (v * f_L * f_E_arr)

    return {
        'Y_E': Y_E,
        'Y_E_bar': Y_E_bar,
        'masses_GeV': m_E,
    }
