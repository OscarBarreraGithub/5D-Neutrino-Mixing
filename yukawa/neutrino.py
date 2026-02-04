"""Neutrino Yukawa computation via seesaw inversion.

Computes the neutrino Yukawa eigenvalues by inverting the Type-I seesaw formula
for warped geometry with UV-localized Majorana mass.

Physics:
    In the universal limit (c_L, c_N, M_N generation-independent), the light
    neutrino mass from the seesaw is:

        m_{ν_i} = (2 k² v² f²_L f²_N) / ((f^UV_N)² M_N) · Y²_{N_i}

    where:
        k = AdS curvature (~M_Planck)
        v = 174 GeV (electroweak VEV)
        f_L = f_IR(c_L, ε) for lepton doublets
        f_N = f_IR(c_N, ε) for RH neutrinos at IR brane
        f^UV_N = f_UV(c_N, ε) for RH neutrinos at UV brane
        M_N = UV-localized Majorana mass (GeV)
        Y_{N_i} = neutrino Yukawa eigenvalues (dimension [mass]^{-1})

    Inverting:
        Y_{N_i} = sqrt[ m_{ν_i} · (f^UV_N)² · M_N / (2 k² v² f²_L f²_N) ]

    The rescaled (dimensionless) Yukawa is:
        Ȳ_{N_i} = 2k · Y_{N_i}

    For perturbativity, we require |Ȳ| < 4.

Reference:
    Perez & Randall, arXiv:0805.4652, Eq. (6)
"""

from typing import Union, Dict, Tuple
import numpy as np

from .constants import EV_TO_GEV


def compute_neutrino_yukawas(
    f_L: float,
    f_N: float,
    f_N_UV: float,
    M_N: float,
    neutrino_masses_eV: Union[np.ndarray, Tuple[float, float, float]],
    k: float,
    v: float = 174.0
) -> Dict[str, np.ndarray]:
    """
    Compute neutrino Yukawa eigenvalues from the seesaw formula.

    This inverts the Type-I seesaw formula in the universal limit where
    c_L, c_N, and M_N are generation-independent.

    Parameters
    ----------
    f_L : float
        IR overlap factor for lepton doublets (universal).
    f_N : float
        IR overlap factor for RH neutrinos (universal).
        Computed as f_IR(c_N, epsilon).
    f_N_UV : float
        UV overlap factor for RH neutrinos.
        Computed as f_UV(c_N, epsilon).
        This enters the canonical normalization of the Majorana mass.
    M_N : float
        UV-localized Majorana mass scale (GeV).
        Typically M_N ~ M_Pl/10 or similar.
    neutrino_masses_eV : array-like of shape (3,)
        Light neutrino masses (m_1, m_2, m_3) in eV.
        Use neutrinoValues.compute_masses() to get these from oscillation data.
    k : float
        AdS curvature scale (GeV). Typically ~M_Planck.
    v : float, optional
        Electroweak VEV (GeV). Default 174 GeV.

    Returns
    -------
    dict with keys:
        'Y_N' : np.ndarray of shape (3,)
            Neutrino Yukawa eigenvalues (dimension [mass]^{-1}).
        'Y_N_bar' : np.ndarray of shape (3,)
            Rescaled (dimensionless) Yukawas Ȳ_{N_i} = 2k · Y_{N_i}.
        'masses_eV' : np.ndarray of shape (3,)
            Input neutrino masses (for verification).

    Raises
    ------
    ValueError
        If any f-factor or M_N is zero or negative.

    Notes
    -----
    - Neutrino masses are input in eV (standard convention) and converted
      internally to GeV for consistency with other quantities.
    - If any neutrino mass is zero (e.g., lightest neutrino massless),
      the corresponding Yukawa is zero.
    - The formula assumes the universal limit. For non-universal bulk masses,
      the full matrix structure would need to be considered.

    Examples
    --------
    >>> from warpConfig.baseParams import get_warp_params, MPL
    >>> from warpConfig.wavefuncs import f_IR, f_UV
    >>> from neutrinos.neutrinoValues import compute_masses
    >>> params = get_warp_params(Lambda_IR=3000)
    >>> epsilon = params['epsilon']
    >>> f_L = f_IR(0.58, epsilon)
    >>> f_N = f_IR(0.27, epsilon)
    >>> f_N_UV = f_UV(0.27, epsilon)
    >>> m1, m2, m3, _ = compute_masses(0.002, 'normal')
    >>> result = compute_neutrino_yukawas(f_L, f_N, f_N_UV, 1.22e18, (m1, m2, m3), MPL)
    >>> print(result['Y_N_bar'])  # Should be O(0.2) to O(1.0)
    """
    # Convert to numpy array and eV to GeV
    m_nu_eV = np.asarray(neutrino_masses_eV, dtype=float)
    m_nu_GeV = m_nu_eV * EV_TO_GEV

    # Validate inputs
    if f_L <= 0:
        raise ValueError(f"f_L must be positive, got {f_L}")
    if f_N <= 0:
        raise ValueError(f"f_N must be positive, got {f_N}")
    if f_N_UV <= 0:
        raise ValueError(f"f_N_UV must be positive, got {f_N_UV}")
    if M_N <= 0:
        raise ValueError(f"M_N must be positive, got {M_N}")
    if m_nu_eV.shape != (3,):
        raise ValueError(f"neutrino_masses_eV must have shape (3,), got {m_nu_eV.shape}")
    if np.any(m_nu_eV < 0):
        raise ValueError(f"Neutrino masses cannot be negative, got {m_nu_eV}")

    # Compute the seesaw prefactor
    # From: m_ν = (2 k² v² f²_L f²_N) / ((f^UV_N)² M_N) · Y²_N
    # Solve for Y²_N: Y²_N = m_ν · (f^UV_N)² · M_N / (2 k² v² f²_L f²_N)
    prefactor = (f_N_UV**2 * M_N) / (2.0 * k**2 * v**2 * f_L**2 * f_N**2)

    # Compute Y_N (handle zero masses gracefully)
    Y_N_sq = m_nu_GeV * prefactor
    Y_N = np.sqrt(Y_N_sq)

    # Compute rescaled (dimensionless) Yukawas
    Y_N_bar = 2.0 * k * Y_N

    return {
        'Y_N': Y_N,
        'Y_N_bar': Y_N_bar,
        'masses_eV': m_nu_eV,
    }
