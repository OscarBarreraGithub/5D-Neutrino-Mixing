"""Unified Yukawa computation API for the 5D RS neutrino mixing model.

This module provides the main entry point for computing all Yukawa couplings
(charged leptons and neutrinos) from the RS model input parameters.

Usage:
    from yukawa import compute_all_yukawas

    result = compute_all_yukawas(
        Lambda_IR=3000,           # 3 TeV KK scale
        c_L=0.58,                 # Lepton doublet bulk mass
        c_E=[0.75, 0.60, 0.50],   # RH charged lepton bulk masses
        c_N=0.27,                 # RH neutrino bulk mass
        M_N=1.22e18,              # M_Pl/10 Majorana mass
        lightest_nu_mass=0.002,   # 2 meV
        ordering='normal'
    )

    print(result.summary())
    print(f"Perturbative: {result.is_perturbative()}")
"""

from dataclasses import dataclass
from typing import Union, Dict, Tuple, Optional
import numpy as np

from .charged_lepton import compute_charged_lepton_yukawas
from .neutrino import compute_neutrino_yukawas


@dataclass
class YukawaResult:
    """Container for all Yukawa computation results.

    Attributes
    ----------
    Y_E : np.ndarray
        5D charged lepton Yukawa couplings (dimension [mass]^{-1}), shape (3,).
    Y_E_bar : np.ndarray
        Rescaled charged lepton Yukawas Ȳ_E = 2k·Y_E (dimensionless), shape (3,).
    Y_N : np.ndarray
        Neutrino Yukawa eigenvalues (dimension [mass]^{-1}), shape (3,).
    Y_N_bar : np.ndarray
        Rescaled neutrino Yukawas Ȳ_N = 2k·Y_N (dimensionless), shape (3,).
    Y_N_matrix : np.ndarray
        Full neutrino Yukawa matrix in charged-lepton mass basis (3×3).
        Computed as Y_N_matrix = V_PMNS · diag(Y_N).
    f_L : float
        IR overlap factor for lepton doublets.
    f_E : np.ndarray
        IR overlap factors for RH charged leptons, shape (3,).
    f_N : float
        IR overlap factor for RH neutrinos.
    f_N_UV : float
        UV overlap factor for RH neutrinos.
    epsilon : float
        Warp factor ε = Λ_IR/k.
    params : dict
        Input parameters used (for reproducibility).
    """
    # Charged leptons
    Y_E: np.ndarray
    Y_E_bar: np.ndarray

    # Neutrinos
    Y_N: np.ndarray
    Y_N_bar: np.ndarray
    Y_N_matrix: np.ndarray

    # Overlap factors
    f_L: float
    f_E: np.ndarray
    f_N: float
    f_N_UV: float

    # Geometry
    epsilon: float

    # Input parameters
    params: Dict

    def is_perturbative(self, max_Y_bar: float = 4.0) -> bool:
        """Check if all rescaled Yukawas satisfy perturbativity bound.

        Parameters
        ----------
        max_Y_bar : float, optional
            Maximum allowed value for |Ȳ|. Default 4.0.
            This corresponds to requiring at least ~3 KK modes before
            the Yukawa becomes strongly coupled.

        Returns
        -------
        bool
            True if all |Ȳ_E| < max_Y_bar and all |Ȳ_N| < max_Y_bar.
        """
        return (np.all(np.abs(self.Y_E_bar) < max_Y_bar) and
                np.all(np.abs(self.Y_N_bar) < max_Y_bar))

    def summary(self) -> str:
        """Return a formatted summary of the Yukawa computation.

        Returns
        -------
        str
            Multi-line summary string suitable for printing.
        """
        lines = [
            "=" * 60,
            "Yukawa Computation Results",
            "=" * 60,
            "",
            "Input Parameters:",
            f"  Λ_IR = {self.params['Lambda_IR']:.0f} GeV",
            f"  c_L = {self.params['c_L']:.4f}",
            f"  c_E = {self.params['c_E']}",
            f"  c_N = {self.params['c_N']:.4f}",
            f"  M_N = {self.params['M_N']:.2e} GeV",
            f"  lightest ν mass = {self.params['lightest_nu_mass']:.4f} eV",
            f"  ordering = {self.params['ordering']}",
            "",
            "Geometry:",
            f"  ε = Λ/k = {self.epsilon:.2e}",
            "",
            "Overlap Factors:",
            f"  f_L = {self.f_L:.6f}",
            f"  f_E = [{self.f_E[0]:.2e}, {self.f_E[1]:.4f}, {self.f_E[2]:.4f}]",
            f"  f_N = {self.f_N:.4f}",
            f"  f_N^UV = {self.f_N_UV:.6f}",
            "",
            "Charged Lepton Yukawas:",
            f"  Y_E (5D) = {self.Y_E}",
            f"  Ȳ_E (rescaled) = [{self.Y_E_bar[0]:.4f}, {self.Y_E_bar[1]:.4f}, {self.Y_E_bar[2]:.4f}]",
            "",
            "Neutrino Yukawas:",
            f"  Y_N (5D) = {self.Y_N}",
            f"  Ȳ_N (rescaled) = [{self.Y_N_bar[0]:.6f}, {self.Y_N_bar[1]:.6f}, {self.Y_N_bar[2]:.6f}]",
            "",
            f"Perturbative (|Ȳ| < 4): {self.is_perturbative()}",
            "=" * 60,
        ]
        return "\n".join(lines)


def compute_all_yukawas(
    Lambda_IR: float,
    c_L: float,
    c_E: Union[np.ndarray, Tuple[float, float, float]],
    c_N: float,
    M_N: float,
    lightest_nu_mass: float = 0.0,
    ordering: str = 'normal',
    majorana_alpha: float = 0.0,
    majorana_beta: float = 0.0,
    k: Optional[float] = None,
    v: float = 174.0
) -> YukawaResult:
    """
    Compute all Yukawa couplings from RS model parameters.

    This is the main entry point for the Yukawa computation module. It computes:
    - Charged lepton Yukawas by inverting the mass formula
    - Neutrino Yukawa eigenvalues by inverting the seesaw formula
    - Full neutrino Yukawa matrix including PMNS mixing

    Parameters
    ----------
    Lambda_IR : float
        IR/KK scale (GeV). Determines warp factor via ε = Λ_IR/k.
        Typical values: 1000-10000 GeV (1-10 TeV).
    c_L : float
        Bulk mass parameter for lepton doublets (universal).
        c > 0.5 localizes toward UV, c < 0.5 toward IR.
    c_E : array-like of shape (3,)
        Bulk mass parameters for RH charged leptons [c_e, c_μ, c_τ].
    c_N : float
        Bulk mass parameter for RH neutrinos (universal).
    M_N : float
        UV-localized Majorana mass scale (GeV).
        Typical values: 10^13 - 10^18 GeV.
    lightest_nu_mass : float, optional
        Lightest neutrino mass (eV). Default 0.0.
    ordering : str, optional
        Neutrino mass ordering: 'normal' or 'inverted'. Default 'normal'.
    majorana_alpha : float, optional
        Majorana CP phase α (radians). Default 0.0.
    majorana_beta : float, optional
        Majorana CP phase β (radians). Default 0.0.
    k : float, optional
        AdS curvature (GeV). Default is M_Planck ≈ 1.22e19 GeV.
    v : float, optional
        Electroweak VEV (GeV). Default 174 GeV.

    Returns
    -------
    YukawaResult
        Dataclass containing:
        - Y_E, Y_E_bar: Charged lepton Yukawas (5D and rescaled)
        - Y_N, Y_N_bar: Neutrino Yukawa eigenvalues (5D and rescaled)
        - Y_N_matrix: Full neutrino Yukawa matrix with PMNS
        - Overlap factors (f_L, f_E, f_N, f_N_UV)
        - Geometry parameters and input params

    Raises
    ------
    ValueError
        If parameters are invalid (e.g., negative masses, invalid ordering).

    Examples
    --------
    >>> # Benchmark point from Perez-Randall paper
    >>> result = compute_all_yukawas(
    ...     Lambda_IR=3000,
    ...     c_L=0.58,
    ...     c_E=[0.75, 0.60, 0.50],
    ...     c_N=0.27,
    ...     M_N=1.22e18,
    ...     lightest_nu_mass=0.002,
    ...     ordering='normal'
    ... )
    >>> print(f"f_L = {result.f_L:.4f}")  # Should be ~0.016
    >>> print(f"Ȳ_E = {result.Y_E_bar}")  # Should be O(1) to O(few)
    >>> print(f"Ȳ_N = {result.Y_N_bar}")  # Should be O(0.01) to O(0.1)
    """
    # Import dependencies (local import to avoid circular deps)
    from warpConfig.baseParams import get_warp_params, MPL, V_EWSB
    from warpConfig.wavefuncs import f_IR, f_UV
    from neutrinos.neutrinoValues import compute_masses, get_pmns

    # Set defaults
    if k is None:
        k = MPL

    # Compute geometry
    geom_params = get_warp_params(k=k, Lambda_IR=Lambda_IR)
    epsilon = geom_params['epsilon']

    # Compute overlap factors
    c_E_arr = np.asarray(c_E, dtype=float)
    f_L_val = float(f_IR(c_L, epsilon))
    f_E_vals = f_IR(c_E_arr, epsilon)
    f_N_val = float(f_IR(c_N, epsilon))
    f_N_UV_val = float(f_UV(c_N, epsilon))

    # Compute neutrino mass spectrum
    m1, m2, m3, M_sum = compute_masses(lightest_nu_mass, ordering)
    nu_masses_eV = (m1, m2, m3)

    # Compute charged lepton Yukawas
    cl_result = compute_charged_lepton_yukawas(
        f_L=f_L_val,
        f_E=f_E_vals,
        k=k,
        v=v
    )

    # Compute neutrino Yukawas
    nu_result = compute_neutrino_yukawas(
        f_L=f_L_val,
        f_N=f_N_val,
        f_N_UV=f_N_UV_val,
        M_N=M_N,
        neutrino_masses_eV=nu_masses_eV,
        k=k,
        v=v
    )

    # Construct full neutrino Yukawa matrix with PMNS
    # In the charged-lepton mass basis: Y_N → V_PMNS · diag(Y_N)
    V_pmns = get_pmns(ordering, majorana_alpha, majorana_beta)
    Y_N_diag = np.diag(nu_result['Y_N'])
    Y_N_matrix = V_pmns @ Y_N_diag

    # Assemble result
    return YukawaResult(
        Y_E=cl_result['Y_E'],
        Y_E_bar=cl_result['Y_E_bar'],
        Y_N=nu_result['Y_N'],
        Y_N_bar=nu_result['Y_N_bar'],
        Y_N_matrix=Y_N_matrix,
        f_L=f_L_val,
        f_E=f_E_vals,
        f_N=f_N_val,
        f_N_UV=f_N_UV_val,
        epsilon=epsilon,
        params={
            'Lambda_IR': Lambda_IR,
            'c_L': c_L,
            'c_E': c_E_arr.tolist(),
            'c_N': c_N,
            'M_N': M_N,
            'lightest_nu_mass': lightest_nu_mass,
            'ordering': ordering,
            'majorana_alpha': majorana_alpha,
            'majorana_beta': majorana_beta,
            'k': k,
            'v': v,
        }
    )
