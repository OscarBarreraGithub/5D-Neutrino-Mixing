"""Yukawa computation module for the 5D RS neutrino mixing model.

This module provides functions to compute Yukawa couplings from RS model parameters.
The main entry point is `compute_all_yukawas()`.

Example
-------
>>> from yukawa import compute_all_yukawas
>>>
>>> result = compute_all_yukawas(
...     Lambda_IR=3000,           # 3 TeV KK scale
...     c_L=0.58,                 # Lepton doublet bulk mass
...     c_E=[0.75, 0.60, 0.50],   # RH charged lepton bulk masses
...     c_N=0.27,                 # RH neutrino bulk mass
...     M_N=1.22e18,              # M_Pl/10 Majorana mass
...     lightest_nu_mass=0.002,   # 2 meV
...     ordering='normal'
... )
>>>
>>> print(result.summary())
>>> print(f"Perturbative: {result.is_perturbative()}")
"""

from .compute_yukawas import compute_all_yukawas, YukawaResult
from .charged_lepton import compute_charged_lepton_yukawas
from .neutrino import compute_neutrino_yukawas
from .constants import (
    M_ELECTRON,
    M_MUON,
    M_TAU,
    LEPTON_MASSES,
    EV_TO_GEV,
)

__all__ = [
    # Main API
    'compute_all_yukawas',
    'YukawaResult',
    # Individual computation functions
    'compute_charged_lepton_yukawas',
    'compute_neutrino_yukawas',
    # Constants
    'M_ELECTRON',
    'M_MUON',
    'M_TAU',
    'LEPTON_MASSES',
    'EV_TO_GEV',
]
