"""PDG 2024 MS-bar quark mass inputs and RG-evolved targets.

This module owns the PDG 2024 "Quark Masses" review numbers (central values,
symmetric 1-sigma uncertainties, the reference scale at which each mass is
quoted, and the active-flavor count at that reference). It also exposes
helper utilities that RG-evolve those values to a single common
renormalization scale (typically ``mu_common = qcd.constants.M_TOP_MS =
162.5 GeV``) so the quark fitter can score against a *consistent* MS-bar
target table rather than an ad-hoc mix of scales.

The Wilson-coefficient ``alpha_s`` reference scale used for Delta-F=2
matching is **independent** of these mass-target scales. That scale stays
at ``mu = 3 TeV`` (`MODERN_DEFAULT_TARGET_SCALE_GEV`); see
``quarkConstraints/modern/inputs.py``. The two scales are kept orthogonal
on purpose: changing the mass target scale must not change WC running.

References
----------
PDG 2024, "Quark Masses" review.
Chetyrkin, Kuhn, Sturm, Eur.Phys.J. C48 (2006) 107 — threshold matching.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional

import numpy as np

from qcd.mass_running import run_msbar_mass


@dataclass(frozen=True)
class PDGQuarkMass:
    """One PDG 2024 MS-bar quark mass input.

    Attributes
    ----------
    central_GeV
        PDG 2024 central value, in GeV.
    sigma_GeV
        Symmetric one-sigma uncertainty in GeV (averaged from PDG asymmetric
        bounds where applicable).
    mu_ref_GeV
        Reference scale at which ``central_GeV`` is quoted (GeV).
    n_f_at_reference
        Number of active MS-bar flavors at ``mu_ref_GeV`` (PDG convention).
    flavor
        Quark flavor identifier ('u', 'd', 's', 'c', 'b', 't').
    """

    flavor: str
    central_GeV: float
    sigma_GeV: float
    mu_ref_GeV: float
    n_f_at_reference: int

    def __post_init__(self) -> None:
        if self.central_GeV <= 0.0:
            raise ValueError("central_GeV must be positive")
        if self.sigma_GeV <= 0.0:
            raise ValueError("sigma_GeV must be positive")
        if self.mu_ref_GeV <= 0.0:
            raise ValueError("mu_ref_GeV must be positive")
        if self.n_f_at_reference not in (3, 4, 5, 6):
            raise ValueError("n_f_at_reference must be in {3,4,5,6}")
        if self.flavor not in ("u", "d", "s", "c", "b", "t"):
            raise ValueError("flavor must be one of u,d,s,c,b,t")


# PDG 2024 "Quark Masses" review.
# Light quarks (u, d, s) are quoted at mu = 2 GeV in n_f = 4 (PDG convention).
# Heavy quarks (c, b, t) are quoted at their own scales m_q(m_q).
PDG_2024_QUARK_MASSES: Dict[str, PDGQuarkMass] = {
    "u": PDGQuarkMass(
        flavor="u",
        central_GeV=2.16e-3,
        sigma_GeV=0.49e-3,
        mu_ref_GeV=2.0,
        n_f_at_reference=4,
    ),
    "d": PDGQuarkMass(
        flavor="d",
        central_GeV=4.70e-3,
        sigma_GeV=0.43e-3,
        mu_ref_GeV=2.0,
        n_f_at_reference=4,
    ),
    "s": PDGQuarkMass(
        flavor="s",
        central_GeV=93.5e-3,
        sigma_GeV=0.8e-3,
        mu_ref_GeV=2.0,
        n_f_at_reference=4,
    ),
    "c": PDGQuarkMass(
        flavor="c",
        central_GeV=1.273,
        sigma_GeV=0.0046,
        mu_ref_GeV=1.273,
        n_f_at_reference=4,
    ),
    "b": PDGQuarkMass(
        flavor="b",
        central_GeV=4.183,
        sigma_GeV=0.007,
        mu_ref_GeV=4.183,
        n_f_at_reference=5,
    ),
    "t": PDGQuarkMass(
        flavor="t",
        central_GeV=162.5,
        # PDG quotes +2.1/-1.5 GeV; the scalar schema stores the
        # symmetrized one-sigma width.
        sigma_GeV=1.8,
        mu_ref_GeV=162.5,
        n_f_at_reference=6,
    ),
}

PDG_QUARK_MASSES_EDITION = "PDG 2024"


def _flavor_order(flavors: tuple[str, ...] | None) -> tuple[str, ...]:
    if flavors is None:
        return ("u", "d", "s", "c", "b", "t")
    out = tuple(flavors)
    for f in out:
        if f not in PDG_2024_QUARK_MASSES:
            raise KeyError(f"unknown flavor {f!r}")
    return out


def pdg_quark_masses_at_scale(
    mu_GeV: float,
    *,
    flavors: tuple[str, ...] | None = None,
    n_loops: int = 4,
) -> Dict[str, float]:
    """Return PDG 2024 quark masses RG-evolved to ``mu_GeV``.

    Parameters
    ----------
    mu_GeV
        Common target scale (GeV) for the returned MS-bar masses. Must be
        strictly above the charm mass — no path runs below ``m_c``.
    flavors
        Subset of flavors to return. Default returns all six.
    n_loops
        Loop order for the alpha_s and mass-anomalous-dimension running
        (1..4). Default 4-loop.

    Returns
    -------
    dict[str, float]
        Mapping ``flavor -> m_q(mu_GeV)`` in GeV.
    """
    if mu_GeV <= 0.0:
        raise ValueError("mu_GeV must be positive")
    out: Dict[str, float] = {}
    for f in _flavor_order(flavors):
        entry = PDG_2024_QUARK_MASSES[f]
        out[f] = run_msbar_mass(
            m_ref=entry.central_GeV,
            mu_ref=entry.mu_ref_GeV,
            mu_target=mu_GeV,
            n_f_ref=entry.n_f_at_reference,
            n_loops=n_loops,
        )
    return out


def pdg_2sigma_relative_at_scale(
    mu_GeV: float,
    *,
    flavors: tuple[str, ...] | None = None,
) -> Dict[str, float]:
    """Return per-flavor PDG 2sigma relative uncertainty at the reference.

    Note: PDG quotes the relative uncertainty at the reference scale, not at
    ``mu_GeV``. The relative uncertainty is invariant under one-loop
    multiplicative running to leading order, so we keep the relative number
    fixed; the absolute uncertainty at ``mu_GeV`` follows by multiplying by
    ``pdg_quark_masses_at_scale(mu_GeV)[f]``.

    Returns
    -------
    dict[str, float]
        Mapping ``flavor -> 2 * sigma_GeV / central_GeV``.
    """
    out: Dict[str, float] = {}
    for f in _flavor_order(flavors):
        entry = PDG_2024_QUARK_MASSES[f]
        out[f] = 2.0 * entry.sigma_GeV / entry.central_GeV
    return out


def pdg_2sigma_absolute_at_scale(
    mu_GeV: float,
    *,
    flavors: tuple[str, ...] | None = None,
    n_loops: int = 4,
) -> Dict[str, float]:
    """Return per-flavor PDG 2sigma absolute uncertainty at ``mu_GeV``.

    Computed as ``relative_2sigma * m_q(mu_GeV)``.
    """
    rel = pdg_2sigma_relative_at_scale(mu_GeV, flavors=flavors)
    central = pdg_quark_masses_at_scale(mu_GeV, flavors=flavors, n_loops=n_loops)
    return {f: rel[f] * central[f] for f in rel}


def pdg_up_down_arrays_at_scale(
    mu_GeV: float,
    *,
    n_loops: int = 4,
) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(up_masses, down_masses)`` numpy arrays at ``mu_GeV``.

    Order matches the existing repo convention:
        up   = [m_u, m_c, m_t]
        down = [m_d, m_s, m_b]
    """
    masses = pdg_quark_masses_at_scale(mu_GeV, n_loops=n_loops)
    up = np.array([masses["u"], masses["c"], masses["t"]], dtype=float)
    down = np.array([masses["d"], masses["s"], masses["b"]], dtype=float)
    return up, down
