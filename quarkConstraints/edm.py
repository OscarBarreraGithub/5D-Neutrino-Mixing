"""Charged-lepton electric dipole moment helpers.

This module intentionally implements only the observable-side EDM conversion:

    d_l / e = c_l^CP-odd

where ``c_l^CP-odd`` has mass dimension -1 in GeV units.  The conversion to
catalog units is therefore

    d_l [e cm] = c_l^CP-odd [GeV^-1] * hbar c [GeV cm].

NEEDS-HUMAN-PHYSICS
-------------------
The RS contribution to a charged-lepton EDM is a one-loop CP-violating dipole
matching problem involving complex lepton/KK-fermion/Higgs couplings.  Those
couplings are not represented by the current ``ParameterPoint``.  This module
does not invent that matching; callers must supply an explicit low-energy
CP-odd dipole proxy coefficient, and catalog constraints must flag that proxy
path.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Mapping

__all__ = [
    "HBARC_GEV_CM",
    "EDM_COEFFICIENT_CONVENTION_V1",
    "EDM_RS_MATCHING_GAP_V1",
    "ChargedLeptonEDMResult",
    "cp_odd_projection_from_chiral_dipole",
    "edm_e_cm_from_cp_odd_dipole",
    "evaluate_charged_lepton_edm",
    "evaluate_charged_lepton_edm_from_chiral_dipole",
]


# CODATA hbar*c in GeV cm: 1 GeV^-1 = 1.973269804e-14 cm.
HBARC_GEV_CM = 1.973269804e-14

EDM_COEFFICIENT_CONVENTION_V1 = (
    "Low-energy EDM proxy convention: c_CPodd has units GeV^-1 and is defined "
    "by d_l/e = c_CPodd, so d_l[e cm] = c_CPodd[GeV^-1] * hbar*c[GeV cm]."
)

EDM_RS_MATCHING_GAP_V1 = (
    "NEEDS-HUMAN-PHYSICS: full RS charged-lepton EDM matching requires "
    "complex lepton, KK-fermion, Higgs, and electroweak couplings that are "
    "not available on ParameterPoint; this module only converts an explicit "
    "low-energy CP-odd dipole proxy coefficient to |d_l|."
)


@dataclass(frozen=True)
class ChargedLeptonEDMResult:
    """Pure-NP charged-lepton EDM prediction in catalog units."""

    lepton: str
    edm_e_cm: float
    abs_edm_e_cm: float
    sm_edm_e_cm: float
    experimental_limit_e_cm: float
    passes: bool
    ratio_to_limit: float
    cp_odd_dipole_coefficient_gev_inv: float
    hbarc_gev_cm: float
    coefficient_source: str
    diagnostics: Mapping[str, Any]


def _finite_float(value: Any, *, name: str) -> float:
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"{name} must be finite")
    return number


def _positive_finite(value: Any, *, name: str) -> float:
    number = _finite_float(value, name=name)
    if number <= 0.0:
        raise ValueError(f"{name} must be positive")
    return number


def _finite_complex(value: Any, *, name: str) -> complex:
    number = complex(value)
    if not math.isfinite(number.real) or not math.isfinite(number.imag):
        raise ValueError(f"{name} must be finite")
    return number


def cp_odd_projection_from_chiral_dipole(
    chiral_dipole_coefficient_gev_inv: complex,
) -> float:
    """Return the CP-odd projection of a complex chiral EDM coefficient.

    In the common chiral-dipole convention the EDM is proportional to the
    imaginary part of the coefficient.  This helper exposes that projection
    without committing to any RS loop matching model.
    """

    coefficient = _finite_complex(
        chiral_dipole_coefficient_gev_inv,
        name="chiral_dipole_coefficient_gev_inv",
    )
    return float(coefficient.imag)


def edm_e_cm_from_cp_odd_dipole(cp_odd_dipole_coefficient_gev_inv: float) -> float:
    """Convert a CP-odd dipole coefficient in GeV^-1 to an EDM in e cm."""

    coefficient = _finite_float(
        cp_odd_dipole_coefficient_gev_inv,
        name="cp_odd_dipole_coefficient_gev_inv",
    )
    return float(coefficient * HBARC_GEV_CM)


def evaluate_charged_lepton_edm(
    cp_odd_dipole_coefficient_gev_inv: float,
    *,
    experimental_limit_e_cm: float,
    lepton: str,
    sm_edm_e_cm: float = 0.0,
    coefficient_source: str = "explicit low-energy CP-odd dipole proxy",
    diagnostics: Mapping[str, Any] | None = None,
) -> ChargedLeptonEDMResult:
    """Evaluate a pure-NP charged-lepton EDM against an experimental limit."""

    limit = _positive_finite(
        experimental_limit_e_cm,
        name="experimental_limit_e_cm",
    )
    sm = _finite_float(sm_edm_e_cm, name="sm_edm_e_cm")
    coefficient = _finite_float(
        cp_odd_dipole_coefficient_gev_inv,
        name="cp_odd_dipole_coefficient_gev_inv",
    )
    edm = edm_e_cm_from_cp_odd_dipole(coefficient)
    abs_edm = abs(edm)
    ratio = float(abs_edm / limit)
    merged_diagnostics = {
        "coefficient_convention": EDM_COEFFICIENT_CONVENTION_V1,
        "rs_matching_gap": EDM_RS_MATCHING_GAP_V1,
        "sm_edm_policy": (
            "SM charged-lepton EDM is negligible at catalog precision; the "
            "HARD bound is applied to the pure NP contribution."
        ),
        **dict(diagnostics or {}),
    }
    return ChargedLeptonEDMResult(
        lepton=str(lepton),
        edm_e_cm=float(edm),
        abs_edm_e_cm=float(abs_edm),
        sm_edm_e_cm=float(sm),
        experimental_limit_e_cm=float(limit),
        passes=bool(ratio <= 1.0),
        ratio_to_limit=ratio,
        cp_odd_dipole_coefficient_gev_inv=float(coefficient),
        hbarc_gev_cm=float(HBARC_GEV_CM),
        coefficient_source=str(coefficient_source),
        diagnostics=merged_diagnostics,
    )


def evaluate_charged_lepton_edm_from_chiral_dipole(
    chiral_dipole_coefficient_gev_inv: complex,
    *,
    experimental_limit_e_cm: float,
    lepton: str,
    sm_edm_e_cm: float = 0.0,
    coefficient_source: str = "explicit complex chiral dipole proxy",
    diagnostics: Mapping[str, Any] | None = None,
) -> ChargedLeptonEDMResult:
    """Evaluate an EDM using ``Im(C_chiral)`` as the CP-odd coefficient."""

    chiral = _finite_complex(
        chiral_dipole_coefficient_gev_inv,
        name="chiral_dipole_coefficient_gev_inv",
    )
    cp_odd = cp_odd_projection_from_chiral_dipole(chiral)
    merged_diagnostics = {
        "chiral_dipole_coefficient_gev_inv": complex(chiral),
        "cp_odd_projection": "imaginary_part",
        **dict(diagnostics or {}),
    }
    return evaluate_charged_lepton_edm(
        cp_odd,
        experimental_limit_e_cm=experimental_limit_e_cm,
        lepton=lepton,
        sm_edm_e_cm=sm_edm_e_cm,
        coefficient_source=coefficient_source,
        diagnostics=merged_diagnostics,
    )
