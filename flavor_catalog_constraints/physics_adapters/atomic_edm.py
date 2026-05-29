"""Stub adapter for Ra-225 and Xe-129 atomic electric dipole moments.

This adapter deliberately does **not** compute atomic EDMs from RS model
parameters.  A grounded prediction for

    d_Ra, d_Xe

needs two missing ingredients:

1. isotope-specific nuclear Schiff moments and relativistic atomic-structure
   factors translating CP-odd low-energy operators into atomic EDMs; and
2. RS matching onto the relevant CP-odd quark, gluon, and semileptonic
   operators.

The only honest v1 operations are therefore finite-number validation of a
catalogued direct limit and an explicit scalar comparison helper for tests or
future diagnostics.  No Schiff-moment, atomic-structure, or RS CP-odd matching
calculation is performed here.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any

__all__ = [
    "ATOMIC_EDM_STUB_MODEL_V1",
    "ATOMIC_EDM_SCHIFF_ATOMIC_NEEDS_HUMAN_PHYSICS",
    "ATOMIC_EDM_RS_MATCHING_NEEDS_HUMAN_PHYSICS",
    "AtomicEDMDirectLimit",
    "AtomicEDMLimitComparison",
    "atomic_edm_direct_limit",
    "compare_atomic_edm_to_limit",
]


ATOMIC_EDM_STUB_MODEL_V1 = (
    "ra225_xe129_atomic_edm_no_schiff_atomic_response_no_rs_cp_odd_matching_stub_v1"
)

ATOMIC_EDM_SCHIFF_ATOMIC_NEEDS_HUMAN_PHYSICS = (
    "NEEDS-HUMAN-PHYSICS: Ra-225/Xe-129 atomic EDM interpretation requires "
    "isotope-specific nuclear Schiff moments plus relativistic atomic-structure "
    "factors; these nuclear/atomic inputs are not fixed by this catalog."
)

ATOMIC_EDM_RS_MATCHING_NEEDS_HUMAN_PHYSICS = (
    "NEEDS-HUMAN-PHYSICS: rigorous RS matching onto CP-odd quark, gluon, "
    "and semileptonic low-energy operators for Ra-225/Xe-129 is not available "
    "on ParameterPoint; no RS atomic-EDM amplitude is computed."
)


@dataclass(frozen=True)
class AtomicEDMDirectLimit:
    """Validated direct atomic-EDM upper limit in ``e cm``."""

    isotope: str
    observable: str
    experimental_limit_e_cm: float


@dataclass(frozen=True)
class AtomicEDMLimitComparison:
    """Scalar comparison between an atomic EDM value and a direct limit.

    All EDM quantities are in ``e cm``.  ``ratio_to_limit`` is
    ``abs(atomic_edm_e_cm) / experimental_limit_e_cm``.  This helper performs
    no RS matching and no nuclear/atomic response calculation.
    """

    isotope: str
    atomic_edm_e_cm: float
    experimental_limit_e_cm: float
    abs_atomic_edm_e_cm: float
    ratio_to_limit: float
    passes: bool


def _finite_float(value: Any, *, name: str) -> float:
    out = float(value)
    if not math.isfinite(out):
        raise ValueError(f"{name} must be finite, got {value!r}")
    return out


def _positive_float(value: Any, *, name: str) -> float:
    out = _finite_float(value, name=name)
    if out <= 0.0:
        raise ValueError(f"{name} must be positive, got {value!r}")
    return out


def atomic_edm_direct_limit(
    *,
    isotope: str,
    observable: str,
    experimental_limit_e_cm: float,
) -> AtomicEDMDirectLimit:
    """Validate and record a direct atomic-EDM limit.

    The returned value is a catalog reference, not a model prediction.
    """

    return AtomicEDMDirectLimit(
        isotope=str(isotope),
        observable=str(observable),
        experimental_limit_e_cm=_positive_float(
            experimental_limit_e_cm,
            name="experimental_limit_e_cm",
        ),
    )


def compare_atomic_edm_to_limit(
    *,
    atomic_edm_e_cm: float,
    experimental_limit_e_cm: float,
    isotope: str = "atomic",
) -> AtomicEDMLimitComparison:
    """Compare an explicit atomic EDM value with a direct upper limit."""

    value = _finite_float(atomic_edm_e_cm, name="atomic_edm_e_cm")
    limit = _positive_float(
        experimental_limit_e_cm,
        name="experimental_limit_e_cm",
    )
    abs_value = abs(value)
    ratio = abs_value / limit
    return AtomicEDMLimitComparison(
        isotope=str(isotope),
        atomic_edm_e_cm=float(value),
        experimental_limit_e_cm=float(limit),
        abs_atomic_edm_e_cm=float(abs_value),
        ratio_to_limit=float(ratio),
        passes=bool(ratio <= 1.0),
    )
