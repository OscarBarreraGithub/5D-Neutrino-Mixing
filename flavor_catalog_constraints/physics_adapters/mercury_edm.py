"""Stub adapter for the Mercury-199 atomic electric dipole moment.

This adapter deliberately does **not** compute a ``199Hg`` EDM from RS model
parameters.  A grounded prediction for

    d_Hg

needs two missing ingredients:

1. nuclear Schiff-moment response and relativistic atomic-structure factors
   translating CP-odd low-energy operators into the atomic EDM; and
2. RS matching onto the relevant CP-odd quark, gluon, and semileptonic
   operators.

The only honest v1 operation is therefore a finite-number bookkeeping
comparison between the measured central value and the experimental upper
limit loaded by the constraint from the catalog sidecar.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

__all__ = [
    "MERCURY_EDM_STUB_MODEL_V1",
    "MERCURY_EDM_SCHIFF_ATOMIC_NEEDS_HUMAN_PHYSICS",
    "MERCURY_EDM_RS_MATCHING_NEEDS_HUMAN_PHYSICS",
    "MercuryEDMLimitComparison",
    "compare_mercury_edm_measurement_to_limit",
]


MERCURY_EDM_STUB_MODEL_V1 = (
    "mercury_199_edm_no_schiff_atomic_response_no_rs_cp_odd_matching_stub_v1"
)

MERCURY_EDM_SCHIFF_ATOMIC_NEEDS_HUMAN_PHYSICS = (
    "NEEDS-HUMAN-PHYSICS: 199Hg atomic EDM interpretation requires a "
    "nuclear Schiff moment calculation plus relativistic atomic-structure "
    "factors; these nonperturbative nuclear/atomic inputs are not fixed by "
    "this catalog."
)

MERCURY_EDM_RS_MATCHING_NEEDS_HUMAN_PHYSICS = (
    "NEEDS-HUMAN-PHYSICS: rigorous RS matching onto CP-odd quark, gluon, "
    "and semileptonic low-energy operators for 199Hg is not available on "
    "ParameterPoint; no RS mercury-EDM amplitude is computed."
)


@dataclass(frozen=True)
class MercuryEDMLimitComparison:
    """Scalar comparison between a measured central value and a Hg EDM limit.

    All EDM quantities are in ``e cm``.  ``ratio_to_limit`` is
    ``abs(measured_mercury_edm_e_cm) / experimental_limit_e_cm``.  This is
    not an RS prediction; it is only the explicit bookkeeping used by the
    E006 INFO stub.
    """

    measured_mercury_edm_e_cm: float
    experimental_limit_e_cm: float
    measurement_abs_e_cm: float
    ratio_to_limit: float
    passes: bool


def _finite_float(value: float, *, name: str) -> float:
    out = float(value)
    if not math.isfinite(out):
        raise ValueError(f"{name} must be finite, got {value!r}")
    return out


def compare_mercury_edm_measurement_to_limit(
    *,
    measured_mercury_edm_e_cm: float,
    experimental_limit_e_cm: float,
) -> MercuryEDMLimitComparison:
    """Compare the measured Hg EDM central value with the upper limit.

    The function performs no Schiff-moment calculation, no atomic-structure
    calculation, and no RS CP-odd matching.  It exists so the constraint can
    expose real numeric ``budget`` and ``ratio`` fields while keeping the
    missing physics visible.
    """

    measured = _finite_float(
        measured_mercury_edm_e_cm,
        name="measured_mercury_edm_e_cm",
    )
    limit = _finite_float(
        experimental_limit_e_cm,
        name="experimental_limit_e_cm",
    )
    if limit <= 0.0:
        raise ValueError("experimental_limit_e_cm must be positive")

    measurement_abs = abs(measured)
    ratio = measurement_abs / limit
    return MercuryEDMLimitComparison(
        measured_mercury_edm_e_cm=float(measured),
        experimental_limit_e_cm=float(limit),
        measurement_abs_e_cm=float(measurement_abs),
        ratio_to_limit=float(ratio),
        passes=bool(ratio <= 1.0),
    )
