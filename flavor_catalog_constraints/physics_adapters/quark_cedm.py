"""Stub adapter for quark chromo-electric dipole moment bounds.

This adapter deliberately does **not** compute quark cEDMs from RS model
parameters, and it does **not** translate quark cEDMs into neutron or atomic
EDMs.  A grounded E008 prediction needs two missing ingredients:

1. RS CP-odd chromo-dipole matching onto low-energy quark cEDM operators.
2. Hadronic and nuclear matrix elements mapping those operators into
   neutron and diamagnetic-atom EDM observables.

The only honest v1 operations are therefore finite-number validation of the
catalogued reference bounds and an explicit caller-supplied comparison helper
for tests or future diagnostics.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

__all__ = [
    "QUARK_CEDM_STUB_MODEL_V1",
    "QUARK_CEDM_HADRONIC_NEEDS_HUMAN_PHYSICS",
    "QUARK_CEDM_RS_MATCHING_NEEDS_HUMAN_PHYSICS",
    "QuarkCEDMReferenceBound",
    "QuarkCEDMReferenceComparison",
    "quark_cedm_reference_bound",
    "compare_quark_cedm_to_reference_bound",
]


QUARK_CEDM_STUB_MODEL_V1 = (
    "quark_cedm_no_rs_cp_odd_matching_no_hadronic_matrix_elements_stub_v1"
)

QUARK_CEDM_HADRONIC_NEEDS_HUMAN_PHYSICS = (
    "NEEDS-HUMAN-PHYSICS: quark cEDM bounds from neutron and atomic EDMs "
    "require hadronic and nuclear matrix elements, CP-odd source choices, "
    "and operator-basis conventions that are not fixed by this catalog."
)

QUARK_CEDM_RS_MATCHING_NEEDS_HUMAN_PHYSICS = (
    "NEEDS-HUMAN-PHYSICS: rigorous RS CP-odd chromo-dipole matching onto "
    "quark cEDM Wilson coefficients is not available on ParameterPoint; "
    "no RS qCEDM amplitude is computed."
)


@dataclass(frozen=True)
class QuarkCEDMReferenceBound:
    """Validated reference bound on a real quark-cEDM combination.

    ``reference_bound_cm`` is the absolute upper bound on the named cEDM
    combination in cm.  It is a catalog reference value, not a prediction.
    """

    bound_observable: str
    reference_bound_cm: float


@dataclass(frozen=True)
class QuarkCEDMReferenceComparison:
    """Scalar comparison to a reference qCEDM bound.

    This helper compares an explicit caller-supplied cEDM combination to a
    catalogued reference bound.  It performs no RS matching and no hadronic or
    nuclear EDM calculation.
    """

    bound_observable: str
    qcedm_combination_cm: float
    reference_bound_cm: float
    abs_qcedm_combination_cm: float
    ratio_to_bound: float
    passes: bool


def _finite_float(value: Any, *, name: str) -> float:
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"{name} must be finite, got {value!r}")
    return number


def _positive_float(value: Any, *, name: str) -> float:
    number = _finite_float(value, name=name)
    if number <= 0.0:
        raise ValueError(f"{name} must be positive, got {value!r}")
    return number


def quark_cedm_reference_bound(
    *,
    bound_observable: str,
    reference_bound_cm: float,
) -> QuarkCEDMReferenceBound:
    """Validate and record a catalogued qCEDM reference bound."""

    return QuarkCEDMReferenceBound(
        bound_observable=str(bound_observable),
        reference_bound_cm=_positive_float(
            reference_bound_cm,
            name="reference_bound_cm",
        ),
    )


def compare_quark_cedm_to_reference_bound(
    *,
    qcedm_combination_cm: float,
    reference_bound_cm: float,
    bound_observable: str = "quark cEDM combination",
) -> QuarkCEDMReferenceComparison:
    """Compare an explicit qCEDM value with a catalogued reference bound.

    The input ``qcedm_combination_cm`` must already be a low-energy qCEDM
    combination in cm.  This function intentionally does not infer it from RS
    parameters or from neutron/atomic EDMs.
    """

    value = _finite_float(
        qcedm_combination_cm,
        name="qcedm_combination_cm",
    )
    bound = _positive_float(
        reference_bound_cm,
        name="reference_bound_cm",
    )
    abs_value = abs(value)
    ratio = abs_value / bound
    return QuarkCEDMReferenceComparison(
        bound_observable=str(bound_observable),
        qcedm_combination_cm=float(value),
        reference_bound_cm=float(bound),
        abs_qcedm_combination_cm=float(abs_value),
        ratio_to_bound=float(ratio),
        passes=bool(ratio <= 1.0),
    )
