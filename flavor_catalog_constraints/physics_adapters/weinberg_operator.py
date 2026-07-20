"""Stub adapter for the Weinberg CP-odd three-gluon operator.

This adapter deliberately does **not** compute a Weinberg-operator Wilson
coefficient from RS model parameters and does **not** translate that
coefficient into a neutron EDM.  A grounded E009 prediction needs two missing
ingredients:

1. RS CP-odd gluonic matching onto a chosen Weinberg-operator basis and
   threshold/RG convention.
2. Non-perturbative hadronic matrix elements for the neutron EDM response to
   the CP-odd gluonic operator.  Current benchmark translations are highly
   convention-dependent and are not fixed by ``ParameterPoint``.

The only honest v1 operations are therefore finite-number validation of a
catalogued single-source benchmark bound and an explicit caller-supplied
comparison helper for tests or future diagnostics.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

__all__ = [
    "WEINBERG_OPERATOR_STUB_MODEL_V1",
    "WEINBERG_OPERATOR_HADRONIC_NEEDS_HUMAN_PHYSICS",
    "WEINBERG_OPERATOR_RS_GLUONIC_MATCHING_NEEDS_HUMAN_PHYSICS",
    "WeinbergOperatorReferenceBound",
    "WeinbergOperatorReferenceComparison",
    "weinberg_operator_reference_bound",
    "compare_weinberg_coefficient_to_reference_bound",
]


WEINBERG_OPERATOR_STUB_MODEL_V1 = (
    "weinberg_operator_no_rs_cp_odd_gluonic_matching_no_hadronic_matrix_"
    "elements_stub_v1"
)

WEINBERG_OPERATOR_HADRONIC_NEEDS_HUMAN_PHYSICS = (
    "NEEDS-HUMAN-PHYSICS: Weinberg-operator neutron-EDM bounds require "
    "non-perturbative CP-odd gluonic hadronic matrix elements, operator "
    "mixing conventions, and source/cancellation assumptions that are not "
    "fixed by this catalog."
)

WEINBERG_OPERATOR_RS_GLUONIC_MATCHING_NEEDS_HUMAN_PHYSICS = (
    "NEEDS-HUMAN-PHYSICS: rigorous RS matching onto the CP-odd gluonic "
    "Weinberg operator is not available on ParameterPoint; no RS "
    "three-gluon Wilson coefficient is computed."
)


@dataclass(frozen=True)
class WeinbergOperatorReferenceBound:
    """Validated reference bound on a real Weinberg coefficient.

    ``reference_bound_gev_minus2`` is the absolute upper bound in the stated
    operator convention.  It is a catalog reference value, not a prediction.
    """

    bound_observable: str
    reference_bound_gev_minus2: float
    convention: str


@dataclass(frozen=True)
class WeinbergOperatorReferenceComparison:
    """Scalar comparison to a Weinberg reference bound.

    This helper compares an explicit caller-supplied coefficient to a
    catalogued benchmark bound.  It performs no RS matching and no hadronic
    neutron-EDM calculation.
    """

    bound_observable: str
    coefficient_gev_minus2: float
    reference_bound_gev_minus2: float
    abs_coefficient_gev_minus2: float
    ratio_to_bound: float
    passes: bool
    convention: str


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


def weinberg_operator_reference_bound(
    *,
    bound_observable: str,
    reference_bound_gev_minus2: float,
    convention: str,
) -> WeinbergOperatorReferenceBound:
    """Validate and record a catalogued Weinberg-operator reference bound."""

    return WeinbergOperatorReferenceBound(
        bound_observable=str(bound_observable),
        reference_bound_gev_minus2=_positive_float(
            reference_bound_gev_minus2,
            name="reference_bound_gev_minus2",
        ),
        convention=str(convention),
    )


def compare_weinberg_coefficient_to_reference_bound(
    *,
    coefficient_gev_minus2: float,
    reference_bound_gev_minus2: float,
    bound_observable: str = "Weinberg three-gluon coefficient",
    convention: str = "caller-supplied Weinberg convention",
) -> WeinbergOperatorReferenceComparison:
    """Compare an explicit Weinberg coefficient with a reference bound.

    The input ``coefficient_gev_minus2`` must already be a low-energy
    coefficient in the same convention as the reference bound.  This function
    intentionally does not infer it from RS parameters or neutron EDM data.
    """

    value = _finite_float(
        coefficient_gev_minus2,
        name="coefficient_gev_minus2",
    )
    bound = _positive_float(
        reference_bound_gev_minus2,
        name="reference_bound_gev_minus2",
    )
    abs_value = abs(value)
    ratio = abs_value / bound
    return WeinbergOperatorReferenceComparison(
        bound_observable=str(bound_observable),
        coefficient_gev_minus2=float(value),
        reference_bound_gev_minus2=float(bound),
        abs_coefficient_gev_minus2=float(abs_value),
        ratio_to_bound=float(ratio),
        passes=bool(ratio <= 1.0),
        convention=str(convention),
    )
