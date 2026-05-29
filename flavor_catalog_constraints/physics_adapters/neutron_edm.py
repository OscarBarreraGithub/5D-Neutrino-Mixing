"""Stub adapter for the neutron electric dipole moment.

This adapter deliberately does **not** compute a neutron EDM from RS model
parameters.  A grounded prediction for

    d_n

needs two missing ingredients:

1. CP-odd quark-level matching in the RS model, including quark EDMs, quark
   chromo-EDMs, and possibly the Weinberg three-gluon operator.
2. Non-perturbative neutron matrix elements translating those low-energy
   operators into ``d_n``.  Lattice and QCD sum-rule inputs carry large,
   basis-dependent uncertainties and are not part of the current
   ``ParameterPoint`` contract.

The only honest v1 operation is therefore a finite-number bookkeeping
comparison between the measured central value and the experimental upper
limit loaded by the constraint from the catalog sidecar.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

__all__ = [
    "NEUTRON_EDM_STUB_MODEL_V1",
    "NEUTRON_EDM_HADRONIC_NEEDS_HUMAN_PHYSICS",
    "NEUTRON_EDM_RS_MATCHING_NEEDS_HUMAN_PHYSICS",
    "NeutronEDMLimitComparison",
    "compare_neutron_edm_measurement_to_limit",
]


NEUTRON_EDM_STUB_MODEL_V1 = (
    "neutron_edm_no_quark_dipole_matching_no_hadronic_matrix_elements_stub_v1"
)

NEUTRON_EDM_HADRONIC_NEEDS_HUMAN_PHYSICS = (
    "NEEDS-HUMAN-PHYSICS: neutron EDM interpretation requires "
    "non-perturbative matrix elements for quark EDMs, quark chromo-EDMs, "
    "and CP-odd gluonic operators; these lattice/QCD-sum-rule inputs are "
    "not fixed by this catalog."
)

NEUTRON_EDM_RS_MATCHING_NEEDS_HUMAN_PHYSICS = (
    "NEEDS-HUMAN-PHYSICS: rigorous RS matching onto CP-odd quark EDM, "
    "quark chromo-EDM, and Weinberg-operator coefficients is not available "
    "on ParameterPoint; no RS neutron-EDM amplitude is computed."
)


@dataclass(frozen=True)
class NeutronEDMLimitComparison:
    """Scalar comparison between a measured central value and an EDM limit.

    All EDM quantities are in ``e cm``.  ``ratio_to_limit`` is
    ``abs(measured_neutron_edm_e_cm) / experimental_limit_e_cm``.  This is
    not an RS prediction; it is only the explicit bookkeeping used by the
    E004 INFO stub.
    """

    measured_neutron_edm_e_cm: float
    experimental_limit_e_cm: float
    measurement_abs_e_cm: float
    ratio_to_limit: float
    passes: bool


def _finite_float(value: float, *, name: str) -> float:
    out = float(value)
    if not math.isfinite(out):
        raise ValueError(f"{name} must be finite, got {value!r}")
    return out


def compare_neutron_edm_measurement_to_limit(
    *,
    measured_neutron_edm_e_cm: float,
    experimental_limit_e_cm: float,
) -> NeutronEDMLimitComparison:
    """Compare the measured neutron EDM central value with the upper limit.

    The function performs no quark-dipole matching and no hadronic matrix
    element calculation.  It exists so the constraint can expose real numeric
    ``budget`` and ``ratio`` fields while keeping the missing physics visible.
    """

    measured = _finite_float(
        measured_neutron_edm_e_cm,
        name="measured_neutron_edm_e_cm",
    )
    limit = _finite_float(
        experimental_limit_e_cm,
        name="experimental_limit_e_cm",
    )
    if limit <= 0.0:
        raise ValueError("experimental_limit_e_cm must be positive")

    measurement_abs = abs(measured)
    ratio = measurement_abs / limit
    return NeutronEDMLimitComparison(
        measured_neutron_edm_e_cm=float(measured),
        experimental_limit_e_cm=float(limit),
        measurement_abs_e_cm=float(measurement_abs),
        ratio_to_limit=float(ratio),
        passes=bool(ratio <= 1.0),
    )
