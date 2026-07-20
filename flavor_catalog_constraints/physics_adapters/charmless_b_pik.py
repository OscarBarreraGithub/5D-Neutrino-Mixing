"""Stub adapter for charmless nonleptonic ``B -> pi K`` observables.

This adapter deliberately does **not** compute Standard-Model or RS
amplitudes for the four ``B -> pi K`` modes.  The SM interpretation requires
QCDF/SCET or phenomenological hadronic amplitudes with strong phases, and a
rigorous RS contribution requires Delta B = 1 QCD/electroweak penguin
matching, running, and nonleptonic matrix elements that are not present on
``ParameterPoint``.

The only honest v1 operation is a finite-number bookkeeping comparison between
a measured anchor and an explicitly documented "NP room" reported by the
constraint.  It is not a prediction and must remain non-vetoing.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

__all__ = [
    "B_TO_PI_K_NONLEPTONIC_STUB_MODEL_V1",
    "B_TO_PI_K_SM_NEEDS_HUMAN_PHYSICS",
    "B_TO_PI_K_RS_NEEDS_HUMAN_PHYSICS",
    "CharmlessBPiKRoomComparison",
    "compare_b_to_pi_k_np_room_to_measurement",
]


B_TO_PI_K_NONLEPTONIC_STUB_MODEL_V1 = (
    "b_to_pi_k_nonleptonic_no_hadronic_amplitude_stub_v1"
)

B_TO_PI_K_SM_NEEDS_HUMAN_PHYSICS = (
    "NEEDS-HUMAN-PHYSICS: SM B -> pi K amplitudes need QCDF/SCET or a "
    "phenomenological nonleptonic amplitude fit with tree, QCD-penguin, "
    "electroweak-penguin, and strong-phase inputs; no reliable first-principles "
    "SM central value is available in this catalog."
)

B_TO_PI_K_RS_NEEDS_HUMAN_PHYSICS = (
    "NEEDS-HUMAN-PHYSICS: rigorous RS B -> pi K requires Delta B=1 QCD and "
    "electroweak penguin matching, RG evolution, and nonleptonic hadronic "
    "matrix elements; no RS penguin amplitude is computed."
)


@dataclass(frozen=True)
class CharmlessBPiKRoomComparison:
    """Scalar comparison between a measured ``B -> pi K`` anchor and NP room.

    All values are dimensionless.  ``ratio_to_room`` is
    ``abs(measured_observable) / documented_np_room_abs``.  This is only the
    explicit bookkeeping used by the B032 INFO stub; it is not a hadronic
    amplitude calculation.
    """

    measured_observable: float
    experimental_uncertainty: float
    measurement_abs: float
    documented_np_room_abs: float
    ratio_to_room: float
    passes: bool


def _finite_float(value: float, *, name: str) -> float:
    out = float(value)
    if not math.isfinite(out):
        raise ValueError(f"{name} must be finite, got {value!r}")
    return out


def compare_b_to_pi_k_np_room_to_measurement(
    *,
    measured_observable: float,
    experimental_uncertainty: float,
    documented_np_room_abs: float,
) -> CharmlessBPiKRoomComparison:
    """Compare a documented NP room with a measured ``B -> pi K`` anchor.

    This function performs no QCDF/SCET, penguin matching, RG evolution, or
    nonleptonic matrix-element calculation.  It exists so the constraint can
    return real scalar ``budget``/``ratio`` fields while keeping the missing
    SM and RS physics visible.
    """

    measured = _finite_float(measured_observable, name="measured_observable")
    uncertainty = _finite_float(
        experimental_uncertainty,
        name="experimental_uncertainty",
    )
    room = _finite_float(
        documented_np_room_abs,
        name="documented_np_room_abs",
    )
    if uncertainty <= 0.0:
        raise ValueError("experimental_uncertainty must be positive")
    if room <= 0.0:
        raise ValueError("documented_np_room_abs must be positive")

    measurement_abs = abs(measured)
    ratio = measurement_abs / room
    return CharmlessBPiKRoomComparison(
        measured_observable=float(measured),
        experimental_uncertainty=float(uncertainty),
        measurement_abs=float(measurement_abs),
        documented_np_room_abs=float(room),
        ratio_to_room=float(ratio),
        passes=bool(ratio <= 1.0),
    )
