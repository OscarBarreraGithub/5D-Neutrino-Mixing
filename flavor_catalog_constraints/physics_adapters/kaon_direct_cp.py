"""Stub adapter for direct CP violation in kaon ``K -> pi pi`` decays.

This adapter deliberately does **not** compute a Standard-Model or RS
penguin amplitude for

    Re(epsilon'/epsilon).

The SM prediction depends on a cancellation between QCD and electroweak
penguins, the Delta I = 1/2 enhancement, and lattice K -> pi pi matrix
elements.  A grounded RS contribution would require Delta S = 1 penguin and
chromomagnetic matching, RG evolution, and hadronic matrix elements that are
not available on ``ParameterPoint``.  The only honest v1 operation is a small
finite-number bookkeeping comparison between the measured value and the
explicitly reported non-vetoing "NP room".
"""

from __future__ import annotations

from dataclasses import dataclass
import math

__all__ = [
    "KAON_DIRECT_CP_STUB_MODEL_V1",
    "KAON_DIRECT_CP_SM_NEEDS_HUMAN_PHYSICS",
    "KAON_DIRECT_CP_RS_NEEDS_HUMAN_PHYSICS",
    "KaonDirectCPRoomComparison",
    "compare_epsilon_prime_np_room_to_measurement",
]


KAON_DIRECT_CP_STUB_MODEL_V1 = "epsilon_prime_over_epsilon_no_penguin_matching_stub_v1"

KAON_DIRECT_CP_SM_NEEDS_HUMAN_PHYSICS = (
    "NEEDS-HUMAN-PHYSICS: SM Re(epsilon'/epsilon) in K -> pi pi requires "
    "the QCD/electroweak penguin cancellation, Delta I=1/2 enhancement, "
    "isospin-breaking corrections, and lattice K -> pi pi matrix elements; "
    "no veto-grade first-principles SM prediction is implemented here."
)

KAON_DIRECT_CP_RS_NEEDS_HUMAN_PHYSICS = (
    "NEEDS-HUMAN-PHYSICS: rigorous RS Re(epsilon'/epsilon) needs Delta S=1 "
    "penguin/chromomagnetic matching, RG evolution, and K -> pi pi hadronic "
    "matrix elements not carried on ParameterPoint; no RS penguin amplitude "
    "is computed."
)


@dataclass(frozen=True)
class KaonDirectCPRoomComparison:
    """Scalar comparison between the measured value and non-vetoing NP room.

    ``measured_re_epsilon_prime_over_epsilon`` and ``documented_np_room_abs``
    are dimensionless.  ``ratio_to_room`` is
    ``abs(measured_re_epsilon_prime_over_epsilon) / documented_np_room_abs``.
    This is not an SM or RS prediction; it is only explicit bookkeeping used
    by the K003 INFO stub.
    """

    measured_re_epsilon_prime_over_epsilon: float
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


def compare_epsilon_prime_np_room_to_measurement(
    *,
    measured_re_epsilon_prime_over_epsilon: float,
    experimental_uncertainty: float,
    documented_np_room_abs: float,
) -> KaonDirectCPRoomComparison:
    """Compare documented NP room with the measured K003 anchor.

    This function performs no penguin matching or hadronic calculation.  It
    exists so the K003 constraint can return real numeric ``budget``/``ratio``
    fields while making the missing SM and RS physics explicit.
    """

    measured = _finite_float(
        measured_re_epsilon_prime_over_epsilon,
        name="measured_re_epsilon_prime_over_epsilon",
    )
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
    return KaonDirectCPRoomComparison(
        measured_re_epsilon_prime_over_epsilon=float(measured),
        experimental_uncertainty=float(uncertainty),
        measurement_abs=float(measurement_abs),
        documented_np_room_abs=float(room),
        ratio_to_room=float(ratio),
        passes=bool(ratio <= 1.0),
    )
