"""Stub adapter for long-distance radiative kaon branching fractions.

This adapter deliberately does **not** compute a Standard-Model ChPT
amplitude or an RS new-physics amplitude for

    K_L -> pi0 gamma gamma.

The measured rate is dominated by chiral-loop and vector-meson-exchange
hadronic dynamics.  A grounded RS prediction would require Delta S = 1
radiative matching and matching onto the chiral amplitude/counterterms, none
of which is available on ``ParameterPoint``.  The only honest v1 operation is
therefore finite-number bookkeeping for the measured branching fraction and
the explicitly reported non-vetoing "NP room".
"""

from __future__ import annotations

import math
from dataclasses import dataclass

__all__ = [
    "RADIATIVE_KAON_STUB_MODEL_V1",
    "RADIATIVE_KAON_SM_NEEDS_HUMAN_PHYSICS",
    "RADIATIVE_KAON_RS_NEEDS_HUMAN_PHYSICS",
    "RadiativeKaonRoomComparison",
    "compare_kl_pi0gammagamma_np_room_to_measurement",
]


RADIATIVE_KAON_STUB_MODEL_V1 = "kl_pi0_gammagamma_no_chpt_matching_stub_v1"

RADIATIVE_KAON_SM_NEEDS_HUMAN_PHYSICS = (
    "NEEDS-HUMAN-PHYSICS: SM K_L -> pi0 gamma gamma is dominated by ChPT "
    "O(p^4)/O(p^6) chiral loops and vector-meson exchange; no dedicated "
    "ChPT amplitude calculation is implemented in this catalog."
)

RADIATIVE_KAON_RS_NEEDS_HUMAN_PHYSICS = (
    "NEEDS-HUMAN-PHYSICS: rigorous RS K_L -> pi0 gamma gamma needs Delta S=1 "
    "radiative matching, dipole/electroweak/two-photon operator treatment, "
    "and matching onto chiral counterterms not carried on ParameterPoint; no "
    "RS amplitude is computed."
)


@dataclass(frozen=True)
class RadiativeKaonRoomComparison:
    """Scalar comparison between the measured BR and non-vetoing NP room.

    ``measured_branching_fraction`` and ``documented_np_room_abs`` are
    dimensionless branching fractions.  ``ratio_to_room`` is
    ``abs(measured_branching_fraction) / documented_np_room_abs``.  This is
    not a ChPT or RS prediction; it is only explicit bookkeeping used by the
    K013 INFO stub.
    """

    measured_branching_fraction: float
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


def compare_kl_pi0gammagamma_np_room_to_measurement(
    *,
    measured_branching_fraction: float,
    experimental_uncertainty: float,
    documented_np_room_abs: float,
) -> RadiativeKaonRoomComparison:
    """Compare a documented NP room with the measured K013 BR anchor.

    This function performs no chiral or new-physics calculation.  It exists so
    the K013 constraint can return real numeric ``budget``/``ratio`` fields
    while making the missing SM and RS physics explicit.
    """

    measured = _finite_float(
        measured_branching_fraction,
        name="measured_branching_fraction",
    )
    uncertainty = _finite_float(
        experimental_uncertainty,
        name="experimental_uncertainty",
    )
    room = _finite_float(
        documented_np_room_abs,
        name="documented_np_room_abs",
    )
    if measured < 0.0:
        raise ValueError("measured_branching_fraction must be non-negative")
    if uncertainty <= 0.0:
        raise ValueError("experimental_uncertainty must be positive")
    if room <= 0.0:
        raise ValueError("documented_np_room_abs must be positive")

    measurement_abs = abs(measured)
    ratio = measurement_abs / room
    return RadiativeKaonRoomComparison(
        measured_branching_fraction=float(measured),
        experimental_uncertainty=float(uncertainty),
        measurement_abs=float(measurement_abs),
        documented_np_room_abs=float(room),
        ratio_to_room=float(ratio),
        passes=bool(ratio <= 1.0),
    )
