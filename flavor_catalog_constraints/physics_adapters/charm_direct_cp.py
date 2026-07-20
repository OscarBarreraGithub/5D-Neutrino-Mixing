"""Stub adapter for direct CP violation in singly Cabibbo-suppressed charm.

This adapter deliberately does **not** compute a Standard-Model or RS
penguin amplitude for

    Delta A_CP(D0 -> K+ K-, pi+ pi-).

The observable is dominated by long-distance, non-perturbative penguin
matrix elements in the SM, and the catalog ``ParameterPoint`` does not carry
the Delta C = 1 RS penguin matching inputs needed for a grounded NP
prediction.  The only honest v1 operation is therefore a finite-number,
documented comparison between the measured asymmetry and the explicitly
chosen "NP room" reported by the constraint.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

__all__ = [
    "CHARM_DIRECT_CP_STUB_MODEL_V1",
    "CHARM_DIRECT_CP_SM_NEEDS_HUMAN_PHYSICS",
    "CHARM_DIRECT_CP_RS_NEEDS_HUMAN_PHYSICS",
    "CharmDirectCPRoomComparison",
    "compare_delta_acp_np_room_to_measurement",
]


CHARM_DIRECT_CP_STUB_MODEL_V1 = "charm_direct_cp_no_penguin_matching_stub_v1"

CHARM_DIRECT_CP_SM_NEEDS_HUMAN_PHYSICS = (
    "NEEDS-HUMAN-PHYSICS: SM Delta A_CP in D0 -> K+K-, pi+pi- is dominated "
    "by long-distance non-perturbative penguin amplitudes; no reliable "
    "first-principles SM central value is available in this catalog."
)

CHARM_DIRECT_CP_RS_NEEDS_HUMAN_PHYSICS = (
    "NEEDS-HUMAN-PHYSICS: rigorous RS Delta C=1 penguin matching for "
    "D0 -> K+K-, pi+pi- is not available on ParameterPoint; no RS penguin "
    "amplitude is computed."
)


@dataclass(frozen=True)
class CharmDirectCPRoomComparison:
    """Scalar comparison between the measured asymmetry and NP room.

    ``measured_delta_acp`` and ``documented_np_room_abs`` are dimensionless
    asymmetries.  ``ratio_to_room`` is
    ``abs(measured_delta_acp) / documented_np_room_abs``.  This is not a
    prediction; it is only the explicit bookkeeping used by the C003 INFO
    stub.
    """

    measured_delta_acp: float
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


def compare_delta_acp_np_room_to_measurement(
    *,
    measured_delta_acp: float,
    experimental_uncertainty: float,
    documented_np_room_abs: float,
) -> CharmDirectCPRoomComparison:
    """Compare a documented NP room with the measured Delta A_CP anchor.

    This function performs no penguin calculation.  It is intentionally
    small so constraints can report a real numeric budget/ratio while keeping
    the missing SM and RS physics visible to reviewers.
    """

    measured = _finite_float(measured_delta_acp, name="measured_delta_acp")
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
    return CharmDirectCPRoomComparison(
        measured_delta_acp=float(measured),
        experimental_uncertainty=float(uncertainty),
        measurement_abs=float(measurement_abs),
        documented_np_room_abs=float(room),
        ratio_to_room=float(ratio),
        passes=bool(ratio <= 1.0),
    )
