"""Stub adapter for ``B_s -> phi phi`` penguin CP phase bookkeeping.

This adapter deliberately does **not** compute the Standard-Model
``B_s -> phi phi`` helicity amplitudes or an RS ``Delta B = 1`` penguin
amplitude.  The SM interpretation of ``phi_s^{s sbar s}`` needs a
QCDF/SCET or phenomenological nonleptonic amplitude treatment with
hadronic penguin pollution and strong phases.  A rigorous RS prediction
would need ``b -> s sbar s`` penguin matching, RG evolution, and
polarization-dependent nonleptonic matrix elements that are not present on
``ParameterPoint``.

The only honest v1 operation is finite scalar bookkeeping: record the
measured ``phi_s^{s sbar s}`` value and compare its magnitude with an
explicitly documented non-vetoing room chosen by the constraint.  It is not a
prediction and must remain INFO-only unless the missing physics is supplied.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

__all__ = [
    "BS_TO_PHI_PHI_NONLEPTONIC_STUB_MODEL_V1",
    "BS_TO_PHI_PHI_SM_NEEDS_HUMAN_PHYSICS",
    "BS_TO_PHI_PHI_RS_NEEDS_HUMAN_PHYSICS",
    "BsPhiPhiPhaseRoomComparison",
    "compare_bs_phiphi_phase_to_room",
]


BS_TO_PHI_PHI_NONLEPTONIC_STUB_MODEL_V1 = (
    "bs_to_phiphi_penguin_cp_no_hadronic_amplitude_stub_v1"
)

BS_TO_PHI_PHI_SM_NEEDS_HUMAN_PHYSICS = (
    "NEEDS-HUMAN-PHYSICS: SM phi_s^{s sbar s} in B_s -> phi phi needs "
    "QCDF/SCET or a phenomenological nonleptonic amplitude treatment with "
    "polarization-dependent hadronic penguin pollution and strong phases; "
    "no first-principles SM central value is available in this catalog."
)

BS_TO_PHI_PHI_RS_NEEDS_HUMAN_PHYSICS = (
    "NEEDS-HUMAN-PHYSICS: rigorous RS phi_s^{s sbar s} needs Delta B=1 "
    "b -> s sbar s penguin matching, RG evolution, and nonleptonic "
    "B_s -> phi phi matrix elements; no RS penguin amplitude is computed."
)


@dataclass(frozen=True)
class BsPhiPhiPhaseRoomComparison:
    """Scalar comparison between measured ``phi_s^{s sbar s}`` and room.

    All phase values are in radians.  ``ratio_to_room`` is
    ``abs(measured_phi_s) / documented_np_room_abs``.  This is only the
    explicit bookkeeping used by the B034 INFO stub; it is not a hadronic
    amplitude or RS penguin calculation.
    """

    measured_phi_s: float
    experimental_uncertainty: float
    measured_abs: float
    documented_np_room_abs: float
    ratio_to_room: float
    passes: bool


def _finite_float(value: float, *, name: str) -> float:
    out = float(value)
    if not math.isfinite(out):
        raise ValueError(f"{name} must be finite, got {value!r}")
    return out


def compare_bs_phiphi_phase_to_room(
    *,
    measured_phi_s: float,
    experimental_uncertainty: float,
    documented_np_room_abs: float,
) -> BsPhiPhiPhaseRoomComparison:
    """Compare a documented phase room with a measured ``B_s -> phi phi`` phase.

    This function performs no QCDF/SCET, ``Delta B = 1`` penguin matching, RG
    evolution, angular-likelihood treatment, or nonleptonic matrix-element
    calculation.  It exists so the constraint can return real scalar
    ``budget``/``ratio`` fields while keeping the missing SM and RS physics
    visible.
    """

    measured = _finite_float(measured_phi_s, name="measured_phi_s")
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

    measured_abs = abs(measured)
    ratio = measured_abs / room
    return BsPhiPhiPhaseRoomComparison(
        measured_phi_s=float(measured),
        experimental_uncertainty=float(uncertainty),
        measured_abs=float(measured_abs),
        documented_np_room_abs=float(room),
        ratio_to_room=float(ratio),
        passes=bool(ratio <= 1.0),
    )
