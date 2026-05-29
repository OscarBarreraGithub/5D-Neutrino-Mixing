"""Stub adapter for ``B0 -> phi K_S`` penguin CP asymmetry.

This adapter deliberately does **not** compute the Standard-Model
``B -> phi K`` amplitude or an RS ``Delta B = 1`` penguin amplitude.
The SM interpretation of ``S_phi K_S`` needs a hadronic framework such as
QCDF/SCET or a phenomenological nonleptonic amplitude fit, including
subleading penguin pollution and strong phases.  A rigorous RS prediction
would additionally need ``b -> s sbar s`` penguin matching, RG evolution, and
nonleptonic matrix elements that are not present on ``ParameterPoint``.

The only honest v1 operation is finite-number bookkeeping: compare the HFLAV
``S_phi K_S`` measurement with the local HFLAV ``sin(2 beta)`` reference and
report the observed ``Delta S``.  It is not a prediction and must remain
non-vetoing.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

__all__ = [
    "B_TO_PHI_KS_NONLEPTONIC_STUB_MODEL_V1",
    "B_TO_PHI_KS_SM_NEEDS_HUMAN_PHYSICS",
    "B_TO_PHI_KS_RS_NEEDS_HUMAN_PHYSICS",
    "SphiKsReferenceComparison",
    "compare_sphiks_to_sin2beta_reference",
]


B_TO_PHI_KS_NONLEPTONIC_STUB_MODEL_V1 = (
    "b_to_phi_ks_penguin_cp_no_hadronic_amplitude_stub_v1"
)

B_TO_PHI_KS_SM_NEEDS_HUMAN_PHYSICS = (
    "NEEDS-HUMAN-PHYSICS: SM S_phiK_S needs QCDF/SCET or a phenomenological "
    "nonleptonic B -> phi K amplitude treatment with subleading hadronic "
    "penguin pollution and strong phases; sin(2 beta) is only the clean "
    "mixing-phase reference, not a first-principles S_phiK_S prediction."
)

B_TO_PHI_KS_RS_NEEDS_HUMAN_PHYSICS = (
    "NEEDS-HUMAN-PHYSICS: rigorous RS S_phiK_S needs Delta B=1 b -> s sbar s "
    "penguin matching, RG evolution, and nonleptonic matrix elements; no RS "
    "penguin amplitude is computed."
)


@dataclass(frozen=True)
class SphiKsReferenceComparison:
    """Bookkeeping comparison of measured ``S_phiK_S`` with ``sin(2 beta)``.

    ``delta_s`` is ``measured_s_phi_ks - sin2beta_reference``.  The
    ``delta_s_uncertainty`` combines the two quoted HFLAV one-sigma errors in
    quadrature and is used only as an informational scale.  This object carries
    no SM or RS decay-amplitude prediction.
    """

    measured_s_phi_ks: float
    measured_uncertainty: float
    sin2beta_reference: float
    sin2beta_uncertainty: float
    delta_s: float
    delta_s_uncertainty: float
    ratio_to_reference_uncertainty: float
    passes: bool


def _finite_float(value: float, *, name: str) -> float:
    out = float(value)
    if not math.isfinite(out):
        raise ValueError(f"{name} must be finite, got {value!r}")
    return out


def _positive_uncertainty(value: float, *, name: str) -> float:
    out = _finite_float(value, name=name)
    if out <= 0.0:
        raise ValueError(f"{name} must be positive, got {value!r}")
    return out


def compare_sphiks_to_sin2beta_reference(
    *,
    measured_s_phi_ks: float,
    measured_uncertainty: float,
    sin2beta_reference: float,
    sin2beta_uncertainty: float,
) -> SphiKsReferenceComparison:
    """Return the non-vetoing ``Delta S`` bookkeeping comparison.

    This function performs no QCDF/SCET, ``Delta B = 1`` penguin matching, RG
    evolution, or nonleptonic matrix-element calculation.  ``passes`` is only
    the advisory statement ``|Delta S| <= sigma_DeltaS`` for the reported
    reference scale; callers keep the result non-vetoing via INFO severity.
    """

    measured = _finite_float(measured_s_phi_ks, name="measured_s_phi_ks")
    measured_sigma = _positive_uncertainty(
        measured_uncertainty,
        name="measured_uncertainty",
    )
    reference = _finite_float(sin2beta_reference, name="sin2beta_reference")
    reference_sigma = _positive_uncertainty(
        sin2beta_uncertainty,
        name="sin2beta_uncertainty",
    )

    delta = measured - reference
    sigma = math.sqrt(measured_sigma * measured_sigma + reference_sigma * reference_sigma)
    ratio = abs(delta) / sigma
    return SphiKsReferenceComparison(
        measured_s_phi_ks=float(measured),
        measured_uncertainty=float(measured_sigma),
        sin2beta_reference=float(reference),
        sin2beta_uncertainty=float(reference_sigma),
        delta_s=float(delta),
        delta_s_uncertainty=float(sigma),
        ratio_to_reference_uncertainty=float(ratio),
        passes=bool(ratio <= 1.0),
    )
