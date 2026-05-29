"""Small adapter for longitudinal vector-boson-scattering CR011 records.

This module intentionally does not implement a first-principles RS recast.
CR011's active observable is a measured fiducial upper limit on
``pp -> jj W_L^+/- W_L^+/-``.  A real theory prediction requires a
polarization-aware SM/EFT likelihood and a model-specific RS strong-EWSB
matching calculation.  The adapter therefore provides only:

* typed validation of the fiducial cross-section limit;
* optional comparison of a human-supplied fiducial prediction to that limit;
* explicit NEEDS-HUMAN-PHYSICS status strings for both missing sides.

No mass proxy, EFT coefficient proxy, form-factor convention, or RS
``M_KK`` reinterpretation is inferred here.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Any, Mapping

VBS_LONGITUDINAL_SM_EFT_GAP_V1 = (
    "NEEDS-HUMAN-PHYSICS: CR011 SM/EFT side needs a dedicated "
    "polarization-aware VBS likelihood or EFT recast, including operator "
    "basis, coefficient profiling, unitarity treatment, and form-factor "
    "convention. The current adapter only records the measured fiducial "
    "longitudinal-WW upper limit."
)

VBS_LONGITUDINAL_RS_MATCHING_GAP_V1 = (
    "NEEDS-HUMAN-PHYSICS: CR011 RS strong-EWSB matching needs the custodial "
    "gauge/scalar resonance spectrum, widths, VV branching fractions, "
    "interference with SM VBS, acceptance, and the mapping from those "
    "amplitudes to the fiducial longitudinal-WW measurement. No M_KK proxy "
    "or composite-Higgs recast is computed."
)

HUMAN_SUPPLIED_SIGMA_RAW_KEY = "cr011_human_fiducial_sigma_fb"
HUMAN_SUPPLIED_SIGMA_SOURCE_RAW_KEY = "cr011_human_fiducial_sigma_source"


@dataclass(frozen=True)
class VBSFiducialLimit:
    """Measured longitudinal-VBS fiducial cross-section upper limit."""

    process_id: str
    value_fb: float
    cl: str | None = None
    experiment: str | None = None
    source: str | None = None
    source_url: str | None = None
    diagnostics: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        _positive_finite(self.value_fb, "value_fb")


@dataclass(frozen=True)
class VBSHumanPrediction:
    """Externally computed fiducial prediction supplied by a human/recast."""

    fiducial_sigma_fb: float
    source: str | None = None

    def __post_init__(self) -> None:
        _nonnegative_finite(self.fiducial_sigma_fb, "fiducial_sigma_fb")


@dataclass(frozen=True)
class VBSLimitComparison:
    """Advisory comparison of a fiducial prediction with a measured limit."""

    passes: bool
    predicted_fiducial_sigma_fb: float | None
    experimental_limit_fb: float
    budget_fb: float
    ratio_to_budget: float | None
    diagnostics: Mapping[str, Any]


def _positive_finite(value: float, name: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0.0:
        raise ValueError(f"{name} must be positive and finite")
    return number


def _nonnegative_finite(value: float, name: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number < 0.0:
        raise ValueError(f"{name} must be non-negative and finite")
    return number


def human_prediction_from_raw(raw: Any) -> VBSHumanPrediction | None:
    """Return an explicitly human-supplied CR011 prediction, if present.

    The accepted raw field is deliberately named
    ``cr011_human_fiducial_sigma_fb`` to avoid mistaking it for a scan-native
    RS prediction.  This path is for externally computed recasts only; missing
    data means CR011 records the measured limit without evaluating a point.
    """

    if raw is None:
        return None

    value: Any
    source: Any = None
    if isinstance(raw, Mapping):
        if HUMAN_SUPPLIED_SIGMA_RAW_KEY not in raw:
            return None
        value = raw[HUMAN_SUPPLIED_SIGMA_RAW_KEY]
        source = raw.get(HUMAN_SUPPLIED_SIGMA_SOURCE_RAW_KEY)
    else:
        if not hasattr(raw, HUMAN_SUPPLIED_SIGMA_RAW_KEY):
            return None
        value = getattr(raw, HUMAN_SUPPLIED_SIGMA_RAW_KEY)
        source = getattr(raw, HUMAN_SUPPLIED_SIGMA_SOURCE_RAW_KEY, None)

    return VBSHumanPrediction(
        fiducial_sigma_fb=float(value),
        source=None if source is None else str(source),
    )


def compare_vbs_fiducial_limit(
    limit: VBSFiducialLimit,
    prediction: VBSHumanPrediction | None,
) -> VBSLimitComparison:
    """Compare an optional human recast prediction to the measured limit.

    With no prediction, the comparison is intentionally non-excluding:
    ``passes=True`` with ``predicted=None`` and ``ratio=None``.  With a
    supplied fiducial cross section, the ratio is ``sigma_fid / limit``.
    """

    diagnostics = {
        "limit_kind": "fiducial_cross_section_upper_limit",
        "limit_units": "fb",
        "prediction_status": (
            "missing-human-recast" if prediction is None else "human-supplied"
        ),
        "needs_human_physics_sm_eft": VBS_LONGITUDINAL_SM_EFT_GAP_V1,
        "needs_human_physics_rs_matching": VBS_LONGITUDINAL_RS_MATCHING_GAP_V1,
    }
    diagnostics.update(dict(limit.diagnostics))

    if prediction is None:
        return VBSLimitComparison(
            passes=True,
            predicted_fiducial_sigma_fb=None,
            experimental_limit_fb=float(limit.value_fb),
            budget_fb=float(limit.value_fb),
            ratio_to_budget=None,
            diagnostics=diagnostics,
        )

    predicted = float(prediction.fiducial_sigma_fb)
    ratio = float(predicted / limit.value_fb)
    diagnostics.update(
        {
            "human_supplied_prediction_used": True,
            "human_prediction_source": prediction.source,
        }
    )
    return VBSLimitComparison(
        passes=bool(ratio <= 1.0),
        predicted_fiducial_sigma_fb=predicted,
        experimental_limit_fb=float(limit.value_fb),
        budget_fb=float(limit.value_fb),
        ratio_to_budget=ratio,
        diagnostics=diagnostics,
    )


__all__ = [
    "HUMAN_SUPPLIED_SIGMA_RAW_KEY",
    "HUMAN_SUPPLIED_SIGMA_SOURCE_RAW_KEY",
    "VBSFiducialLimit",
    "VBSHumanPrediction",
    "VBSLimitComparison",
    "VBS_LONGITUDINAL_SM_EFT_GAP_V1",
    "VBS_LONGITUDINAL_RS_MATCHING_GAP_V1",
    "compare_vbs_fiducial_limit",
    "human_prediction_from_raw",
]
