"""Reusable collider-resonance recast helpers.

The rigorous observable for a direct resonance search is the predicted
``sigma(pp -> X) * BR(X -> final state)`` folded through acceptance, width,
line-shape, interference, and the experiment's limit curve.  The current
``ParameterPoint`` data used by the catalog constraints only carries the KK
mass scale, not a collider production model.  This module therefore provides a
documented mass-limit proxy for benchmark exclusions and a reusable
``sigma*BR`` comparison path for future recasts.

NEEDS-HUMAN-PHYSICS: the mass-limit proxy assumes that the model point follows
the experimental benchmark signal model whose excluded mass interval is quoted
in the catalog sidecar.  A full RS KK-gluon recast must supply the production
cross section, ttbar branching fraction, total width, interference treatment,
and analysis acceptance.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Any, Mapping

MASS_LOWER_BOUND = "mass_lower_bound"
SIGMA_TIMES_BR_UPPER_LIMIT = "sigma_times_br_upper_limit"

COLLIDER_RESONANCE_MASS_PROXY_ASSUMPTION_V1 = (
    "NEEDS-HUMAN-PHYSICS: direct collider resonance recast v1 uses the "
    "catalogued benchmark mass-exclusion edge as a mass lower bound. It does "
    "not compute sigma(pp->X)*BR(X->final state), width, interference, or "
    "acceptance; those require a collider production model and experiment "
    "limit curve."
)


@dataclass(frozen=True)
class ColliderResonanceLimit:
    """Experimental limit used by a collider-resonance constraint."""

    process_id: str
    resonance: str
    final_state: str
    limit_kind: str
    value: float
    units: str
    cl: str | None = None
    source: str | None = None
    source_url: str | None = None
    limit_type: str | None = None
    benchmark_model: str | None = None
    mass_interval_low: float | None = None
    diagnostics: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.limit_kind not in {MASS_LOWER_BOUND, SIGMA_TIMES_BR_UPPER_LIMIT}:
            raise ValueError(f"unsupported collider resonance limit_kind {self.limit_kind!r}")
        _positive_finite(self.value, "limit value")
        if self.mass_interval_low is not None:
            _positive_finite(self.mass_interval_low, "mass_interval_low")


@dataclass(frozen=True)
class ColliderResonancePrediction:
    """Predicted resonance inputs for a direct-search comparison."""

    resonance: str
    final_state: str
    mass_tev: float
    sigma_times_br: float | None = None
    sigma_times_br_units: str | None = None
    matching_assumption: str = COLLIDER_RESONANCE_MASS_PROXY_ASSUMPTION_V1
    diagnostics: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        _positive_finite(self.mass_tev, "mass_tev")
        if self.sigma_times_br is not None:
            _nonnegative_finite(self.sigma_times_br, "sigma_times_br")


@dataclass(frozen=True)
class ColliderResonanceComparison:
    """Result of comparing one prediction to one experimental limit."""

    passes: bool
    predicted_mass_tev: float
    predicted_sigma_times_br: float | None
    experimental_limit: float
    experimental_units: str
    budget: float
    ratio_to_budget: float
    limit_kind: str
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


def kk_mass_tev_from_m_kk_gev(m_kk_gev: float) -> float:
    """Convert a KK mass scale in GeV to a positive finite TeV mass."""

    return float(_positive_finite(float(m_kk_gev), "m_kk_gev") / 1000.0)


def kk_gluon_prediction_from_m_kk_gev(
    m_kk_gev: float,
    *,
    sigma_times_br: float | None = None,
    sigma_times_br_units: str | None = None,
) -> ColliderResonancePrediction:
    """Build the v1 KK-gluon ``ttbar`` resonance prediction from ``M_KK``."""

    return ColliderResonancePrediction(
        resonance="g_KK^(1)",
        final_state="t tbar",
        mass_tev=kk_mass_tev_from_m_kk_gev(m_kk_gev),
        sigma_times_br=sigma_times_br,
        sigma_times_br_units=sigma_times_br_units,
        diagnostics={
            "m_kk_gev": float(m_kk_gev),
            "sigma_times_br_proxy_available": sigma_times_br is not None,
        },
    )


def mass_from_source_gev(source: Any) -> float:
    """Return ``source.M_KK`` as a positive finite GeV mass."""

    if not hasattr(source, "M_KK"):
        raise AttributeError("collider resonance source must provide M_KK")
    return _positive_finite(getattr(source, "M_KK"), "source.M_KK")


def evaluate_resonance_limit(
    prediction: ColliderResonancePrediction,
    limit: ColliderResonanceLimit,
) -> ColliderResonanceComparison:
    """Compare a resonance prediction to a mass or ``sigma*BR`` limit.

    For a mass lower bound, ``ratio_to_budget = m_limit / m_predicted`` so
    ``ratio <= 1`` passes, matching the catalog ``ConstraintResult`` contract.
    For a ``sigma*BR`` upper limit, ``ratio_to_budget = sigmaBR_pred / limit``.
    """

    if prediction.resonance != limit.resonance:
        raise ValueError(
            f"prediction resonance {prediction.resonance!r} does not match "
            f"limit resonance {limit.resonance!r}"
        )
    if prediction.final_state != limit.final_state:
        raise ValueError(
            f"prediction final_state {prediction.final_state!r} does not match "
            f"limit final_state {limit.final_state!r}"
        )

    if limit.limit_kind == MASS_LOWER_BOUND:
        ratio = float(limit.value / prediction.mass_tev)
        passes = bool(prediction.mass_tev >= limit.value)
        predicted_sigma = prediction.sigma_times_br
    else:
        if prediction.sigma_times_br is None:
            raise ValueError(
                "sigma_times_br prediction is required for a sigma*BR limit"
            )
        ratio = float(prediction.sigma_times_br / limit.value)
        passes = bool(ratio <= 1.0)
        predicted_sigma = float(prediction.sigma_times_br)

    diagnostics = {
        "matching_assumption": prediction.matching_assumption,
        "resonance": prediction.resonance,
        "final_state": prediction.final_state,
        "limit_type": limit.limit_type,
        "benchmark_model": limit.benchmark_model,
        "cl": limit.cl,
        "source": limit.source,
        "source_url": limit.source_url,
        "mass_interval_low_tev": limit.mass_interval_low,
        "sigma_times_br_units": prediction.sigma_times_br_units,
        "limit_diagnostics": dict(limit.diagnostics),
        "prediction_diagnostics": dict(prediction.diagnostics),
    }

    return ColliderResonanceComparison(
        passes=passes,
        predicted_mass_tev=float(prediction.mass_tev),
        predicted_sigma_times_br=predicted_sigma,
        experimental_limit=float(limit.value),
        experimental_units=limit.units,
        budget=float(limit.value),
        ratio_to_budget=ratio,
        limit_kind=limit.limit_kind,
        diagnostics=diagnostics,
    )


__all__ = [
    "MASS_LOWER_BOUND",
    "SIGMA_TIMES_BR_UPPER_LIMIT",
    "COLLIDER_RESONANCE_MASS_PROXY_ASSUMPTION_V1",
    "ColliderResonanceLimit",
    "ColliderResonancePrediction",
    "ColliderResonanceComparison",
    "kk_mass_tev_from_m_kk_gev",
    "kk_gluon_prediction_from_m_kk_gev",
    "mass_from_source_gev",
    "evaluate_resonance_limit",
]
