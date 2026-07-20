"""Statistical adapter for semileptonic CKM inclusive/exclusive tensions.

No semileptonic charged-current RS matching core exists in this repository.
This adapter therefore only owns the model-independent observable-level
statistics: asymmetric uncertainty combination and inclusive-vs-exclusive
pulls.  Constraint modules import this adapter rather than duplicating the
statistical machinery locally.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Mapping, Sequence

SEMILEPTONIC_CKM_RS_MATCHING_GAP = (
    "NEEDS-HUMAN-PHYSICS: rigorous RS matching for inclusive/exclusive "
    "|V_cb| and |V_ub| requires charged-current W/W'/KK electroweak "
    "operators, semileptonic form-factor and inclusive-scheme likelihoods, "
    "and CKM-fit correlations that are not available on ParameterPoint."
)


@dataclass(frozen=True)
class AsymmetricUncertainty:
    """One-sigma uncertainty with independent upward and downward sides."""

    upper: float
    lower: float
    components: Mapping[str, float]

    @property
    def symmetric_average(self) -> float:
        return float(0.5 * (self.upper + self.lower))


@dataclass(frozen=True)
class CKMDetermination:
    """A single inclusive or exclusive CKM magnitude determination."""

    label: str
    value: float
    uncertainty_upper: float
    uncertainty_lower: float
    units: str | None = None
    source: str | None = None


@dataclass(frozen=True)
class InclusiveExclusivePull:
    """Pull for one inclusive/exclusive pair."""

    observable: str
    inclusive: CKMDetermination
    exclusive: CKMDetermination
    difference: float
    combined_sigma: float
    pull_sigma: float
    inclusive_sigma_used: float
    exclusive_sigma_used: float
    budget_sigma: float | None = None

    @property
    def ratio_to_budget(self) -> float | None:
        if self.budget_sigma is None:
            return None
        return float(self.pull_sigma / self.budget_sigma)

    @property
    def passes_budget(self) -> bool | None:
        ratio = self.ratio_to_budget
        return None if ratio is None else bool(ratio <= 1.0)


@dataclass(frozen=True)
class CKMTensionSummary:
    """Aggregate scalar summary for a set of CKM tension pulls."""

    pulls: tuple[InclusiveExclusivePull, ...]
    max_pull: InclusiveExclusivePull
    budget_sigma: float
    ratio_to_budget: float
    passes: bool


def _require_finite_nonnegative(name: str, value: object) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name}={value!r} is not numeric") from exc
    if not math.isfinite(number) or number < 0.0:
        raise ValueError(f"{name} must be finite and nonnegative")
    return number


def _quadrature(values: Sequence[float]) -> float:
    return float(math.sqrt(sum(value * value for value in values)))


def uncertainty_from_components(
    components: Mapping[str, object] | None = None,
    *,
    symmetric: object | None = None,
) -> AsymmetricUncertainty:
    """Build a one-sigma asymmetric uncertainty from YAML components.

    Symmetric entries contribute to both sides.  Keys containing ``plus``
    contribute only to the upward side and keys containing ``minus`` contribute
    only to the downward side.
    """

    if symmetric is not None:
        sigma = _require_finite_nonnegative("symmetric", symmetric)
        if sigma <= 0.0:
            raise ValueError("symmetric uncertainty must be positive")
        return AsymmetricUncertainty(
            upper=float(sigma),
            lower=float(sigma),
            components={"symmetric": float(sigma)},
        )

    if not components:
        raise ValueError("uncertainty components are required")

    upper_terms: list[float] = []
    lower_terms: list[float] = []
    normalized: dict[str, float] = {}
    for key, raw_value in components.items():
        label = str(key)
        value = _require_finite_nonnegative(label, raw_value)
        if value <= 0.0:
            continue
        normalized[label] = value
        lowered = label.lower()
        if "plus" in lowered:
            upper_terms.append(value)
        elif "minus" in lowered:
            lower_terms.append(value)
        else:
            upper_terms.append(value)
            lower_terms.append(value)

    upper = _quadrature(upper_terms)
    lower = _quadrature(lower_terms)
    if upper <= 0.0 or lower <= 0.0:
        raise ValueError("uncertainty components must define both sides")
    return AsymmetricUncertainty(
        upper=float(upper),
        lower=float(lower),
        components=normalized,
    )


def inclusive_exclusive_pull(
    observable: str,
    inclusive: CKMDetermination,
    exclusive: CKMDetermination,
    *,
    budget_sigma: float | None = None,
) -> InclusiveExclusivePull:
    """Return the direction-aware pull between inclusive and exclusive values."""

    if inclusive.units is not None and exclusive.units is not None:
        if inclusive.units != exclusive.units:
            raise ValueError(
                f"{observable}: units differ ({inclusive.units!r} vs {exclusive.units!r})"
            )
    difference = float(inclusive.value - exclusive.value)
    if difference >= 0.0:
        inc_sigma = float(inclusive.uncertainty_lower)
        exc_sigma = float(exclusive.uncertainty_upper)
    else:
        inc_sigma = float(inclusive.uncertainty_upper)
        exc_sigma = float(exclusive.uncertainty_lower)
    if inc_sigma <= 0.0 or exc_sigma <= 0.0:
        raise ValueError(f"{observable}: uncertainties must be positive")
    combined = math.sqrt(inc_sigma * inc_sigma + exc_sigma * exc_sigma)
    if combined <= 0.0:
        raise ValueError(f"{observable}: combined sigma must be positive")
    checked_budget = None
    if budget_sigma is not None:
        checked_budget = float(budget_sigma)
        if checked_budget <= 0.0 or not math.isfinite(checked_budget):
            raise ValueError(f"{observable}: budget_sigma must be finite and positive")
    return InclusiveExclusivePull(
        observable=observable,
        inclusive=inclusive,
        exclusive=exclusive,
        difference=difference,
        combined_sigma=float(combined),
        pull_sigma=float(abs(difference) / combined),
        inclusive_sigma_used=inc_sigma,
        exclusive_sigma_used=exc_sigma,
        budget_sigma=checked_budget,
    )


def summarize_ckm_tensions(
    pulls: Sequence[InclusiveExclusivePull],
    *,
    budget_sigma: float,
) -> CKMTensionSummary:
    """Return the largest pull and its ratio to the supplied sigma budget."""

    if not pulls:
        raise ValueError("at least one pull is required")
    budget = float(budget_sigma)
    if budget <= 0.0 or not math.isfinite(budget):
        raise ValueError("budget_sigma must be finite and positive")
    pull_tuple = tuple(pulls)
    max_pull = max(pull_tuple, key=lambda pull: pull.pull_sigma)
    ratio = float(max_pull.pull_sigma / budget)
    return CKMTensionSummary(
        pulls=pull_tuple,
        max_pull=max_pull,
        budget_sigma=budget,
        ratio_to_budget=ratio,
        passes=bool(ratio <= 1.0),
    )
