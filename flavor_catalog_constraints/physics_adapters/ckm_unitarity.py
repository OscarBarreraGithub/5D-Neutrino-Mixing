"""Small arithmetic adapter for first-row CKM unitarity tests.

This module intentionally contains no RS matching.  It provides only the
data-level arithmetic for

    Delta_CKM = |V_ud|^2 + |V_us|^2 + |V_ub|^2 - 1,

or, when a catalog source already quotes the combined first-row sum, the
equivalent sum-vs-unity comparison.  Model-dependent charged-current RS
matching requires electroweak KK/Z/Z'/W couplings and possible G_F shifts that
are not available on the current ParameterPoint.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

__all__ = [
    "CKMFirstRowUnitarityResult",
    "first_row_sum_from_elements",
    "first_row_uncertainty_from_elements",
    "evaluate_first_row_sum",
    "evaluate_first_row_elements",
]


@dataclass(frozen=True)
class CKMFirstRowUnitarityResult:
    """Typed result for the first-row sum-vs-unity comparison."""

    first_row_sum: float
    target_sum: float
    delta_ckm: float
    uncertainty: float
    pull_sigma: float
    passes_one_sigma: bool


def _finite_float(name: str, value: float) -> float:
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"{name} must be finite")
    return number


def _positive_float(name: str, value: float) -> float:
    number = _finite_float(name, value)
    if number <= 0.0:
        raise ValueError(f"{name} must be positive")
    return number


def first_row_sum_from_elements(vud: float, vus: float, vub: float) -> float:
    """Return ``|V_ud|^2 + |V_us|^2 + |V_ub|^2``."""
    vud_f = _finite_float("vud", vud)
    vus_f = _finite_float("vus", vus)
    vub_f = _finite_float("vub", vub)
    return float(vud_f * vud_f + vus_f * vus_f + vub_f * vub_f)


def first_row_uncertainty_from_elements(
    *,
    vud: float,
    sigma_vud: float,
    vus: float,
    sigma_vus: float,
    vub: float,
    sigma_vub: float,
) -> float:
    """Return the uncorrelated one-sigma uncertainty on the first-row sum."""
    vud_f = _finite_float("vud", vud)
    vus_f = _finite_float("vus", vus)
    vub_f = _finite_float("vub", vub)
    sigma_vud_f = _positive_float("sigma_vud", sigma_vud)
    sigma_vus_f = _positive_float("sigma_vus", sigma_vus)
    sigma_vub_f = _positive_float("sigma_vub", sigma_vub)
    return float(
        math.sqrt(
            (2.0 * vud_f * sigma_vud_f) ** 2
            + (2.0 * vus_f * sigma_vus_f) ** 2
            + (2.0 * vub_f * sigma_vub_f) ** 2
        )
    )


def evaluate_first_row_sum(
    first_row_sum: float,
    uncertainty: float,
    *,
    target_sum: float = 1.0,
) -> CKMFirstRowUnitarityResult:
    """Compare a quoted first-row sum with the SM unitarity target."""
    observed = _finite_float("first_row_sum", first_row_sum)
    target = _finite_float("target_sum", target_sum)
    sigma = _positive_float("uncertainty", uncertainty)
    delta = float(observed - target)
    pull = float(abs(delta) / sigma)
    return CKMFirstRowUnitarityResult(
        first_row_sum=observed,
        target_sum=target,
        delta_ckm=delta,
        uncertainty=sigma,
        pull_sigma=pull,
        passes_one_sigma=bool(pull <= 1.0),
    )


def evaluate_first_row_elements(
    *,
    vud: float,
    sigma_vud: float,
    vus: float,
    sigma_vus: float,
    vub: float,
    sigma_vub: float,
    target_sum: float = 1.0,
) -> CKMFirstRowUnitarityResult:
    """Build and compare the first-row sum from measured CKM magnitudes."""
    first_row_sum = first_row_sum_from_elements(vud, vus, vub)
    uncertainty = first_row_uncertainty_from_elements(
        vud=vud,
        sigma_vud=sigma_vud,
        vus=vus,
        sigma_vus=sigma_vus,
        vub=vub,
        sigma_vub=sigma_vub,
    )
    return evaluate_first_row_sum(
        first_row_sum,
        uncertainty,
        target_sum=target_sum,
    )
