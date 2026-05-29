"""CKM extractions from source-level semileptonic inputs.

This module owns small, model-independent arithmetic used by catalog
constraints that need to turn measured products into CKM magnitudes.  It does
not implement RS charged-current matching.  For K_l3,

    |V_us| = (|V_us| f_+(0)) / f_+(0),

with uncorrelated product, lattice, and optional radiative/isospin uncertainty
components propagated in quadrature.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Mapping

__all__ = [
    "KL3VusExtraction",
    "VusConsistencyResult",
    "extract_vus_from_kl3",
    "vus_consistency_pull",
]


@dataclass(frozen=True)
class KL3VusExtraction:
    """Extracted ``|V_us|`` and propagated uncertainty components."""

    fplus_vus_product: float
    fplus_vus_uncertainty: float
    fplus_zero: float
    fplus_zero_uncertainty: float
    vus: float
    experimental_uncertainty: float
    lattice_uncertainty: float
    radiative_isospin_uncertainty: float
    total_uncertainty: float

    @property
    def uncertainty_components(self) -> Mapping[str, float]:
        return {
            "experimental": float(self.experimental_uncertainty),
            "lattice": float(self.lattice_uncertainty),
            "radiative_isospin": float(self.radiative_isospin_uncertainty),
        }


@dataclass(frozen=True)
class VusConsistencyResult:
    """One-sigma consistency pull for two ``|V_us|`` determinations."""

    extracted_vus: float
    reference_vus: float
    delta_vus: float
    budget: float
    pull_sigma: float
    passes: bool


def _finite_float(name: str, value: object) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name}={value!r} is not numeric") from exc
    if not math.isfinite(number):
        raise ValueError(f"{name} must be finite")
    return number


def _positive_float(name: str, value: object) -> float:
    number = _finite_float(name, value)
    if number <= 0.0:
        raise ValueError(f"{name} must be positive")
    return number


def _nonnegative_float(name: str, value: object) -> float:
    number = _finite_float(name, value)
    if number < 0.0:
        raise ValueError(f"{name} must be nonnegative")
    return number


def extract_vus_from_kl3(
    *,
    fplus_vus_product: float,
    fplus_vus_uncertainty: float,
    fplus_zero: float,
    fplus_zero_uncertainty: float,
    radiative_isospin_uncertainty: float = 0.0,
) -> KL3VusExtraction:
    """Return ``|V_us|`` from the K_l3 product and lattice ``f_+(0)``."""

    product = _positive_float("fplus_vus_product", fplus_vus_product)
    product_sigma = _positive_float("fplus_vus_uncertainty", fplus_vus_uncertainty)
    fplus = _positive_float("fplus_zero", fplus_zero)
    fplus_sigma = _positive_float("fplus_zero_uncertainty", fplus_zero_uncertainty)
    radiative_sigma = _nonnegative_float(
        "radiative_isospin_uncertainty",
        radiative_isospin_uncertainty,
    )

    vus = product / fplus
    experimental_sigma = product_sigma / fplus
    lattice_sigma = product * fplus_sigma / (fplus * fplus)
    total_sigma = math.sqrt(
        experimental_sigma * experimental_sigma
        + lattice_sigma * lattice_sigma
        + radiative_sigma * radiative_sigma
    )
    return KL3VusExtraction(
        fplus_vus_product=float(product),
        fplus_vus_uncertainty=float(product_sigma),
        fplus_zero=float(fplus),
        fplus_zero_uncertainty=float(fplus_sigma),
        vus=float(vus),
        experimental_uncertainty=float(experimental_sigma),
        lattice_uncertainty=float(lattice_sigma),
        radiative_isospin_uncertainty=float(radiative_sigma),
        total_uncertainty=float(total_sigma),
    )


def vus_consistency_pull(
    *,
    extracted_vus: float,
    reference_vus: float,
    budget: float,
) -> VusConsistencyResult:
    """Compare an extracted ``|V_us|`` with a reference value."""

    extracted = _finite_float("extracted_vus", extracted_vus)
    reference = _finite_float("reference_vus", reference_vus)
    sigma = _positive_float("budget", budget)
    delta = extracted - reference
    pull = abs(delta) / sigma
    return VusConsistencyResult(
        extracted_vus=float(extracted),
        reference_vus=float(reference),
        delta_vus=float(delta),
        budget=float(sigma),
        pull_sigma=float(pull),
        passes=bool(pull <= 1.0),
    )
