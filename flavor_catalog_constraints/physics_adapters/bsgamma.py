"""Adapter over :mod:`quarkConstraints.bsgamma`.

This is the catalog boundary for inclusive ``Bbar -> X_s gamma`` and reusable
``b -> s gamma`` C7 dipole machinery.  Constraint modules import this adapter
only; the underlying physics implementation remains isolated in
``quarkConstraints``.
"""

from __future__ import annotations

from quarkConstraints.bsgamma import (
    BSGAMMA_INPUT_BUNDLE_V1,
    BSGAMMA_MODEL_V1,
    BSGAMMA_OPERATOR_CONVENTION,
    BSGAMMA_RS_MATCHING_ASSUMPTION_V1,
    BsgammaBranchingResult,
    BsgammaSMInputs,
    BsgammaWilsonCoefficients,
    branching_fraction_from_c7 as _branching_fraction_from_c7,
    compute_bsgamma_wilsons as _compute_bsgamma_wilsons,
    default_sm_inputs as _default_sm_inputs,
    evaluate_inclusive_bsgamma as _evaluate_inclusive_bsgamma,
)
from quarkConstraints.couplings import QuarkMassBasisCouplings

__all__ = [
    "QuarkMassBasisCouplings",
    "BSGAMMA_MODEL_V1",
    "BSGAMMA_OPERATOR_CONVENTION",
    "BSGAMMA_INPUT_BUNDLE_V1",
    "BSGAMMA_RS_MATCHING_ASSUMPTION_V1",
    "BsgammaSMInputs",
    "BsgammaWilsonCoefficients",
    "BsgammaBranchingResult",
    "bsgamma_default_sm_inputs",
    "bsgamma_branching_fraction_from_c7",
    "bsgamma_wilsons_from_couplings",
    "inclusive_bsgamma_from_couplings",
    "inclusive_bsgamma_sm_branching_fraction",
]


def bsgamma_default_sm_inputs() -> BsgammaSMInputs:
    """Return the default C7-normalized ``b -> s gamma`` input bundle."""

    return _default_sm_inputs()


def bsgamma_branching_fraction_from_c7(
    *,
    c7_np: complex = 0.0j,
    c7p_np: complex = 0.0j,
    sm_branching_fraction: float,
    inputs: BsgammaSMInputs | None = None,
) -> BsgammaBranchingResult:
    """Evaluate inclusive ``Bbar -> X_s gamma`` from explicit C7 shifts."""

    return _branching_fraction_from_c7(
        c7_np=c7_np,
        c7p_np=c7p_np,
        sm_branching_fraction=sm_branching_fraction,
        inputs=inputs,
    )


def bsgamma_wilsons_from_couplings(
    couplings: QuarkMassBasisCouplings,
    *,
    m_kk_gev: float | None = None,
    inputs: BsgammaSMInputs | None = None,
) -> BsgammaWilsonCoefficients:
    """Return the v1 ``b -> s gamma`` C7 proxy for mass-basis couplings."""

    return _compute_bsgamma_wilsons(
        couplings,
        m_kk_gev=m_kk_gev,
        inputs=inputs,
    )


def inclusive_bsgamma_from_couplings(
    couplings: QuarkMassBasisCouplings,
    *,
    sm_branching_fraction: float,
    m_kk_gev: float | None = None,
    inputs: BsgammaSMInputs | None = None,
) -> BsgammaBranchingResult:
    """Evaluate ``BR(Bbar -> X_s gamma)`` from mass-basis couplings."""

    return _evaluate_inclusive_bsgamma(
        couplings,
        sm_branching_fraction=sm_branching_fraction,
        m_kk_gev=m_kk_gev,
        inputs=inputs,
    )


def inclusive_bsgamma_sm_branching_fraction(
    *,
    sm_branching_fraction: float,
    inputs: BsgammaSMInputs | None = None,
) -> BsgammaBranchingResult:
    """Evaluate the SM-limit inclusive ``Bbar -> X_s gamma`` branching fraction."""

    return _evaluate_inclusive_bsgamma(
        None,
        sm_branching_fraction=sm_branching_fraction,
        inputs=inputs,
    )

