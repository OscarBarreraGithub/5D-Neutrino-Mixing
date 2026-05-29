"""Adapter over :mod:`quarkConstraints.rare_b_nunu`.

This is the catalog boundary for the Delta-B=1 ``b -> s nu nubar`` machinery.
Constraint modules import this adapter only; the underlying physics
implementation remains isolated in ``quarkConstraints``.  The shared core
exposes both the ``B+ -> K+ nu nubar`` branching-ratio evaluator and reusable
short-distance ``C_L, C_R, epsilon, eta`` parameters for future ``B -> K* nu
nubar`` constraints.
"""

from __future__ import annotations

from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.rare_b_nunu import (
    RARE_B_NUNU_BPLUS_KPLUS_SM_CITATION,
    RARE_B_NUNU_INPUT_BUNDLE_V1,
    RARE_B_NUNU_KSTAR_ETA_COEFFICIENT,
    RARE_B_NUNU_MODEL_V1,
    RARE_B_NUNU_OPERATOR_CONVENTION,
    RARE_B_NUNU_RS_MATCHING_ASSUMPTION_V1,
    RareBNuNuBranchingResult,
    RareBNuNuCKMFactors,
    RareBNuNuSMInputs,
    RareBNuNuShortDistanceResult,
    RareBNuNuWilsonCoefficients,
    ckm_factors as _ckm_factors,
    compute_rare_b_nunu_wilsons as _compute_rare_b_nunu_wilsons,
    default_sm_inputs as _default_sm_inputs,
    evaluate_bplus_to_kplus_nunu as _evaluate_bplus_to_kplus_nunu,
    g_sm_squared as _g_sm_squared,
    inami_lim_x0 as _inami_lim_x0,
    short_distance_response as _short_distance_response,
    sm_branching_fraction as _sm_branching_fraction,
    sm_inputs_with_bplus_kplus_normalization as _sm_inputs_with_bplus_kplus_normalization,
    x_t_top_function as _x_t_top_function,
)

__all__ = [
    "QuarkMassBasisCouplings",
    "RARE_B_NUNU_MODEL_V1",
    "RARE_B_NUNU_OPERATOR_CONVENTION",
    "RARE_B_NUNU_INPUT_BUNDLE_V1",
    "RARE_B_NUNU_RS_MATCHING_ASSUMPTION_V1",
    "RARE_B_NUNU_BPLUS_KPLUS_SM_CITATION",
    "RARE_B_NUNU_KSTAR_ETA_COEFFICIENT",
    "RareBNuNuSMInputs",
    "RareBNuNuCKMFactors",
    "RareBNuNuWilsonCoefficients",
    "RareBNuNuShortDistanceResult",
    "RareBNuNuBranchingResult",
    "rare_b_nunu_default_sm_inputs",
    "rare_b_nunu_sm_inputs_with_bplus_kplus_normalization",
    "rare_b_nunu_inami_lim_x0",
    "rare_b_nunu_x_t_top_function",
    "rare_b_nunu_ckm_factors",
    "rare_b_nunu_g_sm_squared",
    "rare_b_nunu_short_distance_response",
    "rare_b_nunu_sm_branching_fraction",
    "rare_b_nunu_wilsons_from_couplings",
    "bplus_kplus_nunu_from_couplings",
]


def rare_b_nunu_default_sm_inputs() -> RareBNuNuSMInputs:
    """Return the default rare-B ``b -> s nu nubar`` input bundle."""
    return _default_sm_inputs()


def rare_b_nunu_sm_inputs_with_bplus_kplus_normalization(
    branching_fraction: float,
    *,
    inputs: RareBNuNuSMInputs | None = None,
) -> RareBNuNuSMInputs:
    """Return inputs with the ``B+ -> K+`` SM normalization replaced."""
    return _sm_inputs_with_bplus_kplus_normalization(
        branching_fraction,
        inputs=inputs,
    )


def rare_b_nunu_inami_lim_x0(x_t: float) -> float:
    """Return the leading Inami-Lim ``X_0(x_t)`` function."""
    return _inami_lim_x0(x_t)


def rare_b_nunu_x_t_top_function(inputs: RareBNuNuSMInputs | None = None) -> float:
    """Return the short-distance top function ``X_t``."""
    return _x_t_top_function(inputs)


def rare_b_nunu_ckm_factors(
    inputs: RareBNuNuSMInputs | None = None,
) -> RareBNuNuCKMFactors:
    """Return the CKM factors used by the rare-B core."""
    return _ckm_factors(inputs)


def rare_b_nunu_g_sm_squared(inputs: RareBNuNuSMInputs | None = None) -> float:
    """Return the Buras ``g_SM^2`` normalization."""
    return _g_sm_squared(inputs)


def rare_b_nunu_short_distance_response(
    x_np_left: complex = 0.0j,
    x_np_right: complex = 0.0j,
    *,
    inputs: RareBNuNuSMInputs | None = None,
) -> RareBNuNuShortDistanceResult:
    """Return reusable ``C_L, C_R, epsilon, eta`` response parameters."""
    return _short_distance_response(x_np_left, x_np_right, inputs=inputs)


def rare_b_nunu_sm_branching_fraction(
    inputs: RareBNuNuSMInputs | None = None,
) -> RareBNuNuBranchingResult:
    """Evaluate the SM-limit ``B+ -> K+ nu nubar`` branching fraction."""
    return _sm_branching_fraction(inputs)


def rare_b_nunu_wilsons_from_couplings(
    couplings: QuarkMassBasisCouplings,
    *,
    m_kk_gev: float | None = None,
    inputs: RareBNuNuSMInputs | None = None,
) -> RareBNuNuWilsonCoefficients:
    """Return the v1 ``b -> s nu nubar`` Wilson proxy for mass-basis couplings."""
    return _compute_rare_b_nunu_wilsons(
        couplings,
        m_kk_gev=m_kk_gev,
        inputs=inputs,
    )


def bplus_kplus_nunu_from_couplings(
    couplings: QuarkMassBasisCouplings,
    *,
    m_kk_gev: float | None = None,
    inputs: RareBNuNuSMInputs | None = None,
) -> RareBNuNuBranchingResult:
    """Evaluate ``BR(B+ -> K+ nu nubar)`` from mass-basis couplings."""
    return _evaluate_bplus_to_kplus_nunu(
        couplings,
        m_kk_gev=m_kk_gev,
        inputs=inputs,
    )
