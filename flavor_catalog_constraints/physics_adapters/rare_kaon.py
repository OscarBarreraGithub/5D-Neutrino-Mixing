"""Adapter over :mod:`quarkConstraints.rare_kaon_snd`.

This is the catalog boundary for the new Delta-S=1
``K+ -> pi+ nu nubar`` machinery.  Constraint modules import this adapter
only; the underlying physics implementation remains isolated in
``quarkConstraints``.
"""

from __future__ import annotations

from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.rare_kaon_snd import (
    RARE_KAON_INPUT_BUNDLE_V1,
    RARE_KAON_OPERATOR_CONVENTION,
    RARE_KAON_RS_MATCHING_ASSUMPTION_V1,
    RARE_KAON_SND_MODEL_V1,
    RareKaonBranchingResult,
    RareKaonCKMFactors,
    RareKaonSMInputs,
    RareKaonWilsonCoefficients,
    ckm_factors as _ckm_factors,
    compute_rare_kaon_wilsons as _compute_rare_kaon_wilsons,
    default_sm_inputs as _default_sm_inputs,
    evaluate_kplus_to_piplus_nunu as _evaluate_kplus_to_piplus_nunu,
    g_sm_squared as _g_sm_squared,
    kappa_plus as _kappa_plus,
    sm_branching_fraction as _sm_branching_fraction,
)

__all__ = [
    "QuarkMassBasisCouplings",
    "RARE_KAON_SND_MODEL_V1",
    "RARE_KAON_OPERATOR_CONVENTION",
    "RARE_KAON_INPUT_BUNDLE_V1",
    "RARE_KAON_RS_MATCHING_ASSUMPTION_V1",
    "RareKaonSMInputs",
    "RareKaonCKMFactors",
    "RareKaonWilsonCoefficients",
    "RareKaonBranchingResult",
    "rare_kaon_default_sm_inputs",
    "rare_kaon_ckm_factors",
    "rare_kaon_kappa_plus",
    "rare_kaon_g_sm_squared",
    "rare_kaon_sm_branching_fraction",
    "rare_kaon_wilsons_from_couplings",
    "kplus_piplus_nunu_from_couplings",
]


def rare_kaon_default_sm_inputs() -> RareKaonSMInputs:
    """Return the default rare-kaon SM input bundle."""
    return _default_sm_inputs()


def rare_kaon_ckm_factors(
    inputs: RareKaonSMInputs | None = None,
) -> RareKaonCKMFactors:
    """Return the CKM factors used by the rare-kaon core."""
    return _ckm_factors(inputs)


def rare_kaon_kappa_plus(inputs: RareKaonSMInputs | None = None) -> float:
    """Return the charged-mode hadronic factor ``kappa_+``."""
    return _kappa_plus(inputs)


def rare_kaon_g_sm_squared(inputs: RareKaonSMInputs | None = None) -> float:
    """Return the Buras ``g_SM^2`` normalization."""
    return _g_sm_squared(inputs)


def rare_kaon_sm_branching_fraction(
    inputs: RareKaonSMInputs | None = None,
) -> RareKaonBranchingResult:
    """Evaluate the SM-limit charged rare-kaon branching fraction."""
    return _sm_branching_fraction(inputs)


def rare_kaon_wilsons_from_couplings(
    couplings: QuarkMassBasisCouplings,
    *,
    m_kk_gev: float | None = None,
    inputs: RareKaonSMInputs | None = None,
) -> RareKaonWilsonCoefficients:
    """Return the v1 ``s -> d nu nubar`` Wilson proxy for mass-basis couplings."""
    return _compute_rare_kaon_wilsons(
        couplings,
        m_kk_gev=m_kk_gev,
        inputs=inputs,
    )


def kplus_piplus_nunu_from_couplings(
    couplings: QuarkMassBasisCouplings,
    *,
    m_kk_gev: float | None = None,
    inputs: RareKaonSMInputs | None = None,
) -> RareKaonBranchingResult:
    """Evaluate ``BR(K+ -> pi+ nu nubar)`` from mass-basis couplings."""
    return _evaluate_kplus_to_piplus_nunu(
        couplings,
        m_kk_gev=m_kk_gev,
        inputs=inputs,
    )
