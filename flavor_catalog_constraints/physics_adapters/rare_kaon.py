"""Adapter over :mod:`quarkConstraints.rare_kaon_snd`.

This is the catalog boundary for the Delta-S=1 ``K -> pi nu nubar``
machinery.  Constraint modules import this adapter only; the underlying physics
implementation remains isolated in ``quarkConstraints``.  K005 intentionally
appends the ``K_L -> pi0 nu nubar`` machinery to the K004
``quarkConstraints/rare_kaon_snd.py`` module and wraps it here; this is the
designed location for this physics, not an isolation violation.
"""

from __future__ import annotations

from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.rare_kaon_snd import (
    RARE_KAON_INPUT_BUNDLE_V1,
    RARE_KAON_KAPPA_L_CITATION,
    RARE_KAON_KAPPA_L_REF,
    RARE_KAON_OPERATOR_CONVENTION,
    RARE_KAON_RS_MATCHING_ASSUMPTION_V1,
    RARE_KAON_SND_MODEL_V1,
    RareKaonBranchingResult,
    RareKaonCKMFactors,
    RareKaonNeutralBranchingResult,
    RareKaonSMInputs,
    RareKaonWilsonCoefficients,
    ckm_factors as _ckm_factors,
    compute_rare_kaon_wilsons as _compute_rare_kaon_wilsons,
    default_sm_inputs as _default_sm_inputs,
    evaluate_klong_to_pi0_nunu as _evaluate_klong_to_pi0_nunu,
    evaluate_kplus_to_piplus_nunu as _evaluate_kplus_to_piplus_nunu,
    g_sm_squared as _g_sm_squared,
    kappa_l as _kappa_l,
    kappa_plus as _kappa_plus,
    neutral_sm_branching_fraction as _neutral_sm_branching_fraction,
    sm_branching_fraction as _sm_branching_fraction,
)

__all__ = [
    "QuarkMassBasisCouplings",
    "RARE_KAON_SND_MODEL_V1",
    "RARE_KAON_OPERATOR_CONVENTION",
    "RARE_KAON_INPUT_BUNDLE_V1",
    "RARE_KAON_RS_MATCHING_ASSUMPTION_V1",
    "RARE_KAON_KAPPA_L_REF",
    "RARE_KAON_KAPPA_L_CITATION",
    "RareKaonSMInputs",
    "RareKaonCKMFactors",
    "RareKaonWilsonCoefficients",
    "RareKaonBranchingResult",
    "RareKaonNeutralBranchingResult",
    "rare_kaon_default_sm_inputs",
    "rare_kaon_ckm_factors",
    "rare_kaon_kappa_plus",
    "rare_kaon_kappa_l",
    "rare_kaon_g_sm_squared",
    "rare_kaon_sm_branching_fraction",
    "rare_kaon_neutral_sm_branching_fraction",
    "rare_kaon_wilsons_from_couplings",
    "kplus_piplus_nunu_from_couplings",
    "klong_pi0_nunu_from_couplings",
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


def rare_kaon_kappa_l(inputs: RareKaonSMInputs | None = None) -> float:
    """Return the neutral-mode hadronic factor ``kappa_L``."""
    return _kappa_l(inputs)


def rare_kaon_g_sm_squared(inputs: RareKaonSMInputs | None = None) -> float:
    """Return the Buras ``g_SM^2`` normalization."""
    return _g_sm_squared(inputs)


def rare_kaon_sm_branching_fraction(
    inputs: RareKaonSMInputs | None = None,
) -> RareKaonBranchingResult:
    """Evaluate the SM-limit charged rare-kaon branching fraction."""
    return _sm_branching_fraction(inputs)


def rare_kaon_neutral_sm_branching_fraction(
    inputs: RareKaonSMInputs | None = None,
) -> RareKaonNeutralBranchingResult:
    """Evaluate the SM-limit neutral rare-kaon branching fraction."""
    return _neutral_sm_branching_fraction(inputs)


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


def klong_pi0_nunu_from_couplings(
    couplings: QuarkMassBasisCouplings,
    *,
    m_kk_gev: float | None = None,
    inputs: RareKaonSMInputs | None = None,
) -> RareKaonNeutralBranchingResult:
    """Evaluate ``BR(K_L -> pi0 nu nubar)`` from mass-basis couplings."""
    return _evaluate_klong_to_pi0_nunu(
        couplings,
        m_kk_gev=m_kk_gev,
        inputs=inputs,
    )
