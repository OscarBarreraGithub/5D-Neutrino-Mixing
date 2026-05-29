"""Adapter over :mod:`quarkConstraints.rare_kaon_dilepton`.

This is the catalog boundary for Delta-S=1 ``s -> d l+l-`` short-distance
machinery.  Constraint modules import this adapter only; the underlying physics
implementation remains isolated in ``quarkConstraints``.  The module starts
with the ``K_L -> mu+ mu-`` short-distance observable and is designed as the
shared place for K008/K009/K010/K012 dilepton rare-kaon reuse.
"""

from __future__ import annotations

from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.rare_kaon_dilepton import (
    KLongMuMuShortDistanceResult,
    RARE_KAON_DILEPTON_INPUT_BUNDLE_V1,
    RARE_KAON_DILEPTON_MODEL_V1,
    RARE_KAON_DILEPTON_OPERATOR_CONVENTION,
    RARE_KAON_DILEPTON_PARAMETRIZATION_CITATION,
    RARE_KAON_DILEPTON_RS_MATCHING_ASSUMPTION_V1,
    RareKaonDileptonCKMFactors,
    RareKaonDileptonSMInputs,
    RareKaonDileptonWilsonCoefficients,
    ckm_factors as _ckm_factors,
    compute_rare_kaon_dilepton_wilsons as _compute_rare_kaon_dilepton_wilsons,
    default_sm_inputs as _default_sm_inputs,
    evaluate_klong_mumu_short_distance as _evaluate_klong_mumu_short_distance,
    g_sm_squared as _g_sm_squared,
    kappa_mu as _kappa_mu,
    klong_mumu_short_distance_sm as _klong_mumu_short_distance_sm,
)

__all__ = [
    "QuarkMassBasisCouplings",
    "RARE_KAON_DILEPTON_MODEL_V1",
    "RARE_KAON_DILEPTON_OPERATOR_CONVENTION",
    "RARE_KAON_DILEPTON_INPUT_BUNDLE_V1",
    "RARE_KAON_DILEPTON_RS_MATCHING_ASSUMPTION_V1",
    "RARE_KAON_DILEPTON_PARAMETRIZATION_CITATION",
    "RareKaonDileptonSMInputs",
    "RareKaonDileptonCKMFactors",
    "RareKaonDileptonWilsonCoefficients",
    "KLongMuMuShortDistanceResult",
    "rare_kaon_dilepton_default_sm_inputs",
    "rare_kaon_dilepton_ckm_factors",
    "rare_kaon_dilepton_kappa_mu",
    "rare_kaon_dilepton_g_sm_squared",
    "rare_kaon_dilepton_wilsons_from_couplings",
    "klong_mumu_short_distance_sm",
    "klong_mumu_short_distance_from_couplings",
]


def rare_kaon_dilepton_default_sm_inputs() -> RareKaonDileptonSMInputs:
    """Return the default rare-kaon dilepton SM input bundle."""
    return _default_sm_inputs()


def rare_kaon_dilepton_ckm_factors(
    inputs: RareKaonDileptonSMInputs | None = None,
) -> RareKaonDileptonCKMFactors:
    """Return the CKM factors used by the rare-kaon dilepton core."""
    return _ckm_factors(inputs)


def rare_kaon_dilepton_kappa_mu(
    inputs: RareKaonDileptonSMInputs | None = None,
) -> float:
    """Return the ``K_L -> mu+ mu-`` short-distance normalization."""
    return _kappa_mu(inputs)


def rare_kaon_dilepton_g_sm_squared(
    inputs: RareKaonDileptonSMInputs | None = None,
) -> float:
    """Return the Buras ``g_SM^2`` normalization."""
    return _g_sm_squared(inputs)


def rare_kaon_dilepton_wilsons_from_couplings(
    couplings: QuarkMassBasisCouplings,
    *,
    m_kk_gev: float | None = None,
    inputs: RareKaonDileptonSMInputs | None = None,
) -> RareKaonDileptonWilsonCoefficients:
    """Return the v1 ``s -> d mu+mu-`` Wilson proxy for mass-basis couplings."""
    return _compute_rare_kaon_dilepton_wilsons(
        couplings,
        m_kk_gev=m_kk_gev,
        inputs=inputs,
    )


def klong_mumu_short_distance_sm(
    inputs: RareKaonDileptonSMInputs | None = None,
) -> KLongMuMuShortDistanceResult:
    """Evaluate the SM-limit ``K_L -> mu+ mu-`` short-distance rate."""
    return _klong_mumu_short_distance_sm(inputs)


def klong_mumu_short_distance_from_couplings(
    couplings: QuarkMassBasisCouplings,
    *,
    m_kk_gev: float | None = None,
    inputs: RareKaonDileptonSMInputs | None = None,
) -> KLongMuMuShortDistanceResult:
    """Evaluate ``BR(K_L -> mu+ mu-)_SD`` from mass-basis couplings."""
    return _evaluate_klong_mumu_short_distance(
        couplings,
        m_kk_gev=m_kk_gev,
        inputs=inputs,
    )
