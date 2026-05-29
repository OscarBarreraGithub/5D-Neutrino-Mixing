"""Adapter over :mod:`quarkConstraints.rare_charm_dilepton`.

This is the catalog boundary for shared ``c -> u l+ l-`` rare-charm
machinery.  Constraint modules import this adapter only; the underlying
Hamiltonian convention, D0 leptonic short-distance formula, and documented RS
C9/C10 proxy remain isolated in ``quarkConstraints`` for reuse by C004, C005,
and later semileptonic charm siblings.
"""

from __future__ import annotations

from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.rare_charm_dilepton import (
    RARE_CHARM_DILEPTON_INPUT_BUNDLE_V1,
    RARE_CHARM_DILEPTON_MODEL_V1,
    RARE_CHARM_DILEPTON_OPERATOR_CONVENTION,
    RARE_CHARM_DILEPTON_RS_MATCHING_ASSUMPTION_V1,
    RareCharmDileptonCKMFactors,
    RareCharmDileptonSMInputs,
    RareCharmDileptonWilsonCoefficients,
    RareCharmLeptonInputs,
    RareCharmLeptonicBranchingResult,
    RareCharmMesonInputs,
    ckm_factors as _ckm_factors,
    compute_rare_charm_dilepton_wilsons as _compute_rare_charm_dilepton_wilsons,
    default_sm_inputs as _default_sm_inputs,
    evaluate_d0_to_ll as _evaluate_d0_to_ll,
    sm_branching_fraction as _sm_branching_fraction,
)

__all__ = [
    "QuarkMassBasisCouplings",
    "RARE_CHARM_DILEPTON_MODEL_V1",
    "RARE_CHARM_DILEPTON_OPERATOR_CONVENTION",
    "RARE_CHARM_DILEPTON_INPUT_BUNDLE_V1",
    "RARE_CHARM_DILEPTON_RS_MATCHING_ASSUMPTION_V1",
    "RareCharmLeptonInputs",
    "RareCharmMesonInputs",
    "RareCharmDileptonSMInputs",
    "RareCharmDileptonCKMFactors",
    "RareCharmDileptonWilsonCoefficients",
    "RareCharmLeptonicBranchingResult",
    "rare_charm_dilepton_default_sm_inputs",
    "rare_charm_dilepton_ckm_factors",
    "rare_charm_dilepton_sm_branching_fraction",
    "rare_charm_dilepton_wilsons_from_couplings",
    "d0_ll_from_couplings",
    "d0_mumu_from_couplings",
    "d0_ee_from_couplings",
]


def rare_charm_dilepton_default_sm_inputs() -> RareCharmDileptonSMInputs:
    """Return the default rare-charm dilepton SM input bundle."""
    return _default_sm_inputs()


def rare_charm_dilepton_ckm_factors(
    transition: str = "c_u",
    inputs: RareCharmDileptonSMInputs | None = None,
) -> RareCharmDileptonCKMFactors:
    """Return the CKM factors used by the rare-charm core."""
    return _ckm_factors(transition, inputs)


def rare_charm_dilepton_sm_branching_fraction(
    lepton: str = "mu",
    inputs: RareCharmDileptonSMInputs | None = None,
) -> RareCharmLeptonicBranchingResult:
    """Evaluate the short-distance SM-limit ``D0 -> l+ l-`` branching fraction."""
    return _sm_branching_fraction(lepton, inputs)


def rare_charm_dilepton_wilsons_from_couplings(
    couplings: QuarkMassBasisCouplings,
    *,
    transition: str = "c_u",
    m_kk_gev: float | None = None,
    inputs: RareCharmDileptonSMInputs | None = None,
) -> RareCharmDileptonWilsonCoefficients:
    """Return the v1 ``c -> u l l`` Wilson proxy for mass-basis couplings."""
    return _compute_rare_charm_dilepton_wilsons(
        couplings,
        transition=transition,
        m_kk_gev=m_kk_gev,
        inputs=inputs,
    )


def d0_ll_from_couplings(
    couplings: QuarkMassBasisCouplings,
    *,
    lepton: str,
    m_kk_gev: float | None = None,
    inputs: RareCharmDileptonSMInputs | None = None,
) -> RareCharmLeptonicBranchingResult:
    """Evaluate short-distance ``BR(D0 -> l+ l-)`` from mass-basis couplings."""
    return _evaluate_d0_to_ll(
        couplings,
        lepton=lepton,
        m_kk_gev=m_kk_gev,
        inputs=inputs,
    )


def d0_mumu_from_couplings(
    couplings: QuarkMassBasisCouplings,
    *,
    m_kk_gev: float | None = None,
    inputs: RareCharmDileptonSMInputs | None = None,
) -> RareCharmLeptonicBranchingResult:
    """Evaluate short-distance ``BR(D0 -> mu+ mu-)`` from mass-basis couplings."""
    return d0_ll_from_couplings(
        couplings,
        lepton="mu",
        m_kk_gev=m_kk_gev,
        inputs=inputs,
    )


def d0_ee_from_couplings(
    couplings: QuarkMassBasisCouplings,
    *,
    m_kk_gev: float | None = None,
    inputs: RareCharmDileptonSMInputs | None = None,
) -> RareCharmLeptonicBranchingResult:
    """Evaluate short-distance ``BR(D0 -> e+ e-)`` from mass-basis couplings."""
    return d0_ll_from_couplings(
        couplings,
        lepton="e",
        m_kk_gev=m_kk_gev,
        inputs=inputs,
    )
