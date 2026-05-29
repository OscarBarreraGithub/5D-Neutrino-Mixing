"""Adapter over :mod:`quarkConstraints.rare_b_dilepton`.

This is the catalog boundary for shared ``b -> q l+ l-`` rare-beauty
machinery.  Constraint modules import this adapter only; the underlying
Hamiltonian convention, SM formula, and documented RS C9/C10 proxy remain
isolated in ``quarkConstraints`` for reuse by B005, B006, and later
semileptonic ``b -> s l l`` siblings.
"""

from __future__ import annotations

from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.rare_b_dilepton import (
    RARE_B_DILEPTON_INPUT_BUNDLE_V1,
    RARE_B_DILEPTON_MODEL_V1,
    RARE_B_DILEPTON_OPERATOR_CONVENTION,
    RARE_B_DILEPTON_RS_MATCHING_ASSUMPTION_V1,
    RareBDileptonCKMFactors,
    RareBDileptonMesonInputs,
    RareBDileptonSMInputs,
    RareBDileptonWilsonCoefficients,
    RareBLeptonicBranchingResult,
    ckm_factors as _ckm_factors,
    compute_rare_b_dilepton_wilsons as _compute_rare_b_dilepton_wilsons,
    default_sm_inputs as _default_sm_inputs,
    evaluate_bq_to_mumu as _evaluate_bq_to_mumu,
    sm_branching_fraction as _sm_branching_fraction,
)

__all__ = [
    "QuarkMassBasisCouplings",
    "RARE_B_DILEPTON_MODEL_V1",
    "RARE_B_DILEPTON_OPERATOR_CONVENTION",
    "RARE_B_DILEPTON_INPUT_BUNDLE_V1",
    "RARE_B_DILEPTON_RS_MATCHING_ASSUMPTION_V1",
    "RareBDileptonMesonInputs",
    "RareBDileptonSMInputs",
    "RareBDileptonCKMFactors",
    "RareBDileptonWilsonCoefficients",
    "RareBLeptonicBranchingResult",
    "rare_b_dilepton_default_sm_inputs",
    "rare_b_dilepton_ckm_factors",
    "rare_b_dilepton_sm_branching_fraction",
    "rare_b_dilepton_wilsons_from_couplings",
    "bq_mumu_from_couplings",
    "bs_mumu_from_couplings",
    "bd_mumu_from_couplings",
]


def rare_b_dilepton_default_sm_inputs() -> RareBDileptonSMInputs:
    """Return the default rare-beauty dilepton SM input bundle."""
    return _default_sm_inputs()


def rare_b_dilepton_ckm_factors(
    transition: str = "b_s",
    inputs: RareBDileptonSMInputs | None = None,
) -> RareBDileptonCKMFactors:
    """Return the CKM factors used by the rare-beauty core."""
    return _ckm_factors(transition, inputs)


def rare_b_dilepton_sm_branching_fraction(
    transition: str = "b_s",
    inputs: RareBDileptonSMInputs | None = None,
) -> RareBLeptonicBranchingResult:
    """Evaluate the SM-limit ``B_q -> mu+ mu-`` branching fraction."""
    return _sm_branching_fraction(transition, inputs)


def rare_b_dilepton_wilsons_from_couplings(
    couplings: QuarkMassBasisCouplings,
    *,
    transition: str = "b_s",
    m_kk_gev: float | None = None,
    inputs: RareBDileptonSMInputs | None = None,
) -> RareBDileptonWilsonCoefficients:
    """Return the v1 ``b -> q mu mu`` Wilson proxy for mass-basis couplings."""
    return _compute_rare_b_dilepton_wilsons(
        couplings,
        transition=transition,
        m_kk_gev=m_kk_gev,
        inputs=inputs,
    )


def bq_mumu_from_couplings(
    couplings: QuarkMassBasisCouplings,
    *,
    transition: str,
    m_kk_gev: float | None = None,
    inputs: RareBDileptonSMInputs | None = None,
) -> RareBLeptonicBranchingResult:
    """Evaluate ``BR(B_q -> mu+ mu-)`` from mass-basis couplings."""
    return _evaluate_bq_to_mumu(
        couplings,
        transition=transition,
        m_kk_gev=m_kk_gev,
        inputs=inputs,
    )


def bs_mumu_from_couplings(
    couplings: QuarkMassBasisCouplings,
    *,
    m_kk_gev: float | None = None,
    inputs: RareBDileptonSMInputs | None = None,
) -> RareBLeptonicBranchingResult:
    """Evaluate ``BR(B_s -> mu+ mu-)`` from mass-basis couplings."""
    return bq_mumu_from_couplings(
        couplings,
        transition="b_s",
        m_kk_gev=m_kk_gev,
        inputs=inputs,
    )


def bd_mumu_from_couplings(
    couplings: QuarkMassBasisCouplings,
    *,
    m_kk_gev: float | None = None,
    inputs: RareBDileptonSMInputs | None = None,
) -> RareBLeptonicBranchingResult:
    """Evaluate ``BR(B_d -> mu+ mu-)`` from mass-basis couplings."""
    return bq_mumu_from_couplings(
        couplings,
        transition="b_d",
        m_kk_gev=m_kk_gev,
        inputs=inputs,
    )
