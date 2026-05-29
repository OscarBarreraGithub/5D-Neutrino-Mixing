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
    RARE_B_DILEPTON_EXCLUSIVE_BK_FORM_FACTOR_BUNDLE_V1,
    RARE_B_DILEPTON_EXCLUSIVE_BK_FORM_FACTOR_CITATION,
    RARE_B_DILEPTON_EXCLUSIVE_BK_INPUT_BUNDLE_V1,
    RARE_B_DILEPTON_EXCLUSIVE_BK_LIMITATION_V1,
    RARE_B_DILEPTON_EXCLUSIVE_BK_MODEL_V1,
    RARE_B_DILEPTON_INCLUSIVE_XS_INPUT_BUNDLE_V1,
    RARE_B_DILEPTON_INCLUSIVE_XS_LIMITATION_V1,
    RARE_B_DILEPTON_INCLUSIVE_XS_MODEL_V1,
    RARE_B_DILEPTON_INCLUSIVE_XS_PROXY_THEORY_UNCERTAINTY_FRACTION,
    RARE_B_DILEPTON_INCLUSIVE_XS_PROXY_THEORY_UNCERTAINTY_RATIONALE,
    RARE_B_DILEPTON_INPUT_BUNDLE_V1,
    RARE_B_DILEPTON_MODEL_V1,
    RARE_B_DILEPTON_OPERATOR_CONVENTION,
    RARE_B_DILEPTON_RS_MATCHING_ASSUMPTION_V1,
    RareBInclusiveDileptonBranchingResult,
    RareBInclusiveDileptonInputs,
    RareBToKDileptonBranchingResult,
    RareBToKDileptonInputs,
    RareBToKFormFactorInputs,
    RareBDileptonCKMFactors,
    RareBDileptonMesonInputs,
    RareBDileptonSMInputs,
    RareBDileptonWilsonCoefficients,
    RareBLeptonicBranchingResult,
    b_to_k_fplus as _b_to_k_fplus,
    ckm_factors as _ckm_factors,
    compute_rare_b_dilepton_wilsons as _compute_rare_b_dilepton_wilsons,
    default_b_to_k_dilepton_inputs as _default_b_to_k_dilepton_inputs,
    default_inclusive_b_to_xs_dilepton_inputs as _default_inclusive_b_to_xs_dilepton_inputs,
    default_sm_inputs as _default_sm_inputs,
    evaluate_b_to_k_mumu as _evaluate_b_to_k_mumu,
    evaluate_bq_to_mumu as _evaluate_bq_to_mumu,
    evaluate_inclusive_b_to_xs_mumu as _evaluate_inclusive_b_to_xs_mumu,
    sm_b_to_k_mumu_branching_fraction as _sm_b_to_k_mumu_branching_fraction,
    sm_branching_fraction as _sm_branching_fraction,
    sm_inclusive_b_to_xs_mumu_branching_fraction as _sm_inclusive_b_to_xs_mumu_branching_fraction,
)

__all__ = [
    "QuarkMassBasisCouplings",
    "RARE_B_DILEPTON_MODEL_V1",
    "RARE_B_DILEPTON_OPERATOR_CONVENTION",
    "RARE_B_DILEPTON_INPUT_BUNDLE_V1",
    "RARE_B_DILEPTON_RS_MATCHING_ASSUMPTION_V1",
    "RARE_B_DILEPTON_EXCLUSIVE_BK_MODEL_V1",
    "RARE_B_DILEPTON_EXCLUSIVE_BK_INPUT_BUNDLE_V1",
    "RARE_B_DILEPTON_EXCLUSIVE_BK_FORM_FACTOR_BUNDLE_V1",
    "RARE_B_DILEPTON_EXCLUSIVE_BK_FORM_FACTOR_CITATION",
    "RARE_B_DILEPTON_EXCLUSIVE_BK_LIMITATION_V1",
    "RARE_B_DILEPTON_INCLUSIVE_XS_MODEL_V1",
    "RARE_B_DILEPTON_INCLUSIVE_XS_INPUT_BUNDLE_V1",
    "RARE_B_DILEPTON_INCLUSIVE_XS_LIMITATION_V1",
    "RARE_B_DILEPTON_INCLUSIVE_XS_PROXY_THEORY_UNCERTAINTY_FRACTION",
    "RARE_B_DILEPTON_INCLUSIVE_XS_PROXY_THEORY_UNCERTAINTY_RATIONALE",
    "RareBDileptonMesonInputs",
    "RareBDileptonSMInputs",
    "RareBDileptonCKMFactors",
    "RareBDileptonWilsonCoefficients",
    "RareBLeptonicBranchingResult",
    "RareBToKFormFactorInputs",
    "RareBToKDileptonInputs",
    "RareBToKDileptonBranchingResult",
    "RareBInclusiveDileptonInputs",
    "RareBInclusiveDileptonBranchingResult",
    "rare_b_dilepton_default_sm_inputs",
    "rare_b_dilepton_ckm_factors",
    "rare_b_dilepton_sm_branching_fraction",
    "rare_b_dilepton_wilsons_from_couplings",
    "bq_mumu_from_couplings",
    "bs_mumu_from_couplings",
    "bd_mumu_from_couplings",
    "rare_b_to_k_dilepton_default_inputs",
    "rare_b_to_k_fplus",
    "rare_b_to_k_mumu_sm_branching_fraction",
    "rare_b_to_k_mumu_branching_fraction",
    "bplus_kplus_mumu_from_couplings",
    "bzero_kzero_mumu_from_couplings",
    "rare_b_inclusive_xs_dilepton_default_inputs",
    "rare_b_inclusive_xs_mumu_sm_branching_fraction",
    "rare_b_inclusive_xs_mumu_branching_fraction",
    "inclusive_b_to_xs_mumu_from_couplings",
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


def rare_b_to_k_dilepton_default_inputs() -> RareBToKDileptonInputs:
    """Return the default exclusive ``B -> K mu+ mu-`` input bundle."""
    return _default_b_to_k_dilepton_inputs()


def rare_b_to_k_fplus(
    q2_gev2: float,
    mode: RareBToKFormFactorInputs,
) -> float:
    """Evaluate the adapter-exposed ``B -> K`` form factor ``f_+(q^2)``."""
    return _b_to_k_fplus(q2_gev2, mode)


def rare_b_to_k_mumu_sm_branching_fraction(
    *,
    mode: str = "bplus_kplus",
    q2_min_gev2: float | None = None,
    q2_max_gev2: float | None = None,
    inputs: RareBToKDileptonInputs | None = None,
) -> RareBToKDileptonBranchingResult:
    """Evaluate the SM-limit partial ``BR(B -> K mu+ mu-)``."""
    return _sm_b_to_k_mumu_branching_fraction(
        mode=mode,
        q2_min_gev2=q2_min_gev2,
        q2_max_gev2=q2_max_gev2,
        inputs=inputs,
    )


def rare_b_to_k_mumu_branching_fraction(
    source: QuarkMassBasisCouplings | RareBDileptonWilsonCoefficients | None = None,
    *,
    mode: str = "bplus_kplus",
    q2_min_gev2: float | None = None,
    q2_max_gev2: float | None = None,
    m_kk_gev: float | None = None,
    inputs: RareBToKDileptonInputs | None = None,
) -> RareBToKDileptonBranchingResult:
    """Evaluate exclusive ``BR(B -> K mu+ mu-)`` from Wilsons or couplings."""
    return _evaluate_b_to_k_mumu(
        source,
        mode=mode,
        q2_min_gev2=q2_min_gev2,
        q2_max_gev2=q2_max_gev2,
        m_kk_gev=m_kk_gev,
        inputs=inputs,
    )


def bplus_kplus_mumu_from_couplings(
    couplings: QuarkMassBasisCouplings,
    *,
    q2_min_gev2: float | None = None,
    q2_max_gev2: float | None = None,
    m_kk_gev: float | None = None,
    inputs: RareBToKDileptonInputs | None = None,
) -> RareBToKDileptonBranchingResult:
    """Evaluate ``BR(B+ -> K+ mu+ mu-)`` from mass-basis couplings."""
    return rare_b_to_k_mumu_branching_fraction(
        couplings,
        mode="bplus_kplus",
        q2_min_gev2=q2_min_gev2,
        q2_max_gev2=q2_max_gev2,
        m_kk_gev=m_kk_gev,
        inputs=inputs,
    )


def bzero_kzero_mumu_from_couplings(
    couplings: QuarkMassBasisCouplings,
    *,
    q2_min_gev2: float | None = None,
    q2_max_gev2: float | None = None,
    m_kk_gev: float | None = None,
    inputs: RareBToKDileptonInputs | None = None,
) -> RareBToKDileptonBranchingResult:
    """Evaluate ``BR(B0 -> K0 mu+ mu-)`` from mass-basis couplings."""
    return rare_b_to_k_mumu_branching_fraction(
        couplings,
        mode="bzero_kzero",
        q2_min_gev2=q2_min_gev2,
        q2_max_gev2=q2_max_gev2,
        m_kk_gev=m_kk_gev,
        inputs=inputs,
    )


def rare_b_inclusive_xs_dilepton_default_inputs() -> RareBInclusiveDileptonInputs:
    """Return the default inclusive ``B -> X_s mu+ mu-`` input bundle."""
    return _default_inclusive_b_to_xs_dilepton_inputs()


def rare_b_inclusive_xs_mumu_sm_branching_fraction(
    *,
    sm_branching_fraction: float,
    q2_min_gev2: float,
    q2_max_gev2: float,
    inputs: RareBInclusiveDileptonInputs | None = None,
) -> RareBInclusiveDileptonBranchingResult:
    """Evaluate the SM-limit inclusive ``B -> X_s mu+ mu-`` bin."""
    return _sm_inclusive_b_to_xs_mumu_branching_fraction(
        sm_branching_fraction=sm_branching_fraction,
        q2_min_gev2=q2_min_gev2,
        q2_max_gev2=q2_max_gev2,
        inputs=inputs,
    )


def rare_b_inclusive_xs_mumu_branching_fraction(
    source: QuarkMassBasisCouplings | RareBDileptonWilsonCoefficients | None = None,
    *,
    sm_branching_fraction: float,
    q2_min_gev2: float,
    q2_max_gev2: float,
    m_kk_gev: float | None = None,
    inputs: RareBInclusiveDileptonInputs | None = None,
) -> RareBInclusiveDileptonBranchingResult:
    """Evaluate inclusive ``B -> X_s mu+ mu-`` from Wilsons or couplings."""
    return _evaluate_inclusive_b_to_xs_mumu(
        source,
        sm_branching_fraction=sm_branching_fraction,
        q2_min_gev2=q2_min_gev2,
        q2_max_gev2=q2_max_gev2,
        m_kk_gev=m_kk_gev,
        inputs=inputs,
    )


def inclusive_b_to_xs_mumu_from_couplings(
    couplings: QuarkMassBasisCouplings,
    *,
    sm_branching_fraction: float,
    q2_min_gev2: float,
    q2_max_gev2: float,
    m_kk_gev: float | None = None,
    inputs: RareBInclusiveDileptonInputs | None = None,
) -> RareBInclusiveDileptonBranchingResult:
    """Evaluate inclusive ``B -> X_s mu+ mu-`` from mass-basis couplings."""
    return rare_b_inclusive_xs_mumu_branching_fraction(
        couplings,
        sm_branching_fraction=sm_branching_fraction,
        q2_min_gev2=q2_min_gev2,
        q2_max_gev2=q2_max_gev2,
        m_kk_gev=m_kk_gev,
        inputs=inputs,
    )
