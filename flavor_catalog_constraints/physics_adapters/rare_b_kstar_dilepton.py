"""Adapter over :mod:`quarkConstraints.rare_b_kstar_dilepton`.

Constraint modules import this adapter only.  The underlying K* form-factor
proxy reuses the shared ``quarkConstraints.rare_b_dilepton`` C9/C10 Wilson
matching and adds the vector-mode normalization needed by B019.
"""

from __future__ import annotations

from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.rare_b_kstar_dilepton import (
    RARE_B_DILEPTON_EXCLUSIVE_BKSTAR_FORM_FACTOR_BUNDLE_V1,
    RARE_B_DILEPTON_EXCLUSIVE_BKSTAR_FORM_FACTOR_ASSUMPTION_V1,
    RARE_B_DILEPTON_EXCLUSIVE_BKSTAR_FORM_FACTOR_CITATION,
    RARE_B_DILEPTON_EXCLUSIVE_BKSTAR_INPUT_BUNDLE_V1,
    RARE_B_DILEPTON_EXCLUSIVE_BKSTAR_LIMITATION_V1,
    RARE_B_DILEPTON_EXCLUSIVE_BKSTAR_MODEL_V1,
    RareBToKStarDileptonBranchingResult,
    RareBToKStarDileptonInputs,
    RareBToKStarFormFactorInputs,
    b_to_kstar_form_factors as _b_to_kstar_form_factors,
    default_b_to_kstar_dilepton_inputs as _default_b_to_kstar_dilepton_inputs,
    evaluate_b_to_kstar_mumu as _evaluate_b_to_kstar_mumu,
    sm_b_to_kstar_mumu_branching_fraction as _sm_b_to_kstar_mumu_branching_fraction,
)
from quarkConstraints.rare_b_dilepton import (
    RARE_B_DILEPTON_RS_MATCHING_ASSUMPTION_V1,
    RareBDileptonWilsonCoefficients,
)

__all__ = [
    "QuarkMassBasisCouplings",
    "RareBDileptonWilsonCoefficients",
    "RARE_B_DILEPTON_RS_MATCHING_ASSUMPTION_V1",
    "RARE_B_DILEPTON_EXCLUSIVE_BKSTAR_MODEL_V1",
    "RARE_B_DILEPTON_EXCLUSIVE_BKSTAR_INPUT_BUNDLE_V1",
    "RARE_B_DILEPTON_EXCLUSIVE_BKSTAR_FORM_FACTOR_BUNDLE_V1",
    "RARE_B_DILEPTON_EXCLUSIVE_BKSTAR_FORM_FACTOR_ASSUMPTION_V1",
    "RARE_B_DILEPTON_EXCLUSIVE_BKSTAR_FORM_FACTOR_CITATION",
    "RARE_B_DILEPTON_EXCLUSIVE_BKSTAR_LIMITATION_V1",
    "RareBToKStarFormFactorInputs",
    "RareBToKStarDileptonInputs",
    "RareBToKStarDileptonBranchingResult",
    "rare_b_to_kstar_dilepton_default_inputs",
    "rare_b_to_kstar_form_factors",
    "rare_b_to_kstar_mumu_sm_branching_fraction",
    "rare_b_to_kstar_mumu_branching_fraction",
    "bzero_kstarzero_mumu_from_couplings",
]


def rare_b_to_kstar_dilepton_default_inputs() -> RareBToKStarDileptonInputs:
    """Return the default exclusive ``B -> K* mu+ mu-`` proxy input bundle."""

    return _default_b_to_kstar_dilepton_inputs()


def rare_b_to_kstar_form_factors(
    q2_gev2: float,
    mode: RareBToKStarFormFactorInputs,
) -> dict[str, float]:
    """Evaluate adapter-exposed ``B -> K*`` proxy form factors."""

    return dict(_b_to_kstar_form_factors(q2_gev2, mode))


def rare_b_to_kstar_mumu_sm_branching_fraction(
    *,
    mode: str = "bzero_kstarzero",
    q2_min_gev2: float | None = None,
    q2_max_gev2: float | None = None,
    inputs: RareBToKStarDileptonInputs | None = None,
) -> RareBToKStarDileptonBranchingResult:
    """Evaluate the SM-limit partial ``BR(B -> K* mu+ mu-)`` proxy."""

    return _sm_b_to_kstar_mumu_branching_fraction(
        mode=mode,
        q2_min_gev2=q2_min_gev2,
        q2_max_gev2=q2_max_gev2,
        inputs=inputs,
    )


def rare_b_to_kstar_mumu_branching_fraction(
    source: QuarkMassBasisCouplings | RareBDileptonWilsonCoefficients | None = None,
    *,
    mode: str = "bzero_kstarzero",
    q2_min_gev2: float | None = None,
    q2_max_gev2: float | None = None,
    m_kk_gev: float | None = None,
    inputs: RareBToKStarDileptonInputs | None = None,
) -> RareBToKStarDileptonBranchingResult:
    """Evaluate exclusive ``BR(B -> K* mu+ mu-)`` from Wilsons or couplings."""

    return _evaluate_b_to_kstar_mumu(
        source,
        mode=mode,
        q2_min_gev2=q2_min_gev2,
        q2_max_gev2=q2_max_gev2,
        m_kk_gev=m_kk_gev,
        inputs=inputs,
    )


def bzero_kstarzero_mumu_from_couplings(
    couplings: QuarkMassBasisCouplings,
    *,
    q2_min_gev2: float | None = None,
    q2_max_gev2: float | None = None,
    m_kk_gev: float | None = None,
    inputs: RareBToKStarDileptonInputs | None = None,
) -> RareBToKStarDileptonBranchingResult:
    """Evaluate ``BR(B0 -> K*(892)0 mu+ mu-)`` from mass-basis couplings."""

    return rare_b_to_kstar_mumu_branching_fraction(
        couplings,
        mode="bzero_kstarzero",
        q2_min_gev2=q2_min_gev2,
        q2_max_gev2=q2_max_gev2,
        m_kk_gev=m_kk_gev,
        inputs=inputs,
    )
