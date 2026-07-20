"""Adapter over :mod:`quarkConstraints.rare_b_kstar_dilepton`.

Constraint modules import this adapter only.  The underlying K* form-factor
proxy reuses the shared ``quarkConstraints.rare_b_dilepton`` C9/C10 Wilson
matching and adds the vector-mode normalization needed by B019.
"""

from __future__ import annotations

from dataclasses import replace

from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.rare_b_dilepton import (
    RARE_B_DILEPTON_RS_MATCHING_ASSUMPTION_V1,
    RareBDileptonWilsonCoefficients,
)
from quarkConstraints.rare_b_kstar_dilepton import (
    RARE_B_DILEPTON_EXCLUSIVE_BKSTAR_FORM_FACTOR_ASSUMPTION_V1,
    RARE_B_DILEPTON_EXCLUSIVE_BKSTAR_FORM_FACTOR_BUNDLE_V1,
    RARE_B_DILEPTON_EXCLUSIVE_BKSTAR_FORM_FACTOR_CITATION,
    RARE_B_DILEPTON_EXCLUSIVE_BKSTAR_INPUT_BUNDLE_V1,
    RARE_B_DILEPTON_EXCLUSIVE_BKSTAR_LIMITATION_V1,
    RARE_B_DILEPTON_EXCLUSIVE_BKSTAR_MODEL_V1,
    RareBToKStarDileptonBranchingResult,
    RareBToKStarDileptonInputs,
    RareBToKStarFormFactorInputs,
)
from quarkConstraints.rare_b_kstar_dilepton import (
    b_to_kstar_form_factors as _b_to_kstar_form_factors,
)
from quarkConstraints.rare_b_kstar_dilepton import (
    default_b_to_kstar_dilepton_inputs as _default_b_to_kstar_dilepton_inputs,
)
from quarkConstraints.rare_b_kstar_dilepton import (
    evaluate_b_to_kstar_mumu as _evaluate_b_to_kstar_mumu,
)
from quarkConstraints.rare_b_kstar_dilepton import (
    sm_b_to_kstar_mumu_branching_fraction as _sm_b_to_kstar_mumu_branching_fraction,
)
from quarkConstraints.rs_semileptonic_wilsons import RSSemileptonicWilsonBundle

from .rare_b_meson import (
    rare_b_dilepton_wilsons_from_rs_semileptonic,
    rare_b_rs_semileptonic_vector_diagnostics,
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
    "rare_b_to_kstar_mumu_from_rs_semileptonic_wilsons",
    "bzero_kstarzero_mumu_from_rs_semileptonic_wilsons",
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


def rare_b_to_kstar_mumu_from_rs_semileptonic_wilsons(
    source: RSSemileptonicWilsonBundle,
    *,
    lepton: str = "mu",
    mode: str = "bzero_kstarzero",
    q2_min_gev2: float | None = None,
    q2_max_gev2: float | None = None,
    matching_scale_gev: float | None = None,
    inputs: RareBToKStarDileptonInputs | None = None,
) -> RareBToKStarDileptonBranchingResult:
    """Evaluate ``BR(B -> K* l+ l-)`` from Phase-3a RS Wilsons."""

    coeff = source.b_to_s_ll[lepton]
    wilsons = rare_b_dilepton_wilsons_from_rs_semileptonic(
        source,
        transition="b_s",
        lepton=lepton,
        matching_scale_gev=matching_scale_gev,
    )
    result = _evaluate_b_to_kstar_mumu(
        wilsons,
        mode=mode,
        q2_min_gev2=q2_min_gev2,
        q2_max_gev2=q2_max_gev2,
        inputs=inputs,
    )
    diagnostics = dict(result.diagnostics)
    diagnostics.update(rare_b_rs_semileptonic_vector_diagnostics(coeff))
    diagnostics["matching_assumption"] = coeff.matching_assumption
    return replace(result, diagnostics=diagnostics)


def bzero_kstarzero_mumu_from_rs_semileptonic_wilsons(
    source: RSSemileptonicWilsonBundle,
    *,
    lepton: str = "mu",
    q2_min_gev2: float | None = None,
    q2_max_gev2: float | None = None,
    matching_scale_gev: float | None = None,
    inputs: RareBToKStarDileptonInputs | None = None,
) -> RareBToKStarDileptonBranchingResult:
    """Evaluate ``BR(B0 -> K*(892)0 l+ l-)`` from Phase-3a RS Wilsons."""

    return rare_b_to_kstar_mumu_from_rs_semileptonic_wilsons(
        source,
        lepton=lepton,
        mode="bzero_kstarzero",
        q2_min_gev2=q2_min_gev2,
        q2_max_gev2=q2_max_gev2,
        matching_scale_gev=matching_scale_gev,
        inputs=inputs,
    )
