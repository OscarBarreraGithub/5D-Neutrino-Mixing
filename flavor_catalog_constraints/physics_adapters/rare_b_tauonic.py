"""Tau-mode adapter over :mod:`quarkConstraints.rare_b_dilepton`.

The shared rare-beauty core was introduced for B005/B006 muon modes, but its
leptonic formula is parameterized by the charged-lepton mass in the SM input
bundle.  This adapter reuses that core for ``B_q -> tau+ tau-`` by supplying
the tau mass while leaving the muon and electron entry points untouched.
"""

from __future__ import annotations

from dataclasses import replace

from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.rare_b_dilepton import (
    RARE_B_DILEPTON_MODEL_V1,
    RARE_B_DILEPTON_OPERATOR_CONVENTION,
    RARE_B_DILEPTON_RS_MATCHING_ASSUMPTION_V1,
    RareBDileptonCKMFactors,
    RareBDileptonSMInputs,
    RareBDileptonWilsonCoefficients,
    RareBLeptonicBranchingResult,
    ckm_factors as _ckm_factors,
    default_sm_inputs as _default_sm_inputs,
    evaluate_bq_to_mumu as _evaluate_bq_to_mumu,
    sm_branching_fraction as _sm_branching_fraction,
)
from quarkConstraints.rs_semileptonic_wilsons import RSSemileptonicWilsonBundle

from .rare_b_meson import (
    bq_mumu_from_rs_semileptonic_wilsons as _bq_mumu_from_rs_semileptonic_wilsons,
)

__all__ = [
    "TAU_MASS_GEV",
    "RARE_B_TAUONIC_INPUT_BUNDLE_V1",
    "RARE_B_TAUONIC_MODEL_NOTE_V1",
    "RARE_B_DILEPTON_MODEL_V1",
    "RARE_B_DILEPTON_OPERATOR_CONVENTION",
    "RARE_B_DILEPTON_RS_MATCHING_ASSUMPTION_V1",
    "QuarkMassBasisCouplings",
    "RareBDileptonCKMFactors",
    "RareBDileptonSMInputs",
    "RareBDileptonWilsonCoefficients",
    "RareBLeptonicBranchingResult",
    "rare_b_tauonic_default_sm_inputs",
    "rare_b_tauonic_ckm_factors",
    "rare_b_tauonic_sm_branching_fraction",
    "bq_tautau_from_couplings",
    "bs_tautau_from_couplings",
    "bd_tautau_from_couplings",
    "bq_tautau_from_rs_semileptonic_wilsons",
    "bs_tautau_from_rs_semileptonic_wilsons",
    "bd_tautau_from_rs_semileptonic_wilsons",
]

TAU_MASS_GEV = 1.77686
RARE_B_TAUONIC_INPUT_BUNDLE_V1 = "rare_b_dilepton_tau_inputs_buras_repo_ckm_v1"
RARE_B_TAUONIC_MODEL_NOTE_V1 = (
    "B_q -> tau+tau- reuses the rare_b_dilepton Buras C10 formula with "
    "muon_mass_gev replaced by the tau mass; the field name is inherited "
    "from the B005/B006 core and acts as the charged-lepton mass parameter."
)


def rare_b_tauonic_default_sm_inputs() -> RareBDileptonSMInputs:
    """Return the rare-beauty dilepton input bundle specialized to taus."""
    base = _default_sm_inputs()
    return replace(
        base,
        input_bundle=RARE_B_TAUONIC_INPUT_BUNDLE_V1,
        muon_mass_gev=TAU_MASS_GEV,
    )


def rare_b_tauonic_ckm_factors(
    transition: str = "b_s",
    inputs: RareBDileptonSMInputs | None = None,
) -> RareBDileptonCKMFactors:
    """Return the CKM factors used by the tau-mode rare-beauty core."""
    return _ckm_factors(transition, rare_b_tauonic_default_sm_inputs() if inputs is None else inputs)


def rare_b_tauonic_sm_branching_fraction(
    transition: str = "b_s",
    inputs: RareBDileptonSMInputs | None = None,
) -> RareBLeptonicBranchingResult:
    """Evaluate the SM-limit ``B_q -> tau+ tau-`` branching fraction."""
    return _sm_branching_fraction(
        transition,
        rare_b_tauonic_default_sm_inputs() if inputs is None else inputs,
    )


def bq_tautau_from_couplings(
    couplings: QuarkMassBasisCouplings,
    *,
    transition: str,
    m_kk_gev: float | None = None,
    inputs: RareBDileptonSMInputs | None = None,
) -> RareBLeptonicBranchingResult:
    """Evaluate ``BR(B_q -> tau+ tau-)`` from mass-basis couplings."""
    return _evaluate_bq_to_mumu(
        couplings,
        transition=transition,
        m_kk_gev=m_kk_gev,
        inputs=rare_b_tauonic_default_sm_inputs() if inputs is None else inputs,
    )


def bs_tautau_from_couplings(
    couplings: QuarkMassBasisCouplings,
    *,
    m_kk_gev: float | None = None,
    inputs: RareBDileptonSMInputs | None = None,
) -> RareBLeptonicBranchingResult:
    """Evaluate ``BR(B_s -> tau+ tau-)`` from mass-basis couplings."""
    return bq_tautau_from_couplings(
        couplings,
        transition="b_s",
        m_kk_gev=m_kk_gev,
        inputs=inputs,
    )


def bd_tautau_from_couplings(
    couplings: QuarkMassBasisCouplings,
    *,
    m_kk_gev: float | None = None,
    inputs: RareBDileptonSMInputs | None = None,
) -> RareBLeptonicBranchingResult:
    """Evaluate ``BR(B0 -> tau+ tau-)`` from mass-basis couplings."""
    return bq_tautau_from_couplings(
        couplings,
        transition="b_d",
        m_kk_gev=m_kk_gev,
        inputs=inputs,
    )


def bq_tautau_from_rs_semileptonic_wilsons(
    source: RSSemileptonicWilsonBundle,
    *,
    transition: str,
    matching_scale_gev: float | None = None,
    inputs: RareBDileptonSMInputs | None = None,
) -> RareBLeptonicBranchingResult:
    """Evaluate ``BR(B_q -> tau+ tau-)`` from Phase-3a RS Wilsons."""

    return _bq_mumu_from_rs_semileptonic_wilsons(
        source,
        transition=transition,
        lepton="tau",
        matching_scale_gev=matching_scale_gev,
        inputs=rare_b_tauonic_default_sm_inputs() if inputs is None else inputs,
    )


def bs_tautau_from_rs_semileptonic_wilsons(
    source: RSSemileptonicWilsonBundle,
    *,
    matching_scale_gev: float | None = None,
    inputs: RareBDileptonSMInputs | None = None,
) -> RareBLeptonicBranchingResult:
    """Evaluate ``BR(B_s -> tau+ tau-)`` from Phase-3a RS Wilsons."""

    return bq_tautau_from_rs_semileptonic_wilsons(
        source,
        transition="b_s",
        matching_scale_gev=matching_scale_gev,
        inputs=inputs,
    )


def bd_tautau_from_rs_semileptonic_wilsons(
    source: RSSemileptonicWilsonBundle,
    *,
    matching_scale_gev: float | None = None,
    inputs: RareBDileptonSMInputs | None = None,
) -> RareBLeptonicBranchingResult:
    """Evaluate ``BR(B0 -> tau+ tau-)`` from Phase-3a RS Wilsons."""

    return bq_tautau_from_rs_semileptonic_wilsons(
        source,
        transition="b_d",
        matching_scale_gev=matching_scale_gev,
        inputs=inputs,
    )
