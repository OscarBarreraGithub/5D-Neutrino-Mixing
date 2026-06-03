"""Adapter over baryonic :mod:`quarkConstraints.rare_b_baryon_dilepton`.

This is the catalog boundary for ``Lambda_b -> Lambda mu+ mu-``.  Constraint
modules import this adapter only; the underlying high-q2 baryonic form-factor
integral and the reused rare_b_dilepton C9/C10 RS proxy remain isolated in
``quarkConstraints``.
"""

from __future__ import annotations

from dataclasses import replace

from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.rare_b_dilepton import RareBDileptonWilsonCoefficients
from quarkConstraints.rare_b_baryon_dilepton import (
    RARE_B_BARYONIC_DILEPTON_FORM_FACTOR_BUNDLE_V1,
    RARE_B_BARYONIC_DILEPTON_FORM_FACTOR_CITATION,
    RARE_B_BARYONIC_DILEPTON_INPUT_BUNDLE_V1,
    RARE_B_BARYONIC_DILEPTON_LIMITATION_V1,
    RARE_B_BARYONIC_DILEPTON_MODEL_V1,
    RARE_B_BARYONIC_DILEPTON_PROXY_THEORY_UNCERTAINTY_FRACTION,
    RARE_B_BARYONIC_DILEPTON_PROXY_THEORY_UNCERTAINTY_RATIONALE,
    RareBBaryonDileptonBranchingResult,
    RareBBaryonDileptonInputs,
    RareBBaryonFormFactorInputs,
    default_lambdab_to_lambda_dilepton_inputs as _default_lambdab_to_lambda_dilepton_inputs,
    evaluate_lambdab_to_lambda_mumu as _evaluate_lambdab_to_lambda_mumu,
    lambdab_to_lambda_fplus_fminus as _lambdab_to_lambda_fplus_fminus,
    sm_lambdab_to_lambda_mumu_branching_fraction as _sm_lambdab_to_lambda_mumu_branching_fraction,
)
from quarkConstraints.rs_semileptonic_wilsons import RSSemileptonicWilsonBundle

from .rare_b_meson import (
    rare_b_dilepton_wilsons_from_rs_semileptonic,
    rare_b_rs_semileptonic_vector_diagnostics,
)

__all__ = [
    "QuarkMassBasisCouplings",
    "RARE_B_BARYONIC_DILEPTON_MODEL_V1",
    "RARE_B_BARYONIC_DILEPTON_INPUT_BUNDLE_V1",
    "RARE_B_BARYONIC_DILEPTON_FORM_FACTOR_BUNDLE_V1",
    "RARE_B_BARYONIC_DILEPTON_FORM_FACTOR_CITATION",
    "RARE_B_BARYONIC_DILEPTON_LIMITATION_V1",
    "RARE_B_BARYONIC_DILEPTON_PROXY_THEORY_UNCERTAINTY_FRACTION",
    "RARE_B_BARYONIC_DILEPTON_PROXY_THEORY_UNCERTAINTY_RATIONALE",
    "RareBBaryonFormFactorInputs",
    "RareBBaryonDileptonInputs",
    "RareBBaryonDileptonBranchingResult",
    "RareBDileptonWilsonCoefficients",
    "rare_b_baryon_dilepton_default_inputs",
    "rare_b_baryon_fplus_fminus",
    "lambdab_lambda_mumu_sm_branching_fraction",
    "lambdab_lambda_mumu_branching_fraction",
    "lambdab_lambda_mumu_from_couplings",
    "lambdab_lambda_mumu_from_rs_semileptonic_wilsons",
]


def rare_b_baryon_dilepton_default_inputs() -> RareBBaryonDileptonInputs:
    """Return the default baryonic rare-beauty dilepton input bundle."""

    return _default_lambdab_to_lambda_dilepton_inputs()


def rare_b_baryon_fplus_fminus(
    q2_gev2: float,
    form_factor: RareBBaryonFormFactorInputs | None = None,
) -> tuple[float, float]:
    """Evaluate adapter-exposed ``Lambda_b -> Lambda`` ``F_+`` and ``F_-``."""

    return _lambdab_to_lambda_fplus_fminus(q2_gev2, form_factor)


def lambdab_lambda_mumu_sm_branching_fraction(
    *,
    q2_min_gev2: float | None = None,
    q2_max_gev2: float | None = None,
    inputs: RareBBaryonDileptonInputs | None = None,
) -> RareBBaryonDileptonBranchingResult:
    """Evaluate the SM-limit partial ``Lambda_b -> Lambda mu+ mu-`` BR."""

    return _sm_lambdab_to_lambda_mumu_branching_fraction(
        q2_min_gev2=q2_min_gev2,
        q2_max_gev2=q2_max_gev2,
        inputs=inputs,
    )


def lambdab_lambda_mumu_branching_fraction(
    source: QuarkMassBasisCouplings | RareBDileptonWilsonCoefficients | None = None,
    *,
    q2_min_gev2: float | None = None,
    q2_max_gev2: float | None = None,
    m_kk_gev: float | None = None,
    inputs: RareBBaryonDileptonInputs | None = None,
) -> RareBBaryonDileptonBranchingResult:
    """Evaluate baryonic ``BR(Lambda_b -> Lambda mu+ mu-)`` from source inputs."""

    return _evaluate_lambdab_to_lambda_mumu(
        source,
        q2_min_gev2=q2_min_gev2,
        q2_max_gev2=q2_max_gev2,
        m_kk_gev=m_kk_gev,
        inputs=inputs,
    )


def lambdab_lambda_mumu_from_couplings(
    couplings: QuarkMassBasisCouplings,
    *,
    q2_min_gev2: float | None = None,
    q2_max_gev2: float | None = None,
    m_kk_gev: float | None = None,
    inputs: RareBBaryonDileptonInputs | None = None,
) -> RareBBaryonDileptonBranchingResult:
    """Evaluate ``BR(Lambda_b -> Lambda mu+ mu-)`` from mass-basis couplings."""

    return lambdab_lambda_mumu_branching_fraction(
        couplings,
        q2_min_gev2=q2_min_gev2,
        q2_max_gev2=q2_max_gev2,
        m_kk_gev=m_kk_gev,
        inputs=inputs,
    )


def lambdab_lambda_mumu_from_rs_semileptonic_wilsons(
    source: RSSemileptonicWilsonBundle,
    *,
    lepton: str = "mu",
    q2_min_gev2: float | None = None,
    q2_max_gev2: float | None = None,
    matching_scale_gev: float | None = None,
    inputs: RareBBaryonDileptonInputs | None = None,
) -> RareBBaryonDileptonBranchingResult:
    """Evaluate ``BR(Lambda_b -> Lambda l+ l-)`` from Phase-3a RS Wilsons."""

    coeff = source.b_to_s_ll[lepton]
    wilsons = rare_b_dilepton_wilsons_from_rs_semileptonic(
        source,
        transition="b_s",
        lepton=lepton,
        matching_scale_gev=matching_scale_gev,
    )
    result = _evaluate_lambdab_to_lambda_mumu(
        wilsons,
        q2_min_gev2=q2_min_gev2,
        q2_max_gev2=q2_max_gev2,
        inputs=inputs,
    )
    diagnostics = dict(result.diagnostics)
    diagnostics.update(rare_b_rs_semileptonic_vector_diagnostics(coeff))
    diagnostics["matching_assumption"] = coeff.matching_assumption
    return replace(result, diagnostics=diagnostics)
