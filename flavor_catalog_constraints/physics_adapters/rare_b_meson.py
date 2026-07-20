"""Adapter over :mod:`quarkConstraints.rare_b_dilepton`.

This is the catalog boundary for shared ``b -> q l+ l-`` rare-beauty
machinery.  Constraint modules import this adapter only; the underlying
Hamiltonian convention and SM formula remain isolated in ``quarkConstraints``.
Phase-3a RS semileptonic Wilson bundles are translated here into the existing
rare-B Wilson dataclass at the Wilson-value consumption point.
"""

from __future__ import annotations

import math
from dataclasses import replace

from quarkConstraints import bsgamma as _bsgamma_core
from quarkConstraints import rare_b_dilepton as _rare_b_dilepton_core
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
    RareBDileptonCKMFactors,
    RareBDileptonMesonInputs,
    RareBDileptonSMInputs,
    RareBDileptonWilsonCoefficients,
    RareBInclusiveDileptonBranchingResult,
    RareBInclusiveDileptonInputs,
    RareBLeptonicBranchingResult,
    RareBToKDileptonBranchingResult,
    RareBToKDileptonInputs,
    RareBToKFormFactorInputs,
)
from quarkConstraints.rare_b_dilepton import (
    b_to_k_fplus as _b_to_k_fplus,
)
from quarkConstraints.rare_b_dilepton import (
    ckm_factors as _ckm_factors,
)
from quarkConstraints.rare_b_dilepton import (
    compute_rare_b_dilepton_wilsons as _compute_rare_b_dilepton_wilsons,
)
from quarkConstraints.rare_b_dilepton import (
    default_b_to_k_dilepton_inputs as _default_b_to_k_dilepton_inputs,
)
from quarkConstraints.rare_b_dilepton import (
    default_inclusive_b_to_xs_dilepton_inputs as _default_inclusive_b_to_xs_dilepton_inputs,
)
from quarkConstraints.rare_b_dilepton import (
    default_sm_inputs as _default_sm_inputs,
)
from quarkConstraints.rare_b_dilepton import (
    evaluate_b_to_k_mumu as _evaluate_b_to_k_mumu,
)
from quarkConstraints.rare_b_dilepton import (
    evaluate_bq_to_mumu as _evaluate_bq_to_mumu,
)
from quarkConstraints.rare_b_dilepton import (
    evaluate_inclusive_b_to_xs_mumu as _evaluate_inclusive_b_to_xs_mumu,
)
from quarkConstraints.rare_b_dilepton import (
    sm_b_to_k_mumu_branching_fraction as _sm_b_to_k_mumu_branching_fraction,
)
from quarkConstraints.rare_b_dilepton import (
    sm_branching_fraction as _sm_branching_fraction,
)
from quarkConstraints.rare_b_dilepton import (
    sm_inclusive_b_to_xs_mumu_branching_fraction as _sm_inclusive_b_to_xs_mumu_branching_fraction,
)
from quarkConstraints.rs_semileptonic_wilsons import (
    RS_SEMILEPTONIC_MATCHING_ASSUMPTION_V1,
    RSSemileptonicWilsonBundle,
    RSSemileptonicWilsonCoefficients,
)

RARE_B_RS_SEMILEPTONIC_VECTOR_MATCHING_STATUS_V1 = (
    "Phase-3a light-Z RS semileptonic C9/C10/C9p/C10p Wilsons consumed "
    "additively; no rare_b_dilepton _wilson_prefactor call and no second "
    "1/M_KK^2 factor."
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
    "rare_b_dilepton_wilsons_from_rs_semileptonic",
    "rare_b_rs_semileptonic_vector_diagnostics",
    "bq_mumu_from_couplings",
    "bs_mumu_from_couplings",
    "bd_mumu_from_couplings",
    "bq_mumu_from_rs_semileptonic_wilsons",
    "bs_mumu_from_rs_semileptonic_wilsons",
    "bd_mumu_from_rs_semileptonic_wilsons",
    "rare_b_to_k_dilepton_default_inputs",
    "rare_b_to_k_fplus",
    "rare_b_to_k_mumu_sm_branching_fraction",
    "rare_b_to_k_mumu_branching_fraction",
    "bplus_kplus_mumu_from_couplings",
    "bzero_kzero_mumu_from_couplings",
    "rare_b_to_k_mumu_from_rs_semileptonic_wilsons",
    "bplus_kplus_mumu_from_rs_semileptonic_wilsons",
    "bzero_kzero_mumu_from_rs_semileptonic_wilsons",
    "rare_b_inclusive_xs_dilepton_default_inputs",
    "rare_b_inclusive_xs_mumu_sm_branching_fraction",
    "rare_b_inclusive_xs_mumu_branching_fraction",
    "inclusive_b_to_xs_mumu_from_couplings",
    "inclusive_b_to_xs_mumu_from_rs_semileptonic_wilsons",
    "RARE_B_RS_SEMILEPTONIC_VECTOR_MATCHING_STATUS_V1",
    "RS_SEMILEPTONIC_MATCHING_ASSUMPTION_V1",
]

_TRANSITION_TO_RS_BLOCK = {
    "b_s": "b_to_s_ll",
    "b_d": "b_to_d_ll",
}


def _diagnostic_matching_scale(matching_scale_gev: float | None) -> float:
    if matching_scale_gev is None:
        return 0.0
    number = float(matching_scale_gev)
    if not math.isfinite(number) or number <= 0.0:
        raise ValueError("matching_scale_gev must be positive and finite")
    return number


def _rs_semileptonic_coeff(
    source: RSSemileptonicWilsonBundle,
    *,
    transition: str,
    lepton: str,
) -> RSSemileptonicWilsonCoefficients:
    try:
        block_name = _TRANSITION_TO_RS_BLOCK[transition]
    except KeyError as exc:
        raise ValueError(f"unsupported rare-B transition {transition!r}") from exc
    try:
        coeff = getattr(source, block_name)[lepton]
    except (AttributeError, KeyError, TypeError) as exc:
        raise ValueError(
            f"rs_semileptonic_wilsons.{block_name}[{lepton!r}] is not available"
        ) from exc
    if coeff.transition_key != transition:
        raise ValueError(
            f"{block_name}[{lepton!r}] transition_key={coeff.transition_key!r}, "
            f"expected {transition!r}"
        )
    return coeff


def rare_b_rs_semileptonic_vector_diagnostics(
    coeff: RSSemileptonicWilsonCoefficients,
) -> dict[str, object]:
    """Return diagnostics for a Phase-3a rare-B Wilson block."""

    return {
        "rs_semileptonic_model_label": coeff.model_label,
        "rs_semileptonic_operator_convention": coeff.operator_convention,
        "rs_semileptonic_matching_assumption": coeff.matching_assumption,
        "rs_semileptonic_wilsons_present": True,
        "rs_semileptonic_matching_status": (
            "rs_semileptonic_additive_no_second_1_over_M_KK_squared"
        ),
        "rs_semileptonic_vector_matching_status": (
            RARE_B_RS_SEMILEPTONIC_VECTOR_MATCHING_STATUS_V1
        ),
        "rs_semileptonic_transition_key": coeff.transition_key,
        "rs_semileptonic_lepton_key": coeff.lepton_key,
        "rs_semileptonic_quark_sector": coeff.quark_sector,
        "rs_semileptonic_final_quark_index": int(coeff.final_quark_index),
        "rs_semileptonic_initial_quark_index": int(coeff.initial_quark_index),
        "rs_semileptonic_lambda_ckm_name": coeff.lambda_ckm_name,
        "rs_semileptonic_lambda_ckm": complex(coeff.lambda_ckm),
        "rs_semileptonic_contact_units": coeff.contact_units,
        "rs_semileptonic_contacts": {
            key: complex(value) for key, value in coeff.contacts.items()
        },
        "rs_semileptonic_wilson_coefficients": {
            key: complex(value) for key, value in coeff.wilsons.items()
        },
        "wilson_prefactor_reused": False,
        "second_mkk_suppression_applied": False,
    }


def _tag_rs_result(result, coeff: RSSemileptonicWilsonCoefficients):
    diagnostics = dict(result.diagnostics)
    diagnostics.update(rare_b_rs_semileptonic_vector_diagnostics(coeff))
    diagnostics["matching_assumption"] = coeff.matching_assumption
    diagnostics["c9_c10_proxy_reused"] = False
    diagnostics["c9_c10_rs_semileptonic_rewired"] = True
    return replace(result, diagnostics=diagnostics)


def rare_b_dilepton_wilsons_from_rs_semileptonic(
    source: RSSemileptonicWilsonBundle,
    *,
    transition: str = "b_s",
    lepton: str = "mu",
    matching_scale_gev: float | None = None,
) -> RareBDileptonWilsonCoefficients:
    """Translate Phase-3a RS semileptonic Wilsons to the rare-B core shape."""

    coeff = _rs_semileptonic_coeff(source, transition=transition, lepton=lepton)
    scale = _diagnostic_matching_scale(matching_scale_gev)
    left_contact = complex(coeff.contact_LL + coeff.contact_LR)
    right_contact = complex(coeff.contact_RL + coeff.contact_RR)
    return RareBDileptonWilsonCoefficients(
        model_label=coeff.model_label,
        operator_convention=coeff.operator_convention,
        matching_assumption=coeff.matching_assumption,
        transition_key=transition,
        M_KK=scale,
        matching_scale=scale,
        lambda_t=complex(coeff.lambda_ckm),
        left_qb_coupling=left_contact,
        right_qb_coupling=right_contact,
        left_qb_overlap=0.0j,
        right_qb_overlap=0.0j,
        left_quark_delta=left_contact,
        right_quark_delta=right_contact,
        muon_left_delta=0.0,
        muon_right_delta=0.0,
        muon_vector_delta=0.0,
        muon_axial_delta=0.0,
        c9_np=complex(coeff.c9_np),
        c10_np=complex(coeff.c10_np),
        c9p_np=complex(coeff.c9p_np),
        c10p_np=complex(coeff.c10p_np),
    )


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


def bq_mumu_from_rs_semileptonic_wilsons(
    source: RSSemileptonicWilsonBundle,
    *,
    transition: str,
    lepton: str = "mu",
    matching_scale_gev: float | None = None,
    inputs: RareBDileptonSMInputs | None = None,
) -> RareBLeptonicBranchingResult:
    """Evaluate ``BR(B_q -> l+ l-)`` from Phase-3a RS C9/C10 Wilsons."""

    coeff = _rs_semileptonic_coeff(source, transition=transition, lepton=lepton)
    wilsons = rare_b_dilepton_wilsons_from_rs_semileptonic(
        source,
        transition=transition,
        lepton=lepton,
        matching_scale_gev=matching_scale_gev,
    )
    result = _evaluate_bq_to_mumu(
        wilsons,
        transition=transition,
        inputs=inputs,
    )
    return _tag_rs_result(result, coeff)


def bs_mumu_from_rs_semileptonic_wilsons(
    source: RSSemileptonicWilsonBundle,
    *,
    lepton: str = "mu",
    matching_scale_gev: float | None = None,
    inputs: RareBDileptonSMInputs | None = None,
) -> RareBLeptonicBranchingResult:
    """Evaluate ``BR(B_s -> l+ l-)`` from Phase-3a RS C9/C10 Wilsons."""

    return bq_mumu_from_rs_semileptonic_wilsons(
        source,
        transition="b_s",
        lepton=lepton,
        matching_scale_gev=matching_scale_gev,
        inputs=inputs,
    )


def bd_mumu_from_rs_semileptonic_wilsons(
    source: RSSemileptonicWilsonBundle,
    *,
    lepton: str = "mu",
    matching_scale_gev: float | None = None,
    inputs: RareBDileptonSMInputs | None = None,
) -> RareBLeptonicBranchingResult:
    """Evaluate ``BR(B_d -> l+ l-)`` from Phase-3a RS C9/C10 Wilsons."""

    return bq_mumu_from_rs_semileptonic_wilsons(
        source,
        transition="b_d",
        lepton=lepton,
        matching_scale_gev=matching_scale_gev,
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


def rare_b_to_k_mumu_from_rs_semileptonic_wilsons(
    source: RSSemileptonicWilsonBundle,
    *,
    lepton: str = "mu",
    mode: str = "bplus_kplus",
    q2_min_gev2: float | None = None,
    q2_max_gev2: float | None = None,
    matching_scale_gev: float | None = None,
    inputs: RareBToKDileptonInputs | None = None,
) -> RareBToKDileptonBranchingResult:
    """Evaluate exclusive ``BR(B -> K l+ l-)`` from Phase-3a RS Wilsons."""

    coeff = _rs_semileptonic_coeff(source, transition="b_s", lepton=lepton)
    wilsons = rare_b_dilepton_wilsons_from_rs_semileptonic(
        source,
        transition="b_s",
        lepton=lepton,
        matching_scale_gev=matching_scale_gev,
    )
    result = _evaluate_b_to_k_mumu(
        wilsons,
        mode=mode,
        q2_min_gev2=q2_min_gev2,
        q2_max_gev2=q2_max_gev2,
        inputs=inputs,
    )
    return _tag_rs_result(result, coeff)


def bplus_kplus_mumu_from_rs_semileptonic_wilsons(
    source: RSSemileptonicWilsonBundle,
    *,
    lepton: str = "mu",
    q2_min_gev2: float | None = None,
    q2_max_gev2: float | None = None,
    matching_scale_gev: float | None = None,
    inputs: RareBToKDileptonInputs | None = None,
) -> RareBToKDileptonBranchingResult:
    """Evaluate ``BR(B+ -> K+ l+ l-)`` from Phase-3a RS Wilsons."""

    return rare_b_to_k_mumu_from_rs_semileptonic_wilsons(
        source,
        lepton=lepton,
        mode="bplus_kplus",
        q2_min_gev2=q2_min_gev2,
        q2_max_gev2=q2_max_gev2,
        matching_scale_gev=matching_scale_gev,
        inputs=inputs,
    )


def bzero_kzero_mumu_from_rs_semileptonic_wilsons(
    source: RSSemileptonicWilsonBundle,
    *,
    lepton: str = "mu",
    q2_min_gev2: float | None = None,
    q2_max_gev2: float | None = None,
    matching_scale_gev: float | None = None,
    inputs: RareBToKDileptonInputs | None = None,
) -> RareBToKDileptonBranchingResult:
    """Evaluate ``BR(B0 -> K0 l+ l-)`` from Phase-3a RS Wilsons."""

    return rare_b_to_k_mumu_from_rs_semileptonic_wilsons(
        source,
        lepton=lepton,
        mode="bzero_kzero",
        q2_min_gev2=q2_min_gev2,
        q2_max_gev2=q2_max_gev2,
        matching_scale_gev=matching_scale_gev,
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


def inclusive_b_to_xs_mumu_from_rs_semileptonic_wilsons(
    source: RSSemileptonicWilsonBundle,
    *,
    lepton: str = "mu",
    sm_branching_fraction: float,
    q2_min_gev2: float,
    q2_max_gev2: float,
    matching_scale_gev: float | None = None,
    dipole_couplings: QuarkMassBasisCouplings | None = None,
    inputs: RareBInclusiveDileptonInputs | None = None,
) -> RareBInclusiveDileptonBranchingResult:
    """Evaluate inclusive ``B -> X_s l l`` from RS C9/C10 and optional C7."""

    coeff = _rs_semileptonic_coeff(source, transition="b_s", lepton=lepton)
    dilepton_wilsons = rare_b_dilepton_wilsons_from_rs_semileptonic(
        source,
        transition="b_s",
        lepton=lepton,
        matching_scale_gev=matching_scale_gev,
    )
    if dipole_couplings is None:
        result = _evaluate_inclusive_b_to_xs_mumu(
            dilepton_wilsons,
            sm_branching_fraction=sm_branching_fraction,
            q2_min_gev2=q2_min_gev2,
            q2_max_gev2=q2_max_gev2,
            inputs=inputs,
        )
        return _tag_rs_result(result, coeff)

    p = _default_inclusive_b_to_xs_dilepton_inputs() if inputs is None else inputs
    sm_br = _rare_b_dilepton_core._positive_float(  # noqa: SLF001
        sm_branching_fraction,
        "sm_branching_fraction",
    )
    q2_min, q2_max = _rare_b_dilepton_core._inclusive_xs_q2_bounds(  # noqa: SLF001
        p,
        q2_min_gev2,
        q2_max_gev2,
    )
    factors = _ckm_factors("b_s", p.short_distance_inputs)
    dipole_wilsons = _bsgamma_core.compute_bsgamma_wilsons(
        dipole_couplings,
        m_kk_gev=matching_scale_gev,
        inputs=p.dipole_inputs,
    )

    c9_np = complex(dilepton_wilsons.c9_np)
    c9p_np = complex(dilepton_wilsons.c9p_np)
    c10_np = complex(dilepton_wilsons.c10_np)
    c10p_np = complex(dilepton_wilsons.c10p_np)
    c7_np = complex(dipole_wilsons.c7_np)
    c7p_np = complex(dipole_wilsons.c7p_np)

    c7_total = complex(p.dipole_inputs.c7_sm_eff + c7_np)
    c7p_total = complex(p.dipole_inputs.c7p_sm_eff + c7p_np)
    c9_total = complex(p.c9_sm + c9_np)
    c9p_total = complex(c9p_np)
    c10_total = complex(p.short_distance_inputs.c10_sm + c10_np)
    c10p_total = complex(c10p_np)

    sm_integral = _rare_b_dilepton_core._inclusive_xs_shape_integral(  # noqa: SLF001
        inputs=p,
        q2_min_gev2=q2_min,
        q2_max_gev2=q2_max,
        c7_total=p.dipole_inputs.c7_sm_eff,
        c7p_total=p.dipole_inputs.c7p_sm_eff,
        c9_total=p.c9_sm,
        c9p_total=0.0j,
        c10_total=p.short_distance_inputs.c10_sm,
        c10p_total=0.0j,
    )
    total_integral = _rare_b_dilepton_core._inclusive_xs_shape_integral(  # noqa: SLF001
        inputs=p,
        q2_min_gev2=q2_min,
        q2_max_gev2=q2_max,
        c7_total=c7_total,
        c7p_total=c7p_total,
        c9_total=c9_total,
        c9p_total=c9p_total,
        c10_total=c10_total,
        c10p_total=c10p_total,
    )
    if sm_integral <= 0.0 or not math.isfinite(sm_integral):
        raise ValueError("SM inclusive B -> X_s ll shape integral must be positive")
    if total_integral < 0.0 or not math.isfinite(total_integral):
        raise ValueError("total inclusive B -> X_s ll shape integral must be finite")

    ratio_to_sm = float(total_integral / sm_integral)
    br = float(sm_br * ratio_to_sm)
    q2_mid = 0.5 * (q2_min + q2_max)
    diagnostics = {
        "lambda_wolfenstein": float(factors.lambda_wolfenstein),
        "lambda_t": complex(factors.lambda_t),
        "q2_min_gev2": float(q2_min),
        "q2_max_gev2": float(q2_max),
        "q2_bin_width_gev2": float(q2_max - q2_min),
        "q2_mid_gev2": float(q2_mid),
        "partonic_b_mass_gev": float(p.partonic_b_mass_gev),
        "c7_sm_eff": complex(p.dipole_inputs.c7_sm_eff),
        "c7p_sm_eff": complex(p.dipole_inputs.c7p_sm_eff),
        "c9_sm": float(p.c9_sm),
        "c10_sm": float(p.short_distance_inputs.c10_sm),
        "c7_total": complex(c7_total),
        "c7p_total": complex(c7p_total),
        "c9_total": complex(c9_total),
        "c9p_total": complex(c9p_total),
        "c10_total": complex(c10_total),
        "c10p_total": complex(c10p_total),
        "c7_np": complex(c7_np),
        "c7p_np": complex(c7p_np),
        "c9_np": complex(c9_np),
        "c9p_np": complex(c9p_np),
        "c10_np": complex(c10_np),
        "c10p_np": complex(c10p_np),
        "sm_shape_integral": float(sm_integral),
        "total_shape_integral": float(total_integral),
        "kernel_sm_at_bin_center": (
            _rare_b_dilepton_core._inclusive_xs_partonic_kernel(  # noqa: SLF001
                q2_mid,
                inputs=p,
                c7_total=p.dipole_inputs.c7_sm_eff,
                c7p_total=p.dipole_inputs.c7p_sm_eff,
                c9_total=p.c9_sm,
                c9p_total=0.0j,
                c10_total=p.short_distance_inputs.c10_sm,
                c10p_total=0.0j,
            )
        ),
        "kernel_total_at_bin_center": (
            _rare_b_dilepton_core._inclusive_xs_partonic_kernel(  # noqa: SLF001
                q2_mid,
                inputs=p,
                c7_total=c7_total,
                c7p_total=c7p_total,
                c9_total=c9_total,
                c9p_total=c9p_total,
                c10_total=c10_total,
                c10p_total=c10p_total,
            )
        ),
        "integration_steps": float(p.integration_steps),
        "constants_citation": p.constants_citation,
        "matching_assumption": coeff.matching_assumption,
        "dipole_matching_assumption": _bsgamma_core.BSGAMMA_RS_MATCHING_ASSUMPTION_V1,
        "inclusive_limitations": RARE_B_DILEPTON_INCLUSIVE_XS_LIMITATION_V1,
        "c7_included": True,
        "c9_c10_proxy_reused": False,
        "c9_c10_rs_semileptonic_rewired": True,
        "m_kk_gev": float(dilepton_wilsons.M_KK),
        "matching_scale_gev": float(dilepton_wilsons.matching_scale),
        "dilepton_wilson_coefficients": {
            key: complex(value) for key, value in dilepton_wilsons.wilsons.items()
        },
        "c7_np_matching": complex(dipole_wilsons.c7_np_matching),
        "c7p_np_matching": complex(dipole_wilsons.c7p_np_matching),
        "c8_np_matching": complex(dipole_wilsons.c8_np_matching),
        "c8p_np_matching": complex(dipole_wilsons.c8p_np_matching),
        "c8_np": complex(dipole_wilsons.c8_np),
        "c8p_np": complex(dipole_wilsons.c8p_np),
        "c7_running_from_c7": float(dipole_wilsons.c7_running_from_c7),
        "c7_running_from_c8": float(dipole_wilsons.c7_running_from_c8),
        "c8_running_from_c8": float(dipole_wilsons.c8_running_from_c8),
        "alpha_s_matching_scale": float(dipole_wilsons.alpha_s_matching_scale),
        "alpha_s_low_scale": float(dipole_wilsons.alpha_s_low_scale),
        "dipole_wilson_coefficients": {
            key: complex(value) for key, value in dipole_wilsons.wilsons.items()
        },
    }
    diagnostics.update(rare_b_rs_semileptonic_vector_diagnostics(coeff))

    return RareBInclusiveDileptonBranchingResult(
        model_label=RARE_B_DILEPTON_INCLUSIVE_XS_MODEL_V1,
        input_bundle=p.input_bundle,
        q2_min_gev2=float(q2_min),
        q2_max_gev2=float(q2_max),
        branching_fraction=br,
        sm_branching_fraction=float(sm_br),
        np_shift_branching_fraction=float(br - sm_br),
        ratio_to_sm=ratio_to_sm,
        c7_total=complex(c7_total),
        c7p_total=complex(c7p_total),
        c9_total=complex(c9_total),
        c9p_total=complex(c9p_total),
        c10_total=complex(c10_total),
        c10p_total=complex(c10p_total),
        c7_np=complex(c7_np),
        c7p_np=complex(c7p_np),
        c9_np=complex(c9_np),
        c9p_np=complex(c9p_np),
        c10_np=complex(c10_np),
        c10p_np=complex(c10p_np),
        lambda_t=complex(factors.lambda_t),
        dilepton_wilsons=dilepton_wilsons,
        dipole_wilsons=dipole_wilsons,
        diagnostics=diagnostics,
    )


# Append-only export for constraints that need to reuse the B016 proxy-budget
# convention without importing quarkConstraints directly.
RARE_B_DILEPTON_EXCLUSIVE_BK_PROXY_THEORY_UNCERTAINTY_FRACTION = (
    _rare_b_dilepton_core.RARE_B_DILEPTON_EXCLUSIVE_BK_PROXY_THEORY_UNCERTAINTY_FRACTION
)
RARE_B_DILEPTON_EXCLUSIVE_BK_PROXY_THEORY_UNCERTAINTY_RATIONALE = (
    _rare_b_dilepton_core.RARE_B_DILEPTON_EXCLUSIVE_BK_PROXY_THEORY_UNCERTAINTY_RATIONALE
)

__all__.extend(
    [
        "RARE_B_DILEPTON_EXCLUSIVE_BK_PROXY_THEORY_UNCERTAINTY_FRACTION",
        "RARE_B_DILEPTON_EXCLUSIVE_BK_PROXY_THEORY_UNCERTAINTY_RATIONALE",
    ]
)
