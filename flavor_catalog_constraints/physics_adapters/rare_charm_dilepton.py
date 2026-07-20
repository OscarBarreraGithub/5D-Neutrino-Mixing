"""Adapter over :mod:`quarkConstraints.rare_charm_dilepton`.

This is the catalog boundary for shared ``c -> u l+ l-`` rare-charm
machinery.  Constraint modules import this adapter only; the underlying
Hamiltonian convention and D0 leptonic short-distance formula remain isolated
in ``quarkConstraints``.  Phase-3a RS semileptonic Wilson bundles are
translated here into the existing rare-charm Wilson dataclass at the
Wilson-value consumption point.
"""

from __future__ import annotations

import math
from dataclasses import replace

from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.rare_charm_dilepton import (
    RARE_CHARM_DILEPTON_INPUT_BUNDLE_V1,
    RARE_CHARM_DILEPTON_MODEL_V1,
    RARE_CHARM_DILEPTON_OPERATOR_CONVENTION,
    RARE_CHARM_DILEPTON_RS_MATCHING_ASSUMPTION_V1,
    RareCharmDileptonCKMFactors,
    RareCharmDileptonSMInputs,
    RareCharmDileptonWilsonCoefficients,
    RareCharmLeptonicBranchingResult,
    RareCharmLeptonInputs,
    RareCharmMesonInputs,
)
from quarkConstraints.rare_charm_dilepton import (
    ckm_factors as _ckm_factors,
)
from quarkConstraints.rare_charm_dilepton import (
    compute_rare_charm_dilepton_wilsons as _compute_rare_charm_dilepton_wilsons,
)
from quarkConstraints.rare_charm_dilepton import (
    default_sm_inputs as _default_sm_inputs,
)
from quarkConstraints.rare_charm_dilepton import (
    evaluate_d0_to_ll as _evaluate_d0_to_ll,
)
from quarkConstraints.rare_charm_dilepton import (
    sm_branching_fraction as _sm_branching_fraction,
)
from quarkConstraints.rs_semileptonic_wilsons import (
    RS_SEMILEPTONIC_MATCHING_ASSUMPTION_V1,
    RSSemileptonicWilsonBundle,
    RSSemileptonicWilsonCoefficients,
)

RARE_CHARM_RS_SEMILEPTONIC_VECTOR_MATCHING_STATUS_V1 = (
    "Phase-3a light-Z RS semileptonic c->u C9/C10/C9p/C10p Wilsons "
    "consumed additively at the rare-charm Wilson-value input; no "
    "compute_rare_charm_dilepton_wilsons proxy and no second 1/M_KK^2 factor."
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
    "rare_charm_dilepton_wilsons_from_rs_semileptonic",
    "rare_charm_rs_semileptonic_coeff",
    "rare_charm_rs_semileptonic_vector_diagnostics",
    "d0_ll_from_couplings",
    "d0_mumu_from_couplings",
    "d0_ee_from_couplings",
    "d0_ll_from_rs_semileptonic_wilsons",
    "d0_mumu_from_rs_semileptonic_wilsons",
    "d0_ee_from_rs_semileptonic_wilsons",
    "RARE_CHARM_RS_SEMILEPTONIC_VECTOR_MATCHING_STATUS_V1",
    "RS_SEMILEPTONIC_MATCHING_ASSUMPTION_V1",
]


def _diagnostic_matching_scale(matching_scale_gev: float | None) -> float:
    if matching_scale_gev is None:
        return 0.0
    number = float(matching_scale_gev)
    if not math.isfinite(number) or number <= 0.0:
        raise ValueError("matching_scale_gev must be positive and finite")
    return number


def rare_charm_rs_semileptonic_coeff(
    source: RSSemileptonicWilsonBundle,
    *,
    lepton: str,
) -> RSSemileptonicWilsonCoefficients:
    """Return one Phase-3a ``c_to_u_ll`` same-flavor charged-lepton block."""

    try:
        coeff = source.c_to_u_ll[lepton]
    except (AttributeError, KeyError, TypeError) as exc:
        raise ValueError(
            f"rs_semileptonic_wilsons.c_to_u_ll[{lepton!r}] is not available"
        ) from exc
    if coeff.transition_key != "c_u":
        raise ValueError(
            f"c_to_u_ll[{lepton!r}] transition_key={coeff.transition_key!r}, "
            "expected 'c_u'"
        )
    if (
        coeff.quark_sector != "u"
        or coeff.final_quark_index != 0
        or coeff.initial_quark_index != 1
    ):
        raise ValueError(
            f"c_to_u_ll[{lepton!r}] has inconsistent quark indices/sector for c->u"
        )
    return coeff


def rare_charm_rs_semileptonic_vector_diagnostics(
    coeff: RSSemileptonicWilsonCoefficients,
) -> dict[str, object]:
    """Return diagnostics for a Phase-3a rare-charm Wilson block."""

    return {
        "rs_semileptonic_model_label": coeff.model_label,
        "rs_semileptonic_operator_convention": coeff.operator_convention,
        "rs_semileptonic_matching_assumption": coeff.matching_assumption,
        "rs_semileptonic_wilsons_present": True,
        "rs_semileptonic_matching_status": (
            "rs_semileptonic_additive_no_second_1_over_M_KK_squared"
        ),
        "rs_semileptonic_vector_matching_status": (
            RARE_CHARM_RS_SEMILEPTONIC_VECTOR_MATCHING_STATUS_V1
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
        "rare_charm_proxy_reused": False,
        "c_to_u_ll_rs_semileptonic_rewired": True,
        "short_distance_vector_axial_matching_rigorous": True,
    }


def _tag_rs_result(result, coeff: RSSemileptonicWilsonCoefficients):
    diagnostics = dict(result.diagnostics)
    diagnostics.update(rare_charm_rs_semileptonic_vector_diagnostics(coeff))
    diagnostics["matching_assumption"] = coeff.matching_assumption
    return replace(result, diagnostics=diagnostics)


def rare_charm_dilepton_wilsons_from_rs_semileptonic(
    source: RSSemileptonicWilsonBundle,
    *,
    lepton: str = "mu",
    matching_scale_gev: float | None = None,
) -> RareCharmDileptonWilsonCoefficients:
    """Translate Phase-3a ``c -> u l l`` Wilsons to the rare-charm core shape."""

    coeff = rare_charm_rs_semileptonic_coeff(source, lepton=lepton)
    scale = _diagnostic_matching_scale(matching_scale_gev)
    left_vector_contact = complex(coeff.contact_LL + coeff.contact_LR)
    right_vector_contact = complex(coeff.contact_RL + coeff.contact_RR)
    return RareCharmDileptonWilsonCoefficients(
        model_label=coeff.model_label,
        operator_convention=coeff.operator_convention,
        matching_assumption=coeff.matching_assumption,
        transition_key="c_u",
        M_KK=scale,
        matching_scale=scale,
        lambda_b=complex(coeff.lambda_ckm),
        left_uc_coupling=left_vector_contact,
        right_uc_coupling=right_vector_contact,
        left_uc_overlap=0.0j,
        right_uc_overlap=0.0j,
        left_quark_delta=left_vector_contact,
        right_quark_delta=right_vector_contact,
        lepton_left_delta=0.0,
        lepton_right_delta=0.0,
        lepton_vector_delta=0.0,
        lepton_axial_delta=0.0,
        c9_np=complex(coeff.c9_np),
        c10_np=complex(coeff.c10_np),
        c9p_np=complex(coeff.c9p_np),
        c10p_np=complex(coeff.c10p_np),
    )


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


def d0_ll_from_rs_semileptonic_wilsons(
    source: RSSemileptonicWilsonBundle,
    *,
    lepton: str,
    matching_scale_gev: float | None = None,
    inputs: RareCharmDileptonSMInputs | None = None,
) -> RareCharmLeptonicBranchingResult:
    """Evaluate short-distance ``BR(D0 -> l+ l-)`` from Phase-3a RS Wilsons."""

    coeff = rare_charm_rs_semileptonic_coeff(source, lepton=lepton)
    wilsons = rare_charm_dilepton_wilsons_from_rs_semileptonic(
        source,
        lepton=lepton,
        matching_scale_gev=matching_scale_gev,
    )
    result = _evaluate_d0_to_ll(
        wilsons,
        lepton=lepton,
        inputs=inputs,
    )
    return _tag_rs_result(result, coeff)


def d0_mumu_from_rs_semileptonic_wilsons(
    source: RSSemileptonicWilsonBundle,
    *,
    matching_scale_gev: float | None = None,
    inputs: RareCharmDileptonSMInputs | None = None,
) -> RareCharmLeptonicBranchingResult:
    """Evaluate short-distance ``BR(D0 -> mu+ mu-)`` from Phase-3a RS Wilsons."""

    return d0_ll_from_rs_semileptonic_wilsons(
        source,
        lepton="mu",
        matching_scale_gev=matching_scale_gev,
        inputs=inputs,
    )


def d0_ee_from_rs_semileptonic_wilsons(
    source: RSSemileptonicWilsonBundle,
    *,
    matching_scale_gev: float | None = None,
    inputs: RareCharmDileptonSMInputs | None = None,
) -> RareCharmLeptonicBranchingResult:
    """Evaluate short-distance ``BR(D0 -> e+ e-)`` from Phase-3a RS Wilsons."""

    return d0_ll_from_rs_semileptonic_wilsons(
        source,
        lepton="e",
        matching_scale_gev=matching_scale_gev,
        inputs=inputs,
    )
