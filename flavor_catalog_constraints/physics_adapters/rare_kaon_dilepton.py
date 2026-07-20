"""Adapter over :mod:`quarkConstraints.rare_kaon_dilepton`.

This is the catalog boundary for Delta-S=1 ``s -> d l+l-`` short-distance
machinery.  Constraint modules import this adapter only; the underlying physics
implementation remains isolated in ``quarkConstraints``.  The module starts
with the ``K_L -> mu+ mu-`` short-distance observable and is designed as the
shared place for K008/K009/K010/K012 dilepton rare-kaon reuse.
"""

from __future__ import annotations

import math
from dataclasses import replace

from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.rare_kaon_dilepton import (
    RARE_KAON_DILEPTON_INPUT_BUNDLE_V1,
    RARE_KAON_DILEPTON_MODEL_V1,
    RARE_KAON_DILEPTON_OPERATOR_CONVENTION,
    RARE_KAON_DILEPTON_PARAMETRIZATION_CITATION,
    RARE_KAON_DILEPTON_RS_MATCHING_ASSUMPTION_V1,
    RARE_KAON_PI0EE_INTERFERENCE_LIMITATION_V1,
    RARE_KAON_PI0EE_MODEL_V1,
    RARE_KAON_PI0EE_PARAMETRIZATION_CITATION,
    KLongMuMuShortDistanceResult,
    KLongPi0EEChPTInputs,
    KLongPi0EEResult,
    RareKaonDileptonCKMFactors,
    RareKaonDileptonSMInputs,
    RareKaonDileptonWilsonCoefficients,
)
from quarkConstraints.rare_kaon_dilepton import (
    ckm_factors as _ckm_factors,
)
from quarkConstraints.rare_kaon_dilepton import (
    compute_rare_kaon_dilepton_wilsons as _compute_rare_kaon_dilepton_wilsons,
)
from quarkConstraints.rare_kaon_dilepton import (
    default_sm_inputs as _default_sm_inputs,
)
from quarkConstraints.rare_kaon_dilepton import (
    evaluate_klong_mumu_short_distance as _evaluate_klong_mumu_short_distance,
)
from quarkConstraints.rare_kaon_dilepton import (
    evaluate_klong_pi0ee_direct_cp as _evaluate_klong_pi0ee_direct_cp,
)
from quarkConstraints.rare_kaon_dilepton import (
    g_sm_squared as _g_sm_squared,
)
from quarkConstraints.rare_kaon_dilepton import (
    kappa_mu as _kappa_mu,
)
from quarkConstraints.rare_kaon_dilepton import (
    klong_mumu_short_distance_sm as _klong_mumu_short_distance_sm,
)
from quarkConstraints.rare_kaon_dilepton import (
    klong_pi0ee_direct_cp_sm as _klong_pi0ee_direct_cp_sm,
)
from quarkConstraints.rs_semileptonic_wilsons import (
    RS_SEMILEPTONIC_MATCHING_ASSUMPTION_V1,
    RSSemileptonicWilsonBundle,
    RSSemileptonicWilsonCoefficients,
)

RARE_KAON_RS_SEMILEPTONIC_VECTOR_MATCHING_STATUS_V1 = (
    "Phase-3a light-Z RS semileptonic C9/C10/C9p/C10p Wilsons mapped "
    "additively into the rare-kaon Y/y7 effective inputs; no "
    "compute_rare_kaon_dilepton_wilsons proxy and no second 1/M_KK^2 factor."
)

__all__ = [
    "QuarkMassBasisCouplings",
    "RARE_KAON_DILEPTON_MODEL_V1",
    "RARE_KAON_DILEPTON_OPERATOR_CONVENTION",
    "RARE_KAON_DILEPTON_INPUT_BUNDLE_V1",
    "RARE_KAON_DILEPTON_RS_MATCHING_ASSUMPTION_V1",
    "RARE_KAON_DILEPTON_PARAMETRIZATION_CITATION",
    "RARE_KAON_PI0EE_MODEL_V1",
    "RARE_KAON_PI0EE_PARAMETRIZATION_CITATION",
    "RARE_KAON_PI0EE_INTERFERENCE_LIMITATION_V1",
    "RareKaonDileptonSMInputs",
    "RareKaonDileptonCKMFactors",
    "RareKaonDileptonWilsonCoefficients",
    "KLongMuMuShortDistanceResult",
    "KLongPi0EEChPTInputs",
    "KLongPi0EEResult",
    "rare_kaon_dilepton_default_sm_inputs",
    "rare_kaon_dilepton_ckm_factors",
    "rare_kaon_dilepton_kappa_mu",
    "rare_kaon_dilepton_g_sm_squared",
    "rare_kaon_dilepton_wilsons_from_couplings",
    "klong_mumu_short_distance_sm",
    "klong_mumu_short_distance_from_couplings",
    "rare_kaon_y_wilsons_from_rs_semileptonic",
    "rare_kaon_y7_wilsons_from_rs_semileptonic",
    "rare_kaon_rs_semileptonic_vector_diagnostics",
    "klong_mumu_short_distance_from_rs_semileptonic_wilsons",
    "klong_pi0ee_direct_cp_sm",
    "klong_pi0ee_direct_cp_from_couplings",
    "klong_pi0ee_y7_direct_cp_from_rs_semileptonic_wilsons",
    "RARE_KAON_RS_SEMILEPTONIC_VECTOR_MATCHING_STATUS_V1",
    "RS_SEMILEPTONIC_MATCHING_ASSUMPTION_V1",
]


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
    lepton: str,
) -> RSSemileptonicWilsonCoefficients:
    try:
        coeff = source.s_to_d_ll[lepton]
    except (AttributeError, KeyError, TypeError) as exc:
        raise ValueError(
            f"rs_semileptonic_wilsons.s_to_d_ll[{lepton!r}] is not available"
        ) from exc
    if coeff.transition_key != "s_d":
        raise ValueError(
            f"s_to_d_ll[{lepton!r}] transition_key={coeff.transition_key!r}, "
            "expected 's_d'"
        )
    return coeff


def _kaon_c9_c10_to_y_norm(inputs: RareKaonDileptonSMInputs) -> float:
    return float(
        math.sqrt(2.0)
        * inputs.gf_gev_minus2
        * inputs.alpha_em_mz
        / (math.pi * _g_sm_squared(inputs))
    )


def rare_kaon_rs_semileptonic_vector_diagnostics(
    coeff: RSSemileptonicWilsonCoefficients,
    *,
    inputs: RareKaonDileptonSMInputs | None = None,
) -> dict[str, object]:
    p = rare_kaon_dilepton_default_sm_inputs() if inputs is None else inputs
    return {
        "rs_semileptonic_model_label": coeff.model_label,
        "rs_semileptonic_operator_convention": coeff.operator_convention,
        "rs_semileptonic_matching_assumption": coeff.matching_assumption,
        "rs_semileptonic_wilsons_present": True,
        "rs_semileptonic_matching_status": (
            "rs_semileptonic_additive_no_second_1_over_M_KK_squared"
        ),
        "rs_semileptonic_vector_matching_status": (
            RARE_KAON_RS_SEMILEPTONIC_VECTOR_MATCHING_STATUS_V1
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
        "kaon_c9_c10_to_y_normalization": _kaon_c9_c10_to_y_norm(p),
        "wilson_prefactor_reused": False,
        "second_mkk_suppression_applied": False,
    }


def _tag_rs_result(
    result,
    coeff: RSSemileptonicWilsonCoefficients,
    *,
    inputs: RareKaonDileptonSMInputs,
):
    diagnostics = dict(result.diagnostics)
    diagnostics.update(rare_kaon_rs_semileptonic_vector_diagnostics(coeff, inputs=inputs))
    diagnostics["matching_assumption"] = coeff.matching_assumption
    diagnostics["s_to_d_ll_rs_semileptonic_rewired"] = True
    diagnostics["rare_kaon_proxy_reused"] = False
    diagnostics["c9_c10_to_rare_kaon_effective_inputs"] = True
    if "uses_k008_y7v_y7a_rs_proxy" in diagnostics:
        diagnostics["uses_k008_y7v_y7a_rs_proxy"] = False
        diagnostics["uses_k008_y7v_y7a_core_inputs"] = True
    if "uses_vector_y7v_proxy_only_for_a_s" in diagnostics:
        diagnostics["uses_vector_y7v_proxy_only_for_a_s"] = False
        diagnostics["uses_vector_y7v_only_for_a_s"] = True
    if "uses_k006_muon_sd_wilson_proxy" in diagnostics:
        diagnostics["uses_k006_muon_sd_wilson_proxy"] = False
        diagnostics["uses_k006_muon_sd_y_core_inputs"] = True
    if "lambda_y7v_np_proxy" in diagnostics:
        diagnostics["lambda_y7v_np"] = diagnostics["lambda_y7v_np_proxy"]
    if "lambda_y7a_np_proxy" in diagnostics:
        diagnostics["lambda_y7a_np"] = diagnostics["lambda_y7a_np_proxy"]
    if "a_s_np_proxy" in diagnostics:
        diagnostics["a_s_np_shift"] = diagnostics["a_s_np_proxy"]
    return replace(result, diagnostics=diagnostics)


def rare_kaon_y_wilsons_from_rs_semileptonic(
    source: RSSemileptonicWilsonBundle,
    *,
    lepton: str = "mu",
    matching_scale_gev: float | None = None,
    inputs: RareKaonDileptonSMInputs | None = None,
) -> RareKaonDileptonWilsonCoefficients:
    """Translate Phase-3a ``s -> d l l`` C10/C10p to the kaon ``Y_NP`` slot."""

    p = rare_kaon_dilepton_default_sm_inputs() if inputs is None else inputs
    coeff = _rs_semileptonic_coeff(source, lepton=lepton)
    scale = _diagnostic_matching_scale(matching_scale_gev)
    norm = _kaon_c9_c10_to_y_norm(p)
    lambda_ckm = complex(coeff.lambda_ckm)
    y_np_left = -lambda_ckm * norm * complex(coeff.c10_np)
    y_np_right = -lambda_ckm * norm * complex(coeff.c10p_np)
    left_axial_contact = complex(coeff.contact_LR - coeff.contact_LL)
    right_axial_contact = complex(coeff.contact_RR - coeff.contact_RL)
    return RareKaonDileptonWilsonCoefficients(
        model_label=coeff.model_label,
        operator_convention=coeff.operator_convention,
        matching_assumption=coeff.matching_assumption,
        M_KK=scale,
        matching_scale=scale,
        left_sd_coupling=left_axial_contact,
        right_sd_coupling=right_axial_contact,
        left_sd_overlap=0.0j,
        right_sd_overlap=0.0j,
        left_quark_delta=left_axial_contact,
        right_quark_delta=right_axial_contact,
        muon_axial_delta=0.0,
        y_np_left=complex(y_np_left),
        y_np_right=complex(y_np_right),
    )


def rare_kaon_y7_wilsons_from_rs_semileptonic(
    source: RSSemileptonicWilsonBundle,
    *,
    lepton: str = "e",
    matching_scale_gev: float | None = None,
    inputs: RareKaonDileptonSMInputs | None = None,
) -> KLongPi0EEWilsonCoefficients:
    """Translate Phase-3a C9/C10 Wilsons to kaon ``lambda*y7V/y7A`` slots."""

    p = rare_kaon_dilepton_default_sm_inputs() if inputs is None else inputs
    coeff = _rs_semileptonic_coeff(source, lepton=lepton)
    scale = _diagnostic_matching_scale(matching_scale_gev)
    norm = _kaon_c9_c10_to_y_norm(p)
    lambda_ckm = complex(coeff.lambda_ckm)
    lambda_y7v_np = -lambda_ckm * norm * complex(coeff.c9_np + coeff.c9p_np)
    lambda_y7a_np = -lambda_ckm * norm * complex(coeff.c10_np + coeff.c10p_np)
    left_vector_contact = complex(coeff.contact_LL + coeff.contact_LR)
    right_vector_contact = complex(coeff.contact_RL + coeff.contact_RR)
    quark_vector_contact = left_vector_contact + right_vector_contact
    return KLongPi0EEWilsonCoefficients(
        model_label=coeff.model_label,
        operator_convention=coeff.operator_convention,
        matching_assumption=coeff.matching_assumption,
        M_KK=scale,
        matching_scale=scale,
        low_scale_gev=_rare_kaon_dilepton_core.KLONG_PI0EE_Y7_LOW_SCALE_GEV,
        left_sd_coupling=left_vector_contact,
        right_sd_coupling=right_vector_contact,
        left_sd_overlap=0.0j,
        right_sd_overlap=0.0j,
        left_quark_delta=left_vector_contact,
        right_quark_delta=right_vector_contact,
        quark_vector_delta=quark_vector_contact,
        electron_vector_delta=0.0,
        electron_axial_delta=0.0,
        lambda_y7v_np_proxy=complex(lambda_y7v_np),
        lambda_y7a_np_proxy=complex(lambda_y7a_np),
    )


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


def klong_mumu_short_distance_from_rs_semileptonic_wilsons(
    source: RSSemileptonicWilsonBundle,
    *,
    lepton: str = "mu",
    matching_scale_gev: float | None = None,
    inputs: RareKaonDileptonSMInputs | None = None,
) -> KLongMuMuShortDistanceResult:
    """Evaluate ``BR(K_L -> mu+mu-)_SD`` from Phase-3a RS C10/C10p."""

    p = rare_kaon_dilepton_default_sm_inputs() if inputs is None else inputs
    coeff = _rs_semileptonic_coeff(source, lepton=lepton)
    wilsons = rare_kaon_y_wilsons_from_rs_semileptonic(
        source,
        lepton=lepton,
        matching_scale_gev=matching_scale_gev,
        inputs=p,
    )
    result = _evaluate_klong_mumu_short_distance(wilsons, inputs=p)
    return _tag_rs_result(result, coeff, inputs=p)


def klong_pi0ee_direct_cp_sm(
    chpt_inputs: KLongPi0EEChPTInputs,
    *,
    inputs: RareKaonDileptonSMInputs | None = None,
) -> KLongPi0EEResult:
    """Evaluate the SM-limit ``K_L -> pi0 e+ e-`` direct-CP rate."""
    return _klong_pi0ee_direct_cp_sm(chpt_inputs, inputs=inputs)


def klong_pi0ee_direct_cp_from_couplings(
    couplings: QuarkMassBasisCouplings,
    *,
    chpt_inputs: KLongPi0EEChPTInputs,
    m_kk_gev: float | None = None,
    inputs: RareKaonDileptonSMInputs | None = None,
) -> KLongPi0EEResult:
    """Evaluate ``BR(K_L -> pi0 e+ e-)_direct CP`` from mass-basis couplings."""
    return _evaluate_klong_pi0ee_direct_cp(
        couplings,
        chpt_inputs=chpt_inputs,
        m_kk_gev=m_kk_gev,
        inputs=inputs,
    )


from quarkConstraints import (  # noqa: E402, I001 (append-only section)
    rare_kaon_dilepton as _rare_kaon_dilepton_core,
)

RARE_KAON_PI0EE_MODEL_V2 = _rare_kaon_dilepton_core.RARE_KAON_PI0EE_MODEL_V2
RARE_KAON_PI0EE_Y7_OPERATOR_CONVENTION = (
    _rare_kaon_dilepton_core.RARE_KAON_PI0EE_Y7_OPERATOR_CONVENTION
)
RARE_KAON_PI0EE_RS_MATCHING_ASSUMPTION_V2 = (
    _rare_kaon_dilepton_core.RARE_KAON_PI0EE_RS_MATCHING_ASSUMPTION_V2
)
RARE_KAON_PI0EE_SEMILEPTONIC_RUNNING_DIAGNOSTIC_V1 = (
    _rare_kaon_dilepton_core.RARE_KAON_PI0EE_SEMILEPTONIC_RUNNING_DIAGNOSTIC_V1
)
KLongPi0EEWilsonCoefficients = (
    _rare_kaon_dilepton_core.KLongPi0EEWilsonCoefficients
)


def klong_pi0ee_y7_direct_cp_sm(
    chpt_inputs: KLongPi0EEChPTInputs,
    *,
    inputs: RareKaonDileptonSMInputs | None = None,
) -> KLongPi0EEResult:
    """Evaluate the SM-limit ``K_L -> pi0 e+ e-`` y7V/y7A rate."""
    return _rare_kaon_dilepton_core.klong_pi0ee_y7_direct_cp_sm(
        chpt_inputs,
        inputs=inputs,
    )


def klong_pi0ee_y7_direct_cp_from_couplings(
    couplings: QuarkMassBasisCouplings,
    *,
    chpt_inputs: KLongPi0EEChPTInputs,
    m_kk_gev: float | None = None,
    inputs: RareKaonDileptonSMInputs | None = None,
) -> KLongPi0EEResult:
    """Evaluate K008 with explicit ``y7V/y7A`` Wilson structure."""
    return _rare_kaon_dilepton_core.evaluate_klong_pi0ee_y7_direct_cp(
        couplings,
        chpt_inputs=chpt_inputs,
        m_kk_gev=m_kk_gev,
        inputs=inputs,
    )


def klong_pi0ee_y7_direct_cp_from_rs_semileptonic_wilsons(
    source: RSSemileptonicWilsonBundle,
    *,
    chpt_inputs: KLongPi0EEChPTInputs,
    lepton: str = "e",
    matching_scale_gev: float | None = None,
    inputs: RareKaonDileptonSMInputs | None = None,
) -> KLongPi0EEResult:
    """Evaluate K008 from Phase-3a RS C9/C10 mapped to y7V/y7A."""

    p = rare_kaon_dilepton_default_sm_inputs() if inputs is None else inputs
    coeff = _rs_semileptonic_coeff(source, lepton=lepton)
    wilsons = rare_kaon_y7_wilsons_from_rs_semileptonic(
        source,
        lepton=lepton,
        matching_scale_gev=matching_scale_gev,
        inputs=p,
    )
    result = _rare_kaon_dilepton_core.evaluate_klong_pi0ee_y7_direct_cp(
        wilsons,
        chpt_inputs=chpt_inputs,
        inputs=p,
    )
    return _tag_rs_result(result, coeff, inputs=p)


__all__ += [
    "RARE_KAON_PI0EE_MODEL_V2",
    "RARE_KAON_PI0EE_Y7_OPERATOR_CONVENTION",
    "RARE_KAON_PI0EE_RS_MATCHING_ASSUMPTION_V2",
    "RARE_KAON_PI0EE_SEMILEPTONIC_RUNNING_DIAGNOSTIC_V1",
    "KLongPi0EEWilsonCoefficients",
    "klong_pi0ee_y7_direct_cp_sm",
    "klong_pi0ee_y7_direct_cp_from_couplings",
    "klong_pi0ee_y7_direct_cp_from_rs_semileptonic_wilsons",
]
