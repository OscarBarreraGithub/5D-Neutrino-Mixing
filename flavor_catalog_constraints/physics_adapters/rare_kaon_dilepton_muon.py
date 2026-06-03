"""Muon-mode wrapper for ``K_L -> pi0 mu+ mu-`` rare-kaon dileptons.

This adapter is append-only relative to the K008 electron implementation.  It
reuses the existing rare-kaon-dilepton CKM inputs and y7V/y7A core shape, with
Phase-3a ``rs_semileptonic_wilsons.s_to_d_ll`` filling the short-distance RS
C9/C10/C9p/C10p slots for the muon mode.

NEEDS-HUMAN-PHYSICS: K009's long-distance interference/CPC treatment and
neutral-current pieces beyond Phase-3a light-Z matching remain partial.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Mapping

from quarkConstraints import rare_kaon_dilepton as _core

from .rare_kaon_dilepton import (
    KLongPi0EEChPTInputs,
    KLongPi0EEWilsonCoefficients,
    QuarkMassBasisCouplings,
    RareKaonDileptonSMInputs,
    RARE_KAON_PI0EE_PARAMETRIZATION_CITATION,
    RARE_KAON_PI0EE_SEMILEPTONIC_RUNNING_DIAGNOSTIC_V1,
    _rs_semileptonic_coeff,
    _tag_rs_result,
    rare_kaon_y7_wilsons_from_rs_semileptonic,
    rare_kaon_dilepton_ckm_factors,
    rare_kaon_dilepton_default_sm_inputs,
    rare_kaon_dilepton_g_sm_squared,
)
from quarkConstraints.rs_semileptonic_wilsons import RSSemileptonicWilsonBundle

RARE_KAON_PI0MUMU_MODEL_V1 = "rare_kaon_pi0mumu_y7v_y7a_isu_rs_proxy_v1"
RARE_KAON_PI0MUMU_Y7_OPERATOR_CONVENTION = (
    "Q7V=(sbar gamma_mu(1-gamma5)d)(mubar gamma^mu mu), "
    "Q7A=(sbar gamma_mu(1-gamma5)d)(mubar gamma^mu gamma5 mu); "
    "direct CP uses separate lambda_t*y7V and lambda_t*y7A amplitudes, "
    "interference uses the vector Q7V amplitude only"
)
RARE_KAON_PI0MUMU_RS_MATCHING_ASSUMPTION_V1 = (
    "NEEDS-HUMAN-PHYSICS: K009 maps Phase-3a s -> d mu mu C9/C10 "
    "Wilson shifts into the K008 y7V/y7A core shape. Long-distance "
    "interference/CPC treatment and non-light-Z neutral-current effects "
    "remain outside Phase 3a."
)
RARE_KAON_PI0MUMU_PARAMETRIZATION_CITATION = (
    RARE_KAON_PI0EE_PARAMETRIZATION_CITATION
)
KLONG_PI0MUMU_SM_Y7V_BAR = float(_core.KLONG_PI0EE_SM_Y7V_BAR)
KLONG_PI0MUMU_SM_Y7A_BAR = float(_core.KLONG_PI0EE_SM_Y7A_BAR)
KLONG_PI0MUMU_Y7_LOW_SCALE_GEV = float(_core.KLONG_PI0EE_Y7_LOW_SCALE_GEV)


@dataclass(frozen=True)
class KLongPi0MuMuResult:
    """Branching-ratio prediction for ``K_L -> pi0 mu+ mu-``.

    ``direct_cp_branching_fraction`` is the short-distance quantity used by the
    K009 HARD verdict.  Indirect/interference/CPC totals are diagnostic because
    the interference sign and full long-distance treatment need human physics
    review.
    """

    model_label: str
    input_bundle: str
    direct_cp_branching_fraction: float
    sm_direct_cp_branching_fraction: float
    direct_cp_np_shift_branching_fraction: float
    constructive_total_branching_fraction: float
    destructive_total_branching_fraction: float
    sm_constructive_total_branching_fraction: float
    sm_destructive_total_branching_fraction: float
    indirect_cp_branching_fraction: float
    interference_abs_branching_fraction: float
    cpc_branching_fraction: float
    direct_cp_amplitude_ratio: float
    sm_direct_cp_amplitude_ratio: float
    direct_cp_vector_amplitude_ratio: float
    direct_cp_axial_amplitude_ratio: float
    sm_direct_cp_vector_amplitude_ratio: float
    sm_direct_cp_axial_amplitude_ratio: float
    lambda_y7v_eff: complex
    lambda_y7a_eff: complex
    lambda_y7v_np_proxy: complex
    lambda_y7a_np_proxy: complex
    lambda_wolfenstein: float
    lambda_t: complex
    chpt_inputs: KLongPi0EEChPTInputs
    wilsons: KLongPi0EEWilsonCoefficients | None = None
    diagnostics: Mapping[str, float | complex | str | bool] = field(
        default_factory=dict
    )


def _muon_y7_phase_space_coefficients(
    chpt_inputs: KLongPi0EEChPTInputs,
) -> tuple[float, float]:
    """Return y7-normalized direct and interference coefficients.

    K009.yaml stores the ISU coefficients in the traditional
    ``A = Im(lambda_t) / 1e-4`` form.  The y7V/y7A split uses
    ``A_V = Im(lambda_t*y7V) / 1e-4`` and
    ``A_A = Im(lambda_t*y7A) / 1e-4``.  These rescalings make the SM limit
    reduce to the YAML coefficients while preserving the K008 vector/axial
    structure for NP shifts.
    """

    y7_norm = KLONG_PI0MUMU_SM_Y7V_BAR**2 + KLONG_PI0MUMU_SM_Y7A_BAR**2
    if y7_norm <= 0.0:
        raise ValueError("muon y7 direct normalization must be positive")
    y7v_abs = abs(KLONG_PI0MUMU_SM_Y7V_BAR)
    if y7v_abs <= 0.0:
        raise ValueError("muon y7 vector normalization must be positive")
    return float(chpt_inputs.c_dir / y7_norm), float(chpt_inputs.c_int / y7v_abs)


def _klong_pi0mumu_components_from_y7(
    lambda_y7v_np_proxy: complex,
    lambda_y7a_np_proxy: complex,
    *,
    inputs: RareKaonDileptonSMInputs,
    chpt_inputs: KLongPi0EEChPTInputs,
) -> tuple[
    float,
    float,
    float,
    float,
    float,
    object,
    complex,
    complex,
    float,
    float,
    float,
    float,
    float,
]:
    factors = rare_kaon_dilepton_ckm_factors(inputs)
    lambda_y7v_eff = (
        factors.lambda_t * KLONG_PI0MUMU_SM_Y7V_BAR
        + complex(lambda_y7v_np_proxy)
    )
    lambda_y7a_eff = (
        factors.lambda_t * KLONG_PI0MUMU_SM_Y7A_BAR
        + complex(lambda_y7a_np_proxy)
    )
    vector_ratio = lambda_y7v_eff.imag / 1.0e-4
    axial_ratio = lambda_y7a_eff.imag / 1.0e-4
    c_dir_y7va, c_int_y7v = _muon_y7_phase_space_coefficients(chpt_inputs)
    direct = (
        chpt_inputs.normalization
        * c_dir_y7va
        * (vector_ratio**2 + axial_ratio**2)
    )
    indirect = (
        chpt_inputs.normalization
        * chpt_inputs.c_mix
        * chpt_inputs.a_s_abs**2
    )
    interference_abs = (
        chpt_inputs.normalization
        * c_int_y7v
        * chpt_inputs.a_s_abs
        * abs(vector_ratio)
    )
    cpc = chpt_inputs.normalization * chpt_inputs.c_cpc
    constructive = indirect + direct + cpc + interference_abs
    destructive = max(indirect + direct + cpc - interference_abs, 0.0)
    amplitude_quadrature = math.sqrt(vector_ratio**2 + axial_ratio**2)
    return (
        float(direct),
        float(constructive),
        float(destructive),
        float(indirect),
        float(interference_abs),
        factors,
        complex(lambda_y7v_eff),
        complex(lambda_y7a_eff),
        float(vector_ratio),
        float(axial_ratio),
        float(amplitude_quadrature),
        float(c_dir_y7va),
        float(c_int_y7v),
    )


def klong_pi0mumu_y7_direct_cp_sm(
    chpt_inputs: KLongPi0EEChPTInputs,
    *,
    inputs: RareKaonDileptonSMInputs | None = None,
) -> KLongPi0MuMuResult:
    """Evaluate the SM-limit K009 rate with explicit ``y7V/y7A`` structure."""

    return evaluate_klong_pi0mumu_y7_direct_cp(
        None,
        chpt_inputs=chpt_inputs,
        inputs=inputs,
    )


def klong_pi0mumu_y7_direct_cp_from_couplings(
    couplings: QuarkMassBasisCouplings,
    *,
    chpt_inputs: KLongPi0EEChPTInputs,
    m_kk_gev: float | None = None,
    inputs: RareKaonDileptonSMInputs | None = None,
) -> KLongPi0MuMuResult:
    """Evaluate K009 using the reused K008 ``y7V/y7A`` proxy."""

    return evaluate_klong_pi0mumu_y7_direct_cp(
        couplings,
        chpt_inputs=chpt_inputs,
        m_kk_gev=m_kk_gev,
        inputs=inputs,
    )


def klong_pi0mumu_y7_direct_cp_from_rs_semileptonic_wilsons(
    source: RSSemileptonicWilsonBundle,
    *,
    chpt_inputs: KLongPi0EEChPTInputs,
    lepton: str = "mu",
    matching_scale_gev: float | None = None,
    inputs: RareKaonDileptonSMInputs | None = None,
) -> KLongPi0MuMuResult:
    """Evaluate K009 from Phase-3a RS C9/C10 mapped to y7V/y7A."""

    p = rare_kaon_dilepton_default_sm_inputs() if inputs is None else inputs
    coeff = _rs_semileptonic_coeff(source, lepton=lepton)
    wilsons = rare_kaon_y7_wilsons_from_rs_semileptonic(
        source,
        lepton=lepton,
        matching_scale_gev=matching_scale_gev,
        inputs=p,
    )
    result = evaluate_klong_pi0mumu_y7_direct_cp(
        wilsons,
        chpt_inputs=chpt_inputs,
        inputs=p,
    )
    return _tag_rs_result(result, coeff, inputs=p)


def evaluate_klong_pi0mumu_y7_direct_cp(
    source: QuarkMassBasisCouplings | KLongPi0EEWilsonCoefficients | None = None,
    *,
    chpt_inputs: KLongPi0EEChPTInputs,
    m_kk_gev: float | None = None,
    inputs: RareKaonDileptonSMInputs | None = None,
) -> KLongPi0MuMuResult:
    """Evaluate K009 direct CP plus diagnostic indirect/CPC envelopes."""

    p = rare_kaon_dilepton_default_sm_inputs() if inputs is None else inputs
    wilsons: KLongPi0EEWilsonCoefficients | None
    if source is None:
        wilsons = None
        lambda_y7v_np_proxy = 0.0j
        lambda_y7a_np_proxy = 0.0j
    elif isinstance(source, KLongPi0EEWilsonCoefficients):
        wilsons = source
        lambda_y7v_np_proxy = source.lambda_y7v_np_proxy
        lambda_y7a_np_proxy = source.lambda_y7a_np_proxy
    else:
        wilsons = _core.compute_klong_pi0ee_y7_wilsons(
            source,
            m_kk_gev=m_kk_gev,
            inputs=p,
        )
        lambda_y7v_np_proxy = wilsons.lambda_y7v_np_proxy
        lambda_y7a_np_proxy = wilsons.lambda_y7a_np_proxy

    (
        direct,
        constructive,
        destructive,
        indirect,
        interference_abs,
        factors,
        lambda_y7v_eff,
        lambda_y7a_eff,
        vector_ratio,
        axial_ratio,
        amplitude_quadrature,
        c_dir_y7va,
        c_int_y7v,
    ) = _klong_pi0mumu_components_from_y7(
        lambda_y7v_np_proxy,
        lambda_y7a_np_proxy,
        inputs=p,
        chpt_inputs=chpt_inputs,
    )
    (
        sm_direct,
        sm_constructive,
        sm_destructive,
        _sm_indirect,
        _sm_interference_abs,
        _sm_factors,
        sm_lambda_y7v_eff,
        sm_lambda_y7a_eff,
        sm_vector_ratio,
        sm_axial_ratio,
        sm_amplitude_quadrature,
        _sm_c_dir_y7va,
        _sm_c_int_y7v,
    ) = _klong_pi0mumu_components_from_y7(
        0.0j,
        0.0j,
        inputs=p,
        chpt_inputs=chpt_inputs,
    )
    cpc = chpt_inputs.normalization * chpt_inputs.c_cpc
    diagnostics: dict[str, float | complex | str | bool] = {
        "g_sm_squared_gev_minus2": rare_kaon_dilepton_g_sm_squared(p),
        "matching_assumption": RARE_KAON_PI0MUMU_RS_MATCHING_ASSUMPTION_V1,
        "parametrization_citation": RARE_KAON_PI0MUMU_PARAMETRIZATION_CITATION,
        "interference_cpc_limitation": (
            "NEEDS-HUMAN-PHYSICS: K_L -> pi0 mu+mu- direct/indirect "
            "interference depends on the ChPT sign of a_S and the CPC "
            "two-photon treatment; the HARD-veto quantity is the direct-CP "
            "short-distance term, while indirect/interference/CPC totals are "
            "diagnostic envelopes."
        ),
        "short_distance_direct_cp_only": True,
        "interference_sign_fixed": False,
        "cpc_from_two_photon_fixed": False,
        "uses_y7v_y7a_wilson_structure": True,
        "uses_k008_y7v_y7a_rs_proxy": True,
        "uses_y_function_proxy": False,
        "interference_uses_vector_y7v_only": True,
        "direct_uses_vector_and_axial_y7": True,
        "muon_phase_space_from_catalog_coefficients": True,
        "quark_vector_current_uses_left_plus_right": True,
        "sm_y7v_bar": float(KLONG_PI0MUMU_SM_Y7V_BAR),
        "sm_y7a_bar": float(KLONG_PI0MUMU_SM_Y7A_BAR),
        "c_int_y7v_coefficient": float(c_int_y7v),
        "c_dir_y7va_coefficient": float(c_dir_y7va),
        "catalog_c_int_sm_coefficient": float(chpt_inputs.c_int),
        "catalog_c_dir_sm_coefficient": float(chpt_inputs.c_dir),
        "lambda_y7v_eff": complex(lambda_y7v_eff),
        "lambda_y7a_eff": complex(lambda_y7a_eff),
        "lambda_y7v_np_proxy": complex(lambda_y7v_np_proxy),
        "lambda_y7a_np_proxy": complex(lambda_y7a_np_proxy),
        "direct_cp_vector_amplitude_ratio": float(vector_ratio),
        "direct_cp_axial_amplitude_ratio": float(axial_ratio),
        "sm_lambda_y7v_eff": complex(sm_lambda_y7v_eff),
        "sm_lambda_y7a_eff": complex(sm_lambda_y7a_eff),
        "sm_direct_cp_vector_amplitude_ratio": float(sm_vector_ratio),
        "sm_direct_cp_axial_amplitude_ratio": float(sm_axial_ratio),
        "sm_interference_abs_branching_fraction": float(_sm_interference_abs),
        "semileptonic_qcd_running_applied": False,
        "semileptonic_qcd_running_multiplicative_factor": 1.0,
        "semileptonic_qcd_running_effect_fraction": 0.0,
        "semileptonic_qcd_running_diagnostic": (
            RARE_KAON_PI0EE_SEMILEPTONIC_RUNNING_DIAGNOSTIC_V1
        ),
    }
    if wilsons is not None:
        diagnostics.update(
            {
                "m_kk_gev": float(wilsons.M_KK),
                "matching_scale_gev": float(wilsons.matching_scale),
                "semileptonic_low_scale_gev": float(wilsons.low_scale_gev),
                "left_sd_coupling": complex(wilsons.left_sd_coupling),
                "right_sd_coupling": complex(wilsons.right_sd_coupling),
                "left_sd_overlap": complex(wilsons.left_sd_overlap),
                "right_sd_overlap": complex(wilsons.right_sd_overlap),
                "left_quark_delta": complex(wilsons.left_quark_delta),
                "right_quark_delta": complex(wilsons.right_quark_delta),
                "quark_vector_delta": complex(wilsons.quark_vector_delta),
                "muon_vector_delta_proxy": float(wilsons.electron_vector_delta),
                "muon_axial_delta_proxy": float(wilsons.electron_axial_delta),
            }
        )

    return KLongPi0MuMuResult(
        model_label=RARE_KAON_PI0MUMU_MODEL_V1,
        input_bundle=p.input_bundle,
        direct_cp_branching_fraction=direct,
        sm_direct_cp_branching_fraction=sm_direct,
        direct_cp_np_shift_branching_fraction=float(direct - sm_direct),
        constructive_total_branching_fraction=constructive,
        destructive_total_branching_fraction=destructive,
        sm_constructive_total_branching_fraction=sm_constructive,
        sm_destructive_total_branching_fraction=sm_destructive,
        indirect_cp_branching_fraction=indirect,
        interference_abs_branching_fraction=interference_abs,
        cpc_branching_fraction=float(cpc),
        direct_cp_amplitude_ratio=amplitude_quadrature,
        sm_direct_cp_amplitude_ratio=sm_amplitude_quadrature,
        direct_cp_vector_amplitude_ratio=vector_ratio,
        direct_cp_axial_amplitude_ratio=axial_ratio,
        sm_direct_cp_vector_amplitude_ratio=sm_vector_ratio,
        sm_direct_cp_axial_amplitude_ratio=sm_axial_ratio,
        lambda_y7v_eff=complex(lambda_y7v_eff),
        lambda_y7a_eff=complex(lambda_y7a_eff),
        lambda_y7v_np_proxy=complex(lambda_y7v_np_proxy),
        lambda_y7a_np_proxy=complex(lambda_y7a_np_proxy),
        lambda_wolfenstein=float(factors.lambda_wolfenstein),
        lambda_t=complex(factors.lambda_t),
        chpt_inputs=chpt_inputs,
        wilsons=wilsons,
        diagnostics=diagnostics,
    )


__all__ = [
    "RARE_KAON_PI0MUMU_MODEL_V1",
    "RARE_KAON_PI0MUMU_Y7_OPERATOR_CONVENTION",
    "RARE_KAON_PI0MUMU_RS_MATCHING_ASSUMPTION_V1",
    "RARE_KAON_PI0MUMU_PARAMETRIZATION_CITATION",
    "KLONG_PI0MUMU_SM_Y7V_BAR",
    "KLONG_PI0MUMU_SM_Y7A_BAR",
    "KLONG_PI0MUMU_Y7_LOW_SCALE_GEV",
    "KLongPi0MuMuResult",
    "klong_pi0mumu_y7_direct_cp_sm",
    "klong_pi0mumu_y7_direct_cp_from_couplings",
    "klong_pi0mumu_y7_direct_cp_from_rs_semileptonic_wilsons",
    "evaluate_klong_pi0mumu_y7_direct_cp",
]
