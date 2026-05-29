"""Adapter over :mod:`quarkConstraints.rare_b_nunu`.

This is the catalog boundary for the Delta-B=1 ``b -> s nu nubar`` machinery.
Constraint modules import this adapter only; the underlying physics
implementation remains isolated in ``quarkConstraints``.  The shared core
exposes both the ``B+ -> K+ nu nubar`` branching-ratio evaluator and reusable
short-distance ``C_L, C_R, epsilon, eta`` parameters for future ``B -> K* nu
nubar`` constraints.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Mapping

from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.rare_b_nunu import (
    RARE_B_NUNU_BPLUS_KPLUS_SM_CITATION,
    RARE_B_NUNU_INPUT_BUNDLE_V1,
    RARE_B_NUNU_KSTAR_ETA_COEFFICIENT,
    RARE_B_NUNU_MODEL_V1,
    RARE_B_NUNU_OPERATOR_CONVENTION,
    RARE_B_NUNU_RS_MATCHING_ASSUMPTION_V1,
    RareBNuNuBranchingResult,
    RareBNuNuCKMFactors,
    RareBNuNuSMInputs,
    RareBNuNuShortDistanceResult,
    RareBNuNuWilsonCoefficients,
    ckm_factors as _ckm_factors,
    compute_rare_b_nunu_wilsons as _compute_rare_b_nunu_wilsons,
    default_sm_inputs as _default_sm_inputs,
    evaluate_bplus_to_kplus_nunu as _evaluate_bplus_to_kplus_nunu,
    g_sm_squared as _g_sm_squared,
    inami_lim_x0 as _inami_lim_x0,
    short_distance_response as _short_distance_response,
    sm_branching_fraction as _sm_branching_fraction,
    sm_inputs_with_bplus_kplus_normalization as _sm_inputs_with_bplus_kplus_normalization,
    x_t_top_function as _x_t_top_function,
)

__all__ = [
    "QuarkMassBasisCouplings",
    "RARE_B_NUNU_MODEL_V1",
    "RARE_B_NUNU_OPERATOR_CONVENTION",
    "RARE_B_NUNU_INPUT_BUNDLE_V1",
    "RARE_B_NUNU_RS_MATCHING_ASSUMPTION_V1",
    "RARE_B_NUNU_BPLUS_KPLUS_SM_CITATION",
    "RARE_B_NUNU_KSTAR_ETA_COEFFICIENT",
    "RareBNuNuSMInputs",
    "RareBNuNuCKMFactors",
    "RareBNuNuWilsonCoefficients",
    "RareBNuNuShortDistanceResult",
    "RareBNuNuBranchingResult",
    "RareBKstarNuNuBranchingResult",
    "rare_b_nunu_default_sm_inputs",
    "rare_b_nunu_sm_inputs_with_bplus_kplus_normalization",
    "rare_b_nunu_inami_lim_x0",
    "rare_b_nunu_x_t_top_function",
    "rare_b_nunu_ckm_factors",
    "rare_b_nunu_g_sm_squared",
    "rare_b_nunu_short_distance_response",
    "rare_b_nunu_sm_branching_fraction",
    "rare_b_nunu_wilsons_from_couplings",
    "bplus_kplus_nunu_from_couplings",
    "b_to_kstar_nunu_branching_fraction",
    "b_to_kstar_nunu_from_couplings",
]


@dataclass(frozen=True)
class RareBKstarNuNuBranchingResult:
    """Branching-ratio prediction for ``B -> K* nu nubar``.

    The shared core already owns the short-distance ``epsilon, eta, r_kstar``
    response.  This adapter-level result only supplies the exclusive K* SM
    normalization from the catalog sidecar and applies
    ``BR = BR_SM * r_kstar``.
    """

    model_label: str
    input_bundle: str
    branching_fraction: float
    sm_branching_fraction: float
    np_shift_branching_fraction: float
    ratio_to_sm: float
    x_t: float
    lambda_wolfenstein: float
    lambda_t_bs: complex
    c_l_sm: complex
    c_l_total: complex
    c_r_total: complex
    x_eff_left: complex
    x_eff_right: complex
    epsilon: float
    eta: float
    r_k: float
    r_kstar: float
    wilsons: RareBNuNuWilsonCoefficients | None = None
    diagnostics: Mapping[str, float | complex | str] = field(default_factory=dict)


def rare_b_nunu_default_sm_inputs() -> RareBNuNuSMInputs:
    """Return the default rare-B ``b -> s nu nubar`` input bundle."""
    return _default_sm_inputs()


def rare_b_nunu_sm_inputs_with_bplus_kplus_normalization(
    branching_fraction: float,
    *,
    inputs: RareBNuNuSMInputs | None = None,
) -> RareBNuNuSMInputs:
    """Return inputs with the ``B+ -> K+`` SM normalization replaced."""
    return _sm_inputs_with_bplus_kplus_normalization(
        branching_fraction,
        inputs=inputs,
    )


def rare_b_nunu_inami_lim_x0(x_t: float) -> float:
    """Return the leading Inami-Lim ``X_0(x_t)`` function."""
    return _inami_lim_x0(x_t)


def rare_b_nunu_x_t_top_function(inputs: RareBNuNuSMInputs | None = None) -> float:
    """Return the short-distance top function ``X_t``."""
    return _x_t_top_function(inputs)


def rare_b_nunu_ckm_factors(
    inputs: RareBNuNuSMInputs | None = None,
) -> RareBNuNuCKMFactors:
    """Return the CKM factors used by the rare-B core."""
    return _ckm_factors(inputs)


def rare_b_nunu_g_sm_squared(inputs: RareBNuNuSMInputs | None = None) -> float:
    """Return the Buras ``g_SM^2`` normalization."""
    return _g_sm_squared(inputs)


def rare_b_nunu_short_distance_response(
    x_np_left: complex = 0.0j,
    x_np_right: complex = 0.0j,
    *,
    inputs: RareBNuNuSMInputs | None = None,
) -> RareBNuNuShortDistanceResult:
    """Return reusable ``C_L, C_R, epsilon, eta`` response parameters."""
    return _short_distance_response(x_np_left, x_np_right, inputs=inputs)


def rare_b_nunu_sm_branching_fraction(
    inputs: RareBNuNuSMInputs | None = None,
) -> RareBNuNuBranchingResult:
    """Evaluate the SM-limit ``B+ -> K+ nu nubar`` branching fraction."""
    return _sm_branching_fraction(inputs)


def rare_b_nunu_wilsons_from_couplings(
    couplings: QuarkMassBasisCouplings,
    *,
    m_kk_gev: float | None = None,
    inputs: RareBNuNuSMInputs | None = None,
) -> RareBNuNuWilsonCoefficients:
    """Return the v1 ``b -> s nu nubar`` Wilson proxy for mass-basis couplings."""
    return _compute_rare_b_nunu_wilsons(
        couplings,
        m_kk_gev=m_kk_gev,
        inputs=inputs,
    )


def bplus_kplus_nunu_from_couplings(
    couplings: QuarkMassBasisCouplings,
    *,
    m_kk_gev: float | None = None,
    inputs: RareBNuNuSMInputs | None = None,
) -> RareBNuNuBranchingResult:
    """Evaluate ``BR(B+ -> K+ nu nubar)`` from mass-basis couplings."""
    return _evaluate_bplus_to_kplus_nunu(
        couplings,
        m_kk_gev=m_kk_gev,
        inputs=inputs,
    )


def _positive_finite(value: float, name: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0.0:
        raise ValueError(f"{name} must be positive and finite")
    return number


def b_to_kstar_nunu_branching_fraction(
    source: QuarkMassBasisCouplings | RareBNuNuWilsonCoefficients | None = None,
    *,
    sm_branching_fraction: float,
    m_kk_gev: float | None = None,
    inputs: RareBNuNuSMInputs | None = None,
) -> RareBKstarNuNuBranchingResult:
    """Evaluate ``BR(B -> K* nu nubar)`` using the shared K* response.

    ``sm_branching_fraction`` is supplied by the constraint sidecar.  The
    physics core supplies the Wilson proxy and the vector-mode coefficient
    ``r_kstar = (1 + 1.31 eta) epsilon^2``.
    """
    p = _default_sm_inputs() if inputs is None else inputs
    sm = _positive_finite(sm_branching_fraction, "sm_branching_fraction")
    wilsons: RareBNuNuWilsonCoefficients | None
    if source is None:
        wilsons = None
        x_np_left = 0.0j
        x_np_right = 0.0j
    elif isinstance(source, RareBNuNuWilsonCoefficients):
        wilsons = source
        x_np_left = source.x_np_left
        x_np_right = source.x_np_right
    else:
        wilsons = _compute_rare_b_nunu_wilsons(
            source,
            m_kk_gev=m_kk_gev,
            inputs=p,
        )
        x_np_left = wilsons.x_np_left
        x_np_right = wilsons.x_np_right

    response = _short_distance_response(x_np_left, x_np_right, inputs=p)
    factors = _ckm_factors(p)
    br = float(sm * response.r_kstar)
    diagnostics: dict[str, float | complex | str] = {
        "g_sm_squared_gev_minus2": _g_sm_squared(p),
        "matching_assumption": RARE_B_NUNU_RS_MATCHING_ASSUMPTION_V1,
        "operator_convention": RARE_B_NUNU_OPERATOR_CONVENTION,
        "constants_citation": p.constants_citation,
        "bkstar_total_sm_branching_fraction": sm,
        "bkstar_response_formula": "BR(B -> K* nu nubar) = BR_SM * r_kstar",
        "kstar_eta_coefficient": float(RARE_B_NUNU_KSTAR_ETA_COEFFICIENT),
        "m_t_msbar_gev": float(p.m_t_msbar_gev),
        "m_w_gev": float(p.m_w_gev),
        "eta_x": float(p.eta_x),
        "sin2_theta_w": float(p.sin2_theta_w),
        "x_t": float(response.x_t),
        "c_l_sm": complex(response.c_l_sm),
        "c_l_total": complex(response.c_l_total),
        "c_r_total": complex(response.c_r_total),
        "x_eff_left": complex(response.x_eff_left),
        "x_eff_right": complex(response.x_eff_right),
        "epsilon": float(response.epsilon),
        "eta": float(response.eta),
        "r_k": float(response.r_k),
        "r_kstar": float(response.r_kstar),
    }
    if wilsons is not None:
        diagnostics.update(
            {
                "m_kk_gev": float(wilsons.M_KK),
                "matching_scale_gev": float(wilsons.matching_scale),
                "left_sb_coupling": complex(wilsons.left_sb_coupling),
                "right_sb_coupling": complex(wilsons.right_sb_coupling),
                "left_sb_overlap": complex(wilsons.left_sb_overlap),
                "right_sb_overlap": complex(wilsons.right_sb_overlap),
                "left_quark_delta": complex(wilsons.left_quark_delta),
                "right_quark_delta": complex(wilsons.right_quark_delta),
                "neutrino_delta": float(wilsons.neutrino_delta),
                "x_np_left": complex(wilsons.x_np_left),
                "x_np_right": complex(wilsons.x_np_right),
                "x_np_total": complex(wilsons.x_np_total),
            }
        )

    return RareBKstarNuNuBranchingResult(
        model_label=RARE_B_NUNU_MODEL_V1,
        input_bundle=p.input_bundle,
        branching_fraction=br,
        sm_branching_fraction=sm,
        np_shift_branching_fraction=float(br - sm),
        ratio_to_sm=float(br / sm),
        x_t=float(response.x_t),
        lambda_wolfenstein=float(factors.lambda_wolfenstein),
        lambda_t_bs=complex(response.lambda_t_bs),
        c_l_sm=complex(response.c_l_sm),
        c_l_total=complex(response.c_l_total),
        c_r_total=complex(response.c_r_total),
        x_eff_left=complex(response.x_eff_left),
        x_eff_right=complex(response.x_eff_right),
        epsilon=float(response.epsilon),
        eta=float(response.eta),
        r_k=float(response.r_k),
        r_kstar=float(response.r_kstar),
        wilsons=wilsons,
        diagnostics=diagnostics,
    )


def b_to_kstar_nunu_from_couplings(
    couplings: QuarkMassBasisCouplings,
    *,
    sm_branching_fraction: float,
    m_kk_gev: float | None = None,
    inputs: RareBNuNuSMInputs | None = None,
) -> RareBKstarNuNuBranchingResult:
    """Evaluate ``BR(B -> K* nu nubar)`` from mass-basis couplings."""
    return b_to_kstar_nunu_branching_fraction(
        couplings,
        sm_branching_fraction=sm_branching_fraction,
        m_kk_gev=m_kk_gev,
        inputs=inputs,
    )
