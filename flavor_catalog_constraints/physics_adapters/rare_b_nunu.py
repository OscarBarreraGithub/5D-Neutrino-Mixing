"""Adapter over :mod:`quarkConstraints.rare_b_nunu`.

This is the catalog boundary for the Delta-B=1 ``b -> s nu nubar`` machinery.
Constraint modules import this adapter only; the underlying physics
implementation remains isolated in ``quarkConstraints``.  Production Phase-4d
constraints consume the rigorous Phase-4a ``rs_semileptonic_wilsons``
``b_to_s_nunu`` block and map ``X_NP=C/g_SM^2`` directly into the shared
short-distance core, with no charged-lepton Wilson prefactor and no second
``1/M_KK^2`` factor.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field, replace
from typing import Mapping

import numpy as np

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
    RareBNuNuShortDistanceResult,
    RareBNuNuSMInputs,
    RareBNuNuWilsonCoefficients,
)
from quarkConstraints.rare_b_nunu import (
    ckm_factors as _ckm_factors,
)
from quarkConstraints.rare_b_nunu import (
    compute_rare_b_nunu_wilsons as _compute_rare_b_nunu_wilsons,
)
from quarkConstraints.rare_b_nunu import (
    default_sm_inputs as _default_sm_inputs,
)
from quarkConstraints.rare_b_nunu import (
    evaluate_bplus_to_kplus_nunu as _evaluate_bplus_to_kplus_nunu,
)
from quarkConstraints.rare_b_nunu import (
    g_sm_squared as _g_sm_squared,
)
from quarkConstraints.rare_b_nunu import (
    inami_lim_x0 as _inami_lim_x0,
)
from quarkConstraints.rare_b_nunu import (
    short_distance_response as _short_distance_response,
)
from quarkConstraints.rare_b_nunu import (
    sm_branching_fraction as _sm_branching_fraction,
)
from quarkConstraints.rare_b_nunu import (
    sm_inputs_with_bplus_kplus_normalization as _sm_inputs_with_bplus_kplus_normalization,
)
from quarkConstraints.rare_b_nunu import (
    x_t_top_function as _x_t_top_function,
)
from quarkConstraints.rs_semileptonic_wilsons import (
    RS_SEMILEPTONIC_MATCHING_ASSUMPTION_V1,
    RSNuNuWilsonCoefficients,
    RSSemileptonicWilsonBundle,
)

RARE_B_NUNU_RS_SEMILEPTONIC_MATCHING_STATUS_V1 = (
    "Phase-4d light-Z RS active-neutrino Wilsons consumed as X_NP=C/g_SM^2; "
    "no rare_b_nunu proxy matcher, no charged-lepton _wilson_prefactor, and "
    "no second 1/M_KK^2 factor."
)
_ACTIVE_NU_UNIVERSALITY_ATOL = 1.0e-10

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
    "rare_b_nunu_wilsons_from_rs_semileptonic",
    "rare_b_nunu_rs_semileptonic_diagnostics",
    "bplus_kplus_nunu_from_couplings",
    "bplus_kplus_nunu_from_rs_semileptonic_wilsons",
    "b_to_kstar_nunu_branching_fraction",
    "b_to_kstar_nunu_from_couplings",
    "b_to_kstar_nunu_from_rs_semileptonic_wilsons",
    "RARE_B_NUNU_RS_SEMILEPTONIC_MATCHING_STATUS_V1",
    "RS_SEMILEPTONIC_MATCHING_ASSUMPTION_V1",
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


def _diagnostic_matching_scale(matching_scale_gev: float | None) -> float:
    if matching_scale_gev is None:
        return 0.0
    number = float(matching_scale_gev)
    if not math.isfinite(number) or number <= 0.0:
        raise ValueError("matching_scale_gev must be positive and finite")
    return number


def _rs_b_to_s_nunu_coeff(
    source: RSSemileptonicWilsonBundle,
) -> RSNuNuWilsonCoefficients:
    try:
        coeff = source.b_to_s_nunu
    except AttributeError as exc:
        raise TypeError("rs_semileptonic_wilsons.b_to_s_nunu is not available") from exc
    if coeff is None:
        raise ValueError("rs_semileptonic_wilsons.b_to_s_nunu is absent")
    if coeff.transition_key != "b_s":
        raise ValueError(
            "rs_semileptonic_wilsons.b_to_s_nunu transition_key="
            f"{coeff.transition_key!r}, expected 'b_s'"
        )
    if (
        coeff.quark_sector != "d"
        or int(coeff.final_quark_index) != 1
        or int(coeff.initial_quark_index) != 2
    ):
        raise ValueError("rs_semileptonic_wilsons.b_to_s_nunu has inconsistent quark indices")
    return coeff


def _flavor_universal_scalar(
    matrix: object,
    *,
    name: str,
) -> tuple[complex, float, float]:
    arr = np.asarray(matrix, dtype=np.complex128)
    if arr.shape != (3, 3):
        raise ValueError(f"{name} must have shape (3, 3)")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} contains non-finite entries")
    diag = np.diag(arr)
    scalar = complex(np.trace(arr) / 3.0)
    offdiag = np.array(arr, copy=True)
    np.fill_diagonal(offdiag, 0.0)
    max_offdiag = float(np.max(np.abs(offdiag)))
    max_diag_spread = float(np.max(np.abs(diag - scalar)))
    tolerance = _ACTIVE_NU_UNIVERSALITY_ATOL * max(1.0, abs(scalar))
    if max_offdiag > tolerance or max_diag_spread > tolerance:
        raise ValueError(
            f"{name} is not diagonal-universal enough for the scalar rare_b_nunu core"
        )
    return scalar, max_offdiag, max_diag_spread


def _rs_nunu_scalars(
    coeff: RSNuNuWilsonCoefficients,
) -> tuple[complex, complex, dict[str, float | complex]]:
    x_left, left_offdiag, left_diag_spread = _flavor_universal_scalar(
        coeff.x_np_left,
        name="b_to_s_nunu.x_np_left",
    )
    x_right, right_offdiag, right_diag_spread = _flavor_universal_scalar(
        coeff.x_np_right,
        name="b_to_s_nunu.x_np_right",
    )
    return (
        x_left,
        x_right,
        {
            "rs_semileptonic_x_np_left_trace_over_3": complex(x_left),
            "rs_semileptonic_x_np_right_trace_over_3": complex(x_right),
            "rs_semileptonic_x_np_total_trace_over_3": complex(x_left + x_right),
            "rs_semileptonic_x_np_left_max_offdiag_abs": float(left_offdiag),
            "rs_semileptonic_x_np_right_max_offdiag_abs": float(right_offdiag),
            "rs_semileptonic_x_np_left_max_diag_spread_abs": float(left_diag_spread),
            "rs_semileptonic_x_np_right_max_diag_spread_abs": float(right_diag_spread),
        },
    )


def rare_b_nunu_wilsons_from_rs_semileptonic(
    source: RSSemileptonicWilsonBundle,
    *,
    matching_scale_gev: float | None = None,
) -> RareBNuNuWilsonCoefficients:
    """Translate Phase-4a ``b_to_s_nunu`` contacts to the rare-B core shape."""

    coeff = _rs_b_to_s_nunu_coeff(source)
    x_left, x_right, _ = _rs_nunu_scalars(coeff)
    scale = _diagnostic_matching_scale(matching_scale_gev)
    contact_left = complex(x_left * coeff.g_sm_squared_gev_minus2)
    contact_right = complex(x_right * coeff.g_sm_squared_gev_minus2)
    return RareBNuNuWilsonCoefficients(
        model_label=coeff.model_label,
        operator_convention=coeff.operator_convention,
        matching_assumption=coeff.matching_assumption,
        M_KK=scale,
        matching_scale=scale,
        left_sb_coupling=contact_left,
        right_sb_coupling=contact_right,
        left_sb_overlap=0.0j,
        right_sb_overlap=0.0j,
        left_quark_delta=contact_left,
        right_quark_delta=contact_right,
        neutrino_delta=1.0,
        x_np_left=complex(x_left),
        x_np_right=complex(x_right),
    )


def rare_b_nunu_rs_semileptonic_diagnostics(
    coeff: RSNuNuWilsonCoefficients,
) -> dict[str, object]:
    """Return diagnostics for the Phase-4d ``b -> s nu nubar`` rewire."""

    _, _, scalar_diagnostics = _rs_nunu_scalars(coeff)
    return {
        "rs_semileptonic_model_label": coeff.model_label,
        "rs_semileptonic_operator_convention": coeff.operator_convention,
        "rs_semileptonic_matching_assumption": coeff.matching_assumption,
        "rs_semileptonic_wilsons_present": True,
        "rs_semileptonic_matching_status": (
            "rigorous_tree_light_z_nunu_from_rs_semileptonic_wilsons"
        ),
        "rs_semileptonic_nunu_matching_status": (
            RARE_B_NUNU_RS_SEMILEPTONIC_MATCHING_STATUS_V1
        ),
        "rs_semileptonic_transition_key": coeff.transition_key,
        "rs_semileptonic_quark_sector": coeff.quark_sector,
        "rs_semileptonic_final_quark_index": int(coeff.final_quark_index),
        "rs_semileptonic_initial_quark_index": int(coeff.initial_quark_index),
        "rs_semileptonic_contact_units": coeff.contact_units,
        "rs_semileptonic_g_sm_squared_gev_minus2": float(
            coeff.g_sm_squared_gev_minus2
        ),
        "rs_semileptonic_contacts_trace_over_3": {
            "C_LL": complex(np.trace(coeff.contact_LL) / 3.0),
            "C_RL": complex(np.trace(coeff.contact_RL) / 3.0),
        },
        "rs_semileptonic_wilson_coefficients_trace_over_3": {
            "X_NP_L": scalar_diagnostics["rs_semileptonic_x_np_left_trace_over_3"],
            "X_NP_R": scalar_diagnostics["rs_semileptonic_x_np_right_trace_over_3"],
            "X_NP_total": scalar_diagnostics[
                "rs_semileptonic_x_np_total_trace_over_3"
            ],
        },
        "nunu_mapping": "X_NP=C/g_SM^2",
        "active_neutrino_current": "LH",
        "majorana_dirac_rate_factor": 1.0,
        "wilson_prefactor_reused": False,
        "second_mkk_suppression_applied": False,
        "legacy_one_z_proxy_reused": False,
        **scalar_diagnostics,
    }


def _tag_rs_nunu_result(result, coeff: RSNuNuWilsonCoefficients):
    diagnostics = dict(result.diagnostics)
    diagnostics.update(rare_b_nunu_rs_semileptonic_diagnostics(coeff))
    diagnostics["matching_assumption"] = coeff.matching_assumption
    diagnostics["rs_matching_assumption"] = coeff.matching_assumption
    diagnostics["nunu_rs_semileptonic_rewired"] = True
    return replace(result, diagnostics=diagnostics)


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


def bplus_kplus_nunu_from_rs_semileptonic_wilsons(
    source: RSSemileptonicWilsonBundle,
    *,
    matching_scale_gev: float | None = None,
    inputs: RareBNuNuSMInputs | None = None,
) -> RareBNuNuBranchingResult:
    """Evaluate ``BR(B+ -> K+ nu nubar)`` from Phase-4a RS νν Wilsons."""

    coeff = _rs_b_to_s_nunu_coeff(source)
    wilsons = rare_b_nunu_wilsons_from_rs_semileptonic(
        source,
        matching_scale_gev=matching_scale_gev,
    )
    result = _evaluate_bplus_to_kplus_nunu(wilsons, inputs=inputs)
    return _tag_rs_nunu_result(result, coeff)


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


def b_to_kstar_nunu_from_rs_semileptonic_wilsons(
    source: RSSemileptonicWilsonBundle,
    *,
    sm_branching_fraction: float,
    matching_scale_gev: float | None = None,
    inputs: RareBNuNuSMInputs | None = None,
) -> RareBKstarNuNuBranchingResult:
    """Evaluate ``BR(B -> K* nu nubar)`` from Phase-4a RS νν Wilsons."""

    coeff = _rs_b_to_s_nunu_coeff(source)
    wilsons = rare_b_nunu_wilsons_from_rs_semileptonic(
        source,
        matching_scale_gev=matching_scale_gev,
    )
    result = b_to_kstar_nunu_branching_fraction(
        wilsons,
        sm_branching_fraction=sm_branching_fraction,
        inputs=inputs,
    )
    return _tag_rs_nunu_result(result, coeff)
