"""Adapter over :mod:`quarkConstraints.rare_kaon_snd`.

This is the catalog boundary for the Delta-S=1 ``K -> pi nu nubar``
machinery.  Constraint modules import this adapter only; the underlying physics
implementation remains isolated in ``quarkConstraints``.  K005 intentionally
appends the ``K_L -> pi0 nu nubar`` machinery to the K004
``quarkConstraints/rare_kaon_snd.py`` module and wraps it here; this is the
designed location for this physics, not an isolation violation.  Production
Phase-4d constraints consume the rigorous Phase-4a
``rs_semileptonic_wilsons.s_to_d_nunu`` block and map ``X_NP=C/g_SM^2``
directly into the shared Buras/BGS rare-kaon core.
"""

from __future__ import annotations

import math
from dataclasses import replace

import numpy as np

from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.rare_kaon_snd import (
    RARE_KAON_INPUT_BUNDLE_V1,
    RARE_KAON_KAPPA_L_CITATION,
    RARE_KAON_KAPPA_L_REF,
    RARE_KAON_OPERATOR_CONVENTION,
    RARE_KAON_RS_MATCHING_ASSUMPTION_V1,
    RARE_KAON_SND_MODEL_V1,
    RareKaonBranchingResult,
    RareKaonCKMFactors,
    RareKaonNeutralBranchingResult,
    RareKaonSMInputs,
    RareKaonWilsonCoefficients,
)
from quarkConstraints.rare_kaon_snd import (
    ckm_factors as _ckm_factors,
)
from quarkConstraints.rare_kaon_snd import (
    compute_rare_kaon_wilsons as _compute_rare_kaon_wilsons,
)
from quarkConstraints.rare_kaon_snd import (
    default_sm_inputs as _default_sm_inputs,
)
from quarkConstraints.rare_kaon_snd import (
    evaluate_klong_to_pi0_nunu as _evaluate_klong_to_pi0_nunu,
)
from quarkConstraints.rare_kaon_snd import (
    evaluate_kplus_to_piplus_nunu as _evaluate_kplus_to_piplus_nunu,
)
from quarkConstraints.rare_kaon_snd import (
    g_sm_squared as _g_sm_squared,
)
from quarkConstraints.rare_kaon_snd import (
    kappa_l as _kappa_l,
)
from quarkConstraints.rare_kaon_snd import (
    kappa_plus as _kappa_plus,
)
from quarkConstraints.rare_kaon_snd import (
    neutral_sm_branching_fraction as _neutral_sm_branching_fraction,
)
from quarkConstraints.rare_kaon_snd import (
    sm_branching_fraction as _sm_branching_fraction,
)
from quarkConstraints.rs_semileptonic_wilsons import (
    RS_SEMILEPTONIC_MATCHING_ASSUMPTION_V1,
    RSNuNuWilsonCoefficients,
    RSSemileptonicWilsonBundle,
)

RARE_KAON_NUNU_RS_SEMILEPTONIC_MATCHING_STATUS_V1 = (
    "Phase-4d light-Z RS active-neutrino Wilsons consumed as X_NP=C/g_SM^2; "
    "no rare_kaon_snd proxy matcher, no charged-lepton _wilson_prefactor, and "
    "no second 1/M_KK^2 factor."
)
_ACTIVE_NU_UNIVERSALITY_ATOL = 1.0e-10

__all__ = [
    "QuarkMassBasisCouplings",
    "RARE_KAON_SND_MODEL_V1",
    "RARE_KAON_OPERATOR_CONVENTION",
    "RARE_KAON_INPUT_BUNDLE_V1",
    "RARE_KAON_RS_MATCHING_ASSUMPTION_V1",
    "RARE_KAON_KAPPA_L_REF",
    "RARE_KAON_KAPPA_L_CITATION",
    "RareKaonSMInputs",
    "RareKaonCKMFactors",
    "RareKaonWilsonCoefficients",
    "RareKaonBranchingResult",
    "RareKaonNeutralBranchingResult",
    "rare_kaon_default_sm_inputs",
    "rare_kaon_ckm_factors",
    "rare_kaon_kappa_plus",
    "rare_kaon_kappa_l",
    "rare_kaon_g_sm_squared",
    "rare_kaon_sm_branching_fraction",
    "rare_kaon_neutral_sm_branching_fraction",
    "rare_kaon_wilsons_from_couplings",
    "rare_kaon_wilsons_from_rs_semileptonic",
    "rare_kaon_nunu_rs_semileptonic_diagnostics",
    "kplus_piplus_nunu_from_couplings",
    "kplus_piplus_nunu_from_rs_semileptonic_wilsons",
    "klong_pi0_nunu_from_couplings",
    "klong_pi0_nunu_from_rs_semileptonic_wilsons",
    "RARE_KAON_NUNU_RS_SEMILEPTONIC_MATCHING_STATUS_V1",
    "RS_SEMILEPTONIC_MATCHING_ASSUMPTION_V1",
]


def rare_kaon_default_sm_inputs() -> RareKaonSMInputs:
    """Return the default rare-kaon SM input bundle."""
    return _default_sm_inputs()


def rare_kaon_ckm_factors(
    inputs: RareKaonSMInputs | None = None,
) -> RareKaonCKMFactors:
    """Return the CKM factors used by the rare-kaon core."""
    return _ckm_factors(inputs)


def rare_kaon_kappa_plus(inputs: RareKaonSMInputs | None = None) -> float:
    """Return the charged-mode hadronic factor ``kappa_+``."""
    return _kappa_plus(inputs)


def rare_kaon_kappa_l(inputs: RareKaonSMInputs | None = None) -> float:
    """Return the neutral-mode hadronic factor ``kappa_L``."""
    return _kappa_l(inputs)


def rare_kaon_g_sm_squared(inputs: RareKaonSMInputs | None = None) -> float:
    """Return the Buras ``g_SM^2`` normalization."""
    return _g_sm_squared(inputs)


def rare_kaon_sm_branching_fraction(
    inputs: RareKaonSMInputs | None = None,
) -> RareKaonBranchingResult:
    """Evaluate the SM-limit charged rare-kaon branching fraction."""
    return _sm_branching_fraction(inputs)


def rare_kaon_neutral_sm_branching_fraction(
    inputs: RareKaonSMInputs | None = None,
) -> RareKaonNeutralBranchingResult:
    """Evaluate the SM-limit neutral rare-kaon branching fraction."""
    return _neutral_sm_branching_fraction(inputs)


def rare_kaon_wilsons_from_couplings(
    couplings: QuarkMassBasisCouplings,
    *,
    m_kk_gev: float | None = None,
    inputs: RareKaonSMInputs | None = None,
) -> RareKaonWilsonCoefficients:
    """Return the v1 ``s -> d nu nubar`` Wilson proxy for mass-basis couplings."""
    return _compute_rare_kaon_wilsons(
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


def _rs_s_to_d_nunu_coeff(
    source: RSSemileptonicWilsonBundle,
) -> RSNuNuWilsonCoefficients:
    try:
        coeff = source.s_to_d_nunu
    except AttributeError as exc:
        raise TypeError("rs_semileptonic_wilsons.s_to_d_nunu is not available") from exc
    if coeff is None:
        raise ValueError("rs_semileptonic_wilsons.s_to_d_nunu is absent")
    if coeff.transition_key != "s_d":
        raise ValueError(
            "rs_semileptonic_wilsons.s_to_d_nunu transition_key="
            f"{coeff.transition_key!r}, expected 's_d'"
        )
    if (
        coeff.quark_sector != "d"
        or int(coeff.final_quark_index) != 0
        or int(coeff.initial_quark_index) != 1
    ):
        raise ValueError("rs_semileptonic_wilsons.s_to_d_nunu has inconsistent quark indices")
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
            f"{name} is not diagonal-universal enough for the scalar rare_kaon_snd core"
        )
    return scalar, max_offdiag, max_diag_spread


def _rs_nunu_scalars(
    coeff: RSNuNuWilsonCoefficients,
) -> tuple[complex, complex, dict[str, float | complex]]:
    x_left, left_offdiag, left_diag_spread = _flavor_universal_scalar(
        coeff.x_np_left,
        name="s_to_d_nunu.x_np_left",
    )
    x_right, right_offdiag, right_diag_spread = _flavor_universal_scalar(
        coeff.x_np_right,
        name="s_to_d_nunu.x_np_right",
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


def rare_kaon_wilsons_from_rs_semileptonic(
    source: RSSemileptonicWilsonBundle,
    *,
    matching_scale_gev: float | None = None,
) -> RareKaonWilsonCoefficients:
    """Translate Phase-4a ``s_to_d_nunu`` contacts to the rare-kaon core shape."""

    coeff = _rs_s_to_d_nunu_coeff(source)
    x_left, x_right, _ = _rs_nunu_scalars(coeff)
    scale = _diagnostic_matching_scale(matching_scale_gev)
    contact_left = complex(x_left * coeff.g_sm_squared_gev_minus2)
    contact_right = complex(x_right * coeff.g_sm_squared_gev_minus2)
    return RareKaonWilsonCoefficients(
        model_label=coeff.model_label,
        operator_convention=coeff.operator_convention,
        matching_assumption=coeff.matching_assumption,
        M_KK=scale,
        matching_scale=scale,
        left_sd_coupling=contact_left,
        right_sd_coupling=contact_right,
        left_sd_overlap=0.0j,
        right_sd_overlap=0.0j,
        left_quark_delta=contact_left,
        right_quark_delta=contact_right,
        neutrino_delta=1.0,
        x_np_left=complex(x_left),
        x_np_right=complex(x_right),
    )


def rare_kaon_nunu_rs_semileptonic_diagnostics(
    coeff: RSNuNuWilsonCoefficients,
) -> dict[str, object]:
    """Return diagnostics for the Phase-4d ``s -> d nu nubar`` rewire."""

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
            RARE_KAON_NUNU_RS_SEMILEPTONIC_MATCHING_STATUS_V1
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
    diagnostics.update(rare_kaon_nunu_rs_semileptonic_diagnostics(coeff))
    diagnostics["matching_assumption"] = coeff.matching_assumption
    diagnostics["rs_matching_assumption"] = coeff.matching_assumption
    diagnostics["nunu_rs_semileptonic_rewired"] = True
    return replace(result, diagnostics=diagnostics)


def kplus_piplus_nunu_from_couplings(
    couplings: QuarkMassBasisCouplings,
    *,
    m_kk_gev: float | None = None,
    inputs: RareKaonSMInputs | None = None,
) -> RareKaonBranchingResult:
    """Evaluate ``BR(K+ -> pi+ nu nubar)`` from mass-basis couplings."""
    return _evaluate_kplus_to_piplus_nunu(
        couplings,
        m_kk_gev=m_kk_gev,
        inputs=inputs,
    )


def kplus_piplus_nunu_from_rs_semileptonic_wilsons(
    source: RSSemileptonicWilsonBundle,
    *,
    matching_scale_gev: float | None = None,
    inputs: RareKaonSMInputs | None = None,
) -> RareKaonBranchingResult:
    """Evaluate ``BR(K+ -> pi+ nu nubar)`` from Phase-4a RS νν Wilsons."""

    coeff = _rs_s_to_d_nunu_coeff(source)
    wilsons = rare_kaon_wilsons_from_rs_semileptonic(
        source,
        matching_scale_gev=matching_scale_gev,
    )
    result = _evaluate_kplus_to_piplus_nunu(wilsons, inputs=inputs)
    return _tag_rs_nunu_result(result, coeff)


def klong_pi0_nunu_from_couplings(
    couplings: QuarkMassBasisCouplings,
    *,
    m_kk_gev: float | None = None,
    inputs: RareKaonSMInputs | None = None,
) -> RareKaonNeutralBranchingResult:
    """Evaluate ``BR(K_L -> pi0 nu nubar)`` from mass-basis couplings."""
    return _evaluate_klong_to_pi0_nunu(
        couplings,
        m_kk_gev=m_kk_gev,
        inputs=inputs,
    )


def klong_pi0_nunu_from_rs_semileptonic_wilsons(
    source: RSSemileptonicWilsonBundle,
    *,
    matching_scale_gev: float | None = None,
    inputs: RareKaonSMInputs | None = None,
) -> RareKaonNeutralBranchingResult:
    """Evaluate ``BR(K_L -> pi0 nu nubar)`` from Phase-4a RS νν Wilsons."""

    coeff = _rs_s_to_d_nunu_coeff(source)
    wilsons = rare_kaon_wilsons_from_rs_semileptonic(
        source,
        matching_scale_gev=matching_scale_gev,
    )
    result = _evaluate_klong_to_pi0_nunu(wilsons, inputs=inputs)
    return _tag_rs_nunu_result(result, coeff)
