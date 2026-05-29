"""Adapter over :mod:`quarkConstraints.collider_resonance`.

Constraint modules import this adapter only.  The core collider-resonance
comparison code is kept in ``quarkConstraints.collider_resonance`` so future
collider-RS constraints can reuse the same mass-limit and ``sigma*BR`` paths.
"""

from __future__ import annotations

from typing import Any

from quarkConstraints.collider_resonance import (
    COLLIDER_RESONANCE_MASS_PROXY_ASSUMPTION_V1,
    MASS_LOWER_BOUND,
    SIGMA_TIMES_BR_UPPER_LIMIT,
    ColliderResonanceComparison,
    ColliderResonanceLimit,
    ColliderResonancePrediction,
    evaluate_resonance_limit as _evaluate_resonance_limit,
    kk_gluon_prediction_from_m_kk_gev as _kk_gluon_prediction_from_m_kk_gev,
    kk_mass_tev_from_m_kk_gev as _kk_mass_tev_from_m_kk_gev,
    mass_from_source_gev as _mass_from_source_gev,
)

__all__ = [
    "MASS_LOWER_BOUND",
    "SIGMA_TIMES_BR_UPPER_LIMIT",
    "COLLIDER_RESONANCE_MASS_PROXY_ASSUMPTION_V1",
    "ColliderResonanceLimit",
    "ColliderResonancePrediction",
    "ColliderResonanceComparison",
    "kk_gluon_mass_tev_from_m_kk_gev",
    "kk_gluon_prediction_from_m_kk_gev",
    "resolve_kk_gluon_mass_gev",
    "evaluate_collider_resonance_limit",
]


def kk_gluon_mass_tev_from_m_kk_gev(m_kk_gev: float) -> float:
    """Return the positive finite KK-gluon mass in TeV."""

    return _kk_mass_tev_from_m_kk_gev(m_kk_gev)


def kk_gluon_prediction_from_m_kk_gev(
    m_kk_gev: float,
    *,
    sigma_times_br: float | None = None,
    sigma_times_br_units: str | None = None,
) -> ColliderResonancePrediction:
    """Build the v1 KK-gluon resonance prediction from ``M_KK``."""

    return _kk_gluon_prediction_from_m_kk_gev(
        m_kk_gev,
        sigma_times_br=sigma_times_br,
        sigma_times_br_units=sigma_times_br_units,
    )


def resolve_kk_gluon_mass_gev(
    *,
    mass_extra: Any = None,
    couplings: Any = None,
) -> tuple[float | None, str | None]:
    """Resolve the KK-gluon mass in GeV from declared point extras.

    The explicit ``kk_gluon_mass_gev`` extra wins.  If absent, fall back to
    ``quark_mass_basis_couplings.M_KK`` because the point builder already uses
    that convention for quark-sector KK-gluon couplings.
    """

    if mass_extra is not None:
        mass = float(mass_extra)
        _kk_mass_tev_from_m_kk_gev(mass)
        return mass, "kk_gluon_mass_gev"
    if couplings is not None:
        return float(_mass_from_source_gev(couplings)), "quark_mass_basis_couplings.M_KK"
    return None, None


def evaluate_collider_resonance_limit(
    prediction: ColliderResonancePrediction,
    limit: ColliderResonanceLimit,
) -> ColliderResonanceComparison:
    """Compare a resonance prediction to an experimental limit."""

    return _evaluate_resonance_limit(prediction, limit)


VLQ_PAIR_MASS_PROXY_ASSUMPTION_V1 = (
    "NEEDS-HUMAN-PHYSICS: vector-like-quark pair-production recast v1 uses "
    "m_VLQ = M_KK as a KK-fermion mass proxy and compares it to the "
    "catalogued benchmark mass lower limit. It does not compute "
    "sigma(pp->VLQ VLQ)*BR(VLQ->tW)^2, widths, branching-ratio mixtures, "
    "acceptance, or the experiment's mass-dependent limit curve."
)


def resolve_vlq_mkk_gev(*, couplings: Any = None) -> tuple[float | None, str | None]:
    """Resolve the KK-fermion proxy mass in GeV from declared point extras.

    The current scan point carries the common ``M_KK`` scale on the quark
    mass-basis couplings object.  CR002 interprets that scale as the
    vector-like-quark mass proxy and flags the matching as
    ``NEEDS-HUMAN-PHYSICS`` at the constraint result.
    """

    if couplings is not None:
        return float(_mass_from_source_gev(couplings)), "quark_mass_basis_couplings.M_KK"
    return None, None


def vlq_pair_prediction_from_m_kk_gev(
    m_kk_gev: float,
    *,
    resonance: str = "T_5/3 pair",
    final_state: str = "tW tW",
    sigma_times_br: float | None = None,
    sigma_times_br_units: str | None = None,
) -> ColliderResonancePrediction:
    """Build the documented v1 VLQ pair-production mass-proxy prediction."""

    return ColliderResonancePrediction(
        resonance=resonance,
        final_state=final_state,
        mass_tev=_kk_mass_tev_from_m_kk_gev(m_kk_gev),
        sigma_times_br=sigma_times_br,
        sigma_times_br_units=sigma_times_br_units,
        matching_assumption=VLQ_PAIR_MASS_PROXY_ASSUMPTION_V1,
        diagnostics={
            "m_kk_gev": float(m_kk_gev),
            "m_vlq_proxy_gev": float(m_kk_gev),
            "mass_proxy": "m_VLQ = M_KK",
            "sigma_times_br_proxy_available": sigma_times_br is not None,
        },
    )


__all__.extend(
    [
        "VLQ_PAIR_MASS_PROXY_ASSUMPTION_V1",
        "resolve_vlq_mkk_gev",
        "vlq_pair_prediction_from_m_kk_gev",
    ]
)


VLQ_T_PAIR_MASS_PROXY_ASSUMPTION_V1 = (
    "NEEDS-HUMAN-PHYSICS: charge-2/3 vector-like T pair-production recast v1 "
    "uses m_T = M_KK as a KK-fermion mass proxy and compares it to the "
    "catalogued benchmark mass lower limit. It does not compute "
    "sigma(pp->T Tbar)*BR(T->Wb/Zt/Ht)^2, widths, branching-fraction "
    "mixtures, acceptance, or the experiment's mass-dependent limit curve."
)


def resolve_charge_two_thirds_vlq_mkk_gev(
    *,
    couplings: Any = None,
) -> tuple[float | None, str | None]:
    """Resolve the charge-2/3 top-partner proxy mass in GeV.

    The current scan point carries the common quark-sector ``M_KK`` scale, not
    a dedicated custodial charge-2/3 VLQ spectrum.  CR003 therefore interprets
    ``M_KK`` as the top-partner mass proxy and flags the matching as
    ``NEEDS-HUMAN-PHYSICS`` at the result.
    """

    if couplings is not None:
        return float(_mass_from_source_gev(couplings)), "quark_mass_basis_couplings.M_KK"
    return None, None


def vlq_t_pair_prediction_from_m_kk_gev(
    m_kk_gev: float,
    *,
    resonance: str = "T pair",
    final_state: str = "Wb/Zt/Ht",
    sigma_times_br: float | None = None,
    sigma_times_br_units: str | None = None,
) -> ColliderResonancePrediction:
    """Build the documented v1 charge-2/3 VLQ pair-production mass proxy."""

    return ColliderResonancePrediction(
        resonance=resonance,
        final_state=final_state,
        mass_tev=_kk_mass_tev_from_m_kk_gev(m_kk_gev),
        sigma_times_br=sigma_times_br,
        sigma_times_br_units=sigma_times_br_units,
        matching_assumption=VLQ_T_PAIR_MASS_PROXY_ASSUMPTION_V1,
        diagnostics={
            "m_kk_gev": float(m_kk_gev),
            "m_t_partner_proxy_gev": float(m_kk_gev),
            "mass_proxy": "m_T = M_KK",
            "sigma_times_br_proxy_available": sigma_times_br is not None,
        },
    )


__all__.extend(
    [
        "VLQ_T_PAIR_MASS_PROXY_ASSUMPTION_V1",
        "resolve_charge_two_thirds_vlq_mkk_gev",
        "vlq_t_pair_prediction_from_m_kk_gev",
    ]
)


VLQ_B_PAIR_MASS_PROXY_ASSUMPTION_V1 = (
    "NEEDS-HUMAN-PHYSICS: charge -1/3 vector-like B pair-production recast "
    "v1 uses m_B = M_KK as a KK-fermion mass proxy and compares it to the "
    "catalogued benchmark mass lower limit. It does not compute "
    "sigma(pp->B Bbar)*BR(B->tW/bZ/bH)^2, widths, branching-fraction "
    "mixtures, acceptance, or the experiment's mass-dependent limit curve."
)


def resolve_charge_minus_one_third_vlq_mkk_gev(
    *,
    couplings: Any = None,
) -> tuple[float | None, str | None]:
    """Resolve the charge -1/3 bottom-partner proxy mass in GeV.

    The current scan point carries the common quark-sector ``M_KK`` scale, not
    a dedicated custodial charge -1/3 ``B`` spectrum.  CR004 therefore
    interprets ``M_KK`` as the bottom-partner mass proxy and flags the matching
    as ``NEEDS-HUMAN-PHYSICS`` at the result.
    """

    if couplings is not None:
        return float(_mass_from_source_gev(couplings)), "quark_mass_basis_couplings.M_KK"
    return None, None


def vlq_b_pair_prediction_from_m_kk_gev(
    m_kk_gev: float,
    *,
    resonance: str = "B pair",
    final_state: str = "tW/bZ/bH",
    sigma_times_br: float | None = None,
    sigma_times_br_units: str | None = None,
) -> ColliderResonancePrediction:
    """Build the documented v1 charge -1/3 VLQ pair-production mass proxy."""

    return ColliderResonancePrediction(
        resonance=resonance,
        final_state=final_state,
        mass_tev=_kk_mass_tev_from_m_kk_gev(m_kk_gev),
        sigma_times_br=sigma_times_br,
        sigma_times_br_units=sigma_times_br_units,
        matching_assumption=VLQ_B_PAIR_MASS_PROXY_ASSUMPTION_V1,
        diagnostics={
            "m_kk_gev": float(m_kk_gev),
            "m_b_partner_proxy_gev": float(m_kk_gev),
            "mass_proxy": "m_B = M_KK",
            "sigma_times_br_proxy_available": sigma_times_br is not None,
        },
    )


__all__.extend(
    [
        "VLQ_B_PAIR_MASS_PROXY_ASSUMPTION_V1",
        "resolve_charge_minus_one_third_vlq_mkk_gev",
        "vlq_b_pair_prediction_from_m_kk_gev",
    ]
)


KK_EW_DILEPTON_MASS_PROXY_ASSUMPTION_V1 = (
    "NEEDS-HUMAN-PHYSICS: neutral electroweak KK dilepton recast v1 uses "
    "the supplied kk_ew_mass_gev, or M_KK as a fallback mass proxy, and "
    "compares it to a catalogued benchmark spin-1 dilepton mass lower bound. "
    "It does not compute sigma(pp->gamma_KK/Z_KK)*BR(gamma_KK/Z_KK->ll), "
    "light-quark and charged-lepton couplings, total width, interference, "
    "acceptance, or the experiment's mass-dependent limit curve."
)


def resolve_kk_ew_mass_gev(
    *,
    mass_extra: Any = None,
    couplings: Any = None,
) -> tuple[float | None, str | None]:
    """Resolve the neutral electroweak KK mass in GeV from declared extras.

    The explicit ``kk_ew_mass_gev`` extra wins.  If absent, fall back to the
    common quark-sector ``M_KK`` carried by the mass-basis couplings object;
    the resulting collider interpretation is a documented mass proxy and is
    flagged as ``NEEDS-HUMAN-PHYSICS`` by CR005.
    """

    if mass_extra is not None:
        mass = float(mass_extra)
        _kk_mass_tev_from_m_kk_gev(mass)
        return mass, "kk_ew_mass_gev"
    if couplings is not None:
        return float(_mass_from_source_gev(couplings)), "quark_mass_basis_couplings.M_KK"
    return None, None


def kk_ew_dilepton_prediction_from_m_kk_gev(
    m_kk_gev: float,
    *,
    resonance: str = "(gamma^(1), Z^(1))_KK",
    final_state: str = "ee + mumu",
    sigma_times_br: float | None = None,
    sigma_times_br_units: str | None = None,
) -> ColliderResonancePrediction:
    """Build the documented v1 neutral-EW KK dilepton mass proxy."""

    return ColliderResonancePrediction(
        resonance=resonance,
        final_state=final_state,
        mass_tev=_kk_mass_tev_from_m_kk_gev(m_kk_gev),
        sigma_times_br=sigma_times_br,
        sigma_times_br_units=sigma_times_br_units,
        matching_assumption=KK_EW_DILEPTON_MASS_PROXY_ASSUMPTION_V1,
        diagnostics={
            "m_kk_gev": float(m_kk_gev),
            "m_kk_ew_proxy_gev": float(m_kk_gev),
            "mass_proxy": "m_(gamma/Z KK) = kk_ew_mass_gev or M_KK",
            "sigma_times_br_proxy_available": sigma_times_br is not None,
        },
    )


__all__.extend(
    [
        "KK_EW_DILEPTON_MASS_PROXY_ASSUMPTION_V1",
        "resolve_kk_ew_mass_gev",
        "kk_ew_dilepton_prediction_from_m_kk_gev",
    ]
)
