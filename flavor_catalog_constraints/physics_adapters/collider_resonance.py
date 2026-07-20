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
)
from quarkConstraints.collider_resonance import (
    evaluate_resonance_limit as _evaluate_resonance_limit,
)
from quarkConstraints.collider_resonance import (
    kk_gluon_prediction_from_m_kk_gev as _kk_gluon_prediction_from_m_kk_gev,
)
from quarkConstraints.collider_resonance import (
    kk_mass_tev_from_m_kk_gev as _kk_mass_tev_from_m_kk_gev,
)
from quarkConstraints.collider_resonance import (
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
    "reports an advisory comparison to a catalogued benchmark spin-1 dilepton "
    "mass lower bound. Per M-26, raw SSM edges must not HARD-veto "
    "sqrt(L)-suppressed RS KK electroweak bosons. It does not compute "
    "sigma(pp->gamma_KK/Z_KK)*BR(gamma_KK/Z_KK->ll), light-quark and "
    "charged-lepton couplings, total width, interference, acceptance, or the "
    "experiment's mass-dependent limit curve."
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


KK_CHARGED_CURRENT_MASS_PROXY_ASSUMPTION_V1 = (
    "NEEDS-HUMAN-PHYSICS: charged electroweak KK/W' recast v1 uses the "
    "supplied kk_ew_mass_gev, or M_KK as a fallback mass proxy, and reports "
    "an advisory comparison to a catalogued benchmark W' mass lower bound. Per "
    "M-26, raw SSM edges must not HARD-veto sqrt(L)-suppressed RS KK "
    "electroweak bosons. It does not compute sigma(pp->W_KK)*BR(W_KK->l nu,tb), "
    "light-quark, lepton, or tb couplings, total width, interference, "
    "acceptance, or the experiment's mass-dependent limit curve."
)


def kk_charged_current_prediction_from_m_kk_gev(
    m_kk_gev: float,
    *,
    resonance: str = "W_KK^(1)",
    final_state: str = "ell nu",
    sigma_times_br: float | None = None,
    sigma_times_br_units: str | None = None,
) -> ColliderResonancePrediction:
    """Build the documented v1 charged-current EW KK mass proxy."""

    return ColliderResonancePrediction(
        resonance=resonance,
        final_state=final_state,
        mass_tev=_kk_mass_tev_from_m_kk_gev(m_kk_gev),
        sigma_times_br=sigma_times_br,
        sigma_times_br_units=sigma_times_br_units,
        matching_assumption=KK_CHARGED_CURRENT_MASS_PROXY_ASSUMPTION_V1,
        diagnostics={
            "m_kk_gev": float(m_kk_gev),
            "m_kk_charged_current_proxy_gev": float(m_kk_gev),
            "mass_proxy": "m_WKK = kk_ew_mass_gev or M_KK",
            "sigma_times_br_proxy_available": sigma_times_br is not None,
        },
    )


__all__.extend(
    [
        "KK_CHARGED_CURRENT_MASS_PROXY_ASSUMPTION_V1",
        "kk_charged_current_prediction_from_m_kk_gev",
    ]
)


DY_CONTACT_OPERATOR_PROXY_ASSUMPTION_V1 = (
    "NEEDS-HUMAN-PHYSICS: high-mass Drell-Yan contact-operator recast v1 "
    "uses the supplied kk_ew_mass_gev, kk_gluon_mass_gev, or M_KK as a "
    "Lambda_RS contact-scale proxy and compares it to the catalogued llqq "
    "contact-interaction lower limit. It does not compute "
    "4*pi/Lambda^2 ~= |g_q g_l|/M_V^2, helicity currents, interference sign, "
    "light-quark and charged-lepton couplings, total widths, EFT validity, "
    "PDF/electroweak systematics, or a binned high-mass dilepton likelihood."
)


def resolve_dy_contact_scale_gev(
    *,
    ew_mass_extra: Any = None,
    gluon_mass_extra: Any = None,
    couplings: Any = None,
) -> tuple[float | None, str | None]:
    """Resolve the CR009 ``Lambda_RS`` contact-scale proxy in GeV.

    The neutral electroweak KK mass is the closest declared scale for a
    dilepton contact operator.  If it is absent, the explicit KK-gluon scale
    or the common quark-sector ``M_KK`` are used only as documented RS scale
    proxies and are flagged by the CR009 result.
    """

    if ew_mass_extra is not None:
        mass = float(ew_mass_extra)
        _kk_mass_tev_from_m_kk_gev(mass)
        return mass, "kk_ew_mass_gev"
    if gluon_mass_extra is not None:
        mass = float(gluon_mass_extra)
        _kk_mass_tev_from_m_kk_gev(mass)
        return mass, "kk_gluon_mass_gev"
    if couplings is not None:
        return float(_mass_from_source_gev(couplings)), "quark_mass_basis_couplings.M_KK"
    return None, None


def dy_contact_prediction_from_scale_gev(
    lambda_proxy_gev: float,
    *,
    resonance: str = "llqq contact operator",
    final_state: str = "ee + mumu high-mass tail",
    sigma_times_br: float | None = None,
    sigma_times_br_units: str | None = None,
) -> ColliderResonancePrediction:
    """Build the documented CR009 contact-scale lower-bound proxy."""

    return ColliderResonancePrediction(
        resonance=resonance,
        final_state=final_state,
        mass_tev=_kk_mass_tev_from_m_kk_gev(lambda_proxy_gev),
        sigma_times_br=sigma_times_br,
        sigma_times_br_units=sigma_times_br_units,
        matching_assumption=DY_CONTACT_OPERATOR_PROXY_ASSUMPTION_V1,
        diagnostics={
            "lambda_rs_proxy_gev": float(lambda_proxy_gev),
            "contact_scale_proxy": (
                "Lambda_RS = kk_ew_mass_gev, kk_gluon_mass_gev, or M_KK"
            ),
            "sigma_or_eft_proxy_available": sigma_times_br is not None,
        },
    )


__all__.extend(
    [
        "DY_CONTACT_OPERATOR_PROXY_ASSUMPTION_V1",
        "resolve_dy_contact_scale_gev",
        "dy_contact_prediction_from_scale_gev",
    ]
)


VLQ_TB_DOUBLET_MASS_PROXY_ASSUMPTION_V1 = (
    "NEEDS-HUMAN-PHYSICS: weak-isospin (T,B) vector-like-quark doublet "
    "pair-production recast v1 uses m_T = m_B = M_KK as a common KK-fermion "
    "mass proxy and compares it to the catalogued simultaneous doublet mass "
    "lower limit. It does not compute sigma(pp->T Tbar/B Bbar)*BR mixtures, "
    "relative T/B rates, widths, mass splittings, acceptance, nonstandard "
    "cascade decays, or the experiment's mass-dependent limit surface."
)


def resolve_vlq_tb_doublet_mkk_gev(
    *,
    couplings: Any = None,
) -> tuple[float | None, str | None]:
    """Resolve the common (T,B)-doublet proxy mass in GeV.

    The current scan point carries the common quark-sector ``M_KK`` scale, not
    a dedicated weak-isospin ``(T,B)`` vector-like-quark spectrum.  CR010
    therefore interprets ``M_KK`` as the degenerate ``m_T=m_B`` mass proxy and
    flags the matching as ``NEEDS-HUMAN-PHYSICS`` at the result.
    """

    if couplings is not None:
        return float(_mass_from_source_gev(couplings)), "quark_mass_basis_couplings.M_KK"
    return None, None


def vlq_tb_doublet_prediction_from_m_kk_gev(
    m_kk_gev: float,
    *,
    resonance: str = "(T,B) doublet pair",
    final_state: str = "mixed W/Z/H third-generation final states",
    sigma_times_br: float | None = None,
    sigma_times_br_units: str | None = None,
) -> ColliderResonancePrediction:
    """Build the documented v1 weak-isospin (T,B)-doublet mass proxy."""

    return ColliderResonancePrediction(
        resonance=resonance,
        final_state=final_state,
        mass_tev=_kk_mass_tev_from_m_kk_gev(m_kk_gev),
        sigma_times_br=sigma_times_br,
        sigma_times_br_units=sigma_times_br_units,
        matching_assumption=VLQ_TB_DOUBLET_MASS_PROXY_ASSUMPTION_V1,
        diagnostics={
            "m_kk_gev": float(m_kk_gev),
            "m_t_doublet_proxy_gev": float(m_kk_gev),
            "m_b_doublet_proxy_gev": float(m_kk_gev),
            "mass_proxy": "m_T = m_B = M_KK",
            "sigma_times_br_proxy_available": sigma_times_br is not None,
        },
    )


__all__.extend(
    [
        "VLQ_TB_DOUBLET_MASS_PROXY_ASSUMPTION_V1",
        "resolve_vlq_tb_doublet_mkk_gev",
        "vlq_tb_doublet_prediction_from_m_kk_gev",
    ]
)


KK_DIBOSON_SPIN1_MASS_PROXY_ASSUMPTION_V1 = (
    "NEEDS-HUMAN-PHYSICS: spin-1 diboson resonance recast v1 uses the "
    "supplied kk_ew_mass_gev, or M_KK as a fallback electroweak-vector mass "
    "proxy, and compares it to a catalogued HVT model-B W'->WZ benchmark mass "
    "lower bound. It does not compute sigma(pp->W'/Z'/V_KK)*BR(V->WW/WZ/ZZ), "
    "the VV branching surface, total width, production-mode mixture, "
    "interference, or the experiment's mass-dependent acceptance and limit "
    "curve."
)


def resolve_kk_diboson_spin1_mass_gev(
    *,
    ew_mass_extra: Any = None,
    couplings: Any = None,
) -> tuple[float | None, str | None]:
    """Resolve the CR012 spin-1 diboson mass proxy in GeV.

    The explicit electroweak KK mass is the closest declared scan input for a
    W'/Z'-like diboson resonance.  If it is absent, the common quark-sector
    ``M_KK`` carried by the couplings object is used as a documented fallback
    RS scale proxy and flagged as ``NEEDS-HUMAN-PHYSICS`` by CR012.
    """

    if ew_mass_extra is not None:
        mass = float(ew_mass_extra)
        _kk_mass_tev_from_m_kk_gev(mass)
        return mass, "kk_ew_mass_gev"
    if couplings is not None:
        return float(_mass_from_source_gev(couplings)), "quark_mass_basis_couplings.M_KK"
    return None, None


def kk_diboson_spin1_prediction_from_m_kk_gev(
    m_kk_gev: float,
    *,
    resonance: str = "V_KK^(1) spin-1",
    final_state: str = "WW/WZ/ZZ",
    sigma_times_br: float | None = None,
    sigma_times_br_units: str | None = None,
) -> ColliderResonancePrediction:
    """Build the documented CR012 spin-1 diboson mass-proxy prediction."""

    return ColliderResonancePrediction(
        resonance=resonance,
        final_state=final_state,
        mass_tev=_kk_mass_tev_from_m_kk_gev(m_kk_gev),
        sigma_times_br=sigma_times_br,
        sigma_times_br_units=sigma_times_br_units,
        matching_assumption=KK_DIBOSON_SPIN1_MASS_PROXY_ASSUMPTION_V1,
        diagnostics={
            "m_kk_gev": float(m_kk_gev),
            "m_spin1_diboson_proxy_gev": float(m_kk_gev),
            "mass_proxy": "m_spin1_diboson = kk_ew_mass_gev or M_KK",
            "sigma_times_br_proxy_available": sigma_times_br is not None,
        },
    )


__all__.extend(
    [
        "KK_DIBOSON_SPIN1_MASS_PROXY_ASSUMPTION_V1",
        "resolve_kk_diboson_spin1_mass_gev",
        "kk_diboson_spin1_prediction_from_m_kk_gev",
    ]
)


TOP_PHILIC_VECTOR_FOUR_TOP_PROXY_ASSUMPTION_V1 = (
    "NEEDS-HUMAN-PHYSICS: top-philic vector four-top recast v1 uses the "
    "supplied kk_ew_mass_gev, or M_KK as a fallback top-philic vector mass "
    "proxy, and compares it to the catalogued CMS top-philic Z' mass "
    "exclusion at the Gamma/m=50% benchmark. It does not compute "
    "sigma(pp->ttbar Z')*BR(Z'->ttbar), the width-dependent signal "
    "normalization, top-philic coupling pattern, four-top acceptance, "
    "interference, or the SM four-top background likelihood."
)


def resolve_top_philic_vector_mass_gev(
    *,
    ew_mass_extra: Any = None,
    couplings: Any = None,
) -> tuple[float | None, str | None]:
    """Resolve the CR014 top-philic vector mass proxy in GeV.

    The explicit electroweak KK mass is the closest declared scan input for a
    neutral Z'-like vector.  If it is absent, the common quark-sector ``M_KK``
    carried by the mass-basis couplings object is used as the documented RS
    scale proxy and flagged as ``NEEDS-HUMAN-PHYSICS`` by CR014.
    """

    if ew_mass_extra is not None:
        mass = float(ew_mass_extra)
        _kk_mass_tev_from_m_kk_gev(mass)
        return mass, "kk_ew_mass_gev"
    if couplings is not None:
        return float(_mass_from_source_gev(couplings)), "quark_mass_basis_couplings.M_KK"
    return None, None


def top_philic_vector_four_top_prediction_from_mass_gev(
    mass_gev: float,
    *,
    resonance: str = "top-philic vector mediator Z'",
    final_state: str = "t tbar t tbar two-lepton",
    width_over_mass: float | None = None,
    sigma_times_br: float | None = None,
    sigma_times_br_units: str | None = None,
) -> ColliderResonancePrediction:
    """Build the documented CR014 top-philic vector four-top mass proxy."""

    diagnostics = {
        "m_top_philic_vector_proxy_gev": float(mass_gev),
        "mass_proxy": "m_Zprime_top_philic = kk_ew_mass_gev or M_KK",
        "width_over_mass_benchmark": (
            None if width_over_mass is None else float(width_over_mass)
        ),
        "sigma_times_br_proxy_available": sigma_times_br is not None,
    }
    return ColliderResonancePrediction(
        resonance=resonance,
        final_state=final_state,
        mass_tev=_kk_mass_tev_from_m_kk_gev(mass_gev),
        sigma_times_br=sigma_times_br,
        sigma_times_br_units=sigma_times_br_units,
        matching_assumption=TOP_PHILIC_VECTOR_FOUR_TOP_PROXY_ASSUMPTION_V1,
        diagnostics=diagnostics,
    )


__all__.extend(
    [
        "TOP_PHILIC_VECTOR_FOUR_TOP_PROXY_ASSUMPTION_V1",
        "resolve_top_philic_vector_mass_gev",
        "top_philic_vector_four_top_prediction_from_mass_gev",
    ]
)
