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
