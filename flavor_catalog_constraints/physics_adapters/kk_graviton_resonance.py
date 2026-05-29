"""Adapter for spin-2 KK-graviton direct-search constraints.

This module reuses the existing collider-resonance mass-limit machinery.  It
adds only the CR007-specific mass proxy and keeps the ``quarkConstraints``
core behind the established ``physics_adapters.collider_resonance`` boundary.

NEEDS-HUMAN-PHYSICS: ``ParameterPoint`` does not carry a complete physical
spin-2 graviton spectrum, curvature coupling, total width, branching
fractions, or production cross sections.  The v1 recast therefore maps the
available common KK scale to the geometric IR scale and assumes
``m_GKK^(1) = x_G1 * Lambda_IR`` with the first spin-2 root below.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any

from flavor_catalog_constraints.physics_adapters.collider_resonance import (
    MASS_LOWER_BOUND,
    SIGMA_TIMES_BR_UPPER_LIMIT,
    ColliderResonanceComparison,
    ColliderResonanceLimit,
    ColliderResonancePrediction,
    evaluate_collider_resonance_limit,
    kk_gluon_mass_tev_from_m_kk_gev,
    resolve_kk_ew_mass_gev,
    resolve_kk_gluon_mass_gev,
)
from quarkConstraints.scales import (
    GAUGE_KK_ROOT_NN,
    SPIN2_GRAVITON_KK_ROOT,
    spin2_graviton_mass_from_lambda_ir,
)

KK_GRAVITON_SPIN2_MASS_PROXY_ASSUMPTION_V1 = (
    "NEEDS-HUMAN-PHYSICS: spin-2 KK-graviton recast v1 assumes "
    "m_GKK^(1) = x_G1 * Lambda_IR with x_G1=3.8317059702075125. The "
    "available kk_ew_mass_gev/kk_gluon_mass_gev inputs are first-gauge-mode "
    "masses mapped to Lambda_IR with the repo gauge root, while "
    "quark-sector M_KK is treated as the Lambda_IR bookkeeping scale unless "
    "a positive xi_KK is carried by the couplings object. It does not compute "
    "sigma(pp->G_KK)*BR(G_KK->WW/ZZ/gamma gamma/ll/ttbar), k/Mbar_Pl, "
    "total width, branching fractions, acceptance, interference, or the "
    "experiment's mass-dependent limit curve."
)

KK_GRAVITON_SPIN2_SPECTRUM_ASSUMPTION_V1 = (
    "m_GKK^(1) = x_G1 * Lambda_IR with x_G1 equal to the first spin-2 "
    "graviton root; gauge-mode mass extras are divided by the repo first "
    "gauge root before applying x_G1."
)


@dataclass(frozen=True)
class KKGravitonMassMapping:
    """Resolved spin-2 graviton mass and its spectrum bookkeeping."""

    lambda_ir_gev: float
    m_graviton_gev: float
    mass_source: str
    lambda_ir_source: str
    input_mass_gev: float
    input_mass_role: str
    graviton_root: float = SPIN2_GRAVITON_KK_ROOT
    gauge_root: float = GAUGE_KK_ROOT_NN
    spectrum_assumption: str = KK_GRAVITON_SPIN2_SPECTRUM_ASSUMPTION_V1


def _positive_finite(value: float, name: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0.0:
        raise ValueError(f"{name} must be positive and finite")
    return number


def _lambda_ir_from_gauge_mass(mass_gev: float, source: str) -> KKGravitonMassMapping:
    mass = _positive_finite(mass_gev, source)
    lambda_ir = mass / GAUGE_KK_ROOT_NN
    graviton_mass = spin2_graviton_mass_from_lambda_ir(lambda_ir)
    return KKGravitonMassMapping(
        lambda_ir_gev=float(lambda_ir),
        m_graviton_gev=float(graviton_mass),
        mass_source=source,
        lambda_ir_source=f"{source} / GAUGE_KK_ROOT_NN",
        input_mass_gev=float(mass),
        input_mass_role="first gauge KK mass",
    )


def _lambda_ir_from_couplings(m_kk_gev: float, couplings: Any) -> KKGravitonMassMapping:
    m_kk = _positive_finite(m_kk_gev, "quark_mass_basis_couplings.M_KK")
    raw_xi = getattr(couplings, "xi_KK", None)
    xi_kk = (
        1.0
        if raw_xi is None
        else _positive_finite(raw_xi, "quark_mass_basis_couplings.xi_KK")
    )
    lambda_ir = m_kk / xi_kk
    graviton_mass = spin2_graviton_mass_from_lambda_ir(lambda_ir)
    lambda_ir_source = (
        "quark_mass_basis_couplings.M_KK"
        if raw_xi is None
        else "quark_mass_basis_couplings.M_KK / quark_mass_basis_couplings.xi_KK"
    )
    return KKGravitonMassMapping(
        lambda_ir_gev=float(lambda_ir),
        m_graviton_gev=float(graviton_mass),
        mass_source="quark_mass_basis_couplings.M_KK",
        lambda_ir_source=lambda_ir_source,
        input_mass_gev=float(m_kk),
        input_mass_role=(
            "repo Lambda_IR bookkeeping scale"
            if raw_xi is None
            else "xi_KK-scaled KK scale"
        ),
    )


def resolve_kk_graviton_mass_mapping(
    *,
    ew_mass_extra: Any = None,
    gluon_mass_extra: Any = None,
    couplings: Any = None,
) -> KKGravitonMassMapping | None:
    """Resolve the first spin-2 KK-graviton mass from ``Lambda_IR``.

    The current point contract has no dedicated ``kk_graviton_mass_gev`` key.
    We infer ``Lambda_IR`` from declared extras: explicit first gauge-mode
    masses are divided by ``GAUGE_KK_ROOT_NN``; quark-sector ``M_KK`` is the
    repo bookkeeping ``Lambda_IR`` scale unless its couplings object carries a
    positive ``xi_KK``.  Couplings take precedence over ``kk_gluon_mass_gev``
    because the point builder mirrors ``M_KK`` into that extra.
    """

    if ew_mass_extra is not None:
        mass, _ = resolve_kk_ew_mass_gev(mass_extra=ew_mass_extra)
        return _lambda_ir_from_gauge_mass(mass, "kk_ew_mass_gev")
    if couplings is not None:
        mass, _ = resolve_kk_ew_mass_gev(couplings=couplings)
        return _lambda_ir_from_couplings(mass, couplings)
    if gluon_mass_extra is not None:
        mass, _ = resolve_kk_gluon_mass_gev(mass_extra=gluon_mass_extra)
        return _lambda_ir_from_gauge_mass(mass, "kk_gluon_mass_gev")
    return None


def resolve_kk_graviton_mass_gev(
    *,
    ew_mass_extra: Any = None,
    gluon_mass_extra: Any = None,
    couplings: Any = None,
) -> tuple[float | None, str | None]:
    """Resolve the first spin-2 graviton mass in GeV from declared extras."""

    mapping = resolve_kk_graviton_mass_mapping(
        ew_mass_extra=ew_mass_extra,
        gluon_mass_extra=gluon_mass_extra,
        couplings=couplings,
    )
    if mapping is None:
        return None, None
    return mapping.m_graviton_gev, mapping.mass_source


def kk_graviton_spin2_prediction_from_lambda_ir_gev(
    lambda_ir_gev: float,
    *,
    resonance: str = "G_KK^(1)",
    final_state: str = "WW, ZZ",
    sigma_times_br: float | None = None,
    sigma_times_br_units: str | None = None,
) -> ColliderResonancePrediction:
    """Build the documented v1 spin-2 graviton prediction from ``Lambda_IR``."""

    lambda_ir = _positive_finite(lambda_ir_gev, "Lambda_IR")
    m_graviton_gev = spin2_graviton_mass_from_lambda_ir(lambda_ir)
    mass_tev = kk_gluon_mass_tev_from_m_kk_gev(m_graviton_gev)
    return ColliderResonancePrediction(
        resonance=resonance,
        final_state=final_state,
        mass_tev=mass_tev,
        sigma_times_br=sigma_times_br,
        sigma_times_br_units=sigma_times_br_units,
        matching_assumption=KK_GRAVITON_SPIN2_MASS_PROXY_ASSUMPTION_V1,
        diagnostics={
            "lambda_ir_gev": float(lambda_ir),
            "graviton_spin2_root": float(SPIN2_GRAVITON_KK_ROOT),
            "gauge_kk_root_nn": float(GAUGE_KK_ROOT_NN),
            "m_graviton_gev": float(m_graviton_gev),
            "m_graviton_proxy_gev": float(m_graviton_gev),
            "m_graviton_proxy_tev": float(mass_tev),
            "mass_proxy": "m_GKK = SPIN2_GRAVITON_KK_ROOT * Lambda_IR",
            "spectrum_assumption": KK_GRAVITON_SPIN2_SPECTRUM_ASSUMPTION_V1,
            "spectrum_mapping_status": "NEEDS-HUMAN-PHYSICS",
            "sigma_times_br_proxy_available": sigma_times_br is not None,
        },
    )


def kk_graviton_spin2_prediction_from_m_kk_gev(
    m_kk_gev: float,
    *,
    resonance: str = "G_KK^(1)",
    final_state: str = "WW, ZZ",
    sigma_times_br: float | None = None,
    sigma_times_br_units: str | None = None,
) -> ColliderResonancePrediction:
    """Backward-compatible wrapper treating ``m_kk_gev`` as ``Lambda_IR``."""

    return kk_graviton_spin2_prediction_from_lambda_ir_gev(
        m_kk_gev,
        resonance=resonance,
        final_state=final_state,
        sigma_times_br=sigma_times_br,
        sigma_times_br_units=sigma_times_br_units,
    )


__all__ = [
    "MASS_LOWER_BOUND",
    "SIGMA_TIMES_BR_UPPER_LIMIT",
    "ColliderResonanceComparison",
    "ColliderResonanceLimit",
    "ColliderResonancePrediction",
    "GAUGE_KK_ROOT_NN",
    "KKGravitonMassMapping",
    "KK_GRAVITON_SPIN2_MASS_PROXY_ASSUMPTION_V1",
    "KK_GRAVITON_SPIN2_SPECTRUM_ASSUMPTION_V1",
    "SPIN2_GRAVITON_KK_ROOT",
    "evaluate_collider_resonance_limit",
    "kk_graviton_spin2_prediction_from_lambda_ir_gev",
    "kk_graviton_spin2_prediction_from_m_kk_gev",
    "resolve_kk_graviton_mass_gev",
    "resolve_kk_graviton_mass_mapping",
]
