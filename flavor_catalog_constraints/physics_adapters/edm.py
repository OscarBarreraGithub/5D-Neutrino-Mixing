"""Adapter over :mod:`quarkConstraints.edm`.

This is the catalog boundary for charged-lepton electric dipole moments.
Constraint modules import this adapter only; the observable conversion lives
in ``quarkConstraints.edm``.

NEEDS-HUMAN-PHYSICS
-------------------
The available ``ParameterPoint`` does not contain the complex lepton-sector RS
couplings needed for a one-loop EDM matching calculation.  Until that exists,
this adapter accepts an explicit low-energy CP-odd dipole proxy coefficient
for the requested lepton and flags the result.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Mapping

from quarkConstraints.edm import (
    EDM_COEFFICIENT_CONVENTION_V1,
    EDM_RS_MATCHING_GAP_V1,
    HBARC_GEV_CM,
    ChargedLeptonEDMResult,
)
from quarkConstraints.edm import (
    evaluate_charged_lepton_edm as _evaluate_charged_lepton_edm,
)
from quarkConstraints.edm import (
    evaluate_charged_lepton_edm_from_chiral_dipole as _evaluate_from_chiral_dipole,
)

__all__ = [
    "EDM_COEFFICIENT_CONVENTION_V1",
    "EDM_RS_MATCHING_GAP_V1",
    "HBARC_GEV_CM",
    "LEPTON_EDM_PROXY_ASSUMPTION_V1",
    "LeptonEDMProxyInput",
    "ChargedLeptonEDMResult",
    "lepton_edm_proxy_input",
    "charged_lepton_edm_from_lepton_input",
    "electron_edm_from_lepton_input",
]

LEPTON_EDM_PROXY_ASSUMPTION_V1 = (
    "NEEDS-HUMAN-PHYSICS: charged-lepton EDMs need loop-level RS matching "
    "from complex lepton/KK-fermion/Higgs/electroweak couplings. "
    "ParameterPoint does not currently expose those couplings, so this "
    "adapter uses an explicit low-energy CP-odd dipole proxy coefficient."
)


@dataclass(frozen=True)
class LeptonEDMProxyInput:
    """Explicit low-energy charged-lepton EDM proxy.

    ``cp_odd_dipole_coefficient_gev_inv`` is the coefficient ``c_CPodd`` in
    ``d_l/e = c_CPodd``.  ``m_kk_gev`` is diagnostic only for this direct
    low-energy proxy; no mass scaling is applied here.
    """

    lepton: str
    cp_odd_dipole_coefficient_gev_inv: float
    m_kk_gev: float | None = None
    source: str = "caller-supplied charged-lepton EDM proxy"


def _finite_float(value: Any, *, name: str) -> float:
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"{name} must be finite")
    return number


def _positive_finite_or_none(value: Any, *, name: str) -> float | None:
    if value is None:
        return None
    number = _finite_float(value, name=name)
    if number <= 0.0:
        raise ValueError(f"{name} must be positive")
    return number


def lepton_edm_proxy_input(
    cp_odd_dipole_coefficient_gev_inv: float,
    *,
    lepton: str = "e",
    m_kk_gev: float | None = None,
    source: str = "caller-supplied charged-lepton EDM proxy",
) -> LeptonEDMProxyInput:
    """Build a shape-checked proxy for ``charged_lepton_edm_from_lepton_input``."""

    return LeptonEDMProxyInput(
        lepton=_canonical_lepton(lepton),
        cp_odd_dipole_coefficient_gev_inv=_finite_float(
            cp_odd_dipole_coefficient_gev_inv,
            name="cp_odd_dipole_coefficient_gev_inv",
        ),
        m_kk_gev=_positive_finite_or_none(m_kk_gev, name="m_kk_gev"),
        source=str(source),
    )


def _mapping_value(mapping: Mapping[str, Any], keys: tuple[str, ...]) -> Any:
    for key in keys:
        if key in mapping:
            return mapping[key]
    raise KeyError(f"none of keys {keys!r} present")


def _lepton_prefixes(lepton: str) -> tuple[str, ...]:
    lepton = _canonical_lepton(lepton)
    aliases = {
        "e": ("electron", "e"),
        "mu": ("muon", "mu"),
        "tau": ("tau",),
    }
    return aliases.get(lepton, (lepton,))


def _canonical_lepton(lepton: Any) -> str:
    name = str(lepton).strip().lower()
    aliases = {
        "electron": "e",
        "e": "e",
        "muon": "mu",
        "mu": "mu",
        "tau": "tau",
    }
    return aliases.get(name, name)


def _direct_coefficient_keys(lepton: str) -> tuple[str, ...]:
    keys = [
        "cp_odd_dipole_coefficient_gev_inv",
        "cp_odd_dipole_gev_inv",
        "edm_cp_odd_dipole_coefficient_gev_inv",
        "d_over_e_gev_inv",
    ]
    for prefix in _lepton_prefixes(lepton):
        keys.extend(
            [
                f"{prefix}_cp_odd_dipole_coefficient_gev_inv",
                f"{prefix}_cp_odd_dipole_gev_inv",
                f"{prefix}_edm_cp_odd_dipole_coefficient_gev_inv",
                f"d_{prefix}_over_e_gev_inv",
            ]
        )
    return tuple(keys)


def _chiral_coefficient_keys(lepton: str) -> tuple[str, ...]:
    keys = [
        "chiral_dipole_coefficient_gev_inv",
        "edm_chiral_dipole_coefficient_gev_inv",
    ]
    for prefix in _lepton_prefixes(lepton):
        keys.extend(
            [
                f"{prefix}_chiral_dipole_coefficient_gev_inv",
                f"{prefix}_edm_chiral_dipole_coefficient_gev_inv",
            ]
        )
    return tuple(keys)


def _m_kk_from_mapping(mapping: Mapping[str, Any]) -> float | None:
    for key in ("m_kk_gev", "M_KK", "M_KK_gev", "kk_ew_mass_gev"):
        if key in mapping:
            return _positive_finite_or_none(mapping[key], name=key)
    return None


def _proxy_from_mapping(mapping: Mapping[str, Any], *, lepton: str) -> LeptonEDMProxyInput:
    if "edm_proxy" in mapping:
        return _coerce_proxy(mapping["edm_proxy"], lepton=lepton)
    direct = _mapping_value(mapping, _direct_coefficient_keys(lepton))
    return lepton_edm_proxy_input(
        direct,
        lepton=str(mapping.get("lepton", lepton)),
        m_kk_gev=_m_kk_from_mapping(mapping),
        source=str(mapping.get("source", "mapping charged-lepton EDM proxy")),
    )


def _coerce_proxy(value: Any, *, lepton: str) -> LeptonEDMProxyInput:
    lepton = _canonical_lepton(lepton)
    if isinstance(value, LeptonEDMProxyInput):
        proxy = value
    elif isinstance(value, Mapping):
        proxy = _proxy_from_mapping(value, lepton=lepton)
    else:
        raise TypeError(
            "lepton EDM input must be a LeptonEDMProxyInput or a mapping with "
            "a CP-odd dipole coefficient"
        )

    if _canonical_lepton(proxy.lepton) != lepton:
        raise ValueError(f"EDM proxy lepton {proxy.lepton!r} does not match {lepton!r}")
    return proxy


def _mapping_chiral_coefficient(
    value: Any,
    *,
    lepton: str,
) -> tuple[complex, float | None, str] | None:
    if not isinstance(value, Mapping):
        return None
    if "edm_proxy" in value:
        return _mapping_chiral_coefficient(value["edm_proxy"], lepton=lepton)
    for key in _chiral_coefficient_keys(lepton):
        if key in value:
            return (
                complex(value[key]),
                _m_kk_from_mapping(value),
                str(value.get("source", "mapping complex chiral EDM proxy")),
            )
    return None


def charged_lepton_edm_from_lepton_input(
    lepton_input: Any,
    *,
    lepton: str,
    experimental_limit_e_cm: float,
    sm_edm_e_cm: float = 0.0,
    m_kk_gev: float | None = None,
) -> ChargedLeptonEDMResult:
    """Evaluate a charged-lepton EDM from an explicit lepton proxy input."""

    lepton_name = _canonical_lepton(lepton)
    chiral = _mapping_chiral_coefficient(lepton_input, lepton=lepton_name)
    if chiral is not None:
        coefficient, proxy_m_kk, source = chiral
        effective_m_kk = _positive_finite_or_none(
            proxy_m_kk if m_kk_gev is None else m_kk_gev,
            name="m_kk_gev",
        )
        return _evaluate_from_chiral_dipole(
            coefficient,
            experimental_limit_e_cm=experimental_limit_e_cm,
            lepton=lepton_name,
            sm_edm_e_cm=sm_edm_e_cm,
            coefficient_source=source,
            diagnostics={
                "needs_human_physics": LEPTON_EDM_PROXY_ASSUMPTION_V1,
                "used_proxy": True,
                "input_kind": "complex chiral EDM mapping",
                "proxy_source": source,
                "m_kk_gev": effective_m_kk,
                "m_kk_diagnostic_only": True,
            },
        )

    proxy = _coerce_proxy(lepton_input, lepton=lepton_name)
    effective_m_kk = _positive_finite_or_none(
        proxy.m_kk_gev if m_kk_gev is None else m_kk_gev,
        name="m_kk_gev",
    )
    return _evaluate_charged_lepton_edm(
        proxy.cp_odd_dipole_coefficient_gev_inv,
        experimental_limit_e_cm=experimental_limit_e_cm,
        lepton=lepton_name,
        sm_edm_e_cm=sm_edm_e_cm,
        coefficient_source=proxy.source,
        diagnostics={
            "needs_human_physics": LEPTON_EDM_PROXY_ASSUMPTION_V1,
            "used_proxy": True,
            "input_kind": "LeptonEDMProxyInput",
            "proxy_source": proxy.source,
            "m_kk_gev": effective_m_kk,
            "m_kk_diagnostic_only": True,
            "cp_odd_projection": "direct_cp_odd_coefficient",
        },
    )


def electron_edm_from_lepton_input(
    lepton_input: Any,
    *,
    experimental_limit_e_cm: float,
    sm_edm_e_cm: float = 0.0,
    m_kk_gev: float | None = None,
) -> ChargedLeptonEDMResult:
    """Evaluate ``|d_e|`` from an explicit electron EDM proxy input."""

    return charged_lepton_edm_from_lepton_input(
        lepton_input,
        lepton="e",
        experimental_limit_e_cm=experimental_limit_e_cm,
        sm_edm_e_cm=sm_edm_e_cm,
        m_kk_gev=m_kk_gev,
    )
