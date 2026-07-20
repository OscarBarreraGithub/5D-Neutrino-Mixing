"""Adapter for muonium-antimuonium conversion.

This module is the catalog boundary for L006.  There is no repo-local RS
Delta L=2 four-lepton matching core, so the v1 implementation is a small
low-energy proxy.

NEEDS-HUMAN-PHYSICS: a production RS prediction needs a Lorentz-basis
Delta L_mu = -Delta L_e = 2 four-lepton Wilson coefficient, spin-state
normalization, and magnetic-field treatment.  The v1 proxy accepts either an
explicit PDG-convention effective coupling ratio ``G_C/G_F`` or an explicit
conversion probability and calibrates the coupling-to-probability map with the
catalogued MACS/PDG anchors supplied by the constraint.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Mapping

__all__ = [
    "MUONIUM_CONVERSION_MODEL_V1",
    "MUONIUM_CONVERSION_OPERATOR_CONVENTION",
    "MUONIUM_CONVERSION_PROXY_V1",
    "MuoniumConversionProxyInput",
    "MuoniumConversionResult",
    "muonium_conversion_proxy_input",
    "muonium_conversion_from_lepton_input",
]

MUONIUM_CONVERSION_MODEL_V1 = "muonium_antimuonium_probability_proxy_v1"
MUONIUM_CONVERSION_OPERATOR_CONVENTION = (
    "Input g_mmbar_over_gf is interpreted as the PDG G_C/G_F ratio for the "
    "(V-A)x(V-A) effective Lagrangian.  The conversion probability is "
    "calibrated to the catalog anchors as "
    "P_proxy = P_limit * |g_mmbar_over_gf / g_limit_over_gf|^2.  A direct "
    "conversion_probability input is interpreted as already being in the "
    "MACS/PDG 0.1 T probability convention."
)
MUONIUM_CONVERSION_PROXY_V1 = (
    "NEEDS-HUMAN-PHYSICS: full RS Delta L=2 four-lepton matching is not "
    "available on ParameterPoint.  Caller-supplied low-energy proxy inputs "
    "are compared with the MACS/PDG muonium-antimuonium probability bound; "
    "operator-dependent spin and magnetic-field factors are not predicted."
)


@dataclass(frozen=True)
class MuoniumConversionProxyInput:
    """Explicit low-energy proxy input for muonium conversion.

    Provide exactly one of ``g_mmbar_over_gf`` or ``conversion_probability``.
    ``g_mmbar_over_gf`` is the dimensionless PDG-convention effective
    coupling ratio.  ``conversion_probability`` is a direct pure-NP
    probability in the MACS/PDG 0.1 T convention.
    """

    g_mmbar_over_gf: complex | None = None
    conversion_probability: float | None = None
    source: str = "caller-supplied muonium conversion proxy"


@dataclass(frozen=True)
class MuoniumConversionResult:
    """Pure-NP probability prediction for ``M -> Mbar`` conversion."""

    model_label: str
    conversion_probability: float
    sm_conversion_probability: float
    np_conversion_probability: float
    probability_limit: float
    effective_coupling_over_gf: complex | None
    effective_coupling_limit_over_gf: float
    probability_per_coupling_ratio_squared: float
    ratio_to_limit: float
    passes: bool
    input_kind: str
    used_proxy: bool
    source: str
    diagnostics: Mapping[str, Any]


def muonium_conversion_proxy_input(
    *,
    g_mmbar_over_gf: complex | None = None,
    conversion_probability: float | None = None,
    source: str = "caller-supplied muonium conversion proxy",
) -> MuoniumConversionProxyInput:
    """Build a shape-checked proxy input for ``muonium_conversion_from_lepton_input``."""

    has_coupling = g_mmbar_over_gf is not None
    has_probability = conversion_probability is not None
    if has_coupling == has_probability:
        raise ValueError(
            "provide exactly one of g_mmbar_over_gf or conversion_probability"
        )
    return MuoniumConversionProxyInput(
        g_mmbar_over_gf=(
            None
            if g_mmbar_over_gf is None
            else _finite_complex(g_mmbar_over_gf, "g_mmbar_over_gf")
        ),
        conversion_probability=(
            None
            if conversion_probability is None
            else _bounded_probability(
                conversion_probability,
                "conversion_probability",
            )
        ),
        source=str(source),
    )


def muonium_conversion_from_lepton_input(
    lepton_input: Any,
    *,
    probability_limit: float,
    coupling_limit_over_gf: float,
) -> MuoniumConversionResult:
    """Evaluate the v1 muonium-conversion proxy from a mapping or dataclass."""

    p_limit = _bounded_probability(probability_limit, "probability_limit")
    g_limit = _positive_float(coupling_limit_over_gf, "coupling_limit_over_gf")
    proxy = _coerce_proxy(lepton_input)
    calibration = float(p_limit / (g_limit * g_limit))

    if proxy.g_mmbar_over_gf is not None:
        g_value = complex(proxy.g_mmbar_over_gf)
        probability = float(calibration * abs(g_value) ** 2)
        input_kind = "g_mmbar_over_gf"
        direct_probability_used = False
    elif proxy.conversion_probability is not None:
        g_value = None
        probability = float(proxy.conversion_probability)
        input_kind = "conversion_probability"
        direct_probability_used = True
    else:  # pragma: no cover - guarded by constructors/coercion.
        raise ValueError("muonium conversion proxy has no usable input")

    if not math.isfinite(probability) or probability < 0.0:
        raise ValueError("conversion probability must be a nonnegative finite number")

    ratio = float(probability / p_limit)
    diagnostics = {
        "model_label": MUONIUM_CONVERSION_MODEL_V1,
        "operator_convention": MUONIUM_CONVERSION_OPERATOR_CONVENTION,
        "matching_assumption": MUONIUM_CONVERSION_PROXY_V1,
        "needs_human_physics": MUONIUM_CONVERSION_PROXY_V1,
        "sm_conversion_probability": 0.0,
        "sm_lfv_policy": (
            "Muonium-antimuonium conversion is negligible in the SM for "
            "catalog purposes; L006 is applied as a pure-NP probability bound."
        ),
        "probability_limit": float(p_limit),
        "effective_coupling_limit_over_gf": float(g_limit),
        "probability_per_coupling_ratio_squared": float(calibration),
        "effective_coupling_over_gf": g_value,
        "effective_coupling_abs_over_gf": (
            None if g_value is None else float(abs(g_value))
        ),
        "direct_probability_proxy_used": direct_probability_used,
        "input_kind": input_kind,
        "used_proxy": True,
        "proxy_source": proxy.source,
    }
    return MuoniumConversionResult(
        model_label=MUONIUM_CONVERSION_MODEL_V1,
        conversion_probability=float(probability),
        sm_conversion_probability=0.0,
        np_conversion_probability=float(probability),
        probability_limit=float(p_limit),
        effective_coupling_over_gf=g_value,
        effective_coupling_limit_over_gf=float(g_limit),
        probability_per_coupling_ratio_squared=float(calibration),
        ratio_to_limit=ratio,
        passes=bool(ratio <= 1.0),
        input_kind=input_kind,
        used_proxy=True,
        source=proxy.source,
        diagnostics=diagnostics,
    )


def _coerce_proxy(value: Any) -> MuoniumConversionProxyInput:
    if isinstance(value, MuoniumConversionProxyInput):
        return value
    if isinstance(value, Mapping):
        if "muonium_conversion" in value:
            return _coerce_proxy(value["muonium_conversion"])
        return _proxy_from_mapping(value)
    raise TypeError(
        "lepton input must be a MuoniumConversionProxyInput or a mapping with "
        "g_mmbar_over_gf or conversion_probability"
    )


_COUPLING_KEYS = (
    "g_mmbar_over_gf",
    "G_MMbar_over_GF",
    "G_C_over_G_F",
    "g_c_over_gf",
    "muonium_g_over_gf",
)
_PROBABILITY_KEYS = (
    "conversion_probability",
    "muonium_conversion_probability",
    "p_mmbar",
    "P_MMbar",
    "P_M_Mbar",
)


def _proxy_from_mapping(mapping: Mapping[str, Any]) -> MuoniumConversionProxyInput:
    coupling_key = next((key for key in _COUPLING_KEYS if key in mapping), None)
    probability_key = next((key for key in _PROBABILITY_KEYS if key in mapping), None)
    if coupling_key is None and probability_key is None:
        raise KeyError(
            "mapping must contain one of "
            f"{_COUPLING_KEYS!r} or {_PROBABILITY_KEYS!r}"
        )
    if coupling_key is not None and probability_key is not None:
        raise ValueError(
            "mapping must not provide both effective coupling and direct "
            "probability proxies"
        )
    return muonium_conversion_proxy_input(
        g_mmbar_over_gf=(
            None if coupling_key is None else mapping[coupling_key]
        ),
        conversion_probability=(
            None if probability_key is None else mapping[probability_key]
        ),
        source=str(mapping.get("source", "mapping muonium conversion proxy")),
    )


def _finite_complex(value: Any, name: str) -> complex:
    number = complex(value)
    if not math.isfinite(number.real) or not math.isfinite(number.imag):
        raise ValueError(f"{name} must be finite")
    return number


def _positive_float(value: Any, name: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0.0:
        raise ValueError(f"{name} must be a positive finite number")
    return number


def _bounded_probability(value: Any, name: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number < 0.0 or number > 1.0:
        raise ValueError(f"{name} must be a finite probability in [0, 1]")
    return number
