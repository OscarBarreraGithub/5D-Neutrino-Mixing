"""Adapter for ``Z -> mu tau`` LFV constraints.

This module reuses the shared :mod:`quarkConstraints.zpole_lfv`
off-diagonal-coupling branching-fraction machinery, but supplies a mu-tau
proxy extractor. Constraint modules import this adapter only.

NEEDS-HUMAN-PHYSICS: ``ParameterPoint`` does not carry a rigorous
off-diagonal lepton ``Z mu tau`` effective coupling. The proxy path maps
caller-supplied lepton overlap spurions onto
``delta g_mutau ~ (m_Z/M_KK)^2 overlap_mutau`` and is explicitly
diagnostic-only until a full RS EW/lepton matching object exists.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Mapping

import numpy as np

from quarkConstraints.zpole import ZPoleSMInputs, default_sm_inputs
from quarkConstraints.zpole_lfv import (
    ZPOLE_LFV_MODEL_V1,
    ZPoleLFVBranchingResult,
    sm_total_width_weight as _sm_total_width_weight,
    z_lfv_branching_fraction_from_couplings as _branching_from_couplings,
    z_lfv_effective_coupling_limit as _effective_coupling_limit,
)

__all__ = [
    "ZPOLE_LFV_MUTAU_PROXY_V1",
    "ZPoleSMInputs",
    "ZPoleLFVBranchingResult",
    "ZPoleLFVMuTauProxyInput",
    "ZPoleLFVMuTauCouplingProxy",
    "zpole_lfv_mutau_proxy_input",
    "zpole_lfv_mutau_sm_total_width_weight",
    "zpole_lfv_mutau_effective_coupling_limit",
    "zpole_lfv_mutau_coupling_proxy",
    "zpole_lfv_mutau_branching_fraction_with_proxy",
]

ZPOLE_LFV_MUTAU_PROXY_V1 = (
    "NEEDS-HUMAN-PHYSICS: off-diagonal lepton Z mu-tau couplings are not "
    "available on ParameterPoint; caller-supplied lepton overlap spurions are "
    "mapped to delta g_mutau = proxy_strength * (m_Z/M_KK)^2 * overlap_mutau, "
    "standing in for missing model-dependent RS EW KK/Z/Z' and lepton "
    "neutral-current matching."
)

_DEFAULT_CHARGE_STATE_FACTOR = 2.0


@dataclass(frozen=True)
class ZPoleLFVMuTauProxyInput:
    """Explicit lepton-overlap proxy for the LFV ``Z mu tau`` coupling."""

    left_mutau_overlap: complex
    right_mutau_overlap: complex
    m_kk_gev: float
    source: str = "caller-supplied lepton neutral-current mu-tau proxy"


@dataclass(frozen=True)
class ZPoleLFVMuTauCouplingProxy:
    """Documented proxy mapping mu-tau overlap spurions to ``Z mu tau`` shifts."""

    model_label: str
    matching_assumption: str
    M_KK: float
    matching_scale: float
    left_mutau_overlap: complex
    right_mutau_overlap: complex
    scale_factor: float
    delta_g_left_mutau: complex
    delta_g_right_mutau: complex
    source: str
    diagnostics: Mapping[str, Any]


def zpole_lfv_mutau_proxy_input(
    left_mutau_overlap: complex,
    right_mutau_overlap: complex,
    m_kk_gev: float,
    *,
    source: str = "caller-supplied lepton neutral-current mu-tau proxy",
) -> ZPoleLFVMuTauProxyInput:
    """Build a shape-checked explicit proxy input for T017-like constraints."""

    return ZPoleLFVMuTauProxyInput(
        left_mutau_overlap=complex(left_mutau_overlap),
        right_mutau_overlap=complex(right_mutau_overlap),
        m_kk_gev=_positive_float("m_kk_gev", m_kk_gev),
        source=str(source),
    )


def zpole_lfv_mutau_sm_total_width_weight(
    inputs: ZPoleSMInputs | None = None,
) -> dict[str, float]:
    """Return SM total Z-width weights from the shared Z-pole convention."""

    return _sm_total_width_weight(inputs)


def zpole_lfv_mutau_effective_coupling_limit(
    br_limit: float,
    *,
    inputs: ZPoleSMInputs | None = None,
) -> float:
    """Return the limit on ``sqrt(|delta_g_L|^2 + |delta_g_R|^2)``."""

    return _effective_coupling_limit(br_limit, inputs=inputs)


def zpole_lfv_mutau_coupling_proxy(
    source: object,
    *,
    m_kk_gev: float | None = None,
    inputs: ZPoleSMInputs | None = None,
) -> ZPoleLFVMuTauCouplingProxy:
    """Return the documented lepton-overlap proxy for ``Z mu tau`` couplings."""

    p = default_sm_inputs() if inputs is None else inputs
    proxy_input = _coerce_proxy_input(source, m_kk_gev=m_kk_gev)
    resolved_m_kk = _positive_float(
        "m_kk_gev",
        proxy_input.m_kk_gev if m_kk_gev is None else m_kk_gev,
    )
    scale = float((p.m_z_gev / resolved_m_kk) ** 2)
    left = complex(p.proxy_strength * scale * proxy_input.left_mutau_overlap)
    right = complex(p.proxy_strength * scale * proxy_input.right_mutau_overlap)
    return ZPoleLFVMuTauCouplingProxy(
        model_label=ZPOLE_LFV_MODEL_V1,
        matching_assumption=ZPOLE_LFV_MUTAU_PROXY_V1,
        M_KK=resolved_m_kk,
        matching_scale=resolved_m_kk,
        left_mutau_overlap=complex(proxy_input.left_mutau_overlap),
        right_mutau_overlap=complex(proxy_input.right_mutau_overlap),
        scale_factor=scale,
        delta_g_left_mutau=left,
        delta_g_right_mutau=right,
        source=proxy_input.source,
        diagnostics={
            "m_kk_gev": float(resolved_m_kk),
            "matching_scale_gev": float(resolved_m_kk),
            "m_z_gev": float(p.m_z_gev),
            "proxy_strength": float(p.proxy_strength),
            "scale_factor": float(scale),
            "left_mutau_overlap": complex(proxy_input.left_mutau_overlap),
            "right_mutau_overlap": complex(proxy_input.right_mutau_overlap),
            "delta_g_left_mutau": left,
            "delta_g_right_mutau": right,
            "proxy_source": proxy_input.source,
            "matching_assumption": ZPOLE_LFV_MUTAU_PROXY_V1,
        },
    )


def zpole_lfv_mutau_branching_fraction_with_proxy(
    source: object,
    *,
    br_limit: float | None = None,
    m_kk_gev: float | None = None,
    inputs: ZPoleSMInputs | None = None,
) -> tuple[ZPoleLFVBranchingResult, ZPoleLFVMuTauCouplingProxy]:
    """Evaluate ``BR(Z -> mu tau)`` with the documented mu-tau proxy."""

    p = default_sm_inputs() if inputs is None else inputs
    proxy = zpole_lfv_mutau_coupling_proxy(source, m_kk_gev=m_kk_gev, inputs=p)
    result = _branching_from_couplings(
        delta_g_left=proxy.delta_g_left_mutau,
        delta_g_right=proxy.delta_g_right_mutau,
        initial_flavor="mu",
        final_flavor="tau",
        br_limit=br_limit,
        inputs=p,
        charge_state_factor=_DEFAULT_CHARGE_STATE_FACTOR,
    )
    return result, proxy


def _coerce_proxy_input(
    value: Any,
    *,
    m_kk_gev: float | None,
) -> ZPoleLFVMuTauProxyInput:
    if isinstance(value, ZPoleLFVMuTauProxyInput):
        return value
    if isinstance(value, Mapping):
        return _proxy_from_mapping(value, m_kk_gev=m_kk_gev)
    return _proxy_from_object(value, m_kk_gev=m_kk_gev)


def _proxy_from_mapping(
    mapping: Mapping[str, Any],
    *,
    m_kk_gev: float | None,
) -> ZPoleLFVMuTauProxyInput:
    source = str(mapping.get("source", "mapping lepton neutral-current mu-tau proxy"))
    left, right = _overlaps_from_mapping(mapping)
    return ZPoleLFVMuTauProxyInput(
        left_mutau_overlap=left,
        right_mutau_overlap=right,
        m_kk_gev=_resolve_m_kk_from_mapping(mapping, m_kk_gev),
        source=source,
    )


def _proxy_from_object(value: Any, *, m_kk_gev: float | None) -> ZPoleLFVMuTauProxyInput:
    left = _first_present_attr(value, ("left_mutau_overlap", "left_mutau"))
    right = _first_present_attr(value, ("right_mutau_overlap", "right_mutau"))
    if left is None and right is None:
        left_matrix = _first_present_attr(
            value,
            (
                "left_charged_lepton_overlap",
                "left_lepton_overlap",
                "left_overlap",
            ),
        )
        right_matrix = _first_present_attr(
            value,
            (
                "right_charged_lepton_overlap",
                "right_lepton_overlap",
                "right_overlap",
            ),
        )
        if left_matrix is None and right_matrix is None:
            raise TypeError(
                "lepton LFV mu-tau proxy must provide left/right mu-tau overlap "
                "entries or charged-lepton overlap matrices"
            )
        left = 0.0j if left_matrix is None else _matrix_offdiag(left_matrix, "left", 1, 2)
        right = 0.0j if right_matrix is None else _matrix_offdiag(right_matrix, "right", 1, 2)
    return ZPoleLFVMuTauProxyInput(
        left_mutau_overlap=complex(0.0j if left is None else left),
        right_mutau_overlap=complex(0.0j if right is None else right),
        m_kk_gev=_resolve_m_kk_from_object(value, m_kk_gev),
        source=str(getattr(value, "source", "object lepton neutral-current mu-tau proxy")),
    )


def _overlaps_from_mapping(mapping: Mapping[str, Any]) -> tuple[complex, complex]:
    left = _first_present_key(mapping, ("left_mutau_overlap", "left_mutau"))
    right = _first_present_key(mapping, ("right_mutau_overlap", "right_mutau"))
    if left is not None or right is not None:
        return complex(0.0j if left is None else left), complex(0.0j if right is None else right)

    left_matrix = _first_present_key(
        mapping,
        ("left_charged_lepton_overlap", "left_lepton_overlap", "left_overlap"),
    )
    right_matrix = _first_present_key(
        mapping,
        ("right_charged_lepton_overlap", "right_lepton_overlap", "right_overlap"),
    )
    if left_matrix is None and right_matrix is None:
        raise KeyError(
            "mapping must provide left/right mu-tau overlap entries or "
            "charged-lepton overlap matrices"
        )
    left_value = 0.0j if left_matrix is None else _matrix_offdiag(left_matrix, "left", 1, 2)
    right_value = 0.0j if right_matrix is None else _matrix_offdiag(right_matrix, "right", 1, 2)
    return complex(left_value), complex(right_value)


def _resolve_m_kk_from_mapping(
    mapping: Mapping[str, Any],
    override: float | None,
) -> float:
    if override is not None:
        return _positive_float("m_kk_gev", override)
    value = _first_present_key(mapping, ("m_kk_gev", "M_KK_gev", "M_KK"))
    if value is None:
        raise KeyError("mapping must provide m_kk_gev, M_KK_gev, or M_KK")
    return _positive_float("m_kk_gev", value)


def _resolve_m_kk_from_object(value: Any, override: float | None) -> float:
    if override is not None:
        return _positive_float("m_kk_gev", override)
    for name in ("m_kk_gev", "M_KK_gev", "M_KK"):
        if hasattr(value, name):
            return _positive_float("m_kk_gev", getattr(value, name))
    raise AttributeError("lepton LFV mu-tau proxy object must provide m_kk_gev or M_KK")


def _first_present_key(mapping: Mapping[str, Any], keys: tuple[str, ...]) -> Any:
    for key in keys:
        if key in mapping:
            return mapping[key]
    return None


def _first_present_attr(value: Any, names: tuple[str, ...]) -> Any:
    for name in names:
        if hasattr(value, name):
            return getattr(value, name)
    return None


def _matrix_offdiag(value: Any, name: str, row: int, col: int) -> complex:
    matrix = np.asarray(value, dtype=np.complex128)
    if matrix.shape != (3, 3):
        raise ValueError(f"{name} charged-lepton overlap matrix must have shape (3, 3)")
    if not np.all(np.isfinite(matrix.real)) or not np.all(np.isfinite(matrix.imag)):
        raise ValueError(f"{name} charged-lepton overlap matrix entries must be finite")
    return complex(matrix[row, col])


def _positive_float(name: str, value: object) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0.0:
        raise ValueError(f"{name} must be positive and finite")
    return number
