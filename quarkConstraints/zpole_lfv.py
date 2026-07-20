"""LFV Z-pole branching fractions from effective off-diagonal couplings.

This module extends :mod:`quarkConstraints.zpole` without changing the
diagonal T010/T012 machinery.  It keeps the same effective-coupling
normalization,

    L_Z = g_Z Z_mu lbar_i gamma^mu (g_L P_L + g_R P_R) l_j,

and converts off-diagonal ``Z e mu`` chiral couplings into the charge-summed
branching fraction ``BR(Z -> e+- mu-+)``.  The total width is represented by
the same reusable Z-pole partial-width weights used by the diagonal
pseudo-observable code; for LFV rates the SM contribution is zero.

NEEDS-HUMAN-PHYSICS: the proxy helper maps caller-supplied lepton overlap
non-universality onto ``delta g_emu ~ (m_Z/M_KK)^2 overlap_emu``.  A rigorous
RS calculation needs the EW KK/Z/Z' spectrum, lepton mass-basis neutral-current
couplings, brane kinetic terms, and Z-mixing inputs that are not present on the
current ``ParameterPoint``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Mapping

import numpy as np

from .zpole import (
    ZPOLE_MODEL_V1,
    ZPoleSMInputs,
    default_sm_inputs,
    partial_width_weight,
    sm_couplings,
)

ZPOLE_LFV_MODEL_V1 = "zpole_lfv_effective_offdiagonal_couplings_v1"
ZPOLE_LFV_PROXY_V1 = (
    "NEEDS-HUMAN-PHYSICS: off-diagonal lepton Z couplings are not available "
    "on ParameterPoint; caller-supplied lepton overlap spurions are mapped to "
    "delta g_emu = proxy_strength * (m_Z/M_KK)^2 * overlap_emu, standing in "
    "for missing model-dependent RS EW KK/Z/Z' and lepton neutral-current "
    "matching."
)

_QUARK_FLAVORS: tuple[str, ...] = ("u", "d", "s", "c", "b")
_CHARGED_LEPTONS: tuple[str, ...] = ("e", "mu", "tau")
_NEUTRINO_FLAVORS: tuple[str, ...] = ("nu_e", "nu_mu", "nu_tau")
_DEFAULT_CHARGE_STATE_FACTOR = 2.0


@dataclass(frozen=True)
class ZPoleLFVProxyInput:
    """Explicit lepton-overlap proxy for LFV Z coupling shifts."""

    left_emu_overlap: complex
    right_emu_overlap: complex
    m_kk_gev: float
    source: str = "caller-supplied lepton neutral-current proxy"


@dataclass(frozen=True)
class ZPoleLFVCouplingProxy:
    """Documented proxy mapping lepton overlap spurions to ``Z e mu`` couplings."""

    model_label: str
    matching_assumption: str
    M_KK: float
    matching_scale: float
    left_emu_overlap: complex
    right_emu_overlap: complex
    scale_factor: float
    delta_g_left_emu: complex
    delta_g_right_emu: complex
    source: str
    diagnostics: Mapping[str, Any]


@dataclass(frozen=True)
class ZPoleLFVBranchingResult:
    """Branching-fraction result for ``Z -> l_i l_j`` with ``i != j``."""

    model_label: str
    input_bundle: str
    initial_flavor: str
    final_flavor: str
    branching_fraction: float
    sm_branching_fraction: float
    ratio_to_limit: float | None
    br_limit: float | None
    passes: bool | None
    delta_g_left: complex
    delta_g_right: complex
    coupling_norm: float
    lfv_width_weight: float
    sm_total_width_weight: float
    total_width_weight: float
    charge_state_factor: float
    diagnostics: Mapping[str, Any]


def z_lfv_proxy_input(
    left_emu_overlap: complex,
    right_emu_overlap: complex,
    m_kk_gev: float,
    *,
    source: str = "caller-supplied lepton neutral-current proxy",
) -> ZPoleLFVProxyInput:
    """Build a shape-checked explicit proxy input for T015-like constraints."""

    return ZPoleLFVProxyInput(
        left_emu_overlap=complex(left_emu_overlap),
        right_emu_overlap=complex(right_emu_overlap),
        m_kk_gev=_positive_float("m_kk_gev", m_kk_gev),
        source=str(source),
    )


def sm_total_width_weight(inputs: ZPoleSMInputs | None = None) -> dict[str, float]:
    """Return SM Z-width weights for quarks, charged leptons, and neutrinos."""

    p = default_sm_inputs() if inputs is None else inputs
    weights: dict[str, float] = {}
    for flavor in _QUARK_FLAVORS:
        weights[flavor] = partial_width_weight(
            sm_couplings(flavor, p),
            radiator=p.radiator_for(flavor),
        )
    for flavor in _CHARGED_LEPTONS:
        weights[flavor] = partial_width_weight(sm_couplings(flavor, p))
    for flavor in _NEUTRINO_FLAVORS:
        weights[flavor] = partial_width_weight(sm_couplings("nu", p))
    return weights


def z_lfv_branching_fraction_from_couplings(
    *,
    delta_g_left: complex,
    delta_g_right: complex,
    initial_flavor: str = "e",
    final_flavor: str = "mu",
    br_limit: float | None = None,
    inputs: ZPoleSMInputs | None = None,
    charge_state_factor: float = _DEFAULT_CHARGE_STATE_FACTOR,
) -> ZPoleLFVBranchingResult:
    """Return charge-summed ``BR(Z -> l_i l_j)`` from off-diagonal couplings."""

    p = default_sm_inputs() if inputs is None else inputs
    initial = _canonical_charged_lepton(initial_flavor)
    final = _canonical_charged_lepton(final_flavor)
    if initial == final:
        raise ValueError("LFV Z decay requires distinct charged-lepton flavors")

    charge_factor = _positive_float("charge_state_factor", charge_state_factor)
    left = complex(delta_g_left)
    right = complex(delta_g_right)
    coupling_norm = float(abs(left) ** 2 + abs(right) ** 2)
    if not math.isfinite(coupling_norm):
        raise ValueError("LFV coupling norm must be finite")

    sm_weights = sm_total_width_weight(p)
    sm_total = float(sum(sm_weights.values()))
    if sm_total <= 0.0:
        raise ValueError("SM total width weight must be positive")
    lfv_weight = float(charge_factor * coupling_norm)
    total = float(sm_total + lfv_weight)
    branching_fraction = float(lfv_weight / total)

    limit = None if br_limit is None else _bounded_probability("br_limit", br_limit)
    ratio = None if limit is None else float(branching_fraction / limit)
    passes = None if ratio is None else bool(ratio <= 1.0)

    return ZPoleLFVBranchingResult(
        model_label=ZPOLE_LFV_MODEL_V1,
        input_bundle=p.input_bundle,
        initial_flavor=initial,
        final_flavor=final,
        branching_fraction=branching_fraction,
        sm_branching_fraction=0.0,
        ratio_to_limit=ratio,
        br_limit=limit,
        passes=passes,
        delta_g_left=left,
        delta_g_right=right,
        coupling_norm=coupling_norm,
        lfv_width_weight=lfv_weight,
        sm_total_width_weight=sm_total,
        total_width_weight=total,
        charge_state_factor=charge_factor,
        diagnostics={
            "base_model_label": ZPOLE_MODEL_V1,
            "branching_formula": (
                "BR(Z -> e+- mu-+) = charge_state_factor * "
                "(|delta_g_L|^2 + |delta_g_R|^2) / "
                "(SM total Z width weight + LFV width weight)"
            ),
            "sm_width_weights": dict(sm_weights),
            "sin2_theta_eff": float(p.sin2_theta_eff),
            "m_z_gev": float(p.m_z_gev),
        },
    )


def z_lfv_effective_coupling_limit(
    br_limit: float,
    *,
    inputs: ZPoleSMInputs | None = None,
    charge_state_factor: float = _DEFAULT_CHARGE_STATE_FACTOR,
) -> float:
    """Return the limit on ``sqrt(|delta_g_L|^2 + |delta_g_R|^2)``."""

    p = default_sm_inputs() if inputs is None else inputs
    limit = _bounded_probability("br_limit", br_limit)
    charge_factor = _positive_float("charge_state_factor", charge_state_factor)
    sm_total = float(sum(sm_total_width_weight(p).values()))
    norm_limit = limit * sm_total / (charge_factor * (1.0 - limit))
    return float(math.sqrt(norm_limit))


def z_lfv_coupling_proxy(
    source: Any,
    *,
    m_kk_gev: float | None = None,
    inputs: ZPoleSMInputs | None = None,
) -> ZPoleLFVCouplingProxy:
    """Map a lepton-overlap proxy source onto off-diagonal ``Z e mu`` couplings."""

    p = default_sm_inputs() if inputs is None else inputs
    proxy_input = _coerce_proxy_input(source, m_kk_gev=m_kk_gev)
    resolved_m_kk = _positive_float(
        "m_kk_gev",
        proxy_input.m_kk_gev if m_kk_gev is None else m_kk_gev,
    )
    scale = float((p.m_z_gev / resolved_m_kk) ** 2)
    left = complex(p.proxy_strength * scale * proxy_input.left_emu_overlap)
    right = complex(p.proxy_strength * scale * proxy_input.right_emu_overlap)
    return ZPoleLFVCouplingProxy(
        model_label=ZPOLE_LFV_MODEL_V1,
        matching_assumption=ZPOLE_LFV_PROXY_V1,
        M_KK=resolved_m_kk,
        matching_scale=resolved_m_kk,
        left_emu_overlap=complex(proxy_input.left_emu_overlap),
        right_emu_overlap=complex(proxy_input.right_emu_overlap),
        scale_factor=scale,
        delta_g_left_emu=left,
        delta_g_right_emu=right,
        source=proxy_input.source,
        diagnostics={
            "m_kk_gev": float(resolved_m_kk),
            "matching_scale_gev": float(resolved_m_kk),
            "m_z_gev": float(p.m_z_gev),
            "proxy_strength": float(p.proxy_strength),
            "scale_factor": float(scale),
            "left_emu_overlap": complex(proxy_input.left_emu_overlap),
            "right_emu_overlap": complex(proxy_input.right_emu_overlap),
            "delta_g_left_emu": left,
            "delta_g_right_emu": right,
            "proxy_source": proxy_input.source,
            "matching_assumption": ZPOLE_LFV_PROXY_V1,
        },
    )


def z_lfv_branching_fraction_with_proxy(
    source: Any,
    *,
    br_limit: float | None = None,
    m_kk_gev: float | None = None,
    inputs: ZPoleSMInputs | None = None,
    charge_state_factor: float = _DEFAULT_CHARGE_STATE_FACTOR,
) -> tuple[ZPoleLFVBranchingResult, ZPoleLFVCouplingProxy]:
    """Evaluate ``Z -> e mu`` using the documented lepton-overlap proxy."""

    p = default_sm_inputs() if inputs is None else inputs
    proxy = z_lfv_coupling_proxy(source, m_kk_gev=m_kk_gev, inputs=p)
    result = z_lfv_branching_fraction_from_couplings(
        delta_g_left=proxy.delta_g_left_emu,
        delta_g_right=proxy.delta_g_right_emu,
        br_limit=br_limit,
        inputs=p,
        charge_state_factor=charge_state_factor,
    )
    return result, proxy


def _coerce_proxy_input(
    value: Any,
    *,
    m_kk_gev: float | None,
) -> ZPoleLFVProxyInput:
    if isinstance(value, ZPoleLFVProxyInput):
        return value
    if isinstance(value, Mapping):
        return _proxy_from_mapping(value, m_kk_gev=m_kk_gev)
    return _proxy_from_object(value, m_kk_gev=m_kk_gev)


def _proxy_from_mapping(
    mapping: Mapping[str, Any],
    *,
    m_kk_gev: float | None,
) -> ZPoleLFVProxyInput:
    source = str(mapping.get("source", "mapping lepton neutral-current proxy"))
    left, right = _overlaps_from_mapping(mapping)
    return ZPoleLFVProxyInput(
        left_emu_overlap=left,
        right_emu_overlap=right,
        m_kk_gev=_resolve_m_kk_from_mapping(mapping, m_kk_gev),
        source=source,
    )


def _proxy_from_object(value: Any, *, m_kk_gev: float | None) -> ZPoleLFVProxyInput:
    left = _first_present_attr(value, ("left_emu_overlap", "left_emu"))
    right = _first_present_attr(value, ("right_emu_overlap", "right_emu"))
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
                "lepton LFV proxy must provide left/right e-mu overlap entries "
                "or charged-lepton overlap matrices"
            )
        left = 0.0j if left_matrix is None else _matrix_offdiag(left_matrix, "left", 0, 1)
        right = 0.0j if right_matrix is None else _matrix_offdiag(right_matrix, "right", 0, 1)
    return ZPoleLFVProxyInput(
        left_emu_overlap=complex(0.0j if left is None else left),
        right_emu_overlap=complex(0.0j if right is None else right),
        m_kk_gev=_resolve_m_kk_from_object(value, m_kk_gev),
        source=str(getattr(value, "source", "object lepton neutral-current proxy")),
    )


def _overlaps_from_mapping(mapping: Mapping[str, Any]) -> tuple[complex, complex]:
    left = _first_present_key(mapping, ("left_emu_overlap", "left_emu"))
    right = _first_present_key(mapping, ("right_emu_overlap", "right_emu"))
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
            "mapping must provide left/right e-mu overlap entries or "
            "charged-lepton overlap matrices"
        )
    left_value = 0.0j if left_matrix is None else _matrix_offdiag(left_matrix, "left", 0, 1)
    right_value = 0.0j if right_matrix is None else _matrix_offdiag(right_matrix, "right", 0, 1)
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
    raise AttributeError("lepton LFV proxy object must provide m_kk_gev or M_KK")


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


def _canonical_charged_lepton(flavor: str) -> str:
    aliases = {
        "electron": "e",
        "e-": "e",
        "e+": "e",
        "muon": "mu",
        "mu-": "mu",
        "mu+": "mu",
        "tau-": "tau",
        "tau+": "tau",
    }
    name = aliases.get(str(flavor), str(flavor))
    if name not in _CHARGED_LEPTONS:
        raise ValueError(f"unsupported charged-lepton flavor {flavor!r}")
    return name


def _positive_float(name: str, value: object) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0.0:
        raise ValueError(f"{name} must be positive and finite")
    return number


def _bounded_probability(name: str, value: object) -> float:
    number = float(value)
    if not math.isfinite(number) or not 0.0 < number < 1.0:
        raise ValueError(f"{name} must lie between zero and one")
    return number
