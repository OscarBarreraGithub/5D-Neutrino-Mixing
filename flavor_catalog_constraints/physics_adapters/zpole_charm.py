"""Charm-specific adapter over the shared Z-pole machinery.

The reusable pseudo-observable arithmetic lives in
``flavor_catalog_constraints.physics_adapters.zpole`` / ``quarkConstraints.zpole``.
This adapter adds only the T012-specific glue:

* a charm radiator calibration so the SM-limit ``R_c`` can be aligned to the
  YAML-referenced LEP/SLC value, and
* a documented RS ``Zcc`` coupling-shift proxy.

NEEDS-HUMAN-PHYSICS: the proxy maps available up-sector KK-gluon overlap
non-universality onto ``delta g_c ~ (m_Z/M_KK)^2 Delta overlap``.  A rigorous
RS prediction needs the EW KK/Z/Z' spectrum, custodial representations,
fermion embeddings, brane kinetic terms, and Z-mixing inputs that are not
present on ``ParameterPoint``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, replace
from typing import Mapping

import numpy as np

from flavor_catalog_constraints.physics_adapters.zpole import (
    ZPOLE_MODEL_V1,
    QuarkMassBasisCouplings,
    ZPoleQuarkObservables,
    ZPoleSMInputs,
    zpole_default_sm_inputs,
    zpole_evaluate_quark,
    zpole_partial_width_weight,
    zpole_shifted_couplings,
    zpole_sm_couplings,
)

__all__ = [
    "QuarkMassBasisCouplings",
    "ZPOLE_RS_ZCC_PROXY_V1",
    "ZPoleQuarkObservables",
    "ZPoleSMInputs",
    "ZccCouplingShiftProxy",
    "zpole_inputs_with_charm_radiator",
    "zpole_zcc_coupling_shift_proxy",
    "zpole_evaluate_zcc_with_proxy",
]

ZPOLE_RS_ZCC_PROXY_V1 = (
    "NEEDS-HUMAN-PHYSICS: available quark KK-gluon overlap matrices are used "
    "as charm-vs-up neutral-current non-universality proxies; the shift "
    "delta g_c = proxy_strength * (m_Z/M_KK)^2 * (overlap_c - overlap_u) "
    "stands in for missing model-dependent RS EW KK/Z/Z' and custodial Zcc "
    "matching."
)

_HADRONIC_FLAVORS: tuple[str, ...] = ("u", "d", "s", "c", "b")


@dataclass(frozen=True)
class ZccCouplingShiftProxy:
    """Documented RS proxy for charm Z-coupling shifts."""

    model_label: str
    matching_assumption: str
    M_KK: float
    matching_scale: float
    left_charm_overlap: float
    right_charm_overlap: float
    left_up_overlap: float
    right_up_overlap: float
    left_nonuniversality: float
    right_nonuniversality: float
    scale_factor: float
    delta_g_left_c: float
    delta_g_right_c: float
    diagnostics: Mapping[str, float | str]


def zpole_inputs_with_charm_radiator(
    target_r_c: float,
    inputs: ZPoleSMInputs | None = None,
) -> ZPoleSMInputs:
    """Return inputs whose SM-limit ``R_c`` equals ``target_r_c``."""
    p = zpole_default_sm_inputs() if inputs is None else inputs
    target = _bounded_probability("target_r_c", target_r_c)
    radiators = {flavor: p.radiator_for(flavor) for flavor in _HADRONIC_FLAVORS}
    radiators["c"] = 1.0
    weights = {
        flavor: zpole_partial_width_weight(
            zpole_sm_couplings(flavor, p),
            radiator=radiators[flavor],
        )
        for flavor in _HADRONIC_FLAVORS
    }
    charm_unit = weights["c"]
    other_sum = sum(weight for flavor, weight in weights.items() if flavor != "c")
    charm_radiator = target * other_sum / (charm_unit * (1.0 - target))
    if not math.isfinite(charm_radiator) or charm_radiator <= 0.0:
        raise ValueError("calibrated charm radiator is not positive and finite")
    radiators["c"] = float(charm_radiator)
    return replace(p, quark_radiators=radiators)


def zpole_zcc_coupling_shift_proxy(
    source: QuarkMassBasisCouplings,
    *,
    m_kk_gev: float | None = None,
    inputs: ZPoleSMInputs | None = None,
) -> ZccCouplingShiftProxy:
    """Map available mass-basis overlaps onto a documented ``Zcc`` proxy."""
    p = zpole_default_sm_inputs() if inputs is None else inputs
    resolved_m_kk = _positive_float(
        "m_kk_gev",
        source.M_KK if m_kk_gev is None else m_kk_gev,
    )
    left_overlap = _left_up_overlap(source)
    right_overlap = _overlap_matrix(source, "right_up_overlap", "right_up")
    left_c = _real_diagonal(left_overlap, 1, "left_charm_overlap")
    right_c = _real_diagonal(right_overlap, 1, "right_charm_overlap")
    left_u = _real_diagonal(left_overlap, 0, "left_up_overlap")
    right_u = _real_diagonal(right_overlap, 0, "right_up_overlap")
    left_nonuniversality = float(left_c - left_u)
    right_nonuniversality = float(right_c - right_u)
    scale = float((p.m_z_gev / resolved_m_kk) ** 2)
    delta_left = float(p.proxy_strength * scale * left_nonuniversality)
    delta_right = float(p.proxy_strength * scale * right_nonuniversality)
    return ZccCouplingShiftProxy(
        model_label=ZPOLE_MODEL_V1,
        matching_assumption=ZPOLE_RS_ZCC_PROXY_V1,
        M_KK=resolved_m_kk,
        matching_scale=resolved_m_kk,
        left_charm_overlap=float(left_c),
        right_charm_overlap=float(right_c),
        left_up_overlap=float(left_u),
        right_up_overlap=float(right_u),
        left_nonuniversality=left_nonuniversality,
        right_nonuniversality=right_nonuniversality,
        scale_factor=scale,
        delta_g_left_c=delta_left,
        delta_g_right_c=delta_right,
        diagnostics={
            "m_kk_gev": float(resolved_m_kk),
            "matching_scale_gev": float(resolved_m_kk),
            "m_z_gev": float(p.m_z_gev),
            "proxy_strength": float(p.proxy_strength),
            "scale_factor": scale,
            "left_charm_overlap": float(left_c),
            "right_charm_overlap": float(right_c),
            "left_up_overlap": float(left_u),
            "right_up_overlap": float(right_u),
            "left_nonuniversality": left_nonuniversality,
            "right_nonuniversality": right_nonuniversality,
            "delta_g_left_c": delta_left,
            "delta_g_right_c": delta_right,
            "matching_assumption": ZPOLE_RS_ZCC_PROXY_V1,
        },
    )


def zpole_evaluate_zcc_with_proxy(
    source: QuarkMassBasisCouplings,
    *,
    m_kk_gev: float | None = None,
    inputs: ZPoleSMInputs | None = None,
) -> tuple[ZPoleQuarkObservables, ZccCouplingShiftProxy]:
    """Evaluate ``Z -> c cbar`` pseudo-observables with the RS proxy shift."""
    p = zpole_default_sm_inputs() if inputs is None else inputs
    proxy = zpole_zcc_coupling_shift_proxy(source, m_kk_gev=m_kk_gev, inputs=p)
    shifted_charm = zpole_shifted_couplings(
        zpole_sm_couplings("c", p),
        delta_g_left=proxy.delta_g_left_c,
        delta_g_right=proxy.delta_g_right_c,
    )
    observables = zpole_evaluate_quark("c", {"c": shifted_charm}, inputs=p)
    return observables, proxy


def _left_up_overlap(source: QuarkMassBasisCouplings) -> np.ndarray:
    left_up = getattr(source, "left_up", None)
    if left_up is not None:
        g_s = _positive_float("g_s", getattr(source, "g_s", 1.0))
        return _matrix(source, "left_up") / g_s
    return _overlap_matrix(source, "left_overlap", "left_down")


def _overlap_matrix(
    source: QuarkMassBasisCouplings,
    overlap_name: str,
    coupling_name: str,
) -> np.ndarray:
    overlap = getattr(source, overlap_name, None)
    if overlap is not None:
        return _matrix(source, overlap_name)
    coupling = _matrix(source, coupling_name)
    g_s = _positive_float("g_s", getattr(source, "g_s", 1.0))
    return coupling / g_s


def _matrix(source: object, matrix_name: str) -> np.ndarray:
    matrix = np.asarray(getattr(source, matrix_name), dtype=np.complex128)
    if matrix.shape != (3, 3):
        raise ValueError(f"{matrix_name} must have shape (3, 3)")
    if not np.all(np.isfinite(matrix.real)) or not np.all(np.isfinite(matrix.imag)):
        raise ValueError(f"{matrix_name} entries must be finite")
    return matrix


def _real_diagonal(matrix: np.ndarray, index: int, name: str) -> float:
    value = complex(matrix[index, index])
    if abs(value.imag) > 1.0e-12:
        raise ValueError(f"{name} diagonal overlap must be real")
    if not math.isfinite(value.real):
        raise ValueError(f"{name} diagonal overlap must be finite")
    return float(value.real)


def _positive_float(name: str, value: object) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0.0:
        raise ValueError(f"{name} must be positive and finite")
    return number


def _bounded_probability(name: str, value: float) -> float:
    number = float(value)
    if not math.isfinite(number) or not 0.0 < number < 1.0:
        raise ValueError(f"{name} must lie between zero and one")
    return number
