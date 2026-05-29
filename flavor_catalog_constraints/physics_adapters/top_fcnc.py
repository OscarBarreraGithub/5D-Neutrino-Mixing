"""Adapter over :mod:`quarkConstraints.top_fcnc`.

This is the catalog boundary for top-FCNC two-body branching-ratio machinery.
Constraint modules import this adapter only; the effective-coupling formulas
and the documented RS top-Z proxy remain isolated in ``quarkConstraints``.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
import math
from typing import Mapping

import numpy as np

from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.top_fcnc import (
    TOP_FCNC_GLUON_DIPOLE_CONVENTION,
    TOP_FCNC_HIGGS_SCALAR_CONVENTION,
    TOP_FCNC_INPUT_BUNDLE_V1,
    TOP_FCNC_MODEL_V1,
    TOP_FCNC_PHOTON_DIPOLE_CONVENTION,
    TOP_FCNC_RS_Z_MATCHING_ASSUMPTION_V1,
    TOP_FCNC_Z_VECTOR_CONVENTION,
    TopFCNCBranchingResult,
    TopFCNCSMInputs,
    TopZFCNCProxyCouplings,
    compute_top_z_fcnc_proxy as _compute_top_z_fcnc_proxy,
    default_sm_inputs as _default_sm_inputs,
    evaluate_t_to_q_z as _evaluate_t_to_q_z,
    gluon_dipole_branching_fraction as _gluon_dipole_branching_fraction,
    higgs_scalar_branching_fraction as _higgs_scalar_branching_fraction,
    photon_dipole_branching_fraction as _photon_dipole_branching_fraction,
    top_gluon_dipole_partial_width as _top_gluon_dipole_partial_width,
    top_higgs_scalar_partial_width as _top_higgs_scalar_partial_width,
    top_photon_dipole_partial_width as _top_photon_dipole_partial_width,
    top_z_vector_partial_width as _top_z_vector_partial_width,
    weak_neutral_current_coupling as _weak_neutral_current_coupling,
    z_vector_branching_fraction as _z_vector_branching_fraction,
)

__all__ = [
    "QuarkMassBasisCouplings",
    "TOP_FCNC_MODEL_V1",
    "TOP_FCNC_INPUT_BUNDLE_V1",
    "TOP_FCNC_Z_VECTOR_CONVENTION",
    "TOP_FCNC_PHOTON_DIPOLE_CONVENTION",
    "TOP_FCNC_GLUON_DIPOLE_CONVENTION",
    "TOP_FCNC_HIGGS_SCALAR_CONVENTION",
    "TOP_FCNC_RS_Z_MATCHING_ASSUMPTION_V1",
    "TOP_FCNC_RS_PHOTON_DIPOLE_PROXY_ASSUMPTION_V1",
    "TopFCNCSMInputs",
    "TopZFCNCProxyCouplings",
    "TopPhotonDipoleProxyCouplings",
    "TopFCNCBranchingResult",
    "top_fcnc_default_sm_inputs",
    "top_fcnc_weak_neutral_current_coupling",
    "top_fcnc_z_vector_partial_width",
    "top_fcnc_photon_dipole_partial_width",
    "top_fcnc_gluon_dipole_partial_width",
    "top_fcnc_higgs_scalar_partial_width",
    "top_fcnc_z_vector_branching_fraction",
    "top_fcnc_photon_dipole_branching_fraction",
    "top_fcnc_gluon_dipole_branching_fraction",
    "top_fcnc_higgs_scalar_branching_fraction",
    "top_photon_dipole_proxy_from_couplings",
    "top_z_fcnc_proxy_from_couplings",
    "t_to_q_gamma_from_couplings",
    "t_to_q_z_from_couplings",
]

TOP_FCNC_RS_PHOTON_DIPOLE_PROXY_ASSUMPTION_V1 = (
    "NEEDS-HUMAN-PHYSICS: up-sector quark KK-gluon mass-basis couplings are "
    "used only as neutral-current flavor-overlap proxies; the effective "
    "t-q-gamma photon dipole is lambda_qt=(g_qt/g_s)*(m_t^2/M_KK^2), with "
    "no loop, charge, Yukawa, or collider-recast matching. This is a "
    "dimensionally normalized diagnostic proxy for the full RS electromagnetic "
    "dipole calculation."
)


@dataclass(frozen=True)
class TopPhotonDipoleProxyCouplings:
    """Documented RS-proxy effective ``t-q-gamma`` dipole couplings."""

    model_label: str
    operator_convention: str
    matching_assumption: str
    light_quark: str
    light_up_index: int
    M_KK: float
    matching_scale: float
    left_qt_coupling: complex
    right_qt_coupling: complex
    left_qt_overlap: complex
    right_qt_overlap: complex
    dipole_scale_suppression: float
    dipole_left: complex
    dipole_right: complex

    @property
    def effective_couplings(self) -> Mapping[str, complex]:
        return {
            "lambda_L": complex(self.dipole_left),
            "lambda_R": complex(self.dipole_right),
        }


def top_fcnc_default_sm_inputs() -> TopFCNCSMInputs:
    """Return the default top-FCNC input bundle."""
    return _default_sm_inputs()


def top_fcnc_weak_neutral_current_coupling(
    inputs: TopFCNCSMInputs | None = None,
) -> float:
    """Return ``g/(2 c_W)`` from the top-FCNC input bundle."""
    return _weak_neutral_current_coupling(inputs)


def top_fcnc_z_vector_partial_width(
    vector_left: complex = 0.0j,
    vector_right: complex = 0.0j,
    *,
    inputs: TopFCNCSMInputs | None = None,
) -> float:
    """Return ``Gamma(t -> q Z)`` for vector effective couplings."""
    return _top_z_vector_partial_width(vector_left, vector_right, inputs=inputs)


def top_fcnc_photon_dipole_partial_width(
    dipole_left: complex = 0.0j,
    dipole_right: complex = 0.0j,
    *,
    inputs: TopFCNCSMInputs | None = None,
) -> float:
    """Return ``Gamma(t -> q gamma)`` for photon dipole couplings."""
    return _top_photon_dipole_partial_width(dipole_left, dipole_right, inputs=inputs)


def top_fcnc_gluon_dipole_partial_width(
    dipole_left: complex = 0.0j,
    dipole_right: complex = 0.0j,
    *,
    inputs: TopFCNCSMInputs | None = None,
) -> float:
    """Return ``Gamma(t -> q g)`` for chromomagnetic dipole couplings."""
    return _top_gluon_dipole_partial_width(dipole_left, dipole_right, inputs=inputs)


def top_fcnc_higgs_scalar_partial_width(
    scalar_left: complex = 0.0j,
    scalar_right: complex = 0.0j,
    *,
    inputs: TopFCNCSMInputs | None = None,
) -> float:
    """Return ``Gamma(t -> q h)`` for scalar FCNC couplings."""
    return _top_higgs_scalar_partial_width(scalar_left, scalar_right, inputs=inputs)


def top_fcnc_z_vector_branching_fraction(
    *,
    light_quark: str,
    vector_left: complex = 0.0j,
    vector_right: complex = 0.0j,
    inputs: TopFCNCSMInputs | None = None,
) -> TopFCNCBranchingResult:
    """Evaluate ``BR(t -> q Z)`` from vector effective couplings."""
    return _z_vector_branching_fraction(
        light_quark=light_quark,
        vector_left=vector_left,
        vector_right=vector_right,
        inputs=inputs,
    )


def top_fcnc_photon_dipole_branching_fraction(
    *,
    light_quark: str,
    dipole_left: complex = 0.0j,
    dipole_right: complex = 0.0j,
    inputs: TopFCNCSMInputs | None = None,
) -> TopFCNCBranchingResult:
    """Evaluate ``BR(t -> q gamma)`` from photon dipole couplings."""
    return _photon_dipole_branching_fraction(
        light_quark=light_quark,
        dipole_left=dipole_left,
        dipole_right=dipole_right,
        inputs=inputs,
    )


def top_fcnc_gluon_dipole_branching_fraction(
    *,
    light_quark: str,
    dipole_left: complex = 0.0j,
    dipole_right: complex = 0.0j,
    inputs: TopFCNCSMInputs | None = None,
) -> TopFCNCBranchingResult:
    """Evaluate ``BR(t -> q g)`` from chromomagnetic dipole couplings."""
    return _gluon_dipole_branching_fraction(
        light_quark=light_quark,
        dipole_left=dipole_left,
        dipole_right=dipole_right,
        inputs=inputs,
    )


def top_fcnc_higgs_scalar_branching_fraction(
    *,
    light_quark: str,
    scalar_left: complex = 0.0j,
    scalar_right: complex = 0.0j,
    inputs: TopFCNCSMInputs | None = None,
) -> TopFCNCBranchingResult:
    """Evaluate ``BR(t -> q h)`` from scalar effective couplings."""
    return _higgs_scalar_branching_fraction(
        light_quark=light_quark,
        scalar_left=scalar_left,
        scalar_right=scalar_right,
        inputs=inputs,
    )


def _finite_complex(value: complex, name: str) -> complex:
    number = complex(value)
    if not math.isfinite(number.real) or not math.isfinite(number.imag):
        raise ValueError(f"{name} must be finite")
    return number


def _positive_float(value: object, name: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0.0:
        raise ValueError(f"{name} must be positive and finite")
    return number


def _matrix_entry(source: object, matrix_name: str, i: int, j: int) -> complex:
    matrix = np.asarray(getattr(source, matrix_name), dtype=np.complex128)
    if matrix.shape != (3, 3):
        raise ValueError(f"{matrix_name} must have shape (3, 3)")
    return _finite_complex(complex(matrix[i, j]), f"{matrix_name}[{i},{j}]")


def top_photon_dipole_proxy_from_couplings(
    couplings: QuarkMassBasisCouplings,
    *,
    light_quark: str,
    light_up_index: int,
    m_kk_gev: float | None = None,
    inputs: TopFCNCSMInputs | None = None,
) -> TopPhotonDipoleProxyCouplings:
    """Return the v1 documented ``t-q-gamma`` photon-dipole proxy."""
    if light_up_index not in (0, 1):
        raise ValueError("light_up_index must be 0 (u) or 1 (c)")
    p = top_fcnc_default_sm_inputs() if inputs is None else inputs
    resolved_m_kk = _positive_float(
        getattr(couplings, "M_KK") if m_kk_gev is None else m_kk_gev,
        "m_kk_gev",
    )
    left_qt = _matrix_entry(couplings, "left_up", light_up_index, 2)
    right_qt = _matrix_entry(couplings, "right_up", light_up_index, 2)
    g_s = _positive_float(getattr(couplings, "g_s", 1.0), "g_s")
    left_overlap = left_qt / g_s
    right_overlap = right_qt / g_s
    suppression = (p.m_top_gev / resolved_m_kk) ** 2
    dipole_left = left_overlap * suppression
    dipole_right = right_overlap * suppression
    return TopPhotonDipoleProxyCouplings(
        model_label=TOP_FCNC_MODEL_V1,
        operator_convention=TOP_FCNC_PHOTON_DIPOLE_CONVENTION,
        matching_assumption=TOP_FCNC_RS_PHOTON_DIPOLE_PROXY_ASSUMPTION_V1,
        light_quark=light_quark,
        light_up_index=int(light_up_index),
        M_KK=float(resolved_m_kk),
        matching_scale=float(resolved_m_kk),
        left_qt_coupling=complex(left_qt),
        right_qt_coupling=complex(right_qt),
        left_qt_overlap=complex(left_overlap),
        right_qt_overlap=complex(right_overlap),
        dipole_scale_suppression=float(suppression),
        dipole_left=complex(dipole_left),
        dipole_right=complex(dipole_right),
    )


def top_z_fcnc_proxy_from_couplings(
    couplings: QuarkMassBasisCouplings,
    *,
    light_quark: str,
    light_up_index: int,
    m_kk_gev: float | None = None,
    inputs: TopFCNCSMInputs | None = None,
) -> TopZFCNCProxyCouplings:
    """Return the v1 documented ``t-q-Z`` proxy from mass-basis couplings."""
    return _compute_top_z_fcnc_proxy(
        couplings,
        light_quark=light_quark,
        light_up_index=light_up_index,
        m_kk_gev=m_kk_gev,
        inputs=inputs,
    )


def t_to_q_gamma_from_couplings(
    couplings: QuarkMassBasisCouplings | TopPhotonDipoleProxyCouplings | None = None,
    *,
    light_quark: str,
    light_up_index: int,
    m_kk_gev: float | None = None,
    inputs: TopFCNCSMInputs | None = None,
) -> TopFCNCBranchingResult:
    """Evaluate ``BR(t -> q gamma)`` with the v1 photon-dipole proxy."""
    p = top_fcnc_default_sm_inputs() if inputs is None else inputs
    if couplings is None:
        return top_fcnc_photon_dipole_branching_fraction(
            light_quark=light_quark,
            dipole_left=0.0j,
            dipole_right=0.0j,
            inputs=p,
        )
    if isinstance(couplings, TopPhotonDipoleProxyCouplings):
        proxy = couplings
    else:
        proxy = top_photon_dipole_proxy_from_couplings(
            couplings,
            light_quark=light_quark,
            light_up_index=light_up_index,
            m_kk_gev=m_kk_gev,
            inputs=p,
        )
    result = top_fcnc_photon_dipole_branching_fraction(
        light_quark=light_quark,
        dipole_left=proxy.dipole_left,
        dipole_right=proxy.dipole_right,
        inputs=p,
    )
    diagnostics = dict(result.diagnostics)
    diagnostics.update(
        {
            "matching_assumption": proxy.matching_assumption,
            "m_kk_gev": float(proxy.M_KK),
            "matching_scale_gev": float(proxy.matching_scale),
            "left_qt_coupling": complex(proxy.left_qt_coupling),
            "right_qt_coupling": complex(proxy.right_qt_coupling),
            "left_qt_overlap": complex(proxy.left_qt_overlap),
            "right_qt_overlap": complex(proxy.right_qt_overlap),
            "dipole_scale_suppression": float(proxy.dipole_scale_suppression),
        }
    )
    return replace(result, diagnostics=diagnostics)


def t_to_q_z_from_couplings(
    couplings: QuarkMassBasisCouplings | TopZFCNCProxyCouplings | None = None,
    *,
    light_quark: str,
    light_up_index: int,
    m_kk_gev: float | None = None,
    inputs: TopFCNCSMInputs | None = None,
) -> TopFCNCBranchingResult:
    """Evaluate ``BR(t -> q Z)`` from mass-basis couplings or a proxy object."""
    return _evaluate_t_to_q_z(
        couplings,
        light_quark=light_quark,
        light_up_index=light_up_index,
        m_kk_gev=m_kk_gev,
        inputs=inputs,
    )


TOP_FCNC_RS_GLUON_DIPOLE_PROXY_ASSUMPTION_V1 = (
    "NEEDS-HUMAN-PHYSICS: up-sector quark KK-gluon mass-basis couplings are "
    "used only as chromomagnetic flavor-overlap proxies; the effective "
    "t-q-g dipole is zeta_qt=(g_qt/g_s)*(m_t^2/M_KK^2), with no loop/Yukawa "
    "matching, SMEFT C_uG normalization, QCD running, or cg->t collider-recast "
    "matching. This is a dimensionally normalized diagnostic proxy for the "
    "full RS chromomagnetic dipole calculation."
)


@dataclass(frozen=True)
class TopGluonDipoleProxyCouplings:
    """Documented RS-proxy effective ``t-q-g`` chromomagnetic dipoles."""

    model_label: str
    operator_convention: str
    matching_assumption: str
    light_quark: str
    light_up_index: int
    M_KK: float
    matching_scale: float
    left_qt_coupling: complex
    right_qt_coupling: complex
    left_qt_overlap: complex
    right_qt_overlap: complex
    dipole_scale_suppression: float
    dipole_left: complex
    dipole_right: complex

    @property
    def effective_couplings(self) -> Mapping[str, complex]:
        return {
            "zeta_L": complex(self.dipole_left),
            "zeta_R": complex(self.dipole_right),
        }


def top_gluon_dipole_proxy_from_couplings(
    couplings: QuarkMassBasisCouplings,
    *,
    light_quark: str,
    light_up_index: int,
    m_kk_gev: float | None = None,
    inputs: TopFCNCSMInputs | None = None,
) -> TopGluonDipoleProxyCouplings:
    """Return the v1 documented ``t-q-g`` chromomagnetic-dipole proxy."""
    if light_up_index not in (0, 1):
        raise ValueError("light_up_index must be 0 (u) or 1 (c)")
    p = top_fcnc_default_sm_inputs() if inputs is None else inputs
    resolved_m_kk = _positive_float(
        getattr(couplings, "M_KK") if m_kk_gev is None else m_kk_gev,
        "m_kk_gev",
    )
    left_qt = _matrix_entry(couplings, "left_up", light_up_index, 2)
    right_qt = _matrix_entry(couplings, "right_up", light_up_index, 2)
    g_s = _positive_float(getattr(couplings, "g_s", 1.0), "g_s")
    left_overlap = left_qt / g_s
    right_overlap = right_qt / g_s
    suppression = (p.m_top_gev / resolved_m_kk) ** 2
    dipole_left = left_overlap * suppression
    dipole_right = right_overlap * suppression
    return TopGluonDipoleProxyCouplings(
        model_label=TOP_FCNC_MODEL_V1,
        operator_convention=TOP_FCNC_GLUON_DIPOLE_CONVENTION,
        matching_assumption=TOP_FCNC_RS_GLUON_DIPOLE_PROXY_ASSUMPTION_V1,
        light_quark=light_quark,
        light_up_index=int(light_up_index),
        M_KK=float(resolved_m_kk),
        matching_scale=float(resolved_m_kk),
        left_qt_coupling=complex(left_qt),
        right_qt_coupling=complex(right_qt),
        left_qt_overlap=complex(left_overlap),
        right_qt_overlap=complex(right_overlap),
        dipole_scale_suppression=float(suppression),
        dipole_left=complex(dipole_left),
        dipole_right=complex(dipole_right),
    )


def t_to_q_gluon_from_couplings(
    couplings: QuarkMassBasisCouplings | TopGluonDipoleProxyCouplings | None = None,
    *,
    light_quark: str,
    light_up_index: int,
    m_kk_gev: float | None = None,
    inputs: TopFCNCSMInputs | None = None,
) -> TopFCNCBranchingResult:
    """Evaluate ``BR(t -> q g)`` with the v1 chromomagnetic-dipole proxy."""
    p = top_fcnc_default_sm_inputs() if inputs is None else inputs
    if couplings is None:
        return top_fcnc_gluon_dipole_branching_fraction(
            light_quark=light_quark,
            dipole_left=0.0j,
            dipole_right=0.0j,
            inputs=p,
        )
    if isinstance(couplings, TopGluonDipoleProxyCouplings):
        proxy = couplings
    else:
        proxy = top_gluon_dipole_proxy_from_couplings(
            couplings,
            light_quark=light_quark,
            light_up_index=light_up_index,
            m_kk_gev=m_kk_gev,
            inputs=p,
        )
    result = top_fcnc_gluon_dipole_branching_fraction(
        light_quark=light_quark,
        dipole_left=proxy.dipole_left,
        dipole_right=proxy.dipole_right,
        inputs=p,
    )
    diagnostics = dict(result.diagnostics)
    diagnostics.update(
        {
            "matching_assumption": proxy.matching_assumption,
            "m_kk_gev": float(proxy.M_KK),
            "matching_scale_gev": float(proxy.matching_scale),
            "left_qt_coupling": complex(proxy.left_qt_coupling),
            "right_qt_coupling": complex(proxy.right_qt_coupling),
            "left_qt_overlap": complex(proxy.left_qt_overlap),
            "right_qt_overlap": complex(proxy.right_qt_overlap),
            "dipole_scale_suppression": float(proxy.dipole_scale_suppression),
        }
    )
    return replace(result, diagnostics=diagnostics)


__all__ += [
    "TOP_FCNC_RS_GLUON_DIPOLE_PROXY_ASSUMPTION_V1",
    "TopGluonDipoleProxyCouplings",
    "top_gluon_dipole_proxy_from_couplings",
    "t_to_q_gluon_from_couplings",
]


TOP_FCNC_RS_HIGGS_YUKAWA_PROXY_ASSUMPTION_V1 = (
    "NEEDS-HUMAN-PHYSICS: up-sector quark KK-gluon mass-basis couplings are "
    "used only as scalar-Yukawa flavor-overlap proxies; the effective "
    "t-q-H Yukawas are y_qt=(g_qt/g_s)*(m_t^2/M_KK^2), with no Higgs-sector "
    "Yukawa-misalignment calculation, SMEFT normalization, EW KK/Z/Z' mixing, "
    "or combined tH-production plus ttbar-decay collider recast. This is a "
    "dimensionally normalized diagnostic proxy for the full RS top-Higgs FCNC "
    "calculation."
)


@dataclass(frozen=True)
class TopHiggsScalarProxyCouplings:
    """Documented RS-proxy effective ``t-q-H`` scalar Yukawa couplings."""

    model_label: str
    operator_convention: str
    matching_assumption: str
    light_quark: str
    light_up_index: int
    M_KK: float
    matching_scale: float
    left_qt_coupling: complex
    right_qt_coupling: complex
    left_qt_overlap: complex
    right_qt_overlap: complex
    scalar_yukawa_suppression: float
    scalar_left: complex
    scalar_right: complex

    @property
    def effective_couplings(self) -> Mapping[str, complex]:
        return {
            "y_L": complex(self.scalar_left),
            "y_R": complex(self.scalar_right),
        }


def top_higgs_scalar_proxy_from_couplings(
    couplings: QuarkMassBasisCouplings,
    *,
    light_quark: str,
    light_up_index: int,
    m_kk_gev: float | None = None,
    inputs: TopFCNCSMInputs | None = None,
) -> TopHiggsScalarProxyCouplings:
    """Return the v1 documented ``t-q-H`` scalar-Yukawa proxy."""
    if light_up_index not in (0, 1):
        raise ValueError("light_up_index must be 0 (u) or 1 (c)")
    p = top_fcnc_default_sm_inputs() if inputs is None else inputs
    resolved_m_kk = _positive_float(
        getattr(couplings, "M_KK") if m_kk_gev is None else m_kk_gev,
        "m_kk_gev",
    )
    left_qt = _matrix_entry(couplings, "left_up", light_up_index, 2)
    right_qt = _matrix_entry(couplings, "right_up", light_up_index, 2)
    g_s = _positive_float(getattr(couplings, "g_s", 1.0), "g_s")
    left_overlap = left_qt / g_s
    right_overlap = right_qt / g_s
    suppression = (p.m_top_gev / resolved_m_kk) ** 2
    scalar_left = left_overlap * suppression
    scalar_right = right_overlap * suppression
    return TopHiggsScalarProxyCouplings(
        model_label=TOP_FCNC_MODEL_V1,
        operator_convention=TOP_FCNC_HIGGS_SCALAR_CONVENTION,
        matching_assumption=TOP_FCNC_RS_HIGGS_YUKAWA_PROXY_ASSUMPTION_V1,
        light_quark=light_quark,
        light_up_index=int(light_up_index),
        M_KK=float(resolved_m_kk),
        matching_scale=float(resolved_m_kk),
        left_qt_coupling=complex(left_qt),
        right_qt_coupling=complex(right_qt),
        left_qt_overlap=complex(left_overlap),
        right_qt_overlap=complex(right_overlap),
        scalar_yukawa_suppression=float(suppression),
        scalar_left=complex(scalar_left),
        scalar_right=complex(scalar_right),
    )


def t_to_q_higgs_from_couplings(
    couplings: QuarkMassBasisCouplings | TopHiggsScalarProxyCouplings | None = None,
    *,
    light_quark: str,
    light_up_index: int,
    m_kk_gev: float | None = None,
    inputs: TopFCNCSMInputs | None = None,
) -> TopFCNCBranchingResult:
    """Evaluate ``BR(t -> q H)`` with the v1 scalar-Yukawa proxy."""
    p = top_fcnc_default_sm_inputs() if inputs is None else inputs
    if couplings is None:
        return top_fcnc_higgs_scalar_branching_fraction(
            light_quark=light_quark,
            scalar_left=0.0j,
            scalar_right=0.0j,
            inputs=p,
        )
    if isinstance(couplings, TopHiggsScalarProxyCouplings):
        proxy = couplings
    else:
        proxy = top_higgs_scalar_proxy_from_couplings(
            couplings,
            light_quark=light_quark,
            light_up_index=light_up_index,
            m_kk_gev=m_kk_gev,
            inputs=p,
        )
    result = top_fcnc_higgs_scalar_branching_fraction(
        light_quark=light_quark,
        scalar_left=proxy.scalar_left,
        scalar_right=proxy.scalar_right,
        inputs=p,
    )
    diagnostics = dict(result.diagnostics)
    diagnostics.update(
        {
            "matching_assumption": proxy.matching_assumption,
            "m_kk_gev": float(proxy.M_KK),
            "matching_scale_gev": float(proxy.matching_scale),
            "left_qt_coupling": complex(proxy.left_qt_coupling),
            "right_qt_coupling": complex(proxy.right_qt_coupling),
            "left_qt_overlap": complex(proxy.left_qt_overlap),
            "right_qt_overlap": complex(proxy.right_qt_overlap),
            "scalar_yukawa_suppression": float(proxy.scalar_yukawa_suppression),
        }
    )
    return replace(result, diagnostics=diagnostics)


__all__ += [
    "TOP_FCNC_RS_HIGGS_YUKAWA_PROXY_ASSUMPTION_V1",
    "TopHiggsScalarProxyCouplings",
    "top_higgs_scalar_proxy_from_couplings",
    "t_to_q_higgs_from_couplings",
]
