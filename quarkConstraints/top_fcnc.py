"""Reusable top-FCNC branching-ratio machinery.

Effective-coupling convention
-----------------------------
The v1 formulas neglect the light-quark mass and use these dimensionless
couplings:

``t -> q Z``
    ``L = (g / 2 c_W) Z_mu qbar gamma^mu (X_L P_L + X_R P_R) t + h.c.``

``t -> q gamma``
    ``L = e A_mu qbar i sigma^{mu nu} k_nu / m_t
    (lambda_L P_L + lambda_R P_R) t + h.c.``

``t -> q g``
    Same dipole convention as the photon mode, with ``e`` replaced by
    ``g_s T^a``.

``t -> q h``
    ``L = -h qbar (y_L P_L + y_R P_R) t + h.c.``

RS matching status
------------------
NEEDS-HUMAN-PHYSICS: the current ``ParameterPoint`` carries quark mass-basis
``KK-gluon``-style coupling matrices, not the electroweak KK/Z/Z' tower,
Higgs-sector, dipole, or collider-recast inputs required for a complete
top-FCNC prediction.  The only model-facing helper here is therefore a
documented top-Z proxy: divide the up-sector mass-basis coupling by ``g_s`` to
keep the flavor-overlap structure and convert it to a dimensionless
``t-q-Z`` vector coupling with a single Z-like mixing suppression
``m_Z^2 / M_KK^2``.  This is useful for scan diagnostics, but it is not a
substitute for a human-approved RS-to-SMEFT matching convention.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Mapping

import numpy as np

from .couplings import QuarkMassBasisCouplings

TOP_FCNC_MODEL_V1 = "top_fcnc_effective_couplings_v1"
TOP_FCNC_INPUT_BUNDLE_V1 = "top_fcnc_pdg_like_sm_inputs_v1"
TOP_FCNC_Z_VECTOR_CONVENTION = (
    "L=(g/2cW) Z_mu qbar gamma^mu (X_L P_L + X_R P_R) t + h.c.; "
    "light-quark mass neglected"
)
TOP_FCNC_PHOTON_DIPOLE_CONVENTION = (
    "L=e A_mu qbar i sigma^{mu nu} k_nu/m_t "
    "(lambda_L P_L + lambda_R P_R) t + h.c."
)
TOP_FCNC_GLUON_DIPOLE_CONVENTION = (
    "L=g_s G^a_mu qbar T^a i sigma^{mu nu} k_nu/m_t "
    "(zeta_L P_L + zeta_R P_R) t + h.c."
)
TOP_FCNC_HIGGS_SCALAR_CONVENTION = (
    "L=-h qbar (y_L P_L + y_R P_R) t + h.c.; light-quark mass neglected"
)
TOP_FCNC_RS_Z_MATCHING_ASSUMPTION_V1 = (
    "NEEDS-HUMAN-PHYSICS: up-sector quark KK-gluon mass-basis couplings are "
    "used as neutral-current flavor-overlap proxies; the effective t-q-Z "
    "vector coupling is X_qt=(g_qt/g_s)*(m_Z^2/M_KK^2), standing in for the "
    "full RS EW KK/Z/Z' matching."
)


@dataclass(frozen=True)
class TopFCNCSMInputs:
    """Numerical inputs for two-body top-FCNC widths."""

    input_bundle: str = TOP_FCNC_INPUT_BUNDLE_V1
    m_top_gev: float = 172.76
    m_z_gev: float = 91.1876
    m_h_gev: float = 125.25
    alpha_em_mz: float = 1.0 / 127.952
    alpha_s_mt: float = 0.108
    sin2_theta_w: float = 0.23122
    total_top_width_gev: float = 1.41
    constants_citation: str = (
        "PDG-like electroweak/top inputs; width formulas follow the effective "
        "top-FCNC coupling conventions documented in this module."
    )

    def __post_init__(self) -> None:
        for name in (
            "m_top_gev",
            "m_z_gev",
            "m_h_gev",
            "alpha_em_mz",
            "alpha_s_mt",
            "sin2_theta_w",
            "total_top_width_gev",
        ):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be positive and finite")
        if not 0.0 < float(self.sin2_theta_w) < 1.0:
            raise ValueError("sin2_theta_w must lie between zero and one")
        if float(self.m_z_gev) >= float(self.m_top_gev):
            raise ValueError("m_z_gev must be below m_top_gev")
        if float(self.m_h_gev) >= float(self.m_top_gev):
            raise ValueError("m_h_gev must be below m_top_gev")


@dataclass(frozen=True)
class TopZFCNCProxyCouplings:
    """Documented RS-proxy effective ``t-q-Z`` vector couplings."""

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
    z_mixing_suppression: float
    vector_left: complex
    vector_right: complex

    @property
    def effective_couplings(self) -> Mapping[str, complex]:
        return {
            "X_L": complex(self.vector_left),
            "X_R": complex(self.vector_right),
        }


@dataclass(frozen=True)
class TopFCNCBranchingResult:
    """Branching-ratio result for one top-FCNC two-body decay."""

    model_label: str
    input_bundle: str
    process: str
    light_quark: str
    coupling_kind: str
    branching_fraction: float
    partial_width_gev: float
    total_width_gev: float
    coupling_left: complex
    coupling_right: complex
    proxy: TopZFCNCProxyCouplings | None = None
    diagnostics: Mapping[str, float | complex | str | bool] = field(default_factory=dict)


def default_sm_inputs() -> TopFCNCSMInputs:
    """Return the repo-owned default top-FCNC input bundle."""
    return TopFCNCSMInputs()


def weak_neutral_current_coupling(inputs: TopFCNCSMInputs | None = None) -> float:
    """Return ``g / (2 c_W)`` from the input electroweak constants."""
    p = default_sm_inputs() if inputs is None else inputs
    g_weak = math.sqrt(4.0 * math.pi * p.alpha_em_mz / p.sin2_theta_w)
    cos_theta_w = math.sqrt(1.0 - p.sin2_theta_w)
    return float(g_weak / (2.0 * cos_theta_w))


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


def _massive_vector_phase_space(parent_mass: float, vector_mass: float) -> float:
    rho = (vector_mass / parent_mass) ** 2
    return float((1.0 - rho) ** 2 * (1.0 + 2.0 * rho))


def top_z_vector_partial_width(
    vector_left: complex = 0.0j,
    vector_right: complex = 0.0j,
    *,
    inputs: TopFCNCSMInputs | None = None,
) -> float:
    """Return ``Gamma(t -> q Z)`` for vector ``X_L, X_R`` couplings."""
    p = default_sm_inputs() if inputs is None else inputs
    left = _finite_complex(vector_left, "vector_left")
    right = _finite_complex(vector_right, "vector_right")
    norm = abs(left) ** 2 + abs(right) ** 2
    phase = _massive_vector_phase_space(p.m_top_gev, p.m_z_gev)
    g_z = weak_neutral_current_coupling(p)
    width = (
        g_z**2
        * p.m_top_gev**3
        / (32.0 * math.pi * p.m_z_gev**2)
        * phase
        * norm
    )
    return float(width)


def top_photon_dipole_partial_width(
    dipole_left: complex = 0.0j,
    dipole_right: complex = 0.0j,
    *,
    inputs: TopFCNCSMInputs | None = None,
) -> float:
    """Return ``Gamma(t -> q gamma)`` in the module dipole convention."""
    p = default_sm_inputs() if inputs is None else inputs
    left = _finite_complex(dipole_left, "dipole_left")
    right = _finite_complex(dipole_right, "dipole_right")
    return float(0.5 * p.alpha_em_mz * p.m_top_gev * (abs(left) ** 2 + abs(right) ** 2))


def top_gluon_dipole_partial_width(
    dipole_left: complex = 0.0j,
    dipole_right: complex = 0.0j,
    *,
    inputs: TopFCNCSMInputs | None = None,
) -> float:
    """Return ``Gamma(t -> q g)`` in the module chromomagnetic convention."""
    p = default_sm_inputs() if inputs is None else inputs
    left = _finite_complex(dipole_left, "dipole_left")
    right = _finite_complex(dipole_right, "dipole_right")
    color_factor = 4.0 / 3.0
    return float(
        0.5
        * color_factor
        * p.alpha_s_mt
        * p.m_top_gev
        * (abs(left) ** 2 + abs(right) ** 2)
    )


def top_higgs_scalar_partial_width(
    scalar_left: complex = 0.0j,
    scalar_right: complex = 0.0j,
    *,
    inputs: TopFCNCSMInputs | None = None,
) -> float:
    """Return ``Gamma(t -> q h)`` in the module scalar-coupling convention."""
    p = default_sm_inputs() if inputs is None else inputs
    left = _finite_complex(scalar_left, "scalar_left")
    right = _finite_complex(scalar_right, "scalar_right")
    rho = (p.m_h_gev / p.m_top_gev) ** 2
    return float(
        p.m_top_gev
        / (32.0 * math.pi)
        * (1.0 - rho) ** 2
        * (abs(left) ** 2 + abs(right) ** 2)
    )


def _branching_result(
    *,
    process: str,
    light_quark: str,
    coupling_kind: str,
    partial_width_gev: float,
    coupling_left: complex,
    coupling_right: complex,
    inputs: TopFCNCSMInputs,
    convention: str,
    proxy: TopZFCNCProxyCouplings | None = None,
    extra_diagnostics: Mapping[str, float | complex | str | bool] | None = None,
) -> TopFCNCBranchingResult:
    width = float(partial_width_gev)
    if width < 0.0 or not math.isfinite(width):
        raise ValueError("partial width must be finite and non-negative")
    diagnostics: dict[str, float | complex | str | bool] = {
        "operator_convention": convention,
        "m_top_gev": float(inputs.m_top_gev),
        "total_top_width_gev": float(inputs.total_top_width_gev),
        "constants_citation": inputs.constants_citation,
    }
    if extra_diagnostics is not None:
        diagnostics.update(dict(extra_diagnostics))
    return TopFCNCBranchingResult(
        model_label=TOP_FCNC_MODEL_V1,
        input_bundle=inputs.input_bundle,
        process=process,
        light_quark=light_quark,
        coupling_kind=coupling_kind,
        branching_fraction=float(width / inputs.total_top_width_gev),
        partial_width_gev=width,
        total_width_gev=float(inputs.total_top_width_gev),
        coupling_left=complex(coupling_left),
        coupling_right=complex(coupling_right),
        proxy=proxy,
        diagnostics=diagnostics,
    )


def z_vector_branching_fraction(
    *,
    light_quark: str,
    vector_left: complex = 0.0j,
    vector_right: complex = 0.0j,
    inputs: TopFCNCSMInputs | None = None,
    proxy: TopZFCNCProxyCouplings | None = None,
) -> TopFCNCBranchingResult:
    """Evaluate ``BR(t -> q Z)`` from vector effective couplings."""
    p = default_sm_inputs() if inputs is None else inputs
    width = top_z_vector_partial_width(vector_left, vector_right, inputs=p)
    phase = _massive_vector_phase_space(p.m_top_gev, p.m_z_gev)
    diagnostics: dict[str, float | complex | str | bool] = {
        "m_z_gev": float(p.m_z_gev),
        "sin2_theta_w": float(p.sin2_theta_w),
        "alpha_em_mz": float(p.alpha_em_mz),
        "weak_neutral_current_coupling": weak_neutral_current_coupling(p),
        "massive_vector_phase_space": phase,
        "vector_left": complex(vector_left),
        "vector_right": complex(vector_right),
    }
    if proxy is not None:
        diagnostics.update(
            {
                "matching_assumption": proxy.matching_assumption,
                "m_kk_gev": float(proxy.M_KK),
                "matching_scale_gev": float(proxy.matching_scale),
                "left_qt_coupling": complex(proxy.left_qt_coupling),
                "right_qt_coupling": complex(proxy.right_qt_coupling),
                "left_qt_overlap": complex(proxy.left_qt_overlap),
                "right_qt_overlap": complex(proxy.right_qt_overlap),
                "z_mixing_suppression": float(proxy.z_mixing_suppression),
            }
        )
    return _branching_result(
        process=f"t -> {light_quark} Z",
        light_quark=light_quark,
        coupling_kind="Z vector",
        partial_width_gev=width,
        coupling_left=vector_left,
        coupling_right=vector_right,
        inputs=p,
        convention=TOP_FCNC_Z_VECTOR_CONVENTION,
        proxy=proxy,
        extra_diagnostics=diagnostics,
    )


def photon_dipole_branching_fraction(
    *,
    light_quark: str,
    dipole_left: complex = 0.0j,
    dipole_right: complex = 0.0j,
    inputs: TopFCNCSMInputs | None = None,
) -> TopFCNCBranchingResult:
    """Evaluate ``BR(t -> q gamma)`` from photon dipole couplings."""
    p = default_sm_inputs() if inputs is None else inputs
    width = top_photon_dipole_partial_width(dipole_left, dipole_right, inputs=p)
    return _branching_result(
        process=f"t -> {light_quark} gamma",
        light_quark=light_quark,
        coupling_kind="photon dipole",
        partial_width_gev=width,
        coupling_left=dipole_left,
        coupling_right=dipole_right,
        inputs=p,
        convention=TOP_FCNC_PHOTON_DIPOLE_CONVENTION,
        extra_diagnostics={
            "alpha_em_mz": float(p.alpha_em_mz),
            "dipole_left": complex(dipole_left),
            "dipole_right": complex(dipole_right),
        },
    )


def gluon_dipole_branching_fraction(
    *,
    light_quark: str,
    dipole_left: complex = 0.0j,
    dipole_right: complex = 0.0j,
    inputs: TopFCNCSMInputs | None = None,
) -> TopFCNCBranchingResult:
    """Evaluate ``BR(t -> q g)`` from chromomagnetic dipole couplings."""
    p = default_sm_inputs() if inputs is None else inputs
    width = top_gluon_dipole_partial_width(dipole_left, dipole_right, inputs=p)
    return _branching_result(
        process=f"t -> {light_quark} g",
        light_quark=light_quark,
        coupling_kind="gluon dipole",
        partial_width_gev=width,
        coupling_left=dipole_left,
        coupling_right=dipole_right,
        inputs=p,
        convention=TOP_FCNC_GLUON_DIPOLE_CONVENTION,
        extra_diagnostics={
            "alpha_s_mt": float(p.alpha_s_mt),
            "color_factor": 4.0 / 3.0,
            "dipole_left": complex(dipole_left),
            "dipole_right": complex(dipole_right),
        },
    )


def higgs_scalar_branching_fraction(
    *,
    light_quark: str,
    scalar_left: complex = 0.0j,
    scalar_right: complex = 0.0j,
    inputs: TopFCNCSMInputs | None = None,
) -> TopFCNCBranchingResult:
    """Evaluate ``BR(t -> q h)`` from scalar FCNC couplings."""
    p = default_sm_inputs() if inputs is None else inputs
    width = top_higgs_scalar_partial_width(scalar_left, scalar_right, inputs=p)
    rho = (p.m_h_gev / p.m_top_gev) ** 2
    return _branching_result(
        process=f"t -> {light_quark} h",
        light_quark=light_quark,
        coupling_kind="Higgs scalar",
        partial_width_gev=width,
        coupling_left=scalar_left,
        coupling_right=scalar_right,
        inputs=p,
        convention=TOP_FCNC_HIGGS_SCALAR_CONVENTION,
        extra_diagnostics={
            "m_h_gev": float(p.m_h_gev),
            "higgs_phase_space": float((1.0 - rho) ** 2),
            "scalar_left": complex(scalar_left),
            "scalar_right": complex(scalar_right),
        },
    )


def compute_top_z_fcnc_proxy(
    source: QuarkMassBasisCouplings,
    *,
    light_quark: str,
    light_up_index: int,
    m_kk_gev: float | None = None,
    inputs: TopFCNCSMInputs | None = None,
) -> TopZFCNCProxyCouplings:
    """Match mass-basis up-sector couplings onto the v1 top-Z proxy.

    The flavor index convention follows the existing up-sector Delta-F=2 path:
    ``(0, 2)`` is ``u-t`` and ``(1, 2)`` is ``c-t`` in the coefficient of
    ``qbar gamma_mu t``.
    """
    if light_up_index not in (0, 1):
        raise ValueError("light_up_index must be 0 (u) or 1 (c)")
    p = default_sm_inputs() if inputs is None else inputs
    resolved_m_kk = _positive_float(
        getattr(source, "M_KK") if m_kk_gev is None else m_kk_gev,
        "m_kk_gev",
    )
    left_qt = _matrix_entry(source, "left_up", light_up_index, 2)
    right_qt = _matrix_entry(source, "right_up", light_up_index, 2)
    g_s = _positive_float(getattr(source, "g_s", 1.0), "g_s")
    left_overlap = left_qt / g_s
    right_overlap = right_qt / g_s
    suppression = (p.m_z_gev / resolved_m_kk) ** 2
    vector_left = left_overlap * suppression
    vector_right = right_overlap * suppression
    return TopZFCNCProxyCouplings(
        model_label=TOP_FCNC_MODEL_V1,
        operator_convention=TOP_FCNC_Z_VECTOR_CONVENTION,
        matching_assumption=TOP_FCNC_RS_Z_MATCHING_ASSUMPTION_V1,
        light_quark=light_quark,
        light_up_index=int(light_up_index),
        M_KK=float(resolved_m_kk),
        matching_scale=float(resolved_m_kk),
        left_qt_coupling=complex(left_qt),
        right_qt_coupling=complex(right_qt),
        left_qt_overlap=complex(left_overlap),
        right_qt_overlap=complex(right_overlap),
        z_mixing_suppression=float(suppression),
        vector_left=complex(vector_left),
        vector_right=complex(vector_right),
    )


def evaluate_t_to_q_z(
    source: QuarkMassBasisCouplings | TopZFCNCProxyCouplings | None = None,
    *,
    light_quark: str,
    light_up_index: int,
    m_kk_gev: float | None = None,
    inputs: TopFCNCSMInputs | None = None,
) -> TopFCNCBranchingResult:
    """Evaluate ``BR(t -> q Z)`` in the SM-zero limit or with the v1 proxy."""
    p = default_sm_inputs() if inputs is None else inputs
    if source is None:
        return z_vector_branching_fraction(
            light_quark=light_quark,
            vector_left=0.0j,
            vector_right=0.0j,
            inputs=p,
        )
    if isinstance(source, TopZFCNCProxyCouplings):
        proxy = source
    else:
        proxy = compute_top_z_fcnc_proxy(
            source,
            light_quark=light_quark,
            light_up_index=light_up_index,
            m_kk_gev=m_kk_gev,
            inputs=p,
        )
    return z_vector_branching_fraction(
        light_quark=light_quark,
        vector_left=proxy.vector_left,
        vector_right=proxy.vector_right,
        inputs=p,
        proxy=proxy,
    )


__all__ = [
    "TOP_FCNC_MODEL_V1",
    "TOP_FCNC_INPUT_BUNDLE_V1",
    "TOP_FCNC_Z_VECTOR_CONVENTION",
    "TOP_FCNC_PHOTON_DIPOLE_CONVENTION",
    "TOP_FCNC_GLUON_DIPOLE_CONVENTION",
    "TOP_FCNC_HIGGS_SCALAR_CONVENTION",
    "TOP_FCNC_RS_Z_MATCHING_ASSUMPTION_V1",
    "TopFCNCSMInputs",
    "TopZFCNCProxyCouplings",
    "TopFCNCBranchingResult",
    "default_sm_inputs",
    "weak_neutral_current_coupling",
    "top_z_vector_partial_width",
    "top_photon_dipole_partial_width",
    "top_gluon_dipole_partial_width",
    "top_higgs_scalar_partial_width",
    "z_vector_branching_fraction",
    "photon_dipole_branching_fraction",
    "gluon_dipole_branching_fraction",
    "higgs_scalar_branching_fraction",
    "compute_top_z_fcnc_proxy",
    "evaluate_t_to_q_z",
]
