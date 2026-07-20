"""LFV ``D0 -> e mu`` wrapper over the rare-charm dilepton core.

The shared :mod:`quarkConstraints.rare_charm_dilepton` module owns the
``c -> u l l`` CKM normalization, D0 inputs, and documented Z/KK-penguin
quark matching used by C004/C005/C007.  This module reuses that machinery and
adds only the charged-lepton-flavor-violating glue needed for
``D0 -> e+- mu-+``.

The operator convention is the same rare-charm Hamiltonian,

    H_eff = -4 G_F/sqrt(2) lambda_b alpha/(4 pi)
            [ C9 O9 + C10 O10 + C9' O9' + C10' O10' ],

with LFV lepton bilinears ``e bar gamma_mu mu``.  For unequal leptons the
two-body rate uses the standard pseudoscalar phase-space form

    BR(P -> l_i l_j) = tau G_F^2 alpha^2 f_P^2 m_P^3 |lambda|^2
        / (64 pi^3) sqrt(lambda_ij)
        * { [1-(r_i+r_j)^2] |(r_i-r_j)(C9-C9')|^2
          + [1-(r_i-r_j)^2] |(r_i+r_j)(C10-C10')|^2 }.

For equal lepton masses and ``C9=0`` this reduces to the same-flavor formula
implemented in ``rare_charm_dilepton``.

NEEDS-HUMAN-PHYSICS: a rigorous RS prediction requires the off-diagonal
charged-lepton neutral-current coupling after EW KK/Z/Z' mixing and charged
lepton mass-basis rotations.  That object is not present on ParameterPoint.
The v1 proxy accepts an explicit e-mu overlap spurion and maps it to a
Z-like LFV lepton coupling, while the quark-side matching is reused from the
rare-charm dilepton core.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Mapping

import numpy as np

from .couplings import QuarkMassBasisCouplings
from .rare_charm_dilepton import (
    RARE_CHARM_DILEPTON_INPUT_BUNDLE_V1,
    RARE_CHARM_DILEPTON_MODEL_V1,
    RARE_CHARM_DILEPTON_OPERATOR_CONVENTION,
    RARE_CHARM_DILEPTON_RS_MATCHING_ASSUMPTION_V1,
    RareCharmDileptonSMInputs,
    RareCharmDileptonWilsonCoefficients,
    compute_rare_charm_dilepton_wilsons,
    default_sm_inputs,
)

RARE_CHARM_LFV_DILEPTON_MODEL_V1 = "rare_charm_lfv_d0_emu_proxy_v1"
RARE_CHARM_LFV_DILEPTON_OPERATOR_CONVENTION = (
    RARE_CHARM_DILEPTON_OPERATOR_CONVENTION
    + "; LFV extension uses e-mu lepton bilinears and the charge-summed "
    "D0 -> e+- mu-+ branching fraction"
)
RARE_CHARM_LFV_DILEPTON_PROXY_V1 = (
    "NEEDS-HUMAN-PHYSICS: off-diagonal e-mu lepton neutral-current couplings "
    "are not available on ParameterPoint; caller-supplied lepton overlap "
    "spurions are mapped to Z-like LFV couplings delta_l_L/R = g_Z * "
    "overlap_emu, standing in for the missing RS EW KK/Z/Z' and charged-"
    "lepton mass-basis matching."
)
_DEFAULT_CHARGE_STATE_FACTOR = 2.0


@dataclass(frozen=True)
class RareCharmLFVLeptonProxyInput:
    """Explicit e-mu lepton-overlap proxy for LFV neutral currents."""

    left_emu_overlap: complex
    right_emu_overlap: complex
    m_kk_gev: float
    source: str = "caller-supplied rare-charm e-mu lepton proxy"


@dataclass(frozen=True)
class RareCharmLFVLeptonCouplingProxy:
    """Mapped Z-like e-mu lepton coupling used by the LFV charm proxy."""

    model_label: str
    matching_assumption: str
    M_KK: float
    matching_scale: float
    left_emu_overlap: complex
    right_emu_overlap: complex
    lepton_left_delta_emu: complex
    lepton_right_delta_emu: complex
    lepton_vector_delta_emu: complex
    lepton_axial_delta_emu: complex
    source: str
    diagnostics: Mapping[str, Any]


@dataclass(frozen=True)
class RareCharmLFVWilsonCoefficients:
    """LFV ``c -> u e mu`` Wilson proxy from quark and lepton spurions."""

    model_label: str
    base_model_label: str
    operator_convention: str
    matching_assumption: str
    transition_key: str
    M_KK: float
    matching_scale: float
    lambda_b: complex
    left_uc_coupling: complex
    right_uc_coupling: complex
    left_uc_overlap: complex
    right_uc_overlap: complex
    left_quark_delta: complex
    right_quark_delta: complex
    lepton_left_delta_emu: complex
    lepton_right_delta_emu: complex
    lepton_vector_delta_emu: complex
    lepton_axial_delta_emu: complex
    c9_lfv_np: complex
    c10_lfv_np: complex
    c9p_lfv_np: complex
    c10p_lfv_np: complex
    base_same_flavor_wilsons: RareCharmDileptonWilsonCoefficients

    @property
    def c9_leptonic_lfv_np(self) -> complex:
        """Combination entering the vector part: ``C9 - C9'``."""

        return complex(self.c9_lfv_np - self.c9p_lfv_np)

    @property
    def c10_leptonic_lfv_np(self) -> complex:
        """Combination entering the axial part: ``C10 - C10'``."""

        return complex(self.c10_lfv_np - self.c10p_lfv_np)

    @property
    def wilsons(self) -> Mapping[str, complex]:
        return {
            "C9_LFV_NP": complex(self.c9_lfv_np),
            "C10_LFV_NP": complex(self.c10_lfv_np),
            "C9p_LFV_NP": complex(self.c9p_lfv_np),
            "C10p_LFV_NP": complex(self.c10p_lfv_np),
        }


@dataclass(frozen=True)
class RareCharmLFVBranchingResult:
    """Short-distance ``D0 -> e+- mu-+`` branching-fraction prediction."""

    model_label: str
    input_bundle: str
    transition_key: str
    branching_fraction: float
    sm_branching_fraction: float
    np_shift_branching_fraction: float
    c9_lfv_combination: complex
    c10_lfv_combination: complex
    lambda_b: complex
    electron_mass_gev: float
    muon_mass_gev: float
    phase_space_lambda_sqrt: float
    charge_state_factor: float
    wilsons: RareCharmLFVWilsonCoefficients | None = None
    diagnostics: Mapping[str, Any] = field(default_factory=dict)


def rare_charm_lfv_proxy_input(
    left_emu_overlap: complex,
    right_emu_overlap: complex,
    m_kk_gev: float,
    *,
    source: str = "caller-supplied rare-charm e-mu lepton proxy",
) -> RareCharmLFVLeptonProxyInput:
    """Build a shape-checked LFV lepton proxy input."""

    return RareCharmLFVLeptonProxyInput(
        left_emu_overlap=_finite_complex(left_emu_overlap, "left_emu_overlap"),
        right_emu_overlap=_finite_complex(right_emu_overlap, "right_emu_overlap"),
        m_kk_gev=_positive_float(m_kk_gev, "m_kk_gev"),
        source=str(source),
    )


def rare_charm_lfv_lepton_coupling_proxy(
    source: Any,
    *,
    m_kk_gev: float | None = None,
    inputs: RareCharmDileptonSMInputs | None = None,
) -> RareCharmLFVLeptonCouplingProxy:
    """Return the documented Z-like e-mu lepton-coupling proxy."""

    p = default_sm_inputs() if inputs is None else inputs
    resolved_m_kk = None if m_kk_gev is None else _positive_float(m_kk_gev, "m_kk_gev")
    proxy_input = _coerce_proxy_input(source, m_kk_gev=resolved_m_kk)
    if resolved_m_kk is None:
        resolved_m_kk = proxy_input.m_kk_gev
    g_weak = math.sqrt(4.0 * math.pi * p.alpha_em_mz / p.sin2_theta_w)
    g_z = g_weak / math.sqrt(1.0 - p.sin2_theta_w)
    left = complex(g_z * proxy_input.left_emu_overlap)
    right = complex(g_z * proxy_input.right_emu_overlap)
    vector = complex(left + right)
    axial = complex(right - left)
    return RareCharmLFVLeptonCouplingProxy(
        model_label=RARE_CHARM_LFV_DILEPTON_MODEL_V1,
        matching_assumption=RARE_CHARM_LFV_DILEPTON_PROXY_V1,
        M_KK=float(resolved_m_kk),
        matching_scale=float(resolved_m_kk),
        left_emu_overlap=complex(proxy_input.left_emu_overlap),
        right_emu_overlap=complex(proxy_input.right_emu_overlap),
        lepton_left_delta_emu=left,
        lepton_right_delta_emu=right,
        lepton_vector_delta_emu=vector,
        lepton_axial_delta_emu=axial,
        source=proxy_input.source,
        diagnostics={
            "m_kk_gev": float(resolved_m_kk),
            "matching_scale_gev": float(resolved_m_kk),
            "g_z": float(g_z),
            "left_emu_overlap": complex(proxy_input.left_emu_overlap),
            "right_emu_overlap": complex(proxy_input.right_emu_overlap),
            "lepton_left_delta_emu": left,
            "lepton_right_delta_emu": right,
            "lepton_vector_delta_emu": vector,
            "lepton_axial_delta_emu": axial,
            "proxy_source": proxy_input.source,
            "matching_assumption": RARE_CHARM_LFV_DILEPTON_PROXY_V1,
        },
    )


def compute_rare_charm_lfv_wilsons(
    quark_source: QuarkMassBasisCouplings,
    lepton_source: Any,
    *,
    transition: str = "c_u",
    m_kk_gev: float | None = None,
    inputs: RareCharmDileptonSMInputs | None = None,
) -> RareCharmLFVWilsonCoefficients:
    """Match quark and LFV lepton spurions onto ``c -> u e mu`` Wilsons."""

    if transition != "c_u":
        raise ValueError(f"unsupported rare-charm LFV transition {transition!r}")
    p = default_sm_inputs() if inputs is None else inputs
    base = compute_rare_charm_dilepton_wilsons(
        quark_source,
        transition=transition,
        m_kk_gev=m_kk_gev,
        inputs=p,
    )
    proxy = rare_charm_lfv_lepton_coupling_proxy(
        lepton_source,
        m_kk_gev=base.M_KK,
        inputs=p,
    )
    vector_scale = _safe_ratio(
        proxy.lepton_vector_delta_emu,
        base.lepton_vector_delta,
        "lepton_vector_delta",
    )
    axial_scale = _safe_ratio(
        proxy.lepton_axial_delta_emu,
        base.lepton_axial_delta,
        "lepton_axial_delta",
    )
    return RareCharmLFVWilsonCoefficients(
        model_label=RARE_CHARM_LFV_DILEPTON_MODEL_V1,
        base_model_label=RARE_CHARM_DILEPTON_MODEL_V1,
        operator_convention=RARE_CHARM_LFV_DILEPTON_OPERATOR_CONVENTION,
        matching_assumption=RARE_CHARM_LFV_DILEPTON_PROXY_V1,
        transition_key=transition,
        M_KK=float(base.M_KK),
        matching_scale=float(base.matching_scale),
        lambda_b=complex(base.lambda_b),
        left_uc_coupling=complex(base.left_uc_coupling),
        right_uc_coupling=complex(base.right_uc_coupling),
        left_uc_overlap=complex(base.left_uc_overlap),
        right_uc_overlap=complex(base.right_uc_overlap),
        left_quark_delta=complex(base.left_quark_delta),
        right_quark_delta=complex(base.right_quark_delta),
        lepton_left_delta_emu=complex(proxy.lepton_left_delta_emu),
        lepton_right_delta_emu=complex(proxy.lepton_right_delta_emu),
        lepton_vector_delta_emu=complex(proxy.lepton_vector_delta_emu),
        lepton_axial_delta_emu=complex(proxy.lepton_axial_delta_emu),
        c9_lfv_np=complex(base.c9_np * vector_scale),
        c10_lfv_np=complex(base.c10_np * axial_scale),
        c9p_lfv_np=complex(base.c9p_np * vector_scale),
        c10p_lfv_np=complex(base.c10p_np * axial_scale),
        base_same_flavor_wilsons=base,
    )


def rare_charm_lfv_sm_branching_fraction(
    inputs: RareCharmDileptonSMInputs | None = None,
) -> RareCharmLFVBranchingResult:
    """Return the catalog SM-limit ``D0 -> e mu`` rate, zero for LFV."""

    p = default_sm_inputs() if inputs is None else inputs
    return _branching_from_wilsons(None, inputs=p)


def evaluate_d0_to_emu(
    quark_source: QuarkMassBasisCouplings | RareCharmLFVWilsonCoefficients,
    lepton_source: Any | None = None,
    *,
    m_kk_gev: float | None = None,
    inputs: RareCharmDileptonSMInputs | None = None,
    charge_state_factor: float = _DEFAULT_CHARGE_STATE_FACTOR,
) -> RareCharmLFVBranchingResult:
    """Evaluate the charge-summed short-distance ``BR(D0 -> e+- mu-+)``."""

    p = default_sm_inputs() if inputs is None else inputs
    if isinstance(quark_source, RareCharmLFVWilsonCoefficients):
        wilsons = quark_source
    else:
        if lepton_source is None:
            raise TypeError("lepton_source is required for rare-charm LFV matching")
        wilsons = compute_rare_charm_lfv_wilsons(
            quark_source,
            lepton_source,
            m_kk_gev=m_kk_gev,
            inputs=p,
        )
    return _branching_from_wilsons(
        wilsons,
        inputs=p,
        charge_state_factor=charge_state_factor,
    )


def _branching_from_wilsons(
    wilsons: RareCharmLFVWilsonCoefficients | None,
    *,
    inputs: RareCharmDileptonSMInputs,
    charge_state_factor: float = _DEFAULT_CHARGE_STATE_FACTOR,
) -> RareCharmLFVBranchingResult:
    charge_factor = _positive_float(charge_state_factor, "charge_state_factor")
    meson = inputs.d0
    electron = inputs.electron
    muon = inputs.muon
    r_e = electron.mass_gev / meson.meson_mass_gev
    r_mu = muon.mass_gev / meson.meson_mass_gev
    lambda_phase = _kallen_dimensionless(r_e * r_e, r_mu * r_mu)
    if lambda_phase <= 0.0:
        raise ValueError("D0 -> e mu phase space is closed")
    sqrt_lambda = math.sqrt(lambda_phase)

    if wilsons is None:
        c9_combo = 0.0j
        c10_combo = 0.0j
        lambda_b = 0.0j
    else:
        c9_combo = wilsons.c9_leptonic_lfv_np
        c10_combo = wilsons.c10_leptonic_lfv_np
        lambda_b = wilsons.lambda_b

    mass_sum = r_e + r_mu
    mass_diff = r_mu - r_e
    vector_term = (1.0 - mass_sum * mass_sum) * abs(mass_diff * c9_combo) ** 2
    axial_term = (1.0 - mass_diff * mass_diff) * abs(mass_sum * c10_combo) ** 2
    tau = meson.lifetime_ps * 1.0e-12 / inputs.hbar_gev_s
    prefactor = (
        tau
        * inputs.gf_gev_minus2**2
        * inputs.alpha_em_mz**2
        / (64.0 * math.pi**3)
        * meson.decay_constant_gev**2
        * meson.meson_mass_gev**3
        * abs(lambda_b) ** 2
        * sqrt_lambda
    )
    branching_one_charge = float(prefactor * (vector_term + axial_term))
    branching = float(charge_factor * branching_one_charge)
    diagnostics: dict[str, Any] = {
        "base_model_label": RARE_CHARM_DILEPTON_MODEL_V1,
        "base_matching_assumption": RARE_CHARM_DILEPTON_RS_MATCHING_ASSUMPTION_V1,
        "matching_assumption": RARE_CHARM_LFV_DILEPTON_PROXY_V1,
        "operator_convention": RARE_CHARM_LFV_DILEPTON_OPERATOR_CONVENTION,
        "input_bundle": RARE_CHARM_DILEPTON_INPUT_BUNDLE_V1,
        "sm_branching_fraction": 0.0,
        "sm_lfv_policy": "Charged-LFV D0 -> e mu has zero SM rate for catalog purposes.",
        "electron_mass_gev": float(electron.mass_gev),
        "muon_mass_gev": float(muon.mass_gev),
        "r_e": float(r_e),
        "r_mu": float(r_mu),
        "phase_space_lambda": float(lambda_phase),
        "phase_space_lambda_sqrt": float(sqrt_lambda),
        "mass_sum_ratio": float(mass_sum),
        "mass_difference_ratio": float(mass_diff),
        "charge_state_factor": float(charge_factor),
        "branching_fraction_one_charge": float(branching_one_charge),
        "vector_term": float(vector_term),
        "axial_term": float(axial_term),
        "c9_lfv_combination": complex(c9_combo),
        "c10_lfv_combination": complex(c10_combo),
        "lambda_b": complex(lambda_b),
    }
    if wilsons is not None:
        diagnostics.update(
            {
                "m_kk_gev": float(wilsons.M_KK),
                "matching_scale_gev": float(wilsons.matching_scale),
                "left_uc_coupling": complex(wilsons.left_uc_coupling),
                "right_uc_coupling": complex(wilsons.right_uc_coupling),
                "left_uc_overlap": complex(wilsons.left_uc_overlap),
                "right_uc_overlap": complex(wilsons.right_uc_overlap),
                "left_quark_delta": complex(wilsons.left_quark_delta),
                "right_quark_delta": complex(wilsons.right_quark_delta),
                "lepton_left_delta_emu": complex(wilsons.lepton_left_delta_emu),
                "lepton_right_delta_emu": complex(wilsons.lepton_right_delta_emu),
                "lepton_vector_delta_emu": complex(wilsons.lepton_vector_delta_emu),
                "lepton_axial_delta_emu": complex(wilsons.lepton_axial_delta_emu),
                "c9_lfv_np": complex(wilsons.c9_lfv_np),
                "c10_lfv_np": complex(wilsons.c10_lfv_np),
                "c9p_lfv_np": complex(wilsons.c9p_lfv_np),
                "c10p_lfv_np": complex(wilsons.c10p_lfv_np),
            }
        )

    return RareCharmLFVBranchingResult(
        model_label=RARE_CHARM_LFV_DILEPTON_MODEL_V1,
        input_bundle=RARE_CHARM_DILEPTON_INPUT_BUNDLE_V1,
        transition_key="c_u",
        branching_fraction=branching,
        sm_branching_fraction=0.0,
        np_shift_branching_fraction=branching,
        c9_lfv_combination=complex(c9_combo),
        c10_lfv_combination=complex(c10_combo),
        lambda_b=complex(lambda_b),
        electron_mass_gev=float(electron.mass_gev),
        muon_mass_gev=float(muon.mass_gev),
        phase_space_lambda_sqrt=float(sqrt_lambda),
        charge_state_factor=float(charge_factor),
        wilsons=wilsons,
        diagnostics=diagnostics,
    )


def _coerce_proxy_input(
    value: Any,
    *,
    m_kk_gev: float | None,
) -> RareCharmLFVLeptonProxyInput:
    if isinstance(value, RareCharmLFVLeptonProxyInput):
        return rare_charm_lfv_proxy_input(
            value.left_emu_overlap,
            value.right_emu_overlap,
            value.m_kk_gev if m_kk_gev is None else m_kk_gev,
            source=value.source,
        )
    if isinstance(value, Mapping):
        return _proxy_from_mapping(value, m_kk_gev=m_kk_gev)
    return _proxy_from_object(value, m_kk_gev=m_kk_gev)


def _proxy_from_mapping(
    mapping: Mapping[str, Any],
    *,
    m_kk_gev: float | None,
) -> RareCharmLFVLeptonProxyInput:
    left, right = _overlaps_from_mapping(mapping)
    return rare_charm_lfv_proxy_input(
        left,
        right,
        _resolve_m_kk_from_mapping(mapping, m_kk_gev),
        source=str(mapping.get("source", "mapping rare-charm e-mu lepton proxy")),
    )


def _proxy_from_object(
    value: Any,
    *,
    m_kk_gev: float | None,
) -> RareCharmLFVLeptonProxyInput:
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
                "rare-charm LFV lepton proxy must provide left/right e-mu "
                "overlap entries or charged-lepton overlap matrices"
            )
        left = 0.0j if left_matrix is None else _matrix_offdiag(left_matrix, "left", 0, 1)
        right = 0.0j if right_matrix is None else _matrix_offdiag(right_matrix, "right", 0, 1)
    return rare_charm_lfv_proxy_input(
        0.0j if left is None else left,
        0.0j if right is None else right,
        _resolve_m_kk_from_object(value, m_kk_gev),
        source=str(getattr(value, "source", "object rare-charm e-mu lepton proxy")),
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
        return _positive_float(override, "m_kk_gev")
    value = _first_present_key(mapping, ("m_kk_gev", "M_KK_gev", "M_KK"))
    if value is None:
        raise KeyError("mapping must provide m_kk_gev, M_KK_gev, or M_KK")
    return _positive_float(value, "m_kk_gev")


def _resolve_m_kk_from_object(value: Any, override: float | None) -> float:
    if override is not None:
        return _positive_float(override, "m_kk_gev")
    for name in ("m_kk_gev", "M_KK_gev", "M_KK"):
        if hasattr(value, name):
            return _positive_float(getattr(value, name), "m_kk_gev")
    raise AttributeError("rare-charm LFV lepton proxy object must provide m_kk_gev or M_KK")


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


def _finite_complex(value: Any, name: str) -> complex:
    out = complex(value)
    if not math.isfinite(out.real) or not math.isfinite(out.imag):
        raise ValueError(f"{name} must be finite")
    return out


def _positive_float(value: Any, name: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0.0:
        raise ValueError(f"{name} must be positive and finite")
    return number


def _safe_ratio(numerator: complex, denominator: float, name: str) -> complex:
    denom = float(denominator)
    if not math.isfinite(denom) or denom == 0.0:
        raise ValueError(f"{name} denominator must be finite and non-zero")
    return complex(numerator / denom)


def _kallen_dimensionless(x: float, y: float) -> float:
    return float(1.0 + x * x + y * y - 2.0 * (x + y + x * y))


__all__ = [
    "RARE_CHARM_LFV_DILEPTON_MODEL_V1",
    "RARE_CHARM_LFV_DILEPTON_OPERATOR_CONVENTION",
    "RARE_CHARM_LFV_DILEPTON_PROXY_V1",
    "RareCharmLFVLeptonProxyInput",
    "RareCharmLFVLeptonCouplingProxy",
    "RareCharmLFVWilsonCoefficients",
    "RareCharmLFVBranchingResult",
    "rare_charm_lfv_proxy_input",
    "rare_charm_lfv_lepton_coupling_proxy",
    "compute_rare_charm_lfv_wilsons",
    "rare_charm_lfv_sm_branching_fraction",
    "evaluate_d0_to_emu",
]
