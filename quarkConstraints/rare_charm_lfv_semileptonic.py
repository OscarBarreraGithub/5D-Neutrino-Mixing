"""LFV ``D+ -> pi+ e mu`` wrapper over rare-charm semileptonic machinery.

The D-to-pi form factors, CKM normalization, and rare-charm Hamiltonian
convention are reused from :mod:`quarkConstraints.rare_charm_semileptonic`.
The charged-LFV ``c -> u e mu`` Wilson proxy is reused from
:mod:`quarkConstraints.rare_charm_lfv_dilepton`.

The rate implemented here is a full-kinematic-q2, short-distance vector/axial
proxy for one charge mode.  It evaluates the unequal-lepton phase space with a
direct spin-summed lepton-tensor contraction and uses the same normalization as
the C007 ``D+ -> pi+ mu+ mu-`` smooth short-distance formula.  Scalar/tensor
operators, resonance amplitudes, and LHCb window/acceptance effects are not
included.

NEEDS-HUMAN-PHYSICS: a rigorous RS prediction requires the off-diagonal
charged-lepton neutral-current coupling after EW KK/Z/Z' mixing and charged
lepton mass-basis rotations.  That object is not present on ParameterPoint.
The v1 proxy accepts an explicit e-mu overlap spurion and maps it to the same
Z-like LFV lepton coupling used by C006, while the D-to-pi form factor and
quark-side matching are reused from C007/C004.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from functools import lru_cache
from typing import Any, Mapping

import numpy as np

from .couplings import QuarkMassBasisCouplings
from .rare_charm_dilepton import (
    RARE_CHARM_DILEPTON_RS_MATCHING_ASSUMPTION_V1,
    ckm_factors,
)
from .rare_charm_lfv_dilepton import (
    RARE_CHARM_LFV_DILEPTON_OPERATOR_CONVENTION,
    RARE_CHARM_LFV_DILEPTON_PROXY_V1,
    RareCharmLFVWilsonCoefficients,
    compute_rare_charm_lfv_wilsons,
)
from .rare_charm_semileptonic import (
    RARE_CHARM_DTOPI_MUMU_FORM_FACTOR_MODEL_V1,
    RARE_CHARM_DTOPI_MUMU_INPUT_BUNDLE_V1,
    RARE_CHARM_DTOPI_MUMU_PARAMETRIZATION_CITATION,
    RARE_CHARM_DTOPI_MUMU_RESONANCE_LIMITATION_V1,
    RareCharmDToPiMuMuBranchingResult,
    RareCharmDToPiMuMuInputs,
    default_dtopi_mumu_inputs,
    dtopi_fplus,
    dtopi_fzero,
)

RARE_CHARM_DTOPI_EMU_MODEL_V1 = "rare_charm_dtopi_emu_lfv_form_factor_proxy_v1"
RARE_CHARM_DTOPI_EMU_OPERATOR_CONVENTION = (
    RARE_CHARM_LFV_DILEPTON_OPERATOR_CONVENTION
    + "; D->pi LFV semileptonic rate uses C9_LFV+C9p_LFV and "
    "C10_LFV+C10p_LFV with unequal e/mu phase space."
)
RARE_CHARM_DTOPI_EMU_PARAMETRIZATION_CITATION = (
    RARE_CHARM_DTOPI_MUMU_PARAMETRIZATION_CITATION
    + "; unequal-lepton e-mu phase space evaluated by direct vector/axial "
    "lepton-tensor contraction in the dilepton rest frame."
)
RARE_CHARM_DTOPI_EMU_PROXY_V1 = (
    "NEEDS-HUMAN-PHYSICS: D+ -> pi+ e mu v1 reuses the C007 D->pi "
    "form-factor normalization and the C006 e-mu LFV Z-like coupling proxy; "
    "scalar/tensor operators, resonance amplitudes, charge-mode-specific "
    "lepton matching, and LHCb acceptance windows are not available."
)
RARE_CHARM_DTOPI_EMU_Q2_TREATMENT_V1 = (
    "full_kinematic_q2_smooth_short_distance_lfv_proxy_no_lhcb_window_or_"
    "resonance_acceptance"
)
_ANGULAR_QUADRATURE_POINTS = 48
_METRIC = np.diag([1.0, -1.0, -1.0, -1.0])


@dataclass(frozen=True)
class RareCharmDToPiLFVBranchingResult:
    """Short-distance branching-ratio prediction for one ``D+ -> pi+ e mu`` mode."""

    model_label: str
    input_bundle: str
    transition_key: str
    charge_mode: str
    branching_fraction: float
    sm_branching_fraction: float
    np_shift_branching_fraction: float
    q2_min_gev2: float
    q2_max_gev2: float
    electron_mass_gev: float
    muon_mass_gev: float
    c9_lfv_semileptonic_np: complex
    c10_lfv_semileptonic_np: complex
    c9_lfv_np: complex
    c10_lfv_np: complex
    c9p_lfv_np: complex
    c10p_lfv_np: complex
    lambda_b: complex
    wilsons: RareCharmLFVWilsonCoefficients | None = None
    diagnostics: Mapping[str, Any] = field(default_factory=dict)


def default_dtopi_emu_inputs() -> RareCharmDToPiMuMuInputs:
    """Return the shared ``D+ -> pi+`` short-distance input bundle."""

    return default_dtopi_mumu_inputs()


def dtopi_emu_q2_range(
    inputs: RareCharmDToPiMuMuInputs | None = None,
) -> tuple[float, float]:
    """Return the kinematic ``q2`` range for ``D+ -> pi+ e mu``."""

    p = default_dtopi_emu_inputs() if inputs is None else inputs
    m_e = p.rare_charm.electron.mass_gev
    m_mu = p.rare_charm.muon.mass_gev
    q2_min = (m_e + m_mu) ** 2
    q2_max = (p.dplus_mass_gev - p.piplus_mass_gev) ** 2
    if q2_max <= q2_min:
        raise ValueError("D+ -> pi+ e mu phase space is closed")
    return float(q2_min), float(q2_max)


def dtopi_emu_differential_branching_fraction(
    q2_gev2: float,
    *,
    c9_lfv_semileptonic: complex,
    c10_lfv_semileptonic: complex,
    inputs: RareCharmDToPiMuMuInputs | None = None,
    charge_mode: str = "eplus_muminus",
) -> float:
    """Evaluate ``dBR_SD(D+ -> pi+ e mu) / dq2`` for one LFV charge mode."""

    p = default_dtopi_emu_inputs() if inputs is None else inputs
    q2 = float(q2_gev2)
    q2_min, q2_max = dtopi_emu_q2_range(p)
    if q2 < q2_min or q2 > q2_max:
        raise ValueError(f"q2={q2} outside D+ -> pi+ e mu physical range")

    if charge_mode == "eplus_muminus":
        m_minus = p.rare_charm.muon.mass_gev
        m_plus = p.rare_charm.electron.mass_gev
    elif charge_mode == "eminus_muplus":
        m_minus = p.rare_charm.electron.mass_gev
        m_plus = p.rare_charm.muon.mass_gev
    else:
        raise ValueError(f"unsupported D+ -> pi+ e mu charge mode {charge_mode!r}")

    m_d = p.dplus_mass_gev
    m_pi = p.piplus_mass_gev
    m_d2 = m_d * m_d
    m_pi2 = m_pi * m_pi
    lam_had = max(0.0, _kallen(m_d2, m_pi2, q2))
    lam_lep = max(0.0, _kallen(q2, m_minus * m_minus, m_plus * m_plus))
    if lam_had <= 0.0 or lam_lep <= 0.0:
        return 0.0

    factors = ckm_factors("c_u", p.rare_charm)
    gamma0 = (
        p.rare_charm.gf_gev_minus2**2
        * p.rare_charm.alpha_em_mz**2
        * abs(factors.lambda_b) ** 2
        / ((4.0 * math.pi) ** 5 * m_d**3)
    )
    tau = p.dplus_lifetime_ps * 1.0e-12 / p.rare_charm.hbar_gev_s
    beta_lfv = math.sqrt(lam_lep) / q2
    angular = _angular_tensor_integral(
        q2,
        m_minus=m_minus,
        m_plus=m_plus,
        c9_lfv_semileptonic=complex(c9_lfv_semileptonic),
        c10_lfv_semileptonic=complex(c10_lfv_semileptonic),
        inputs=p,
    )
    rate = tau * gamma0 * math.sqrt(lam_had) * beta_lfv * angular / 4.0
    if rate < 0.0 and abs(rate) < 1.0e-30:
        return 0.0
    if rate < 0.0:
        raise ValueError(f"negative dBR/dq2={rate} at q2={q2}")
    return float(rate)


def dtopi_emu_sm(
    inputs: RareCharmDToPiMuMuInputs | None = None,
) -> RareCharmDToPiLFVBranchingResult:
    """Return the catalog SM-limit ``D+ -> pi+ e mu`` rate, zero for LFV."""

    return evaluate_dplus_to_piplus_emu(None, inputs=inputs)


def evaluate_dplus_to_piplus_emu(
    quark_source: QuarkMassBasisCouplings | RareCharmLFVWilsonCoefficients | None = None,
    lepton_source: Any | None = None,
    *,
    m_kk_gev: float | None = None,
    inputs: RareCharmDToPiMuMuInputs | None = None,
    charge_mode: str = "eplus_muminus",
) -> RareCharmDToPiLFVBranchingResult:
    """Evaluate smooth ``BR_SD(D+ -> pi+ e mu)`` with the v1 LFV proxy."""

    p = default_dtopi_emu_inputs() if inputs is None else inputs
    if quark_source is None:
        wilsons = None
        c9_np = c10_np = c9p_np = c10p_np = 0.0j
        lambda_b = 0.0j
    elif isinstance(quark_source, RareCharmLFVWilsonCoefficients):
        wilsons = quark_source
        c9_np = wilsons.c9_lfv_np
        c10_np = wilsons.c10_lfv_np
        c9p_np = wilsons.c9p_lfv_np
        c10p_np = wilsons.c10p_lfv_np
        lambda_b = wilsons.lambda_b
    else:
        if lepton_source is None:
            raise TypeError("lepton_source is required for D+ -> pi+ e mu LFV matching")
        wilsons = compute_rare_charm_lfv_wilsons(
            quark_source,
            lepton_source,
            transition="c_u",
            m_kk_gev=m_kk_gev,
            inputs=p.rare_charm,
        )
        c9_np = wilsons.c9_lfv_np
        c10_np = wilsons.c10_lfv_np
        c9p_np = wilsons.c9p_lfv_np
        c10p_np = wilsons.c10p_lfv_np
        lambda_b = wilsons.lambda_b

    c9_semileptonic = complex(c9_np + c9p_np)
    c10_semileptonic = complex(c10_np + c10p_np)
    branching = _integrate_branching_fraction(
        c9_lfv_semileptonic=c9_semileptonic,
        c10_lfv_semileptonic=c10_semileptonic,
        inputs=p,
        charge_mode=charge_mode,
    )
    q2_min, q2_max = dtopi_emu_q2_range(p)
    diagnostics: dict[str, Any] = {
        "base_dtopi_input_bundle": RARE_CHARM_DTOPI_MUMU_INPUT_BUNDLE_V1,
        "base_form_factor_model": RARE_CHARM_DTOPI_MUMU_FORM_FACTOR_MODEL_V1,
        "base_matching_assumption": RARE_CHARM_DILEPTON_RS_MATCHING_ASSUMPTION_V1,
        "matching_assumption": RARE_CHARM_DTOPI_EMU_PROXY_V1,
        "lfv_lepton_matching_assumption": RARE_CHARM_LFV_DILEPTON_PROXY_V1,
        "operator_convention": RARE_CHARM_DTOPI_EMU_OPERATOR_CONVENTION,
        "parametrization_citation": RARE_CHARM_DTOPI_EMU_PARAMETRIZATION_CITATION,
        "resonance_limitation": RARE_CHARM_DTOPI_MUMU_RESONANCE_LIMITATION_V1,
        "q2_treatment": RARE_CHARM_DTOPI_EMU_Q2_TREATMENT_V1,
        "short_distance_only": True,
        "sm_branching_fraction": 0.0,
        "sm_lfv_policy": "Charged-LFV D+ -> pi+ e mu has zero SM rate for catalog purposes.",
        "charge_mode": charge_mode,
        "charge_mode_specific_lepton_matching_available": False,
        "charge_mode_prediction_policy": (
            "The v1 e-mu proxy is not orientation-specific; the same one-mode "
            "prediction is compared to both C008 charge-mode limits."
        ),
        "semileptonic_primed_combination_is_plus": True,
        "scalar_tensor_operators_included": False,
        "resonance_amplitudes_included": False,
        "lhcb_window_acceptance_applied": False,
        "q2_min_gev2": float(q2_min),
        "q2_max_gev2": float(q2_max),
        "quadrature_points": float(p.quadrature_points),
        "angular_quadrature_points": float(_ANGULAR_QUADRATURE_POINTS),
        "dplus_mass_gev": float(p.dplus_mass_gev),
        "piplus_mass_gev": float(p.piplus_mass_gev),
        "dplus_lifetime_ps": float(p.dplus_lifetime_ps),
        "electron_mass_gev": float(p.rare_charm.electron.mass_gev),
        "muon_mass_gev": float(p.rare_charm.muon.mass_gev),
        "lambda_b": complex(lambda_b),
        "c9_lfv_semileptonic_np": complex(c9_semileptonic),
        "c10_lfv_semileptonic_np": complex(c10_semileptonic),
        "form_factor_model": p.form_factor.model_label,
        "form_factor_source": p.form_factor.source,
        "fplus_0": float(p.form_factor.fplus_0),
        "fplus_shape_a": float(p.form_factor.fplus_shape_a),
        "fzero_shape_b": float(p.form_factor.fzero_shape_b),
        "form_factor_pole_mass_gev": float(p.form_factor.pole_mass_gev),
        "constants_citation": p.constants_citation,
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

    return RareCharmDToPiLFVBranchingResult(
        model_label=RARE_CHARM_DTOPI_EMU_MODEL_V1,
        input_bundle=p.input_bundle,
        transition_key="c_u",
        charge_mode=charge_mode,
        branching_fraction=float(branching),
        sm_branching_fraction=0.0,
        np_shift_branching_fraction=float(branching),
        q2_min_gev2=float(q2_min),
        q2_max_gev2=float(q2_max),
        electron_mass_gev=float(p.rare_charm.electron.mass_gev),
        muon_mass_gev=float(p.rare_charm.muon.mass_gev),
        c9_lfv_semileptonic_np=complex(c9_semileptonic),
        c10_lfv_semileptonic_np=complex(c10_semileptonic),
        c9_lfv_np=complex(c9_np),
        c10_lfv_np=complex(c10_np),
        c9p_lfv_np=complex(c9p_np),
        c10p_lfv_np=complex(c10p_np),
        lambda_b=complex(lambda_b),
        wilsons=wilsons,
        diagnostics=diagnostics,
    )


def _integrate_branching_fraction(
    *,
    c9_lfv_semileptonic: complex,
    c10_lfv_semileptonic: complex,
    inputs: RareCharmDToPiMuMuInputs,
    charge_mode: str,
) -> float:
    q2_min, q2_max = dtopi_emu_q2_range(inputs)
    nodes, weights = np.polynomial.legendre.leggauss(int(inputs.quadrature_points))
    mid = 0.5 * (q2_max + q2_min)
    half_width = 0.5 * (q2_max - q2_min)
    total = 0.0
    for node, weight in zip(nodes, weights):
        q2 = mid + half_width * float(node)
        total += float(weight) * dtopi_emu_differential_branching_fraction(
            q2,
            c9_lfv_semileptonic=c9_lfv_semileptonic,
            c10_lfv_semileptonic=c10_lfv_semileptonic,
            inputs=inputs,
            charge_mode=charge_mode,
        )
    return float(half_width * total)


def _angular_tensor_integral(
    q2: float,
    *,
    m_minus: float,
    m_plus: float,
    c9_lfv_semileptonic: complex,
    c10_lfv_semileptonic: complex,
    inputs: RareCharmDToPiMuMuInputs,
) -> float:
    nodes, weights = _legendre_nodes(_ANGULAR_QUADRATURE_POINTS)
    total = 0.0
    for cos_theta, weight in zip(nodes, weights):
        vector, axial = _tensor_contractions(
            q2,
            float(cos_theta),
            m_minus=m_minus,
            m_plus=m_plus,
            inputs=inputs,
        )
        total += float(weight) * (
            abs(c9_lfv_semileptonic) ** 2 * vector
            + abs(c10_lfv_semileptonic) ** 2 * axial
        )
    if total < 0.0 and abs(total) < 1.0e-24:
        return 0.0
    if total < 0.0:
        raise ValueError(f"negative angular tensor contraction {total}")
    return float(total)


def _tensor_contractions(
    q2: float,
    cos_theta: float,
    *,
    m_minus: float,
    m_plus: float,
    inputs: RareCharmDToPiMuMuInputs,
) -> tuple[float, float]:
    root_q2 = math.sqrt(q2)
    m_d = inputs.dplus_mass_gev
    m_pi = inputs.piplus_mass_gev
    m_d2 = m_d * m_d
    m_pi2 = m_pi * m_pi
    delta = m_d2 - m_pi2
    lam_had = max(0.0, _kallen(m_d2, m_pi2, q2))
    lam_lep = max(0.0, _kallen(q2, m_minus * m_minus, m_plus * m_plus))
    p_had = math.sqrt(lam_had) / (2.0 * root_q2)
    p_lep = math.sqrt(lam_lep) / (2.0 * root_q2)
    e_d = (m_d2 - m_pi2 + q2) / (2.0 * root_q2)
    e_pi = (m_d2 - m_pi2 - q2) / (2.0 * root_q2)
    e_minus = (q2 + m_minus * m_minus - m_plus * m_plus) / (2.0 * root_q2)
    e_plus = (q2 + m_plus * m_plus - m_minus * m_minus) / (2.0 * root_q2)
    sin_theta = math.sqrt(max(0.0, 1.0 - cos_theta * cos_theta))

    p_d = np.array([e_d, 0.0, 0.0, p_had], dtype=np.complex128)
    p_pi = np.array([e_pi, 0.0, 0.0, p_had], dtype=np.complex128)
    q_vec = np.array([root_q2, 0.0, 0.0, 0.0], dtype=np.complex128)
    p_minus = np.array(
        [e_minus, p_lep * sin_theta, 0.0, p_lep * cos_theta],
        dtype=np.complex128,
    )
    p_plus = np.array(
        [e_plus, -p_lep * sin_theta, 0.0, -p_lep * cos_theta],
        dtype=np.complex128,
    )
    fplus = dtopi_fplus(q2, inputs)
    fzero = dtopi_fzero(q2, inputs)
    hadronic = fplus * (p_d + p_pi - delta / q2 * q_vec) + fzero * delta / q2 * q_vec
    hadronic_cov = _METRIC @ hadronic
    p_dot = float(np.real(p_minus @ _METRIC @ p_plus))
    vector_tensor = _lepton_tensor(p_minus, p_plus, p_dot + m_minus * m_plus)
    axial_tensor = _lepton_tensor(p_minus, p_plus, p_dot - m_minus * m_plus)
    vector = _contract_hadronic_tensor(hadronic_cov, vector_tensor)
    axial = _contract_hadronic_tensor(hadronic_cov, axial_tensor)
    return float(vector), float(axial)


def _lepton_tensor(
    p_minus: np.ndarray,
    p_plus: np.ndarray,
    metric_term: float,
) -> np.ndarray:
    tensor = np.zeros((4, 4), dtype=np.complex128)
    for mu in range(4):
        for nu in range(4):
            tensor[mu, nu] = 4.0 * (
                p_minus[mu] * p_plus[nu]
                + p_minus[nu] * p_plus[mu]
                - _METRIC[mu, nu] * metric_term
            )
    return tensor


def _contract_hadronic_tensor(hadronic_cov: np.ndarray, tensor: np.ndarray) -> float:
    value = 0.0j
    for mu in range(4):
        for nu in range(4):
            value += hadronic_cov[mu] * np.conjugate(hadronic_cov[nu]) * tensor[mu, nu]
    out = float(np.real(value))
    if out < 0.0 and abs(out) < 1.0e-24:
        return 0.0
    if out < 0.0:
        raise ValueError(f"negative lepton-tensor contraction {out}")
    return out


@lru_cache(maxsize=None)
def _legendre_nodes(points: int) -> tuple[tuple[float, ...], tuple[float, ...]]:
    nodes, weights = np.polynomial.legendre.leggauss(int(points))
    return tuple(float(node) for node in nodes), tuple(float(weight) for weight in weights)


def _kallen(a: float, b: float, c: float) -> float:
    return float(a * a + b * b + c * c - 2.0 * (a * b + a * c + b * c))


__all__ = [
    "RARE_CHARM_DTOPI_EMU_MODEL_V1",
    "RARE_CHARM_DTOPI_EMU_OPERATOR_CONVENTION",
    "RARE_CHARM_DTOPI_EMU_PARAMETRIZATION_CITATION",
    "RARE_CHARM_DTOPI_EMU_PROXY_V1",
    "RARE_CHARM_DTOPI_EMU_Q2_TREATMENT_V1",
    "RareCharmDToPiLFVBranchingResult",
    "RareCharmDToPiMuMuBranchingResult",
    "RareCharmDToPiMuMuInputs",
    "default_dtopi_emu_inputs",
    "dtopi_emu_q2_range",
    "dtopi_emu_differential_branching_fraction",
    "dtopi_emu_sm",
    "evaluate_dplus_to_piplus_emu",
]
