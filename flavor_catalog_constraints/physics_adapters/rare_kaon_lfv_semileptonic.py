"""LFV ``K+ -> pi+ e mu`` semileptonic rare-kaon adapter.

This adapter composes the shared ``s -> d l l`` rare-kaon Wilson proxy with
the K019 e-mu LFV lepton-coupling pattern, then evaluates a full three-body
``K -> pi`` vector-current rate.  The quark-side matching is deliberately the
same documented Z-like proxy used by the rare-kaon dilepton path, but the
pseudoscalar-to-pseudoscalar hadronic current uses the vector combination
``left + right`` rather than the K_L two-body ``left - right`` combination.

NEEDS-HUMAN-PHYSICS: a rigorous RS prediction requires off-diagonal physical
charged-lepton neutral-current couplings after EW KK/Z/Z' mixing and
charged-lepton mass-basis rotations.  ``ParameterPoint`` does not carry those
couplings, so this v1 path accepts an explicit e-mu lepton proxy and reports
that limitation in diagnostics.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache
import math
from typing import Any, Mapping

import numpy as np

from .rare_kaon_dilepton import (
    QuarkMassBasisCouplings,
    RARE_KAON_DILEPTON_INPUT_BUNDLE_V1,
    RARE_KAON_DILEPTON_MODEL_V1,
    RARE_KAON_DILEPTON_RS_MATCHING_ASSUMPTION_V1,
    RareKaonDileptonSMInputs,
    RareKaonDileptonWilsonCoefficients,
    rare_kaon_dilepton_ckm_factors,
    rare_kaon_dilepton_default_sm_inputs,
    rare_kaon_dilepton_g_sm_squared,
    rare_kaon_dilepton_wilsons_from_couplings,
)
from .rare_kaon_lfv_dilepton import (
    RARE_KAON_LFV_DILEPTON_PARAMETRIZATION_CITATION,
    RARE_KAON_LFV_DILEPTON_PROXY_V1,
    RareKaonLFVLeptonCouplingProxy,
    RareKaonLFVLeptonProxyInput,
    rare_kaon_lfv_lepton_coupling_proxy,
    rare_kaon_lfv_proxy_input,
)

RARE_KAON_KTOPI_EMU_MODEL_V1 = "rare_kaon_kplus_piplus_emu_lfv_form_factor_proxy_v1"
RARE_KAON_KTOPI_EMU_INPUT_BUNDLE_V1 = (
    "rare_kaon_kplus_piplus_emu_full_q2_form_factor_inputs_v1"
)
RARE_KAON_KTOPI_EMU_FORM_FACTOR_MODEL_V1 = (
    "K_to_pi_fplus_fzero_single_pole_fplus0_FLAG2024_v1"
)
RARE_KAON_KTOPI_EMU_OPERATOR_CONVENTION = (
    "H_eff=(G_F/sqrt(2))*alpha/(2*pi*sin2thetaW) "
    "[YV_LFV (sbar gamma_mu d)(ebar gamma^mu mu) + "
    "YA_LFV (sbar gamma_mu d)(ebar gamma^mu gamma5 mu)]; "
    "K->pi uses the quark vector current, so left/right s-d proxy pieces "
    "enter as left+right."
)
RARE_KAON_KTOPI_EMU_PARAMETRIZATION_CITATION = (
    RARE_KAON_LFV_DILEPTON_PARAMETRIZATION_CITATION
    + "; K->pi vector-current three-body phase space with f_+(q2)"
)
RARE_KAON_KTOPI_EMU_PROXY_V1 = (
    "NEEDS-HUMAN-PHYSICS: K+ -> pi+ e mu v1 reuses the rare-kaon s->d "
    "neutral-current proxy and the K019 e-mu lepton spurion. Scalar/tensor "
    "operators, charge-orientation-specific lepton matching, and a complete "
    "RS EW KK/Z/Z' tower are not available on ParameterPoint."
)
RARE_KAON_KTOPI_EMU_Q2_TREATMENT_V1 = (
    "full_kinematic_q2_short_distance_lfv_proxy_with_k_to_pi_form_factor"
)

_METRIC = np.diag([1.0, -1.0, -1.0, -1.0])
_ANGULAR_QUADRATURE_POINTS = 48
_DEFAULT_CHARGE_STATE_FACTOR = 1.0


@dataclass(frozen=True)
class RareKaonKToPiFormFactorInputs:
    """K-to-pi vector-current form-factor proxy inputs."""

    model_label: str = RARE_KAON_KTOPI_EMU_FORM_FACTOR_MODEL_V1
    fplus_0: float = 0.9698
    fzero_0: float = 0.9698
    vector_pole_mass_gev: float = 0.89166
    scalar_pole_mass_gev: float = 1.43
    source: str = (
        "FLAG/PDG-era K_l3 f_+(0)=0.9698 and single-pole K*(892) vector "
        "shape; scalar-pole f_0 included only for the unequal-lepton "
        "q_mu component of the same vector current."
    )

    def __post_init__(self) -> None:
        for name in (
            "fplus_0",
            "fzero_0",
            "vector_pole_mass_gev",
            "scalar_pole_mass_gev",
        ):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be positive and finite")


@dataclass(frozen=True)
class RareKaonKToPiLFVInputs:
    """Numerical inputs for the ``K+ -> pi+ e mu`` LFV proxy."""

    input_bundle: str = RARE_KAON_KTOPI_EMU_INPUT_BUNDLE_V1
    rare_kaon: RareKaonDileptonSMInputs = field(
        default_factory=rare_kaon_dilepton_default_sm_inputs
    )
    form_factor: RareKaonKToPiFormFactorInputs = field(
        default_factory=RareKaonKToPiFormFactorInputs
    )
    kplus_mass_gev: float = 0.493677
    piplus_mass_gev: float = 0.13957039
    kplus_lifetime_s: float = 1.2380e-8
    hbar_gev_s: float = 6.582119569e-25
    electron_mass_gev: float = 0.00051099895
    muon_mass_gev: float = 0.1056583745
    charge_state_factor: float = _DEFAULT_CHARGE_STATE_FACTOR
    quadrature_points: int = 160
    constants_citation: str = (
        "PDG-era K+/pi+/charged-lepton masses and K+ lifetime; G_F, alpha, "
        "sin2thetaW, and CKM from the shared rare_kaon_dilepton input bundle"
    )

    def __post_init__(self) -> None:
        for name in (
            "kplus_mass_gev",
            "piplus_mass_gev",
            "kplus_lifetime_s",
            "hbar_gev_s",
            "electron_mass_gev",
            "muon_mass_gev",
            "charge_state_factor",
        ):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be positive and finite")
        if int(self.quadrature_points) < 16:
            raise ValueError("quadrature_points must be at least 16")
        if self.kplus_mass_gev <= (
            self.piplus_mass_gev + self.electron_mass_gev + self.muon_mass_gev
        ):
            raise ValueError("K+ -> pi+ e mu phase space is closed")


@dataclass(frozen=True)
class RareKaonKToPiLFVWilsonCoefficients:
    """LFV ``s -> d e mu`` Wilson proxy for the ``K -> pi`` vector current."""

    model_label: str
    base_model_label: str
    operator_convention: str
    matching_assumption: str
    M_KK: float
    matching_scale: float
    left_sd_coupling: complex
    right_sd_coupling: complex
    left_sd_overlap: complex
    right_sd_overlap: complex
    left_quark_delta: complex
    right_quark_delta: complex
    quark_vector_delta: complex
    lepton_left_delta_emu: complex
    lepton_right_delta_emu: complex
    lepton_vector_delta_emu: complex
    lepton_axial_delta_emu: complex
    y_vector_lfv: complex
    y_axial_lfv: complex
    base_same_flavor_wilsons: RareKaonDileptonWilsonCoefficients
    lepton_proxy: RareKaonLFVLeptonCouplingProxy

    @property
    def wilsons(self) -> Mapping[str, complex]:
        return {
            "YV_LFV_semileptonic": complex(self.y_vector_lfv),
            "YA_LFV_semileptonic": complex(self.y_axial_lfv),
        }


@dataclass(frozen=True)
class RareKaonKToPiLFVBranchingResult:
    """Short-distance branching-ratio prediction for one LFV charge mode."""

    model_label: str
    input_bundle: str
    charge_mode: str
    branching_fraction: float
    sm_branching_fraction: float
    np_shift_branching_fraction: float
    q2_min_gev2: float
    q2_max_gev2: float
    electron_mass_gev: float
    muon_mass_gev: float
    y_vector_lfv: complex
    y_axial_lfv: complex
    lambda_wolfenstein: float
    wilsons: RareKaonKToPiLFVWilsonCoefficients | None = None
    diagnostics: Mapping[str, Any] = field(default_factory=dict)


def rare_kaon_ktopi_emu_default_inputs() -> RareKaonKToPiLFVInputs:
    """Return the default ``K+ -> pi+ e mu`` LFV input bundle."""

    return RareKaonKToPiLFVInputs()


def rare_kaon_ktopi_emu_proxy_input(
    left_emu_overlap: complex,
    right_emu_overlap: complex,
    m_kk_gev: float,
    *,
    source: str = "caller-supplied rare-kaon K->pi e-mu lepton proxy",
) -> RareKaonLFVLeptonProxyInput:
    """Build a proxy input accepted by the K019/K020 LFV adapters."""

    return rare_kaon_lfv_proxy_input(
        left_emu_overlap,
        right_emu_overlap,
        m_kk_gev,
        source=source,
    )


def rare_kaon_ktopi_emu_lepton_coupling_proxy(
    source: object,
    *,
    m_kk_gev: float | None = None,
    inputs: RareKaonKToPiLFVInputs | None = None,
) -> RareKaonLFVLeptonCouplingProxy:
    """Return the documented Z-like e-mu lepton-coupling proxy."""

    p = rare_kaon_ktopi_emu_default_inputs() if inputs is None else inputs
    return rare_kaon_lfv_lepton_coupling_proxy(
        source,
        m_kk_gev=m_kk_gev,
        inputs=p,
    )


def rare_kaon_ktopi_fplus(
    q2_gev2: float,
    inputs: RareKaonKToPiLFVInputs | None = None,
) -> float:
    """Return the proxy ``K -> pi`` vector form factor ``f_+(q2)``."""

    p = rare_kaon_ktopi_emu_default_inputs() if inputs is None else inputs
    q2 = float(q2_gev2)
    ff = p.form_factor
    denominator = 1.0 - q2 / ff.vector_pole_mass_gev**2
    if denominator <= 0.0 or not math.isfinite(denominator):
        raise ValueError(f"invalid K->pi f_+ denominator at q2={q2_gev2}")
    return float(ff.fplus_0 / denominator)


def rare_kaon_ktopi_fzero(
    q2_gev2: float,
    inputs: RareKaonKToPiLFVInputs | None = None,
) -> float:
    """Return the proxy ``K -> pi`` scalar form factor ``f_0(q2)``."""

    p = rare_kaon_ktopi_emu_default_inputs() if inputs is None else inputs
    q2 = float(q2_gev2)
    ff = p.form_factor
    denominator = 1.0 - q2 / ff.scalar_pole_mass_gev**2
    if denominator <= 0.0 or not math.isfinite(denominator):
        raise ValueError(f"invalid K->pi f_0 denominator at q2={q2_gev2}")
    return float(ff.fzero_0 / denominator)


def rare_kaon_ktopi_emu_q2_range(
    inputs: RareKaonKToPiLFVInputs | None = None,
) -> tuple[float, float]:
    """Return the kinematic ``q2`` range for ``K+ -> pi+ e mu``."""

    p = rare_kaon_ktopi_emu_default_inputs() if inputs is None else inputs
    q2_min = (p.electron_mass_gev + p.muon_mass_gev) ** 2
    q2_max = (p.kplus_mass_gev - p.piplus_mass_gev) ** 2
    if q2_max <= q2_min:
        raise ValueError("K+ -> pi+ e mu phase space is closed")
    return float(q2_min), float(q2_max)


def rare_kaon_ktopi_emu_wilsons_from_couplings(
    quark_couplings: QuarkMassBasisCouplings,
    lepton_couplings: object,
    *,
    m_kk_gev: float | None = None,
    inputs: RareKaonKToPiLFVInputs | None = None,
) -> RareKaonKToPiLFVWilsonCoefficients:
    """Return the v1 LFV Wilson proxy for the ``K -> pi`` vector current."""

    p = rare_kaon_ktopi_emu_default_inputs() if inputs is None else inputs
    base = rare_kaon_dilepton_wilsons_from_couplings(
        quark_couplings,
        m_kk_gev=m_kk_gev,
        inputs=p.rare_kaon,
    )
    proxy = rare_kaon_lfv_lepton_coupling_proxy(
        lepton_couplings,
        m_kk_gev=base.M_KK,
        inputs=p,
    )
    quark_vector_delta = complex(base.left_quark_delta + base.right_quark_delta)
    normalizer = rare_kaon_dilepton_g_sm_squared(p.rare_kaon) * base.M_KK**2
    y_vector = quark_vector_delta * proxy.lepton_vector_delta_emu / normalizer
    y_axial = quark_vector_delta * proxy.lepton_axial_delta_emu / normalizer
    return RareKaonKToPiLFVWilsonCoefficients(
        model_label=RARE_KAON_KTOPI_EMU_MODEL_V1,
        base_model_label=RARE_KAON_DILEPTON_MODEL_V1,
        operator_convention=RARE_KAON_KTOPI_EMU_OPERATOR_CONVENTION,
        matching_assumption=RARE_KAON_KTOPI_EMU_PROXY_V1,
        M_KK=float(base.M_KK),
        matching_scale=float(base.matching_scale),
        left_sd_coupling=complex(base.left_sd_coupling),
        right_sd_coupling=complex(base.right_sd_coupling),
        left_sd_overlap=complex(base.left_sd_overlap),
        right_sd_overlap=complex(base.right_sd_overlap),
        left_quark_delta=complex(base.left_quark_delta),
        right_quark_delta=complex(base.right_quark_delta),
        quark_vector_delta=quark_vector_delta,
        lepton_left_delta_emu=complex(proxy.lepton_left_delta_emu),
        lepton_right_delta_emu=complex(proxy.lepton_right_delta_emu),
        lepton_vector_delta_emu=complex(proxy.lepton_vector_delta_emu),
        lepton_axial_delta_emu=complex(proxy.lepton_axial_delta_emu),
        y_vector_lfv=complex(y_vector),
        y_axial_lfv=complex(y_axial),
        base_same_flavor_wilsons=base,
        lepton_proxy=proxy,
    )


def rare_kaon_ktopi_emu_differential_branching_fraction(
    q2_gev2: float,
    *,
    y_vector_lfv: complex,
    y_axial_lfv: complex,
    inputs: RareKaonKToPiLFVInputs | None = None,
    charge_mode: str = "muplus_eminus",
) -> float:
    """Evaluate ``dBR_SD(K+ -> pi+ e mu) / dq2`` for one LFV charge mode."""

    p = rare_kaon_ktopi_emu_default_inputs() if inputs is None else inputs
    q2 = float(q2_gev2)
    q2_min, q2_max = rare_kaon_ktopi_emu_q2_range(p)
    if q2 < q2_min or q2 > q2_max:
        raise ValueError(f"q2={q2} outside K+ -> pi+ e mu physical range")

    if charge_mode == "muplus_eminus":
        m_minus = p.electron_mass_gev
        m_plus = p.muon_mass_gev
    elif charge_mode == "muminus_eplus":
        m_minus = p.muon_mass_gev
        m_plus = p.electron_mass_gev
    else:
        raise ValueError(f"unsupported K+ -> pi+ e mu charge mode {charge_mode!r}")

    m_k = p.kplus_mass_gev
    m_pi = p.piplus_mass_gev
    lam_had = max(0.0, _kallen(m_k * m_k, m_pi * m_pi, q2))
    lam_lep = max(0.0, _kallen(q2, m_minus * m_minus, m_plus * m_plus))
    if lam_had <= 0.0 or lam_lep <= 0.0:
        return 0.0

    angular = _angular_tensor_integral(
        q2,
        m_minus=m_minus,
        m_plus=m_plus,
        y_vector_lfv=complex(y_vector_lfv),
        y_axial_lfv=complex(y_axial_lfv),
        inputs=p,
    )
    c0 = _effective_hamiltonian_prefactor(p.rare_kaon)
    tau = p.kplus_lifetime_s / p.hbar_gev_s
    phase_space = math.sqrt(lam_had) * math.sqrt(lam_lep) / (
        512.0 * math.pi**3 * m_k**3 * q2
    )
    rate = p.charge_state_factor * tau * c0 * c0 * phase_space * angular
    if rate < 0.0 and abs(rate) < 1.0e-30:
        return 0.0
    if rate < 0.0:
        raise ValueError(f"negative dBR/dq2={rate} at q2={q2}")
    return float(rate)


def kplus_piplus_emu_sm(
    inputs: RareKaonKToPiLFVInputs | None = None,
    *,
    charge_mode: str = "muplus_eminus",
) -> RareKaonKToPiLFVBranchingResult:
    """Return the catalog SM-limit ``K+ -> pi+ e mu`` rate, zero for LFV."""

    p = rare_kaon_ktopi_emu_default_inputs() if inputs is None else inputs
    return _branching_from_wilsons(None, inputs=p, charge_mode=charge_mode)


def kplus_piplus_emu_from_couplings(
    quark_couplings: QuarkMassBasisCouplings | RareKaonKToPiLFVWilsonCoefficients,
    lepton_couplings: object | None = None,
    *,
    m_kk_gev: float | None = None,
    inputs: RareKaonKToPiLFVInputs | None = None,
    charge_mode: str = "muplus_eminus",
) -> RareKaonKToPiLFVBranchingResult:
    """Evaluate smooth ``BR_SD(K+ -> pi+ e mu)`` from proxy couplings."""

    p = rare_kaon_ktopi_emu_default_inputs() if inputs is None else inputs
    if isinstance(quark_couplings, RareKaonKToPiLFVWilsonCoefficients):
        wilsons = quark_couplings
    else:
        if lepton_couplings is None:
            raise TypeError("lepton_couplings is required for K+ -> pi+ e mu LFV matching")
        wilsons = rare_kaon_ktopi_emu_wilsons_from_couplings(
            quark_couplings,
            lepton_couplings,
            m_kk_gev=m_kk_gev,
            inputs=p,
        )
    return _branching_from_wilsons(wilsons, inputs=p, charge_mode=charge_mode)


def _branching_from_wilsons(
    wilsons: RareKaonKToPiLFVWilsonCoefficients | None,
    *,
    inputs: RareKaonKToPiLFVInputs,
    charge_mode: str,
) -> RareKaonKToPiLFVBranchingResult:
    if wilsons is None:
        y_vector = 0.0j
        y_axial = 0.0j
    else:
        y_vector = complex(wilsons.y_vector_lfv)
        y_axial = complex(wilsons.y_axial_lfv)

    q2_min, q2_max = rare_kaon_ktopi_emu_q2_range(inputs)
    branching = _integrate_branching_fraction(
        y_vector_lfv=y_vector,
        y_axial_lfv=y_axial,
        inputs=inputs,
        charge_mode=charge_mode,
    )
    factors = rare_kaon_dilepton_ckm_factors(inputs.rare_kaon)
    q2_mid = 0.5 * (q2_min + q2_max)
    diagnostics: dict[str, Any] = {
        "base_model_label": RARE_KAON_DILEPTON_MODEL_V1,
        "base_input_bundle": RARE_KAON_DILEPTON_INPUT_BUNDLE_V1,
        "base_matching_assumption": RARE_KAON_DILEPTON_RS_MATCHING_ASSUMPTION_V1,
        "lfv_lepton_matching_assumption": RARE_KAON_LFV_DILEPTON_PROXY_V1,
        "matching_assumption": RARE_KAON_KTOPI_EMU_PROXY_V1,
        "operator_convention": RARE_KAON_KTOPI_EMU_OPERATOR_CONVENTION,
        "parametrization_citation": RARE_KAON_KTOPI_EMU_PARAMETRIZATION_CITATION,
        "q2_treatment": RARE_KAON_KTOPI_EMU_Q2_TREATMENT_V1,
        "short_distance_only": True,
        "sm_branching_fraction": 0.0,
        "sm_lfv_policy": (
            "K+ -> pi+ e mu is charged-LFV and has zero SM rate for catalog "
            "purposes."
        ),
        "charge_mode": charge_mode,
        "charge_state_factor": float(inputs.charge_state_factor),
        "charge_conjugate_modes_included": False,
        "charge_mode_specific_lepton_matching_available": False,
        "quark_vector_current_uses_left_plus_right": True,
        "scalar_tensor_operators_included": False,
        "q2_min_gev2": float(q2_min),
        "q2_max_gev2": float(q2_max),
        "quadrature_points": float(inputs.quadrature_points),
        "angular_quadrature_points": float(_ANGULAR_QUADRATURE_POINTS),
        "effective_hamiltonian_prefactor_gev_minus2": float(
            _effective_hamiltonian_prefactor(inputs.rare_kaon)
        ),
        "g_sm_squared_gev_minus2": float(
            rare_kaon_dilepton_g_sm_squared(inputs.rare_kaon)
        ),
        "lambda_wolfenstein": float(factors.lambda_wolfenstein),
        "lambda_c": complex(factors.lambda_c),
        "lambda_t": complex(factors.lambda_t),
        "kplus_mass_gev": float(inputs.kplus_mass_gev),
        "piplus_mass_gev": float(inputs.piplus_mass_gev),
        "kplus_lifetime_s": float(inputs.kplus_lifetime_s),
        "electron_mass_gev": float(inputs.electron_mass_gev),
        "muon_mass_gev": float(inputs.muon_mass_gev),
        "form_factor_model": inputs.form_factor.model_label,
        "form_factor_source": inputs.form_factor.source,
        "fplus_0": float(inputs.form_factor.fplus_0),
        "fzero_0": float(inputs.form_factor.fzero_0),
        "form_factor_vector_pole_mass_gev": float(
            inputs.form_factor.vector_pole_mass_gev
        ),
        "form_factor_scalar_pole_mass_gev": float(
            inputs.form_factor.scalar_pole_mass_gev
        ),
        "fplus_q2_min": rare_kaon_ktopi_fplus(q2_min, inputs),
        "fplus_q2_mid": rare_kaon_ktopi_fplus(q2_mid, inputs),
        "fplus_q2_max": rare_kaon_ktopi_fplus(q2_max, inputs),
        "constants_citation": inputs.constants_citation,
        "y_vector_lfv": complex(y_vector),
        "y_axial_lfv": complex(y_axial),
    }
    if wilsons is not None:
        diagnostics.update(
            {
                "m_kk_gev": float(wilsons.M_KK),
                "matching_scale_gev": float(wilsons.matching_scale),
                "left_sd_coupling": complex(wilsons.left_sd_coupling),
                "right_sd_coupling": complex(wilsons.right_sd_coupling),
                "left_sd_overlap": complex(wilsons.left_sd_overlap),
                "right_sd_overlap": complex(wilsons.right_sd_overlap),
                "left_quark_delta": complex(wilsons.left_quark_delta),
                "right_quark_delta": complex(wilsons.right_quark_delta),
                "quark_vector_delta": complex(wilsons.quark_vector_delta),
                "lepton_left_delta_emu": complex(wilsons.lepton_left_delta_emu),
                "lepton_right_delta_emu": complex(wilsons.lepton_right_delta_emu),
                "lepton_vector_delta_emu": complex(wilsons.lepton_vector_delta_emu),
                "lepton_axial_delta_emu": complex(wilsons.lepton_axial_delta_emu),
                "base_muon_axial_delta": float(
                    wilsons.base_same_flavor_wilsons.muon_axial_delta
                ),
                "base_y_np_left": complex(
                    wilsons.base_same_flavor_wilsons.y_np_left
                ),
                "base_y_np_right": complex(
                    wilsons.base_same_flavor_wilsons.y_np_right
                ),
                "proxy_source": wilsons.lepton_proxy.source,
            }
        )

    return RareKaonKToPiLFVBranchingResult(
        model_label=RARE_KAON_KTOPI_EMU_MODEL_V1,
        input_bundle=inputs.input_bundle,
        charge_mode=charge_mode,
        branching_fraction=float(branching),
        sm_branching_fraction=0.0,
        np_shift_branching_fraction=float(branching),
        q2_min_gev2=float(q2_min),
        q2_max_gev2=float(q2_max),
        electron_mass_gev=float(inputs.electron_mass_gev),
        muon_mass_gev=float(inputs.muon_mass_gev),
        y_vector_lfv=complex(y_vector),
        y_axial_lfv=complex(y_axial),
        lambda_wolfenstein=float(factors.lambda_wolfenstein),
        wilsons=wilsons,
        diagnostics=diagnostics,
    )


def _integrate_branching_fraction(
    *,
    y_vector_lfv: complex,
    y_axial_lfv: complex,
    inputs: RareKaonKToPiLFVInputs,
    charge_mode: str,
) -> float:
    q2_min, q2_max = rare_kaon_ktopi_emu_q2_range(inputs)
    nodes, weights = np.polynomial.legendre.leggauss(int(inputs.quadrature_points))
    mid = 0.5 * (q2_max + q2_min)
    half_width = 0.5 * (q2_max - q2_min)
    total = 0.0
    for node, weight in zip(nodes, weights):
        q2 = mid + half_width * float(node)
        total += float(weight) * rare_kaon_ktopi_emu_differential_branching_fraction(
            q2,
            y_vector_lfv=y_vector_lfv,
            y_axial_lfv=y_axial_lfv,
            inputs=inputs,
            charge_mode=charge_mode,
        )
    return float(half_width * total)


def _angular_tensor_integral(
    q2: float,
    *,
    m_minus: float,
    m_plus: float,
    y_vector_lfv: complex,
    y_axial_lfv: complex,
    inputs: RareKaonKToPiLFVInputs,
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
            abs(y_vector_lfv) ** 2 * vector + abs(y_axial_lfv) ** 2 * axial
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
    inputs: RareKaonKToPiLFVInputs,
) -> tuple[float, float]:
    root_q2 = math.sqrt(q2)
    m_k = inputs.kplus_mass_gev
    m_pi = inputs.piplus_mass_gev
    m_k2 = m_k * m_k
    m_pi2 = m_pi * m_pi
    delta = m_k2 - m_pi2
    lam_had = max(0.0, _kallen(m_k2, m_pi2, q2))
    lam_lep = max(0.0, _kallen(q2, m_minus * m_minus, m_plus * m_plus))
    p_had = math.sqrt(lam_had) / (2.0 * root_q2)
    p_lep = math.sqrt(lam_lep) / (2.0 * root_q2)
    e_k = (m_k2 - m_pi2 + q2) / (2.0 * root_q2)
    e_pi = (m_k2 - m_pi2 - q2) / (2.0 * root_q2)
    e_minus = (q2 + m_minus * m_minus - m_plus * m_plus) / (2.0 * root_q2)
    e_plus = (q2 + m_plus * m_plus - m_minus * m_minus) / (2.0 * root_q2)
    sin_theta = math.sqrt(max(0.0, 1.0 - cos_theta * cos_theta))

    p_k = np.array([e_k, 0.0, 0.0, p_had], dtype=np.complex128)
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
    fplus = rare_kaon_ktopi_fplus(q2, inputs)
    fzero = rare_kaon_ktopi_fzero(q2, inputs)
    hadronic = fplus * (p_k + p_pi - delta / q2 * q_vec) + fzero * delta / q2 * q_vec
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


def _effective_hamiltonian_prefactor(inputs: RareKaonDileptonSMInputs) -> float:
    return float(
        inputs.gf_gev_minus2
        / math.sqrt(2.0)
        * inputs.alpha_em_mz
        / (2.0 * math.pi * inputs.sin2_theta_w)
    )


def _kallen(a: float, b: float, c: float) -> float:
    return float(a * a + b * b + c * c - 2.0 * (a * b + a * c + b * c))


__all__ = [
    "QuarkMassBasisCouplings",
    "RareKaonKToPiFormFactorInputs",
    "RareKaonKToPiLFVInputs",
    "RareKaonKToPiLFVWilsonCoefficients",
    "RareKaonKToPiLFVBranchingResult",
    "RareKaonLFVLeptonProxyInput",
    "RareKaonLFVLeptonCouplingProxy",
    "RARE_KAON_KTOPI_EMU_MODEL_V1",
    "RARE_KAON_KTOPI_EMU_INPUT_BUNDLE_V1",
    "RARE_KAON_KTOPI_EMU_FORM_FACTOR_MODEL_V1",
    "RARE_KAON_KTOPI_EMU_OPERATOR_CONVENTION",
    "RARE_KAON_KTOPI_EMU_PARAMETRIZATION_CITATION",
    "RARE_KAON_KTOPI_EMU_PROXY_V1",
    "RARE_KAON_KTOPI_EMU_Q2_TREATMENT_V1",
    "RARE_KAON_LFV_DILEPTON_PROXY_V1",
    "rare_kaon_ktopi_emu_default_inputs",
    "rare_kaon_ktopi_emu_proxy_input",
    "rare_kaon_ktopi_emu_lepton_coupling_proxy",
    "rare_kaon_ktopi_fplus",
    "rare_kaon_ktopi_fzero",
    "rare_kaon_ktopi_emu_q2_range",
    "rare_kaon_ktopi_emu_wilsons_from_couplings",
    "rare_kaon_ktopi_emu_differential_branching_fraction",
    "kplus_piplus_emu_sm",
    "kplus_piplus_emu_from_couplings",
]
