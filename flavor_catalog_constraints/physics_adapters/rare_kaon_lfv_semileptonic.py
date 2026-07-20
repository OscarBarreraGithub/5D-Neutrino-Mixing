"""LFV ``K+ -> pi+ e mu`` semileptonic rare-kaon adapter.

The production Phase-4c path reads the rigorous off-diagonal ``e mu`` block
from ``rs_semileptonic_wilsons.lfv_llqq`` and maps its C9/C10 inputs into the
rare-kaon K-to-pi Y normalization without reusing the old one-boson proxy.  The
older explicit-spurion helpers remain as legacy compatibility entry points
only.  The pseudoscalar-to-pseudoscalar hadronic current uses the vector
combination ``left + right`` rather than the K_L two-body ``left - right``
combination.

With the current diagonal charged-lepton fit, the Phase-4a LFV block is
rigorously zero at tree level, so the evaluated LFV rate is zero and
non-vetoing.  Nonzero tree-level LFV requires non-diagonal lepton structure;
loop-induced LFV is deferred.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field, replace
from functools import lru_cache
from typing import Any, Mapping

import numpy as np

from quarkConstraints.rs_semileptonic_wilsons import RSSemileptonicWilsonBundle

from .rare_kaon_dilepton import (
    RARE_KAON_DILEPTON_INPUT_BUNDLE_V1,
    RARE_KAON_DILEPTON_MODEL_V1,
    RARE_KAON_DILEPTON_RS_MATCHING_ASSUMPTION_V1,
    QuarkMassBasisCouplings,
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
    RARE_KAON_LFV_TREE_LEVEL_NOTE_V1,
    RareKaonLFVLeptonCouplingProxy,
    RareKaonLFVLeptonProxyInput,
    rare_kaon_lfv_coeff_from_rs_semileptonic,
    rare_kaon_lfv_lepton_coupling_proxy,
    rare_kaon_lfv_proxy_input,
    rare_kaon_lfv_rs_semileptonic_diagnostics,
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
    "LEGACY-PROXY: K+ -> pi+ e mu explicit-spurion compatibility path. The "
    "Phase-4c K020 production path uses rs_semileptonic_wilsons.lfv_llqq; "
    "scalar/tensor operators and charge-orientation-specific loop effects "
    "remain deferred."
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
    base_same_flavor_wilsons: RareKaonDileptonWilsonCoefficients | None
    lepton_proxy: RareKaonLFVLeptonCouplingProxy | None

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


def rare_kaon_ktopi_emu_wilsons_from_rs_semileptonic(
    source: RSSemileptonicWilsonBundle,
    *,
    matching_scale_gev: float | None = None,
    inputs: RareKaonKToPiLFVInputs | None = None,
    lepton_pair: str = "e_mu",
) -> RareKaonKToPiLFVWilsonCoefficients:
    """Map Phase-4a LFV C9/C10 Wilsons into the K->pi LFV Y inputs."""

    p = rare_kaon_ktopi_emu_default_inputs() if inputs is None else inputs
    coeff = rare_kaon_lfv_coeff_from_rs_semileptonic(
        source,
        transition="s_to_d",
        lepton_pair=lepton_pair,
    )
    scale = _diagnostic_matching_scale(matching_scale_gev)
    norm = _kaon_lfv_c9_c10_to_y_norm(p)
    lambda_ckm = complex(coeff.lambda_ckm)
    left_vector_contact = complex(coeff.contact_LL + coeff.contact_LR)
    right_vector_contact = complex(coeff.contact_RL + coeff.contact_RR)
    quark_vector_contact = left_vector_contact + right_vector_contact
    y_vector = -lambda_ckm * norm * complex(
        coeff.c9_lfv_np + coeff.c9p_lfv_np
    )
    y_axial = -lambda_ckm * norm * complex(
        coeff.c10_lfv_np + coeff.c10p_lfv_np
    )
    return RareKaonKToPiLFVWilsonCoefficients(
        model_label=coeff.model_label,
        base_model_label=RARE_KAON_DILEPTON_MODEL_V1,
        operator_convention=coeff.operator_convention,
        matching_assumption=coeff.matching_assumption,
        M_KK=scale,
        matching_scale=scale,
        left_sd_coupling=left_vector_contact,
        right_sd_coupling=right_vector_contact,
        left_sd_overlap=0.0j,
        right_sd_overlap=0.0j,
        left_quark_delta=left_vector_contact,
        right_quark_delta=right_vector_contact,
        quark_vector_delta=quark_vector_contact,
        lepton_left_delta_emu=complex(coeff.contact_LL + coeff.contact_RL),
        lepton_right_delta_emu=complex(coeff.contact_LR + coeff.contact_RR),
        lepton_vector_delta_emu=complex(
            coeff.contact_LL + coeff.contact_LR + coeff.contact_RL + coeff.contact_RR
        ),
        lepton_axial_delta_emu=complex(
            coeff.contact_LR + coeff.contact_RR - coeff.contact_LL - coeff.contact_RL
        ),
        y_vector_lfv=complex(y_vector),
        y_axial_lfv=complex(y_axial),
        base_same_flavor_wilsons=None,
        lepton_proxy=None,
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


def kplus_piplus_emu_from_rs_semileptonic_wilsons(
    source: RSSemileptonicWilsonBundle,
    *,
    matching_scale_gev: float | None = None,
    inputs: RareKaonKToPiLFVInputs | None = None,
    charge_mode: str = "muplus_eminus",
    lepton_pair: str = "e_mu",
) -> RareKaonKToPiLFVBranchingResult:
    """Evaluate ``K+ -> pi+ e mu`` from Phase-4a LFV llqq Wilsons."""

    p = rare_kaon_ktopi_emu_default_inputs() if inputs is None else inputs
    coeff = rare_kaon_lfv_coeff_from_rs_semileptonic(
        source,
        transition="s_to_d",
        lepton_pair=lepton_pair,
    )
    wilsons = rare_kaon_ktopi_emu_wilsons_from_rs_semileptonic(
        source,
        matching_scale_gev=matching_scale_gev,
        inputs=p,
        lepton_pair=lepton_pair,
    )
    result = _branching_from_wilsons(wilsons, inputs=p, charge_mode=charge_mode)
    diagnostics = dict(result.diagnostics)
    diagnostics.update(
        rare_kaon_lfv_rs_semileptonic_diagnostics(coeff, source=source)
    )
    diagnostics["base_matching_assumption"] = coeff.matching_assumption
    diagnostics["lfv_lepton_matching_assumption"] = coeff.matching_assumption
    diagnostics["matching_assumption"] = coeff.matching_assumption
    diagnostics["s_to_d_lfv_llqq_rs_semileptonic_rewired"] = True
    diagnostics["rare_kaon_lfv_proxy_reused"] = False
    diagnostics["c9_c10_to_rare_kaon_effective_inputs"] = True
    diagnostics["lfv_tree_level_note"] = RARE_KAON_LFV_TREE_LEVEL_NOTE_V1
    return replace(result, diagnostics=diagnostics)


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
        "lfv_lepton_matching_assumption": (
            RARE_KAON_LFV_DILEPTON_PROXY_V1
            if wilsons is None
            else wilsons.matching_assumption
        ),
        "matching_assumption": (
            RARE_KAON_KTOPI_EMU_PROXY_V1
            if wilsons is None
            else wilsons.matching_assumption
        ),
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
            }
        )
        if wilsons.base_same_flavor_wilsons is not None:
            diagnostics.update(
                {
                    "base_muon_axial_delta": float(
                        wilsons.base_same_flavor_wilsons.muon_axial_delta
                    ),
                    "base_y_np_left": complex(
                        wilsons.base_same_flavor_wilsons.y_np_left
                    ),
                    "base_y_np_right": complex(
                        wilsons.base_same_flavor_wilsons.y_np_right
                    ),
                }
            )
        if wilsons.lepton_proxy is not None:
            diagnostics["proxy_source"] = wilsons.lepton_proxy.source

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


def _diagnostic_matching_scale(matching_scale_gev: float | None) -> float:
    if matching_scale_gev is None:
        return 0.0
    number = float(matching_scale_gev)
    if not math.isfinite(number) or number <= 0.0:
        raise ValueError("matching_scale_gev must be positive and finite")
    return number


def _kaon_lfv_c9_c10_to_y_norm(inputs: RareKaonKToPiLFVInputs) -> float:
    return float(
        math.sqrt(2.0)
        * inputs.rare_kaon.gf_gev_minus2
        * inputs.rare_kaon.alpha_em_mz
        / (math.pi * rare_kaon_dilepton_g_sm_squared(inputs.rare_kaon))
    )


def _kallen(a: float, b: float, c: float) -> float:
    return float(a * a + b * b + c * c - 2.0 * (a * b + a * c + b * c))


RARE_KAON_KL_PI0_EMU_MODEL_V1 = "rare_kaon_klong_pi0_emu_lfv_form_factor_proxy_v1"
RARE_KAON_KL_PI0_EMU_INPUT_BUNDLE_V1 = (
    "rare_kaon_klong_pi0_emu_full_q2_neutral_form_factor_inputs_v1"
)
RARE_KAON_KL_PI0_EMU_FORM_FACTOR_MODEL_V1 = (
    "K0_to_pi0_fplus_fzero_single_pole_fplus0_FLAG2024_v1"
)
RARE_KAON_KL_PI0_EMU_PARAMETRIZATION_CITATION = (
    RARE_KAON_KTOPI_EMU_PARAMETRIZATION_CITATION
    + "; neutral K_L -> pi0 mode with standard neutral K_l3 f_+(q2)"
)
RARE_KAON_KL_PI0_EMU_PROXY_V1 = (
    "LEGACY-PROXY: K_L -> pi0 e+- mu-+ explicit-spurion compatibility path. "
    "The Phase-4c K021 production path uses rs_semileptonic_wilsons.lfv_llqq; "
    "scalar/tensor operators, charge-orientation-specific loop effects, and "
    "full sd/ds CP-phase matching remain deferred."
)
RARE_KAON_KL_PI0_EMU_Q2_TREATMENT_V1 = (
    "full_kinematic_q2_short_distance_lfv_proxy_with_klong_to_pi0_form_factor"
)


def rare_kaon_klong_pi0_emu_default_inputs() -> RareKaonKToPiLFVInputs:
    """Return the neutral-mode ``K_L -> pi0 e+- mu-+`` LFV input bundle."""

    return RareKaonKToPiLFVInputs(
        input_bundle=RARE_KAON_KL_PI0_EMU_INPUT_BUNDLE_V1,
        form_factor=RareKaonKToPiFormFactorInputs(
            model_label=RARE_KAON_KL_PI0_EMU_FORM_FACTOR_MODEL_V1,
            fplus_0=0.9698,
            fzero_0=0.9698,
            vector_pole_mass_gev=0.89166,
            scalar_pole_mass_gev=1.43,
            source=(
                "FLAG/PDG-era neutral K_l3 f_+(0)=0.9698 with the same "
                "single-pole K*(892) vector shape used by K020; K_L is "
                "treated with the standard neutral-mode semileptonic "
                "normalization for the K->pi vector current."
            ),
        ),
        kplus_mass_gev=0.497611,
        piplus_mass_gev=0.1349768,
        kplus_lifetime_s=5.116e-8,
        charge_state_factor=2.0,
        constants_citation=(
            "PDG-era K_L/K0 mass, pi0 mass, and K_L lifetime; neutral K_l3 "
            "f_+(0) normalization; G_F, alpha, sin2thetaW, and CKM from the "
            "shared rare_kaon_dilepton input bundle. charge_state_factor=2 "
            "matches the K021 YAML limit summed over e+- mu-+ charge states."
        ),
    )


def _neutral_klong_pi0_result(
    result: RareKaonKToPiLFVBranchingResult,
    *,
    inputs: RareKaonKToPiLFVInputs,
) -> RareKaonKToPiLFVBranchingResult:
    diagnostics = dict(result.diagnostics)
    diagnostics.update(
        {
            "base_model_label": RARE_KAON_KTOPI_EMU_MODEL_V1,
            "model_label": RARE_KAON_KL_PI0_EMU_MODEL_V1,
            "input_bundle": RARE_KAON_KL_PI0_EMU_INPUT_BUNDLE_V1,
            "matching_assumption": (
                RARE_KAON_KL_PI0_EMU_PROXY_V1
                if result.wilsons is None
                else result.wilsons.matching_assumption
            ),
            "parametrization_citation": RARE_KAON_KL_PI0_EMU_PARAMETRIZATION_CITATION,
            "q2_treatment": RARE_KAON_KL_PI0_EMU_Q2_TREATMENT_V1,
            "sm_lfv_policy": (
                "K_L -> pi0 e mu is charged-LFV and has zero SM rate for "
                "catalog purposes."
            ),
            "neutral_mode": True,
            "klong_cp_eigenstate_proxy": (
                "K_L is treated with the standard neutral K->pi semileptonic "
                "form-factor normalization; rigorous sd/ds CP-phase matching "
                "is not available on ParameterPoint."
            ),
            "charge_conjugate_modes_included": True,
            "charge_state_factor": float(inputs.charge_state_factor),
            "summed_charge_states": "e+- mu-+",
            "klong_mass_gev": float(inputs.kplus_mass_gev),
            "pi0_mass_gev": float(inputs.piplus_mass_gev),
            "klong_lifetime_s": float(inputs.kplus_lifetime_s),
            "legacy_input_field_names": (
                "RareKaonKToPiLFVInputs keeps kplus_mass_gev, piplus_mass_gev, "
                "and kplus_lifetime_s field names for K020 compatibility; "
                "K021 fills them with K_L, pi0, and K_L lifetime values."
            ),
        }
    )
    return RareKaonKToPiLFVBranchingResult(
        model_label=RARE_KAON_KL_PI0_EMU_MODEL_V1,
        input_bundle=inputs.input_bundle,
        charge_mode=result.charge_mode,
        branching_fraction=float(result.branching_fraction),
        sm_branching_fraction=0.0,
        np_shift_branching_fraction=float(result.np_shift_branching_fraction),
        q2_min_gev2=float(result.q2_min_gev2),
        q2_max_gev2=float(result.q2_max_gev2),
        electron_mass_gev=float(result.electron_mass_gev),
        muon_mass_gev=float(result.muon_mass_gev),
        y_vector_lfv=complex(result.y_vector_lfv),
        y_axial_lfv=complex(result.y_axial_lfv),
        lambda_wolfenstein=float(result.lambda_wolfenstein),
        wilsons=result.wilsons,
        diagnostics=diagnostics,
    )


def klong_pi0_emu_sm(
    inputs: RareKaonKToPiLFVInputs | None = None,
    *,
    charge_mode: str = "muplus_eminus",
) -> RareKaonKToPiLFVBranchingResult:
    """Return the catalog SM-limit ``K_L -> pi0 e+- mu-+`` rate, zero for LFV."""

    p = rare_kaon_klong_pi0_emu_default_inputs() if inputs is None else inputs
    return _neutral_klong_pi0_result(
        _branching_from_wilsons(None, inputs=p, charge_mode=charge_mode),
        inputs=p,
    )


def klong_pi0_emu_from_couplings(
    quark_couplings: QuarkMassBasisCouplings | RareKaonKToPiLFVWilsonCoefficients,
    lepton_couplings: object | None = None,
    *,
    m_kk_gev: float | None = None,
    inputs: RareKaonKToPiLFVInputs | None = None,
    charge_mode: str = "muplus_eminus",
) -> RareKaonKToPiLFVBranchingResult:
    """Evaluate charge-summed ``BR_SD(K_L -> pi0 e+- mu-+)`` from proxy couplings."""

    p = rare_kaon_klong_pi0_emu_default_inputs() if inputs is None else inputs
    if isinstance(quark_couplings, RareKaonKToPiLFVWilsonCoefficients):
        wilsons = quark_couplings
    else:
        if lepton_couplings is None:
            raise TypeError("lepton_couplings is required for K_L -> pi0 e mu LFV matching")
        wilsons = rare_kaon_ktopi_emu_wilsons_from_couplings(
            quark_couplings,
            lepton_couplings,
            m_kk_gev=m_kk_gev,
            inputs=p,
        )
    return _neutral_klong_pi0_result(
        _branching_from_wilsons(wilsons, inputs=p, charge_mode=charge_mode),
        inputs=p,
    )


def klong_pi0_emu_from_rs_semileptonic_wilsons(
    source: RSSemileptonicWilsonBundle,
    *,
    matching_scale_gev: float | None = None,
    inputs: RareKaonKToPiLFVInputs | None = None,
    charge_mode: str = "muplus_eminus",
    lepton_pair: str = "e_mu",
) -> RareKaonKToPiLFVBranchingResult:
    """Evaluate ``K_L -> pi0 e mu`` from Phase-4a LFV llqq Wilsons."""

    p = rare_kaon_klong_pi0_emu_default_inputs() if inputs is None else inputs
    base = kplus_piplus_emu_from_rs_semileptonic_wilsons(
        source,
        matching_scale_gev=matching_scale_gev,
        inputs=p,
        charge_mode=charge_mode,
        lepton_pair=lepton_pair,
    )
    return _neutral_klong_pi0_result(base, inputs=p)


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
    "RARE_KAON_KL_PI0_EMU_MODEL_V1",
    "RARE_KAON_KL_PI0_EMU_INPUT_BUNDLE_V1",
    "RARE_KAON_KL_PI0_EMU_FORM_FACTOR_MODEL_V1",
    "RARE_KAON_KL_PI0_EMU_PARAMETRIZATION_CITATION",
    "RARE_KAON_KL_PI0_EMU_PROXY_V1",
    "RARE_KAON_KL_PI0_EMU_Q2_TREATMENT_V1",
    "RARE_KAON_LFV_DILEPTON_PROXY_V1",
    "RARE_KAON_LFV_TREE_LEVEL_NOTE_V1",
    "rare_kaon_ktopi_emu_default_inputs",
    "rare_kaon_klong_pi0_emu_default_inputs",
    "rare_kaon_ktopi_emu_proxy_input",
    "rare_kaon_ktopi_emu_lepton_coupling_proxy",
    "rare_kaon_ktopi_fplus",
    "rare_kaon_ktopi_fzero",
    "rare_kaon_ktopi_emu_q2_range",
    "rare_kaon_ktopi_emu_wilsons_from_couplings",
    "rare_kaon_ktopi_emu_wilsons_from_rs_semileptonic",
    "rare_kaon_ktopi_emu_differential_branching_fraction",
    "kplus_piplus_emu_sm",
    "kplus_piplus_emu_from_couplings",
    "kplus_piplus_emu_from_rs_semileptonic_wilsons",
    "klong_pi0_emu_sm",
    "klong_pi0_emu_from_couplings",
    "klong_pi0_emu_from_rs_semileptonic_wilsons",
]
