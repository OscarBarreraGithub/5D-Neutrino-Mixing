"""LFV ``K_L -> e mu`` wrapper over rare-kaon dilepton machinery.

The production Phase-4c path reads the rigorous off-diagonal ``e mu`` block
from ``rs_semileptonic_wilsons.lfv_llqq`` and maps its C9/C10 inputs into the
rare-kaon Y normalization without reusing the old one-boson proxy.  The older
explicit-spurion helpers remain as legacy compatibility entry points only.

The LFV two-body rate is normalized to the K006 ``K_L -> mu+mu-`` convention.
For equal muon masses and a single axial coefficient it reduces to
``kappa_mu * |Y_NP / lambda^5|^2``.  For the charge-summed unequal-lepton mode,

    BR(K_L -> e+- mu-+) = 2 kappa_mu / lambda^10
        * sqrt(lambda_e_mu) / sqrt(1 - 4 r_mu^2)
        * { [1 - (r_e + r_mu)^2] |(r_mu - r_e) Y_V|^2
          + [1 - (r_mu - r_e)^2] |(r_e + r_mu) Y_A|^2 }
        / (4 r_mu^2).

With the current diagonal charged-lepton fit, the Phase-4a LFV block is
rigorously zero at tree level, so the evaluated LFV rate is zero and
non-vetoing.  Nonzero tree-level LFV requires non-diagonal lepton structure;
loop-induced LFV is deferred.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
import math
from typing import Any, Mapping

from quarkConstraints.deltaf2 import M_K

from .rare_kaon_dilepton import (
    QuarkMassBasisCouplings,
    RARE_KAON_DILEPTON_INPUT_BUNDLE_V1,
    RARE_KAON_DILEPTON_MODEL_V1,
    RARE_KAON_DILEPTON_OPERATOR_CONVENTION,
    RARE_KAON_DILEPTON_PARAMETRIZATION_CITATION,
    RARE_KAON_DILEPTON_RS_MATCHING_ASSUMPTION_V1,
    RareKaonDileptonSMInputs,
    RareKaonDileptonWilsonCoefficients,
    rare_kaon_dilepton_ckm_factors,
    rare_kaon_dilepton_default_sm_inputs,
    rare_kaon_dilepton_g_sm_squared,
    rare_kaon_dilepton_kappa_mu,
    rare_kaon_dilepton_wilsons_from_couplings,
)
from quarkConstraints.rs_semileptonic_wilsons import (
    RSLFVSemileptonicWilsonCoefficients,
    RSSemileptonicWilsonBundle,
)

RARE_KAON_LFV_DILEPTON_MODEL_V1 = "rare_kaon_lfv_klong_emu_proxy_v1"
RARE_KAON_LFV_DILEPTON_OPERATOR_CONVENTION = (
    RARE_KAON_DILEPTON_OPERATOR_CONVENTION
    + "; LFV extension uses e-mu lepton bilinears and the charge-summed "
    "K_L -> e+- mu-+ branching fraction"
)
RARE_KAON_LFV_DILEPTON_PROXY_V1 = (
    "LEGACY-PROXY: caller-supplied e-mu lepton overlap spurions are mapped "
    "to Z-like LFV couplings. The Phase-4c K019/K020/K021 production path "
    "does not use this proxy; it reads rs_semileptonic_wilsons.lfv_llqq."
)
RARE_KAON_LFV_TREE_LEVEL_NOTE_V1 = (
    "tree-level LFV rigorous from Phase-4a lfv_llqq light-Z contacts "
    "(=0 for the diagonal charged-lepton fit); nonzero only with "
    "non-diagonal lepton structure / loop-induced LFV deferred"
)
RARE_KAON_LFV_RS_MATCHING_STATUS_V1 = (
    "rs_semileptonic_lfv_llqq_additive_no_wilson_prefactor_reuse_no_second_"
    "1_over_M_KK_squared"
)
RARE_KAON_LFV_DILEPTON_PARAMETRIZATION_CITATION = (
    RARE_KAON_DILEPTON_PARAMETRIZATION_CITATION
    + "; unequal-lepton pseudoscalar two-body phase-space normalization"
)
RARE_KAON_LFV_DILEPTON_INPUT_BUNDLE_V1 = (
    "rare_kaon_lfv_klong_emu_inputs_k006_normalized_v1"
)
_DEFAULT_CHARGE_STATE_FACTOR = 2.0


@dataclass(frozen=True)
class RareKaonLFVDileptonSMInputs:
    """Numerical inputs for the K019 LFV two-body proxy."""

    input_bundle: str = RARE_KAON_LFV_DILEPTON_INPUT_BUNDLE_V1
    rare_kaon: RareKaonDileptonSMInputs = field(
        default_factory=rare_kaon_dilepton_default_sm_inputs
    )
    kaon_mass_gev: float = M_K
    electron_mass_gev: float = 0.00051099895
    muon_mass_gev: float = 0.1056583745
    charge_state_factor: float = _DEFAULT_CHARGE_STATE_FACTOR
    constants_citation: str = (
        "PDG charged-lepton masses and K0 mass mirrored from "
        "quarkConstraints.deltaf2; K006 kappa_mu normalization"
    )

    def __post_init__(self) -> None:
        for name in (
            "kaon_mass_gev",
            "electron_mass_gev",
            "muon_mass_gev",
            "charge_state_factor",
        ):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be positive and finite")


@dataclass(frozen=True)
class RareKaonLFVLeptonProxyInput:
    """Explicit e-mu lepton-overlap proxy for LFV neutral currents."""

    left_emu_overlap: complex
    right_emu_overlap: complex
    m_kk_gev: float
    source: str = "caller-supplied rare-kaon e-mu lepton proxy"


@dataclass(frozen=True)
class RareKaonLFVLeptonCouplingProxy:
    """Mapped Z-like e-mu lepton coupling used by the K019 proxy."""

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
class RareKaonLFVWilsonCoefficients:
    """LFV ``s -> d e mu`` Wilson proxy from quark and lepton spurions."""

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
    lepton_left_delta_emu: complex
    lepton_right_delta_emu: complex
    lepton_vector_delta_emu: complex
    lepton_axial_delta_emu: complex
    y_vector_left: complex
    y_vector_right: complex
    y_axial_left: complex
    y_axial_right: complex
    base_same_flavor_wilsons: RareKaonDileptonWilsonCoefficients | None

    @property
    def y_vector_lfv(self) -> complex:
        """Vector LFV combination entering the two-body rate."""

        return complex(self.y_vector_left - self.y_vector_right)

    @property
    def y_axial_lfv(self) -> complex:
        """Axial LFV combination entering the two-body rate."""

        return complex(self.y_axial_left - self.y_axial_right)

    @property
    def wilsons(self) -> Mapping[str, complex]:
        return {
            "YV_LFV_L": complex(self.y_vector_left),
            "YV_LFV_R": complex(self.y_vector_right),
            "YA_LFV_L": complex(self.y_axial_left),
            "YA_LFV_R": complex(self.y_axial_right),
            "YV_LFV_total": self.y_vector_lfv,
            "YA_LFV_total": self.y_axial_lfv,
        }


@dataclass(frozen=True)
class RareKaonLFVBranchingResult:
    """Short-distance ``K_L -> e+- mu-+`` branching-fraction prediction."""

    model_label: str
    input_bundle: str
    branching_fraction: float
    sm_branching_fraction: float
    np_shift_branching_fraction: float
    y_vector_lfv: complex
    y_axial_lfv: complex
    lambda_wolfenstein: float
    electron_mass_gev: float
    muon_mass_gev: float
    phase_space_lambda_sqrt: float
    charge_state_factor: float
    wilsons: RareKaonLFVWilsonCoefficients | None = None
    diagnostics: Mapping[str, Any] = field(default_factory=dict)


def rare_kaon_lfv_default_sm_inputs() -> RareKaonLFVDileptonSMInputs:
    """Return the K019 LFV input bundle."""

    return RareKaonLFVDileptonSMInputs()


def rare_kaon_lfv_proxy_input(
    left_emu_overlap: complex,
    right_emu_overlap: complex,
    m_kk_gev: float,
    *,
    source: str = "caller-supplied rare-kaon e-mu lepton proxy",
) -> RareKaonLFVLeptonProxyInput:
    """Build a shape-checked LFV lepton proxy input."""

    return RareKaonLFVLeptonProxyInput(
        left_emu_overlap=_finite_complex(left_emu_overlap, "left_emu_overlap"),
        right_emu_overlap=_finite_complex(right_emu_overlap, "right_emu_overlap"),
        m_kk_gev=_positive_float(m_kk_gev, "m_kk_gev"),
        source=str(source),
    )


def rare_kaon_lfv_lepton_coupling_proxy(
    source: Any,
    *,
    m_kk_gev: float | None = None,
    inputs: RareKaonLFVDileptonSMInputs | None = None,
) -> RareKaonLFVLeptonCouplingProxy:
    """Return the documented Z-like e-mu lepton-coupling proxy."""

    p = rare_kaon_lfv_default_sm_inputs() if inputs is None else inputs
    resolved_m_kk = None if m_kk_gev is None else _positive_float(m_kk_gev, "m_kk_gev")
    proxy_input = _coerce_proxy_input(source, m_kk_gev=resolved_m_kk)
    if resolved_m_kk is None:
        resolved_m_kk = proxy_input.m_kk_gev
    g_z = _weak_z_coupling(p.rare_kaon)
    left = complex(g_z * proxy_input.left_emu_overlap)
    right = complex(g_z * proxy_input.right_emu_overlap)
    vector = complex(left + right)
    axial = complex(right - left)
    return RareKaonLFVLeptonCouplingProxy(
        model_label=RARE_KAON_LFV_DILEPTON_MODEL_V1,
        matching_assumption=RARE_KAON_LFV_DILEPTON_PROXY_V1,
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
            "matching_assumption": RARE_KAON_LFV_DILEPTON_PROXY_V1,
        },
    )


def rare_kaon_lfv_wilsons_from_couplings(
    quark_couplings: QuarkMassBasisCouplings,
    lepton_couplings: Any,
    *,
    m_kk_gev: float | None = None,
    inputs: RareKaonLFVDileptonSMInputs | None = None,
) -> RareKaonLFVWilsonCoefficients:
    """Return the v1 ``s -> d e mu`` Wilson proxy."""

    p = rare_kaon_lfv_default_sm_inputs() if inputs is None else inputs
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
    vector_scale = _safe_ratio(
        proxy.lepton_vector_delta_emu,
        base.muon_axial_delta,
        "muon_axial_delta",
    )
    axial_scale = _safe_ratio(
        proxy.lepton_axial_delta_emu,
        base.muon_axial_delta,
        "muon_axial_delta",
    )
    return RareKaonLFVWilsonCoefficients(
        model_label=RARE_KAON_LFV_DILEPTON_MODEL_V1,
        base_model_label=RARE_KAON_DILEPTON_MODEL_V1,
        operator_convention=RARE_KAON_LFV_DILEPTON_OPERATOR_CONVENTION,
        matching_assumption=RARE_KAON_LFV_DILEPTON_PROXY_V1,
        M_KK=float(base.M_KK),
        matching_scale=float(base.matching_scale),
        left_sd_coupling=complex(base.left_sd_coupling),
        right_sd_coupling=complex(base.right_sd_coupling),
        left_sd_overlap=complex(base.left_sd_overlap),
        right_sd_overlap=complex(base.right_sd_overlap),
        left_quark_delta=complex(base.left_quark_delta),
        right_quark_delta=complex(base.right_quark_delta),
        lepton_left_delta_emu=complex(proxy.lepton_left_delta_emu),
        lepton_right_delta_emu=complex(proxy.lepton_right_delta_emu),
        lepton_vector_delta_emu=complex(proxy.lepton_vector_delta_emu),
        lepton_axial_delta_emu=complex(proxy.lepton_axial_delta_emu),
        y_vector_left=complex(base.y_np_left * vector_scale),
        y_vector_right=complex(base.y_np_right * vector_scale),
        y_axial_left=complex(base.y_np_left * axial_scale),
        y_axial_right=complex(base.y_np_right * axial_scale),
        base_same_flavor_wilsons=base,
    )


def rare_kaon_lfv_coeff_from_rs_semileptonic(
    source: RSSemileptonicWilsonBundle,
    *,
    transition: str = "s_to_d",
    lepton_pair: str = "e_mu",
) -> RSLFVSemileptonicWilsonCoefficients:
    """Return the Phase-4a LFV ``s -> d e mu`` Wilson block."""

    try:
        coeff = source.lfv_llqq[transition][lepton_pair]
    except (AttributeError, KeyError, TypeError) as exc:
        raise ValueError(
            "rs_semileptonic_wilsons.lfv_llqq"
            f"[{transition!r}][{lepton_pair!r}] is not available"
        ) from exc
    if coeff.transition_key != "s_d":
        raise ValueError(
            f"lfv_llqq[{transition!r}][{lepton_pair!r}] "
            f"transition_key={coeff.transition_key!r}, expected 's_d'"
        )
    if coeff.lepton_pair_key != lepton_pair:
        raise ValueError(
            f"LFV lepton pair {coeff.lepton_pair_key!r} does not match "
            f"{lepton_pair!r}"
        )
    return coeff


def rare_kaon_lfv_rs_semileptonic_diagnostics(
    coeff: RSLFVSemileptonicWilsonCoefficients,
    *,
    source: RSSemileptonicWilsonBundle | None = None,
    inputs: RareKaonLFVDileptonSMInputs | None = None,
) -> dict[str, Any]:
    """Diagnostics common to the rigorous LFV rare-kaon rewire."""

    p = rare_kaon_lfv_default_sm_inputs() if inputs is None else inputs
    return {
        "rs_semileptonic_wilsons_present": True,
        "rs_semileptonic_model_label": coeff.model_label,
        "rs_semileptonic_operator_convention": coeff.operator_convention,
        "rs_semileptonic_matching_assumption": coeff.matching_assumption,
        "rs_semileptonic_matching_status": RARE_KAON_LFV_RS_MATCHING_STATUS_V1,
        "tree_level_matching_status": (
            "rigorous_tree_light_z_lfv_llqq_from_rs_semileptonic_wilsons"
        ),
        "lfv_tree_level_note": RARE_KAON_LFV_TREE_LEVEL_NOTE_V1,
        "loop_lfv_status": "loop_induced_lfv_deferred",
        "rs_semileptonic_transition_key": coeff.transition_key,
        "rs_semileptonic_lepton_pair_key": coeff.lepton_pair_key,
        "rs_semileptonic_quark_sector": coeff.quark_sector,
        "rs_semileptonic_final_quark_index": int(coeff.final_quark_index),
        "rs_semileptonic_initial_quark_index": int(coeff.initial_quark_index),
        "rs_semileptonic_final_lepton_index": int(coeff.final_lepton_index),
        "rs_semileptonic_initial_lepton_index": int(coeff.initial_lepton_index),
        "rs_semileptonic_lambda_ckm_name": coeff.lambda_ckm_name,
        "rs_semileptonic_lambda_ckm": complex(coeff.lambda_ckm),
        "rs_semileptonic_contact_units": coeff.contact_units,
        "rs_semileptonic_contacts": {
            key: complex(value) for key, value in coeff.contacts.items()
        },
        "rs_semileptonic_wilson_coefficients": {
            key: complex(value) for key, value in coeff.wilsons.items()
        },
        "kaon_c9_c10_to_y_normalization": _kaon_lfv_c9_c10_to_y_norm(p),
        "lfv_contact_factorized": False,
        "wilson_prefactor_reused": False,
        "second_mkk_suppression_applied": False,
        "includes_heavy_neutral_exchange": (
            None if source is None else bool(source.includes_heavy_neutral_exchange)
        ),
        "includes_heavy_neutral_lepton": (
            None if source is None else bool(source.includes_heavy_neutral_lepton)
        ),
    }


def rare_kaon_lfv_wilsons_from_rs_semileptonic(
    source: RSSemileptonicWilsonBundle,
    *,
    matching_scale_gev: float | None = None,
    inputs: RareKaonLFVDileptonSMInputs | None = None,
    lepton_pair: str = "e_mu",
) -> RareKaonLFVWilsonCoefficients:
    """Map Phase-4a LFV C9/C10 Wilsons into the K019 Y inputs."""

    p = rare_kaon_lfv_default_sm_inputs() if inputs is None else inputs
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
    left_axial_contact = complex(coeff.contact_LR - coeff.contact_LL)
    right_axial_contact = complex(coeff.contact_RR - coeff.contact_RL)
    return RareKaonLFVWilsonCoefficients(
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
        lepton_left_delta_emu=complex(coeff.contact_LL + coeff.contact_RL),
        lepton_right_delta_emu=complex(coeff.contact_LR + coeff.contact_RR),
        lepton_vector_delta_emu=complex(
            coeff.contact_LL + coeff.contact_LR + coeff.contact_RL + coeff.contact_RR
        ),
        lepton_axial_delta_emu=complex(
            coeff.contact_LR + coeff.contact_RR - coeff.contact_LL - coeff.contact_RL
        ),
        y_vector_left=complex(-lambda_ckm * norm * complex(coeff.c9_lfv_np)),
        y_vector_right=complex(-lambda_ckm * norm * complex(coeff.c9p_lfv_np)),
        y_axial_left=complex(-lambda_ckm * norm * complex(coeff.c10_lfv_np)),
        y_axial_right=complex(-lambda_ckm * norm * complex(coeff.c10p_lfv_np)),
        base_same_flavor_wilsons=None,
    )


def klong_emu_from_rs_semileptonic_wilsons(
    source: RSSemileptonicWilsonBundle,
    *,
    matching_scale_gev: float | None = None,
    inputs: RareKaonLFVDileptonSMInputs | None = None,
    lepton_pair: str = "e_mu",
) -> RareKaonLFVBranchingResult:
    """Evaluate ``K_L -> e mu`` from Phase-4a LFV llqq Wilsons."""

    p = rare_kaon_lfv_default_sm_inputs() if inputs is None else inputs
    coeff = rare_kaon_lfv_coeff_from_rs_semileptonic(
        source,
        transition="s_to_d",
        lepton_pair=lepton_pair,
    )
    wilsons = rare_kaon_lfv_wilsons_from_rs_semileptonic(
        source,
        matching_scale_gev=matching_scale_gev,
        inputs=p,
        lepton_pair=lepton_pair,
    )
    result = _branching_from_wilsons(wilsons, inputs=p)
    diagnostics = dict(result.diagnostics)
    diagnostics.update(
        rare_kaon_lfv_rs_semileptonic_diagnostics(coeff, source=source, inputs=p)
    )
    diagnostics["base_matching_assumption"] = coeff.matching_assumption
    diagnostics["matching_assumption"] = coeff.matching_assumption
    diagnostics["s_to_d_lfv_llqq_rs_semileptonic_rewired"] = True
    diagnostics["rare_kaon_lfv_proxy_reused"] = False
    diagnostics["c9_c10_to_rare_kaon_effective_inputs"] = True
    return replace(result, diagnostics=diagnostics)


def rare_kaon_lfv_sm_branching_fraction(
    inputs: RareKaonLFVDileptonSMInputs | None = None,
) -> RareKaonLFVBranchingResult:
    """Return the catalog SM-limit ``K_L -> e mu`` rate, zero for LFV."""

    p = rare_kaon_lfv_default_sm_inputs() if inputs is None else inputs
    return _branching_from_wilsons(None, inputs=p)


def klong_emu_from_couplings(
    quark_couplings: QuarkMassBasisCouplings | RareKaonLFVWilsonCoefficients,
    lepton_couplings: Any | None = None,
    *,
    m_kk_gev: float | None = None,
    inputs: RareKaonLFVDileptonSMInputs | None = None,
) -> RareKaonLFVBranchingResult:
    """Evaluate charge-summed short-distance ``BR(K_L -> e+- mu-+)``."""

    p = rare_kaon_lfv_default_sm_inputs() if inputs is None else inputs
    if isinstance(quark_couplings, RareKaonLFVWilsonCoefficients):
        wilsons = quark_couplings
    else:
        if lepton_couplings is None:
            raise TypeError("lepton_couplings is required for rare-kaon LFV matching")
        wilsons = rare_kaon_lfv_wilsons_from_couplings(
            quark_couplings,
            lepton_couplings,
            m_kk_gev=m_kk_gev,
            inputs=p,
        )
    return _branching_from_wilsons(wilsons, inputs=p)


def _branching_from_wilsons(
    wilsons: RareKaonLFVWilsonCoefficients | None,
    *,
    inputs: RareKaonLFVDileptonSMInputs,
) -> RareKaonLFVBranchingResult:
    m_k = float(inputs.kaon_mass_gev)
    m_e = float(inputs.electron_mass_gev)
    m_mu = float(inputs.muon_mass_gev)
    r_e = m_e / m_k
    r_mu = m_mu / m_k
    phase_lambda = _kallen_dimensionless(r_e * r_e, r_mu * r_mu)
    if phase_lambda <= 0.0:
        raise ValueError("K_L -> e mu phase space is closed")
    sqrt_lambda = math.sqrt(phase_lambda)
    beta_mumu_sq = 1.0 - 4.0 * r_mu * r_mu
    if beta_mumu_sq <= 0.0:
        raise ValueError("K_L -> mu mu normalization phase space is closed")
    beta_mumu = math.sqrt(beta_mumu_sq)

    if wilsons is None:
        y_vector = 0.0j
        y_axial = 0.0j
    else:
        y_vector = wilsons.y_vector_lfv
        y_axial = wilsons.y_axial_lfv

    factors = rare_kaon_dilepton_ckm_factors(inputs.rare_kaon)
    lam = float(factors.lambda_wolfenstein)
    mass_sum = r_e + r_mu
    mass_diff = r_mu - r_e
    vector_term = (1.0 - mass_sum * mass_sum) * abs(mass_diff * y_vector) ** 2
    axial_term = (1.0 - mass_diff * mass_diff) * abs(mass_sum * y_axial) ** 2
    kappa = rare_kaon_dilepton_kappa_mu(inputs.rare_kaon)
    equal_muon_normalization = 4.0 * r_mu * r_mu
    branching = float(
        inputs.charge_state_factor
        * kappa
        * (sqrt_lambda / beta_mumu)
        * (vector_term + axial_term)
        / (equal_muon_normalization * lam**10)
    )

    diagnostics: dict[str, Any] = {
        "base_model_label": RARE_KAON_DILEPTON_MODEL_V1,
        "base_matching_assumption": RARE_KAON_DILEPTON_RS_MATCHING_ASSUMPTION_V1,
        "matching_assumption": (
            RARE_KAON_LFV_DILEPTON_PROXY_V1
            if wilsons is None
            else wilsons.matching_assumption
        ),
        "operator_convention": RARE_KAON_LFV_DILEPTON_OPERATOR_CONVENTION,
        "input_bundle": inputs.input_bundle,
        "base_input_bundle": RARE_KAON_DILEPTON_INPUT_BUNDLE_V1,
        "sm_branching_fraction": 0.0,
        "sm_lfv_policy": (
            "K_L -> e mu is charged-LFV and has zero SM rate for catalog purposes."
        ),
        "kaon_mass_gev": float(m_k),
        "electron_mass_gev": float(m_e),
        "muon_mass_gev": float(m_mu),
        "r_e": float(r_e),
        "r_mu": float(r_mu),
        "phase_space_lambda": float(phase_lambda),
        "phase_space_lambda_sqrt": float(sqrt_lambda),
        "k006_mumu_phase_space_beta": float(beta_mumu),
        "mass_sum_ratio": float(mass_sum),
        "mass_difference_ratio": float(mass_diff),
        "charge_state_factor": float(inputs.charge_state_factor),
        "vector_term": float(vector_term),
        "axial_term": float(axial_term),
        "kappa_mu": float(kappa),
        "lambda_wolfenstein": float(lam),
        "normalization_formula": (
            "charge_state_factor*kappa_mu*(sqrt(lambda_emu)/beta_mumu)"
            "*(vector_term+axial_term)/(4*r_mu^2*lambda^10)"
        ),
        "parametrization_citation": (
            RARE_KAON_LFV_DILEPTON_PARAMETRIZATION_CITATION
        ),
        "y_vector_lfv": complex(y_vector),
        "y_axial_lfv": complex(y_axial),
        "lambda_c": complex(factors.lambda_c),
        "lambda_t": complex(factors.lambda_t),
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
                "lepton_left_delta_emu": complex(wilsons.lepton_left_delta_emu),
                "lepton_right_delta_emu": complex(wilsons.lepton_right_delta_emu),
                "lepton_vector_delta_emu": complex(wilsons.lepton_vector_delta_emu),
                "lepton_axial_delta_emu": complex(wilsons.lepton_axial_delta_emu),
                "y_vector_left": complex(wilsons.y_vector_left),
                "y_vector_right": complex(wilsons.y_vector_right),
                "y_axial_left": complex(wilsons.y_axial_left),
                "y_axial_right": complex(wilsons.y_axial_right),
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
                    "base_y_np_total": complex(
                        wilsons.base_same_flavor_wilsons.y_np_total
                    ),
                }
            )

    return RareKaonLFVBranchingResult(
        model_label=RARE_KAON_LFV_DILEPTON_MODEL_V1,
        input_bundle=inputs.input_bundle,
        branching_fraction=branching,
        sm_branching_fraction=0.0,
        np_shift_branching_fraction=branching,
        y_vector_lfv=complex(y_vector),
        y_axial_lfv=complex(y_axial),
        lambda_wolfenstein=float(lam),
        electron_mass_gev=float(m_e),
        muon_mass_gev=float(m_mu),
        phase_space_lambda_sqrt=float(sqrt_lambda),
        charge_state_factor=float(inputs.charge_state_factor),
        wilsons=wilsons,
        diagnostics=diagnostics,
    )


def _weak_z_coupling(inputs: RareKaonDileptonSMInputs) -> float:
    g_weak = math.sqrt(4.0 * math.pi * inputs.alpha_em_mz / inputs.sin2_theta_w)
    cos_theta_w = math.sqrt(1.0 - inputs.sin2_theta_w)
    return float(g_weak / cos_theta_w)


def _diagnostic_matching_scale(matching_scale_gev: float | None) -> float:
    if matching_scale_gev is None:
        return 0.0
    return _positive_float(matching_scale_gev, "matching_scale_gev")


def _kaon_lfv_c9_c10_to_y_norm(inputs: RareKaonLFVDileptonSMInputs) -> float:
    return float(
        math.sqrt(2.0)
        * inputs.rare_kaon.gf_gev_minus2
        * inputs.rare_kaon.alpha_em_mz
        / (math.pi * rare_kaon_dilepton_g_sm_squared(inputs.rare_kaon))
    )


def _positive_float(value: object, name: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0.0:
        raise ValueError(f"{name} must be positive and finite")
    return number


def _finite_complex(value: object, name: str) -> complex:
    number = complex(value)
    if not math.isfinite(number.real) or not math.isfinite(number.imag):
        raise ValueError(f"{name} must be finite")
    return number


def _coerce_proxy_input(
    source: Any,
    *,
    m_kk_gev: float | None,
) -> RareKaonLFVLeptonProxyInput:
    if isinstance(source, RareKaonLFVLeptonProxyInput):
        if m_kk_gev is None:
            return source
        return replace(source, m_kk_gev=_positive_float(m_kk_gev, "m_kk_gev"))

    if isinstance(source, Mapping):
        left = source.get("left_emu_overlap", source.get("left_em_overlap"))
        right = source.get("right_emu_overlap", source.get("right_em_overlap"))
        source_label = source.get("source", "mapping rare-kaon e-mu lepton proxy")
        source_mkk = source.get("m_kk_gev", source.get("M_KK", m_kk_gev))
    else:
        left = getattr(source, "left_emu_overlap", None)
        right = getattr(source, "right_emu_overlap", None)
        source_label = getattr(
            source,
            "source",
            f"{type(source).__name__} rare-kaon e-mu lepton proxy",
        )
        source_mkk = getattr(source, "m_kk_gev", getattr(source, "M_KK", m_kk_gev))

    if left is None or right is None:
        raise TypeError(
            "lepton LFV proxy must provide left_emu_overlap and right_emu_overlap"
        )
    if source_mkk is None:
        raise TypeError("lepton LFV proxy needs m_kk_gev or an override")
    return rare_kaon_lfv_proxy_input(
        left,
        right,
        _positive_float(source_mkk, "m_kk_gev"),
        source=str(source_label),
    )


def _safe_ratio(numerator: complex, denominator: complex | float, name: str) -> complex:
    denom = complex(denominator)
    if abs(denom) <= 0.0:
        raise ValueError(f"{name} normalization is zero")
    return complex(numerator) / denom


def _kallen_dimensionless(x: float, y: float) -> float:
    return float(1.0 + x * x + y * y - 2.0 * (x + y + x * y))


__all__ = [
    "QuarkMassBasisCouplings",
    "RareKaonLFVDileptonSMInputs",
    "RareKaonLFVLeptonProxyInput",
    "RareKaonLFVLeptonCouplingProxy",
    "RareKaonLFVWilsonCoefficients",
    "RareKaonLFVBranchingResult",
    "RARE_KAON_LFV_DILEPTON_MODEL_V1",
    "RARE_KAON_LFV_DILEPTON_OPERATOR_CONVENTION",
    "RARE_KAON_LFV_DILEPTON_PROXY_V1",
    "RARE_KAON_LFV_TREE_LEVEL_NOTE_V1",
    "RARE_KAON_LFV_RS_MATCHING_STATUS_V1",
    "RARE_KAON_LFV_DILEPTON_PARAMETRIZATION_CITATION",
    "rare_kaon_lfv_default_sm_inputs",
    "rare_kaon_lfv_proxy_input",
    "rare_kaon_lfv_lepton_coupling_proxy",
    "rare_kaon_lfv_wilsons_from_couplings",
    "rare_kaon_lfv_coeff_from_rs_semileptonic",
    "rare_kaon_lfv_rs_semileptonic_diagnostics",
    "rare_kaon_lfv_wilsons_from_rs_semileptonic",
    "rare_kaon_lfv_sm_branching_fraction",
    "klong_emu_from_couplings",
    "klong_emu_from_rs_semileptonic_wilsons",
]
