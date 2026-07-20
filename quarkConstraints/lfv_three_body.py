"""Charged-LFV three-body lepton decays from dipole and contact amplitudes.

This module provides a small reusable low-energy observable for
``l_i -> 3 l_j`` modes such as ``mu -> 3e`` and, with different flavor labels,
``tau -> 3mu`` or ``tau -> 3e``.

The implemented v1 convention is

    BR(l_i -> 3 l_j) =
        BR(l_i -> l_j gamma) * alpha/(3 pi)
            * (log(m_i^2 / m_j^2) - 11/4)
        + BR(l_i -> l_j nu nubar)
            * [2 (|G_LL|^2 + |G_RR|^2) + |G_LR|^2 + |G_RL|^2
               + I_dipole-contact],

where ``G_AB = G_AB^Z + G_AB^box`` are dimensionless vector-contact
amplitudes in the Hamiltonian convention

    H_eff = 4 G_F/sqrt(2) G_AB
            (lbar_j gamma_mu P_A l_i)(lbar_j gamma^mu P_B l_j).

The Z-penguin proxy maps caller-supplied off-diagonal lepton overlap spurions
onto ``delta g_L/R(l_j l_i) = (m_Z/M_KK)^2 overlap_L/R`` and then integrates
out the SM Z with electron/tau/muon diagonal chiral couplings.  Box amplitudes
are accepted as explicit dimensionless ``G_AB^box`` inputs in the same
convention.

NEEDS-HUMAN-PHYSICS: a production RS prediction requires the EW KK/Z/Z' tower,
charged-lepton mass-basis neutral-current couplings, loop-level dipole
matching, and four-lepton box matching.  The v1 contact path is a documented
proxy for missing lepton-sector RS matching, not a full model calculation.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Mapping

import numpy as np

LFV_THREE_BODY_MODEL_V1 = "lfv_three_body_dipole_z_box_proxy_v1"
LFV_THREE_BODY_INPUT_BUNDLE_V1 = "lfv_three_body_sm_inputs_2026_v1"
LFV_THREE_BODY_OPERATOR_CONVENTION = (
    "BR(l_i -> 3 l_j) = BR(l_i -> l_j gamma) * alpha/(3*pi) * "
    "(log(m_i^2/m_j^2)-11/4) + BR(l_i -> l_j nu nubar) * "
    "[2(|G_LL|^2+|G_RR|^2)+|G_LR|^2+|G_RL|^2 "
    "+ I_dipole-contact]; "
    "G_AB are dimensionless vector-contact amplitudes normalized to "
    "4 G_F/sqrt(2)."
)
LFV_THREE_BODY_DIPOLE_CONTACT_INTERFERENCE_CONVENTION = (
    "Dipole-contact interference uses chiral dipole amplitudes A_L,A_R "
    "normalized by BR(l_i -> l_j gamma)=384*pi^2*(|A_L|^2+|A_R|^2) "
    "for muons and by 384*pi^2*BR(l_i -> l_j nu nubar) for taus.  "
    "Kuno-Okada give I_dc=-8*e*Re[A_R*(2G_LL+G_LR)^* + "
    "A_L*(2G_RR+G_RL)^*] for the width normalized to "
    "G_F^2*m_i^5/(192*pi^3).  If only the parent dipole branching fraction "
    "is known, the relative phase/chiral split is NEEDS-HUMAN-PHYSICS and "
    "the reported HARD prediction uses the constructive |8e| sign-envelope "
    "upper bound; the destructive lower envelope is kept in diagnostics."
)
LFV_THREE_BODY_PROXY_V1 = (
    "NEEDS-HUMAN-PHYSICS: off-diagonal charged-lepton neutral-current and "
    "four-lepton box couplings are not available on ParameterPoint; caller "
    "proxy inputs map lepton overlaps to Z-penguin contact amplitudes "
    "G_AB^Z = 2 delta_g_A(l_j l_i) g_B(l_j) with "
    "delta_g_A=(m_Z/M_KK)^2 overlap_A, while box_* inputs are explicit "
    "dimensionless G_AB^box amplitudes in the same low-energy convention."
)

_CHARGED_LEPTON_MASSES_GEV = {
    "e": 0.00051099895000,
    "mu": 0.1056583755,
    "tau": 1.77686,
}
_CHARGED_LEPTONS = frozenset(_CHARGED_LEPTON_MASSES_GEV)
TAU_TO_E_NUNU_BRANCHING_FRACTION = 0.1782
TAU_TO_MU_NUNU_BRANCHING_FRACTION = 0.1739
PDG_TAU_LEPTONIC_BRANCHING_FRACTION_SOURCE = (
    "PDG tau branching-fraction review: B(tau->e nu nubar)=0.1782 and "
    "B(tau->mu nu nubar)=0.1739."
)
_LEPTONIC_NORMALIZATION_BRANCHING_FRACTIONS = {
    ("mu", "e"): 1.0,
    ("tau", "e"): TAU_TO_E_NUNU_BRANCHING_FRACTION,
    ("tau", "mu"): TAU_TO_MU_NUNU_BRANCHING_FRACTION,
}


@dataclass(frozen=True)
class LFVThreeBodySMInputs:
    """Numerical inputs for the low-energy ``l_i -> 3 l_j`` formula."""

    input_bundle: str = LFV_THREE_BODY_INPUT_BUNDLE_V1
    alpha_em: float = 1.0 / 137.035999084
    sin2_theta_w: float = 0.23122
    m_z_gev: float = 91.1876
    proxy_strength: float = 1.0
    charged_lepton_masses_gev: Mapping[str, float] = field(
        default_factory=lambda: dict(_CHARGED_LEPTON_MASSES_GEV)
    )
    leptonic_normalization_branching_fractions: Mapping[tuple[str, str], float] = field(
        default_factory=lambda: dict(_LEPTONIC_NORMALIZATION_BRANCHING_FRACTIONS)
    )

    def charged_lepton_mass(self, flavor: str) -> float:
        key = _canonical_charged_lepton(flavor)
        return float(self.charged_lepton_masses_gev[key])

    def z_chiral_couplings(self, flavor: str) -> tuple[float, float]:
        """Return SM diagonal ``Z l l`` chiral couplings ``(g_L, g_R)``."""

        _canonical_charged_lepton(flavor)
        return float(-0.5 + self.sin2_theta_w), float(self.sin2_theta_w)

    def leptonic_normalization_branching_fraction(
        self,
        initial_flavor: str,
        final_flavor: str,
    ) -> float:
        """Return ``BR(l_i -> l_j nu nubar)`` converting KO widths to BRs."""

        key = (
            _canonical_charged_lepton(initial_flavor),
            _canonical_charged_lepton(final_flavor),
        )
        if key not in self.leptonic_normalization_branching_fractions:
            raise ValueError(
                "missing leptonic normalization branching fraction for "
                f"{key[0]} -> {key[1]} nu nubar"
            )
        number = float(self.leptonic_normalization_branching_fractions[key])
        if not math.isfinite(number) or not 0.0 < number <= 1.0:
            raise ValueError(
                "leptonic normalization branching fraction must lie in (0, 1]"
            )
        return number


@dataclass(frozen=True)
class LFVThreeBodyContactProxyInput:
    """Explicit proxy source for LFV Z-penguin and box contact amplitudes."""

    initial_flavor: str
    final_flavor: str
    left_lfv_overlap: complex
    right_lfv_overlap: complex
    m_kk_gev: float
    box_ll: complex = 0.0j
    box_lr: complex = 0.0j
    box_rl: complex = 0.0j
    box_rr: complex = 0.0j
    source: str = "caller-supplied lfv three-body contact proxy"


@dataclass(frozen=True)
class LFVThreeBodyContactAmplitudes:
    """Dimensionless vector-contact amplitudes for ``l_i -> 3 l_j``."""

    model_label: str
    matching_assumption: str
    operator_convention: str
    initial_flavor: str
    final_flavor: str
    M_KK: float
    matching_scale: float
    scale_factor: float
    left_lfv_overlap: complex
    right_lfv_overlap: complex
    delta_g_left_lfv: complex
    delta_g_right_lfv: complex
    z_ll: complex
    z_lr: complex
    z_rl: complex
    z_rr: complex
    box_ll: complex
    box_lr: complex
    box_rl: complex
    box_rr: complex
    total_ll: complex
    total_lr: complex
    total_rl: complex
    total_rr: complex
    source: str
    diagnostics: Mapping[str, Any]


@dataclass(frozen=True)
class LFVThreeBodyBranchingResult:
    """Branching-fraction prediction for ``l_i -> 3 l_j``."""

    model_label: str
    input_bundle: str
    initial_flavor: str
    final_flavor: str
    branching_fraction: float
    sm_branching_fraction: float
    np_shift_branching_fraction: float
    dipole_component: float
    z_penguin_component: float
    box_component: float
    z_box_interference_component: float
    contact_component: float
    dipole_conversion_factor: float
    dipole_parent_branching_fraction: float
    dipole_contact_interference_component: float
    dipole_contact_interference_lower: float
    dipole_contact_interference_upper: float
    dipole_contact_interference_treatment: str
    dipole_amplitude_left: complex | None
    dipole_amplitude_right: complex | None
    ratio_to_limit: float | None
    br_limit: float | None
    passes: bool | None
    contact_amplitudes: LFVThreeBodyContactAmplitudes
    diagnostics: Mapping[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class _DipoleContactInterference:
    component: float
    lower: float
    upper: float
    treatment: str
    amplitude_left: complex | None
    amplitude_right: complex | None
    amplitude_norm: float
    chiral_combination_left: complex
    chiral_combination_right: complex


def default_sm_inputs() -> LFVThreeBodySMInputs:
    """Return the default low-energy input bundle."""

    return LFVThreeBodySMInputs()


def lfv_three_body_proxy_input(
    left_lfv_overlap: complex,
    right_lfv_overlap: complex,
    m_kk_gev: float,
    *,
    initial_flavor: str = "mu",
    final_flavor: str = "e",
    box_ll: complex = 0.0j,
    box_lr: complex = 0.0j,
    box_rl: complex = 0.0j,
    box_rr: complex = 0.0j,
    source: str = "caller-supplied lfv three-body contact proxy",
) -> LFVThreeBodyContactProxyInput:
    """Build a shape-checked contact proxy input."""

    initial = _canonical_charged_lepton(initial_flavor)
    final = _canonical_charged_lepton(final_flavor)
    if initial == final:
        raise ValueError("LFV three-body decay requires distinct lepton flavors")
    return LFVThreeBodyContactProxyInput(
        initial_flavor=initial,
        final_flavor=final,
        left_lfv_overlap=_finite_complex(left_lfv_overlap, "left_lfv_overlap"),
        right_lfv_overlap=_finite_complex(right_lfv_overlap, "right_lfv_overlap"),
        m_kk_gev=_positive_float(m_kk_gev, "m_kk_gev"),
        box_ll=_finite_complex(box_ll, "box_ll"),
        box_lr=_finite_complex(box_lr, "box_lr"),
        box_rl=_finite_complex(box_rl, "box_rl"),
        box_rr=_finite_complex(box_rr, "box_rr"),
        source=str(source),
    )


def lfv_three_body_has_contact_proxy(source: Any) -> bool:
    """Return ``True`` when ``source`` carries contact-proxy information."""

    if isinstance(source, LFVThreeBodyContactProxyInput):
        return True
    if isinstance(source, Mapping):
        return any(key in source for key in _CONTACT_MAPPING_KEYS)
    return any(hasattr(source, name) for name in _CONTACT_ATTR_NAMES)


def lfv_three_body_contact_amplitudes(
    source: Any,
    *,
    initial_flavor: str = "mu",
    final_flavor: str = "e",
    m_kk_gev: float | None = None,
    inputs: LFVThreeBodySMInputs | None = None,
) -> LFVThreeBodyContactAmplitudes:
    """Map a contact proxy source onto low-energy vector amplitudes."""

    p = default_sm_inputs() if inputs is None else inputs
    proxy = _coerce_contact_proxy(
        source,
        initial_flavor=initial_flavor,
        final_flavor=final_flavor,
        m_kk_gev=m_kk_gev,
    )
    scale = float(p.proxy_strength * (p.m_z_gev / proxy.m_kk_gev) ** 2)
    delta_left = complex(scale * proxy.left_lfv_overlap)
    delta_right = complex(scale * proxy.right_lfv_overlap)
    g_left, g_right = p.z_chiral_couplings(proxy.final_flavor)
    z_ll = complex(2.0 * delta_left * g_left)
    z_lr = complex(2.0 * delta_left * g_right)
    z_rl = complex(2.0 * delta_right * g_left)
    z_rr = complex(2.0 * delta_right * g_right)

    total_ll = complex(z_ll + proxy.box_ll)
    total_lr = complex(z_lr + proxy.box_lr)
    total_rl = complex(z_rl + proxy.box_rl)
    total_rr = complex(z_rr + proxy.box_rr)
    diagnostics = {
        "m_kk_gev": float(proxy.m_kk_gev),
        "matching_scale_gev": float(proxy.m_kk_gev),
        "m_z_gev": float(p.m_z_gev),
        "sin2_theta_w": float(p.sin2_theta_w),
        "z_final_lepton_g_left": float(g_left),
        "z_final_lepton_g_right": float(g_right),
        "proxy_strength": float(p.proxy_strength),
        "scale_factor": float(scale),
        "left_lfv_overlap": complex(proxy.left_lfv_overlap),
        "right_lfv_overlap": complex(proxy.right_lfv_overlap),
        "delta_g_left_lfv": delta_left,
        "delta_g_right_lfv": delta_right,
        "z_ll": z_ll,
        "z_lr": z_lr,
        "z_rl": z_rl,
        "z_rr": z_rr,
        "box_ll": complex(proxy.box_ll),
        "box_lr": complex(proxy.box_lr),
        "box_rl": complex(proxy.box_rl),
        "box_rr": complex(proxy.box_rr),
        "total_ll": total_ll,
        "total_lr": total_lr,
        "total_rl": total_rl,
        "total_rr": total_rr,
        "proxy_source": proxy.source,
        "matching_assumption": LFV_THREE_BODY_PROXY_V1,
    }
    return LFVThreeBodyContactAmplitudes(
        model_label=LFV_THREE_BODY_MODEL_V1,
        matching_assumption=LFV_THREE_BODY_PROXY_V1,
        operator_convention=LFV_THREE_BODY_OPERATOR_CONVENTION,
        initial_flavor=proxy.initial_flavor,
        final_flavor=proxy.final_flavor,
        M_KK=float(proxy.m_kk_gev),
        matching_scale=float(proxy.m_kk_gev),
        scale_factor=float(scale),
        left_lfv_overlap=complex(proxy.left_lfv_overlap),
        right_lfv_overlap=complex(proxy.right_lfv_overlap),
        delta_g_left_lfv=delta_left,
        delta_g_right_lfv=delta_right,
        z_ll=z_ll,
        z_lr=z_lr,
        z_rl=z_rl,
        z_rr=z_rr,
        box_ll=complex(proxy.box_ll),
        box_lr=complex(proxy.box_lr),
        box_rl=complex(proxy.box_rl),
        box_rr=complex(proxy.box_rr),
        total_ll=total_ll,
        total_lr=total_lr,
        total_rl=total_rl,
        total_rr=total_rr,
        source=proxy.source,
        diagnostics=diagnostics,
    )


def lfv_three_body_from_components(
    *,
    dipole_parent_branching_fraction: float,
    contact_amplitudes: LFVThreeBodyContactAmplitudes,
    br_limit: float | None = None,
    inputs: LFVThreeBodySMInputs | None = None,
    dipole_amplitude_left: complex | None = None,
    dipole_amplitude_right: complex | None = None,
) -> LFVThreeBodyBranchingResult:
    """Evaluate ``BR(l_i -> 3 l_j)`` from dipole BR and contact amplitudes."""

    p = default_sm_inputs() if inputs is None else inputs
    dipole_br = _nonnegative_float(
        dipole_parent_branching_fraction,
        "dipole_parent_branching_fraction",
    )
    initial = contact_amplitudes.initial_flavor
    final = contact_amplitudes.final_flavor
    initial_mass = p.charged_lepton_mass(initial)
    final_mass = p.charged_lepton_mass(final)
    if initial_mass <= final_mass:
        raise ValueError("parent lepton must be heavier than final lepton")
    leptonic_normalization_br = p.leptonic_normalization_branching_fraction(
        initial,
        final,
    )

    log_enhancement = math.log((initial_mass / final_mass) ** 2) - 11.0 / 4.0
    dipole_factor = float(p.alpha_em / (3.0 * math.pi) * log_enhancement)
    if dipole_factor < 0.0:
        raise ValueError("dipole conversion factor is negative")

    dipole_component = float(dipole_br * dipole_factor)
    z_width_component = _vector_contact_rate(
        contact_amplitudes.z_ll,
        contact_amplitudes.z_lr,
        contact_amplitudes.z_rl,
        contact_amplitudes.z_rr,
    )
    box_width_component = _vector_contact_rate(
        contact_amplitudes.box_ll,
        contact_amplitudes.box_lr,
        contact_amplitudes.box_rl,
        contact_amplitudes.box_rr,
    )
    contact_width_component = _vector_contact_rate(
        contact_amplitudes.total_ll,
        contact_amplitudes.total_lr,
        contact_amplitudes.total_rl,
        contact_amplitudes.total_rr,
    )
    z_component = float(leptonic_normalization_br * z_width_component)
    box_component = float(leptonic_normalization_br * box_width_component)
    contact_component = float(leptonic_normalization_br * contact_width_component)
    dipole_contact = _dipole_contact_interference(
        dipole_parent_branching_fraction=dipole_br,
        total_ll=contact_amplitudes.total_ll,
        total_lr=contact_amplitudes.total_lr,
        total_rl=contact_amplitudes.total_rl,
        total_rr=contact_amplitudes.total_rr,
        alpha_em=p.alpha_em,
        leptonic_normalization_branching_fraction=leptonic_normalization_br,
        dipole_amplitude_left=dipole_amplitude_left,
        dipole_amplitude_right=dipole_amplitude_right,
    )
    interference = float(contact_component - z_component - box_component)
    branching = float(
        dipole_component + contact_component + dipole_contact.component
    )
    limit = None if br_limit is None else _bounded_probability(br_limit, "br_limit")
    ratio = None if limit is None else float(branching / limit)
    passes = None if ratio is None else bool(ratio <= 1.0)

    diagnostics = {
        "model_label": LFV_THREE_BODY_MODEL_V1,
        "input_bundle": p.input_bundle,
        "operator_convention": LFV_THREE_BODY_OPERATOR_CONVENTION,
        "matching_assumption": LFV_THREE_BODY_PROXY_V1,
        "sm_branching_fraction": 0.0,
        "sm_lfv_policy": (
            "Charged-LFV l_i -> 3 l_j has negligible SM rate for catalog "
            "purposes; the bound is applied to the pure NP prediction."
        ),
        "alpha_em": float(p.alpha_em),
        "initial_lepton_mass_gev": float(initial_mass),
        "final_lepton_mass_gev": float(final_mass),
        "leptonic_normalization_branching_fraction": float(
            leptonic_normalization_br
        ),
        "leptonic_normalization_source": (
            "muon normalization is unity; tau normalization uses "
            f"{PDG_TAU_LEPTONIC_BRANCHING_FRACTION_SOURCE}"
        ),
        "dipole_log_enhancement": float(log_enhancement),
        "dipole_conversion_factor": float(dipole_factor),
        "dipole_parent_branching_fraction": float(dipole_br),
        "dipole_component": float(dipole_component),
        "z_penguin_component": float(z_component),
        "box_component": float(box_component),
        "z_box_interference_component": float(interference),
        "contact_component": float(contact_component),
        "z_penguin_width_component": float(z_width_component),
        "box_width_component": float(box_width_component),
        "z_box_interference_width_component": float(
            contact_width_component - z_width_component - box_width_component
        ),
        "contact_width_component": float(contact_width_component),
        "dipole_contact_interference_component": float(dipole_contact.component),
        "dipole_contact_interference_lower": float(dipole_contact.lower),
        "dipole_contact_interference_upper": float(dipole_contact.upper),
        "dipole_contact_interference_treatment": dipole_contact.treatment,
        "dipole_contact_interference_convention": (
            LFV_THREE_BODY_DIPOLE_CONTACT_INTERFERENCE_CONVENTION
        ),
        "dipole_contact_interference_needs_human_physics": (
            "NEEDS-HUMAN-PHYSICS" in dipole_contact.treatment
        ),
        "dipole_chiral_amplitude_norm": float(dipole_contact.amplitude_norm),
        "dipole_amplitude_left": dipole_contact.amplitude_left,
        "dipole_amplitude_right": dipole_contact.amplitude_right,
        "dipole_contact_chiral_combination_left": (
            dipole_contact.chiral_combination_left
        ),
        "dipole_contact_chiral_combination_right": (
            dipole_contact.chiral_combination_right
        ),
        "branching_formula": (
            "BR = BR(l_i -> l_j gamma) * alpha/(3*pi) * "
            "(log(m_i^2/m_j^2)-11/4) + "
            "BR(l_i -> l_j nu nubar) * "
            "[2(|G_LL|^2+|G_RR|^2)+|G_LR|^2+|G_RL|^2 with "
            "G_AB=G_AB^Z+G_AB^box plus "
            "I_dc=-8*e*Re[A_R*(2G_LL+G_LR)^* + "
            "A_L*(2G_RR+G_RL)^*]] (Kuno-Okada)."
        ),
        **dict(contact_amplitudes.diagnostics),
    }
    return LFVThreeBodyBranchingResult(
        model_label=LFV_THREE_BODY_MODEL_V1,
        input_bundle=p.input_bundle,
        initial_flavor=initial,
        final_flavor=final,
        branching_fraction=branching,
        sm_branching_fraction=0.0,
        np_shift_branching_fraction=branching,
        dipole_component=float(dipole_component),
        z_penguin_component=float(z_component),
        box_component=float(box_component),
        z_box_interference_component=float(interference),
        contact_component=float(contact_component),
        dipole_conversion_factor=float(dipole_factor),
        dipole_parent_branching_fraction=float(dipole_br),
        dipole_contact_interference_component=float(dipole_contact.component),
        dipole_contact_interference_lower=float(dipole_contact.lower),
        dipole_contact_interference_upper=float(dipole_contact.upper),
        dipole_contact_interference_treatment=dipole_contact.treatment,
        dipole_amplitude_left=dipole_contact.amplitude_left,
        dipole_amplitude_right=dipole_contact.amplitude_right,
        ratio_to_limit=ratio,
        br_limit=limit,
        passes=passes,
        contact_amplitudes=contact_amplitudes,
        diagnostics=diagnostics,
    )


def _coerce_contact_proxy(
    source: Any,
    *,
    initial_flavor: str,
    final_flavor: str,
    m_kk_gev: float | None,
) -> LFVThreeBodyContactProxyInput:
    if isinstance(source, LFVThreeBodyContactProxyInput):
        return lfv_three_body_proxy_input(
            source.left_lfv_overlap,
            source.right_lfv_overlap,
            source.m_kk_gev if m_kk_gev is None else m_kk_gev,
            initial_flavor=source.initial_flavor,
            final_flavor=source.final_flavor,
            box_ll=source.box_ll,
            box_lr=source.box_lr,
            box_rl=source.box_rl,
            box_rr=source.box_rr,
            source=source.source,
        )
    if isinstance(source, Mapping):
        return _proxy_from_mapping(
            source,
            initial_flavor=initial_flavor,
            final_flavor=final_flavor,
            m_kk_gev=m_kk_gev,
        )
    return _proxy_from_object(
        source,
        initial_flavor=initial_flavor,
        final_flavor=final_flavor,
        m_kk_gev=m_kk_gev,
    )


def _proxy_from_mapping(
    mapping: Mapping[str, Any],
    *,
    initial_flavor: str,
    final_flavor: str,
    m_kk_gev: float | None,
) -> LFVThreeBodyContactProxyInput:
    initial = _canonical_charged_lepton(
        _first_present_key(mapping, ("initial_flavor", "parent_flavor")) or initial_flavor
    )
    final = _canonical_charged_lepton(
        _first_present_key(mapping, ("final_flavor", "daughter_flavor")) or final_flavor
    )
    left, right = _overlaps_from_mapping(mapping)
    return lfv_three_body_proxy_input(
        left,
        right,
        _resolve_m_kk_from_mapping(mapping, m_kk_gev),
        initial_flavor=initial,
        final_flavor=final,
        box_ll=_optional_complex_from_mapping(mapping, _BOX_LL_KEYS, "box_ll"),
        box_lr=_optional_complex_from_mapping(mapping, _BOX_LR_KEYS, "box_lr"),
        box_rl=_optional_complex_from_mapping(mapping, _BOX_RL_KEYS, "box_rl"),
        box_rr=_optional_complex_from_mapping(mapping, _BOX_RR_KEYS, "box_rr"),
        source=str(mapping.get("source", "mapping lfv three-body contact proxy")),
    )


def _proxy_from_object(
    value: Any,
    *,
    initial_flavor: str,
    final_flavor: str,
    m_kk_gev: float | None,
) -> LFVThreeBodyContactProxyInput:
    initial = _canonical_charged_lepton(
        _first_present_attr(value, ("initial_flavor", "parent_flavor")) or initial_flavor
    )
    final = _canonical_charged_lepton(
        _first_present_attr(value, ("final_flavor", "daughter_flavor")) or final_flavor
    )
    left = _first_present_attr(
        value,
        ("left_lfv_overlap", "left_emu_overlap", "left_emu"),
    )
    right = _first_present_attr(
        value,
        ("right_lfv_overlap", "right_emu_overlap", "right_emu"),
    )
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
        if left_matrix is not None:
            left = _matrix_offdiag(left_matrix, "left", final, initial)
        if right_matrix is not None:
            right = _matrix_offdiag(right_matrix, "right", final, initial)
    return lfv_three_body_proxy_input(
        0.0j if left is None else left,
        0.0j if right is None else right,
        _resolve_m_kk_from_object(value, m_kk_gev),
        initial_flavor=initial,
        final_flavor=final,
        box_ll=_optional_complex_from_object(value, _BOX_LL_KEYS, "box_ll"),
        box_lr=_optional_complex_from_object(value, _BOX_LR_KEYS, "box_lr"),
        box_rl=_optional_complex_from_object(value, _BOX_RL_KEYS, "box_rl"),
        box_rr=_optional_complex_from_object(value, _BOX_RR_KEYS, "box_rr"),
        source=str(getattr(value, "source", "object lfv three-body contact proxy")),
    )


def _overlaps_from_mapping(mapping: Mapping[str, Any]) -> tuple[complex, complex]:
    left = _first_present_key(
        mapping,
        ("left_lfv_overlap", "left_emu_overlap", "left_emu"),
    )
    right = _first_present_key(
        mapping,
        ("right_lfv_overlap", "right_emu_overlap", "right_emu"),
    )
    if left is not None or right is not None:
        return (
            _finite_complex(0.0j if left is None else left, "left_lfv_overlap"),
            _finite_complex(0.0j if right is None else right, "right_lfv_overlap"),
        )

    left_matrix = _first_present_key(
        mapping,
        ("left_charged_lepton_overlap", "left_lepton_overlap", "left_overlap"),
    )
    right_matrix = _first_present_key(
        mapping,
        ("right_charged_lepton_overlap", "right_lepton_overlap", "right_overlap"),
    )
    initial = _canonical_charged_lepton(
        _first_present_key(mapping, ("initial_flavor", "parent_flavor")) or "mu"
    )
    final = _canonical_charged_lepton(
        _first_present_key(mapping, ("final_flavor", "daughter_flavor")) or "e"
    )
    left_value = 0.0j
    right_value = 0.0j
    if left_matrix is not None:
        left_value = _matrix_offdiag(left_matrix, "left", final, initial)
    if right_matrix is not None:
        right_value = _matrix_offdiag(right_matrix, "right", final, initial)
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
    raise AttributeError("lfv three-body contact proxy must provide m_kk_gev or M_KK")


def _optional_complex_from_mapping(
    mapping: Mapping[str, Any],
    keys: tuple[str, ...],
    name: str,
) -> complex:
    value = _first_present_key(mapping, keys)
    return 0.0j if value is None else _finite_complex(value, name)


def _optional_complex_from_object(value: Any, keys: tuple[str, ...], name: str) -> complex:
    item = _first_present_attr(value, keys)
    return 0.0j if item is None else _finite_complex(item, name)


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


def _matrix_offdiag(value: Any, name: str, final: str, initial: str) -> complex:
    matrix = np.asarray(value, dtype=np.complex128)
    if matrix.shape != (3, 3):
        raise ValueError(f"{name} charged-lepton overlap matrix must have shape (3, 3)")
    if not np.all(np.isfinite(matrix.real)) or not np.all(np.isfinite(matrix.imag)):
        raise ValueError(f"{name} charged-lepton overlap matrix entries must be finite")
    return complex(matrix[_flavor_index(final), _flavor_index(initial)])


def _vector_contact_rate(ll: complex, lr: complex, rl: complex, rr: complex) -> float:
    return float(2.0 * (abs(ll) ** 2 + abs(rr) ** 2) + abs(lr) ** 2 + abs(rl) ** 2)


def _dipole_contact_interference(
    *,
    dipole_parent_branching_fraction: float,
    total_ll: complex,
    total_lr: complex,
    total_rl: complex,
    total_rr: complex,
    alpha_em: float,
    leptonic_normalization_branching_fraction: float,
    dipole_amplitude_left: complex | None,
    dipole_amplitude_right: complex | None,
) -> _DipoleContactInterference:
    chiral_left = complex(2.0 * total_rr + total_rl)
    chiral_right = complex(2.0 * total_ll + total_lr)
    leptonic_br = float(leptonic_normalization_branching_fraction)
    if not math.isfinite(leptonic_br) or not 0.0 < leptonic_br <= 1.0:
        raise ValueError("leptonic normalization branching fraction must lie in (0, 1]")
    dipole_norm = math.sqrt(
        float(dipole_parent_branching_fraction)
        / (384.0 * math.pi**2 * leptonic_br)
    )
    electric_charge = math.sqrt(4.0 * math.pi * float(alpha_em))
    # Kuno-Okada RMP 73 (2001), Eq. (2.14): -8e times this chiral pairing.
    coefficient = -8.0 * electric_charge
    has_explicit_amplitudes = (
        dipole_amplitude_left is not None or dipole_amplitude_right is not None
    )

    if has_explicit_amplitudes:
        left = _finite_complex(
            0.0j if dipole_amplitude_left is None else dipole_amplitude_left,
            "dipole_amplitude_left",
        )
        right = _finite_complex(
            0.0j if dipole_amplitude_right is None else dipole_amplitude_right,
            "dipole_amplitude_right",
        )
        component = float(
            leptonic_br
            * coefficient
            * (right * chiral_right.conjugate() + left * chiral_left.conjugate()).real
        )
        return _DipoleContactInterference(
            component=component,
            lower=component,
            upper=component,
            treatment="explicit_chiral_dipole_amplitudes",
            amplitude_left=left,
            amplitude_right=right,
            amplitude_norm=float(math.sqrt(abs(left) ** 2 + abs(right) ** 2)),
            chiral_combination_left=chiral_left,
            chiral_combination_right=chiral_right,
        )

    envelope = float(
        abs(coefficient)
        * leptonic_br
        * dipole_norm
        * math.sqrt(abs(chiral_left) ** 2 + abs(chiral_right) ** 2)
    )
    treatment = (
        "not_applicable_zero_dipole_or_contact"
        if envelope == 0.0
        else "constructive_sign_envelope_NEEDS-HUMAN-PHYSICS"
    )
    return _DipoleContactInterference(
        component=envelope,
        lower=-envelope,
        upper=envelope,
        treatment=treatment,
        amplitude_left=None,
        amplitude_right=None,
        amplitude_norm=float(dipole_norm),
        chiral_combination_left=chiral_left,
        chiral_combination_right=chiral_right,
    )


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
    key = aliases.get(str(flavor), str(flavor))
    if key not in _CHARGED_LEPTONS:
        raise ValueError(f"unsupported charged-lepton flavor {flavor!r}")
    return key


def _flavor_index(flavor: str) -> int:
    return {"e": 0, "mu": 1, "tau": 2}[_canonical_charged_lepton(flavor)]


def _positive_float(value: Any, name: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0.0:
        raise ValueError(f"{name} must be positive and finite")
    return number


def _nonnegative_float(value: Any, name: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number < 0.0:
        raise ValueError(f"{name} must be nonnegative and finite")
    return number


def _bounded_probability(value: Any, name: str) -> float:
    number = float(value)
    if not math.isfinite(number) or not 0.0 < number < 1.0:
        raise ValueError(f"{name} must lie between zero and one")
    return number


def _finite_complex(value: Any, name: str) -> complex:
    number = complex(value)
    if not math.isfinite(number.real) or not math.isfinite(number.imag):
        raise ValueError(f"{name} must be finite")
    return number


_BOX_LL_KEYS = ("box_ll", "box_left_left", "box_L_L")
_BOX_LR_KEYS = ("box_lr", "box_left_right", "box_L_R")
_BOX_RL_KEYS = ("box_rl", "box_right_left", "box_R_L")
_BOX_RR_KEYS = ("box_rr", "box_right_right", "box_R_R")
_CONTACT_MAPPING_KEYS = frozenset(
    {
        "left_lfv_overlap",
        "right_lfv_overlap",
        "left_emu_overlap",
        "right_emu_overlap",
        "left_emu",
        "right_emu",
        "left_charged_lepton_overlap",
        "right_charged_lepton_overlap",
        "left_lepton_overlap",
        "right_lepton_overlap",
        "left_overlap",
        "right_overlap",
        *_BOX_LL_KEYS,
        *_BOX_LR_KEYS,
        *_BOX_RL_KEYS,
        *_BOX_RR_KEYS,
    }
)
_CONTACT_ATTR_NAMES = tuple(_CONTACT_MAPPING_KEYS)


__all__ = [
    "LFV_THREE_BODY_MODEL_V1",
    "LFV_THREE_BODY_INPUT_BUNDLE_V1",
    "LFV_THREE_BODY_OPERATOR_CONVENTION",
    "LFV_THREE_BODY_DIPOLE_CONTACT_INTERFERENCE_CONVENTION",
    "LFV_THREE_BODY_PROXY_V1",
    "LFVThreeBodySMInputs",
    "LFVThreeBodyContactProxyInput",
    "LFVThreeBodyContactAmplitudes",
    "LFVThreeBodyBranchingResult",
    "default_sm_inputs",
    "lfv_three_body_proxy_input",
    "lfv_three_body_has_contact_proxy",
    "lfv_three_body_contact_amplitudes",
    "lfv_three_body_from_components",
]
