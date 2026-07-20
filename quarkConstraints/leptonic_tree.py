"""Tree-level charged-current leptonic meson decays.

The Standard Model rate for a charged pseudoscalar meson is

    BR(P+ -> l+ nu_l) =
        G_F^2 / (8 pi) * m_P * m_l^2
        * (1 - m_l^2 / m_P^2)^2
        * f_P^2 * |V_ij|^2 * tau_P.

The module is intentionally generic: the default input bundle is the
``B+ -> tau+ nu_tau`` mode needed by B009, while callers can instantiate
``LeptonicTreeInputs`` for other charged leptonic modes.

RS matching status
------------------
NEEDS-HUMAN-PHYSICS: the current ``ParameterPoint`` does not carry the
charged-current W/W' tower, charged-Higgs, lepton, or neutrino couplings
needed for a rigorous RS prediction.  The v1 new-physics helper is a
documented charged-current proxy only: it rescales the SM amplitude by a
unit-normalized real shift

    A_total / A_SM = 1 + proxy_strength * (m_P^2 / M_KK^2).

The sign, flavor structure, lepton coupling, scalar helicity factor, and
interference convention are therefore not physics-complete inputs; they are
made explicit in diagnostics for catalog-level scans.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Mapping

from .deltaf2 import F_BD

LEPTONIC_TREE_MODEL_V1 = "leptonic_tree_charged_current_proxy_v1"
LEPTONIC_TREE_INPUT_BUNDLE_V1 = "leptonic_tree_bplus_taunu_repo_defaults_v1"
LEPTONIC_TREE_BPLUS_TAUNU_YAML_INPUT_BUNDLE_V1 = (
    "leptonic_tree_bplus_taunu_b009_yaml_anchors_v1"
)
LEPTONIC_TREE_KPLUS_EMUNU_YAML_INPUT_BUNDLE_V1 = (
    "leptonic_tree_kplus_enu_over_munu_k017_yaml_anchors_v1"
)
LEPTONIC_TREE_OPERATOR_CONVENTION = (
    "BR(P+ -> l+ nu)=G_F^2 m_P m_l^2 f_P^2 |V_ij|^2 tau_P "
    "(1-m_l^2/m_P^2)^2/(8 pi)"
)
LEPTONIC_TREE_LFU_RATIO_OPERATOR_CONVENTION = (
    "R(P+; l1/l2)=Gamma(P+ -> l1+ nu)/Gamma(P+ -> l2+ nu), "
    "computed from the charged-pseudoscalar tree kernel with a caller-supplied "
    "radiative multiplier"
)
LEPTONIC_TREE_CHARGED_CURRENT_PROXY_ASSUMPTION_V1 = (
    "NEEDS-HUMAN-PHYSICS: no charged-current W/W', charged-Higgs, lepton, "
    "or neutrino matching inputs are available on ParameterPoint; v1 uses a "
    "unit-normalized amplitude proxy A/A_SM=1+(m_P^2/M_KK^2)."
)
LEPTONIC_TREE_LFU_RATIO_PROXY_ASSUMPTION_V1 = (
    "NEEDS-HUMAN-PHYSICS: full R_K matching needs lepton-nonuniversal "
    "charged-current W/W', charged-Higgs, heavy-neutrino, lepton-profile, and "
    "radiative-convention inputs that are not available on ParameterPoint; v1 "
    "uses a documented electron-channel-only amplitude proxy "
    "A_e/A_e_SM=1+(m_K^2/M_KK^2), A_mu/A_mu_SM=1."
)


@dataclass(frozen=True)
class LeptonicTreeInputs:
    """Inputs for one charged pseudoscalar leptonic decay."""

    input_bundle: str = LEPTONIC_TREE_INPUT_BUNDLE_V1
    mode_key: str = "bplus_taunu"
    display_name: str = "B+ -> tau+ nu_tau"
    meson_mass_gev: float = 5.27934
    lepton_mass_gev: float = 1.77686
    decay_constant_gev: float = F_BD
    ckm_abs: float = 0.00368
    lifetime_ps: float = 1.638
    gf_gev_minus2: float = 1.1663787e-5
    hbar_gev_s: float = 6.582119569e-25
    proxy_strength: float = 1.0
    constants_citation: str = (
        "PDG-era B+/tau lifetime and masses; f_B from quarkConstraints.deltaf2 "
        "FLAG-style F_BD=0.190 GeV; repo CKM target |V_ub|=0.00368. "
        "Live B009 constraints override f_B and |V_ub| from YAML anchors."
    )

    def __post_init__(self) -> None:
        for name in (
            "meson_mass_gev",
            "lepton_mass_gev",
            "decay_constant_gev",
            "ckm_abs",
            "lifetime_ps",
            "gf_gev_minus2",
            "hbar_gev_s",
        ):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be positive and finite")
        if float(self.lepton_mass_gev) >= float(self.meson_mass_gev):
            raise ValueError("lepton_mass_gev must be smaller than meson_mass_gev")
        if not math.isfinite(float(self.proxy_strength)):
            raise ValueError("proxy_strength must be finite")


@dataclass(frozen=True)
class LeptonicTreeBranchingResult:
    """Branching-fraction prediction for a charged leptonic meson mode."""

    model_label: str
    input_bundle: str
    mode_key: str
    branching_fraction: float
    sm_branching_fraction: float
    np_shift_branching_fraction: float
    charged_current_np_amplitude_ratio: complex
    amplitude_multiplier: float
    M_KK: float | None = None
    matching_scale: float | None = None
    diagnostics: Mapping[str, float | complex | str | bool | None] = field(
        default_factory=dict
    )


@dataclass(frozen=True)
class LeptonicTreeRatioInputs:
    """Inputs for one charged-pseudoscalar leptonic LFU ratio."""

    input_bundle: str = LEPTONIC_TREE_KPLUS_EMUNU_YAML_INPUT_BUNDLE_V1
    mode_key: str = "kplus_enu_over_munu"
    display_name: str = "K+ -> e+ nu over K+ -> mu+ nu"
    numerator_mode_key: str = "kplus_enu"
    denominator_mode_key: str = "kplus_munu"
    meson_mass_gev: float = 0.493677
    numerator_lepton_mass_gev: float = 0.00051099895
    denominator_lepton_mass_gev: float = 0.1056583755
    radiative_correction_multiplier: float = 1.0
    proxy_strength: float = 1.0
    constants_citation: str = (
        "PDG-era charged-kaon, electron, and muon masses. The radiative "
        "multiplier is supplied by the caller; K017 derives it from its "
        "Cirigliano-Rosell SM anchor."
    )

    def __post_init__(self) -> None:
        for name in (
            "meson_mass_gev",
            "numerator_lepton_mass_gev",
            "denominator_lepton_mass_gev",
            "radiative_correction_multiplier",
        ):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be positive and finite")
        if float(self.numerator_lepton_mass_gev) >= float(self.meson_mass_gev):
            raise ValueError("numerator_lepton_mass_gev must be below meson_mass_gev")
        if float(self.denominator_lepton_mass_gev) >= float(self.meson_mass_gev):
            raise ValueError("denominator_lepton_mass_gev must be below meson_mass_gev")
        if not math.isfinite(float(self.proxy_strength)):
            raise ValueError("proxy_strength must be finite")


@dataclass(frozen=True)
class LeptonicTreeRatioResult:
    """Prediction for a charged-leptonic LFU ratio."""

    model_label: str
    input_bundle: str
    mode_key: str
    ratio: float
    sm_ratio: float
    np_shift_ratio: float
    numerator_np_amplitude_ratio: complex
    denominator_np_amplitude_ratio: complex
    numerator_amplitude_multiplier: float
    denominator_amplitude_multiplier: float
    M_KK: float | None = None
    matching_scale: float | None = None
    diagnostics: Mapping[str, float | complex | str | bool | None] = field(
        default_factory=dict
    )


def default_bplus_tau_nu_inputs() -> LeptonicTreeInputs:
    """Return the default ``B+ -> tau+ nu_tau`` input bundle."""
    return LeptonicTreeInputs()


def bplus_tau_nu_inputs_from_theory_anchors(
    *,
    decay_constant_gev: float,
    ckm_abs: float,
    constants_citation: str,
    input_bundle: str = LEPTONIC_TREE_BPLUS_TAUNU_YAML_INPUT_BUNDLE_V1,
) -> LeptonicTreeInputs:
    """Return ``B+ -> tau+ nu_tau`` inputs with YAML-backed theory values."""
    return LeptonicTreeInputs(
        input_bundle=input_bundle,
        decay_constant_gev=float(decay_constant_gev),
        ckm_abs=float(ckm_abs),
        constants_citation=str(constants_citation),
    )


def default_kplus_enu_over_munu_inputs() -> LeptonicTreeRatioInputs:
    """Return default ``R_K = Gamma(K -> e nu)/Gamma(K -> mu nu)`` inputs."""
    return LeptonicTreeRatioInputs()


def _leptonic_lfu_width_kernel(*, meson_mass_gev: float, lepton_mass_gev: float) -> float:
    meson_mass = float(meson_mass_gev)
    lepton_mass = float(lepton_mass_gev)
    if not math.isfinite(meson_mass) or meson_mass <= 0.0:
        raise ValueError("meson_mass_gev must be positive and finite")
    if not math.isfinite(lepton_mass) or lepton_mass <= 0.0:
        raise ValueError("lepton_mass_gev must be positive and finite")
    if lepton_mass >= meson_mass:
        raise ValueError("lepton_mass_gev must be below meson_mass_gev")
    phase = (1.0 - (lepton_mass / meson_mass) ** 2) ** 2
    return float(lepton_mass * lepton_mass * phase)


def leptonic_lfu_tree_ratio_without_radiation(
    inputs: LeptonicTreeRatioInputs | None = None,
) -> float:
    """Return the helicity-suppressed tree ratio before radiative corrections."""
    p = default_kplus_enu_over_munu_inputs() if inputs is None else inputs
    numerator = _leptonic_lfu_width_kernel(
        meson_mass_gev=p.meson_mass_gev,
        lepton_mass_gev=p.numerator_lepton_mass_gev,
    )
    denominator = _leptonic_lfu_width_kernel(
        meson_mass_gev=p.meson_mass_gev,
        lepton_mass_gev=p.denominator_lepton_mass_gev,
    )
    return float(numerator / denominator)


def kplus_enu_over_munu_inputs_from_sm_ratio_anchor(
    *,
    sm_ratio: float,
    constants_citation: str,
    input_bundle: str = LEPTONIC_TREE_KPLUS_EMUNU_YAML_INPUT_BUNDLE_V1,
) -> LeptonicTreeRatioInputs:
    """Return K LFU inputs with the radiative multiplier fixed by an SM anchor."""
    sm_value = float(sm_ratio)
    if not math.isfinite(sm_value) or sm_value <= 0.0:
        raise ValueError("sm_ratio must be positive and finite")
    base = LeptonicTreeRatioInputs(
        input_bundle=input_bundle,
        radiative_correction_multiplier=1.0,
        constants_citation=str(constants_citation),
    )
    tree_ratio = leptonic_lfu_tree_ratio_without_radiation(base)
    return LeptonicTreeRatioInputs(
        input_bundle=input_bundle,
        radiative_correction_multiplier=float(sm_value / tree_ratio),
        constants_citation=str(constants_citation),
    )


def lifetime_ps_to_gev_inverse(tau_ps: float, hbar_gev_s: float) -> float:
    """Convert a lifetime in ps to natural units."""
    tau = float(tau_ps)
    hbar = float(hbar_gev_s)
    if not math.isfinite(tau) or tau <= 0.0:
        raise ValueError("tau_ps must be positive and finite")
    if not math.isfinite(hbar) or hbar <= 0.0:
        raise ValueError("hbar_gev_s must be positive and finite")
    return float(tau * 1.0e-12 / hbar)


def sm_branching_fraction(inputs: LeptonicTreeInputs | None = None) -> LeptonicTreeBranchingResult:
    """Evaluate the SM tree-level branching fraction."""
    p = default_bplus_tau_nu_inputs() if inputs is None else inputs
    tau_gev_inv = lifetime_ps_to_gev_inverse(p.lifetime_ps, p.hbar_gev_s)
    phase = (1.0 - (p.lepton_mass_gev / p.meson_mass_gev) ** 2) ** 2
    br = (
        p.gf_gev_minus2**2
        / (8.0 * math.pi)
        * p.meson_mass_gev
        * p.lepton_mass_gev**2
        * phase
        * p.decay_constant_gev**2
        * p.ckm_abs**2
        * tau_gev_inv
    )
    diagnostics: dict[str, float | complex | str | bool | None] = {
        "operator_convention": LEPTONIC_TREE_OPERATOR_CONVENTION,
        "matching_assumption": "SM tree level; no new-physics amplitude applied",
        "gf_gev_minus2": float(p.gf_gev_minus2),
        "meson_mass_gev": float(p.meson_mass_gev),
        "lepton_mass_gev": float(p.lepton_mass_gev),
        "decay_constant_gev": float(p.decay_constant_gev),
        "ckm_abs": float(p.ckm_abs),
        "lifetime_ps": float(p.lifetime_ps),
        "lifetime_gev_inverse": float(tau_gev_inv),
        "helicity_phase_space": float(phase),
        "constants_citation": p.constants_citation,
    }
    return LeptonicTreeBranchingResult(
        model_label=LEPTONIC_TREE_MODEL_V1,
        input_bundle=p.input_bundle,
        mode_key=p.mode_key,
        branching_fraction=float(br),
        sm_branching_fraction=float(br),
        np_shift_branching_fraction=0.0,
        charged_current_np_amplitude_ratio=0.0j,
        amplitude_multiplier=1.0,
        diagnostics=diagnostics,
    )


def sm_leptonic_lfu_ratio(
    inputs: LeptonicTreeRatioInputs | None = None,
) -> LeptonicTreeRatioResult:
    """Evaluate the SM charged-leptonic LFU ratio."""
    p = default_kplus_enu_over_munu_inputs() if inputs is None else inputs
    tree_ratio = leptonic_lfu_tree_ratio_without_radiation(p)
    sm_ratio = float(tree_ratio * p.radiative_correction_multiplier)
    diagnostics: dict[str, float | complex | str | bool | None] = {
        "operator_convention": LEPTONIC_TREE_LFU_RATIO_OPERATOR_CONVENTION,
        "matching_assumption": "SM tree level with caller-supplied radiative multiplier",
        "meson_mass_gev": float(p.meson_mass_gev),
        "numerator_lepton_mass_gev": float(p.numerator_lepton_mass_gev),
        "denominator_lepton_mass_gev": float(p.denominator_lepton_mass_gev),
        "tree_ratio_without_radiation": float(tree_ratio),
        "radiative_correction_multiplier": float(p.radiative_correction_multiplier),
        "constants_citation": p.constants_citation,
        "numerator_mode_key": p.numerator_mode_key,
        "denominator_mode_key": p.denominator_mode_key,
    }
    return LeptonicTreeRatioResult(
        model_label=LEPTONIC_TREE_MODEL_V1,
        input_bundle=p.input_bundle,
        mode_key=p.mode_key,
        ratio=sm_ratio,
        sm_ratio=sm_ratio,
        np_shift_ratio=0.0,
        numerator_np_amplitude_ratio=0.0j,
        denominator_np_amplitude_ratio=0.0j,
        numerator_amplitude_multiplier=1.0,
        denominator_amplitude_multiplier=1.0,
        diagnostics=diagnostics,
    )


def charged_current_proxy_amplitude_ratio(
    *,
    meson_mass_gev: float,
    m_kk_gev: float,
    proxy_strength: float = 1.0,
) -> complex:
    """Return the v1 charged-current NP amplitude ratio."""
    meson_mass = float(meson_mass_gev)
    m_kk = float(m_kk_gev)
    strength = float(proxy_strength)
    if not math.isfinite(meson_mass) or meson_mass <= 0.0:
        raise ValueError("meson_mass_gev must be positive and finite")
    if not math.isfinite(m_kk) or m_kk <= 0.0:
        raise ValueError("m_kk_gev must be positive and finite")
    if not math.isfinite(strength):
        raise ValueError("proxy_strength must be finite")
    return complex(strength * (meson_mass / m_kk) ** 2)


def _checked_complex_ratio(value: complex, *, field_name: str) -> complex:
    ratio = complex(value)
    if not math.isfinite(ratio.real) or not math.isfinite(ratio.imag):
        raise ValueError(f"{field_name} must be finite")
    return ratio


def evaluate_leptonic_branching_fraction(
    *,
    m_kk_gev: float | None = None,
    np_amplitude_ratio: complex | None = None,
    inputs: LeptonicTreeInputs | None = None,
) -> LeptonicTreeBranchingResult:
    """Evaluate a charged leptonic branching fraction.

    If ``np_amplitude_ratio`` is supplied it is used directly.  Otherwise,
    a provided ``m_kk_gev`` activates the documented charged-current proxy.
    With neither input, the result is the SM limit.
    """
    p = default_bplus_tau_nu_inputs() if inputs is None else inputs
    sm = sm_branching_fraction(p)
    if np_amplitude_ratio is None:
        if m_kk_gev is None:
            amplitude_ratio = 0.0j
        else:
            amplitude_ratio = charged_current_proxy_amplitude_ratio(
                meson_mass_gev=p.meson_mass_gev,
                m_kk_gev=m_kk_gev,
                proxy_strength=p.proxy_strength,
            )
    else:
        amplitude_ratio = complex(np_amplitude_ratio)
        if not math.isfinite(amplitude_ratio.real) or not math.isfinite(
            amplitude_ratio.imag
        ):
            raise ValueError("np_amplitude_ratio must be finite")

    multiplier = abs(1.0 + amplitude_ratio) ** 2
    total = float(sm.sm_branching_fraction * multiplier)
    m_kk = None if m_kk_gev is None else float(m_kk_gev)
    diagnostics = dict(sm.diagnostics)
    diagnostics.update(
        {
            "matching_assumption": (
                LEPTONIC_TREE_CHARGED_CURRENT_PROXY_ASSUMPTION_V1
                if (m_kk_gev is not None or np_amplitude_ratio is not None)
                else "SM tree level; no new-physics amplitude applied"
            ),
            "m_kk_gev": m_kk,
            "matching_scale_gev": m_kk,
            "proxy_strength": float(p.proxy_strength),
            "charged_current_np_amplitude_ratio": complex(amplitude_ratio),
            "amplitude_multiplier": float(multiplier),
            "np_shift_branching_fraction": float(total - sm.sm_branching_fraction),
            "uses_charged_current_proxy": bool(
                np_amplitude_ratio is None and m_kk_gev is not None
            ),
            "needs_human_physics": LEPTONIC_TREE_CHARGED_CURRENT_PROXY_ASSUMPTION_V1,
        }
    )
    return LeptonicTreeBranchingResult(
        model_label=LEPTONIC_TREE_MODEL_V1,
        input_bundle=p.input_bundle,
        mode_key=p.mode_key,
        branching_fraction=total,
        sm_branching_fraction=float(sm.sm_branching_fraction),
        np_shift_branching_fraction=float(total - sm.sm_branching_fraction),
        charged_current_np_amplitude_ratio=complex(amplitude_ratio),
        amplitude_multiplier=float(multiplier),
        M_KK=m_kk,
        matching_scale=m_kk,
        diagnostics=diagnostics,
    )


def evaluate_leptonic_lfu_ratio(
    *,
    m_kk_gev: float | None = None,
    numerator_np_amplitude_ratio: complex | None = None,
    denominator_np_amplitude_ratio: complex | None = None,
    inputs: LeptonicTreeRatioInputs | None = None,
) -> LeptonicTreeRatioResult:
    """Evaluate a charged-leptonic LFU ratio in the SM or v1 LNU proxy model."""
    p = default_kplus_enu_over_munu_inputs() if inputs is None else inputs
    sm = sm_leptonic_lfu_ratio(p)

    if numerator_np_amplitude_ratio is None and denominator_np_amplitude_ratio is None:
        if m_kk_gev is None:
            numerator_ratio = 0.0j
            denominator_ratio = 0.0j
        else:
            numerator_ratio = charged_current_proxy_amplitude_ratio(
                meson_mass_gev=p.meson_mass_gev,
                m_kk_gev=m_kk_gev,
                proxy_strength=p.proxy_strength,
            )
            denominator_ratio = 0.0j
    else:
        numerator_ratio = (
            0.0j
            if numerator_np_amplitude_ratio is None
            else _checked_complex_ratio(
                numerator_np_amplitude_ratio,
                field_name="numerator_np_amplitude_ratio",
            )
        )
        denominator_ratio = (
            0.0j
            if denominator_np_amplitude_ratio is None
            else _checked_complex_ratio(
                denominator_np_amplitude_ratio,
                field_name="denominator_np_amplitude_ratio",
            )
        )

    numerator_multiplier = abs(1.0 + numerator_ratio) ** 2
    denominator_multiplier = abs(1.0 + denominator_ratio) ** 2
    if denominator_multiplier <= 0.0:
        raise ValueError("denominator amplitude multiplier must be positive")

    total = float(sm.sm_ratio * numerator_multiplier / denominator_multiplier)
    m_kk = None if m_kk_gev is None else float(m_kk_gev)
    diagnostics = dict(sm.diagnostics)
    diagnostics.update(
        {
            "matching_assumption": (
                LEPTONIC_TREE_LFU_RATIO_PROXY_ASSUMPTION_V1
                if (
                    m_kk_gev is not None
                    or numerator_np_amplitude_ratio is not None
                    or denominator_np_amplitude_ratio is not None
                )
                else "SM tree level; no lepton-nonuniversal amplitude applied"
            ),
            "m_kk_gev": m_kk,
            "matching_scale_gev": m_kk,
            "proxy_strength": float(p.proxy_strength),
            "numerator_np_amplitude_ratio": complex(numerator_ratio),
            "denominator_np_amplitude_ratio": complex(denominator_ratio),
            "numerator_amplitude_multiplier": float(numerator_multiplier),
            "denominator_amplitude_multiplier": float(denominator_multiplier),
            "np_shift_ratio": float(total - sm.sm_ratio),
            "uses_lepton_nonuniversal_proxy": bool(
                numerator_np_amplitude_ratio is None
                and denominator_np_amplitude_ratio is None
                and m_kk_gev is not None
            ),
            "needs_human_physics": LEPTONIC_TREE_LFU_RATIO_PROXY_ASSUMPTION_V1,
        }
    )
    return LeptonicTreeRatioResult(
        model_label=LEPTONIC_TREE_MODEL_V1,
        input_bundle=p.input_bundle,
        mode_key=p.mode_key,
        ratio=total,
        sm_ratio=float(sm.sm_ratio),
        np_shift_ratio=float(total - sm.sm_ratio),
        numerator_np_amplitude_ratio=complex(numerator_ratio),
        denominator_np_amplitude_ratio=complex(denominator_ratio),
        numerator_amplitude_multiplier=float(numerator_multiplier),
        denominator_amplitude_multiplier=float(denominator_multiplier),
        M_KK=m_kk,
        matching_scale=m_kk,
        diagnostics=diagnostics,
    )
