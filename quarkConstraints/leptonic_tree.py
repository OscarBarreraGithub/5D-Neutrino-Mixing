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

from dataclasses import dataclass, field
import math
from typing import Mapping

from .deltaf2 import F_BD

LEPTONIC_TREE_MODEL_V1 = "leptonic_tree_charged_current_proxy_v1"
LEPTONIC_TREE_INPUT_BUNDLE_V1 = "leptonic_tree_bplus_taunu_repo_defaults_v1"
LEPTONIC_TREE_BPLUS_TAUNU_YAML_INPUT_BUNDLE_V1 = (
    "leptonic_tree_bplus_taunu_b009_yaml_anchors_v1"
)
LEPTONIC_TREE_OPERATOR_CONVENTION = (
    "BR(P+ -> l+ nu)=G_F^2 m_P m_l^2 f_P^2 |V_ij|^2 tau_P "
    "(1-m_l^2/m_P^2)^2/(8 pi)"
)
LEPTONIC_TREE_CHARGED_CURRENT_PROXY_ASSUMPTION_V1 = (
    "NEEDS-HUMAN-PHYSICS: no charged-current W/W', charged-Higgs, lepton, "
    "or neutrino matching inputs are available on ParameterPoint; v1 uses a "
    "unit-normalized amplitude proxy A/A_SM=1+(m_P^2/M_KK^2)."
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
