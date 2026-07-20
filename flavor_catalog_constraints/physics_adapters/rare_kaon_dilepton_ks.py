"""K_S -> pi0 e+e- wrapper over the rare-kaon dilepton machinery.

This adapter is append-only relative to the K008/K009 ``K -> pi0 l+l-`` paths.
It maps Phase-3a ``rs_semileptonic_wilsons.s_to_d_ll`` C9/C9p into the
existing K008 y7V slot as a documented shift of the CP-conserving vector
``a_S`` amplitude:

    BR(K_S -> pi0 e+e-) = 5.2e-9 |a_S^eff|^2,
    a_S^eff(s_a) = s_a |a_S| + Re(lambda_y7V^NP) / 1e-4,
    s_a = +/- 1.

The coefficient is the Buchalla-D'Ambrosio-Isidori /
Isidori-Smith-Unterdorfer electron-mode ChPT approximation for the full
kinematic region.  The constraint supplies ``|a_S|`` from the K010 full-rate
anchor, so the SM-limit baseline is measurement-normalized.  The measured
``|a_S|`` leaves the sign free; diagnostics expose both signs, and K010 applies
the conservative sign envelope when forming the HARD verdict.

NEEDS-HUMAN-PHYSICS: a complete K010 prediction needs the long-distance K_S
form factor plus EW KK/Z/Z'/photon-penguin and electron-sector inputs beyond
the Phase-3a light-Z semileptonic Wilsons.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Mapping

from quarkConstraints import rare_kaon_dilepton as _core
from quarkConstraints.rs_semileptonic_wilsons import RSSemileptonicWilsonBundle

from .rare_kaon_dilepton import (
    RARE_KAON_PI0EE_RS_MATCHING_ASSUMPTION_V2,
    RARE_KAON_PI0EE_SEMILEPTONIC_RUNNING_DIAGNOSTIC_V1,
    KLongPi0EEWilsonCoefficients,
    QuarkMassBasisCouplings,
    RareKaonDileptonSMInputs,
    _rs_semileptonic_coeff,
    _tag_rs_result,
    rare_kaon_dilepton_default_sm_inputs,
    rare_kaon_dilepton_g_sm_squared,
    rare_kaon_y7_wilsons_from_rs_semileptonic,
)

RARE_KAON_KS_PI0EE_MODEL_V1 = "rare_kaon_ks_pi0ee_a_s_bdi_y7v_rs_proxy_v1"
RARE_KAON_KS_PI0EE_PARAMETRIZATION_CITATION = (
    "Buchalla, D'Ambrosio, Isidori Nucl. Phys. B672 (2003) 387, "
    "arXiv:hep-ph/0308008; Isidori, Smith, Unterdorfer Eur. Phys. J. C36 "
    "(2004) 57, arXiv:hep-ph/0404127"
)
RARE_KAON_KS_PI0EE_RS_MATCHING_ASSUMPTION_V1 = (
    "NEEDS-HUMAN-PHYSICS: K010 uses the a_S-driven K_S -> pi0 e+e- "
    "full-rate formula, while Phase-3a supplies only the semileptonic "
    "light-Z C9/C9p shift mapped into the vector y7V input. Full EW "
    "KK/Z/Z'/photon-penguin, electron-sector, and K_S form-factor matching "
    "remain outside Phase 3a."
)
KSHORT_PI0EE_A_S_BRANCHING_COEFFICIENT = 5.2e-9
KSHORT_PI0EE_A_S_PROXY_NORMALIZATION = 1.0e-4


@dataclass(frozen=True)
class KShortPi0EEChPTInputs:
    """ChPT input bundle for the K_S electron-mode a_S term."""

    rate_coefficient: float
    a_s_abs: float
    citation: str = RARE_KAON_KS_PI0EE_PARAMETRIZATION_CITATION

    def __post_init__(self) -> None:
        for name in ("rate_coefficient", "a_s_abs"):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be positive and finite")


@dataclass(frozen=True)
class KShortPi0EEResult:
    """Branching-ratio prediction for K_S -> pi0 e+e-.

    ``branching_fraction`` keeps the positive-a_S branch for compatibility.
    K010's HARD verdict uses the explicit positive/negative branch fields. The
    NP shift is a documented proxy, not a complete RS matching.
    """

    model_label: str
    input_bundle: str
    branching_fraction: float
    sm_branching_fraction: float
    np_shift_branching_fraction: float
    branching_fraction_positive_a_s: float
    branching_fraction_negative_a_s: float
    rate_coefficient: float
    a_s_abs: float
    a_s_effective: float
    a_s_effective_positive: float
    a_s_effective_negative: float
    a_s_np_proxy: float
    a_s_np_proxy_complex: complex
    lambda_y7v_np_proxy: complex
    lambda_y7a_np_proxy: complex
    chpt_inputs: KShortPi0EEChPTInputs
    wilsons: KLongPi0EEWilsonCoefficients | None = None
    diagnostics: Mapping[str, float | complex | str | bool] = field(
        default_factory=dict
    )


def kshort_pi0ee_a_s_sm(
    chpt_inputs: KShortPi0EEChPTInputs,
    *,
    inputs: RareKaonDileptonSMInputs | None = None,
) -> KShortPi0EEResult:
    """Evaluate the SM-limit K010 a_S-driven full-region rate."""

    return evaluate_kshort_pi0ee_a_s(None, chpt_inputs=chpt_inputs, inputs=inputs)


def kshort_pi0ee_a_s_from_couplings(
    couplings: QuarkMassBasisCouplings,
    *,
    chpt_inputs: KShortPi0EEChPTInputs,
    m_kk_gev: float | None = None,
    inputs: RareKaonDileptonSMInputs | None = None,
) -> KShortPi0EEResult:
    """Evaluate K010 using the documented K008 y7V proxy shift."""

    return evaluate_kshort_pi0ee_a_s(
        couplings,
        chpt_inputs=chpt_inputs,
        m_kk_gev=m_kk_gev,
        inputs=inputs,
    )


def kshort_pi0ee_a_s_from_rs_semileptonic_wilsons(
    source: RSSemileptonicWilsonBundle,
    *,
    chpt_inputs: KShortPi0EEChPTInputs,
    lepton: str = "e",
    matching_scale_gev: float | None = None,
    inputs: RareKaonDileptonSMInputs | None = None,
) -> KShortPi0EEResult:
    """Evaluate K010 from Phase-3a RS C9/C10 mapped to the y7V a_S shift."""

    p = rare_kaon_dilepton_default_sm_inputs() if inputs is None else inputs
    coeff = _rs_semileptonic_coeff(source, lepton=lepton)
    wilsons = rare_kaon_y7_wilsons_from_rs_semileptonic(
        source,
        lepton=lepton,
        matching_scale_gev=matching_scale_gev,
        inputs=p,
    )
    result = evaluate_kshort_pi0ee_a_s(
        wilsons,
        chpt_inputs=chpt_inputs,
        inputs=p,
    )
    return _tag_rs_result(result, coeff, inputs=p)


def evaluate_kshort_pi0ee_a_s(
    source: QuarkMassBasisCouplings | KLongPi0EEWilsonCoefficients | None = None,
    *,
    chpt_inputs: KShortPi0EEChPTInputs,
    m_kk_gev: float | None = None,
    inputs: RareKaonDileptonSMInputs | None = None,
) -> KShortPi0EEResult:
    """Evaluate the a_S-driven K_S -> pi0 e+e- rate plus RS proxy shift."""

    p = rare_kaon_dilepton_default_sm_inputs() if inputs is None else inputs
    wilsons: KLongPi0EEWilsonCoefficients | None
    if source is None:
        wilsons = None
        lambda_y7v_np_proxy = 0.0j
        lambda_y7a_np_proxy = 0.0j
    elif isinstance(source, KLongPi0EEWilsonCoefficients):
        wilsons = source
        lambda_y7v_np_proxy = source.lambda_y7v_np_proxy
        lambda_y7a_np_proxy = source.lambda_y7a_np_proxy
    else:
        wilsons = _core.compute_klong_pi0ee_y7_wilsons(
            source,
            m_kk_gev=m_kk_gev,
            inputs=p,
        )
        lambda_y7v_np_proxy = wilsons.lambda_y7v_np_proxy
        lambda_y7a_np_proxy = wilsons.lambda_y7a_np_proxy

    a_s_np_proxy_complex = (
        complex(lambda_y7v_np_proxy) / KSHORT_PI0EE_A_S_PROXY_NORMALIZATION
    )
    a_s_np_proxy = float(a_s_np_proxy_complex.real)
    a_s_effective_positive = float(chpt_inputs.a_s_abs + a_s_np_proxy)
    a_s_effective_negative = float(-chpt_inputs.a_s_abs + a_s_np_proxy)
    sm_branching_fraction = float(
        chpt_inputs.rate_coefficient * chpt_inputs.a_s_abs**2
    )
    branching_fraction_positive_a_s = float(
        chpt_inputs.rate_coefficient * a_s_effective_positive**2
    )
    branching_fraction_negative_a_s = float(
        chpt_inputs.rate_coefficient * a_s_effective_negative**2
    )
    branching_fraction = branching_fraction_positive_a_s

    diagnostics: dict[str, float | complex | str | bool] = {
        "g_sm_squared_gev_minus2": rare_kaon_dilepton_g_sm_squared(p),
        "matching_assumption": RARE_KAON_KS_PI0EE_RS_MATCHING_ASSUMPTION_V1,
        "k008_matching_assumption_reused": RARE_KAON_PI0EE_RS_MATCHING_ASSUMPTION_V2,
        "parametrization_citation": RARE_KAON_KS_PI0EE_PARAMETRIZATION_CITATION,
        "ks_a_s_driven_full_rate": True,
        "cp_conserving_ks_mode": True,
        "direct_cp_kl_formula_used": False,
        "short_distance_direct_cp_only": False,
        "uses_k008_y7v_y7a_rs_proxy": True,
        "uses_vector_y7v_proxy_only_for_a_s": True,
        "axial_y7a_not_used_for_ks_a_s_proxy": True,
        "a_s_sign_fixed": False,
        "a_s_sign_ambiguity_evaluated": True,
        "a_s_sign_branches": "both",
        "rate_coefficient_per_a_s_squared": float(chpt_inputs.rate_coefficient),
        "a_s_abs": float(chpt_inputs.a_s_abs),
        "a_s_effective": float(a_s_effective_positive),
        "a_s_effective_positive": float(a_s_effective_positive),
        "a_s_effective_negative": float(a_s_effective_negative),
        "a_s_np_proxy": float(a_s_np_proxy),
        "a_s_np_proxy_complex": complex(a_s_np_proxy_complex),
        "branching_fraction_positive_a_s": float(branching_fraction_positive_a_s),
        "branching_fraction_negative_a_s": float(branching_fraction_negative_a_s),
        "np_shift_branching_fraction_positive_a_s": float(
            branching_fraction_positive_a_s - sm_branching_fraction
        ),
        "np_shift_branching_fraction_negative_a_s": float(
            branching_fraction_negative_a_s - sm_branching_fraction
        ),
        "lambda_y7v_np_proxy": complex(lambda_y7v_np_proxy),
        "lambda_y7a_np_proxy": complex(lambda_y7a_np_proxy),
        "semileptonic_qcd_running_applied": False,
        "semileptonic_qcd_running_multiplicative_factor": 1.0,
        "semileptonic_qcd_running_effect_fraction": 0.0,
        "semileptonic_qcd_running_diagnostic": (
            RARE_KAON_PI0EE_SEMILEPTONIC_RUNNING_DIAGNOSTIC_V1
        ),
    }
    if wilsons is not None:
        diagnostics.update(
            {
                "m_kk_gev": float(wilsons.M_KK),
                "matching_scale_gev": float(wilsons.matching_scale),
                "semileptonic_low_scale_gev": float(wilsons.low_scale_gev),
                "left_sd_coupling": complex(wilsons.left_sd_coupling),
                "right_sd_coupling": complex(wilsons.right_sd_coupling),
                "left_sd_overlap": complex(wilsons.left_sd_overlap),
                "right_sd_overlap": complex(wilsons.right_sd_overlap),
                "left_quark_delta": complex(wilsons.left_quark_delta),
                "right_quark_delta": complex(wilsons.right_quark_delta),
                "quark_vector_delta": complex(wilsons.quark_vector_delta),
                "electron_vector_delta_proxy": float(wilsons.electron_vector_delta),
                "electron_axial_delta_proxy": float(wilsons.electron_axial_delta),
            }
        )

    return KShortPi0EEResult(
        model_label=RARE_KAON_KS_PI0EE_MODEL_V1,
        input_bundle=p.input_bundle,
        branching_fraction=branching_fraction,
        sm_branching_fraction=sm_branching_fraction,
        np_shift_branching_fraction=float(branching_fraction - sm_branching_fraction),
        branching_fraction_positive_a_s=float(branching_fraction_positive_a_s),
        branching_fraction_negative_a_s=float(branching_fraction_negative_a_s),
        rate_coefficient=float(chpt_inputs.rate_coefficient),
        a_s_abs=float(chpt_inputs.a_s_abs),
        a_s_effective=float(a_s_effective_positive),
        a_s_effective_positive=float(a_s_effective_positive),
        a_s_effective_negative=float(a_s_effective_negative),
        a_s_np_proxy=float(a_s_np_proxy),
        a_s_np_proxy_complex=complex(a_s_np_proxy_complex),
        lambda_y7v_np_proxy=complex(lambda_y7v_np_proxy),
        lambda_y7a_np_proxy=complex(lambda_y7a_np_proxy),
        chpt_inputs=chpt_inputs,
        wilsons=wilsons,
        diagnostics=diagnostics,
    )
