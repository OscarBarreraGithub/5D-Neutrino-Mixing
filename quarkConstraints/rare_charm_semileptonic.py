"""Short-distance ``D+ -> pi+ mu+ mu-`` rare-charm machinery.

This module deliberately reuses :mod:`quarkConstraints.rare_charm_dilepton`
for the shared ``c -> u l+ l-`` Wilson proxy introduced for C004/C005.  The
new code here is only the exclusive ``D -> pi`` short-distance rate:

    dBR/dq2 = 2 tau_D [a_mu(q2) + c_mu(q2) / 3],

with the vector/axial pieces of the standard ``D -> pi l l`` distribution and
the Becirevic-Kaidalov-style ``D -> pi`` form-factor normalization used in
rare-charm phenomenology.  Resonance amplitudes are not modeled here.

NEEDS-HUMAN-PHYSICS: the RS contribution is still the same documented
Z/KK-penguin proxy from the rare-charm dilepton core.  The current
``ParameterPoint`` does not carry the full electroweak KK/Z/Z', charged-lepton,
scalar/tensor, resonance, or experimental dimuon-window information required
for a complete ``D+ -> pi+ mu+ mu-`` RS recast.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Mapping

import numpy as np

from .couplings import QuarkMassBasisCouplings
from .rare_charm_dilepton import (
    RARE_CHARM_DILEPTON_RS_MATCHING_ASSUMPTION_V1,
    RareCharmDileptonSMInputs,
    RareCharmDileptonWilsonCoefficients,
    ckm_factors,
    compute_rare_charm_dilepton_wilsons,
)
from .rare_charm_dilepton import (
    default_sm_inputs as rare_charm_dilepton_default_sm_inputs,
)

RARE_CHARM_DTOPI_MUMU_MODEL_V1 = "rare_charm_dtopi_mumu_sd_form_factor_proxy_v1"
RARE_CHARM_DTOPI_MUMU_INPUT_BUNDLE_V1 = (
    "rare_charm_dtopi_mumu_sd_inputs_form_factor_v1"
)
RARE_CHARM_DTOPI_MUMU_FORM_FACTOR_MODEL_V1 = (
    "D_to_pi_fplus_fzero_BK_shape_fplus0_0p67_v1"
)
RARE_CHARM_DTOPI_MUMU_OPERATOR_CONVENTION = (
    "H_eff=-4 G_F/sqrt(2) lambda_b alpha/(4 pi) "
    "[C9 O9 + C10 O10 + C9' O9' + C10' O10']; "
    "D->pi semileptonic vector/axial rate uses C9+C9' and C10+C10'."
)
RARE_CHARM_DTOPI_MUMU_PARAMETRIZATION_CITATION = (
    "D+ -> pi+ l+l- vector/axial distribution and D->pi form factors as in "
    "Sahoo et al. EPJC 77 (2017) 344, Eqs. 43-50; rare-charm "
    "long-distance caveat from de Boer-Hiller Phys. Rev. D93 (2016) 074001"
)
RARE_CHARM_DTOPI_MUMU_RESONANCE_LIMITATION_V1 = (
    "Long-distance phi/rho/omega resonance amplitudes and the experiment's "
    "exact nonresonant dimuon-window acceptance are not modeled; v1 evaluates "
    "the smooth short-distance full-q2 proxy and the constraint layer compares "
    "that proxy to the catalogued C007 upper-limit budget."
)
RARE_CHARM_DTOPI_MUMU_Q2_TREATMENT_V1 = (
    "full_kinematic_q2_smooth_short_distance_proxy_no_lhcb_nonresonant_window_"
    "or_resonance_acceptance"
)


@dataclass(frozen=True)
class RareCharmDToPiFormFactorInputs:
    """Exclusive ``D -> pi`` form-factor normalization for the SD proxy."""

    model_label: str = RARE_CHARM_DTOPI_MUMU_FORM_FACTOR_MODEL_V1
    fplus_0: float = 0.67
    fplus_shape_a: float = 0.28
    fzero_shape_b: float = 1.27
    pole_mass_gev: float = 1.90
    source: str = (
        "Sahoo et al. EPJC 77 (2017) 344 quoting D -> pi l nu inputs: "
        "f_+(0)=0.67(3), m_pole=1.90(8) GeV, a=0.28(14), b=1.27(17)"
    )

    def __post_init__(self) -> None:
        for name in ("fplus_0", "fzero_shape_b", "pole_mass_gev"):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be positive and finite")
        if not math.isfinite(float(self.fplus_shape_a)):
            raise ValueError("fplus_shape_a must be finite")


@dataclass(frozen=True)
class RareCharmDToPiMuMuInputs:
    """Numerical inputs for the ``D+ -> pi+ mu+ mu-`` SD proxy."""

    input_bundle: str = RARE_CHARM_DTOPI_MUMU_INPUT_BUNDLE_V1
    rare_charm: RareCharmDileptonSMInputs = field(
        default_factory=rare_charm_dilepton_default_sm_inputs
    )
    form_factor: RareCharmDToPiFormFactorInputs = field(
        default_factory=RareCharmDToPiFormFactorInputs
    )
    dplus_mass_gev: float = 1.86965
    piplus_mass_gev: float = 0.13957039
    dplus_lifetime_ps: float = 1.033
    quadrature_points: int = 160
    constants_citation: str = (
        "PDG-era D+/pi+ masses and tau_D+; muon, G_F, alpha_e, hbar and CKM "
        "from the shared rare_charm_dilepton input bundle"
    )

    def __post_init__(self) -> None:
        for name in ("dplus_mass_gev", "piplus_mass_gev", "dplus_lifetime_ps"):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be positive and finite")
        if self.dplus_mass_gev <= self.piplus_mass_gev + 2.0 * self.rare_charm.muon.mass_gev:
            raise ValueError("D+ -> pi+ mu+ mu- phase space is closed")
        if int(self.quadrature_points) < 16:
            raise ValueError("quadrature_points must be at least 16")


@dataclass(frozen=True)
class RareCharmDToPiMuMuBranchingResult:
    """Short-distance branching-ratio prediction for ``D+ -> pi+ mu+ mu-``."""

    model_label: str
    input_bundle: str
    transition_key: str
    lepton_key: str
    branching_fraction: float
    sm_branching_fraction: float
    np_shift_branching_fraction: float
    q2_min_gev2: float
    q2_max_gev2: float
    c9_semileptonic_np: complex
    c10_semileptonic_np: complex
    c9_np: complex
    c10_np: complex
    c9p_np: complex
    c10p_np: complex
    lambda_b: complex
    wilsons: RareCharmDileptonWilsonCoefficients | None = None
    diagnostics: Mapping[str, float | complex | str | bool] = field(
        default_factory=dict
    )


def default_dtopi_mumu_inputs() -> RareCharmDToPiMuMuInputs:
    """Return the default ``D+ -> pi+ mu+ mu-`` SD input bundle."""
    return RareCharmDToPiMuMuInputs()


def _positive_float(value: object, name: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0.0:
        raise ValueError(f"{name} must be positive and finite")
    return number


def _tau_ps_to_gev_inverse(tau_ps: float, hbar_gev_s: float) -> float:
    return float(tau_ps * 1.0e-12 / hbar_gev_s)


def _kallen(a: float, b: float, c: float) -> float:
    return float(a * a + b * b + c * c - 2.0 * (a * b + a * c + b * c))


def dtopi_fplus(q2_gev2: float, inputs: RareCharmDToPiMuMuInputs | None = None) -> float:
    """Return ``f_+(q2)`` for the repo-owned ``D -> pi`` proxy shape."""
    p = default_dtopi_mumu_inputs() if inputs is None else inputs
    ff = p.form_factor
    x = float(q2_gev2) / ff.pole_mass_gev**2
    denominator = (1.0 - x) * (1.0 - ff.fplus_shape_a * x)
    if denominator <= 0.0 or not math.isfinite(denominator):
        raise ValueError(f"invalid D->pi f_+ denominator at q2={q2_gev2}")
    return float(ff.fplus_0 / denominator)


def dtopi_fzero(q2_gev2: float, inputs: RareCharmDToPiMuMuInputs | None = None) -> float:
    """Return ``f_0(q2)`` for the repo-owned ``D -> pi`` proxy shape."""
    p = default_dtopi_mumu_inputs() if inputs is None else inputs
    ff = p.form_factor
    x = float(q2_gev2) / ff.pole_mass_gev**2
    denominator = 1.0 - x / ff.fzero_shape_b
    if denominator <= 0.0 or not math.isfinite(denominator):
        raise ValueError(f"invalid D->pi f_0 denominator at q2={q2_gev2}")
    return float(ff.fplus_0 / denominator)


def dtopi_mumu_q2_range(
    inputs: RareCharmDToPiMuMuInputs | None = None,
) -> tuple[float, float]:
    """Return the kinematic dimuon ``q2`` range for ``D+ -> pi+ mu+ mu-``."""
    p = default_dtopi_mumu_inputs() if inputs is None else inputs
    m_mu = p.rare_charm.muon.mass_gev
    return (
        float(4.0 * m_mu * m_mu),
        float((p.dplus_mass_gev - p.piplus_mass_gev) ** 2),
    )


def dtopi_mumu_differential_branching_fraction(
    q2_gev2: float,
    *,
    c9_semileptonic: complex,
    c10_semileptonic: complex,
    inputs: RareCharmDToPiMuMuInputs | None = None,
) -> float:
    """Evaluate ``dBR(D+ -> pi+ mu+ mu-)_SD / dq2`` at one ``q2`` point."""
    p = default_dtopi_mumu_inputs() if inputs is None else inputs
    q2 = float(q2_gev2)
    q2_min, q2_max = dtopi_mumu_q2_range(p)
    if q2 < q2_min or q2 > q2_max:
        raise ValueError(f"q2={q2} outside D+ -> pi+ mu+ mu- physical range")

    m_d = p.dplus_mass_gev
    m_pi = p.piplus_mass_gev
    m_mu = p.rare_charm.muon.mass_gev
    m_d2 = m_d * m_d
    m_pi2 = m_pi * m_pi
    lam = max(0.0, _kallen(m_d2, m_pi2, q2))
    beta2 = max(0.0, 1.0 - 4.0 * m_mu * m_mu / q2)
    sqrt_lam = math.sqrt(lam)
    beta = math.sqrt(beta2)

    fplus = dtopi_fplus(q2, p)
    fzero = dtopi_fzero(q2, p)
    vector_amplitude = complex(fplus * c9_semileptonic)
    c10_total = complex(c10_semileptonic)
    axial_amplitude = complex(fplus * c10_total)
    pseudoscalar_amplitude = complex(
        -m_mu
        * (
            fplus
            - (m_d2 - m_pi2) / q2 * (fzero - fplus)
        )
        * c10_total
    )

    factors = ckm_factors("c_u", p.rare_charm)
    gamma0 = (
        p.rare_charm.gf_gev_minus2**2
        * p.rare_charm.alpha_em_mz**2
        * abs(factors.lambda_b) ** 2
        / ((4.0 * math.pi) ** 5 * m_d**3)
    )
    a_term = gamma0 * sqrt_lam * beta * (
        2.0 * q2 * abs(pseudoscalar_amplitude) ** 2
        + 0.5 * lam * (
            abs(axial_amplitude) ** 2 + abs(vector_amplitude) ** 2
        )
        + 4.0
        * m_mu
        * (m_d2 - m_pi2 + q2)
        * (axial_amplitude * pseudoscalar_amplitude.conjugate()).real
        + 8.0 * m_mu * m_mu * m_d2 * abs(axial_amplitude) ** 2
    )
    c_term = gamma0 * sqrt_lam * beta * (
        -0.5
        * lam
        * beta2
        * (abs(vector_amplitude) ** 2 + abs(axial_amplitude) ** 2)
    )
    tau = _tau_ps_to_gev_inverse(
        p.dplus_lifetime_ps,
        p.rare_charm.hbar_gev_s,
    )
    rate = 2.0 * tau * (a_term + c_term / 3.0)
    if rate < 0.0 and abs(rate) < 1.0e-30:
        return 0.0
    if rate < 0.0:
        raise ValueError(f"negative dBR/dq2={rate} at q2={q2}")
    return float(rate)


def _integrate_branching_fraction(
    *,
    c9_semileptonic: complex,
    c10_semileptonic: complex,
    inputs: RareCharmDToPiMuMuInputs,
) -> float:
    q2_min, q2_max = dtopi_mumu_q2_range(inputs)
    nodes, weights = np.polynomial.legendre.leggauss(int(inputs.quadrature_points))
    mid = 0.5 * (q2_max + q2_min)
    half_width = 0.5 * (q2_max - q2_min)
    total = 0.0
    for node, weight in zip(nodes, weights):
        q2 = mid + half_width * float(node)
        total += float(weight) * dtopi_mumu_differential_branching_fraction(
            q2,
            c9_semileptonic=c9_semileptonic,
            c10_semileptonic=c10_semileptonic,
            inputs=inputs,
        )
    return float(half_width * total)


def dtopi_mumu_sm(
    inputs: RareCharmDToPiMuMuInputs | None = None,
) -> RareCharmDToPiMuMuBranchingResult:
    """Evaluate the zero-Wilson SM-limit of the smooth SD proxy."""
    return evaluate_dplus_to_piplus_mumu(None, inputs=inputs)


def evaluate_dplus_to_piplus_mumu(
    source: QuarkMassBasisCouplings | RareCharmDileptonWilsonCoefficients | None = None,
    *,
    m_kk_gev: float | None = None,
    inputs: RareCharmDToPiMuMuInputs | None = None,
) -> RareCharmDToPiMuMuBranchingResult:
    """Evaluate smooth ``BR_SD(D+ -> pi+ mu+ mu-)`` with the v1 RS proxy."""
    p = default_dtopi_mumu_inputs() if inputs is None else inputs
    wilsons: RareCharmDileptonWilsonCoefficients | None
    if source is None:
        wilsons = None
        c9_np = c10_np = c9p_np = c10p_np = 0.0j
    elif isinstance(source, RareCharmDileptonWilsonCoefficients):
        wilsons = source
        c9_np = source.c9_np
        c10_np = source.c10_np
        c9p_np = source.c9p_np
        c10p_np = source.c10p_np
    else:
        wilsons = compute_rare_charm_dilepton_wilsons(
            source,
            transition="c_u",
            m_kk_gev=m_kk_gev,
            inputs=p.rare_charm,
        )
        c9_np = wilsons.c9_np
        c10_np = wilsons.c10_np
        c9p_np = wilsons.c9p_np
        c10p_np = wilsons.c10p_np

    c9_semileptonic = complex(c9_np + c9p_np)
    c10_semileptonic = complex(p.rare_charm.c10_sm + c10_np + c10p_np)
    branching = _integrate_branching_fraction(
        c9_semileptonic=c9_semileptonic,
        c10_semileptonic=c10_semileptonic,
        inputs=p,
    )
    sm_branching = _integrate_branching_fraction(
        c9_semileptonic=0.0j,
        c10_semileptonic=complex(p.rare_charm.c10_sm),
        inputs=p,
    )
    factors = ckm_factors("c_u", p.rare_charm)
    q2_min, q2_max = dtopi_mumu_q2_range(p)
    diagnostics: dict[str, float | complex | str | bool] = {
        "lambda_wolfenstein": factors.lambda_wolfenstein,
        "lambda_b": factors.lambda_b,
        "c10_sm": float(p.rare_charm.c10_sm),
        "c9_semileptonic_np": complex(c9_semileptonic),
        "c10_semileptonic_np": complex(c10_semileptonic - p.rare_charm.c10_sm),
        "semileptonic_primed_combination_is_plus": True,
        "leptonic_c10_minus_c10p_not_used_for_dtopi": True,
        "matching_assumption": RARE_CHARM_DILEPTON_RS_MATCHING_ASSUMPTION_V1,
        "operator_convention": RARE_CHARM_DTOPI_MUMU_OPERATOR_CONVENTION,
        "parametrization_citation": RARE_CHARM_DTOPI_MUMU_PARAMETRIZATION_CITATION,
        "resonance_limitation": RARE_CHARM_DTOPI_MUMU_RESONANCE_LIMITATION_V1,
        "q2_treatment": RARE_CHARM_DTOPI_MUMU_Q2_TREATMENT_V1,
        "lhcb_nonresonant_dimuon_window_applied": False,
        "short_distance_only": True,
        "resonance_amplitudes_included": False,
        "q2_min_gev2": float(q2_min),
        "q2_max_gev2": float(q2_max),
        "quadrature_points": float(p.quadrature_points),
        "dplus_mass_gev": float(p.dplus_mass_gev),
        "piplus_mass_gev": float(p.piplus_mass_gev),
        "dplus_lifetime_ps": float(p.dplus_lifetime_ps),
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
                "lepton_left_delta": float(wilsons.lepton_left_delta),
                "lepton_right_delta": float(wilsons.lepton_right_delta),
                "lepton_vector_delta": float(wilsons.lepton_vector_delta),
                "lepton_axial_delta": float(wilsons.lepton_axial_delta),
                "c9_np": complex(wilsons.c9_np),
                "c10_np": complex(wilsons.c10_np),
                "c9p_np": complex(wilsons.c9p_np),
                "c10p_np": complex(wilsons.c10p_np),
            }
        )

    return RareCharmDToPiMuMuBranchingResult(
        model_label=RARE_CHARM_DTOPI_MUMU_MODEL_V1,
        input_bundle=p.input_bundle,
        transition_key="c_u",
        lepton_key="mu",
        branching_fraction=float(branching),
        sm_branching_fraction=float(sm_branching),
        np_shift_branching_fraction=float(branching - sm_branching),
        q2_min_gev2=float(q2_min),
        q2_max_gev2=float(q2_max),
        c9_semileptonic_np=complex(c9_semileptonic),
        c10_semileptonic_np=complex(c10_semileptonic - p.rare_charm.c10_sm),
        c9_np=complex(c9_np),
        c10_np=complex(c10_np),
        c9p_np=complex(c9p_np),
        c10p_np=complex(c10p_np),
        lambda_b=factors.lambda_b,
        wilsons=wilsons,
        diagnostics=diagnostics,
    )


__all__ = [
    "RARE_CHARM_DTOPI_MUMU_MODEL_V1",
    "RARE_CHARM_DTOPI_MUMU_INPUT_BUNDLE_V1",
    "RARE_CHARM_DTOPI_MUMU_FORM_FACTOR_MODEL_V1",
    "RARE_CHARM_DTOPI_MUMU_OPERATOR_CONVENTION",
    "RARE_CHARM_DTOPI_MUMU_PARAMETRIZATION_CITATION",
    "RARE_CHARM_DTOPI_MUMU_RESONANCE_LIMITATION_V1",
    "RARE_CHARM_DTOPI_MUMU_Q2_TREATMENT_V1",
    "RareCharmDToPiFormFactorInputs",
    "RareCharmDToPiMuMuInputs",
    "RareCharmDToPiMuMuBranchingResult",
    "default_dtopi_mumu_inputs",
    "dtopi_fplus",
    "dtopi_fzero",
    "dtopi_mumu_q2_range",
    "dtopi_mumu_differential_branching_fraction",
    "dtopi_mumu_sm",
    "evaluate_dplus_to_piplus_mumu",
]
