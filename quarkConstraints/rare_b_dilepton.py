"""Short-distance ``b -> q l+ l-`` machinery for rare beauty decays.

Physics convention
------------------
The leptonic ``B_q -> mu+ mu-`` rate is evaluated in the standard Buras
``b -> q l l`` Hamiltonian convention

    H_eff = -4 G_F/sqrt(2) lambda_t alpha_e/(4 pi)
            [ C9 O9 + C10 O10 + C9' O9' + C10' O10' + ... ],

where ``lambda_t = V_tb V_tq^*`` and
``O10 = (qbar gamma_mu P_L b)(mubar gamma^mu gamma5 mu)``.  For the pure
leptonic rate, the SM and the v1 prediction are dominated by
``C10 - C10'``; ``C9`` is computed and returned for shared ``b -> q l l``
reuse, but it does not contribute to ``B_q -> mu+ mu-``.

The SM branching fraction is

    tau_Bq G_F^2 alpha_e^2 f_Bq^2 m_Bq m_mu^2 |lambda_t|^2
    / (16 pi^3) * sqrt(1 - 4 m_mu^2 / m_Bq^2) * |C10_SM|^2,

with the standard time-integration factor
``(1 + A_DeltaGamma y_q) / (1 - y_q^2)``.  ``A_DeltaGamma`` is computed from
the complex scalar and pseudoscalar amplitudes.  The pure-C10 SM-like limit has
``A_DeltaGamma = +1`` and gives ``BR(B_s -> mu mu) ~= 3.65e-9`` with the repo
CKM target, compatible with the B005 Czaja-Misiak 2024 anchor ``3.64(12)e-9``.

RS matching assumption
----------------------
NEEDS-HUMAN-PHYSICS: the current ``ParameterPoint`` carries quark mass-basis
``KK-gluon``-style coupling matrices, but not the electroweak KK/Z/Z' tower,
muon localization, Higgs/radion, or scalar/pseudoscalar couplings required for
a model-complete RS ``b -> q l l`` prediction.  The v1 matching below is a
documented Z/KK-penguin proxy: divide the supplied quark coupling by ``g_s`` to
keep the flavor-overlap structure, couple that overlap to one Z-like neutral
boson at ``M_KK``, use SM-Z muon vector/axial charges, and match the result
onto ``C9``/``C10`` in the Hamiltonian above.  Scalar and pseudoscalar Wilsons
are present in the data model for future reuse but are zero in this proxy.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Mapping

import numpy as np

from .couplings import QuarkMassBasisCouplings
from .deltaf2 import F_BD, F_BS, M_BD, M_BS, M_B_QUARK, M_D_QUARK_BD, M_S_QUARK_BS
from .model import RotationParameters, ckm_like_unitary

RARE_B_DILEPTON_MODEL_V1 = "rare_b_dilepton_buras_c10_rs_proxy_v1"
RARE_B_DILEPTON_OPERATOR_CONVENTION = (
    "H_eff=-4 G_F/sqrt(2) lambda_t alpha/(4 pi) "
    "[C9 O9 + C10 O10 + C9' O9' + C10' O10']; "
    "O10=(qbar gamma_mu P_L b)(mubar gamma^mu gamma5 mu)"
)
RARE_B_DILEPTON_INPUT_BUNDLE_V1 = "rare_b_dilepton_sm_inputs_buras_repo_ckm_v1"
RARE_B_DILEPTON_RS_MATCHING_ASSUMPTION_V1 = (
    "NEEDS-HUMAN-PHYSICS: quark KK-gluon mass-basis couplings are used as "
    "neutral-current flavor-overlap proxies; a single Z-like EW KK/Z-penguin "
    "boson with SM-Z muon vector/axial charges stands in for the full RS "
    "EW KK/Z/Z', lepton, Higgs/radion, scalar and pseudoscalar matching."
)


@dataclass(frozen=True)
class RareBDileptonMesonInputs:
    """Hadronic and lifetime inputs for one ``B_q -> mu+ mu-`` channel."""

    transition_key: str
    display_name: str
    light_down_index: int
    meson_mass_gev: float
    decay_constant_gev: float
    lifetime_ps: float
    light_quark_mass_gev: float
    width_difference_y: float = 0.0
    source: str = "PDG/FLAG constants mirrored from quarkConstraints.deltaf2"

    def __post_init__(self) -> None:
        if self.transition_key not in {"b_s", "b_d"}:
            raise ValueError("transition_key must be 'b_s' or 'b_d'")
        if self.light_down_index not in {0, 1}:
            raise ValueError("light_down_index must select d=0 or s=1")
        for name in (
            "meson_mass_gev",
            "decay_constant_gev",
            "lifetime_ps",
            "light_quark_mass_gev",
        ):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be positive and finite")
        if abs(float(self.width_difference_y)) >= 1.0:
            raise ValueError("width_difference_y must satisfy |y| < 1")


@dataclass(frozen=True)
class RareBDileptonSMInputs:
    """Numerical inputs for the shared rare-beauty dilepton module."""

    input_bundle: str = RARE_B_DILEPTON_INPUT_BUNDLE_V1
    gf_gev_minus2: float = 1.1663787e-5
    alpha_em_mz: float = 1.0 / 127.952
    sin2_theta_w: float = 0.23122
    c10_sm: float = -4.103
    muon_mass_gev: float = 0.1056583745
    hbar_gev_s: float = 6.582119569e-25
    theta12: float = 0.2274
    theta13: float = 0.00368
    theta23: float = 0.0415
    delta: float = 1.196
    bs: RareBDileptonMesonInputs = field(
        default_factory=lambda: RareBDileptonMesonInputs(
            transition_key="b_s",
            display_name="B_s -> mu+ mu-",
            light_down_index=1,
            meson_mass_gev=M_BS,
            decay_constant_gev=F_BS,
            lifetime_ps=1.515,
            light_quark_mass_gev=M_S_QUARK_BS,
            width_difference_y=0.0645,
            source=(
                "F_Bs and m_Bs from quarkConstraints.deltaf2 FLAG/PDG; "
                "tau_Bs and y_s PDG/HFLAV-era defaults"
            ),
        )
    )
    bd: RareBDileptonMesonInputs = field(
        default_factory=lambda: RareBDileptonMesonInputs(
            transition_key="b_d",
            display_name="B_d -> mu+ mu-",
            light_down_index=0,
            meson_mass_gev=M_BD,
            decay_constant_gev=F_BD,
            lifetime_ps=1.519,
            light_quark_mass_gev=M_D_QUARK_BD,
            width_difference_y=0.0,
            source=(
                "F_Bd and m_Bd from quarkConstraints.deltaf2 FLAG/PDG; "
                "tau_Bd PDG-era default"
            ),
        )
    )
    constants_citation: str = (
        "Buras b->sll Hamiltonian convention; Czaja-Misiak arXiv:2407.03810 "
        "for B_s -> mu mu SM validation; repo CKM target "
        "quarkConstraints.modern.inputs.ModernDefaultCKMTarget"
    )

    def __post_init__(self) -> None:
        for name in (
            "gf_gev_minus2",
            "alpha_em_mz",
            "sin2_theta_w",
            "muon_mass_gev",
            "hbar_gev_s",
        ):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be positive and finite")
        if not 0.0 < float(self.sin2_theta_w) < 1.0:
            raise ValueError("sin2_theta_w must lie between zero and one")
        if not math.isfinite(float(self.c10_sm)) or float(self.c10_sm) == 0.0:
            raise ValueError("c10_sm must be finite and non-zero")

    def meson(self, transition: str) -> RareBDileptonMesonInputs:
        """Return the meson-input block for ``transition``."""
        if transition == "b_s":
            return self.bs
        if transition == "b_d":
            return self.bd
        raise ValueError(f"unsupported rare-b transition {transition!r}")


@dataclass(frozen=True)
class RareBDileptonCKMFactors:
    """CKM factors entering ``b -> q l l``."""

    transition_key: str
    lambda_wolfenstein: float
    lambda_t: complex
    matrix: tuple[tuple[complex, ...], ...]


@dataclass(frozen=True)
class RareBDileptonWilsonCoefficients:
    """Leading ``b -> q mu+ mu-`` Wilson proxy from mass-basis couplings."""

    model_label: str
    operator_convention: str
    matching_assumption: str
    transition_key: str
    M_KK: float
    matching_scale: float
    lambda_t: complex
    left_qb_coupling: complex
    right_qb_coupling: complex
    left_qb_overlap: complex
    right_qb_overlap: complex
    left_quark_delta: complex
    right_quark_delta: complex
    muon_left_delta: float
    muon_right_delta: float
    muon_vector_delta: float
    muon_axial_delta: float
    c9_np: complex
    c10_np: complex
    c9p_np: complex
    c10p_np: complex
    cs_np: complex = 0.0j
    csp_np: complex = 0.0j
    cp_np: complex = 0.0j
    cpp_np: complex = 0.0j

    @property
    def c9_effective_np(self) -> complex:
        return complex(self.c9_np + self.c9p_np)

    @property
    def c10_leptonic_np(self) -> complex:
        """Combination entering ``B_q -> mu mu``: ``C10 - C10'``."""
        return complex(self.c10_np - self.c10p_np)

    @property
    def wilsons(self) -> Mapping[str, complex]:
        return {
            "C9_NP": complex(self.c9_np),
            "C10_NP": complex(self.c10_np),
            "C9p_NP": complex(self.c9p_np),
            "C10p_NP": complex(self.c10p_np),
            "CS_NP": complex(self.cs_np),
            "CSp_NP": complex(self.csp_np),
            "CP_NP": complex(self.cp_np),
            "CPp_NP": complex(self.cpp_np),
        }


@dataclass(frozen=True)
class RareBLeptonicBranchingResult:
    """Branching-ratio prediction for ``B_q -> mu+ mu-``."""

    model_label: str
    input_bundle: str
    transition_key: str
    branching_fraction: float
    sm_branching_fraction: float
    np_shift_branching_fraction: float
    prompt_branching_fraction: float
    sm_prompt_branching_fraction: float
    time_integrated_factor: float
    a_delta_gamma: float
    c10_total: complex
    c10_sm: float
    c10_leptonic_np: complex
    c9_effective_np: complex
    scalar_amplitude: complex
    pseudoscalar_amplitude: complex
    lambda_t: complex
    wilsons: RareBDileptonWilsonCoefficients | None = None
    diagnostics: Mapping[str, float | complex | str | bool] = field(default_factory=dict)


def default_sm_inputs() -> RareBDileptonSMInputs:
    """Return the repo-owned default rare-beauty dilepton input bundle."""
    return RareBDileptonSMInputs()


def ckm_factors(
    transition: str = "b_s",
    inputs: RareBDileptonSMInputs | None = None,
) -> RareBDileptonCKMFactors:
    """Return ``lambda_t = V_tb V_tq^*`` for ``b -> q``."""
    p = default_sm_inputs() if inputs is None else inputs
    meson = p.meson(transition)
    matrix_np = ckm_like_unitary(
        RotationParameters(
            theta12=p.theta12,
            theta13=p.theta13,
            theta23=p.theta23,
            delta=p.delta,
        )
    )
    lam = float(abs(matrix_np[0, 1]))
    lambda_t = complex(matrix_np[2, 2] * np.conjugate(matrix_np[2, meson.light_down_index]))
    matrix = tuple(tuple(complex(entry) for entry in row) for row in matrix_np)
    if abs(lambda_t) <= 0.0:
        raise ValueError(f"{transition}: lambda_t must be non-zero")
    return RareBDileptonCKMFactors(
        transition_key=transition,
        lambda_wolfenstein=lam,
        lambda_t=lambda_t,
        matrix=matrix,
    )


def _matrix_entry(source: object, matrix_name: str, i: int, j: int) -> complex:
    matrix = np.asarray(getattr(source, matrix_name), dtype=np.complex128)
    if matrix.shape != (3, 3):
        raise ValueError(f"{matrix_name} must have shape (3, 3)")
    value = complex(matrix[i, j])
    if not math.isfinite(value.real) or not math.isfinite(value.imag):
        raise ValueError(f"{matrix_name}[{i},{j}] must be finite")
    return value


def _positive_float(value: object, name: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0.0:
        raise ValueError(f"{name} must be positive and finite")
    return number


def _weak_couplings(
    inputs: RareBDileptonSMInputs,
) -> tuple[float, float, float, float, float]:
    g_weak = math.sqrt(4.0 * math.pi * inputs.alpha_em_mz / inputs.sin2_theta_w)
    cos_theta_w = math.sqrt(1.0 - inputs.sin2_theta_w)
    g_z = g_weak / cos_theta_w
    mu_left = g_z * (-0.5 + inputs.sin2_theta_w)
    mu_right = g_z * inputs.sin2_theta_w
    mu_vector = mu_left + mu_right
    mu_axial = mu_right - mu_left
    neutral_delta = g_z / 2.0
    return (
        float(neutral_delta),
        float(mu_left),
        float(mu_right),
        float(mu_vector),
        float(mu_axial),
    )


def _wilson_prefactor(
    *,
    lambda_t: complex,
    m_kk_gev: float,
    inputs: RareBDileptonSMInputs,
) -> complex:
    return complex(
        -math.pi
        / (
            math.sqrt(2.0)
            * inputs.gf_gev_minus2
            * inputs.alpha_em_mz
            * lambda_t
            * m_kk_gev**2
        )
    )


def compute_rare_b_dilepton_wilsons(
    source: QuarkMassBasisCouplings,
    *,
    transition: str = "b_s",
    m_kk_gev: float | None = None,
    inputs: RareBDileptonSMInputs | None = None,
) -> RareBDileptonWilsonCoefficients:
    """Match mass-basis couplings onto the v1 ``b -> q mu mu`` proxy."""
    p = default_sm_inputs() if inputs is None else inputs
    meson = p.meson(transition)
    resolved_m_kk = _positive_float(
        getattr(source, "M_KK") if m_kk_gev is None else m_kk_gev,
        "m_kk_gev",
    )
    q_idx = meson.light_down_index
    b_idx = 2
    left_qb = _matrix_entry(source, "left_down", q_idx, b_idx)
    right_qb = _matrix_entry(source, "right_down", q_idx, b_idx)
    g_s = _positive_float(getattr(source, "g_s", 1.0), "g_s")
    left_overlap = left_qb / g_s
    right_overlap = right_qb / g_s

    factors = ckm_factors(transition, p)
    neutral_delta, mu_left, mu_right, mu_vector, mu_axial = _weak_couplings(p)
    left_quark_delta = neutral_delta * left_overlap
    right_quark_delta = neutral_delta * right_overlap
    prefactor = _wilson_prefactor(
        lambda_t=factors.lambda_t,
        m_kk_gev=resolved_m_kk,
        inputs=p,
    )
    c9_np = prefactor * left_quark_delta * mu_vector
    c10_np = prefactor * left_quark_delta * mu_axial
    c9p_np = prefactor * right_quark_delta * mu_vector
    c10p_np = prefactor * right_quark_delta * mu_axial

    return RareBDileptonWilsonCoefficients(
        model_label=RARE_B_DILEPTON_MODEL_V1,
        operator_convention=RARE_B_DILEPTON_OPERATOR_CONVENTION,
        matching_assumption=RARE_B_DILEPTON_RS_MATCHING_ASSUMPTION_V1,
        transition_key=transition,
        M_KK=resolved_m_kk,
        matching_scale=resolved_m_kk,
        lambda_t=factors.lambda_t,
        left_qb_coupling=left_qb,
        right_qb_coupling=right_qb,
        left_qb_overlap=complex(left_overlap),
        right_qb_overlap=complex(right_overlap),
        left_quark_delta=complex(left_quark_delta),
        right_quark_delta=complex(right_quark_delta),
        muon_left_delta=mu_left,
        muon_right_delta=mu_right,
        muon_vector_delta=mu_vector,
        muon_axial_delta=mu_axial,
        c9_np=complex(c9_np),
        c10_np=complex(c10_np),
        c9p_np=complex(c9p_np),
        c10p_np=complex(c10p_np),
    )


def _tau_ps_to_gev_inverse(tau_ps: float, hbar_gev_s: float) -> float:
    return float(tau_ps * 1.0e-12 / hbar_gev_s)


def _a_delta_gamma_from_amplitudes(
    *,
    scalar_amplitude: complex,
    pseudoscalar_amplitude: complex,
) -> float:
    """Return amplitude-dependent ``A_DeltaGamma`` with no extra mixing phase."""
    scalar = complex(scalar_amplitude)
    pseudoscalar = complex(pseudoscalar_amplitude)
    denominator = abs(pseudoscalar) ** 2 + abs(scalar) ** 2
    if not math.isfinite(denominator):
        raise ValueError("B_q -> mu mu amplitude norm must be finite")
    if denominator <= 0.0:
        return 1.0
    value = ((pseudoscalar * pseudoscalar).real - (scalar * scalar).real) / denominator
    return float(min(1.0, max(-1.0, value)))


def _time_integrated_factor(
    meson: RareBDileptonMesonInputs,
    *,
    a_delta_gamma: float,
) -> float:
    y = float(meson.width_difference_y)
    return float((1.0 + float(a_delta_gamma) * y) / (1.0 - y * y))


def _branching_from_wilson_values(
    *,
    transition: str,
    c10_np: complex,
    c10p_np: complex,
    cs_np: complex,
    csp_np: complex,
    cp_np: complex,
    cpp_np: complex,
    inputs: RareBDileptonSMInputs,
) -> tuple[float, float, float, float, complex, complex, complex, RareBDileptonCKMFactors]:
    meson = inputs.meson(transition)
    factors = ckm_factors(transition, inputs)
    m_b = M_B_QUARK
    m_q = meson.light_quark_mass_gev
    m_mu = inputs.muon_mass_gev
    m_bq = meson.meson_mass_gev
    beta2 = 1.0 - 4.0 * m_mu * m_mu / (m_bq * m_bq)
    if beta2 <= 0.0:
        raise ValueError(f"{transition}: leptonic phase-space beta^2 <= 0")
    beta = math.sqrt(beta2)
    scalar_scale = m_bq * m_bq / (2.0 * m_mu) * (m_b / (m_b + m_q))

    c10_total = complex(inputs.c10_sm + c10_np)
    pseudoscalar_amplitude = c10_total - complex(c10p_np) + scalar_scale * (
        complex(cp_np) - complex(cpp_np)
    )
    scalar_amplitude = beta * scalar_scale * (complex(cs_np) - complex(csp_np))
    tau = _tau_ps_to_gev_inverse(meson.lifetime_ps, inputs.hbar_gev_s)
    prefactor = (
        tau
        * inputs.gf_gev_minus2**2
        * inputs.alpha_em_mz**2
        / (16.0 * math.pi**3)
        * meson.decay_constant_gev**2
        * m_bq
        * m_mu**2
        * abs(factors.lambda_t) ** 2
        * beta
    )
    prompt = prefactor * (
        abs(pseudoscalar_amplitude) ** 2 + abs(scalar_amplitude) ** 2
    )
    a_delta_gamma = _a_delta_gamma_from_amplitudes(
        scalar_amplitude=scalar_amplitude,
        pseudoscalar_amplitude=pseudoscalar_amplitude,
    )
    time_factor = _time_integrated_factor(meson, a_delta_gamma=a_delta_gamma)
    return (
        float(prompt * time_factor),
        float(prompt),
        float(time_factor),
        float(a_delta_gamma),
        c10_total,
        scalar_amplitude,
        pseudoscalar_amplitude,
        factors,
    )


def sm_branching_fraction(
    transition: str = "b_s",
    inputs: RareBDileptonSMInputs | None = None,
) -> RareBLeptonicBranchingResult:
    """Evaluate the SM-limit ``B_q -> mu+ mu-`` branching fraction."""
    p = default_sm_inputs() if inputs is None else inputs
    (
        br,
        prompt,
        time_factor,
        a_delta_gamma,
        c10_total,
        scalar_amp,
        pseudoscalar_amp,
        factors,
    ) = _branching_from_wilson_values(
        transition=transition,
        c10_np=0.0j,
        c10p_np=0.0j,
        cs_np=0.0j,
        csp_np=0.0j,
        cp_np=0.0j,
        cpp_np=0.0j,
        inputs=p,
    )
    return RareBLeptonicBranchingResult(
        model_label=RARE_B_DILEPTON_MODEL_V1,
        input_bundle=p.input_bundle,
        transition_key=transition,
        branching_fraction=br,
        sm_branching_fraction=br,
        np_shift_branching_fraction=0.0,
        prompt_branching_fraction=prompt,
        sm_prompt_branching_fraction=prompt,
        time_integrated_factor=time_factor,
        a_delta_gamma=a_delta_gamma,
        c10_total=c10_total,
        c10_sm=float(p.c10_sm),
        c10_leptonic_np=0.0j,
        c9_effective_np=0.0j,
        scalar_amplitude=scalar_amp,
        pseudoscalar_amplitude=pseudoscalar_amp,
        lambda_t=factors.lambda_t,
        diagnostics={
            "lambda_wolfenstein": factors.lambda_wolfenstein,
            "lambda_t": factors.lambda_t,
            "c10_sm": float(p.c10_sm),
            "a_delta_gamma": float(a_delta_gamma),
            "time_integrated_factor": time_factor,
            "constants_citation": p.constants_citation,
        },
    )


def evaluate_bq_to_mumu(
    source: QuarkMassBasisCouplings | RareBDileptonWilsonCoefficients | None = None,
    *,
    transition: str = "b_s",
    m_kk_gev: float | None = None,
    inputs: RareBDileptonSMInputs | None = None,
) -> RareBLeptonicBranchingResult:
    """Evaluate ``BR(B_q -> mu+ mu-)`` in the SM or with the v1 RS proxy."""
    p = default_sm_inputs() if inputs is None else inputs
    wilsons: RareBDileptonWilsonCoefficients | None
    if source is None:
        wilsons = None
        c10_np = c10p_np = cs_np = csp_np = cp_np = cpp_np = 0.0j
        c9_effective = 0.0j
    elif isinstance(source, RareBDileptonWilsonCoefficients):
        wilsons = source
        c10_np = source.c10_np
        c10p_np = source.c10p_np
        cs_np = source.cs_np
        csp_np = source.csp_np
        cp_np = source.cp_np
        cpp_np = source.cpp_np
        c9_effective = source.c9_effective_np
    else:
        wilsons = compute_rare_b_dilepton_wilsons(
            source,
            transition=transition,
            m_kk_gev=m_kk_gev,
            inputs=p,
        )
        c10_np = wilsons.c10_np
        c10p_np = wilsons.c10p_np
        cs_np = wilsons.cs_np
        csp_np = wilsons.csp_np
        cp_np = wilsons.cp_np
        cpp_np = wilsons.cpp_np
        c9_effective = wilsons.c9_effective_np

    (
        br,
        prompt,
        time_factor,
        a_delta_gamma,
        c10_total,
        scalar_amp,
        pseudoscalar_amp,
        factors,
    ) = _branching_from_wilson_values(
        transition=transition,
        c10_np=c10_np,
        c10p_np=c10p_np,
        cs_np=cs_np,
        csp_np=csp_np,
        cp_np=cp_np,
        cpp_np=cpp_np,
        inputs=p,
    )
    sm = sm_branching_fraction(transition, p)
    diagnostics: dict[str, float | complex | str | bool] = {
        "lambda_wolfenstein": factors.lambda_wolfenstein,
        "lambda_t": factors.lambda_t,
        "c10_sm": float(p.c10_sm),
        "c10_total": complex(c10_total),
        "c10_leptonic_np": complex(c10_np - c10p_np),
        "c9_effective_np": complex(c9_effective),
        "c9_does_not_enter_leptonic_rate": True,
        "scalar_amplitude": complex(scalar_amp),
        "pseudoscalar_amplitude": complex(pseudoscalar_amp),
        "prompt_branching_fraction": float(prompt),
        "a_delta_gamma": float(a_delta_gamma),
        "time_integrated_factor": float(time_factor),
        "matching_assumption": RARE_B_DILEPTON_RS_MATCHING_ASSUMPTION_V1,
    }
    if wilsons is not None:
        diagnostics.update(
            {
                "m_kk_gev": float(wilsons.M_KK),
                "matching_scale_gev": float(wilsons.matching_scale),
                "left_qb_coupling": complex(wilsons.left_qb_coupling),
                "right_qb_coupling": complex(wilsons.right_qb_coupling),
                "left_qb_overlap": complex(wilsons.left_qb_overlap),
                "right_qb_overlap": complex(wilsons.right_qb_overlap),
                "left_quark_delta": complex(wilsons.left_quark_delta),
                "right_quark_delta": complex(wilsons.right_quark_delta),
                "muon_left_delta": float(wilsons.muon_left_delta),
                "muon_right_delta": float(wilsons.muon_right_delta),
                "muon_vector_delta": float(wilsons.muon_vector_delta),
                "muon_axial_delta": float(wilsons.muon_axial_delta),
                "c9_np": complex(wilsons.c9_np),
                "c10_np": complex(wilsons.c10_np),
                "c9p_np": complex(wilsons.c9p_np),
                "c10p_np": complex(wilsons.c10p_np),
            }
        )

    return RareBLeptonicBranchingResult(
        model_label=RARE_B_DILEPTON_MODEL_V1,
        input_bundle=p.input_bundle,
        transition_key=transition,
        branching_fraction=br,
        sm_branching_fraction=float(sm.branching_fraction),
        np_shift_branching_fraction=float(br - sm.branching_fraction),
        prompt_branching_fraction=prompt,
        sm_prompt_branching_fraction=float(sm.prompt_branching_fraction),
        time_integrated_factor=time_factor,
        a_delta_gamma=a_delta_gamma,
        c10_total=c10_total,
        c10_sm=float(p.c10_sm),
        c10_leptonic_np=complex(c10_np - c10p_np),
        c9_effective_np=complex(c9_effective),
        scalar_amplitude=scalar_amp,
        pseudoscalar_amplitude=pseudoscalar_amp,
        lambda_t=factors.lambda_t,
        wilsons=wilsons,
        diagnostics=diagnostics,
    )


__all__ = [
    "RARE_B_DILEPTON_MODEL_V1",
    "RARE_B_DILEPTON_OPERATOR_CONVENTION",
    "RARE_B_DILEPTON_INPUT_BUNDLE_V1",
    "RARE_B_DILEPTON_RS_MATCHING_ASSUMPTION_V1",
    "RareBDileptonMesonInputs",
    "RareBDileptonSMInputs",
    "RareBDileptonCKMFactors",
    "RareBDileptonWilsonCoefficients",
    "RareBLeptonicBranchingResult",
    "default_sm_inputs",
    "ckm_factors",
    "compute_rare_b_dilepton_wilsons",
    "sm_branching_fraction",
    "evaluate_bq_to_mumu",
]
