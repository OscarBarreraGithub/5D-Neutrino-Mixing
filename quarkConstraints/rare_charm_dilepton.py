"""Short-distance ``c -> u l+ l-`` machinery for rare charm decays.

Physics convention
------------------
The neutral leptonic ``D0 -> l+ l-`` short-distance rate is evaluated in the
standard rare-charm Hamiltonian normalization

    H_eff = -4 G_F/sqrt(2) lambda_b alpha_e/(4 pi)
            [ C9 O9 + C10 O10 + C9' O9' + C10' O10' + ... ],

where ``lambda_b = V_cb^* V_ub`` and
``O10 = (ubar gamma_mu P_L c)(lbar gamma^mu gamma5 l)``.  For the pure
leptonic rate, ``C9`` is diagnostic-only and the short-distance amplitude is
controlled by ``C10 - C10'`` plus optional scalar/pseudoscalar coefficients
kept in the data model for later reuse.

Charm caveat
------------
``D0 -> mu+ mu-`` is long-distance dominated in the Standard Model.  This
module deliberately owns only the short-distance piece that can be constrained
cleanly by the experimental upper bound; long-distance context is carried by
the catalog constraint sidecar.

RS matching assumption
----------------------
NEEDS-HUMAN-PHYSICS: the current ``ParameterPoint`` carries quark mass-basis
``KK-gluon``-style coupling matrices, but not the electroweak KK/Z/Z' tower,
charged-lepton localization, Higgs/radion, or scalar/pseudoscalar couplings
required for a model-complete RS ``c -> u l l`` prediction.  The v1 matching
below is a documented Z/KK-penguin proxy: divide the supplied quark coupling by
``g_s`` to keep the flavor-overlap structure, couple that overlap to one
Z-like neutral boson at ``M_KK``, use SM-Z charged-lepton vector/axial charges,
and match the result onto ``C9``/``C10`` in the Hamiltonian above.  The chiral
lepton currents are decomposed as ``lbar gamma_mu P_{L,R} l = (V +/- A)/2``,
giving the tree-level proxy normalization ``+pi Delta_uc Delta_{V,A}/(2
sqrt(2) G_F alpha lambda_b M_KK^2)`` in this module's Lagrangian convention.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Mapping

import numpy as np

from .couplings import QuarkMassBasisCouplings
from .deltaf2 import F_D, M_C_QUARK, M_D0, M_U_QUARK
from .model import RotationParameters, ckm_like_unitary

RARE_CHARM_DILEPTON_MODEL_V1 = "rare_charm_dilepton_cu_c10_rs_proxy_v1"
RARE_CHARM_DILEPTON_OPERATOR_CONVENTION = (
    "H_eff=-4 G_F/sqrt(2) lambda_b alpha/(4 pi) "
    "[C9 O9 + C10 O10 + C9' O9' + C10' O10']; "
    "lambda_b=V_cb^* V_ub; "
    "O10=(ubar gamma_mu P_L c)(lbar gamma^mu gamma5 l)"
)
RARE_CHARM_DILEPTON_INPUT_BUNDLE_V1 = "rare_charm_dilepton_sm_inputs_charm_sd_v1"
RARE_CHARM_DILEPTON_RS_MATCHING_ASSUMPTION_V1 = (
    "NEEDS-HUMAN-PHYSICS: quark KK-gluon mass-basis couplings are used as "
    "neutral-current flavor-overlap proxies; a single Z-like EW KK/Z-penguin "
    "boson with SM-Z charged-lepton vector/axial charges is WET-matched with "
    "+pi/(2 sqrt(2) G_F alpha lambda_b M_KK^2); it stands in for the full "
    "RS EW KK/Z/Z', lepton, Higgs/radion, scalar and pseudoscalar c->u l l "
    "matching."
)


@dataclass(frozen=True)
class RareCharmLeptonInputs:
    """Charged-lepton inputs for one ``D0 -> l+ l-`` channel."""

    lepton_key: str
    display_name: str
    mass_gev: float
    source: str

    def __post_init__(self) -> None:
        if self.lepton_key not in {"e", "mu"}:
            raise ValueError("lepton_key must be 'e' or 'mu'")
        value = float(self.mass_gev)
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError("lepton mass must be positive and finite")


@dataclass(frozen=True)
class RareCharmMesonInputs:
    """Hadronic and lifetime inputs for a neutral charm leptonic mode."""

    transition_key: str
    display_name: str
    light_up_index: int
    charm_index: int
    meson_mass_gev: float
    decay_constant_gev: float
    lifetime_ps: float
    charm_quark_mass_gev: float
    light_quark_mass_gev: float
    source: str

    def __post_init__(self) -> None:
        if self.transition_key != "c_u":
            raise ValueError("only the c_u transition is currently supported")
        if self.light_up_index != 0 or self.charm_index != 1:
            raise ValueError("c_u must use up index 0 and charm index 1")
        for name in (
            "meson_mass_gev",
            "decay_constant_gev",
            "lifetime_ps",
            "charm_quark_mass_gev",
            "light_quark_mass_gev",
        ):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be positive and finite")


@dataclass(frozen=True)
class RareCharmDileptonSMInputs:
    """Numerical inputs for the shared rare-charm dilepton module."""

    input_bundle: str = RARE_CHARM_DILEPTON_INPUT_BUNDLE_V1
    gf_gev_minus2: float = 1.1663787e-5
    alpha_em_mz: float = 1.0 / 127.952
    sin2_theta_w: float = 0.23122
    c10_sm: float = 0.0
    hbar_gev_s: float = 6.582119569e-25
    theta12: float = 0.2274
    theta13: float = 0.00368
    theta23: float = 0.0415
    delta: float = 1.196
    d0: RareCharmMesonInputs = field(
        default_factory=lambda: RareCharmMesonInputs(
            transition_key="c_u",
            display_name="D0 -> l+ l-",
            light_up_index=0,
            charm_index=1,
            meson_mass_gev=M_D0,
            decay_constant_gev=F_D,
            lifetime_ps=0.4101,
            charm_quark_mass_gev=M_C_QUARK,
            light_quark_mass_gev=M_U_QUARK,
            source=(
                "F_D, m_D0, m_c, m_u from quarkConstraints.deltaf2 "
                "FLAG/PDG constants; tau_D0 PDG-era default"
            ),
        )
    )
    muon: RareCharmLeptonInputs = field(
        default_factory=lambda: RareCharmLeptonInputs(
            lepton_key="mu",
            display_name="mu+ mu-",
            mass_gev=0.1056583745,
            source="PDG charged-lepton mass",
        )
    )
    electron: RareCharmLeptonInputs = field(
        default_factory=lambda: RareCharmLeptonInputs(
            lepton_key="e",
            display_name="e+ e-",
            mass_gev=0.00051099895,
            source="PDG charged-lepton mass",
        )
    )
    constants_citation: str = (
        "Rare-charm c->u l l Hamiltonian convention; Burdman-Golowich-"
        "Hewett-Pakvasa long-distance caveat; repo CKM target "
        "quarkConstraints.modern.inputs.ModernDefaultCKMTarget"
    )

    def __post_init__(self) -> None:
        for name in (
            "gf_gev_minus2",
            "alpha_em_mz",
            "sin2_theta_w",
            "hbar_gev_s",
        ):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be positive and finite")
        if not 0.0 < float(self.sin2_theta_w) < 1.0:
            raise ValueError("sin2_theta_w must lie between zero and one")
        if not math.isfinite(float(self.c10_sm)):
            raise ValueError("c10_sm must be finite")

    def lepton(self, lepton: str) -> RareCharmLeptonInputs:
        """Return the charged-lepton input block for ``lepton``."""
        if lepton in {"mu", "muon"}:
            return self.muon
        if lepton in {"e", "electron"}:
            return self.electron
        raise ValueError(f"unsupported rare-charm lepton {lepton!r}")


@dataclass(frozen=True)
class RareCharmDileptonCKMFactors:
    """CKM factors entering ``c -> u l l``."""

    transition_key: str
    lambda_wolfenstein: float
    lambda_b: complex
    matrix: tuple[tuple[complex, ...], ...]


@dataclass(frozen=True)
class RareCharmDileptonWilsonCoefficients:
    """Leading ``c -> u l+ l-`` Wilson proxy from mass-basis couplings."""

    model_label: str
    operator_convention: str
    matching_assumption: str
    transition_key: str
    M_KK: float
    matching_scale: float
    lambda_b: complex
    left_uc_coupling: complex
    right_uc_coupling: complex
    left_uc_overlap: complex
    right_uc_overlap: complex
    left_quark_delta: complex
    right_quark_delta: complex
    lepton_left_delta: float
    lepton_right_delta: float
    lepton_vector_delta: float
    lepton_axial_delta: float
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
        """Combination entering ``D0 -> l l``: ``C10 - C10'``."""
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
class RareCharmLeptonicBranchingResult:
    """Short-distance branching-ratio prediction for ``D0 -> l+ l-``."""

    model_label: str
    input_bundle: str
    transition_key: str
    lepton_key: str
    branching_fraction: float
    sm_branching_fraction: float
    np_shift_branching_fraction: float
    c10_total: complex
    c10_sm: float
    c10_leptonic_np: complex
    c9_effective_np: complex
    scalar_amplitude: complex
    pseudoscalar_amplitude: complex
    lambda_b: complex
    wilsons: RareCharmDileptonWilsonCoefficients | None = None
    diagnostics: Mapping[str, float | complex | str | bool] = field(default_factory=dict)


def default_sm_inputs() -> RareCharmDileptonSMInputs:
    """Return the repo-owned default rare-charm dilepton input bundle."""
    return RareCharmDileptonSMInputs()


def ckm_factors(
    transition: str = "c_u",
    inputs: RareCharmDileptonSMInputs | None = None,
) -> RareCharmDileptonCKMFactors:
    """Return ``lambda_b = V_cb^* V_ub`` for ``c -> u``."""
    if transition != "c_u":
        raise ValueError(f"unsupported rare-charm transition {transition!r}")
    p = default_sm_inputs() if inputs is None else inputs
    matrix_np = ckm_like_unitary(
        RotationParameters(
            theta12=p.theta12,
            theta13=p.theta13,
            theta23=p.theta23,
            delta=p.delta,
        )
    )
    lam = float(abs(matrix_np[0, 1]))
    lambda_b = complex(np.conjugate(matrix_np[1, 2]) * matrix_np[0, 2])
    matrix = tuple(tuple(complex(entry) for entry in row) for row in matrix_np)
    if abs(lambda_b) <= 0.0:
        raise ValueError("c_u: lambda_b must be non-zero")
    return RareCharmDileptonCKMFactors(
        transition_key=transition,
        lambda_wolfenstein=lam,
        lambda_b=lambda_b,
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
    inputs: RareCharmDileptonSMInputs,
) -> tuple[float, float, float, float, float]:
    g_weak = math.sqrt(4.0 * math.pi * inputs.alpha_em_mz / inputs.sin2_theta_w)
    cos_theta_w = math.sqrt(1.0 - inputs.sin2_theta_w)
    g_z = g_weak / cos_theta_w
    lepton_left = g_z * (-0.5 + inputs.sin2_theta_w)
    lepton_right = g_z * inputs.sin2_theta_w
    lepton_vector = lepton_left + lepton_right
    lepton_axial = lepton_right - lepton_left
    up_neutral_delta = g_z / 2.0
    return (
        float(up_neutral_delta),
        float(lepton_left),
        float(lepton_right),
        float(lepton_vector),
        float(lepton_axial),
    )


def _wilson_prefactor(
    *,
    lambda_b: complex,
    m_kk_gev: float,
    inputs: RareCharmDileptonSMInputs,
) -> complex:
    return complex(
        math.pi
        / (
            2.0
            * math.sqrt(2.0)
            * inputs.gf_gev_minus2
            * inputs.alpha_em_mz
            * lambda_b
            * m_kk_gev**2
        )
    )


def compute_rare_charm_dilepton_wilsons(
    source: QuarkMassBasisCouplings,
    *,
    transition: str = "c_u",
    m_kk_gev: float | None = None,
    inputs: RareCharmDileptonSMInputs | None = None,
) -> RareCharmDileptonWilsonCoefficients:
    """Match mass-basis couplings onto the v1 ``c -> u l l`` proxy."""
    if transition != "c_u":
        raise ValueError(f"unsupported rare-charm transition {transition!r}")
    p = default_sm_inputs() if inputs is None else inputs
    meson = p.d0
    resolved_m_kk = _positive_float(
        getattr(source, "M_KK") if m_kk_gev is None else m_kk_gev,
        "m_kk_gev",
    )
    u_idx = meson.light_up_index
    c_idx = meson.charm_index
    left_uc = _matrix_entry(source, "left_up", u_idx, c_idx)
    right_uc = _matrix_entry(source, "right_up", u_idx, c_idx)
    g_s = _positive_float(getattr(source, "g_s", 1.0), "g_s")
    left_overlap = left_uc / g_s
    right_overlap = right_uc / g_s

    factors = ckm_factors(transition, p)
    neutral_delta, lep_left, lep_right, lep_vector, lep_axial = _weak_couplings(p)
    left_quark_delta = neutral_delta * left_overlap
    right_quark_delta = neutral_delta * right_overlap
    prefactor = _wilson_prefactor(
        lambda_b=factors.lambda_b,
        m_kk_gev=resolved_m_kk,
        inputs=p,
    )
    c9_np = prefactor * left_quark_delta * lep_vector
    c10_np = prefactor * left_quark_delta * lep_axial
    c9p_np = prefactor * right_quark_delta * lep_vector
    c10p_np = prefactor * right_quark_delta * lep_axial

    return RareCharmDileptonWilsonCoefficients(
        model_label=RARE_CHARM_DILEPTON_MODEL_V1,
        operator_convention=RARE_CHARM_DILEPTON_OPERATOR_CONVENTION,
        matching_assumption=RARE_CHARM_DILEPTON_RS_MATCHING_ASSUMPTION_V1,
        transition_key=transition,
        M_KK=resolved_m_kk,
        matching_scale=resolved_m_kk,
        lambda_b=factors.lambda_b,
        left_uc_coupling=left_uc,
        right_uc_coupling=right_uc,
        left_uc_overlap=complex(left_overlap),
        right_uc_overlap=complex(right_overlap),
        left_quark_delta=complex(left_quark_delta),
        right_quark_delta=complex(right_quark_delta),
        lepton_left_delta=lep_left,
        lepton_right_delta=lep_right,
        lepton_vector_delta=lep_vector,
        lepton_axial_delta=lep_axial,
        c9_np=complex(c9_np),
        c10_np=complex(c10_np),
        c9p_np=complex(c9p_np),
        c10p_np=complex(c10p_np),
    )


def _tau_ps_to_gev_inverse(tau_ps: float, hbar_gev_s: float) -> float:
    return float(tau_ps * 1.0e-12 / hbar_gev_s)


def _branching_from_wilson_values(
    *,
    lepton: str,
    c10_np: complex,
    c10p_np: complex,
    cs_np: complex,
    csp_np: complex,
    cp_np: complex,
    cpp_np: complex,
    inputs: RareCharmDileptonSMInputs,
) -> tuple[float, complex, complex, complex, RareCharmDileptonCKMFactors]:
    meson = inputs.d0
    lep = inputs.lepton(lepton)
    factors = ckm_factors("c_u", inputs)
    m_c = meson.charm_quark_mass_gev
    m_u = meson.light_quark_mass_gev
    m_l = lep.mass_gev
    m_d = meson.meson_mass_gev
    beta2 = 1.0 - 4.0 * m_l * m_l / (m_d * m_d)
    if beta2 <= 0.0:
        raise ValueError(f"{lepton}: leptonic phase-space beta^2 <= 0")
    beta = math.sqrt(beta2)
    scalar_scale = m_d * m_d / (2.0 * m_l) * (m_c / (m_c + m_u))

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
        * m_d
        * m_l**2
        * abs(factors.lambda_b) ** 2
        * beta
    )
    branching = prefactor * (
        abs(pseudoscalar_amplitude) ** 2 + abs(scalar_amplitude) ** 2
    )
    return (
        float(branching),
        c10_total,
        scalar_amplitude,
        pseudoscalar_amplitude,
        factors,
    )


def sm_branching_fraction(
    lepton: str = "mu",
    inputs: RareCharmDileptonSMInputs | None = None,
) -> RareCharmLeptonicBranchingResult:
    """Evaluate the short-distance SM-limit ``D0 -> l+ l-`` branching fraction."""
    p = default_sm_inputs() if inputs is None else inputs
    br, c10_total, scalar_amp, pseudoscalar_amp, factors = _branching_from_wilson_values(
        lepton=lepton,
        c10_np=0.0j,
        c10p_np=0.0j,
        cs_np=0.0j,
        csp_np=0.0j,
        cp_np=0.0j,
        cpp_np=0.0j,
        inputs=p,
    )
    lep = p.lepton(lepton)
    return RareCharmLeptonicBranchingResult(
        model_label=RARE_CHARM_DILEPTON_MODEL_V1,
        input_bundle=p.input_bundle,
        transition_key="c_u",
        lepton_key=lep.lepton_key,
        branching_fraction=br,
        sm_branching_fraction=br,
        np_shift_branching_fraction=0.0,
        c10_total=c10_total,
        c10_sm=float(p.c10_sm),
        c10_leptonic_np=0.0j,
        c9_effective_np=0.0j,
        scalar_amplitude=scalar_amp,
        pseudoscalar_amplitude=pseudoscalar_amp,
        lambda_b=factors.lambda_b,
        diagnostics={
            "lambda_wolfenstein": factors.lambda_wolfenstein,
            "lambda_b": factors.lambda_b,
            "c10_sm": float(p.c10_sm),
            "sm_short_distance_policy": (
                "No nonzero short-distance SM C10 anchor is provided in the "
                "C004 catalog sidecar; the clean constrained piece is evaluated "
                "with C10_SM=0 and the long-distance SM context is carried by "
                "the constraint diagnostics."
            ),
            "constants_citation": p.constants_citation,
        },
    )


def evaluate_d0_to_ll(
    source: QuarkMassBasisCouplings | RareCharmDileptonWilsonCoefficients | None = None,
    *,
    lepton: str = "mu",
    m_kk_gev: float | None = None,
    inputs: RareCharmDileptonSMInputs | None = None,
) -> RareCharmLeptonicBranchingResult:
    """Evaluate the short-distance ``BR(D0 -> l+ l-)`` with the v1 RS proxy."""
    p = default_sm_inputs() if inputs is None else inputs
    wilsons: RareCharmDileptonWilsonCoefficients | None
    if source is None:
        wilsons = None
        c10_np = c10p_np = cs_np = csp_np = cp_np = cpp_np = 0.0j
        c9_effective = 0.0j
    elif isinstance(source, RareCharmDileptonWilsonCoefficients):
        wilsons = source
        c10_np = source.c10_np
        c10p_np = source.c10p_np
        cs_np = source.cs_np
        csp_np = source.csp_np
        cp_np = source.cp_np
        cpp_np = source.cpp_np
        c9_effective = source.c9_effective_np
    else:
        wilsons = compute_rare_charm_dilepton_wilsons(
            source,
            transition="c_u",
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

    br, c10_total, scalar_amp, pseudoscalar_amp, factors = _branching_from_wilson_values(
        lepton=lepton,
        c10_np=c10_np,
        c10p_np=c10p_np,
        cs_np=cs_np,
        csp_np=csp_np,
        cp_np=cp_np,
        cpp_np=cpp_np,
        inputs=p,
    )
    sm = sm_branching_fraction(lepton, p)
    diagnostics: dict[str, float | complex | str | bool] = {
        "lambda_wolfenstein": factors.lambda_wolfenstein,
        "lambda_b": factors.lambda_b,
        "c10_sm": float(p.c10_sm),
        "c10_total": complex(c10_total),
        "c10_leptonic_np": complex(c10_np - c10p_np),
        "c9_effective_np": complex(c9_effective),
        "c9_does_not_enter_leptonic_rate": True,
        "scalar_amplitude": complex(scalar_amp),
        "pseudoscalar_amplitude": complex(pseudoscalar_amp),
        "matching_assumption": RARE_CHARM_DILEPTON_RS_MATCHING_ASSUMPTION_V1,
        "charm_long_distance_dominated": True,
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

    lep = p.lepton(lepton)
    return RareCharmLeptonicBranchingResult(
        model_label=RARE_CHARM_DILEPTON_MODEL_V1,
        input_bundle=p.input_bundle,
        transition_key="c_u",
        lepton_key=lep.lepton_key,
        branching_fraction=br,
        sm_branching_fraction=float(sm.branching_fraction),
        np_shift_branching_fraction=float(br - sm.branching_fraction),
        c10_total=c10_total,
        c10_sm=float(p.c10_sm),
        c10_leptonic_np=complex(c10_np - c10p_np),
        c9_effective_np=complex(c9_effective),
        scalar_amplitude=scalar_amp,
        pseudoscalar_amplitude=pseudoscalar_amp,
        lambda_b=factors.lambda_b,
        wilsons=wilsons,
        diagnostics=diagnostics,
    )


__all__ = [
    "RARE_CHARM_DILEPTON_MODEL_V1",
    "RARE_CHARM_DILEPTON_OPERATOR_CONVENTION",
    "RARE_CHARM_DILEPTON_INPUT_BUNDLE_V1",
    "RARE_CHARM_DILEPTON_RS_MATCHING_ASSUMPTION_V1",
    "RareCharmLeptonInputs",
    "RareCharmMesonInputs",
    "RareCharmDileptonSMInputs",
    "RareCharmDileptonCKMFactors",
    "RareCharmDileptonWilsonCoefficients",
    "RareCharmLeptonicBranchingResult",
    "default_sm_inputs",
    "ckm_factors",
    "compute_rare_charm_dilepton_wilsons",
    "sm_branching_fraction",
    "evaluate_d0_to_ll",
]
