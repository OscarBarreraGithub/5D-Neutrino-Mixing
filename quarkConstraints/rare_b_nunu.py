"""Short-distance ``b -> s nu nubar`` machinery.

Physics convention
------------------
The Standard Model short-distance coefficient is written in the common
``b -> s nu nubar`` convention

    C_L^SM = - X_t / sin^2(theta_W),    C_R^SM = 0,

with ``X_t`` evaluated from the Inami-Lim top function ``X_0(x_t)`` and a
small QCD/electroweak multiplier.  The exclusive ``B+ -> K+ nu nubar``
normalization is supplied as an input branching fraction.  For the charged
``B+ -> K+`` mode, the HPQCD long-distance term is not modified by
short-distance new physics; only the short-distance remainder is rescaled by

    R_K = |C_L + C_R|^2 / |C_L^SM|^2.

The charged prediction is therefore

    BR = BR_LD + (BR_SM_total - BR_LD) R_K.

This separation is intentional: form-factor integration lives in the
catalogued SM normalization, the charged-B long-distance piece remains fixed,
and this module owns the reusable ``b -> s nu nubar`` short-distance response.
The same ``C_L, C_R, epsilon, eta`` response is designed for reuse by
``B -> K* nu nubar`` constraints.

RS matching assumption
----------------------
NEEDS-HUMAN-PHYSICS: the current ParameterPoint carries quark mass-basis
``KK-gluon``-style coupling matrices, but not the electroweak KK/Z/Z' tower or
lepton/neutrino localization couplings needed for a model-complete RS
prediction.  The v1 matching below is therefore a documented leading proxy:
divide the supplied quark coupling by its ``g_s`` to keep the flavor-overlap
structure, couple that overlap to one Z-like neutral boson with quark and
neutrino couplings ``g/(2 c_W)``, and match

    X_NP = Delta_sb Delta_nu / (g_SM^2 M_KK^2),
    g_SM^2 = 4 G_F^2 M_W^2 / (2 pi^2).

Both left- and right-handed down-quark currents are kept.  For the
pseudoscalar ``B -> K`` mode they enter through ``C_L + C_R``.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
import math
from typing import Mapping

import numpy as np

from qcd.constants import M_TOP_MS

from .couplings import QuarkMassBasisCouplings
from .model import RotationParameters, ckm_like_unitary

RARE_B_NUNU_MODEL_V1 = "rare_b_nunu_xt_rs_proxy_v1"
RARE_B_NUNU_OPERATOR_CONVENTION = (
    "O_LR=(sbar gamma_mu P_{L,R} b)(nubar gamma^mu P_L nu)"
)
RARE_B_NUNU_INPUT_BUNDLE_V1 = "rare_b_nunu_sm_inputs_xt_hpqcd_repo_ckm_v1"
RARE_B_NUNU_RS_MATCHING_ASSUMPTION_V1 = (
    "NEEDS-HUMAN-PHYSICS: quark KK-gluon mass-basis couplings are used as "
    "neutral-current flavor-overlap proxies; a single Z-like boson with "
    "Delta_q=Delta_nu=g/(2 c_W) stands in for the full RS EW KK/Z/Z' and "
    "neutrino-sector matching."
)
RARE_B_NUNU_BPLUS_KPLUS_SM_CITATION = (
    "HPQCD 2023 B -> K nu nubar SM normalization, arXiv:2207.13371; "
    "catalog constraints should override this from their YAML sidecar"
)
RARE_B_NUNU_BPLUS_KPLUS_LONG_DISTANCE_CITATION = (
    "HPQCD 2023 charged-mode long-distance contribution "
    "BR_LD(B+ -> K+ nu nubar)=6.09(53)e-7, arXiv:2207.13371"
)
RARE_B_NUNU_KSTAR_ETA_COEFFICIENT = 1.31


@dataclass(frozen=True)
class RareBNuNuSMInputs:
    """Numerical inputs for the reusable ``b -> s nu nubar`` response."""

    input_bundle: str = RARE_B_NUNU_INPUT_BUNDLE_V1
    br_bplus_kplus_sm: float = 5.58e-6
    br_bplus_kplus_long_distance: float = 6.09e-7
    m_t_msbar_gev: float = M_TOP_MS
    m_w_gev: float = 80.379
    eta_x: float = 0.994
    sin2_theta_w: float = 0.23122
    alpha_em_mz: float = 1.0 / 127.952
    gf_gev_minus2: float = 1.1663787e-5
    theta12: float = 0.2274
    theta13: float = 0.00368
    theta23: float = 0.0415
    delta: float = 1.196
    constants_citation: str = (
        "Inami-Lim top function with eta_X=0.994 and m_t(m_t)=162.5 GeV; "
        "electroweak constants follow the rare-kaon input bundle in this repo; "
        "repo CKM target quarkConstraints.modern.inputs.ModernDefaultCKMTarget"
    )
    sm_normalization_citation: str = RARE_B_NUNU_BPLUS_KPLUS_SM_CITATION
    long_distance_citation: str = RARE_B_NUNU_BPLUS_KPLUS_LONG_DISTANCE_CITATION

    def __post_init__(self) -> None:
        for name in (
            "br_bplus_kplus_sm",
            "br_bplus_kplus_long_distance",
            "m_t_msbar_gev",
            "m_w_gev",
            "eta_x",
            "sin2_theta_w",
            "alpha_em_mz",
            "gf_gev_minus2",
        ):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be positive and finite")
        if float(self.br_bplus_kplus_long_distance) >= float(self.br_bplus_kplus_sm):
            raise ValueError(
                "br_bplus_kplus_long_distance must be smaller than "
                "br_bplus_kplus_sm"
            )
        if not 0.0 < float(self.sin2_theta_w) < 1.0:
            raise ValueError("sin2_theta_w must lie between zero and one")


@dataclass(frozen=True)
class RareBNuNuCKMFactors:
    """CKM factors entering the ``b -> s nu nubar`` Hamiltonian."""

    lambda_wolfenstein: float
    lambda_t_bs: complex
    matrix: tuple[tuple[complex, ...], ...]


@dataclass(frozen=True)
class RareBNuNuWilsonCoefficients:
    """Leading ``b -> s nu nubar`` Wilson proxy from mass-basis couplings."""

    model_label: str
    operator_convention: str
    matching_assumption: str
    M_KK: float
    matching_scale: float
    left_sb_coupling: complex
    right_sb_coupling: complex
    left_sb_overlap: complex
    right_sb_overlap: complex
    left_quark_delta: complex
    right_quark_delta: complex
    neutrino_delta: float
    x_np_left: complex
    x_np_right: complex

    @property
    def x_np_total(self) -> complex:
        return complex(self.x_np_left + self.x_np_right)

    @property
    def wilsons(self) -> Mapping[str, complex]:
        return {
            "X_NP_L": complex(self.x_np_left),
            "X_NP_R": complex(self.x_np_right),
            "X_NP_total": self.x_np_total,
        }


@dataclass(frozen=True)
class RareBNuNuShortDistanceResult:
    """Reusable short-distance response for ``b -> s nu nubar`` modes."""

    model_label: str
    input_bundle: str
    x_t: float
    lambda_t_bs: complex
    c_l_sm: complex
    c_l_total: complex
    c_r_total: complex
    x_eff_left: complex
    x_eff_right: complex
    epsilon: float
    eta: float
    r_k: float
    r_kstar: float


@dataclass(frozen=True)
class RareBNuNuBranchingResult:
    """Branching-ratio prediction for ``B+ -> K+ nu nubar``."""

    model_label: str
    input_bundle: str
    branching_fraction: float
    sm_branching_fraction: float
    np_shift_branching_fraction: float
    ratio_to_sm: float
    x_t: float
    lambda_wolfenstein: float
    lambda_t_bs: complex
    c_l_sm: complex
    c_l_total: complex
    c_r_total: complex
    x_eff_left: complex
    x_eff_right: complex
    epsilon: float
    eta: float
    r_k: float
    r_kstar: float
    wilsons: RareBNuNuWilsonCoefficients | None = None
    diagnostics: Mapping[str, float | complex | str] = field(default_factory=dict)


def default_sm_inputs() -> RareBNuNuSMInputs:
    """Return the repo-owned default rare-B ``b -> s nu nubar`` input bundle."""
    return RareBNuNuSMInputs()


def sm_inputs_with_bplus_kplus_normalization(
    branching_fraction: float,
    *,
    inputs: RareBNuNuSMInputs | None = None,
) -> RareBNuNuSMInputs:
    """Return inputs with the ``B+ -> K+`` SM normalization replaced."""
    number = _positive_float(branching_fraction, "branching_fraction")
    p = default_sm_inputs() if inputs is None else inputs
    return replace(p, br_bplus_kplus_sm=number)


def inami_lim_x0(x_t: float) -> float:
    """Return the leading Inami-Lim ``X_0(x_t)`` function."""
    x = _positive_float(x_t, "x_t")
    if math.isclose(x, 1.0, rel_tol=0.0, abs_tol=1.0e-12):
        raise ValueError("x_t must not be one")
    return float(
        x
        / 8.0
        * (
            (x + 2.0) / (x - 1.0)
            + (3.0 * x - 6.0) / ((x - 1.0) ** 2) * math.log(x)
        )
    )


def x_t_top_function(inputs: RareBNuNuSMInputs | None = None) -> float:
    """Return ``X_t`` from ``X_0((m_t/M_W)^2)`` times ``eta_X``."""
    p = default_sm_inputs() if inputs is None else inputs
    x = (float(p.m_t_msbar_gev) / float(p.m_w_gev)) ** 2
    return float(float(p.eta_x) * inami_lim_x0(x))


def ckm_factors(inputs: RareBNuNuSMInputs | None = None) -> RareBNuNuCKMFactors:
    """Return ``lambda_t = V_ts^* V_tb`` from the repo CKM target."""
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
    lambda_t_bs = complex(np.conjugate(matrix_np[2, 1]) * matrix_np[2, 2])
    matrix = tuple(tuple(complex(entry) for entry in row) for row in matrix_np)
    return RareBNuNuCKMFactors(
        lambda_wolfenstein=lam,
        lambda_t_bs=lambda_t_bs,
        matrix=matrix,
    )


def g_sm_squared(inputs: RareBNuNuSMInputs | None = None) -> float:
    """Return the Buras ``g_SM^2`` normalization in GeV^-2."""
    p = default_sm_inputs() if inputs is None else inputs
    return float(4.0 * p.gf_gev_minus2**2 * p.m_w_gev**2 / (2.0 * math.pi**2))


def short_distance_response(
    x_np_left: complex = 0.0j,
    x_np_right: complex = 0.0j,
    *,
    inputs: RareBNuNuSMInputs | None = None,
) -> RareBNuNuShortDistanceResult:
    """Return reusable ``C_L, C_R, epsilon, eta`` response parameters."""
    p = default_sm_inputs() if inputs is None else inputs
    factors = ckm_factors(p)
    x_t = x_t_top_function(p)
    if abs(factors.lambda_t_bs) <= 0.0:
        raise ValueError("lambda_t_bs must be nonzero")

    x_eff_left = factors.lambda_t_bs * x_t + complex(x_np_left)
    x_eff_right = complex(x_np_right)
    c_l_sm = complex(-x_t / p.sin2_theta_w)
    c_l_total = -x_eff_left / (factors.lambda_t_bs * p.sin2_theta_w)
    c_r_total = -x_eff_right / (factors.lambda_t_bs * p.sin2_theta_w)

    denom = abs(c_l_total) ** 2 + abs(c_r_total) ** 2
    if denom <= 0.0:
        epsilon = 0.0
        eta = 0.0
    else:
        epsilon = math.sqrt(denom / (abs(c_l_sm) ** 2))
        eta = -float((c_l_total * c_r_total.conjugate()).real) / denom

    r_k = abs(c_l_total + c_r_total) ** 2 / (abs(c_l_sm) ** 2)
    r_kstar = (1.0 + RARE_B_NUNU_KSTAR_ETA_COEFFICIENT * eta) * epsilon**2

    return RareBNuNuShortDistanceResult(
        model_label=RARE_B_NUNU_MODEL_V1,
        input_bundle=p.input_bundle,
        x_t=float(x_t),
        lambda_t_bs=factors.lambda_t_bs,
        c_l_sm=complex(c_l_sm),
        c_l_total=complex(c_l_total),
        c_r_total=complex(c_r_total),
        x_eff_left=complex(x_eff_left),
        x_eff_right=complex(x_eff_right),
        epsilon=float(epsilon),
        eta=float(eta),
        r_k=float(r_k),
        r_kstar=float(r_kstar),
    )


def sm_branching_fraction(
    inputs: RareBNuNuSMInputs | None = None,
) -> RareBNuNuBranchingResult:
    """Evaluate the SM-limit ``B+ -> K+ nu nubar`` branching fraction."""
    return evaluate_bplus_to_kplus_nunu(None, inputs=inputs)


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


def _weak_neutral_current_coupling(inputs: RareBNuNuSMInputs) -> float:
    g_weak = math.sqrt(4.0 * math.pi * inputs.alpha_em_mz / inputs.sin2_theta_w)
    cos_theta_w = math.sqrt(1.0 - inputs.sin2_theta_w)
    return float(g_weak / (2.0 * cos_theta_w))


def compute_rare_b_nunu_wilsons(
    source: QuarkMassBasisCouplings,
    *,
    m_kk_gev: float | None = None,
    inputs: RareBNuNuSMInputs | None = None,
) -> RareBNuNuWilsonCoefficients:
    """Match mass-basis couplings onto the v1 ``b -> s nu nubar`` proxy.

    The flavor index convention follows the existing B_s Delta-F=2 path:
    the ``(1, 2)`` down-sector entry is reported as the ``s-b`` coupling.
    """
    p = default_sm_inputs() if inputs is None else inputs
    resolved_m_kk = _positive_float(
        getattr(source, "M_KK") if m_kk_gev is None else m_kk_gev,
        "m_kk_gev",
    )
    left_sb = _matrix_entry(source, "left_down", 1, 2)
    right_sb = _matrix_entry(source, "right_down", 1, 2)
    g_s = _positive_float(getattr(source, "g_s", 1.0), "g_s")
    left_overlap = left_sb / g_s
    right_overlap = right_sb / g_s

    neutral_delta = _weak_neutral_current_coupling(p)
    left_quark_delta = neutral_delta * left_overlap
    right_quark_delta = neutral_delta * right_overlap
    neutrino_delta = neutral_delta
    normalizer = g_sm_squared(p) * resolved_m_kk**2
    x_np_left = left_quark_delta * neutrino_delta / normalizer
    x_np_right = right_quark_delta * neutrino_delta / normalizer

    return RareBNuNuWilsonCoefficients(
        model_label=RARE_B_NUNU_MODEL_V1,
        operator_convention=RARE_B_NUNU_OPERATOR_CONVENTION,
        matching_assumption=RARE_B_NUNU_RS_MATCHING_ASSUMPTION_V1,
        M_KK=resolved_m_kk,
        matching_scale=resolved_m_kk,
        left_sb_coupling=left_sb,
        right_sb_coupling=right_sb,
        left_sb_overlap=complex(left_overlap),
        right_sb_overlap=complex(right_overlap),
        left_quark_delta=complex(left_quark_delta),
        right_quark_delta=complex(right_quark_delta),
        neutrino_delta=float(neutrino_delta),
        x_np_left=complex(x_np_left),
        x_np_right=complex(x_np_right),
    )


def evaluate_bplus_to_kplus_nunu(
    source: QuarkMassBasisCouplings | RareBNuNuWilsonCoefficients | None = None,
    *,
    m_kk_gev: float | None = None,
    inputs: RareBNuNuSMInputs | None = None,
) -> RareBNuNuBranchingResult:
    """Evaluate ``BR(B+ -> K+ nu nubar)`` in the SM or with the v1 RS proxy."""
    p = default_sm_inputs() if inputs is None else inputs
    wilsons: RareBNuNuWilsonCoefficients | None
    if source is None:
        wilsons = None
        x_np_left = 0.0j
        x_np_right = 0.0j
    elif isinstance(source, RareBNuNuWilsonCoefficients):
        wilsons = source
        x_np_left = source.x_np_left
        x_np_right = source.x_np_right
    else:
        wilsons = compute_rare_b_nunu_wilsons(source, m_kk_gev=m_kk_gev, inputs=p)
        x_np_left = wilsons.x_np_left
        x_np_right = wilsons.x_np_right

    response = short_distance_response(x_np_left, x_np_right, inputs=p)
    sm = float(p.br_bplus_kplus_sm)
    br_long_distance = float(p.br_bplus_kplus_long_distance)
    br_short_distance_sm = float(sm - br_long_distance)
    br = float(br_long_distance + br_short_distance_sm * response.r_k)
    factors = ckm_factors(p)
    diagnostics: dict[str, float | complex | str] = {
        "g_sm_squared_gev_minus2": g_sm_squared(p),
        "matching_assumption": RARE_B_NUNU_RS_MATCHING_ASSUMPTION_V1,
        "operator_convention": RARE_B_NUNU_OPERATOR_CONVENTION,
        "sm_normalization_citation": p.sm_normalization_citation,
        "long_distance_citation": p.long_distance_citation,
        "constants_citation": p.constants_citation,
        "bplus_kplus_total_sm_branching_fraction": sm,
        "bplus_kplus_long_distance_branching_fraction": br_long_distance,
        "bplus_kplus_short_distance_sm_branching_fraction": br_short_distance_sm,
        "m_t_msbar_gev": float(p.m_t_msbar_gev),
        "m_w_gev": float(p.m_w_gev),
        "eta_x": float(p.eta_x),
        "sin2_theta_w": float(p.sin2_theta_w),
        "x_t": float(response.x_t),
        "c_l_sm": complex(response.c_l_sm),
        "c_l_total": complex(response.c_l_total),
        "c_r_total": complex(response.c_r_total),
        "x_eff_left": complex(response.x_eff_left),
        "x_eff_right": complex(response.x_eff_right),
        "epsilon": float(response.epsilon),
        "eta": float(response.eta),
        "r_k": float(response.r_k),
        "r_kstar": float(response.r_kstar),
    }
    if wilsons is not None:
        diagnostics.update(
            {
                "m_kk_gev": float(wilsons.M_KK),
                "matching_scale_gev": float(wilsons.matching_scale),
                "left_sb_coupling": complex(wilsons.left_sb_coupling),
                "right_sb_coupling": complex(wilsons.right_sb_coupling),
                "left_sb_overlap": complex(wilsons.left_sb_overlap),
                "right_sb_overlap": complex(wilsons.right_sb_overlap),
                "left_quark_delta": complex(wilsons.left_quark_delta),
                "right_quark_delta": complex(wilsons.right_quark_delta),
                "neutrino_delta": float(wilsons.neutrino_delta),
                "x_np_left": complex(wilsons.x_np_left),
                "x_np_right": complex(wilsons.x_np_right),
                "x_np_total": complex(wilsons.x_np_total),
            }
        )

    return RareBNuNuBranchingResult(
        model_label=RARE_B_NUNU_MODEL_V1,
        input_bundle=p.input_bundle,
        branching_fraction=br,
        sm_branching_fraction=sm,
        np_shift_branching_fraction=float(br - sm),
        ratio_to_sm=float(br / sm),
        x_t=float(response.x_t),
        lambda_wolfenstein=float(factors.lambda_wolfenstein),
        lambda_t_bs=complex(response.lambda_t_bs),
        c_l_sm=complex(response.c_l_sm),
        c_l_total=complex(response.c_l_total),
        c_r_total=complex(response.c_r_total),
        x_eff_left=complex(response.x_eff_left),
        x_eff_right=complex(response.x_eff_right),
        epsilon=float(response.epsilon),
        eta=float(response.eta),
        r_k=float(response.r_k),
        r_kstar=float(response.r_kstar),
        wilsons=wilsons,
        diagnostics=diagnostics,
    )


__all__ = [
    "RARE_B_NUNU_MODEL_V1",
    "RARE_B_NUNU_OPERATOR_CONVENTION",
    "RARE_B_NUNU_INPUT_BUNDLE_V1",
    "RARE_B_NUNU_RS_MATCHING_ASSUMPTION_V1",
    "RARE_B_NUNU_BPLUS_KPLUS_SM_CITATION",
    "RARE_B_NUNU_BPLUS_KPLUS_LONG_DISTANCE_CITATION",
    "RARE_B_NUNU_KSTAR_ETA_COEFFICIENT",
    "RareBNuNuSMInputs",
    "RareBNuNuCKMFactors",
    "RareBNuNuWilsonCoefficients",
    "RareBNuNuShortDistanceResult",
    "RareBNuNuBranchingResult",
    "default_sm_inputs",
    "sm_inputs_with_bplus_kplus_normalization",
    "inami_lim_x0",
    "x_t_top_function",
    "ckm_factors",
    "g_sm_squared",
    "short_distance_response",
    "sm_branching_fraction",
    "compute_rare_b_nunu_wilsons",
    "evaluate_bplus_to_kplus_nunu",
]
