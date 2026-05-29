"""Short-distance ``s -> d nu nubar`` machinery for ``K+ -> pi+ nu nubar``.

Physics convention
------------------
The SM branching fraction is evaluated with the standard Buras/Brod-Gorbahn-
Stamou charged-mode parametrization

    BR(K+ -> pi+ nu nubar) = kappa_+ (1 + Delta_EM) [
        (Im(lambda_t X_eff) / lambda^5)^2
        + (Re(lambda_c) P_c / lambda + Re(lambda_t X_eff) / lambda^5)^2
    ],

where ``lambda_q = V_qs^* V_qd`` and the repo's fixed CKM target is used.
The default constants are ``kappa_+ = 5.173e-11 (lambda/0.225)^8``,
``P_c = 0.404``, ``X_t = 1.481``, and ``1 + Delta_EM = 0.997``.  These are
the common numerical inputs used in the Buras et al. status review
(JHEP 11 (2015) 033, arXiv:1503.02693).  The K004 catalog sidecar uses the
Buras-Venturini 2022 SM value as its main anchor (arXiv:2203.10099), while the
Brod-Gorbahn-Stamou 2021 value is kept as a validation reference
(arXiv:2105.02868).

RS matching assumption
----------------------
NEEDS-HUMAN-PHYSICS: the current ParameterPoint carries quark mass-basis
``KK-gluon``-style coupling matrices, but not the electroweak KK/Z/Z' tower or
neutrino-localization couplings needed for a model-complete RS prediction.
The v1 matching below is therefore a documented leading proxy: divide the
supplied quark coupling by its ``g_s`` to keep the flavor-overlap structure,
couple that overlap to one Z-like neutral boson with quark and neutrino
couplings ``g/(2 c_W)``, and match

    X_NP = Delta_sd Delta_nu / (g_SM^2 M_KK^2),
    g_SM^2 = 4 G_F^2 M_W^2 / (2 pi^2).

Both left- and right-handed down-quark currents enter the charged-mode
pseudoscalar matrix element with the same sign in this proxy.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Mapping

import numpy as np

from .couplings import QuarkMassBasisCouplings
from .model import RotationParameters, ckm_like_unitary

RARE_KAON_SND_MODEL_V1 = "rare_kaon_snd_buras_bgs_rs_proxy_v1"
RARE_KAON_OPERATOR_CONVENTION = (
    "O_LR=(sbar gamma_mu P_{L,R} d)(nubar gamma^mu P_L nu)"
)
RARE_KAON_INPUT_BUNDLE_V1 = "rare_kaon_sm_inputs_buras_bv_deltaem_repo_ckm_v1"
RARE_KAON_RS_MATCHING_ASSUMPTION_V1 = (
    "NEEDS-HUMAN-PHYSICS: quark KK-gluon mass-basis couplings are used as "
    "neutral-current flavor-overlap proxies; a single Z-like boson with "
    "Delta_q=Delta_nu=g/(2 c_W) stands in for the full RS EW KK/Z/Z' and "
    "neutrino-sector matching."
)


@dataclass(frozen=True)
class RareKaonSMInputs:
    """Numerical inputs for the charged rare-kaon SM parametrization."""

    input_bundle: str = RARE_KAON_INPUT_BUNDLE_V1
    kappa_plus_ref: float = 5.173e-11
    kappa_lambda_ref: float = 0.225
    p_c: float = 0.404
    x_t: float = 1.481
    delta_em_correction: float = 0.997
    theta12: float = 0.2274
    theta13: float = 0.00368
    theta23: float = 0.0415
    delta: float = 1.196
    gf_gev_minus2: float = 1.1663787e-5
    m_w_gev: float = 80.379
    alpha_em_mz: float = 1.0 / 127.952
    sin2_theta_w: float = 0.23122
    constants_citation: str = (
        "Buras et al. JHEP 11 (2015) 033, arXiv:1503.02693; "
        "Buras-Venturini arXiv:2203.10099; "
        "Brod-Gorbahn-Stamou validation reference arXiv:2105.02868; "
        "repo CKM target quarkConstraints.modern.inputs.ModernDefaultCKMTarget"
    )

    def __post_init__(self) -> None:
        for name in (
            "kappa_plus_ref",
            "kappa_lambda_ref",
            "p_c",
            "x_t",
            "delta_em_correction",
            "gf_gev_minus2",
            "m_w_gev",
            "alpha_em_mz",
            "sin2_theta_w",
        ):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be positive and finite")
        if not 0.0 < float(self.sin2_theta_w) < 1.0:
            raise ValueError("sin2_theta_w must lie between zero and one")


@dataclass(frozen=True)
class RareKaonCKMFactors:
    """CKM factors entering the ``K+ -> pi+ nu nubar`` formula."""

    lambda_wolfenstein: float
    lambda_c: complex
    lambda_t: complex
    matrix: tuple[tuple[complex, ...], ...]


@dataclass(frozen=True)
class RareKaonWilsonCoefficients:
    """Leading ``s -> d nu nubar`` Wilson proxy from mass-basis couplings."""

    model_label: str
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
class RareKaonBranchingResult:
    """Branching-ratio prediction for ``K+ -> pi+ nu nubar``."""

    model_label: str
    input_bundle: str
    branching_fraction: float
    sm_branching_fraction: float
    np_shift_branching_fraction: float
    kappa_plus: float
    p_c: float
    x_t: float
    delta_em_correction: float
    x_eff_top: complex
    x_np_total: complex
    lambda_wolfenstein: float
    lambda_c: complex
    lambda_t: complex
    wilsons: RareKaonWilsonCoefficients | None = None
    diagnostics: Mapping[str, float | complex | str] = field(default_factory=dict)


def default_sm_inputs() -> RareKaonSMInputs:
    """Return the repo-owned default rare-kaon SM input bundle."""
    return RareKaonSMInputs()


def ckm_factors(inputs: RareKaonSMInputs | None = None) -> RareKaonCKMFactors:
    """Return ``lambda``, ``lambda_c``, and ``lambda_t`` from the repo CKM target."""
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
    lambda_c = complex(np.conjugate(matrix_np[1, 1]) * matrix_np[1, 0])
    lambda_t = complex(np.conjugate(matrix_np[2, 1]) * matrix_np[2, 0])
    matrix = tuple(tuple(complex(entry) for entry in row) for row in matrix_np)
    return RareKaonCKMFactors(
        lambda_wolfenstein=lam,
        lambda_c=lambda_c,
        lambda_t=lambda_t,
        matrix=matrix,
    )


def kappa_plus(inputs: RareKaonSMInputs | None = None) -> float:
    """Return ``kappa_+`` with the standard ``(lambda/0.225)^8`` rescaling."""
    p = default_sm_inputs() if inputs is None else inputs
    lam = ckm_factors(p).lambda_wolfenstein
    return float(p.kappa_plus_ref * (lam / p.kappa_lambda_ref) ** 8)


def _branching_fraction_from_x_np(
    x_np_total: complex,
    *,
    inputs: RareKaonSMInputs,
) -> tuple[float, float, RareKaonCKMFactors, complex]:
    factors = ckm_factors(inputs)
    lam = factors.lambda_wolfenstein
    x_eff_top = factors.lambda_t * inputs.x_t + complex(x_np_total)
    imag_term = (x_eff_top.imag / lam**5) ** 2
    real_term = (
        factors.lambda_c.real / lam * inputs.p_c
        + x_eff_top.real / lam**5
    ) ** 2
    kappa = kappa_plus(inputs)
    br = inputs.delta_em_correction * kappa * (imag_term + real_term)
    return float(br), float(kappa), factors, x_eff_top


def sm_branching_fraction(inputs: RareKaonSMInputs | None = None) -> RareKaonBranchingResult:
    """Evaluate the SM-limit charged rare-kaon branching fraction."""
    p = default_sm_inputs() if inputs is None else inputs
    br, kappa, factors, x_eff_top = _branching_fraction_from_x_np(0.0j, inputs=p)
    return RareKaonBranchingResult(
        model_label=RARE_KAON_SND_MODEL_V1,
        input_bundle=p.input_bundle,
        branching_fraction=br,
        sm_branching_fraction=br,
        np_shift_branching_fraction=0.0,
        kappa_plus=kappa,
        p_c=float(p.p_c),
        x_t=float(p.x_t),
        delta_em_correction=float(p.delta_em_correction),
        x_eff_top=x_eff_top,
        x_np_total=0.0j,
        lambda_wolfenstein=factors.lambda_wolfenstein,
        lambda_c=factors.lambda_c,
        lambda_t=factors.lambda_t,
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


def _weak_neutral_current_coupling(inputs: RareKaonSMInputs) -> float:
    g_weak = math.sqrt(4.0 * math.pi * inputs.alpha_em_mz / inputs.sin2_theta_w)
    cos_theta_w = math.sqrt(1.0 - inputs.sin2_theta_w)
    return float(g_weak / (2.0 * cos_theta_w))


def g_sm_squared(inputs: RareKaonSMInputs | None = None) -> float:
    """Return the Buras ``g_SM^2`` normalization in GeV^-2."""
    p = default_sm_inputs() if inputs is None else inputs
    return float(4.0 * p.gf_gev_minus2**2 * p.m_w_gev**2 / (2.0 * math.pi**2))


def compute_rare_kaon_wilsons(
    source: QuarkMassBasisCouplings,
    *,
    m_kk_gev: float | None = None,
    inputs: RareKaonSMInputs | None = None,
) -> RareKaonWilsonCoefficients:
    """Match mass-basis couplings onto the v1 ``s -> d nu nubar`` proxy.

    The flavor index convention follows the existing kaon Delta-F=2 path:
    the ``(0, 1)`` down-sector entry is reported as the ``s-d`` coupling.
    """
    p = default_sm_inputs() if inputs is None else inputs
    resolved_m_kk = _positive_float(
        getattr(source, "M_KK") if m_kk_gev is None else m_kk_gev,
        "m_kk_gev",
    )
    left_sd = _matrix_entry(source, "left_down", 0, 1)
    right_sd = _matrix_entry(source, "right_down", 0, 1)
    g_s = _positive_float(getattr(source, "g_s", 1.0), "g_s")
    left_overlap = left_sd / g_s
    right_overlap = right_sd / g_s

    neutral_delta = _weak_neutral_current_coupling(p)
    left_quark_delta = neutral_delta * left_overlap
    right_quark_delta = neutral_delta * right_overlap
    neutrino_delta = neutral_delta
    normalizer = g_sm_squared(p) * resolved_m_kk**2
    x_np_left = left_quark_delta * neutrino_delta / normalizer
    x_np_right = right_quark_delta * neutrino_delta / normalizer

    return RareKaonWilsonCoefficients(
        model_label=RARE_KAON_SND_MODEL_V1,
        operator_convention=RARE_KAON_OPERATOR_CONVENTION,
        matching_assumption=RARE_KAON_RS_MATCHING_ASSUMPTION_V1,
        M_KK=resolved_m_kk,
        matching_scale=resolved_m_kk,
        left_sd_coupling=left_sd,
        right_sd_coupling=right_sd,
        left_sd_overlap=complex(left_overlap),
        right_sd_overlap=complex(right_overlap),
        left_quark_delta=complex(left_quark_delta),
        right_quark_delta=complex(right_quark_delta),
        neutrino_delta=neutrino_delta,
        x_np_left=complex(x_np_left),
        x_np_right=complex(x_np_right),
    )


def evaluate_kplus_to_piplus_nunu(
    source: QuarkMassBasisCouplings | RareKaonWilsonCoefficients | None = None,
    *,
    m_kk_gev: float | None = None,
    inputs: RareKaonSMInputs | None = None,
) -> RareKaonBranchingResult:
    """Evaluate ``BR(K+ -> pi+ nu nubar)`` in the SM or with the v1 RS proxy."""
    p = default_sm_inputs() if inputs is None else inputs
    wilsons: RareKaonWilsonCoefficients | None
    if source is None:
        wilsons = None
        x_np_total = 0.0j
    elif isinstance(source, RareKaonWilsonCoefficients):
        wilsons = source
        x_np_total = source.x_np_total
    else:
        wilsons = compute_rare_kaon_wilsons(source, m_kk_gev=m_kk_gev, inputs=p)
        x_np_total = wilsons.x_np_total

    br, kappa, factors, x_eff_top = _branching_fraction_from_x_np(
        x_np_total,
        inputs=p,
    )
    sm = sm_branching_fraction(p).branching_fraction
    diagnostics: dict[str, float | complex | str] = {
        "g_sm_squared_gev_minus2": g_sm_squared(p),
        "delta_em_correction": float(p.delta_em_correction),
        "matching_assumption": RARE_KAON_RS_MATCHING_ASSUMPTION_V1,
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
                "neutrino_delta": float(wilsons.neutrino_delta),
                "x_np_left": complex(wilsons.x_np_left),
                "x_np_right": complex(wilsons.x_np_right),
                "x_np_total": complex(wilsons.x_np_total),
            }
        )

    return RareKaonBranchingResult(
        model_label=RARE_KAON_SND_MODEL_V1,
        input_bundle=p.input_bundle,
        branching_fraction=br,
        sm_branching_fraction=sm,
        np_shift_branching_fraction=float(br - sm),
        kappa_plus=kappa,
        p_c=float(p.p_c),
        x_t=float(p.x_t),
        delta_em_correction=float(p.delta_em_correction),
        x_eff_top=x_eff_top,
        x_np_total=complex(x_np_total),
        lambda_wolfenstein=factors.lambda_wolfenstein,
        lambda_c=factors.lambda_c,
        lambda_t=factors.lambda_t,
        wilsons=wilsons,
        diagnostics=diagnostics,
    )


__all__ = [
    "RARE_KAON_SND_MODEL_V1",
    "RARE_KAON_OPERATOR_CONVENTION",
    "RARE_KAON_INPUT_BUNDLE_V1",
    "RARE_KAON_RS_MATCHING_ASSUMPTION_V1",
    "RareKaonSMInputs",
    "RareKaonCKMFactors",
    "RareKaonWilsonCoefficients",
    "RareKaonBranchingResult",
    "default_sm_inputs",
    "ckm_factors",
    "kappa_plus",
    "g_sm_squared",
    "sm_branching_fraction",
    "compute_rare_kaon_wilsons",
    "evaluate_kplus_to_piplus_nunu",
]
