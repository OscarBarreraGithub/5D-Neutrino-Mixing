"""Short-distance ``s -> d l+l-`` machinery for rare kaon decays.

Physics convention
------------------
This module starts with the reusable short-distance piece needed for
``K_L -> mu+ mu-`` and is intentionally separate from
``quarkConstraints.rare_kaon_snd`` (the ``K -> pi nu nubar`` machinery).
The ``K_L -> mu+ mu-`` observable is long-distance dominated; only the
short-distance dispersive component is evaluated here.  Following the
Buras/Isidori convention,

    BR(K_L -> mu+ mu-)_SD = kappa_mu [
        Re(Y_eff) / lambda^5 + Re(lambda_c) P_c(Y) / lambda
    ]^2,

with ``Y_eff = lambda_t (Y_L - Y_R)``.  The SM default uses
``kappa_mu = 2.01e-9 (lambda/0.2252)^8``, ``P_c(Y)=0.115``, and
``Y_L^SM ~= 0.94``, giving the expected ``O(0.8e-9)`` short-distance
branching fraction with the repo CKM target.  The normalization and
right-handed sign follow Buras, Buttazzo, Knegjens, JHEP 11 (2015) 166,
arXiv:1507.08672.  The YAML SM theory anchor uses the Gorbahn-Haisch NNLO
charm result, Phys. Rev. Lett. 97 (2006) 122002, arXiv:hep-ph/0605203,
``BR(K_L -> mu+ mu-)_SD = (0.79 +/- 0.12)e-9``; the conservative
short-distance-bound interpretation follows Isidori and Unterdorfer,
JHEP 01 (2004) 009, arXiv:hep-ph/0311084.

RS matching assumption
----------------------
NEEDS-HUMAN-PHYSICS: the current ParameterPoint carries quark mass-basis
``KK-gluon``-style coupling matrices, but not the electroweak KK/Z/Z' tower or
muon axial couplings needed for a model-complete RS prediction.  The v1
matching below is therefore a documented leading proxy: divide the supplied
quark coupling by its ``g_s`` to keep the flavor-overlap structure, couple that
overlap to one Z-like neutral boson with quark and muon axial couplings
``g/(2 c_W)``, and match

    Y_NP = Delta_A^mumu (Delta_L^sd - Delta_R^sd)
           / (g_SM^2 M_KK^2),
    g_SM^2 = 4 G_F^2 M_W^2 / (2 pi^2).

The relative minus sign between left- and right-handed quark currents is the
standard axial-lepton ``Y_L - Y_R`` convention for ``K_L -> mu+ mu-``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Mapping

import numpy as np

from .couplings import QuarkMassBasisCouplings
from .model import RotationParameters, ckm_like_unitary

RARE_KAON_DILEPTON_MODEL_V1 = "rare_kaon_dilepton_buras_isidori_rs_proxy_v1"
RARE_KAON_DILEPTON_OPERATOR_CONVENTION = (
    "Y_eff=lambda_t(Y_L-Y_R), O_A=(sbar gamma_mu P_{L,R} d)"
    "(mubar gamma^mu gamma5 mu)"
)
RARE_KAON_DILEPTON_INPUT_BUNDLE_V1 = (
    "rare_kaon_dilepton_sm_inputs_buras_isidori_repo_ckm_v1"
)
RARE_KAON_DILEPTON_RS_MATCHING_ASSUMPTION_V1 = (
    "NEEDS-HUMAN-PHYSICS: quark KK-gluon mass-basis couplings are used as "
    "neutral-current flavor-overlap proxies; a single Z-like boson with "
    "Delta_q=Delta_A^mumu=g/(2 c_W) stands in for the full RS EW KK/Z/Z' "
    "and lepton-sector matching."
)
RARE_KAON_DILEPTON_PARAMETRIZATION_CITATION = (
    "Buras, Buttazzo, Knegjens JHEP 11 (2015) 166, arXiv:1507.08672; "
    "Gorbahn-Haisch Phys. Rev. Lett. 97 (2006) 122002, "
    "arXiv:hep-ph/0605203; "
    "Isidori-Unterdorfer JHEP 01 (2004) 009, arXiv:hep-ph/0311084"
)


@dataclass(frozen=True)
class RareKaonDileptonSMInputs:
    """Numerical inputs for the ``K_L -> mu+ mu-`` SD parametrization."""

    input_bundle: str = RARE_KAON_DILEPTON_INPUT_BUNDLE_V1
    kappa_mu_ref: float = 2.01e-9
    kappa_lambda_ref: float = 0.2252
    p_c_y: float = 0.115
    y_t: float = 0.94
    theta12: float = 0.2274
    theta13: float = 0.00368
    theta23: float = 0.0415
    delta: float = 1.196
    gf_gev_minus2: float = 1.1663787e-5
    m_w_gev: float = 80.379
    alpha_em_mz: float = 1.0 / 127.952
    sin2_theta_w: float = 0.23122
    constants_citation: str = RARE_KAON_DILEPTON_PARAMETRIZATION_CITATION

    def __post_init__(self) -> None:
        for name in (
            "kappa_mu_ref",
            "kappa_lambda_ref",
            "p_c_y",
            "y_t",
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
class RareKaonDileptonCKMFactors:
    """CKM factors entering short-distance ``s -> d l+l-`` formulae."""

    lambda_wolfenstein: float
    lambda_c: complex
    lambda_t: complex
    matrix: tuple[tuple[complex, ...], ...]


@dataclass(frozen=True)
class RareKaonDileptonWilsonCoefficients:
    """Leading ``s -> d mu+mu-`` Wilson proxy from mass-basis couplings."""

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
    muon_axial_delta: float
    y_np_left: complex
    y_np_right: complex

    @property
    def y_np_total(self) -> complex:
        return complex(self.y_np_left - self.y_np_right)

    @property
    def wilsons(self) -> Mapping[str, complex]:
        return {
            "Y_NP_L": complex(self.y_np_left),
            "Y_NP_R": complex(self.y_np_right),
            "Y_NP_total": self.y_np_total,
        }


@dataclass(frozen=True)
class KLongMuMuShortDistanceResult:
    """Short-distance branching-ratio prediction for ``K_L -> mu+ mu-``."""

    model_label: str
    input_bundle: str
    branching_fraction: float
    sm_branching_fraction: float
    np_shift_branching_fraction: float
    kappa_mu: float
    p_c_y: float
    y_t: float
    y_eff: complex
    y_np_total: complex
    short_distance_amplitude: float
    lambda_wolfenstein: float
    lambda_c: complex
    lambda_t: complex
    wilsons: RareKaonDileptonWilsonCoefficients | None = None
    diagnostics: Mapping[str, float | complex | str | bool] = field(
        default_factory=dict
    )


def default_sm_inputs() -> RareKaonDileptonSMInputs:
    """Return the repo-owned default rare-kaon dilepton SM input bundle."""
    return RareKaonDileptonSMInputs()


def ckm_factors(
    inputs: RareKaonDileptonSMInputs | None = None,
) -> RareKaonDileptonCKMFactors:
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
    return RareKaonDileptonCKMFactors(
        lambda_wolfenstein=lam,
        lambda_c=lambda_c,
        lambda_t=lambda_t,
        matrix=matrix,
    )


def kappa_mu(inputs: RareKaonDileptonSMInputs | None = None) -> float:
    """Return ``kappa_mu`` with the Buras ``(lambda/0.2252)^8`` rescaling."""
    p = default_sm_inputs() if inputs is None else inputs
    lam = ckm_factors(p).lambda_wolfenstein
    return float(p.kappa_mu_ref * (lam / p.kappa_lambda_ref) ** 8)


def _positive_float(value: object, name: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0.0:
        raise ValueError(f"{name} must be positive and finite")
    return number


def _matrix_entry(source: object, matrix_name: str, i: int, j: int) -> complex:
    matrix = np.asarray(getattr(source, matrix_name), dtype=np.complex128)
    if matrix.shape != (3, 3):
        raise ValueError(f"{matrix_name} must have shape (3, 3)")
    value = complex(matrix[i, j])
    if not math.isfinite(value.real) or not math.isfinite(value.imag):
        raise ValueError(f"{matrix_name}[{i},{j}] must be finite")
    return value


def _weak_neutral_current_coupling(inputs: RareKaonDileptonSMInputs) -> float:
    g_weak = math.sqrt(4.0 * math.pi * inputs.alpha_em_mz / inputs.sin2_theta_w)
    cos_theta_w = math.sqrt(1.0 - inputs.sin2_theta_w)
    return float(g_weak / (2.0 * cos_theta_w))


def g_sm_squared(inputs: RareKaonDileptonSMInputs | None = None) -> float:
    """Return the Buras ``g_SM^2`` normalization in GeV^-2."""
    p = default_sm_inputs() if inputs is None else inputs
    return float(4.0 * p.gf_gev_minus2**2 * p.m_w_gev**2 / (2.0 * math.pi**2))


def compute_rare_kaon_dilepton_wilsons(
    source: QuarkMassBasisCouplings,
    *,
    m_kk_gev: float | None = None,
    inputs: RareKaonDileptonSMInputs | None = None,
) -> RareKaonDileptonWilsonCoefficients:
    """Match mass-basis couplings onto the v1 ``s -> d mu+mu-`` proxy.

    The flavor index convention follows the existing kaon Delta-F=2 and rare
    kaon paths: the ``(0, 1)`` down-sector entry is reported as the ``s-d``
    coupling.
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
    muon_axial_delta = neutral_delta
    normalizer = g_sm_squared(p) * resolved_m_kk**2
    y_np_left = left_quark_delta * muon_axial_delta / normalizer
    y_np_right = right_quark_delta * muon_axial_delta / normalizer

    return RareKaonDileptonWilsonCoefficients(
        model_label=RARE_KAON_DILEPTON_MODEL_V1,
        operator_convention=RARE_KAON_DILEPTON_OPERATOR_CONVENTION,
        matching_assumption=RARE_KAON_DILEPTON_RS_MATCHING_ASSUMPTION_V1,
        M_KK=resolved_m_kk,
        matching_scale=resolved_m_kk,
        left_sd_coupling=left_sd,
        right_sd_coupling=right_sd,
        left_sd_overlap=complex(left_overlap),
        right_sd_overlap=complex(right_overlap),
        left_quark_delta=complex(left_quark_delta),
        right_quark_delta=complex(right_quark_delta),
        muon_axial_delta=muon_axial_delta,
        y_np_left=complex(y_np_left),
        y_np_right=complex(y_np_right),
    )


def _klong_mumu_sd_from_y_np(
    y_np_total: complex,
    *,
    inputs: RareKaonDileptonSMInputs,
) -> tuple[float, float, float, RareKaonDileptonCKMFactors, complex]:
    factors = ckm_factors(inputs)
    lam = factors.lambda_wolfenstein
    y_eff = factors.lambda_t * inputs.y_t + complex(y_np_total)
    amplitude = y_eff.real / lam**5 + factors.lambda_c.real / lam * inputs.p_c_y
    kappa = kappa_mu(inputs)
    br = kappa * amplitude**2
    return float(br), float(kappa), float(amplitude), factors, y_eff


def klong_mumu_short_distance_sm(
    inputs: RareKaonDileptonSMInputs | None = None,
) -> KLongMuMuShortDistanceResult:
    """Evaluate the SM-limit ``K_L -> mu+ mu-`` short-distance branching fraction."""
    return evaluate_klong_mumu_short_distance(None, inputs=inputs)


def evaluate_klong_mumu_short_distance(
    source: QuarkMassBasisCouplings | RareKaonDileptonWilsonCoefficients | None = None,
    *,
    m_kk_gev: float | None = None,
    inputs: RareKaonDileptonSMInputs | None = None,
) -> KLongMuMuShortDistanceResult:
    """Evaluate ``BR(K_L -> mu+ mu-)_SD`` in the SM or with the v1 RS proxy."""
    p = default_sm_inputs() if inputs is None else inputs
    wilsons: RareKaonDileptonWilsonCoefficients | None
    if source is None:
        wilsons = None
        y_np_total = 0.0j
    elif isinstance(source, RareKaonDileptonWilsonCoefficients):
        wilsons = source
        y_np_total = source.y_np_total
    else:
        wilsons = compute_rare_kaon_dilepton_wilsons(
            source,
            m_kk_gev=m_kk_gev,
            inputs=p,
        )
        y_np_total = wilsons.y_np_total

    br, kappa, amplitude, factors, y_eff = _klong_mumu_sd_from_y_np(
        y_np_total,
        inputs=p,
    )
    sm, _, sm_amplitude, _, sm_y_eff = _klong_mumu_sd_from_y_np(0.0j, inputs=p)
    diagnostics: dict[str, float | complex | str | bool] = {
        "g_sm_squared_gev_minus2": g_sm_squared(p),
        "matching_assumption": RARE_KAON_DILEPTON_RS_MATCHING_ASSUMPTION_V1,
        "parametrization_citation": RARE_KAON_DILEPTON_PARAMETRIZATION_CITATION,
        "short_distance_only": True,
        "long_distance_dominated_total_rate": True,
        "right_handed_enters_with_minus_sign": True,
        "sm_short_distance_amplitude": float(sm_amplitude),
        "sm_y_eff": complex(sm_y_eff),
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
                "muon_axial_delta": float(wilsons.muon_axial_delta),
                "y_np_left": complex(wilsons.y_np_left),
                "y_np_right": complex(wilsons.y_np_right),
                "y_np_total": complex(wilsons.y_np_total),
            }
        )

    return KLongMuMuShortDistanceResult(
        model_label=RARE_KAON_DILEPTON_MODEL_V1,
        input_bundle=p.input_bundle,
        branching_fraction=br,
        sm_branching_fraction=sm,
        np_shift_branching_fraction=float(br - sm),
        kappa_mu=kappa,
        p_c_y=float(p.p_c_y),
        y_t=float(p.y_t),
        y_eff=complex(y_eff),
        y_np_total=complex(y_np_total),
        short_distance_amplitude=amplitude,
        lambda_wolfenstein=factors.lambda_wolfenstein,
        lambda_c=factors.lambda_c,
        lambda_t=factors.lambda_t,
        wilsons=wilsons,
        diagnostics=diagnostics,
    )


RARE_KAON_PI0EE_MODEL_V1 = "rare_kaon_pi0ee_isidori_smith_unterdorfer_rs_proxy_v1"
RARE_KAON_PI0EE_PARAMETRIZATION_CITATION = (
    "Isidori, Smith, Unterdorfer Eur. Phys. J. C 36 (2004) 57, "
    "arXiv:hep-ph/0404127; Buras-Mescia-Smith rare-kaon CPV "
    "decomposition convention"
)
RARE_KAON_PI0EE_INTERFERENCE_LIMITATION_V1 = (
    "NEEDS-HUMAN-PHYSICS: K_L -> pi0 e+e- direct/indirect interference "
    "depends on the ChPT sign of a_S and on long-distance inputs not present "
    "on ParameterPoint; the rigorous HARD-veto quantity is the direct-CP "
    "short-distance term, while indirect/interference/CPC totals are reported "
    "as diagnostic envelopes."
)


@dataclass(frozen=True)
class KLongPi0EEChPTInputs:
    """ChPT coefficient bundle for ``K_L -> pi0 e+ e-``.

    The coefficients are intentionally supplied by the catalog constraint from
    ``K008.yaml`` rather than read here.  With central values, the
    Isidori-Smith-Unterdorfer form is

        BR = [C_mix |a_S|^2 +/- C_int |a_S| A + C_dir A^2 + C_CPC] N,

    where ``A = Im(lambda_t) / 1e-4`` in the SM.  For the RS proxy we replace
    the direct-CP short-distance amplitude by
    ``A_eff = Im(lambda_t Y_t + Y_NP) / (Y_t 1e-4)`` so the SM limit reduces
    exactly to the published parametrization.
    """

    normalization: float
    c_mix: float
    c_int: float
    c_dir: float
    c_cpc: float
    a_s_abs: float
    sm_im_lambda_t_over_1e_minus4: float
    citation: str = RARE_KAON_PI0EE_PARAMETRIZATION_CITATION

    def __post_init__(self) -> None:
        for name in (
            "normalization",
            "c_mix",
            "c_int",
            "c_dir",
            "a_s_abs",
            "sm_im_lambda_t_over_1e_minus4",
        ):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be positive and finite")
        c_cpc = float(self.c_cpc)
        if not math.isfinite(c_cpc) or c_cpc < 0.0:
            raise ValueError("c_cpc must be finite and non-negative")


@dataclass(frozen=True)
class KLongPi0EEResult:
    """Branching-ratio prediction for ``K_L -> pi0 e+ e-``.

    ``direct_cp_branching_fraction`` is the short-distance quantity suitable
    for the K008 HARD verdict.  The constructive/destructive totals are the
    central ISU decomposition with both interference signs and are diagnostic
    only until the sign and CPC treatment are fixed by a human physics review.
    """

    model_label: str
    input_bundle: str
    direct_cp_branching_fraction: float
    sm_direct_cp_branching_fraction: float
    direct_cp_np_shift_branching_fraction: float
    constructive_total_branching_fraction: float
    destructive_total_branching_fraction: float
    sm_constructive_total_branching_fraction: float
    sm_destructive_total_branching_fraction: float
    indirect_cp_branching_fraction: float
    interference_abs_branching_fraction: float
    cpc_branching_fraction: float
    direct_cp_amplitude_ratio: float
    sm_direct_cp_amplitude_ratio: float
    y_eff: complex
    y_np_total: complex
    lambda_wolfenstein: float
    lambda_t: complex
    y_t: float
    chpt_inputs: KLongPi0EEChPTInputs
    wilsons: RareKaonDileptonWilsonCoefficients | None = None
    diagnostics: Mapping[str, float | complex | str | bool] = field(
        default_factory=dict
    )


def _klong_pi0ee_components_from_y_np(
    y_np_total: complex,
    *,
    inputs: RareKaonDileptonSMInputs,
    chpt_inputs: KLongPi0EEChPTInputs,
) -> tuple[
    float,
    float,
    float,
    float,
    float,
    RareKaonDileptonCKMFactors,
    complex,
    float,
]:
    factors = ckm_factors(inputs)
    y_eff = factors.lambda_t * inputs.y_t + complex(y_np_total)
    amplitude_ratio = y_eff.imag / (inputs.y_t * 1.0e-4)
    direct = chpt_inputs.normalization * chpt_inputs.c_dir * amplitude_ratio**2
    indirect = (
        chpt_inputs.normalization
        * chpt_inputs.c_mix
        * chpt_inputs.a_s_abs**2
    )
    interference_abs = (
        chpt_inputs.normalization
        * chpt_inputs.c_int
        * chpt_inputs.a_s_abs
        * abs(amplitude_ratio)
    )
    cpc = chpt_inputs.normalization * chpt_inputs.c_cpc
    constructive = indirect + direct + cpc + interference_abs
    destructive = max(indirect + direct + cpc - interference_abs, 0.0)
    return (
        float(direct),
        float(constructive),
        float(destructive),
        float(indirect),
        float(interference_abs),
        factors,
        complex(y_eff),
        float(amplitude_ratio),
    )


def klong_pi0ee_direct_cp_sm(
    chpt_inputs: KLongPi0EEChPTInputs,
    *,
    inputs: RareKaonDileptonSMInputs | None = None,
) -> KLongPi0EEResult:
    """Evaluate the SM-limit ``K_L -> pi0 e+ e-`` direct-CP rate."""
    return evaluate_klong_pi0ee_direct_cp(None, chpt_inputs=chpt_inputs, inputs=inputs)


def evaluate_klong_pi0ee_direct_cp(
    source: QuarkMassBasisCouplings | RareKaonDileptonWilsonCoefficients | None = None,
    *,
    chpt_inputs: KLongPi0EEChPTInputs,
    m_kk_gev: float | None = None,
    inputs: RareKaonDileptonSMInputs | None = None,
) -> KLongPi0EEResult:
    """Evaluate ``K_L -> pi0 e+ e-`` with the reused v1 ``Y_NP`` proxy."""
    p = default_sm_inputs() if inputs is None else inputs
    wilsons: RareKaonDileptonWilsonCoefficients | None
    if source is None:
        wilsons = None
        y_np_total = 0.0j
    elif isinstance(source, RareKaonDileptonWilsonCoefficients):
        wilsons = source
        y_np_total = source.y_np_total
    else:
        wilsons = compute_rare_kaon_dilepton_wilsons(
            source,
            m_kk_gev=m_kk_gev,
            inputs=p,
        )
        y_np_total = wilsons.y_np_total

    (
        direct,
        constructive,
        destructive,
        indirect,
        interference_abs,
        factors,
        y_eff,
        amplitude_ratio,
    ) = _klong_pi0ee_components_from_y_np(
        y_np_total,
        inputs=p,
        chpt_inputs=chpt_inputs,
    )
    (
        sm_direct,
        sm_constructive,
        sm_destructive,
        _sm_indirect,
        _sm_interference_abs,
        _sm_factors,
        _sm_y_eff,
        sm_amplitude_ratio,
    ) = _klong_pi0ee_components_from_y_np(
        0.0j,
        inputs=p,
        chpt_inputs=chpt_inputs,
    )
    cpc = chpt_inputs.normalization * chpt_inputs.c_cpc
    diagnostics: dict[str, float | complex | str | bool] = {
        "g_sm_squared_gev_minus2": g_sm_squared(p),
        "matching_assumption": RARE_KAON_DILEPTON_RS_MATCHING_ASSUMPTION_V1,
        "parametrization_citation": RARE_KAON_PI0EE_PARAMETRIZATION_CITATION,
        "interference_cpc_limitation": RARE_KAON_PI0EE_INTERFERENCE_LIMITATION_V1,
        "short_distance_direct_cp_only": True,
        "interference_sign_fixed": False,
        "cpc_from_two_photon_fixed": False,
        "uses_y_function_proxy": True,
        "right_handed_enters_with_minus_sign": True,
        "sm_y_eff": complex(_sm_y_eff),
        "sm_interference_abs_branching_fraction": float(_sm_interference_abs),
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
                "electron_axial_delta_proxy": float(wilsons.muon_axial_delta),
                "y_np_left": complex(wilsons.y_np_left),
                "y_np_right": complex(wilsons.y_np_right),
                "y_np_total": complex(wilsons.y_np_total),
            }
        )

    return KLongPi0EEResult(
        model_label=RARE_KAON_PI0EE_MODEL_V1,
        input_bundle=p.input_bundle,
        direct_cp_branching_fraction=direct,
        sm_direct_cp_branching_fraction=sm_direct,
        direct_cp_np_shift_branching_fraction=float(direct - sm_direct),
        constructive_total_branching_fraction=constructive,
        destructive_total_branching_fraction=destructive,
        sm_constructive_total_branching_fraction=sm_constructive,
        sm_destructive_total_branching_fraction=sm_destructive,
        indirect_cp_branching_fraction=indirect,
        interference_abs_branching_fraction=interference_abs,
        cpc_branching_fraction=float(cpc),
        direct_cp_amplitude_ratio=amplitude_ratio,
        sm_direct_cp_amplitude_ratio=sm_amplitude_ratio,
        y_eff=complex(y_eff),
        y_np_total=complex(y_np_total),
        lambda_wolfenstein=factors.lambda_wolfenstein,
        lambda_t=factors.lambda_t,
        y_t=float(p.y_t),
        chpt_inputs=chpt_inputs,
        wilsons=wilsons,
        diagnostics=diagnostics,
    )


__all__ = [
    "RARE_KAON_DILEPTON_MODEL_V1",
    "RARE_KAON_DILEPTON_OPERATOR_CONVENTION",
    "RARE_KAON_DILEPTON_INPUT_BUNDLE_V1",
    "RARE_KAON_DILEPTON_RS_MATCHING_ASSUMPTION_V1",
    "RARE_KAON_DILEPTON_PARAMETRIZATION_CITATION",
    "RARE_KAON_PI0EE_MODEL_V1",
    "RARE_KAON_PI0EE_PARAMETRIZATION_CITATION",
    "RARE_KAON_PI0EE_INTERFERENCE_LIMITATION_V1",
    "RareKaonDileptonSMInputs",
    "RareKaonDileptonCKMFactors",
    "RareKaonDileptonWilsonCoefficients",
    "KLongMuMuShortDistanceResult",
    "KLongPi0EEChPTInputs",
    "KLongPi0EEResult",
    "default_sm_inputs",
    "ckm_factors",
    "kappa_mu",
    "g_sm_squared",
    "compute_rare_kaon_dilepton_wilsons",
    "klong_mumu_short_distance_sm",
    "evaluate_klong_mumu_short_distance",
    "klong_pi0ee_direct_cp_sm",
    "evaluate_klong_pi0ee_direct_cp",
]


RARE_KAON_PI0EE_MODEL_V2 = "rare_kaon_pi0ee_y7v_y7a_isu_rs_proxy_v2"
RARE_KAON_PI0EE_Y7_OPERATOR_CONVENTION = (
    "Q7V=(sbar gamma_mu(1-gamma5)d)(ebar gamma^mu e), "
    "Q7A=(sbar gamma_mu(1-gamma5)d)(ebar gamma^mu gamma5 e); "
    "direct CP uses separate lambda_t*y7V and lambda_t*y7A amplitudes, "
    "interference uses the vector Q7V amplitude only"
)
RARE_KAON_PI0EE_RS_MATCHING_ASSUMPTION_V2 = (
    "NEEDS-HUMAN-PHYSICS: K_L -> pi0 e+e- now uses the Wilson-level "
    "y7V/y7A ISU structure, but the RS y7V/y7A values are proxy-matched "
    "from quark KK-gluon mass-basis flavor overlaps and one Z-like neutral "
    "boson. A complete EW KK/Z/Z' plus electron-sector matching is not "
    "available on ParameterPoint."
)
RARE_KAON_PI0EE_SEMILEPTONIC_RUNNING_DIAGNOSTIC_V1 = (
    "No additional multiplicative QCD running is applied to the NP "
    "semileptonic Q7V/Q7A proxy between M_KK and the low scale: the "
    "color-singlet vector/axial lepton currents have LO QCD factor 1.0 in "
    "this restricted semileptonic matching. The SM y7V/y7A inputs are the "
    "low-scale ISU values, with y7V quoted for mu=0.8-1.2 GeV; four-quark "
    "mixing into Q7V is not generated by this proxy and remains part of the "
    "NEEDS-HUMAN-PHYSICS matching limitation."
)
KLONG_PI0EE_SM_Y7V_BAR = 0.73
KLONG_PI0EE_SM_Y7A_BAR = -0.68
KLONG_PI0EE_C_INT_Y7V_COEFFICIENT = 8.91
KLONG_PI0EE_C_DIR_Y7VA_COEFFICIENT = 2.67
KLONG_PI0EE_Y7_LOW_SCALE_GEV = 1.0


@dataclass(frozen=True)
class KLongPi0EEWilsonCoefficients:
    """Wilson-level ``K_L -> pi0 e+e-`` proxy coefficients for K008.

    ``lambda_y7v_np_proxy`` and ``lambda_y7a_np_proxy`` are NP shifts in the
    same CP-odd amplitude slots as ``lambda_t * ybar_7V`` and
    ``lambda_t * ybar_7A``.  They are not claimed to be complete RS Wilson
    coefficients; the point is to preserve the correct ISU vector/axial
    structure while making the missing matching explicit.
    """

    model_label: str
    operator_convention: str
    matching_assumption: str
    M_KK: float
    matching_scale: float
    low_scale_gev: float
    left_sd_coupling: complex
    right_sd_coupling: complex
    left_sd_overlap: complex
    right_sd_overlap: complex
    left_quark_delta: complex
    right_quark_delta: complex
    quark_vector_delta: complex
    electron_vector_delta: float
    electron_axial_delta: float
    lambda_y7v_np_proxy: complex
    lambda_y7a_np_proxy: complex
    sm_y7v_bar: float = KLONG_PI0EE_SM_Y7V_BAR
    sm_y7a_bar: float = KLONG_PI0EE_SM_Y7A_BAR

    @property
    def wilsons(self) -> Mapping[str, complex | float]:
        return {
            "lambda_y7v_np_proxy": complex(self.lambda_y7v_np_proxy),
            "lambda_y7a_np_proxy": complex(self.lambda_y7a_np_proxy),
            "sm_y7v_bar": float(self.sm_y7v_bar),
            "sm_y7a_bar": float(self.sm_y7a_bar),
        }


def _electron_neutral_current_deltas(
    inputs: RareKaonDileptonSMInputs,
) -> tuple[float, float]:
    neutral_delta = _weak_neutral_current_coupling(inputs)
    electron_vector_delta = neutral_delta * (-1.0 + 4.0 * inputs.sin2_theta_w)
    electron_axial_delta = -neutral_delta
    return float(electron_vector_delta), float(electron_axial_delta)


def compute_klong_pi0ee_y7_wilsons(
    source: QuarkMassBasisCouplings,
    *,
    m_kk_gev: float | None = None,
    inputs: RareKaonDileptonSMInputs | None = None,
) -> KLongPi0EEWilsonCoefficients:
    """Match mass-basis couplings onto the K008 ``y7V/y7A`` proxy.

    For the pseudoscalar-to-pseudoscalar ``K -> pi`` matrix element the quark
    vector current is the relevant proxy structure, so left- and right-handed
    flavor-changing overlaps enter as a sum here.  This differs from the
    ``K_L -> mu+mu-`` axial-lepton ``Y_L - Y_R`` convention and is intentionally
    isolated to the appended K008 path.
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
    quark_vector_delta = left_quark_delta + right_quark_delta
    electron_vector_delta, electron_axial_delta = _electron_neutral_current_deltas(p)
    normalizer = g_sm_squared(p) * resolved_m_kk**2
    lambda_y7v_np_proxy = quark_vector_delta * electron_vector_delta / normalizer
    lambda_y7a_np_proxy = quark_vector_delta * electron_axial_delta / normalizer

    return KLongPi0EEWilsonCoefficients(
        model_label=RARE_KAON_PI0EE_MODEL_V2,
        operator_convention=RARE_KAON_PI0EE_Y7_OPERATOR_CONVENTION,
        matching_assumption=RARE_KAON_PI0EE_RS_MATCHING_ASSUMPTION_V2,
        M_KK=resolved_m_kk,
        matching_scale=resolved_m_kk,
        low_scale_gev=KLONG_PI0EE_Y7_LOW_SCALE_GEV,
        left_sd_coupling=left_sd,
        right_sd_coupling=right_sd,
        left_sd_overlap=complex(left_overlap),
        right_sd_overlap=complex(right_overlap),
        left_quark_delta=complex(left_quark_delta),
        right_quark_delta=complex(right_quark_delta),
        quark_vector_delta=complex(quark_vector_delta),
        electron_vector_delta=electron_vector_delta,
        electron_axial_delta=electron_axial_delta,
        lambda_y7v_np_proxy=complex(lambda_y7v_np_proxy),
        lambda_y7a_np_proxy=complex(lambda_y7a_np_proxy),
    )


def _klong_pi0ee_components_from_y7(
    lambda_y7v_np_proxy: complex,
    lambda_y7a_np_proxy: complex,
    *,
    inputs: RareKaonDileptonSMInputs,
    chpt_inputs: KLongPi0EEChPTInputs,
) -> tuple[
    float,
    float,
    float,
    float,
    float,
    RareKaonDileptonCKMFactors,
    complex,
    complex,
    float,
    float,
    float,
]:
    factors = ckm_factors(inputs)
    lambda_y7v_eff = (
        factors.lambda_t * KLONG_PI0EE_SM_Y7V_BAR + complex(lambda_y7v_np_proxy)
    )
    lambda_y7a_eff = (
        factors.lambda_t * KLONG_PI0EE_SM_Y7A_BAR + complex(lambda_y7a_np_proxy)
    )
    vector_ratio = lambda_y7v_eff.imag / 1.0e-4
    axial_ratio = lambda_y7a_eff.imag / 1.0e-4
    direct = (
        chpt_inputs.normalization
        * KLONG_PI0EE_C_DIR_Y7VA_COEFFICIENT
        * (vector_ratio**2 + axial_ratio**2)
    )
    indirect = (
        chpt_inputs.normalization
        * chpt_inputs.c_mix
        * chpt_inputs.a_s_abs**2
    )
    interference_abs = (
        chpt_inputs.normalization
        * KLONG_PI0EE_C_INT_Y7V_COEFFICIENT
        * chpt_inputs.a_s_abs
        * abs(vector_ratio)
    )
    cpc = chpt_inputs.normalization * chpt_inputs.c_cpc
    constructive = indirect + direct + cpc + interference_abs
    destructive = max(indirect + direct + cpc - interference_abs, 0.0)
    amplitude_quadrature = math.sqrt(vector_ratio**2 + axial_ratio**2)
    return (
        float(direct),
        float(constructive),
        float(destructive),
        float(indirect),
        float(interference_abs),
        factors,
        complex(lambda_y7v_eff),
        complex(lambda_y7a_eff),
        float(vector_ratio),
        float(axial_ratio),
        float(amplitude_quadrature),
    )


def klong_pi0ee_y7_direct_cp_sm(
    chpt_inputs: KLongPi0EEChPTInputs,
    *,
    inputs: RareKaonDileptonSMInputs | None = None,
) -> KLongPi0EEResult:
    """Evaluate the SM-limit K008 rate with explicit ``y7V/y7A`` structure."""
    return evaluate_klong_pi0ee_y7_direct_cp(
        None,
        chpt_inputs=chpt_inputs,
        inputs=inputs,
    )


def evaluate_klong_pi0ee_y7_direct_cp(
    source: QuarkMassBasisCouplings | KLongPi0EEWilsonCoefficients | None = None,
    *,
    chpt_inputs: KLongPi0EEChPTInputs,
    m_kk_gev: float | None = None,
    inputs: RareKaonDileptonSMInputs | None = None,
) -> KLongPi0EEResult:
    """Evaluate K008 using ISU ``y7V/y7A`` direct CP and vector interference."""
    p = default_sm_inputs() if inputs is None else inputs
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
        wilsons = compute_klong_pi0ee_y7_wilsons(
            source,
            m_kk_gev=m_kk_gev,
            inputs=p,
        )
        lambda_y7v_np_proxy = wilsons.lambda_y7v_np_proxy
        lambda_y7a_np_proxy = wilsons.lambda_y7a_np_proxy

    (
        direct,
        constructive,
        destructive,
        indirect,
        interference_abs,
        factors,
        lambda_y7v_eff,
        lambda_y7a_eff,
        vector_ratio,
        axial_ratio,
        amplitude_quadrature,
    ) = _klong_pi0ee_components_from_y7(
        lambda_y7v_np_proxy,
        lambda_y7a_np_proxy,
        inputs=p,
        chpt_inputs=chpt_inputs,
    )
    (
        sm_direct,
        sm_constructive,
        sm_destructive,
        _sm_indirect,
        _sm_interference_abs,
        _sm_factors,
        sm_lambda_y7v_eff,
        sm_lambda_y7a_eff,
        sm_vector_ratio,
        sm_axial_ratio,
        sm_amplitude_quadrature,
    ) = _klong_pi0ee_components_from_y7(
        0.0j,
        0.0j,
        inputs=p,
        chpt_inputs=chpt_inputs,
    )
    cpc = chpt_inputs.normalization * chpt_inputs.c_cpc
    diagnostics: dict[str, float | complex | str | bool] = {
        "g_sm_squared_gev_minus2": g_sm_squared(p),
        "matching_assumption": RARE_KAON_PI0EE_RS_MATCHING_ASSUMPTION_V2,
        "parametrization_citation": RARE_KAON_PI0EE_PARAMETRIZATION_CITATION,
        "interference_cpc_limitation": RARE_KAON_PI0EE_INTERFERENCE_LIMITATION_V1,
        "short_distance_direct_cp_only": True,
        "interference_sign_fixed": False,
        "cpc_from_two_photon_fixed": False,
        "uses_y_function_proxy": False,
        "uses_y7v_y7a_wilson_structure": True,
        "interference_uses_vector_y7v_only": True,
        "direct_uses_vector_and_axial_y7": True,
        "quark_vector_current_uses_left_plus_right": True,
        "legacy_single_axial_y_rescale_removed_for_k008": True,
        "sm_y7v_bar": float(KLONG_PI0EE_SM_Y7V_BAR),
        "sm_y7a_bar": float(KLONG_PI0EE_SM_Y7A_BAR),
        "c_int_y7v_coefficient": float(KLONG_PI0EE_C_INT_Y7V_COEFFICIENT),
        "c_dir_y7va_coefficient": float(KLONG_PI0EE_C_DIR_Y7VA_COEFFICIENT),
        "catalog_c_int_sm_coefficient": float(chpt_inputs.c_int),
        "catalog_c_dir_sm_coefficient": float(chpt_inputs.c_dir),
        "lambda_y7v_eff": complex(lambda_y7v_eff),
        "lambda_y7a_eff": complex(lambda_y7a_eff),
        "lambda_y7v_np_proxy": complex(lambda_y7v_np_proxy),
        "lambda_y7a_np_proxy": complex(lambda_y7a_np_proxy),
        "direct_cp_vector_amplitude_ratio": float(vector_ratio),
        "direct_cp_axial_amplitude_ratio": float(axial_ratio),
        "sm_lambda_y7v_eff": complex(sm_lambda_y7v_eff),
        "sm_lambda_y7a_eff": complex(sm_lambda_y7a_eff),
        "sm_direct_cp_vector_amplitude_ratio": float(sm_vector_ratio),
        "sm_direct_cp_axial_amplitude_ratio": float(sm_axial_ratio),
        "sm_interference_abs_branching_fraction": float(_sm_interference_abs),
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

    return KLongPi0EEResult(
        model_label=RARE_KAON_PI0EE_MODEL_V2,
        input_bundle=p.input_bundle,
        direct_cp_branching_fraction=direct,
        sm_direct_cp_branching_fraction=sm_direct,
        direct_cp_np_shift_branching_fraction=float(direct - sm_direct),
        constructive_total_branching_fraction=constructive,
        destructive_total_branching_fraction=destructive,
        sm_constructive_total_branching_fraction=sm_constructive,
        sm_destructive_total_branching_fraction=sm_destructive,
        indirect_cp_branching_fraction=indirect,
        interference_abs_branching_fraction=interference_abs,
        cpc_branching_fraction=float(cpc),
        direct_cp_amplitude_ratio=amplitude_quadrature,
        sm_direct_cp_amplitude_ratio=sm_amplitude_quadrature,
        y_eff=complex(lambda_y7a_eff),
        y_np_total=complex(lambda_y7a_np_proxy),
        lambda_wolfenstein=factors.lambda_wolfenstein,
        lambda_t=factors.lambda_t,
        y_t=float(p.y_t),
        chpt_inputs=chpt_inputs,
        wilsons=wilsons,
        diagnostics=diagnostics,
    )


__all__ += [
    "RARE_KAON_PI0EE_MODEL_V2",
    "RARE_KAON_PI0EE_Y7_OPERATOR_CONVENTION",
    "RARE_KAON_PI0EE_RS_MATCHING_ASSUMPTION_V2",
    "RARE_KAON_PI0EE_SEMILEPTONIC_RUNNING_DIAGNOSTIC_V1",
    "KLONG_PI0EE_SM_Y7V_BAR",
    "KLONG_PI0EE_SM_Y7A_BAR",
    "KLONG_PI0EE_C_INT_Y7V_COEFFICIENT",
    "KLONG_PI0EE_C_DIR_Y7VA_COEFFICIENT",
    "KLONG_PI0EE_Y7_LOW_SCALE_GEV",
    "KLongPi0EEWilsonCoefficients",
    "compute_klong_pi0ee_y7_wilsons",
    "klong_pi0ee_y7_direct_cp_sm",
    "evaluate_klong_pi0ee_y7_direct_cp",
]
