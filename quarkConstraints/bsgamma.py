"""Shared ``b -> s gamma`` dipole machinery.

Physics convention
------------------
The inclusive ``Bbar -> X_s gamma`` rate is represented by the leading
electromagnetic dipole scaling

    BR = BR_SM * (|C7_SM + C7_NP|^2 + |C7p_SM + C7p_NP|^2)
              / (|C7_SM|^2 + |C7p_SM|^2),

where the caller supplies the NNLO SM branching fraction for the desired
photon-energy cut.  The default effective SM coefficient is
``C7_SM(mu_b) = -0.304``.  New-physics dipoles matched at ``M_KK`` are evolved
to ``mu_b`` before entering this low-scale rate; this module does not hardcode
any experimental or SM branching-fraction anchor.

RS matching assumption
----------------------
NEEDS-HUMAN-PHYSICS: the current ParameterPoint carries quark mass-basis
``KK-gluon``-style coupling matrices, but not the Higgs/Goldstone, KK-fermion,
charged-current, or chromomagnetic matching inputs needed for a model-complete
RS prediction of ``C7`` and ``C8``.  The v1 matching below is therefore a
documented dipole-coefficient proxy: divide the supplied b-s quark coupling by
``g_s`` to keep the flavor-overlap structure, scale it as
``(3 TeV / M_KK)^2``, and identify the left-handed b-s overlap with high-scale
``C7_NP``/``C8_NP`` and the right-handed overlap with high-scale
``C7p_NP``/``C8p_NP``.  The proxy normalizations are order-one coefficients set
to 1 by default and exposed in diagnostics.  The high-scale dipoles are then
run to ``mu_b`` with the standard leading-log ``C7``-``C8`` mixing.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Mapping

import numpy as np

from .couplings import QuarkMassBasisCouplings
from .qcd_running import run_alpha_s

BSGAMMA_MODEL_V1 = "bsgamma_c7_dipole_rs_proxy_v1"
BSGAMMA_OPERATOR_CONVENTION = (
    "O7=(e/16 pi^2) m_b (sbar_L sigma^{mu nu} b_R) F_{mu nu}; "
    "O7p is the L<->R flipped operator"
)
BSGAMMA_INPUT_BUNDLE_V1 = "bsgamma_c7eff_mub_rate_normalized_v1"
BSGAMMA_RS_MATCHING_ASSUMPTION_V1 = (
    "NEEDS-HUMAN-PHYSICS: quark KK-gluon mass-basis couplings are used as "
    "b-s flavor-overlap proxies for loop-level RS dipoles; full KK fermion, "
    "Higgs/Goldstone, W/charged-current, and finite C7/C8 matching is not "
    "present. The proxy C7/C8 dipoles are evolved with leading-log QCD running."
)
BSGAMMA_LL_RUNNING_LABEL = "bsgamma_c7_c8_ll_qcd_running_v1"

_GAMMA_77_0 = 32.0 / 3.0
_GAMMA_88_0 = 28.0 / 3.0
_C8_TO_C7_MIXING = 8.0 / 3.0


@dataclass(frozen=True)
class BsgammaSMInputs:
    """Numerical inputs for the C7-normalized inclusive-rate proxy."""

    input_bundle: str = BSGAMMA_INPUT_BUNDLE_V1
    c7_sm_eff: complex = -0.304 + 0.0j
    c7p_sm_eff: complex = 0.0j
    reference_scale_gev: float = 3000.0
    low_scale_gev: float = 4.8
    photon_energy_cut_gev: float = 1.6
    c7_proxy_normalization: float = 1.0
    c8_proxy_normalization: float = 1.0
    alpha_s_mz: float = 0.1179
    m_z_gev: float = 91.1876
    m_t_gev: float = 163.5
    m_b_gev: float = 4.18
    m_c_gev: float = 1.27
    constants_citation: str = (
        "C7_eff(mu_b) representative low-scale SM value used for a normalized "
        "dipole-rate proxy; branching-fraction normalization is supplied by "
        "the catalog anchor, e.g. Misiak-Rehman-Steinhauser arXiv:2002.01548."
    )

    def __post_init__(self) -> None:
        for name in (
            "reference_scale_gev",
            "low_scale_gev",
            "photon_energy_cut_gev",
            "alpha_s_mz",
            "m_z_gev",
            "m_t_gev",
            "m_b_gev",
            "m_c_gev",
        ):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be positive and finite")
        if not (self.m_c_gev < self.m_b_gev < self.m_z_gev < self.m_t_gev):
            raise ValueError("thresholds must satisfy m_c < m_b < m_z < m_t")
        for name in ("c7_proxy_normalization", "c8_proxy_normalization"):
            norm = float(getattr(self, name))
            if not math.isfinite(norm):
                raise ValueError(f"{name} must be finite")
        denominator = abs(complex(self.c7_sm_eff)) ** 2 + abs(
            complex(self.c7p_sm_eff)
        ) ** 2
        if denominator <= 0.0 or not math.isfinite(denominator):
            raise ValueError("SM C7 denominator must be positive and finite")


@dataclass(frozen=True)
class BsgammaWilsonCoefficients:
    """Leading b-s dipole proxy from mass-basis couplings."""

    model_label: str
    operator_convention: str
    matching_assumption: str
    M_KK: float
    matching_scale: float
    low_scale: float
    left_bs_coupling: complex
    right_bs_coupling: complex
    left_bs_overlap: complex
    right_bs_overlap: complex
    proxy_scale_factor: float
    c7_np_matching: complex
    c7p_np_matching: complex
    c8_np_matching: complex
    c8p_np_matching: complex
    c7_np: complex
    c7p_np: complex
    c8_np: complex
    c8p_np: complex
    c7_running_from_c7: float
    c7_running_from_c8: float
    c8_running_from_c8: float
    alpha_s_matching_scale: float
    alpha_s_low_scale: float

    @property
    def wilsons(self) -> Mapping[str, complex]:
        return {
            "C7_NP": complex(self.c7_np),
            "C7p_NP": complex(self.c7p_np),
            "C8_NP": complex(self.c8_np),
            "C8p_NP": complex(self.c8p_np),
            "C7_NP_MKK": complex(self.c7_np_matching),
            "C7p_NP_MKK": complex(self.c7p_np_matching),
            "C8_NP_MKK": complex(self.c8_np_matching),
            "C8p_NP_MKK": complex(self.c8p_np_matching),
        }


@dataclass(frozen=True)
class BsgammaBranchingResult:
    """Branching-ratio prediction for inclusive ``Bbar -> X_s gamma``."""

    model_label: str
    input_bundle: str
    branching_fraction: float
    sm_branching_fraction: float
    np_shift_branching_fraction: float
    ratio_to_sm: float
    c7_sm_eff: complex
    c7p_sm_eff: complex
    c7_total: complex
    c7p_total: complex
    c7_np: complex
    c7p_np: complex
    wilsons: BsgammaWilsonCoefficients | None = None
    diagnostics: Mapping[str, float | complex | str] = field(default_factory=dict)


def default_sm_inputs() -> BsgammaSMInputs:
    """Return the default C7-normalized inclusive-rate input bundle."""

    return BsgammaSMInputs()


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


def _dipole_power(c7: complex, c7p: complex) -> float:
    power = abs(complex(c7)) ** 2 + abs(complex(c7p)) ** 2
    if not math.isfinite(power):
        raise ValueError("dipole power must be finite")
    return float(power)


def _beta0(nf: int) -> float:
    return (33.0 - 2.0 * nf) / 3.0


def _nf_for_bsgamma_scale(mu: float, inputs: BsgammaSMInputs) -> int:
    if mu > inputs.m_t_gev:
        return 6
    if mu > inputs.m_b_gev:
        return 5
    if mu > inputs.m_c_gev:
        return 4
    return 3


def _alpha_s(mu: float, inputs: BsgammaSMInputs) -> float:
    return run_alpha_s(
        mu,
        alpha_s_mz=float(inputs.alpha_s_mz),
        m_z=float(inputs.m_z_gev),
        m_t=float(inputs.m_t_gev),
        m_b=float(inputs.m_b_gev),
        m_c=float(inputs.m_c_gev),
    )


def _segment_c7_c8_running(
    *,
    mu_high: float,
    mu_low: float,
    nf: int,
    inputs: BsgammaSMInputs,
) -> tuple[float, float, float]:
    alpha_high = _alpha_s(mu_high, inputs)
    alpha_low = _alpha_s(mu_low, inputs)
    eta = alpha_high / alpha_low
    b0 = _beta0(nf)
    c7_from_c7 = eta ** (_GAMMA_77_0 / (2.0 * b0))
    c8_from_c8 = eta ** (_GAMMA_88_0 / (2.0 * b0))
    c7_from_c8 = _C8_TO_C7_MIXING * (c8_from_c8 - c7_from_c7)
    return float(c7_from_c7), float(c7_from_c8), float(c8_from_c8)


def _ll_bsgamma_running_coefficients(
    *,
    matching_scale: float,
    low_scale: float,
    inputs: BsgammaSMInputs,
) -> tuple[float, float, float]:
    if matching_scale <= low_scale:
        raise ValueError("matching_scale must be greater than low_scale")

    thresholds = sorted(
        [
            scale
            for scale in (inputs.m_t_gev, inputs.m_b_gev, inputs.m_c_gev)
            if low_scale < scale < matching_scale
        ],
        reverse=True,
    )
    boundaries = [matching_scale, *thresholds, low_scale]
    u77_total = 1.0
    u78_total = 0.0
    u88_total = 1.0

    for mu_high, mu_low in zip(boundaries, boundaries[1:]):
        nf = _nf_for_bsgamma_scale(mu_high, inputs)
        u77, u78, u88 = _segment_c7_c8_running(
            mu_high=mu_high,
            mu_low=mu_low,
            nf=nf,
            inputs=inputs,
        )
        u78_total = u77 * u78_total + u78 * u88_total
        u77_total *= u77
        u88_total *= u88

    return float(u77_total), float(u78_total), float(u88_total)


def _run_c7_c8_to_low_scale(
    *,
    c7_matching: complex,
    c8_matching: complex,
    matching_scale: float,
    inputs: BsgammaSMInputs,
) -> tuple[complex, complex, float, float, float]:
    u77, u78, u88 = _ll_bsgamma_running_coefficients(
        matching_scale=matching_scale,
        low_scale=float(inputs.low_scale_gev),
        inputs=inputs,
    )
    c7_low = u77 * complex(c7_matching) + u78 * complex(c8_matching)
    c8_low = u88 * complex(c8_matching)
    return complex(c7_low), complex(c8_low), u77, u78, u88


def _branching_from_dipoles(
    *,
    c7_np: complex,
    c7p_np: complex,
    sm_branching_fraction: float,
    inputs: BsgammaSMInputs,
) -> tuple[float, float, complex, complex]:
    sm_br = _positive_float(sm_branching_fraction, "sm_branching_fraction")
    c7_total = complex(inputs.c7_sm_eff) + complex(c7_np)
    c7p_total = complex(inputs.c7p_sm_eff) + complex(c7p_np)
    sm_power = _dipole_power(inputs.c7_sm_eff, inputs.c7p_sm_eff)
    total_power = _dipole_power(c7_total, c7p_total)
    ratio_to_sm = total_power / sm_power
    return float(sm_br * ratio_to_sm), float(ratio_to_sm), c7_total, c7p_total


def compute_bsgamma_wilsons(
    source: QuarkMassBasisCouplings,
    *,
    m_kk_gev: float | None = None,
    inputs: BsgammaSMInputs | None = None,
) -> BsgammaWilsonCoefficients:
    """Match mass-basis couplings onto the v1 ``b -> s gamma`` C7 proxy.

    The down-sector flavor index convention follows the B_s Delta-F=2 path:
    the ``(1, 2)`` entry is the s-b coupling relevant for ``b -> s``.
    """

    p = default_sm_inputs() if inputs is None else inputs
    resolved_m_kk = _positive_float(
        getattr(source, "M_KK") if m_kk_gev is None else m_kk_gev,
        "m_kk_gev",
    )
    left_bs = _matrix_entry(source, "left_down", 1, 2)
    right_bs = _matrix_entry(source, "right_down", 1, 2)
    g_s = _positive_float(getattr(source, "g_s", 1.0), "g_s")
    left_overlap = left_bs / g_s
    right_overlap = right_bs / g_s
    proxy_scale = float(
        p.c7_proxy_normalization * (p.reference_scale_gev / resolved_m_kk) ** 2
    )
    c7_matching = proxy_scale * left_overlap
    c7p_matching = proxy_scale * right_overlap
    c8_matching = float(p.c8_proxy_normalization) * proxy_scale * left_overlap
    c8p_matching = float(p.c8_proxy_normalization) * proxy_scale * right_overlap
    c7_np, c8_np, u77, u78, u88 = _run_c7_c8_to_low_scale(
        c7_matching=c7_matching,
        c8_matching=c8_matching,
        matching_scale=resolved_m_kk,
        inputs=p,
    )
    c7p_np, c8p_np, _u77p, _u78p, _u88p = _run_c7_c8_to_low_scale(
        c7_matching=c7p_matching,
        c8_matching=c8p_matching,
        matching_scale=resolved_m_kk,
        inputs=p,
    )

    return BsgammaWilsonCoefficients(
        model_label=BSGAMMA_MODEL_V1,
        operator_convention=BSGAMMA_OPERATOR_CONVENTION,
        matching_assumption=BSGAMMA_RS_MATCHING_ASSUMPTION_V1,
        M_KK=resolved_m_kk,
        matching_scale=resolved_m_kk,
        low_scale=float(p.low_scale_gev),
        left_bs_coupling=left_bs,
        right_bs_coupling=right_bs,
        left_bs_overlap=complex(left_overlap),
        right_bs_overlap=complex(right_overlap),
        proxy_scale_factor=proxy_scale,
        c7_np_matching=complex(c7_matching),
        c7p_np_matching=complex(c7p_matching),
        c8_np_matching=complex(c8_matching),
        c8p_np_matching=complex(c8p_matching),
        c7_np=complex(c7_np),
        c7p_np=complex(c7p_np),
        c8_np=complex(c8_np),
        c8p_np=complex(c8p_np),
        c7_running_from_c7=u77,
        c7_running_from_c8=u78,
        c8_running_from_c8=u88,
        alpha_s_matching_scale=float(_alpha_s(resolved_m_kk, p)),
        alpha_s_low_scale=float(_alpha_s(float(p.low_scale_gev), p)),
    )


def branching_fraction_from_c7(
    *,
    c7_np: complex = 0.0j,
    c7p_np: complex = 0.0j,
    sm_branching_fraction: float,
    inputs: BsgammaSMInputs | None = None,
) -> BsgammaBranchingResult:
    """Evaluate the inclusive branching fraction from explicit C7 shifts."""

    p = default_sm_inputs() if inputs is None else inputs
    br, ratio_to_sm, c7_total, c7p_total = _branching_from_dipoles(
        c7_np=complex(c7_np),
        c7p_np=complex(c7p_np),
        sm_branching_fraction=sm_branching_fraction,
        inputs=p,
    )
    sm_br = _positive_float(sm_branching_fraction, "sm_branching_fraction")
    diagnostics: dict[str, float | complex | str] = {
        "branching_formula": (
            "BR = BR_SM * (|C7_SM + C7_NP|^2 + |C7p_SM + C7p_NP|^2) "
            "/ (|C7_SM|^2 + |C7p_SM|^2)"
        ),
        "photon_energy_cut_gev": float(p.photon_energy_cut_gev),
        "reference_scale_gev": float(p.reference_scale_gev),
        "low_scale_gev": float(p.low_scale_gev),
        "c7_proxy_normalization": float(p.c7_proxy_normalization),
        "c8_proxy_normalization": float(p.c8_proxy_normalization),
        "dipole_power_sm": _dipole_power(p.c7_sm_eff, p.c7p_sm_eff),
        "dipole_power_total": _dipole_power(c7_total, c7p_total),
        "constants_citation": p.constants_citation,
        "wilson_scale_assumption": (
            "c7_np and c7p_np are low-scale coefficients at mu_b"
        ),
    }
    return BsgammaBranchingResult(
        model_label=BSGAMMA_MODEL_V1,
        input_bundle=p.input_bundle,
        branching_fraction=br,
        sm_branching_fraction=sm_br,
        np_shift_branching_fraction=float(br - sm_br),
        ratio_to_sm=ratio_to_sm,
        c7_sm_eff=complex(p.c7_sm_eff),
        c7p_sm_eff=complex(p.c7p_sm_eff),
        c7_total=c7_total,
        c7p_total=c7p_total,
        c7_np=complex(c7_np),
        c7p_np=complex(c7p_np),
        diagnostics=diagnostics,
    )


def evaluate_inclusive_bsgamma(
    source: QuarkMassBasisCouplings | BsgammaWilsonCoefficients | None = None,
    *,
    sm_branching_fraction: float,
    m_kk_gev: float | None = None,
    inputs: BsgammaSMInputs | None = None,
) -> BsgammaBranchingResult:
    """Evaluate inclusive ``Bbar -> X_s gamma`` in the SM or v1 C7 proxy."""

    p = default_sm_inputs() if inputs is None else inputs
    wilsons: BsgammaWilsonCoefficients | None
    if source is None:
        wilsons = None
        c7_np = 0.0j
        c7p_np = 0.0j
    elif isinstance(source, BsgammaWilsonCoefficients):
        wilsons = source
        c7_np = source.c7_np
        c7p_np = source.c7p_np
    else:
        wilsons = compute_bsgamma_wilsons(source, m_kk_gev=m_kk_gev, inputs=p)
        c7_np = wilsons.c7_np
        c7p_np = wilsons.c7p_np

    result = branching_fraction_from_c7(
        c7_np=c7_np,
        c7p_np=c7p_np,
        sm_branching_fraction=sm_branching_fraction,
        inputs=p,
    )
    diagnostics = dict(result.diagnostics)
    diagnostics["matching_assumption"] = BSGAMMA_RS_MATCHING_ASSUMPTION_V1
    if wilsons is not None:
        diagnostics.update(
            {
                "m_kk_gev": float(wilsons.M_KK),
                "matching_scale_gev": float(wilsons.matching_scale),
                "low_scale_gev": float(wilsons.low_scale),
                "left_bs_coupling": complex(wilsons.left_bs_coupling),
                "right_bs_coupling": complex(wilsons.right_bs_coupling),
                "left_bs_overlap": complex(wilsons.left_bs_overlap),
                "right_bs_overlap": complex(wilsons.right_bs_overlap),
                "proxy_scale_factor": float(wilsons.proxy_scale_factor),
                "c7_np_matching": complex(wilsons.c7_np_matching),
                "c7p_np_matching": complex(wilsons.c7p_np_matching),
                "c8_np_matching": complex(wilsons.c8_np_matching),
                "c8p_np_matching": complex(wilsons.c8p_np_matching),
                "c7_np": complex(wilsons.c7_np),
                "c7p_np": complex(wilsons.c7p_np),
                "c8_np": complex(wilsons.c8_np),
                "c8p_np": complex(wilsons.c8p_np),
                "bsgamma_rg_running_applied": True,
                "bsgamma_rg_running_label": BSGAMMA_LL_RUNNING_LABEL,
                "c7_running_from_c7": float(wilsons.c7_running_from_c7),
                "c7_running_from_c8": float(wilsons.c7_running_from_c8),
                "c8_running_from_c8": float(wilsons.c8_running_from_c8),
                "alpha_s_matching_scale": float(wilsons.alpha_s_matching_scale),
                "alpha_s_low_scale": float(wilsons.alpha_s_low_scale),
            }
        )
    return BsgammaBranchingResult(
        model_label=result.model_label,
        input_bundle=result.input_bundle,
        branching_fraction=result.branching_fraction,
        sm_branching_fraction=result.sm_branching_fraction,
        np_shift_branching_fraction=result.np_shift_branching_fraction,
        ratio_to_sm=result.ratio_to_sm,
        c7_sm_eff=result.c7_sm_eff,
        c7p_sm_eff=result.c7p_sm_eff,
        c7_total=result.c7_total,
        c7p_total=result.c7p_total,
        c7_np=result.c7_np,
        c7p_np=result.c7p_np,
        wilsons=wilsons,
        diagnostics=diagnostics,
    )
