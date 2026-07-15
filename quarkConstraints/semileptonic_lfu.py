"""Semileptonic charged-current LFU-ratio proxy machinery.

Physics convention
------------------
This module owns reusable catalog-level machinery for ``b -> c tau nu`` LFU
ratios such as ``R(D)`` and ``R(D*)``.  The rigorous SM normalization is not
computed here from form factors; catalog constraints must pass the SM ratio
loaded from their YAML sidecar.  The core then applies a documented
amplitude-level charged-current stress proxy,

    R_X = R_X^SM |1 + C_tau^proxy|^2,

with

    C_tau^proxy = xi_cb^scalar m_b m_tau / (2 sqrt(2) M_KK^2).

``xi_cb^scalar`` is built from the current quark mass-basis coupling object by
dividing out ``g_s`` and multiplying the charm-left and bottom-right diagonal
overlap proxies.  This is intentionally a proxy, not an RS prediction.

RS matching assumption
----------------------
NEEDS-HUMAN-PHYSICS: a model-complete ``R(D)`` or ``R(D*)`` prediction needs
the charged electroweak KK/W' tower, charged-Higgs or leptoquark sector,
lepton and neutrino localization/couplings, the full ``b -> c tau nu`` WET
operator basis, and mode-specific form-factor integration.  Those inputs are
not present on ``ParameterPoint``.  The v1 proxy keeps only the expected
chirality-enhanced charged-current scaling ``m_b m_tau / M_KK^2`` and a
dimensionless quark-overlap stress factor.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Mapping

import numpy as np

from .couplings import QuarkMassBasisCouplings

SEMILEPTONIC_LFU_MODEL_V1 = "semileptonic_lfu_rd_charged_current_proxy_v1"
SEMILEPTONIC_LFU_OPERATOR_CONVENTION = (
    "amplitude proxy for (cbar P_L b)(taubar P_L nu_tau) normalized to "
    "the dimensionless m_b*m_tau/M_KK^2 charged-current stress scaling"
)
SEMILEPTONIC_LFU_INPUT_BUNDLE_V1 = "semileptonic_lfu_yaml_sm_ratio_proxy_inputs_v1"
SEMILEPTONIC_LFU_RS_MATCHING_ASSUMPTION_V1 = (
    "NEEDS-HUMAN-PHYSICS: quark KK-gluon mass-basis couplings are used only "
    "as charged-current flavor-overlap stress proxies; full W/W'/charged-"
    "Higgs/leptoquark, lepton, neutrino, and form-factor matching is absent."
)


@dataclass(frozen=True)
class SemileptonicLFUInputs:
    """Numerical inputs for semileptonic LFU ratio proxy evaluation."""

    input_bundle: str = SEMILEPTONIC_LFU_INPUT_BUNDLE_V1
    mode: str = "B->D"
    sm_lfu_ratio: float = 1.0
    bottom_mass_gev: float = 4.183
    tau_mass_gev: float = 1.77686
    gf_gev_minus2: float = 1.1663787e-5
    constants_citation: str = (
        "m_b(m_b) from PDG 2024 quark masses as used in "
        "quarkConstraints.pdg_quark_masses; tau mass and G_F use standard "
        "PDG electroweak constants. Constraint modules override sm_lfu_ratio "
        "from their YAML sidecar."
    )

    def __post_init__(self) -> None:
        for name in (
            "sm_lfu_ratio",
            "bottom_mass_gev",
            "tau_mass_gev",
            "gf_gev_minus2",
        ):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be positive and finite")


@dataclass(frozen=True)
class SemileptonicLFUWilsonProxy:
    """Leading charged-current LFU proxy from mass-basis couplings."""

    model_label: str
    operator_convention: str
    matching_assumption: str
    M_KK: float
    matching_scale: float
    charm_left_coupling: complex
    bottom_right_coupling: complex
    charm_left_overlap: complex
    bottom_right_overlap: complex
    scalar_overlap_proxy: complex
    scalar_amplitude_shift: complex

    @property
    def wilsons(self) -> Mapping[str, complex]:
        return {"C_tau_scalar_proxy": complex(self.scalar_amplitude_shift)}


@dataclass(frozen=True)
class SemileptonicLFUResult:
    """Prediction for one charged-current LFU ratio."""

    model_label: str
    input_bundle: str
    mode: str
    ratio: float
    sm_ratio: float
    np_shift_ratio: float
    response_factor: float
    scalar_amplitude_shift: complex
    wilsons: SemileptonicLFUWilsonProxy | None = None
    diagnostics: Mapping[str, float | complex | str | tuple[int, int]] = field(
        default_factory=dict
    )


def default_inputs() -> SemileptonicLFUInputs:
    """Return a generic input bundle; constraints should override ``sm_lfu_ratio``."""

    return SemileptonicLFUInputs()


def inputs_with_sm_ratio(
    sm_lfu_ratio: float,
    *,
    mode: str = "B->D",
    inputs: SemileptonicLFUInputs | None = None,
) -> SemileptonicLFUInputs:
    """Return inputs with the catalog-loaded SM LFU ratio installed."""

    base = default_inputs() if inputs is None else inputs
    return SemileptonicLFUInputs(
        input_bundle=base.input_bundle,
        mode=mode,
        sm_lfu_ratio=float(sm_lfu_ratio),
        bottom_mass_gev=float(base.bottom_mass_gev),
        tau_mass_gev=float(base.tau_mass_gev),
        gf_gev_minus2=float(base.gf_gev_minus2),
        constants_citation=base.constants_citation,
    )


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


def compute_semileptonic_lfu_wilson_proxy(
    source: QuarkMassBasisCouplings,
    *,
    m_kk_gev: float | None = None,
    inputs: SemileptonicLFUInputs | None = None,
) -> SemileptonicLFUWilsonProxy:
    """Return the v1 charged-current scalar LFU proxy for ``b -> c tau nu``.

    The quark index convention is ``c = up-sector index 1`` and ``b =
    down-sector index 2``.  The scalar overlap uses the charm-left and
    bottom-right diagonal overlap proxies because the intended stress direction
    is charged-Higgs/leptoquark-like and chirality enhanced by ``m_b m_tau``.
    """

    p = default_inputs() if inputs is None else inputs
    resolved_m_kk = _positive_float(
        getattr(source, "M_KK") if m_kk_gev is None else m_kk_gev,
        "m_kk_gev",
    )
    g_s = _positive_float(getattr(source, "g_s", 1.0), "g_s")
    charm_left = _matrix_entry(source, "left_up", 1, 1)
    bottom_right = _matrix_entry(source, "right_down", 2, 2)
    charm_left_overlap = charm_left / g_s
    bottom_right_overlap = bottom_right / g_s
    scalar_overlap = charm_left_overlap * bottom_right_overlap
    # Report 17/M-33: m_b*m_tau/M_KK^2 is already dimensionless; no extra 1/G_F.
    normalizer = 2.0 * math.sqrt(2.0) * resolved_m_kk**2
    scalar_shift = (
        scalar_overlap * p.bottom_mass_gev * p.tau_mass_gev / normalizer
    )

    return SemileptonicLFUWilsonProxy(
        model_label=SEMILEPTONIC_LFU_MODEL_V1,
        operator_convention=SEMILEPTONIC_LFU_OPERATOR_CONVENTION,
        matching_assumption=SEMILEPTONIC_LFU_RS_MATCHING_ASSUMPTION_V1,
        M_KK=float(resolved_m_kk),
        matching_scale=float(resolved_m_kk),
        charm_left_coupling=complex(charm_left),
        bottom_right_coupling=complex(bottom_right),
        charm_left_overlap=complex(charm_left_overlap),
        bottom_right_overlap=complex(bottom_right_overlap),
        scalar_overlap_proxy=complex(scalar_overlap),
        scalar_amplitude_shift=complex(scalar_shift),
    )


def ratio_from_scalar_shift(
    scalar_amplitude_shift: complex = 0.0j,
    *,
    inputs: SemileptonicLFUInputs | None = None,
) -> SemileptonicLFUResult:
    """Evaluate ``R_X = R_X^SM |1 + C_tau^proxy|^2``."""

    p = default_inputs() if inputs is None else inputs
    response = abs(1.0 + complex(scalar_amplitude_shift)) ** 2
    ratio = float(p.sm_lfu_ratio * response)
    diagnostics: dict[str, float | complex | str | tuple[int, int]] = {
        "matching_assumption": SEMILEPTONIC_LFU_RS_MATCHING_ASSUMPTION_V1,
        "operator_convention": SEMILEPTONIC_LFU_OPERATOR_CONVENTION,
        "constants_citation": p.constants_citation,
        "bottom_mass_gev": float(p.bottom_mass_gev),
        "tau_mass_gev": float(p.tau_mass_gev),
        "gf_gev_minus2": float(p.gf_gev_minus2),
        "response_factor": float(response),
        "scalar_amplitude_shift": complex(scalar_amplitude_shift),
    }
    return SemileptonicLFUResult(
        model_label=SEMILEPTONIC_LFU_MODEL_V1,
        input_bundle=p.input_bundle,
        mode=p.mode,
        ratio=ratio,
        sm_ratio=float(p.sm_lfu_ratio),
        np_shift_ratio=float(ratio - p.sm_lfu_ratio),
        response_factor=float(response),
        scalar_amplitude_shift=complex(scalar_amplitude_shift),
        wilsons=None,
        diagnostics=diagnostics,
    )


def sm_lfu_ratio(
    inputs: SemileptonicLFUInputs | None = None,
) -> SemileptonicLFUResult:
    """Evaluate the SM-limit LFU ratio from the supplied input bundle."""

    return ratio_from_scalar_shift(0.0j, inputs=inputs)


def evaluate_rd_lfu_ratio(
    source: QuarkMassBasisCouplings | SemileptonicLFUWilsonProxy | None = None,
    *,
    m_kk_gev: float | None = None,
    inputs: SemileptonicLFUInputs | None = None,
) -> SemileptonicLFUResult:
    """Evaluate ``R(D)`` or a sibling LFU ratio with the v1 proxy."""

    p = default_inputs() if inputs is None else inputs
    wilsons: SemileptonicLFUWilsonProxy | None
    if source is None:
        wilsons = None
        scalar_shift = 0.0j
    elif isinstance(source, SemileptonicLFUWilsonProxy):
        wilsons = source
        scalar_shift = source.scalar_amplitude_shift
    else:
        wilsons = compute_semileptonic_lfu_wilson_proxy(
            source,
            m_kk_gev=m_kk_gev,
            inputs=p,
        )
        scalar_shift = wilsons.scalar_amplitude_shift

    result = ratio_from_scalar_shift(scalar_shift, inputs=p)
    diagnostics = dict(result.diagnostics)
    if wilsons is not None:
        diagnostics.update(
            {
                "m_kk_gev": float(wilsons.M_KK),
                "matching_scale_gev": float(wilsons.matching_scale),
                "charm_left_coupling": complex(wilsons.charm_left_coupling),
                "bottom_right_coupling": complex(wilsons.bottom_right_coupling),
                "charm_left_overlap": complex(wilsons.charm_left_overlap),
                "bottom_right_overlap": complex(wilsons.bottom_right_overlap),
                "scalar_overlap_proxy": complex(wilsons.scalar_overlap_proxy),
                "up_down_sector_indices": (1, 2),
            }
        )
    return SemileptonicLFUResult(
        model_label=result.model_label,
        input_bundle=result.input_bundle,
        mode=result.mode,
        ratio=float(result.ratio),
        sm_ratio=float(result.sm_ratio),
        np_shift_ratio=float(result.np_shift_ratio),
        response_factor=float(result.response_factor),
        scalar_amplitude_shift=complex(result.scalar_amplitude_shift),
        wilsons=wilsons,
        diagnostics=diagnostics,
    )


__all__ = [
    "SEMILEPTONIC_LFU_MODEL_V1",
    "SEMILEPTONIC_LFU_OPERATOR_CONVENTION",
    "SEMILEPTONIC_LFU_INPUT_BUNDLE_V1",
    "SEMILEPTONIC_LFU_RS_MATCHING_ASSUMPTION_V1",
    "SemileptonicLFUInputs",
    "SemileptonicLFUWilsonProxy",
    "SemileptonicLFUResult",
    "default_inputs",
    "inputs_with_sm_ratio",
    "compute_semileptonic_lfu_wilson_proxy",
    "ratio_from_scalar_shift",
    "sm_lfu_ratio",
    "evaluate_rd_lfu_ratio",
]
