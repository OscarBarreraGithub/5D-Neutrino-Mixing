"""Adapter over :mod:`quarkConstraints.bsgamma`.

This is the catalog boundary for inclusive ``Bbar -> X_s gamma`` and reusable
``b -> s gamma`` C7 dipole machinery.  Constraint modules import this adapter
only; the underlying physics implementation remains isolated in
``quarkConstraints``.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from quarkConstraints.bsgamma import (
    BSGAMMA_INPUT_BUNDLE_V1,
    BSGAMMA_MODEL_V1,
    BSGAMMA_OPERATOR_CONVENTION,
    BSGAMMA_RS_MATCHING_ASSUMPTION_V1,
    BsgammaBranchingResult,
    BsgammaSMInputs,
    BsgammaWilsonCoefficients,
    branching_fraction_from_c7 as _branching_fraction_from_c7,
    compute_bsgamma_wilsons as _compute_bsgamma_wilsons,
    default_sm_inputs as _default_sm_inputs,
    evaluate_inclusive_bsgamma as _evaluate_inclusive_bsgamma,
)
from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.modern.inputs import ModernDefaultCKMTarget
from quarkConstraints.model import RotationParameters, ckm_like_unitary

__all__ = [
    "QuarkMassBasisCouplings",
    "BSGAMMA_MODEL_V1",
    "BSGAMMA_OPERATOR_CONVENTION",
    "BSGAMMA_INPUT_BUNDLE_V1",
    "BSGAMMA_RS_MATCHING_ASSUMPTION_V1",
    "BsgammaSMInputs",
    "BsgammaWilsonCoefficients",
    "BsgammaBranchingResult",
    "BdgammaCKMFactors",
    "bdgamma_ckm_factors",
    "bdgamma_wilsons_from_couplings",
    "bsgamma_default_sm_inputs",
    "bsgamma_branching_fraction_from_c7",
    "bsgamma_wilsons_from_couplings",
    "exclusive_bdrhogamma_from_c7",
    "exclusive_bdrhogamma_from_couplings",
    "exclusive_bdrhogamma_sm_branching_fraction",
    "exclusive_btokstargamma_from_c7",
    "exclusive_btokstargamma_from_couplings",
    "exclusive_btokstargamma_sm_branching_fraction",
    "inclusive_bsgamma_from_couplings",
    "inclusive_bsgamma_sm_branching_fraction",
]


@dataclass(frozen=True)
class BdgammaCKMFactors:
    """Repo-default CKM normalization factors for ``b -> d gamma``."""

    input_bundle: str
    lambda_t_bd: complex
    lambda_t_bs: complex
    abs_vtd_over_vts: float
    abs_lambda_t_bd_over_lambda_t_bs: float
    ckm_power_suppression: float
    matrix: tuple[tuple[complex, ...], ...]


def _ckm_matrix_tuple(matrix: np.ndarray) -> tuple[tuple[complex, ...], ...]:
    return tuple(tuple(complex(value) for value in row) for row in matrix)


def bdgamma_ckm_factors() -> BdgammaCKMFactors:
    """Return ``b -> d`` CKM factors relative to the shared ``b -> s`` setup."""

    target = ModernDefaultCKMTarget()
    matrix = ckm_like_unitary(
        RotationParameters(
            theta12=float(target.theta12),
            theta13=float(target.theta13),
            theta23=float(target.theta23),
            delta=float(target.delta),
        )
    )
    lambda_t_bd = complex(np.conjugate(matrix[2, 0]) * matrix[2, 2])
    lambda_t_bs = complex(np.conjugate(matrix[2, 1]) * matrix[2, 2])
    abs_vts = abs(complex(matrix[2, 1]))
    if abs_vts <= 0.0 or abs(lambda_t_bs) <= 0.0:
        raise ValueError("repo CKM target has zero V_ts/lambda_t_bs")
    abs_vtd_over_vts = float(abs(complex(matrix[2, 0])) / abs_vts)
    abs_lambda_ratio = float(abs(lambda_t_bd) / abs(lambda_t_bs))
    if not math.isfinite(abs_vtd_over_vts) or not math.isfinite(abs_lambda_ratio):
        raise ValueError("b -> d CKM ratio must be finite")
    return BdgammaCKMFactors(
        input_bundle=target.target_id,
        lambda_t_bd=lambda_t_bd,
        lambda_t_bs=lambda_t_bs,
        abs_vtd_over_vts=abs_vtd_over_vts,
        abs_lambda_t_bd_over_lambda_t_bs=abs_lambda_ratio,
        ckm_power_suppression=float(abs_lambda_ratio * abs_lambda_ratio),
        matrix=_ckm_matrix_tuple(matrix),
    )


def _matrix(source: object, matrix_name: str) -> np.ndarray:
    matrix = np.asarray(getattr(source, matrix_name), dtype=np.complex128)
    if matrix.shape != (3, 3):
        raise ValueError(f"{matrix_name} must have shape (3, 3)")
    return matrix


def _bd_source_as_bs_source(source: QuarkMassBasisCouplings) -> QuarkMassBasisCouplings:
    """Map the d-b coupling slot into the existing b-s C7 proxy input slot."""

    left_down = _matrix(source, "left_down").copy()
    right_down = _matrix(source, "right_down").copy()
    left_bd = complex(left_down[0, 2])
    right_bd = complex(right_down[0, 2])
    if not (
        math.isfinite(left_bd.real)
        and math.isfinite(left_bd.imag)
        and math.isfinite(right_bd.real)
        and math.isfinite(right_bd.imag)
    ):
        raise ValueError("b-d dipole proxy couplings must be finite")

    left_down[1, 2] = left_bd
    left_down[2, 1] = np.conjugate(left_bd)
    right_down[1, 2] = right_bd
    right_down[2, 1] = np.conjugate(right_bd)

    return QuarkMassBasisCouplings(
        M_KK=float(getattr(source, "M_KK")),
        xi_KK=float(getattr(source, "xi_KK")),
        alpha_s=float(getattr(source, "alpha_s")),
        g_s=float(getattr(source, "g_s")),
        left_overlap=_matrix(source, "left_overlap").copy(),
        right_up_overlap=_matrix(source, "right_up_overlap").copy(),
        right_down_overlap=_matrix(source, "right_down_overlap").copy(),
        left_up=_matrix(source, "left_up").copy(),
        left_down=left_down,
        right_up=_matrix(source, "right_up").copy(),
        right_down=right_down,
    )


def bdgamma_wilsons_from_couplings(
    couplings: QuarkMassBasisCouplings,
    *,
    m_kk_gev: float | None = None,
    inputs: BsgammaSMInputs | None = None,
) -> BsgammaWilsonCoefficients:
    """Return the v1 ``b -> d gamma`` C7 proxy using the shared LL machinery.

    The underlying core only exposes a ``b -> s`` entry point.  This adapter
    selects the down-sector ``(0, 2)`` d-b slot and maps it into the core's
    s-b slot before invoking the unchanged C7/C8 matching and running code.
    """

    return _compute_bsgamma_wilsons(
        _bd_source_as_bs_source(couplings),
        m_kk_gev=m_kk_gev,
        inputs=inputs,
    )


def bsgamma_default_sm_inputs() -> BsgammaSMInputs:
    """Return the default C7-normalized ``b -> s gamma`` input bundle."""

    return _default_sm_inputs()


def bsgamma_branching_fraction_from_c7(
    *,
    c7_np: complex = 0.0j,
    c7p_np: complex = 0.0j,
    sm_branching_fraction: float,
    inputs: BsgammaSMInputs | None = None,
) -> BsgammaBranchingResult:
    """Evaluate inclusive ``Bbar -> X_s gamma`` from explicit C7 shifts."""

    return _branching_fraction_from_c7(
        c7_np=c7_np,
        c7p_np=c7p_np,
        sm_branching_fraction=sm_branching_fraction,
        inputs=inputs,
    )


def bsgamma_wilsons_from_couplings(
    couplings: QuarkMassBasisCouplings,
    *,
    m_kk_gev: float | None = None,
    inputs: BsgammaSMInputs | None = None,
) -> BsgammaWilsonCoefficients:
    """Return the v1 ``b -> s gamma`` C7 proxy for mass-basis couplings."""

    return _compute_bsgamma_wilsons(
        couplings,
        m_kk_gev=m_kk_gev,
        inputs=inputs,
    )


def inclusive_bsgamma_from_couplings(
    couplings: QuarkMassBasisCouplings,
    *,
    sm_branching_fraction: float,
    m_kk_gev: float | None = None,
    inputs: BsgammaSMInputs | None = None,
) -> BsgammaBranchingResult:
    """Evaluate ``BR(Bbar -> X_s gamma)`` from mass-basis couplings."""

    return _evaluate_inclusive_bsgamma(
        couplings,
        sm_branching_fraction=sm_branching_fraction,
        m_kk_gev=m_kk_gev,
        inputs=inputs,
    )


def inclusive_bsgamma_sm_branching_fraction(
    *,
    sm_branching_fraction: float,
    inputs: BsgammaSMInputs | None = None,
) -> BsgammaBranchingResult:
    """Evaluate the SM-limit inclusive ``Bbar -> X_s gamma`` branching fraction."""

    return _evaluate_inclusive_bsgamma(
        None,
        sm_branching_fraction=sm_branching_fraction,
        inputs=inputs,
    )


def _with_exclusive_bdrhogamma_metadata(
    result: BsgammaBranchingResult,
    *,
    normalization_label: str,
    normalization_source: str | None,
    mode_label: str,
) -> BsgammaBranchingResult:
    ckm = bdgamma_ckm_factors()
    diagnostics = dict(result.diagnostics)
    diagnostics.update(
        {
            "exclusive_mode": str(mode_label),
            "flavor_transition": "b -> d gamma",
            "exclusive_normalization_label": str(normalization_label),
            "exclusive_normalization_source": (
                "" if normalization_source is None else str(normalization_source)
            ),
            "exclusive_normalization_branching_fraction": float(
                result.sm_branching_fraction
            ),
            "exclusive_branching_formula": (
                f"BR({mode_label}) = BR_norm^(b->d) * "
                "(|C7_SM + C7_NP^(d)|^2 + |C7p_SM + C7p_NP^(d)|^2) / "
                "(|C7_SM|^2 + |C7p_SM|^2)"
            ),
            "exclusive_normalization_note": (
                "The caller-supplied B014 exclusive normalization is the "
                "CKM-suppressed b -> d branching-fraction scale; the adapter "
                "records |V_td/V_ts|^2 for comparison with b -> s gamma but "
                "does not reinterpret a b -> s normalization as b -> d."
            ),
            "bdgamma_adapter_note": (
                "The d-b down-sector coupling entry (0, 2) is mapped into the "
                "existing bsgamma s-b slot before invoking the unchanged C7/C8 "
                "leading-log running machinery."
            ),
            "ckm_input_bundle": ckm.input_bundle,
            "lambda_t_bd": complex(ckm.lambda_t_bd),
            "lambda_t_bs": complex(ckm.lambda_t_bs),
            "abs_vtd_over_vts": float(ckm.abs_vtd_over_vts),
            "abs_lambda_t_bd_over_lambda_t_bs": float(
                ckm.abs_lambda_t_bd_over_lambda_t_bs
            ),
            "ckm_power_suppression_vtd_over_vts_squared": float(
                ckm.ckm_power_suppression
            ),
            "ckm_normalization": (
                "b -> d gamma rates are normalized with lambda_t^d = "
                "V_td^* V_tb and are |V_td/V_ts|^2-suppressed relative to "
                "the b -> s gamma normalization at fixed hadronic response."
            ),
        }
    )
    if result.wilsons is not None:
        diagnostics.update(
            {
                "left_bd_coupling": complex(result.wilsons.left_bs_coupling),
                "right_bd_coupling": complex(result.wilsons.right_bs_coupling),
                "left_bd_overlap": complex(result.wilsons.left_bs_overlap),
                "right_bd_overlap": complex(result.wilsons.right_bs_overlap),
                "bsgamma_wilson_field_alias": (
                    "BsgammaWilsonCoefficients left_bs/right_bs fields carry "
                    "the adapter-selected b-d coupling for this result."
                ),
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
        wilsons=result.wilsons,
        diagnostics=diagnostics,
    )


def exclusive_bdrhogamma_from_c7(
    *,
    c7_np: complex = 0.0j,
    c7p_np: complex = 0.0j,
    exclusive_sm_branching_fraction: float,
    inputs: BsgammaSMInputs | None = None,
    normalization_label: str = "B -> rho/omega gamma exclusive normalization",
    normalization_source: str | None = None,
    mode_label: str = "B -> rho/omega gamma",
) -> BsgammaBranchingResult:
    """Evaluate exclusive ``b -> d gamma`` from explicit low-scale C7 shifts."""

    result = _branching_fraction_from_c7(
        c7_np=c7_np,
        c7p_np=c7p_np,
        sm_branching_fraction=exclusive_sm_branching_fraction,
        inputs=inputs,
    )
    return _with_exclusive_bdrhogamma_metadata(
        result,
        normalization_label=normalization_label,
        normalization_source=normalization_source,
        mode_label=mode_label,
    )


def exclusive_bdrhogamma_from_couplings(
    couplings: QuarkMassBasisCouplings,
    *,
    exclusive_sm_branching_fraction: float,
    m_kk_gev: float | None = None,
    inputs: BsgammaSMInputs | None = None,
    normalization_label: str = "B -> rho/omega gamma exclusive normalization",
    normalization_source: str | None = None,
    mode_label: str = "B -> rho/omega gamma",
) -> BsgammaBranchingResult:
    """Evaluate exclusive ``b -> d gamma`` from mass-basis d-b couplings."""

    p = _default_sm_inputs() if inputs is None else inputs
    wilsons = bdgamma_wilsons_from_couplings(
        couplings,
        m_kk_gev=m_kk_gev,
        inputs=p,
    )
    result = _evaluate_inclusive_bsgamma(
        wilsons,
        sm_branching_fraction=exclusive_sm_branching_fraction,
        inputs=p,
    )
    return _with_exclusive_bdrhogamma_metadata(
        result,
        normalization_label=normalization_label,
        normalization_source=normalization_source,
        mode_label=mode_label,
    )


def exclusive_bdrhogamma_sm_branching_fraction(
    *,
    exclusive_sm_branching_fraction: float,
    inputs: BsgammaSMInputs | None = None,
    normalization_label: str = "B -> rho/omega gamma exclusive normalization",
    normalization_source: str | None = None,
    mode_label: str = "B -> rho/omega gamma",
) -> BsgammaBranchingResult:
    """Evaluate the no-NP exclusive ``b -> d gamma`` normalization."""

    result = _evaluate_inclusive_bsgamma(
        None,
        sm_branching_fraction=exclusive_sm_branching_fraction,
        inputs=inputs,
    )
    return _with_exclusive_bdrhogamma_metadata(
        result,
        normalization_label=normalization_label,
        normalization_source=normalization_source,
        mode_label=mode_label,
    )


def _with_exclusive_btokstargamma_metadata(
    result: BsgammaBranchingResult,
    *,
    normalization_label: str,
    normalization_source: str | None,
) -> BsgammaBranchingResult:
    diagnostics = dict(result.diagnostics)
    diagnostics.update(
        {
            "exclusive_mode": "B -> K*(892) gamma",
            "exclusive_normalization_label": str(normalization_label),
            "exclusive_normalization_source": (
                "" if normalization_source is None else str(normalization_source)
            ),
            "exclusive_normalization_branching_fraction": float(
                result.sm_branching_fraction
            ),
            "exclusive_branching_formula": (
                "BR(B -> K* gamma) = BR_norm * "
                "(|C7_SM + C7_NP|^2 + |C7p_SM + C7p_NP|^2) / "
                "(|C7_SM|^2 + |C7p_SM|^2)"
            ),
            "exclusive_normalization_note": (
                "The exclusive form-factor and normalization dependence is "
                "absorbed into the caller-supplied BR_norm; C7/C8 running is "
                "the shared b -> s gamma machinery."
            ),
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
        wilsons=result.wilsons,
        diagnostics=diagnostics,
    )


def exclusive_btokstargamma_from_c7(
    *,
    c7_np: complex = 0.0j,
    c7p_np: complex = 0.0j,
    exclusive_sm_branching_fraction: float,
    inputs: BsgammaSMInputs | None = None,
    normalization_label: str = "B -> K*(892) gamma exclusive normalization",
    normalization_source: str | None = None,
) -> BsgammaBranchingResult:
    """Evaluate exclusive ``B -> K*(892) gamma`` from explicit C7 shifts."""

    result = _branching_fraction_from_c7(
        c7_np=c7_np,
        c7p_np=c7p_np,
        sm_branching_fraction=exclusive_sm_branching_fraction,
        inputs=inputs,
    )
    return _with_exclusive_btokstargamma_metadata(
        result,
        normalization_label=normalization_label,
        normalization_source=normalization_source,
    )


def exclusive_btokstargamma_from_couplings(
    couplings: QuarkMassBasisCouplings,
    *,
    exclusive_sm_branching_fraction: float,
    m_kk_gev: float | None = None,
    inputs: BsgammaSMInputs | None = None,
    normalization_label: str = "B -> K*(892) gamma exclusive normalization",
    normalization_source: str | None = None,
) -> BsgammaBranchingResult:
    """Evaluate exclusive ``B -> K*(892) gamma`` from mass-basis couplings."""

    result = _evaluate_inclusive_bsgamma(
        couplings,
        sm_branching_fraction=exclusive_sm_branching_fraction,
        m_kk_gev=m_kk_gev,
        inputs=inputs,
    )
    return _with_exclusive_btokstargamma_metadata(
        result,
        normalization_label=normalization_label,
        normalization_source=normalization_source,
    )


def exclusive_btokstargamma_sm_branching_fraction(
    *,
    exclusive_sm_branching_fraction: float,
    inputs: BsgammaSMInputs | None = None,
    normalization_label: str = "B -> K*(892) gamma exclusive normalization",
    normalization_source: str | None = None,
) -> BsgammaBranchingResult:
    """Evaluate the SM-limit exclusive ``B -> K*(892) gamma`` branching fraction."""

    result = _evaluate_inclusive_bsgamma(
        None,
        sm_branching_fraction=exclusive_sm_branching_fraction,
        inputs=inputs,
    )
    return _with_exclusive_btokstargamma_metadata(
        result,
        normalization_label=normalization_label,
        normalization_source=normalization_source,
    )


def _with_exclusive_bsphigamma_metadata(
    result: BsgammaBranchingResult,
    *,
    normalization_label: str,
    normalization_source: str | None,
) -> BsgammaBranchingResult:
    diagnostics = dict(result.diagnostics)
    diagnostics.update(
        {
            "exclusive_mode": "B_s0 -> phi(1020) gamma",
            "exclusive_normalization_label": str(normalization_label),
            "exclusive_normalization_source": (
                "" if normalization_source is None else str(normalization_source)
            ),
            "exclusive_normalization_branching_fraction": float(
                result.sm_branching_fraction
            ),
            "exclusive_branching_formula": (
                "BR(B_s -> phi gamma) = BR_norm * "
                "(|C7_SM + C7_NP|^2 + |C7p_SM + C7p_NP|^2) / "
                "(|C7_SM|^2 + |C7p_SM|^2)"
            ),
            "exclusive_normalization_note": (
                "The B_s -> phi gamma exclusive form-factor and normalization "
                "dependence is absorbed into the caller-supplied BR_norm; "
                "C7/C8 running is the shared b -> s gamma machinery."
            ),
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
        wilsons=result.wilsons,
        diagnostics=diagnostics,
    )


def exclusive_bsphigamma_from_c7(
    *,
    c7_np: complex = 0.0j,
    c7p_np: complex = 0.0j,
    exclusive_sm_branching_fraction: float,
    inputs: BsgammaSMInputs | None = None,
    normalization_label: str = "B_s0 -> phi(1020) gamma exclusive normalization",
    normalization_source: str | None = None,
) -> BsgammaBranchingResult:
    """Evaluate exclusive ``B_s0 -> phi(1020) gamma`` from explicit C7 shifts."""

    result = _branching_fraction_from_c7(
        c7_np=c7_np,
        c7p_np=c7p_np,
        sm_branching_fraction=exclusive_sm_branching_fraction,
        inputs=inputs,
    )
    return _with_exclusive_bsphigamma_metadata(
        result,
        normalization_label=normalization_label,
        normalization_source=normalization_source,
    )


def exclusive_bsphigamma_from_couplings(
    couplings: QuarkMassBasisCouplings,
    *,
    exclusive_sm_branching_fraction: float,
    m_kk_gev: float | None = None,
    inputs: BsgammaSMInputs | None = None,
    normalization_label: str = "B_s0 -> phi(1020) gamma exclusive normalization",
    normalization_source: str | None = None,
) -> BsgammaBranchingResult:
    """Evaluate exclusive ``B_s0 -> phi(1020) gamma`` from mass-basis couplings."""

    result = _evaluate_inclusive_bsgamma(
        couplings,
        sm_branching_fraction=exclusive_sm_branching_fraction,
        m_kk_gev=m_kk_gev,
        inputs=inputs,
    )
    return _with_exclusive_bsphigamma_metadata(
        result,
        normalization_label=normalization_label,
        normalization_source=normalization_source,
    )


def exclusive_bsphigamma_sm_branching_fraction(
    *,
    exclusive_sm_branching_fraction: float,
    inputs: BsgammaSMInputs | None = None,
    normalization_label: str = "B_s0 -> phi(1020) gamma exclusive normalization",
    normalization_source: str | None = None,
) -> BsgammaBranchingResult:
    """Evaluate the SM-limit exclusive ``B_s0 -> phi(1020) gamma`` rate."""

    result = _evaluate_inclusive_bsgamma(
        None,
        sm_branching_fraction=exclusive_sm_branching_fraction,
        inputs=inputs,
    )
    return _with_exclusive_bsphigamma_metadata(
        result,
        normalization_label=normalization_label,
        normalization_source=normalization_source,
    )
