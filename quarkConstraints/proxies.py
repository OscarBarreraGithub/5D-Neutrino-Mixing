"""Fast flavor proxies and matrix-level alignment diagnostics."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np

from .benchmarks import benchmark_spurion_input, default_quark_targets, default_spurion_seed
from .fit import QuarkFitResult, evaluate_quark_fit, fit_quark_sector
from .model import QuarkBulkState, QuarkSpurionPoint, derive_bulk_state
from .scales import DEFAULT_QUARK_XI_KK, default_quark_m_kk_from_lambda_ir


@dataclass(frozen=True)
class AlignmentDiagnostics:
    """Matrix-level quark-basis diagnostics for one quark-sector point.

    The ``*_offdiag_ratio_in_q_basis`` entries are off-diagonal Frobenius
    fractions, so smaller values mean stronger alignment. The older
    ``alignment`` naming is kept for compatibility, but these are best read as
    misalignment fractions.
    """

    down_commutator_norm: float
    up_commutator_norm: float
    down_offdiag_ratio_in_q_basis: float
    up_offdiag_ratio_in_q_basis: float

    @property
    def down_misalignment_in_q_basis(self) -> float:
        return self.down_offdiag_ratio_in_q_basis

    @property
    def up_misalignment_in_q_basis(self) -> float:
        return self.up_offdiag_ratio_in_q_basis

    @property
    def down_to_up_misalignment_ratio(self) -> float:
        return float(self.down_offdiag_ratio_in_q_basis / max(self.up_offdiag_ratio_in_q_basis, 1e-30))


@dataclass(frozen=True)
class ProxySummary:
    """Fast proxy summary for one quark-sector point."""

    h_rs_proxy: float
    f_q3: float
    m_kk: float
    q_overlap_hierarchy: np.ndarray
    r_value: float
    diagnostics: AlignmentDiagnostics


def _commutator_norm(a: np.ndarray, b: np.ndarray) -> float:
    """Return a normalized Frobenius norm of the commutator."""
    commutator = a @ b - b @ a
    denom = np.linalg.norm(a, ord="fro") * np.linalg.norm(b, ord="fro")
    if denom <= 0.0:
        return 0.0
    return float(np.linalg.norm(commutator, ord="fro") / denom)


def _offdiag_ratio(matrix: np.ndarray) -> float:
    """Return the off-diagonal Frobenius fraction of a matrix.

    Smaller values indicate stronger alignment with the chosen basis.
    """
    arr = np.asarray(matrix, dtype=np.complex128).copy()
    np.fill_diagonal(arr, 0.0)
    denom = np.linalg.norm(matrix, ord="fro")
    if denom <= 0.0:
        return 0.0
    return float(np.linalg.norm(arr, ord="fro") / denom)


def compute_alignment_diagnostics(
    point: QuarkSpurionPoint,
    *,
    bulk_state: QuarkBulkState | None = None,
) -> AlignmentDiagnostics:
    """Compute matrix-level alignment diagnostics for one MFV point."""
    state = derive_bulk_state(point) if bulk_state is None else bulk_state
    down_left = point.Y_d @ point.Y_d.conjugate().T
    up_left = point.Y_u @ point.Y_u.conjugate().T
    down_in_q_basis = state.rotation_Q.conjugate().T @ down_left @ state.rotation_Q
    up_in_q_basis = state.rotation_Q.conjugate().T @ up_left @ state.rotation_Q
    return AlignmentDiagnostics(
        down_commutator_norm=_commutator_norm(state.C_Q, down_left),
        up_commutator_norm=_commutator_norm(state.C_Q, up_left),
        down_offdiag_ratio_in_q_basis=_offdiag_ratio(down_in_q_basis),
        up_offdiag_ratio_in_q_basis=_offdiag_ratio(up_in_q_basis),
    )


def compute_proxy_summary(
    fit_result: QuarkFitResult,
    *,
    m_kk: float | None = None,
    xi_KK: float = DEFAULT_QUARK_XI_KK,
) -> ProxySummary:
    """Compute the paper-inspired down-sector proxy summary.

    The `h_RS` quantity here is a proxy, not a full observable prediction.
    When `m_kk` is omitted, this helper follows the repo convention
    `M_KK ≡ Lambda_IR` for internal comparisons. That is a bookkeeping choice,
    not a literal identification of the physical first KK mass, which remains
    sector dependent in the broader repo conventions.
    """
    state = fit_result.bulk_state
    kk_mass = float(
        m_kk
        if m_kk is not None
        else default_quark_m_kk_from_lambda_ir(state.point.Lambda_IR, xi_KK=xi_KK)
    )
    f_q3 = float(np.max(state.F_Q))
    q_overlap_hierarchy = np.asarray(state.F_Q / max(f_q3, 1e-30), dtype=float)
    h_rs_proxy = 0.5 * (3000.0 / kk_mass) ** 2 * (f_q3 / 0.3) ** 4
    diagnostics = compute_alignment_diagnostics(state.point, bulk_state=state)
    return ProxySummary(
        h_rs_proxy=float(h_rs_proxy),
        f_q3=f_q3,
        m_kk=kk_mass,
        q_overlap_hierarchy=q_overlap_hierarchy,
        r_value=float(state.point.r),
        diagnostics=diagnostics,
    )


def summarize_flavor_diagnostics(
    fit_result: QuarkFitResult,
    *,
    m_kk: float | None = None,
    xi_KK: float = DEFAULT_QUARK_XI_KK,
) -> ProxySummary:
    """Backward-compatible alias used by scripts and validation helpers."""
    return compute_proxy_summary(fit_result, m_kk=m_kk, xi_KK=xi_KK)


def alignment_diagnostics(fit_result: QuarkFitResult) -> AlignmentDiagnostics:
    """Compatibility wrapper returning the alignment diagnostics for one fit."""
    return compute_alignment_diagnostics(fit_result.point, bulk_state=fit_result.bulk_state)


def h_rs_proxy(
    fit_result: QuarkFitResult,
    *,
    m_kk: float | None = None,
    xi_KK: float = DEFAULT_QUARK_XI_KK,
) -> float:
    """Compatibility wrapper returning the paper-style proxy value only."""
    return compute_proxy_summary(fit_result, m_kk=m_kk, xi_KK=xi_KK).h_rs_proxy


def suppression_summary(
    fit_result: QuarkFitResult,
    *,
    m_kk: float | None = None,
    xi_KK: float = DEFAULT_QUARK_XI_KK,
    include_legacy_alignment_aliases: bool = False,
) -> dict[str, float]:
    """Compatibility wrapper returning a flat summary dict.

    By default this helper returns only correctly named misalignment quantities.
    Set ``include_legacy_alignment_aliases=True`` to also emit the old
    ``alignment`` keys as aliases for those same misalignment fractions.
    """
    summary = compute_proxy_summary(fit_result, m_kk=m_kk, xi_KK=xi_KK)
    down_misalignment = float(summary.diagnostics.down_misalignment_in_q_basis)
    up_misalignment = float(summary.diagnostics.up_misalignment_in_q_basis)
    misalignment_ratio = float(summary.diagnostics.down_to_up_misalignment_ratio)
    out = {
        "h_rs_proxy": float(summary.h_rs_proxy),
        "down_misalignment": down_misalignment,
        "up_misalignment": up_misalignment,
        "down_to_up_misalignment_ratio": misalignment_ratio,
        "f_q_hierarchy": float(np.min(summary.q_overlap_hierarchy)),
        "M_KK": float(summary.m_kk),
        "r": float(summary.r_value),
    }
    if include_legacy_alignment_aliases:
        out.update(
            {
                "down_alignment": down_misalignment,
                "up_alignment": up_misalignment,
                "down_to_up_alignment_ratio": misalignment_ratio,
            }
        )
    return out


def sweep_r_proxy_summary(
    r_values: Sequence[float],
    *,
    Lambda_IR: float = 3000.0,
    refit: bool = True,
    xi_KK: float = DEFAULT_QUARK_XI_KK,
) -> list[ProxySummary]:
    """Evaluate the benchmark while varying the down-sector protection dial."""
    targets = default_quark_targets()
    seed = default_spurion_seed()
    summaries: list[ProxySummary] = []
    for r_value in np.asarray(r_values, dtype=float):
        if refit:
            solution = fit_quark_sector(
                targets,
                r=float(r_value),
                Lambda_IR=Lambda_IR,
                seed=seed,
                overall_scale=seed.overall_scale,
                max_nfev=120,
            )
            seed = solution.seed
            fit = solution.result
        else:
            fit = fit_quark_sector(
                targets,
                r=float(r_value),
                Lambda_IR=Lambda_IR,
                seed=seed,
                overall_scale=seed.overall_scale,
                max_nfev=1,
                fit_orientation=False,
            ).result
        summaries.append(compute_proxy_summary(fit, xi_KK=xi_KK))
    return summaries


def sweep_r(
    point: QuarkSpurionPoint,
    r_values: Sequence[float],
) -> list[dict[str, float]]:
    """Compatibility wrapper returning flat rows over an r sweep.

    This sweep keeps the supplied spurion matrices fixed and varies only ``r``
    and the derived bulk-state observables.
    """
    rows: list[dict[str, float]] = []
    for r_value in np.asarray(r_values, dtype=float):
        swept_point = QuarkSpurionPoint(
            Y_u=point.Y_u,
            Y_d=point.Y_d,
            r=float(r_value),
            Lambda_IR=point.Lambda_IR,
            k=point.k,
            v=point.v,
            label=f"{point.label}-r-{r_value:.6g}",
            metadata=dict(point.metadata),
        )
        fit = evaluate_quark_fit(swept_point)
        rows.append(suppression_summary(fit))
    return rows
