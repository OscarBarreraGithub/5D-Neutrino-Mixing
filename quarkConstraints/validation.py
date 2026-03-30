"""Reusable validation data for plots, notebooks, and tests."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Sequence

import numpy as np

from .benchmarks import QuarkTargets, default_quark_targets, default_spurion_seed
from .couplings import compute_quark_kk_gluon_couplings
from .deltaf2 import evaluate_delta_f2_constraints
from .fit import QuarkFitSolution, fit_quark_sector
from .model import BulkMassMap
from .proxies import summarize_flavor_diagnostics
from .scales import DEFAULT_QUARK_XI_KK, default_quark_m_kk_from_lambda_ir


@dataclass(frozen=True)
class BenchmarkFitSummary:
    """Compact benchmark summary for plot and notebook validation."""

    mass_residual_norm: float
    ckm_residual_norm: float
    m_kk: float
    down_proxy: float
    alignment_metric: float
    up_alignment_metric: float
    deltaf2_max_ratio: float
    epsilon_k_ratio: float
    b_d_ratio: float
    b_s_ratio: float
    d_ratio: float
    passes_mass: bool
    passes_ckm: bool
    passes_proxy: bool
    passes_alignment: bool
    passes_deltaf2: bool

    @property
    def passes_all(self) -> bool:
        return (
            self.passes_mass
            and self.passes_ckm
            and self.passes_proxy
            and self.passes_alignment
            and self.passes_deltaf2
        )


def benchmark_solution(
    *,
    r: float = 0.25,
    overall_scale: float | None = None,
    targets: QuarkTargets | None = None,
    bulk_map: BulkMassMap | None = None,
    max_nfev: int = 150,
) -> QuarkFitSolution:
    """Return a deterministic benchmark fit solution."""
    targets = default_quark_targets() if targets is None else targets
    seed = default_spurion_seed()
    scale = seed.overall_scale if overall_scale is None else overall_scale
    return fit_quark_sector(
        targets,
        r=r,
        overall_scale=scale,
        seed=seed,
        bulk_map=bulk_map,
        max_nfev=max_nfev,
    )


def benchmark_fit_summary(
    solution: QuarkFitSolution,
    *,
    xi_KK: float = DEFAULT_QUARK_XI_KK,
) -> BenchmarkFitSummary:
    """Return the benchmark summary used by validation plots."""
    result = solution.result
    m_kk = default_quark_m_kk_from_lambda_ir(result.point.Lambda_IR, xi_KK=xi_KK)
    diagnostics = summarize_flavor_diagnostics(result, m_kk=m_kk)
    deltaf2 = evaluate_delta_f2_constraints(
        compute_quark_kk_gluon_couplings(result, M_KK=m_kk, xi_KK=xi_KK),
        M_KK=m_kk,
    )
    deltaf2_by_system = deltaf2.by_system
    mass_residual_norm = float(
        np.linalg.norm(np.concatenate([result.mass_residuals_up, result.mass_residuals_down]))
    )
    ckm_residual_norm = float(np.linalg.norm(result.ckm_residuals))
    down_alignment = float(diagnostics.diagnostics.down_offdiag_ratio_in_q_basis)
    up_alignment = float(diagnostics.diagnostics.up_offdiag_ratio_in_q_basis)
    alignment_ratio = float(down_alignment / max(up_alignment, 1e-30))
    return BenchmarkFitSummary(
        mass_residual_norm=mass_residual_norm,
        ckm_residual_norm=ckm_residual_norm,
        m_kk=m_kk,
        down_proxy=float(diagnostics.h_rs_proxy),
        alignment_metric=alignment_ratio,
        up_alignment_metric=up_alignment,
        deltaf2_max_ratio=deltaf2.max_ratio_to_bound,
        epsilon_k_ratio=deltaf2_by_system["K"].ratio_to_bound,
        b_d_ratio=deltaf2_by_system["B_d"].ratio_to_bound,
        b_s_ratio=deltaf2_by_system["B_s"].ratio_to_bound,
        d_ratio=deltaf2_by_system["D"].ratio_to_bound,
        passes_mass=mass_residual_norm < 5.0e-3,
        passes_ckm=ckm_residual_norm < 5.0e-3,
        passes_proxy=diagnostics.h_rs_proxy < 1.0,
        passes_alignment=alignment_ratio < 6.0,
        passes_deltaf2=deltaf2.passes_all,
    )


def benchmark_plot_data(
    *,
    r: float = 0.25,
    overall_scale: float | None = None,
    max_nfev: int = 150,
    xi_KK: float = DEFAULT_QUARK_XI_KK,
) -> Dict[str, np.ndarray]:
    """Return arrays used by the benchmark notebook panels."""
    solution = benchmark_solution(r=r, overall_scale=overall_scale, max_nfev=max_nfev)
    result = solution.result
    m_kk = default_quark_m_kk_from_lambda_ir(result.point.Lambda_IR, xi_KK=xi_KK)
    diagnostics = summarize_flavor_diagnostics(result, m_kk=m_kk)
    deltaf2 = evaluate_delta_f2_constraints(
        compute_quark_kk_gluon_couplings(result, M_KK=m_kk, xi_KK=xi_KK),
        M_KK=m_kk,
    )
    deltaf2_by_system = deltaf2.by_system
    targets = default_quark_targets()
    return {
        "target_masses_up": targets.up_masses,
        "target_masses_down": targets.down_masses,
        "fit_masses_up": result.masses_up,
        "fit_masses_down": result.masses_down,
        "target_ckm_observables": targets.ckm_observables,
        "fit_ckm_observables": result.ckm_observables,
        "c_Q": result.state.c_Q,
        "c_u": result.state.c_u,
        "c_d": result.state.c_d,
        "F_Q": result.state.F_Q,
        "F_u": result.state.F_u,
        "F_d": result.state.F_d,
        "M_KK": np.array([m_kk], dtype=float),
        "mass_residual_norm": np.array(
            [np.linalg.norm(np.concatenate([result.mass_residuals_up, result.mass_residuals_down]))],
            dtype=float,
        ),
        "ckm_residual_norm": np.array([np.linalg.norm(result.ckm_residuals)], dtype=float),
        "h_rs_proxy": np.array([diagnostics.h_rs_proxy], dtype=float),
        "down_alignment": np.array(
            [diagnostics.diagnostics.down_offdiag_ratio_in_q_basis],
            dtype=float,
        ),
        "up_alignment": np.array(
            [diagnostics.diagnostics.up_offdiag_ratio_in_q_basis],
            dtype=float,
        ),
        "alignment_ratio": np.array(
            [
                diagnostics.diagnostics.down_offdiag_ratio_in_q_basis
                / max(diagnostics.diagnostics.up_offdiag_ratio_in_q_basis, 1e-30)
            ],
            dtype=float,
        ),
        "deltaf2_max_ratio": np.array([deltaf2.max_ratio_to_bound], dtype=float),
        "epsilon_k_ratio": np.array([deltaf2_by_system["K"].ratio_to_bound], dtype=float),
        "b_d_ratio": np.array([deltaf2_by_system["B_d"].ratio_to_bound], dtype=float),
        "b_s_ratio": np.array([deltaf2_by_system["B_s"].ratio_to_bound], dtype=float),
        "d_ratio": np.array([deltaf2_by_system["D"].ratio_to_bound], dtype=float),
    }


def r_sweep_plot_data(
    r_values: Sequence[float],
    *,
    overall_scale: float | None = None,
    targets: QuarkTargets | None = None,
    bulk_map: BulkMassMap | None = None,
    max_nfev: int = 120,
    xi_KK: float = DEFAULT_QUARK_XI_KK,
) -> Dict[str, np.ndarray]:
    """Return `r`-sweep arrays for alignment and proxy validation."""
    targets = default_quark_targets() if targets is None else targets
    seed = default_spurion_seed()
    scale = seed.overall_scale if overall_scale is None else overall_scale
    r_arr = np.sort(np.asarray(r_values, dtype=float))
    if r_arr.ndim != 1 or r_arr.size == 0:
        raise ValueError("r_values must be a non-empty 1D array")

    scores = []
    down_proxy = []
    alignment_metric = []
    up_alignment = []
    alignment_ratio = []
    epsilon_k_ratio = []
    b_d_ratio = []
    b_s_ratio = []
    d_ratio = []
    c_q = []
    c_u = []
    c_d = []
    f_q3 = []
    current_seed = seed

    for r_value in r_arr:
        solution = fit_quark_sector(
            targets,
            r=float(r_value),
            overall_scale=scale,
            seed=current_seed,
            bulk_map=bulk_map,
            max_nfev=max_nfev,
        )
        current_seed = solution.seed
        result = solution.result
        m_kk = default_quark_m_kk_from_lambda_ir(result.point.Lambda_IR, xi_KK=xi_KK)
        diagnostics = summarize_flavor_diagnostics(result, m_kk=m_kk)
        deltaf2 = evaluate_delta_f2_constraints(
            compute_quark_kk_gluon_couplings(result, M_KK=m_kk, xi_KK=xi_KK),
            M_KK=m_kk,
        )
        deltaf2_by_system = deltaf2.by_system
        scores.append(result.score)
        down_proxy.append(diagnostics.h_rs_proxy)
        alignment_metric.append(diagnostics.diagnostics.down_offdiag_ratio_in_q_basis)
        up_alignment.append(diagnostics.diagnostics.up_offdiag_ratio_in_q_basis)
        alignment_ratio.append(
            diagnostics.diagnostics.down_offdiag_ratio_in_q_basis
            / max(diagnostics.diagnostics.up_offdiag_ratio_in_q_basis, 1e-30)
        )
        epsilon_k_ratio.append(deltaf2_by_system["K"].ratio_to_bound)
        b_d_ratio.append(deltaf2_by_system["B_d"].ratio_to_bound)
        b_s_ratio.append(deltaf2_by_system["B_s"].ratio_to_bound)
        d_ratio.append(deltaf2_by_system["D"].ratio_to_bound)
        c_q.append(result.state.c_Q)
        c_u.append(result.state.c_u)
        c_d.append(result.state.c_d)
        f_q3.append(diagnostics.f_q3)

    return {
        "r_values": r_arr,
        "score": np.asarray(scores, dtype=float),
        "down_proxy": np.asarray(down_proxy, dtype=float),
        "alignment_metric": np.asarray(alignment_metric, dtype=float),
        "up_alignment": np.asarray(up_alignment, dtype=float),
        "alignment_ratio": np.asarray(alignment_ratio, dtype=float),
        "epsilon_k_ratio": np.asarray(epsilon_k_ratio, dtype=float),
        "b_d_ratio": np.asarray(b_d_ratio, dtype=float),
        "b_s_ratio": np.asarray(b_s_ratio, dtype=float),
        "d_ratio": np.asarray(d_ratio, dtype=float),
        "c_Q": np.asarray(c_q, dtype=float),
        "c_u": np.asarray(c_u, dtype=float),
        "c_d": np.asarray(c_d, dtype=float),
        "F_Q3": np.asarray(f_q3, dtype=float),
    }


def is_monotonic_decreasing(values: Sequence[float], *, atol: float = 1e-12) -> bool:
    """Return whether the array is monotonically non-increasing."""
    arr = np.asarray(values, dtype=float)
    return bool(np.all(np.diff(arr) <= atol))


def is_monotonic_increasing(values: Sequence[float], *, atol: float = 1e-12) -> bool:
    """Return whether the array is monotonically non-decreasing."""
    arr = np.asarray(values, dtype=float)
    return bool(np.all(np.diff(arr) >= -atol))


def r_sweep_trends_ok(data: Dict[str, np.ndarray]) -> bool:
    """Return whether the `r`-sweep shows the expected suppression trend."""
    return (
        is_monotonic_increasing(data["down_proxy"])
        and is_monotonic_increasing(data["alignment_ratio"])
        and is_monotonic_increasing(data["epsilon_k_ratio"])
        and is_monotonic_increasing(data["b_d_ratio"])
        and is_monotonic_increasing(data["b_s_ratio"])
    )
