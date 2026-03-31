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
from .scales import (
    DEFAULT_QUARK_BENCHMARK_XI_KK,
    DEFAULT_QUARK_BENCHMARK_H_RS_MAX,
    DEFAULT_QUARK_MAX_MISALIGNMENT_STRESS,
    DEFAULT_QUARK_PAPER_H_RS_MAX,
    default_quark_m_kk_from_lambda_ir,
)


@dataclass(frozen=True)
class BenchmarkFitSummary:
    """Compact benchmark summary for plot and notebook validation."""

    mass_residual_norm: float
    ckm_residual_norm: float
    m_kk: float
    xi_kk: float
    down_proxy: float
    proxy_limit: float
    paper_proxy_target: float
    down_misalignment: float
    up_misalignment: float
    down_to_up_misalignment_ratio: float
    misalignment_limit: float
    deltaf2_max_ratio: float
    epsilon_k_ratio: float
    b_d_ratio: float
    b_s_ratio: float
    d_ratio: float
    passes_mass: bool
    passes_ckm: bool
    passes_proxy: bool
    passes_paper_proxy: bool
    passes_misalignment: bool
    passes_deltaf2: bool

    @property
    def alignment_metric(self) -> float:
        return self.down_misalignment

    @property
    def up_alignment_metric(self) -> float:
        return self.up_misalignment

    @property
    def alignment_ratio(self) -> float:
        return self.down_to_up_misalignment_ratio

    @property
    def passes_alignment(self) -> bool:
        return self.passes_misalignment

    @property
    def passes_all(self) -> bool:
        return (
            self.passes_mass
            and self.passes_ckm
            and self.passes_proxy
            and self.passes_misalignment
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
    xi_KK: float = DEFAULT_QUARK_BENCHMARK_XI_KK,
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
    down_misalignment = float(diagnostics.diagnostics.down_misalignment_in_q_basis)
    up_misalignment = float(diagnostics.diagnostics.up_misalignment_in_q_basis)
    misalignment_ratio = float(diagnostics.diagnostics.down_to_up_misalignment_ratio)
    return BenchmarkFitSummary(
        mass_residual_norm=mass_residual_norm,
        ckm_residual_norm=ckm_residual_norm,
        m_kk=m_kk,
        xi_kk=float(xi_KK),
        down_proxy=float(diagnostics.h_rs_proxy),
        proxy_limit=DEFAULT_QUARK_BENCHMARK_H_RS_MAX,
        paper_proxy_target=DEFAULT_QUARK_PAPER_H_RS_MAX,
        down_misalignment=down_misalignment,
        up_misalignment=up_misalignment,
        down_to_up_misalignment_ratio=misalignment_ratio,
        misalignment_limit=DEFAULT_QUARK_MAX_MISALIGNMENT_STRESS,
        deltaf2_max_ratio=deltaf2.max_ratio_to_bound,
        epsilon_k_ratio=deltaf2_by_system["K"].ratio_to_bound,
        b_d_ratio=deltaf2_by_system["B_d"].ratio_to_bound,
        b_s_ratio=deltaf2_by_system["B_s"].ratio_to_bound,
        d_ratio=deltaf2_by_system["D"].ratio_to_bound,
        passes_mass=mass_residual_norm < 5.0e-3,
        passes_ckm=ckm_residual_norm < 5.0e-3,
        passes_proxy=diagnostics.h_rs_proxy < DEFAULT_QUARK_BENCHMARK_H_RS_MAX,
        passes_paper_proxy=diagnostics.h_rs_proxy < DEFAULT_QUARK_PAPER_H_RS_MAX,
        passes_misalignment=misalignment_ratio < DEFAULT_QUARK_MAX_MISALIGNMENT_STRESS,
        passes_deltaf2=deltaf2.passes_all,
    )


def benchmark_plot_data(
    *,
    r: float = 0.25,
    overall_scale: float | None = None,
    max_nfev: int = 150,
    xi_KK: float = DEFAULT_QUARK_BENCHMARK_XI_KK,
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
        "xi_KK": np.array([xi_KK], dtype=float),
        "M_KK": np.array([m_kk], dtype=float),
        "mass_residual_norm": np.array(
            [np.linalg.norm(np.concatenate([result.mass_residuals_up, result.mass_residuals_down]))],
            dtype=float,
        ),
        "ckm_residual_norm": np.array([np.linalg.norm(result.ckm_residuals)], dtype=float),
        "benchmark_proxy_limit": np.array([DEFAULT_QUARK_BENCHMARK_H_RS_MAX], dtype=float),
        "proxy_limit": np.array([DEFAULT_QUARK_BENCHMARK_H_RS_MAX], dtype=float),
        "paper_proxy_target": np.array([DEFAULT_QUARK_PAPER_H_RS_MAX], dtype=float),
        "h_rs_proxy": np.array([diagnostics.h_rs_proxy], dtype=float),
        "down_misalignment": np.array(
            [diagnostics.diagnostics.down_misalignment_in_q_basis],
            dtype=float,
        ),
        "up_misalignment": np.array(
            [diagnostics.diagnostics.up_misalignment_in_q_basis],
            dtype=float,
        ),
        "down_to_up_misalignment_ratio": np.array(
            [
                diagnostics.diagnostics.down_to_up_misalignment_ratio
            ],
            dtype=float,
        ),
        "down_alignment": np.array(
            [diagnostics.diagnostics.down_misalignment_in_q_basis],
            dtype=float,
        ),
        "up_alignment": np.array(
            [diagnostics.diagnostics.up_misalignment_in_q_basis],
            dtype=float,
        ),
        "alignment_ratio": np.array(
            [diagnostics.diagnostics.down_to_up_misalignment_ratio],
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
    xi_KK: float = DEFAULT_QUARK_BENCHMARK_XI_KK,
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
    down_misalignment = []
    up_misalignment = []
    misalignment_ratio = []
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
        down_misalignment.append(diagnostics.diagnostics.down_misalignment_in_q_basis)
        up_misalignment.append(diagnostics.diagnostics.up_misalignment_in_q_basis)
        misalignment_ratio.append(diagnostics.diagnostics.down_to_up_misalignment_ratio)
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
        "down_misalignment": np.asarray(down_misalignment, dtype=float),
        "up_misalignment": np.asarray(up_misalignment, dtype=float),
        "down_to_up_misalignment_ratio": np.asarray(misalignment_ratio, dtype=float),
        "alignment_metric": np.asarray(down_misalignment, dtype=float),
        "up_alignment": np.asarray(up_misalignment, dtype=float),
        "alignment_ratio": np.asarray(misalignment_ratio, dtype=float),
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
    misalignment_key = (
        "down_to_up_misalignment_ratio"
        if "down_to_up_misalignment_ratio" in data
        else "alignment_ratio"
    )
    return (
        is_monotonic_increasing(data["down_proxy"])
        and is_monotonic_increasing(data[misalignment_key])
        and is_monotonic_increasing(data["epsilon_k_ratio"])
        and is_monotonic_increasing(data["b_d_ratio"])
        and is_monotonic_increasing(data["b_s_ratio"])
    )
