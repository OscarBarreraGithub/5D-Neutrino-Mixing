"""Tests for quark-sector validation and plot-data helpers."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from quarkConstraints.validation import (
    benchmark_fit_summary,
    benchmark_plot_data,
    benchmark_solution,
    r_sweep_plot_data,
    r_sweep_trends_ok,
)


def test_benchmark_fit_summary_passes_for_default_solution():
    """The benchmark fit should satisfy the validation gates."""
    solution = benchmark_solution(max_nfev=120)
    summary = benchmark_fit_summary(solution)

    assert summary.passes_mass
    assert summary.passes_ckm
    assert summary.passes_proxy
    assert not summary.passes_paper_proxy
    assert summary.passes_misalignment
    assert summary.passes_deltaf2
    assert summary.passes_all
    assert summary.down_proxy < summary.proxy_limit
    assert summary.down_proxy > summary.paper_proxy_target
    assert summary.down_to_up_misalignment_ratio < summary.misalignment_limit


def test_benchmark_plot_data_has_expected_shapes():
    """Benchmark plot data should expose spectra and residual scalars."""
    data = benchmark_plot_data(max_nfev=120)

    assert data["target_masses_up"].shape == (3,)
    assert data["target_masses_down"].shape == (3,)
    assert data["fit_masses_up"].shape == (3,)
    assert data["fit_masses_down"].shape == (3,)
    assert data["c_Q"].shape == (3,)
    assert data["c_u"].shape == (3,)
    assert data["c_d"].shape == (3,)
    assert data["F_Q"].shape == (3,)
    assert data["xi_KK"].shape == (1,)
    assert data["M_KK"].shape == (1,)
    assert data["mass_residual_norm"].shape == (1,)
    assert data["ckm_residual_norm"].shape == (1,)
    assert data["benchmark_proxy_limit"].shape == (1,)
    assert data["down_misalignment"].shape == (1,)
    assert data["up_misalignment"].shape == (1,)
    assert data["down_to_up_misalignment_ratio"].shape == (1,)
    assert data["alignment_ratio"].shape == (1,)
    assert data["deltaf2_max_ratio"].shape == (1,)
    assert data["epsilon_k_ratio"].shape == (1,)
    assert data["b_d_ratio"].shape == (1,)
    assert data["b_s_ratio"].shape == (1,)
    assert data["d_ratio"].shape == (1,)


def test_r_sweep_plot_data_shows_expected_suppression_trend():
    """The r sweep should show smaller down-sector proxy and alignment at small r."""
    data = r_sweep_plot_data([0.05, 0.1, 0.25, 0.4, 1.0], max_nfev=100)

    assert np.all(np.diff(data["r_values"]) > 0.0)
    assert data["c_Q"].shape == (5, 3)
    assert data["c_u"].shape == (5, 3)
    assert data["c_d"].shape == (5, 3)
    assert data["down_misalignment"].shape == (5,)
    assert data["up_misalignment"].shape == (5,)
    assert data["down_to_up_misalignment_ratio"].shape == (5,)
    assert data["alignment_ratio"].shape == (5,)
    assert data["epsilon_k_ratio"].shape == (5,)
    assert data["b_d_ratio"].shape == (5,)
    assert data["b_s_ratio"].shape == (5,)
    assert data["d_ratio"].shape == (5,)
    assert r_sweep_trends_ok(data)
