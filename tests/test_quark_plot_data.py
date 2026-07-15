"""Tests for quark-sector validation and plot-data helpers."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from quarkConstraints.validation import (
    benchmark_fit_summary,
    benchmark_plot_data,
    benchmark_solution,
    bulk_mass_map_comparison_data,
    r_sweep_plot_data,
    r_sweep_trends_ok,
)


def test_benchmark_fit_summary_reports_default_solution_gates():
    """The benchmark fit summary should expose the current validation gates."""
    solution = benchmark_solution(max_nfev=120)
    summary = benchmark_fit_summary(solution)

    assert summary.passes_mass
    assert summary.passes_ckm
    assert summary.passes_proxy
    assert summary.passes_misalignment
    # The Phase 2 hadronic-input and Wilson-RG audits make the default point an
    # epsilon_K failure, while the non-DeltaF2 fit-quality gates still pass.
    assert not summary.passes_deltaf2
    assert not summary.passes_all
    # Re-pinned after M-15: validation summaries use the fitted normalized seed
    # representative rather than the old CKM-rebuilt reported seed.
    assert np.isclose(summary.deltaf2_max_ratio, 4.950085995246542)
    assert summary.passes_proxy == (summary.down_proxy < summary.proxy_limit)
    assert summary.passes_paper_proxy == (summary.down_proxy < summary.paper_proxy_target)
    assert summary.proxy_limit >= summary.paper_proxy_target
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
    assert data["paper_proxy_target"].shape == (1,)
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
    """The small-r sweep should show the expected suppression trend."""
    data = r_sweep_plot_data([0.0, 0.05, 0.1, 0.25, 0.4], max_nfev=100)

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


def test_bulk_mass_map_comparison_data_has_expected_shapes_and_window():
    """The bulk-mass comparison helper should expose a stable sample cloud."""
    data = bulk_mass_map_comparison_data(
        r_values=[0.0, 0.25, 1.0],
        overall_scale_values=[1.5, 3.0],
    )

    sample_count = 3 * 2 * 9
    assert data["eig_samples"].shape == (sample_count,)
    assert data["c_sigmoid_samples"].shape == (sample_count,)
    assert data["c_affine_samples"].shape == (sample_count,)
    assert data["sample_r_values"].shape == (sample_count,)
    assert data["sample_overall_scale_values"].shape == (sample_count,)
    assert data["sector_ids"].shape == (sample_count,)
    assert data["generation_ids"].shape == (sample_count,)
    assert data["eig_Q_samples"].shape == (3 * 2 * 3,)
    assert data["eig_u_samples"].shape == (3 * 2 * 3,)
    assert data["eig_d_samples"].shape == (3 * 2 * 3,)
    assert np.isclose(data["c_sigmoid_grid"][0], data["c_uv"][0])
    assert np.isclose(data["c_affine_grid"][0], data["c_uv"][0])
    assert np.all(data["c_sigmoid_samples"] <= data["c_uv"][0] + 1.0e-12)
    assert np.all(data["c_sigmoid_samples"] >= data["c_ir"][0] - 1.0e-12)


def test_bulk_mass_map_comparison_shows_tangent_match_and_large_lambda_saturation():
    """The sigmoid should agree with the affine tangent at small lambda and saturate later."""
    data = bulk_mass_map_comparison_data(
        r_values=[0.0, 0.25, 1.0],
        overall_scale_values=[0.5, 6.0],
    )

    delta = np.abs(data["c_sigmoid_samples"] - data["c_affine_samples"])
    small_lambda = data["eig_samples"] < 1.0e-2
    large_lambda = data["eig_samples"] > 1.0

    assert np.any(small_lambda)
    assert np.any(large_lambda)
    assert np.max(delta[small_lambda]) < 1.0e-3
    assert np.any(data["c_affine_samples"] < data["c_ir"][0])
    assert float(data["affine_below_c_ir_fraction"][0]) > 0.0
    assert np.median(delta[large_lambda]) > 5.0e-2
