"""Tests for the quark-sector proxy and alignment diagnostics."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from quarkConstraints.benchmarks import benchmark_spurion_input, default_quark_targets
from quarkConstraints.fit import fit_quark_sector
from quarkConstraints.model import QuarkSpurionPoint
from quarkConstraints.proxies import summarize_flavor_diagnostics, suppression_summary, sweep_r, sweep_r_proxy_summary


def test_proxy_summary_is_positive_and_down_sector_misalignment_stays_controlled():
    """The fitted benchmark should keep the down-sector misalignment stress O(1)."""
    result = fit_quark_sector(default_quark_targets(), overall_scale=3.0, max_nfev=120).result
    summary = summarize_flavor_diagnostics(result)
    ratio = summary.diagnostics.down_to_up_misalignment_ratio

    assert summary.h_rs_proxy > 0.0
    assert ratio < 5.0


def test_proxy_summary_respects_explicit_mkk_override():
    result = fit_quark_sector(default_quark_targets(), overall_scale=3.0, max_nfev=120).result
    default_summary = summarize_flavor_diagnostics(result)
    overridden = summarize_flavor_diagnostics(result, m_kk=6000.0)

    assert np.isclose(default_summary.m_kk, result.point.Lambda_IR)
    assert np.isclose(overridden.m_kk, 6000.0)
    assert np.isclose(overridden.h_rs_proxy / default_summary.h_rs_proxy, 0.25)


def test_r_sweep_proxy_summary_shows_down_sector_suppression_for_small_r():
    """The down-sector proxy and relative misalignment stress should shrink as r decreases."""
    sweep = sweep_r_proxy_summary([0.05, 0.1, 0.25, 0.4, 1.0])
    down_proxy = np.array([item.h_rs_proxy for item in sweep], dtype=float)
    alignment_ratio = np.array(
        [item.diagnostics.down_to_up_misalignment_ratio for item in sweep],
        dtype=float,
    )

    assert down_proxy[0] <= down_proxy[-1]
    assert alignment_ratio[0] <= alignment_ratio[-1]


def test_suppression_summary_uses_misalignment_keys_by_default():
    result = fit_quark_sector(default_quark_targets(), overall_scale=3.0, max_nfev=120).result
    summary = suppression_summary(result)

    assert "down_misalignment" in summary
    assert "up_misalignment" in summary
    assert "down_to_up_misalignment_ratio" in summary
    assert "down_alignment" not in summary


def test_sweep_r_tracks_the_supplied_point_rather_than_the_benchmark_seed():
    point_a = benchmark_spurion_input(r=0.25)
    point_b = QuarkSpurionPoint(
        Y_u=1.1 * point_a.Y_u,
        Y_d=0.9 * point_a.Y_d,
        r=point_a.r,
        Lambda_IR=point_a.Lambda_IR,
        k=point_a.k,
        v=point_a.v,
        label="modified-point",
    )
    rows_a = sweep_r(point_a, [0.1, 0.3])
    rows_b = sweep_r(point_b, [0.1, 0.3])

    assert rows_a != rows_b
    assert rows_a[0]["r"] == 0.1
    assert rows_b[1]["r"] == 0.3
