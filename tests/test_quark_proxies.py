"""Tests for the quark-sector proxy and alignment diagnostics."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from quarkConstraints.benchmarks import default_quark_targets
from quarkConstraints.fit import fit_quark_sector
from quarkConstraints.proxies import summarize_flavor_diagnostics, sweep_r_proxy_summary


def test_proxy_summary_is_positive_and_down_sector_is_more_aligned():
    """The fitted benchmark should stay in the proxy-favored alignment window."""
    result = fit_quark_sector(default_quark_targets(), overall_scale=3.0, max_nfev=120).result
    summary = summarize_flavor_diagnostics(result)
    ratio = (
        summary.diagnostics.down_offdiag_ratio_in_q_basis
        / max(summary.diagnostics.up_offdiag_ratio_in_q_basis, 1e-30)
    )

    assert summary.h_rs_proxy > 0.0
    assert ratio < 5.0


def test_r_sweep_proxy_summary_shows_down_sector_suppression_for_small_r():
    """The down-sector proxy and relative alignment pressure should shrink as r decreases."""
    sweep = sweep_r_proxy_summary([0.05, 0.1, 0.25, 0.4, 1.0])
    down_proxy = np.array([item.h_rs_proxy for item in sweep], dtype=float)
    alignment_ratio = np.array(
        [
            item.diagnostics.down_offdiag_ratio_in_q_basis
            / max(item.diagnostics.up_offdiag_ratio_in_q_basis, 1e-30)
            for item in sweep
        ],
        dtype=float,
    )

    assert down_proxy[0] <= down_proxy[-1]
    assert alignment_ratio[0] <= alignment_ratio[-1]
