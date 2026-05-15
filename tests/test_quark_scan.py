"""Tests for the quark-sector scan wrapper."""

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from quarkConstraints.scan import QuarkScanConfig, run_quark_scan


def test_quark_scan_returns_rows_and_writes_csv(tmp_path):
    """A minimal quark scan should return rows and emit the documented schema."""
    csv_path = tmp_path / "quark_scan.csv"
    config = QuarkScanConfig(
        r_values=[0.1, 0.25],
        overall_scale_values=[3.0],
        Lambda_IR_values=[3000.0],
        record_git_metadata=False,
        max_nfev=80,
    )
    rows = run_quark_scan(config, output_csv=str(csv_path), progress_every=0)

    assert len(rows) == 2
    assert "fit_score" in rows[0]
    assert "M_KK" in rows[0]
    assert "xi_KK" in rows[0]
    assert "proxy_h_rs" in rows[0]
    assert "deltaf2_passes" in rows[0]
    assert "epsilon_k_ratio" in rows[0]
    assert "b_d_mix_ratio" in rows[0]
    assert "b_s_mix_ratio" in rows[0]
    assert "d_mix_ratio" in rows[0]
    assert "passes_all" in rows[0]

    with open(csv_path, encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        file_rows = list(reader)

    assert len(file_rows) == 2
    assert "alignment_ratio" in file_rows[0]
    assert "deltaf2_max_ratio" in file_rows[0]
    assert "fit_parameterization" in file_rows[0]


def test_quark_scan_threads_explicit_xi_kk_into_mkk():
    config = QuarkScanConfig(
        r_values=[0.25],
        overall_scale_values=[3.0],
        Lambda_IR_values=[3000.0],
        xi_KK=2.0,
        record_git_metadata=False,
        max_nfev=80,
    )
    row = run_quark_scan(config, progress_every=0)[0]

    assert row["xi_KK"] == 2.0
    assert row["M_KK"] == 6000.0


def test_quark_scan_rejects_points_that_fail_the_repo_proxy_gate():
    # Use a tightened proxy gate so the test does not rely on the exact
    # numerical proxy of r=0.4 (which depends on the target spectrum).
    config = QuarkScanConfig(
        r_values=[0.4],
        overall_scale_values=[3.0],
        Lambda_IR_values=[3000.0],
        record_git_metadata=False,
        max_nfev=100,
        max_proxy_h_rs=0.5,
    )
    row = run_quark_scan(config, progress_every=0)[0]

    assert row["proxy_h_rs"] > config.max_proxy_h_rs
    assert row["passes_all"] is False
    assert "proxy_h_rs" in row["reject_reason"]


def test_quark_scan_rejects_unimplemented_rng_seed_global():
    try:
        QuarkScanConfig(rng_seed_global=123)
    except ValueError as exc:
        assert "not yet supported" in str(exc)
    else:
        raise AssertionError("rng_seed_global should be rejected until stochastic seeding exists")
