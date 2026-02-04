"""Tests for scanParams module."""

import csv as csv_mod

import numpy as np

from scanParams import ScanConfig, run_scan


def test_benchmark_point_in_scan():
    """The Perez-Randall benchmark should appear in a 1-point scan."""
    config = ScanConfig(
        c_L_values=np.array([0.58]),
        c_N_values=np.array([0.27]),
        c_E_fixed=[0.75, 0.60, 0.50],
        Lambda_IR=3000.0,
        M_N=1.22e18,
        lightest_nu_mass=0.002,
        ordering='normal',
    )
    results = run_scan(config, progress_every=0)
    assert len(results) == 1
    row = results[0]
    # Perturbative (max |Y_bar| < 4 for default; Y_E_bar_3 ~ 5.4 exceeds)
    assert not row['perturbative']
    assert np.isclose(row['Y_E_bar_1'], 2.94, rtol=0.05)
    assert np.isclose(row['Y_N_bar_3'], 1.024, rtol=0.05)


def test_scan_csv_output(tmp_path):
    """CSV output should have correct headers and row count."""
    config = ScanConfig(
        c_L_values=np.array([0.55, 0.60]),
        c_N_values=np.array([0.25, 0.30]),
        c_E_fixed=[0.75, 0.60, 0.50],
    )
    csv_path = str(tmp_path / "test_scan.csv")
    results = run_scan(config, output_csv=csv_path, progress_every=0)
    assert len(results) == 4  # 2 x 2 grid

    with open(csv_path) as f:
        reader = csv_mod.DictReader(f)
        rows = list(reader)
    assert len(rows) == 4
    assert 'passes_all' in rows[0]
    assert 'reject_reason' in rows[0]


def test_extra_filter():
    """Custom filter should be applied and recorded."""
    config = ScanConfig(
        c_L_values=np.array([0.58]),
        c_N_values=np.array([0.27]),
        c_E_fixed=[0.75, 0.60, 0.50],
    )

    def always_reject(result):
        return (False, "custom_reject")

    results = run_scan(config, extra_filters=[always_reject], progress_every=0)
    assert not results[0]['passes_all']
    assert 'custom_reject' in results[0]['reject_reason']
