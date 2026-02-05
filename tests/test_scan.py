"""Tests for scanParams module."""

import csv as csv_mod

import numpy as np
import pytest

from flavorConstraints import coefficient_from_br_limit
from scanParams import AnarchyConfig, ScanConfig, run_scan


def _benchmark_config(**overrides):
    """Build a 1-point benchmark config with deterministic metadata."""
    base = dict(
        record_git_metadata=False,
        rng_seed_global=12345,
        Lambda_IR_values=np.array([3000.0]),
        c_L_values=np.array([0.58]),
        c_N_values=np.array([0.27]),
        c_E_fixed=[0.75, 0.60, 0.50],
        MN_mode="fixed_ratio",
        MN_over_k=1.22e18 / 1.2209e19,
        lightest_nu_mass_values=np.array([0.002]),
    )
    base.update(overrides)
    return ScanConfig(**base)


def test_benchmark_point_in_scan():
    """Perez-Randall-like benchmark point should evaluate consistently."""
    config = _benchmark_config()
    results = run_scan(config, progress_every=0)
    assert len(results) == 1

    row = results[0]
    # Perturbative (max |Y_bar| < 4 for default; Y_E_bar_3 ~ 5.4 exceeds)
    assert not row["perturbative"]
    assert np.isclose(row["Y_E_bar_1"], 2.94, rtol=0.05)
    assert np.isclose(row["Y_N_bar_3"], 1.024, rtol=0.05)
    assert np.isclose(row["M_N"], 1.22e18, rtol=1e-12)
    assert np.isclose(row["M_KK"], 3000.0)
    assert np.isclose(
        row["lfv_C"],
        coefficient_from_br_limit(config.br_limit, prefactor=config.prefac_br),
        rtol=1e-12,
    )


def test_scan_csv_output(tmp_path):
    """CSV output should have correct headers and row count."""
    config = ScanConfig(
        record_git_metadata=False,
        Lambda_IR_values=np.array([3000.0]),
        c_L_values=np.array([0.55, 0.60]),
        c_N_values=np.array([0.25, 0.30]),
        c_E_fixed=[0.75, 0.60, 0.50],
        MN_mode="fixed_ratio",
        MN_over_k=0.1,
        lightest_nu_mass_values=np.array([0.002]),
    )
    csv_path = str(tmp_path / "test_scan.csv")
    results = run_scan(config, output_csv=csv_path, progress_every=0)
    assert len(results) == 4  # 2 x 2 grid

    with open(csv_path, encoding="utf-8") as handle:
        reader = csv_mod.DictReader(handle)
        rows = list(reader)
    assert len(rows) == 4
    assert "passes_all" in rows[0]
    assert "reject_reason" in rows[0]
    assert "lfv_model" in rows[0]
    assert "anarchy_score" in rows[0]


def test_extra_filter():
    """Custom filter should be applied and recorded."""
    config = _benchmark_config()

    def always_reject(result):
        return (False, "custom_reject")

    results = run_scan(config, extra_filters=[always_reject], progress_every=0)
    assert not results[0]["passes_all"]
    assert "custom_reject" in results[0]["reject_reason"]


@pytest.mark.parametrize(
    ("kwargs", "expected_msg"),
    [
        ({"c_E_fixed": None, "c_E_grid": None}, "Either c_E_fixed or c_E_grid"),
        ({"c_E_fixed": [0.7, 0.6]}, "c_E_fixed must contain exactly 3 values"),
        (
            {"c_E_grid": [np.array([0.7]), np.array([0.6])]},
            "c_E_grid must contain exactly 3 arrays",
        ),
        ({"MN_mode": "scan_ratio", "MN_over_k_values": None}, "MN_over_k_values must be provided"),
        ({"ordering": "inverted"}, "supports only ordering='normal'"),
    ],
)
def test_scan_config_validates_inputs(kwargs, expected_msg):
    """ScanConfig should reject malformed inputs."""
    base = dict(
        Lambda_IR_values=np.array([3000.0]),
        c_L_values=np.array([0.58]),
        c_N_values=np.array([0.27]),
        lightest_nu_mass_values=np.array([0.002]),
        MN_mode="fixed_ratio",
        MN_over_k=0.1,
    )
    base.update(kwargs)
    with pytest.raises(ValueError, match=expected_msg):
        ScanConfig(**base)


def test_scan_config_accepts_c_e_grid():
    """ScanConfig should accept a 3-array c_E_grid."""
    config = ScanConfig(
        Lambda_IR_values=np.array([3000.0]),
        c_L_values=np.array([0.58]),
        c_N_values=np.array([0.27]),
        lightest_nu_mass_values=np.array([0.002]),
        MN_mode="fixed_ratio",
        MN_over_k=0.1,
        c_E_fixed=None,
        c_E_grid=[
            np.array([0.70, 0.72]),
            np.array([0.55]),
            np.array([0.45]),
        ],
    )
    assert config.total_points == 2


def test_sorts_c_e_descending_for_fixed_point():
    """Fixed c_E values should be sorted descending by default."""
    config = _benchmark_config(c_E_fixed=[0.50, 0.80, 0.60])
    row = run_scan(config, progress_every=0)[0]
    assert (row["c_E1"], row["c_E2"], row["c_E3"]) == (0.8, 0.6, 0.5)


def test_lfv_uses_xi_kk_mapping():
    """LFV RHS should use M_KK = xi_KK * Lambda_IR."""
    config = _benchmark_config(xi_KK=2.0)
    row = run_scan(config, progress_every=0)[0]

    c_val = coefficient_from_br_limit(config.br_limit, prefactor=config.prefac_br)
    assert np.isclose(row["M_KK"], 6000.0)
    assert np.isclose(row["lfv_rhs"], c_val * (6000.0 / config.lfv_reference_scale) ** 2)


def test_anarchy_scoring_is_deterministic_for_fixed_seed():
    """Anarchy score should be reproducible with a fixed global seed."""
    config = _benchmark_config(anarchy=AnarchyConfig(), rng_seed_global=777)

    row1 = run_scan(config, progress_every=0)[0]
    row2 = run_scan(config, progress_every=0)[0]

    assert np.isclose(row1["anarchy_score"], row2["anarchy_score"])
    assert np.isclose(row1["anarchy_yN_overall"], row2["anarchy_yN_overall"])


def test_anarchy_min_score_filter_rejects_low_score_points():
    """Points below anarchy_min_score should fail with explicit reason."""
    config = _benchmark_config(
        anarchy=AnarchyConfig(),
        anarchy_min_score=0.0,
        rng_seed_global=42,
    )
    row = run_scan(config, progress_every=0)[0]
    assert not row["passes_all"]
    assert "anarchy_score" in row["reject_reason"]
