"""Tests for scan post-processing/reclassification."""

import numpy as np

from scanParams import AnarchyConfig, ReclassifyConfig, ScanConfig, classify_row, run_scan


def _benchmark_row():
    config = ScanConfig(
        record_git_metadata=False,
        Lambda_IR_values=np.array([3000.0]),
        c_L_values=np.array([0.58]),
        c_N_values=np.array([0.27]),
        c_E_fixed=[0.75, 0.60, 0.50],
        MN_mode="fixed_ratio",
        MN_over_k=1.22e18 / 1.2209e19,
        lightest_nu_mass_values=np.array([0.002]),
    )
    return run_scan(config, progress_every=0)[0]


def test_classify_row_reproduces_default_core_flags():
    row = _benchmark_row()
    rec = classify_row(row, ReclassifyConfig())

    assert rec["reclass_perturbative"] == row["perturbative"]
    assert rec["reclass_natural"] == row["natural"]
    assert rec["reclass_lfv_passes"] == row["lfv_passes"]


def test_classify_row_accepts_csv_style_string_values():
    row = _benchmark_row()
    row_as_strings = {k: str(v) for k, v in row.items()}
    rec = classify_row(row_as_strings, ReclassifyConfig())

    assert isinstance(rec["reclass_perturbative"], bool)
    assert isinstance(rec["reclass_natural"], bool)
    assert isinstance(rec["reclass_lfv_passes"], bool)


def test_reclassify_allows_posthoc_anarchy_threshold_changes():
    row = _benchmark_row()
    base_kwargs = dict(
        max_Y_bar=10.0,
        naturalness_range=(1e-6, 10.0),
        require_lfv=False,
        anarchy=AnarchyConfig(),
    )

    rec_loose = classify_row(row, ReclassifyConfig(**base_kwargs, anarchy_min_score=-100.0))
    rec_tight = classify_row(row, ReclassifyConfig(**base_kwargs, anarchy_min_score=0.0))

    assert rec_loose["reclass_passes_all"]
    assert not rec_tight["reclass_passes_all"]
    assert "anarchy_score" in rec_tight["reclass_reject_reason"]
