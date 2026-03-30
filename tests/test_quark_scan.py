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
    assert "proxy_h_rs" in rows[0]
    assert "passes_all" in rows[0]

    with open(csv_path, encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        file_rows = list(reader)

    assert len(file_rows) == 2
    assert "alignment_ratio" in file_rows[0]
