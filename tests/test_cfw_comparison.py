from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts" / "rs_anarchy_cfw_comparison.py"
FOLLOWUP_SUMMARY = REPO / "scan_outputs" / "followup_crossings_summary.json"


def test_cfw_comparison_default_headline_matches_followup_summary(tmp_path):
    """The comparison driver default must preserve the signed-off headline.

    The full RUNA draws file is several GB, so this test builds a one-row
    fixture whose perturbative M_KK^min is calibrated from the current
    followup summary.  Running the driver with its defaults must map that row
    back to the same g_s^*=3 central value.
    """
    followup = json.loads(FOLLOWUP_SUMMARY.read_text())
    expected = followup["runs"]["RUNA"]["p50_gs_star_3_TeV"]
    gs_pert = followup["convention"]["g_s_pert"]
    gs_star = followup["convention"]["g_s_star"]
    p50_pert = expected * gs_pert / gs_star

    run_dir = tmp_path / "run"
    run_dir.mkdir()
    row = {
        "M_KK_GeV": 1000.0,
        "ok": True,
        "passes_pdg": True,
        "up_log_max": 0.0,
        "down_log_max": 0.0,
        "ckm_log_max": 0.0,
        "j_log": 0.0,
        "deltaf2_ratios": {
            "epsilon_K": p50_pert * p50_pert,
            "Delta_m_K": 0.0,
            "Delta_m_Bd": 0.0,
            "Delta_m_Bs": 0.0,
            "Delta_m_D0": 0.0,
        },
        "max_ratio": p50_pert * p50_pert,
    }
    (run_dir / "draws.jsonl").write_text(json.dumps(row) + "\n")

    summary_out = tmp_path / "summary.json"
    subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--run",
            str(run_dir),
            "--out-dir",
            str(tmp_path / "figures"),
            "--summary-out",
            str(summary_out),
            "--no-plot",
        ],
        cwd=REPO,
        check=True,
    )

    summary = json.loads(summary_out.read_text())
    got = summary["curves"]["post_audit_default"]["p50_gs_star_TeV"]
    assert summary["gs_star"] == pytest.approx(gs_star, rel=0, abs=0)
    assert got == pytest.approx(expected, rel=1e-12)
