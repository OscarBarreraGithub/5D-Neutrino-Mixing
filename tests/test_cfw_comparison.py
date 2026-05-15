from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts" / "rs_anarchy_cfw_comparison.py"
FOLLOWUP_SUMMARY = REPO / "scan_outputs" / "followup_crossings_summary.json"

sys.path.insert(0, str(REPO))
from scripts import rs_anarchy_cfw_comparison as cfw_script  # noqa: E402


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


def test_cfw_matched_projection_is_factor_2p2_at_matched_conventions(tmp_path):
    """The CFW comparison is not a percent-level match to the 21 TeV marker."""
    expected_gs3 = 23.37157
    matched_budget = abs(cfw_script.EPSILON_K_EXP - 1.81e-3)
    matched_scale = cfw_script.EPSILON_K_BUDGET_DEFAULT / matched_budget
    p50_pert = expected_gs3 * cfw_script.G_S_PERT / cfw_script.DEFAULT_GS_STAR

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
            "epsilon_K": (p50_pert * p50_pert) / matched_scale,
            "Delta_m_K": 0.0,
            "Delta_m_Bd": 0.0,
            "Delta_m_Bs": 0.0,
            "Delta_m_D0": 0.0,
        },
        "max_ratio": (p50_pert * p50_pert) / matched_scale,
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
            "--eps-k-sm",
            "1.81e-3",
            "--pdg-relative-tolerance",
            "0.30",
            "--summary-out",
            str(summary_out),
            "--no-plot",
        ],
        cwd=REPO,
        check=True,
    )

    summary = json.loads(summary_out.read_text())
    got_gs3 = summary["curves"]["cfw_matched"]["p50_gs_star_TeV"]
    got_gs6 = got_gs3 * (cfw_script.CFW_RS_GS6_TEV / cfw_script.CFW_RS_GS3_TEV)
    gs3_ratio = got_gs3 / cfw_script.CFW_RS_GS3_TEV
    gs6_ratio = got_gs6 / cfw_script.CFW_RS_GS6_TEV

    assert got_gs3 == pytest.approx(expected_gs3, rel=1e-6)
    assert gs3_ratio == pytest.approx(2.23, rel=0.01)
    assert gs6_ratio == pytest.approx(gs3_ratio, rel=1e-12)
    assert gs3_ratio > 2.0
