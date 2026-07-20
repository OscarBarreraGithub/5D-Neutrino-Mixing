"""Tests for the quark-sector scan wrapper."""

import csv
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from quarkConstraints.fit import QuarkFitSeed
from quarkConstraints.model import RotationParameters
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
    with pytest.warns(RuntimeWarning, match="xi_KK=1.0"):
        rows = run_quark_scan(config, output_csv=str(csv_path), progress_every=0)

    assert len(rows) == 2
    assert "fit_score" in rows[0]
    assert "M_KK" in rows[0]
    assert "m_kk_physical_gev" in rows[0]
    assert "lambda_ir_gev" in rows[0]
    assert "xi_KK" in rows[0]
    assert "mass_convention_id" in rows[0]
    assert "coupling_policy_id" in rows[0]
    assert "g_s_4d" in rows[0]
    assert "g_eff" in rows[0]
    assert "g_s_multiplier" in rows[0]
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
    assert "mass_convention_id" in file_rows[0]


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


def test_quark_scan_threads_epsilon_k_budget_override(monkeypatch):
    seen = {}

    def fake_evaluate_delta_f2_constraints(*args, **kwargs):
        seen["epsilon_k_np_budget_override"] = kwargs.get(
            "epsilon_k_np_budget_override"
        )
        by_system = {
            system: SimpleNamespace(ratio_to_bound=0.0, passes=True)
            for system in ("K", "B_d", "B_s", "D")
        }
        return SimpleNamespace(
            by_system=by_system,
            operator_convention="test_deltaf2",
            input_bundle="test_inputs",
            passes_all=True,
            max_ratio_to_bound=0.0,
        )

    monkeypatch.setattr(
        "quarkConstraints.scan.evaluate_delta_f2_constraints",
        fake_evaluate_delta_f2_constraints,
    )
    config = QuarkScanConfig(
        r_values=[0.25],
        overall_scale_values=[3.0],
        Lambda_IR_values=[3000.0],
        record_git_metadata=False,
        max_nfev=80,
        epsilon_k_np_budget_override=3.0e-4,
    )

    with pytest.warns(RuntimeWarning, match="xi_KK=1.0"):
        rows = run_quark_scan(config, progress_every=0)

    assert len(rows) == 1
    assert seen["epsilon_k_np_budget_override"] == 3.0e-4


def test_quark_scan_chained_seed_does_not_reapply_overall_scale(monkeypatch):
    seen_overall_scales = []

    seed = QuarkFitSeed(
        up_singular_values=np.ones(3),
        down_singular_values=np.ones(3),
        overall_scale=1.0,
        up_left=RotationParameters(),
        up_right=RotationParameters(),
        down_left=RotationParameters(),
        down_right=RotationParameters(),
    )
    state = SimpleNamespace(
        c_Q=np.ones(3),
        c_u=np.ones(3),
        c_d=np.ones(3),
        F_Q=np.ones(3),
        F_u=np.ones(3),
        F_d=np.ones(3),
    )
    result = SimpleNamespace(
        score=0.0,
        residual_norm=0.0,
        mass_residuals_up=np.zeros(3),
        mass_residuals_down=np.zeros(3),
        ckm_residuals=np.zeros(4),
        state=state,
        masses_up=np.ones(3),
        masses_down=np.ones(3),
        ckm_observables=np.zeros(4),
        point=SimpleNamespace(Y_u=np.eye(3), Y_d=np.eye(3)),
    )

    def fake_fit_quark_sector(*_args, **kwargs):
        seen_overall_scales.append(kwargs["overall_scale"])
        return SimpleNamespace(
            seed=seed,
            result=result,
            success=True,
            message="ok",
            nfev=1,
            initial_score=0.0,
        )

    def fake_evaluate_delta_f2_constraints(*_args, **_kwargs):
        by_system = {
            system: SimpleNamespace(ratio_to_bound=0.0, passes=True)
            for system in ("K", "B_d", "B_s", "D")
        }
        return SimpleNamespace(
            by_system=by_system,
            operator_convention="test_deltaf2",
            input_bundle="test_inputs",
            passes_all=True,
            max_ratio_to_bound=0.0,
        )

    monkeypatch.setattr("quarkConstraints.scan.fit_quark_sector", fake_fit_quark_sector)
    monkeypatch.setattr(
        "quarkConstraints.scan.summarize_flavor_diagnostics",
        lambda *_args, **_kwargs: SimpleNamespace(
            h_rs_proxy=0.0,
            diagnostics=SimpleNamespace(
                down_offdiag_ratio_in_q_basis=1.0,
                up_offdiag_ratio_in_q_basis=1.0,
            ),
        ),
    )
    monkeypatch.setattr(
        "quarkConstraints.scan.compute_quark_kk_gluon_couplings",
        lambda *_args, **_kwargs: SimpleNamespace(
            coupling_policy_id="test_policy",
            g_s_4d=1.0,
            g_eff=1.0,
            g_s_multiplier=1.0,
        ),
    )
    monkeypatch.setattr(
        "quarkConstraints.scan.evaluate_delta_f2_constraints",
        fake_evaluate_delta_f2_constraints,
    )
    config = QuarkScanConfig(
        r_values=[0.1, 0.25],
        overall_scale_values=[3.0],
        Lambda_IR_values=[3000.0],
        record_git_metadata=False,
        max_nfev=1,
    )

    with pytest.warns(RuntimeWarning, match="xi_KK=1.0"):
        rows = run_quark_scan(config, progress_every=0)

    assert len(rows) == 2
    assert seen_overall_scales == [3.0, None]


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
    with pytest.warns(RuntimeWarning, match="xi_KK=1.0"):
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
