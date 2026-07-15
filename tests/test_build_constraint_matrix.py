from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest

from quarkConstraints.scales import KKGluonMassConventionError

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = REPO_ROOT / "scripts" / "build_constraint_matrix.py"


def _load_matrix_builder():
    spec = importlib.util.spec_from_file_location("build_constraint_matrix", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _write_row(tmp_path: Path, params: dict[str, object]) -> Path:
    path = tmp_path / "tile-00000.jsonl"
    row = {
        "params": params,
        "fit_diagnostics": {"success": True},
        "constraints": {
            "K001": {
                "passes": True,
                "ratio": 0.25,
            }
        },
    }
    path.write_text(json.dumps(row) + "\n", encoding="utf-8")
    return path


def test_build_constraint_matrix_writes_typed_mass_columns(tmp_path):
    builder = _load_matrix_builder()
    path = _write_row(
        tmp_path,
        {
            "M_KK": 7000.0,
            "m_kk_physical_gev": 7000.0,
            "Lambda_IR": 3500.0,
            "lambda_ir_gev": 3500.0,
            "xi_KK": 2.0,
            "mass_convention_id": "kk_gluon_m1_physical.xi_kk_times_lambda_ir.v1",
            "quark_fit_r": 0.05,
        },
    )

    records, seen_ids, n_total = builder.build_records([str(path)], collider_gev=5500.0)

    assert n_total == 1
    assert seen_ids == {"K001"}
    assert records[0]["M_KK_GeV"] == pytest.approx(7000.0)
    assert records[0]["m_kk_physical_gev"] == pytest.approx(7000.0)
    assert records[0]["lambda_ir_gev"] == pytest.approx(3500.0)
    assert records[0]["xi_KK"] == pytest.approx(2.0)
    assert records[0]["mass_convention_id"] == (
        "kk_gluon_m1_physical.xi_kk_times_lambda_ir.v1"
    )
    assert records[0]["pass_COLLIDER"] is True


def test_build_constraint_matrix_rejects_geometric_labeled_as_physical_tile(tmp_path):
    builder = _load_matrix_builder()
    path = _write_row(
        tmp_path,
        {
            "M_KK": 7000.0,
            "m_kk_physical_gev": 7000.0,
            "Lambda_IR": 7000.0,
            "lambda_ir_gev": 7000.0,
            "xi_KK": 2.0,
            "mass_convention_id": "kk_gluon_m1_physical.xi_kk_times_lambda_ir.v1",
        },
    )

    with pytest.raises(KKGluonMassConventionError, match="mass convention mismatch"):
        builder.build_records([str(path)], collider_gev=5500.0)


def test_build_constraint_matrix_rejects_xi_one_under_physical_label(tmp_path):
    builder = _load_matrix_builder()
    path = _write_row(
        tmp_path,
        {
            "M_KK": 7000.0,
            "m_kk_physical_gev": 7000.0,
            "Lambda_IR": 7000.0,
            "lambda_ir_gev": 7000.0,
            "xi_KK": 1.0,
            "mass_convention_id": "kk_gluon_m1_physical.xi_kk_times_lambda_ir.v1",
        },
    )

    with pytest.raises(KKGluonMassConventionError, match="physical KK-gluon"):
        builder.build_records([str(path)], collider_gev=5500.0)
