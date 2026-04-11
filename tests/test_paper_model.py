"""Tests for the paper-facing 0710.1869 benchmark/model helpers."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from quarkConstraints.paper_0710_1869.inputs import default_paper_0710_1869_table_i_inputs
from quarkConstraints.paper_0710_1869.model import (
    PAPER_0710_1869_EQ3_RELATION_ID,
    PAPER_0710_1869_MODEL_SCHEMA_ID,
    Paper07101869BenchmarkSector,
    Paper07101869DiagonalCMatrices,
    Paper07101869Eq3Example,
    Paper07101869Eq3ResidualSummary,
    Paper07101869TableIBenchmark,
    Paper07101869V5KMParameters,
    build_diagonal_c_matrix,
    build_table_i_benchmark_from_inputs,
    build_table_i_diagonal_c_matrices,
    default_paper_0710_1869_eq3_example,
    default_paper_0710_1869_table_i_benchmark,
    evaluate_default_paper_0710_1869_eq3_consistency,
    evaluate_eq3_diagonal_consistency,
    paper_v5km_matrix,
)
from quarkConstraints.paper_0710_1869.validation import module_has_forbidden_import

REPO_ROOT = Path(__file__).resolve().parents[1]
MODEL_MODULE_PATH = REPO_ROOT / "quarkConstraints" / "paper_0710_1869" / "model.py"


def test_default_table_i_benchmark_preserves_quoted_sector_values():
    benchmark = default_paper_0710_1869_table_i_benchmark()

    assert benchmark.schema_id == PAPER_0710_1869_MODEL_SCHEMA_ID
    assert benchmark.label == "table_i_benchmark"
    assert benchmark.q_sector.c_eigenvalues == pytest.approx((0.64, 0.59, 0.46))
    assert benchmark.u_sector.c_eigenvalues == pytest.approx((0.68, 0.53, -0.06))
    assert benchmark.d_sector.c_eigenvalues == pytest.approx((0.65, 0.60, 0.58))
    np.testing.assert_allclose(benchmark.q_sector.f_vector, [2.0e-3, 1.0e-2, 2.0e-1])
    np.testing.assert_allclose(benchmark.u_sector.f_vector, [7.0e-4, 6.0e-2, 8.0e-1])
    np.testing.assert_allclose(benchmark.d_sector.f_vector, [2.0e-3, 8.0e-3, 2.0e-2])


def test_diagonal_c_helpers_build_expected_shapes_and_entries():
    sector = Paper07101869BenchmarkSector(
        label="Q",
        c_eigenvalues=(0.64, 0.59, 0.46),
        f_eigenvalues=(2.0e-3, 1.0e-2, 2.0e-1),
    )
    direct = build_diagonal_c_matrix(sector.c_eigenvalues)
    bundled = build_table_i_diagonal_c_matrices()

    assert direct.shape == (3, 3)
    np.testing.assert_allclose(np.diag(direct), [0.64, 0.59, 0.46])
    np.testing.assert_allclose(direct - np.diag(np.diag(direct)), 0.0, atol=1.0e-15)
    assert bundled.C_Q.shape == (3, 3)
    assert bundled.C_u.shape == (3, 3)
    assert bundled.C_d.shape == (3, 3)
    np.testing.assert_allclose(np.diag(bundled.C_u), [0.68, 0.53, -0.06])
    np.testing.assert_allclose(np.diag(bundled.C_d), [0.65, 0.60, 0.58])


def test_model_benchmark_is_built_from_sourced_table_i_inputs() -> None:
    sourced = default_paper_0710_1869_table_i_inputs()
    benchmark = build_table_i_benchmark_from_inputs(sourced)

    assert isinstance(benchmark, Paper07101869TableIBenchmark)
    assert benchmark.q_sector.c_eigenvalues == pytest.approx((0.64, 0.59, 0.46))
    assert benchmark.u_sector.f_eigenvalues == pytest.approx((7.0e-4, 6.0e-2, 8.0e-1))


def test_paper_v5km_matrix_for_quoted_example_is_unitary():
    example = default_paper_0710_1869_eq3_example()
    v5km = paper_v5km_matrix(example.v5km_parameters)

    assert v5km.shape == (3, 3)
    assert example.relation_id == PAPER_0710_1869_EQ3_RELATION_ID
    assert example.v5km_parameters.theta12 == pytest.approx(np.deg2rad(115.0))
    assert example.v5km_parameters.theta23 == pytest.approx(np.deg2rad(65.0))
    assert example.v5km_parameters.theta13 == pytest.approx(np.deg2rad(70.0))
    np.testing.assert_allclose(v5km.conjugate().T @ v5km, np.eye(3), atol=1.0e-12)
    assert abs(np.linalg.det(v5km)) == pytest.approx(1.0, abs=1.0e-12)


def test_default_eq3_residual_summary_is_bounded_and_not_claimed_exact():
    summary = evaluate_default_paper_0710_1869_eq3_consistency()
    expected_residual = np.array([0.12379469, -0.01376433, -0.16003036])

    assert isinstance(summary, Paper07101869Eq3ResidualSummary)
    assert summary.exact_agreement is False
    np.testing.assert_allclose(summary.lhs_diag, [0.64, 0.59, 0.46], atol=1.0e-12)
    np.testing.assert_allclose(summary.residual, expected_residual, atol=1.0e-8)
    assert summary.max_rhs_imag_abs < 1.0e-12
    assert summary.max_abs_residual == pytest.approx(np.max(np.abs(expected_residual)))
    assert 1.0e-3 < summary.max_abs_residual < 0.25
    payload = summary.as_dict()
    assert payload["exact_agreement"] is False
    assert payload["relation_id"] == PAPER_0710_1869_EQ3_RELATION_ID


def test_eq3_residual_summary_can_recognize_an_exact_constructed_case():
    summary = evaluate_eq3_diagonal_consistency(
        c_q_eigenvalues=(1.0, 2.0, 3.0),
        c_u_eigenvalues=(9.0, 8.0, 7.0),
        c_d_eigenvalues=(1.0, 2.0, 3.0),
        a=1.0,
        r=0.0,
        v5km=np.eye(3, dtype=np.complex128),
        tolerance=1.0e-15,
    )

    assert summary.exact_agreement is True
    np.testing.assert_allclose(summary.rhs_diag, [1.0, 2.0, 3.0], atol=1.0e-15)
    np.testing.assert_allclose(summary.residual, [0.0, 0.0, 0.0], atol=1.0e-15)
    assert summary.max_abs_residual == pytest.approx(0.0, abs=1.0e-15)
    assert summary.l2_residual == pytest.approx(0.0, abs=1.0e-15)


@pytest.mark.parametrize(
    ("factory", "expected_message"),
    [
        (
            lambda: Paper07101869BenchmarkSector(
                label="",
                c_eigenvalues=(0.64, 0.59, 0.46),
                f_eigenvalues=(2.0e-3, 1.0e-2, 2.0e-1),
            ),
            "label",
        ),
        (
            lambda: Paper07101869BenchmarkSector(
                label="Q",
                c_eigenvalues=(0.64, 0.59),
                f_eigenvalues=(2.0e-3, 1.0e-2, 2.0e-1),
            ),
            "c_eigenvalues",
        ),
        (
            lambda: Paper07101869BenchmarkSector(
                label="Q",
                c_eigenvalues=(0.64, np.inf, 0.46),
                f_eigenvalues=(2.0e-3, 1.0e-2, 2.0e-1),
            ),
            "c_eigenvalues",
        ),
        (
            lambda: Paper07101869TableIBenchmark(schema_id="wrong.schema"),
            "schema_id",
        ),
        (
            lambda: Paper07101869TableIBenchmark(q_sector="bad"),
            "q_sector",
        ),
        (
            lambda: Paper07101869TableIBenchmark(u_sector="bad"),
            "u_sector",
        ),
        (
            lambda: Paper07101869TableIBenchmark(d_sector="bad"),
            "d_sector",
        ),
        (
            lambda: Paper07101869DiagonalCMatrices(
                C_Q=np.eye(3, dtype=np.complex128),
                C_u=np.array(
                    [[1.0, 1.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]],
                    dtype=np.complex128,
                ),
                C_d=np.eye(3, dtype=np.complex128),
            ),
            "C_u",
        ),
        (
            lambda: Paper07101869DiagonalCMatrices(
                C_Q=np.eye(3, dtype=np.complex128),
                C_u=np.eye(3, dtype=np.complex128),
                C_d=np.diag([1.0 + 1.0j, 2.0, 3.0]).astype(np.complex128),
            ),
            "C_d diagonal entries must be real",
        ),
        (
            lambda: Paper07101869V5KMParameters(theta12=float("nan")),
            "theta12",
        ),
        (
            lambda: Paper07101869Eq3Example(a=0.0),
            "a",
        ),
        (
            lambda: Paper07101869Eq3Example(v5km_parameters="bad"),
            "v5km_parameters",
        ),
        (
            lambda: paper_v5km_matrix(np.eye(3, dtype=np.complex128)),
            "params",
        ),
        (
            lambda: evaluate_eq3_diagonal_consistency(
                c_q_eigenvalues=(1.0, 2.0, 3.0),
                c_u_eigenvalues=(1.0, 2.0, 3.0),
                c_d_eigenvalues=(1.0, 2.0, 3.0),
                a=1.0,
                r=0.0,
                v5km=np.array(
                    [
                        [1.0, 0.0, 0.0],
                        [0.0, np.nan, 0.0],
                        [0.0, 0.0, 1.0],
                    ],
                    dtype=np.complex128,
                ),
            ),
            "finite",
        ),
        (
            lambda: evaluate_eq3_diagonal_consistency(
                c_q_eigenvalues=(1.0, 2.0, 3.0),
                c_u_eigenvalues=(1.0, 2.0, 3.0),
                c_d_eigenvalues=(1.0, 2.0, 3.0),
                a=1.0,
                r=0.0,
                v5km=np.eye(2, dtype=np.complex128),
            ),
            "v5km",
        ),
        (
            lambda: Paper07101869Eq3ResidualSummary(
                lhs_diag=np.array([1.0, 2.0, 3.0]),
                rhs_diag=np.array([1.0, 2.0, 3.0]),
                residual=np.array([0.0, 0.0, 0.0]),
                a=1.0,
                r=0.0,
                max_rhs_imag_abs=-1.0,
                max_abs_residual=0.0,
                l2_residual=0.0,
                exact_agreement=True,
                tolerance=1.0e-12,
            ),
            "max_rhs_imag_abs",
        ),
        (
            lambda: evaluate_eq3_diagonal_consistency(
                c_q_eigenvalues=(1.0, 2.0, 3.0),
                c_u_eigenvalues=(1.0, 2.0, 3.0),
                c_d_eigenvalues=(1.0, 2.0, 3.0),
                a=1.0,
                r=0.0,
                v5km=np.eye(3, dtype=np.complex128),
                tolerance=0.0,
            ),
            "tolerance",
        ),
    ],
)
def test_model_helpers_reject_invalid_contracts(factory, expected_message):
    with pytest.raises(ValueError, match=expected_message):
        factory()


def test_model_module_does_not_import_repo_v1_model_or_fit():
    assert MODEL_MODULE_PATH.exists()
    assert not module_has_forbidden_import(
        MODEL_MODULE_PATH,
        {
            "quarkConstraints.benchmarks",
            "quarkConstraints.couplings",
            "quarkConstraints.deltaf2",
            "quarkConstraints.fit",
            "quarkConstraints.model",
            "quarkConstraints.proxies",
            "quarkConstraints.scan",
        },
    )


def test_importing_model_module_does_not_load_repo_v1_runtime_modules() -> None:
    script = """
import importlib
import json
import sys

importlib.import_module("quarkConstraints.paper_0710_1869.model")
forbidden = sorted(
    name
    for name in sys.modules
    if name in {
        "quarkConstraints.benchmarks",
        "quarkConstraints.couplings",
        "quarkConstraints.deltaf2",
        "quarkConstraints.fit",
        "quarkConstraints.model",
        "quarkConstraints.proxies",
        "quarkConstraints.scan",
    }
)
print(json.dumps(forbidden))
"""
    completed = subprocess.run(
        [sys.executable, "-c", script],
        check=True,
        capture_output=True,
        text=True,
        cwd=str(REPO_ROOT),
    )
    assert completed.stdout.strip() == "[]"
