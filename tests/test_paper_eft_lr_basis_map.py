"""Tests for the frozen paper<->BMU LR basis map in the paper-owned RG slice."""

from __future__ import annotations

import importlib
from collections.abc import Sequence

import pytest

EXPECTED_CONTRACT_ID = "paper_q4q5_to_bmu_lr_basis.map_frozen.audit_ready.v3"
EXPECTED_STATUS_ID = (
    "paper_q4q5_to_bmu_lr_basis.map_frozen."
    "lr_rg_active.custom_lr_hadronic_active."
    "custom_lr_only_observable_active.custom_combined_observable_active."
    "default_export_lr_capable.v7"
)
EXPECTED_PAPER_OPERATOR_ORDER = ("Q4_LR", "Q5_LR")
EXPECTED_BMU_OPERATOR_ORDER = ("Q1_LR_BMU", "Q2_LR_BMU")
EXPECTED_PAPER_TO_BMU_OPERATOR_MAP = ((0.0, -2.0), (1.0, 0.0))
EXPECTED_BMU_TO_PAPER_OPERATOR_MAP = ((0.0, 1.0), (-0.5, 0.0))
EXPECTED_PAPER_TO_BMU_WILSON_MAP = ((0.0, -0.5), (1.0, 0.0))
EXPECTED_BMU_TO_PAPER_WILSON_MAP = ((0.0, 1.0), (-2.0, 0.0))


def _load_rg_module():
    return importlib.import_module("quarkConstraints.paper_0710_1869.eft_deltaf2.rg")


def _real_matrix_2x2(
    value: Sequence[Sequence[complex]],
) -> tuple[tuple[float, float], tuple[float, float]]:
    assert len(value) == 2
    rows: list[tuple[float, float]] = []
    for row in value:
        assert len(row) == 2
        rows.append(tuple(float(complex(entry).real) for entry in row))
    return (rows[0], rows[1])


def _matrix_multiply(
    left: tuple[tuple[float, float], tuple[float, float]],
    right: tuple[tuple[float, float], tuple[float, float]],
) -> tuple[tuple[float, float], tuple[float, float]]:
    return (
        (
            (left[0][0] * right[0][0]) + (left[0][1] * right[1][0]),
            (left[0][0] * right[0][1]) + (left[0][1] * right[1][1]),
        ),
        (
            (left[1][0] * right[0][0]) + (left[1][1] * right[1][0]),
            (left[1][0] * right[0][1]) + (left[1][1] * right[1][1]),
        ),
    )


def _matrix_vector_multiply(
    matrix: tuple[tuple[float, float], tuple[float, float]],
    vector: tuple[float, float],
) -> tuple[float, float]:
    return (
        (matrix[0][0] * vector[0]) + (matrix[0][1] * vector[1]),
        (matrix[1][0] * vector[0]) + (matrix[1][1] * vector[1]),
    )


def test_lr_basis_contract_exposes_exact_frozen_map_metadata() -> None:
    rg_module = _load_rg_module()
    contract = rg_module.default_paper_0710_1869_deltaf2_lr_basis_contract()

    assert contract.contract_id == EXPECTED_CONTRACT_ID
    assert contract.status_id == EXPECTED_STATUS_ID
    assert ".lr_rg_active." in contract.status_id
    assert ".custom_lr_hadronic_active." in contract.status_id
    assert ".custom_lr_only_observable_active." in contract.status_id
    assert ".custom_combined_observable_active." in contract.status_id
    assert contract.status_id.endswith(".default_export_lr_capable.v7")
    assert tuple(contract.paper_operator_order) == EXPECTED_PAPER_OPERATOR_ORDER
    assert tuple(contract.bmu_lr_operator_order) == EXPECTED_BMU_OPERATOR_ORDER
    assert contract.mapping_matrix_frozen is True
    assert contract.lr_running_activated is True
    assert tuple(
        tuple(float(value) for value in row) for row in contract.paper_to_bmu_operator_map_matrix
    ) == (
        EXPECTED_PAPER_TO_BMU_OPERATOR_MAP
    )
    assert tuple(
        tuple(float(value) for value in row) for row in contract.bmu_to_paper_operator_map_matrix
    ) == (
        EXPECTED_BMU_TO_PAPER_OPERATOR_MAP
    )
    assert tuple(
        tuple(float(value) for value in row) for row in contract.paper_to_bmu_wilson_map_matrix
    ) == (
        EXPECTED_PAPER_TO_BMU_WILSON_MAP
    )
    assert tuple(
        tuple(float(value) for value in row) for row in contract.bmu_to_paper_wilson_map_matrix
    ) == (
        EXPECTED_BMU_TO_PAPER_WILSON_MAP
    )


def test_public_lr_map_helper_matrices_match_the_frozen_contract() -> None:
    rg_module = _load_rg_module()

    assert _real_matrix_2x2(rg_module.paper_lr_to_bmu_lr_operator_map_matrix()) == (
        EXPECTED_PAPER_TO_BMU_OPERATOR_MAP
    )
    assert _real_matrix_2x2(rg_module.bmu_lr_to_paper_lr_operator_map_matrix()) == (
        EXPECTED_BMU_TO_PAPER_OPERATOR_MAP
    )
    assert _real_matrix_2x2(rg_module.paper_lr_to_bmu_lr_wilson_map_matrix()) == (
        EXPECTED_PAPER_TO_BMU_WILSON_MAP
    )
    assert _real_matrix_2x2(rg_module.bmu_lr_to_paper_lr_wilson_map_matrix()) == (
        EXPECTED_BMU_TO_PAPER_WILSON_MAP
    )


def test_lr_map_round_trip_identity_and_basis_vectors() -> None:
    rg_module = _load_rg_module()

    operator_forward = _real_matrix_2x2(rg_module.paper_lr_to_bmu_lr_operator_map_matrix())
    operator_inverse = _real_matrix_2x2(rg_module.bmu_lr_to_paper_lr_operator_map_matrix())
    wilson_forward = _real_matrix_2x2(rg_module.paper_lr_to_bmu_lr_wilson_map_matrix())
    wilson_inverse = _real_matrix_2x2(rg_module.bmu_lr_to_paper_lr_wilson_map_matrix())

    assert _matrix_multiply(operator_inverse, operator_forward) == (
        (1.0, 0.0),
        (0.0, 1.0),
    )
    assert _matrix_multiply(operator_forward, operator_inverse) == (
        (1.0, 0.0),
        (0.0, 1.0),
    )
    assert _matrix_multiply(wilson_inverse, wilson_forward) == (
        (1.0, 0.0),
        (0.0, 1.0),
    )
    assert _matrix_multiply(wilson_forward, wilson_inverse) == (
        (1.0, 0.0),
        (0.0, 1.0),
    )

    assert _matrix_vector_multiply(operator_forward, (1.0, 0.0)) == (0.0, 1.0)
    assert _matrix_vector_multiply(operator_forward, (0.0, 1.0)) == (-2.0, 0.0)
    assert _matrix_vector_multiply(wilson_forward, (1.0, 0.0)) == (0.0, 1.0)
    assert _matrix_vector_multiply(wilson_forward, (0.0, 1.0)) == (-0.5, 0.0)

    assert rg_module.map_paper_lr_wilsons_to_bmu_lr((1.0 + 0.0j, 0.0 + 0.0j)) == (
        0.0 + 0.0j,
        1.0 + 0.0j,
    )
    assert rg_module.map_paper_lr_wilsons_to_bmu_lr((0.0 + 0.0j, 1.0 + 0.0j)) == (
        -0.5 + 0.0j,
        0.0 + 0.0j,
    )
    assert rg_module.map_bmu_lr_wilsons_to_paper_lr((1.0 + 0.0j, 0.0 + 0.0j)) == (
        0.0 + 0.0j,
        -2.0 + 0.0j,
    )
    assert rg_module.map_bmu_lr_wilsons_to_paper_lr((0.0 + 0.0j, 1.0 + 0.0j)) == (
        1.0 + 0.0j,
        0.0 + 0.0j,
    )
    assert rg_module.map_bmu_lr_wilsons_to_paper_lr(
        rg_module.map_paper_lr_wilsons_to_bmu_lr((1.25 - 0.5j, -0.75 + 0.25j))
    ) == pytest.approx((1.25 - 0.5j, -0.75 + 0.25j))
