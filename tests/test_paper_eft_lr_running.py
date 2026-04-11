"""Acceptance tests for the LR-RG-1 paper-owned LR evolution slice."""

from __future__ import annotations

import importlib

import pytest

PROBE_INPUT_PAPER_LR = (1.25 - 0.5j, -0.75 + 0.25j)


def _load_rg_module():
    return importlib.import_module("quarkConstraints.paper_0710_1869.eft_deltaf2.rg")


def _load_operators_module():
    return importlib.import_module("quarkConstraints.paper_0710_1869.eft_deltaf2.operators")


def _complex_matrix_multiply(
    left: tuple[tuple[complex, complex], tuple[complex, complex]],
    right: tuple[tuple[complex, complex], tuple[complex, complex]],
) -> tuple[tuple[complex, complex], tuple[complex, complex]]:
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


def _complex_matrix_vector_multiply(
    matrix: tuple[tuple[complex, complex], tuple[complex, complex]],
    vector: tuple[complex, complex],
) -> tuple[complex, complex]:
    return (
        (matrix[0][0] * vector[0]) + (matrix[0][1] * vector[1]),
        (matrix[1][0] * vector[0]) + (matrix[1][1] * vector[1]),
    )


def _assert_complex_pair_close(
    observed: tuple[complex, complex],
    expected: tuple[complex, complex],
    *,
    abs_tol: float = 1.0e-12,
) -> None:
    for lhs, rhs in zip(observed, expected, strict=True):
        assert lhs.real == pytest.approx(rhs.real, rel=0.0, abs=abs_tol)
        assert lhs.imag == pytest.approx(rhs.imag, rel=0.0, abs=abs_tol)


def _assert_complex_matrix_close(
    observed: tuple[tuple[complex, complex], tuple[complex, complex]],
    expected: tuple[tuple[complex, complex], tuple[complex, complex]],
    *,
    abs_tol: float = 1.0e-12,
) -> None:
    for observed_row, expected_row in zip(observed, expected, strict=True):
        _assert_complex_pair_close(observed_row, expected_row, abs_tol=abs_tol)


def _lr_block(
    paper_basis_evolution_matrix: tuple[tuple[complex, ...], ...],
) -> tuple[tuple[complex, complex], tuple[complex, complex]]:
    return (
        (
            complex(paper_basis_evolution_matrix[2][2]),
            complex(paper_basis_evolution_matrix[2][3]),
        ),
        (
            complex(paper_basis_evolution_matrix[3][2]),
            complex(paper_basis_evolution_matrix[3][3]),
        ),
    )


def _toy_wilsons(*, matching_scale_GeV: float, q4_lr: complex, q5_lr: complex):
    operators_module = _load_operators_module()
    return operators_module.Paper07101869DeltaF2WilsonCoefficients(
        contract=operators_module.default_paper_0710_1869_deltaf2_wilson_contract(),
        benchmark_id="toy_lr_running_probe",
        scale_label=f"toy_mu_{matching_scale_GeV:g}",
        system_id="kaon",
        sector_id="down",
        generations=(0, 1),
        matching_scale_GeV=float(matching_scale_GeV),
        propagator_mass_GeV=3000.0,
        left_coupling=0.0,
        right_coupling=0.0,
        q1_vll=0.0,
        q1_vrr=0.0,
        q4_lr=q4_lr,
        q5_lr=q5_lr,
    )


def test_paper_lr_evolution_block_equals_winv_u_bmu_w() -> None:
    rg_module = _load_rg_module()
    evolution = rg_module.compute_deltaf2_lo_evolution_matrix(3000.0, 2.0)
    expected = _complex_matrix_multiply(
        rg_module.bmu_lr_to_paper_lr_wilson_map_matrix(),
        _complex_matrix_multiply(
            evolution.bmu_lr_evolution_matrix,
            rg_module.paper_lr_to_bmu_lr_wilson_map_matrix(),
        ),
    )

    _assert_complex_matrix_close(_lr_block(evolution.paper_basis_evolution_matrix), expected)


def test_thresholded_lr_segment_product_matches_bmu_total_and_maps_back() -> None:
    rg_module = _load_rg_module()
    evolution = rg_module.compute_deltaf2_lo_evolution_matrix(3000.0, 2.0)

    segment_product = (
        (1.0 + 0.0j, 0.0 + 0.0j),
        (0.0 + 0.0j, 1.0 + 0.0j),
    )
    for segment in evolution.segments:
        segment_product = _complex_matrix_multiply(segment.bmu_lr_matrix, segment_product)

    _assert_complex_matrix_close(segment_product, evolution.bmu_lr_evolution_matrix)

    mapped_input = rg_module.map_paper_lr_wilsons_to_bmu_lr(PROBE_INPUT_PAPER_LR)
    mapped_back = rg_module.map_bmu_lr_wilsons_to_paper_lr(
        _complex_matrix_vector_multiply(segment_product, mapped_input)
    )
    expected = _complex_matrix_vector_multiply(
        _lr_block(evolution.paper_basis_evolution_matrix),
        PROBE_INPUT_PAPER_LR,
    )
    _assert_complex_pair_close(expected, mapped_back)


def test_public_lr_running_matches_conjugated_bmu_result() -> None:
    rg_module = _load_rg_module()
    evolution = rg_module.compute_deltaf2_lo_evolution_matrix(3000.0, 2.0)
    expected = rg_module.map_bmu_lr_wilsons_to_paper_lr(
        _complex_matrix_vector_multiply(
            evolution.bmu_lr_evolution_matrix,
            rg_module.map_paper_lr_wilsons_to_bmu_lr(PROBE_INPUT_PAPER_LR),
        )
    )

    public_result = rg_module.evolve_deltaf2_wilsons_lo(
        _toy_wilsons(
            matching_scale_GeV=3000.0,
            q4_lr=PROBE_INPUT_PAPER_LR[0],
            q5_lr=PROBE_INPUT_PAPER_LR[1],
        ),
        mu_low_GeV=2.0,
    )
    observed = (
        complex(public_result.coefficients["Q4_LR"]),
        complex(public_result.coefficients["Q5_LR"]),
    )

    _assert_complex_pair_close(observed, expected)
