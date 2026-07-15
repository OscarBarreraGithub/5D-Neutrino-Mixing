"""Tests for the paper-facing 0710.1869 benchmark/model helpers."""

from __future__ import annotations

from dataclasses import replace
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from quarkConstraints.paper_0710_1869.benchmarks import (
    Paper07101869BenchmarkSpurionSeed,
    build_paper_0710_1869_seeded_structural_reference_point,
    build_paper_0710_1869_seeded_physical_point,
    default_paper_0710_1869_pr1_benchmark,
)
from quarkConstraints.paper_0710_1869.conventions import (
    PAPER_0710_1869_PROFILE_DERIVATION_POLICY_ID,
    PAPER_0710_1869_SEED_TO_PROFILE_MAPPING_POLICY_ID,
    PAPER_0710_1869_UNIVERSAL_TERM_COEFFICIENT_POLICY_ID,
)
from quarkConstraints.paper_0710_1869.inputs import default_paper_0710_1869_table_i_inputs
from quarkConstraints.paper_0710_1869.inputs import (
    PAPER_0710_1869_AFFINE_BULK_MASS_LEADING_TERM_COEFFICIENT,
    PAPER_0710_1869_AFFINE_BULK_MASS_UNIVERSAL_OFFSET,
    PAPER_0710_1869_PHYSICAL_SEED_TO_PROFILE_CONTRACT_SCHEMA_ID,
    default_paper_0710_1869_physical_seed_to_profile_contract,
)
from quarkConstraints.paper_0710_1869.model import (
    PAPER_0710_1869_EQ3_RELATION_ID,
    PAPER_0710_1869_MODEL_SCHEMA_ID,
    PAPER_0710_1869_NONCANONICAL_CLAIM_LEVEL_ID,
    PAPER_0710_1869_DERIVED_PROFILE_INPUT_POLICY_ID,
    PAPER_0710_1869_PHYSICAL_BENCHMARK_SEED_CONSTRUCTION_ID,
    PAPER_0710_1869_PHYSICAL_BULK_STATE_SCHEMA_ID,
    PAPER_0710_1869_PHYSICAL_CLAIM_LEVEL_ID,
    PAPER_0710_1869_PHYSICAL_CORE_CONTRACT_ID,
    PAPER_0710_1869_PHYSICAL_POINT_KIND_ID,
    PAPER_0710_1869_PHYSICAL_POINT_SCHEMA_ID,
    PAPER_0710_1869_POINT_DERIVED_PHYSICAL_STATUS_ID,
    PAPER_0710_1869_NO_BULK_MASS_MAP_POLICY_ID,
    PAPER_0710_1869_REFERENCE_ONLY_PROFILE_INPUT_POLICY_ID,
    PAPER_0710_1869_STRUCTURAL_REFERENCE_STATUS_ID,
    Paper07101869PhysicalBulkState,
    Paper07101869PhysicalPoint,
    Paper07101869BenchmarkSector,
    Paper07101869DiagonalCMatrices,
    Paper07101869Eq3Example,
    Paper07101869Eq3ResidualSummary,
    Paper07101869TableIBenchmark,
    Paper07101869RotationParameters,
    Paper07101869V5KMParameters,
    build_diagonal_c_matrix,
    build_table_i_benchmark_from_inputs,
    build_table_i_diagonal_c_matrices,
    default_paper_0710_1869_eq3_example,
    default_paper_0710_1869_table_i_benchmark,
    evaluate_default_paper_0710_1869_eq3_consistency,
    derive_paper_0710_1869_bulk_state,
    derive_paper_0710_1869_physical_bulk_state,
    evaluate_eq3_diagonal_consistency,
    paper_v5km_matrix,
    _ordered_hermitian_spectrum,
)
from quarkConstraints.paper_0710_1869.validation import module_has_forbidden_import
from warpConfig.wavefuncs import f_IR

REPO_ROOT = Path(__file__).resolve().parents[1]
MODEL_MODULE_PATH = REPO_ROOT / "quarkConstraints" / "paper_0710_1869" / "model.py"


def _physical_sector_policy(contract, sector_id: str):
    for policy in contract.universal_term_policy.sector_policies:
        if policy.sector_id == sector_id:
            return policy
    raise AssertionError(f"missing physical sector policy {sector_id!r}")


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


def _seed() -> Paper07101869BenchmarkSpurionSeed:
    return Paper07101869BenchmarkSpurionSeed(
        up_singular_values=(0.45, 1.2, 3.4),
        down_singular_values=(0.15, 0.55, 1.7),
        overall_scale=0.8,
        up_left=Paper07101869RotationParameters.from_degrees(
            theta12_deg=19.0,
            theta13_deg=4.0,
            theta23_deg=11.0,
            delta=0.3,
        ),
        up_right=Paper07101869RotationParameters.from_degrees(
            theta12_deg=13.0,
            theta13_deg=7.0,
            theta23_deg=5.0,
            delta=0.1,
        ),
        down_left=Paper07101869RotationParameters.from_degrees(
            theta12_deg=7.0,
            theta13_deg=2.0,
            theta23_deg=9.0,
            delta=0.5,
        ),
        down_right=Paper07101869RotationParameters.from_degrees(
            theta12_deg=17.0,
            theta13_deg=3.0,
            theta23_deg=6.0,
            delta=0.2,
        ),
        notes="test-seed",
    )


def _seeded_point():
    benchmark = default_paper_0710_1869_pr1_benchmark()
    return build_paper_0710_1869_seeded_structural_reference_point(
        benchmark,
        _seed(),
        metadata={"test_case": "model_acceptance"},
    )


def _physical_point():
    benchmark = default_paper_0710_1869_pr1_benchmark()
    return build_paper_0710_1869_seeded_physical_point(
        benchmark,
        _seed(),
        metadata={"test_case": "model_acceptance"},
    )


def _degenerate_physical_point():
    benchmark = default_paper_0710_1869_pr1_benchmark()
    degenerate_seed = Paper07101869BenchmarkSpurionSeed(
        up_singular_values=(0.75, 0.75, 1.5),
        down_singular_values=(0.25, 0.25, 0.9),
        overall_scale=0.8,
        notes="degenerate-seed",
    )
    return build_paper_0710_1869_seeded_physical_point(
        benchmark,
        degenerate_seed,
        metadata={"test_case": "degenerate_physical_acceptance"},
    )


def test_seeded_structural_reference_state_separates_reference_aliases_from_derived_profiles():
    point = _seeded_point()
    bulk_state = derive_paper_0710_1869_bulk_state(point)

    assert point.profile_input_policy_id == PAPER_0710_1869_REFERENCE_ONLY_PROFILE_INPUT_POLICY_ID
    assert bulk_state.claim_level_id == PAPER_0710_1869_NONCANONICAL_CLAIM_LEVEL_ID
    assert bulk_state.bulk_mass_map_policy_id == PAPER_0710_1869_NO_BULK_MASS_MAP_POLICY_ID
    assert bulk_state.structural_reference_status_id == PAPER_0710_1869_STRUCTURAL_REFERENCE_STATUS_ID
    assert bulk_state.profile_consistency_Q.policy_id == PAPER_0710_1869_REFERENCE_ONLY_PROFILE_INPUT_POLICY_ID
    assert bulk_state.profile_consistency_u.policy_id == PAPER_0710_1869_REFERENCE_ONLY_PROFILE_INPUT_POLICY_ID
    assert bulk_state.profile_consistency_d.policy_id == PAPER_0710_1869_REFERENCE_ONLY_PROFILE_INPUT_POLICY_ID

    np.testing.assert_allclose(bulk_state.reference_diagnostics_Q.quoted_c, bulk_state.c_Q)
    np.testing.assert_allclose(bulk_state.reference_diagnostics_u.quoted_c, bulk_state.c_u)
    np.testing.assert_allclose(bulk_state.reference_diagnostics_d.quoted_c, bulk_state.c_d)
    np.testing.assert_allclose(bulk_state.reference_diagnostics_Q.quoted_f, bulk_state.F_Q)
    np.testing.assert_allclose(bulk_state.reference_diagnostics_u.quoted_f, bulk_state.F_u)
    np.testing.assert_allclose(bulk_state.reference_diagnostics_d.quoted_f, bulk_state.F_d)
    np.testing.assert_allclose(
        bulk_state.reference_diagnostics_Q.geometry_derived_f_from_quoted_c,
        bulk_state.derived_F_Q,
    )
    np.testing.assert_allclose(
        bulk_state.reference_diagnostics_u.geometry_derived_f_from_quoted_c,
        bulk_state.derived_F_u,
    )
    np.testing.assert_allclose(
        bulk_state.reference_diagnostics_d.geometry_derived_f_from_quoted_c,
        bulk_state.derived_F_d,
    )
    assert not np.allclose(bulk_state.F_Q, bulk_state.derived_F_Q)
    assert not np.allclose(bulk_state.F_u, bulk_state.derived_F_u)
    assert not np.allclose(bulk_state.F_d, bulk_state.derived_F_d)


def test_geometry_derived_profiles_are_deterministic_under_frozen_reference_contract():
    point = _seeded_point()
    relabeled_point = replace(
        point,
        q_sector=replace(point.q_sector, f_eigenvalues=(9.0e-4, 1.7e-2, 3.5e-1)),
        u_sector=replace(point.u_sector, f_eigenvalues=(1.3e-3, 7.1e-2, 9.2e-1)),
        d_sector=replace(point.d_sector, f_eigenvalues=(3.1e-3, 9.4e-3, 3.4e-2)),
    )

    base_state = derive_paper_0710_1869_bulk_state(point)
    relabeled_state = derive_paper_0710_1869_bulk_state(relabeled_point)

    np.testing.assert_allclose(base_state.derived_F_Q, relabeled_state.derived_F_Q)
    np.testing.assert_allclose(base_state.derived_F_u, relabeled_state.derived_F_u)
    np.testing.assert_allclose(base_state.derived_F_d, relabeled_state.derived_F_d)
    assert not np.allclose(base_state.F_Q, relabeled_state.F_Q)
    assert not np.allclose(base_state.F_u, relabeled_state.F_u)
    assert not np.allclose(base_state.F_d, relabeled_state.F_d)


def test_physical_contract_is_exactly_frozen():
    contract = default_paper_0710_1869_physical_seed_to_profile_contract()

    assert contract.schema_id == PAPER_0710_1869_PHYSICAL_SEED_TO_PROFILE_CONTRACT_SCHEMA_ID
    assert contract.mapping_policy.policy_id == PAPER_0710_1869_SEED_TO_PROFILE_MAPPING_POLICY_ID
    assert (
        contract.mapping_policy.profile_derivation_policy_id
        == PAPER_0710_1869_PROFILE_DERIVATION_POLICY_ID
    )
    assert (
        contract.universal_term_policy.policy_id
        == PAPER_0710_1869_UNIVERSAL_TERM_COEFFICIENT_POLICY_ID
    )
    assert not contract.mapping_policy.uses_hidden_bulk_mass_map_surrogate
    for sector_id in ("Q", "u", "d"):
        policy = _physical_sector_policy(contract, sector_id)
        assert policy.leading_term_coefficient == pytest.approx(
            PAPER_0710_1869_AFFINE_BULK_MASS_LEADING_TERM_COEFFICIENT
        )
        assert policy.universal_offset == pytest.approx(
            PAPER_0710_1869_AFFINE_BULK_MASS_UNIVERSAL_OFFSET
        )


def test_physical_point_and_bulk_state_use_point_derived_frozen_contract_and_stay_distinct_from_reference_only_path():
    physical_point = _physical_point()
    physical_bulk_state = derive_paper_0710_1869_physical_bulk_state(physical_point)
    physical_bulk_state_again = derive_paper_0710_1869_physical_bulk_state(_physical_point())
    reference_point = _seeded_point()
    reference_bulk_state = derive_paper_0710_1869_bulk_state(reference_point)

    assert isinstance(physical_point, Paper07101869PhysicalPoint)
    assert isinstance(physical_bulk_state, Paper07101869PhysicalBulkState)
    assert physical_point.schema_id == PAPER_0710_1869_PHYSICAL_POINT_SCHEMA_ID
    assert physical_bulk_state.schema_id == PAPER_0710_1869_PHYSICAL_BULK_STATE_SCHEMA_ID
    assert physical_point.point_kind_id == PAPER_0710_1869_PHYSICAL_POINT_KIND_ID
    assert physical_point.claim_level_id == PAPER_0710_1869_PHYSICAL_CLAIM_LEVEL_ID
    assert physical_point.profile_input_policy_id == PAPER_0710_1869_DERIVED_PROFILE_INPUT_POLICY_ID
    assert physical_point.physical_profile_status_id == PAPER_0710_1869_POINT_DERIVED_PHYSICAL_STATUS_ID
    assert physical_point.construction_id == PAPER_0710_1869_PHYSICAL_BENCHMARK_SEED_CONSTRUCTION_ID
    assert physical_point.physical_contract == default_paper_0710_1869_physical_seed_to_profile_contract()
    assert physical_bulk_state.core_contract_id == PAPER_0710_1869_PHYSICAL_CORE_CONTRACT_ID
    assert physical_bulk_state.point_kind_id == PAPER_0710_1869_PHYSICAL_POINT_KIND_ID
    assert physical_bulk_state.claim_level_id == PAPER_0710_1869_PHYSICAL_CLAIM_LEVEL_ID
    assert physical_bulk_state.physical_profile_status_id == PAPER_0710_1869_POINT_DERIVED_PHYSICAL_STATUS_ID
    assert physical_bulk_state.bulk_mass_map_policy_id == PAPER_0710_1869_SEED_TO_PROFILE_MAPPING_POLICY_ID
    assert physical_bulk_state.profile_input_policy_id == PAPER_0710_1869_DERIVED_PROFILE_INPUT_POLICY_ID
    assert physical_bulk_state.physical_contract == default_paper_0710_1869_physical_seed_to_profile_contract()

    q_policy = _physical_sector_policy(physical_bulk_state.physical_contract, "Q")
    u_policy = _physical_sector_policy(physical_bulk_state.physical_contract, "u")
    d_policy = _physical_sector_policy(physical_bulk_state.physical_contract, "d")
    np.testing.assert_allclose(
        physical_bulk_state.c_Q,
        (q_policy.leading_term_coefficient * physical_bulk_state.eig_Q)
        + q_policy.universal_offset,
        atol=1.0e-12,
    )
    np.testing.assert_allclose(
        physical_bulk_state.c_u,
        (u_policy.leading_term_coefficient * physical_bulk_state.eig_u)
        + u_policy.universal_offset,
        atol=1.0e-12,
    )
    np.testing.assert_allclose(
        physical_bulk_state.c_d,
        (d_policy.leading_term_coefficient * physical_bulk_state.eig_d)
        + d_policy.universal_offset,
        atol=1.0e-12,
    )
    np.testing.assert_allclose(
        physical_bulk_state.F_Q,
        np.asarray(f_IR(physical_bulk_state.c_Q, physical_bulk_state.epsilon), dtype=float),
        atol=1.0e-12,
    )
    np.testing.assert_allclose(
        physical_bulk_state.F_u,
        np.asarray(f_IR(physical_bulk_state.c_u, physical_bulk_state.epsilon), dtype=float),
        atol=1.0e-12,
    )
    np.testing.assert_allclose(
        physical_bulk_state.F_d,
        np.asarray(f_IR(physical_bulk_state.c_d, physical_bulk_state.epsilon), dtype=float),
        atol=1.0e-12,
    )
    np.testing.assert_allclose(physical_bulk_state.c_Q, physical_bulk_state_again.c_Q, atol=1.0e-12)
    np.testing.assert_allclose(physical_bulk_state.F_Q, physical_bulk_state_again.F_Q, atol=1.0e-12)
    np.testing.assert_allclose(physical_bulk_state.c_u, physical_bulk_state_again.c_u, atol=1.0e-12)
    np.testing.assert_allclose(physical_bulk_state.F_u, physical_bulk_state_again.F_u, atol=1.0e-12)
    np.testing.assert_allclose(physical_bulk_state.c_d, physical_bulk_state_again.c_d, atol=1.0e-12)
    np.testing.assert_allclose(physical_bulk_state.F_d, physical_bulk_state_again.F_d, atol=1.0e-12)
    assert not np.allclose(physical_bulk_state.c_Q, reference_bulk_state.c_Q)
    assert not np.allclose(physical_bulk_state.c_u, reference_bulk_state.c_u)
    assert not np.allclose(physical_bulk_state.c_d, reference_bulk_state.c_d)
    assert not np.allclose(physical_bulk_state.F_Q, reference_bulk_state.F_Q)
    assert not np.allclose(physical_bulk_state.F_u, reference_bulk_state.F_u)
    assert not np.allclose(physical_bulk_state.F_d, reference_bulk_state.F_d)
    assert physical_point.metadata["benchmark_reference_kind"] == "physical_qs1_seed_to_profile"
    assert physical_point.metadata["physical_contract_schema_id"] == PAPER_0710_1869_PHYSICAL_SEED_TO_PROFILE_CONTRACT_SCHEMA_ID
    assert physical_point.metadata["physical_mapping_policy_id"] == PAPER_0710_1869_SEED_TO_PROFILE_MAPPING_POLICY_ID
    assert physical_point.metadata["physical_universal_term_policy_id"] == PAPER_0710_1869_UNIVERSAL_TERM_COEFFICIENT_POLICY_ID


def test_physical_point_builder_rejects_non_frozen_contract_variants():
    benchmark = default_paper_0710_1869_pr1_benchmark()
    contract = default_paper_0710_1869_physical_seed_to_profile_contract()
    mutated_contract = replace(
        contract,
        mapping_policy=replace(contract.mapping_policy, notes="mutated mapping policy"),
    )

    with pytest.raises(ValueError, match="physical_contract.mapping_policy must match the exact frozen QS1 default"):
        build_paper_0710_1869_seeded_physical_point(
            benchmark,
            _seed(),
            physical_contract=mutated_contract,
        )


def test_physical_bulk_state_rejects_forged_payload_and_mutated_nested_contracts():
    physical_point = _physical_point()
    physical_state = derive_paper_0710_1869_physical_bulk_state(physical_point)
    mutated_contract = replace(
        physical_state.physical_contract,
        universal_term_policy=replace(
            physical_state.physical_contract.universal_term_policy,
            notes="mutated universal-term policy",
        ),
    )

    with pytest.raises(ValueError, match="point-derived physical payload"):
        replace(
            physical_state,
            c_Q=physical_state.c_Q + np.array([1.0e-8, 0.0, 0.0]),
        )

    with pytest.raises(ValueError, match="physical bulk-state contract must match point.physical_contract"):
        replace(physical_state, physical_contract=mutated_contract)


def test_exact_degenerate_physical_bulk_state_canonicalization_remains_deterministic():
    physical_point = _degenerate_physical_point()
    first = derive_paper_0710_1869_physical_bulk_state(physical_point)
    second = derive_paper_0710_1869_physical_bulk_state(physical_point)

    np.testing.assert_allclose(first.eig_Q[:2], first.eig_Q[0], atol=1.0e-12)
    np.testing.assert_allclose(first.eig_u[:2], first.eig_u[0], atol=1.0e-12)
    np.testing.assert_allclose(first.eig_d[:2], first.eig_d[0], atol=1.0e-12)
    np.testing.assert_allclose(first.c_Q, second.c_Q, atol=1.0e-12)
    np.testing.assert_allclose(first.c_u, second.c_u, atol=1.0e-12)
    np.testing.assert_allclose(first.c_d, second.c_d, atol=1.0e-12)
    np.testing.assert_allclose(first.F_Q, second.F_Q, atol=1.0e-12)
    np.testing.assert_allclose(first.F_u, second.F_u, atol=1.0e-12)
    np.testing.assert_allclose(first.F_d, second.F_d, atol=1.0e-12)
    np.testing.assert_allclose(first.rotation_Q, second.rotation_Q, atol=1.0e-12)
    np.testing.assert_allclose(first.rotation_u, second.rotation_u, atol=1.0e-12)
    np.testing.assert_allclose(first.rotation_d, second.rotation_d, atol=1.0e-12)
    np.testing.assert_allclose(first.Y_u_bulk_basis, second.Y_u_bulk_basis, atol=1.0e-12)
    np.testing.assert_allclose(first.Y_d_bulk_basis, second.Y_d_bulk_basis, atol=1.0e-12)


def test_near_degenerate_physical_eigensystems_are_rejected_instead_of_using_runtime_dependent_basis():
    rotation = np.array(
        [
            [1.0 / np.sqrt(2.0), 1.0 / np.sqrt(2.0), 0.0],
            [-1.0 / np.sqrt(2.0), 1.0 / np.sqrt(2.0), 0.0],
            [0.0, 0.0, 1.0],
        ],
        dtype=np.complex128,
    )
    eigenvalues = np.diag([1.0e6, 1.0e6 + 5.0e-7, 2.0e6]).astype(np.complex128)
    matrix = rotation @ eigenvalues @ rotation.conjugate().T

    with pytest.raises(ValueError, match="near-degenerate Hermitian spectra are not supported"):
        _ordered_hermitian_spectrum(matrix)


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


def test_eq3_residual_summary_rejects_nonunitary_raw_v5km_arrays():
    with pytest.raises(ValueError, match="unitary"):
        evaluate_eq3_diagonal_consistency(
            c_q_eigenvalues=(1.0, 2.0, 3.0),
            c_u_eigenvalues=(9.0, 8.0, 7.0),
            c_d_eigenvalues=(1.0, 2.0, 3.0),
            a=1.0,
            r=0.0,
            v5km=np.diag([2.0, 1.0, 1.0]).astype(np.complex128),
        )


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
