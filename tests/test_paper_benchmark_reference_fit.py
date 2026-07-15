from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

import quarkConstraints.paper_0710_1869 as paper_pkg
from quarkConstraints.paper_0710_1869.benchmarks import (
    Paper07101869BenchmarkSpurionSeed,
    build_paper_0710_1869_benchmark_point,
    build_paper_0710_1869_seeded_physical_point,
    build_paper_0710_1869_seeded_structural_reference_point,
    default_paper_0710_1869_pr1_benchmark,
)
from quarkConstraints.paper_0710_1869.fit import (
    PAPER_0710_1869_BENCHMARK_REFERENCE_ATTACHMENT_CONVENTION_ID,
    PAPER_0710_1869_BENCHMARK_REFERENCE_CONTRACT_ID,
    PAPER_0710_1869_BENCHMARK_REFERENCE_PROFILE_SOURCE_ID,
    PAPER_0710_1869_BENCHMARK_REFERENCE_RESULT_SCHEMA_ID,
    PAPER_0710_1869_BENCHMARK_REFERENCE_TARGETS_SCHEMA_ID,
    Paper07101869BenchmarkReferenceMassProbeResult,
    Paper07101869BenchmarkReferenceTargets,
    build_paper_0710_1869_benchmark_reference_mass_matrices,
    evaluate_paper_0710_1869_benchmark_reference_mass_probe,
    paper_benchmark_reference_mass_probe_observables,
    paper_benchmark_reference_mass_probe_residuals,
)
from quarkConstraints.paper_0710_1869.model import (
    PAPER_0710_1869_BENCHMARK_SEED_CONSTRUCTION_ID,
    PAPER_0710_1869_DERIVED_PROFILE_INPUT_POLICY_ID,
    PAPER_0710_1869_NONCANONICAL_CLAIM_LEVEL_ID,
    PAPER_0710_1869_NO_BULK_MASS_MAP_POLICY_ID,
    PAPER_0710_1869_PHYSICAL_BENCHMARK_SEED_CONSTRUCTION_ID,
    PAPER_0710_1869_PHYSICAL_BULK_STATE_SCHEMA_ID,
    PAPER_0710_1869_PHYSICAL_CLAIM_LEVEL_ID,
    PAPER_0710_1869_PHYSICAL_CORE_CONTRACT_ID,
    PAPER_0710_1869_PHYSICAL_POINT_KIND_ID,
    PAPER_0710_1869_PHYSICAL_POINT_SCHEMA_ID,
    PAPER_0710_1869_POINT_DERIVED_PHYSICAL_STATUS_ID,
    PAPER_0710_1869_QUOTED_PROFILE_INPUT_POLICY_ID,
    PAPER_0710_1869_REFERENCE_ONLY_PROFILE_INPUT_POLICY_ID,
    PAPER_0710_1869_STRUCTURAL_CORE_CONTRACT_ID,
    PAPER_0710_1869_STRUCTURAL_REFERENCE_STATUS_ID,
    Paper07101869BenchmarkSector,
    Paper07101869RotationParameters,
    derive_paper_0710_1869_physical_bulk_state,
    derive_paper_0710_1869_bulk_state,
)
from quarkConstraints.paper_0710_1869.conventions import (
    PAPER_0710_1869_PROFILE_DERIVATION_POLICY_ID,
    PAPER_0710_1869_SEED_TO_PROFILE_MAPPING_POLICY_ID,
    PAPER_0710_1869_UNIVERSAL_TERM_COEFFICIENT_POLICY_ID,
)
from quarkConstraints.paper_0710_1869.inputs import (
    PAPER_0710_1869_AFFINE_BULK_MASS_LEADING_TERM_COEFFICIENT,
    PAPER_0710_1869_AFFINE_BULK_MASS_UNIVERSAL_OFFSET,
    PAPER_0710_1869_PHYSICAL_SEED_TO_PROFILE_CONTRACT_SCHEMA_ID,
    default_paper_0710_1869_physical_seed_to_profile_contract,
)


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


@pytest.fixture
def seeded_point():
    return build_paper_0710_1869_seeded_structural_reference_point(
        default_paper_0710_1869_pr1_benchmark(),
        _seed(),
        metadata={"test_case": "benchmark_reference_acceptance"},
    )


@pytest.fixture
def bulk_state(seeded_point):
    return derive_paper_0710_1869_bulk_state(seeded_point)


@pytest.fixture
def probe_result(seeded_point):
    return evaluate_paper_0710_1869_benchmark_reference_mass_probe(seeded_point)


def _off_diagonal(matrix: np.ndarray) -> np.ndarray:
    return matrix - np.diag(np.diag(matrix))


def _permute_reference_pairs(
    sector: Paper07101869BenchmarkSector,
    order: tuple[int, int, int],
) -> Paper07101869BenchmarkSector:
    c_vector = np.asarray(sector.c_eigenvalues, dtype=float)
    f_vector = np.asarray(sector.f_eigenvalues, dtype=float)
    return Paper07101869BenchmarkSector(
        label=sector.label,
        c_eigenvalues=tuple(c_vector[list(order)]),
        f_eigenvalues=tuple(f_vector[list(order)]),
    )


def test_seeded_points_and_bulk_states_carry_structural_only_markers(
    seeded_point,
    bulk_state,
):
    assert seeded_point.claim_level_id == PAPER_0710_1869_NONCANONICAL_CLAIM_LEVEL_ID
    assert seeded_point.bulk_mass_map_policy_id == PAPER_0710_1869_NO_BULK_MASS_MAP_POLICY_ID
    assert seeded_point.profile_input_policy_id == PAPER_0710_1869_REFERENCE_ONLY_PROFILE_INPUT_POLICY_ID
    assert seeded_point.structural_reference_status_id == PAPER_0710_1869_STRUCTURAL_REFERENCE_STATUS_ID
    assert seeded_point.construction_id == PAPER_0710_1869_BENCHMARK_SEED_CONSTRUCTION_ID
    assert seeded_point.metadata["benchmark_status"] == "sourced_structural_only"
    assert seeded_point.metadata["benchmark_reference_kind"] == "structural_only"
    assert seeded_point.metadata["test_case"] == "benchmark_reference_acceptance"

    assert bulk_state.core_contract_id == PAPER_0710_1869_STRUCTURAL_CORE_CONTRACT_ID
    assert bulk_state.claim_level_id == seeded_point.claim_level_id
    assert bulk_state.bulk_mass_map_policy_id == seeded_point.bulk_mass_map_policy_id
    assert bulk_state.structural_reference_status_id == seeded_point.structural_reference_status_id
    assert bulk_state.profile_consistency_Q.policy_id == PAPER_0710_1869_REFERENCE_ONLY_PROFILE_INPUT_POLICY_ID
    assert bulk_state.profile_consistency_u.policy_id == PAPER_0710_1869_REFERENCE_ONLY_PROFILE_INPUT_POLICY_ID
    assert bulk_state.profile_consistency_d.policy_id == PAPER_0710_1869_REFERENCE_ONLY_PROFILE_INPUT_POLICY_ID
    assert bulk_state.reference_eq3_summary is not None


@pytest.mark.parametrize(
    "stale_policy_id",
    [
        PAPER_0710_1869_QUOTED_PROFILE_INPUT_POLICY_ID,
        PAPER_0710_1869_DERIVED_PROFILE_INPUT_POLICY_ID,
    ],
)
def test_seeded_points_and_benchmark_builders_reject_stale_profile_policy_ids(
    seeded_point,
    stale_policy_id,
):
    with pytest.raises(ValueError, match="profile_input_policy_id"):
        replace(seeded_point, profile_input_policy_id=stale_policy_id)

    with pytest.raises(ValueError, match="diagnostics-only Table I reference profile policy"):
        build_paper_0710_1869_seeded_structural_reference_point(
            default_paper_0710_1869_pr1_benchmark(),
            _seed(),
            profile_input_policy_id=stale_policy_id,
        )


def test_benchmark_builder_is_honest_about_seeded_structural_reference_output():
    point = build_paper_0710_1869_benchmark_point(
        default_paper_0710_1869_pr1_benchmark(),
        _seed(),
    )

    assert point.construction_id == PAPER_0710_1869_BENCHMARK_SEED_CONSTRUCTION_ID
    assert point.claim_level_id == PAPER_0710_1869_NONCANONICAL_CLAIM_LEVEL_ID
    assert point.bulk_mass_map_policy_id == PAPER_0710_1869_NO_BULK_MASS_MAP_POLICY_ID
    assert point.structural_reference_status_id == PAPER_0710_1869_STRUCTURAL_REFERENCE_STATUS_ID
    assert "custom seeded structural reference point" in (point.notes or "")
    assert "not a sourced exact paper benchmark" in (point.notes or "")
    assert "both left- and right-handed flavor rotations" in (point.notes or "")
    assert point.metadata["benchmark_id"] == "pr1.table_i_eq3_example.v1"
    assert point.metadata["benchmark_status"] == "sourced_structural_only"
    assert point.metadata["benchmark_reference_kind"] == "structural_only"


def test_structural_benchmark_builder_rejects_conflicting_canonical_metadata_and_keeps_extras():
    benchmark = default_paper_0710_1869_pr1_benchmark()

    with pytest.raises(
        ValueError,
        match="metadata\\['benchmark_reference_kind'\\] must not override the frozen canonical provenance value",
    ):
        build_paper_0710_1869_seeded_structural_reference_point(
            benchmark,
            _seed(),
            metadata={"benchmark_reference_kind": "forged_override"},
        )

    point = build_paper_0710_1869_seeded_structural_reference_point(
        benchmark,
        _seed(),
        metadata={
            "benchmark_reference_kind": "structural_only",
            "test_case": "structural_metadata_acceptance",
        },
    )

    assert point.metadata["benchmark_reference_kind"] == "structural_only"
    assert point.metadata["test_case"] == "structural_metadata_acceptance"


def test_physical_benchmark_builder_uses_exact_frozen_contract_and_stays_separate_from_reference_only_path():
    benchmark = default_paper_0710_1869_pr1_benchmark()
    physical_point = build_paper_0710_1869_seeded_physical_point(
        benchmark,
        _seed(),
        metadata={"test_case": "physical_benchmark_acceptance"},
    )
    physical_state = derive_paper_0710_1869_physical_bulk_state(physical_point)
    structural_state = derive_paper_0710_1869_bulk_state(
        build_paper_0710_1869_seeded_structural_reference_point(
            benchmark,
            _seed(),
            metadata={"test_case": "structural_benchmark_acceptance"},
        )
    )

    assert physical_point.construction_id == PAPER_0710_1869_PHYSICAL_BENCHMARK_SEED_CONSTRUCTION_ID
    assert physical_point.schema_id == PAPER_0710_1869_PHYSICAL_POINT_SCHEMA_ID
    assert physical_state.schema_id == PAPER_0710_1869_PHYSICAL_BULK_STATE_SCHEMA_ID
    assert physical_point.point_kind_id == PAPER_0710_1869_PHYSICAL_POINT_KIND_ID
    assert physical_point.claim_level_id == PAPER_0710_1869_PHYSICAL_CLAIM_LEVEL_ID
    assert physical_point.profile_input_policy_id == PAPER_0710_1869_DERIVED_PROFILE_INPUT_POLICY_ID
    assert physical_point.physical_profile_status_id == PAPER_0710_1869_POINT_DERIVED_PHYSICAL_STATUS_ID
    assert physical_point.bulk_mass_map_policy_id == PAPER_0710_1869_SEED_TO_PROFILE_MAPPING_POLICY_ID
    assert physical_point.physical_contract == default_paper_0710_1869_physical_seed_to_profile_contract()
    assert physical_point.metadata["benchmark_reference_kind"] == "physical_qs1_seed_to_profile"
    assert physical_point.metadata["physical_contract_schema_id"] == PAPER_0710_1869_PHYSICAL_SEED_TO_PROFILE_CONTRACT_SCHEMA_ID
    assert physical_point.metadata["physical_mapping_policy_id"] == PAPER_0710_1869_SEED_TO_PROFILE_MAPPING_POLICY_ID
    assert physical_point.metadata["physical_universal_term_policy_id"] == PAPER_0710_1869_UNIVERSAL_TERM_COEFFICIENT_POLICY_ID
    assert physical_state.core_contract_id == PAPER_0710_1869_PHYSICAL_CORE_CONTRACT_ID
    assert physical_state.point_kind_id == PAPER_0710_1869_PHYSICAL_POINT_KIND_ID
    assert physical_state.claim_level_id == PAPER_0710_1869_PHYSICAL_CLAIM_LEVEL_ID
    assert physical_state.profile_input_policy_id == PAPER_0710_1869_DERIVED_PROFILE_INPUT_POLICY_ID
    assert physical_state.physical_profile_status_id == PAPER_0710_1869_POINT_DERIVED_PHYSICAL_STATUS_ID
    affine_policies = {
        policy.sector_id: policy
        for policy in physical_point.physical_contract.universal_term_policy.sector_policies
    }
    for sector_id, c_values, eig_values in (
        ("Q", physical_state.c_Q, physical_state.eig_Q),
        ("u", physical_state.c_u, physical_state.eig_u),
        ("d", physical_state.c_d, physical_state.eig_d),
    ):
        policy = affine_policies[sector_id]
        assert policy.leading_term_coefficient == (
            PAPER_0710_1869_AFFINE_BULK_MASS_LEADING_TERM_COEFFICIENT
        )
        assert policy.universal_offset == PAPER_0710_1869_AFFINE_BULK_MASS_UNIVERSAL_OFFSET
        # C-6: physical seeded profiles use the negative-slope affine map c = -lambda + 0.6.
        np.testing.assert_allclose(
            c_values,
            policy.leading_term_coefficient * eig_values + policy.universal_offset,
            atol=1.0e-12,
        )
    assert not np.allclose(physical_state.c_Q, structural_state.c_Q)
    assert not np.allclose(physical_state.c_u, structural_state.c_u)
    assert not np.allclose(physical_state.c_d, structural_state.c_d)
    assert not np.allclose(physical_state.F_Q, structural_state.F_Q)
    assert not np.allclose(physical_state.F_u, structural_state.F_u)
    assert not np.allclose(physical_state.F_d, structural_state.F_d)


def test_physical_benchmark_builder_rejects_conflicting_canonical_metadata_and_keeps_extras():
    benchmark = default_paper_0710_1869_pr1_benchmark()

    with pytest.raises(
        ValueError,
        match="metadata\\['physical_mapping_policy_id'\\] must not override the frozen canonical provenance value",
    ):
        build_paper_0710_1869_seeded_physical_point(
            benchmark,
            _seed(),
            metadata={"physical_mapping_policy_id": "forged_override"},
        )

    point = build_paper_0710_1869_seeded_physical_point(
        benchmark,
        _seed(),
        metadata={
            "physical_mapping_policy_id": PAPER_0710_1869_SEED_TO_PROFILE_MAPPING_POLICY_ID,
            "test_case": "physical_metadata_acceptance",
        },
    )

    assert point.metadata["physical_mapping_policy_id"] == PAPER_0710_1869_SEED_TO_PROFILE_MAPPING_POLICY_ID
    assert point.metadata["test_case"] == "physical_metadata_acceptance"


def test_physical_benchmark_builder_rejects_mutated_frozen_contract():
    benchmark = default_paper_0710_1869_pr1_benchmark()
    contract = default_paper_0710_1869_physical_seed_to_profile_contract()
    mutated_contract = replace(
        contract,
        mapping_policy=replace(
            contract.mapping_policy,
            notes="mutated seed-to-profile mapping policy",
        ),
    )

    with pytest.raises(
        ValueError,
        match="physical_contract.mapping_policy must match the exact frozen QS1 default",
    ):
        build_paper_0710_1869_seeded_physical_point(
            benchmark,
            _seed(),
            physical_contract=mutated_contract,
        )


def test_reference_only_fit_surface_rejects_qs1_physical_point_and_state():
    benchmark = default_paper_0710_1869_pr1_benchmark()
    physical_point = build_paper_0710_1869_seeded_physical_point(
        benchmark,
        _seed(),
    )
    physical_state = derive_paper_0710_1869_physical_bulk_state(physical_point)

    with pytest.raises(ValueError, match="point must be a Paper07101869Point"):
        evaluate_paper_0710_1869_benchmark_reference_mass_probe(physical_point)

    with pytest.raises(ValueError, match="bulk_state must be a Paper07101869BulkState"):
        build_paper_0710_1869_benchmark_reference_mass_matrices(physical_state)


def test_bulk_state_diagonalization_and_bulk_basis_invariants(seeded_point, bulk_state):
    np.testing.assert_allclose(bulk_state.C_Q, bulk_state.C_Q.conjugate().T, atol=1.0e-12)
    np.testing.assert_allclose(bulk_state.C_u, bulk_state.C_u.conjugate().T, atol=1.0e-12)
    np.testing.assert_allclose(bulk_state.C_d, bulk_state.C_d.conjugate().T, atol=1.0e-12)

    diagonalized_q = bulk_state.rotation_Q.conjugate().T @ bulk_state.C_Q @ bulk_state.rotation_Q
    diagonalized_u = bulk_state.rotation_u.conjugate().T @ bulk_state.C_u @ bulk_state.rotation_u
    diagonalized_d = bulk_state.rotation_d.conjugate().T @ bulk_state.C_d @ bulk_state.rotation_d

    np.testing.assert_allclose(_off_diagonal(diagonalized_q), 0.0, atol=1.0e-10)
    np.testing.assert_allclose(_off_diagonal(diagonalized_u), 0.0, atol=1.0e-10)
    np.testing.assert_allclose(_off_diagonal(diagonalized_d), 0.0, atol=1.0e-10)
    np.testing.assert_allclose(np.diag(diagonalized_q).real, bulk_state.eig_Q, atol=1.0e-10)
    np.testing.assert_allclose(np.diag(diagonalized_u).real, bulk_state.eig_u, atol=1.0e-10)
    np.testing.assert_allclose(np.diag(diagonalized_d).real, bulk_state.eig_d, atol=1.0e-10)

    np.testing.assert_allclose(
        bulk_state.Y_u_bulk_basis,
        bulk_state.rotation_Q.conjugate().T @ seeded_point.Y_u @ bulk_state.rotation_u,
        atol=1.0e-12,
    )
    np.testing.assert_allclose(
        bulk_state.Y_d_bulk_basis,
        bulk_state.rotation_Q.conjugate().T @ seeded_point.Y_d @ bulk_state.rotation_d,
        atol=1.0e-12,
    )


def test_structural_core_keeps_table_i_aliases_as_reference_data_only(bulk_state):
    np.testing.assert_allclose(
        bulk_state.reference_diagnostics_Q.quoted_c,
        bulk_state.c_Q,
        atol=1.0e-12,
    )
    np.testing.assert_allclose(
        bulk_state.reference_diagnostics_u.quoted_c,
        bulk_state.c_u,
        atol=1.0e-12,
    )
    np.testing.assert_allclose(
        bulk_state.reference_diagnostics_d.quoted_c,
        bulk_state.c_d,
        atol=1.0e-12,
    )
    np.testing.assert_allclose(
        bulk_state.reference_diagnostics_u.structural_eigenvalues,
        bulk_state.eig_u,
        atol=1.0e-12,
    )
    assert np.all(bulk_state.eig_u >= -1.0e-12)
    assert bulk_state.c_u[2] < 0.0
    assert not np.allclose(bulk_state.eig_u, bulk_state.c_u, atol=1.0e-6)
    assert bulk_state.reference_diagnostics_u.max_abs_delta_c > 1.0e-3

    np.testing.assert_allclose(bulk_state.profile_consistency_Q.quoted_f, bulk_state.F_Q, atol=1.0e-12)
    np.testing.assert_allclose(bulk_state.profile_consistency_u.quoted_f, bulk_state.F_u, atol=1.0e-12)
    np.testing.assert_allclose(bulk_state.profile_consistency_d.quoted_f, bulk_state.F_d, atol=1.0e-12)
    np.testing.assert_allclose(
        bulk_state.reference_diagnostics_Q.geometry_derived_f_from_quoted_c,
        bulk_state.derived_F_Q,
        atol=1.0e-12,
    )
    np.testing.assert_allclose(
        bulk_state.reference_diagnostics_u.geometry_derived_f_from_quoted_c,
        bulk_state.derived_F_u,
        atol=1.0e-12,
    )
    np.testing.assert_allclose(
        bulk_state.reference_diagnostics_d.geometry_derived_f_from_quoted_c,
        bulk_state.derived_F_d,
        atol=1.0e-12,
    )


def test_benchmark_reference_mass_probe_uses_attached_reference_profile_aliases(bulk_state):
    M_u, M_d = build_paper_0710_1869_benchmark_reference_mass_matrices(bulk_state)
    attached_f_q = np.array([2.0e-3, 1.0e-2, 2.0e-1], dtype=float)
    attached_f_u = np.array([7.0e-4, 6.0e-2, 8.0e-1], dtype=float)
    attached_f_d = np.array([2.0e-3, 8.0e-3, 2.0e-2], dtype=float)
    reference_M_u = (
        2.0
        * bulk_state.point.v
        * np.diag(attached_f_q)
        @ bulk_state.Y_u_bulk_basis
        @ np.diag(attached_f_u)
    )
    reference_M_d = (
        2.0
        * bulk_state.point.v
        * np.diag(attached_f_q)
        @ bulk_state.Y_d_bulk_basis
        @ np.diag(attached_f_d)
    )
    derived_M_u = (
        2.0
        * bulk_state.point.v
        * np.diag(bulk_state.derived_F_Q)
        @ bulk_state.Y_u_bulk_basis
        @ np.diag(bulk_state.derived_F_u)
    )

    np.testing.assert_allclose(M_u, reference_M_u, atol=1.0e-12)
    np.testing.assert_allclose(M_d, reference_M_d, atol=1.0e-12)
    assert np.max(np.abs(M_u - derived_M_u)) > 1.0e-8


def test_benchmark_reference_mass_probe_returns_full_diagonalization_outputs(
    bulk_state,
    probe_result,
):
    M_u, M_d = build_paper_0710_1869_benchmark_reference_mass_matrices(bulk_state)
    observables = paper_benchmark_reference_mass_probe_observables(M_u, M_d)

    np.testing.assert_allclose(probe_result.M_u, M_u, atol=1.0e-12)
    np.testing.assert_allclose(probe_result.M_d, M_d, atol=1.0e-12)
    np.testing.assert_allclose(probe_result.masses_up, observables["masses_up"], atol=1.0e-12)
    np.testing.assert_allclose(
        probe_result.masses_down,
        observables["masses_down"],
        atol=1.0e-12,
    )
    np.testing.assert_allclose(probe_result.ckm, observables["ckm"], atol=1.0e-12)
    np.testing.assert_allclose(
        probe_result.U_L_u.conjugate().T @ probe_result.M_u @ probe_result.U_R_u,
        np.diag(probe_result.masses_up),
        atol=1.0e-10,
    )
    np.testing.assert_allclose(
        probe_result.U_L_d.conjugate().T @ probe_result.M_d @ probe_result.U_R_d,
        np.diag(probe_result.masses_down),
        atol=1.0e-10,
    )
    np.testing.assert_allclose(
        probe_result.ckm,
        probe_result.U_L_u.conjugate().T @ probe_result.U_L_d,
        atol=1.0e-12,
    )
    np.testing.assert_allclose(
        probe_result.U_L_u.conjugate().T @ probe_result.U_L_u,
        np.eye(3),
        atol=1.0e-12,
    )
    np.testing.assert_allclose(
        probe_result.U_L_d.conjugate().T @ probe_result.U_L_d,
        np.eye(3),
        atol=1.0e-12,
    )


def test_probe_attachment_convention_is_machine_readable_and_permutation_stable(
    seeded_point,
    probe_result,
):
    targets = Paper07101869BenchmarkReferenceTargets(
        up_masses=probe_result.masses_up,
        down_masses=probe_result.masses_down,
        ckm=probe_result.ckm,
    )
    assert (
        targets.attachment_convention_id
        == PAPER_0710_1869_BENCHMARK_REFERENCE_ATTACHMENT_CONVENTION_ID
    )
    assert (
        probe_result.attachment_convention_id
        == PAPER_0710_1869_BENCHMARK_REFERENCE_ATTACHMENT_CONVENTION_ID
    )

    permuted_point = replace(
        seeded_point,
        q_sector=_permute_reference_pairs(seeded_point.q_sector, (1, 2, 0)),
        u_sector=_permute_reference_pairs(seeded_point.u_sector, (2, 0, 1)),
        d_sector=_permute_reference_pairs(seeded_point.d_sector, (2, 1, 0)),
    )
    permuted_result = evaluate_paper_0710_1869_benchmark_reference_mass_probe(permuted_point)

    assert (
        permuted_result.attachment_convention_id
        == PAPER_0710_1869_BENCHMARK_REFERENCE_ATTACHMENT_CONVENTION_ID
    )
    np.testing.assert_allclose(permuted_result.M_u, probe_result.M_u, atol=1.0e-12)
    np.testing.assert_allclose(permuted_result.M_d, probe_result.M_d, atol=1.0e-12)
    np.testing.assert_allclose(permuted_result.masses_up, probe_result.masses_up, atol=1.0e-12)
    np.testing.assert_allclose(
        permuted_result.masses_down,
        probe_result.masses_down,
        atol=1.0e-12,
    )
    np.testing.assert_allclose(permuted_result.ckm, probe_result.ckm, atol=1.0e-12)


def test_benchmark_reference_targets_reject_unsorted_masses_and_nonunitary_ckm(probe_result):
    with pytest.raises(ValueError, match="strictly ascending"):
        Paper07101869BenchmarkReferenceTargets(
            up_masses=(probe_result.masses_up[2], probe_result.masses_up[1], probe_result.masses_up[0]),
            down_masses=probe_result.masses_down,
            ckm=probe_result.ckm,
        )

    bad_ckm = probe_result.ckm.copy()
    bad_ckm[0, 0] += 0.05
    with pytest.raises(ValueError, match="unitary"):
        Paper07101869BenchmarkReferenceTargets(
            up_masses=probe_result.masses_up,
            down_masses=probe_result.masses_down,
            ckm=bad_ckm,
        )


def test_benchmark_reference_probe_self_target_has_zero_score_and_validates_invariants(
    seeded_point,
    probe_result,
):
    targets = Paper07101869BenchmarkReferenceTargets(
        up_masses=probe_result.masses_up,
        down_masses=probe_result.masses_down,
        ckm=probe_result.ckm,
        label="self_target",
    )
    result = evaluate_paper_0710_1869_benchmark_reference_mass_probe(seeded_point, targets)

    assert isinstance(result, Paper07101869BenchmarkReferenceMassProbeResult)
    assert result.schema_id == PAPER_0710_1869_BENCHMARK_REFERENCE_RESULT_SCHEMA_ID
    assert result.target_schema_id == PAPER_0710_1869_BENCHMARK_REFERENCE_TARGETS_SCHEMA_ID
    assert result.target_label == "self_target"
    assert result.probe_contract_id == PAPER_0710_1869_BENCHMARK_REFERENCE_CONTRACT_ID
    assert result.profile_source_id == PAPER_0710_1869_BENCHMARK_REFERENCE_PROFILE_SOURCE_ID
    assert (
        result.attachment_convention_id
        == PAPER_0710_1869_BENCHMARK_REFERENCE_ATTACHMENT_CONVENTION_ID
    )
    assert result.total_score == pytest.approx(0.0, abs=1.0e-12)
    np.testing.assert_allclose(result.mass_residuals_up, 0.0, atol=1.0e-12)
    np.testing.assert_allclose(result.mass_residuals_down, 0.0, atol=1.0e-12)
    np.testing.assert_allclose(result.ckm_residuals, 0.0, atol=1.0e-12)
    assert result.max_ckm_unitarity_residual == pytest.approx(0.0, abs=1.0e-12)
    assert result.max_rotation_unitarity_residual == pytest.approx(0.0, abs=1.0e-12)
    assert result.max_reconstruction_residual == pytest.approx(0.0, abs=1.0e-10)
    assert "benchmark-reference mass probe" in result.summary()
    assert "not a seed-derived physical bulk-profile fit" in result.summary()


def test_benchmark_reference_result_rejects_nonunitary_ckm_and_bad_reconstruction(probe_result):
    bad_ckm = probe_result.ckm.copy()
    bad_ckm[0, 0] += 0.05
    with pytest.raises(ValueError, match="ckm must be unitary"):
        replace(probe_result, ckm=bad_ckm)

    bad_matrix = probe_result.M_u.copy()
    bad_matrix[0, 0] += 1.0
    with pytest.raises(ValueError, match="SVD reconstruction"):
        replace(probe_result, M_u=bad_matrix)


def test_benchmark_reference_residual_helper_rejects_nonunitary_ckm(probe_result):
    targets = Paper07101869BenchmarkReferenceTargets(
        up_masses=probe_result.masses_up,
        down_masses=probe_result.masses_down,
        ckm=probe_result.ckm,
    )
    bad_ckm = probe_result.ckm.copy()
    bad_ckm[0, 0] += 0.05

    with pytest.raises(ValueError, match="ckm must be unitary"):
        paper_benchmark_reference_mass_probe_residuals(
            probe_result.masses_up,
            probe_result.masses_down,
            bad_ckm,
            targets,
        )


def test_benchmark_reference_result_rejects_unitary_ckm_inconsistent_with_left_rotations(
    probe_result,
):
    inconsistent_ckm = probe_result.ckm @ np.diag(
        np.array([np.exp(0.02j), 1.0, 1.0], dtype=np.complex128)
    )

    with pytest.raises(ValueError, match=r"U_L_u\^dagger U_L_d"):
        replace(probe_result, ckm=inconsistent_ckm)


def test_package_surface_exports_only_canonical_benchmark_reference_probe(seeded_point):
    assert not hasattr(paper_pkg, "Paper07101869Targets")
    assert not hasattr(paper_pkg, "Paper07101869FitResult")
    assert not hasattr(paper_pkg, "evaluate_paper_0710_1869_fit")
    assert not hasattr(paper_pkg, "PAPER_0710_1869_QUOTED_PROFILE_INPUT_POLICY_ID")
    assert not hasattr(paper_pkg, "PAPER_0710_1869_DERIVED_PROFILE_INPUT_POLICY_ID")

    assert (
        paper_pkg.Paper07101869BenchmarkReferenceTargets
        is Paper07101869BenchmarkReferenceTargets
    )
    assert (
        paper_pkg.Paper07101869BenchmarkReferenceMassProbeResult
        is Paper07101869BenchmarkReferenceMassProbeResult
    )
    assert (
        paper_pkg.PAPER_0710_1869_BENCHMARK_REFERENCE_TARGETS_SCHEMA_ID
        == PAPER_0710_1869_BENCHMARK_REFERENCE_TARGETS_SCHEMA_ID
    )
    assert (
        paper_pkg.PAPER_0710_1869_BENCHMARK_REFERENCE_RESULT_SCHEMA_ID
        == PAPER_0710_1869_BENCHMARK_REFERENCE_RESULT_SCHEMA_ID
    )
    assert (
        paper_pkg.PAPER_0710_1869_BENCHMARK_REFERENCE_ATTACHMENT_CONVENTION_ID
        == PAPER_0710_1869_BENCHMARK_REFERENCE_ATTACHMENT_CONVENTION_ID
    )
    assert (
        paper_pkg.PAPER_0710_1869_REFERENCE_ONLY_PROFILE_INPUT_POLICY_ID
        == PAPER_0710_1869_REFERENCE_ONLY_PROFILE_INPUT_POLICY_ID
    )

    result = paper_pkg.evaluate_paper_0710_1869_benchmark_reference_mass_probe(seeded_point)
    assert isinstance(result, Paper07101869BenchmarkReferenceMassProbeResult)
    assert result.probe_contract_id == PAPER_0710_1869_BENCHMARK_REFERENCE_CONTRACT_ID
    assert result.profile_source_id == PAPER_0710_1869_BENCHMARK_REFERENCE_PROFILE_SOURCE_ID
    assert (
        result.attachment_convention_id
        == PAPER_0710_1869_BENCHMARK_REFERENCE_ATTACHMENT_CONVENTION_ID
    )
    assert "benchmark-reference mass probe" in result.summary()
