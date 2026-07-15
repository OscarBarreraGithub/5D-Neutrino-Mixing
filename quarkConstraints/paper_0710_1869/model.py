"""Paper-facing model and exact-point helpers for the dedicated 0710.1869 mode."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from warpConfig.baseParams import DEFAULT_LAMBDA_IR, MPL, V_EWSB, get_warp_params
from warpConfig.wavefuncs import f_IR

from .conventions import (
    PAPER_0710_1869_MODE_ID,
    PAPER_0710_1869_PAPER_ID,
    PAPER_0710_1869_PROFILE_DERIVATION_POLICY_ID,
    PAPER_0710_1869_SEED_TO_PROFILE_MAPPING_POLICY_ID,
    PAPER_0710_1869_UNIVERSAL_TERM_COEFFICIENT_POLICY_ID,
)
from .inputs import (
    Paper07101869Eq3Example as Paper07101869Eq3Inputs,
)
from .inputs import (
    Paper07101869AffineBulkMassSectorPolicy,
    Paper07101869PhysicalSeedToProfileContract,
    Paper07101869TableIInputs,
    default_paper_0710_1869_table_i_inputs,
)
from .inputs import (
    default_paper_0710_1869_eq3_example as default_paper_0710_1869_eq3_inputs,
    default_paper_0710_1869_physical_seed_to_profile_contract,
)
from .validation import (
    require_known_schema_id,
    require_member,
    require_nonempty_identifier,
    require_positive_finite,
)

PAPER_0710_1869_MODEL_SCHEMA_ID = "quarkConstraints.paper_0710_1869.model.v1"
PAPER_0710_1869_POINT_SCHEMA_ID = "quarkConstraints.paper_0710_1869.point.v1"
PAPER_0710_1869_BULK_STATE_SCHEMA_ID = "quarkConstraints.paper_0710_1869.bulk_state.v1"
PAPER_0710_1869_PHYSICAL_POINT_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.physical_point.v1"
)
PAPER_0710_1869_PHYSICAL_BULK_STATE_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.physical_bulk_state.v1"
)
PAPER_0710_1869_PROFILE_CONSISTENCY_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.profile_consistency.v1"
)
PAPER_0710_1869_REFERENCE_DIAGNOSTICS_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.reference_diagnostics.v1"
)
PAPER_0710_1869_QUOTED_PROFILE_INPUT_POLICY_ID = "profile_inputs.quoted_table_i.primary.v1"
PAPER_0710_1869_DERIVED_PROFILE_INPUT_POLICY_ID = (
    "profile_inputs.derived_from_c_and_geometry.v1"
)
PAPER_0710_1869_REFERENCE_ONLY_PROFILE_INPUT_POLICY_ID = (
    "profile_inputs.table_i_reference.diagnostics_only.v1"
)
PAPER_0710_1869_NO_BULK_MASS_MAP_POLICY_ID = (
    "quarkConstraints.paper_0710_1869.bulk_mass_map.none.structural_only.v1"
)
PAPER_0710_1869_CUSTOM_POINT_KIND_ID = (
    "quarkConstraints.paper_0710_1869.point_kind.custom_seeded_structural_reference.v1"
)
PAPER_0710_1869_NONCANONICAL_CLAIM_LEVEL_ID = (
    "quarkConstraints.paper_0710_1869.claim.custom_seeded_structural_reference_only.v1"
)
PAPER_0710_1869_STRUCTURAL_REFERENCE_STATUS_ID = (
    "quarkConstraints.paper_0710_1869.structural_reference.table_i_reference_only.v1"
)
PAPER_0710_1869_STRUCTURAL_CORE_CONTRACT_ID = (
    "quarkConstraints.paper_0710_1869.contract.seeded_structural_mfv_core.v1"
)
PAPER_0710_1869_EXPLICIT_SEED_CONSTRUCTION_ID = (
    "quarkConstraints.paper_0710_1869.point_construction.explicit_spurion_seed.v1"
)
PAPER_0710_1869_BENCHMARK_SEED_CONSTRUCTION_ID = (
    "quarkConstraints.paper_0710_1869.point_construction."
    "structural_benchmark_plus_explicit_spurion_seed.v1"
)
PAPER_0710_1869_PHYSICAL_POINT_KIND_ID = (
    "quarkConstraints.paper_0710_1869.point_kind.custom_seeded_physical_qs1.v1"
)
PAPER_0710_1869_PHYSICAL_CLAIM_LEVEL_ID = (
    "quarkConstraints.paper_0710_1869.claim.qs1_seed_to_profile_physical_state_only.v1"
)
PAPER_0710_1869_POINT_DERIVED_PHYSICAL_STATUS_ID = (
    "quarkConstraints.paper_0710_1869.physical_status.point_derived_seed_to_profile.v1"
)
PAPER_0710_1869_PHYSICAL_CORE_CONTRACT_ID = (
    "quarkConstraints.paper_0710_1869.contract.seeded_physical_qs1_core.v1"
)
PAPER_0710_1869_PHYSICAL_SEED_CONSTRUCTION_ID = (
    "quarkConstraints.paper_0710_1869.point_construction.explicit_spurion_seed_physical_qs1.v1"
)
PAPER_0710_1869_PHYSICAL_BENCHMARK_SEED_CONSTRUCTION_ID = (
    "quarkConstraints.paper_0710_1869.point_construction."
    "structural_benchmark_plus_explicit_spurion_seed_physical_qs1.v1"
)
PAPER_0710_1869_EQ3_RELATION_ID = "diag(C_Q)=a*diag(r*V5KM^dagger*C_u*V5KM + C_d)"
PAPER_0710_1869_TABLE_I_REFERENCE = "Table I, arXiv:0710.1869v1"
PAPER_0710_1869_EQ3_REFERENCE = "Eq. (3), arXiv:0710.1869v1"
PAPER_0710_1869_EQ3_UNITARITY_TOLERANCE = 1.0e-10
PAPER_0710_1869_PHYSICAL_UNITARITY_TOLERANCE = 1.0e-10
PAPER_0710_1869_PHYSICAL_STATE_CONSISTENCY_ATOL = 1.0e-12
PAPER_0710_1869_EIGENVALUE_DEGENERACY_ATOL = 1.0e-12
PAPER_0710_1869_EIGENVECTOR_COMPONENT_TOL = 1.0e-12

_ALLOWED_PROFILE_INPUT_POLICY_IDS = (
    PAPER_0710_1869_REFERENCE_ONLY_PROFILE_INPUT_POLICY_ID,
)
_ALLOWED_BULK_MASS_MAP_POLICY_IDS = (PAPER_0710_1869_NO_BULK_MASS_MAP_POLICY_ID,)
_ALLOWED_POINT_KIND_IDS = (PAPER_0710_1869_CUSTOM_POINT_KIND_ID,)
_ALLOWED_CLAIM_LEVEL_IDS = (PAPER_0710_1869_NONCANONICAL_CLAIM_LEVEL_ID,)
_ALLOWED_STRUCTURAL_REFERENCE_STATUS_IDS = (PAPER_0710_1869_STRUCTURAL_REFERENCE_STATUS_ID,)
_ALLOWED_CONSTRUCTION_IDS = (
    PAPER_0710_1869_EXPLICIT_SEED_CONSTRUCTION_ID,
    PAPER_0710_1869_BENCHMARK_SEED_CONSTRUCTION_ID,
)


def _as_real_triplet(
    name: str, values: tuple[float, float, float] | list[float] | np.ndarray
) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.shape != (3,):
        raise ValueError(f"{name} must have shape (3,)")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must contain only finite values")
    return arr.astype(float, copy=True)


def _as_positive_real_triplet(
    name: str, values: tuple[float, float, float] | list[float] | np.ndarray
) -> np.ndarray:
    arr = _as_real_triplet(name, values)
    if np.any(arr <= 0.0):
        raise ValueError(f"{name} must contain only positive values")
    return arr


def _as_nonnegative_real_triplet(
    name: str, values: tuple[float, float, float] | list[float] | np.ndarray
) -> np.ndarray:
    arr = _as_real_triplet(name, values)
    if np.any(arr < 0.0):
        raise ValueError(f"{name} must contain only non-negative values")
    return arr


def _as_complex_matrix3(name: str, values: np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=np.complex128)
    if arr.shape != (3, 3):
        raise ValueError(f"{name} must have shape (3, 3)")
    if not np.all(np.isfinite(arr.real)) or not np.all(np.isfinite(arr.imag)):
        raise ValueError(f"{name} must contain only finite entries")
    return arr


def _unitarity_residual(matrix: np.ndarray) -> float:
    arr = _as_complex_matrix3("matrix", matrix)
    identity = np.eye(arr.shape[0], dtype=np.complex128)
    return float(np.max(np.abs(arr.conjugate().T @ arr - identity)))


def _require_unitary_matrix3(name: str, values: np.ndarray, *, tolerance: float) -> np.ndarray:
    arr = _as_complex_matrix3(name, values)
    residual = _unitarity_residual(arr)
    if residual > tolerance:
        raise ValueError(
            f"{name} must be unitary within tolerance {tolerance:.3e}; "
            f"got residual {residual:.3e}"
        )
    return arr


def _require_finite(name: str, value: float) -> float:
    numeric = float(value)
    if not np.isfinite(numeric):
        raise ValueError(f"{name} must be finite")
    return numeric


def _require_nonnegative_finite(name: str, value: float) -> float:
    numeric = float(value)
    if not np.isfinite(numeric) or numeric < 0.0:
        raise ValueError(f"{name} must be a non-negative finite float")
    return numeric


def _canonicalize_vector_phase(
    vector: np.ndarray, *, tolerance: float = PAPER_0710_1869_EIGENVECTOR_COMPONENT_TOL
) -> np.ndarray:
    arr = np.asarray(vector, dtype=np.complex128).copy()
    for component in arr:
        magnitude = float(abs(component))
        if magnitude > tolerance:
            return arr / (component / magnitude)
    return arr


def _canonical_basis_from_projector(
    projector: np.ndarray,
    *,
    dimension: int,
    tolerance: float = PAPER_0710_1869_EIGENVECTOR_COMPONENT_TOL,
) -> np.ndarray:
    seeds = np.eye(projector.shape[0], dtype=np.complex128)
    basis: list[np.ndarray] = []
    for seed in seeds.T:
        candidate = projector @ seed
        for existing in basis:
            candidate = candidate - existing * np.vdot(existing, candidate)
        norm = float(np.linalg.norm(candidate))
        if norm <= tolerance:
            continue
        basis.append(_canonicalize_vector_phase(candidate / norm, tolerance=tolerance))
        if len(basis) == dimension:
            break
    if len(basis) != dimension:
        raise ValueError("could not construct a canonical basis for the degenerate eigenspace")
    canonical_basis = np.column_stack(basis)
    q_matrix, _ = np.linalg.qr(canonical_basis)
    for column in range(q_matrix.shape[1]):
        q_matrix[:, column] = _canonicalize_vector_phase(
            q_matrix[:, column],
            tolerance=tolerance,
        )
    return q_matrix


def _degenerate_cluster_slices(
    eigenvalues: np.ndarray, *, tolerance: float = PAPER_0710_1869_EIGENVALUE_DEGENERACY_ATOL
) -> list[tuple[int, int]]:
    values = np.asarray(eigenvalues, dtype=float)
    if values.shape != (3,):
        raise ValueError("eigenvalues must have shape (3,) for clustering")
    scale = max(1.0, float(np.max(np.abs(values))))
    threshold = tolerance * scale
    slices: list[tuple[int, int]] = []
    start = 0
    while start < values.size:
        stop = start + 1
        while stop < values.size and abs(values[stop] - values[stop - 1]) <= threshold:
            stop += 1
        slices.append((start, stop))
        start = stop
    return slices


def _max_abs_delta(actual: np.ndarray, expected: np.ndarray) -> float:
    return float(np.max(np.abs(np.asarray(actual) - np.asarray(expected))))


def _require_allclose(
    name: str,
    actual: np.ndarray | float,
    expected: np.ndarray | float,
    *,
    atol: float = PAPER_0710_1869_PHYSICAL_STATE_CONSISTENCY_ATOL,
) -> None:
    actual_arr = np.asarray(actual)
    expected_arr = np.asarray(expected)
    if actual_arr.shape != expected_arr.shape:
        raise ValueError(
            f"{name} must match the point-derived physical payload shape {expected_arr.shape}; "
            f"got {actual_arr.shape}"
        )
    if not np.allclose(actual_arr, expected_arr, atol=atol, rtol=0.0):
        raise ValueError(
            f"{name} must match the point-derived physical payload within tolerance "
            f"{atol:.3e}; got max abs delta {_max_abs_delta(actual_arr, expected_arr):.3e}"
        )


def _ordered_hermitian_spectrum(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    hermitian = 0.5 * (matrix + matrix.conjugate().T)
    values, vectors = np.linalg.eigh(hermitian)
    order = np.argsort(values.real, kind="stable")
    ordered_values = np.real_if_close(values[order]).astype(float)
    ordered_vectors = np.asarray(vectors[:, order], dtype=np.complex128)
    for start, stop in _degenerate_cluster_slices(ordered_values):
        if stop - start == 1:
            ordered_vectors[:, start] = _canonicalize_vector_phase(ordered_vectors[:, start])
            continue
        cluster_values = ordered_values[start:stop]
        cluster_vectors = ordered_vectors[:, start:stop]
        projector = cluster_vectors @ cluster_vectors.conjugate().T
        canonical_cluster = _canonical_basis_from_projector(
            projector,
            dimension=stop - start,
        )
        cluster_matrix = canonical_cluster.conjugate().T @ hermitian @ canonical_cluster
        cluster_expected = np.diag(cluster_values)
        if np.allclose(cluster_matrix, cluster_expected, atol=1.0e-10, rtol=0.0):
            ordered_vectors[:, start:stop] = canonical_cluster
            continue
        raise ValueError(
            "near-degenerate Hermitian spectra are not supported by the frozen paper "
            "contract because the eigensystem basis would be runtime-dependent; "
            f"cluster eigenvalues were {cluster_values.tolist()}"
        )
    ordered_vectors = _require_unitary_matrix3(
        "ordered_eigenvectors",
        ordered_vectors,
        tolerance=PAPER_0710_1869_PHYSICAL_UNITARITY_TOLERANCE,
    )
    return ordered_values, ordered_vectors


def _require_physical_sector_policy(
    *,
    name: str,
    policy: Paper07101869AffineBulkMassSectorPolicy,
    expected_sector_id: str,
) -> Paper07101869AffineBulkMassSectorPolicy:
    if not isinstance(policy, Paper07101869AffineBulkMassSectorPolicy):
        raise ValueError(f"{name} must be a Paper07101869AffineBulkMassSectorPolicy")
    if policy.sector_id != expected_sector_id:
        raise ValueError(f"{name}.sector_id must be exactly {expected_sector_id!r}")
    if policy.leading_term_coefficient >= 0.0:
        raise ValueError(f"{name}.leading_term_coefficient must be negative")
    if not np.isfinite(float(policy.universal_offset)):
        raise ValueError(f"{name}.universal_offset must be finite")
    return policy


def _require_exact_frozen_physical_contract(
    contract: Paper07101869PhysicalSeedToProfileContract,
) -> Paper07101869PhysicalSeedToProfileContract:
    if not isinstance(contract, Paper07101869PhysicalSeedToProfileContract):
        raise ValueError("physical_contract must be a Paper07101869PhysicalSeedToProfileContract")

    if contract.mapping_policy.policy_id != PAPER_0710_1869_SEED_TO_PROFILE_MAPPING_POLICY_ID:
        raise ValueError("physical_contract must use the exact frozen seed-to-profile policy id")
    if (
        contract.universal_term_policy.policy_id
        != PAPER_0710_1869_UNIVERSAL_TERM_COEFFICIENT_POLICY_ID
    ):
        raise ValueError(
            "physical_contract must use the exact frozen universal-term/coefficient policy id"
        )
    if (
        contract.mapping_policy.profile_derivation_policy_id
        != PAPER_0710_1869_PROFILE_DERIVATION_POLICY_ID
    ):
        raise ValueError(
            "physical_contract must use the exact frozen profile-derivation policy id"
        )
    if contract.mapping_policy.uses_hidden_bulk_mass_map_surrogate:
        raise ValueError("physical_contract must not use a hidden BulkMassMap surrogate")

    for index, sector_id in enumerate(("Q", "u", "d")):
        _require_physical_sector_policy(
            name=f"physical_contract.universal_term_policy.sector_policies[{index}]",
            policy=contract.universal_term_policy.sector_policies[index],
            expected_sector_id=sector_id,
        )

    default_mapping_policy = (
        default_paper_0710_1869_physical_seed_to_profile_contract().mapping_policy
    )
    if contract.mapping_policy != default_mapping_policy:
        raise ValueError("physical_contract.mapping_policy must match the exact frozen QS1 default")
    return contract


@dataclass(frozen=True)
class Paper07101869BenchmarkSector:
    """One Table I sector with benchmark bulk-mass and profile eigenvalues."""

    label: str
    c_eigenvalues: tuple[float, float, float]
    f_eigenvalues: tuple[float, float, float]

    def __post_init__(self) -> None:
        object.__setattr__(self, "label", require_nonempty_identifier("label", self.label))
        object.__setattr__(
            self,
            "c_eigenvalues",
            tuple(_as_real_triplet("c_eigenvalues", self.c_eigenvalues)),
        )
        object.__setattr__(
            self,
            "f_eigenvalues",
            tuple(_as_positive_real_triplet("f_eigenvalues", self.f_eigenvalues)),
        )

    @property
    def c_vector(self) -> np.ndarray:
        return _as_real_triplet("c_eigenvalues", self.c_eigenvalues)

    @property
    def f_vector(self) -> np.ndarray:
        return _as_positive_real_triplet("f_eigenvalues", self.f_eigenvalues)

    def build_diagonal_c_matrix(self) -> np.ndarray:
        """Return the 3x3 diagonal bulk-mass matrix for this sector."""
        return np.diag(self.c_vector).astype(np.complex128)


@dataclass(frozen=True)
class Paper07101869TableIBenchmark:
    """Canonical Table I benchmark eigenvalues from the paper."""

    schema_id: str = PAPER_0710_1869_MODEL_SCHEMA_ID
    label: str = "table_i_benchmark"
    reference: str = PAPER_0710_1869_TABLE_I_REFERENCE
    q_sector: Paper07101869BenchmarkSector = field(
        default_factory=lambda: Paper07101869BenchmarkSector(
            label="Q",
            c_eigenvalues=(0.64, 0.59, 0.46),
            f_eigenvalues=(2.0e-3, 1.0e-2, 2.0e-1),
        )
    )
    u_sector: Paper07101869BenchmarkSector = field(
        default_factory=lambda: Paper07101869BenchmarkSector(
            label="u",
            c_eigenvalues=(0.68, 0.53, -0.06),
            f_eigenvalues=(7.0e-4, 6.0e-2, 8.0e-1),
        )
    )
    d_sector: Paper07101869BenchmarkSector = field(
        default_factory=lambda: Paper07101869BenchmarkSector(
            label="d",
            c_eigenvalues=(0.65, 0.60, 0.58),
            f_eigenvalues=(2.0e-3, 8.0e-3, 2.0e-2),
        )
    )

    def __post_init__(self) -> None:
        if self.schema_id != PAPER_0710_1869_MODEL_SCHEMA_ID:
            raise ValueError(f"schema_id must be exactly {PAPER_0710_1869_MODEL_SCHEMA_ID!r}")
        object.__setattr__(self, "label", require_nonempty_identifier("label", self.label))
        object.__setattr__(
            self, "reference", require_nonempty_identifier("reference", self.reference)
        )
        if not isinstance(self.q_sector, Paper07101869BenchmarkSector):
            raise ValueError("q_sector must be a Paper07101869BenchmarkSector")
        if not isinstance(self.u_sector, Paper07101869BenchmarkSector):
            raise ValueError("u_sector must be a Paper07101869BenchmarkSector")
        if not isinstance(self.d_sector, Paper07101869BenchmarkSector):
            raise ValueError("d_sector must be a Paper07101869BenchmarkSector")


@dataclass(frozen=True)
class Paper07101869DiagonalCMatrices:
    """Diagonal paper-side bulk-mass matrices used in the Eq. (3) check."""

    C_Q: np.ndarray
    C_u: np.ndarray
    C_d: np.ndarray

    def __post_init__(self) -> None:
        object.__setattr__(self, "C_Q", _as_complex_matrix3("C_Q", self.C_Q))
        object.__setattr__(self, "C_u", _as_complex_matrix3("C_u", self.C_u))
        object.__setattr__(self, "C_d", _as_complex_matrix3("C_d", self.C_d))
        for name in ("C_Q", "C_u", "C_d"):
            matrix = getattr(self, name)
            diagonal_projection = np.diag(np.diag(matrix))
            if not np.allclose(matrix, diagonal_projection, atol=1.0e-12):
                raise ValueError(f"{name} must be diagonal")
            if not np.allclose(np.diag(matrix).imag, 0.0, atol=1.0e-12):
                raise ValueError(f"{name} diagonal entries must be real")


@dataclass(frozen=True)
class Paper07101869RotationParameters:
    """Paper-owned CKM-like unitary parameterization."""

    theta12: float = 0.0
    theta13: float = 0.0
    theta23: float = 0.0
    delta: float = 0.0

    def __post_init__(self) -> None:
        object.__setattr__(self, "theta12", _require_finite("theta12", self.theta12))
        object.__setattr__(self, "theta13", _require_finite("theta13", self.theta13))
        object.__setattr__(self, "theta23", _require_finite("theta23", self.theta23))
        object.__setattr__(self, "delta", _require_finite("delta", self.delta))

    @classmethod
    def from_degrees(
        cls,
        theta12_deg: float = 0.0,
        theta13_deg: float = 0.0,
        theta23_deg: float = 0.0,
        delta: float = 0.0,
    ) -> "Paper07101869RotationParameters":
        return cls(
            theta12=np.deg2rad(theta12_deg),
            theta13=np.deg2rad(theta13_deg),
            theta23=np.deg2rad(theta23_deg),
            delta=delta,
        )


@dataclass(frozen=True)
class Paper07101869V5KMParameters:
    """CKM-like angle/phase parameterization for the paper-side ``V5KM``."""

    theta12: float = 0.0
    theta23: float = 0.0
    theta13: float = 0.0
    delta: float = 0.0

    def __post_init__(self) -> None:
        object.__setattr__(self, "theta12", _require_finite("theta12", self.theta12))
        object.__setattr__(self, "theta23", _require_finite("theta23", self.theta23))
        object.__setattr__(self, "theta13", _require_finite("theta13", self.theta13))
        object.__setattr__(self, "delta", _require_finite("delta", self.delta))

    @classmethod
    def from_degrees(
        cls,
        *,
        theta12_deg: float,
        theta23_deg: float,
        theta13_deg: float,
        delta: float,
    ) -> "Paper07101869V5KMParameters":
        return cls(
            theta12=np.deg2rad(theta12_deg),
            theta23=np.deg2rad(theta23_deg),
            theta13=np.deg2rad(theta13_deg),
            delta=delta,
        )


@dataclass(frozen=True)
class Paper07101869ProfileConsistency:
    """Diagnostic comparison between quoted and geometry-derived profile overlaps."""

    quoted_f: np.ndarray
    derived_f: np.ndarray
    delta_f: np.ndarray
    max_abs_delta: float
    policy_id: str
    schema_id: str = PAPER_0710_1869_PROFILE_CONSISTENCY_SCHEMA_ID

    def __post_init__(self) -> None:
        object.__setattr__(self, "quoted_f", _as_positive_real_triplet("quoted_f", self.quoted_f))
        object.__setattr__(self, "derived_f", _as_real_triplet("derived_f", self.derived_f))
        object.__setattr__(self, "delta_f", _as_real_triplet("delta_f", self.delta_f))
        object.__setattr__(
            self,
            "max_abs_delta",
            _require_nonnegative_finite("max_abs_delta", self.max_abs_delta),
        )
        object.__setattr__(
            self,
            "policy_id",
            require_member("policy_id", self.policy_id, _ALLOWED_PROFILE_INPUT_POLICY_IDS),
        )
        object.__setattr__(
            self,
            "schema_id",
            require_known_schema_id(
                "schema_id",
                self.schema_id,
                expected=PAPER_0710_1869_PROFILE_CONSISTENCY_SCHEMA_ID,
            ),
        )

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_id": self.schema_id,
            "policy_id": self.policy_id,
            "quoted_f": self.quoted_f.tolist(),
            "derived_f": self.derived_f.tolist(),
            "delta_f": self.delta_f.tolist(),
            "max_abs_delta": self.max_abs_delta,
        }


@dataclass(frozen=True)
class Paper07101869BenchmarkReferenceDiagnostics:
    """Reference-only comparison between seeded structural data and quoted Table I values."""

    quoted_c: np.ndarray
    structural_eigenvalues: np.ndarray
    delta_c: np.ndarray
    max_abs_delta_c: float
    quoted_f: np.ndarray
    geometry_derived_f_from_quoted_c: np.ndarray
    delta_f: np.ndarray
    max_abs_delta_f: float
    schema_id: str = PAPER_0710_1869_REFERENCE_DIAGNOSTICS_SCHEMA_ID

    def __post_init__(self) -> None:
        object.__setattr__(self, "quoted_c", _as_real_triplet("quoted_c", self.quoted_c))
        object.__setattr__(
            self,
            "structural_eigenvalues",
            _as_real_triplet("structural_eigenvalues", self.structural_eigenvalues),
        )
        object.__setattr__(self, "delta_c", _as_real_triplet("delta_c", self.delta_c))
        object.__setattr__(self, "quoted_f", _as_positive_real_triplet("quoted_f", self.quoted_f))
        object.__setattr__(
            self,
            "geometry_derived_f_from_quoted_c",
            _as_real_triplet(
                "geometry_derived_f_from_quoted_c",
                self.geometry_derived_f_from_quoted_c,
            ),
        )
        object.__setattr__(self, "delta_f", _as_real_triplet("delta_f", self.delta_f))
        object.__setattr__(
            self,
            "max_abs_delta_c",
            _require_nonnegative_finite("max_abs_delta_c", self.max_abs_delta_c),
        )
        object.__setattr__(
            self,
            "max_abs_delta_f",
            _require_nonnegative_finite("max_abs_delta_f", self.max_abs_delta_f),
        )
        object.__setattr__(
            self,
            "schema_id",
            require_known_schema_id(
                "schema_id",
                self.schema_id,
                expected=PAPER_0710_1869_REFERENCE_DIAGNOSTICS_SCHEMA_ID,
            ),
        )

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_id": self.schema_id,
            "quoted_c": self.quoted_c.tolist(),
            "structural_eigenvalues": self.structural_eigenvalues.tolist(),
            "delta_c": self.delta_c.tolist(),
            "max_abs_delta_c": self.max_abs_delta_c,
            "quoted_f": self.quoted_f.tolist(),
            "geometry_derived_f_from_quoted_c": self.geometry_derived_f_from_quoted_c.tolist(),
            "delta_f": self.delta_f.tolist(),
            "max_abs_delta_f": self.max_abs_delta_f,
        }


@dataclass(frozen=True)
class Paper07101869Point:
    """Custom seeded quark spurion point with structural reference metadata only."""

    Y_u: np.ndarray
    Y_d: np.ndarray
    q_sector: Paper07101869BenchmarkSector
    u_sector: Paper07101869BenchmarkSector
    d_sector: Paper07101869BenchmarkSector
    r: float
    Lambda_IR: float = DEFAULT_LAMBDA_IR
    k: float = MPL
    v: float = V_EWSB
    label: str = "custom_seeded_structural_reference"
    profile_input_policy_id: str = PAPER_0710_1869_REFERENCE_ONLY_PROFILE_INPUT_POLICY_ID
    point_kind_id: str = PAPER_0710_1869_CUSTOM_POINT_KIND_ID
    claim_level_id: str = PAPER_0710_1869_NONCANONICAL_CLAIM_LEVEL_ID
    bulk_mass_map_policy_id: str = PAPER_0710_1869_NO_BULK_MASS_MAP_POLICY_ID
    structural_reference_status_id: str = PAPER_0710_1869_STRUCTURAL_REFERENCE_STATUS_ID
    construction_id: str = PAPER_0710_1869_EXPLICIT_SEED_CONSTRUCTION_ID
    metadata: Mapping[str, Any] = field(default_factory=dict)
    notes: str | None = None
    schema_id: str = PAPER_0710_1869_POINT_SCHEMA_ID
    mode_id: str = PAPER_0710_1869_MODE_ID
    paper_id: str = PAPER_0710_1869_PAPER_ID

    def __post_init__(self) -> None:
        object.__setattr__(self, "Y_u", _as_complex_matrix3("Y_u", self.Y_u))
        object.__setattr__(self, "Y_d", _as_complex_matrix3("Y_d", self.Y_d))
        if not isinstance(self.q_sector, Paper07101869BenchmarkSector):
            raise ValueError("q_sector must be a Paper07101869BenchmarkSector")
        if not isinstance(self.u_sector, Paper07101869BenchmarkSector):
            raise ValueError("u_sector must be a Paper07101869BenchmarkSector")
        if not isinstance(self.d_sector, Paper07101869BenchmarkSector):
            raise ValueError("d_sector must be a Paper07101869BenchmarkSector")
        object.__setattr__(self, "r", _require_nonnegative_finite("r", self.r))
        object.__setattr__(self, "Lambda_IR", require_positive_finite("Lambda_IR", self.Lambda_IR))
        object.__setattr__(self, "k", require_positive_finite("k", self.k))
        object.__setattr__(self, "v", require_positive_finite("v", self.v))
        object.__setattr__(self, "label", require_nonempty_identifier("label", self.label))
        object.__setattr__(
            self,
            "profile_input_policy_id",
            require_member(
                "profile_input_policy_id",
                self.profile_input_policy_id,
                _ALLOWED_PROFILE_INPUT_POLICY_IDS,
            ),
        )
        object.__setattr__(
            self,
            "point_kind_id",
            require_member("point_kind_id", self.point_kind_id, _ALLOWED_POINT_KIND_IDS),
        )
        object.__setattr__(
            self,
            "claim_level_id",
            require_member("claim_level_id", self.claim_level_id, _ALLOWED_CLAIM_LEVEL_IDS),
        )
        object.__setattr__(
            self,
            "bulk_mass_map_policy_id",
            require_member(
                "bulk_mass_map_policy_id",
                self.bulk_mass_map_policy_id,
                _ALLOWED_BULK_MASS_MAP_POLICY_IDS,
            ),
        )
        object.__setattr__(
            self,
            "structural_reference_status_id",
            require_member(
                "structural_reference_status_id",
                self.structural_reference_status_id,
                _ALLOWED_STRUCTURAL_REFERENCE_STATUS_IDS,
            ),
        )
        object.__setattr__(
            self,
            "construction_id",
            require_member("construction_id", self.construction_id, _ALLOWED_CONSTRUCTION_IDS),
        )
        if self.metadata is None:
            object.__setattr__(self, "metadata", {})
        elif not isinstance(self.metadata, Mapping):
            raise ValueError("metadata must be a mapping")
        else:
            object.__setattr__(self, "metadata", dict(self.metadata))
        if self.notes is not None:
            object.__setattr__(self, "notes", require_nonempty_identifier("notes", self.notes))
        object.__setattr__(
            self,
            "schema_id",
            require_known_schema_id(
                "schema_id",
                self.schema_id,
                expected=PAPER_0710_1869_POINT_SCHEMA_ID,
            ),
        )
        object.__setattr__(self, "mode_id", require_nonempty_identifier("mode_id", self.mode_id))
        object.__setattr__(self, "paper_id", require_nonempty_identifier("paper_id", self.paper_id))
        require_member("mode_id", self.mode_id, (PAPER_0710_1869_MODE_ID,))
        require_member("paper_id", self.paper_id, (PAPER_0710_1869_PAPER_ID,))


@dataclass(frozen=True)
class Paper07101869BulkState:
    """Structural MFV eigensystem state plus benchmark-reference diagnostics."""

    point: Paper07101869Point
    epsilon: float
    C_Q: np.ndarray
    C_u: np.ndarray
    C_d: np.ndarray
    # Legacy aliases retained for downstream compatibility. These are quoted
    # Table I reference eigenvalues only, not point-derived physical bulk masses.
    c_Q: np.ndarray
    c_u: np.ndarray
    c_d: np.ndarray
    # Legacy aliases retained for downstream compatibility. These are quoted
    # Table I reference profiles only, not active inputs used by the core state.
    F_Q: np.ndarray
    F_u: np.ndarray
    F_d: np.ndarray
    # Geometry-derived profiles from the quoted Table I reference c values.
    derived_F_Q: np.ndarray
    derived_F_u: np.ndarray
    derived_F_d: np.ndarray
    profile_consistency_Q: Paper07101869ProfileConsistency
    profile_consistency_u: Paper07101869ProfileConsistency
    profile_consistency_d: Paper07101869ProfileConsistency
    reference_diagnostics_Q: Paper07101869BenchmarkReferenceDiagnostics
    reference_diagnostics_u: Paper07101869BenchmarkReferenceDiagnostics
    reference_diagnostics_d: Paper07101869BenchmarkReferenceDiagnostics
    reference_eq3_summary: Paper07101869Eq3ResidualSummary | None
    rotation_Q: np.ndarray
    rotation_u: np.ndarray
    rotation_d: np.ndarray
    eig_Q: np.ndarray
    eig_u: np.ndarray
    eig_d: np.ndarray
    Y_u_bulk_basis: np.ndarray
    Y_d_bulk_basis: np.ndarray
    core_contract_id: str = PAPER_0710_1869_STRUCTURAL_CORE_CONTRACT_ID
    point_kind_id: str = PAPER_0710_1869_CUSTOM_POINT_KIND_ID
    claim_level_id: str = PAPER_0710_1869_NONCANONICAL_CLAIM_LEVEL_ID
    bulk_mass_map_policy_id: str = PAPER_0710_1869_NO_BULK_MASS_MAP_POLICY_ID
    structural_reference_status_id: str = PAPER_0710_1869_STRUCTURAL_REFERENCE_STATUS_ID
    schema_id: str = PAPER_0710_1869_BULK_STATE_SCHEMA_ID

    def __post_init__(self) -> None:
        if not isinstance(self.point, Paper07101869Point):
            raise ValueError("point must be a Paper07101869Point")
        object.__setattr__(self, "epsilon", require_positive_finite("epsilon", self.epsilon))
        for name in (
            "C_Q",
            "C_u",
            "C_d",
            "rotation_Q",
            "rotation_u",
            "rotation_d",
            "Y_u_bulk_basis",
            "Y_d_bulk_basis",
        ):
            object.__setattr__(self, name, _as_complex_matrix3(name, getattr(self, name)))
        for name in (
            "c_Q",
            "c_u",
            "c_d",
            "F_Q",
            "F_u",
            "F_d",
            "derived_F_Q",
            "derived_F_u",
            "derived_F_d",
            "eig_Q",
            "eig_u",
            "eig_d",
        ):
            object.__setattr__(self, name, _as_real_triplet(name, getattr(self, name)))
        for name in (
            "profile_consistency_Q",
            "profile_consistency_u",
            "profile_consistency_d",
        ):
            if not isinstance(getattr(self, name), Paper07101869ProfileConsistency):
                raise ValueError(f"{name} must be a Paper07101869ProfileConsistency")
        for name in (
            "reference_diagnostics_Q",
            "reference_diagnostics_u",
            "reference_diagnostics_d",
        ):
            if not isinstance(getattr(self, name), Paper07101869BenchmarkReferenceDiagnostics):
                raise ValueError(f"{name} must be a Paper07101869BenchmarkReferenceDiagnostics")
        if self.reference_eq3_summary is not None and not isinstance(
            self.reference_eq3_summary, Paper07101869Eq3ResidualSummary
        ):
            raise ValueError("reference_eq3_summary must be a Paper07101869Eq3ResidualSummary")
        object.__setattr__(
            self,
            "core_contract_id",
            require_nonempty_identifier("core_contract_id", self.core_contract_id),
        )
        object.__setattr__(
            self,
            "point_kind_id",
            require_member("point_kind_id", self.point_kind_id, _ALLOWED_POINT_KIND_IDS),
        )
        object.__setattr__(
            self,
            "claim_level_id",
            require_member("claim_level_id", self.claim_level_id, _ALLOWED_CLAIM_LEVEL_IDS),
        )
        object.__setattr__(
            self,
            "bulk_mass_map_policy_id",
            require_member(
                "bulk_mass_map_policy_id",
                self.bulk_mass_map_policy_id,
                _ALLOWED_BULK_MASS_MAP_POLICY_IDS,
            ),
        )
        object.__setattr__(
            self,
            "structural_reference_status_id",
            require_member(
                "structural_reference_status_id",
                self.structural_reference_status_id,
                _ALLOWED_STRUCTURAL_REFERENCE_STATUS_IDS,
            ),
        )
        if self.point.point_kind_id != self.point_kind_id:
            raise ValueError("bulk-state point_kind_id must match point.point_kind_id")
        if self.point.claim_level_id != self.claim_level_id:
            raise ValueError("bulk-state claim_level_id must match point.claim_level_id")
        if self.point.bulk_mass_map_policy_id != self.bulk_mass_map_policy_id:
            raise ValueError(
                "bulk-state bulk_mass_map_policy_id must match point.bulk_mass_map_policy_id"
            )
        if self.point.structural_reference_status_id != self.structural_reference_status_id:
            raise ValueError(
                "bulk-state structural_reference_status_id must match "
                "point.structural_reference_status_id"
            )
        if self.core_contract_id != PAPER_0710_1869_STRUCTURAL_CORE_CONTRACT_ID:
            raise ValueError(
                "core_contract_id must identify the seeded structural MFV core contract"
            )
        for name, diagnostic, quoted_c, quoted_f, derived_f in (
            (
                "Q",
                self.reference_diagnostics_Q,
                self.c_Q,
                self.F_Q,
                self.derived_F_Q,
            ),
            (
                "u",
                self.reference_diagnostics_u,
                self.c_u,
                self.F_u,
                self.derived_F_u,
            ),
            (
                "d",
                self.reference_diagnostics_d,
                self.c_d,
                self.F_d,
                self.derived_F_d,
            ),
        ):
            if not np.allclose(diagnostic.quoted_c, quoted_c, atol=1.0e-12):
                raise ValueError(f"{name}-sector c alias must match quoted reference values")
            if not np.allclose(diagnostic.quoted_f, quoted_f, atol=1.0e-12):
                raise ValueError(f"{name}-sector F alias must match quoted reference values")
            if not np.allclose(
                diagnostic.geometry_derived_f_from_quoted_c,
                derived_f,
                atol=1.0e-12,
            ):
                raise ValueError(
                    f"{name}-sector derived_F alias must match geometry-derived reference values"
                )
        object.__setattr__(
            self,
            "schema_id",
            require_known_schema_id(
                "schema_id",
                self.schema_id,
                expected=PAPER_0710_1869_BULK_STATE_SCHEMA_ID,
            ),
        )


@dataclass(frozen=True)
class Paper07101869PhysicalPoint:
    """Custom seeded point for the hard-gated QS1 physical seed-to-profile bridge."""

    Y_u: np.ndarray
    Y_d: np.ndarray
    physical_contract: Paper07101869PhysicalSeedToProfileContract = field(
        default_factory=default_paper_0710_1869_physical_seed_to_profile_contract
    )
    r: float = 0.0
    Lambda_IR: float = DEFAULT_LAMBDA_IR
    k: float = MPL
    v: float = V_EWSB
    label: str = "custom_seeded_physical_qs1"
    profile_input_policy_id: str = PAPER_0710_1869_DERIVED_PROFILE_INPUT_POLICY_ID
    point_kind_id: str = PAPER_0710_1869_PHYSICAL_POINT_KIND_ID
    claim_level_id: str = PAPER_0710_1869_PHYSICAL_CLAIM_LEVEL_ID
    bulk_mass_map_policy_id: str = PAPER_0710_1869_SEED_TO_PROFILE_MAPPING_POLICY_ID
    physical_profile_status_id: str = PAPER_0710_1869_POINT_DERIVED_PHYSICAL_STATUS_ID
    construction_id: str = PAPER_0710_1869_PHYSICAL_SEED_CONSTRUCTION_ID
    metadata: Mapping[str, Any] = field(default_factory=dict)
    notes: str | None = None
    schema_id: str = PAPER_0710_1869_PHYSICAL_POINT_SCHEMA_ID
    mode_id: str = PAPER_0710_1869_MODE_ID
    paper_id: str = PAPER_0710_1869_PAPER_ID

    def __post_init__(self) -> None:
        object.__setattr__(self, "Y_u", _as_complex_matrix3("Y_u", self.Y_u))
        object.__setattr__(self, "Y_d", _as_complex_matrix3("Y_d", self.Y_d))
        object.__setattr__(
            self,
            "physical_contract",
            _require_exact_frozen_physical_contract(self.physical_contract),
        )
        object.__setattr__(self, "r", _require_nonnegative_finite("r", self.r))
        object.__setattr__(self, "Lambda_IR", require_positive_finite("Lambda_IR", self.Lambda_IR))
        object.__setattr__(self, "k", require_positive_finite("k", self.k))
        object.__setattr__(self, "v", require_positive_finite("v", self.v))
        object.__setattr__(self, "label", require_nonempty_identifier("label", self.label))
        object.__setattr__(
            self,
            "profile_input_policy_id",
            require_member(
                "profile_input_policy_id",
                self.profile_input_policy_id,
                (PAPER_0710_1869_DERIVED_PROFILE_INPUT_POLICY_ID,),
            ),
        )
        object.__setattr__(
            self,
            "point_kind_id",
            require_member(
                "point_kind_id",
                self.point_kind_id,
                (PAPER_0710_1869_PHYSICAL_POINT_KIND_ID,),
            ),
        )
        object.__setattr__(
            self,
            "claim_level_id",
            require_member(
                "claim_level_id",
                self.claim_level_id,
                (PAPER_0710_1869_PHYSICAL_CLAIM_LEVEL_ID,),
            ),
        )
        object.__setattr__(
            self,
            "bulk_mass_map_policy_id",
            require_member(
                "bulk_mass_map_policy_id",
                self.bulk_mass_map_policy_id,
                (PAPER_0710_1869_SEED_TO_PROFILE_MAPPING_POLICY_ID,),
            ),
        )
        object.__setattr__(
            self,
            "physical_profile_status_id",
            require_member(
                "physical_profile_status_id",
                self.physical_profile_status_id,
                (PAPER_0710_1869_POINT_DERIVED_PHYSICAL_STATUS_ID,),
            ),
        )
        object.__setattr__(
            self,
            "construction_id",
            require_member(
                "construction_id",
                self.construction_id,
                (
                    PAPER_0710_1869_PHYSICAL_SEED_CONSTRUCTION_ID,
                    PAPER_0710_1869_PHYSICAL_BENCHMARK_SEED_CONSTRUCTION_ID,
                ),
            ),
        )
        if self.metadata is None:
            object.__setattr__(self, "metadata", {})
        elif not isinstance(self.metadata, Mapping):
            raise ValueError("metadata must be a mapping")
        else:
            object.__setattr__(self, "metadata", dict(self.metadata))
        if self.notes is not None:
            object.__setattr__(self, "notes", require_nonempty_identifier("notes", self.notes))
        object.__setattr__(
            self,
            "schema_id",
            require_known_schema_id(
                "schema_id",
                self.schema_id,
                expected=PAPER_0710_1869_PHYSICAL_POINT_SCHEMA_ID,
            ),
        )
        object.__setattr__(self, "mode_id", require_nonempty_identifier("mode_id", self.mode_id))
        object.__setattr__(self, "paper_id", require_nonempty_identifier("paper_id", self.paper_id))
        require_member("mode_id", self.mode_id, (PAPER_0710_1869_MODE_ID,))
        require_member("paper_id", self.paper_id, (PAPER_0710_1869_PAPER_ID,))
        if self.bulk_mass_map_policy_id != self.physical_contract.mapping_policy.policy_id:
            raise ValueError(
                "bulk_mass_map_policy_id must match the exact frozen seed-to-profile policy id"
            )


@dataclass(frozen=True)
class Paper07101869PhysicalBulkState:
    """Point-derived QS1 physical state with affine c and geometry-derived F values."""

    point: Paper07101869PhysicalPoint
    epsilon: float
    C_Q: np.ndarray
    C_u: np.ndarray
    C_d: np.ndarray
    c_Q: np.ndarray
    c_u: np.ndarray
    c_d: np.ndarray
    F_Q: np.ndarray
    F_u: np.ndarray
    F_d: np.ndarray
    rotation_Q: np.ndarray
    rotation_u: np.ndarray
    rotation_d: np.ndarray
    eig_Q: np.ndarray
    eig_u: np.ndarray
    eig_d: np.ndarray
    Y_u_bulk_basis: np.ndarray
    Y_d_bulk_basis: np.ndarray
    physical_contract: Paper07101869PhysicalSeedToProfileContract = field(
        default_factory=default_paper_0710_1869_physical_seed_to_profile_contract
    )
    core_contract_id: str = PAPER_0710_1869_PHYSICAL_CORE_CONTRACT_ID
    point_kind_id: str = PAPER_0710_1869_PHYSICAL_POINT_KIND_ID
    claim_level_id: str = PAPER_0710_1869_PHYSICAL_CLAIM_LEVEL_ID
    bulk_mass_map_policy_id: str = PAPER_0710_1869_SEED_TO_PROFILE_MAPPING_POLICY_ID
    profile_input_policy_id: str = PAPER_0710_1869_DERIVED_PROFILE_INPUT_POLICY_ID
    physical_profile_status_id: str = PAPER_0710_1869_POINT_DERIVED_PHYSICAL_STATUS_ID
    schema_id: str = PAPER_0710_1869_PHYSICAL_BULK_STATE_SCHEMA_ID

    def __post_init__(self) -> None:
        if not isinstance(self.point, Paper07101869PhysicalPoint):
            raise ValueError("point must be a Paper07101869PhysicalPoint")
        object.__setattr__(self, "epsilon", require_positive_finite("epsilon", self.epsilon))
        for name in (
            "C_Q",
            "C_u",
            "C_d",
            "rotation_Q",
            "rotation_u",
            "rotation_d",
            "Y_u_bulk_basis",
            "Y_d_bulk_basis",
        ):
            object.__setattr__(self, name, _as_complex_matrix3(name, getattr(self, name)))
        for name in ("c_Q", "c_u", "c_d", "eig_Q", "eig_u", "eig_d"):
            object.__setattr__(self, name, _as_real_triplet(name, getattr(self, name)))
        for name in ("F_Q", "F_u", "F_d"):
            object.__setattr__(self, name, _as_positive_real_triplet(name, getattr(self, name)))
        object.__setattr__(
            self,
            "physical_contract",
            _require_exact_frozen_physical_contract(self.physical_contract),
        )
        object.__setattr__(
            self,
            "core_contract_id",
            require_member(
                "core_contract_id",
                self.core_contract_id,
                (PAPER_0710_1869_PHYSICAL_CORE_CONTRACT_ID,),
            ),
        )
        object.__setattr__(
            self,
            "point_kind_id",
            require_member(
                "point_kind_id",
                self.point_kind_id,
                (PAPER_0710_1869_PHYSICAL_POINT_KIND_ID,),
            ),
        )
        object.__setattr__(
            self,
            "claim_level_id",
            require_member(
                "claim_level_id",
                self.claim_level_id,
                (PAPER_0710_1869_PHYSICAL_CLAIM_LEVEL_ID,),
            ),
        )
        object.__setattr__(
            self,
            "bulk_mass_map_policy_id",
            require_member(
                "bulk_mass_map_policy_id",
                self.bulk_mass_map_policy_id,
                (PAPER_0710_1869_SEED_TO_PROFILE_MAPPING_POLICY_ID,),
            ),
        )
        object.__setattr__(
            self,
            "profile_input_policy_id",
            require_member(
                "profile_input_policy_id",
                self.profile_input_policy_id,
                (PAPER_0710_1869_DERIVED_PROFILE_INPUT_POLICY_ID,),
            ),
        )
        object.__setattr__(
            self,
            "physical_profile_status_id",
            require_member(
                "physical_profile_status_id",
                self.physical_profile_status_id,
                (PAPER_0710_1869_POINT_DERIVED_PHYSICAL_STATUS_ID,),
            ),
        )
        if self.point.point_kind_id != self.point_kind_id:
            raise ValueError("physical bulk-state point_kind_id must match point.point_kind_id")
        if self.point.claim_level_id != self.claim_level_id:
            raise ValueError("physical bulk-state claim_level_id must match point.claim_level_id")
        if self.point.bulk_mass_map_policy_id != self.bulk_mass_map_policy_id:
            raise ValueError(
                "physical bulk-state bulk_mass_map_policy_id must match point.bulk_mass_map_policy_id"
            )
        if self.point.profile_input_policy_id != self.profile_input_policy_id:
            raise ValueError(
                "physical bulk-state profile_input_policy_id must match point.profile_input_policy_id"
            )
        if self.point.physical_profile_status_id != self.physical_profile_status_id:
            raise ValueError(
                "physical bulk-state physical_profile_status_id must match "
                "point.physical_profile_status_id"
            )
        if self.point.physical_contract != self.physical_contract:
            raise ValueError("physical bulk-state contract must match point.physical_contract")
        if self.bulk_mass_map_policy_id != self.physical_contract.mapping_policy.policy_id:
            raise ValueError(
                "physical bulk-state bulk_mass_map_policy_id must match "
                "physical_contract.mapping_policy.policy_id"
            )
        for name in ("rotation_Q", "rotation_u", "rotation_d"):
            object.__setattr__(
                self,
                name,
                _require_unitary_matrix3(
                    name,
                    getattr(self, name),
                    tolerance=PAPER_0710_1869_PHYSICAL_UNITARITY_TOLERANCE,
                ),
            )
        expected_payload = _build_point_derived_physical_bulk_payload(self.point)
        if self.physical_contract != expected_payload["physical_contract"]:
            raise ValueError("physical bulk-state contract must match the point-derived payload")
        _require_allclose("epsilon", self.epsilon, expected_payload["epsilon"])
        for name in (
            "C_Q",
            "C_u",
            "C_d",
            "c_Q",
            "c_u",
            "c_d",
            "F_Q",
            "F_u",
            "F_d",
            "rotation_Q",
            "rotation_u",
            "rotation_d",
            "eig_Q",
            "eig_u",
            "eig_d",
            "Y_u_bulk_basis",
            "Y_d_bulk_basis",
        ):
            _require_allclose(name, getattr(self, name), expected_payload[name])
        object.__setattr__(
            self,
            "schema_id",
            require_known_schema_id(
                "schema_id",
                self.schema_id,
                expected=PAPER_0710_1869_PHYSICAL_BULK_STATE_SCHEMA_ID,
            ),
        )

@dataclass(frozen=True)
class Paper07101869Eq3Example:
    """Quoted Eq. (3) example from the paper."""

    relation_id: str = PAPER_0710_1869_EQ3_RELATION_ID
    reference: str = PAPER_0710_1869_EQ3_REFERENCE
    a: float = 0.8
    r: float = 0.3
    v5km_parameters: Paper07101869V5KMParameters = field(
        default_factory=lambda: Paper07101869V5KMParameters.from_degrees(
            theta12_deg=115.0,
            theta23_deg=65.0,
            theta13_deg=70.0,
            delta=0.6,
        )
    )

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "relation_id", require_nonempty_identifier("relation_id", self.relation_id)
        )
        object.__setattr__(
            self, "reference", require_nonempty_identifier("reference", self.reference)
        )
        object.__setattr__(self, "a", require_positive_finite("a", self.a))
        object.__setattr__(self, "r", _require_finite("r", self.r))
        if not isinstance(self.v5km_parameters, Paper07101869V5KMParameters):
            raise ValueError("v5km_parameters must be a Paper07101869V5KMParameters")


@dataclass(frozen=True)
class Paper07101869Eq3ResidualSummary:
    """Residual summary for the paper-side Eq. (3) consistency check."""

    lhs_diag: np.ndarray
    rhs_diag: np.ndarray
    residual: np.ndarray
    a: float
    r: float
    max_rhs_imag_abs: float
    max_abs_residual: float
    l2_residual: float
    exact_agreement: bool
    tolerance: float
    relation_id: str = PAPER_0710_1869_EQ3_RELATION_ID

    def __post_init__(self) -> None:
        object.__setattr__(self, "lhs_diag", _as_real_triplet("lhs_diag", self.lhs_diag))
        object.__setattr__(self, "rhs_diag", _as_real_triplet("rhs_diag", self.rhs_diag))
        object.__setattr__(self, "residual", _as_real_triplet("residual", self.residual))
        object.__setattr__(self, "a", require_positive_finite("a", self.a))
        object.__setattr__(self, "r", _require_finite("r", self.r))
        object.__setattr__(
            self,
            "max_rhs_imag_abs",
            _require_nonnegative_finite("max_rhs_imag_abs", self.max_rhs_imag_abs),
        )
        object.__setattr__(
            self,
            "max_abs_residual",
            _require_nonnegative_finite("max_abs_residual", self.max_abs_residual),
        )
        object.__setattr__(
            self,
            "l2_residual",
            _require_nonnegative_finite("l2_residual", self.l2_residual),
        )
        object.__setattr__(
            self, "tolerance", require_positive_finite("tolerance", self.tolerance)
        )
        object.__setattr__(
            self, "relation_id", require_nonempty_identifier("relation_id", self.relation_id)
        )

    def as_dict(self) -> dict[str, float | bool | list[float] | str]:
        """Return a stable summary representation for tests and scripts."""
        return {
            "relation_id": self.relation_id,
            "lhs_diag": self.lhs_diag.tolist(),
            "rhs_diag": self.rhs_diag.tolist(),
            "residual": self.residual.tolist(),
            "a": self.a,
            "r": self.r,
            "max_rhs_imag_abs": self.max_rhs_imag_abs,
            "max_abs_residual": self.max_abs_residual,
            "l2_residual": self.l2_residual,
            "exact_agreement": self.exact_agreement,
            "tolerance": self.tolerance,
        }


def build_diagonal_c_matrix(
    eigenvalues: tuple[float, float, float] | list[float] | np.ndarray,
) -> np.ndarray:
    """Build a diagonal 3x3 ``C`` matrix from benchmark eigenvalues."""
    return np.diag(_as_real_triplet("eigenvalues", eigenvalues)).astype(np.complex128)


def build_table_i_diagonal_c_matrices(
    benchmark: Paper07101869TableIBenchmark | None = None,
) -> Paper07101869DiagonalCMatrices:
    """Return the paper-side diagonal ``C_Q``, ``C_u``, and ``C_d`` matrices."""
    resolved = default_paper_0710_1869_table_i_benchmark() if benchmark is None else benchmark
    return Paper07101869DiagonalCMatrices(
        C_Q=resolved.q_sector.build_diagonal_c_matrix(),
        C_u=resolved.u_sector.build_diagonal_c_matrix(),
        C_d=resolved.d_sector.build_diagonal_c_matrix(),
    )


def build_table_i_benchmark_from_inputs(
    table_inputs: Paper07101869TableIInputs,
) -> Paper07101869TableIBenchmark:
    """Convert the sourced Table I input contract into model-side benchmark sectors."""
    if not isinstance(table_inputs, Paper07101869TableIInputs):
        raise ValueError("table_inputs must be a Paper07101869TableIInputs")

    sector_map = {
        sector.sector_id: Paper07101869BenchmarkSector(
            label=sector.sector_id,
            c_eigenvalues=tuple(entry.c_value for entry in sector.entries),
            f_eigenvalues=tuple(entry.f_value for entry in sector.entries),
        )
        for sector in table_inputs.sectors
    }
    return Paper07101869TableIBenchmark(
        q_sector=sector_map["Q"],
        u_sector=sector_map["u"],
        d_sector=sector_map["d"],
    )


def paper_ckm_like_unitary(params: Paper07101869RotationParameters) -> np.ndarray:
    """Return a CKM-like unitary matrix from paper-owned angles and phase."""
    if not isinstance(params, Paper07101869RotationParameters):
        raise ValueError("params must be a Paper07101869RotationParameters")
    s12, s13, s23 = np.sin([params.theta12, params.theta13, params.theta23])
    c12, c13, c23 = np.cos([params.theta12, params.theta13, params.theta23])
    phase = np.exp(-1j * params.delta)
    phase_conj = np.conjugate(phase)
    return np.array(
        [
            [c12 * c13, s12 * c13, s13 * phase],
            [
                -s12 * c23 - c12 * s23 * s13 * phase_conj,
                c12 * c23 - s12 * s23 * s13 * phase_conj,
                s23 * c13,
            ],
            [
                s12 * s23 - c12 * c23 * s13 * phase_conj,
                -c12 * s23 - s12 * c23 * s13 * phase_conj,
                c23 * c13,
            ],
        ],
        dtype=np.complex128,
    )


def paper_v5km_matrix(params: Paper07101869V5KMParameters) -> np.ndarray:
    """Return the paper-side CKM-like ``V5KM`` matrix."""
    if not isinstance(params, Paper07101869V5KMParameters):
        raise ValueError("params must be a Paper07101869V5KMParameters")
    return paper_ckm_like_unitary(
        Paper07101869RotationParameters(
            theta12=params.theta12,
            theta13=params.theta13,
            theta23=params.theta23,
            delta=params.delta,
        )
    )


def build_paper_0710_1869_spurion_matrix(
    singular_values: np.ndarray | list[float] | tuple[float, ...],
    *,
    overall_scale: float = 1.0,
    left_rotation: Paper07101869RotationParameters | None = None,
    right_rotation: Paper07101869RotationParameters | None = None,
) -> np.ndarray:
    """Construct a paper-owned Yukawa spurion from singular values and rotations."""
    sigma = _as_nonnegative_real_triplet("singular_values", singular_values)
    resolved_scale = require_positive_finite("overall_scale", overall_scale)
    U_left = paper_ckm_like_unitary(left_rotation or Paper07101869RotationParameters())
    U_right = paper_ckm_like_unitary(right_rotation or Paper07101869RotationParameters())
    return resolved_scale * U_left @ np.diag(sigma) @ U_right.conjugate().T


def build_paper_0710_1869_point_from_singular_values(
    *,
    up_singular_values: np.ndarray | list[float] | tuple[float, ...],
    down_singular_values: np.ndarray | list[float] | tuple[float, ...],
    q_sector: Paper07101869BenchmarkSector,
    u_sector: Paper07101869BenchmarkSector,
    d_sector: Paper07101869BenchmarkSector,
    profile_input_policy_id: str = PAPER_0710_1869_REFERENCE_ONLY_PROFILE_INPUT_POLICY_ID,
    overall_scale: float = 1.0,
    r: float,
    up_left: Paper07101869RotationParameters | None = None,
    up_right: Paper07101869RotationParameters | None = None,
    down_left: Paper07101869RotationParameters | None = None,
    down_right: Paper07101869RotationParameters | None = None,
    Lambda_IR: float = DEFAULT_LAMBDA_IR,
    k: float = MPL,
    v: float = V_EWSB,
    label: str = "custom_seeded_structural_reference",
    construction_id: str = PAPER_0710_1869_EXPLICIT_SEED_CONSTRUCTION_ID,
    metadata: Mapping[str, Any] | None = None,
    notes: str | None = None,
) -> Paper07101869Point:
    """Build a custom seeded point with diagnostics-only Table I references."""
    Y_u = build_paper_0710_1869_spurion_matrix(
        up_singular_values,
        overall_scale=overall_scale,
        left_rotation=up_left,
        right_rotation=up_right,
    )
    Y_d = build_paper_0710_1869_spurion_matrix(
        down_singular_values,
        overall_scale=overall_scale,
        left_rotation=down_left,
        right_rotation=down_right,
    )
    return Paper07101869Point(
        Y_u=Y_u,
        Y_d=Y_d,
        q_sector=q_sector,
        u_sector=u_sector,
        d_sector=d_sector,
        r=r,
        Lambda_IR=Lambda_IR,
        k=k,
        v=v,
        label=label,
        profile_input_policy_id=profile_input_policy_id,
        construction_id=construction_id,
        metadata=dict(metadata or {}),
        notes=notes,
    )


def build_paper_0710_1869_physical_point_from_singular_values(
    *,
    up_singular_values: np.ndarray | list[float] | tuple[float, ...],
    down_singular_values: np.ndarray | list[float] | tuple[float, ...],
    physical_contract: Paper07101869PhysicalSeedToProfileContract | None = None,
    overall_scale: float = 1.0,
    r: float,
    up_left: Paper07101869RotationParameters | None = None,
    up_right: Paper07101869RotationParameters | None = None,
    down_left: Paper07101869RotationParameters | None = None,
    down_right: Paper07101869RotationParameters | None = None,
    Lambda_IR: float = DEFAULT_LAMBDA_IR,
    k: float = MPL,
    v: float = V_EWSB,
    label: str = "custom_seeded_physical_qs1",
    construction_id: str = PAPER_0710_1869_PHYSICAL_SEED_CONSTRUCTION_ID,
    metadata: Mapping[str, Any] | None = None,
    notes: str | None = None,
) -> Paper07101869PhysicalPoint:
    """Build a hard-gated QS1 physical point from seeded Yukawa spurions only."""
    resolved_contract = _require_exact_frozen_physical_contract(
        default_paper_0710_1869_physical_seed_to_profile_contract()
        if physical_contract is None
        else physical_contract
    )
    Y_u = build_paper_0710_1869_spurion_matrix(
        up_singular_values,
        overall_scale=overall_scale,
        left_rotation=up_left,
        right_rotation=up_right,
    )
    Y_d = build_paper_0710_1869_spurion_matrix(
        down_singular_values,
        overall_scale=overall_scale,
        left_rotation=down_left,
        right_rotation=down_right,
    )
    return Paper07101869PhysicalPoint(
        Y_u=Y_u,
        Y_d=Y_d,
        physical_contract=resolved_contract,
        r=r,
        Lambda_IR=Lambda_IR,
        k=k,
        v=v,
        label=label,
        construction_id=construction_id,
        metadata=dict(metadata or {}),
        notes=notes,
    )


def _build_profile_consistency(
    *,
    quoted_f: np.ndarray,
    derived_f: np.ndarray,
    policy_id: str,
) -> Paper07101869ProfileConsistency:
    delta = derived_f - quoted_f
    return Paper07101869ProfileConsistency(
        quoted_f=quoted_f,
        derived_f=derived_f,
        delta_f=delta,
        max_abs_delta=float(np.max(np.abs(delta))),
        policy_id=policy_id,
    )


def _build_reference_diagnostics(
    *,
    quoted_c: np.ndarray,
    structural_eigenvalues: np.ndarray,
    quoted_f: np.ndarray,
    geometry_derived_f_from_quoted_c: np.ndarray,
) -> Paper07101869BenchmarkReferenceDiagnostics:
    delta_c = structural_eigenvalues - quoted_c
    delta_f = geometry_derived_f_from_quoted_c - quoted_f
    return Paper07101869BenchmarkReferenceDiagnostics(
        quoted_c=quoted_c,
        structural_eigenvalues=structural_eigenvalues,
        delta_c=delta_c,
        max_abs_delta_c=float(np.max(np.abs(delta_c))),
        quoted_f=quoted_f,
        geometry_derived_f_from_quoted_c=geometry_derived_f_from_quoted_c,
        delta_f=delta_f,
        max_abs_delta_f=float(np.max(np.abs(delta_f))),
    )


def _sector_policy_by_id(
    contract: Paper07101869PhysicalSeedToProfileContract,
    sector_id: str,
) -> Paper07101869AffineBulkMassSectorPolicy:
    resolved_contract = _require_exact_frozen_physical_contract(contract)
    for policy in resolved_contract.universal_term_policy.sector_policies:
        if policy.sector_id == sector_id:
            return policy
    raise ValueError(f"physical contract is missing sector policy {sector_id!r}")


def _build_point_derived_physical_sector(
    *,
    matrix: np.ndarray,
    epsilon: float,
    policy: Paper07101869AffineBulkMassSectorPolicy,
    c_name: str,
    f_name: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    eig, rotation = _ordered_hermitian_spectrum(matrix)
    c_values = _as_real_triplet(
        c_name,
        policy.leading_term_coefficient * eig + policy.universal_offset,
    )
    # RESIDUAL(C-6): exact Table-I affine coefficients pending paper 0710.1869.
    F_values = _as_positive_real_triplet(
        f_name,
        np.asarray(f_IR(c_values, epsilon), dtype=float),
    )
    return eig, rotation, c_values, F_values


def _build_point_derived_physical_bulk_payload(
    point: Paper07101869PhysicalPoint,
) -> dict[str, Any]:
    if not isinstance(point, Paper07101869PhysicalPoint):
        raise ValueError("point must be a Paper07101869PhysicalPoint")

    params = get_warp_params(k=point.k, Lambda_IR=point.Lambda_IR)
    epsilon = float(params["epsilon"])
    physical_contract = _require_exact_frozen_physical_contract(point.physical_contract)

    C_u = point.Y_u.conjugate().T @ point.Y_u
    C_d = point.Y_d.conjugate().T @ point.Y_d
    C_Q = point.r * (point.Y_u @ point.Y_u.conjugate().T) + (
        point.Y_d @ point.Y_d.conjugate().T
    )

    eig_Q, rotation_Q, c_Q, F_Q = _build_point_derived_physical_sector(
        matrix=C_Q,
        epsilon=epsilon,
        policy=_sector_policy_by_id(physical_contract, "Q"),
        c_name="c_Q",
        f_name="F_Q",
    )
    eig_u, rotation_u, c_u, F_u = _build_point_derived_physical_sector(
        matrix=C_u,
        epsilon=epsilon,
        policy=_sector_policy_by_id(physical_contract, "u"),
        c_name="c_u",
        f_name="F_u",
    )
    eig_d, rotation_d, c_d, F_d = _build_point_derived_physical_sector(
        matrix=C_d,
        epsilon=epsilon,
        policy=_sector_policy_by_id(physical_contract, "d"),
        c_name="c_d",
        f_name="F_d",
    )

    return {
        "epsilon": epsilon,
        "C_Q": C_Q,
        "C_u": C_u,
        "C_d": C_d,
        "c_Q": c_Q,
        "c_u": c_u,
        "c_d": c_d,
        "F_Q": F_Q,
        "F_u": F_u,
        "F_d": F_d,
        "rotation_Q": rotation_Q,
        "rotation_u": rotation_u,
        "rotation_d": rotation_d,
        "eig_Q": eig_Q,
        "eig_u": eig_u,
        "eig_d": eig_d,
        "Y_u_bulk_basis": rotation_Q.conjugate().T @ point.Y_u @ rotation_u,
        "Y_d_bulk_basis": rotation_Q.conjugate().T @ point.Y_d @ rotation_d,
        "physical_contract": physical_contract,
    }


def _reference_eq3_summary_from_metadata(
    point: Paper07101869Point,
) -> Paper07101869Eq3ResidualSummary | None:
    metadata = dict(point.metadata)
    payload = metadata.get("benchmark_eq3_example")
    if not isinstance(payload, Mapping):
        return None

    required_keys = (
        "a",
        "r",
        "theta12_deg",
        "theta23_deg",
        "theta13_deg",
        "delta",
    )
    if any(key not in payload for key in required_keys):
        return None

    return evaluate_eq3_diagonal_consistency(
        c_q_eigenvalues=point.q_sector.c_vector,
        c_u_eigenvalues=point.u_sector.c_vector,
        c_d_eigenvalues=point.d_sector.c_vector,
        a=float(payload["a"]),
        r=float(payload["r"]),
        v5km=Paper07101869V5KMParameters.from_degrees(
            theta12_deg=float(payload["theta12_deg"]),
            theta23_deg=float(payload["theta23_deg"]),
            theta13_deg=float(payload["theta13_deg"]),
            delta=float(payload["delta"]),
        ),
    )


def derive_paper_0710_1869_bulk_state(point: Paper07101869Point) -> Paper07101869BulkState:
    """Construct the seeded structural MFV state and reference diagnostics."""
    if not isinstance(point, Paper07101869Point):
        raise ValueError("point must be a Paper07101869Point")

    params = get_warp_params(k=point.k, Lambda_IR=point.Lambda_IR)
    epsilon = float(params["epsilon"])

    C_u = point.Y_u.conjugate().T @ point.Y_u
    C_d = point.Y_d.conjugate().T @ point.Y_d
    C_Q = point.r * (point.Y_u @ point.Y_u.conjugate().T) + (
        point.Y_d @ point.Y_d.conjugate().T
    )

    eig_Q, rotation_Q = _ordered_hermitian_spectrum(C_Q)
    eig_u, rotation_u = _ordered_hermitian_spectrum(C_u)
    eig_d, rotation_d = _ordered_hermitian_spectrum(C_d)

    reference_c_Q = point.q_sector.c_vector
    reference_c_u = point.u_sector.c_vector
    reference_c_d = point.d_sector.c_vector

    quoted_F_Q = point.q_sector.f_vector
    quoted_F_u = point.u_sector.f_vector
    quoted_F_d = point.d_sector.f_vector
    derived_F_Q = _as_real_triplet(
        "derived_F_Q",
        np.asarray(f_IR(reference_c_Q, epsilon), dtype=float),
    )
    derived_F_u = _as_real_triplet(
        "derived_F_u",
        np.asarray(f_IR(reference_c_u, epsilon), dtype=float),
    )
    derived_F_d = _as_real_triplet(
        "derived_F_d",
        np.asarray(f_IR(reference_c_d, epsilon), dtype=float),
    )

    profile_consistency_Q = _build_profile_consistency(
        quoted_f=quoted_F_Q,
        derived_f=derived_F_Q,
        policy_id=point.profile_input_policy_id,
    )
    profile_consistency_u = _build_profile_consistency(
        quoted_f=quoted_F_u,
        derived_f=derived_F_u,
        policy_id=point.profile_input_policy_id,
    )
    profile_consistency_d = _build_profile_consistency(
        quoted_f=quoted_F_d,
        derived_f=derived_F_d,
        policy_id=point.profile_input_policy_id,
    )
    reference_diagnostics_Q = _build_reference_diagnostics(
        quoted_c=reference_c_Q,
        structural_eigenvalues=eig_Q,
        quoted_f=quoted_F_Q,
        geometry_derived_f_from_quoted_c=derived_F_Q,
    )
    reference_diagnostics_u = _build_reference_diagnostics(
        quoted_c=reference_c_u,
        structural_eigenvalues=eig_u,
        quoted_f=quoted_F_u,
        geometry_derived_f_from_quoted_c=derived_F_u,
    )
    reference_diagnostics_d = _build_reference_diagnostics(
        quoted_c=reference_c_d,
        structural_eigenvalues=eig_d,
        quoted_f=quoted_F_d,
        geometry_derived_f_from_quoted_c=derived_F_d,
    )

    Y_u_bulk_basis = rotation_Q.conjugate().T @ point.Y_u @ rotation_u
    Y_d_bulk_basis = rotation_Q.conjugate().T @ point.Y_d @ rotation_d

    return Paper07101869BulkState(
        point=point,
        epsilon=epsilon,
        C_Q=C_Q,
        C_u=C_u,
        C_d=C_d,
        c_Q=reference_c_Q,
        c_u=reference_c_u,
        c_d=reference_c_d,
        F_Q=quoted_F_Q,
        F_u=quoted_F_u,
        F_d=quoted_F_d,
        derived_F_Q=derived_F_Q,
        derived_F_u=derived_F_u,
        derived_F_d=derived_F_d,
        profile_consistency_Q=profile_consistency_Q,
        profile_consistency_u=profile_consistency_u,
        profile_consistency_d=profile_consistency_d,
        reference_diagnostics_Q=reference_diagnostics_Q,
        reference_diagnostics_u=reference_diagnostics_u,
        reference_diagnostics_d=reference_diagnostics_d,
        reference_eq3_summary=_reference_eq3_summary_from_metadata(point),
        rotation_Q=rotation_Q,
        rotation_u=rotation_u,
        rotation_d=rotation_d,
        eig_Q=eig_Q,
        eig_u=eig_u,
        eig_d=eig_d,
        Y_u_bulk_basis=Y_u_bulk_basis,
        Y_d_bulk_basis=Y_d_bulk_basis,
    )


def derive_paper_0710_1869_physical_bulk_state(
    point: Paper07101869PhysicalPoint,
) -> Paper07101869PhysicalBulkState:
    """Construct the hard-gated QS1 physical bulk state from seeded spurions."""
    expected_payload = _build_point_derived_physical_bulk_payload(point)

    return Paper07101869PhysicalBulkState(
        point=point,
        epsilon=float(expected_payload["epsilon"]),
        C_Q=np.asarray(expected_payload["C_Q"], dtype=np.complex128),
        C_u=np.asarray(expected_payload["C_u"], dtype=np.complex128),
        C_d=np.asarray(expected_payload["C_d"], dtype=np.complex128),
        c_Q=np.asarray(expected_payload["c_Q"], dtype=float),
        c_u=np.asarray(expected_payload["c_u"], dtype=float),
        c_d=np.asarray(expected_payload["c_d"], dtype=float),
        F_Q=np.asarray(expected_payload["F_Q"], dtype=float),
        F_u=np.asarray(expected_payload["F_u"], dtype=float),
        F_d=np.asarray(expected_payload["F_d"], dtype=float),
        rotation_Q=np.asarray(expected_payload["rotation_Q"], dtype=np.complex128),
        rotation_u=np.asarray(expected_payload["rotation_u"], dtype=np.complex128),
        rotation_d=np.asarray(expected_payload["rotation_d"], dtype=np.complex128),
        eig_Q=np.asarray(expected_payload["eig_Q"], dtype=float),
        eig_u=np.asarray(expected_payload["eig_u"], dtype=float),
        eig_d=np.asarray(expected_payload["eig_d"], dtype=float),
        Y_u_bulk_basis=np.asarray(expected_payload["Y_u_bulk_basis"], dtype=np.complex128),
        Y_d_bulk_basis=np.asarray(expected_payload["Y_d_bulk_basis"], dtype=np.complex128),
        physical_contract=expected_payload["physical_contract"],
    )


def evaluate_eq3_diagonal_consistency(
    *,
    c_q_eigenvalues: tuple[float, float, float] | list[float] | np.ndarray,
    c_u_eigenvalues: tuple[float, float, float] | list[float] | np.ndarray,
    c_d_eigenvalues: tuple[float, float, float] | list[float] | np.ndarray,
    a: float,
    r: float,
    v5km: np.ndarray | Paper07101869V5KMParameters,
    tolerance: float = 1.0e-12,
) -> Paper07101869Eq3ResidualSummary:
    """Evaluate the Eq. (3) diagonal residual for one paper-side point."""
    resolved_a = require_positive_finite("a", a)
    resolved_r = _require_finite("r", r)
    lhs_diag = _as_real_triplet("c_q_eigenvalues", c_q_eigenvalues)
    c_u = build_diagonal_c_matrix(c_u_eigenvalues)
    c_d = build_diagonal_c_matrix(c_d_eigenvalues)
    resolved_v5km = (
        paper_v5km_matrix(v5km)
        if isinstance(v5km, Paper07101869V5KMParameters)
        else _require_unitary_matrix3(
            "v5km",
            v5km,
            tolerance=PAPER_0710_1869_EQ3_UNITARITY_TOLERANCE,
        )
    )
    rhs_full = resolved_r * resolved_v5km.conjugate().T @ c_u @ resolved_v5km + c_d
    rhs_diag_complex = np.diag(rhs_full)
    max_rhs_imag_abs = float(np.max(np.abs(rhs_diag_complex.imag)))
    rhs_diag = rhs_diag_complex.real.astype(float)
    residual = lhs_diag - (resolved_a * rhs_diag)
    max_abs_residual = float(np.max(np.abs(residual)))
    l2_residual = float(np.linalg.norm(residual))
    resolved_tolerance = require_positive_finite("tolerance", tolerance)
    return Paper07101869Eq3ResidualSummary(
        lhs_diag=lhs_diag,
        rhs_diag=rhs_diag,
        residual=residual,
        a=resolved_a,
        r=resolved_r,
        max_rhs_imag_abs=max_rhs_imag_abs,
        max_abs_residual=max_abs_residual,
        l2_residual=l2_residual,
        exact_agreement=max(max_abs_residual, max_rhs_imag_abs) <= resolved_tolerance,
        tolerance=resolved_tolerance,
    )


def default_paper_0710_1869_table_i_benchmark() -> Paper07101869TableIBenchmark:
    """Return the canonical Table I benchmark bundle."""
    return build_table_i_benchmark_from_inputs(default_paper_0710_1869_table_i_inputs())


def default_paper_0710_1869_eq3_example() -> Paper07101869Eq3Example:
    """Return the quoted Eq. (3) angle/phase example."""
    sourced_example = default_paper_0710_1869_eq3_inputs()
    if not isinstance(sourced_example, Paper07101869Eq3Inputs):
        raise ValueError("default_paper_0710_1869_eq3_inputs() returned the wrong type")
    return Paper07101869Eq3Example(
        a=sourced_example.a,
        r=sourced_example.r,
        v5km_parameters=Paper07101869V5KMParameters.from_degrees(
            theta12_deg=sourced_example.theta12_deg,
            theta23_deg=sourced_example.theta23_deg,
            theta13_deg=sourced_example.theta13_deg,
            delta=sourced_example.delta,
        ),
    )


def evaluate_default_paper_0710_1869_eq3_consistency(
    *,
    tolerance: float = 1.0e-12,
) -> Paper07101869Eq3ResidualSummary:
    """Evaluate the quoted Eq. (3) example against the Table I eigenvalues."""
    benchmark = default_paper_0710_1869_table_i_benchmark()
    example = default_paper_0710_1869_eq3_example()
    return evaluate_eq3_diagonal_consistency(
        c_q_eigenvalues=benchmark.q_sector.c_eigenvalues,
        c_u_eigenvalues=benchmark.u_sector.c_eigenvalues,
        c_d_eigenvalues=benchmark.d_sector.c_eigenvalues,
        a=example.a,
        r=example.r,
        v5km=example.v5km_parameters,
        tolerance=tolerance,
    )


__all__ = [
    "PAPER_0710_1869_BULK_STATE_SCHEMA_ID",
    "PAPER_0710_1869_PHYSICAL_BULK_STATE_SCHEMA_ID",
    "PAPER_0710_1869_PHYSICAL_POINT_SCHEMA_ID",
    "PAPER_0710_1869_EQ3_REFERENCE",
    "PAPER_0710_1869_EQ3_RELATION_ID",
    "PAPER_0710_1869_MODEL_SCHEMA_ID",
    "PAPER_0710_1869_NO_BULK_MASS_MAP_POLICY_ID",
    "PAPER_0710_1869_PHYSICAL_BENCHMARK_SEED_CONSTRUCTION_ID",
    "PAPER_0710_1869_PHYSICAL_CLAIM_LEVEL_ID",
    "PAPER_0710_1869_PHYSICAL_CORE_CONTRACT_ID",
    "PAPER_0710_1869_PHYSICAL_POINT_KIND_ID",
    "PAPER_0710_1869_PHYSICAL_SEED_CONSTRUCTION_ID",
    "PAPER_0710_1869_POINT_SCHEMA_ID",
    "PAPER_0710_1869_POINT_DERIVED_PHYSICAL_STATUS_ID",
    "PAPER_0710_1869_PROFILE_CONSISTENCY_SCHEMA_ID",
    "PAPER_0710_1869_REFERENCE_DIAGNOSTICS_SCHEMA_ID",
    "PAPER_0710_1869_REFERENCE_ONLY_PROFILE_INPUT_POLICY_ID",
    "PAPER_0710_1869_STRUCTURAL_CORE_CONTRACT_ID",
    "PAPER_0710_1869_TABLE_I_REFERENCE",
    "Paper07101869BenchmarkSector",
    "Paper07101869BenchmarkReferenceDiagnostics",
    "Paper07101869BulkState",
    "Paper07101869DiagonalCMatrices",
    "Paper07101869Eq3Example",
    "Paper07101869Eq3ResidualSummary",
    "Paper07101869PhysicalBulkState",
    "Paper07101869PhysicalPoint",
    "Paper07101869Point",
    "Paper07101869ProfileConsistency",
    "Paper07101869RotationParameters",
    "Paper07101869TableIBenchmark",
    "Paper07101869V5KMParameters",
    "build_diagonal_c_matrix",
    "build_paper_0710_1869_physical_point_from_singular_values",
    "build_paper_0710_1869_point_from_singular_values",
    "build_paper_0710_1869_spurion_matrix",
    "build_table_i_benchmark_from_inputs",
    "build_table_i_diagonal_c_matrices",
    "default_paper_0710_1869_eq3_example",
    "default_paper_0710_1869_table_i_benchmark",
    "derive_paper_0710_1869_bulk_state",
    "derive_paper_0710_1869_physical_bulk_state",
    "evaluate_default_paper_0710_1869_eq3_consistency",
    "evaluate_eq3_diagonal_consistency",
    "paper_ckm_like_unitary",
    "paper_v5km_matrix",
]
