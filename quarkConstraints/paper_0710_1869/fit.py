"""Mass/CKM probe helpers for the dedicated 0710.1869 mode."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from diagonalization.diag import SVD

from .model import (
    PAPER_0710_1869_DERIVED_PROFILE_INPUT_POLICY_ID,
    PAPER_0710_1869_NONCANONICAL_CLAIM_LEVEL_ID,
    PAPER_0710_1869_NO_BULK_MASS_MAP_POLICY_ID,
    PAPER_0710_1869_PHYSICAL_CLAIM_LEVEL_ID,
    PAPER_0710_1869_PHYSICAL_CORE_CONTRACT_ID,
    PAPER_0710_1869_POINT_DERIVED_PHYSICAL_STATUS_ID,
    PAPER_0710_1869_SEED_TO_PROFILE_MAPPING_POLICY_ID,
    PAPER_0710_1869_STRUCTURAL_REFERENCE_STATUS_ID,
    Paper07101869BulkState,
    Paper07101869PhysicalBulkState,
    Paper07101869PhysicalPoint,
    Paper07101869Point,
    derive_paper_0710_1869_bulk_state,
    derive_paper_0710_1869_physical_bulk_state,
)
from .validation import (
    require_known_schema_id,
    require_nonempty_identifier,
    require_positive_finite,
)

PAPER_0710_1869_BENCHMARK_REFERENCE_TARGETS_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.benchmark_reference_targets.v1"
)
PAPER_0710_1869_BENCHMARK_REFERENCE_RESULT_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.benchmark_reference_result.v1"
)
PAPER_0710_1869_BENCHMARK_REFERENCE_CONTRACT_ID = (
    "quarkConstraints.paper_0710_1869.contract.benchmark_reference_mass_probe.v1"
)
PAPER_0710_1869_BENCHMARK_REFERENCE_ATTACHMENT_CONVENTION_ID = (
    "quarkConstraints.paper_0710_1869.reference_attachment."
    "rank_order_by_quoted_c_descending.v1"
)
PAPER_0710_1869_BENCHMARK_REFERENCE_PROFILE_SOURCE_ID = (
    "quarkConstraints.paper_0710_1869.profile_source."
    "rank_order_attached_table_i_c_f_reference_pairs_by_descending_c.v1"
)
PAPER_0710_1869_BENCHMARK_REFERENCE_PROVENANCE_MODE_ID = (
    "quarkConstraints.paper_0710_1869.provenance."
    "rank_order_attach_table_i_alias_pairs_by_descending_c_to_sorted_structural_basis.v1"
)
PAPER_0710_1869_BENCHMARK_REFERENCE_TARGET_KIND_ID = (
    "quarkConstraints.paper_0710_1869.target_kind.benchmark_reference_mass_probe.v1"
)
PAPER_0710_1869_BENCHMARK_REFERENCE_RESULT_KIND_ID = (
    "quarkConstraints.paper_0710_1869.result_kind.benchmark_reference_mass_probe.v1"
)
PAPER_0710_1869_BENCHMARK_REFERENCE_UNITARITY_TOLERANCE = 1.0e-10
PAPER_0710_1869_BENCHMARK_REFERENCE_RECONSTRUCTION_TOLERANCE = 1.0e-10
PAPER_0710_1869_PHYSICAL_MASS_CKM_TARGETS_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.physical_mass_ckm_targets.v1"
)
PAPER_0710_1869_PHYSICAL_MASS_CKM_RESULT_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.physical_mass_ckm_result.v1"
)
PAPER_0710_1869_PHYSICAL_MASS_CKM_CONTRACT_ID = (
    "quarkConstraints.paper_0710_1869.contract.physical_mass_ckm_probe.v1"
)
PAPER_0710_1869_PHYSICAL_MASS_CKM_PROFILE_SOURCE_ID = (
    "quarkConstraints.paper_0710_1869.profile_source."
    "point_derived_geometry_from_affine_bulk_masses.v1"
)
PAPER_0710_1869_PHYSICAL_MASS_CKM_PROVENANCE_MODE_ID = (
    "quarkConstraints.paper_0710_1869.provenance."
    "qs1_point_derived_seed_to_profile_mass_ckm_probe.v1"
)
PAPER_0710_1869_PHYSICAL_MASS_CKM_TARGET_KIND_ID = (
    "quarkConstraints.paper_0710_1869.target_kind.qs1_physical_mass_ckm_probe.v1"
)
PAPER_0710_1869_PHYSICAL_MASS_CKM_RESULT_KIND_ID = (
    "quarkConstraints.paper_0710_1869.result_kind.qs1_physical_mass_ckm_probe.v1"
)
PAPER_0710_1869_PHYSICAL_MASS_CKM_UNITARITY_TOLERANCE = 1.0e-10
PAPER_0710_1869_PHYSICAL_MASS_CKM_RECONSTRUCTION_TOLERANCE = 1.0e-10

# Compatibility aliases retained for the current package export surface.
PAPER_0710_1869_TARGETS_SCHEMA_ID = PAPER_0710_1869_BENCHMARK_REFERENCE_TARGETS_SCHEMA_ID
PAPER_0710_1869_FIT_RESULT_SCHEMA_ID = PAPER_0710_1869_BENCHMARK_REFERENCE_RESULT_SCHEMA_ID


def _as_real_vector(
    name: str,
    values: np.ndarray | list[float] | tuple[float, ...],
    expected_shape: tuple[int, ...],
) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.shape != expected_shape:
        raise ValueError(f"{name} must have shape {expected_shape}")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must contain only finite values")
    return arr.astype(float, copy=True)


def _as_positive_ascending_real_triplet(
    name: str,
    values: np.ndarray | list[float] | tuple[float, ...],
) -> np.ndarray:
    arr = _as_real_vector(name, values, (3,))
    if np.any(arr <= 0.0):
        raise ValueError(f"{name} must be strictly positive")
    if np.any(np.diff(arr) <= 0.0):
        raise ValueError(f"{name} must be strictly ascending")
    return arr


def _as_nonnegative_ascending_real_triplet(
    name: str,
    values: np.ndarray | list[float] | tuple[float, ...],
) -> np.ndarray:
    arr = _as_real_vector(name, values, (3,))
    if np.any(arr < 0.0):
        raise ValueError(f"{name} must be non-negative")
    if np.any(np.diff(arr) < 0.0):
        raise ValueError(f"{name} must be ascending")
    return arr


def _as_complex_matrix(name: str, values: np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=np.complex128)
    if arr.shape != (3, 3):
        raise ValueError(f"{name} must have shape (3, 3)")
    if not np.all(np.isfinite(arr.real)) or not np.all(np.isfinite(arr.imag)):
        raise ValueError(f"{name} must contain only finite entries")
    return arr


def _as_positive_real_triplet(
    name: str,
    values: np.ndarray | list[float] | tuple[float, ...],
) -> np.ndarray:
    arr = _as_real_vector(name, values, (3,))
    if np.any(arr <= 0.0):
        raise ValueError(f"{name} must be strictly positive")
    return arr


def _require_exact_identifier(name: str, value: str, *, expected: str) -> str:
    normalized = require_nonempty_identifier(name, value)
    if normalized != expected:
        raise ValueError(f"{name} must be exactly {expected!r}")
    return normalized


def _unitarity_residual(matrix: np.ndarray) -> float:
    arr = _as_complex_matrix("matrix", matrix)
    identity = np.eye(arr.shape[0], dtype=np.complex128)
    return float(np.max(np.abs(arr.conjugate().T @ arr - identity)))


def _require_unitary(name: str, matrix: np.ndarray, *, tolerance: float) -> np.ndarray:
    arr = _as_complex_matrix(name, matrix)
    residual = _unitarity_residual(arr)
    if residual > tolerance:
        raise ValueError(
            f"{name} must be unitary within tolerance {tolerance:.3e}; "
            f"got residual {residual:.3e}"
        )
    return arr


def _max_svd_reconstruction_residual(
    matrix: np.ndarray,
    left_rotation: np.ndarray,
    singular_values: np.ndarray,
    right_rotation: np.ndarray,
) -> float:
    reconstructed = (
        left_rotation
        @ np.diag(_as_nonnegative_ascending_real_triplet("singular_values", singular_values))
        @ right_rotation.conjugate().T
    )
    return float(np.max(np.abs(_as_complex_matrix("matrix", matrix) - reconstructed)))


def _ordered_dirac_svd(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    left_rotation, singular_values, right_rotation = SVD(_as_complex_matrix("matrix", matrix))
    order = np.argsort(singular_values)
    return (
        left_rotation[:, order],
        np.asarray(singular_values[order], dtype=float),
        right_rotation[:, order],
    )


def _require_structural_reference_bulk_state(
    bulk_state: Paper07101869BulkState,
) -> Paper07101869BulkState:
    if not isinstance(bulk_state, Paper07101869BulkState):
        raise ValueError("bulk_state must be a Paper07101869BulkState")
    if bulk_state.claim_level_id != PAPER_0710_1869_NONCANONICAL_CLAIM_LEVEL_ID:
        raise ValueError("bulk_state must carry the noncanonical structural-reference claim level")
    if bulk_state.bulk_mass_map_policy_id != PAPER_0710_1869_NO_BULK_MASS_MAP_POLICY_ID:
        raise ValueError("bulk_state must carry the structural-only bulk-mass map policy")
    if bulk_state.structural_reference_status_id != PAPER_0710_1869_STRUCTURAL_REFERENCE_STATUS_ID:
        raise ValueError("bulk_state must carry the structural-reference status marker")
    return bulk_state


def _require_physical_bulk_state(
    bulk_state: Paper07101869PhysicalBulkState,
) -> Paper07101869PhysicalBulkState:
    if not isinstance(bulk_state, Paper07101869PhysicalBulkState):
        raise ValueError("bulk_state must be a Paper07101869PhysicalBulkState")
    if bulk_state.claim_level_id != PAPER_0710_1869_PHYSICAL_CLAIM_LEVEL_ID:
        raise ValueError("bulk_state must carry the QS1 physical claim level")
    if bulk_state.bulk_mass_map_policy_id != PAPER_0710_1869_SEED_TO_PROFILE_MAPPING_POLICY_ID:
        raise ValueError("bulk_state must carry the frozen seed-to-profile mapping policy")
    if bulk_state.profile_input_policy_id != PAPER_0710_1869_DERIVED_PROFILE_INPUT_POLICY_ID:
        raise ValueError("bulk_state must carry the derived-profile input policy marker")
    if bulk_state.physical_profile_status_id != PAPER_0710_1869_POINT_DERIVED_PHYSICAL_STATUS_ID:
        raise ValueError("bulk_state must carry the point-derived physical status marker")
    if bulk_state.core_contract_id != PAPER_0710_1869_PHYSICAL_CORE_CONTRACT_ID:
        raise ValueError("bulk_state must carry the seeded physical QS1 core contract id")
    return bulk_state


def _require_targets_align_with_bulk_state(
    targets: "Paper07101869BenchmarkReferenceTargets",
    bulk_state: Paper07101869BulkState,
) -> None:
    if targets.claim_level_id != bulk_state.claim_level_id:
        raise ValueError("targets.claim_level_id must match bulk_state.claim_level_id")
    if targets.bulk_mass_map_policy_id != bulk_state.bulk_mass_map_policy_id:
        raise ValueError(
            "targets.bulk_mass_map_policy_id must match bulk_state.bulk_mass_map_policy_id"
        )
    if targets.structural_reference_status_id != bulk_state.structural_reference_status_id:
        raise ValueError(
            "targets.structural_reference_status_id must match "
            "bulk_state.structural_reference_status_id"
        )


def _require_targets_align_with_physical_bulk_state(
    targets: "Paper07101869PhysicalMassCKMTargets",
    bulk_state: Paper07101869PhysicalBulkState,
) -> None:
    if targets.claim_level_id != bulk_state.claim_level_id:
        raise ValueError("targets.claim_level_id must match bulk_state.claim_level_id")
    if targets.bulk_mass_map_policy_id != bulk_state.bulk_mass_map_policy_id:
        raise ValueError(
            "targets.bulk_mass_map_policy_id must match bulk_state.bulk_mass_map_policy_id"
        )
    if targets.profile_input_policy_id != bulk_state.profile_input_policy_id:
        raise ValueError(
            "targets.profile_input_policy_id must match bulk_state.profile_input_policy_id"
        )
    if targets.physical_profile_status_id != bulk_state.physical_profile_status_id:
        raise ValueError(
            "targets.physical_profile_status_id must match "
            "bulk_state.physical_profile_status_id"
        )
    if targets.core_contract_id != bulk_state.core_contract_id:
        raise ValueError("targets.core_contract_id must match bulk_state.core_contract_id")


def _rank_order_attach_reference_profiles(
    *,
    sector_name: str,
    quoted_c: np.ndarray,
    quoted_f: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Attach quoted Table I reference pairs to the sorted structural basis.

    The structural eigensystems in ``Paper07101869BulkState`` are sorted by
    ascending structural eigenvalue. This benchmark-reference probe is
    intentionally noncanonical and uses one explicit deterministic convention:
    sort each quoted Table I ``(c_i, f_i)`` pair by quoted ``c_i`` descending,
    preserve each quoted pair, and attach the reordered ``f_i`` aliases
    slot-by-slot onto that sorted structural basis.
    """

    quoted_c_arr = _as_real_vector(f"{sector_name}_quoted_c", quoted_c, (3,))
    quoted_f_arr = _as_positive_real_triplet(f"{sector_name}_quoted_f", quoted_f)
    order = np.lexsort((np.arange(quoted_c_arr.size), -quoted_c_arr))
    return quoted_c_arr[order], quoted_f_arr[order]


def paper_jarlskog_invariant(ckm: np.ndarray) -> float:
    """Return the CKM Jarlskog invariant for benchmark-reference probes."""
    v = _as_complex_matrix("ckm", ckm)
    return float(np.imag(v[0, 1] * v[1, 2] * np.conjugate(v[0, 2]) * np.conjugate(v[1, 1])))


def paper_ckm_observables(ckm: np.ndarray) -> np.ndarray:
    """Return the reduced CKM observable set for benchmark-reference probes."""
    v = _as_complex_matrix("ckm", ckm)
    return np.array(
        [abs(v[0, 1]), abs(v[1, 2]), abs(v[0, 2]), paper_jarlskog_invariant(v)],
        dtype=float,
    )


def _mass_probe_observables(M_u: np.ndarray, M_d: np.ndarray) -> dict[str, np.ndarray]:
    U_L_u, masses_up, U_R_u = _ordered_dirac_svd(M_u)
    U_L_d, masses_down, U_R_d = _ordered_dirac_svd(M_d)
    ckm = U_L_u.conjugate().T @ U_L_d
    return {
        "U_L_u": U_L_u,
        "U_R_u": U_R_u,
        "U_L_d": U_L_d,
        "U_R_d": U_R_d,
        "masses_up": masses_up,
        "masses_down": masses_down,
        "ckm": ckm,
    }


def _mass_probe_residuals(
    *,
    masses_up: np.ndarray,
    masses_down: np.ndarray,
    ckm: np.ndarray,
    target_up_masses: np.ndarray,
    target_down_masses: np.ndarray,
    target_ckm: np.ndarray,
    ckm_unitarity_tolerance: float,
) -> dict[str, np.ndarray | float]:
    masses_up_arr = _as_nonnegative_ascending_real_triplet("masses_up", masses_up)
    masses_down_arr = _as_nonnegative_ascending_real_triplet("masses_down", masses_down)
    ckm_arr = _require_unitary("ckm", ckm, tolerance=ckm_unitarity_tolerance)

    up_log_residuals = np.log(np.maximum(masses_up_arr, 1.0e-30) / target_up_masses)
    down_log_residuals = np.log(np.maximum(masses_down_arr, 1.0e-30) / target_down_masses)
    ckm_abs_residuals = np.abs(np.abs(ckm_arr) - np.abs(target_ckm))
    target_ckm_observables = paper_ckm_observables(target_ckm)
    ckm_observable_residuals = (
        paper_ckm_observables(ckm_arr) - target_ckm_observables
    ) / np.maximum(np.abs(target_ckm_observables), 1.0e-8)

    total_score = float(
        np.sqrt(
            np.mean(
                np.concatenate(
                    [up_log_residuals, down_log_residuals, ckm_observable_residuals]
                )
                ** 2
            )
        )
    )
    return {
        "up_log_residuals": up_log_residuals,
        "down_log_residuals": down_log_residuals,
        "ckm_abs_residuals": ckm_abs_residuals,
        "ckm_observable_residuals": ckm_observable_residuals,
        "total_score": total_score,
    }


@dataclass(frozen=True)
class Paper07101869BenchmarkReferenceTargets:
    """Targets for scoring a noncanonical benchmark-reference mass probe."""

    up_masses: np.ndarray
    down_masses: np.ndarray
    ckm: np.ndarray
    label: str = "benchmark_reference_targets"
    claim_level_id: str = PAPER_0710_1869_NONCANONICAL_CLAIM_LEVEL_ID
    bulk_mass_map_policy_id: str = PAPER_0710_1869_NO_BULK_MASS_MAP_POLICY_ID
    structural_reference_status_id: str = PAPER_0710_1869_STRUCTURAL_REFERENCE_STATUS_ID
    target_kind_id: str = PAPER_0710_1869_BENCHMARK_REFERENCE_TARGET_KIND_ID
    probe_contract_id: str = PAPER_0710_1869_BENCHMARK_REFERENCE_CONTRACT_ID
    attachment_convention_id: str = PAPER_0710_1869_BENCHMARK_REFERENCE_ATTACHMENT_CONVENTION_ID
    profile_source_id: str = PAPER_0710_1869_BENCHMARK_REFERENCE_PROFILE_SOURCE_ID
    provenance_mode_id: str = PAPER_0710_1869_BENCHMARK_REFERENCE_PROVENANCE_MODE_ID
    ckm_unitarity_tolerance: float = PAPER_0710_1869_BENCHMARK_REFERENCE_UNITARITY_TOLERANCE
    schema_id: str = PAPER_0710_1869_BENCHMARK_REFERENCE_TARGETS_SCHEMA_ID

    def __post_init__(self) -> None:
        tolerance = require_positive_finite(
            "ckm_unitarity_tolerance",
            self.ckm_unitarity_tolerance,
        )
        object.__setattr__(
            self,
            "up_masses",
            _as_positive_ascending_real_triplet("up_masses", self.up_masses),
        )
        object.__setattr__(
            self,
            "down_masses",
            _as_positive_ascending_real_triplet("down_masses", self.down_masses),
        )
        object.__setattr__(
            self,
            "ckm",
            _require_unitary("ckm", self.ckm, tolerance=tolerance),
        )
        object.__setattr__(self, "label", require_nonempty_identifier("label", self.label))
        object.__setattr__(
            self,
            "claim_level_id",
            _require_exact_identifier(
                "claim_level_id",
                self.claim_level_id,
                expected=PAPER_0710_1869_NONCANONICAL_CLAIM_LEVEL_ID,
            ),
        )
        object.__setattr__(
            self,
            "bulk_mass_map_policy_id",
            _require_exact_identifier(
                "bulk_mass_map_policy_id",
                self.bulk_mass_map_policy_id,
                expected=PAPER_0710_1869_NO_BULK_MASS_MAP_POLICY_ID,
            ),
        )
        object.__setattr__(
            self,
            "structural_reference_status_id",
            _require_exact_identifier(
                "structural_reference_status_id",
                self.structural_reference_status_id,
                expected=PAPER_0710_1869_STRUCTURAL_REFERENCE_STATUS_ID,
            ),
        )
        object.__setattr__(
            self,
            "target_kind_id",
            _require_exact_identifier(
                "target_kind_id",
                self.target_kind_id,
                expected=PAPER_0710_1869_BENCHMARK_REFERENCE_TARGET_KIND_ID,
            ),
        )
        object.__setattr__(
            self,
            "probe_contract_id",
            _require_exact_identifier(
                "probe_contract_id",
                self.probe_contract_id,
                expected=PAPER_0710_1869_BENCHMARK_REFERENCE_CONTRACT_ID,
            ),
        )
        object.__setattr__(
            self,
            "attachment_convention_id",
            _require_exact_identifier(
                "attachment_convention_id",
                self.attachment_convention_id,
                expected=PAPER_0710_1869_BENCHMARK_REFERENCE_ATTACHMENT_CONVENTION_ID,
            ),
        )
        object.__setattr__(
            self,
            "profile_source_id",
            _require_exact_identifier(
                "profile_source_id",
                self.profile_source_id,
                expected=PAPER_0710_1869_BENCHMARK_REFERENCE_PROFILE_SOURCE_ID,
            ),
        )
        object.__setattr__(
            self,
            "provenance_mode_id",
            _require_exact_identifier(
                "provenance_mode_id",
                self.provenance_mode_id,
                expected=PAPER_0710_1869_BENCHMARK_REFERENCE_PROVENANCE_MODE_ID,
            ),
        )
        object.__setattr__(self, "ckm_unitarity_tolerance", tolerance)
        object.__setattr__(
            self,
            "schema_id",
            require_known_schema_id(
                "schema_id",
                self.schema_id,
                expected=PAPER_0710_1869_BENCHMARK_REFERENCE_TARGETS_SCHEMA_ID,
            ),
        )

    @property
    def abs_ckm(self) -> np.ndarray:
        return np.abs(self.ckm)

    @property
    def ckm_observables(self) -> np.ndarray:
        return paper_ckm_observables(self.ckm)

    @property
    def max_ckm_unitarity_residual(self) -> float:
        return _unitarity_residual(self.ckm)


@dataclass(frozen=True)
class Paper07101869BenchmarkReferenceMassProbeResult:
    """Mass/CKM probe built from rank-order attached Table I reference profiles."""

    bulk_state: Paper07101869BulkState
    M_u: np.ndarray
    M_d: np.ndarray
    U_L_u: np.ndarray
    U_R_u: np.ndarray
    U_L_d: np.ndarray
    U_R_d: np.ndarray
    masses_up: np.ndarray
    masses_down: np.ndarray
    ckm: np.ndarray
    up_log_residuals: np.ndarray | None = None
    down_log_residuals: np.ndarray | None = None
    ckm_abs_residuals: np.ndarray | None = None
    ckm_observable_residuals: np.ndarray | None = None
    total_score: float | None = None
    target_label: str | None = None
    target_schema_id: str | None = None
    claim_level_id: str = PAPER_0710_1869_NONCANONICAL_CLAIM_LEVEL_ID
    bulk_mass_map_policy_id: str = PAPER_0710_1869_NO_BULK_MASS_MAP_POLICY_ID
    structural_reference_status_id: str = PAPER_0710_1869_STRUCTURAL_REFERENCE_STATUS_ID
    result_kind_id: str = PAPER_0710_1869_BENCHMARK_REFERENCE_RESULT_KIND_ID
    probe_contract_id: str = PAPER_0710_1869_BENCHMARK_REFERENCE_CONTRACT_ID
    attachment_convention_id: str = PAPER_0710_1869_BENCHMARK_REFERENCE_ATTACHMENT_CONVENTION_ID
    profile_source_id: str = PAPER_0710_1869_BENCHMARK_REFERENCE_PROFILE_SOURCE_ID
    provenance_mode_id: str = PAPER_0710_1869_BENCHMARK_REFERENCE_PROVENANCE_MODE_ID
    validation_tolerance: float = PAPER_0710_1869_BENCHMARK_REFERENCE_UNITARITY_TOLERANCE
    max_ckm_unitarity_residual: float | None = None
    max_rotation_unitarity_residual: float | None = None
    max_reconstruction_residual: float | None = None
    schema_id: str = PAPER_0710_1869_BENCHMARK_REFERENCE_RESULT_SCHEMA_ID

    def __post_init__(self) -> None:
        _require_structural_reference_bulk_state(self.bulk_state)
        tolerance = require_positive_finite("validation_tolerance", self.validation_tolerance)
        object.__setattr__(self, "M_u", _as_complex_matrix("M_u", self.M_u))
        object.__setattr__(self, "M_d", _as_complex_matrix("M_d", self.M_d))
        object.__setattr__(
            self,
            "U_L_u",
            _require_unitary("U_L_u", self.U_L_u, tolerance=tolerance),
        )
        object.__setattr__(
            self,
            "U_R_u",
            _require_unitary("U_R_u", self.U_R_u, tolerance=tolerance),
        )
        object.__setattr__(
            self,
            "U_L_d",
            _require_unitary("U_L_d", self.U_L_d, tolerance=tolerance),
        )
        object.__setattr__(
            self,
            "U_R_d",
            _require_unitary("U_R_d", self.U_R_d, tolerance=tolerance),
        )
        object.__setattr__(
            self,
            "masses_up",
            _as_nonnegative_ascending_real_triplet("masses_up", self.masses_up),
        )
        object.__setattr__(
            self,
            "masses_down",
            _as_nonnegative_ascending_real_triplet("masses_down", self.masses_down),
        )
        object.__setattr__(
            self,
            "ckm",
            _require_unitary("ckm", self.ckm, tolerance=tolerance),
        )
        if self.up_log_residuals is not None:
            object.__setattr__(
                self,
                "up_log_residuals",
                _as_real_vector("up_log_residuals", self.up_log_residuals, (3,)),
            )
        if self.down_log_residuals is not None:
            object.__setattr__(
                self,
                "down_log_residuals",
                _as_real_vector("down_log_residuals", self.down_log_residuals, (3,)),
            )
        if self.ckm_abs_residuals is not None:
            object.__setattr__(
                self,
                "ckm_abs_residuals",
                _as_real_vector("ckm_abs_residuals", self.ckm_abs_residuals, (3, 3)),
            )
        if self.ckm_observable_residuals is not None:
            object.__setattr__(
                self,
                "ckm_observable_residuals",
                _as_real_vector("ckm_observable_residuals", self.ckm_observable_residuals, (4,)),
            )
        if self.total_score is not None:
            score = float(self.total_score)
            if not np.isfinite(score):
                raise ValueError("total_score must be finite")
            object.__setattr__(self, "total_score", score)
        if self.target_label is not None:
            object.__setattr__(
                self,
                "target_label",
                require_nonempty_identifier("target_label", self.target_label),
            )
        if self.target_schema_id is not None:
            object.__setattr__(
                self,
                "target_schema_id",
                require_known_schema_id(
                    "target_schema_id",
                    self.target_schema_id,
                    expected=PAPER_0710_1869_BENCHMARK_REFERENCE_TARGETS_SCHEMA_ID,
                ),
            )
        object.__setattr__(
            self,
            "claim_level_id",
            _require_exact_identifier(
                "claim_level_id",
                self.claim_level_id,
                expected=PAPER_0710_1869_NONCANONICAL_CLAIM_LEVEL_ID,
            ),
        )
        object.__setattr__(
            self,
            "bulk_mass_map_policy_id",
            _require_exact_identifier(
                "bulk_mass_map_policy_id",
                self.bulk_mass_map_policy_id,
                expected=PAPER_0710_1869_NO_BULK_MASS_MAP_POLICY_ID,
            ),
        )
        object.__setattr__(
            self,
            "structural_reference_status_id",
            _require_exact_identifier(
                "structural_reference_status_id",
                self.structural_reference_status_id,
                expected=PAPER_0710_1869_STRUCTURAL_REFERENCE_STATUS_ID,
            ),
        )
        object.__setattr__(
            self,
            "result_kind_id",
            _require_exact_identifier(
                "result_kind_id",
                self.result_kind_id,
                expected=PAPER_0710_1869_BENCHMARK_REFERENCE_RESULT_KIND_ID,
            ),
        )
        object.__setattr__(
            self,
            "probe_contract_id",
            _require_exact_identifier(
                "probe_contract_id",
                self.probe_contract_id,
                expected=PAPER_0710_1869_BENCHMARK_REFERENCE_CONTRACT_ID,
            ),
        )
        object.__setattr__(
            self,
            "attachment_convention_id",
            _require_exact_identifier(
                "attachment_convention_id",
                self.attachment_convention_id,
                expected=PAPER_0710_1869_BENCHMARK_REFERENCE_ATTACHMENT_CONVENTION_ID,
            ),
        )
        object.__setattr__(
            self,
            "profile_source_id",
            _require_exact_identifier(
                "profile_source_id",
                self.profile_source_id,
                expected=PAPER_0710_1869_BENCHMARK_REFERENCE_PROFILE_SOURCE_ID,
            ),
        )
        object.__setattr__(
            self,
            "provenance_mode_id",
            _require_exact_identifier(
                "provenance_mode_id",
                self.provenance_mode_id,
                expected=PAPER_0710_1869_BENCHMARK_REFERENCE_PROVENANCE_MODE_ID,
            ),
        )
        object.__setattr__(self, "validation_tolerance", tolerance)
        object.__setattr__(
            self,
            "schema_id",
            require_known_schema_id(
                "schema_id",
                self.schema_id,
                expected=PAPER_0710_1869_BENCHMARK_REFERENCE_RESULT_SCHEMA_ID,
            ),
        )
        if self.claim_level_id != self.bulk_state.claim_level_id:
            raise ValueError("result claim_level_id must match bulk_state.claim_level_id")
        if self.bulk_mass_map_policy_id != self.bulk_state.bulk_mass_map_policy_id:
            raise ValueError(
                "result bulk_mass_map_policy_id must match bulk_state.bulk_mass_map_policy_id"
            )
        if self.structural_reference_status_id != self.bulk_state.structural_reference_status_id:
            raise ValueError(
                "result structural_reference_status_id must match "
                "bulk_state.structural_reference_status_id"
            )

        ckm_from_left_rotations = self.U_L_u.conjugate().T @ self.U_L_d
        ckm_consistency_residual = float(np.max(np.abs(self.ckm - ckm_from_left_rotations)))
        if ckm_consistency_residual > tolerance:
            raise ValueError(
                "ckm must match U_L_u^dagger U_L_d within validation_tolerance; "
                f"got residual {ckm_consistency_residual:.3e}"
            )

        computed_ckm_residual = _unitarity_residual(self.ckm)
        computed_rotation_residual = float(
            max(
                _unitarity_residual(self.U_L_u),
                _unitarity_residual(self.U_R_u),
                _unitarity_residual(self.U_L_d),
                _unitarity_residual(self.U_R_d),
            )
        )
        computed_reconstruction_residual = float(
            max(
                _max_svd_reconstruction_residual(
                    self.M_u,
                    self.U_L_u,
                    self.masses_up,
                    self.U_R_u,
                ),
                _max_svd_reconstruction_residual(
                    self.M_d,
                    self.U_L_d,
                    self.masses_down,
                    self.U_R_d,
                ),
            )
        )

        object.__setattr__(self, "max_ckm_unitarity_residual", computed_ckm_residual)
        object.__setattr__(self, "max_rotation_unitarity_residual", computed_rotation_residual)
        object.__setattr__(self, "max_reconstruction_residual", computed_reconstruction_residual)

        if computed_reconstruction_residual > tolerance:
            raise ValueError(
                "SVD reconstruction must hold within validation_tolerance; "
                f"got residual {computed_reconstruction_residual:.3e}"
            )

    @property
    def state(self) -> Paper07101869BulkState:
        return self.bulk_state

    @property
    def point(self) -> Paper07101869Point:
        return self.bulk_state.point

    @property
    def score(self) -> float:
        return float(0.0 if self.total_score is None else self.total_score)

    @property
    def mass_residuals_up(self) -> np.ndarray:
        return np.zeros(3) if self.up_log_residuals is None else self.up_log_residuals

    @property
    def mass_residuals_down(self) -> np.ndarray:
        return np.zeros(3) if self.down_log_residuals is None else self.down_log_residuals

    @property
    def ckm_residuals(self) -> np.ndarray:
        if self.ckm_observable_residuals is not None:
            return self.ckm_observable_residuals
        if self.ckm_abs_residuals is None:
            return np.zeros(4)
        return np.array(
            [
                self.ckm_abs_residuals[0, 1],
                self.ckm_abs_residuals[1, 2],
                self.ckm_abs_residuals[0, 2],
                0.0,
            ],
            dtype=float,
        )

    @property
    def ckm_matrix(self) -> np.ndarray:
        return self.ckm

    @property
    def ckm_observables(self) -> np.ndarray:
        return paper_ckm_observables(self.ckm)

    @property
    def residual_norm(self) -> float:
        return float(
            np.linalg.norm(
                np.concatenate(
                    [self.mass_residuals_up, self.mass_residuals_down, self.ckm_residuals]
                )
            )
        )

    def summary(self) -> str:
        _, attached_f_q = _rank_order_attach_reference_profiles(
            sector_name="Q",
            quoted_c=self.bulk_state.c_Q,
            quoted_f=self.bulk_state.F_Q,
        )
        _, attached_f_u = _rank_order_attach_reference_profiles(
            sector_name="u",
            quoted_c=self.bulk_state.c_u,
            quoted_f=self.bulk_state.F_u,
        )
        _, attached_f_d = _rank_order_attach_reference_profiles(
            sector_name="d",
            quoted_c=self.bulk_state.c_d,
            quoted_f=self.bulk_state.F_d,
        )
        lines = [
            "=" * 60,
            f"Paper 0710.1869 benchmark-reference mass probe: {self.point.label}",
            "=" * 60,
            "Uses a noncanonical rank-order attachment of quoted Table I (c_i, f_i) "
            "reference pairs sorted by quoted c descending onto the sorted "
            "structural basis; not a seed-derived physical bulk-profile fit.",
            f"r = {self.point.r:.4f}",
            f"epsilon = {self.bulk_state.epsilon:.3e}",
            f"c_Q ref = {np.array2string(self.bulk_state.c_Q, precision=4)}",
            f"c_u ref = {np.array2string(self.bulk_state.c_u, precision=4)}",
            f"c_d ref = {np.array2string(self.bulk_state.c_d, precision=4)}",
            f"F_Q attached(ref,c-desc) = {np.array2string(attached_f_q, precision=5)}",
            f"F_u attached(ref,c-desc) = {np.array2string(attached_f_u, precision=5)}",
            f"F_d attached(ref,c-desc) = {np.array2string(attached_f_d, precision=5)}",
            f"m_u probe = {np.array2string(self.masses_up, precision=5)}",
            f"m_d probe = {np.array2string(self.masses_down, precision=5)}",
            f"CKM obs = {np.array2string(self.ckm_observables, precision=6)}",
            f"attachment convention = {self.attachment_convention_id}",
            f"profile source = {self.profile_source_id}",
        ]
        if self.total_score is not None:
            lines.append(f"probe score = {self.total_score:.6e}")
        return "\n".join(lines)


@dataclass(frozen=True)
class Paper07101869PhysicalMassCKMTargets:
    """Targets for the hard-gated QS1 physical masses/CKM probe."""

    up_masses: np.ndarray
    down_masses: np.ndarray
    ckm: np.ndarray
    label: str = "physical_mass_ckm_targets"
    claim_level_id: str = PAPER_0710_1869_PHYSICAL_CLAIM_LEVEL_ID
    bulk_mass_map_policy_id: str = PAPER_0710_1869_SEED_TO_PROFILE_MAPPING_POLICY_ID
    profile_input_policy_id: str = PAPER_0710_1869_DERIVED_PROFILE_INPUT_POLICY_ID
    physical_profile_status_id: str = PAPER_0710_1869_POINT_DERIVED_PHYSICAL_STATUS_ID
    core_contract_id: str = PAPER_0710_1869_PHYSICAL_CORE_CONTRACT_ID
    target_kind_id: str = PAPER_0710_1869_PHYSICAL_MASS_CKM_TARGET_KIND_ID
    probe_contract_id: str = PAPER_0710_1869_PHYSICAL_MASS_CKM_CONTRACT_ID
    profile_source_id: str = PAPER_0710_1869_PHYSICAL_MASS_CKM_PROFILE_SOURCE_ID
    provenance_mode_id: str = PAPER_0710_1869_PHYSICAL_MASS_CKM_PROVENANCE_MODE_ID
    ckm_unitarity_tolerance: float = PAPER_0710_1869_PHYSICAL_MASS_CKM_UNITARITY_TOLERANCE
    schema_id: str = PAPER_0710_1869_PHYSICAL_MASS_CKM_TARGETS_SCHEMA_ID

    def __post_init__(self) -> None:
        tolerance = require_positive_finite(
            "ckm_unitarity_tolerance",
            self.ckm_unitarity_tolerance,
        )
        object.__setattr__(
            self,
            "up_masses",
            _as_positive_ascending_real_triplet("up_masses", self.up_masses),
        )
        object.__setattr__(
            self,
            "down_masses",
            _as_positive_ascending_real_triplet("down_masses", self.down_masses),
        )
        object.__setattr__(
            self,
            "ckm",
            _require_unitary("ckm", self.ckm, tolerance=tolerance),
        )
        object.__setattr__(self, "label", require_nonempty_identifier("label", self.label))
        object.__setattr__(
            self,
            "claim_level_id",
            _require_exact_identifier(
                "claim_level_id",
                self.claim_level_id,
                expected=PAPER_0710_1869_PHYSICAL_CLAIM_LEVEL_ID,
            ),
        )
        object.__setattr__(
            self,
            "bulk_mass_map_policy_id",
            _require_exact_identifier(
                "bulk_mass_map_policy_id",
                self.bulk_mass_map_policy_id,
                expected=PAPER_0710_1869_SEED_TO_PROFILE_MAPPING_POLICY_ID,
            ),
        )
        object.__setattr__(
            self,
            "profile_input_policy_id",
            _require_exact_identifier(
                "profile_input_policy_id",
                self.profile_input_policy_id,
                expected=PAPER_0710_1869_DERIVED_PROFILE_INPUT_POLICY_ID,
            ),
        )
        object.__setattr__(
            self,
            "physical_profile_status_id",
            _require_exact_identifier(
                "physical_profile_status_id",
                self.physical_profile_status_id,
                expected=PAPER_0710_1869_POINT_DERIVED_PHYSICAL_STATUS_ID,
            ),
        )
        object.__setattr__(
            self,
            "core_contract_id",
            _require_exact_identifier(
                "core_contract_id",
                self.core_contract_id,
                expected=PAPER_0710_1869_PHYSICAL_CORE_CONTRACT_ID,
            ),
        )
        object.__setattr__(
            self,
            "target_kind_id",
            _require_exact_identifier(
                "target_kind_id",
                self.target_kind_id,
                expected=PAPER_0710_1869_PHYSICAL_MASS_CKM_TARGET_KIND_ID,
            ),
        )
        object.__setattr__(
            self,
            "probe_contract_id",
            _require_exact_identifier(
                "probe_contract_id",
                self.probe_contract_id,
                expected=PAPER_0710_1869_PHYSICAL_MASS_CKM_CONTRACT_ID,
            ),
        )
        object.__setattr__(
            self,
            "profile_source_id",
            _require_exact_identifier(
                "profile_source_id",
                self.profile_source_id,
                expected=PAPER_0710_1869_PHYSICAL_MASS_CKM_PROFILE_SOURCE_ID,
            ),
        )
        object.__setattr__(
            self,
            "provenance_mode_id",
            _require_exact_identifier(
                "provenance_mode_id",
                self.provenance_mode_id,
                expected=PAPER_0710_1869_PHYSICAL_MASS_CKM_PROVENANCE_MODE_ID,
            ),
        )
        object.__setattr__(self, "ckm_unitarity_tolerance", tolerance)
        object.__setattr__(
            self,
            "schema_id",
            require_known_schema_id(
                "schema_id",
                self.schema_id,
                expected=PAPER_0710_1869_PHYSICAL_MASS_CKM_TARGETS_SCHEMA_ID,
            ),
        )

    @property
    def abs_ckm(self) -> np.ndarray:
        return np.abs(self.ckm)

    @property
    def ckm_observables(self) -> np.ndarray:
        return paper_ckm_observables(self.ckm)

    @property
    def max_ckm_unitarity_residual(self) -> float:
        return _unitarity_residual(self.ckm)


@dataclass(frozen=True)
class Paper07101869PhysicalMassCKMProbeResult:
    """QS1 physical masses/CKM probe built from point-derived physical profiles."""

    bulk_state: Paper07101869PhysicalBulkState
    M_u: np.ndarray
    M_d: np.ndarray
    U_L_u: np.ndarray
    U_R_u: np.ndarray
    U_L_d: np.ndarray
    U_R_d: np.ndarray
    masses_up: np.ndarray
    masses_down: np.ndarray
    ckm: np.ndarray
    up_log_residuals: np.ndarray | None = None
    down_log_residuals: np.ndarray | None = None
    ckm_abs_residuals: np.ndarray | None = None
    ckm_observable_residuals: np.ndarray | None = None
    total_score: float | None = None
    target_label: str | None = None
    target_schema_id: str | None = None
    claim_level_id: str = PAPER_0710_1869_PHYSICAL_CLAIM_LEVEL_ID
    bulk_mass_map_policy_id: str = PAPER_0710_1869_SEED_TO_PROFILE_MAPPING_POLICY_ID
    profile_input_policy_id: str = PAPER_0710_1869_DERIVED_PROFILE_INPUT_POLICY_ID
    physical_profile_status_id: str = PAPER_0710_1869_POINT_DERIVED_PHYSICAL_STATUS_ID
    core_contract_id: str = PAPER_0710_1869_PHYSICAL_CORE_CONTRACT_ID
    result_kind_id: str = PAPER_0710_1869_PHYSICAL_MASS_CKM_RESULT_KIND_ID
    probe_contract_id: str = PAPER_0710_1869_PHYSICAL_MASS_CKM_CONTRACT_ID
    profile_source_id: str = PAPER_0710_1869_PHYSICAL_MASS_CKM_PROFILE_SOURCE_ID
    provenance_mode_id: str = PAPER_0710_1869_PHYSICAL_MASS_CKM_PROVENANCE_MODE_ID
    validation_tolerance: float = PAPER_0710_1869_PHYSICAL_MASS_CKM_UNITARITY_TOLERANCE
    max_ckm_unitarity_residual: float | None = None
    max_rotation_unitarity_residual: float | None = None
    max_reconstruction_residual: float | None = None
    schema_id: str = PAPER_0710_1869_PHYSICAL_MASS_CKM_RESULT_SCHEMA_ID

    def __post_init__(self) -> None:
        _require_physical_bulk_state(self.bulk_state)
        tolerance = require_positive_finite("validation_tolerance", self.validation_tolerance)
        object.__setattr__(self, "M_u", _as_complex_matrix("M_u", self.M_u))
        object.__setattr__(self, "M_d", _as_complex_matrix("M_d", self.M_d))
        object.__setattr__(
            self,
            "U_L_u",
            _require_unitary("U_L_u", self.U_L_u, tolerance=tolerance),
        )
        object.__setattr__(
            self,
            "U_R_u",
            _require_unitary("U_R_u", self.U_R_u, tolerance=tolerance),
        )
        object.__setattr__(
            self,
            "U_L_d",
            _require_unitary("U_L_d", self.U_L_d, tolerance=tolerance),
        )
        object.__setattr__(
            self,
            "U_R_d",
            _require_unitary("U_R_d", self.U_R_d, tolerance=tolerance),
        )
        object.__setattr__(
            self,
            "masses_up",
            _as_nonnegative_ascending_real_triplet("masses_up", self.masses_up),
        )
        object.__setattr__(
            self,
            "masses_down",
            _as_nonnegative_ascending_real_triplet("masses_down", self.masses_down),
        )
        object.__setattr__(
            self,
            "ckm",
            _require_unitary("ckm", self.ckm, tolerance=tolerance),
        )
        if self.up_log_residuals is not None:
            object.__setattr__(
                self,
                "up_log_residuals",
                _as_real_vector("up_log_residuals", self.up_log_residuals, (3,)),
            )
        if self.down_log_residuals is not None:
            object.__setattr__(
                self,
                "down_log_residuals",
                _as_real_vector("down_log_residuals", self.down_log_residuals, (3,)),
            )
        if self.ckm_abs_residuals is not None:
            object.__setattr__(
                self,
                "ckm_abs_residuals",
                _as_real_vector("ckm_abs_residuals", self.ckm_abs_residuals, (3, 3)),
            )
        if self.ckm_observable_residuals is not None:
            object.__setattr__(
                self,
                "ckm_observable_residuals",
                _as_real_vector("ckm_observable_residuals", self.ckm_observable_residuals, (4,)),
            )
        if self.total_score is not None:
            score = float(self.total_score)
            if not np.isfinite(score):
                raise ValueError("total_score must be finite")
            object.__setattr__(self, "total_score", score)
        if self.target_label is not None:
            object.__setattr__(
                self,
                "target_label",
                require_nonempty_identifier("target_label", self.target_label),
            )
        if self.target_schema_id is not None:
            object.__setattr__(
                self,
                "target_schema_id",
                require_known_schema_id(
                    "target_schema_id",
                    self.target_schema_id,
                    expected=PAPER_0710_1869_PHYSICAL_MASS_CKM_TARGETS_SCHEMA_ID,
                ),
            )
        object.__setattr__(
            self,
            "claim_level_id",
            _require_exact_identifier(
                "claim_level_id",
                self.claim_level_id,
                expected=PAPER_0710_1869_PHYSICAL_CLAIM_LEVEL_ID,
            ),
        )
        object.__setattr__(
            self,
            "bulk_mass_map_policy_id",
            _require_exact_identifier(
                "bulk_mass_map_policy_id",
                self.bulk_mass_map_policy_id,
                expected=PAPER_0710_1869_SEED_TO_PROFILE_MAPPING_POLICY_ID,
            ),
        )
        object.__setattr__(
            self,
            "profile_input_policy_id",
            _require_exact_identifier(
                "profile_input_policy_id",
                self.profile_input_policy_id,
                expected=PAPER_0710_1869_DERIVED_PROFILE_INPUT_POLICY_ID,
            ),
        )
        object.__setattr__(
            self,
            "physical_profile_status_id",
            _require_exact_identifier(
                "physical_profile_status_id",
                self.physical_profile_status_id,
                expected=PAPER_0710_1869_POINT_DERIVED_PHYSICAL_STATUS_ID,
            ),
        )
        object.__setattr__(
            self,
            "core_contract_id",
            _require_exact_identifier(
                "core_contract_id",
                self.core_contract_id,
                expected=PAPER_0710_1869_PHYSICAL_CORE_CONTRACT_ID,
            ),
        )
        object.__setattr__(
            self,
            "result_kind_id",
            _require_exact_identifier(
                "result_kind_id",
                self.result_kind_id,
                expected=PAPER_0710_1869_PHYSICAL_MASS_CKM_RESULT_KIND_ID,
            ),
        )
        object.__setattr__(
            self,
            "probe_contract_id",
            _require_exact_identifier(
                "probe_contract_id",
                self.probe_contract_id,
                expected=PAPER_0710_1869_PHYSICAL_MASS_CKM_CONTRACT_ID,
            ),
        )
        object.__setattr__(
            self,
            "profile_source_id",
            _require_exact_identifier(
                "profile_source_id",
                self.profile_source_id,
                expected=PAPER_0710_1869_PHYSICAL_MASS_CKM_PROFILE_SOURCE_ID,
            ),
        )
        object.__setattr__(
            self,
            "provenance_mode_id",
            _require_exact_identifier(
                "provenance_mode_id",
                self.provenance_mode_id,
                expected=PAPER_0710_1869_PHYSICAL_MASS_CKM_PROVENANCE_MODE_ID,
            ),
        )
        object.__setattr__(self, "validation_tolerance", tolerance)
        object.__setattr__(
            self,
            "schema_id",
            require_known_schema_id(
                "schema_id",
                self.schema_id,
                expected=PAPER_0710_1869_PHYSICAL_MASS_CKM_RESULT_SCHEMA_ID,
            ),
        )
        if self.claim_level_id != self.bulk_state.claim_level_id:
            raise ValueError("result claim_level_id must match bulk_state.claim_level_id")
        if self.bulk_mass_map_policy_id != self.bulk_state.bulk_mass_map_policy_id:
            raise ValueError(
                "result bulk_mass_map_policy_id must match bulk_state.bulk_mass_map_policy_id"
            )
        if self.profile_input_policy_id != self.bulk_state.profile_input_policy_id:
            raise ValueError(
                "result profile_input_policy_id must match bulk_state.profile_input_policy_id"
            )
        if self.physical_profile_status_id != self.bulk_state.physical_profile_status_id:
            raise ValueError(
                "result physical_profile_status_id must match "
                "bulk_state.physical_profile_status_id"
            )
        if self.core_contract_id != self.bulk_state.core_contract_id:
            raise ValueError("result core_contract_id must match bulk_state.core_contract_id")

        ckm_from_left_rotations = self.U_L_u.conjugate().T @ self.U_L_d
        ckm_consistency_residual = float(np.max(np.abs(self.ckm - ckm_from_left_rotations)))
        if ckm_consistency_residual > tolerance:
            raise ValueError(
                "ckm must match U_L_u^dagger U_L_d within validation_tolerance; "
                f"got residual {ckm_consistency_residual:.3e}"
            )

        computed_ckm_residual = _unitarity_residual(self.ckm)
        computed_rotation_residual = float(
            max(
                _unitarity_residual(self.U_L_u),
                _unitarity_residual(self.U_R_u),
                _unitarity_residual(self.U_L_d),
                _unitarity_residual(self.U_R_d),
            )
        )
        computed_reconstruction_residual = float(
            max(
                _max_svd_reconstruction_residual(
                    self.M_u,
                    self.U_L_u,
                    self.masses_up,
                    self.U_R_u,
                ),
                _max_svd_reconstruction_residual(
                    self.M_d,
                    self.U_L_d,
                    self.masses_down,
                    self.U_R_d,
                ),
            )
        )

        object.__setattr__(self, "max_ckm_unitarity_residual", computed_ckm_residual)
        object.__setattr__(self, "max_rotation_unitarity_residual", computed_rotation_residual)
        object.__setattr__(self, "max_reconstruction_residual", computed_reconstruction_residual)

        if computed_reconstruction_residual > tolerance:
            raise ValueError(
                "SVD reconstruction must hold within validation_tolerance; "
                f"got residual {computed_reconstruction_residual:.3e}"
            )

    @property
    def state(self) -> Paper07101869PhysicalBulkState:
        return self.bulk_state

    @property
    def point(self) -> Paper07101869PhysicalPoint:
        return self.bulk_state.point

    @property
    def score(self) -> float:
        return float(0.0 if self.total_score is None else self.total_score)

    @property
    def mass_residuals_up(self) -> np.ndarray:
        return np.zeros(3) if self.up_log_residuals is None else self.up_log_residuals

    @property
    def mass_residuals_down(self) -> np.ndarray:
        return np.zeros(3) if self.down_log_residuals is None else self.down_log_residuals

    @property
    def ckm_residuals(self) -> np.ndarray:
        if self.ckm_observable_residuals is not None:
            return self.ckm_observable_residuals
        if self.ckm_abs_residuals is None:
            return np.zeros(4)
        return np.array(
            [
                self.ckm_abs_residuals[0, 1],
                self.ckm_abs_residuals[1, 2],
                self.ckm_abs_residuals[0, 2],
                0.0,
            ],
            dtype=float,
        )

    @property
    def ckm_matrix(self) -> np.ndarray:
        return self.ckm

    @property
    def ckm_observables(self) -> np.ndarray:
        return paper_ckm_observables(self.ckm)

    @property
    def residual_norm(self) -> float:
        return float(
            np.linalg.norm(
                np.concatenate(
                    [self.mass_residuals_up, self.mass_residuals_down, self.ckm_residuals]
                )
            )
        )

    def summary(self) -> str:
        lines = [
            "=" * 60,
            f"Paper 0710.1869 QS1 physical mass/CKM probe: {self.point.label}",
            "=" * 60,
            "Uses the hard-gated point-derived QS1 physical seed-to-profile path. "
            "Masses and CKM come from full SVD diagonalization of point-derived mass "
            "matrices. This remains a QS1-only mass/CKM diagnostic surface, not a "
            "coupling, matching, or observable evaluation.",
            f"r = {self.point.r:.4f}",
            f"epsilon = {self.bulk_state.epsilon:.3e}",
            f"c_Q derived = {np.array2string(self.bulk_state.c_Q, precision=4)}",
            f"c_u derived = {np.array2string(self.bulk_state.c_u, precision=4)}",
            f"c_d derived = {np.array2string(self.bulk_state.c_d, precision=4)}",
            f"F_Q derived = {np.array2string(self.bulk_state.F_Q, precision=5)}",
            f"F_u derived = {np.array2string(self.bulk_state.F_u, precision=5)}",
            f"F_d derived = {np.array2string(self.bulk_state.F_d, precision=5)}",
            f"m_u phys = {np.array2string(self.masses_up, precision=5)}",
            f"m_d phys = {np.array2string(self.masses_down, precision=5)}",
            f"CKM obs = {np.array2string(self.ckm_observables, precision=6)}",
            f"profile source = {self.profile_source_id}",
            f"provenance mode = {self.provenance_mode_id}",
        ]
        if self.total_score is not None:
            lines.append(f"probe score = {self.total_score:.6e}")
        return "\n".join(lines)


def build_paper_0710_1869_benchmark_reference_mass_matrices(
    bulk_state: Paper07101869BulkState,
) -> tuple[np.ndarray, np.ndarray]:
    """Build mass matrices from noncanonically attached Table I profile aliases."""
    resolved_bulk_state = _require_structural_reference_bulk_state(bulk_state)
    _, attached_f_q = _rank_order_attach_reference_profiles(
        sector_name="Q",
        quoted_c=resolved_bulk_state.c_Q,
        quoted_f=resolved_bulk_state.F_Q,
    )
    _, attached_f_u = _rank_order_attach_reference_profiles(
        sector_name="u",
        quoted_c=resolved_bulk_state.c_u,
        quoted_f=resolved_bulk_state.F_u,
    )
    _, attached_f_d = _rank_order_attach_reference_profiles(
        sector_name="d",
        quoted_c=resolved_bulk_state.c_d,
        quoted_f=resolved_bulk_state.F_d,
    )
    F_Q = np.diag(attached_f_q)
    F_u = np.diag(attached_f_u)
    F_d = np.diag(attached_f_d)
    prefactor = 2.0 * resolved_bulk_state.point.v
    M_u = prefactor * F_Q @ resolved_bulk_state.Y_u_bulk_basis @ F_u
    M_d = prefactor * F_Q @ resolved_bulk_state.Y_d_bulk_basis @ F_d
    return M_u, M_d


def paper_benchmark_reference_mass_probe_observables(
    M_u: np.ndarray,
    M_d: np.ndarray,
) -> dict[str, np.ndarray]:
    """Extract masses and CKM data from a benchmark-reference mass probe."""
    return _mass_probe_observables(M_u, M_d)


def paper_benchmark_reference_mass_probe_residuals(
    masses_up: np.ndarray,
    masses_down: np.ndarray,
    ckm: np.ndarray,
    targets: Paper07101869BenchmarkReferenceTargets,
) -> dict[str, np.ndarray | float]:
    """Compute residuals against benchmark-reference targets."""
    if not isinstance(targets, Paper07101869BenchmarkReferenceTargets):
        raise ValueError("targets must be a Paper07101869BenchmarkReferenceTargets")
    return _mass_probe_residuals(
        masses_up=masses_up,
        masses_down=masses_down,
        ckm=ckm,
        target_up_masses=targets.up_masses,
        target_down_masses=targets.down_masses,
        target_ckm=targets.ckm,
        ckm_unitarity_tolerance=PAPER_0710_1869_BENCHMARK_REFERENCE_UNITARITY_TOLERANCE,
    )


def build_paper_0710_1869_physical_mass_matrices(
    bulk_state: Paper07101869PhysicalBulkState,
) -> tuple[np.ndarray, np.ndarray]:
    """Build mass matrices from the point-derived QS1 physical bulk state."""
    resolved_bulk_state = _require_physical_bulk_state(bulk_state)
    F_Q = np.diag(resolved_bulk_state.F_Q)
    F_u = np.diag(resolved_bulk_state.F_u)
    F_d = np.diag(resolved_bulk_state.F_d)
    prefactor = 2.0 * resolved_bulk_state.point.v
    M_u = prefactor * F_Q @ resolved_bulk_state.Y_u_bulk_basis @ F_u
    M_d = prefactor * F_Q @ resolved_bulk_state.Y_d_bulk_basis @ F_d
    return M_u, M_d


def paper_physical_mass_ckm_probe_observables(
    M_u: np.ndarray,
    M_d: np.ndarray,
) -> dict[str, np.ndarray]:
    """Extract masses and CKM data from a QS1 physical mass probe."""
    return _mass_probe_observables(M_u, M_d)


def paper_physical_mass_ckm_probe_residuals(
    masses_up: np.ndarray,
    masses_down: np.ndarray,
    ckm: np.ndarray,
    targets: Paper07101869PhysicalMassCKMTargets,
) -> dict[str, np.ndarray | float]:
    """Compute residuals against QS1 physical masses/CKM targets."""
    if not isinstance(targets, Paper07101869PhysicalMassCKMTargets):
        raise ValueError("targets must be a Paper07101869PhysicalMassCKMTargets")
    return _mass_probe_residuals(
        masses_up=masses_up,
        masses_down=masses_down,
        ckm=ckm,
        target_up_masses=targets.up_masses,
        target_down_masses=targets.down_masses,
        target_ckm=targets.ckm,
        ckm_unitarity_tolerance=PAPER_0710_1869_PHYSICAL_MASS_CKM_UNITARITY_TOLERANCE,
    )


def evaluate_paper_0710_1869_physical_mass_ckm_probe(
    point: Paper07101869PhysicalPoint,
    targets: Paper07101869PhysicalMassCKMTargets | None = None,
) -> Paper07101869PhysicalMassCKMProbeResult:
    """Evaluate the QS1 physical masses/CKM probe for one physical point."""
    if not isinstance(point, Paper07101869PhysicalPoint):
        raise ValueError("point must be a Paper07101869PhysicalPoint")
    bulk_state = derive_paper_0710_1869_physical_bulk_state(point)
    _require_physical_bulk_state(bulk_state)
    if targets is not None:
        _require_targets_align_with_physical_bulk_state(targets, bulk_state)

    M_u, M_d = build_paper_0710_1869_physical_mass_matrices(bulk_state)
    observables = paper_physical_mass_ckm_probe_observables(M_u, M_d)
    residuals: dict[str, np.ndarray | float] = {}
    target_label: str | None = None
    target_schema_id: str | None = None
    if targets is not None:
        residuals = paper_physical_mass_ckm_probe_residuals(
            observables["masses_up"],
            observables["masses_down"],
            observables["ckm"],
            targets,
        )
        target_label = targets.label
        target_schema_id = targets.schema_id

    return Paper07101869PhysicalMassCKMProbeResult(
        bulk_state=bulk_state,
        M_u=M_u,
        M_d=M_d,
        U_L_u=observables["U_L_u"],
        U_R_u=observables["U_R_u"],
        U_L_d=observables["U_L_d"],
        U_R_d=observables["U_R_d"],
        masses_up=observables["masses_up"],
        masses_down=observables["masses_down"],
        ckm=observables["ckm"],
        up_log_residuals=residuals.get("up_log_residuals"),  # type: ignore[arg-type]
        down_log_residuals=residuals.get("down_log_residuals"),  # type: ignore[arg-type]
        ckm_abs_residuals=residuals.get("ckm_abs_residuals"),  # type: ignore[arg-type]
        ckm_observable_residuals=residuals.get("ckm_observable_residuals"),  # type: ignore[arg-type]
        total_score=residuals.get("total_score"),  # type: ignore[arg-type]
        target_label=target_label,
        target_schema_id=target_schema_id,
    )


def evaluate_paper_0710_1869_benchmark_reference_mass_probe(
    point: Paper07101869Point,
    targets: Paper07101869BenchmarkReferenceTargets | None = None,
) -> Paper07101869BenchmarkReferenceMassProbeResult:
    """Evaluate a noncanonical benchmark-reference mass probe for one structural point."""
    if not isinstance(point, Paper07101869Point):
        raise ValueError("point must be a Paper07101869Point")
    bulk_state = derive_paper_0710_1869_bulk_state(point)
    _require_structural_reference_bulk_state(bulk_state)
    if targets is not None:
        _require_targets_align_with_bulk_state(targets, bulk_state)

    M_u, M_d = build_paper_0710_1869_benchmark_reference_mass_matrices(bulk_state)
    observables = paper_benchmark_reference_mass_probe_observables(M_u, M_d)
    residuals: dict[str, np.ndarray | float] = {}
    target_label: str | None = None
    target_schema_id: str | None = None
    if targets is not None:
        residuals = paper_benchmark_reference_mass_probe_residuals(
            observables["masses_up"],
            observables["masses_down"],
            observables["ckm"],
            targets,
        )
        target_label = targets.label
        target_schema_id = targets.schema_id

    return Paper07101869BenchmarkReferenceMassProbeResult(
        bulk_state=bulk_state,
        M_u=M_u,
        M_d=M_d,
        U_L_u=observables["U_L_u"],
        U_R_u=observables["U_R_u"],
        U_L_d=observables["U_L_d"],
        U_R_d=observables["U_R_d"],
        masses_up=observables["masses_up"],
        masses_down=observables["masses_down"],
        ckm=observables["ckm"],
        up_log_residuals=residuals.get("up_log_residuals"),  # type: ignore[arg-type]
        down_log_residuals=residuals.get("down_log_residuals"),  # type: ignore[arg-type]
        ckm_abs_residuals=residuals.get("ckm_abs_residuals"),  # type: ignore[arg-type]
        ckm_observable_residuals=residuals.get("ckm_observable_residuals"),  # type: ignore[arg-type]
        total_score=residuals.get("total_score"),  # type: ignore[arg-type]
        target_label=target_label,
        target_schema_id=target_schema_id,
    )


# Compatibility shims for the current package export surface.
Paper07101869Targets = Paper07101869BenchmarkReferenceTargets
Paper07101869FitResult = Paper07101869BenchmarkReferenceMassProbeResult


def build_paper_0710_1869_mass_matrices(
    bulk_state: Paper07101869BulkState,
) -> tuple[np.ndarray, np.ndarray]:
    """Compatibility shim for benchmark-reference mass matrices."""
    return build_paper_0710_1869_benchmark_reference_mass_matrices(bulk_state)


def paper_mass_matrix_observables(M_u: np.ndarray, M_d: np.ndarray) -> dict[str, np.ndarray]:
    """Compatibility shim for benchmark-reference mass-probe observables."""
    return paper_benchmark_reference_mass_probe_observables(M_u, M_d)


def paper_fit_residuals(
    masses_up: np.ndarray,
    masses_down: np.ndarray,
    ckm: np.ndarray,
    targets: Paper07101869BenchmarkReferenceTargets,
) -> dict[str, np.ndarray | float]:
    """Compatibility shim for benchmark-reference mass-probe residuals."""
    return paper_benchmark_reference_mass_probe_residuals(masses_up, masses_down, ckm, targets)


def evaluate_paper_0710_1869_fit(
    point: Paper07101869Point,
    targets: Paper07101869BenchmarkReferenceTargets | None = None,
) -> Paper07101869BenchmarkReferenceMassProbeResult:
    """Compatibility shim for the benchmark-reference mass probe."""
    return evaluate_paper_0710_1869_benchmark_reference_mass_probe(point, targets)


__all__ = [
    "PAPER_0710_1869_BENCHMARK_REFERENCE_ATTACHMENT_CONVENTION_ID",
    "PAPER_0710_1869_BENCHMARK_REFERENCE_CONTRACT_ID",
    "PAPER_0710_1869_BENCHMARK_REFERENCE_PROFILE_SOURCE_ID",
    "PAPER_0710_1869_BENCHMARK_REFERENCE_PROVENANCE_MODE_ID",
    "PAPER_0710_1869_BENCHMARK_REFERENCE_RECONSTRUCTION_TOLERANCE",
    "PAPER_0710_1869_BENCHMARK_REFERENCE_RESULT_KIND_ID",
    "PAPER_0710_1869_BENCHMARK_REFERENCE_RESULT_SCHEMA_ID",
    "PAPER_0710_1869_BENCHMARK_REFERENCE_TARGET_KIND_ID",
    "PAPER_0710_1869_BENCHMARK_REFERENCE_TARGETS_SCHEMA_ID",
    "PAPER_0710_1869_BENCHMARK_REFERENCE_UNITARITY_TOLERANCE",
    "PAPER_0710_1869_FIT_RESULT_SCHEMA_ID",
    "PAPER_0710_1869_PHYSICAL_MASS_CKM_CONTRACT_ID",
    "PAPER_0710_1869_PHYSICAL_MASS_CKM_PROFILE_SOURCE_ID",
    "PAPER_0710_1869_PHYSICAL_MASS_CKM_PROVENANCE_MODE_ID",
    "PAPER_0710_1869_PHYSICAL_MASS_CKM_RECONSTRUCTION_TOLERANCE",
    "PAPER_0710_1869_PHYSICAL_MASS_CKM_RESULT_KIND_ID",
    "PAPER_0710_1869_PHYSICAL_MASS_CKM_RESULT_SCHEMA_ID",
    "PAPER_0710_1869_PHYSICAL_MASS_CKM_TARGET_KIND_ID",
    "PAPER_0710_1869_PHYSICAL_MASS_CKM_TARGETS_SCHEMA_ID",
    "PAPER_0710_1869_PHYSICAL_MASS_CKM_UNITARITY_TOLERANCE",
    "PAPER_0710_1869_TARGETS_SCHEMA_ID",
    "Paper07101869BenchmarkReferenceMassProbeResult",
    "Paper07101869BenchmarkReferenceTargets",
    "Paper07101869FitResult",
    "Paper07101869PhysicalMassCKMProbeResult",
    "Paper07101869PhysicalMassCKMTargets",
    "Paper07101869Targets",
    "build_paper_0710_1869_benchmark_reference_mass_matrices",
    "build_paper_0710_1869_mass_matrices",
    "build_paper_0710_1869_physical_mass_matrices",
    "evaluate_paper_0710_1869_benchmark_reference_mass_probe",
    "evaluate_paper_0710_1869_fit",
    "evaluate_paper_0710_1869_physical_mass_ckm_probe",
    "paper_benchmark_reference_mass_probe_observables",
    "paper_benchmark_reference_mass_probe_residuals",
    "paper_ckm_observables",
    "paper_fit_residuals",
    "paper_physical_mass_ckm_probe_observables",
    "paper_physical_mass_ckm_probe_residuals",
    "paper_jarlskog_invariant",
    "paper_mass_matrix_observables",
]
