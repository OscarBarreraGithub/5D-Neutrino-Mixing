"""First paper-owned KK-gluon coupling layer for ``paper_0710_1869``."""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass

import numpy as np

from .benchmarks import (
    Paper07101869Benchmark,
    default_paper_0710_1869_pr1_benchmark,
)
from .conventions import PAPER_0710_1869_MODE_ID, PAPER_0710_1869_PAPER_ID
from .couplings import (
    PAPER_0710_1869_DIMENSIONLESS_MATRIX_POLICY_ID,
    PAPER_0710_1869_EFFECTIVE_SCALE_POLICY_ID,
    PAPER_0710_1869_PROPAGATOR_MASS_RULE_ID,
    PAPER_0710_1869_SINGLE_MODE_SCALE_POLICY_ID,
    PAPER_0710_1869_UNIVERSAL_SUBTRACTION_POLICY_ID,
    Paper07101869CouplingContract,
    Paper07101869GaugeCouplingNormalization,
    default_paper_0710_1869_coupling_contract,
    evaluate_paper_0710_1869_gauge_coupling,
)
from .fit import (
    build_paper_0710_1869_physical_mass_matrices,
    paper_physical_mass_ckm_probe_observables,
)
from .model import (
    PAPER_0710_1869_POINT_DERIVED_PHYSICAL_STATUS_ID,
    Paper07101869PhysicalBulkState,
    Paper07101869V5KMParameters,
    build_table_i_benchmark_from_inputs,
    paper_v5km_matrix,
)
from .scales import Paper07101869ScalePoint, default_paper_0710_1869_scales
from .validation import (
    require_known_schema_id,
    require_member,
    require_nonempty_identifier,
    require_positive_finite,
)

PAPER_0710_1869_KK_GLUON_MATRIX_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.kk_gluon_matrix.v1"
)
PAPER_0710_1869_KK_GLUON_SCHEMA_ID = "quarkConstraints.paper_0710_1869.kk_gluon.v1"
PAPER_0710_1869_KK_GLUON_SUMMARY_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.kk_gluon_summary.v1"
)
PAPER_0710_1869_LEFT_UP_ALIGNED_BASIS_ID = "left_handed.up_aligned_identity.v1"
PAPER_0710_1869_LEFT_DOWN_ALIGNED_BASIS_ID = (
    "left_handed.down_aligned_using_v5km_example.v1"
)
PAPER_0710_1869_PHYSICAL_LEFT_UP_ALIGNED_BASIS_ID = (
    "left_handed.up_aligned_using_qs1_mass_probe_u_l_u.v1"
)
PAPER_0710_1869_PHYSICAL_LEFT_DOWN_ALIGNED_BASIS_ID = (
    "left_handed.down_aligned_using_qs1_mass_probe_u_l_d.v1"
)
PAPER_0710_1869_RIGHT_DIAGONAL_BASIS_ID = "right_handed.table_i_diagonal_basis.v1"
PAPER_0710_1869_PHYSICAL_RIGHT_UP_ALIGNED_BASIS_ID = (
    "right_handed.up_aligned_using_qs1_mass_probe_u_r_u.v1"
)
PAPER_0710_1869_PHYSICAL_RIGHT_DOWN_ALIGNED_BASIS_ID = (
    "right_handed.down_aligned_using_qs1_mass_probe_u_r_d.v1"
)
PAPER_0710_1869_LEFT_HANDED_BASIS_POLICY_ID = (
    "left_handed.paired_up_and_down_aligned_outputs.v1"
)
PAPER_0710_1869_PHYSICAL_LEFT_HANDED_BASIS_POLICY_ID = (
    "left_handed.paired_up_and_down_aligned_outputs.from_qs1_mass_probe.v1"
)
PAPER_0710_1869_RIGHT_HANDED_BASIS_POLICY_ID = (
    "right_handed.diagonal_only_until_paper_owned_rotations_exist.v1"
)
PAPER_0710_1869_PHYSICAL_RIGHT_HANDED_BASIS_POLICY_ID = (
    "right_handed.paired_up_and_down_aligned_outputs.from_qs1_mass_probe.v1"
)
PAPER_0710_1869_PHYSICAL_KK_GLUON_STATUS_ID = PAPER_0710_1869_POINT_DERIVED_PHYSICAL_STATUS_ID


def _as_complex_matrix3(name: str, values: np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=np.complex128)
    if arr.shape != (3, 3):
        raise ValueError(f"{name} must have shape (3, 3)")
    if not np.all(np.isfinite(arr.real)) or not np.all(np.isfinite(arr.imag)):
        raise ValueError(f"{name} must contain only finite entries")
    return arr


def _as_real_triplet(
    name: str, values: tuple[float, float, float] | list[float] | np.ndarray
) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.shape != (3,):
        raise ValueError(f"{name} must have shape (3,)")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must contain only finite values")
    return arr.astype(float, copy=True)


def _hermitian(matrix: np.ndarray) -> np.ndarray:
    arr = _as_complex_matrix3("matrix", matrix)
    return 0.5 * (arr + arr.conjugate().T)


def _diagonal_profile_kernel(
    name: str, profile_values: tuple[float, float, float] | list[float] | np.ndarray
) -> np.ndarray:
    vector = _as_real_triplet(name, profile_values)
    if np.any(vector <= 0.0):
        raise ValueError(f"{name} must contain positive profile eigenvalues")
    return np.diag(vector**2).astype(np.complex128)


def _universal_component(matrix: np.ndarray) -> float:
    arr = _as_complex_matrix3("matrix", matrix)
    trace = np.trace(arr)
    if not math.isclose(float(trace.imag), 0.0, abs_tol=1.0e-12):
        raise ValueError("trace of a flavor matrix must be real")
    return float(trace.real / 3.0)


def _subtract_universal_piece(matrix: np.ndarray) -> tuple[np.ndarray, float]:
    arr = _as_complex_matrix3("matrix", matrix)
    universal_component = _universal_component(arr)
    shifted = _hermitian(arr - (universal_component * np.eye(3, dtype=np.complex128)))
    return shifted, universal_component


def _offdiag_fro_norm(matrix: np.ndarray) -> float:
    arr = _as_complex_matrix3("matrix", matrix).copy()
    np.fill_diagonal(arr, 0.0)
    return float(np.linalg.norm(arr, ord="fro"))


def _max_abs_imag(matrix: np.ndarray) -> float:
    arr = _as_complex_matrix3("matrix", matrix)
    return float(np.max(np.abs(arr.imag)))


def _real_diag(matrix: np.ndarray) -> tuple[float, float, float]:
    arr = _as_complex_matrix3("matrix", matrix)
    diagonal = np.diag(arr)
    if not np.allclose(diagonal.imag, 0.0, atol=1.0e-12):
        raise ValueError("diagonal entries must be real")
    return tuple(float(value) for value in diagonal.real)


def _matrix_to_rows(matrix: np.ndarray) -> list[list[list[float]]]:
    arr = _as_complex_matrix3("matrix", matrix)
    return [
        [[float(value.real), float(value.imag)] for value in row]
        for row in arr
    ]


@dataclass(frozen=True)
class Paper07101869KKGluonFlavorMatrix:
    """One benchmark-facing KK-gluon flavor matrix block."""

    label: str
    sector_id: str
    chiral_id: str
    basis_id: str
    g_s_mu_gs: float
    rs_volume_L: float
    rs_volume_sqrt_2L: float
    universal_component: float
    raw_dimensionless: np.ndarray
    universal_subtracted_dimensionless: np.ndarray
    raw_gs_normalized: np.ndarray
    universal_subtracted_gs_normalized: np.ndarray
    schema_id: str = PAPER_0710_1869_KK_GLUON_MATRIX_SCHEMA_ID

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "schema_id",
            require_known_schema_id(
                "schema_id",
                self.schema_id,
                expected=PAPER_0710_1869_KK_GLUON_MATRIX_SCHEMA_ID,
            ),
        )
        object.__setattr__(self, "label", require_nonempty_identifier("label", self.label))
        object.__setattr__(self, "basis_id", require_nonempty_identifier("basis_id", self.basis_id))
        require_member("sector_id", self.sector_id, ("Q", "u", "d"))
        require_member("chiral_id", self.chiral_id, ("left", "right"))
        object.__setattr__(self, "g_s_mu_gs", require_positive_finite("g_s_mu_gs", self.g_s_mu_gs))
        object.__setattr__(
            self,
            "rs_volume_L",
            require_positive_finite("rs_volume_L", self.rs_volume_L),
        )
        object.__setattr__(
            self,
            "rs_volume_sqrt_2L",
            require_positive_finite("rs_volume_sqrt_2L", self.rs_volume_sqrt_2L),
        )
        if not math.isclose(
            self.rs_volume_sqrt_2L,
            math.sqrt(2.0 * self.rs_volume_L),
            rel_tol=1.0e-12,
            abs_tol=0.0,
        ):
            raise ValueError("rs_volume_sqrt_2L must equal sqrt(2*rs_volume_L)")
        if not math.isfinite(float(self.universal_component)):
            raise ValueError("universal_component must be finite")
        object.__setattr__(self, "raw_dimensionless", _hermitian(self.raw_dimensionless))
        object.__setattr__(
            self,
            "universal_subtracted_dimensionless",
            _hermitian(self.universal_subtracted_dimensionless),
        )
        object.__setattr__(self, "raw_gs_normalized", _hermitian(self.raw_gs_normalized))
        object.__setattr__(
            self,
            "universal_subtracted_gs_normalized",
            _hermitian(self.universal_subtracted_gs_normalized),
        )

        expected_subtracted = self.raw_dimensionless - (
            self.universal_component * np.eye(3, dtype=np.complex128)
        )
        if not np.allclose(
            self.universal_subtracted_dimensionless,
            expected_subtracted,
            atol=1.0e-12,
        ):
            raise ValueError(
                "universal_subtracted_dimensionless must subtract the universal identity piece"
            )
        expected_raw_gs = self.g_s_mu_gs * (
            self.rs_volume_sqrt_2L * self.raw_dimensionless
            - (np.eye(3, dtype=np.complex128) / self.rs_volume_sqrt_2L)
        )
        if not np.allclose(self.raw_gs_normalized, expected_raw_gs, atol=1.0e-12):
            raise ValueError(
                "raw_gs_normalized must equal "
                "g_s_mu_gs * (sqrt(2L) * raw_dimensionless - I / sqrt(2L))"
            )
        expected_subtracted_gs = (
            self.g_s_mu_gs
            * self.rs_volume_sqrt_2L
            * self.universal_subtracted_dimensionless
        )
        if not np.allclose(
            self.universal_subtracted_gs_normalized,
            expected_subtracted_gs,
            atol=1.0e-12,
        ):
            raise ValueError(
                "universal_subtracted_gs_normalized must equal "
                "g_s_mu_gs * sqrt(2L) * universal_subtracted_dimensionless"
            )
        if not math.isclose(
            float(np.trace(self.universal_subtracted_dimensionless).real),
            0.0,
            abs_tol=1.0e-12,
        ):
            raise ValueError("universal_subtracted_dimensionless must be traceless")

    def summary(self) -> "Paper07101869KKGluonMatrixSummary":
        """Return a deterministic scalar summary for this matrix block."""
        return Paper07101869KKGluonMatrixSummary(
            label=self.label,
            sector_id=self.sector_id,
            chiral_id=self.chiral_id,
            basis_id=self.basis_id,
            universal_component=float(self.universal_component),
            raw_diag=_real_diag(self.raw_dimensionless),
            universal_subtracted_diag=_real_diag(self.universal_subtracted_dimensionless),
            raw_trace=float(np.trace(self.raw_dimensionless).real),
            raw_offdiag_fro_norm=_offdiag_fro_norm(self.raw_dimensionless),
            universal_subtracted_offdiag_fro_norm=_offdiag_fro_norm(
                self.universal_subtracted_dimensionless
            ),
            raw_gs_normalized_offdiag_fro_norm=_offdiag_fro_norm(self.raw_gs_normalized),
            universal_subtracted_gs_normalized_offdiag_fro_norm=_offdiag_fro_norm(
                self.universal_subtracted_gs_normalized
            ),
            max_abs_imag_raw=_max_abs_imag(self.raw_dimensionless),
        )

    def as_dict(self) -> dict[str, object]:
        """Return a deterministic representation including full matrices."""
        return {
            "schema_id": self.schema_id,
            "label": self.label,
            "sector_id": self.sector_id,
            "chiral_id": self.chiral_id,
            "basis_id": self.basis_id,
            "g_s_mu_gs": self.g_s_mu_gs,
            "rs_volume_L": self.rs_volume_L,
            "rs_volume_sqrt_2L": self.rs_volume_sqrt_2L,
            "universal_component": self.universal_component,
            "raw_dimensionless": _matrix_to_rows(self.raw_dimensionless),
            "universal_subtracted_dimensionless": _matrix_to_rows(
                self.universal_subtracted_dimensionless
            ),
            "raw_gs_normalized": _matrix_to_rows(self.raw_gs_normalized),
            "universal_subtracted_gs_normalized": _matrix_to_rows(
                self.universal_subtracted_gs_normalized
            ),
            "summary": self.summary().as_dict(),
        }


@dataclass(frozen=True)
class Paper07101869KKGluonMatrixSummary:
    """Scalar summary for one KK-gluon matrix block."""

    label: str
    sector_id: str
    chiral_id: str
    basis_id: str
    universal_component: float
    raw_diag: tuple[float, float, float]
    universal_subtracted_diag: tuple[float, float, float]
    raw_trace: float
    raw_offdiag_fro_norm: float
    universal_subtracted_offdiag_fro_norm: float
    raw_gs_normalized_offdiag_fro_norm: float
    universal_subtracted_gs_normalized_offdiag_fro_norm: float
    max_abs_imag_raw: float

    def __post_init__(self) -> None:
        object.__setattr__(self, "label", require_nonempty_identifier("label", self.label))
        object.__setattr__(self, "basis_id", require_nonempty_identifier("basis_id", self.basis_id))
        require_member("sector_id", self.sector_id, ("Q", "u", "d"))
        require_member("chiral_id", self.chiral_id, ("left", "right"))
        for field_name in (
            "universal_component",
            "raw_trace",
            "raw_offdiag_fro_norm",
            "universal_subtracted_offdiag_fro_norm",
            "raw_gs_normalized_offdiag_fro_norm",
            "universal_subtracted_gs_normalized_offdiag_fro_norm",
            "max_abs_imag_raw",
        ):
            numeric = float(getattr(self, field_name))
            if not math.isfinite(numeric):
                raise ValueError(f"{field_name} must be finite")
            object.__setattr__(self, field_name, numeric)
        object.__setattr__(self, "raw_diag", _real_triplet_tuple("raw_diag", self.raw_diag))
        object.__setattr__(
            self,
            "universal_subtracted_diag",
            _real_triplet_tuple("universal_subtracted_diag", self.universal_subtracted_diag),
        )

    def as_dict(self) -> dict[str, object]:
        return {
            "label": self.label,
            "sector_id": self.sector_id,
            "chiral_id": self.chiral_id,
            "basis_id": self.basis_id,
            "universal_component": self.universal_component,
            "raw_diag": list(self.raw_diag),
            "universal_subtracted_diag": list(self.universal_subtracted_diag),
            "raw_trace": self.raw_trace,
            "raw_offdiag_fro_norm": self.raw_offdiag_fro_norm,
            "universal_subtracted_offdiag_fro_norm": self.universal_subtracted_offdiag_fro_norm,
            "raw_gs_normalized_offdiag_fro_norm": self.raw_gs_normalized_offdiag_fro_norm,
            "universal_subtracted_gs_normalized_offdiag_fro_norm": (
                self.universal_subtracted_gs_normalized_offdiag_fro_norm
            ),
            "max_abs_imag_raw": self.max_abs_imag_raw,
        }


@dataclass(frozen=True)
class Paper07101869KKGluonCouplings:
    """Deterministic KK-gluon benchmark couplings from frozen paper inputs."""

    contract: Paper07101869CouplingContract
    normalization: Paper07101869GaugeCouplingNormalization
    left_up_aligned: Paper07101869KKGluonFlavorMatrix
    left_down_aligned: Paper07101869KKGluonFlavorMatrix
    right_up: Paper07101869KKGluonFlavorMatrix
    right_down: Paper07101869KKGluonFlavorMatrix
    benchmark_id: str
    benchmark_label: str
    benchmark_status: str
    theta12_deg: float
    theta23_deg: float
    theta13_deg: float
    delta: float
    schema_id: str = PAPER_0710_1869_KK_GLUON_SCHEMA_ID
    mode_id: str = PAPER_0710_1869_MODE_ID
    paper_id: str = PAPER_0710_1869_PAPER_ID
    left_handed_basis_policy_id: str = PAPER_0710_1869_LEFT_HANDED_BASIS_POLICY_ID
    right_handed_basis_policy_id: str = PAPER_0710_1869_RIGHT_HANDED_BASIS_POLICY_ID
    notes: str = (
        "Q-sector couplings are exposed in paired left-handed benchmark bases. "
        "Right-handed sectors remain diagonal until paper-owned rotations are frozen."
    )

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "schema_id",
            require_known_schema_id(
                "schema_id",
                self.schema_id,
                expected=PAPER_0710_1869_KK_GLUON_SCHEMA_ID,
            ),
        )
        object.__setattr__(self, "mode_id", require_nonempty_identifier("mode_id", self.mode_id))
        object.__setattr__(self, "paper_id", require_nonempty_identifier("paper_id", self.paper_id))
        require_member("mode_id", self.mode_id, (PAPER_0710_1869_MODE_ID,))
        require_member("paper_id", self.paper_id, (PAPER_0710_1869_PAPER_ID,))
        if not isinstance(self.contract, Paper07101869CouplingContract):
            raise ValueError("contract must be a Paper07101869CouplingContract")
        if not isinstance(self.normalization, Paper07101869GaugeCouplingNormalization):
            raise ValueError(
                "normalization must be a Paper07101869GaugeCouplingNormalization"
            )
        for field_name in (
            "benchmark_id",
            "benchmark_label",
            "benchmark_status",
            "left_handed_basis_policy_id",
            "right_handed_basis_policy_id",
            "notes",
        ):
            object.__setattr__(
                self,
                field_name,
                require_nonempty_identifier(field_name, getattr(self, field_name)),
            )
        require_member(
            "benchmark_status",
            self.benchmark_status,
            _allowed_kk_gluon_status_ids(),
        )
        for field_name in ("theta12_deg", "theta23_deg", "theta13_deg", "delta"):
            numeric = float(getattr(self, field_name))
            if not math.isfinite(numeric):
                raise ValueError(f"{field_name} must be finite")
            object.__setattr__(self, field_name, numeric)
        for field_name in (
            "left_up_aligned",
            "left_down_aligned",
            "right_up",
            "right_down",
        ):
            if not isinstance(getattr(self, field_name), Paper07101869KKGluonFlavorMatrix):
                raise ValueError(f"{field_name} must be a Paper07101869KKGluonFlavorMatrix")

    def as_dict(self) -> dict[str, object]:
        """Return a deterministic representation of the KK-gluon couplings object."""
        return {
            "schema_id": self.schema_id,
            "mode_id": self.mode_id,
            "paper_id": self.paper_id,
            "benchmark_id": self.benchmark_id,
            "benchmark_label": self.benchmark_label,
            "benchmark_status": self.benchmark_status,
            "theta12_deg": self.theta12_deg,
            "theta23_deg": self.theta23_deg,
            "theta13_deg": self.theta13_deg,
            "delta": self.delta,
            "left_handed_basis_policy_id": self.left_handed_basis_policy_id,
            "right_handed_basis_policy_id": self.right_handed_basis_policy_id,
            "contract": {
                "schema_id": self.contract.schema_id,
                "dimensionless_matrix_policy_id": self.contract.dimensionless_matrix_policy_id,
                "mu_gs_semantics_id": self.contract.mu_gs_semantics_id,
                "universal_subtraction_policy_id": self.contract.universal_subtraction_policy_id,
                "propagator_mass_rule_id": self.contract.propagator_mass_rule_id,
                "single_mode_scale_policy_id": self.contract.single_mode_scale_policy_id,
                "effective_scale_policy_id": self.contract.effective_scale_policy_id,
            },
            "normalization": {
                "schema_id": self.normalization.schema_id,
                "contract_schema_id": self.normalization.contract_schema_id,
                "scale_label": self.normalization.scale_label,
                "mu_gs_GeV": self.normalization.mu_gs_GeV,
                "mu_match_GeV": self.normalization.mu_match_GeV,
                "m_g1_GeV": self.normalization.m_g1_GeV,
                "propagator_mass_GeV": self.normalization.propagator_mass_GeV,
                "alpha_s_mu_gs": self.normalization.alpha_s_mu_gs,
                "g_s_mu_gs": self.normalization.g_s_mu_gs,
                "kk_gluon_normalization_id": self.normalization.kk_gluon_normalization_id,
                "mu_gs_semantics_id": self.normalization.mu_gs_semantics_id,
                "scale_policy_id": self.normalization.scale_policy_id,
                "propagator_mass_rule_id": self.normalization.propagator_mass_rule_id,
            },
            "left_up_aligned": self.left_up_aligned.as_dict(),
            "left_down_aligned": self.left_down_aligned.as_dict(),
            "right_up": self.right_up.as_dict(),
            "right_down": self.right_down.as_dict(),
            "notes": self.notes,
            "summary": self.summary().as_dict(),
        }

    def summary(self) -> "Paper07101869KKGluonBenchmarkSummary":
        """Return the deterministic benchmark summary."""
        return Paper07101869KKGluonBenchmarkSummary(
            benchmark_id=self.benchmark_id,
            benchmark_label=self.benchmark_label,
            benchmark_status=self.benchmark_status,
            scale_label=self.normalization.scale_label,
            scale_policy_id=self.normalization.scale_policy_id,
            kk_gluon_normalization_id=self.normalization.kk_gluon_normalization_id,
            dimensionless_matrix_policy_id=self.contract.dimensionless_matrix_policy_id,
            mu_gs_semantics_id=self.contract.mu_gs_semantics_id,
            universal_subtraction_policy_id=self.contract.universal_subtraction_policy_id,
            propagator_mass_rule_id=self.contract.propagator_mass_rule_id,
            left_handed_basis_policy_id=self.left_handed_basis_policy_id,
            right_handed_basis_policy_id=self.right_handed_basis_policy_id,
            mu_gs_GeV=self.normalization.mu_gs_GeV,
            mu_match_GeV=self.normalization.mu_match_GeV,
            m_g1_GeV=self.normalization.m_g1_GeV,
            propagator_mass_GeV=self.normalization.propagator_mass_GeV,
            alpha_s_mu_gs=self.normalization.alpha_s_mu_gs,
            g_s_mu_gs=self.normalization.g_s_mu_gs,
            theta12_deg=self.theta12_deg,
            theta23_deg=self.theta23_deg,
            theta13_deg=self.theta13_deg,
            delta=self.delta,
            matrix_summaries=(
                self.left_up_aligned.summary(),
                self.left_down_aligned.summary(),
                self.right_up.summary(),
                self.right_down.summary(),
            ),
            notes=self.notes,
        )


@dataclass(frozen=True)
class Paper07101869PhysicalKKGluonCouplings(Paper07101869KKGluonCouplings):
    """Physical-side KK-gluon wrapper built from a point-derived bulk state."""

    physical_profile_status_id: str = PAPER_0710_1869_PHYSICAL_KK_GLUON_STATUS_ID
    notes: str = (
        "Physical-side KK-gluon adapter sourced from the point-derived QS1 bulk state. "
        "Benchmark-facing schema compatibility is retained, but the left-handed Q-sector "
        "bases are the physical QS1 SVD outputs U_L_u and U_L_d. Right-handed u/d blocks "
        "are likewise built in the physical QS1 SVD bases U_R_u and U_R_d."
    )

    def __post_init__(self) -> None:
        super().__post_init__()
        object.__setattr__(
            self,
            "physical_profile_status_id",
            require_member(
                "physical_profile_status_id",
                self.physical_profile_status_id,
                (PAPER_0710_1869_PHYSICAL_KK_GLUON_STATUS_ID,),
            ),
        )
        require_member(
            "benchmark_status",
            self.benchmark_status,
            (PAPER_0710_1869_PHYSICAL_KK_GLUON_STATUS_ID,),
        )
        require_member(
            "left_handed_basis_policy_id",
            self.left_handed_basis_policy_id,
            (PAPER_0710_1869_PHYSICAL_LEFT_HANDED_BASIS_POLICY_ID,),
        )
        require_member(
            "right_handed_basis_policy_id",
            self.right_handed_basis_policy_id,
            (PAPER_0710_1869_PHYSICAL_RIGHT_HANDED_BASIS_POLICY_ID,),
        )

    def summary(self) -> "Paper07101869PhysicalKKGluonBenchmarkSummary":
        """Return the deterministic summary for the physical adapter wrapper."""
        return Paper07101869PhysicalKKGluonBenchmarkSummary(
            benchmark_id=self.benchmark_id,
            benchmark_label=self.benchmark_label,
            benchmark_status=self.benchmark_status,
            scale_label=self.normalization.scale_label,
            scale_policy_id=self.normalization.scale_policy_id,
            kk_gluon_normalization_id=self.normalization.kk_gluon_normalization_id,
            dimensionless_matrix_policy_id=self.contract.dimensionless_matrix_policy_id,
            mu_gs_semantics_id=self.contract.mu_gs_semantics_id,
            universal_subtraction_policy_id=self.contract.universal_subtraction_policy_id,
            propagator_mass_rule_id=self.contract.propagator_mass_rule_id,
            left_handed_basis_policy_id=self.left_handed_basis_policy_id,
            right_handed_basis_policy_id=self.right_handed_basis_policy_id,
            mu_gs_GeV=self.normalization.mu_gs_GeV,
            mu_match_GeV=self.normalization.mu_match_GeV,
            m_g1_GeV=self.normalization.m_g1_GeV,
            propagator_mass_GeV=self.normalization.propagator_mass_GeV,
            alpha_s_mu_gs=self.normalization.alpha_s_mu_gs,
            g_s_mu_gs=self.normalization.g_s_mu_gs,
            theta12_deg=self.theta12_deg,
            theta23_deg=self.theta23_deg,
            theta13_deg=self.theta13_deg,
            delta=self.delta,
            matrix_summaries=(
                self.left_up_aligned.summary(),
                self.left_down_aligned.summary(),
                self.right_up.summary(),
                self.right_down.summary(),
            ),
            notes=self.notes,
            physical_profile_status_id=self.physical_profile_status_id,
        )

    def as_dict(self) -> dict[str, object]:
        payload = super().as_dict()
        payload["physical_profile_status_id"] = self.physical_profile_status_id
        return payload


@dataclass(frozen=True)
class Paper07101869KKGluonBenchmarkSummary:
    """Stable summary payload for the benchmark-facing KK-gluon layer."""

    benchmark_id: str
    benchmark_label: str
    benchmark_status: str
    scale_label: str
    scale_policy_id: str
    kk_gluon_normalization_id: str
    dimensionless_matrix_policy_id: str
    mu_gs_semantics_id: str
    universal_subtraction_policy_id: str
    propagator_mass_rule_id: str
    left_handed_basis_policy_id: str
    right_handed_basis_policy_id: str
    mu_gs_GeV: float
    mu_match_GeV: float
    m_g1_GeV: float
    propagator_mass_GeV: float
    alpha_s_mu_gs: float
    g_s_mu_gs: float
    theta12_deg: float
    theta23_deg: float
    theta13_deg: float
    delta: float
    matrix_summaries: tuple[Paper07101869KKGluonMatrixSummary, ...]
    notes: str
    schema_id: str = PAPER_0710_1869_KK_GLUON_SUMMARY_SCHEMA_ID
    mode_id: str = PAPER_0710_1869_MODE_ID
    paper_id: str = PAPER_0710_1869_PAPER_ID

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "schema_id",
            require_known_schema_id(
                "schema_id",
                self.schema_id,
                expected=PAPER_0710_1869_KK_GLUON_SUMMARY_SCHEMA_ID,
            ),
        )
        object.__setattr__(self, "mode_id", require_nonempty_identifier("mode_id", self.mode_id))
        object.__setattr__(self, "paper_id", require_nonempty_identifier("paper_id", self.paper_id))
        require_member("mode_id", self.mode_id, (PAPER_0710_1869_MODE_ID,))
        require_member("paper_id", self.paper_id, (PAPER_0710_1869_PAPER_ID,))
        require_member(
            "benchmark_status",
            self.benchmark_status,
            _allowed_kk_gluon_status_ids(),
        )
        require_member(
            "scale_policy_id",
            self.scale_policy_id,
            (
                PAPER_0710_1869_SINGLE_MODE_SCALE_POLICY_ID,
                PAPER_0710_1869_EFFECTIVE_SCALE_POLICY_ID,
            ),
        )
        require_member(
            "dimensionless_matrix_policy_id",
            self.dimensionless_matrix_policy_id,
            (PAPER_0710_1869_DIMENSIONLESS_MATRIX_POLICY_ID,),
        )
        require_member(
            "universal_subtraction_policy_id",
            self.universal_subtraction_policy_id,
            (PAPER_0710_1869_UNIVERSAL_SUBTRACTION_POLICY_ID,),
        )
        require_member(
            "propagator_mass_rule_id",
            self.propagator_mass_rule_id,
            (PAPER_0710_1869_PROPAGATOR_MASS_RULE_ID,),
        )
        for field_name in (
            "benchmark_id",
            "benchmark_label",
            "scale_label",
            "kk_gluon_normalization_id",
            "mu_gs_semantics_id",
            "left_handed_basis_policy_id",
            "right_handed_basis_policy_id",
            "notes",
        ):
            object.__setattr__(
                self,
                field_name,
                require_nonempty_identifier(field_name, getattr(self, field_name)),
            )
        for field_name in (
            "mu_gs_GeV",
            "mu_match_GeV",
            "m_g1_GeV",
            "propagator_mass_GeV",
            "alpha_s_mu_gs",
            "g_s_mu_gs",
            "theta12_deg",
            "theta23_deg",
            "theta13_deg",
            "delta",
        ):
            numeric = float(getattr(self, field_name))
            if field_name not in {"theta12_deg", "theta23_deg", "theta13_deg", "delta"}:
                numeric = require_positive_finite(field_name, numeric)
            elif not math.isfinite(numeric):
                raise ValueError(f"{field_name} must be finite")
            object.__setattr__(self, field_name, numeric)
        if len(self.matrix_summaries) != 4:
            raise ValueError("matrix_summaries must contain exactly four matrix summaries")
        if not all(
            isinstance(item, Paper07101869KKGluonMatrixSummary)
            for item in self.matrix_summaries
        ):
            raise ValueError(
                "matrix_summaries must contain only Paper07101869KKGluonMatrixSummary items"
            )

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_id": self.schema_id,
            "mode_id": self.mode_id,
            "paper_id": self.paper_id,
            "benchmark_id": self.benchmark_id,
            "benchmark_label": self.benchmark_label,
            "benchmark_status": self.benchmark_status,
            "scale_label": self.scale_label,
            "scale_policy_id": self.scale_policy_id,
            "kk_gluon_normalization_id": self.kk_gluon_normalization_id,
            "dimensionless_matrix_policy_id": self.dimensionless_matrix_policy_id,
            "mu_gs_semantics_id": self.mu_gs_semantics_id,
            "universal_subtraction_policy_id": self.universal_subtraction_policy_id,
            "propagator_mass_rule_id": self.propagator_mass_rule_id,
            "left_handed_basis_policy_id": self.left_handed_basis_policy_id,
            "right_handed_basis_policy_id": self.right_handed_basis_policy_id,
            "mu_gs_GeV": self.mu_gs_GeV,
            "mu_match_GeV": self.mu_match_GeV,
            "m_g1_GeV": self.m_g1_GeV,
            "propagator_mass_GeV": self.propagator_mass_GeV,
            "alpha_s_mu_gs": self.alpha_s_mu_gs,
            "g_s_mu_gs": self.g_s_mu_gs,
            "theta12_deg": self.theta12_deg,
            "theta23_deg": self.theta23_deg,
            "theta13_deg": self.theta13_deg,
            "delta": self.delta,
            "matrix_summaries": [item.as_dict() for item in self.matrix_summaries],
            "notes": self.notes,
        }


@dataclass(frozen=True)
class Paper07101869PhysicalKKGluonBenchmarkSummary(
    Paper07101869KKGluonBenchmarkSummary
):
    """Physical-side summary wrapper for the KK-gluon adapter."""

    physical_profile_status_id: str = PAPER_0710_1869_PHYSICAL_KK_GLUON_STATUS_ID

    def __post_init__(self) -> None:
        super().__post_init__()
        object.__setattr__(
            self,
            "physical_profile_status_id",
            require_member(
                "physical_profile_status_id",
                self.physical_profile_status_id,
                (PAPER_0710_1869_PHYSICAL_KK_GLUON_STATUS_ID,),
            ),
        )
        require_member(
            "benchmark_status",
            self.benchmark_status,
            (PAPER_0710_1869_PHYSICAL_KK_GLUON_STATUS_ID,),
        )

    def as_dict(self) -> dict[str, object]:
        payload = super().as_dict()
        payload["physical_profile_status_id"] = self.physical_profile_status_id
        return payload


def _real_triplet_tuple(
    name: str, values: tuple[float, float, float] | list[float] | np.ndarray
) -> tuple[float, float, float]:
    arr = _as_real_triplet(name, values)
    return tuple(float(value) for value in arr)


def _allowed_kk_gluon_status_ids() -> tuple[str, ...]:
    return (
        "sourced_structural_only",
        PAPER_0710_1869_PHYSICAL_KK_GLUON_STATUS_ID,
    )


def _physical_benchmark_metadata(
    physical_bulk_state: Paper07101869PhysicalBulkState,
) -> tuple[str, str, float, float, float, float]:
    metadata = physical_bulk_state.point.metadata
    benchmark_id = physical_bulk_state.point.label
    benchmark_label = physical_bulk_state.point.label
    theta12_deg = 0.0
    theta23_deg = 0.0
    theta13_deg = 0.0
    delta = 0.0

    if isinstance(metadata, Mapping):
        raw_benchmark_id = metadata.get("benchmark_id")
        raw_benchmark_label = metadata.get("benchmark_label")
        if isinstance(raw_benchmark_id, str) and raw_benchmark_id:
            benchmark_id = raw_benchmark_id
        if isinstance(raw_benchmark_label, str) and raw_benchmark_label:
            benchmark_label = raw_benchmark_label
        raw_eq3 = metadata.get("benchmark_eq3_example")
        if isinstance(raw_eq3, Mapping):
            theta12_deg = float(raw_eq3.get("theta12_deg", theta12_deg))
            theta23_deg = float(raw_eq3.get("theta23_deg", theta23_deg))
            theta13_deg = float(raw_eq3.get("theta13_deg", theta13_deg))
            delta = float(raw_eq3.get("delta", delta))
    return benchmark_id, benchmark_label, theta12_deg, theta23_deg, theta13_deg, delta


def _build_matrix_block(
    *,
    label: str,
    sector_id: str,
    chiral_id: str,
    basis_id: str,
    raw_dimensionless: np.ndarray,
    normalization: Paper07101869GaugeCouplingNormalization,
) -> Paper07101869KKGluonFlavorMatrix:
    raw = _hermitian(raw_dimensionless)
    subtracted, universal_component = _subtract_universal_piece(raw)
    g_s = normalization.g_s_mu_gs
    sqrt_2L = normalization.rs_volume_sqrt_2L
    return Paper07101869KKGluonFlavorMatrix(
        label=label,
        sector_id=sector_id,
        chiral_id=chiral_id,
        basis_id=basis_id,
        g_s_mu_gs=g_s,
        rs_volume_L=normalization.rs_volume_L,
        rs_volume_sqrt_2L=sqrt_2L,
        universal_component=universal_component,
        raw_dimensionless=raw,
        universal_subtracted_dimensionless=subtracted,
        raw_gs_normalized=_hermitian(
            g_s * (sqrt_2L * raw - (np.eye(3, dtype=np.complex128) / sqrt_2L))
        ),
        universal_subtracted_gs_normalized=_hermitian(g_s * sqrt_2L * subtracted),
    )


def _build_physical_matrix_block(
    *,
    label: str,
    sector_id: str,
    chiral_id: str,
    basis_id: str,
    profiles: tuple[float, float, float] | list[float] | np.ndarray,
    rotation: list[list[float]] | list[list[complex]] | np.ndarray,
    normalization: Paper07101869GaugeCouplingNormalization,
) -> Paper07101869KKGluonFlavorMatrix:
    profile_vector = _as_real_triplet("profiles", profiles)
    rotation_matrix = _as_complex_matrix3("rotation", np.asarray(rotation, dtype=np.complex128))
    raw_dimensionless = _hermitian(
        rotation_matrix.conjugate().T
        @ np.diag(profile_vector**2).astype(np.complex128)
        @ rotation_matrix
    )
    return _build_matrix_block(
        label=label,
        sector_id=sector_id,
        chiral_id=chiral_id,
        basis_id=basis_id,
        raw_dimensionless=raw_dimensionless,
        normalization=normalization,
    )


def build_paper_0710_1869_kk_gluon_couplings(
    benchmark: Paper07101869Benchmark | None = None,
    *,
    scale_point: Paper07101869ScalePoint | None = None,
    contract: Paper07101869CouplingContract | None = None,
) -> Paper07101869KKGluonCouplings:
    """Build the first deterministic KK-gluon benchmark layer from paper inputs."""
    resolved_benchmark = (
        default_paper_0710_1869_pr1_benchmark() if benchmark is None else benchmark
    )
    if not isinstance(resolved_benchmark, Paper07101869Benchmark):
        raise ValueError("benchmark must be a Paper07101869Benchmark")

    resolved_scale = default_paper_0710_1869_scales() if scale_point is None else scale_point
    if not isinstance(resolved_scale, Paper07101869ScalePoint):
        raise ValueError("scale_point must be a Paper07101869ScalePoint")

    resolved_contract = (
        default_paper_0710_1869_coupling_contract() if contract is None else contract
    )
    if not isinstance(resolved_contract, Paper07101869CouplingContract):
        raise ValueError("contract must be a Paper07101869CouplingContract")

    normalization = evaluate_paper_0710_1869_gauge_coupling(
        resolved_scale,
        contract=resolved_contract,
    )
    table_benchmark = build_table_i_benchmark_from_inputs(resolved_benchmark.table_i_inputs)
    v5km_parameters = Paper07101869V5KMParameters.from_degrees(
        theta12_deg=resolved_benchmark.eq3_example.theta12_deg,
        theta23_deg=resolved_benchmark.eq3_example.theta23_deg,
        theta13_deg=resolved_benchmark.eq3_example.theta13_deg,
        delta=resolved_benchmark.eq3_example.delta,
    )
    v5km = paper_v5km_matrix(v5km_parameters)

    q_kernel = _diagonal_profile_kernel(
        "q_sector.f_eigenvalues",
        table_benchmark.q_sector.f_eigenvalues,
    )
    u_kernel = _diagonal_profile_kernel(
        "u_sector.f_eigenvalues",
        table_benchmark.u_sector.f_eigenvalues,
    )
    d_kernel = _diagonal_profile_kernel(
        "d_sector.f_eigenvalues",
        table_benchmark.d_sector.f_eigenvalues,
    )

    left_up_aligned = _build_matrix_block(
        label="left_up_aligned",
        sector_id="Q",
        chiral_id="left",
        basis_id=PAPER_0710_1869_LEFT_UP_ALIGNED_BASIS_ID,
        raw_dimensionless=q_kernel,
        normalization=normalization,
    )
    left_down_aligned = _build_matrix_block(
        label="left_down_aligned",
        sector_id="Q",
        chiral_id="left",
        basis_id=PAPER_0710_1869_LEFT_DOWN_ALIGNED_BASIS_ID,
        raw_dimensionless=_hermitian(v5km.conjugate().T @ q_kernel @ v5km),
        normalization=normalization,
    )
    right_up = _build_matrix_block(
        label="right_up_diagonal",
        sector_id="u",
        chiral_id="right",
        basis_id=PAPER_0710_1869_RIGHT_DIAGONAL_BASIS_ID,
        raw_dimensionless=u_kernel,
        normalization=normalization,
    )
    # RESIDUAL(C-2): default RH-down alignment model choice pending paper 0710.1869.
    right_down = _build_matrix_block(
        label="right_down_diagonal",
        sector_id="d",
        chiral_id="right",
        basis_id=PAPER_0710_1869_RIGHT_DIAGONAL_BASIS_ID,
        raw_dimensionless=d_kernel,
        normalization=normalization,
    )

    return Paper07101869KKGluonCouplings(
        contract=resolved_contract,
        normalization=normalization,
        left_up_aligned=left_up_aligned,
        left_down_aligned=left_down_aligned,
        right_up=right_up,
        right_down=right_down,
        benchmark_id=resolved_benchmark.benchmark_id,
        benchmark_label=resolved_benchmark.label,
        benchmark_status=resolved_benchmark.status,
        theta12_deg=resolved_benchmark.eq3_example.theta12_deg,
        theta23_deg=resolved_benchmark.eq3_example.theta23_deg,
        theta13_deg=resolved_benchmark.eq3_example.theta13_deg,
        delta=resolved_benchmark.eq3_example.delta,
    )


def build_paper_0710_1869_kk_gluon_couplings_from_physical_bulk_state(
    physical_bulk_state: Paper07101869PhysicalBulkState,
    *,
    scale_point: Paper07101869ScalePoint | None = None,
    contract: Paper07101869CouplingContract | None = None,
) -> Paper07101869PhysicalKKGluonCouplings:
    """Build the honest QS2 physical KK-gluon adapter from a QS1 physical bulk state."""
    if not isinstance(physical_bulk_state, Paper07101869PhysicalBulkState):
        raise ValueError("physical_bulk_state must be a Paper07101869PhysicalBulkState")
    if (
        physical_bulk_state.physical_profile_status_id
        != PAPER_0710_1869_PHYSICAL_KK_GLUON_STATUS_ID
    ):
        raise ValueError(
            "physical_bulk_state must carry the exact point-derived QS1 physical status"
        )

    resolved_scale = default_paper_0710_1869_scales() if scale_point is None else scale_point
    if not isinstance(resolved_scale, Paper07101869ScalePoint):
        raise ValueError("scale_point must be a Paper07101869ScalePoint")

    resolved_contract = (
        default_paper_0710_1869_coupling_contract() if contract is None else contract
    )
    if not isinstance(resolved_contract, Paper07101869CouplingContract):
        raise ValueError("contract must be a Paper07101869CouplingContract")

    normalization = evaluate_paper_0710_1869_gauge_coupling(
        resolved_scale,
        contract=resolved_contract,
    )
    M_u, M_d = build_paper_0710_1869_physical_mass_matrices(physical_bulk_state)
    observables = paper_physical_mass_ckm_probe_observables(M_u, M_d)

    left_up_aligned = _build_physical_matrix_block(
        label="left_up_aligned",
        sector_id="Q",
        chiral_id="left",
        basis_id=PAPER_0710_1869_PHYSICAL_LEFT_UP_ALIGNED_BASIS_ID,
        profiles=physical_bulk_state.F_Q,
        rotation=observables["U_L_u"],
        normalization=normalization,
    )
    left_down_aligned = _build_physical_matrix_block(
        label="left_down_aligned",
        sector_id="Q",
        chiral_id="left",
        basis_id=PAPER_0710_1869_PHYSICAL_LEFT_DOWN_ALIGNED_BASIS_ID,
        profiles=physical_bulk_state.F_Q,
        rotation=observables["U_L_d"],
        normalization=normalization,
    )
    right_up = _build_physical_matrix_block(
        label="right_up_aligned",
        sector_id="u",
        chiral_id="right",
        basis_id=PAPER_0710_1869_PHYSICAL_RIGHT_UP_ALIGNED_BASIS_ID,
        profiles=physical_bulk_state.F_u,
        rotation=observables["U_R_u"],
        normalization=normalization,
    )
    right_down = _build_physical_matrix_block(
        label="right_down_aligned",
        sector_id="d",
        chiral_id="right",
        basis_id=PAPER_0710_1869_PHYSICAL_RIGHT_DOWN_ALIGNED_BASIS_ID,
        profiles=physical_bulk_state.F_d,
        rotation=observables["U_R_d"],
        normalization=normalization,
    )
    benchmark_id, benchmark_label, theta12_deg, theta23_deg, theta13_deg, delta = (
        _physical_benchmark_metadata(physical_bulk_state)
    )

    return Paper07101869PhysicalKKGluonCouplings(
        contract=resolved_contract,
        normalization=normalization,
        left_up_aligned=left_up_aligned,
        left_down_aligned=left_down_aligned,
        right_up=right_up,
        right_down=right_down,
        benchmark_id=benchmark_id,
        benchmark_label=benchmark_label,
        benchmark_status=physical_bulk_state.physical_profile_status_id,
        theta12_deg=theta12_deg,
        theta23_deg=theta23_deg,
        theta13_deg=theta13_deg,
        delta=delta,
        left_handed_basis_policy_id=PAPER_0710_1869_PHYSICAL_LEFT_HANDED_BASIS_POLICY_ID,
        right_handed_basis_policy_id=PAPER_0710_1869_PHYSICAL_RIGHT_HANDED_BASIS_POLICY_ID,
        physical_profile_status_id=physical_bulk_state.physical_profile_status_id,
    )


def default_paper_0710_1869_kk_gluon_couplings() -> Paper07101869KKGluonCouplings:
    """Return the default deterministic KK-gluon benchmark layer."""
    return build_paper_0710_1869_kk_gluon_couplings()


def build_paper_0710_1869_kk_gluon_benchmark_summary(
    benchmark: Paper07101869Benchmark | None = None,
    *,
    scale_point: Paper07101869ScalePoint | None = None,
    contract: Paper07101869CouplingContract | None = None,
) -> Paper07101869KKGluonBenchmarkSummary:
    """Build the deterministic benchmark summary for the KK-gluon layer."""
    return build_paper_0710_1869_kk_gluon_couplings(
        benchmark,
        scale_point=scale_point,
        contract=contract,
    ).summary()


def default_paper_0710_1869_kk_gluon_benchmark_summary() -> Paper07101869KKGluonBenchmarkSummary:
    """Return the default benchmark summary for the KK-gluon layer."""
    return build_paper_0710_1869_kk_gluon_benchmark_summary()


__all__ = [
    "PAPER_0710_1869_KK_GLUON_MATRIX_SCHEMA_ID",
    "PAPER_0710_1869_KK_GLUON_SCHEMA_ID",
    "PAPER_0710_1869_KK_GLUON_SUMMARY_SCHEMA_ID",
    "PAPER_0710_1869_LEFT_DOWN_ALIGNED_BASIS_ID",
    "PAPER_0710_1869_LEFT_HANDED_BASIS_POLICY_ID",
    "PAPER_0710_1869_LEFT_UP_ALIGNED_BASIS_ID",
    "PAPER_0710_1869_PHYSICAL_KK_GLUON_STATUS_ID",
    "PAPER_0710_1869_PHYSICAL_LEFT_DOWN_ALIGNED_BASIS_ID",
    "PAPER_0710_1869_PHYSICAL_LEFT_HANDED_BASIS_POLICY_ID",
    "PAPER_0710_1869_PHYSICAL_LEFT_UP_ALIGNED_BASIS_ID",
    "PAPER_0710_1869_PHYSICAL_RIGHT_DOWN_ALIGNED_BASIS_ID",
    "PAPER_0710_1869_PHYSICAL_RIGHT_HANDED_BASIS_POLICY_ID",
    "PAPER_0710_1869_PHYSICAL_RIGHT_UP_ALIGNED_BASIS_ID",
    "PAPER_0710_1869_RIGHT_DIAGONAL_BASIS_ID",
    "PAPER_0710_1869_RIGHT_HANDED_BASIS_POLICY_ID",
    "Paper07101869KKGluonBenchmarkSummary",
    "Paper07101869KKGluonCouplings",
    "Paper07101869KKGluonFlavorMatrix",
    "Paper07101869KKGluonMatrixSummary",
    "Paper07101869PhysicalKKGluonBenchmarkSummary",
    "Paper07101869PhysicalKKGluonCouplings",
    "build_paper_0710_1869_kk_gluon_benchmark_summary",
    "build_paper_0710_1869_kk_gluon_couplings",
    "build_paper_0710_1869_kk_gluon_couplings_from_physical_bulk_state",
    "default_paper_0710_1869_kk_gluon_benchmark_summary",
    "default_paper_0710_1869_kk_gluon_couplings",
]
