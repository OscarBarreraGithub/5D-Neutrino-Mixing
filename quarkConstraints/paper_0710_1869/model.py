"""Paper-facing benchmark helpers for the 0710.1869 Eq. (3) relation."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .inputs import (
    Paper07101869Eq3Example as Paper07101869Eq3Inputs,
)
from .inputs import (
    Paper07101869TableIInputs,
    default_paper_0710_1869_table_i_inputs,
)
from .inputs import (
    default_paper_0710_1869_eq3_example as default_paper_0710_1869_eq3_inputs,
)
from .validation import require_nonempty_identifier, require_positive_finite

PAPER_0710_1869_MODEL_SCHEMA_ID = "quarkConstraints.paper_0710_1869.model.v1"
PAPER_0710_1869_EQ3_RELATION_ID = "diag(C_Q)=a*diag(r*V5KM^dagger*C_u*V5KM + C_d)"
PAPER_0710_1869_TABLE_I_REFERENCE = "Table I, arXiv:0710.1869v1"
PAPER_0710_1869_EQ3_REFERENCE = "Eq. (3), arXiv:0710.1869v1"


def _as_real_triplet(
    name: str, values: tuple[float, float, float] | list[float] | np.ndarray
) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.shape != (3,):
        raise ValueError(f"{name} must have shape (3,)")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must contain only finite values")
    return arr.astype(float, copy=True)


def _as_complex_matrix3(name: str, values: np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=np.complex128)
    if arr.shape != (3, 3):
        raise ValueError(f"{name} must have shape (3, 3)")
    if not np.all(np.isfinite(arr.real)) or not np.all(np.isfinite(arr.imag)):
        raise ValueError(f"{name} must contain only finite entries")
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
            tuple(_as_real_triplet("f_eigenvalues", self.f_eigenvalues)),
        )

    @property
    def c_vector(self) -> np.ndarray:
        return _as_real_triplet("c_eigenvalues", self.c_eigenvalues)

    @property
    def f_vector(self) -> np.ndarray:
        return _as_real_triplet("f_eigenvalues", self.f_eigenvalues)

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


def paper_v5km_matrix(params: Paper07101869V5KMParameters) -> np.ndarray:
    """Return the paper-side CKM-like ``V5KM`` matrix."""
    if not isinstance(params, Paper07101869V5KMParameters):
        raise ValueError("params must be a Paper07101869V5KMParameters")
    s12, s23, s13 = np.sin([params.theta12, params.theta23, params.theta13])
    c12, c23, c13 = np.cos([params.theta12, params.theta23, params.theta13])
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
        else _as_complex_matrix3("v5km", v5km)
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
    "PAPER_0710_1869_EQ3_REFERENCE",
    "PAPER_0710_1869_EQ3_RELATION_ID",
    "PAPER_0710_1869_MODEL_SCHEMA_ID",
    "PAPER_0710_1869_TABLE_I_REFERENCE",
    "Paper07101869BenchmarkSector",
    "Paper07101869DiagonalCMatrices",
    "Paper07101869Eq3Example",
    "Paper07101869Eq3ResidualSummary",
    "Paper07101869TableIBenchmark",
    "Paper07101869V5KMParameters",
    "build_diagonal_c_matrix",
    "build_table_i_benchmark_from_inputs",
    "build_table_i_diagonal_c_matrices",
    "default_paper_0710_1869_eq3_example",
    "default_paper_0710_1869_table_i_benchmark",
    "evaluate_default_paper_0710_1869_eq3_consistency",
    "evaluate_eq3_diagonal_consistency",
    "paper_v5km_matrix",
]
