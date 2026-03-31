"""Core MFV quark-sector model primitives.

This module keeps the quark-sector implementation MFV-native:

- the primary flavor inputs are the Yukawa spurions ``Y_u`` and ``Y_d``
- bulk masses are derived from ``C_Q``, ``C_u``, and ``C_d``
- zero-mode overlaps reuse ``warpConfig.wavefuncs.f_IR``
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, Mapping, Optional

import numpy as np

from warpConfig.baseParams import DEFAULT_LAMBDA_IR, MPL, V_EWSB, get_warp_params
from warpConfig.wavefuncs import f_IR


def _as_real_vector(name: str, values: np.ndarray | list[float] | tuple[float, ...]) -> np.ndarray:
    """Return a length-3 real vector."""
    arr = np.asarray(values, dtype=float)
    if arr.shape != (3,):
        raise ValueError(f"{name} must have shape (3,)")
    return arr


def _as_complex_matrix(name: str, values: np.ndarray | list[list[complex]]) -> np.ndarray:
    """Return a complex 3x3 matrix."""
    arr = np.asarray(values, dtype=np.complex128)
    if arr.shape != (3, 3):
        raise ValueError(f"{name} must have shape (3, 3)")
    return arr


@dataclass(frozen=True)
class RotationParameters:
    """CKM-like unitary parameterization."""

    theta12: float = 0.0
    theta13: float = 0.0
    theta23: float = 0.0
    delta: float = 0.0


@dataclass(frozen=True)
class BulkMassMap:
    """Bounded map from positive spurion eigenvalues to RS bulk masses."""

    c_uv: float = 0.72
    c_ir: float = 0.30
    eigen_scale: float = 1.0

    def map_eigenvalues(self, eigenvalues: np.ndarray) -> np.ndarray:
        """Map non-negative eigenvalues onto the physical bulk-mass window."""
        vals = np.asarray(eigenvalues, dtype=float)
        if np.any(vals < -1e-12):
            raise ValueError("bulk-mass mapping requires non-negative eigenvalues")
        vals = np.maximum(vals, 0.0)
        scaled = vals / (vals + self.eigen_scale)
        return self.c_uv - (self.c_uv - self.c_ir) * scaled


@dataclass(frozen=True)
class QuarkSpurionPoint:
    """MFV-native quark-sector input point."""

    Y_u: np.ndarray
    Y_d: np.ndarray
    r: float
    Lambda_IR: float = DEFAULT_LAMBDA_IR
    k: float = MPL
    v: float = V_EWSB
    label: str = "anonymous"
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(self, "Y_u", _as_complex_matrix("Y_u", self.Y_u))
        object.__setattr__(self, "Y_d", _as_complex_matrix("Y_d", self.Y_d))
        if not np.isfinite(self.r):
            raise ValueError("r must be finite")
        if self.r < 0.0:
            raise ValueError("r must be non-negative in the current MFV implementation")
        if self.Lambda_IR <= 0.0:
            raise ValueError("Lambda_IR must be positive")
        if self.k <= 0.0:
            raise ValueError("k must be positive")
        if self.v <= 0.0:
            raise ValueError("v must be positive")


@dataclass(frozen=True)
class QuarkBulkState:
    """Derived bulk-mass and overlap data."""

    point: QuarkSpurionPoint
    bulk_mass_map: BulkMassMap
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


def ckm_like_unitary(params: RotationParameters) -> np.ndarray:
    """Return a CKM-like unitary matrix from three angles and one phase."""
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


def unitary_angle_summary(U: np.ndarray) -> Dict[str, float]:
    """Return a compact angle summary from a CKM-like unitary."""
    U = _as_complex_matrix("U", U)
    theta13 = float(np.arcsin(np.clip(np.abs(U[0, 2]), 0.0, 1.0)))
    theta12 = float(np.arctan2(np.abs(U[0, 1]), max(np.abs(U[0, 0]), 1e-15)))
    theta23 = float(np.arctan2(np.abs(U[1, 2]), max(np.abs(U[2, 2]), 1e-15)))
    return {
        "theta12": theta12,
        "theta13": theta13,
        "theta23": theta23,
    }


def build_spurion_matrix(
    singular_values: np.ndarray | list[float] | tuple[float, ...],
    *,
    overall_scale: float = 1.0,
    left_rotation: Optional[RotationParameters] = None,
    right_rotation: Optional[RotationParameters] = None,
) -> np.ndarray:
    """Construct a Yukawa spurion from singular values and unitary rotations."""
    sigma = _as_real_vector("singular_values", singular_values)
    if np.any(sigma < 0.0):
        raise ValueError("singular values must be non-negative")
    U_left = ckm_like_unitary(left_rotation or RotationParameters())
    U_right = ckm_like_unitary(right_rotation or RotationParameters())
    return overall_scale * U_left @ np.diag(sigma) @ U_right.conj().T


def build_mfv_point_from_singular_values(
    *,
    up_singular_values: np.ndarray | list[float] | tuple[float, ...],
    down_singular_values: np.ndarray | list[float] | tuple[float, ...],
    overall_scale: float = 1.0,
    r: float,
    up_left: Optional[RotationParameters] = None,
    up_right: Optional[RotationParameters] = None,
    down_left: Optional[RotationParameters] = None,
    down_right: Optional[RotationParameters] = None,
    Lambda_IR: float = DEFAULT_LAMBDA_IR,
    k: float = MPL,
    v: float = V_EWSB,
    label: str = "constructed",
    metadata: Optional[Mapping[str, Any]] = None,
) -> QuarkSpurionPoint:
    """Build a deterministic MFV point from singular values plus rotations."""
    Y_u = build_spurion_matrix(
        up_singular_values,
        overall_scale=overall_scale,
        left_rotation=up_left,
        right_rotation=up_right,
    )
    Y_d = build_spurion_matrix(
        down_singular_values,
        overall_scale=overall_scale,
        left_rotation=down_left,
        right_rotation=down_right,
    )
    return QuarkSpurionPoint(
        Y_u=Y_u,
        Y_d=Y_d,
        r=r,
        Lambda_IR=Lambda_IR,
        k=k,
        v=v,
        label=label,
        metadata=dict(metadata or {}),
    )


def _ordered_hermitian_spectrum(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return eigenvalues/eigenvectors ordered from light to heavy generation."""
    vals, vecs = np.linalg.eigh(matrix)
    # Light-to-heavy generation labeling: smallest spurion eigenvalue first,
    # largest third-generation eigenvalue last.
    idx = np.argsort(vals.real)
    vals = np.real_if_close(vals[idx]).astype(float)
    vecs = vecs[:, idx]
    return vals, vecs


def derive_bulk_state(
    point: QuarkSpurionPoint,
    *,
    bulk_mass_map: Optional[BulkMassMap] = None,
) -> QuarkBulkState:
    """Construct derived bulk masses and overlaps from the MFV spurions."""
    bulk_mass_map = bulk_mass_map or BulkMassMap()
    params = get_warp_params(k=point.k, Lambda_IR=point.Lambda_IR)
    epsilon = float(params["epsilon"])

    C_u = point.Y_u.conj().T @ point.Y_u
    C_d = point.Y_d.conj().T @ point.Y_d
    C_Q = point.r * (point.Y_u @ point.Y_u.conj().T) + (point.Y_d @ point.Y_d.conj().T)

    eig_Q, rotation_Q = _ordered_hermitian_spectrum(C_Q)
    eig_u, rotation_u = _ordered_hermitian_spectrum(C_u)
    eig_d, rotation_d = _ordered_hermitian_spectrum(C_d)

    c_Q = bulk_mass_map.map_eigenvalues(eig_Q)
    c_u = bulk_mass_map.map_eigenvalues(eig_u)
    c_d = bulk_mass_map.map_eigenvalues(eig_d)

    F_Q = np.asarray(f_IR(c_Q, epsilon), dtype=float)
    F_u = np.asarray(f_IR(c_u, epsilon), dtype=float)
    F_d = np.asarray(f_IR(c_d, epsilon), dtype=float)

    Y_u_bulk_basis = rotation_Q.conj().T @ point.Y_u @ rotation_u
    Y_d_bulk_basis = rotation_Q.conj().T @ point.Y_d @ rotation_d

    return QuarkBulkState(
        point=point,
        bulk_mass_map=bulk_mass_map,
        epsilon=epsilon,
        C_Q=C_Q,
        C_u=C_u,
        C_d=C_d,
        c_Q=c_Q,
        c_u=c_u,
        c_d=c_d,
        F_Q=F_Q,
        F_u=F_u,
        F_d=F_d,
        rotation_Q=rotation_Q,
        rotation_u=rotation_u,
        rotation_d=rotation_d,
        eig_Q=eig_Q,
        eig_u=eig_u,
        eig_d=eig_d,
        Y_u_bulk_basis=Y_u_bulk_basis,
        Y_d_bulk_basis=Y_d_bulk_basis,
    )


def spurion_svd_summary(point: QuarkSpurionPoint) -> Dict[str, Dict[str, Any]]:
    """Return singular-value and mixing summaries for the two spurions."""
    out: Dict[str, Dict[str, Any]] = {}
    for label, matrix in (("up", point.Y_u), ("down", point.Y_d)):
        U_l, singular_values, U_r_h = np.linalg.svd(matrix, full_matrices=True)
        out[label] = {
            "singular_values": singular_values,
            "left_angles": unitary_angle_summary(U_l),
            "right_angles": unitary_angle_summary(U_r_h.conj().T),
        }
    return out
