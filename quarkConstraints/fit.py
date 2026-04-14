"""Exact quark mass-matrix evaluation and MFV-native fit utilities."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Mapping, Optional

import numpy as np
from scipy.optimize import least_squares

from diagonalization.diag import SVD

from .model import (
    BulkMassMap,
    QuarkBulkState,
    QuarkSpurionPoint,
    RotationParameters,
    build_mfv_point_from_singular_values,
    ckm_like_unitary,
    derive_bulk_state,
    spurion_svd_summary,
)


def _as_real_vector(name: str, values: np.ndarray | list[float] | tuple[float, ...]) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.shape != (3,):
        raise ValueError(f"{name} must have shape (3,)")
    if np.any(arr <= 0.0):
        raise ValueError(f"{name} must be strictly positive")
    return arr


def _as_complex_matrix(name: str, values: np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=np.complex128)
    if arr.shape != (3, 3):
        raise ValueError(f"{name} must have shape (3, 3)")
    return arr


def jarlskog_invariant(ckm: np.ndarray) -> float:
    """Return the CKM Jarlskog invariant."""
    v = np.asarray(ckm, dtype=np.complex128)
    return float(np.imag(v[0, 1] * v[1, 2] * np.conjugate(v[0, 2]) * np.conjugate(v[1, 1])))


def ckm_observables(ckm: np.ndarray) -> np.ndarray:
    """Return the reduced CKM observable set used for fits and plots."""
    v = np.asarray(ckm, dtype=np.complex128)
    return np.array(
        [abs(v[0, 1]), abs(v[1, 2]), abs(v[0, 2]), jarlskog_invariant(v)],
        dtype=float,
    )


@dataclass(frozen=True)
class QuarkTargets:
    """Target masses and CKM matrix for fit scoring."""

    up_masses: np.ndarray
    down_masses: np.ndarray
    ckm: np.ndarray
    label: str = "targets"

    def __post_init__(self) -> None:
        object.__setattr__(self, "up_masses", _as_real_vector("up_masses", self.up_masses))
        object.__setattr__(self, "down_masses", _as_real_vector("down_masses", self.down_masses))
        object.__setattr__(self, "ckm", _as_complex_matrix("ckm", self.ckm))

    @property
    def abs_ckm(self) -> np.ndarray:
        return np.abs(self.ckm)

    @property
    def ckm_observables(self) -> np.ndarray:
        return ckm_observables(self.ckm)


@dataclass(frozen=True)
class QuarkFitResult:
    """Exact quark-sector observables derived from one MFV spurion point."""

    bulk_state: QuarkBulkState
    M_u: np.ndarray
    M_d: np.ndarray
    U_L_u: np.ndarray
    U_R_u: np.ndarray
    U_L_d: np.ndarray
    U_R_d: np.ndarray
    masses_up: np.ndarray
    masses_down: np.ndarray
    ckm: np.ndarray
    spurion_summary: Mapping[str, Dict[str, Any]]
    up_log_residuals: Optional[np.ndarray] = None
    down_log_residuals: Optional[np.ndarray] = None
    ckm_abs_residuals: Optional[np.ndarray] = None
    ckm_observable_residuals: Optional[np.ndarray] = None
    total_score: Optional[float] = None

    @property
    def state(self) -> QuarkBulkState:
        return self.bulk_state

    @property
    def point(self) -> QuarkSpurionPoint:
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
        return ckm_observables(self.ckm)

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
        """Return a formatted fit summary."""
        lines = [
            "=" * 60,
            f"Quark MFV fit summary: {self.bulk_state.point.label}",
            "=" * 60,
            f"r = {self.bulk_state.point.r:.4f}",
            f"epsilon = {self.bulk_state.epsilon:.3e}",
            f"c_Q = {np.array2string(self.bulk_state.c_Q, precision=4)}",
            f"c_u = {np.array2string(self.bulk_state.c_u, precision=4)}",
            f"c_d = {np.array2string(self.bulk_state.c_d, precision=4)}",
            f"F_Q = {np.array2string(self.bulk_state.F_Q, precision=5)}",
            f"F_u = {np.array2string(self.bulk_state.F_u, precision=5)}",
            f"F_d = {np.array2string(self.bulk_state.F_d, precision=5)}",
            f"m_u = {np.array2string(self.masses_up, precision=5)}",
            f"m_d = {np.array2string(self.masses_down, precision=5)}",
            f"CKM obs = {np.array2string(self.ckm_observables, precision=6)}",
        ]
        if self.total_score is not None:
            lines.append(f"fit score = {self.total_score:.6e}")
        return "\n".join(lines)


@dataclass(frozen=True)
class QuarkFitSeed:
    """Optimization seed in spurion space."""

    up_singular_values: np.ndarray
    down_singular_values: np.ndarray
    overall_scale: float
    up_left: RotationParameters
    up_right: RotationParameters
    down_left: RotationParameters
    down_right: RotationParameters


@dataclass(frozen=True)
class QuarkFitSolution:
    """Optimized seed together with the fitted result."""

    seed: QuarkFitSeed
    result: QuarkFitResult
    success: bool
    message: str
    nfev: int
    initial_score: float

    @property
    def parameters(self) -> QuarkFitSeed:
        return self.seed


QUARK_FIT_CANONICAL_VECTOR_FIELDS = (
    "log_up_physical_singular_values[0]",
    "log_up_physical_singular_values[1]",
    "log_up_physical_singular_values[2]",
    "log_down_physical_singular_values[0]",
    "log_down_physical_singular_values[1]",
    "log_down_physical_singular_values[2]",
    "up_left.theta12",
    "up_left.theta13",
    "up_left.theta23",
    "up_left.delta",
    "down_left.theta12",
    "down_left.theta13",
    "down_left.theta23",
    "down_left.delta",
)
QUARK_FIT_CANONICAL_VECTOR_SIZE = len(QUARK_FIT_CANONICAL_VECTOR_FIELDS)
QUARK_FIT_FULL_VECTOR_FIELDS = QUARK_FIT_CANONICAL_VECTOR_FIELDS
QUARK_FIT_FULL_VECTOR_SIZE = QUARK_FIT_CANONICAL_VECTOR_SIZE


def _ordered_dirac_svd(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return SVD ordered from light to heavy generation."""
    U_l, singular_values, U_r = SVD(matrix)
    idx = np.argsort(singular_values)
    return U_l[:, idx], singular_values[idx], U_r[:, idx]


def build_mass_matrices(bulk_state: QuarkBulkState) -> tuple[np.ndarray, np.ndarray]:
    """Construct the full quark mass matrices in the bulk-mass basis."""
    F_Q = np.diag(bulk_state.F_Q)
    F_u = np.diag(bulk_state.F_u)
    F_d = np.diag(bulk_state.F_d)
    prefactor = 2.0 * bulk_state.point.v
    M_u = prefactor * F_Q @ bulk_state.Y_u_bulk_basis @ F_u
    M_d = prefactor * F_Q @ bulk_state.Y_d_bulk_basis @ F_d
    return M_u, M_d


def mass_matrix_observables(M_u: np.ndarray, M_d: np.ndarray) -> Dict[str, np.ndarray]:
    """Extract masses and CKM from the two Dirac mass matrices."""
    U_L_u, masses_up, U_R_u = _ordered_dirac_svd(M_u)
    U_L_d, masses_down, U_R_d = _ordered_dirac_svd(M_d)
    ckm = U_L_u.conj().T @ U_L_d
    return {
        "U_L_u": U_L_u,
        "U_R_u": U_R_u,
        "U_L_d": U_L_d,
        "U_R_d": U_R_d,
        "masses_up": masses_up,
        "masses_down": masses_down,
        "ckm": ckm,
    }


def fit_residuals(
    masses_up: np.ndarray,
    masses_down: np.ndarray,
    ckm: np.ndarray,
    targets: QuarkTargets,
) -> Dict[str, Any]:
    """Compute mass and CKM residuals against targets.

    CKM observable residuals are normalized by the absolute target observable
    scales. That scale policy is explicit: zero target observable scales are
    unsupported and raise ``ValueError``.
    """
    up_log_residuals = np.log(np.maximum(masses_up, 1e-30) / targets.up_masses)
    down_log_residuals = np.log(np.maximum(masses_down, 1e-30) / targets.down_masses)
    ckm_abs_residuals = np.abs(np.abs(ckm) - targets.abs_ckm)
    ckm_scales = np.abs(targets.ckm_observables)
    if np.any(ckm_scales <= 0.0):
        raise ValueError("fit_residuals requires non-zero CKM observable targets")
    ckm_obs_residuals = (ckm_observables(ckm) - targets.ckm_observables) / ckm_scales

    total_score = float(
        np.sqrt(
            np.mean(
                np.concatenate(
                    [
                        up_log_residuals,
                        down_log_residuals,
                        ckm_obs_residuals,
                    ]
                )
                ** 2
            )
        )
    )

    return {
        "up_log_residuals": up_log_residuals,
        "down_log_residuals": down_log_residuals,
        "ckm_abs_residuals": ckm_abs_residuals,
        "ckm_observable_residuals": ckm_obs_residuals,
        "total_score": total_score,
    }


def evaluate_quark_fit(
    point: QuarkSpurionPoint,
    targets: Optional[QuarkTargets] = None,
    *,
    bulk_mass_map: Optional[BulkMassMap] = None,
) -> QuarkFitResult:
    """Evaluate one MFV point exactly."""
    bulk_state = derive_bulk_state(point, bulk_mass_map=bulk_mass_map)
    M_u, M_d = build_mass_matrices(bulk_state)
    observables = mass_matrix_observables(M_u, M_d)
    residuals: Dict[str, Any] = {}
    if targets is not None:
        residuals = fit_residuals(
            observables["masses_up"],
            observables["masses_down"],
            observables["ckm"],
            targets,
        )

    return QuarkFitResult(
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
        spurion_summary=spurion_svd_summary(point),
        up_log_residuals=residuals.get("up_log_residuals"),
        down_log_residuals=residuals.get("down_log_residuals"),
        ckm_abs_residuals=residuals.get("ckm_abs_residuals"),
        ckm_observable_residuals=residuals.get("ckm_observable_residuals"),
        total_score=residuals.get("total_score"),
    )


def _seed_from_benchmark_module(seed: Any) -> QuarkFitSeed:
    """Normalize a benchmark seed object into a fit seed."""
    return QuarkFitSeed(
        up_singular_values=np.asarray(seed.up_singular_values, dtype=float),
        down_singular_values=np.asarray(seed.down_singular_values, dtype=float),
        overall_scale=float(seed.overall_scale),
        up_left=seed.up_left,
        up_right=seed.up_right,
        down_left=seed.down_left,
        down_right=seed.down_right,
    )


def _as_fit_seed(seed: Any) -> QuarkFitSeed:
    """Normalize a seed-like object into the local fit seed container."""
    if isinstance(seed, QuarkFitSeed):
        return seed
    return _seed_from_benchmark_module(seed)


def _canonicalize_fit_seed(seed: Any, *, overall_scale: float | None = None) -> QuarkFitSeed:
    """Return the canonical quotient-chart representative for a fit seed.

    The canonical representative absorbs the common scale into the physical
    singular spectra, wraps the left rotations into the fixed domain, and
    fixes the right rotations to the identity representative.
    """
    fit_seed = _as_fit_seed(seed)
    effective_scale = fit_seed.overall_scale if overall_scale is None else float(overall_scale)
    if not np.isfinite(effective_scale) or effective_scale <= 0.0:
        raise ValueError("overall_scale must be finite and strictly positive")
    up_physical = _require_positive_physical_spectrum(
        "up physical singular values",
        effective_scale * fit_seed.up_singular_values,
    )
    down_physical = _require_positive_physical_spectrum(
        "down physical singular values",
        effective_scale * fit_seed.down_singular_values,
    )
    return QuarkFitSeed(
        up_singular_values=up_physical,
        down_singular_values=down_physical,
        overall_scale=1.0,
        up_left=_canonicalize_rotation(fit_seed.up_left),
        up_right=RotationParameters(),
        down_left=_canonicalize_rotation(fit_seed.down_left),
        down_right=RotationParameters(),
    )


def _rotation_to_vector(rotation: RotationParameters) -> np.ndarray:
    return np.array(
        [rotation.theta12, rotation.theta13, rotation.theta23, rotation.delta],
        dtype=float,
    )


def _rotation_from_vector(values: np.ndarray) -> RotationParameters:
    return RotationParameters(*np.asarray(values, dtype=float))


def _require_positive_physical_spectrum(name: str, values: np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.shape != (3,):
        raise ValueError(f"{name} must have shape (3,)")
    if np.any(~np.isfinite(arr)) or np.any(arr <= 0.0):
        raise ValueError(f"{name} must be strictly positive")
    return arr


def _canonicalize_angle(angle: float) -> float:
    """Map an angle into the half-open interval [-pi, pi)."""
    return float(((float(angle) + np.pi) % (2.0 * np.pi)) - np.pi)


def _canonicalize_rotation(rotation: RotationParameters) -> RotationParameters:
    """Return the canonical representative for a left-rotation block.

    All four coordinates are wrapped into the same half-open interval so the
    quotient chart has no exact `2π` translation redundancy.
    """
    return RotationParameters(
        theta12=_canonicalize_angle(rotation.theta12),
        theta13=_canonicalize_angle(rotation.theta13),
        theta23=_canonicalize_angle(rotation.theta23),
        delta=_canonicalize_angle(rotation.delta),
    )


def _rotation_from_unitary(matrix: np.ndarray) -> RotationParameters:
    """Extract a canonical ``RotationParameters`` representative from a unitary."""
    U = np.asarray(matrix, dtype=np.complex128)
    if U.shape != (3, 3):
        raise ValueError("rotation unitary must have shape (3, 3)")
    theta13 = float(np.arcsin(np.clip(np.abs(U[0, 2]), 0.0, 1.0)))
    theta12 = float(np.arctan2(np.abs(U[0, 1]), max(np.abs(U[0, 0]), 1e-15)))
    theta23 = float(np.arctan2(np.abs(U[1, 2]), max(np.abs(U[2, 2]), 1e-15)))
    quartet = U[0, 0] * U[1, 2] * np.conjugate(U[0, 2]) * np.conjugate(U[2, 2])
    delta = float(_canonicalize_angle(np.angle(quartet)))
    return _canonicalize_rotation(
        RotationParameters(theta12=theta12, theta13=theta13, theta23=theta23, delta=delta)
    )


def _rotation_from_ckm_observables(observables: np.ndarray) -> RotationParameters:
    """Build a canonical CKM-like representative from observable magnitudes.

    The fit score only depends on the reduced CKM observable set, so the
    returned seed is reported in a fixed common-left gauge with the up-sector
    left rotation set to the identity representative and the down-sector left
    rotation rebuilt from the fitted observables.
    """
    obs = np.round(np.asarray(observables, dtype=float), 5)
    if obs.shape != (4,):
        raise ValueError("CKM observables must have shape (4,)")
    vus, vcb, vub, jarlskog = obs
    s13 = float(np.clip(abs(vub), 0.0, 1.0))
    c13 = float(np.sqrt(max(1.0 - s13 * s13, 0.0)))
    s12 = float(np.clip(abs(vus) / max(c13, 1e-15), 0.0, 1.0))
    s23 = float(np.clip(abs(vcb) / max(c13, 1e-15), 0.0, 1.0))
    c12 = float(np.sqrt(max(1.0 - s12 * s12, 0.0)))
    c23 = float(np.sqrt(max(1.0 - s23 * s23, 0.0)))
    denominator = max(c12 * c13 * c13 * c23 * s12 * s13 * s23, 1e-30)
    delta = float(np.arcsin(np.clip(jarlskog / denominator, -1.0, 1.0)))
    return _canonicalize_rotation(
        RotationParameters(
            theta12=float(np.arcsin(s12)),
            theta13=float(np.arcsin(s13)),
            theta23=float(np.arcsin(s23)),
            delta=delta,
        )
    )


def _sort_physical_spectrum(
    singular_values: np.ndarray,
    rotation: RotationParameters,
) -> tuple[np.ndarray, RotationParameters]:
    """Sort a physical spectrum and permute the matching left unitary."""
    spectrum = np.asarray(singular_values, dtype=float)
    if spectrum.shape != (3,):
        raise ValueError("physical singular values must have shape (3,)")
    order = np.argsort(spectrum, kind="mergesort")
    sorted_spectrum = spectrum[order]
    rotated = ckm_like_unitary(rotation)[:, order]
    return sorted_spectrum, _rotation_from_unitary(rotated)


def _canonicalize_reported_seed(
    seed: QuarkFitSeed,
    *,
    ckm_observables: np.ndarray,
) -> QuarkFitSeed:
    """Canonicalize a fitted seed into the frozen reported representative."""
    up_singular_values = np.sort(
        _require_positive_physical_spectrum("up physical singular values", seed.up_singular_values)
    )
    down_singular_values = np.sort(
        _require_positive_physical_spectrum("down physical singular values", seed.down_singular_values)
    )
    up_singular_values = np.round(up_singular_values, 5)
    down_singular_values = np.round(down_singular_values, 5)
    down_left = _rotation_from_ckm_observables(ckm_observables)
    return QuarkFitSeed(
        up_singular_values=up_singular_values,
        down_singular_values=down_singular_values,
        overall_scale=1.0,
        up_left=RotationParameters(),
        up_right=RotationParameters(),
        down_left=RotationParameters(
            theta12=float(np.round(down_left.theta12, 5)),
            theta13=float(np.round(down_left.theta13, 5)),
            theta23=float(np.round(down_left.theta23, 5)),
            delta=float(np.round(down_left.delta, 5)),
        ),
        down_right=RotationParameters(),
    )


def encode_quark_fit_canonical_vector(seed: Any) -> np.ndarray:
    """Encode the quotient-chart fit vector.

    Frozen order:
    log physical up singular values (light to heavy), log physical down
    singular values (light to heavy), then the left-rotation blocks
    ``up_left`` and ``down_left`` in ``(theta12, theta13, theta23, delta)``
    order.

    The common overall scale is absorbed into the singular values. Right
    rotations and the external ``r`` dial remain outside this chart.
    """
    fit_seed = _as_fit_seed(seed)
    up_physical = _require_positive_physical_spectrum(
        "up physical singular values",
        fit_seed.overall_scale * fit_seed.up_singular_values,
    )
    down_physical = _require_positive_physical_spectrum(
        "down physical singular values",
        fit_seed.overall_scale * fit_seed.down_singular_values,
    )
    values = [
        *np.log(up_physical),
        *np.log(down_physical),
        *_rotation_to_vector(_canonicalize_rotation(fit_seed.up_left)),
        *_rotation_to_vector(_canonicalize_rotation(fit_seed.down_left)),
    ]
    vector = np.asarray(values, dtype=float)
    if vector.shape != (QUARK_FIT_CANONICAL_VECTOR_SIZE,):
        raise AssertionError("unexpected quark canonical-vector shape")
    return vector


def decode_quark_fit_canonical_vector(vector: np.ndarray) -> QuarkFitSeed:
    """Decode the quotient-chart fit vector back into a fit seed.

    The decoded representative is normalized with ``overall_scale = 1`` and
    identity right rotations.
    """
    arr = np.asarray(vector, dtype=float)
    if arr.shape != (QUARK_FIT_CANONICAL_VECTOR_SIZE,):
        raise ValueError(
            f"canonical quark fit vector must have shape ({QUARK_FIT_CANONICAL_VECTOR_SIZE},)"
        )
    return QuarkFitSeed(
        up_singular_values=np.exp(arr[:3]),
        down_singular_values=np.exp(arr[3:6]),
        overall_scale=1.0,
        up_left=_canonicalize_rotation(_rotation_from_vector(arr[6:10])),
        up_right=RotationParameters(),
        down_left=_canonicalize_rotation(_rotation_from_vector(arr[10:14])),
        down_right=RotationParameters(),
    )


def encode_quark_fit_full_vector(seed: Any) -> np.ndarray:
    """Compatibility alias for the canonical quotient-chart vector."""
    return encode_quark_fit_canonical_vector(seed)


def decode_quark_fit_full_vector(
    vector: np.ndarray,
) -> QuarkFitSeed:
    """Compatibility alias for the canonical quotient-chart decoder."""
    return decode_quark_fit_canonical_vector(vector)


def _encode_seed(seed: QuarkFitSeed, fit_orientation: bool) -> np.ndarray:
    """Encode a fit seed into an optimization vector."""
    values = [*np.log(seed.up_singular_values), *np.log(seed.down_singular_values)]
    if fit_orientation:
        values.extend(
            [
                _canonicalize_rotation(seed.up_left).theta12,
                _canonicalize_rotation(seed.up_left).theta13,
                _canonicalize_rotation(seed.up_left).theta23,
                _canonicalize_rotation(seed.up_left).delta,
                _canonicalize_rotation(seed.down_left).theta12,
                _canonicalize_rotation(seed.down_left).theta13,
                _canonicalize_rotation(seed.down_left).theta23,
                _canonicalize_rotation(seed.down_left).delta,
            ]
        )
    return np.asarray(values, dtype=float)


def _decode_seed(
    vector: np.ndarray,
    template: QuarkFitSeed,
    fit_orientation: bool,
) -> QuarkFitSeed:
    """Decode an optimization vector back into a fit seed."""
    arr = np.asarray(vector, dtype=float)
    up_singular_values = np.exp(arr[:3])
    down_singular_values = np.exp(arr[3:6])
    if fit_orientation:
        up_left = _canonicalize_rotation(RotationParameters(*arr[6:10]))
        down_left = _canonicalize_rotation(RotationParameters(*arr[10:14]))
    else:
        up_left = _canonicalize_rotation(template.up_left)
        down_left = _canonicalize_rotation(template.down_left)
    return QuarkFitSeed(
        up_singular_values=up_singular_values,
        down_singular_values=down_singular_values,
        overall_scale=1.0,
        up_left=up_left,
        up_right=RotationParameters(),
        down_left=down_left,
        down_right=RotationParameters(),
    )


def _seed_to_point(
    seed: QuarkFitSeed,
    *,
    r: float,
    Lambda_IR: float,
    k: float,
) -> QuarkSpurionPoint:
    """Build a quark MFV point from a fit seed."""
    return build_mfv_point_from_singular_values(
        up_singular_values=seed.up_singular_values,
        down_singular_values=seed.down_singular_values,
        overall_scale=seed.overall_scale,
        r=r,
        up_left=seed.up_left,
        up_right=seed.up_right,
        down_left=seed.down_left,
        down_right=seed.down_right,
        Lambda_IR=Lambda_IR,
        k=k,
        label="fitted_mfv_point",
    )


def fit_quark_sector(
    targets: QuarkTargets,
    *,
    r: float = 0.25,
    overall_scale: float | None = None,
    seed: Any | None = None,
    k: float = 1.2209e19,
    Lambda_IR: float = 3000.0,
    bulk_map: Optional[BulkMassMap] = None,
    max_nfev: int = 200,
    fit_orientation: bool = True,
) -> QuarkFitSolution:
    """Fit a compact spurion seed to target masses and CKM data.

    The optimizer state is the canonical quotient chart:
    absorbed-scale physical singular values plus left-rotation coordinates.
    When ``fit_orientation`` is false, the left rotations are held fixed at the
    canonicalized template representative and only the physical spectra are
    searched.
    """
    if seed is None:
        from .benchmarks import default_spurion_seed

        seed = default_spurion_seed()
    template = _canonicalize_fit_seed(seed, overall_scale=overall_scale)

    initial_point = _seed_to_point(template, r=r, Lambda_IR=Lambda_IR, k=k)
    initial_result = evaluate_quark_fit(initial_point, targets, bulk_mass_map=bulk_map)

    def residual_vector(vector: np.ndarray) -> np.ndarray:
        try:
            candidate = _decode_seed(vector, template, fit_orientation)
            point = _seed_to_point(candidate, r=r, Lambda_IR=Lambda_IR, k=k)
            result = evaluate_quark_fit(point, targets, bulk_mass_map=bulk_map)
            return np.concatenate(
                [result.mass_residuals_up, result.mass_residuals_down, result.ckm_residuals]
            )
        except (ValueError, np.linalg.LinAlgError):
            # Return large residuals for unphysical parameter regions
            # (e.g., negative bulk-mass eigenvalues at extreme r/scale)
            return np.full(10, 1e6)

    optimum = least_squares(
        residual_vector,
        x0=_encode_seed(template, fit_orientation),
        max_nfev=max_nfev,
    )
    raw_best_seed = _decode_seed(optimum.x, template, fit_orientation)
    best_point = _seed_to_point(raw_best_seed, r=r, Lambda_IR=Lambda_IR, k=k)
    best_result = evaluate_quark_fit(best_point, targets, bulk_mass_map=bulk_map)
    best_seed = _canonicalize_reported_seed(raw_best_seed, ckm_observables=best_result.ckm_observables)
    return QuarkFitSolution(
        seed=best_seed,
        result=best_result,
        success=bool(optimum.success),
        message=str(optimum.message),
        nfev=int(optimum.nfev),
        initial_score=initial_result.score,
    )
