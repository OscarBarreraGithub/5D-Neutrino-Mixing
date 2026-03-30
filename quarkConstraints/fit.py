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
                np.concatenate([self.mass_residuals_up, self.mass_residuals_down, self.ckm_residuals])
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
    """Compute mass and CKM residuals against targets."""
    up_log_residuals = np.log(np.maximum(masses_up, 1e-30) / targets.up_masses)
    down_log_residuals = np.log(np.maximum(masses_down, 1e-30) / targets.down_masses)
    ckm_abs_residuals = np.abs(np.abs(ckm) - targets.abs_ckm)
    ckm_obs_residuals = (
        ckm_observables(ckm) - targets.ckm_observables
    ) / np.maximum(targets.ckm_observables, 1e-8)

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


def _encode_seed(seed: QuarkFitSeed, fit_orientation: bool) -> np.ndarray:
    """Encode a fit seed into an optimization vector."""
    values = [*np.log(seed.up_singular_values), *np.log(seed.down_singular_values)]
    if fit_orientation:
        values.extend(
            [
                seed.up_left.theta12,
                seed.up_left.theta13,
                seed.up_left.theta23,
                seed.up_left.delta,
                seed.down_left.theta12,
                seed.down_left.theta13,
                seed.down_left.theta23,
                seed.down_left.delta,
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
        up_left = RotationParameters(*arr[6:10])
        down_left = RotationParameters(*arr[10:14])
    else:
        up_left = template.up_left
        down_left = template.down_left
    return QuarkFitSeed(
        up_singular_values=up_singular_values,
        down_singular_values=down_singular_values,
        overall_scale=template.overall_scale,
        up_left=up_left,
        up_right=template.up_right,
        down_left=down_left,
        down_right=template.down_right,
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
    """Fit a compact spurion seed to target masses and CKM data."""
    if seed is None:
        from .benchmarks import default_spurion_seed

        seed = default_spurion_seed()
    template = _seed_from_benchmark_module(seed)
    if overall_scale is not None and overall_scale != template.overall_scale:
        template = QuarkFitSeed(
            up_singular_values=template.up_singular_values,
            down_singular_values=template.down_singular_values,
            overall_scale=float(overall_scale),
            up_left=template.up_left,
            up_right=template.up_right,
            down_left=template.down_left,
            down_right=template.down_right,
        )

    initial_point = _seed_to_point(template, r=r, Lambda_IR=Lambda_IR, k=k)
    initial_result = evaluate_quark_fit(initial_point, targets, bulk_mass_map=bulk_map)

    def residual_vector(vector: np.ndarray) -> np.ndarray:
        candidate = _decode_seed(vector, template, fit_orientation)
        point = _seed_to_point(candidate, r=r, Lambda_IR=Lambda_IR, k=k)
        result = evaluate_quark_fit(point, targets, bulk_mass_map=bulk_map)
        return np.concatenate(
            [result.mass_residuals_up, result.mass_residuals_down, result.ckm_residuals]
        )

    optimum = least_squares(
        residual_vector,
        x0=_encode_seed(template, fit_orientation),
        max_nfev=max_nfev,
    )
    best_seed = _decode_seed(optimum.x, template, fit_orientation)
    best_point = _seed_to_point(best_seed, r=r, Lambda_IR=Lambda_IR, k=k)
    best_result = evaluate_quark_fit(best_point, targets, bulk_mass_map=bulk_map)
    return QuarkFitSolution(
        seed=best_seed,
        result=best_result,
        success=bool(optimum.success),
        message=str(optimum.message),
        nfev=int(optimum.nfev),
        initial_score=initial_result.score,
    )
