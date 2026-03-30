"""Deterministic benchmarks and fit seeds for the quark-sector module."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .fit import QuarkFitResult, QuarkFitSolution, QuarkTargets, evaluate_quark_fit, fit_quark_sector
from .model import (
    RotationParameters,
    QuarkSpurionPoint,
    build_mfv_point_from_singular_values,
    ckm_like_unitary,
)


def _as_real_vector(name: str, values: np.ndarray | list[float] | tuple[float, ...]) -> np.ndarray:
    """Return a validated positive real three-vector."""
    arr = np.asarray(values, dtype=float)
    if arr.shape != (3,):
        raise ValueError(f"{name} must have shape (3,)")
    if np.any(arr <= 0.0):
        raise ValueError(f"{name} must be strictly positive")
    return arr


@dataclass(frozen=True)
class SpurionSeed:
    """Compact fit seed in MFV spurion space."""

    up_singular_values: np.ndarray
    down_singular_values: np.ndarray
    overall_scale: float = 3.0
    up_left: RotationParameters = field(default_factory=RotationParameters)
    up_right: RotationParameters = field(default_factory=RotationParameters)
    down_left: RotationParameters = field(default_factory=RotationParameters)
    down_right: RotationParameters = field(default_factory=RotationParameters)

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "up_singular_values",
            _as_real_vector("up_singular_values", self.up_singular_values),
        )
        object.__setattr__(
            self,
            "down_singular_values",
            _as_real_vector("down_singular_values", self.down_singular_values),
        )
        if self.overall_scale <= 0.0:
            raise ValueError("overall_scale must be positive")


@dataclass(frozen=True)
class QuarkBenchmark:
    """Repo-local deterministic benchmark."""

    name: str
    point: QuarkSpurionPoint
    targets: QuarkTargets
    notes: str


@dataclass(frozen=True)
class EvaluatedBenchmark:
    """Benchmark point together with its evaluated result."""

    benchmark: QuarkBenchmark
    result: QuarkFitResult

    @property
    def point(self) -> QuarkSpurionPoint:
        return self.benchmark.point

    @property
    def targets(self) -> QuarkTargets:
        return self.benchmark.targets


def default_quark_targets() -> QuarkTargets:
    """Return a rough SM-like target set for exploratory fits."""
    ckm = ckm_like_unitary(
        RotationParameters(theta12=0.2274, theta13=0.00368, theta23=0.0415, delta=1.196)
    )
    return QuarkTargets(
        up_masses=np.array([0.0013, 0.62, 172.0], dtype=float),
        down_masses=np.array([0.0028, 0.057, 2.86], dtype=float),
        ckm=ckm,
        label="rough-sm-like",
    )


def rough_sm_targets() -> QuarkTargets:
    """Backward-compatible alias for the default exploratory targets."""
    return default_quark_targets()


def default_spurion_seed() -> SpurionSeed:
    """Return the deterministic seed used by fits, scans, and plots."""
    return SpurionSeed(
        up_singular_values=np.array([0.003, 0.12, 1.0], dtype=float),
        down_singular_values=np.array([0.006, 0.035, 0.45], dtype=float),
        overall_scale=2.8,
        up_left=RotationParameters(theta12=0.02, theta13=0.002, theta23=0.03, delta=0.8),
        up_right=RotationParameters(theta12=0.04, theta13=0.01, theta23=0.02, delta=0.4),
        down_left=RotationParameters(theta12=0.23, theta13=0.0035, theta23=0.041, delta=1.2),
        down_right=RotationParameters(theta12=0.06, theta13=0.015, theta23=0.03, delta=0.5),
    )


def benchmark_spurion_input(
    *,
    r: float = 0.25,
    Lambda_IR: float = 3000.0,
    seed: SpurionSeed | None = None,
) -> QuarkSpurionPoint:
    """Build the default benchmark spurion point with an explicit `r` dial."""
    seed = default_spurion_seed() if seed is None else seed
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
        label="repo-local-mfv-benchmark",
        metadata={"preferred_r_window": (0.1, 0.4), "overall_y": seed.overall_scale},
    )


def default_quark_benchmark() -> QuarkBenchmark:
    """Return the repo-local regression benchmark."""
    point = benchmark_spurion_input()
    fit = evaluate_quark_fit(point)
    targets = QuarkTargets(
        up_masses=fit.masses_up,
        down_masses=fit.masses_down,
        ckm=fit.ckm,
        label="repo-local-mfv-benchmark",
    )
    return QuarkBenchmark(
        name="repo-local-mfv-benchmark",
        point=point,
        targets=targets,
        notes=(
            "Deterministic repo-local regression point. "
            "Use default_quark_targets() when fitting to rough SM-like data."
        ),
    )


def evaluate_default_benchmark() -> EvaluatedBenchmark:
    """Evaluate the stored deterministic benchmark against its own targets."""
    benchmark = default_quark_benchmark()
    result = evaluate_quark_fit(benchmark.point, targets=benchmark.targets)
    return EvaluatedBenchmark(benchmark=benchmark, result=result)


def solve_default_benchmark(
    *,
    r: float = 0.25,
    overall_scale: float | None = None,
    max_nfev: int = 150,
) -> QuarkFitSolution:
    """Solve the rough SM-like fit starting from the deterministic seed."""
    seed = default_spurion_seed()
    scale = seed.overall_scale if overall_scale is None else overall_scale
    return fit_quark_sector(
        default_quark_targets(),
        r=r,
        overall_scale=scale,
        seed=seed,
        max_nfev=max_nfev,
    )
