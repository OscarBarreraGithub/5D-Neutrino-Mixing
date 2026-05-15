"""Deterministic benchmarks and fit seeds for the quark-sector module."""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field

import numpy as np

from qcd.constants import M_TOP_MS

from .fit import (
    QuarkFitResult,
    QuarkFitSolution,
    QuarkTargets,
    evaluate_quark_fit,
    fit_quark_sector,
)
from .model import (
    QuarkSpurionPoint,
    RotationParameters,
    build_mfv_point_from_singular_values,
    ckm_like_unitary,
)
from .pdg_quark_masses import (
    PDG_QUARK_MASSES_EDITION,
    pdg_2sigma_relative_at_scale,
    pdg_up_down_arrays_at_scale,
)
from .scales import DEFAULT_QUARK_FIT_SCALE_GEV, DEFAULT_QUARK_TARGET_SCALE_GEV

_TARGET_CKM = ckm_like_unitary(
    RotationParameters(theta12=0.2274, theta13=0.00368, theta23=0.0415, delta=1.196)
)
# Legacy v1 fixed-scale target bundle (mu = 3 TeV; ad-hoc round numbers).
# Retained to support the legacy ``rough_sm_targets`` compatibility path and
# the synthetic regression test fixtures. Production code now consumes
# ``_FIXED_SCALE_TARGETS_PDG2024_MT_V1`` below.
_FIXED_SCALE_TARGETS_MU_3TEV_V1 = {
    "scale_GeV": DEFAULT_QUARK_TARGET_SCALE_GEV,
    "provenance": "legacy v1 fixed-scale target bundle (mu = 3 TeV, ad-hoc rounds)",
    "up_masses": np.array([0.0013, 0.62, 172.0], dtype=float),
    "down_masses": np.array([0.0028, 0.057, 2.86], dtype=float),
    "ckm": _TARGET_CKM,
}


def _build_pdg2024_mt_target_bundle() -> dict[str, object]:
    """Build the PDG 2024 MS-bar target bundle at ``mu = m_t(m_t)``.

    The light/charm/bottom inputs are RG-evolved from their PDG 2024
    reference scales to ``mu_common = qcd.constants.M_TOP_MS = 163.5 GeV``
    using :func:`qcd.mass_running.run_msbar_mass`. Top is run from
    ``m_t(m_t) = 162.5`` to 163.5 GeV; the effect is per-mille.

    The Wilson-coefficient ``alpha_s`` reference scale is **independent**
    of this mass-target scale; see ``modern/inputs.py``.
    """
    up, down = pdg_up_down_arrays_at_scale(DEFAULT_QUARK_FIT_SCALE_GEV)
    rel_2sigma = pdg_2sigma_relative_at_scale(DEFAULT_QUARK_FIT_SCALE_GEV)
    return {
        "scale_GeV": DEFAULT_QUARK_FIT_SCALE_GEV,
        "provenance": (
            "PDG 2024 MS-bar quark masses RG-evolved to mu = m_t(m_t) = "
            f"{DEFAULT_QUARK_FIT_SCALE_GEV} GeV via qcd.mass_running."
        ),
        "edition": PDG_QUARK_MASSES_EDITION,
        "up_masses": up,
        "down_masses": down,
        "ckm": _TARGET_CKM,
        "up_2sigma_relative": np.array(
            [rel_2sigma["u"], rel_2sigma["c"], rel_2sigma["t"]],
            dtype=float,
        ),
        "down_2sigma_relative": np.array(
            [rel_2sigma["d"], rel_2sigma["s"], rel_2sigma["b"]],
            dtype=float,
        ),
        "label": "pdg-2024-msbar-mu-mt-v1",
    }


# Repo-owned PDG-2024 fixed-scale target bundle at mu_common = m_t(m_t).
# Computed at import time so the values are deterministic and unit-tested.
_FIXED_SCALE_TARGETS_PDG2024_MT_V1 = _build_pdg2024_mt_target_bundle()


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
    """Return the default PDG-2024 quark target table at ``mu = m_t(m_t)``.

    The mass targets are PDG 2024 MS-bar values RG-evolved to
    ``DEFAULT_QUARK_FIT_SCALE_GEV = qcd.constants.M_TOP_MS = 163.5 GeV``;
    the CKM target is the PDG 2024 central unitary. The Wilson-coefficient
    ``alpha_s`` reference scale stays at 3 TeV (see ``modern/inputs.py``);
    the two scales are kept orthogonal.
    """
    return QuarkTargets(
        up_masses=_FIXED_SCALE_TARGETS_PDG2024_MT_V1["up_masses"],
        down_masses=_FIXED_SCALE_TARGETS_PDG2024_MT_V1["down_masses"],
        ckm=_FIXED_SCALE_TARGETS_PDG2024_MT_V1["ckm"],
        label="pdg-2024-msbar-mu-mt-v1",
    )


def rough_sm_targets() -> QuarkTargets:
    """Compatibility helper for the old rough-SM-like target nomenclature.

    Deprecated: prefer ``default_quark_targets()``. The shape and ordering
    are preserved.
    """
    warnings.warn(
        "rough_sm_targets() is deprecated; use default_quark_targets() "
        "for the PDG-2024 MS-bar target bundle at mu = m_t(m_t).",
        DeprecationWarning,
        stacklevel=2,
    )
    targets = default_quark_targets()
    return QuarkTargets(
        up_masses=targets.up_masses,
        down_masses=targets.down_masses,
        ckm=targets.ckm,
        label="rough-sm-like-compatibility",
    )


def default_spurion_seed() -> SpurionSeed:
    """Return the deterministic seed used by fits, scans, and plots."""
    return SpurionSeed(
        up_singular_values=np.array(
            [
                0.1434280953329283,
                0.33368791612279205,
                1.2394579588587887,
            ],
            dtype=float,
        ),
        down_singular_values=np.array(
            [
                0.22903978232967957,
                0.16350920983304943,
                0.28800013126033523,
            ],
            dtype=float,
        ),
        overall_scale=2.8,
        up_left=RotationParameters(
            theta12=-2.4040939384155906,
            theta13=-2.8383004089408406,
            theta23=-0.20007460664685306,
            delta=2.05348624211683,
        ),
        up_right=RotationParameters(),
        down_left=RotationParameters(
            theta12=3.1151255159110542,
            theta13=2.755960572425262,
            theta23=0.287658452173571,
            delta=-0.2529736730248784,
        ),
        down_right=RotationParameters(),
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
        metadata={
            "preferred_r_window": (0.1, 0.4),
            "overall_y": seed.overall_scale,
            "default_target_scale_GeV": DEFAULT_QUARK_FIT_SCALE_GEV,
            "default_target_label": default_quark_targets().label,
            "kk_scale_convention": "repo-default-xi_KK=1.0-so-M_KK-equals-Lambda_IR",
        },
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
            "Use default_quark_targets() for the PDG-2024 MS-bar target "
            "bundle RG-evolved to mu = m_t(m_t) = 163.5 GeV."
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
    """Solve the fixed-scale fit starting from the deterministic seed."""
    seed = default_spurion_seed()
    scale = seed.overall_scale if overall_scale is None else overall_scale
    return fit_quark_sector(
        default_quark_targets(),
        r=r,
        overall_scale=scale,
        seed=seed,
        max_nfev=max_nfev,
    )
