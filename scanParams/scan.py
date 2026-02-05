"""Grid-scan driver for the RS lepton-sector parameter space.

This scanner evaluates ``compute_all_yukawas()`` over configurable grids,
applies perturbativity / naturalness / mu-to-e-gamma filters, and writes
results to CSV for later inspection.

Conventions implemented here follow ``scanParams/THEORY_PRIORS.md``:
- explicit ``k`` input (fixed by default),
- ``M_N`` represented by ``MN_over_k = M_N / k``,
- explicit ``M_KK = xi_KK * Lambda_IR`` mapping for LFV.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple
import csv
import itertools
import subprocess
import time

import numpy as np

from flavorConstraints import (
    PREFAC_BR,
    check_mu_to_e_gamma,
    coefficient_from_br_limit,
)
from warpConfig.baseParams import MPL
from warpConfig.wavefuncs import f_IR
from yukawa import YukawaResult, compute_all_yukawas
from .anarchy import AnarchyConfig, sample_anarchy_state


# Published MEG II bound (2025): BR(mu -> e gamma) < 1.5e-13 (90% CL).
BR_LIMIT_MEGII_2025 = 1.5e-13


# ── CSV schema ───────────────────────────────────────────────────────────

CSV_COLUMNS = [
    "sample_index",
    "git_commit",
    "dirty_tree",
    "rng_seed_global",
    "rng_seed_sample",
    "lfv_model",
    "br_limit",
    "prefac_br",
    "lfv_C",
    "lfv_reference_scale",
    "xi_KK",
    "max_fL_ratio",
    "delta_cL_max_symmetric",
    "delta_cL_max_one_sided",
    "delta_cL_max_symmetric_over_cL_pct",
    "delta_cL_max_one_sided_over_cL_pct",
    "k",
    "MN_mode",
    "MN_over_k",
    "M_N",
    "Lambda_IR",
    "M_KK",
    "lightest_nu_mass",
    "ordering",
    "majorana_alpha",
    "majorana_beta",
    "c_L",
    "c_N",
    "c_E1",
    "c_E2",
    "c_E3",
    "Y_E_bar_1",
    "Y_E_bar_2",
    "Y_E_bar_3",
    "Y_N_bar_1",
    "Y_N_bar_2",
    "Y_N_bar_3",
    "f_L",
    "f_N",
    "f_N_UV",
    "max_Y_bar",
    "perturbative",
    "natural",
    "lfv_passes",
    "lfv_lhs",
    "lfv_rhs",
    "lfv_ratio",
    "passes_all",
    "reject_reason",
    "anarchy_enabled",
    "anarchy_score",
    "anarchy_band_penalty",
    "anarchy_condition_penalty",
    "anarchy_w_band",
    "anarchy_w_cond",
    "anarchy_w_fit",
    "anarchy_yN_overall",
]


# ── Helpers ───────────────────────────────────────────────────────────────


def _as_1d_float_array(name: str, values: Sequence[float] | np.ndarray) -> np.ndarray:
    """Convert input into a non-empty 1D float array."""
    arr = np.asarray(values, dtype=float)
    if arr.ndim == 0:
        arr = arr.reshape(1)
    if arr.ndim != 1:
        raise ValueError(f"{name} must be 1D")
    if arr.size == 0:
        raise ValueError(f"{name} must not be empty")
    return arr


def _run_git(args: List[str], repo_root: Path) -> str:
    """Run a git command and return stripped stdout (empty string on failure)."""
    proc = subprocess.run(
        ["git", *args],
        cwd=str(repo_root),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        return ""
    return proc.stdout.strip()


def _resolve_git_metadata(enabled: bool) -> Tuple[str, Optional[bool]]:
    """Return (git_commit, dirty_tree) metadata."""
    if not enabled:
        return "disabled", None

    repo_root = Path(__file__).resolve().parents[1]
    git_commit = _run_git(["rev-parse", "HEAD"], repo_root)
    if not git_commit:
        return "unknown", None

    status = _run_git(["status", "--porcelain"], repo_root)
    dirty_tree = len(status) > 0
    return git_commit, dirty_tree


def _sample_seed(global_seed: Optional[int], sample_index: int) -> Optional[int]:
    """Derive a deterministic per-sample seed from a global seed."""
    if global_seed is None:
        return None
    seed_seq = np.random.SeedSequence([int(global_seed), int(sample_index)])
    return int(seed_seq.generate_state(1, dtype=np.uint64)[0])


def _ratio_from_delta(c0: float, delta: float, epsilon: float, mode: str) -> float:
    """Return f_IR ratio at a given delta for a selected degeneracy mode."""
    if mode == "symmetric":
        num = float(f_IR(c0 - delta, epsilon))
        den = float(f_IR(c0 + delta, epsilon))
        return num / den
    if mode == "one_sided":
        num = float(f_IR(c0 - delta, epsilon))
        den = float(f_IR(c0, epsilon))
        return num / den
    raise ValueError(f"Unknown mode: {mode}")


def _solve_delta_c_for_ratio(
    c0: float,
    epsilon: float,
    max_fL_ratio: float,
    mode: str,
) -> float:
    """Solve for delta c such that the selected ratio equals max_fL_ratio."""
    # Keep search away from c -> 0.5 singular neighborhood.
    hi = min(0.05, max(1e-6, c0 - 0.500001))
    lo = 0.0

    # Expand bracket if needed (should be rare in scanner ranges).
    ratio_hi = _ratio_from_delta(c0, hi, epsilon, mode)
    if ratio_hi < max_fL_ratio:
        for _ in range(8):
            candidate = hi * 2.0
            if candidate >= 0.2:
                break
            hi = candidate
            ratio_hi = _ratio_from_delta(c0, hi, epsilon, mode)
            if ratio_hi >= max_fL_ratio:
                break
        if ratio_hi < max_fL_ratio:
            return hi

    for _ in range(80):
        mid = 0.5 * (lo + hi)
        ratio_mid = _ratio_from_delta(c0, mid, epsilon, mode)
        if ratio_mid > max_fL_ratio:
            hi = mid
        else:
            lo = mid
    return lo


def _derive_cL_degeneracy_metadata(c_L0: float, k: float, Lambda_IR: float, max_fL_ratio: float) -> Dict[str, float]:
    """Compute derived delta-c tolerances from a wavefunction-ratio prior."""
    epsilon = float(Lambda_IR / k)
    d_sym = _solve_delta_c_for_ratio(c_L0, epsilon, max_fL_ratio=max_fL_ratio, mode="symmetric")
    d_one = _solve_delta_c_for_ratio(c_L0, epsilon, max_fL_ratio=max_fL_ratio, mode="one_sided")
    return {
        "delta_cL_max_symmetric": d_sym,
        "delta_cL_max_one_sided": d_one,
        "delta_cL_max_symmetric_over_cL_pct": 100.0 * d_sym / c_L0 if c_L0 != 0 else np.nan,
        "delta_cL_max_one_sided_over_cL_pct": 100.0 * d_one / c_L0 if c_L0 != 0 else np.nan,
    }


# ── Configuration ────────────────────────────────────────────────────────


@dataclass
class ScanConfig:
    """Configuration for an RS lepton-sector parameter scan.

    Notes
    -----
    - ``k`` is explicit and fixed by default.
    - ``M_N`` is derived from ``MN_over_k``.
    - ``M_KK`` used in LFV checks is ``xi_KK * Lambda_IR``.
    """

    # Geometry / fixed physics
    k: float = MPL
    ordering: str = "normal"
    majorana_alpha: float = 0.0
    majorana_beta: float = 0.0

    # Scan grids
    Lambda_IR_values: np.ndarray = field(
        default_factory=lambda: np.array([3000.0], dtype=float)
    )
    c_L_values: np.ndarray = field(
        default_factory=lambda: np.linspace(0.52, 0.70, 10)
    )
    c_N_values: np.ndarray = field(
        default_factory=lambda: np.linspace(0.15, 0.45, 10)
    )
    lightest_nu_mass_values: np.ndarray = field(
        default_factory=lambda: np.array([0.002], dtype=float)
    )

    # Charged-lepton bulk masses
    c_E_fixed: Optional[List[float]] = field(default_factory=lambda: [0.75, 0.60, 0.50])
    c_E_grid: Optional[List[np.ndarray]] = None
    sort_c_E_descending: bool = True

    # UV Majorana control
    MN_mode: str = "fixed_ratio"  # "fixed_ratio" or "scan_ratio"
    MN_over_k: float = 0.1
    MN_over_k_values: Optional[np.ndarray] = None

    # LFV convention
    lfv_model: str = "mu_to_e_gamma_nda_v1"
    br_limit: float = BR_LIMIT_MEGII_2025
    prefac_br: float = PREFAC_BR
    lfv_reference_scale: float = 3000.0
    xi_KK: float = 1.0

    # Filter thresholds
    max_Y_bar: float = 4.0
    naturalness_range: Tuple[float, float] = (0.1, 4.0)
    max_fL_ratio: float = 1.1

    # Metadata controls
    rng_seed_global: Optional[int] = None
    record_git_metadata: bool = True

    # Optional anarchic-prior sampling/scoring
    anarchy: Optional[AnarchyConfig] = None
    anarchy_min_score: Optional[float] = None

    # Internal caches (built in __post_init__)
    _c_E_points: List[Tuple[float, float, float]] = field(init=False, repr=False)
    _MN_over_k_points: np.ndarray = field(init=False, repr=False)

    def __post_init__(self):
        if self.ordering.strip().lower() != "normal":
            raise ValueError("v1 scanner currently supports only ordering='normal'")

        if self.k <= 0:
            raise ValueError("k must be positive")
        if self.xi_KK <= 0:
            raise ValueError("xi_KK must be positive")
        if self.br_limit <= 0:
            raise ValueError("br_limit must be positive")
        if self.prefac_br <= 0:
            raise ValueError("prefac_br must be positive")
        if self.max_fL_ratio <= 1.0:
            raise ValueError("max_fL_ratio must be > 1")

        self.Lambda_IR_values = _as_1d_float_array("Lambda_IR_values", self.Lambda_IR_values)
        self.c_L_values = _as_1d_float_array("c_L_values", self.c_L_values)
        self.c_N_values = _as_1d_float_array("c_N_values", self.c_N_values)
        self.lightest_nu_mass_values = _as_1d_float_array(
            "lightest_nu_mass_values", self.lightest_nu_mass_values
        )

        if np.any(self.Lambda_IR_values <= 0):
            raise ValueError("Lambda_IR_values must be positive")
        if np.any(self.lightest_nu_mass_values < 0):
            raise ValueError("lightest_nu_mass_values must be non-negative")
        if self.anarchy_min_score is not None and not np.isfinite(self.anarchy_min_score):
            raise ValueError("anarchy_min_score must be finite when provided")

        if self.MN_mode not in {"fixed_ratio", "scan_ratio"}:
            raise ValueError("MN_mode must be 'fixed_ratio' or 'scan_ratio'")

        if self.MN_mode == "fixed_ratio":
            if self.MN_over_k <= 0:
                raise ValueError("MN_over_k must be positive for MN_mode='fixed_ratio'")
            self._MN_over_k_points = np.array([float(self.MN_over_k)], dtype=float)
        else:
            if self.MN_over_k_values is None:
                raise ValueError("MN_over_k_values must be provided for MN_mode='scan_ratio'")
            self._MN_over_k_points = _as_1d_float_array(
                "MN_over_k_values", self.MN_over_k_values
            )
            if np.any(self._MN_over_k_points <= 0):
                raise ValueError("MN_over_k_values must be positive")

        # Build charged-lepton c_E points.
        if self.c_E_grid is not None:
            if len(self.c_E_grid) != 3:
                raise ValueError(f"c_E_grid must contain exactly 3 arrays, got {len(self.c_E_grid)}")
            c1 = _as_1d_float_array("c_E_grid[0]", self.c_E_grid[0])
            c2 = _as_1d_float_array("c_E_grid[1]", self.c_E_grid[1])
            c3 = _as_1d_float_array("c_E_grid[2]", self.c_E_grid[2])
            raw_points = itertools.product(c1, c2, c3)
            points = [tuple(map(float, p)) for p in raw_points]
        else:
            if self.c_E_fixed is None:
                raise ValueError("Either c_E_fixed or c_E_grid must be provided")
            if len(self.c_E_fixed) != 3:
                raise ValueError(f"c_E_fixed must contain exactly 3 values, got {len(self.c_E_fixed)}")
            points = [tuple(map(float, self.c_E_fixed))]

        if self.sort_c_E_descending:
            self._c_E_points = [tuple(sorted(p, reverse=True)) for p in points]
        else:
            self._c_E_points = points

    @property
    def c_E_points(self) -> List[Tuple[float, float, float]]:
        """Materialized scan points for (c_E1, c_E2, c_E3)."""
        return self._c_E_points

    @property
    def MN_over_k_points(self) -> np.ndarray:
        """Materialized scan points for MN_over_k."""
        return self._MN_over_k_points

    @property
    def total_points(self) -> int:
        """Total number of scan points."""
        return (
            len(self.Lambda_IR_values)
            * len(self.c_L_values)
            * len(self.c_N_values)
            * len(self.c_E_points)
            * len(self.MN_over_k_points)
            * len(self.lightest_nu_mass_values)
        )


# ── Single-point evaluation ──────────────────────────────────────────────


def _evaluate_point(
    sample_index: int,
    Lambda_IR: float,
    c_L: float,
    c_N: float,
    c_E: Tuple[float, float, float],
    MN_over_k: float,
    lightest_nu_mass: float,
    config: ScanConfig,
    lfv_C: float,
    git_commit: str,
    dirty_tree: Optional[bool],
    rng_seed_sample: Optional[int],
    extra_filters: List[Callable[[YukawaResult], Tuple[bool, str]]],
) -> Dict[str, Any]:
    """Evaluate one parameter point and return a row dict."""
    M_N = float(MN_over_k * config.k)
    M_KK = float(config.xi_KK * Lambda_IR)
    cL_degeneracy = _derive_cL_degeneracy_metadata(
        c_L0=c_L,
        k=config.k,
        Lambda_IR=Lambda_IR,
        max_fL_ratio=config.max_fL_ratio,
    )

    row: Dict[str, Any] = {
        "sample_index": sample_index,
        "git_commit": git_commit,
        "dirty_tree": dirty_tree,
        "rng_seed_global": config.rng_seed_global,
        "rng_seed_sample": rng_seed_sample,
        "lfv_model": config.lfv_model,
        "br_limit": config.br_limit,
        "prefac_br": config.prefac_br,
        "lfv_C": lfv_C,
        "lfv_reference_scale": config.lfv_reference_scale,
        "xi_KK": config.xi_KK,
        "max_fL_ratio": config.max_fL_ratio,
        "delta_cL_max_symmetric": cL_degeneracy["delta_cL_max_symmetric"],
        "delta_cL_max_one_sided": cL_degeneracy["delta_cL_max_one_sided"],
        "delta_cL_max_symmetric_over_cL_pct": cL_degeneracy["delta_cL_max_symmetric_over_cL_pct"],
        "delta_cL_max_one_sided_over_cL_pct": cL_degeneracy["delta_cL_max_one_sided_over_cL_pct"],
        "k": config.k,
        "MN_mode": config.MN_mode,
        "MN_over_k": MN_over_k,
        "M_N": M_N,
        "Lambda_IR": Lambda_IR,
        "M_KK": M_KK,
        "lightest_nu_mass": lightest_nu_mass,
        "ordering": config.ordering,
        "majorana_alpha": config.majorana_alpha,
        "majorana_beta": config.majorana_beta,
        "c_L": c_L,
        "c_N": c_N,
        "c_E1": c_E[0],
        "c_E2": c_E[1],
        "c_E3": c_E[2],
    }

    row.update(
        {
            "Y_E_bar_1": np.nan,
            "Y_E_bar_2": np.nan,
            "Y_E_bar_3": np.nan,
            "Y_N_bar_1": np.nan,
            "Y_N_bar_2": np.nan,
            "Y_N_bar_3": np.nan,
            "f_L": np.nan,
            "f_N": np.nan,
            "f_N_UV": np.nan,
            "max_Y_bar": np.nan,
            "perturbative": False,
            "natural": False,
            "lfv_passes": False,
            "lfv_lhs": np.nan,
            "lfv_rhs": np.nan,
            "lfv_ratio": np.nan,
            "passes_all": False,
            "reject_reason": "",
            "anarchy_enabled": config.anarchy is not None,
            "anarchy_score": np.nan,
            "anarchy_band_penalty": np.nan,
            "anarchy_condition_penalty": np.nan,
            "anarchy_w_band": np.nan,
            "anarchy_w_cond": np.nan,
            "anarchy_w_fit": np.nan,
            "anarchy_yN_overall": np.nan,
        }
    )

    try:
        result = compute_all_yukawas(
            Lambda_IR=Lambda_IR,
            c_L=c_L,
            c_E=c_E,
            c_N=c_N,
            M_N=M_N,
            lightest_nu_mass=lightest_nu_mass,
            ordering=config.ordering,
            majorana_alpha=config.majorana_alpha,
            majorana_beta=config.majorana_beta,
            k=config.k,
        )
    except Exception as exc:
        row["reject_reason"] = f"exception:{exc}"
        return row

    all_y = np.concatenate([np.abs(result.Y_E_bar), np.abs(result.Y_N_bar)])
    max_y = float(np.max(all_y))
    min_y = float(np.min(all_y))

    row.update(
        {
            "Y_E_bar_1": result.Y_E_bar[0],
            "Y_E_bar_2": result.Y_E_bar[1],
            "Y_E_bar_3": result.Y_E_bar[2],
            "Y_N_bar_1": result.Y_N_bar[0],
            "Y_N_bar_2": result.Y_N_bar[1],
            "Y_N_bar_3": result.Y_N_bar[2],
            "f_L": result.f_L,
            "f_N": result.f_N,
            "f_N_UV": result.f_N_UV,
            "max_Y_bar": max_y,
        }
    )

    reasons: List[str] = []

    # 1) Perturbativity
    is_perturbative = result.is_perturbative(config.max_Y_bar)
    row["perturbative"] = is_perturbative
    if not is_perturbative:
        reasons.append("perturbativity")

    # 2) Naturalness
    lo, hi = config.naturalness_range
    is_natural = (min_y >= lo) and (max_y <= hi)
    row["natural"] = is_natural
    if not is_natural:
        reasons.append("naturalness")

    # 3) mu -> e gamma
    lfv = check_mu_to_e_gamma(
        result,
        C=lfv_C,
        reference_scale=config.lfv_reference_scale,
        M_KK_override=M_KK,
    )
    row["lfv_passes"] = bool(lfv["passes"])
    row["lfv_lhs"] = float(lfv["lhs"])
    row["lfv_rhs"] = float(lfv["rhs"])
    row["lfv_ratio"] = float(lfv["ratio"])
    if not row["lfv_passes"]:
        reasons.append("mu_to_e_gamma")

    # 4) User-supplied filters
    for filt in extra_filters:
        ok, label = filt(result)
        if not ok:
            reasons.append(label)

    # 5) Optional anarchic-prior score.
    if config.anarchy is not None:
        if rng_seed_sample is None:
            rng = np.random.default_rng()
        else:
            rng = np.random.default_rng(rng_seed_sample)
        anarchy_state = sample_anarchy_state(rng=rng, config=config.anarchy)
        row["anarchy_score"] = anarchy_state["score"]
        row["anarchy_band_penalty"] = anarchy_state["band_penalty"]
        row["anarchy_condition_penalty"] = anarchy_state["condition_penalty"]
        row["anarchy_w_band"] = anarchy_state["w_band"]
        row["anarchy_w_cond"] = anarchy_state["w_cond"]
        row["anarchy_w_fit"] = anarchy_state["w_fit"]
        row["anarchy_yN_overall"] = anarchy_state["yN_overall"]
        if config.anarchy_min_score is not None and row["anarchy_score"] < config.anarchy_min_score:
            reasons.append("anarchy_score")

    row["reject_reason"] = ";".join(reasons)
    row["passes_all"] = len(reasons) == 0
    return row


# ── Main scan loop ───────────────────────────────────────────────────────


def run_scan(
    config: ScanConfig,
    output_csv: Optional[str] = None,
    extra_filters: Optional[List[Callable[[YukawaResult], Tuple[bool, str]]]] = None,
    progress_every: int = 100,
) -> List[Dict[str, Any]]:
    """Run a parameter scan over the configured grid.

    Parameters
    ----------
    config : ScanConfig
        Scan configuration.
    output_csv : str or None
        Path to write CSV results. If None, results are only returned.
    extra_filters : list of callables or None
        Each callable takes ``YukawaResult`` and returns ``(passes, label)``.
    progress_every : int
        Print progress every N points. Set to 0 to suppress.

    Returns
    -------
    list of dict
        One row per scan point.
    """
    if extra_filters is None:
        extra_filters = []

    lfv_C = coefficient_from_br_limit(config.br_limit, prefactor=config.prefac_br)
    git_commit, dirty_tree = _resolve_git_metadata(config.record_git_metadata)

    total = config.total_points
    results: List[Dict[str, Any]] = []
    n_pass = 0
    t_start = time.time()

    scan_iter = itertools.product(
        config.Lambda_IR_values,
        config.c_L_values,
        config.c_N_values,
        config.c_E_points,
        config.MN_over_k_points,
        config.lightest_nu_mass_values,
    )

    for sample_index, (Lambda_IR, c_L, c_N, c_E, MN_over_k, lightest_nu_mass) in enumerate(scan_iter):
        seed = _sample_seed(config.rng_seed_global, sample_index)

        row = _evaluate_point(
            sample_index=sample_index,
            Lambda_IR=float(Lambda_IR),
            c_L=float(c_L),
            c_N=float(c_N),
            c_E=c_E,
            MN_over_k=float(MN_over_k),
            lightest_nu_mass=float(lightest_nu_mass),
            config=config,
            lfv_C=lfv_C,
            git_commit=git_commit,
            dirty_tree=dirty_tree,
            rng_seed_sample=seed,
            extra_filters=extra_filters,
        )
        results.append(row)

        if row["passes_all"]:
            n_pass += 1

        if progress_every > 0 and (sample_index + 1) % progress_every == 0:
            elapsed = time.time() - t_start
            rate = (sample_index + 1) / elapsed if elapsed > 0 else 0
            print(
                f"  [{sample_index + 1}/{total}]  accepted: {n_pass}  ({rate:.0f} pts/s)"
            )

    elapsed = time.time() - t_start
    print(
        f"Scan complete: {total} points in {elapsed:.1f}s, "
        f"{n_pass} accepted ({100 * n_pass / max(total, 1):.1f}%)"
    )

    if output_csv is not None:
        with open(output_csv, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=CSV_COLUMNS)
            writer.writeheader()
            writer.writerows(results)
        print(f"Results written to {output_csv}")

    return results
