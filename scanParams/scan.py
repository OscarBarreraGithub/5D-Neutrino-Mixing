"""Grid-scan driver for the RS lepton-sector parameter space.

Loops over bulk mass parameters (c_L, c_N, c_E) at fixed geometry and
Majorana mass, evaluates Yukawa couplings via ``compute_all_yukawas()``,
applies perturbativity / naturalness / mu-to-e-gamma filters, and writes
accepted and rejected points to CSV for later inspection.
"""

from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional, Tuple
import csv
import itertools
import time

import numpy as np

from yukawa import compute_all_yukawas, YukawaResult
from flavorConstraints import check_mu_to_e_gamma


# ── CSV schema ───────────────────────────────────────────────────────────

CSV_COLUMNS = [
    'c_L', 'c_N', 'c_E1', 'c_E2', 'c_E3',
    'Lambda_IR', 'M_N', 'lightest_nu_mass', 'ordering',
    'Y_E_bar_1', 'Y_E_bar_2', 'Y_E_bar_3',
    'Y_N_bar_1', 'Y_N_bar_2', 'Y_N_bar_3',
    'f_L', 'f_N', 'f_N_UV', 'max_Y_bar',
    'perturbative', 'natural', 'lfv_passes', 'lfv_ratio',
    'passes_all', 'reject_reason',
]


# ── Configuration ────────────────────────────────────────────────────────

@dataclass
class ScanConfig:
    """Configuration for a parameter space scan.

    Parameters
    ----------
    Lambda_IR : float
        Fixed KK scale in GeV.
    M_N : float
        Fixed UV Majorana mass in GeV.
    lightest_nu_mass : float
        Lightest neutrino mass in eV.
    ordering : str
        ``'normal'`` or ``'inverted'``.
    c_L_values, c_N_values : array-like
        1-D grids of bulk mass parameters to scan.
    c_E_fixed : list of 3 floats or None
        Fixed charged-lepton bulk masses (used when *c_E_grid* is None).
    c_E_grid : list of 3 arrays or None
        If provided, independent grids for [c_E1, c_E2, c_E3].
    max_Y_bar : float
        Perturbativity bound on rescaled Yukawas.
    naturalness_range : (float, float)
        Accepted range for every |Y_bar| entry.
    lfv_C : float
        Coefficient for the mu-to-e-gamma bound.
    lfv_reference_scale : float
        Reference KK scale for the LFV bound (GeV).
    """

    # Fixed physics
    Lambda_IR: float = 3000.0
    M_N: float = 1.22e18
    lightest_nu_mass: float = 0.002
    ordering: str = 'normal'

    # Scan grids
    c_L_values: np.ndarray = field(
        default_factory=lambda: np.linspace(0.50, 0.70, 11))
    c_N_values: np.ndarray = field(
        default_factory=lambda: np.linspace(0.20, 0.50, 11))

    # Charged-lepton bulk masses
    c_E_fixed: Optional[List[float]] = field(
        default_factory=lambda: [0.75, 0.60, 0.50])
    c_E_grid: Optional[List[np.ndarray]] = None

    # Filter thresholds
    max_Y_bar: float = 4.0
    naturalness_range: Tuple[float, float] = (0.1, 4.0)
    lfv_C: float = 0.02
    lfv_reference_scale: float = 3000.0

    def __post_init__(self):
        if self.c_E_grid is not None:
            if len(self.c_E_grid) != 3:
                raise ValueError(
                    f"c_E_grid must contain exactly 3 arrays, got {len(self.c_E_grid)}")
        elif self.c_E_fixed is None:
            raise ValueError("Either c_E_fixed or c_E_grid must be provided")
        elif len(self.c_E_fixed) != 3:
            raise ValueError(
                f"c_E_fixed must contain exactly 3 values, got {len(self.c_E_fixed)}")

    @property
    def total_points(self) -> int:
        """Total number of grid points."""
        n = len(self.c_L_values) * len(self.c_N_values)
        if self.c_E_grid is not None:
            for arr in self.c_E_grid:
                n *= len(arr)
        return n


# ── Single-point evaluation ──────────────────────────────────────────────

def _evaluate_point(
    c_L: float,
    c_N: float,
    c_E: List[float],
    config: ScanConfig,
    extra_filters: List[Callable[[YukawaResult], Tuple[bool, str]]],
) -> Dict[str, Any]:
    """Evaluate one parameter point and return a row dict."""
    row: Dict[str, Any] = {
        'c_L': c_L, 'c_N': c_N,
        'c_E1': c_E[0], 'c_E2': c_E[1], 'c_E3': c_E[2],
        'Lambda_IR': config.Lambda_IR, 'M_N': config.M_N,
        'lightest_nu_mass': config.lightest_nu_mass,
        'ordering': config.ordering,
    }

    nan_defaults: Dict[str, Any] = {
        'Y_E_bar_1': np.nan, 'Y_E_bar_2': np.nan, 'Y_E_bar_3': np.nan,
        'Y_N_bar_1': np.nan, 'Y_N_bar_2': np.nan, 'Y_N_bar_3': np.nan,
        'f_L': np.nan, 'f_N': np.nan, 'f_N_UV': np.nan, 'max_Y_bar': np.nan,
        'perturbative': False, 'natural': False,
        'lfv_passes': False, 'lfv_ratio': np.nan,
        'passes_all': False, 'reject_reason': '',
    }

    try:
        result = compute_all_yukawas(
            Lambda_IR=config.Lambda_IR,
            c_L=c_L,
            c_E=c_E,
            c_N=c_N,
            M_N=config.M_N,
            lightest_nu_mass=config.lightest_nu_mass,
            ordering=config.ordering,
        )
    except Exception as exc:
        row.update(nan_defaults)
        row['reject_reason'] = f"exception:{exc}"
        return row

    # Extract outputs
    all_Y = np.concatenate([np.abs(result.Y_E_bar), np.abs(result.Y_N_bar)])
    max_y = float(np.max(all_Y))
    min_y = float(np.min(all_Y))

    row.update({
        'Y_E_bar_1': result.Y_E_bar[0],
        'Y_E_bar_2': result.Y_E_bar[1],
        'Y_E_bar_3': result.Y_E_bar[2],
        'Y_N_bar_1': result.Y_N_bar[0],
        'Y_N_bar_2': result.Y_N_bar[1],
        'Y_N_bar_3': result.Y_N_bar[2],
        'f_L': result.f_L,
        'f_N': result.f_N,
        'f_N_UV': result.f_N_UV,
        'max_Y_bar': max_y,
    })

    # ── Filters ──
    reasons: List[str] = []

    # 1. Perturbativity
    is_pert = result.is_perturbative(config.max_Y_bar)
    row['perturbative'] = is_pert
    if not is_pert:
        reasons.append('perturbativity')

    # 2. Naturalness
    lo, hi = config.naturalness_range
    is_natural = (min_y >= lo) and (max_y <= hi)
    row['natural'] = is_natural
    if not is_natural:
        reasons.append('naturalness')

    # 3. μ→eγ
    lfv = check_mu_to_e_gamma(
        result, C=config.lfv_C, reference_scale=config.lfv_reference_scale)
    row['lfv_passes'] = lfv['passes']
    row['lfv_ratio'] = lfv['ratio']
    if not lfv['passes']:
        reasons.append('mu_to_e_gamma')

    # 4. User-supplied filters
    for filt in extra_filters:
        ok, label = filt(result)
        if not ok:
            reasons.append(label)

    row['reject_reason'] = ';'.join(reasons)
    row['passes_all'] = len(reasons) == 0
    return row


# ── Main scan loop ───────────────────────────────────────────────────────

def run_scan(
    config: ScanConfig,
    output_csv: Optional[str] = None,
    extra_filters: Optional[List[Callable[[YukawaResult], Tuple[bool, str]]]] = None,
    progress_every: int = 100,
) -> List[Dict[str, Any]]:
    """Run a grid scan over the parameter space.

    Parameters
    ----------
    config : ScanConfig
        Scan configuration.
    output_csv : str or None
        Path to write CSV results.  If None, results are only returned.
    extra_filters : list of callables or None
        Each takes a *YukawaResult* and returns ``(passes, label)``.
    progress_every : int
        Print progress every *N* points.  Set to 0 to suppress.

    Returns
    -------
    list of dict
        One dict per grid point (columns match *CSV_COLUMNS*).
    """
    if extra_filters is None:
        extra_filters = []

    # Build c_E iterator
    if config.c_E_grid is not None:
        c_E_iter = list(itertools.product(
            config.c_E_grid[0], config.c_E_grid[1], config.c_E_grid[2]))
    else:
        c_E_iter = [tuple(config.c_E_fixed)]

    total = config.total_points
    results: List[Dict[str, Any]] = []
    n_pass = 0
    t_start = time.time()

    for idx, (c_L, c_N, c_E_tuple) in enumerate(
        itertools.product(config.c_L_values, config.c_N_values, c_E_iter)
    ):
        row = _evaluate_point(c_L, c_N, list(c_E_tuple), config, extra_filters)
        results.append(row)
        if row['passes_all']:
            n_pass += 1

        if progress_every > 0 and (idx + 1) % progress_every == 0:
            elapsed = time.time() - t_start
            rate = (idx + 1) / elapsed if elapsed > 0 else 0
            print(f"  [{idx + 1}/{total}]  accepted: {n_pass}  "
                  f"({rate:.0f} pts/s)")

    elapsed = time.time() - t_start
    print(f"Scan complete: {total} points in {elapsed:.1f}s, "
          f"{n_pass} accepted ({100 * n_pass / max(total, 1):.1f}%)")

    if output_csv is not None:
        with open(output_csv, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
            writer.writeheader()
            writer.writerows(results)
        print(f"Results written to {output_csv}")

    return results
