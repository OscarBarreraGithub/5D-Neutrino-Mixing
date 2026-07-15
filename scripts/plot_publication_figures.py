#!/usr/bin/env python3
"""Generate publication-quality figures from quark-sector scan results.

Produces two figures:

1. **System-by-system exclusion boundaries** in the (r, M_KK) plane with
   individual contours for each meson system and the combined exclusion.

2. **SUPERSEDED/CORRECTED 2007-vs-modern notice** replacing the former
   post-hoc rescaled lower-bound comparison.

Usage
-----
    python scripts/plot_publication_figures.py <results.jsonl> \
        --output-dir results/figures/quark
"""

from __future__ import annotations

import argparse
import json
import sys
import textwrap
from collections import defaultdict
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

# NOTE: no top-level `matplotlib.use(...)` call. Setting a backend at import
# time hijacks the backend of any notebook that imports this module for its
# helper functions (e.g. reproduce_fig1_fig2.ipynb uses %matplotlib inline).
# The CLI entry point in main() sets Agg explicitly before doing any plotting.
import matplotlib  # noqa: F401
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.interpolate import griddata  # noqa: E402

from quarkConstraints.scales import GAUGE_KK_ROOT_NN  # noqa: E402

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

SYSTEM_IDS = ("epsilon_K", "K", "B_d", "B_s", "D0")

SYSTEM_LABELS: dict[str, str] = {
    "epsilon_K": r"$\varepsilon_K$",
    "K": r"$\Delta M_K$",
    "B_d": r"$\Delta M_{B_d}$",
    "B_s": r"$\Delta M_{B_s}$",
    "D0": r"$\Delta M_{D^0}$",
}

SYSTEM_CONTOUR_STYLES: dict[str, dict[str, Any]] = {
    "epsilon_K": {"color": "#7B2D8E", "linestyle": "solid", "linewidth": 2},
    "K": {"color": "#2166AC", "linestyle": "dashed", "linewidth": 2},
    "B_d": {"color": "#D95F02", "linestyle": "dashdot", "linewidth": 2},
    "B_s": {"color": "#C0392B", "linestyle": "dotted", "linewidth": 2},
    "D0": {"color": "#1B7837", "linestyle": (0, (3, 1, 1, 1, 1, 1)), "linewidth": 2},
}

SUPERSEDED_2007_COMPARISON_TITLE = (
    "SUPERSEDED/CORRECTED: 2007-vs-modern all-system rescale"
)

SUPERSEDED_2007_COMPARISON_NOTICE = (
    "The former 2007-vs-modern lower-bound comparison is retracted. Modern scan "
    "rows store ratios against hadronic |M12| budgets in GeV for B_d, B_s, and "
    "D0. The old rescale used legacy dimensionless operator-weight bounds for "
    "those systems, so it multiplied incompatible quantities. The prior 9.4x "
    "D0 tightening claim and the shaded tightening band were units artifacts."
)

DEFAULT_PUBLICATION_XI_KK = GAUGE_KK_ROOT_NN


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "results_file",
        type=Path,
        help="Path to merged scan results JSONL file.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help=(
            "Directory for output figures.  "
            "Defaults to results/figures/quark relative to the repo root."
        ),
    )
    parser.add_argument(
        "--xi-kk",
        type=float,
        default=DEFAULT_PUBLICATION_XI_KK,
        help=(
            "Physical KK-gluon mass convention used on the plot axis, "
            "m_g^(1) = xi_kk * Lambda_IR. "
            f"Default: {DEFAULT_PUBLICATION_XI_KK:.6f}."
        ),
    )
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def _load_results(path: Path) -> list[dict[str, Any]]:
    """Load JSONL results, skipping lines that fail to parse."""
    records: list[dict[str, Any]] = []
    with open(path) as fh:
        for line_no, line in enumerate(fh, 1):
            line = line.strip()
            if not line:
                continue
            try:
                record = json.loads(line)
            except json.JSONDecodeError as exc:
                print(f"WARNING: skipping line {line_no}: {exc}", file=sys.stderr)
                continue
            records.append(record)
    return records


def _extract_arrays(records: list[dict[str, Any]]) -> dict[str, Any]:
    """Extract numpy arrays from loaded records.

    Points that did not converge are dropped.  Supports both
    ``fit_converged`` and ``fit_success`` field names.
    """
    r_vals: list[float] = []
    mkk_vals: list[float] = []
    scale_vals: list[float] = []
    xi_vals: list[float] = []
    accepted_vals: list[bool] = []
    ratios_list: list[dict[str, float]] = []
    failing_list: list[list[str]] = []

    for rec in records:
        converged = rec.get("fit_converged", rec.get("fit_success", True))
        if not converged:
            continue
        # Skip points with poor fit quality (score > 0.1 means quark
        # masses/CKM are not reproduced, so Wilson coefficients are wrong)
        fit_score = rec.get("fit_score", 0.0)
        if fit_score > 0.1:
            continue

        r_vals.append(float(rec["r"]))
        mkk = float(rec.get("M_KK", rec["Lambda_IR"]))
        mkk_vals.append(mkk)
        scale_vals.append(float(rec["overall_scale"]))
        xi_vals.append(float(rec.get("xi_KK", 1.0)))
        accepted_vals.append(bool(rec["accepted"]))
        ratios_list.append(
            {str(k): float(v) for k, v in rec.get("ratio_to_bound_by_system", {}).items()}
        )
        failing_list.append(
            [str(s) for s in rec.get("failing_system_ids", [])]
        )

    return {
        "r": np.asarray(r_vals),
        "M_KK": np.asarray(mkk_vals),
        "overall_scale": np.asarray(scale_vals),
        "xi_KK": np.asarray(xi_vals),
        "accepted": np.asarray(accepted_vals, dtype=bool),
        "ratios": ratios_list,
        "failing": failing_list,
    }


def _apply_publication_convention(
    data: dict[str, Any],
    *,
    xi_kk: float,
) -> dict[str, Any]:
    """Relabel the saved scan into a paper-like KK-gluon convention.

    The stored `repo_v1` scan rows use the bookkeeping convention
    `M_KK = xi_scan * Lambda_IR`.  Publication plots should instead be
    expressed in terms of a physical first KK-gluon mass
    `m_g^(1) = xi_kk * Lambda_IR`.

    Because the tree-level `Delta F = 2` ratios scale as `1 / M^2`, the saved
    ratios can be converted post hoc without re-running the fit.
    """
    if xi_kk <= 0.0:
        raise ValueError("xi_kk must be positive")

    scan_mkk = np.asarray(data["M_KK"], dtype=float)
    scan_xi = np.asarray(data["xi_KK"], dtype=float)

    mkk_out: list[float] = []
    ratios_out: list[dict[str, float]] = []
    accepted_out: list[bool] = []

    for idx, scan_mass in enumerate(scan_mkk):
        xi_scan = float(scan_xi[idx])
        axis_scale = xi_kk / xi_scan
        mkk_out.append(float(scan_mass) * axis_scale)

        # Only convert the mass convention. The saved scan already carries the
        # coupling choice used at evaluation time.
        ratio_scale = 1.0 / (axis_scale**2)

        ratios = {
            system_id: float(data["ratios"][idx].get(system_id, 0.0)) * ratio_scale
            for system_id in SYSTEM_IDS
        }
        ratios_out.append(ratios)
        accepted_out.append(max(ratios.values()) <= 1.0)

    out = dict(data)
    out["M_KK"] = np.asarray(mkk_out, dtype=float)
    out["ratios"] = ratios_out
    out["accepted"] = np.asarray(accepted_out, dtype=bool)
    out["publication_convention"] = {
        "xi_kk": float(xi_kk),
        "scan_coupling_mode": "saved_scan",
    }
    return out


# ---------------------------------------------------------------------------
# Publication styling
# ---------------------------------------------------------------------------

def _configure_style() -> None:
    """Set up clean publication-quality matplotlib styling."""
    plt.rcParams.update({
        "figure.dpi": 150,
        "font.size": 12,
        "font.family": "serif",
        "mathtext.fontset": "cm",
        "text.usetex": False,
        "axes.titlesize": 14,
        "axes.labelsize": 14,
        "legend.fontsize": 11,
        "legend.framealpha": 0.92,
        "legend.edgecolor": "0.7",
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
        "axes.linewidth": 1.0,
        "xtick.major.width": 0.9,
        "ytick.major.width": 0.9,
        "xtick.minor.width": 0.5,
        "ytick.minor.width": 0.5,
        "xtick.major.size": 5,
        "ytick.major.size": 5,
        "xtick.minor.size": 3,
        "ytick.minor.size": 3,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "lines.linewidth": 1.8,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.08,
    })


# ---------------------------------------------------------------------------
# Aggregation helpers
# ---------------------------------------------------------------------------

def _aggregate_best_scale(
    data: dict[str, Any],
) -> tuple[np.ndarray, np.ndarray, dict[str, np.ndarray], np.ndarray]:
    """At each unique (r, M_KK) grid point, pick one common best-fit point.

    The earlier plotting code took an independent envelope over
    `overall_scale` for each system.  That was overly optimistic because the
    five exclusion contours no longer came from one physical MFV point.
    Publication figures instead use the single overall_scale that minimizes the
    combined max(all 5 ratios), and then read off every system from that same
    point.

    Returns
    -------
    r_agg, mkk_agg : 1-D arrays of unique (r, M_KK) values.
    ratios_agg : dict mapping system_id -> 1-D array of ratio values evaluated
        at the combined best point.
    max_ratio_agg : 1-D array of the minimized combined max-ratio values.
    """
    r = data["r"]
    mkk = data["M_KK"]
    ratios_list = data["ratios"]

    grid: dict[tuple[float, float], list[int]] = defaultdict(list)
    for i in range(len(r)):
        grid[(r[i], mkk[i])].append(i)

    r_out: list[float] = []
    mkk_out: list[float] = []
    sys_ratios: dict[str, list[float]] = {s: [] for s in SYSTEM_IDS}
    max_ratios: list[float] = []

    for (r_val, mkk_val), indices in grid.items():
        # Filter to indices that have non-empty ratio dicts
        valid_indices = [idx for idx in indices if ratios_list[idx]]
        if not valid_indices:
            continue

        r_out.append(r_val)
        mkk_out.append(mkk_val)

        best_idx = min(
            valid_indices,
            key=lambda idx: max(ratios_list[idx].get(s, 0.0) for s in SYSTEM_IDS),
        )
        best_ratios = ratios_list[best_idx]
        max_ratios.append(max(best_ratios.get(s, 0.0) for s in SYSTEM_IDS))

        # Every system is read off from the same combined best point.
        for s in SYSTEM_IDS:
            sys_ratios[s].append(best_ratios.get(s, 0.0))

    r_agg = np.asarray(r_out)
    mkk_agg = np.asarray(mkk_out)
    ratios_agg = {s: np.asarray(vals) for s, vals in sys_ratios.items()}
    max_ratio_agg = np.asarray(max_ratios)

    return r_agg, mkk_agg, ratios_agg, max_ratio_agg


# ---------------------------------------------------------------------------
# Figure 1: System-by-system exclusion boundaries
# ---------------------------------------------------------------------------

def build_exclusion_boundaries_figure(data: dict[str, Any]) -> plt.Figure | None:
    """Build fig 1 (per-system + combined exclusion contours) and return the Figure.

    Returns ``None`` if there are too few grid points to interpolate.
    """
    r_agg, mkk_agg, ratios_agg, max_ratio_agg = _aggregate_best_scale(data)

    if len(r_agg) < 4:
        print("WARNING: too few grid points for contour interpolation.", file=sys.stderr)
        return None

    # Convert M_KK from GeV to TeV for plotting
    mkk_tev = mkk_agg / 1e3

    # Build a fine regular grid in log space
    log_r = np.log10(r_agg)
    log_mkk = np.log10(mkk_tev)

    n_grid = 300
    r_lin = np.linspace(log_r.min(), log_r.max(), n_grid)
    mkk_lin = np.linspace(log_mkk.min(), log_mkk.max(), n_grid)
    R_grid, M_grid = np.meshgrid(r_lin, mkk_lin)
    points = np.column_stack([log_r, log_mkk])

    # Interpolate each system's ratio onto the fine grid
    interp_ratios: dict[str, np.ndarray] = {}
    for sys_id in SYSTEM_IDS:
        vals = np.log10(np.maximum(ratios_agg[sys_id], 1e-10))
        interp_ratios[sys_id] = griddata(
            points, vals, (R_grid, M_grid), method="cubic", fill_value=np.nan
        )

    # Combined: max ratio across systems
    max_ratio_log = np.log10(np.maximum(max_ratio_agg, 1e-10))
    interp_max = griddata(
        points, max_ratio_log, (R_grid, M_grid), method="cubic", fill_value=np.nan
    )

    # --- Plot ---
    fig, ax = plt.subplots(figsize=(8, 6))

    # Shade allowed (above combined contour) and excluded (below)
    # Use contourf on the combined max-ratio field
    # ratio < 1 => log < 0 => allowed; ratio > 1 => log > 0 => excluded
    allowed_mask = np.where(np.isnan(interp_max), np.nan, interp_max)

    # Shading: excluded region in light red, allowed in light green
    ax.contourf(
        R_grid, M_grid, allowed_mask,
        levels=[-10, 0.0],
        colors=["#FFCCCC"],  # light red for excluded (ratio > 1 below contour)
        alpha=0.4,
        extend="min",
    )
    ax.contourf(
        R_grid, M_grid, allowed_mask,
        levels=[0.0, 10],
        colors=["#FFCCCC"],  # will be overwritten
        alpha=0.0,
    )
    # The logic: in the (log_r, log_mkk) plane, higher M_KK => smaller ratio
    # => allowed is ABOVE the contour.  contourf with levels [-10, 0] shades
    # where log(max_ratio) < 0, i.e. the allowed region.
    # Let's redo this more carefully.
    #
    # At fixed r, increasing M_KK reduces the ratio (KK states decouple).
    # So the ALLOWED region (ratio < 1, log < 0) is at high M_KK (upper part).
    # The EXCLUDED region (ratio > 1, log > 0) is at low M_KK (lower part).

    # Clear previous contourf and redo
    ax.clear()

    # Allowed: log(max_ratio) < 0
    ax.contourf(
        R_grid, M_grid, allowed_mask,
        levels=[-10, 0.0],
        colors=["#C8F7C5"],  # light green
        alpha=0.45,
    )
    # Excluded: log(max_ratio) > 0
    ax.contourf(
        R_grid, M_grid, allowed_mask,
        levels=[0.0, 10],
        colors=["#F5B7B1"],  # light pink/red
        alpha=0.45,
    )

    # Individual system contours at ratio = 1 (log10 = 0)
    for sys_id in SYSTEM_IDS:
        style = SYSTEM_CONTOUR_STYLES[sys_id]
        cs = ax.contour(
            R_grid, M_grid, interp_ratios[sys_id],
            levels=[0.0],
            colors=[style["color"]],
            linewidths=[style["linewidth"]],
            linestyles=[style["linestyle"]],
        )
        # Invisible artist for the legend
        ax.plot(
            [], [],
            color=style["color"],
            linestyle=style["linestyle"],
            linewidth=style["linewidth"],
            label=SYSTEM_LABELS[sys_id],
        )

    # The green/red shading already marks the combined allowed/excluded
    # boundary, so no additional black contour line is needed.

    # Region labels
    # Find a point safely in the allowed region (high M_KK, mid r)
    mid_r = 0.5 * (log_r.min() + log_r.max())
    high_mkk = log_mkk.max() - 0.15 * (log_mkk.max() - log_mkk.min())
    low_mkk = log_mkk.min() + 0.25 * (log_mkk.max() - log_mkk.min())

    ax.text(
        mid_r, high_mkk, "Allowed",
        fontsize=16, fontweight="bold", color="#1B7837",
        ha="center", va="center", alpha=0.8,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="none", alpha=0.6),
    )
    ax.text(
        mid_r, low_mkk, "Excluded",
        fontsize=16, fontweight="bold", color="#C0392B",
        ha="center", va="center", alpha=0.8,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="none", alpha=0.6),
    )

    # Convert log-scale tick labels back to real values
    ax.set_xlabel(r"$r$", fontsize=14)
    ax.set_ylabel(r"$m_{g^{(1)}}$ [TeV]", fontsize=14)

    # Custom tick formatter for log-scale axes
    _set_log_ticks(ax, "x", log_r.min(), log_r.max())
    _set_log_ticks(ax, "y", log_mkk.min(), log_mkk.max())

    ax.legend(
        loc="upper left",
        framealpha=0.92,
        edgecolor="0.7",
        fontsize=10,
        ncol=2,
    )

    fig.tight_layout()
    return fig


def _plot_exclusion_boundaries(data: dict[str, Any], output_dir: Path) -> list[Path]:
    """Build fig 1 and save PDF + PNG to ``output_dir``."""
    fig = build_exclusion_boundaries_figure(data)
    if fig is None:
        return []
    saved: list[Path] = []
    for ext in ("pdf", "png"):
        out = output_dir / f"fig1_exclusion_boundaries.{ext}"
        fig.savefig(out, dpi=300, format=ext)
        saved.append(out)
    plt.close(fig)
    return saved


def _set_log_ticks(ax: plt.Axes, axis: str, vmin: float, vmax: float) -> None:
    """Set tick labels that show real values on a linearly-spaced log10 axis."""
    from matplotlib.ticker import FuncFormatter

    def _fmt(val: float, _pos: Any) -> str:
        real = 10**val
        if real >= 1.0:
            if real == int(real):
                return f"{int(real)}"
            return f"{real:.1f}"
        if real >= 0.01:
            return f"{real:.2f}"
        return f"{real:.1e}"

    if axis == "x":
        ax.xaxis.set_major_formatter(FuncFormatter(_fmt))
    else:
        ax.yaxis.set_major_formatter(FuncFormatter(_fmt))


# ---------------------------------------------------------------------------
# Figure 2: M_KK lower bound -- 2007 vs modern
# ---------------------------------------------------------------------------

def _compute_min_mkk_by_r(
    data: dict[str, Any],
    *,
    use_2007: bool = False,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """For each unique r value, find the minimum M_KK at which a point passes.

    The former ``use_2007=True`` mode is disabled by C-7: it mixed hadronic
    |M12|/budget ratios with legacy dimensionless operator-weight bounds.

    Returns
    -------
    r_vals : sorted unique r values
    min_mkk : minimum passing M_KK (GeV) for each r (NaN if none passes)
    has_passing : boolean mask -- True where at least one point passes
    """
    r = data["r"]
    mkk = data["M_KK"]
    ratios_list = data["ratios"]
    if use_2007:
        raise ValueError(
            "2007 rescaling is superseded: no same-convention B_d/B_s/D0 "
            "hadronic 2007 budgets are available in the scan rows."
        )

    # Group by r
    by_r: dict[float, list[tuple[float, dict[str, float]]]] = defaultdict(list)
    for i in range(len(r)):
        by_r[r[i]].append((mkk[i], ratios_list[i]))

    sorted_r = np.sort(np.array(list(by_r.keys())))
    min_mkk = np.full_like(sorted_r, np.nan)
    has_passing = np.zeros_like(sorted_r, dtype=bool)

    for j, r_val in enumerate(sorted_r):
        pts = by_r[r_val]
        best_mkk = float("inf")
        found = False
        for mkk_val, ratios in pts:
            if not ratios:
                continue

            modern_ratios = {s: ratios.get(s, 0.0) for s in SYSTEM_IDS}

            # Aggregate over overall_scale: accept if max ratio <= 1
            max_r = max(modern_ratios.values())
            if max_r <= 1.0 and mkk_val < best_mkk:
                best_mkk = mkk_val
                found = True

        if found:
            min_mkk[j] = best_mkk
            has_passing[j] = True

    return sorted_r, min_mkk, has_passing


def build_mkk_bound_comparison_figure(data: dict[str, Any]) -> plt.Figure:
    """Build fig 2 as a superseded/corrected banner."""
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.axis("off")
    ax.text(
        0.5,
        0.72,
        SUPERSEDED_2007_COMPARISON_TITLE,
        ha="center",
        va="center",
        fontsize=14,
        fontweight="bold",
    )
    ax.text(
        0.5,
        0.45,
        textwrap.fill(SUPERSEDED_2007_COMPARISON_NOTICE, width=82),
        ha="center",
        va="center",
        fontsize=10.5,
        linespacing=1.35,
    )
    ax.text(
        0.5,
        0.18,
        "This figure intentionally contains no 2007-vs-modern lower-bound curve.",
        ha="center",
        va="center",
        fontsize=10,
        style="italic",
    )
    fig.tight_layout()
    return fig


def _plot_mkk_bound_comparison(data: dict[str, Any], output_dir: Path) -> list[Path]:
    """Build fig 2 and save PDF + PNG to ``output_dir``."""
    fig = build_mkk_bound_comparison_figure(data)
    saved: list[Path] = []
    for ext in ("pdf", "png"):
        out = output_dir / f"fig2_mkk_bound_2007_vs_modern.{ext}"
        fig.savefig(out, dpi=300, format=ext)
        saved.append(out)
    plt.close(fig)
    return saved


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    matplotlib.use("Agg")
    args = _parse_args()
    _configure_style()

    repo_root = Path(__file__).resolve().parents[1]
    output_dir = args.output_dir or (repo_root / "results" / "figures" / "quark")
    output_dir.mkdir(parents=True, exist_ok=True)

    results_path: Path = args.results_file
    if not results_path.is_file():
        print(f"ERROR: results file not found: {results_path}", file=sys.stderr)
        return 1

    print(f"Loading results from {results_path} ...")
    records = _load_results(results_path)
    if not records:
        print("ERROR: no valid records found in results file.", file=sys.stderr)
        return 1

    data = _extract_arrays(records)
    data = _apply_publication_convention(
        data,
        xi_kk=float(args.xi_kk),
    )
    n_total = len(data["r"])
    n_accepted = int(np.sum(data["accepted"]))
    print(f"Loaded {n_total} converged points ({n_accepted} accepted).")
    print(
        "Publication convention: "
        f"m_g^(1) = {args.xi_kk:.6f} * Lambda_IR; "
        "saved scan couplings left unchanged."
    )

    if n_total == 0:
        print("ERROR: no converged points to plot.", file=sys.stderr)
        return 1

    # --- Figure 1 ---
    print("Plotting Figure 1: system-by-system exclusion boundaries ...")
    saved_1 = _plot_exclusion_boundaries(data, output_dir)
    for p in saved_1:
        print(f"  Saved: {p}")

    # --- Figure 2 ---
    print("Writing Figure 2 superseded/corrected 2007-vs-modern notice ...")
    saved_2 = _plot_mkk_bound_comparison(data, output_dir)
    for p in saved_2:
        print(f"  Saved: {p}")

    all_saved = saved_1 + saved_2
    print(f"\nDone. {len(all_saved)} files written to {output_dir}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
