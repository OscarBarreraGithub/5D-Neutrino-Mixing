#!/usr/bin/env python3
"""Generate publication-quality figures from quark-sector scan results.

Produces two figures:

1. **System-by-system exclusion boundaries** in the (r, M_KK) plane with
   individual contours for each meson system and the combined exclusion.

2. **M_KK lower bound: 2007 vs modern** showing how constraints have
   tightened since the original Fitzpatrick-Perez-Randall paper.

Usage
-----
    python scripts/plot_publication_figures.py <results.jsonl> \
        --output-dir results/figures/publication
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import matplotlib
matplotlib.use("Agg")
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

# 2007 vs modern bound rescaling: ratio = bound_modern / bound_2007
# If < 1, modern bound is tighter.
BOUND_RATIOS: dict[str, float] = {
    "epsilon_K": 0.70,
    "K": 0.70,
    "B_d": 0.67,
    "B_s": 0.37,
    "D0": 0.106,
}

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
            "Defaults to results/figures/publication relative to the repo root."
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

def _plot_exclusion_boundaries(data: dict[str, Any], output_dir: Path) -> list[Path]:
    """Draw per-system and combined exclusion contours in (r, M_KK)."""
    r_agg, mkk_agg, ratios_agg, max_ratio_agg = _aggregate_best_scale(data)

    if len(r_agg) < 4:
        print("WARNING: too few grid points for contour interpolation.", file=sys.stderr)
        return []

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

    If *use_2007* is True, rescale ratios to 2007-era bounds before checking.

    Returns
    -------
    r_vals : sorted unique r values
    min_mkk : minimum passing M_KK (GeV) for each r (NaN if none passes)
    has_passing : boolean mask -- True where at least one point passes
    """
    r = data["r"]
    mkk = data["M_KK"]
    ratios_list = data["ratios"]

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

            if use_2007:
                rescaled = {s: ratios.get(s, 0.0) * BOUND_RATIOS.get(s, 1.0) for s in SYSTEM_IDS}
            else:
                rescaled = {s: ratios.get(s, 0.0) for s in SYSTEM_IDS}

            # Aggregate over overall_scale: accept if max ratio <= 1
            max_r = max(rescaled.values())
            if max_r <= 1.0 and mkk_val < best_mkk:
                best_mkk = mkk_val
                found = True

        if found:
            min_mkk[j] = best_mkk
            has_passing[j] = True

    return sorted_r, min_mkk, has_passing


def _plot_mkk_bound_comparison(data: dict[str, Any], output_dir: Path) -> list[Path]:
    """Figure 2: M_KK lower bound comparing 2007 and modern constraints."""
    r_modern, mkk_modern, has_modern = _compute_min_mkk_by_r(data, use_2007=False)
    r_2007, mkk_2007, has_2007 = _compute_min_mkk_by_r(data, use_2007=True)

    # Both arrays share the same r grid by construction
    assert np.array_equal(r_modern, r_2007)
    r_vals = r_modern

    # Convert GeV -> TeV
    mkk_modern_tev = mkk_modern / 1e3
    mkk_2007_tev = mkk_2007 / 1e3

    fig, ax = plt.subplots(figsize=(8, 5))

    # Plot modern bound
    m_mask = has_modern
    ax.plot(
        r_vals[m_mask], mkk_modern_tev[m_mask],
        marker="s", ms=5, lw=2.2,
        color="#2166AC", markerfacecolor="#2166AC",
        markeredgecolor="#1A4E7A", markeredgewidth=0.8,
        label=r"Modern (2024+) bound",
        zorder=4,
    )

    # Plot 2007 bound
    a_mask = has_2007
    ax.plot(
        r_vals[a_mask], mkk_2007_tev[a_mask],
        marker="o", ms=5, lw=2.2,
        color="#C0392B", markerfacecolor="#C0392B",
        markeredgecolor="#922B21", markeredgewidth=0.8,
        linestyle="dashed",
        label=r"2007 (FPR) bound",
        zorder=4,
    )

    # Shade the region between the two curves where both have data
    both_mask = has_modern & has_2007
    if np.any(both_mask):
        r_both = r_vals[both_mask]
        mkk_mod_both = mkk_modern_tev[both_mask]
        mkk_07_both = mkk_2007_tev[both_mask]
        ax.fill_between(
            r_both,
            mkk_07_both,
            mkk_mod_both,
            alpha=0.25,
            color="#9B59B6",
            label="Constraint tightening since 2007",
            zorder=2,
        )

    # Mark fully excluded r values (no passing point at any M_KK) with arrows
    excluded_modern = ~has_modern
    if np.any(excluded_modern):
        y_top = 20.0  # TeV -- mark as > 20 TeV
        ax.scatter(
            r_vals[excluded_modern],
            np.full(np.sum(excluded_modern), y_top),
            marker="^", s=60, c="#2166AC", edgecolors="#1A4E7A",
            linewidths=0.8, zorder=5,
        )

    excluded_2007 = ~has_2007
    if np.any(excluded_2007):
        y_top = 20.0
        ax.scatter(
            r_vals[excluded_2007],
            np.full(np.sum(excluded_2007), y_top),
            marker="^", s=60, c="#C0392B", edgecolors="#922B21",
            linewidths=0.8, zorder=5,
        )

    # (LHC direct search line removed for cleaner presentation)

    ax.set_xscale("log")
    ax.set_xlabel(r"$r$", fontsize=14)
    ax.set_ylabel(r"Minimum $m_{g^{(1)}}$ [TeV]", fontsize=14)

    ax.legend(
        loc="best",
        framealpha=0.92,
        edgecolor="0.7",
        fontsize=10,
    )
    ax.grid(True, which="both", alpha=0.2, linewidth=0.5)

    # Minor ticks
    ax.minorticks_on()

    fig.tight_layout()

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
    args = _parse_args()
    _configure_style()

    repo_root = Path(__file__).resolve().parents[1]
    output_dir = args.output_dir or (repo_root / "results" / "figures" / "publication")
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
    print("Plotting Figure 2: M_KK lower bound -- 2007 vs modern ...")
    saved_2 = _plot_mkk_bound_comparison(data, output_dir)
    for p in saved_2:
        print(f"  Saved: {p}")

    all_saved = saved_1 + saved_2
    print(f"\nDone. {len(all_saved)} files written to {output_dir}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
