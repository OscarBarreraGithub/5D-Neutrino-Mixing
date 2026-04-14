#!/usr/bin/env python3
"""Visualize quark-sector allowed Yukawa parameter space from scan results.

Reads merged scan results from a JSONL file (one JSON object per line) and
produces publication-quality figures showing the allowed region in the
(r, M_KK) plane, binding constraints, exclusion curves, and Yukawa structure.

Usage
-----
    python scripts/plot_allowed_region.py <results.jsonl> [--output-dir results/figures]
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.colors as mcolors  # noqa: E402
import numpy as np  # noqa: E402

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

SYSTEM_IDS = ("epsilon_K", "K", "B_d", "B_s", "D0")

SYSTEM_COLORS: dict[str, str] = {
    "epsilon_K": "#8e44ad",  # purple
    "K": "#2980b9",          # blue
    "B_d": "#e67e22",        # orange
    "B_s": "#c0392b",        # red
    "D0": "#27ae60",         # green
}

SYSTEM_LABELS: dict[str, str] = {
    "epsilon_K": r"$\epsilon_K$",
    "K": r"$K^0\text{-}\bar{K}^0$",
    "B_d": r"$B_d\text{-}\bar{B}_d$",
    "B_s": r"$B_s\text{-}\bar{B}_s$",
    "D0": r"$D^0\text{-}\bar{D}^0$",
}


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
            "Defaults to results/figures/ relative to the repo root."
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

    Handles both the user-facing field names from the prompt
    (``fit_converged``) and the actual codebase field names
    (``fit_success``).  Points that did not converge are dropped.
    """
    r_vals: list[float] = []
    mkk_vals: list[float] = []
    scale_vals: list[float] = []
    accepted_vals: list[bool] = []
    ratios_list: list[dict[str, float]] = []
    failing_list: list[list[str]] = []

    for rec in records:
        # Convergence: support both field names
        converged = rec.get("fit_converged", rec.get("fit_success", True))
        if not converged:
            continue

        r_vals.append(float(rec["r"]))
        # M_KK = Lambda_IR when xi_KK = 1; prefer explicit M_KK if present
        mkk = float(rec.get("M_KK", rec["Lambda_IR"]))
        mkk_vals.append(mkk)
        scale_vals.append(float(rec["overall_scale"]))
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
        "accepted": np.asarray(accepted_vals, dtype=bool),
        "ratios": ratios_list,
        "failing": failing_list,
    }


# ---------------------------------------------------------------------------
# Style
# ---------------------------------------------------------------------------

def _configure_style() -> None:
    plt.rcParams.update({
        "figure.dpi": 150,
        "font.size": 11,
        "font.family": "serif",
        "mathtext.fontset": "cm",
        "axes.titlesize": 14,
        "axes.labelsize": 13,
        "legend.fontsize": 10,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "xtick.minor.width": 0.5,
        "ytick.minor.width": 0.5,
        "lines.linewidth": 1.5,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.05,
    })


# ---------------------------------------------------------------------------
# Figure 1: Allowed region in (r, M_KK) plane
# ---------------------------------------------------------------------------

def _plot_allowed_region(data: dict[str, Any], output_dir: Path) -> Path:
    """Allowed / rejected scatter in the (r, M_KK) plane.

    One panel per unique ``overall_scale`` value.
    """
    r = data["r"]
    mkk = data["M_KK"] / 1e3  # GeV -> TeV
    accepted = data["accepted"]
    scales = data["overall_scale"]
    unique_scales = np.sort(np.unique(scales))

    n_panels = len(unique_scales)
    ncols = min(n_panels, 3)
    nrows = int(np.ceil(n_panels / ncols))

    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(5.5 * ncols, 4.5 * nrows),
        squeeze=False,
    )

    for idx, s_val in enumerate(unique_scales):
        ax = axes[idx // ncols, idx % ncols]
        mask = scales == s_val

        rej = mask & ~accepted
        acc = mask & accepted

        ax.scatter(
            r[rej], mkk[rej],
            c="#e74c3c", marker="x", s=20, alpha=0.5, linewidths=0.8,
            label="Rejected", rasterized=True,
        )
        ax.scatter(
            r[acc], mkk[acc],
            c="#27ae60", marker="o", s=22, alpha=0.7, edgecolors="none",
            label="Accepted", rasterized=True,
        )

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(r"$r$")
        ax.set_ylabel(r"$M_{\mathrm{KK}}$ [TeV]")
        ax.set_title(rf"overall\_scale $= {s_val:.2f}$")
        ax.legend(loc="upper left", framealpha=0.9, markerscale=1.3)
        ax.grid(True, which="both", alpha=0.25, linewidth=0.5)

    # Hide unused panels
    for idx in range(n_panels, nrows * ncols):
        axes[idx // ncols, idx % ncols].set_visible(False)

    fig.suptitle(
        r"Allowed Region in $(r,\, M_{\mathrm{KK}})$ Plane",
        fontsize=15, y=1.02,
    )
    fig.tight_layout()
    out = output_dir / "allowed_region_r_mkk.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return out


# ---------------------------------------------------------------------------
# Figure 2: Binding constraint map
# ---------------------------------------------------------------------------

def _binding_constraint(ratios: dict[str, float]) -> str | None:
    """Return the system ID with the highest ratio_to_bound."""
    if not ratios:
        return None
    return max(ratios, key=ratios.get)  # type: ignore[arg-type]


def _plot_binding_constraint(data: dict[str, Any], output_dir: Path) -> Path:
    """Color-code by the binding (highest-ratio) constraint system."""
    r = data["r"]
    mkk = data["M_KK"] / 1e3
    scales = data["overall_scale"]
    ratios_list = data["ratios"]
    unique_scales = np.sort(np.unique(scales))

    n_panels = len(unique_scales)
    ncols = min(n_panels, 3)
    nrows = int(np.ceil(n_panels / ncols))

    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(5.5 * ncols, 4.5 * nrows),
        squeeze=False,
    )

    for idx, s_val in enumerate(unique_scales):
        ax = axes[idx // ncols, idx % ncols]
        mask = np.where(scales == s_val)[0]

        plotted_systems: set[str] = set()
        for i in mask:
            binding = _binding_constraint(ratios_list[i])
            if binding is None:
                continue
            color = SYSTEM_COLORS.get(binding, "#7f8c8d")
            label = SYSTEM_LABELS.get(binding, binding) if binding not in plotted_systems else None
            ax.scatter(
                r[i], mkk[i],
                c=color, marker="o", s=18, alpha=0.6, edgecolors="none",
                label=label, rasterized=True,
            )
            plotted_systems.add(binding)

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(r"$r$")
        ax.set_ylabel(r"$M_{\mathrm{KK}}$ [TeV]")
        ax.set_title(rf"overall\_scale $= {s_val:.2f}$")
        ax.legend(loc="upper left", framealpha=0.9, fontsize=8, markerscale=1.5)
        ax.grid(True, which="both", alpha=0.25, linewidth=0.5)

    for idx in range(n_panels, nrows * ncols):
        axes[idx // ncols, idx % ncols].set_visible(False)

    fig.suptitle(
        r"Binding $\Delta F=2$ Constraint in $(r,\, M_{\mathrm{KK}})$ Plane",
        fontsize=15, y=1.02,
    )
    fig.tight_layout()
    out = output_dir / "binding_constraint_map.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return out


# ---------------------------------------------------------------------------
# Figure 3: Exclusion curve — minimum M_KK vs r
# ---------------------------------------------------------------------------

def _plot_exclusion_curve(data: dict[str, Any], output_dir: Path) -> Path:
    """For each r, find the minimum M_KK at which any point is accepted."""
    r = data["r"]
    mkk = data["M_KK"] / 1e3
    accepted = data["accepted"]

    unique_r = np.sort(np.unique(r))
    min_mkk: list[float] = []
    excluded_r: list[float] = []
    included_r: list[float] = []

    for r_val in unique_r:
        mask = (r == r_val) & accepted
        if np.any(mask):
            min_mkk.append(float(np.min(mkk[mask])))
            included_r.append(r_val)
        else:
            excluded_r.append(r_val)

    fig, ax = plt.subplots(figsize=(7, 5))

    if included_r:
        ax.plot(
            included_r, min_mkk,
            marker="s", ms=7, lw=2.0,
            color="#2c3e50", markerfacecolor="#2980b9",
            markeredgecolor="#2c3e50", markeredgewidth=1.0,
            label=r"Min allowed $M_{\mathrm{KK}}$",
            zorder=3,
        )
        ax.fill_between(
            included_r, min_mkk, [ax.get_ylim()[0] if ax.get_ylim()[0] > 0 else 0.1] * len(min_mkk),
            alpha=0.12, color="#e74c3c", label="Excluded",
        )

    if excluded_r:
        y_top = max(min_mkk) * 1.5 if min_mkk else 10.0
        for xr in excluded_r:
            ax.axvline(xr, color="#e74c3c", ls=":", lw=1.0, alpha=0.6)
        ax.scatter(
            excluded_r, [y_top] * len(excluded_r),
            marker="v", c="#e74c3c", s=60, zorder=4,
            label="Fully excluded (no passing point)",
        )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$r$")
    ax.set_ylabel(r"Minimum $M_{\mathrm{KK}}$ [TeV]")
    ax.set_title(r"Exclusion Curve: Minimum $M_{\mathrm{KK}}$ vs $r$")
    ax.legend(loc="upper right", framealpha=0.95)
    ax.grid(True, which="both", alpha=0.25, linewidth=0.5)

    fig.tight_layout()
    out = output_dir / "exclusion_curve_min_mkk_vs_r.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return out


# ---------------------------------------------------------------------------
# Figure 4: Yukawa structure vs M_KK
# ---------------------------------------------------------------------------

def _plot_yukawa_vs_mkk(data: dict[str, Any], output_dir: Path) -> Path:
    """For accepted points, show overall_scale vs M_KK, colored by r."""
    accepted = data["accepted"]
    if not np.any(accepted):
        # Nothing to plot — create a placeholder
        fig, ax = plt.subplots(figsize=(7, 5))
        ax.text(
            0.5, 0.5, "No accepted points",
            transform=ax.transAxes, ha="center", va="center", fontsize=14,
        )
        ax.set_title("Yukawa Scale vs $M_{\\mathrm{KK}}$ (no data)")
        out = output_dir / "yukawa_vs_mkk.png"
        fig.savefig(out, dpi=200)
        plt.close(fig)
        return out

    r_acc = data["r"][accepted]
    mkk_acc = data["M_KK"][accepted] / 1e3
    scale_acc = data["overall_scale"][accepted]

    fig, ax = plt.subplots(figsize=(7, 5))

    log_r = np.log10(r_acc)
    norm = mcolors.Normalize(vmin=log_r.min(), vmax=log_r.max())
    cmap = plt.cm.viridis  # type: ignore[attr-defined]

    sc = ax.scatter(
        scale_acc, mkk_acc,
        c=log_r, cmap=cmap, norm=norm,
        s=25, alpha=0.7, edgecolors="none",
        rasterized=True,
    )
    cbar = fig.colorbar(sc, ax=ax, pad=0.02)
    cbar.set_label(r"$\log_{10}(r)$")

    ax.set_xlabel("Overall Yukawa scale")
    ax.set_ylabel(r"$M_{\mathrm{KK}}$ [TeV]")
    ax.set_yscale("log")
    ax.set_title(r"Accepted Points: Yukawa Scale vs $M_{\mathrm{KK}}$")
    ax.grid(True, which="both", alpha=0.25, linewidth=0.5)

    fig.tight_layout()
    out = output_dir / "yukawa_vs_mkk.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return out


# ---------------------------------------------------------------------------
# Figure 5: Heat map of worst constraint ratio
# ---------------------------------------------------------------------------

def _max_ratio(ratios: dict[str, float]) -> float:
    """Return the maximum ratio_to_bound across all systems."""
    if not ratios:
        return 0.0
    return max(ratios.values())


def _plot_heatmap(data: dict[str, Any], output_dir: Path) -> Path:
    """2D histogram of log10(max ratio_to_bound) in (r, M_KK) plane."""
    r = data["r"]
    mkk = data["M_KK"] / 1e3
    scales = data["overall_scale"]
    ratios_list = data["ratios"]
    unique_scales = np.sort(np.unique(scales))

    max_ratios = np.array([_max_ratio(rt) for rt in ratios_list])

    n_panels = len(unique_scales)
    ncols = min(n_panels, 3)
    nrows = int(np.ceil(n_panels / ncols))

    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(5.5 * ncols, 4.5 * nrows),
        squeeze=False,
    )

    for idx, s_val in enumerate(unique_scales):
        ax = axes[idx // ncols, idx % ncols]
        mask = scales == s_val

        r_sub = r[mask]
        mkk_sub = mkk[mask]
        mr_sub = max_ratios[mask]

        if len(r_sub) == 0:
            ax.set_visible(False)
            continue

        # Use log10 of the ratio for color; clamp small values
        log_ratio = np.log10(np.maximum(mr_sub, 1e-6))

        # Build regular grid via pivot-style binning
        unique_r_sub = np.sort(np.unique(r_sub))
        unique_mkk_sub = np.sort(np.unique(mkk_sub))

        if len(unique_r_sub) > 1 and len(unique_mkk_sub) > 1:
            # Grid-based approach: build 2D array keyed on (r, mkk)
            grid = np.full((len(unique_mkk_sub), len(unique_r_sub)), np.nan)
            r_idx_map = {v: i for i, v in enumerate(unique_r_sub)}
            mkk_idx_map = {v: i for i, v in enumerate(unique_mkk_sub)}
            for i_pt in range(len(r_sub)):
                ri = r_idx_map.get(r_sub[i_pt])
                mi = mkk_idx_map.get(mkk_sub[i_pt])
                if ri is not None and mi is not None:
                    # Keep worst (maximum) ratio if multiple points at same location
                    existing = grid[mi, ri]
                    if np.isnan(existing) or log_ratio[i_pt] > existing:
                        grid[mi, ri] = log_ratio[i_pt]

            # For pcolormesh we need bin edges on a log scale
            log_r_edges = _bin_edges_log(unique_r_sub)
            log_mkk_edges = _bin_edges_log(unique_mkk_sub)

            vmin = np.nanmin(grid) if np.any(~np.isnan(grid)) else -3.0
            vmax = np.nanmax(grid) if np.any(~np.isnan(grid)) else 1.0

            pcm = ax.pcolormesh(
                log_r_edges, log_mkk_edges, grid,
                cmap="RdYlGn_r", vmin=vmin, vmax=vmax,
                shading="flat", rasterized=True,
            )
            cbar = fig.colorbar(pcm, ax=ax, pad=0.02)
            cbar.set_label(r"$\log_{10}(\max\, \mathrm{ratio})$")

            # Contour at ratio = 1 (log10 = 0) if data spans this threshold
            if vmin < 0.0 < vmax:
                # Contour needs cell-center coordinates
                r_centers = unique_r_sub
                mkk_centers = unique_mkk_sub
                rr, mm = np.meshgrid(np.log10(r_centers), np.log10(mkk_centers))
                with np.errstate(invalid="ignore"):
                    ax.contour(
                        rr, mm, grid,
                        levels=[0.0], colors="white", linewidths=2.0,
                    )

            ax.set_xlabel(r"$\log_{10}(r)$")
            ax.set_ylabel(r"$\log_{10}(M_{\mathrm{KK}}$ [TeV]$)$")
        else:
            # Fallback: simple scatter
            sc = ax.scatter(
                np.log10(r_sub), np.log10(mkk_sub),
                c=log_ratio, cmap="RdYlGn_r", s=30, edgecolors="none",
                rasterized=True,
            )
            cbar = fig.colorbar(sc, ax=ax, pad=0.02)
            cbar.set_label(r"$\log_{10}(\max\, \mathrm{ratio})$")
            ax.set_xlabel(r"$\log_{10}(r)$")
            ax.set_ylabel(r"$\log_{10}(M_{\mathrm{KK}}$ [TeV]$)$")

        ax.set_title(rf"overall\_scale $= {s_val:.2f}$")

    for idx in range(n_panels, nrows * ncols):
        axes[idx // ncols, idx % ncols].set_visible(False)

    fig.suptitle(
        r"Worst $\Delta F=2$ Ratio Heat Map",
        fontsize=15, y=1.02,
    )
    fig.tight_layout()
    out = output_dir / "worst_ratio_heatmap.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return out


def _bin_edges_log(centers: np.ndarray) -> np.ndarray:
    """Compute pcolormesh bin edges in log-space from sorted center values.

    Returns an array of length ``len(centers) + 1`` in log10 coordinates.
    """
    log_c = np.log10(centers)
    edges = np.empty(len(log_c) + 1)
    if len(log_c) == 1:
        edges[0] = log_c[0] - 0.5
        edges[1] = log_c[0] + 0.5
    else:
        half_gaps = np.diff(log_c) / 2.0
        edges[0] = log_c[0] - half_gaps[0]
        edges[-1] = log_c[-1] + half_gaps[-1]
        edges[1:-1] = log_c[:-1] + half_gaps
    return edges


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    args = _parse_args()
    _configure_style()

    repo_root = Path(__file__).resolve().parents[1]
    output_dir = args.output_dir or (repo_root / "results" / "figures")
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
    n_total = len(data["r"])
    n_accepted = int(np.sum(data["accepted"]))
    print(f"Loaded {n_total} converged points ({n_accepted} accepted).")

    if n_total == 0:
        print("ERROR: no converged points to plot.", file=sys.stderr)
        return 1

    # Generate all figures
    saved: list[Path] = []

    print("Plotting Figure 1: allowed region ...")
    saved.append(_plot_allowed_region(data, output_dir))

    print("Plotting Figure 2: binding constraint map ...")
    saved.append(_plot_binding_constraint(data, output_dir))

    print("Plotting Figure 3: exclusion curve ...")
    saved.append(_plot_exclusion_curve(data, output_dir))

    print("Plotting Figure 4: Yukawa scale vs M_KK ...")
    saved.append(_plot_yukawa_vs_mkk(data, output_dir))

    print("Plotting Figure 5: worst-ratio heatmap ...")
    saved.append(_plot_heatmap(data, output_dir))

    print()
    print("Saved figures:")
    for p in saved:
        try:
            rel = p.relative_to(repo_root)
        except ValueError:
            rel = p
        print(f"  {rel}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
