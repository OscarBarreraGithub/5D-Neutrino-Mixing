#!/usr/bin/env python3
"""Comprehensive physics analysis and plotting for quark-sector scan results.

Reads merged JSONL results from the quark-sector scan and produces deep
physics analysis plots examining constraint dominance, CP violation
structure, exclusion contours, and Yukawa scale dependence.

Usage
-----
python scripts/analyze_quark_scan.py results.jsonl [--compare results2.jsonl] [--output-dir results/figures/analysis]
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import warnings
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker
import numpy as np
from matplotlib.patches import Patch
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

SYSTEM_IDS = ["epsilon_K", "K", "B_d", "B_s", "D0"]
SYSTEM_COLORS = {
    "epsilon_K": "#7B2FBE",  # purple
    "K": "#2176FF",          # blue
    "B_d": "#FF8C00",        # orange
    "B_s": "#D7263D",        # red
    "D0": "#2CA02C",         # green
}
SYSTEM_LABELS = {
    "epsilon_K": r"$\epsilon_K$",
    "K": r"$\Delta M_K$",
    "B_d": r"$\Delta M_{B_d}$",
    "B_s": r"$\Delta M_{B_s}$",
    "D0": r"$\Delta M_{D^0}$",
}
SYSTEM_LINESTYLES = {
    "epsilon_K": "-",
    "K": "--",
    "B_d": "-.",
    "B_s": ":",
    "D0": (0, (3, 1, 1, 1, 1, 1)),
}

# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------


def load_jsonl(path: str | Path) -> list[dict[str, Any]]:
    """Load a JSONL file, skipping blank lines."""
    rows: list[dict[str, Any]] = []
    with open(path, encoding="utf-8") as f:
        for lineno, line in enumerate(f, start=1):
            stripped = line.strip()
            if not stripped:
                continue
            try:
                row = json.loads(stripped)
            except json.JSONDecodeError as exc:
                print(f"WARNING: skipping malformed JSON on line {lineno}: {exc}",
                      file=sys.stderr)
                continue
            rows.append(row)
    return rows


# ---------------------------------------------------------------------------
# Data wrangling helpers
# ---------------------------------------------------------------------------


def extract_arrays(rows: list[dict[str, Any]]) -> dict[str, np.ndarray]:
    """Convert list-of-dicts into parallel numpy arrays."""
    n = len(rows)
    r = np.empty(n)
    mkk = np.empty(n)
    scale = np.empty(n)
    accepted = np.empty(n, dtype=bool)
    fit_score = np.empty(n)
    ratios = {sid: np.full(n, np.nan) for sid in SYSTEM_IDS}

    for i, row in enumerate(rows):
        r[i] = row["r"]
        mkk[i] = row.get("M_KK", row.get("Lambda_IR", np.nan))
        scale[i] = row.get("overall_scale", 1.0)
        accepted[i] = bool(row.get("accepted", False))
        fit_score[i] = row.get("fit_score", np.nan)
        rb = row.get("ratio_to_bound_by_system", {})
        for sid in SYSTEM_IDS:
            val = rb.get(sid, np.nan)
            if val is not None:
                ratios[sid][i] = float(val)

    out = dict(r=r, mkk=mkk, scale=scale, accepted=accepted,
               fit_score=fit_score, **{f"ratio_{sid}": ratios[sid]
                                       for sid in SYSTEM_IDS})
    # Filter out points with poor fit quality (score > 0.1 means quark
    # masses/CKM are not reproduced, so Wilson coefficients are unreliable)
    good = fit_score < 0.1
    if not np.all(good):
        n_bad = int((~good).sum())
        print(f"  Filtered {n_bad} point(s) with fit_score > 0.1")
        out = {k: v[good] for k, v in out.items()}
    return out


def aggregate_over_scale(data: dict[str, np.ndarray],
                         agg: str = "min_max_ratio") -> dict[str, np.ndarray]:
    """Aggregate scan results to unique (r, M_KK) grid points.

    For each unique (r, M_KK) pair, selects the overall_scale that gives
    the smallest max-ratio (the 'best' point). Returns one entry per
    unique (r, M_KK).
    """
    r_vals = data["r"]
    mkk_vals = data["mkk"]

    # Build an index keyed on (r, mkk) rounded to avoid floating noise
    grid_map: dict[tuple[float, float], list[int]] = {}
    for i in range(len(r_vals)):
        key = (round(r_vals[i], 12), round(mkk_vals[i], 6))
        grid_map.setdefault(key, []).append(i)

    n_unique = len(grid_map)
    out: dict[str, np.ndarray] = {
        "r": np.empty(n_unique),
        "mkk": np.empty(n_unique),
        "scale": np.empty(n_unique),
        "accepted": np.empty(n_unique, dtype=bool),
        "fit_score": np.empty(n_unique),
    }
    for sid in SYSTEM_IDS:
        out[f"ratio_{sid}"] = np.empty(n_unique)

    for j, ((rv, mv), indices) in enumerate(sorted(grid_map.items())):
        # For each (r, mkk), take the scale that minimises max-ratio
        best_idx = None
        best_max_ratio = np.inf
        for idx in indices:
            max_ratio = max(data[f"ratio_{sid}"][idx] for sid in SYSTEM_IDS
                           if np.isfinite(data[f"ratio_{sid}"][idx]))
            if max_ratio < best_max_ratio:
                best_max_ratio = max_ratio
                best_idx = idx
        if best_idx is None:
            best_idx = indices[0]
        out["r"][j] = data["r"][best_idx]
        out["mkk"][j] = data["mkk"][best_idx]
        out["scale"][j] = data["scale"][best_idx]
        out["accepted"][j] = data["accepted"][best_idx]
        out["fit_score"][j] = data["fit_score"][best_idx]
        for sid in SYSTEM_IDS:
            out[f"ratio_{sid}"][j] = data[f"ratio_{sid}"][best_idx]

    return out


def _make_grid(log_r: np.ndarray, log_mkk: np.ndarray,
               n_r: int = 200, n_mkk: int = 200):
    """Return meshgrid arrays for interpolation."""
    r_min, r_max = log_r.min(), log_r.max()
    m_min, m_max = log_mkk.min(), log_mkk.max()
    pad_r = 0.05 * (r_max - r_min) if r_max > r_min else 0.1
    pad_m = 0.05 * (m_max - m_min) if m_max > m_min else 0.1
    ri = np.linspace(r_min - pad_r, r_max + pad_r, n_r)
    mi = np.linspace(m_min - pad_m, m_max + pad_m, n_mkk)
    return np.meshgrid(ri, mi)


def _interp(log_r: np.ndarray, log_mkk: np.ndarray, values: np.ndarray,
            grid_r: np.ndarray, grid_m: np.ndarray,
            method: str = "cubic") -> np.ndarray:
    """Interpolate scattered data onto a regular grid with fallback."""
    finite = np.isfinite(values)
    if finite.sum() < 4:
        method = "nearest"
    pts = np.column_stack([log_r[finite], log_mkk[finite]])
    vals = values[finite]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = griddata(pts, vals, (grid_r, grid_m), method=method)
    # Fill remaining NaN with nearest-neighbor
    mask_nan = np.isnan(result)
    if mask_nan.any():
        fill = griddata(pts, vals, (grid_r, grid_m), method="nearest")
        result[mask_nan] = fill[mask_nan]
    return result


# ---------------------------------------------------------------------------
# Figure 1: Kaon vs D0 constraint dominance map
# ---------------------------------------------------------------------------


def figure1_dominance(agg: dict[str, np.ndarray], output_dir: Path) -> None:
    """Kaon vs D0 constraint dominance map in the (r, M_KK) plane."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    fig.suptitle("Figure 1: Kaon vs D$^0$ Constraint Dominance", fontsize=15,
                 fontweight="bold")

    log_r = np.log10(agg["r"])
    log_mkk = np.log10(agg["mkk"] / 1e3)  # TeV
    grid_r, grid_m = _make_grid(log_r, log_mkk)

    # --- Panel A: dominant system ---
    ax = axes[0, 0]
    ax.set_title("(A) Binding constraint", fontsize=12)
    # Assign integer label to the system with highest ratio
    sys_index = {sid: k for k, sid in enumerate(SYSTEM_IDS)}
    dominant = np.empty(len(log_r))
    for i in range(len(log_r)):
        best_sys = max(SYSTEM_IDS,
                       key=lambda s: agg[f"ratio_{s}"][i]
                       if np.isfinite(agg[f"ratio_{s}"][i]) else -np.inf)
        dominant[i] = sys_index[best_sys]

    dom_grid = _interp(log_r, log_mkk, dominant, grid_r, grid_m,
                       method="nearest")
    cmap = mcolors.ListedColormap([SYSTEM_COLORS[s] for s in SYSTEM_IDS])
    bounds = np.arange(-0.5, len(SYSTEM_IDS) + 0.5, 1)
    norm = mcolors.BoundaryNorm(bounds, cmap.N)
    im = ax.pcolormesh(grid_r, grid_m, dom_grid, cmap=cmap, norm=norm,
                       shading="auto", rasterized=True)
    patches = [Patch(facecolor=SYSTEM_COLORS[s], label=SYSTEM_LABELS[s])
               for s in SYSTEM_IDS]
    ax.legend(handles=patches, fontsize=8, loc="upper left",
              framealpha=0.9)
    ax.set_xlabel(r"$\log_{10} r$")
    ax.set_ylabel(r"$\log_{10}(M_{\mathrm{KK}}/\mathrm{TeV})$")

    # --- Panel B: kaon-only exclusion ---
    ax = axes[0, 1]
    ax.set_title(r"(B) Kaon-only exclusion ($\epsilon_K + \Delta M_K$)",
                 fontsize=11)
    kaon_max = np.maximum(
        np.where(np.isfinite(agg["ratio_epsilon_K"]), agg["ratio_epsilon_K"], 0),
        np.where(np.isfinite(agg["ratio_K"]), agg["ratio_K"], 0),
    )
    kaon_grid = _interp(log_r, log_mkk, np.log10(np.clip(kaon_max, 1e-6, None)),
                        grid_r, grid_m)
    kaon_grid_smooth = gaussian_filter(kaon_grid, sigma=1.5)
    cf = ax.contourf(grid_r, grid_m, kaon_grid_smooth,
                     levels=np.linspace(-4, 2, 25),
                     cmap="RdYlGn_r", extend="both")
    ax.contour(grid_r, grid_m, kaon_grid_smooth, levels=[0.0],
               colors="k", linewidths=2.5)
    plt.colorbar(cf, ax=ax, label=r"$\log_{10}(\max\;r_{\mathrm{kaon}})$")
    ax.set_xlabel(r"$\log_{10} r$")
    ax.set_ylabel(r"$\log_{10}(M_{\mathrm{KK}}/\mathrm{TeV})$")

    # --- Panel C: D0-only exclusion ---
    ax = axes[1, 0]
    ax.set_title(r"(C) $D^0$-only exclusion", fontsize=12)
    d0_ratio = np.where(np.isfinite(agg["ratio_D0"]), agg["ratio_D0"], 1e-6)
    d0_grid = _interp(log_r, log_mkk, np.log10(np.clip(d0_ratio, 1e-6, None)),
                      grid_r, grid_m)
    d0_grid_smooth = gaussian_filter(d0_grid, sigma=1.5)
    cf = ax.contourf(grid_r, grid_m, d0_grid_smooth,
                     levels=np.linspace(-4, 2, 25),
                     cmap="RdYlGn_r", extend="both")
    ax.contour(grid_r, grid_m, d0_grid_smooth, levels=[0.0],
               colors="k", linewidths=2.5)
    plt.colorbar(cf, ax=ax, label=r"$\log_{10}(r_{D^0})$")
    ax.set_xlabel(r"$\log_{10} r$")
    ax.set_ylabel(r"$\log_{10}(M_{\mathrm{KK}}/\mathrm{TeV})$")

    # --- Panel D: combined exclusion ---
    ax = axes[1, 1]
    ax.set_title("(D) Combined exclusion (all 5 systems)", fontsize=12)
    all_max = np.zeros(len(log_r))
    for sid in SYSTEM_IDS:
        vals = agg[f"ratio_{sid}"]
        all_max = np.maximum(all_max,
                             np.where(np.isfinite(vals), vals, 0))
    all_grid = _interp(log_r, log_mkk, np.log10(np.clip(all_max, 1e-6, None)),
                       grid_r, grid_m)
    all_grid_smooth = gaussian_filter(all_grid, sigma=1.5)
    cf = ax.contourf(grid_r, grid_m, all_grid_smooth,
                     levels=np.linspace(-4, 2, 25),
                     cmap="RdYlGn_r", extend="both")
    ax.contour(grid_r, grid_m, all_grid_smooth, levels=[0.0],
               colors="k", linewidths=2.5)
    plt.colorbar(cf, ax=ax, label=r"$\log_{10}(\max\;r_{\mathrm{all}})$")
    ax.set_xlabel(r"$\log_{10} r$")
    ax.set_ylabel(r"$\log_{10}(M_{\mathrm{KK}}/\mathrm{TeV})$")

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    outpath = output_dir / "fig1_dominance_map.pdf"
    fig.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {outpath}")


# ---------------------------------------------------------------------------
# Figure 2: CP violation structure
# ---------------------------------------------------------------------------


def figure2_cp_violation(agg: dict[str, np.ndarray], output_dir: Path) -> None:
    """CP violation structure in the (r, M_KK) plane."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    fig.suptitle("Figure 2: CP Violation Structure", fontsize=15,
                 fontweight="bold")

    log_r = np.log10(agg["r"])
    log_mkk = np.log10(agg["mkk"] / 1e3)
    grid_r, grid_m = _make_grid(log_r, log_mkk)

    # --- Panel A: CP phase proxy = epsilon_K / Delta M_K ratio ---
    ax = axes[0]
    ax.set_title(r"(A) CP phase proxy: $r_{\epsilon_K} / r_{\Delta M_K}$",
                 fontsize=11)
    ek = agg["ratio_epsilon_K"]
    dk = agg["ratio_K"]
    # Avoid division by zero
    safe_dk = np.where((np.isfinite(dk)) & (dk > 0), dk, np.nan)
    cp_proxy = np.where(np.isfinite(ek) & np.isfinite(safe_dk),
                        ek / safe_dk, np.nan)
    log_cp = np.log10(np.clip(cp_proxy, 1e-6, 1e6))
    cp_grid = _interp(log_r, log_mkk, log_cp, grid_r, grid_m)
    cp_grid_smooth = gaussian_filter(np.nan_to_num(cp_grid, nan=0.0), sigma=1.5)
    cf = ax.contourf(grid_r, grid_m, cp_grid_smooth,
                     levels=np.linspace(-3, 3, 25), cmap="coolwarm",
                     extend="both")
    plt.colorbar(cf, ax=ax,
                 label=r"$\log_{10}(r_{\epsilon_K}/r_{\Delta M_K})$")
    ax.set_xlabel(r"$\log_{10} r$")
    ax.set_ylabel(r"$\log_{10}(M_{\mathrm{KK}}/\mathrm{TeV})$")

    # --- Panel B: epsilon_K contour ---
    ax = axes[1]
    ax.set_title(r"(B) $\epsilon_K$ ratio-to-bound", fontsize=12)
    ek_log = np.log10(np.clip(agg["ratio_epsilon_K"], 1e-8, None))
    ek_grid = _interp(log_r, log_mkk, ek_log, grid_r, grid_m)
    ek_grid_smooth = gaussian_filter(ek_grid, sigma=1.5)
    cf = ax.contourf(grid_r, grid_m, ek_grid_smooth,
                     levels=np.linspace(-6, 2, 30),
                     cmap="magma_r", extend="both")
    cs = ax.contour(grid_r, grid_m, ek_grid_smooth, levels=[0.0],
                    colors="lime", linewidths=2.5)
    ax.clabel(cs, fmt=r"$r=1$", fontsize=9)
    plt.colorbar(cf, ax=ax,
                 label=r"$\log_{10}(r_{\epsilon_K})$")
    ax.set_xlabel(r"$\log_{10} r$")
    ax.set_ylabel(r"$\log_{10}(M_{\mathrm{KK}}/\mathrm{TeV})$")

    # --- Panel C: Delta M_K contour ---
    ax = axes[2]
    ax.set_title(r"(C) $\Delta M_K$ ratio-to-bound", fontsize=12)
    dk_log = np.log10(np.clip(agg["ratio_K"], 1e-8, None))
    dk_grid = _interp(log_r, log_mkk, dk_log, grid_r, grid_m)
    dk_grid_smooth = gaussian_filter(dk_grid, sigma=1.5)
    cf = ax.contourf(grid_r, grid_m, dk_grid_smooth,
                     levels=np.linspace(-6, 2, 30),
                     cmap="magma_r", extend="both")
    cs = ax.contour(grid_r, grid_m, dk_grid_smooth, levels=[0.0],
                    colors="lime", linewidths=2.5)
    ax.clabel(cs, fmt=r"$r=1$", fontsize=9)
    plt.colorbar(cf, ax=ax,
                 label=r"$\log_{10}(r_{\Delta M_K})$")
    ax.set_xlabel(r"$\log_{10} r$")
    ax.set_ylabel(r"$\log_{10}(M_{\mathrm{KK}}/\mathrm{TeV})$")

    fig.tight_layout(rect=[0, 0, 1, 0.93])
    outpath = output_dir / "fig2_cp_violation.pdf"
    fig.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {outpath}")


# ---------------------------------------------------------------------------
# Figure 3: Smooth exclusion contours (max ratio)
# ---------------------------------------------------------------------------


def figure3_exclusion_contours(data: dict[str, np.ndarray],
                               agg: dict[str, np.ndarray],
                               output_dir: Path) -> None:
    """Smooth interpolated exclusion contours in (log r, log M_KK/TeV)."""

    # --- Per-scale panels + combined best-scale panel ---
    unique_scales = np.unique(data["scale"])
    n_panels = len(unique_scales) + 1
    ncols = min(n_panels, 3)
    nrows = (n_panels + ncols - 1) // ncols
    fig, axes_flat = plt.subplots(nrows, ncols, figsize=(6 * ncols, 5 * nrows),
                                  squeeze=False)
    axes = axes_flat.flatten()
    fig.suptitle("Figure 3: Exclusion Contours (max ratio-to-bound)",
                 fontsize=15, fontweight="bold")

    # Diverging colormap: green (allowed) to red (excluded)
    cmap = mcolors.TwoSlopeNorm(vcenter=0.0, vmin=-2.0, vmax=2.0)

    for panel_idx, s_val in enumerate(unique_scales):
        ax = axes[panel_idx]
        mask = data["scale"] == s_val
        lr = np.log10(data["r"][mask])
        lm = np.log10(data["mkk"][mask] / 1e3)
        max_rat = np.zeros(mask.sum())
        for sid in SYSTEM_IDS:
            v = data[f"ratio_{sid}"][mask]
            max_rat = np.maximum(max_rat, np.where(np.isfinite(v), v, 0))
        log_max = np.log10(np.clip(max_rat, 1e-6, None))

        gr, gm = _make_grid(lr, lm)
        Z = _interp(lr, lm, log_max, gr, gm)
        Z_smooth = gaussian_filter(Z, sigma=1.5)

        cf = ax.contourf(gr, gm, Z_smooth, levels=np.linspace(-2, 2, 40),
                         cmap="RdYlGn_r", extend="both")
        ax.contour(gr, gm, Z_smooth, levels=[0.0], colors="k",
                   linewidths=3, linestyles="-")
        plt.colorbar(cf, ax=ax, label=r"$\log_{10}(\max\;r)$")
        ax.set_title(f"overall_scale = {s_val:.2g}", fontsize=11)
        ax.set_xlabel(r"$\log_{10} r$")
        ax.set_ylabel(r"$\log_{10}(M_{\mathrm{KK}}/\mathrm{TeV})$")

    # Combined panel: best scale
    ax = axes[n_panels - 1]
    log_r = np.log10(agg["r"])
    log_mkk = np.log10(agg["mkk"] / 1e3)
    max_rat_agg = np.zeros(len(log_r))
    for sid in SYSTEM_IDS:
        v = agg[f"ratio_{sid}"]
        max_rat_agg = np.maximum(max_rat_agg, np.where(np.isfinite(v), v, 0))
    log_max_agg = np.log10(np.clip(max_rat_agg, 1e-6, None))

    gr, gm = _make_grid(log_r, log_mkk)
    Z = _interp(log_r, log_mkk, log_max_agg, gr, gm)
    Z_smooth = gaussian_filter(Z, sigma=1.5)

    cf = ax.contourf(gr, gm, Z_smooth, levels=np.linspace(-2, 2, 40),
                     cmap="RdYlGn_r", extend="both")
    ax.contour(gr, gm, Z_smooth, levels=[0.0], colors="k",
               linewidths=3, linestyles="-")
    plt.colorbar(cf, ax=ax, label=r"$\log_{10}(\max\;r)$")
    ax.set_title("Best overall_scale", fontsize=11, fontweight="bold")
    ax.set_xlabel(r"$\log_{10} r$")
    ax.set_ylabel(r"$\log_{10}(M_{\mathrm{KK}}/\mathrm{TeV})$")

    # Hide unused panels
    for k in range(n_panels, len(axes)):
        axes[k].set_visible(False)

    fig.tight_layout(rect=[0, 0, 1, 0.94])
    outpath = output_dir / "fig3_exclusion_contours.pdf"
    fig.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {outpath}")


# ---------------------------------------------------------------------------
# Figure 4: System-by-system exclusion boundaries
# ---------------------------------------------------------------------------


def figure4_system_boundaries(agg: dict[str, np.ndarray],
                              output_dir: Path) -> None:
    """Overlay exclusion contours (ratio=1) for each meson system."""
    fig, ax = plt.subplots(figsize=(10, 8))
    fig.suptitle("Figure 4: System-by-System Exclusion Boundaries",
                 fontsize=15, fontweight="bold")

    log_r = np.log10(agg["r"])
    log_mkk = np.log10(agg["mkk"] / 1e3)
    gr, gm = _make_grid(log_r, log_mkk)

    legend_handles = []
    for sid in SYSTEM_IDS:
        vals = agg[f"ratio_{sid}"]
        log_vals = np.log10(np.clip(vals, 1e-8, None))
        Z = _interp(log_r, log_mkk, log_vals, gr, gm)
        Z_smooth = gaussian_filter(Z, sigma=1.5)
        ax.contour(gr, gm, Z_smooth, levels=[0.0],
                   colors=[SYSTEM_COLORS[sid]],
                   linewidths=2.5,
                   linestyles=[SYSTEM_LINESTYLES[sid]])
        # Use a proxy Line2D for the legend
        legend_handles.append(
            plt.Line2D([], [], color=SYSTEM_COLORS[sid],
                       linewidth=2.5, linestyle=SYSTEM_LINESTYLES[sid],
                       label=SYSTEM_LABELS[sid])
        )

    # Shade the combined allowed region lightly
    max_rat = np.zeros(len(log_r))
    for sid in SYSTEM_IDS:
        v = agg[f"ratio_{sid}"]
        max_rat = np.maximum(max_rat, np.where(np.isfinite(v), v, 0))
    log_max = np.log10(np.clip(max_rat, 1e-6, None))
    Z_all = _interp(log_r, log_mkk, log_max, gr, gm)
    Z_all_smooth = gaussian_filter(Z_all, sigma=1.5)
    ax.contourf(gr, gm, Z_all_smooth, levels=[-10, 0.0],
                colors=["#d4edda"], alpha=0.3)
    ax.contour(gr, gm, Z_all_smooth, levels=[0.0],
               colors="k", linewidths=3, linestyles="-")
    legend_handles.append(
        plt.Line2D([], [], color="k", linewidth=3, linestyle="-",
                   label="Combined")
    )

    ax.legend(handles=legend_handles, fontsize=11, loc="upper left",
              framealpha=0.9)
    ax.set_xlabel(r"$\log_{10} r$", fontsize=13)
    ax.set_ylabel(r"$\log_{10}(M_{\mathrm{KK}}/\mathrm{TeV})$", fontsize=13)
    ax.grid(True, alpha=0.3)

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    outpath = output_dir / "fig4_system_boundaries.pdf"
    fig.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {outpath}")


# ---------------------------------------------------------------------------
# Figure 5: M_KK lower bound vs r with error band
# ---------------------------------------------------------------------------


def figure5_mkk_bound_vs_r(data: dict[str, np.ndarray],
                           output_dir: Path) -> None:
    """M_KK lower bound vs r, with band from overall_scale variation."""
    fig, ax = plt.subplots(figsize=(10, 7))
    fig.suptitle(r"Figure 5: $M_{\mathrm{KK}}$ Lower Bound vs $r$",
                 fontsize=15, fontweight="bold")

    unique_r = np.sort(np.unique(np.round(data["r"], 12)))
    unique_scales = np.sort(np.unique(data["scale"]))

    # For each r and each scale, find the minimum M_KK that is accepted.
    # If none accepted, find where max_ratio crosses 1 by interpolation.
    bound_by_scale: dict[float, np.ndarray] = {}
    for s_val in unique_scales:
        bounds = np.full(len(unique_r), np.nan)
        for j, rv in enumerate(unique_r):
            mask = (np.round(data["r"], 12) == rv) & (data["scale"] == s_val)
            if not mask.any():
                continue
            mkk_vals = data["mkk"][mask]
            max_rats = np.zeros(mask.sum())
            for sid in SYSTEM_IDS:
                v = data[f"ratio_{sid}"][mask]
                max_rats = np.maximum(max_rats,
                                      np.where(np.isfinite(v), v, 0))
            # Sort by M_KK
            order = np.argsort(mkk_vals)
            mkk_sorted = mkk_vals[order]
            rat_sorted = max_rats[order]

            # Find the lowest M_KK where ratio <= 1
            accepted_mask = rat_sorted <= 1.0
            if accepted_mask.any():
                bounds[j] = mkk_sorted[accepted_mask][0]
            else:
                # Interpolate: find where ratio crosses 1 if it decreases
                # with M_KK (as expected physically)
                if len(rat_sorted) >= 2:
                    for k in range(len(rat_sorted) - 1):
                        if rat_sorted[k] > 1.0 > rat_sorted[k + 1]:
                            # Linear interpolation in log space
                            frac = ((np.log10(1.0) - np.log10(rat_sorted[k]))
                                    / (np.log10(rat_sorted[k + 1])
                                       - np.log10(rat_sorted[k])))
                            bounds[j] = 10 ** (np.log10(mkk_sorted[k])
                                               + frac * (np.log10(mkk_sorted[k + 1])
                                                         - np.log10(mkk_sorted[k])))
                            break
        bound_by_scale[s_val] = bounds

    # Stack all scale bounds
    all_bounds = np.array(list(bound_by_scale.values()))  # (n_scales, n_r)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        median_bound = np.nanmedian(all_bounds, axis=0)
        min_bound = np.nanmin(all_bounds, axis=0)
        max_bound = np.nanmax(all_bounds, axis=0)

    valid = np.isfinite(median_bound)
    if valid.any():
        ax.fill_between(unique_r[valid], min_bound[valid] / 1e3,
                        max_bound[valid] / 1e3,
                        alpha=0.25, color="steelblue",
                        label="Scale variation band")
        ax.plot(unique_r[valid], median_bound[valid] / 1e3,
                "o-", color="steelblue", linewidth=2, markersize=5,
                label="Median bound")

        # Also plot individual scales
        for s_val, bnd in bound_by_scale.items():
            v = np.isfinite(bnd)
            if v.any():
                ax.plot(unique_r[v], bnd[v] / 1e3, "x--", alpha=0.4,
                        label=f"scale={s_val:.2g}", markersize=4)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$r$", fontsize=13)
    ax.set_ylabel(r"$M_{\mathrm{KK}}^{\min}$ [TeV]", fontsize=13)
    ax.legend(fontsize=9, loc="best", framealpha=0.9)
    ax.grid(True, which="both", alpha=0.3)

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    outpath = output_dir / "fig5_mkk_bound_vs_r.pdf"
    fig.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {outpath}")


# ---------------------------------------------------------------------------
# Figure 6: 2007-era vs modern bounds (comparison)
# ---------------------------------------------------------------------------


def figure6_comparison(agg1: dict[str, np.ndarray],
                       agg2: dict[str, np.ndarray],
                       data1: dict[str, np.ndarray],
                       data2: dict[str, np.ndarray],
                       label1: str, label2: str,
                       output_dir: Path) -> None:
    """Compare two scan results: exclusion contours, M_KK bounds, scatter."""
    fig, axes = plt.subplots(1, 3, figsize=(20, 7))
    fig.suptitle("Figure 6: Comparison of Constraint Inputs",
                 fontsize=15, fontweight="bold")

    # Helper to compute max ratio
    def _max_ratio(d: dict[str, np.ndarray]) -> np.ndarray:
        mr = np.zeros(len(d["r"]))
        for sid in SYSTEM_IDS:
            v = d[f"ratio_{sid}"]
            mr = np.maximum(mr, np.where(np.isfinite(v), v, 0))
        return mr

    # --- Panel A: exclusion contours ---
    ax = axes[0]
    ax.set_title("(A) Exclusion contours", fontsize=12)

    comp_handles = []
    for agg, ls, lbl, col in [(agg1, "-", label1, "navy"),
                               (agg2, "--", label2, "firebrick")]:
        lr = np.log10(agg["r"])
        lm = np.log10(agg["mkk"] / 1e3)
        mr = _max_ratio(agg)
        log_mr = np.log10(np.clip(mr, 1e-6, None))
        gr, gm = _make_grid(lr, lm)
        Z = _interp(lr, lm, log_mr, gr, gm)
        Z_smooth = gaussian_filter(Z, sigma=1.5)
        ax.contour(gr, gm, Z_smooth, levels=[0.0],
                   colors=[col], linewidths=2.5, linestyles=[ls])
        comp_handles.append(
            plt.Line2D([], [], color=col, linewidth=2.5, linestyle=ls,
                       label=lbl)
        )

    ax.legend(handles=comp_handles, fontsize=10)
    ax.set_xlabel(r"$\log_{10} r$", fontsize=12)
    ax.set_ylabel(r"$\log_{10}(M_{\mathrm{KK}}/\mathrm{TeV})$", fontsize=12)
    ax.grid(True, alpha=0.3)

    # --- Panel B: ratio of M_KK lower bounds ---
    ax = axes[1]
    ax.set_title(r"(B) $M_{\mathrm{KK}}^{\min}$ ratio", fontsize=12)

    for dat, lbl, col in [(data1, label1, "navy"), (data2, label2, "firebrick")]:
        unique_r = np.sort(np.unique(np.round(dat["r"], 12)))
        bounds = np.full(len(unique_r), np.nan)
        for j, rv in enumerate(unique_r):
            mask = np.round(dat["r"], 12) == rv
            if not mask.any():
                continue
            mkk_vals = dat["mkk"][mask]
            mr = _max_ratio({k: dat[k][mask] for k in dat})
            order = np.argsort(mkk_vals)
            acc = mr[order] <= 1.0
            if acc.any():
                bounds[j] = mkk_vals[order][acc][0]
        valid = np.isfinite(bounds)
        if valid.any():
            ax.plot(unique_r[valid], bounds[valid] / 1e3, "o-",
                    color=col, linewidth=2, label=lbl)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$r$", fontsize=12)
    ax.set_ylabel(r"$M_{\mathrm{KK}}^{\min}$ [TeV]", fontsize=12)
    ax.legend(fontsize=10)
    ax.grid(True, which="both", alpha=0.3)

    # --- Panel C: scatter overlay ---
    ax = axes[2]
    ax.set_title("(C) Yukawa scale vs $M_{\\mathrm{KK}}$", fontsize=12)

    for dat, lbl, col, marker in [(data1, label1, "navy", "o"),
                                   (data2, label2, "firebrick", "s")]:
        mr = _max_ratio(dat)
        accepted = mr <= 1.0
        ax.scatter(dat["mkk"][~accepted] / 1e3, dat["scale"][~accepted],
                   c=col, alpha=0.15, s=8, marker=marker)
        ax.scatter(dat["mkk"][accepted] / 1e3, dat["scale"][accepted],
                   c=col, alpha=0.8, s=25, marker=marker, edgecolors="k",
                   linewidths=0.5, label=f"{lbl} (accepted)")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$M_{\mathrm{KK}}$ [TeV]", fontsize=12)
    ax.set_ylabel("overall_scale", fontsize=12)
    ax.legend(fontsize=9)
    ax.grid(True, which="both", alpha=0.3)

    fig.tight_layout(rect=[0, 0, 1, 0.94])
    outpath = output_dir / "fig6_comparison.pdf"
    fig.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {outpath}")


# ---------------------------------------------------------------------------
# Figure 7: Yukawa scale dependence
# ---------------------------------------------------------------------------


def figure7_yukawa_scale(data: dict[str, np.ndarray],
                         output_dir: Path) -> None:
    """Acceptance fraction vs overall_scale, binned by M_KK."""
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    fig.suptitle("Figure 7: Yukawa Scale Dependence of Acceptance",
                 fontsize=15, fontweight="bold")

    unique_scales = np.sort(np.unique(data["scale"]))
    unique_mkk = np.sort(np.unique(np.round(data["mkk"], 6)))

    # Max ratio for each point
    max_rat = np.zeros(len(data["r"]))
    for sid in SYSTEM_IDS:
        v = data[f"ratio_{sid}"]
        max_rat = np.maximum(max_rat, np.where(np.isfinite(v), v, 0))

    # --- Panel A: acceptance fraction vs scale, one curve per M_KK bin ---
    ax = axes[0]
    ax.set_title("(A) Acceptance fraction vs Yukawa scale", fontsize=12)

    if len(unique_scales) > 1:
        n_mkk_bins = min(6, len(unique_mkk))
        mkk_bins = np.array_split(unique_mkk, n_mkk_bins)
        cmap_lines = plt.cm.viridis(np.linspace(0.1, 0.9, n_mkk_bins))

        for bi, mkk_group in enumerate(mkk_bins):
            acc_frac = []
            for sv in unique_scales:
                mask = ((data["scale"] == sv)
                        & np.isin(np.round(data["mkk"], 6), mkk_group))
                if mask.sum() == 0:
                    acc_frac.append(np.nan)
                else:
                    acc_frac.append(np.mean(max_rat[mask] <= 1.0))
            lo = mkk_group[0] / 1e3
            hi = mkk_group[-1] / 1e3
            ax.plot(unique_scales, acc_frac, "o-", color=cmap_lines[bi],
                    linewidth=2, markersize=5,
                    label=f"$M_{{KK}}$: {lo:.1f}-{hi:.1f} TeV")
    else:
        # Only one scale: show acceptance vs M_KK instead
        acc_frac_by_mkk = []
        for mv in unique_mkk:
            mask = np.round(data["mkk"], 6) == mv
            if mask.sum() == 0:
                acc_frac_by_mkk.append(np.nan)
            else:
                acc_frac_by_mkk.append(np.mean(max_rat[mask] <= 1.0))
        ax.plot(unique_mkk / 1e3, acc_frac_by_mkk, "o-", color="steelblue",
                linewidth=2, markersize=6,
                label="Acceptance (single scale)")
        ax.set_xlabel(r"$M_{\mathrm{KK}}$ [TeV]", fontsize=12)
        ax.set_xscale("log")

    if len(unique_scales) > 1:
        ax.set_xlabel("overall_scale", fontsize=12)
        ax.set_xscale("log")
    ax.set_ylabel("Acceptance fraction", fontsize=12)
    ax.set_ylim(-0.05, 1.05)
    # Only show legend when there are labeled artists
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(fontsize=8, loc="best", framealpha=0.9)
    ax.grid(True, alpha=0.3)

    # --- Panel B: heatmap of acceptance in (scale, M_KK) ---
    ax = axes[1]
    ax.set_title("(B) Acceptance heatmap", fontsize=12)

    if len(unique_scales) > 1 and len(unique_mkk) > 1:
        acc_map = np.full((len(unique_mkk), len(unique_scales)), np.nan)
        for im, mv in enumerate(unique_mkk):
            for js, sv in enumerate(unique_scales):
                mask = ((np.round(data["mkk"], 6) == mv)
                        & (data["scale"] == sv))
                if mask.sum() > 0:
                    acc_map[im, js] = np.mean(max_rat[mask] <= 1.0)
        im_plot = ax.pcolormesh(
            np.log10(unique_scales),
            np.log10(unique_mkk / 1e3),
            acc_map, cmap="RdYlGn", vmin=0, vmax=1, shading="auto")
        plt.colorbar(im_plot, ax=ax, label="Acceptance fraction")
        ax.set_xlabel(r"$\log_{10}$(overall_scale)", fontsize=12)
        ax.set_ylabel(r"$\log_{10}(M_{\mathrm{KK}}/\mathrm{TeV})$", fontsize=12)
    elif len(unique_mkk) > 1:
        # Single scale: show max_ratio distribution vs M_KK
        for mv in unique_mkk:
            mask = np.round(data["mkk"], 6) == mv
            if mask.sum() > 0:
                ax.scatter(np.full(mask.sum(), mv / 1e3), max_rat[mask],
                           s=15, alpha=0.5)
        ax.axhline(1.0, color="k", linewidth=1.5, linestyle="--",
                   label="ratio = 1")
        ax.set_xlabel(r"$M_{\mathrm{KK}}$ [TeV]", fontsize=12)
        ax.set_ylabel("max ratio-to-bound", fontsize=12)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.legend(fontsize=10)
    else:
        ax.text(0.5, 0.5, "Insufficient variation\nfor heatmap",
                transform=ax.transAxes, ha="center", va="center",
                fontsize=14, color="gray")
    ax.grid(True, alpha=0.3)

    fig.tight_layout(rect=[0, 0, 1, 0.94])
    outpath = output_dir / "fig7_yukawa_scale.pdf"
    fig.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {outpath}")


# ---------------------------------------------------------------------------
# Summary statistics
# ---------------------------------------------------------------------------


def print_summary(data: dict[str, np.ndarray], label: str) -> None:
    """Print summary statistics to stdout."""
    n = len(data["r"])
    max_rat = np.zeros(n)
    for sid in SYSTEM_IDS:
        v = data[f"ratio_{sid}"]
        max_rat = np.maximum(max_rat, np.where(np.isfinite(v), v, 0))
    n_acc = int(np.sum(max_rat <= 1.0))

    print(f"\n{'='*60}")
    print(f"Summary for: {label}")
    print(f"{'='*60}")
    print(f"  Total points:    {n}")
    print(f"  Accepted:        {n_acc} ({100*n_acc/n:.1f}%)" if n > 0 else "  No data")
    print(f"  r range:         [{data['r'].min():.4g}, {data['r'].max():.4g}]")
    print(f"  M_KK range:      [{data['mkk'].min():.0f}, {data['mkk'].max():.0f}] GeV")
    print(f"                   [{data['mkk'].min()/1e3:.2f}, {data['mkk'].max()/1e3:.2f}] TeV")
    print(f"  Scale range:     [{data['scale'].min():.3g}, {data['scale'].max():.3g}]")

    # Dominant constraint statistics
    dom_count = {sid: 0 for sid in SYSTEM_IDS}
    for i in range(n):
        best = max(SYSTEM_IDS,
                   key=lambda s: data[f"ratio_{s}"][i]
                   if np.isfinite(data[f"ratio_{s}"][i]) else -np.inf)
        dom_count[best] += 1
    print("  Dominant constraint frequency:")
    for sid in SYSTEM_IDS:
        print(f"    {SYSTEM_LABELS[sid]:>15s}: {dom_count[sid]:5d} "
              f"({100*dom_count[sid]/n:.1f}%)" if n > 0 else "")

    # Per-system ratio statistics
    print("  Ratio-to-bound statistics (median [min, max]):")
    for sid in SYSTEM_IDS:
        v = data[f"ratio_{sid}"]
        fin = v[np.isfinite(v)]
        if len(fin) > 0:
            print(f"    {SYSTEM_LABELS[sid]:>15s}: {np.median(fin):.4g} "
                  f"[{fin.min():.4g}, {fin.max():.4g}]")
    print()


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Physics analysis plots for quark-sector scan results.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("results", type=str,
                        help="Path to the primary merged JSONL results file.")
    parser.add_argument("--compare", type=str, default=None,
                        help="Path to a second JSONL for comparison (Figure 6).")
    parser.add_argument("--output-dir", type=str,
                        default="results/figures/analysis",
                        help="Output directory for figures (default: "
                             "results/figures/analysis).")
    parser.add_argument("--label1", type=str, default=None,
                        help="Label for the primary dataset (default: filename).")
    parser.add_argument("--label2", type=str, default=None,
                        help="Label for the comparison dataset (default: filename).")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load primary dataset
    print(f"Loading {args.results} ...")
    rows1 = load_jsonl(args.results)
    if not rows1:
        print("ERROR: no data rows in primary results file.", file=sys.stderr)
        sys.exit(1)
    data1 = extract_arrays(rows1)
    agg1 = aggregate_over_scale(data1)
    label1 = args.label1 or Path(args.results).stem

    print_summary(data1, label1)

    # Generate figures 1-5, 7
    print("Generating figures ...")
    figure1_dominance(agg1, output_dir)
    figure2_cp_violation(agg1, output_dir)
    figure3_exclusion_contours(data1, agg1, output_dir)
    figure4_system_boundaries(agg1, output_dir)
    figure5_mkk_bound_vs_r(data1, output_dir)
    figure7_yukawa_scale(data1, output_dir)

    # Figure 6: comparison (if second file provided)
    if args.compare:
        print(f"\nLoading comparison file {args.compare} ...")
        rows2 = load_jsonl(args.compare)
        if not rows2:
            print("WARNING: no data in comparison file; skipping Figure 6.",
                  file=sys.stderr)
        else:
            data2 = extract_arrays(rows2)
            agg2 = aggregate_over_scale(data2)
            label2 = args.label2 or Path(args.compare).stem
            print_summary(data2, label2)
            figure6_comparison(agg1, agg2, data1, data2, label1, label2,
                               output_dir)

    print(f"\nAll figures saved to {output_dir.resolve()}")


if __name__ == "__main__":
    main()
