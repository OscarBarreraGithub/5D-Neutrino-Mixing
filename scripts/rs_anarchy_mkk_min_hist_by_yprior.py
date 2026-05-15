"""Run B plotter: M_KK^min histogram overlaid for several Y-prior variants.

Reads several rs_anarchy run directories (each with its own draws.jsonl
and tile_summary.json), extracts M_KK_min for PDG-passing draws, and
overlays the resulting distributions. Typical usage:

    python scripts/rs_anarchy_mkk_min_hist_by_yprior.py \
        --run baseline=scan_outputs/rs_anarchy_20260507T030811 \
        --run narrow_uniform=scan_outputs/rs_anarchy_runB_narrow_uniform_<TS> \
        --run wide_uniform=scan_outputs/rs_anarchy_runB_wide_uniform_<TS> \
        --run gaussian_3sigma=scan_outputs/rs_anarchy_runB_gaussian_3sigma_<TS>
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Dict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter

REPO = Path(__file__).resolve().parents[1]
DEFAULT_OUT = REPO / "results/figures/quark"
MKK_TICKS = [1, 2, 3, 5, 10, 20, 30, 50, 100, 200]
PALETTE = {
    "baseline":         "#1f77b4",
    "narrow_uniform":   "#2ca02c",
    "wide_uniform":     "#9467bd",
    "gaussian_3sigma":  "#d62728",
}


def _load_pdg_passing(draws_path: Path) -> np.ndarray:
    out = []
    with draws_path.open() as fh:
        for line in fh:
            try:
                row = json.loads(line)
            except json.JSONDecodeError:
                continue
            if not row.get("ok") or not row.get("passes_pdg"):
                continue
            mr = row.get("max_ratio")
            if mr is None or not np.isfinite(mr) or mr <= 0:
                continue
            mkk_t_tev = float(row["M_KK_GeV"]) / 1000.0
            out.append(mkk_t_tev * math.sqrt(float(mr)))
    return np.asarray(out)


def _parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run", action="append", default=[],
                   help="<tag>=<run_dir>; repeatable.")
    p.add_argument("--out-dir", type=str, default=str(DEFAULT_OUT))
    p.add_argument("--label-tag", type=str, default="")
    return p.parse_args()


def _parse_run_specs(items) -> Dict[str, Path]:
    if not items:
        sys.stderr.write("ERROR: at least one --run TAG=DIR required.\n")
        sys.exit(1)
    out: Dict[str, Path] = {}
    for it in items:
        if "=" not in it:
            sys.stderr.write(f"ERROR: bad --run spec: {it}\n")
            sys.exit(1)
        tag, _, path = it.partition("=")
        rd = Path(path)
        if not (rd / "draws.jsonl").exists():
            sys.stderr.write(f"ERROR: draws.jsonl not found in {rd}\n")
            sys.exit(1)
        out[tag] = rd
    return out


def main():
    args = _parse_args()
    runs = _parse_run_specs(args.run)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    per_run: Dict[str, np.ndarray] = {}
    for tag, run_dir in runs.items():
        draws = run_dir / "draws.jsonl"
        print(f"loading {tag} -> {draws}", flush=True)
        per_run[tag] = _load_pdg_passing(draws)
        print(f"  {tag}: {per_run[tag].size:,} PDG-passing draws")

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.7))
    bins = np.logspace(np.log10(0.3), np.log10(200.0), 70)

    # Panel 1: density
    ax = axes[0]
    for tag, arr in per_run.items():
        if arr.size == 0:
            continue
        color = PALETTE.get(tag)
        ax.hist(arr, bins=bins, density=True, histtype="step", linewidth=2.0,
                color=color, label=f"{tag} ({arr.size:,})")
    ax.set_xscale("log")
    ax.xaxis.set_major_locator(FixedLocator(MKK_TICKS))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:g}"))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel(r"$M_{\mathrm{KK}}^{\min}$ per draw  [TeV]")
    ax.set_ylabel("probability density")
    ax.set_title("Per-draw $M_{\\mathrm{KK}}^{\\min}$ vs Y-prior")
    ax.legend(fontsize=10, loc="upper right", title="Y prior")
    ax.grid(True, which="both", alpha=0.25, linestyle=":")

    # Panel 2: cumulative
    ax = axes[1]
    for tag, arr in per_run.items():
        if arr.size == 0:
            continue
        s = np.sort(arr)
        cdf = (np.arange(1, s.size + 1) / s.size) * 100.0
        color = PALETTE.get(tag)
        ax.plot(s, cdf, linewidth=2.0, color=color,
                label=f"{tag} ({arr.size:,})")
    ax.axhline(50.0, color="k", linewidth=0.7, alpha=0.5, linestyle=":")
    ax.axhline(95.0, color="k", linewidth=0.7, alpha=0.5, linestyle=":")
    ax.set_xscale("log")
    ax.xaxis.set_major_locator(FixedLocator(MKK_TICKS))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:g}"))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel(r"$M_{\mathrm{KK}}$  [TeV]")
    ax.set_ylabel(r"fraction with $M_{\mathrm{KK}}^{\min}\leq M_{\mathrm{KK}}$ [%]")
    ax.set_title("Cumulative distribution")
    ax.set_xlim(0.5, 250.0)
    ax.set_ylim(0, 102)
    ax.legend(fontsize=10, loc="lower right", title="Y prior")
    ax.grid(True, which="both", alpha=0.25, linestyle=":")

    suptitle = "RS-anarchy ensemble: $M_{\\mathrm{KK}}^{\\min}$ vs Y-prior choice"
    if args.label_tag:
        suptitle = f"{suptitle}  [{args.label_tag}]"
    fig.suptitle(suptitle, fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])

    stem = "rs_anarchy_mkk_min_hist_by_yprior"
    if args.label_tag:
        stem = f"{stem}_{args.label_tag}"
    for ext in ("pdf", "png"):
        path = out_dir / f"{stem}.{ext}"
        fig.savefig(path, dpi=200)
        print(f"wrote {path}")

    print()
    print("=== Median / p95 (TeV) ===")
    for tag, arr in per_run.items():
        if arr.size == 0:
            continue
        print(f"  {tag:18s}: median {np.percentile(arr, 50):6.2f}  "
              f"p95 {np.percentile(arr, 95):6.2f}   (n={arr.size:,})")


if __name__ == "__main__":
    main()
