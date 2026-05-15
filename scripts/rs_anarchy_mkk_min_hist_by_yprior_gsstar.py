"""Y-prior sensitivity histogram in g_s* = 3 convention only.

Differs from rs_anarchy_mkk_min_hist_by_yprior.py only in that the
x-axis is rescaled to g_s* = 3 (CFW convention) — the on-disk M_KK^min
values are multiplied by 3 / 1.05 ~= 2.857. No perturbative-axis display.

Default writes:
  results/figures/quark/rs_anarchy_mkk_min_hist_by_yprior_gsstar.{pdf,png}
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter

REPO = Path(__file__).resolve().parents[1]
DEFAULT_OUT = REPO / "results/figures/quark"
MKK_TICKS = [3, 5, 10, 20, 30, 50, 100, 200, 500]

G_S_PERT = 1.05
G_S_STAR = 3.0
GS_RESCALE_MKK = G_S_STAR / G_S_PERT

PALETTE = {
    "baseline":         "C0",
    "narrow_uniform":   "C2",
    "wide_uniform":     "C3",
    "gaussian_3sigma":  "C4",
}


def _load_pdg_passing(draws_path: Path) -> np.ndarray:
    """Return per-draw M_KK^min in TeV (perturbative; rescale at plot time)."""
    out = []
    with draws_path.open() as fh:
        for line in fh:
            try:
                row = json.loads(line)
            except json.JSONDecodeError:
                continue
            if not row.get("ok") or not row.get("passes_pdg"):
                continue
            mr = float(row.get("max_ratio", float("nan")))
            if not np.isfinite(mr) or mr <= 0:
                continue
            mkk_t = float(row["M_KK_GeV"]) / 1000.0
            out.append(mkk_t * float(np.sqrt(mr)))
    return np.asarray(out, dtype=float)


def _parse_run_specs(specs: List[str]) -> Dict[str, Path]:
    runs: Dict[str, Path] = {}
    for s in specs:
        if "=" not in s:
            sys.stderr.write(f"ERROR: --run argument must be TAG=DIR, got: {s}\n")
            sys.exit(2)
        tag, raw = s.split("=", 1)
        d = Path(raw)
        if not d.exists():
            sys.stderr.write(f"ERROR: run dir does not exist: {d}\n")
            sys.exit(1)
        runs[tag] = d
    return runs


def _parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run", action="append", default=[],
                   help="TAG=DIR; pass multiple --run flags")
    p.add_argument("--out-dir", type=str, default=str(DEFAULT_OUT))
    p.add_argument("--label-tag", type=str, default="")
    return p.parse_args()


def main():
    args = _parse_args()
    if not args.run:
        sys.stderr.write("ERROR: at least one --run TAG=DIR is required.\n")
        sys.exit(2)
    runs = _parse_run_specs(args.run)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    per_run: Dict[str, np.ndarray] = {}
    for tag, run_dir in runs.items():
        draws = run_dir / "draws.jsonl"
        print(f"loading {tag} -> {draws}", flush=True)
        arr = _load_pdg_passing(draws) * GS_RESCALE_MKK  # rescale once, here
        per_run[tag] = arr
        print(f"  {tag}: {arr.size:,} PDG-passing draws  "
              f"(median {np.median(arr) if arr.size else float('nan'):.2f} TeV at g_s*=3)")

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.7))
    bins = np.logspace(np.log10(1.0), np.log10(500.0), 70)

    # Panel 1: density
    ax = axes[0]
    for tag, arr in per_run.items():
        if arr.size == 0:
            continue
        ax.hist(arr, bins=bins, density=True, histtype="step", linewidth=2.0,
                color=PALETTE.get(tag), label=f"{tag} ({arr.size:,})")
    ax.set_xscale("log")
    ax.xaxis.set_major_locator(FixedLocator(MKK_TICKS))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:g}"))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlim(2.0, 500.0)
    ax.set_xlabel(r"$M_{\mathrm{KK}}^{\min}$  [TeV]   ($g_s^\star = 3$)",
                  fontsize=12)
    ax.set_ylabel("probability density", fontsize=12)
    ax.set_title(r"Per-draw $M_{\mathrm{KK}}^{\min}$ vs Y-prior  ($g_s^\star = 3$)",
                 fontsize=12)
    ax.legend(fontsize=10, loc="upper right", title="Y prior")
    ax.grid(True, which="both", alpha=0.25, linestyle=":")

    # Panel 2: cumulative
    ax = axes[1]
    for tag, arr in per_run.items():
        if arr.size == 0:
            continue
        s = np.sort(arr)
        cdf = (np.arange(1, s.size + 1) / s.size) * 100.0
        ax.plot(s, cdf, linewidth=2.0, color=PALETTE.get(tag),
                label=f"{tag} ({arr.size:,})")
    ax.axhline(50.0, color="k", linewidth=0.7, alpha=0.5, linestyle=":")
    ax.axhline(95.0, color="k", linewidth=0.7, alpha=0.5, linestyle=":")
    ax.set_xscale("log")
    ax.xaxis.set_major_locator(FixedLocator(MKK_TICKS))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:g}"))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlim(2.0, 500.0)
    ax.set_ylim(0, 102)
    ax.set_xlabel(r"$M_{\mathrm{KK}}$  [TeV]   ($g_s^\star = 3$)", fontsize=12)
    ax.set_ylabel(r"fraction with $M_{\mathrm{KK}}^{\min}\leq M_{\mathrm{KK}}$ [%]",
                  fontsize=12)
    ax.set_title(r"Cumulative distribution  ($g_s^\star = 3$)", fontsize=12)
    ax.legend(fontsize=10, loc="lower right", title="Y prior")
    ax.grid(True, which="both", alpha=0.25, linestyle=":")

    fig.tight_layout()
    stem = "rs_anarchy_mkk_min_hist_by_yprior_gsstar"
    if args.label_tag:
        stem = f"{stem}_{args.label_tag}"
    for ext in ("pdf", "png"):
        path = out_dir / f"{stem}.{ext}"
        fig.savefig(path, dpi=200)
        print(f"  wrote {path}")


if __name__ == "__main__":
    main()
