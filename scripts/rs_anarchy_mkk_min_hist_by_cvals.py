"""Run 3 plotter: M_KK^min histogram overlaid for several c-value patterns.

Reads several rs_anarchy run directories (each with its own draws.jsonl),
extracts M_KK_min for PDG-passing draws, and overlays the resulting
distributions.

Usage:
    python scripts/rs_anarchy_mkk_min_hist_by_cvals.py \
        --run baseline=scan_outputs/rs_anarchy_run3_baseline_<TS> \
        --run qtop_shifted=scan_outputs/rs_anarchy_run3_qtop_shifted_<TS> \
        --run moreUV=scan_outputs/rs_anarchy_run3_moreUV_<TS> \
        --run moreIR=scan_outputs/rs_anarchy_run3_moreIR_<TS> \
        --out-dir results/figures/quark
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
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from quarkConstraints.finite_stats import wilson_upper_limit

DEFAULT_OUT = REPO / "results/figures/quark"
MKK_TICKS = [1, 2, 3, 5, 10, 20, 30, 50, 100, 200]
PALETTE = {
    "baseline":       "#1f77b4",
    "qtop_shifted":   "#d62728",
    "moreUV":         "#2ca02c",
    "moreIR":         "#ff7f0e",
}


def _format_count_short(n: int) -> str:
    if n >= 1_000_000:
        return f"{n / 1_000_000:.1f}M"
    if n >= 1_000:
        return f"{n / 1_000:.0f}k"
    return str(n)


def _format_sci_tex(x: float) -> str:
    mantissa, _, exponent = f"{x:.1e}".partition("e")
    return rf"{mantissa}\times10^{{{int(exponent)}}}"


def _zero_pass_label(run_dir: Path) -> str:
    summary_path = run_dir / "tile_summary.json"
    if not summary_path.exists():
        return "0 PDG-passes"
    data = json.loads(summary_path.read_text())
    n_total = sum(int(tile["n_draws"]) for tile in data["tiles"])
    n_pass = sum(int(tile["n_pdg_pass"]) for tile in data["tiles"])
    if n_pass != 0:
        return f"{n_pass:,} PDG-passes"
    p_ul = wilson_upper_limit(0, n_total)
    return (
        rf"N={_format_count_short(n_total)}, "
        rf"$p\leq{_format_sci_tex(p_ul)}$ $z=1.92$ UL"
    )


def _summary_has_zero_pdg_passes(run_dir: Path) -> bool:
    summary_path = run_dir / "tile_summary.json"
    if not summary_path.exists():
        return False
    data = json.loads(summary_path.read_text())
    return sum(int(tile["n_pdg_pass"]) for tile in data["tiles"]) == 0


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
    p.add_argument(
        "--run", action="append", default=[],
        help="<tag>=<run_dir> mapping. Repeatable. The run_dir must contain a "
             "draws.jsonl file.",
    )
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
            sys.stderr.write(f"ERROR: bad --run spec (expected TAG=DIR): {it}\n")
            sys.exit(1)
        tag, _, path = it.partition("=")
        rd = Path(path)
        draws = rd / "draws.jsonl"
        if not draws.exists():
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
        if _summary_has_zero_pdg_passes(run_dir):
            print(f"loading {tag} -> {run_dir / 'tile_summary.json'}", flush=True)
            print(f"  {tag}: 0 PDG-passing draws")
            per_run[tag] = np.asarray([])
            continue
        draws = run_dir / "draws.jsonl"
        print(f"loading {tag} -> {draws}", flush=True)
        arr = _load_pdg_passing(draws)
        print(f"  {tag}: {arr.size:,} PDG-passing draws")
        per_run[tag] = arr

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.7))
    bins = np.logspace(np.log10(0.3), np.log10(200.0), 70)

    # Panel 1: density
    ax = axes[0]
    empty_tags = [tag for tag, arr in per_run.items() if arr.size == 0]
    for tag, arr in per_run.items():
        if arr.size == 0:
            continue
        color = PALETTE.get(tag, None)
        ax.hist(arr, bins=bins, density=True, histtype="step", linewidth=2.0,
                color=color, label=f"{tag} ({arr.size:,})")
    # Add a "no PDG-passes" entry per tag with empty arr so the legend reflects this
    for tag in empty_tags:
        color = PALETTE.get(tag, None)
        ax.plot([], [], color=color, linewidth=2.0, linestyle=":",
                label=f"{tag} ({_zero_pass_label(runs[tag])})")
    ax.set_xscale("log")
    ax.xaxis.set_major_locator(FixedLocator(MKK_TICKS))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:g}"))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel(r"$M_{\mathrm{KK}}^{\min}$ per draw  [TeV]")
    ax.set_ylabel("probability density")
    ax.set_title("Per-draw $M_{\\mathrm{KK}}^{\\min}$ vs c-value pattern")
    ax.legend(fontsize=10, loc="upper right", title="c-pattern")
    ax.grid(True, which="both", alpha=0.25, linestyle=":")

    # Panel 2: cumulative
    ax = axes[1]
    for tag, arr in per_run.items():
        if arr.size == 0:
            continue
        s = np.sort(arr)
        cdf = (np.arange(1, s.size + 1) / s.size) * 100.0
        color = PALETTE.get(tag, None)
        ax.plot(s, cdf, linewidth=2.0, color=color,
                label=f"{tag} ({arr.size:,})")
    for tag in empty_tags:
        color = PALETTE.get(tag, None)
        ax.plot([], [], color=color, linewidth=2.0, linestyle=":",
                label=f"{tag} ({_zero_pass_label(runs[tag])})")
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
    ax.legend(fontsize=10, loc="lower right", title="c-pattern")
    ax.grid(True, which="both", alpha=0.25, linestyle=":")

    suptitle = "RS-anarchy ensemble: $M_{\\mathrm{KK}}^{\\min}$ vs c-value pattern"
    if args.label_tag:
        suptitle = f"{suptitle}  [{args.label_tag}]"
    fig.suptitle(suptitle, fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])

    stem = "rs_anarchy_mkk_min_hist_by_cvals"
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
        print(f"  {tag:14s}: median {np.percentile(arr, 50):6.2f}  "
              f"p95 {np.percentile(arr, 95):6.2f}   (n={arr.size:,})")


if __name__ == "__main__":
    main()
