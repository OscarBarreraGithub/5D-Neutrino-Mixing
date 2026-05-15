"""Run 1 plotter: M_KK^min histogram split by PDG-match tightness.

Re-tabulates an existing draws.jsonl into three sub-ensembles by tightening
the PDG-match factor (the original scan ran with factor-of-3 mass + CKM
tolerance and factor-of-5 J tolerance; this script re-derives the
factor-1.5 and factor-2 sub-ensembles using the per-draw log-residual
fields already stored in draws.jsonl).

A draw passes "factor X" iff
    max(up_log_max, down_log_max, ckm_log_max) < ln(X)   (mass + CKM)
    j_log < ln(5*X/3)                                   (J tracks at 5/3*X)

The per-draw M_KK^min is M_KK_tile * sqrt(max_ratio).

Default outputs:
  results/figures/quark/rs_anarchy_mkk_min_hist_by_pdg_tightness.{pdf,png}
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DRAWS = REPO / "scan_outputs/rs_anarchy_20260507T030811/draws.jsonl"
DEFAULT_OUT = REPO / "results/figures/quark"

TIGHTNESS_FACTORS = [1.5, 2.0, 3.0]
TIGHTNESS_COLORS = {1.5: "C0", 2.0: "C2", 3.0: "C3"}
MKK_TICKS = [1, 2, 3, 5, 10, 20, 30, 50, 100, 200]


def _passes_tightness(row: dict, factor: float) -> bool:
    if not row.get("ok"):
        return False
    up_lm = row.get("up_log_max")
    dn_lm = row.get("down_log_max")
    ckm_lm = row.get("ckm_log_max")
    j_l = row.get("j_log")
    if up_lm is None or dn_lm is None or ckm_lm is None or j_l is None:
        return False
    log_mass = math.log(factor)
    log_j = math.log(5.0 * factor / 3.0)
    if max(float(up_lm), float(dn_lm), float(ckm_lm)) >= log_mass:
        return False
    if float(j_l) >= log_j:
        return False
    return True


def _load_per_factor(draws_path: Path):
    """Single streaming pass over draws.jsonl. Returns dict keyed by
    tightness factor with arrays of per-draw M_KK^min (TeV)."""
    out = {f: [] for f in TIGHTNESS_FACTORS}
    n_total = 0
    with draws_path.open() as fh:
        for line in fh:
            try:
                row = json.loads(line)
            except json.JSONDecodeError:
                continue
            n_total += 1
            if not row.get("ok"):
                continue
            mr = row.get("max_ratio")
            if mr is None or not np.isfinite(mr) or mr <= 0:
                continue
            mkk_t_tev = float(row["M_KK_GeV"]) / 1000.0
            mkk_min_tev = mkk_t_tev * math.sqrt(float(mr))
            for f in TIGHTNESS_FACTORS:
                if _passes_tightness(row, f):
                    out[f].append(mkk_min_tev)
    return n_total, {f: np.asarray(v) for f, v in out.items()}


def _parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--draws", type=str, default=str(DEFAULT_DRAWS))
    p.add_argument("--out-dir", type=str, default=str(DEFAULT_OUT))
    p.add_argument("--label-tag", type=str, default="")
    return p.parse_args()


def main():
    args = _parse_args()
    draws_path = Path(args.draws)
    if not draws_path.exists():
        sys.stderr.write(f"ERROR: draws file not found: {draws_path}\n")
        sys.exit(1)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"loading {draws_path} ...", flush=True)
    n_total, per_factor = _load_per_factor(draws_path)
    print(f"  {n_total:,} rows total")
    for f in TIGHTNESS_FACTORS:
        print(f"  factor-{f}: {per_factor[f].size:,} draws pass")

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.7))

    # ---- Panel 1: differential probability density ----
    ax = axes[0]
    bins = np.logspace(np.log10(0.3), np.log10(200.0), 70)
    for f in TIGHTNESS_FACTORS:
        arr = per_factor[f]
        if arr.size == 0:
            continue
        ax.hist(
            arr, bins=bins, density=True, histtype="step",
            linewidth=2.0, color=TIGHTNESS_COLORS[f],
            label=f"factor {f:g} ({arr.size:,})",
        )
    ax.set_xscale("log")
    ax.xaxis.set_major_locator(FixedLocator(MKK_TICKS))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:g}"))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel(r"$M_{\mathrm{KK}}^{\min}$ per draw  [TeV]   (perturbative $g_s\!\approx\!1.05$)")
    ax.set_ylabel("probability density")
    ax.set_title("Per-draw $M_{\\mathrm{KK}}^{\\min}$ vs PDG-match tightness")
    ax.legend(fontsize=10, loc="upper right",
              title="PDG mass/CKM tolerance")
    ax.grid(True, which="both", alpha=0.25, linestyle=":")

    # ---- Panel 2: cumulative ----
    ax = axes[1]
    for f in TIGHTNESS_FACTORS:
        arr = per_factor[f]
        if arr.size == 0:
            continue
        s = np.sort(arr)
        cdf = (np.arange(1, s.size + 1) / s.size) * 100.0
        ax.plot(s, cdf, color=TIGHTNESS_COLORS[f], linewidth=2.0,
                label=f"factor {f:g} ({arr.size:,})")
    ax.axhline(50.0, color="k", linewidth=0.7, alpha=0.5, linestyle=":")
    ax.axhline(95.0, color="k", linewidth=0.7, alpha=0.5, linestyle=":")
    ax.set_xscale("log")
    ax.xaxis.set_major_locator(FixedLocator(MKK_TICKS))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:g}"))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel(r"$M_{\mathrm{KK}}$  [TeV]")
    ax.set_ylabel(r"fraction of accepted draws with $M_{\mathrm{KK}}^{\min}\leq M_{\mathrm{KK}}$ [%]")
    ax.set_title("Cumulative distribution by tightness")
    ax.set_xlim(0.5, 250.0)
    ax.set_ylim(0, 102)
    ax.legend(fontsize=10, loc="lower right",
              title="PDG mass/CKM tolerance")
    ax.grid(True, which="both", alpha=0.25, linestyle=":")

    suptitle = ("RS-anarchy ensemble: $M_{\\mathrm{KK}}^{\\min}$ split by PDG "
                "match-tightness")
    if args.label_tag:
        suptitle = f"{suptitle}  [{args.label_tag}]"
    fig.suptitle(suptitle, fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])

    stem = "rs_anarchy_mkk_min_hist_by_pdg_tightness"
    if args.label_tag:
        stem = f"{stem}_{args.label_tag}"
    for ext in ("pdf", "png"):
        path = out_dir / f"{stem}.{ext}"
        fig.savefig(path, dpi=200)
        print(f"wrote {path}")

    print()
    print("=== Median (TeV) per tightness factor ===")
    for f in TIGHTNESS_FACTORS:
        arr = per_factor[f]
        if arr.size == 0:
            continue
        print(f"  factor {f:g}: median {np.percentile(arr, 50):6.2f}, "
              f"p95 {np.percentile(arr, 95):6.2f}    (n={arr.size:,})")


if __name__ == "__main__":
    main()
