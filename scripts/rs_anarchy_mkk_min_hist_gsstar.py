"""Single-panel M_KK^min histogram for the RS-anarchy ensemble in the
g_s* = 3 (CFW / strong-coupling) convention only. PDG-passing draws,
split by binding observable.

Differs from rs_anarchy_mkk_min_hist.py in that:
  - x-axis is M_KK^min in g_s* = 3 convention (i.e. on-disk perturbative
    M_KK_min rescaled by 3/1.05 ~= 2.857)
  - no secondary axis, no cumulative panel
  - 50% / 95% percentile crossings annotated as vertical guides

Default writes:
  results/figures/quark/rs_anarchy_mkk_min_hist_gsstar.{pdf,png}
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DRAWS = REPO / "scan_outputs/rs_anarchy_runA_20260507T065040/draws.jsonl"
DEFAULT_OUT = REPO / "results/figures/quark"

SYSTEMS = ["epsilon_K", "Delta_m_K", "Delta_m_Bd", "Delta_m_Bs", "Delta_m_D0"]
SYSTEM_LABELS = {
    "epsilon_K":  r"$\varepsilon_K$",
    "Delta_m_K":  r"$\Delta m_K$",
    "Delta_m_Bd": r"$\Delta m_{B_d}$",
    "Delta_m_Bs": r"$\Delta m_{B_s}$",
    "Delta_m_D0": r"$\Delta m_{D^0}$",
}
SYSTEM_COLORS = {
    "epsilon_K":  "C3",
    "Delta_m_K":  "C0",
    "Delta_m_Bd": "C2",
    "Delta_m_Bs": "C4",
    "Delta_m_D0": "C1",
}

G_S_PERT = 1.05
G_S_STAR = 3.0
GS_RESCALE_MKK = G_S_STAR / G_S_PERT  # ~ 2.857


def load_pdg_passing(draws_path: Path):
    records = []
    with draws_path.open() as fh:
        for line in fh:
            try:
                row = json.loads(line)
            except json.JSONDecodeError:
                continue
            if not row.get("ok") or not row.get("passes_pdg"):
                continue
            mkk_t_tev = float(row["M_KK_GeV"]) / 1000.0
            ratios = row["deltaf2_ratios"]
            best_sys = None
            best_r = -np.inf
            for s in SYSTEMS:
                r = ratios.get(s)
                if r is None or not np.isfinite(r):
                    continue
                if r > best_r:
                    best_r = r
                    best_sys = s
            if best_sys is None or best_r <= 0:
                continue
            records.append((mkk_t_tev * float(np.sqrt(best_r)), best_sys))
    mkk_min = np.array([r[0] for r in records])
    binding = np.array([r[1] for r in records])
    return mkk_min, binding


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

    print(f"loading PDG-passing draws from {draws_path} ...", flush=True)
    mkk_min_pert, binding = load_pdg_passing(draws_path)
    if mkk_min_pert.size == 0:
        sys.stderr.write("ERROR: zero usable PDG-passing draws.\n")
        sys.exit(1)

    # Rescale to g_s* = 3 convention
    mkk_min = mkk_min_pert * GS_RESCALE_MKK
    print(f"  {mkk_min.size:,} PDG-passing draws "
          f"(scaled by {GS_RESCALE_MKK:.3f} to g_s* = {G_S_STAR:.0f})")

    fig, ax = plt.subplots(figsize=(8.5, 5.0))
    bins = np.logspace(np.log10(1.0), np.log10(500.0), 70)

    ax.hist(mkk_min, bins=bins, color="lightgrey", alpha=0.75,
            label=f"all PDG-passing ({mkk_min.size:,})", zorder=0)
    for s in SYSTEMS:
        sel = binding == s
        if not sel.any():
            continue
        ax.hist(mkk_min[sel], bins=bins, histtype="step", linewidth=1.8,
                color=SYSTEM_COLORS[s],
                label=f"binding: {SYSTEM_LABELS[s]} ({sel.sum():,})")

    # 50% and 95% percentile guides (computed on the rescaled distribution)
    p50 = float(np.percentile(mkk_min, 50))
    p95 = float(np.percentile(mkk_min, 95))
    for pct, val, color in [(50, p50, "k"), (95, p95, "k")]:
        ax.axvline(val, color=color, linewidth=0.9, linestyle="--", alpha=0.6)
        ax.text(val * 1.04, ax.get_ylim()[1] * 0.92 if ax.get_ylim()[1] > 0 else 1.0,
                f"{pct}% : {val:.1f} TeV", fontsize=9, color=color,
                rotation=0, va="top", ha="left")

    ax.set_xscale("log")
    mkk_ticks = [3, 5, 10, 20, 30, 50, 100, 200, 500]
    ax.xaxis.set_major_locator(FixedLocator(mkk_ticks))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:g}"))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlim(2.0, 500.0)
    ax.set_xlabel(r"$M_{\mathrm{KK}}^{\min}$  [TeV]   ($g_s^\star = 3$)",
                  fontsize=12)
    ax.set_ylabel("# PDG-passing draws", fontsize=12)
    ax.set_title(r"Per-draw $M_{\mathrm{KK}}^{\min}$ distribution"
                 r" (RS anarchy, $g_s^\star = 3$ convention)",
                 fontsize=12)
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(True, which="both", alpha=0.25, linestyle=":")
    fig.tight_layout()

    stem = "rs_anarchy_mkk_min_hist_gsstar"
    if args.label_tag:
        stem = f"{stem}_{args.label_tag}"
    for ext in ("pdf", "png"):
        path = out_dir / f"{stem}.{ext}"
        fig.savefig(path, dpi=200)
        print(f"  wrote {path}")
    print(f"  median (g_s* = 3): {p50:.2f} TeV")
    print(f"  95%    (g_s* = 3): {p95:.2f} TeV")


if __name__ == "__main__":
    main()
