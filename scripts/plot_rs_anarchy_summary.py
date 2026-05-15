"""RS-anarchy ensemble summary plots: PDG-pass fraction, max-ratio
percentiles, and per-system pass rate.

Default reads:
  scan_outputs/rs_anarchy_20260507T030811/tile_summary.json
Default writes (one file per figure, .pdf and .png each):
  results/figures/quark/rs_anarchy_pdg_pass_fraction.{pdf,png}
  results/figures/quark/rs_anarchy_max_ratio_vs_mkk.{pdf,png}
  results/figures/quark/rs_anarchy_per_system_pass_rate.{pdf,png}
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
DEFAULT_RUN = REPO / "scan_outputs" / "rs_anarchy_20260507T030811"
DEFAULT_FIG = REPO / "results" / "figures" / "quark"

MKK_TICKS_TEV = [3, 5, 7, 10, 15, 20, 30, 50]
RATIO_TICKS = [0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0]


def plain_log_axis(ax, axis: str, ticks):
    """Set explicit tick positions and plain-number labels on a log-scale axis."""
    a = ax.xaxis if axis == "x" else ax.yaxis
    a.set_major_locator(FixedLocator(ticks))
    a.set_major_formatter(FuncFormatter(lambda x, _: f"{x:g}"))
    a.set_minor_formatter(NullFormatter())


def _parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    # Accept either a directory (--summary <dir>/tile_summary.json or just <dir>)
    # or a draws JSONL plus its summary explicitly.
    p.add_argument("--summary", type=str,
                   default=str(DEFAULT_RUN / "tile_summary.json"),
                   help="Path to tile_summary.json")
    p.add_argument("--draws", type=str,
                   default=str(DEFAULT_RUN / "draws.jsonl"),
                   help="(unused here, accepted for orchestrator parity).")
    p.add_argument("--out-dir", type=str, default=str(DEFAULT_FIG),
                   help="Directory to write figures into.")
    p.add_argument("--label-tag", type=str, default="",
                   help="If set, append '_<label-tag>' to each figure stem.")
    return p.parse_args()


def main():
    args = _parse_args()
    summary_path = Path(args.summary)
    if not summary_path.exists():
        sys.stderr.write(f"ERROR: summary file not found: {summary_path}\n")
        sys.exit(1)
    fig_dir = Path(args.out_dir)
    fig_dir.mkdir(parents=True, exist_ok=True)
    label_tag = args.label_tag

    def stem(name: str) -> str:
        return f"{name}_{label_tag}" if label_tag else name

    with summary_path.open() as f:
        summary = json.load(f)
    tiles = sorted(summary["tiles"], key=lambda t: t["M_KK_GeV"])
    mkks = np.array([t["M_KK_GeV"] / 1000 for t in tiles])

    # -- Figure 1: PDG-pass fraction vs M_KK --
    fig, ax = plt.subplots(figsize=(8, 5))
    fracs = np.array([t["pdg_pass_fraction"] for t in tiles])
    ax.plot(mkks, fracs * 100, "o-", linewidth=2.4, markersize=9, color="C0")
    ax.set_xscale("log")
    plain_log_axis(ax, "x", MKK_TICKS_TEV)
    ax.set_xlabel(r"$M_{\mathrm{KK}}$ [TeV]", fontsize=12)
    ax.set_ylabel("PDG-pass fraction (anarchic Y) [%]", fontsize=12)
    title = r"RS flavor anarchy --- PDG-pass rate vs $M_{\mathrm{KK}}$"
    if label_tag:
        title = f"{title} [{label_tag}]"
    ax.set_title(title, fontsize=13)
    ax.grid(True, alpha=0.3, which="both")
    ax.set_ylim(0, max(fracs) * 110)
    fig.tight_layout()
    fig.savefig(fig_dir / f"{stem('rs_anarchy_pdg_pass_fraction')}.png",
                dpi=200, bbox_inches="tight")
    fig.savefig(fig_dir / f"{stem('rs_anarchy_pdg_pass_fraction')}.pdf",
                bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {stem('rs_anarchy_pdg_pass_fraction')}.{{pdf,png}}")

    # -- Figure 2: max-ratio percentiles among PDG-passing draws --
    p05 = np.array([t["max_ratio_percentiles_pdg"]["p05"] for t in tiles])
    p25 = np.array([t["max_ratio_percentiles_pdg"]["p25"] for t in tiles])
    p50 = np.array([t["max_ratio_percentiles_pdg"]["p50"] for t in tiles])
    p75 = np.array([t["max_ratio_percentiles_pdg"]["p75"] for t in tiles])
    p95 = np.array([t["max_ratio_percentiles_pdg"]["p95"] for t in tiles])

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(mkks, p50, "o-", color="C3", linewidth=2.5, label="median", markersize=9)
    ax.fill_between(mkks, p25, p75, color="C3", alpha=0.30, label="25%--75%")
    ax.fill_between(mkks, p05, p95, color="C3", alpha=0.12, label="5%--95%")
    ax.axhline(1.0, color="k", linestyle="--", alpha=0.6, label=r"$\Delta F=2$ bound")
    ax.set_xscale("log")
    ax.set_yscale("log")
    plain_log_axis(ax, "x", MKK_TICKS_TEV)
    plain_log_axis(ax, "y", RATIO_TICKS)
    ax.set_xlabel(r"$M_{\mathrm{KK}}$ [TeV]", fontsize=12)
    ax.set_ylabel(r"max $\Delta F=2$ ratio (PDG-passing draws)", fontsize=12)
    title2 = r"RS flavor anarchy --- typical FCNC ratio vs $M_{\mathrm{KK}}$"
    if label_tag:
        title2 = f"{title2} [{label_tag}]"
    ax.set_title(title2, fontsize=13)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(fig_dir / f"{stem('rs_anarchy_max_ratio_vs_mkk')}.png",
                dpi=200, bbox_inches="tight")
    fig.savefig(fig_dir / f"{stem('rs_anarchy_max_ratio_vs_mkk')}.pdf",
                bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {stem('rs_anarchy_max_ratio_vs_mkk')}.{{pdf,png}}")

    # -- Figure 3: per-system pass rate vs M_KK --
    fig, ax = plt.subplots(figsize=(8, 5))
    systems = ["epsilon_K", "Delta_m_K", "Delta_m_Bd", "Delta_m_Bs", "Delta_m_D0"]
    labels = [r"$\varepsilon_K$", r"$\Delta M_K$", r"$\Delta M_{B_d}$",
              r"$\Delta M_{B_s}$", r"$\Delta M_{D^0}$"]
    colors_sys = ["#7B2D8E", "#2166AC", "#D95F02", "#C0392B", "#1B7837"]
    for s, l, c in zip(systems, labels, colors_sys):
        fr = np.array([t["per_system_pass_fraction_pdg"][s] for t in tiles]) * 100
        ax.plot(mkks, fr, "o-", color=c, label=l, linewidth=2, markersize=8)
    ax.axhline(50, color="k", linestyle=":", alpha=0.5)
    ax.set_xscale("log")
    plain_log_axis(ax, "x", MKK_TICKS_TEV)
    ax.set_xlabel(r"$M_{\mathrm{KK}}$ [TeV]", fontsize=12)
    ax.set_ylabel("Per-system pass rate among PDG-pass draws [%]", fontsize=12)
    title3 = r"RS flavor anarchy --- which $\Delta F=2$ system binds, vs $M_{\mathrm{KK}}$"
    if label_tag:
        title3 = f"{title3} [{label_tag}]"
    ax.set_title(title3, fontsize=12)
    ax.legend(fontsize=10, loc="lower right")
    ax.grid(True, alpha=0.3, which="both")
    ax.set_ylim(0, 105)
    fig.tight_layout()
    fig.savefig(fig_dir / f"{stem('rs_anarchy_per_system_pass_rate')}.png",
                dpi=200, bbox_inches="tight")
    fig.savefig(fig_dir / f"{stem('rs_anarchy_per_system_pass_rate')}.pdf",
                bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {stem('rs_anarchy_per_system_pass_rate')}.{{pdf,png}}")


if __name__ == "__main__":
    main()
