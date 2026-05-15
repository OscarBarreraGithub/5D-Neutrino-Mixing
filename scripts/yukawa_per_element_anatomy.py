"""Anatomy of the envelope-fit Yukawa entries: are the off-diagonals all
suppressed by the same factor (a single 'r'), or are they individually
tuned?

Outputs:
  results/figures/quark/yukawa_per_element_anatomy.{pdf,png}
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter

REPO = Path(__file__).resolve().parents[1]
ENV_CSV = REPO / "scan_outputs/dense_20260506T141321/derived/accepted_points_with_yukawas.csv"
OUT = REPO / "results/figures/quark"


def load_yukawa_magnitudes(df: pd.DataFrame) -> dict:
    """Return |Y_f_ij| arrays for f in 'ud' and (i,j) in 1..3."""
    out = {}
    for f in ("u", "d"):
        for i in (1, 2, 3):
            for j in (1, 2, 3):
                mag = np.hypot(df[f"Y_{f}_{i}{j}_re"].values,
                               df[f"Y_{f}_{i}{j}_im"].values)
                out[(f, i, j)] = mag
    return out


def main():
    df = pd.read_csv(ENV_CSV)
    print(f"loaded {len(df):,} accepted points")
    Y = load_yukawa_magnitudes(df)

    # ============================================================
    # Panel A: per-element median |Y_ij| as 3x3 heatmaps for u, d
    # ============================================================
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 9.0),
                              gridspec_kw={"width_ratios": [1.0, 1.4]})

    for row, f in enumerate("ud"):
        med = np.zeros((3, 3))
        for i in (1, 2, 3):
            for j in (1, 2, 3):
                med[i - 1, j - 1] = np.median(Y[(f, i, j)])
        # log10 for readability
        ax = axes[row, 0]
        im = ax.imshow(np.log10(med), cmap="RdBu_r", vmin=-3, vmax=0.7,
                       aspect="equal")
        for i in range(3):
            for j in range(3):
                ax.text(j, i, f"{med[i, j]:.3g}", ha="center", va="center",
                        fontsize=11,
                        color="white" if abs(np.log10(med[i, j])) > 1.0 else "black")
        ax.set_xticks([0, 1, 2]); ax.set_yticks([0, 1, 2])
        ax.set_xticklabels(["1", "2", "3"]); ax.set_yticklabels(["1", "2", "3"])
        ax.set_xlabel("j (right)"); ax.set_ylabel("i (left)")
        ax.set_title(fr"Median $|Y_{{{f},ij}}|$  (envelope fit)")
        plt.colorbar(im, ax=ax, label=r"$\log_{10}\,|Y_{ij}|$",
                     fraction=0.046, pad=0.04)

    # ============================================================
    # Panel B: distribution of each off-diagonal element's |Y|, stacked
    # ============================================================
    OFFDIAG = [(i, j) for i in (1, 2, 3) for j in (1, 2, 3) if i != j]

    bins = np.logspace(-4.5, 1.0, 80)
    for row, f in enumerate("ud"):
        ax = axes[row, 1]
        for (i, j), color in zip(OFFDIAG, ["C0","C1","C2","C3","C4","C5"]):
            arr = Y[(f, i, j)]
            ax.hist(arr, bins=bins, density=True, histtype="step",
                    linewidth=1.7, color=color,
                    label=fr"$|Y_{{{f},{i}{j}}}|$  (med {np.median(arr):.3g})")
        ax.set_xscale("log")
        ax.xaxis.set_major_locator(FixedLocator([1e-4, 1e-3, 1e-2, 0.1, 1.0, 10.0]))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:g}"))
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.set_xlim(1e-4, 6.0)
        ax.set_xlabel(r"$|Y_{f,ij}|$")
        ax.set_ylabel("probability density")
        ax.set_title(fr"Per-element off-diagonal distribution: $Y_{f}$")
        ax.legend(fontsize=8.5, loc="upper left", ncol=2)
        ax.grid(True, which="both", alpha=0.25, linestyle=":")

    fig.suptitle("Envelope-fit Yukawa anatomy: is there a single off-diagonal "
                 r"suppression factor $r$?", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.97])

    out_png = OUT / "yukawa_per_element_anatomy.png"
    out_pdf = OUT / "yukawa_per_element_anatomy.pdf"
    fig.savefig(out_png, dpi=180)
    fig.savefig(out_pdf)
    print(f"wrote {out_png}")
    print(f"wrote {out_pdf}")

    # ============================================================
    # Diagnostic: per-point spread of log10|Y_off|.
    # If a single 'r' suppresses all off-diagonals, std should be tiny.
    # If each off-diagonal is independently tuned, std should be O(1).
    # ============================================================
    n = len(df)
    log_off_u = np.zeros((n, 6))
    log_off_d = np.zeros((n, 6))
    for k, (i, j) in enumerate(OFFDIAG):
        log_off_u[:, k] = np.log10(Y[("u", i, j)])
        log_off_d[:, k] = np.log10(Y[("d", i, j)])
    std_per_point_u = np.std(log_off_u, axis=1)
    std_per_point_d = np.std(log_off_d, axis=1)

    fig2, ax2 = plt.subplots(figsize=(8.0, 4.5))
    ax2.hist(std_per_point_u, bins=60, density=True, histtype="step",
             linewidth=2.0, color="C3",
             label=fr"$Y_u$ off-diagonals (median spread "
                   fr"{np.median(std_per_point_u):.2f} dex)")
    ax2.hist(std_per_point_d, bins=60, density=True, histtype="step",
             linewidth=2.0, color="C0",
             label=fr"$Y_d$ off-diagonals (median spread "
                   fr"{np.median(std_per_point_d):.2f} dex)")
    ax2.axvline(0.1, color="grey", linestyle=":", alpha=0.7)
    ax2.text(0.11, ax2.get_ylim()[1] * 0.9 if ax2.get_ylim()[1] > 0 else 1,
             "0.1 dex (single-$r$ would be here)",
             fontsize=9, color="grey")
    ax2.set_xlabel(r"per-point std$_{i\neq j}[\log_{10} |Y_{ij}|]$  [dex]",
                   fontsize=11)
    ax2.set_ylabel("probability density")
    ax2.set_title("Spread of off-diagonal Yukawa magnitudes within a single fit point",
                  fontsize=12)
    ax2.legend(fontsize=10)
    ax2.grid(True, which="both", alpha=0.25, linestyle=":")
    fig2.tight_layout()
    out2_png = OUT / "yukawa_offdiag_spread.png"
    out2_pdf = OUT / "yukawa_offdiag_spread.pdf"
    fig2.savefig(out2_png, dpi=180)
    fig2.savefig(out2_pdf)
    print(f"wrote {out2_png}")
    print(f"wrote {out2_pdf}")

    # Print numerical summary
    print()
    print("=== Per-element median |Y_ij| ===")
    for f in "ud":
        print(f"\n  Y_{f}:")
        for i in (1, 2, 3):
            row = "  ".join(f"{np.median(Y[(f, i, j)]):8.4f}" for j in (1, 2, 3))
            print(f"   i={i}:  {row}")
    print()
    print(f"=== Per-point spread of log10|Y_off| (single-r would be ~0) ===")
    print(f"  Y_u: median {np.median(std_per_point_u):.3f} dex, "
          f"95% {np.percentile(std_per_point_u, 95):.3f} dex")
    print(f"  Y_d: median {np.median(std_per_point_d):.3f} dex, "
          f"95% {np.percentile(std_per_point_d, 95):.3f} dex")
    print()
    print("Interpretation: if median spread << 0.3 dex, the off-diagonals are")
    print("co-suppressed (single-r). If >= 0.5 dex, they're independent.")


if __name__ == "__main__":
    main()
