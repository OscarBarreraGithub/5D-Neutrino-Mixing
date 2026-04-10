#!/usr/bin/env python3
"""Plot the quark-sector kaon proxy ratio over an ``r`` sweep."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from quarkConstraints.validation import r_sweep_plot_data

FIG_DIR = REPO_ROOT / "results" / "figures"

DEFAULT_R_VALUES = (0.05, 0.1, 0.25, 0.4, 1.0)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--r-values",
        type=float,
        nargs="+",
        default=list(DEFAULT_R_VALUES),
        help="Monotonic r values to evaluate.",
    )
    parser.add_argument(
        "--overall-scale",
        type=float,
        default=None,
        help="Optional fixed overall spurion scale passed into the fit.",
    )
    parser.add_argument(
        "--xi-kk",
        type=float,
        default=1.0,
        help="Explicit KK-scale convention factor.",
    )
    parser.add_argument(
        "--max-nfev",
        type=int,
        default=100,
        help="Maximum least-squares evaluations per sweep point.",
    )
    return parser.parse_args()


def _configure_style() -> None:
    plt.rcParams.update(
        {
            "figure.dpi": 150,
            "font.size": 11,
            "axes.titlesize": 14,
            "axes.labelsize": 12,
            "legend.fontsize": 10,
        }
    )


def _plot_epsilon_k_vs_r(data: dict[str, np.ndarray], path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    r_values = data["r_values"]
    epsilon_k_ratio = data["epsilon_k_ratio"]

    ax.plot(
        r_values,
        epsilon_k_ratio,
        marker="o",
        lw=2.5,
        ms=8,
        color="#c0392b",
        markeredgecolor="#7b241c",
        markeredgewidth=1.2,
        label=r"$\epsilon_K$ ratio",
    )
    ax.axhline(1.0, color="black", ls="--", lw=1.5, alpha=0.8, label="repo pass/fail = 1")
    ax.fill_between(r_values, epsilon_k_ratio, 1.0, where=epsilon_k_ratio <= 1.0, color="#d5f5e3", alpha=0.45)
    ax.fill_between(r_values, epsilon_k_ratio, 1.0, where=epsilon_k_ratio > 1.0, color="#fadbd8", alpha=0.45)

    for x_val, y_val in zip(r_values, epsilon_k_ratio):
        ax.annotate(f"{y_val:.2f}", (x_val, y_val), textcoords="offset points", xytext=(0, 8), ha="center")

    ax.set_title(r"$K^0-\bar{K}^0$ Proxy Constraint vs $r$")
    ax.set_xlabel(r"$r$")
    ax.set_ylabel(r"$\epsilon_K$ ratio to repo bound")
    ax.set_yscale("log")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="upper left", framealpha=0.95)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    args = _parse_args()
    _configure_style()
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    data = r_sweep_plot_data(
        args.r_values,
        overall_scale=args.overall_scale,
        max_nfev=args.max_nfev,
        xi_KK=args.xi_kk,
    )

    epsilon_path = FIG_DIR / "quark_epsilon_k_vs_r.png"
    _plot_epsilon_k_vs_r(data, epsilon_path)

    print("Saved:", epsilon_path.relative_to(REPO_ROOT))
    print("r_values =", data["r_values"])
    print("epsilon_k_ratio =", data["epsilon_k_ratio"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
