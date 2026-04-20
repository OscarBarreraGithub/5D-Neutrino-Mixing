#!/usr/bin/env python3
"""Compare the repo sigmoid bulk-mass map to its affine alpha-style surrogate."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

FIG_DIR = REPO_ROOT / "results" / "figures"
DEFAULT_R_VALUES = (0.0, 0.02, 0.05, 0.1, 0.25, 0.5, 1.0, 2.0)
DEFAULT_OVERALL_SCALE_VALUES = (1.5, 2.8, 4.0, 6.0)
SECTOR_COLORS = {"Q": "#1f77b4", "u": "#d62728", "d": "#2ca02c"}


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--r-values",
        type=float,
        nargs="+",
        default=list(DEFAULT_R_VALUES),
        help="Non-negative r values used to sample the C_Q spectrum.",
    )
    parser.add_argument(
        "--overall-scale-values",
        type=float,
        nargs="+",
        default=list(DEFAULT_OVERALL_SCALE_VALUES),
        help="Positive overall spurion scales used to sample the C matrices.",
    )
    parser.add_argument(
        "--lambda-ir",
        type=float,
        default=3000.0,
        help="IR scale in GeV used in the comparison state construction.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=FIG_DIR / "quark_bulk_mass_map_comparison.png",
        help="Output image path.",
    )
    return parser.parse_args()


def _configure_style() -> None:
    plt.rcParams.update(
        {
            "figure.dpi": 150,
            "font.size": 11,
            "axes.titlesize": 13,
            "axes.labelsize": 12,
            "legend.fontsize": 10,
        }
    )


def _quantile_summary(name: str, values: np.ndarray) -> str:
    q0, q50, q90, q99, q100 = np.quantile(values, [0.0, 0.5, 0.9, 0.99, 1.0])
    return (
        f"{name}: min={q0:.4g}, p50={q50:.4g}, "
        f"p90={q90:.4g}, p99={q99:.4g}, max={q100:.4g}"
    )


def _plot_curve_panel(ax: plt.Axes, data: dict[str, np.ndarray]) -> None:
    positive = data["lambda_grid"] > 0.0
    ax.plot(
        data["lambda_grid"][positive],
        data["c_sigmoid_grid"][positive],
        lw=2.5,
        color="#1f77b4",
        label="sigmoid map",
    )
    ax.plot(
        data["lambda_grid"][positive],
        data["c_affine_grid"][positive],
        lw=2.0,
        ls="--",
        color="#d62728",
        label="affine alpha-style surrogate",
    )
    ax.axhline(float(data["c_uv"][0]), color="black", lw=1.0, alpha=0.35)
    ax.axhline(float(data["c_ir"][0]), color="black", lw=1.0, alpha=0.35)
    ax.fill_between(
        data["lambda_grid"][positive],
        float(data["c_ir"][0]),
        float(data["c_uv"][0]),
        color="#d6eaf8",
        alpha=0.35,
    )
    ax.set_xscale("log")
    ax.set_xlabel(r"spurion eigenvalue $\lambda$")
    ax.set_ylabel(r"bulk mass parameter $c$")
    ax.set_title(r"Map shape: $c(\lambda)$")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="upper right", framealpha=0.95)


def _plot_sample_panel(ax: plt.Axes, data: dict[str, np.ndarray]) -> None:
    positive = data["lambda_grid"] > 0.0
    for sector_id in ("Q", "u", "d"):
        ax.scatter(
            data[f"eig_{sector_id}_samples"],
            data[f"c_{sector_id}_sigmoid_samples"],
            s=28,
            alpha=0.75,
            color=SECTOR_COLORS[sector_id],
            label=rf"${sector_id}$ samples",
        )
    ax.plot(
        data["lambda_grid"][positive],
        data["c_sigmoid_grid"][positive],
        lw=2.2,
        color="#1f77b4",
    )
    ax.plot(
        data["lambda_grid"][positive],
        data["c_affine_grid"][positive],
        lw=1.8,
        ls="--",
        color="#d62728",
    )
    ax.axhline(float(data["c_ir"][0]), color="black", lw=1.0, alpha=0.35)
    ax.axhline(float(data["c_uv"][0]), color="black", lw=1.0, alpha=0.35)
    ax.set_xscale("log")
    ax.set_xlabel(r"sampled spurion eigenvalue $\lambda$")
    ax.set_ylabel(r"sigmoid-mapped $c$")
    ax.set_title("Observed eigenvalue cloud")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="lower left", framealpha=0.95)


def _save_plot(data: dict[str, np.ndarray], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.8))
    _plot_curve_panel(axes[0], data)
    _plot_sample_panel(axes[1], data)
    fig.tight_layout()
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    from quarkConstraints.validation import bulk_mass_map_comparison_data

    args = _parse_args()
    _configure_style()

    data = bulk_mass_map_comparison_data(
        r_values=args.r_values,
        overall_scale_values=args.overall_scale_values,
        Lambda_IR=args.lambda_ir,
    )
    _save_plot(data, args.output)

    affine_alpha = float(data["affine_alpha"][0])
    affine_beta = float(data["affine_beta"][0])
    delta = np.abs(data["c_sigmoid_samples"] - data["c_affine_samples"])

    try:
        relative_output = args.output.relative_to(REPO_ROOT)
    except ValueError:
        relative_output = args.output

    print("Saved:", relative_output)
    print(f"Affine surrogate: c = {affine_alpha:.4f} {affine_beta:+.4f} * lambda")
    print(_quantile_summary("lambda samples", data["eig_samples"]))
    print(_quantile_summary("sigmoid c samples", data["c_sigmoid_samples"]))
    print(_quantile_summary("affine c samples", data["c_affine_samples"]))
    print(_quantile_summary("window-clipped affine c", data["c_affine_clipped_samples"]))
    print(
        "Affine below c_ir fraction = "
        f"{float(data['affine_below_c_ir_fraction'][0]):.3f}"
    )
    print(f"Median |c_sigmoid - c_affine| = {float(np.median(delta)):.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
