"""Re-tabulate the RS-anarchy ensemble under varying pass-gates and
g_s conventions, to make the M_KK bound a band rather than a single number.

Default outputs:
  results/figures/quark/rs_anarchy_gate_sensitivity.{pdf,png}
  scan_outputs/rs_anarchy_20260507T030811/gate_sensitivity.json
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
DEFAULT_DRAWS = REPO / "scan_outputs/rs_anarchy_20260507T030811/draws.jsonl"
DEFAULT_SUMMARY = REPO / "scan_outputs/rs_anarchy_20260507T030811/tile_summary.json"
DEFAULT_OUT_DIR = REPO / "results/figures/quark"


# Perturbative g_s used by the code (alpha_s(M_KK) -> g_s ~ 1.05 at M_KK ~ 3 TeV).
# Wilson coefficients scale as g_s^2; ratio_to_bound therefore scales as g_s^2.
G_S_PERT = 1.05
G_S_STAR = 3.0
RATIO_RESCALE_GS_STAR = (G_S_STAR / G_S_PERT) ** 2  # ~ 8.16


def load_per_tile(draws_path: Path):
    """Stream draws.jsonl, group max_ratio by (M_KK, passes_pdg)."""
    by_tile_pdg = {}
    by_tile_all = {}
    with draws_path.open() as fh:
        for line in fh:
            try:
                row = json.loads(line)
            except json.JSONDecodeError:
                continue
            if not row.get("ok"):
                continue
            mkk = float(row["M_KK_GeV"])
            mr = float(row["max_ratio"])
            by_tile_all.setdefault(mkk, []).append(mr)
            if row.get("passes_pdg"):
                by_tile_pdg.setdefault(mkk, []).append(mr)
    for d in (by_tile_pdg, by_tile_all):
        for k in d:
            d[k] = np.asarray(d[k], dtype=float)
    return by_tile_pdg, by_tile_all


def acceptance_fraction(ratios: np.ndarray, threshold: float) -> float:
    if ratios.size == 0:
        return float("nan")
    return float(np.mean(ratios < threshold))


def _parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--draws", type=str, default=str(DEFAULT_DRAWS))
    p.add_argument("--summary", type=str, default=str(DEFAULT_SUMMARY),
                   help="(unused; accepted for orchestrator parity).")
    p.add_argument("--out-dir", type=str, default=str(DEFAULT_OUT_DIR))
    p.add_argument("--label-tag", type=str, default="",
                   help="If set, append '_<label-tag>' to file stems.")
    return p.parse_args()


def main():
    args = _parse_args()
    draws_path = Path(args.draws)
    if not draws_path.exists():
        sys.stderr.write(f"ERROR: draws file not found: {draws_path}\n")
        sys.exit(1)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_json = draws_path.parent / (
        f"gate_sensitivity_{args.label_tag}.json"
        if args.label_tag else "gate_sensitivity.json"
    )

    print("loading draws ...", flush=True)
    by_tile_pdg, _ = load_per_tile(draws_path)
    if not by_tile_pdg:
        sys.stderr.write("ERROR: no PDG-passing draws found.\n")
        sys.exit(1)
    mkks = np.array(sorted(by_tile_pdg.keys()))
    print(f"  {len(mkks)} tiles, {sum(r.size for r in by_tile_pdg.values()):,} PDG-passing draws")

    plot_gates = [
        ("ratio < 1.0  (saturate NP budget)", 1.0),
        ("ratio < 0.5  (half NP budget)", 0.5),
        ("ratio < 0.3  (Buras / UTfit style)", 0.3),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.5), sharey=True)

    for ax, (g_s_label, rescale) in zip(
        axes, [("perturbative $g_s\\!\\approx\\!1.05$", 1.0),
               (f"$g_s^\\star = {G_S_STAR:.0f}$ (rescale $\\times {RATIO_RESCALE_GS_STAR:.1f}$)",
                RATIO_RESCALE_GS_STAR)]):
        for label, thresh in plot_gates:
            fracs = np.array([
                acceptance_fraction(by_tile_pdg[m] * rescale, thresh)
                for m in mkks
            ])
            ax.plot(mkks / 1000.0, fracs, marker="o", label=label, linewidth=1.8)
        ax.set_xscale("log")
        mkk_ticks = [3, 5, 7, 10, 15, 20, 30, 50, 100, 200]
        ax.xaxis.set_major_locator(FixedLocator(mkk_ticks))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:g}"))
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.set_xlabel(r"$M_{\mathrm{KK}}$  [TeV]")
        ax.set_title(g_s_label)
        ax.set_ylim(0.0, 1.02)
        ax.grid(True, which="both", alpha=0.25, linestyle=":")
        ax.axhline(0.5, color="k", linewidth=0.5, alpha=0.4)
        ax.axhline(0.95, color="k", linewidth=0.5, alpha=0.4)
        ax.legend(fontsize=8.5, loc="lower right")

    axes[0].set_ylabel("fraction of PDG-passing draws with all FCNC ratios below gate")
    suptitle = "Gate sensitivity of the RS-anarchy $M_{\\mathrm{KK}}$ bound"
    if args.label_tag:
        suptitle = f"{suptitle}  [{args.label_tag}]"
    fig.suptitle(suptitle, fontsize=12)
    fig.tight_layout()

    stem = "rs_anarchy_gate_sensitivity"
    if args.label_tag:
        stem = f"{stem}_{args.label_tag}"
    for ext in ("pdf", "png"):
        fig.savefig(out_dir / f"{stem}.{ext}", dpi=200)
        print(f"wrote {out_dir}/{stem}.{ext}")

    out = {"g_s_pert": G_S_PERT, "g_s_star": G_S_STAR, "rescale": RATIO_RESCALE_GS_STAR,
           "tiles_TeV": (mkks / 1000.0).tolist(), "crossings": {}}

    def crossing(mkks_tev: np.ndarray, fracs: np.ndarray, target: float) -> float:
        for i in range(len(fracs) - 1):
            f0, f1 = fracs[i], fracs[i + 1]
            if (f0 - target) * (f1 - target) <= 0 and f0 != f1:
                t = (target - f0) / (f1 - f0)
                lm = np.log10(mkks_tev[i]) + t * (np.log10(mkks_tev[i + 1]) - np.log10(mkks_tev[i]))
                return float(10.0 ** lm)
        return float("nan")

    for g_s_label, rescale in [("g_s_pert", 1.0), ("g_s_star_3", RATIO_RESCALE_GS_STAR)]:
        out["crossings"][g_s_label] = {}
        for label, thresh in plot_gates:
            fracs = np.array([
                acceptance_fraction(by_tile_pdg[m] * rescale, thresh)
                for m in mkks
            ])
            out["crossings"][g_s_label][label] = {
                "M_KK_at_50pct_TeV": crossing(mkks / 1000.0, fracs, 0.5),
                "M_KK_at_95pct_TeV": crossing(mkks / 1000.0, fracs, 0.95),
            }

    out_json.write_text(json.dumps(out, indent=2))
    print(f"wrote {out_json}")
    print()
    print("=== Crossings (M_KK [TeV] at given acceptance threshold) ===")
    for g_s_label, by_gate in out["crossings"].items():
        print(f"\n{g_s_label}:")
        for gate_label, vals in by_gate.items():
            print(f"  {gate_label:40s}  50% -> {vals['M_KK_at_50pct_TeV']:6.2f}    95% -> {vals['M_KK_at_95pct_TeV']:6.2f}")


if __name__ == "__main__":
    main()
