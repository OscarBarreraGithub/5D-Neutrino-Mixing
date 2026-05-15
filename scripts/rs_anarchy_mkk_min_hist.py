"""Per-draw M_KK^min histogram for the RS-anarchy ensemble, restricted to
PDG-passing draws and split by which Delta-F=2 system binds.

For each PDG-passing draw at tile M_KK_t with per-system ratios r_X, the
M_KK at which that draw would just satisfy the FCNC bound for system X is
  M_KK_min_X = M_KK_t * sqrt(r_X)
because all five ratios scale as 1/M_KK^2 (the overlap factors only depend
on epsilon = M_KK / (xi*k) and that dependence is sub-percent across our
M_KK range). The per-draw M_KK^min is then max_X M_KK_min_X, and the
"binding system" is argmax_X r_X.

Pooling across tiles is justified because the rescaled M_KK_min collapses
all tiles onto the same distribution (modulo the negligible epsilon
dependence in the overlap factors).

Default outputs:
  results/figures/quark/rs_anarchy_mkk_min_hist.{pdf,png}
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
GS_RESCALE_MKK = G_S_STAR / G_S_PERT  # ~ 2.86; M_KK at fixed acceptance scales as g_s


def load_pdg_passing(draws_path: Path):
    """Stream draws.jsonl, return (M_KK_min_TeV, binding_system) arrays."""
    records = []
    n_total = 0
    n_pdg = 0
    with draws_path.open() as fh:
        for line in fh:
            try:
                row = json.loads(line)
            except json.JSONDecodeError:
                continue
            n_total += 1
            if not row.get("ok"):
                continue
            if not row.get("passes_pdg"):
                continue
            n_pdg += 1

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
            mkk_min_tev = mkk_t_tev * float(np.sqrt(best_r))
            records.append((mkk_min_tev, best_sys))
    print(f"  {n_total:,} draws total; {n_pdg:,} PDG-passing; {len(records):,} usable.")
    mkk_min = np.array([r[0] for r in records])
    binding = np.array([r[1] for r in records])
    return mkk_min, binding


def _parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--draws", type=str, default=str(DEFAULT_DRAWS),
                   help="Path to draws.jsonl")
    p.add_argument("--out-dir", type=str, default=str(DEFAULT_OUT),
                   help="Directory to write the figure(s) into.")
    p.add_argument("--label-tag", type=str, default="",
                   help="If set, append '_<label-tag>' to the output filename "
                        "stem to disambiguate runs.")
    return p.parse_args()


def main():
    args = _parse_args()
    draws_path = Path(args.draws)
    out_dir = Path(args.out_dir)
    if not draws_path.exists():
        sys.stderr.write(f"ERROR: draws file not found: {draws_path}\n")
        sys.exit(1)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"loading PDG-passing draws from {draws_path} ...", flush=True)
    mkk_min_pert, binding = load_pdg_passing(draws_path)
    if mkk_min_pert.size == 0:
        sys.stderr.write("ERROR: zero usable PDG-passing draws.\n")
        sys.exit(1)

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.7))

    # ---- Panel 1: differential histogram, perturbative g_s, split by binding ----
    ax = axes[0]
    bins = np.logspace(np.log10(0.3), np.log10(200.0), 70)
    ax.hist(mkk_min_pert, bins=bins, color="lightgrey", alpha=0.7,
            label=f"all PDG-passing ({len(mkk_min_pert):,})", zorder=0)
    for s in SYSTEMS:
        sel = binding == s
        if not sel.any():
            continue
        ax.hist(mkk_min_pert[sel], bins=bins, histtype="step", linewidth=1.8,
                color=SYSTEM_COLORS[s],
                label=f"binding: {SYSTEM_LABELS[s]} ({sel.sum():,})")
    ax.set_xscale("log")
    mkk_ticks = [1, 2, 3, 5, 10, 20, 30, 50, 100, 200]
    ax.xaxis.set_major_locator(FixedLocator(mkk_ticks))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:g}"))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel(r"$M_{\mathrm{KK}}^{\min}$ per draw  [TeV]   (perturbative $g_s\!\approx\!1.05$)")
    ax.set_ylabel("# PDG-passing draws")
    ax.set_title(r"Per-draw $M_{\mathrm{KK}}^{\min}$ distribution, by binding observable")
    ax.legend(fontsize=8.8, loc="upper right")
    ax.grid(True, which="both", alpha=0.25, linestyle=":")

    secax = ax.secondary_xaxis("top", functions=(lambda x: x * GS_RESCALE_MKK,
                                                 lambda x: x / GS_RESCALE_MKK))
    secax.set_xlabel(r"$M_{\mathrm{KK}}^{\min}$ at $g_s^\star = 3$  [TeV]", fontsize=10)

    # ---- Panel 2: cumulative, both conventions ----
    ax = axes[1]
    sorted_pert = np.sort(mkk_min_pert)
    cdf = (np.arange(1, len(sorted_pert) + 1) / len(sorted_pert)) * 100.0

    ax.plot(sorted_pert, cdf, color="C0", linewidth=2.0,
            label=r"perturbative $g_s\approx 1.05$")
    ax.plot(sorted_pert * GS_RESCALE_MKK, cdf, color="C3", linewidth=2.0,
            linestyle="--", label=r"strong-coupling $g_s^\star = 3$")
    ax.axhline(50.0, color="k", linewidth=0.7, alpha=0.5, linestyle=":")
    ax.axhline(95.0, color="k", linewidth=0.7, alpha=0.5, linestyle=":")
    ax.text(0.4, 51, "50%", fontsize=8.5, color="grey")
    ax.text(0.4, 96, "95%", fontsize=8.5, color="grey")

    for cross_label, frac in [("50%", 0.50), ("95%", 0.95)]:
        idx = int(frac * len(sorted_pert))
        m_pert = sorted_pert[idx]
        m_gss  = m_pert * GS_RESCALE_MKK
        ax.scatter([m_pert], [frac * 100.0], color="C0", zorder=5, s=24)
        ax.scatter([m_gss], [frac * 100.0], color="C3", zorder=5, s=24)
        ax.annotate(f"{m_pert:.1f} TeV", xy=(m_pert, frac*100), xytext=(6, -12),
                    textcoords="offset points", fontsize=8, color="C0")
        ax.annotate(f"{m_gss:.1f} TeV", xy=(m_gss, frac*100), xytext=(6, 6),
                    textcoords="offset points", fontsize=8, color="C3")

    ax.set_xscale("log")
    mkk_ticks2 = [1, 2, 3, 5, 10, 20, 30, 50, 100, 200]
    ax.xaxis.set_major_locator(FixedLocator(mkk_ticks2))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:g}"))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel(r"$M_{\mathrm{KK}}$  [TeV]")
    ax.set_ylabel(r"fraction of PDG-passing draws with $M_{\mathrm{KK}}^{\min}\leq M_{\mathrm{KK}}$ [%]")
    ax.set_title("Cumulative: fraction allowed at given $M_{\mathrm{KK}}$")
    ax.set_xlim(0.5, 250.0)
    ax.set_ylim(0, 102)
    ax.legend(fontsize=9, loc="lower right")
    ax.grid(True, which="both", alpha=0.25, linestyle=":")

    suptitle = "RS-anarchy ensemble: $M_{\\mathrm{KK}}^{\\min}$ histogram (PDG-passing draws only)"
    if args.label_tag:
        suptitle = f"{suptitle}  [{args.label_tag}]"
    fig.suptitle(suptitle, fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])

    stem = "rs_anarchy_mkk_min_hist"
    if args.label_tag:
        stem = f"{stem}_{args.label_tag}"
    for ext in ("pdf", "png"):
        path = out_dir / f"{stem}.{ext}"
        fig.savefig(path, dpi=200)
        print("wrote", path)

    print()
    print("=== M_KK_min percentiles (TeV) ===")
    for q in [5, 25, 50, 75, 90, 95, 99]:
        v = np.percentile(mkk_min_pert, q)
        print(f"  p{q:>2}: pert {v:7.2f}    g_s*=3 {v*GS_RESCALE_MKK:7.2f}")
    print()
    print("=== Binding-system fractions ===")
    for s in SYSTEMS:
        f = float((binding == s).mean())
        print(f"  {s:12s}: {f*100:5.2f}%")


if __name__ == "__main__":
    main()
