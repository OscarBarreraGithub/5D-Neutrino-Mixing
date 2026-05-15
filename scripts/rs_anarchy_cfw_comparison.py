"""Run C plotter: CFW-style comparison of the M_KK^min CDF.

Reads a run directory (with the tightened PDG gate + tighter Y-floor),
computes M_KK^min at the perturbative coupling, and overlays the
g_s* = 4 rescaling (M_KK -> M_KK * sqrt((4/1.05)^2) ~ x3.81). Marks
the CFW abstract benchmarks at 21 TeV ("generic anarchic") and 33 TeV
("PGB Higgs") as vertical dashed reference lines.

Usage:
    python scripts/rs_anarchy_cfw_comparison.py \
        --run scan_outputs/rs_anarchy_runC_<TS> \
        --out-dir results/figures/quark
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
DEFAULT_OUT = REPO / "results/figures/quark"
MKK_TICKS = [1, 2, 3, 5, 10, 20, 30, 50, 100, 200]

G_S_PERT = 1.05
G_S_STAR_CFW = 4.0
GS_RESCALE_MKK = G_S_STAR_CFW / G_S_PERT  # ~ 3.81

CFW_GENERIC_TEV = 21.0
CFW_PGB_TEV = 33.0


def _load_pdg_passing(draws_path: Path, require_pdg: bool = True) -> np.ndarray:
    out = []
    with draws_path.open() as fh:
        for line in fh:
            try:
                row = json.loads(line)
            except json.JSONDecodeError:
                continue
            if not row.get("ok"):
                continue
            if require_pdg and not row.get("passes_pdg"):
                continue
            mr = row.get("max_ratio")
            if mr is None or not np.isfinite(mr) or mr <= 0:
                continue
            mkk_t_tev = float(row["M_KK_GeV"]) / 1000.0
            out.append(mkk_t_tev * math.sqrt(float(mr)))
    return np.asarray(out)


def _parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run", type=str, required=True,
                   help="Run directory containing draws.jsonl")
    p.add_argument("--out-dir", type=str, default=str(DEFAULT_OUT))
    p.add_argument("--label-tag", type=str, default="")
    return p.parse_args()


def main():
    args = _parse_args()
    run_dir = Path(args.run)
    draws = run_dir / "draws.jsonl"
    if not draws.exists():
        sys.stderr.write(f"ERROR: draws.jsonl not found in {run_dir}\n")
        sys.exit(1)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"loading {draws} ...", flush=True)
    arr_pert = _load_pdg_passing(draws, require_pdg=True)
    relaxed = False
    if arr_pert.size == 0:
        sys.stderr.write(
            "WARNING: no PDG-passing draws under tightened gate; "
            "falling back to all 'ok' draws (Delta F=2 ratios still computed).\n"
        )
        arr_pert = _load_pdg_passing(draws, require_pdg=False)
        relaxed = True
    if arr_pert.size == 0:
        sys.stderr.write("ERROR: no usable draws at all.\n")
        sys.exit(1)
    print(f"  {arr_pert.size:,} draws used (relaxed={relaxed})")
    arr_gss = arr_pert * GS_RESCALE_MKK

    fig, ax = plt.subplots(figsize=(8.5, 5.5))

    s_pert = np.sort(arr_pert)
    cdf = (np.arange(1, s_pert.size + 1) / s_pert.size) * 100.0
    ax.plot(s_pert, cdf, color="C0", linewidth=2.4,
            label=r"perturbative $g_s\!\approx\!1.05$")
    ax.plot(s_pert * GS_RESCALE_MKK, cdf, color="C3", linewidth=2.4,
            linestyle="--", label=fr"$g_s^\star = {G_S_STAR_CFW:.0f}$ (CFW-like)")

    # CFW abstract benchmarks
    ax.axvline(CFW_GENERIC_TEV, color="grey", linewidth=1.2, linestyle=":",
               label=fr"CFW generic anarchic ({CFW_GENERIC_TEV:.0f} TeV)")
    ax.axvline(CFW_PGB_TEV, color="black", linewidth=1.2, linestyle="--",
               label=fr"CFW PGB Higgs ({CFW_PGB_TEV:.0f} TeV)")

    ax.axhline(50.0, color="k", linewidth=0.7, alpha=0.4, linestyle=":")
    ax.axhline(95.0, color="k", linewidth=0.7, alpha=0.4, linestyle=":")

    ax.set_xscale("log")
    ax.xaxis.set_major_locator(FixedLocator(MKK_TICKS))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:g}"))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel(r"$M_{\mathrm{KK}}$  [TeV]")
    ax.set_ylabel(r"fraction of PDG-passing draws with $M_{\mathrm{KK}}^{\min}\leq M_{\mathrm{KK}}$ [%]")
    if relaxed:
        title = ("RS-anarchy CFW-like comparison "
                 "(tightened gate yields 0 PDG passes; using all OK draws)")
    else:
        title = ("RS-anarchy CFW-like comparison "
                 "(tightened PDG gate, Y-floor 0.5)")
    if args.label_tag:
        title = f"{title}  [{args.label_tag}]"
    ax.set_title(title, fontsize=11)
    ax.set_xlim(0.5, 250.0)
    ax.set_ylim(0, 102)
    ax.legend(fontsize=9, loc="lower right")
    ax.grid(True, which="both", alpha=0.25, linestyle=":")

    fig.tight_layout()

    stem = "rs_anarchy_cfw_comparison"
    if args.label_tag:
        stem = f"{stem}_{args.label_tag}"
    for ext in ("pdf", "png"):
        path = out_dir / f"{stem}.{ext}"
        fig.savefig(path, dpi=200)
        print(f"wrote {path}")

    print()
    print("=== M_KK^min percentiles (TeV) ===")
    for q in [5, 25, 50, 75, 90, 95, 99]:
        v_pert = np.percentile(arr_pert, q)
        v_gss = v_pert * GS_RESCALE_MKK
        print(f"  p{q:>2}: pert {v_pert:7.2f}    g_s*={G_S_STAR_CFW:.0f} {v_gss:7.2f}")


if __name__ == "__main__":
    main()
