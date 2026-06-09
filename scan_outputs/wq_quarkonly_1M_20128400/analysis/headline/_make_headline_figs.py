#!/usr/bin/env python
"""Headline figures for the WQ quark-only 1M r x M_KK scan.

Read-only on the scan cache; writes only PNGs (+README) into this directory.
Run: source ~/.bashrc && conda activate ising_bootstrap && python _make_headline_figs.py
"""
import os
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

HERE = os.path.dirname(os.path.abspath(__file__))
RUN = os.path.abspath(os.path.join(HERE, "..", ".."))
CACHE = os.path.join(RUN, "wq_quarkonly_cache.parquet")
OUT = HERE

SUBTITLE = ("quark-only, non-custodial RS  |  Z$\\to b\\bar b$ floor "
            "$\\sim$25--30 TeV is the non-custodial RS $Zb_L$ bound")

# Physics names for constraint IDs that veto anything
NAMES = {
    "T010": r"Z$\to b\bar b$  ($R_b^0$, $A_b$)",
    "T011": r"Z$\to b\bar b$ asym. ($A_{FB}^{0,b}$, $A_b$)",
    "CR001": r"KK gluon $g_{KK}\to t\bar t$ (collider)",
    "EW001": r"$S,T,U$ oblique (EWPT)",
    "CR012": r"KK $V^{(1)}\to WW/WZ/ZZ$ (collider)",
    "CR013": r"KK graviton $G^{(1)}\to\gamma\gamma$ (collider)",
    "CR002": r"KK $T_{5/3}$ pair (collider)",
    "CR003": r"KK $T$ pair (collider)",
    "CR004": r"KK $B$ pair (collider)",
    "CR008": r"KK $T$ singlet pair (collider)",
    "CR010": r"KK $(T,B)$ doublet pair (collider)",
    "B012": r"BR($B^0\to K^{*0}\gamma$) radiative",
    "B011": r"BR($\bar B\to X_s\gamma$) radiative",
    "B013": r"BR($B_s\to\phi\gamma$) radiative",
    "B003": r"$\Delta m_s$ ($B_s$ mixing)",
    "B004": r"$\phi_s$ ($B_s\to J/\psi\phi$)",
    "B001": r"$\Delta m_d$ ($B_d$ mixing)",
    "K001": r"$\epsilon_K$ (kaon CPV)",
    "K002": r"$\Delta m_K$ (kaon mixing)",
    "C001": r"$D^0$-$\bar D^0$ mixing",
}

RIG = plt.cm.tab10(3)   # red-ish for rigorous
PRX = plt.cm.tab10(0)   # blue for proxy


def titles(fig, ax, title):
    """Place a bold title + grey subtitle cleanly above the axes."""
    ax.set_title("")
    fig.suptitle(title, fontsize=17, fontweight="bold", y=0.985)
    fig.text(0.5, 0.93, SUBTITLE, ha="center", va="top",
             fontsize=11.5, color="0.35")


plt.rcParams.update({
    "font.size": 15, "axes.titlesize": 17, "axes.labelsize": 15,
    "xtick.labelsize": 13, "ytick.labelsize": 13, "legend.fontsize": 13,
    "figure.dpi": 130, "savefig.dpi": 150, "axes.grid": True,
    "grid.alpha": 0.3,
})


def load():
    pf = pq.ParquetFile(CACHE)
    names = [f.name for f in pf.schema_arrow]
    cids = sorted({n.split("__")[0] for n in names if n.endswith("__pass")})
    cols = (["M_KK", "r", "skipped", "survives_strict",
             "cQ1", "cQ2", "cQ3", "cu1", "cu2", "cu3", "cd1", "cd2", "cd3",
             "up_sv1", "up_sv2", "up_sv3", "dn_sv1", "dn_sv2", "dn_sv3"]
            + [c + "__pass" for c in cids] + [c + "__tag" for c in cids])
    df = pq.read_table(CACHE, columns=cols).to_pandas()
    ns = df[~df.skipped].copy()
    ns["TeV"] = ns.M_KK / 1000.0

    def tagof(c):
        v = ns.loc[ns[c + "__tag"] != "", c + "__tag"]
        return v.mode().iloc[0] if len(v) else ""
    tags = {c: tagof(c) for c in cids}
    return ns, cids, tags


def fig1_survival(ns):
    """Survival fraction vs M_KK (rigorous strict)."""
    grid = sorted(ns.TeV.unique())
    grid = [g for g in grid if g >= 4]
    fig, ax = plt.subplots(figsize=(10, 6.2))

    # overall (pooled over r) but plotted per-r so we never silently pool
    rs = sorted(ns.r.unique())
    cmap = plt.cm.viridis(np.linspace(0.1, 0.9, len(rs)))
    for r, col in zip(rs, cmap):
        sub = ns[ns.r == r]
        y = [sub.loc[sub.TeV == g, "survives_strict"].mean() for g in grid]
        ax.plot(grid, y, "-o", color=col, lw=2.4, ms=7,
                label=fr"$r={r:g}$")

    ax.axvspan(20, 30, color="orange", alpha=0.16, zorder=0)
    ax.annotate("rigorous $M_{KK}$ floor\n(survival 0$\\to$1)",
                xy=(25, 0.5), xytext=(33, 0.42), fontsize=13,
                ha="left", color="darkorange",
                arrowprops=dict(arrowstyle="->", color="darkorange", lw=1.8))
    ax.set_xscale("log")
    ax.set_xticks(grid)
    ax.set_xticklabels([f"{g:g}" for g in grid])
    ax.set_xlim(4, 52)
    ax.set_ylim(-0.03, 1.05)
    ax.set_xlabel(r"$M_{KK}$  [TeV]")
    ax.set_ylabel("rigorous strict survival fraction")
    ax.legend(title="fit draw $r$", loc="center left", framealpha=0.9)
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    titles(fig, ax, "Survival vs KK scale: sharp Z$\\to b\\bar b$ floor near 25--30 TeV")
    p = os.path.join(OUT, "fig1_survival_vs_MKK.png")
    fig.savefig(p); plt.close(fig)
    return p


def fig2_reach(ns, cids, tags):
    """Money plot: per-constraint reach = highest M_KK with >50% veto."""
    grid = sorted(ns.TeV.unique())
    reach = {}
    for c in cids:
        if tags[c] == "partial" or c not in NAMES:
            continue
        rr = 0.0
        for g in grid:
            sub = ns[ns.TeV == g]
            if (~sub[c + "__pass"]).mean() > 0.5:
                rr = g
        if rr > 0:
            reach[c] = (rr, tags[c])
    items = sorted(reach.items(), key=lambda x: x[1][0])  # ascending for barh
    labels = [NAMES[c] for c, _ in items]
    vals = [v[0] for _, v in items]
    cols = [RIG if v[1] == "rigorous" else PRX for _, v in items]

    fig, ax = plt.subplots(figsize=(11.5, 7.5))
    y = np.arange(len(items))
    ax.barh(y, vals, color=cols, edgecolor="black", lw=0.7, height=0.72)
    for yi, v, (c, (rr, t)) in zip(y, vals, items):
        # T010 grid reach is 20 TeV but true >50% crossing is 20-30 TeV
        txt = f"{v:g} TeV" + ("  (crossing 20--30)" if c == "T010" else "")
        ax.text(v + 0.6, yi, txt, va="center", fontsize=12.5)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlabel(r"reach: highest $M_{KK}$ that still vetoes $>$50% of points  [TeV]")
    ax.set_xlim(0, 36)
    ax.grid(axis="y", visible=False)
    leg = [Patch(facecolor=RIG, edgecolor="black", label="rigorous"),
           Patch(facecolor=PRX, edgecolor="black", label="proxy")]
    ax.legend(handles=leg, loc="lower right", framealpha=0.95)
    fig.tight_layout(rect=[0, 0, 1, 0.91])
    titles(fig, ax, "Which constraint sets the KK floor?  Z$\\to b\\bar b$ dominates by ~4$\\times$")
    p = os.path.join(OUT, "fig2_constraint_reach_ranking.png")
    fig.savefig(p); plt.close(fig)
    return p, reach


def fig3_zbb_curve(ns):
    """Z->bb veto vs M_KK with eps_K and CR001 for contrast."""
    grid = [g for g in sorted(ns.TeV.unique()) if g >= 4]
    fig, ax = plt.subplots(figsize=(10, 6.2))
    curves = [("T010__pass", NAMES["T010"], plt.cm.tab10(3), 3.2, "o"),
              ("CR001__pass", NAMES["CR001"], plt.cm.tab10(0), 2.2, "s"),
              ("K001__pass", NAMES["K001"], plt.cm.tab10(2), 2.2, "^")]
    for col, lab, c, lw, mk in curves:
        y = [(~ns.loc[ns.TeV == g, col]).mean() for g in grid]
        ax.plot(grid, y, "-" + mk, color=c, lw=lw, ms=7, label=lab)
    ax.axhline(0.5, color="0.5", ls="--", lw=1.2)
    ax.axvspan(20, 30, color="orange", alpha=0.14, zorder=0)
    ax.annotate("Z$\\to b\\bar b$ veto crosses 50%\nbetween 20 and 30 TeV",
                xy=(22, 0.5), xytext=(7, 0.30), fontsize=12.5, color="darkred",
                arrowprops=dict(arrowstyle="->", color="darkred", lw=1.7))
    ax.set_xscale("log")
    ax.set_xticks(grid); ax.set_xticklabels([f"{g:g}" for g in grid])
    ax.set_xlim(4, 52); ax.set_ylim(-0.03, 1.05)
    ax.set_xlabel(r"$M_{KK}$  [TeV]")
    ax.set_ylabel("fraction of points vetoed")
    ax.legend(loc="center left", framealpha=0.92)
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    titles(fig, ax, "Z$\\to b\\bar b$ veto dwarfs $\\epsilon_K$ and the KK-gluon collider bound")
    p = os.path.join(OUT, "fig3_Zbb_veto_vs_MKK.png")
    fig.savefig(p); plt.close(fig)
    return p


def fig4_bulkmass(ns):
    """Bulk-mass localization: c_Q, c_u, c_d per generation."""
    fig, ax = plt.subplots(figsize=(12, 6.6))
    groups = [("cQ", r"$c_Q$ (LH doublet)", plt.cm.tab10(0)),
              ("cu", r"$c_u$ (RH up)", plt.cm.tab10(1)),
              ("cd", r"$c_d$ (RH down)", plt.cm.tab10(2))]
    positions, data, fcolors = [], [], []
    xticklabels = []
    base = 0
    for prefix, lab, col in groups:
        for g in (1, 2, 3):
            data.append(ns[f"{prefix}{g}"].values)
            positions.append(base + g)
            fcolors.append(col)
            xticklabels.append(f"gen{g}")
        base += 4
    vp = ax.violinplot(data, positions=positions, widths=0.8,
                       showmeans=False, showextrema=False)
    for body, c in zip(vp["bodies"], fcolors):
        body.set_facecolor(c); body.set_alpha(0.65); body.set_edgecolor("black")
    bp = ax.boxplot(data, positions=positions, widths=0.18, showfliers=False,
                    patch_artist=True, medianprops=dict(color="black", lw=1.6))
    for box in bp["boxes"]:
        box.set_facecolor("white"); box.set_alpha=0.9
    ax.set_xlim(0, positions[-1] + 1)
    ax.axhline(0.5, color="red", ls="--", lw=1.8)
    ax.text(0.3, 0.505, "c = 1/2   (IR$\\leftrightarrow$UV crossover)",
            color="red", va="bottom", ha="left", fontsize=12.5)
    ax.set_xticks(positions)
    ax.set_xticklabels(xticklabels)
    for (prefix, lab, col), gx in zip(groups, [2, 6, 10]):
        ax.text(gx, 0.95, lab, transform=ax.get_xaxis_transform(),
                ha="center", va="top", fontsize=14, color=col, fontweight="bold")
    ax.set_ylabel("bulk mass parameter $c = M_5/k$")
    ax.grid(axis="x", visible=False)
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    titles(fig, ax, "Geometric flavor hierarchy: light quarks UV ($c>1/2$), 3rd gen IR ($c<1/2$)")
    p = os.path.join(OUT, "fig4_bulk_mass_localization.png")
    fig.savefig(p); plt.close(fig)
    return p


def fig5_yukawa(ns):
    """Fitted up/down Yukawa singular-value hierarchy."""
    fig, ax = plt.subplots(figsize=(10.5, 6.2))
    cols = [("up_sv1", "up_sv2", "up_sv3", "up-type", plt.cm.tab10(1)),
            ("dn_sv1", "dn_sv2", "dn_sv3", "down-type", plt.cm.tab10(0))]
    positions, data, fcolors, xtl = [], [], [], []
    base = 0
    for s1, s2, s3, lab, col in cols:
        for sv, gl in zip((s1, s2, s3), ("SV1", "SV2", "SV3")):
            v = ns[sv].values
            v = v[(v > 0) & np.isfinite(v)]
            data.append(np.log10(v))
            positions.append(base + len(xtl) % 3 + 1)
            fcolors.append(col); xtl.append(gl)
        base += 4
    positions = [1, 2, 3, 5, 6, 7]
    vp = ax.violinplot(data, positions=positions, widths=0.8,
                       showmeans=False, showextrema=False)
    for body, c in zip(vp["bodies"], fcolors):
        body.set_facecolor(c); body.set_alpha(0.65); body.set_edgecolor("black")
    bp = ax.boxplot(data, positions=positions, widths=0.18, showfliers=False,
                    patch_artist=True, medianprops=dict(color="black", lw=1.6))
    for box in bp["boxes"]:
        box.set_facecolor("white")
    ax.set_xticks(positions)
    ax.set_xticklabels(["SV1", "SV2", "SV3", "SV1", "SV2", "SV3"])
    ax.text(2, 0.95, "up-type", transform=ax.get_xaxis_transform(),
            ha="center", va="top", fontsize=14, color=plt.cm.tab10(1), fontweight="bold")
    ax.text(6, 0.95, "down-type", transform=ax.get_xaxis_transform(),
            ha="center", va="top", fontsize=14, color=plt.cm.tab10(0), fontweight="bold")
    ax.set_ylabel(r"$\log_{10}$ fitted Yukawa singular value")
    ax.grid(axis="x", visible=False)
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    titles(fig, ax, "Fitted Yukawa singular values span the quark-mass hierarchy")
    p = os.path.join(OUT, "fig5_yukawa_hierarchy.png")
    fig.savefig(p); plt.close(fig)
    return p


def main():
    ns, cids, tags = load()
    print(f"loaded {len(ns):,} non-skipped rows")
    p1 = fig1_survival(ns)
    p2, reach = fig2_reach(ns, cids, tags)
    p3 = fig3_zbb_curve(ns)
    p4 = fig4_bulkmass(ns)
    p5 = fig5_yukawa(ns)
    print("REACH RANKING:")
    for c, (rr, t) in sorted(reach.items(), key=lambda x: -x[1][0]):
        print(f"  {c:6s} {rr:5.0f} TeV  [{t}]  {NAMES[c]}")
    for p in (p1, p2, p3, p4, p5):
        sz = os.path.getsize(p)
        print(f"WROTE {p}  ({sz} bytes)")


if __name__ == "__main__":
    main()
