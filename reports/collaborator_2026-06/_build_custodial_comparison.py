#!/usr/bin/env python3
r"""Custodial ours-vs-literature comparison figures (LEADING custodial physics only).

Builds:
  1. custodial_decoupling_ladder.png -- the headline: per-constraint KK-mass floor,
     MINIMAL vs CUSTODIAL, ours vs literature. The visual point: custodial moves the
     ELECTROWEAK floors (S,T,U and Z->bb) LEFT, but the eps_K wall does NOT move.
  2. custodial_epsK_wall_flat.png -- eps_K floor vs model: minimal (Bauer/CFW/Archer)
     and custodial (Blanke, full DeltaF=2) literature anchors + our Lane-A band, all on
     one axis, showing the wall is flat under custodialization.

CONVENTION: x-axis is the PHYSICAL first KK gauge/gluon mass [TeV] (= our M_KK =
x1*Lambda_IR, x1=2.4487 ~ the literature M_g(1)/m_1^gauge/M_1). Literature numbers are
read in that same physical-mass convention; a residual O(1) convention spread across
papers is absorbed into the plotted BANDS and noted. Custodial ONLY touches S,T,U and
Z->bb; the flavor constraints (eps_K, Dm_K, Dm_Bd/Bs, D0, S_psiphi) are custodial-blind
(KK-gluon, EW-singlet) -- see docs/CUSTODIAL_PROVENANCE.md "CRITICAL" section: the
eps_K-blindness is REAL (Blanke full custodial calc), not a proxy artifact.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Patch
from matplotlib.lines import Line2D

REPO = Path("/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing")
OUT = REPO / "reports" / "collaborator_2026-06" / "figures_solo"
OUT.mkdir(parents=True, exist_ok=True)

mpl.rcParams.update({
    "figure.dpi": 150, "savefig.dpi": 150, "font.size": 11,
    "axes.titlesize": 12, "axes.labelsize": 12,
    "legend.fontsize": 8.5, "xtick.labelsize": 10, "ytick.labelsize": 10,
})

# ---------------------------------------------------------------------------
# Floors in PHYSICAL first-KK-mass [TeV].  (min, max) = plotted band; cite = source.
# OURS: from our scoreboard / C1 table / floor summary (LANE A anarchic where flavor).
# LIT : from the literature catalog (docs/CUSTODIAL_LITERATURE_CATALOG.md).
# ---------------------------------------------------------------------------
ROWS = [
    dict(name=r"$\varepsilon_K$  (anarchic, KK-gluon $C_4$)", kind="flavor",
         ours_min=(9, 24), ours_cust=(9, 24),
         lit_min=(10, 30), lit_cust=(18, 30),
         note="custodial-BLIND (EW-singlet); Blanke full custodial calc confirms"),
    dict(name=r"$S,T,U$  (oblique)", kind="ew",
         ours_min=(16, 20), ours_cust=(6, 6),
         lit_min=(10, 10), lit_cust=(2.5, 3.0),
         note="custodial protects T (not S); S sets the residual floor"),
    dict(name=r"$Z\to b\bar b$", kind="ew",
         ours_min=(5, 5), ours_cust=(2.5, 3.0),
         lit_min=(8.75, 8.75), lit_cust=(0.8, 1.5),
         note="P_LR protects g_L^b (b_L bidoublet)"),
    dict(name=r"$\Delta m_K,\,\Delta m_{B_{d,s}},\,D^0,\,S_{\psi\phi}$"+"\n(other flavor)",
         kind="flavor",
         ours_min=(2, 7), ours_cust=(2, 7),
         lit_min=(2, 7), lit_cust=(2, 7),
         note="custodial-BLIND; sub-dominant to eps_K (not binding)"),
]

# ============================================================ FIG 1: ladder
fig, ax = plt.subplots(figsize=(10.6, 6.4))
yvals = np.arange(len(ROWS))[::-1]   # top row first
C_OURS, C_LIT = "#1f77b4", "#d62728"

def band_marker(ax, y, rng, color, marker, dy, label=None):
    lo, hi = rng
    if hi > lo:
        ax.plot([lo, hi], [y+dy, y+dy], color=color, lw=5, alpha=0.30, solid_capstyle="round")
    xc = np.sqrt(lo*hi)
    ax.plot([xc], [y+dy], marker=marker, color=color, ms=10, mec="k", mew=0.6,
            zorder=5, label=label)
    return xc

for row, y in zip(ROWS, yvals):
    ew = row["kind"] == "ew"
    # OURS minimal -> custodial
    x_om = band_marker(ax, y, row["ours_min"], C_OURS, "o", +0.16)
    x_oc = band_marker(ax, y, row["ours_cust"], C_OURS, "s", +0.16)
    # LIT minimal -> custodial
    x_lm = band_marker(ax, y, row["lit_min"], C_LIT, "o", -0.16)
    x_lc = band_marker(ax, y, row["lit_cust"], C_LIT, "s", -0.16)
    # arrows minimal(circle) -> custodial(square)
    for (xa, xb, yy, col) in [(x_om, x_oc, y+0.16, C_OURS), (x_lm, x_lc, y-0.16, C_LIT)]:
        if ew and xb < xa*0.97:
            ax.add_patch(FancyArrowPatch((xa, yy), (xb, yy), arrowstyle="-|>",
                         mutation_scale=13, color=col, lw=1.6, zorder=4, alpha=0.85))
        else:
            # flavor: no shift -> draw an "=" marker between
            ax.text(np.sqrt(xa*xb), yy, "=", color=col, fontsize=12, ha="center",
                    va="center", fontweight="bold", zorder=6)
    # right-edge note
    ax.text(46, y, row["note"], fontsize=7.6, va="center", ha="left", color="0.25")

# shaded "custodial brings EW into reach" zone (2-7 TeV)
ax.axvspan(2, 7, color="green", alpha=0.06, zorder=0)
ax.text(np.sqrt(2*7), -0.45, "few-TeV (reachable)", color="green",
        fontsize=7.5, ha="center", va="center", alpha=0.85)
ax.axvspan(18, 40, color="grey", alpha=0.08, zorder=0)
ax.text(26, -0.45, r"$\varepsilon_K$ / existence wall", color="0.4",
        fontsize=7.5, ha="center", va="center")

ax.set_yticks(yvals)
ax.set_yticklabels([r["name"] for r in ROWS])
ax.set_xscale("log")
ax.set_xlim(1, 95)
ax.set_xticks([1, 2, 3, 5, 7, 10, 20, 30])
ax.set_xticklabels(["1", "2", "3", "5", "7", "10", "20", "30"])
ax.set_xlabel(r"KK-scale floor  $M_{\rm KK}$ (physical first KK mass) [TeV]")
ax.set_ylim(-0.6, len(ROWS)-0.4)
ax.set_title("Custodial decoupling ladder: custodial moves the EW floors LEFT, "
             "the $\\varepsilon_K$ wall does NOT move\n"
             "(circle = minimal $\\to$ square = custodial; arrow = shift, '=' = custodial-blind)")
ax.grid(axis="x", alpha=0.25, which="both")
handles = [
    Line2D([], [], marker="o", ls="", color=C_OURS, mec="k", label="ours (minimal)"),
    Line2D([], [], marker="s", ls="", color=C_OURS, mec="k", label="ours (custodial)"),
    Line2D([], [], marker="o", ls="", color=C_LIT, mec="k", label="literature (minimal)"),
    Line2D([], [], marker="s", ls="", color=C_LIT, mec="k", label="literature (custodial)"),
]
ax.legend(handles=handles, loc="lower right", ncol=2, framealpha=0.95)
fig.subplots_adjust(right=0.74)
p = OUT / "custodial_decoupling_ladder.png"
fig.savefig(p, dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"[saved] {p.name} ({p.stat().st_size//1024} KB)")

# ============================================================ FIG 2: eps_K wall flat
fig, ax = plt.subplots(figsize=(8.8, 5.2))
anchors = [
    ("Bauer II 2009 (minimal, S1 10%)",        10, "min",  "#ff7f0e"),
    ("CFW 2008 (minimal/RS anarchic)",          21, "min",  "#ff7f0e"),
    ("CFW 2008 (pGB/GHU)",                      33, "min",  "#ff7f0e"),
    ("Archer 2012 (minimal, tuning band)",  (10,30),"min",  "#ff7f0e"),
    ("Blanke 2008 (CUSTODIAL model, dT<20)",    18, "cust", "#1f77b4"),
    ("Blanke 2008 (CUSTODIAL model, dT<10)",    30, "cust", "#1f77b4"),
    ("OURS Lane-A (anarchic, |Vcb| band)",  (9,24), "ours", "#2ca02c"),
]
for i, (lab, val, grp, col) in enumerate(anchors):
    y = len(anchors)-1-i
    if isinstance(val, tuple):
        ax.plot([val[0], val[1]], [y, y], color=col, lw=6, alpha=0.4, solid_capstyle="round")
        ax.plot([np.sqrt(val[0]*val[1])], [y], "D", color=col, ms=9, mec="k", mew=0.6)
    else:
        ax.plot([val], [y], "D", color=col, ms=10, mec="k", mew=0.6)
    ax.text(0.93, y, lab, ha="right", va="center", fontsize=8.4, transform=ax.get_yaxis_transform())
ax.axvspan(9, 33, color="grey", alpha=0.10)
ax.text(np.sqrt(9*33), -0.45, "the wall: ~9-33 TeV in BOTH minimal AND custodial",
        ha="center", va="center", fontsize=8.5, color="0.3")
ax.set_xscale("log")
ax.set_xlim(5, 60)
ax.set_xticks([5, 7, 10, 20, 30, 50]); ax.set_xticklabels(["5","7","10","20","30","50"])
ax.set_yticks([])
ax.set_ylim(-0.6, len(anchors)-0.3)
ax.set_xlabel(r"$\varepsilon_K$ floor  $M_{\rm KK}$ (physical KK-gluon mass) [TeV]")
ax.set_title(r"The $\varepsilon_K$ wall is FLAT under custodialization"+"\n"
             r"(orange=minimal, blue=full custodial model [Blanke], green=ours) "
             r"$\Rightarrow$ custodial is an EW cure, not a flavor cure")
ax.grid(axis="x", alpha=0.25, which="both")
p2 = OUT / "custodial_epsK_wall_flat.png"
fig.savefig(p2, dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"[saved] {p2.name} ({p2.stat().st_size//1024} KB)")
print("DONE (custodial comparison)")
