#!/usr/bin/env python3
r"""Stitch paper-crop + our-figure into labelled side-by-side composites for the
custodial comparison report (custodial_comparison.tex).

Mirrors the style of notebooks/_build_sidebyside_vs_papers.py: the ACTUAL
published figure (cropped PNG from the paper PDF) gets a "PAPER" subtitle, our
fixed-code reproduction gets an "OURS" subtitle, both normalised to a common
panel height so neither is squashed.

Outputs (into reports/collaborator_2026-06/figures_solo/):
  sxs_STU.png          CGHNP Fig.4 (paper, minimal)   | our solo_ST_recentered
  sxs_Zbb_top.png      CGHNP Fig.8 (paper, minimal)   | our solo_Zbb_gLgR
  sxs_Zbb_cpsw.png     CPSW 0607106 Fig.1 (protected) | CPSW 0701055 Fig.1 (loop)
  sxs_epsK_min.png     Bauer Fig.4-S1 + CFW Fig.1 (paper, minimal) | our solo_epsK_cloud
  sxs_epsK_cust.png    Blanke Fig.9 (paper, CUSTODIAL) | our custodial_epsK_wall_flat
"""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np

REPO = Path("/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing")
PAP = REPO / "references" / "paper_figures"
SOLO = REPO / "reports" / "collaborator_2026-06" / "figures_solo"
OUT = SOLO

mpl.rcParams.update({"figure.dpi": 160, "savefig.dpi": 160})

LBL_KW = dict(fontsize=12, fontweight="bold", ha="center", va="top")
PAPER_C = "#7a3b12"   # warm brown for "PAPER"
OURS_C = "#145a96"    # accent blue for "OURS"


def _imshow(ax, path, title, tcolor):
    img = mpimg.imread(str(path))
    ax.imshow(img)
    ax.axis("off")
    ax.set_title(title, color=tcolor, fontsize=12, fontweight="bold", pad=6)


def crop_top_left_quadrant(path, out_path):
    """Bauer Fig.4 is a 2x2 (S1..S4); keep only the S1 (top-left) panel."""
    img = mpimg.imread(str(path))
    h, w = img.shape[:2]
    # take ~54% height to include the S1 panel's bottom axis + legend (the
    # panels sit slightly above the 50% line), and the full left column width.
    sub = img[: int(0.55 * h), : int(0.52 * w)]
    plt.imsave(str(out_path), sub)
    return out_path


def two_up(left, ltitle, right, rtitle, out_name, figsize=(13.0, 5.4),
           lcolor=PAPER_C, rcolor=OURS_C, suptitle=None):
    fig, (axL, axR) = plt.subplots(1, 2, figsize=figsize)
    _imshow(axL, left, ltitle, lcolor)
    _imshow(axR, right, rtitle, rcolor)
    if suptitle:
        fig.suptitle(suptitle, fontsize=12.5, y=0.995)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.90 if not suptitle else 0.86,
                        bottom=0.01, wspace=0.03)
    p = OUT / out_name
    fig.savefig(p, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print(f"[saved] {p.name}  ({p.stat().st_size//1024} KB)")


def stack_two_up(left_top, lt_t, right_top, rt_t,
                 left_bot, lb_t, right_bot, rb_t,
                 out_name, figsize=(13.0, 8.6), suptitle=None):
    """2 rows x 2 cols: paper(L) ours(R) on each row (or two paper panels)."""
    fig, axs = plt.subplots(2, 2, figsize=figsize)
    _imshow(axs[0, 0], left_top, lt_t, PAPER_C)
    _imshow(axs[0, 1], right_top, rt_t, OURS_C)
    _imshow(axs[1, 0], left_bot, lb_t, PAPER_C)
    _imshow(axs[1, 1], right_bot, rb_t, PAPER_C)
    if suptitle:
        fig.suptitle(suptitle, fontsize=12.5, y=0.998)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.92, bottom=0.01,
                        wspace=0.03, hspace=0.14)
    p = OUT / out_name
    fig.savefig(p, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print(f"[saved] {p.name}  ({p.stat().st_size//1024} KB)")


# ============================================================ (a) S,T,U
two_up(
    PAP / "casagrande_0807.4937_fig4_ST.png",
    "PAPER — CGHNP 0807.4937 Fig. 4  (minimal RS, $S$–$T$)",
    SOLO / "solo_ST_recentered.png",
    "OURS — minimal + custodial $S$–$T$ trajectories",
    "sxs_STU.png", figsize=(13.6, 5.2),
)

# ============================================================ (b) Z->bb top row
two_up(
    PAP / "casagrande_0807.4937_fig8_gLgR.png",
    "PAPER — CGHNP 0807.4937 Fig. 8  (minimal RS, $g_L^b$–$g_R^b$)",
    SOLO / "solo_Zbb_gLgR.png",
    "OURS — $(g_L^b,g_R^b)$ stripe + $Z$-pole CL ellipses",
    "sxs_Zbb_top.png", figsize=(13.6, 5.0),
)

# ============================================================ (b) Z->bb CPSW row
two_up(
    PAP / "carena_0607106_fig1_dgbL_vs_c.png",
    "PAPER — CPSW 0607106 Fig. 1  ($\\delta g_{bL}/g$ vs $c_q$: bidoublet protected)",
    PAP / "carena_0701055_fig1_T_dgbL.png",
    "PAPER — CPSW 0701055 Fig. 1  ($\\delta g_{bL}/g$ vs $\\Delta T$: loop undoes it)",
    "sxs_Zbb_cpsw.png", figsize=(13.6, 4.8),
    rcolor=PAPER_C,
)

# ============================================================ (c) eps_K minimal row
s1 = crop_top_left_quadrant(PAP / "bauer_0912.1625_fig4_epsK.png",
                            OUT / "_bauer_fig4_S1.png")
stack_two_up(
    s1, "PAPER — Bauer 0912.1625 Fig. 4 (S1)  (minimal anarchic $\\varepsilon_K$)",
    SOLO / "solo_epsK_cloud.png", "OURS — Lane-A anarchic $|\\varepsilon_K|$ cloud",
    PAP / "csaki_falkowski_weiler_0804.1954_fig1_epsK.png",
    "PAPER — CFW 0804.1954 Fig. 1  (minimal anarchic $\\mathrm{Im}\\,\\Lambda_{LR}$ vs $m_G$)",
    SOLO / "custodial_epsK_wall_flat.png",
    "OURS — $\\varepsilon_K$ wall (minimal + Blanke custodial + ours)",
    "sxs_epsK_min.png", figsize=(13.6, 9.2),
)

# ============================================================ (c) eps_K custodial row
two_up(
    PAP / "blanke_0809.1073_fig9_epsK_MKK.png",
    "PAPER — Blanke 0809.1073 Fig. 9  (CUSTODIAL model: $\\varepsilon_K$ tuning vs $M_{KK}$)",
    SOLO / "custodial_epsK_wall_flat.png",
    "OURS — $\\varepsilon_K$ floor flat under custodialization",
    "sxs_epsK_cust.png", figsize=(13.6, 5.0),
)

print("DONE (custodial side-by-side composites)")
