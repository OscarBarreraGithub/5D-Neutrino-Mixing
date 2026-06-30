#!/usr/bin/env python3
"""Stitch the 2x2-per-observable composites for the custodial comparison report.

For each observable we lay out FOUR panels on a strict 2x2 grid:

         |   PAPER             |   OURS
  -------+---------------------+---------------------
  no     |  paper (no cust.)   |  ours  (no cust.)
  cust.  |                     |
  -------+---------------------+---------------------
  cust-  |  paper (cust.)      |  ours  (cust.)
  odial  |                     |

Column headers (PAPER / OURS) sit across the top in bold; row labels
(no custodial / custodial) run down the left edge in bold.  Each panel keeps its
native aspect ratio (centered in its cell) so neither the paper crop nor our
render is squashed.

Outputs (into figures_solo/):  cmp2x2_epsK.png, cmp2x2_STU.png, cmp2x2_Zbb.png
"""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

REPO = Path("/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing")
SOLO = REPO / "reports" / "collaborator_2026-06" / "figures_solo"
PAPER = REPO / "references" / "paper_figures"

ACCENT = "#145a96"   # OURS column / blue
WARN = "#964614"     # PAPER column / orange-brown

# Each observable: (output name, [ (path, cell-row, cell-col, short-tag) ... ]).
# row 0 = "no custodial", row 1 = "custodial"; col 0 = PAPER, col 1 = OURS.
OBSERVABLES = {
    "cmp2x2_epsK.png": dict(
        title=r"$\varepsilon_K$  (neutral-kaon CP violation)",
        panels=[
            (PAPER / "bauer_0912.1625_fig4_epsK.png", 0, 0,
             "Bauer 0912.1625 Fig.4 (S1-S4 scenarios) — minimal anarchic"),
            (SOLO / "solo_epsK_cloud.png", 0, 1,
             "ours — Lane-A anarchic cloud (cf. S1)"),
            (PAPER / "blanke_0809.1073_fig9_epsK_MKK.png", 1, 0,
             r"Blanke 0809.1073 Fig.9 — FULL custodial ($\Delta_{BG}$ vs $M_{KK}$, not a cloud)"),
            (SOLO / "ours_epsK_custodial.png", 1, 1,
             "ours — SAME cloud (custodial-blind)"),
        ],
    ),
    "cmp2x2_STU.png": dict(
        title=r"$S,T,U$  (electroweak oblique)",
        panels=[
            (PAPER / "casagrande_0807.4937_fig4_ST.png", 0, 0,
             "CGHNP 0807.4937 Fig.4 — minimal RS $S$-$T$"),
            (SOLO / "ours_ST_minimal.png", 0, 1,
             "ours — minimal $S$-$T$ (shoots up in $T$)"),
            (PAPER / "carena_0701055_fig1_T_dgbL.png", 1, 0,
             r"CPSW 0701055 Fig.1 — custodial $T\to0$ ($\delta g_{bL}$ vs $\Delta T$ plane)"),
            (SOLO / "ours_ST_custodial.png", 1, 1,
             r"ours — custodial $P_{LR}$ ($T$ hugs 0)"),
        ],
    ),
    "cmp2x2_Zbb.png": dict(
        title=r"$Z\to b\bar b$  ($b$-quark $Z$ couplings)",
        panels=[
            (PAPER / "casagrande_0807.4937_fig8_gLgR.png", 0, 0,
             "CGHNP 0807.4937 Fig.8 — minimal $g_L$-$g_R$"),
            (SOLO / "ours_Zbb_minimal.png", 0, 1,
             "ours — minimal stripe (sweeps from SM)"),
            (PAPER / "carena_0607106_fig1_dgbL_vs_c.png", 1, 0,
             r"CPSW 0607106 Fig.1 — bidoublet $P_{LR}$ protected ($\delta g_{bL}/g_{bL}$ vs $c_q$)"),
            (SOLO / "ours_Zbb_custodial.png", 1, 1,
             "ours — custodial $P_{LR}$ (collapses to SM)"),
        ],
    ),
}

# Grid geometry (figure fractions).  Left/top margins host the row/column labels.
LEFT_MARGIN = 0.050      # row-label strip
TOP_MARGIN = 0.092       # title + column-header strip
RIGHT_PAD = 0.012
BOT_PAD = 0.012
GUTTER_X = 0.020         # space between the two columns
GUTTER_Y = 0.052         # space between the two rows (room for per-panel sub-tag)
TAG_H = 0.030            # height reserved under each panel for its source tag


def cell_box(row, col):
    """Return (x0, y0, w, h) of the drawing area for grid cell (row, col)."""
    avail_w = 1.0 - LEFT_MARGIN - RIGHT_PAD
    avail_h = 1.0 - TOP_MARGIN - BOT_PAD
    cw = (avail_w - GUTTER_X) / 2.0
    ch = (avail_h - GUTTER_Y) / 2.0
    x0 = LEFT_MARGIN + col * (cw + GUTTER_X)
    # row 0 is the TOP row -> highest y.
    y_top = 1.0 - TOP_MARGIN
    y0 = y_top - (row + 1) * ch - row * GUTTER_Y
    return x0, y0, cw, ch


def place_image(fig, path, cell, tag):
    """Place an image centered in its cell (native aspect), with a tag below."""
    x0, y0, cw, ch = cell
    # Reserve a strip at the bottom of the cell for the source tag.
    img_y0 = y0 + TAG_H
    img_h = ch - TAG_H
    img = plt.imread(str(path))
    ih, iw = img.shape[0], img.shape[1]
    ar = iw / ih
    # Fit image into (cw, img_h) preserving aspect ratio, centered.
    cell_ar = cw / img_h
    if ar > cell_ar:          # image is wider -> width-limited
        draw_w = cw
        draw_h = cw / ar
    else:                     # image is taller -> height-limited
        draw_h = img_h
        draw_w = img_h * ar
    dx0 = x0 + (cw - draw_w) / 2.0
    dy0 = img_y0 + (img_h - draw_h) / 2.0
    ax = fig.add_axes([dx0, dy0, draw_w, draw_h])
    ax.imshow(img)
    ax.axis("off")
    # source tag centered under the image
    fig.text(x0 + cw / 2.0, y0 + TAG_H * 0.45, tag, ha="center", va="center",
             fontsize=8.0, color="0.25", style="italic")


def build(out_name, spec):
    fig = plt.figure(figsize=(13.6, 12.4), dpi=150)
    fig.patch.set_facecolor("white")

    # ---- title across the very top ----
    fig.text(0.5 + LEFT_MARGIN / 2.0, 0.992, "Custodial 2$\\times$2:  " + spec["title"],
             ha="center", va="top", fontsize=17, fontweight="bold", color="#222222")

    # ---- column headers (PAPER / OURS) ----
    for col, (label, color) in enumerate([("PAPER", WARN), ("OURS", ACCENT)]):
        x0, _, cw, _ = cell_box(0, col)
        fig.text(x0 + cw / 2.0, 1.0 - TOP_MARGIN + 0.022, label,
                 ha="center", va="center", fontsize=18, fontweight="bold",
                 color="white",
                 bbox=dict(boxstyle="round,pad=0.5", fc=color, ec="none"))

    # ---- row labels (no custodial / custodial) down the left edge ----
    for row, label in enumerate(["no custodial", "custodial"]):
        _, y0, _, ch = cell_box(row, 0)
        fig.text(LEFT_MARGIN * 0.42, y0 + ch / 2.0, label,
                 ha="center", va="center", rotation=90,
                 fontsize=15, fontweight="bold", color="#333333")

    # ---- faint frames around each row to reinforce the no-cust / cust split ----
    for row, edge in enumerate(["#bcbcbc", "#bcbcbc"]):
        x0a, y0a, cwa, cha = cell_box(row, 0)
        x0b, _, cwb, _ = cell_box(row, 1)
        full_w = (x0b + cwb) - x0a
        rect = FancyBboxPatch(
            (x0a - 0.006, y0a - 0.004), full_w + 0.012, cha + 0.008,
            boxstyle="round,pad=0.004", transform=fig.transFigure,
            fc="none", ec=edge, lw=1.0, zorder=0)
        fig.add_artist(rect)

    # ---- place the four panels ----
    for path, row, col, tag in spec["panels"]:
        assert path.exists(), f"missing panel: {path}"
        place_image(fig, path, cell_box(row, col), tag)

    out = SOLO / out_name
    fig.savefig(out, dpi=150, facecolor="white")
    plt.close(fig)
    print(f"[saved] {out_name}  ({out.stat().st_size // 1024} KB)")


if __name__ == "__main__":
    for name, spec in OBSERVABLES.items():
        build(name, spec)
    print("DONE (cmp2x2 composites)")
