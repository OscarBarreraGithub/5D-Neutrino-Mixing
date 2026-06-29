#!/usr/bin/env python3
"""Regenerate the side-by-side figures/fig_epsK_cloud.png for report.tex.

The previous side-by-side (built by notebooks/_build_anarchic_reproduction_vs_papers.py)
used the OLD 2-colour, separately-subsampled construction, which read as a banded
(column-per-tile) cloud and did NOT carry the faithful Bauer grey/blue/orange Z->bb
colouring.  This compositor stitches the Bauer paper crop (left) with the NOW-dense,
de-banded single-ensemble solo panel (right), so the report side-by-side matches the
notes.pdf solo panel and the dense data exactly.
"""
from pathlib import Path
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

REPO = Path("/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing")
REPORT = REPO / "reports" / "collaborator_2026-06"
PAPER = REPO / "references" / "paper_figures" / "bauer_0912.1625_fig4_epsK.png"
SOLO = REPORT / "figures_solo" / "solo_epsK_cloud.png"
OUT = REPORT / "figures" / "fig_epsK_cloud.png"

img_paper = mpimg.imread(PAPER)
img_solo = mpimg.imread(SOLO)

fig, axes = plt.subplots(1, 2, figsize=(14, 5.6))
axes[0].imshow(img_paper)
axes[0].axis("off")
axes[0].set_title("Bauer 0912.1625 Fig. 4 (S1-S4): $|\\varepsilon_K|$ vs $M_{KK}$",
                  fontsize=10)
axes[1].imshow(img_solo)
axes[1].axis("off")
fig.tight_layout()
fig.savefig(OUT, dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"[saved] {OUT} ({OUT.stat().st_size//1024} KB)")
