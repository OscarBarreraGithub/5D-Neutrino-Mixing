#!/usr/bin/env python
"""Plots for the Yukawa perturbation study (RS-FLAVOR-ALIGNMENT-2026-07).
Reads perturb_study.npz written by scripts/yukawa_perturbation_study.py."""
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

RUN = Path(".orchestration/runs/RS-FLAVOR-ALIGNMENT-2026-07")
d = np.load(RUN / "perturb_study.npz", allow_pickle=True)
COL = {"flat_typical": "0.55", "flat_tuned": "crimson",
       "nelson_barr": "tab:green", "u2": "tab:blue"}
LAB = {"flat_typical": "flat anarchy (magnitude-suppressed)",
       "flat_tuned": "flat anarchy (phase-tuned)",
       "nelson_barr": "Nelson-Barr (real $Y_d$)",
       "u2": "rank-one / $U(2)$"}

# ---- FIG 1: multiplicative-noise breakdown ----
if any(k.startswith("noise_") for k in d.files):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 5))
    for c in ["flat_typical", "flat_tuned", "nelson_barr", "u2"]:
        k = f"noise_{c}"
        if k not in d.files: continue
        arr = d[k]  # (nsigma, 3): sigma, passfrac, median_ratio
        s, pf, mr = arr[:, 0], arr[:, 1], arr[:, 2]
        ax1.semilogx(s, pf, "o-", color=COL[c], label=LAB[c])
        ax2.loglog(s, mr, "o-", color=COL[c], label=LAB[c])
    ax1.axhline(0.5, color="k", ls=":", lw=0.8)
    ax1.set_xlabel(r"noise level $\sigma$  ($Y_{ij}\to Y_{ij}(1+\sigma z_{ij})$)")
    ax1.set_ylabel(r"fraction still passing $\varepsilon_K$")
    ax1.set_title("How much multiplicative noise before $\\varepsilon_K$ breaks")
    ax1.set_ylim(-0.03, 1.05); ax1.grid(alpha=0.3, which="both"); ax1.legend(fontsize=8)
    ax2.axhline(1.0, color="k", ls=":", lw=0.8); ax2.text(3e-4, 1.2, "$\\varepsilon_K$ bound", fontsize=8)
    ax2.set_xlabel(r"noise level $\sigma$"); ax2.set_ylabel(r"median $\varepsilon_K/\mathrm{bound}$ (how badly it fails)")
    ax2.set_title("Magnitude of failure vs noise")
    ax2.grid(alpha=0.3, which="both"); ax2.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(RUN / "fig_perturb_noise.png", dpi=140)
    print("wrote fig_perturb_noise.png")

# ---- FIG 2: gradient magnitude / tuning radius per class ----
if any(k.startswith("grad_") for k in d.files):
    fig, ax = plt.subplots(figsize=(7, 5))
    classes = [c for c in ["flat_typical", "flat_tuned", "nelson_barr", "u2"] if f"grad_{c}" in d.files]
    gnorms = [np.linalg.norm(d[f"grad_{c}"]) for c in classes]
    ax.bar(range(len(classes)), gnorms, color=[COL[c] for c in classes])
    ax.set_yscale("log")
    ax.set_xticks(range(len(classes))); ax.set_xticklabels([LAB[c] for c in classes], rotation=20, ha="right", fontsize=8)
    ax.set_ylabel(r"$|\nabla_{Y_d}\,\mathrm{Im}\,C_4|$   (inverse tuning radius)")
    ax.set_title("Local sensitivity of $\\varepsilon_K$ to the Yukawas\n(large = tuned/fragile, small = protected)")
    ax.grid(alpha=0.3, axis="y", which="both")
    fig.tight_layout(); fig.savefig(RUN / "fig_perturb_gradient.png", dpi=140)
    print("wrote fig_perturb_gradient.png")

# ---- FIG 3: gradient-field PCA (dangerous-subspace dimension) ----
if any(k.startswith("pca_") for k in d.files):
    fig, ax = plt.subplots(figsize=(7.5, 5))
    for c in ["flat_typical", "nelson_barr", "u2"]:
        k = f"pca_{c}"
        if k not in d.files: continue
        sv = d[k]; v = sv**2; v = v / v.sum()
        peff = (sv**2).sum()**2 / ((sv**2) @ (sv**2))
        ax.plot(range(1, len(v) + 1), v, "o-", color=COL[c],
                label=f"{LAB[c]}  (eff. dim {peff:.1f})")
    ax.set_xlabel("principal component"); ax.set_ylabel("fraction of gradient variance")
    ax.set_title("Directions $\\varepsilon_K$ is sensitive to, across each ensemble\n"
                 "anarchy spreads over many directions; $U(2)$ confines to a few")
    ax.set_yscale("log"); ax.grid(alpha=0.3, which="both"); ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(RUN / "fig_perturb_pca.png", dpi=140)
    print("wrote fig_perturb_pca.png")
print("DONE")
