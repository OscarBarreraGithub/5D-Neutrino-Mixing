#!/usr/bin/env python3
"""Standalone 'ours-only' S2 (aligned) anarchic eps_K panel -> solo_epsK_cloud_S2.png.

Companion to ``_render_solo_anarchic.py`` (which renders the S1 "standard" panel).
This renders the SAME three-colour eps_K cloud for Bauer scenario S2 "aligned"
(0912.1625 Sec 5.2): anarchic Yukawas with Y_max=3 BUT a COMMON right-handed
down-type bulk mass c_d (U(3) flavour symmetry on the RH down sector), which
suppresses the dF=2 corrections to eps_K.  The construction is byte-identical to
the S1 panel except it reads the S2 parquets, so the only physics difference is
the down-sector alignment baked into the scan (scripts/anarchic_bauer_s1.py
common_cd=True for S2).

Colour scheme (Bauer Fig. 4 caption, verbatim):
  grey   = "without Z->bb"  = the full mass+CKM-consistent ensemble (passes_pdg),
                              no Z->bb cut;
  blue   = "with Z->bb"     = passes_pdg AND consistent with both Z->bb couplings
                              at 99% CL (Z99=2.5758);
  orange = "with |eps_K|"   = passes_pdg AND Z->bb AND |eps_K| in [1.2,3.2]e-3.

As in the S1 panel, the WHOLE figure (grey backdrop + 5/50/95% quantiles +
blue/orange overlays) is restricted to the passes_pdg (mass+CKM-reproducing)
subset, because Bauer plots only the ensemble that passes his chi^2/3sigma gate.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

REPO = Path("/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing")
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
OUT = REPO / "reports" / "collaborator_2026-06" / "figures_solo"
OUT.mkdir(parents=True, exist_ok=True)
ARDIR = REPO / "scan_outputs" / "anarchic_reproduction"

mpl.rcParams.update({
    "figure.dpi": 150, "savefig.dpi": 150, "font.size": 11,
    "axes.titlesize": 12, "axes.labelsize": 12,
    "legend.fontsize": 9, "xtick.labelsize": 10, "ytick.labelsize": 10,
})

# experimental inputs (verbatim, same as the S1 panel)
SM_EPS_K = 2.161e-3
WIN_PAPER = (1.2e-3, 3.2e-3)

# --- S2 main ensemble (eps_K cloud + quantiles), gated to passes_pdg ----------
df = pd.read_parquet(ARDIR / "anarchic_bauer_S2.parquet")
if "eps_k_np" not in df.columns and "eps_K_np" in df.columns:
    df["eps_k_np"] = df["eps_K_np"].values
_N_RAW = len(df)
if "passes_pdg" in df.columns:
    df = df[df["passes_pdg"].astype(bool)].copy()
print(f"[epsK_cloud_S2] S2 ensemble: {_N_RAW} -> {len(df)} mass+CKM-consistent "
      f"(passes_pdg, {len(df)/max(_N_RAW,1):.1%})")
if "c_Q3" in df.columns:
    print(f"[epsK_cloud_S2] gated c_Q3 median = {df['c_Q3'].median():+.3f} "
          f"(Bauer Table 1 S2: -0.24 +/- 0.43 -> +0.24 repo convention)")

# Realize self-consistent |eps_K| with a random NP phase (same recipe as the S1
# hybrid in _render_solo_anarchic.py: rng(11) on eps_K_np), so the grey backdrop
# and cyan quantiles come from the FAT S2 ensemble.
df["eps_k_total"] = np.abs(
    SM_EPS_K + df["eps_k_np"].values
    * np.exp(1j * np.random.default_rng(11).uniform(0, 2 * np.pi, len(df))))
tiles = np.array(sorted(df.M_KK_TeV.unique()))


def quantile_curves(frame, col, tiles_, qs=(5, 50, 95)):
    return {q: np.array([np.percentile(frame[frame.M_KK_TeV == m][col], q)
                         for m in tiles_]) for q in qs}


qc = quantile_curves(df, "eps_k_total", tiles)

# --- S2 Z->bb companion (blue/orange) -- gated to passes_pdg ------------------
ZBB = ARDIR / "anarchic_bauer_s2_zbb.parquet"
have_zbb = ZBB.exists()
if have_zbb:
    dz = pd.read_parquet(ZBB)
    if "eps_k_np" not in dz.columns and "eps_K_np" in dz.columns:
        dz["eps_k_np"] = dz["eps_K_np"].values
    phz = np.random.default_rng(11).uniform(0, 2 * np.pi, len(dz))
    dz["eps_k_total"] = np.abs(SM_EPS_K + dz["eps_k_np"].values * np.exp(1j * phz))
    if "passes_pdg" in dz.columns:
        dz_acc = dz[dz["passes_pdg"].astype(bool)].copy()
    else:
        dz_acc = dz
    dz_pass = dz_acc[dz_acc["passes_Zbb"].astype(bool)].copy()
    dz_in = ((dz_pass.eps_k_total >= WIN_PAPER[0]) &
             (dz_pass.eps_k_total <= WIN_PAPER[1]))
    print(f"[epsK_cloud_S2] Z->bb companion: {len(dz)} draws "
          f"({len(dz_acc)} mass+CKM-consistent), {len(dz_pass)} ALSO pass Z->bb "
          f"(99% CL, blue), {int(dz_in.sum())} pass BOTH (orange)")
else:
    print("[epsK_cloud_S2] WARNING: no S2 Z->bb companion parquet; "
          "blue/orange overlays will be omitted.")


def _xjit(mkk, seed, w=0.11):
    r = np.random.default_rng(seed).uniform(-w, w, len(mkk))
    return np.asarray(mkk, float) * (1.0 + r)


fig, ax = plt.subplots(figsize=(8.2, 6.2))
ax.axhspan(WIN_PAPER[0], WIN_PAPER[1], color="green", alpha=0.10, zorder=0)
# GREY backdrop: full mass+CKM-consistent S2 cloud (no Z->bb cut) -- explicit
# 7000/tile subsample (rng 23), MIRRORING the S1 hybrid styling exactly.
_GREY_PER_TILE = 7000
_grey_parts, _rng_sub = [], np.random.default_rng(23)
for _m in tiles:
    _gg = df[df.M_KK_TeV == _m]
    if len(_gg) > _GREY_PER_TILE:
        _gg = _gg.iloc[_rng_sub.choice(len(_gg), _GREY_PER_TILE, replace=False)]
    _grey_parts.append(_gg)
df_grey = pd.concat(_grey_parts, ignore_index=True)
ax.scatter(_xjit(df_grey.M_KK_TeV.values, 3), df_grey.eps_k_total.values, s=5,
           c="0.62", alpha=0.35, rasterized=True, zorder=1, linewidths=0,
           label=r"without $Z\to b\bar b$ (full ensemble)")
if have_zbb and len(dz_pass):
    xjb = _xjit(dz_pass.M_KK_TeV.values, 5)
    ax.scatter(xjb, dz_pass.eps_k_total.values, s=7, c="#1f4fff", alpha=0.6,
               rasterized=True, zorder=3, linewidths=0,
               label=r"with $Z\to b\bar b$ (99\% CL)")
    ax.scatter(xjb[dz_in.values], dz_pass.eps_k_total.values[dz_in.values], s=11,
               c="#ff8c00", alpha=0.95, rasterized=True, zorder=4, linewidths=0,
               label=r"with $|\varepsilon_K|$ (both, $[1.2,3.2]\times10^{-3}$)")
for q, ls, lab in [(95, "-", "95%"), (50, "--", "median"), (5, "-", "5%")]:
    ax.plot(tiles, qc[q], color="#00cccc", ls=ls, lw=1.8, zorder=5,
            label=f"{lab} quantile")
ax.axhline(WIN_PAPER[0], color="green", lw=0.6, ls="--", alpha=0.5, zorder=2)
ax.axhline(WIN_PAPER[1], color="green", lw=0.6, ls="--", alpha=0.5, zorder=2)
ax.fill_between([], [], [], color="green", alpha=0.10,
                label=r"exp. band $[1.2,3.2]\times10^{-3}$")
ax.axhline(SM_EPS_K, color="k", lw=0.8, ls=":")
ax.text(1.05, SM_EPS_K * 1.15, "SM", fontsize=8)
ax.set_yscale("log")
ax.set_xlim(1, 10.5)
ax.set_ylim(1e-7, 1e3)
ax.set_xlabel(r"$M_{\rm KK}$ [TeV]")
ax.set_ylabel(r"$|\varepsilon_K|$")
ax.set_title("Anarchic-forward $|\\varepsilon_K|$ cloud (Bauer-S2 ALIGNED: common $c_d$)\n"
             r"grey = without $Z\to b\bar b$, blue = with $Z\to b\bar b$, "
             r"orange = with $|\varepsilon_K|$ (cf.\ 0912.1625 Fig.\ 4, S2)")
ax.legend(loc="upper right", fontsize=7.5, markerscale=2.2, framealpha=0.92, ncol=1)
ax.grid(alpha=0.22, which="both")
p = OUT / "solo_epsK_cloud_S2.png"
fig.savefig(p, dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"[saved] solo_epsK_cloud_S2.png  ({p.stat().st_size // 1024} KB)")
