#!/usr/bin/env python3
"""Standalone 'ours-only' anarchic-reproduction panels (solo PNGs).

Re-renders OUR panels from notebooks/_build_anarchic_reproduction_vs_papers.py
as their OWN matplotlib figures (no cropped paper image). Physics/data unchanged;
only the figure is standalone with full labels/legend/colorbar.

Covers:  solo_epsK_cloud, solo_consistency, solo_D0_funnel (DECLUTTERED), solo_ReIm_M12.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from matplotlib.colors import LogNorm

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
OURS_CMAP = "viridis"

def save(fig, name):
    p = OUT / name
    fig.savefig(p, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[saved] {name}  ({p.stat().st_size//1024} KB)")

# ---------------------------------------------------------------------------
# Load: Bauer-S1-matched ensemble + complex-M12 ensemble (verbatim from build)
# ---------------------------------------------------------------------------
df = pd.read_parquet(ARDIR / "anarchic_bauer_S1.parquet")
if "eps_k_np" not in df.columns and "eps_K_np" in df.columns:
    df["eps_k_np"] = df["eps_K_np"].values

DATA_CM12 = (ARDIR / "anarchic_bauer_s1_m12.parquet")
if not DATA_CM12.exists():
    DATA_CM12 = ARDIR / "anarchic_complex_m12.parquet"
dc = pd.read_parquet(DATA_CM12)
tiles_c = np.array(sorted(dc.M_KK_TeV.unique()))

# experimental inputs (verbatim)
SM_EPS_K     = 2.161e-3
EXP_EPS_NOW  = 2.228e-3
WIN_PAPER = (1.2e-3, 3.2e-3)
BUDGET_NOW = abs(EXP_EPS_NOW - SM_EPS_K)

# SM mixing amplitudes (verbatim from build) for the M12 ratio planes
import quarkConstraints.deltaf2 as _d2
BETA_D = np.radians(22.0); BETA_S = np.radians(-1.0)
M12_SM_BS = (_d2.DELTA_M_BS_SM / 2.0) * np.exp(2j*BETA_S)
M12_SM_K  = (_d2.DELTA_M_K / 2.0) + 1j*(_d2.EPSILON_K_SM*np.sqrt(2.0)*_d2.DELTA_M_K/_d2.KAPPA_EPSILON)
HBAR_GEVS = 6.582119569e-25; TAU_D0 = 410.3e-15
GAMMA_D = HBAR_GEVS / TAU_D0
X12_EXP_D = 1.00e-2
SIN_PHI_BOUND_D = 0.0022/0.012

def cm12(df_, sys_, mkk):
    g = df_[df_.M_KK_TeV == mkk]
    return g[f"re_m12_{sys_}"].values + 1j*g[f"im_m12_{sys_}"].values

# ============================================================ 3. solo_epsK_cloud
# (verbatim physics from build FIGURE 1)
rng = np.random.default_rng(7)
phi = rng.uniform(0, 2*np.pi, len(df))
df["eps_k_total"] = np.abs(SM_EPS_K + df["eps_k_np"].values * np.exp(1j*phi))
tiles = np.array(sorted(df.M_KK_TeV.unique()))
def quantile_curves(col, qs=(5, 50, 95)):
    return {q: np.array([np.percentile(df[df.M_KK_TeV == m][col], q) for m in tiles]) for q in qs}
qc = quantile_curves("eps_k_total")

fig, ax = plt.subplots(figsize=(8.2, 6.2))
sub = pd.concat([g.sample(min(6000, len(g)), random_state=1)
                 for _, g in df.groupby("M_KK_TeV")], ignore_index=True)
cons = (sub.eps_k_total >= WIN_PAPER[0]) & (sub.eps_k_total <= WIN_PAPER[1])
xj = sub.M_KK_TeV.values * (1 + np.random.default_rng(3).uniform(-0.15, 0.15, len(sub)))
ax.scatter(xj[~cons], sub.eps_k_total.values[~cons], s=2, c="0.6", alpha=0.22,
           rasterized=True, label="anarchic draw (grey)")
ax.scatter(xj[cons], sub.eps_k_total.values[cons], s=2, c="orange", alpha=0.5,
           rasterized=True, label=r"consistent $|\varepsilon_K|$ (orange)")
for q, ls, lab in [(95, "-", "95%"), (50, "--", "median"), (5, "-", "5%")]:
    ax.plot(tiles, qc[q], color="navy", ls=ls, lw=1.8, zorder=5,
            label=f"{lab} quantile")
ax.axhspan(WIN_PAPER[0], WIN_PAPER[1], color="green", alpha=0.12, zorder=0,
           label=r"exp. band $[1.2,3.2]\times10^{-3}$")
ax.axhline(SM_EPS_K, color="k", lw=0.8, ls=":")
ax.text(1.05, SM_EPS_K*1.15, "SM", fontsize=8)
ax.set_yscale("log"); ax.set_xlim(1, 10.5); ax.set_ylim(1e-7, 1e3)
ax.set_xlabel(r"$M_{\rm KK}$ [TeV]"); ax.set_ylabel(r"$|\varepsilon_K|$")
ax.set_title("Anarchic-forward $|\\varepsilon_K|$ fat cloud (Bauer-S1-matched)\n"
             "+ 5/50/95% quantile curves + exp. band")
ax.legend(loc="upper right", fontsize=8, markerscale=4)
ax.grid(alpha=0.25, which="both")
save(fig, "solo_epsK_cloud.png")

# ============================================================ 4. solo_consistency
# (verbatim physics from build FIGURE 3)
df["pass_eps_K_current"] = df["eps_k_np"].values <= BUDGET_NOW
frac_paper = [((df[df.M_KK_TeV == m].eps_k_total >= WIN_PAPER[0]) &
               (df[df.M_KK_TeV == m].eps_k_total <= WIN_PAPER[1])).mean()*100 for m in tiles]
frac_now = [df[df.M_KK_TeV == m].pass_eps_K_current.mean()*100 for m in tiles]

fig, ax = plt.subplots(figsize=(8.0, 6.0))
ax.plot(tiles, frac_paper, "o-", color="orange", lw=2, ms=5,
        label=r"paper-era window $[1.2,3.2]\times10^{-3}$")
ax.plot(tiles, frac_now, "s-", color="crimson", lw=2, ms=5,
        label=r"current gate (NP $\leq$ |exp$-$SM| $\approx6.8\times10^{-5}$)")
ax.axhline(10, color="grey", ls=":", lw=0.8); ax.text(9.0, 11.5, "10%", fontsize=8, color="grey")
ax.set_xlim(1, 10.5); ax.set_ylim(0, 100)
ax.set_xlabel(r"$M_{\rm KK}$ [TeV]")
ax.set_ylabel(r"% of anarchic draws consistent with $|\varepsilon_K|$")
ax.set_title("Percent $|\\varepsilon_K|$-consistent vs $M_{\\rm KK}$\n(paper-era window vs current strict gate)")
ax.legend(fontsize=9, loc="upper left"); ax.grid(alpha=0.3)
save(fig, "solo_consistency.png")

# ============================================================ 5. solo_D0_funnel  (DECLUTTERED)
# Physics identical to build FIGURE 9 (Gedalia funnel: y = x12^NP/x12 = 2|M12^NP|/(Gamma_D x12_exp);
# x = sin(arg M12^NP); funnel = min(1, sinphi_bound/|x|); GMFV/LMFV bands).
# DECLUTTER: instead of 3 overlaid colored clouds, show ONE 2D hexbin density of the
# pooled D cloud (all complex-M12 tiles) for an immediately-readable distribution, plus
# 3 thin median-y guide lines (one per representative M_KK) so the M_KK trend is still shown.
xs = np.linspace(-1, 1, 600)
ymax = np.minimum(1.0, SIN_PHI_BOUND_D/np.maximum(np.abs(xs), 1e-3))
GMFV_TOP = 0.12
LMFV_X0 = 0.85

# DECLUTTER decision: pooling all tiles flattens the cloud into a uniform blanket
# that hides both the funnel and the M_KK trend. Instead show ONE clean density at a
# single representative scale (3 TeV = Bauer S1), so the funnel/bands read clearly;
# then overlay median-y markers for 1/3/10 TeV so the "cloud sinks as M_KK rises"
# trend is still explicit. Physics (x12^NP/x12, sin2sigma, funnel) unchanged.
MKK_REP = 3.0
m12_rep = cm12(dc, "D", MKK_REP)
y_rep = 2*np.abs(m12_rep)/(GAMMA_D*X12_EXP_D)
x_rep = np.sin(np.angle(m12_rep))
inwin = y_rep <= 1.0

fig, ax = plt.subplots(figsize=(9.0, 6.6))
# allowed funnel + bands (drawn UNDER the density; opaque so they show through)
ax.fill_between(xs, 0, ymax, color="#c9cdec", alpha=0.85, zorder=0)
ax.plot(xs, ymax, color="#3a3a7a", lw=1.8, zorder=4)
ax.axhspan(0, GMFV_TOP, color="#e9e520", alpha=0.85, zorder=0.5)
ax.add_patch(Rectangle((LMFV_X0, 0), 1.0-LMFV_X0, GMFV_TOP, facecolor="#d98b9e",
                       edgecolor="none", alpha=0.9, zorder=0.6))
ax.text(-0.95, GMFV_TOP*0.5, "GMFV", fontsize=9, va="center", zorder=5)
ax.text(LMFV_X0-0.02, GMFV_TOP*1.35, "LMFV", fontsize=8, va="center", ha="right", zorder=5)
# ONE density representation of the D cloud at 3 TeV (greyscale hexbin so the
# coloured funnel/bands stay legible underneath); log counts.
hb = ax.hexbin(x_rep[inwin], y_rep[inwin], gridsize=50, cmap="Greys", bins="log",
               mincnt=1, extent=(-1, 1, 0, 1.0), zorder=2, alpha=0.80)
cb = fig.colorbar(hb, ax=ax, pad=0.02)
cb.set_label(r"our $D$ draws / bin at $M_{\rm KK}=3$ TeV (log count)")
# median-y guide lines for 1/3/10 TeV -> the M_KK trend (cloud sinks into funnel)
med_cols = {1.0: "tab:red", 3.0: "tab:orange", 10.0: "tab:green"}
for mkk, c in med_cols.items():
    if mkk not in tiles_c:
        continue
    m12 = cm12(dc, "D", mkk)
    ymed = float(np.median(2*np.abs(m12)/(GAMMA_D*X12_EXP_D)))
    ax.axhline(ymed, color=c, lw=2.2, ls="--", zorder=3)
    ax.text(0.985, min(ymed*1.02, 0.96), f"median $y$, {mkk:g} TeV = {ymed:.3f}",
            color="k", fontsize=8, ha="right", va="bottom", zorder=6,
            bbox=dict(boxstyle="round,pad=0.18", fc=c, ec="none", alpha=0.85))
ax.set_xlim(-1, 1); ax.set_ylim(0, 1.0)
ax.set_xlabel(r"$\sin 2\sigma_D = \sin(\arg M_{12}^{\rm NP})$")
ax.set_ylabel(r"$x_{12}^{\rm NP}/x_{12}$")
ax.set_title("$D^0$-mixing: anarchic cloud on the Gedalia funnel\n"
             r"(density at 3 TeV; median-$y$ trend for $M_{\rm KK}=1,3,10$ TeV)")
# explicit legend explaining every element
handles = [
    Line2D([], [], marker="s", ls="", color="#c9cdec", ms=11, label="allowed region (funnel interior)"),
    Line2D([], [], color="#3a3a7a", lw=2, label=r"funnel boundary ($|y\sin2\sigma|\leq0.18$, $y\leq1$)"),
    Line2D([], [], marker="s", ls="", color="#e9e520", ms=11, label="GMFV band"),
    Line2D([], [], marker="s", ls="", color="#d98b9e", ms=11, label="LMFV box"),
    Line2D([], [], marker="h", ls="", color="0.35", ms=10, label=r"our $D$ cloud (3 TeV density)"),
    Line2D([], [], color="tab:orange", lw=2, ls="--", label=r"median $y$ vs $M_{\rm KK}$ (sinks as $M_{\rm KK}\uparrow$)"),
]
ax.legend(handles=handles, loc="upper left", fontsize=8, framealpha=0.95)
save(fig, "solo_D0_funnel.png")

# ============================================================ 6. solo_ReIm_M12
# (verbatim physics from build FIGURE 10: density-coloured scatter of |Re/Re_SM|,|Im/Im_SM|)
MKK_BQ = 3.0
def density_scatter(ax, x, y, xext, yext, cmap=OURS_CMAP):
    lx, ly = np.log10(x), np.log10(y)
    H, xe, ye = np.histogram2d(lx, ly, bins=80, range=[np.log10(xext), np.log10(yext)])
    ix = np.clip(np.searchsorted(xe, lx) - 1, 0, H.shape[0]-1)
    iy = np.clip(np.searchsorted(ye, ly) - 1, 0, H.shape[1]-1)
    cval = np.clip(H[ix, iy], 1, None)
    order = np.argsort(cval)
    sc = ax.scatter(x[order], y[order], c=cval[order], s=2.5, cmap=cmap,
                    norm=LogNorm(vmin=1, vmax=max(cval.max(), 10)),
                    alpha=0.85, edgecolor="none", rasterized=True)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(*xext); ax.set_ylim(*yext)
    ax.axhline(1, color="grey", lw=0.7); ax.axvline(1, color="grey", lw=0.7)
    return sc

fig, (axK, axS) = plt.subplots(1, 2, figsize=(13.0, 5.6))
# Kaon panel
m12_K = cm12(dc, "K", MKK_BQ)
reK = np.abs(m12_K.real/M12_SM_K.real); imK = np.abs(m12_K.imag/M12_SM_K.imag)
scK = density_scatter(axK, reK, imK, (1e-3, 1e3), (1e-5, 1e5))
axK.set_xlabel(r"$|\mathrm{Re}(M_{12}^K)_{\rm KK}/\mathrm{Re}(M_{12}^K)_{\rm SM}|$")
axK.set_ylabel(r"$|\mathrm{Im}(M_{12}^K)_{\rm KK}/\mathrm{Im}(M_{12}^K)_{\rm SM}|$")
axK.set_title(r"Kaon $M_{12}$ plane ($M_{\rm KK}=3$ TeV)")
cbK = fig.colorbar(scK, ax=axK, pad=0.02); cbK.set_label("point count")
# Bs panel
m12_Bs = cm12(dc, "Bs", MKK_BQ)
reBs = np.abs(m12_Bs.real/abs(M12_SM_BS)); imBs = np.abs(m12_Bs.imag/abs(M12_SM_BS))
scS = density_scatter(axS, reBs, imBs, (1e-5, 1e3), (1e-5, 1e3))
axS.set_xlabel(r"$|\mathrm{Re}(M_{12}^s)_{\rm KK}/(M_{12}^s)_{\rm SM}|$")
axS.set_ylabel(r"$|\mathrm{Im}(M_{12}^s)_{\rm KK}/(M_{12}^s)_{\rm SM}|$")
axS.set_title(r"$B_s$ $M_{12}$ plane ($M_{\rm KK}=3$ TeV)")
cbS = fig.colorbar(scS, ax=axS, pad=0.02); cbS.set_label("point count")
fig.suptitle(r"NP (KK) contribution to $M_{12}$, normalized to SM — kaon (left), $B_s$ (right)",
             y=1.01, fontsize=12)
save(fig, "solo_ReIm_M12.png")

print("DONE (anarchic group)")
