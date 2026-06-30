#!/usr/bin/env python3
"""OURS panels for the 2x2-per-observable custodial comparison report.

For each of the three observables (eps_K, S-T-U oblique, Z->bb) this builds the
TWO "ours" panels (no-custodial + custodial) that sit in the right column of the
2x2 grid; the LEFT column is the cropped published figure (paper_figures/).

Outputs (all into figures_solo/):
  eps_K  no-cust  : solo_epsK_cloud.png        (REUSED as-is; not rebuilt here)
  eps_K  cust     : ours_epsK_custodial.png    (same cloud + BOLD "IDENTICAL" annotation)
  S,T,U  no-cust  : ours_ST_minimal.png        (S-T scatter, minimal shoots up in T)
  S,T,U  cust     : ours_ST_custodial.png      (same axes, custodial hugs T=0)
  Z->bb  no-cust  : ours_Zbb_minimal.png       (g_L-g_R minimal stripe + CL ellipses)
  Z->bb  cust     : ours_Zbb_custodial.png     (custodial P_LR proxy collapses to SM)

Physics is lifted verbatim from the existing renderers:
  - S-T ellipse helper + custodial trajectory: _build_recentered_st.py
  - Z->bb CL contours: _render_solo_sidebyside.py section 2
  - eps_K cloud: reuse solo_epsK_cloud.png (figures_solo/), re-render w/ annotation.

The custodial Z->bb P_LR proxy (per CLAUDE.md / the report brief):
    g_L_cust = SM_GL_B + (1/L)*(g_L_b - SM_GL_B)   (b_L shift suppressed by 1/L)
    g_R_cust = SM_GR_B                              (g_R^b zeroed: b_R in (1,1)_{-1/3})
with L = 35 (RS volume).  These COLLAPSE onto SM (~L x closer than minimal).
"""
from __future__ import annotations

import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset

REPO = Path("/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing")
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
OUT = REPO / "reports" / "collaborator_2026-06" / "figures_solo"
OUT.mkdir(parents=True, exist_ok=True)

from quarkConstraints.oblique_stu import (  # noqa: E402
    ObliqueSTFit,
    evaluate_rs_oblique_proxy,
    CHI2_2DOF_95,
)

mpl.rcParams.update({
    "figure.dpi": 150, "savefig.dpi": 150, "font.size": 11,
    "axes.grid": True, "grid.alpha": 0.25, "axes.titlesize": 12,
    "axes.labelsize": 12, "legend.fontsize": 9, "legend.framealpha": 0.9,
})

# RS volume (logarithm of the warp hierarchy) used by the P_LR proxy.
RS_VOLUME_L = 35.0
# SM Z b-bbar couplings (brief-specified reference point).
SM_GL_B, SM_GR_B = -0.42114, 0.077420
# EW001 / PDG-2025 U-fixed oblique fit (recentered ellipse).
FIT_UFIXED = dict(s=0.026, t=0.047, ss=0.075, st=0.066, rho=0.90)
S_COEFFICIENT = 30.0  # EW001.yaml warped S coefficient


def save(fig, name):
    p = OUT / name
    fig.savefig(p, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[saved] {name}  ({p.stat().st_size // 1024} KB)")


def to_ObliqueSTFit(fit: dict) -> ObliqueSTFit:
    return ObliqueSTFit(
        s_central=fit["s"], t_central=fit["t"],
        sigma_s=fit["ss"], sigma_t=fit["st"], rho_st=fit["rho"],
        u_fixed=0.0, chi2_budget=CHI2_2DOF_95, confidence_level=0.95,
    )


def st_ellipse_patches(fit):
    """68/95/99% S-T contours (mirrors _build_recentered_st.py)."""
    cov = np.array([[fit["ss"]**2, fit["rho"]*fit["ss"]*fit["st"]],
                    [fit["rho"]*fit["ss"]*fit["st"], fit["st"]**2]])
    vals, vecs = np.linalg.eigh(cov)
    ang = np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0]))
    out = []
    for nsig in (np.sqrt(2.279), np.sqrt(5.991), np.sqrt(9.210)):
        w, h = 2*nsig*np.sqrt(vals)
        out.append((fit["s"], fit["t"], w, h, ang))
    return out


def add_st_ellipses(ax, fit, color, ls="-", lw=1.7, alpha=1.0):
    for (cx, cy, w, h, ang) in st_ellipse_patches(fit):
        ax.add_patch(Ellipse((cx, cy), w, h, angle=ang, fill=False,
                             edgecolor=color, ls=ls, lw=lw, alpha=alpha))


# ===========================================================================
# Load the minimal-RS side-by-side points (S_pred, T_pred, g_L_b, g_R_b, M_KK).
# ===========================================================================
dfc = pd.read_parquet(REPO / "scan_outputs" / "sidebyside_points.parquet")
if "replay_error" in dfc.columns:
    dfc = dfc[dfc["replay_error"].isna()].copy()
dfc["M_KK_TeV"] = dfc["M_KK_GeV"] / 1000.0
dfc = dfc.dropna(subset=["S_pred", "T_pred"]).sort_values("M_KK_TeV").reset_index(drop=True)

# Custodial oblique proxy over the SAME M_KK grid (REUSE _build_recentered_st method).
m_grid_gev = dfc["M_KK_GeV"].to_numpy()
fit_obj = to_ObliqueSTFit(FIT_UFIXED)
min_S, min_T, cust_S, cust_T = [], [], [], []
for m in m_grid_gev:
    pm = evaluate_rs_oblique_proxy(m_kk_gev=float(m), fit=fit_obj,
                                   s_coefficient=S_COEFFICIENT, ew_model="minimal_rs").prediction
    pc = evaluate_rs_oblique_proxy(m_kk_gev=float(m), fit=fit_obj,
                                   s_coefficient=S_COEFFICIENT, ew_model="custodial_rs_plr").prediction
    min_S.append(pm.s); min_T.append(pm.t)
    cust_S.append(pc.s); cust_T.append(pc.t)
min_S, min_T = np.array(min_S), np.array(min_T)
cust_S, cust_T = np.array(cust_S), np.array(cust_T)

# Validation: minimal proxy reproduces parquet S_pred/T_pred (<1%); custodial T < minimal T.
ds = np.abs(min_S - dfc["S_pred"].to_numpy()) / np.abs(dfc["S_pred"].to_numpy()) * 100
dt = np.abs(min_T - dfc["T_pred"].to_numpy()) / np.abs(dfc["T_pred"].to_numpy()) * 100
assert ds.max() < 1.0 and dt.max() < 1.0, "minimal proxy disagrees w/ parquet >1%"
# Custodial T must genuinely HUG ZERO, not merely be "< minimal" (a large negative
# T would pass that weak check).  The CGHNP Eq.(153) proxy gives a fixed ratio
# T_cust/T_min = -1/(2 L^2) (sign flip + 1/L^2 suppression).  Assert BOTH the
# correct sign/ratio AND that |T_cust| is small in absolute oblique units.
EXPECTED_TCUST_OVER_TMIN = -1.0 / (2.0 * RS_VOLUME_L**2)   # ~ -4.08e-4 for L=35
ratio_T = cust_T / min_T
assert np.allclose(ratio_T, EXPECTED_TCUST_OVER_TMIN, rtol=1e-3), (
    f"custodial T/minimal T = [{ratio_T.min():.2e},{ratio_T.max():.2e}] "
    f"!= expected -1/(2L^2)={EXPECTED_TCUST_OVER_TMIN:.2e}")
assert np.abs(cust_T).max() < 0.02, (
    f"custodial |T| must hug 0 (<0.02); got max|T|={np.abs(cust_T).max():.4f}")
assert np.all(cust_T < 0), "custodial T should be small-NEGATIVE (sign flip vs minimal)"
print(f"[validate] minimal proxy vs parquet: max|dS|={ds.max():.3f}% max|dT|={dt.max():.3f}%")
print(f"[validate] custodial T/minimal T = {ratio_T.mean():.3e} "
      f"(expected {EXPECTED_TCUST_OVER_TMIN:.3e}); max|T_cust|={np.abs(cust_T).max():.4f} (hugs 0)")
# Custodial S must be LEFT UNCHANGED (S is unprotected by P_LR).
assert np.allclose(cust_S, min_S, rtol=1e-9), "custodial S must equal minimal S (S unprotected)"
print(f"[validate] custodial S == minimal S (unprotected): max|dS|={np.abs(cust_S-min_S).max():.2e}")

# Shared S-T axes so the two panels are directly comparable.  The minimal-RS
# trajectory shoots up to T ~ O(10-25) at the lowest M_KK, so NO finite linear
# y-range can contain it; we cap the view at T_VIEW_TOP and add an explicit
# off-frame annotation (below) stating where the trajectory actually goes.  The
# cap is set ABOVE the visible point cluster so nothing is jammed against the
# top spine (the Round-1 clipping bug).
ST_XLIM = (-0.4, 0.6)
ST_YLIM = (-0.4, 0.86)   # headroom above the highest in-frame point (T~0.755)
T_VIEW_TOP = ST_YLIM[1]
T_MAX_DATA = float(dfc["T_pred"].max())   # ~25; reported in the off-frame note
order = np.argsort(dfc["M_KK_TeV"].to_numpy())
MKK_TeV = dfc["M_KK_TeV"].to_numpy()
# Shared colour normalisation for M_KK across both S-T panels.
ST_NORM = mpl.colors.LogNorm(vmin=MKK_TeV.min(), vmax=MKK_TeV.max())


def _st_base(ax):
    add_st_ellipses(ax, FIT_UFIXED, "#1f77b4", "-", lw=1.7)
    ax.plot(0, 0, "k+", ms=12, mew=2, zorder=8)
    ax.annotate("SM", (0, 0), textcoords="offset points", xytext=(7, 6), fontweight="bold")
    ax.set_xlim(*ST_XLIM); ax.set_ylim(*ST_YLIM)
    ax.axhline(0, color="k", lw=0.4); ax.axvline(0, color="k", lw=0.4)
    ax.set_xlabel("S"); ax.set_ylabel("T")


# ----------------------------- OURS S,T,U no-cust (minimal) -----------------------------
fig, ax = plt.subplots(figsize=(7.4, 6.6))
_st_base(ax)
ax.plot(dfc["S_pred"].to_numpy()[order], dfc["T_pred"].to_numpy()[order],
        "-", color="navy", lw=1.0, alpha=0.5, zorder=4)
sc = ax.scatter(dfc["S_pred"], dfc["T_pred"], c=MKK_TeV, cmap="plasma",
                norm=ST_NORM, s=14, alpha=0.9, zorder=5)
for mtev in (1, 2, 3, 5, 10):
    i = int((dfc["M_KK_TeV"] - mtev).abs().argmin())
    ax.annotate(f"{mtev}", (dfc["S_pred"].iloc[i], dfc["T_pred"].iloc[i]),
                textcoords="offset points", xytext=(5, -3), fontsize=8, color="navy")
# "Shoots up" annotation sits in the open RIGHT half of the panel (well clear of
# both the vertical trajectory at S~0.03 AND the relocated lower-right legend).
# Anchor its arrow to an IN-FRAME high-T point: order[2] is at T~25 (off-frame),
# which made the arrow target invisible in Round 1.  Pick the trajectory point
# closest to T=0.55 inside the view box instead.
_in_view = (dfc["T_pred"] < T_VIEW_TOP) & (dfc["S_pred"] < ST_XLIM[1])
_anchor_i = int((dfc["T_pred"][_in_view] - 0.55).abs().idxmin())
ax.annotate("minimal RS:\nshoots UP in $T$\n(volume-enhanced\n$T$ problem)",
            xy=(dfc["S_pred"].iloc[_anchor_i], dfc["T_pred"].iloc[_anchor_i]),
            xytext=(0.255, 0.50), fontsize=10, color="navy", fontweight="bold",
            ha="left", va="center",
            arrowprops=dict(arrowstyle="->", color="navy", lw=1.4))
# Explicit off-frame note: the trajectory does NOT stop at the top spine.
ax.annotate(rf"$\uparrow$ continues to $T\!\approx\!{T_MAX_DATA:.0f}$"
            "\nat $M_{\\rm KK}\\!=\\!1$ TeV (off-frame)",
            xy=(0.035, T_VIEW_TOP), xytext=(0.105, T_VIEW_TOP - 0.02),
            fontsize=8.5, color="navy", va="top", ha="left",
            arrowprops=dict(arrowstyle="-|>", color="navy", lw=1.3))
cb = fig.colorbar(sc, ax=ax, pad=0.02); cb.set_label(r"$M_{\rm KK}$ [TeV]")
ax.set_title("OURS — minimal RS $S$-$T$ (no custodial)\n"
             "+ recentered PDG-2025 $U$-fixed ellipse", fontsize=11)
ax.legend(handles=[
    Line2D([], [], color="#1f77b4", label="PDG-2025 $U$-fixed (68/95/99%)"),
    Line2D([], [], color="navy", marker="o", ms=5, label="minimal RS (ours)"),
    Line2D([], [], marker="+", ls="", color="k", label="SM"),
], loc="lower right", fontsize=8.5)
save(fig, "ours_ST_minimal.png")

# ----------------------------- OURS S,T,U cust (custodial) -----------------------------
fig, ax = plt.subplots(figsize=(7.4, 6.6))
_st_base(ax)
# faint minimal trajectory for reference (so the contrast is explicit)
ax.plot(dfc["S_pred"].to_numpy()[order], dfc["T_pred"].to_numpy()[order],
        "-", color="0.6", lw=0.9, alpha=0.45, zorder=3)
ax.scatter(dfc["S_pred"], dfc["T_pred"], c="0.7", s=6, alpha=0.35, zorder=3)
ax.plot(cust_S[order], cust_T[order], "-", color="#d62728", lw=2.2, alpha=0.95, zorder=6)
sc = ax.scatter(cust_S, cust_T, c=MKK_TeV, cmap="plasma", norm=ST_NORM,
                marker="D", s=24, edgecolor="k", lw=0.3, alpha=0.95, zorder=7)
for mtev in (1, 2, 3, 5, 10):
    i = int((dfc["M_KK_TeV"] - mtev).abs().argmin())
    ax.annotate(f"{mtev}", (cust_S[i], cust_T[i]),
                textcoords="offset points", xytext=(5, 5), fontsize=8, color="#d62728")
ax.annotate(r"custodial $P_{LR}$: $T\to0$"
            "\n(hugs the $S$-axis;\n$\\Delta T\\propto+L\\to-1/L$)",
            xy=(cust_S[order[40]], cust_T[order[40]]),
            xytext=(0.20, 0.30), fontsize=10, color="#d62728", fontweight="bold",
            arrowprops=dict(arrowstyle="->", color="#d62728", lw=1.4))
cb = fig.colorbar(sc, ax=ax, pad=0.02); cb.set_label(r"$M_{\rm KK}$ [TeV]")
ax.set_title("OURS — custodial $P_{LR}$ $S$-$T$\n"
             "$T$ collapses to $\\sim0$; residual set by unprotected $S$", fontsize=11)
ax.legend(handles=[
    Line2D([], [], color="#1f77b4", label="PDG-2025 $U$-fixed (68/95/99%)"),
    Line2D([], [], color="0.6", marker="o", ms=4, label="minimal RS (faint ref.)"),
    Line2D([], [], color="#d62728", marker="D", ms=6, label="custodial $P_{LR}$ (ours)"),
    Line2D([], [], marker="+", ls="", color="k", label="SM"),
], loc="lower right", fontsize=8.5)  # lower-right is empty (T<0); clears faint ref pts
save(fig, "ours_ST_custodial.png")


# ===========================================================================
# Z->bb panels: (g_L^b, g_R^b) plane with Z-pole CL contours + SM point.
# CL machinery copied verbatim from _render_solo_sidebyside.py section 2.
# ===========================================================================
ETA_QCD, ETA_QED, Z_B = 0.9954, 0.9997, 0.997e-3
GL_U, GR_U = 0.34674, -0.15470
GL_D, GR_D = -0.42434, 0.077345
A_E = 0.1462


def Rb_AB_AFB(gL, gR):
    num_light = 4.0*((GL_U**2+GR_U**2)+(GL_D**2+GR_D**2))
    denom_b = ETA_QCD*ETA_QED*((1-6*Z_B)*(gL-gR)**2+(gL+gR)**2)
    Rb = 1.0/(1.0+num_light/denom_b)
    r = (gL+gR)/(gL-gR)
    Ab = (2.0*np.sqrt(1-4*Z_B)*r)/(1-4*Z_B+(1+2*Z_B)*r**2)
    Afb = 0.75*A_E*Ab
    return Rb, Ab, Afb


obs = np.array([0.21629, 0.923, 0.0992]); sig = np.array([0.00066, 0.020, 0.0016])
corr = np.array([[1, -0.08, -0.10], [-0.08, 1, 0.06], [-0.10, 0.06, 1]])
cov_inv = np.linalg.inv(np.outer(sig, sig)*corr)
ZBB_XLIM = (-0.44, -0.40); ZBB_YLIM = (0.05, 0.12)
gL_ax = np.linspace(*ZBB_XLIM, 320); gR_ax = np.linspace(*ZBB_YLIM, 320)
GLg, GRg = np.meshgrid(gL_ax, gR_ax)
pred = np.stack(Rb_AB_AFB(GLg, GRg), axis=-1)
d = pred - obs
chi2 = np.einsum("...i,ij,...j->...", d, cov_inv, d)
LEV = [2.279, 5.991, 9.210]

# minimal stripe = the in-window sidebyside draws (same selection as solo_Zbb_gLgR)
win = ((dfc["g_L_b"] >= ZBB_XLIM[0]) & (dfc["g_L_b"] <= ZBB_XLIM[1]) &
       (dfc["g_R_b"] >= ZBB_YLIM[0]) & (dfc["g_R_b"] <= ZBB_YLIM[1]))
dpl = dfc[win].sort_values("g_L_b").copy()

# custodial P_LR proxy points (brief-specified): b_L shift / L, b_R -> SM (zeroed).
gL_cust = SM_GL_B + (1.0/RS_VOLUME_L)*(dpl["g_L_b"].to_numpy() - SM_GL_B)
gR_cust = np.full(len(dpl), SM_GR_B)
d_min = np.hypot(dpl["g_L_b"].to_numpy() - SM_GL_B, dpl["g_R_b"].to_numpy() - SM_GR_B)
d_cust = np.hypot(gL_cust - SM_GL_B, gR_cust - SM_GR_B)
# The P_LR proxy divides the b_L deviation by L=35 and zeros the (already-tiny)
# b_R deviation, so the mean distance-to-SM ratio is L=35 in the limit g_R^b==SM
# for minimal, and somewhat MORE once the minimal g_R^b spread is included.  A
# bare ">20" tolerated a much weaker collapse; require the ratio to actually be
# consistent with the claimed ~Lx (i.e. >= ~0.9*L), matching the on-figure "35x".
assert RS_VOLUME_L == 35.0  # keep the on-figure "35x" annotation honest
ratio = d_min.mean() / max(d_cust.mean(), 1e-12)
assert ratio >= 0.9 * RS_VOLUME_L, (
    f"custodial Zbb collapse must be ~Lx (>= {0.9*RS_VOLUME_L:.0f}); got {ratio:.1f}x")
print(f"[validate] Z->bb: minimal mean |.-SM|={d_min.mean():.5f}  "
      f"custodial mean={d_cust.mean():.6f}  -> custodial {ratio:.1f}x closer to SM "
      f"(claim 'L={RS_VOLUME_L:.0f}x')")


def _zbb_base(ax):
    ax.contourf(GLg, GRg, chi2, levels=[0]+LEV,
                colors=["#ffe14d", "#ff9f1c", "#e23b3b"], alpha=0.9)
    ax.contour(GLg, GRg, chi2, levels=LEV, colors="k", linewidths=0.4, alpha=0.4)
    ax.plot(SM_GL_B, SM_GR_B, "ko", ms=7, zorder=8)
    ax.annotate("SM", (SM_GL_B, SM_GR_B), textcoords="offset points",
                xytext=(6, -13), fontweight="bold")
    ax.set_xlim(*ZBB_XLIM); ax.set_ylim(*ZBB_YLIM)
    ax.set_xlabel(r"$g_L^b$"); ax.set_ylabel(r"$g_R^b$")


ZBB_NORM = mpl.colors.LogNorm(vmin=dpl["M_KK_TeV"].min(), vmax=dpl["M_KK_TeV"].max())

# ----------------------------- OURS Z->bb no-cust (minimal stripe) -----------------------------
fig, ax = plt.subplots(figsize=(7.6, 6.4))
_zbb_base(ax)
ax.plot(dpl["g_L_b"], dpl["g_R_b"], "-", color="#3b528b", lw=2.2, alpha=0.7, zorder=5)
sc = ax.scatter(dpl["g_L_b"], dpl["g_R_b"], c=dpl["M_KK_TeV"], cmap="viridis",
                norm=ZBB_NORM, s=10, edgecolor="none", alpha=0.9, zorder=6)
ax.annotate("minimal RS:\nstripe sweeps AWAY\nfrom SM as $M_{\\rm KK}$ drops",
            xy=(dpl["g_L_b"].iloc[len(dpl)//2], dpl["g_R_b"].iloc[len(dpl)//2]),
            xytext=(-0.418, 0.112), fontsize=9.5, color="#22305f", fontweight="bold",
            ha="left", arrowprops=dict(arrowstyle="->", color="#22305f", lw=1.3))
cb = fig.colorbar(sc, ax=ax, pad=0.02); cb.set_label(r"$M_{\rm KK}$ [TeV]")
ax.set_title("OURS — minimal RS $(g_L^b,g_R^b)$ stripe (no custodial)\n"
             "+ $Z$-pole 68/95/99% CL ellipses", fontsize=11)
ax.legend(handles=[
    Line2D([], [], marker="s", ls="", color="#ffe14d", label="68% CL"),
    Line2D([], [], marker="s", ls="", color="#ff9f1c", label="95% CL"),
    Line2D([], [], marker="s", ls="", color="#e23b3b", label="99% CL"),
    Line2D([], [], marker="o", ls="", color="k", label="SM"),
    Line2D([], [], marker="o", ls="", color="#3b528b", label="minimal RS (ours)"),
], loc="lower right", fontsize=8.3)
save(fig, "ours_Zbb_minimal.png")

# ----------------------------- OURS Z->bb cust (P_LR proxy collapses to SM) -----------------------------
fig, ax = plt.subplots(figsize=(7.6, 6.4))
# Draw the base WITHOUT the SM marker here, then re-draw SM as a HOLLOW ring at
# LOW zorder so the collapsed custodial diamonds (on top) stay visible — the
# Round-1 bug was the solid black SM dot masking the very points it collapses to.
ax.contourf(GLg, GRg, chi2, levels=[0]+LEV,
            colors=["#ffe14d", "#ff9f1c", "#e23b3b"], alpha=0.9)
ax.contour(GLg, GRg, chi2, levels=LEV, colors="k", linewidths=0.4, alpha=0.4)
# Moderate zoom: enough to show the custodial clump hugging SM AND the minimal
# stripe peeling away, with NO cluttering inset.
ax.set_xlim(-0.4240, -0.4110); ax.set_ylim(0.0650, 0.1120)
ax.set_xlabel(r"$g_L^b$"); ax.set_ylabel(r"$g_R^b$")
# minimal stripe sweeping AWAY from SM (clear grey) + a direction arrow.
ax.plot(dpl["g_L_b"], dpl["g_R_b"], "-", color="0.55", lw=1.6, alpha=0.6, zorder=4)
ax.scatter(dpl["g_L_b"], dpl["g_R_b"], c="0.62", s=16, alpha=0.55, zorder=4)
ax.annotate("", xy=(-0.4122, SM_GR_B), xytext=(-0.4180, SM_GR_B),
            arrowprops=dict(arrowstyle="-|>", color="0.45", lw=1.6))
ax.text(-0.4150, SM_GR_B - 0.0030, "minimal RS sweeps away",
        color="0.35", fontsize=9.5, ha="center", va="top")
# SM as a hollow ring UNDER the custodial diamonds (so they stay visible).
ax.plot(SM_GL_B, SM_GR_B, marker="o", ms=14, mfc="none", mec="k", mew=1.8, zorder=6)
ax.annotate("SM", (SM_GL_B, SM_GR_B), textcoords="offset points",
            xytext=(-23, -3), fontweight="bold", fontsize=11, zorder=9)
# custodial diamonds: a tight bright clump hugging SM.
sc = ax.scatter(gL_cust, gR_cust, c=dpl["M_KK_TeV"].to_numpy(), cmap="viridis",
                norm=ZBB_NORM, s=30, edgecolor="k", lw=0.4, alpha=0.97, zorder=7,
                marker="D")
# one clean label in the open upper-right, no crossing arrow.
ax.text(0.975, 0.96,
        "custodial $P_{LR}$:\n$b$ couplings collapse onto SM\n"
        f"(${ratio:.0f}\\times$ closer than minimal)",
        transform=ax.transAxes, ha="right", va="top", color="#1b5e20",
        fontweight="bold", fontsize=10.5,
        bbox=dict(boxstyle="round,pad=0.4", fc="#eaf6ea", ec="#1b5e20", alpha=0.95))

cb = fig.colorbar(sc, ax=ax, pad=0.02); cb.set_label(r"$M_{\rm KK}$ [TeV]")
ax.set_title("OURS — custodial $P_{LR}$ proxy $(g_L^b,g_R^b)$\n"
             "$b_L$ shift $\\div L$, $b_R\\in(1,1)_{-1/3}\\Rightarrow\\delta g_R^b=0$", fontsize=11)
ax.legend(handles=[
    Line2D([], [], marker="s", ls="", color="#ffe14d", label="68% CL"),
    Line2D([], [], marker="s", ls="", color="#ff9f1c", label="95% CL"),
    Line2D([], [], marker="s", ls="", color="#e23b3b", label="99% CL"),
    Line2D([], [], marker="o", ls="", color="k", label="SM"),
    Line2D([], [], marker="o", ls="", color="0.6", label="minimal (faint ref.)"),
    Line2D([], [], marker="D", ls="", color="#3b528b", label="custodial $P_{LR}$ (ours)"),
], loc="upper left", fontsize=8.0)
save(fig, "ours_Zbb_custodial.png")


# ===========================================================================
# eps_K custodial panel: reuse the existing solo_epsK_cloud.png exactly, but add
# a BOLD annotation that the cloud is UNCHANGED under custodial (KK-gluon C_4 is
# an EW singlet -> blind).  Build by re-loading the PNG as an image and overlaying.
# (The point is that it looks IDENTICAL to solo_epsK_cloud.png + the annotation.)
# ===========================================================================
src_cloud = OUT / "solo_epsK_cloud.png"
assert src_cloud.exists(), f"missing {src_cloud} (run _render_solo_anarchic.py first)"
img = plt.imread(str(src_cloud))
h, w = img.shape[0], img.shape[1]
fig = plt.figure(figsize=(w/150.0, h/150.0))
ax = fig.add_axes([0, 0, 1, 1]); ax.imshow(img); ax.axis("off")
# BOLD banner: custodial leaves the eps_K cloud identical.  Placed in the OPEN
# lower band of the cloud (|eps_K| ~ 1e-5..1e-6 is empty across all M_KK), so it
# does NOT cover the title, the top-right legend, or the dense data band near the
# SM line (the Round-1 bug had it at the very top, over the title+legend).
ax.text(0.52, 0.20,
        "custodial: IDENTICAL cloud\n(KK-gluon $C_4$ is an EW singlet $\\Rightarrow$ blind)",
        transform=ax.transAxes, ha="center", va="center",
        fontsize=15, fontweight="bold", color="#7a2d2d",
        bbox=dict(boxstyle="round,pad=0.5", fc="#fff3f0", ec="#7a2d2d", lw=2.0, alpha=0.96))
fig.savefig(OUT / "ours_epsK_custodial.png", dpi=150)
plt.close(fig)
print(f"[saved] ours_epsK_custodial.png  "
      f"({(OUT / 'ours_epsK_custodial.png').stat().st_size // 1024} KB)  "
      f"[same cloud as solo_epsK_cloud.png + IDENTICAL banner]")

# Sanity: confirm the eps_K 'ours' no-cust panel exists (reused as-is).
assert (OUT / "solo_epsK_cloud.png").exists()
print("[ok] eps_K no-cust panel reused as-is: solo_epsK_cloud.png")

print("\nDONE (ours custodial panels)")
print(f"  custodial S-T hugs T=0 : minimal T up to {min_T.max():.1f}, "
      f"custodial |T|<{np.abs(cust_T).max():.3f}")
print(f"  custodial Z->bb collapses to SM: {ratio:.0f}x closer than minimal stripe")
print(f"  eps_K ours unchanged    : ours_epsK_custodial == solo_epsK_cloud + annotation")
