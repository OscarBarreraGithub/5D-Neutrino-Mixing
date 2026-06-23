#!/usr/bin/env python3
"""Standalone 'ours-only' reproduction panels (solo PNGs) for the notes report.

This re-renders OUR panels as their OWN matplotlib figures (no cropped paper
image, no side-by-side). The data + plotting logic are lifted verbatim from
notebooks/_build_sidebyside_vs_papers.py (STU, Zbb) and the existence/census
cells (anarchic_reproduction_vs_papers.ipynb @99664c9 + _constraint_explorer_src.py).
Physics unchanged; only the figure is now standalone with full labels/legend/colorbar.

Covers:  solo_STU, solo_Zbb_gLgR, solo_existence_vs_typical, solo_census.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.lines import Line2D

REPO = Path("/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing")
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
OUT = REPO / "reports" / "collaborator_2026-06" / "figures_solo"
OUT.mkdir(parents=True, exist_ok=True)

mpl.rcParams.update({
    "figure.dpi": 150, "savefig.dpi": 150, "font.size": 11,
    "axes.grid": True, "grid.alpha": 0.25, "axes.titlesize": 12,
    "axes.labelsize": 12, "legend.fontsize": 9, "legend.framealpha": 0.9,
})

def save(fig, name):
    p = OUT / name
    fig.savefig(p, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[saved] {name}  ({p.stat().st_size//1024} KB)")

# ---------------------------------------------------------------------------
# Continuous-M_KK points (STU, Zbb): same load as _build_sidebyside_vs_papers.py
# ---------------------------------------------------------------------------
dfc = pd.read_parquet(REPO / "scan_outputs" / "sidebyside_points.parquet")
if "replay_error" in dfc.columns:
    dfc = dfc[dfc["replay_error"].isna()].copy()
dfc["M_KK_TeV"] = dfc["M_KK_GeV"] / 1000.0
dfc = dfc.sort_values("M_KK_TeV").reset_index(drop=True)
SM_GL_B, SM_GR_B = -0.42114, 0.077420

# ============================================================ 1. solo_STU
# (verbatim physics from _build_sidebyside_vs_papers.py section 1)
FIT_OURS  = dict(s=0.026, t=0.047, ss=0.075, st=0.066, rho=0.90)   # EW001 / PDG-2025
FIT_PAPER = dict(s=0.07,  t=0.16,  ss=0.10,  st=0.10,  rho=0.85)   # CGHNP 2008

def ellipse_patches(fit, color):
    cov = np.array([[fit["ss"]**2, fit["rho"]*fit["ss"]*fit["st"]],
                    [fit["rho"]*fit["ss"]*fit["st"], fit["st"]**2]])
    vals, vecs = np.linalg.eigh(cov)
    ang = np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0]))
    out = []
    for nsig, ls in zip([np.sqrt(2.279), np.sqrt(5.991), np.sqrt(9.210)],
                        ["-", "--", ":"]):
        w, h = 2*nsig*np.sqrt(vals)
        out.append(Ellipse((fit["s"], fit["t"]), w, h, angle=ang,
                           fill=False, edgecolor=color, ls=ls, lw=1.6))
    return out

fig, ax = plt.subplots(figsize=(7.6, 6.6))
for e in ellipse_patches(FIT_OURS, "#1f77b4"):  ax.add_patch(e)
for e in ellipse_patches(FIT_PAPER, "#888888"): ax.add_patch(e)
ax.plot(0, 0, "k+", ms=11, mew=2, zorder=6)
ax.annotate("SM", (0, 0), textcoords="offset points", xytext=(6, 6))
sub = dfc.dropna(subset=["S_pred", "T_pred"])
order = np.argsort(sub["M_KK_TeV"].values)
ax.plot(sub["S_pred"].values[order], sub["T_pred"].values[order],
        "-", color="navy", lw=1.0, alpha=0.5, zorder=4)
sc = ax.scatter(sub["S_pred"], sub["T_pred"], c=sub["M_KK_TeV"], cmap="plasma",
                norm=mpl.colors.LogNorm(), s=12, alpha=0.85, zorder=5)
for mtev in (1, 2, 3, 5, 10):
    row = sub.iloc[(sub["M_KK_TeV"] - mtev).abs().argmin()]
    ax.annotate(f"{mtev} TeV", (row.S_pred, row.T_pred),
                textcoords="offset points", xytext=(6, -2), fontsize=8, color="navy")
cb = fig.colorbar(sc, ax=ax, pad=0.02); cb.set_label(r"$M_{\rm KK}$ [TeV]")
ax.set_xlim(-0.4, 0.6); ax.set_ylim(-0.4, 0.6)
ax.axhline(0, color="k", lw=0.4); ax.axvline(0, color="k", lw=0.4)
ax.set_xlabel("S"); ax.set_ylabel("T")
ax.set_title("RS minimal $S$–$T$ proxy (continuous $M_{\\rm KK}$)\n+ live EW001 / PDG-2025 ellipse")
handles = [Line2D([], [], color="#1f77b4", label="EW001 / PDG-2025 (68/95/99%)"),
           Line2D([], [], color="#888888", label="paper CGHNP-2008 (faint ref.)"),
           Line2D([], [], marker="+", ls="", color="k", label="SM")]
ax.legend(handles=handles, loc="upper left")
save(fig, "solo_STU.png")

# ============================================================ 2. solo_Zbb_gLgR
# (verbatim physics from _build_sidebyside_vs_papers.py section 3)
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
corr = np.array([[1,-0.08,-0.10],[-0.08,1,0.06],[-0.10,0.06,1]])
cov_inv = np.linalg.inv(np.outer(sig, sig)*corr)
gL_ax = np.linspace(-0.44, -0.40, 320); gR_ax = np.linspace(0.05, 0.12, 320)
GLg, GRg = np.meshgrid(gL_ax, gR_ax)
pred = np.stack(Rb_AB_AFB(GLg, GRg), axis=-1)
d = pred - obs
chi2 = np.einsum("...i,ij,...j->...", d, cov_inv, d)
LEV = [2.279, 5.991, 9.210]

fig, ax = plt.subplots(figsize=(7.8, 6.4))
ax.contourf(GLg, GRg, chi2, levels=[0]+LEV,
            colors=["#ffe14d", "#ff9f1c", "#e23b3b"], alpha=0.9)
ax.contour(GLg, GRg, chi2, levels=LEV, colors="k", linewidths=0.4, alpha=0.4)
ax.plot(SM_GL_B, SM_GR_B, "ko", ms=6, zorder=6)
ax.annotate("SM", (SM_GL_B, SM_GR_B), textcoords="offset points", xytext=(6, -12))
win = ((dfc["g_L_b"] >= -0.44) & (dfc["g_L_b"] <= -0.40) &
       (dfc["g_R_b"] >= 0.05) & (dfc["g_R_b"] <= 0.12))
dpl = dfc[win].sort_values("g_L_b")
ax.plot(dpl["g_L_b"], dpl["g_R_b"], "-", color="#3b528b", lw=2.2, alpha=0.7, zorder=4)
sc = ax.scatter(dpl["g_L_b"], dpl["g_R_b"], c=dpl["M_KK_TeV"], cmap="viridis",
                norm=mpl.colors.LogNorm(), s=8, edgecolor="none", alpha=0.9, zorder=5)
cb = fig.colorbar(sc, ax=ax, pad=0.02); cb.set_label(r"$M_{\rm KK}$ [TeV]")
ax.set_xlim(-0.44, -0.40); ax.set_ylim(0.05, 0.12)
ax.set_xlabel(r"$g_L^b$"); ax.set_ylabel(r"$g_R^b$")
ax.set_title("RS $(g_L^b,g_R^b)$ stripe (continuous $M_{\\rm KK}$)\n+ Z-pole CL ellipses (paper Eq. 173 inputs)")
handles = [Line2D([], [], marker="s", ls="", color="#ffe14d", label="68% CL"),
           Line2D([], [], marker="s", ls="", color="#ff9f1c", label="95% CL"),
           Line2D([], [], marker="s", ls="", color="#e23b3b", label="99% CL"),
           Line2D([], [], marker="o", ls="", color="k", label="SM ref"),
           Line2D([], [], marker="o", ls="", color="#3b528b", label="our RS draws")]
ax.legend(handles=handles, loc="upper left")
save(fig, "solo_Zbb_gLgR.png")

# ============================================================ 7. solo_existence_vs_typical
# (verbatim physics from anarchic_reproduction_vs_papers.ipynb @99664c9, cell 24)
dc = pd.read_parquet(REPO / "scan_outputs" / "anarchic_reproduction" / "anarchic_bauer_s1_m12.parquet")
df = pd.read_parquet(REPO / "scan_outputs" / "anarchic_reproduction" / "anarchic_draws.parquet")
sb = pd.read_parquet(REPO / "scan_outputs" / "sidebyside_points.parquet")
sb["M_KK_TeV"] = sb["M_KK_GeV"] / 1000.0  # ratio_EW001 / ratio_T010 are keyed by tile

flav = [("ratio_eps_K", r"$\varepsilon_K$",  "C3"),
        ("ratio_dm_K",  r"$\Delta m_K$",      "C0"),
        ("ratio_B_d",   r"$\Delta m_{B_d}$",  "C2"),
        ("ratio_B_s",   r"$\Delta m_{B_s}$",  "C4"),
        ("ratio_D",     r"$\Delta m_{D^0}$",  "C1")]
flav_tiles = np.array(sorted(set(dc.M_KK_TeV.unique()) | set(df.M_KK_TeV.unique())))

def flav_min_med(col):
    mins, meds, xs = [], [], []
    for m in flav_tiles:
        parts = []
        if col in dc.columns:
            parts.append(dc.loc[dc.M_KK_TeV == m, col])
        if col in df.columns:
            parts.append(df.loc[df.M_KK_TeV == m, col])
        v = pd.concat(parts) if parts else pd.Series([], dtype=float)
        if len(v) == 0:
            continue
        xs.append(m); mins.append(float(v.min())); meds.append(float(v.median()))
    return np.array(xs), np.array(mins), np.array(meds)

ir_bins = np.linspace(1.0, 10.0, 19)
sb["_bin"] = pd.cut(sb.M_KK_TeV, ir_bins)
def irr_min_med(col):
    g = sb.groupby("_bin", observed=True).agg(mid=("M_KK_TeV", "mean"),
                                              mn=(col, "min"), md=(col, "median"))
    return g["mid"].values, g["mn"].values, g["md"].values

def cross1(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    for i in range(1, len(x)):
        if y[i-1] > 1.0 >= y[i]:
            lx0, lx1 = np.log(x[i-1]), np.log(x[i])
            ly0, ly1 = np.log(y[i-1]), np.log(y[i])
            return float(np.exp(lx0 + (0.0 - ly0) * (lx1 - lx0) / (ly1 - ly0)))
    return None

xS, mnS, mdS = irr_min_med("ratio_EW001")
A_min = mnS[-1] * xS[-1]**4
A_med = mdS[-1] * xS[-1]**4
M_ext = np.linspace(xS[-1], 30.0, 200)
STU_min_ext = A_min / M_ext**4
STU_med_ext = A_med / M_ext**4
stu_floor_min = (A_min)**0.25
stu_floor_med = (A_med)**0.25
xZ, mnZ, mdZ = irr_min_med("ratio_T010")

fig, ax = plt.subplots(figsize=(10.5, 7.0))
crossings = {}
for col, lab, c in flav:
    x, mn, md = flav_min_med(col)
    ax.plot(x, md, "-",  color=c, lw=1.8, marker="o", ms=4, label=f"{lab}  (median)")
    ax.plot(x, mn, "--", color=c, lw=1.1, alpha=0.65)
    crossings[lab + " min"] = cross1(x, mn)
    crossings[lab + " median"] = cross1(x, md)
ax.plot(xZ, mdZ, "-",  color="k", lw=2.4, marker="s", ms=4, label=r"$Z\to b\bar b$ (median, irreducible)")
ax.plot(xZ, mnZ, "--", color="k", lw=1.3, alpha=0.8, label=r"$Z\to b\bar b$ (min/existence)")
ax.fill_between(xZ, mnZ, mdZ, color="k", alpha=0.12)
crossings["Z->bb min"] = cross1(xZ, mnZ)
crossings["Z->bb median"] = cross1(xZ, mdZ)
ax.plot(xS, mdS, "-",  color="tab:purple", lw=2.8, marker="D", ms=4,
        label=r"$S,T,U$ (median, irreducible)")
ax.plot(xS, mnS, "--", color="tab:purple", lw=1.6, alpha=0.85,
        label=r"$S,T,U$ (min/existence)")
ax.fill_between(xS, mnS, mdS, color="tab:purple", alpha=0.12)
ax.plot(M_ext, STU_med_ext, "-",  color="tab:purple", lw=1.6, alpha=0.45)
ax.plot(M_ext, STU_min_ext, "--", color="tab:purple", lw=1.2, alpha=0.45)
ax.text(13.5, STU_med_ext[np.argmin(np.abs(M_ext-13.5))]*1.3,
        r"$\propto 1/M_{\rm KK}^4$ (extrapolated)", color="tab:purple", fontsize=8, rotation=-32)
crossings["S,T,U min (existence)"] = stu_floor_min
crossings["S,T,U median (typical)"] = stu_floor_med
ax.axhline(1.0, color="r", lw=1.6, ls="-", zorder=1)
ax.text(1.05, 1.4, "exclusion  (NP/bound = 1)", color="r", fontsize=9)
for M, c, txt, dy in [(stu_floor_med, "tab:purple", "", 1), (20.0, "C3", "", 1)]:
    ax.axvline(M if c == "tab:purple" else 30.0, color=c, ls=":", lw=1.3, alpha=0.7)
ax.scatter([stu_floor_med], [1.0], color="tab:purple", s=90, zorder=6, edgecolor="k")
ax.annotate(f"S,T,U median crosses 1\n@ {stu_floor_med:.0f} TeV  (existence floor)",
            xy=(stu_floor_med, 1.0), xytext=(21, 5.5), color="tab:purple", fontsize=8.5,
            ha="left", arrowprops=dict(arrowstyle="->", color="tab:purple"))
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlim(1.0, 32); ax.set_ylim(1e-5, 1e5)
ax.set_xticks([1, 2, 3, 5, 10, 18, 20, 30])
ax.set_xticklabels(["1", "2", "3", "5", "10", "18", "20", "30"])
ax.set_xlabel(r"$M_{\rm KK}$  [TeV]")
ax.set_ylabel("NP / bound  (per-$M_{\\rm KK}$ min = existence,  median = typical)")
ax.set_title("EXISTENCE (min, best/tuned point) vs TYPICAL (median) constraint floor")
ax.grid(alpha=0.3, which="both")
head = ("Combined EXISTENCE floor $\\approx$18-20 TeV, set by $S,T,U$ (irreducible);\n"
        "flavor tunes away (min crashes below 1 at $\\sim$1 TeV).\n"
        "Typical (median) combined floor $\\sim$30 TeV, set by $\\varepsilon_K$.")
ax.text(0.015, 0.022, head, transform=ax.transAxes, fontsize=9.5, va="bottom", ha="left",
        bbox=dict(boxstyle="round,pad=0.5", fc="#fffbe6", ec="0.4", alpha=0.95))
ax.legend(ncol=2, fontsize=7.6, loc="upper right", framealpha=0.9)
save(fig, "solo_existence_vs_typical.png")
print("  existence/typical crossings:", {k: (round(v,1) if isinstance(v,float) else v)
                                          for k, v in crossings.items()})

# ============================================================ 8. solo_census
# (verbatim physics from _constraint_explorer_src.py: veto_fraction_by_mkk + per_constraint_floor)
cm_raw = pd.read_parquet(REPO / "scan_outputs" / "fix100k_minimal_20260622T080053" / "constraint_matrix.parquet")
cm = cm_raw[cm_raw["evaluated"].fillna(False).astype(bool)].copy()
LEPTON_IDS = ["L001"]
ALL_IDS = [c[len("pass_"):] for c in cm.columns if c.startswith("pass_")]
PRESENT_IDS = [cid for cid in ALL_IDS if cid not in LEPTON_IDS and not cm[f"pass_{cid}"].isna().all()]
LABELS = {
    "K001": r"$\varepsilon_K$", "C001": r"$D^0$ mixing", "C002": r"$D^0$ CPV",
    "B001": r"$\Delta m_d$", "B003": r"$\Delta m_s$", "L001": r"$\mu\to e\gamma$",
    "T010": r"$Z\to b\bar b$", "T011": r"$Z\to b\bar b$ (aux)", "EW001": "S,T,U",
    "COLLIDER": "Collider $M_{\\rm KK}$ cut",
}
MKK_TILES = np.array(sorted(cm["M_KK_TeV"].unique()))

def veto_fraction_by_mkk(cid):
    col = f"pass_{cid}"
    passes = cm[col].fillna(True).to_numpy(dtype=bool) if col in cm.columns else np.ones(len(cm), bool)
    out = []
    for m in MKK_TILES:
        sel = cm["M_KK_TeV"].to_numpy() == m
        n = sel.sum()
        out.append((1.0 - passes[sel].mean()) if n else np.nan)
    return np.array(out)

PER_CONSTRAINT_VETO_THRESHOLD = 0.50
def per_constraint_floor(cid, threshold=PER_CONSTRAINT_VETO_THRESHOLD):
    v = veto_fraction_by_mkk(cid)
    ok = np.where(v < threshold)[0]
    return float(MKK_TILES[ok[0]]) if len(ok) else None

# Physics constraints = the present quark/EW IDs WITHOUT the artificial COLLIDER
# cut (a hard M_KK>=5.5 TeV step, not a flavor/EW veto). COLLIDER is drawn
# separately as a faint dashed reference. (COLLIDER is already inside PRESENT_IDS
# via pass_COLLIDER; appending it again would double-draw and hide K001, whose
# veto curve coincides with it on this evaluated subset.) Dedup + rank by floor.
PHYS_IDS = [cid for cid in dict.fromkeys(PRESENT_IDS) if cid != "COLLIDER"]
ranked = sorted(
    [(cid, per_constraint_floor(cid)) for cid in PHYS_IDS],
    key=lambda kv: -(kv[1] if (kv[1] is not None) else 1e9))

fig, ax = plt.subplots(figsize=(9.6, 6.8))
# faint reference: artificial collider M_KK cut (step, not a physics veto)
if "COLLIDER" in PRESENT_IDS:
    ax.plot(MKK_TILES, veto_fraction_by_mkk("COLLIDER"), "--", color="0.6", lw=1.4,
            alpha=0.7, label=r"collider $M_{\rm KK}$ cut (ref., not a flavor/EW veto)")
for cid, floor in ranked:
    v = veto_fraction_by_mkk(cid)
    flab = f"{LABELS.get(cid, cid)}  (floor {floor:g} TeV)" if floor is not None \
        else f"{LABELS.get(cid, cid)}  (never $\\geq$50%)"
    ax.plot(MKK_TILES, v, "o-", lw=1.9, alpha=0.9, ms=5, label=flab)
ax.axhline(0.5, ls=":", color="grey", lw=1.0)
ax.text(MKK_TILES[0]*1.02, 0.52, "50% veto (floor threshold)", fontsize=8, color="grey")
ax.set_xscale("log")
ax.set_xticks([1, 2, 3, 5, 7, 10, 15, 20, 30, 50])
ax.set_xticklabels(["1", "2", "3", "5", "7", "10", "15", "20", "30", "50"])
ax.set_xlabel(r"$M_{\rm KK}$  [TeV]")
ax.set_ylabel("individual veto fraction  (1 $-$ pass rate)")
ax.set_ylim(-0.02, 1.02)
ax.set_title("Constraint census — our 100k minimal-model scan\n"
             "per-constraint veto fraction vs $M_{\\rm KK}$ (legend ordered by floor)")
ax.legend(fontsize=8, ncol=2, loc="upper right")
ax.grid(alpha=0.3, which="both")
save(fig, "solo_census.png")
print("  census floor ranking:", [(cid, f) for cid, f in ranked])

print("DONE (sidebyside group)")
