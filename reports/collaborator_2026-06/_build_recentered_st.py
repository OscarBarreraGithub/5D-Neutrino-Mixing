#!/usr/bin/env python3
"""STEP 1 of the RS-vs-SM flavor note: recentered S-T figure + C1 oblique floor table.

Plan-only post-processing (NO new scan).  Builds:

  DELIVERABLE A  figures_solo/solo_ST_recentered.png
    - PDG-2025 U-fixed ellipse (solid blue, EW001 anchor)
    - CGHNP-2008 ellipse (faint grey reference)
    - PDG-2025 U-FREE ellipse (dashed), rho(S,T)=0.91 from EW001 anchor file
    - MINIMAL RS S-T trajectory (parquet S_pred/T_pred, colored by M_KK)
    - CUSTODIAL RS S-T trajectory (proxy ew_model='custodial_rs_plr', same M_KK grid)
    - SM at origin

  DELIVERABLE B  oblique_floor_table_C1.json (+ printed markdown)
    2x3 table {minimal_rs, custodial_rs_plr} x {CGHNP-2008, PDG-2025 U-fixed,
    PDG-2025 U-free}: the largest M_KK (TeV) at which the proxy is EXCLUDED at 95%
    (ratio_to_budget crosses 1, passes flips False->True as M_KK rises).

The ellipse helper + FIT dicts mirror _render_solo_sidebyside.py.  The S coefficient
(c_S = 30.0) and proxy conventions are the EW001 catalog anchor
(flavor_catalog/processes/top_higgs_ew/EW001.yaml).
"""
from __future__ import annotations

import json
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
JSON_OUT = REPO / "reports" / "collaborator_2026-06" / "oblique_floor_table_C1.json"

from quarkConstraints.oblique_stu import (  # noqa: E402
    ObliqueSTFit,
    evaluate_rs_oblique_proxy,
    CHI2_2DOF_95,
    DEFAULT_HIGGS_VEV_GEV,
    DEFAULT_SIN2_THETA_W,
    DEFAULT_RS_VOLUME_LOG,
)

mpl.rcParams.update({
    "figure.dpi": 150, "savefig.dpi": 150, "font.size": 11,
    "axes.grid": True, "grid.alpha": 0.25, "axes.titlesize": 12,
    "axes.labelsize": 12, "legend.fontsize": 9, "legend.framealpha": 0.9,
})

# ---------------------------------------------------------------------------
# EW001 anchor / proxy conventions (s_coefficient = 30.0 from EW001.yaml)
# ---------------------------------------------------------------------------
S_COEFFICIENT = 30.0  # flavor_catalog/processes/top_higgs_ew/EW001.yaml warped S coeff

# Three fit ellipses.  ss=sigma_S, st=sigma_T, rho=rho(S,T).
# FIT_UFIXED + FIT_PAPER mirror _render_solo_sidebyside.py exactly.
# FIT_UFREE: PDG-2025 S/T/U-floating fit, rho(S,T)=0.91 read from EW001
#   anchor file (flavor_catalog/references/EW001/pdg_2025_electroweak_stu.txt:27).
FIT_UFIXED = dict(s=0.026, t=0.047, ss=0.075, st=0.066, rho=0.90)   # EW001 / PDG-2025 U-fixed
FIT_PAPER  = dict(s=0.07,  t=0.16,  ss=0.10,  st=0.10,  rho=0.85)   # CGHNP-2008
FIT_UFREE  = dict(s=0.021, t=0.04,  ss=0.096, st=0.12,  rho=0.91)   # PDG-2025 U-free


def to_ObliqueSTFit(fit: dict) -> ObliqueSTFit:
    return ObliqueSTFit(
        s_central=fit["s"], t_central=fit["t"],
        sigma_s=fit["ss"], sigma_t=fit["st"], rho_st=fit["rho"],
        u_fixed=0.0, chi2_budget=CHI2_2DOF_95, confidence_level=0.95,
    )


def ellipse_patches(fit, color):
    """68/95/99% contours (mirrors _render_solo_sidebyside.py)."""
    cov = np.array([[fit["ss"]**2, fit["rho"]*fit["ss"]*fit["st"]],
                    [fit["rho"]*fit["ss"]*fit["st"], fit["st"]**2]])
    vals, vecs = np.linalg.eigh(cov)
    ang = np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0]))
    out = []
    for nsig in (np.sqrt(2.279), np.sqrt(5.991), np.sqrt(9.210)):
        w, h = 2*nsig*np.sqrt(vals)
        out.append((fit["s"], fit["t"], w, h, ang))
    return out


def add_ellipses(ax, fit, color, ls, lw=1.6, alpha=1.0):
    for (cx, cy, w, h, ang) in ellipse_patches(fit, color):
        ax.add_patch(Ellipse((cx, cy), w, h, angle=ang, fill=False,
                             edgecolor=color, ls=ls, lw=lw, alpha=alpha))


# ---------------------------------------------------------------------------
# Load minimal RS trajectory (parquet) and build the custodial one on same grid
# ---------------------------------------------------------------------------
dfc = pd.read_parquet(REPO / "scan_outputs" / "sidebyside_points.parquet")
if "replay_error" in dfc.columns:
    dfc = dfc[dfc["replay_error"].isna()].copy()
dfc["M_KK_TeV"] = dfc["M_KK_GeV"] / 1000.0
dfc = dfc.dropna(subset=["S_pred", "T_pred"]).sort_values("M_KK_TeV").reset_index(drop=True)

# Custodial proxy S,T over the SAME M_KK grid as the minimal trajectory.
m_grid_gev = dfc["M_KK_GeV"].to_numpy()
cust_S, cust_T = [], []
for m in m_grid_gev:
    cmp_ = evaluate_rs_oblique_proxy(
        m_kk_gev=float(m), fit=to_ObliqueSTFit(FIT_UFIXED),
        s_coefficient=S_COEFFICIENT, ew_model="custodial_rs_plr",
    )
    cust_S.append(cmp_.prediction.s)
    cust_T.append(cmp_.prediction.t)
cust_S = np.array(cust_S)
cust_T = np.array(cust_T)

# VALIDATION 1: minimal proxy S,T must match parquet S_pred/T_pred to ~1%.
min_S, min_T = [], []
for m in m_grid_gev:
    cmp_ = evaluate_rs_oblique_proxy(
        m_kk_gev=float(m), fit=to_ObliqueSTFit(FIT_UFIXED),
        s_coefficient=S_COEFFICIENT, ew_model="minimal_rs",
    )
    min_S.append(cmp_.prediction.s)
    min_T.append(cmp_.prediction.t)
min_S = np.array(min_S)
min_T = np.array(min_T)
ds_pct = np.abs(min_S - dfc["S_pred"].to_numpy()) / np.abs(dfc["S_pred"].to_numpy()) * 100
dt_pct = np.abs(min_T - dfc["T_pred"].to_numpy()) / np.abs(dfc["T_pred"].to_numpy()) * 100
print(f"[validate] minimal proxy vs parquet: max |dS|={ds_pct.max():.3f}%  max |dT|={dt_pct.max():.3f}%")
assert ds_pct.max() < 1.0 and dt_pct.max() < 1.0, "minimal proxy disagrees with parquet >1% -> wrong s_coefficient"

# VALIDATION 2: custodial T < minimal T at every M_KK.
assert np.all(cust_T < min_T), "custodial T must be < minimal T at every M_KK"
print(f"[validate] custodial T < minimal T at all {len(cust_T)} grid points: OK "
      f"(custodial T spans {cust_T.min():.4f}..{cust_T.max():.4f}, "
      f"minimal T spans {min_T.min():.4f}..{min_T.max():.4f})")

# ---------------------------------------------------------------------------
# DELIVERABLE A: figure
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8.2, 7.0))

add_ellipses(ax, FIT_PAPER,  "#888888", "-",  lw=1.3, alpha=0.55)  # 2008 faint grey
add_ellipses(ax, FIT_UFIXED, "#1f77b4", "-",  lw=1.6)              # 2025 U-fixed solid blue
add_ellipses(ax, FIT_UFREE,  "#2ca02c", "--", lw=1.5)             # 2025 U-free dashed green

ax.plot(0, 0, "k+", ms=12, mew=2, zorder=8)
ax.annotate("SM", (0, 0), textcoords="offset points", xytext=(7, 6))

# Minimal RS trajectory (parquet), colored by M_KK.
order = np.argsort(dfc["M_KK_TeV"].to_numpy())
ax.plot(dfc["S_pred"].to_numpy()[order], dfc["T_pred"].to_numpy()[order],
        "-", color="navy", lw=1.0, alpha=0.5, zorder=4)
sc = ax.scatter(dfc["S_pred"], dfc["T_pred"], c=dfc["M_KK_TeV"], cmap="plasma",
                norm=mpl.colors.LogNorm(), s=12, alpha=0.85, zorder=5)

# Custodial RS trajectory: same M_KK grid, distinct line/markers.
ax.plot(cust_S[order], cust_T[order], "-", color="#d62728", lw=1.8,
        alpha=0.9, zorder=6)
ax.scatter(cust_S[order][::40], cust_T[order][::40], marker="D", s=22,
           facecolor="#d62728", edgecolor="k", lw=0.4, zorder=7)

# M_KK annotations on both trajectories.
for mtev in (1, 2, 3, 5, 10):
    i = int((dfc["M_KK_TeV"] - mtev).abs().argmin())
    ax.annotate(f"{mtev}", (dfc["S_pred"].iloc[i], dfc["T_pred"].iloc[i]),
                textcoords="offset points", xytext=(5, -3), fontsize=7.5, color="navy")
    ax.annotate(f"{mtev}", (cust_S[i], cust_T[i]),
                textcoords="offset points", xytext=(5, 3), fontsize=7.5, color="#d62728")

cb = fig.colorbar(sc, ax=ax, pad=0.02)
cb.set_label(r"$M_{\rm KK}$ [TeV]  (minimal trajectory color)")
ax.set_xlim(-0.4, 0.6)
ax.set_ylim(-0.4, 0.7)
ax.axhline(0, color="k", lw=0.4)
ax.axvline(0, color="k", lw=0.4)
ax.set_xlabel("S")
ax.set_ylabel("T")
ax.set_title(
    "RS minimal vs custodial $S$-$T$ trajectory + recentered EW fit\n"
    "2008$\\to$2025 recentering tightens the ellipse; custodial $P_{LR}$ relieves "
    "$T$\nand pulls the trajectory back toward the recentered ellipse",
    fontsize=10.5,
)
handles = [
    Line2D([], [], color="#1f77b4", ls="-",  label="PDG-2025 U-fixed (EW001)"),
    Line2D([], [], color="#2ca02c", ls="--", label="PDG-2025 U-free"),
    Line2D([], [], color="#888888", ls="-",  alpha=0.6, label="CGHNP-2008 (ref.)"),
    Line2D([], [], color="navy",    ls="-",  marker="o", ms=4,
           label="minimal RS (parquet)"),
    Line2D([], [], color="#d62728", ls="-",  marker="D", ms=5,
           label="custodial RS $P_{LR}$ (proxy)"),
    Line2D([], [], marker="+", ls="", color="k", label="SM"),
]
ax.legend(handles=handles, loc="upper left", fontsize=8.5)
fig.savefig(OUT / "solo_ST_recentered.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"[saved] solo_ST_recentered.png  "
      f"({(OUT / 'solo_ST_recentered.png').stat().st_size // 1024} KB)")

# ---------------------------------------------------------------------------
# DELIVERABLE B: C1 oblique floor table
# ---------------------------------------------------------------------------
FITS = {
    "CGHNP-2008":       FIT_PAPER,
    "PDG-2025 U-fixed": FIT_UFIXED,
    "PDG-2025 U-free":  FIT_UFREE,
}
MODELS = ("minimal_rs", "custodial_rs_plr")

# Fine M_KK grid 1-40 TeV.  Floor = largest M_KK still EXCLUDED at 95%
# (ratio_to_budget > 1, passes=False).  As M_KK rises the prediction shrinks
# toward the SM, so passes flips False->True once; the floor is the crossing.
M_GRID_TEV = np.linspace(1.0, 40.0, 39001)  # 1 GeV steps


def oblique_floor_tev(ew_model: str, fit: dict) -> float | None:
    obj = to_ObliqueSTFit(fit)
    excluded = np.array([
        not evaluate_rs_oblique_proxy(
            m_kk_gev=float(m * 1000.0), fit=obj,
            s_coefficient=S_COEFFICIENT, ew_model=ew_model,
        ).passes
        for m in M_GRID_TEV
    ])
    idx = np.where(excluded)[0]
    if len(idx) == 0:
        return None  # never excluded on this grid (floor < 1 TeV)
    return float(M_GRID_TEV[idx[-1]])  # largest excluded M_KK


table: dict[str, dict[str, float | None]] = {}
for model in MODELS:
    table[model] = {}
    for fit_name, fit in FITS.items():
        table[model][fit_name] = oblique_floor_tev(model, fit)

# Printed markdown table.
fit_names = list(FITS.keys())
print("\n### C1 oblique M_KK floor (TeV) — largest M_KK excluded at 95%\n")
hdr = "| ew_model \\ fit | " + " | ".join(fit_names) + " |"
sep = "|" + "---|" * (len(fit_names) + 1)
print(hdr)
print(sep)
for model in MODELS:
    cells = []
    for fn in fit_names:
        v = table[model][fn]
        cells.append(f"{v:.2f}" if v is not None else "<1")
    print(f"| {model} | " + " | ".join(cells) + " |")

# JSON dump.
payload = {
    "description": "C1 oblique M_KK floor (TeV) = largest M_KK excluded at 95% by the "
                   "RS oblique proxy (ratio_to_budget>1, passes=False).",
    "s_coefficient": S_COEFFICIENT,
    "s_coefficient_source": "flavor_catalog/processes/top_higgs_ew/EW001.yaml (warped S coeff)",
    "proxy_defaults": {
        "higgs_vev_gev": DEFAULT_HIGGS_VEV_GEV,
        "sin2_theta_w": DEFAULT_SIN2_THETA_W,
        "rs_volume_log": DEFAULT_RS_VOLUME_LOG,
    },
    "chi2_budget_95_2dof": CHI2_2DOF_95,
    "m_grid_tev": {"min": 1.0, "max": 40.0, "n": int(M_GRID_TEV.size)},
    "fits": {
        name: {"s": f["s"], "t": f["t"], "sigma_s": f["ss"], "sigma_t": f["st"],
               "rho_st": f["rho"]}
        for name, f in FITS.items()
    },
    "floors_tev": {model: dict(table[model]) for model in MODELS},
}
JSON_OUT.write_text(json.dumps(payload, indent=2))
print(f"\n[saved] {JSON_OUT.name}")

# ---------------------------------------------------------------------------
# Narrative checks (recentering raises minimal floor; custodial relieves EW).
# ---------------------------------------------------------------------------
mf_2008 = table["minimal_rs"]["CGHNP-2008"]
mf_2025f = table["minimal_rs"]["PDG-2025 U-fixed"]
mf_2025u = table["minimal_rs"]["PDG-2025 U-free"]
print("\n--- narrative checks ---")
if mf_2008 is not None and mf_2025f is not None:
    if mf_2025f > mf_2008:
        print(f"OK recentering raises minimal floor: "
              f"U-fixed {mf_2025f:.2f} > 2008 {mf_2008:.2f} TeV")
    else:
        print(f"*** FLAG *** recentering does NOT raise minimal floor: "
              f"U-fixed {mf_2025f:.2f} <= 2008 {mf_2008:.2f} TeV — narrative BROKEN")
for fn in fit_names:
    mv = table["minimal_rs"][fn]
    cv = table["custodial_rs_plr"][fn]
    mvs = f"{mv:.2f}" if mv is not None else "<1"
    cvs = f"{cv:.2f}" if cv is not None else "<1"
    rel = "OK" if (mv is not None and cv is not None and cv < mv) else "CHECK"
    print(f"[{fn}] custodial {cvs} vs minimal {mvs} TeV -> {rel}")

print("\nDONE (recentered ST / C1 floor)")
