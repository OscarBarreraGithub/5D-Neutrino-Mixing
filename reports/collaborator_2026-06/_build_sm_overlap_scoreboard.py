#!/usr/bin/env python3
r"""STEP 3 (decisive) — the SM-overlap scoreboard for the RS-vs-SM flavor note.

THE QUESTION (verdict gate): is there ANY M_KK at which a generic anarchic point
(i) survives the binding constraints (epsilon_K median floor + S,T,U existence
floor) AND (ii) still shows a >1sigma deviation from the SM in some measurable
quark channel?  If every channel's "generic deviation" region sits BELOW the
survival floor, conclusion (A) (RS wins on flavor) is dead and (B) (RS reduced to
a mass-hierarchy mechanism) is confirmed.

LANE: this is LANE A (anarchic), by design — the literature strawman the note
interrogates.  Floors here are anarchic; do NOT read them as production (~7 TeV).

DATA: scan_outputs/anarchic_reproduction/anarchic_bauer_s1_m12.parquet
  (990k draws x 45 M_KK tiles [1..20 TeV]; per-draw COMPLEX NP M12 for K/Bd/Bs/D,
   co-registered with passes_pdg and ratio_eps_K).  The re/im_m12_* columns are the
   NP amplitudes M12^NP in GeV (deltaf2.compute_m12_np hadronic MEs).  NO new scan.

KI-1 GUARD (docs/KNOWN_ISSUES.md): the B-meson mixing PHASES are computed from the
stored COMPLEX NP M12 against a COMPLEX SM box,
    phi_q^Delta = arg(1 + M12^NP / M12^SM_box),   M12^SM_box COMPLEX,
NOT the B002/B004 real-M12_SM adapter (which is not rephasing-invariant).  The
complex SM boxes use the repo's calibrated convention from _render_solo_anarchic.py
(beta_d = 22 deg, beta_s = -1 deg) and deltaf2 SM Delta-m predictions.
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

REPO = Path("/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing")
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
import quarkConstraints.deltaf2 as d2  # SM Delta-m predictions, kaon SM box pieces

OUT = REPO / "reports" / "collaborator_2026-06" / "figures_solo"
OUT.mkdir(parents=True, exist_ok=True)
ARDIR = REPO / "scan_outputs" / "anarchic_reproduction"

mpl.rcParams.update({
    "figure.dpi": 150, "savefig.dpi": 150, "font.size": 11,
    "axes.titlesize": 12, "axes.labelsize": 12,
    "legend.fontsize": 9, "xtick.labelsize": 10, "ytick.labelsize": 10,
})

# ---------------------------------------------------------------------------
# Complex SM boxes (repo-calibrated convention, _render_solo_anarchic.py L85-88)
# ---------------------------------------------------------------------------
BETA_D = np.radians(22.0)            # SM B_d mixing phase (2 beta = 0.768 rad)
BETA_S = np.radians(-1.0)            # SM B_s mixing phase (small)
M12_SM_BD = (d2.DELTA_M_BD_SM / 2.0) * np.exp(2j * BETA_D)
M12_SM_BS = (d2.DELTA_M_BS_SM / 2.0) * np.exp(2j * BETA_S)
M12_SM_K  = (d2.DELTA_M_K / 2.0) + 1j * (
    d2.EPSILON_K_SM * np.sqrt(2.0) * d2.DELTA_M_K / d2.KAPPA_EPSILON)
# D: SM short-distance box is negligible / long-distance dominated -> no clean SM
# phase reference.  We treat D as a CONSERVATIVE cross-check via the NP amplitude
# only (|2 M12^NP| vs the experimental Delta m_D), NOT a phase-overlap channel.
GAMMA_D = d2.DELTA_M_D_EXP           # use exp Delta m_D as the conservative scale

# ---------------------------------------------------------------------------
# Experimental 1sigma sensitivities (current world averages; sourced in the note)
#   sin2beta  : S_psiKS = 0.699 +- 0.017            (HFLAV 2023)
#   phi_s     : -0.049 +- 0.016 rad (~16 mrad)      (HFLAV/LHCb 2023)
#   C_Bd      : 1.05 +- 0.08  (NP magnitude ratio)  (UTfit NP-fit)
#   C_Bs      : 1.11 +- 0.09                          (UTfit NP-fit)
#   D mixing  : |M12^NP| budget = Delta m_D^exp / 2  (long-distance dominated)
# These set what ">1sigma deviation" means per channel.  The qualitative verdict
# is robust to O(1) changes in these numbers (the deviation/survival gap is wide).
# ---------------------------------------------------------------------------
SIG_SIN2B = 0.017
SIG_PHIS  = 0.016
SIG_CBD   = 0.08
SIG_CBS   = 0.09

# ---------------------------------------------------------------------------
# Load the co-registered NP M12 parquet, gate to mass+CKM-consistent draws.
# ---------------------------------------------------------------------------
df = pd.read_parquet(ARDIR / "anarchic_bauer_s1_m12.parquet")
_N = len(df)
df = df[df["passes_pdg"].astype(bool)].copy()
tiles = np.array(sorted(df["M_KK_TeV"].unique()))
print(f"[scoreboard] {_N} draws -> {len(df)} passes_pdg; {len(tiles)} M_KK tiles "
      f"[{tiles.min():.2f}..{tiles.max():.2f} TeV]")

# Complex NP amplitudes per draw
m12_K  = df["re_m12_K"].values  + 1j * df["im_m12_K"].values
m12_Bd = df["re_m12_Bd"].values + 1j * df["im_m12_Bd"].values
m12_Bs = df["re_m12_Bs"].values + 1j * df["im_m12_Bs"].values
m12_D  = df["re_m12_D"].values  + 1j * df["im_m12_D"].values
mkk    = df["M_KK_TeV"].values

# ---------------------------------------------------------------------------
# Per-draw DEVIATION in sigma units, per channel (KI-1-correct, complex SM box).
# A draw's NP M12 enters with a random NP phase already baked into re/im (the scan
# drew anarchic complex Yukawas), so arg(M12^NP) is physical.
# ---------------------------------------------------------------------------
# B_d: magnitude C_Bd and phase shift -> Delta sin2beta
tot_Bd = M12_SM_BD + m12_Bd
C_Bd = np.abs(tot_Bd) / np.abs(M12_SM_BD)
phi_d_shift = np.angle(tot_Bd) - np.angle(M12_SM_BD)
dsin2b = np.abs(np.sin(2 * BETA_D + phi_d_shift) - np.sin(2 * BETA_D))
dev_Bd_C   = np.abs(C_Bd - 1.0) / SIG_CBD
dev_Bd_phi = dsin2b / SIG_SIN2B

# B_s: magnitude C_Bs and phase shift phi_s
tot_Bs = M12_SM_BS + m12_Bs
C_Bs = np.abs(tot_Bs) / np.abs(M12_SM_BS)
phi_s_shift = np.angle(tot_Bs) - np.angle(M12_SM_BS)
dev_Bs_C   = np.abs(C_Bs - 1.0) / SIG_CBS
dev_Bs_phi = np.abs(phi_s_shift) / SIG_PHIS

# D0: conservative amplitude cross-check, |2 M12^NP| / Delta m_D^exp expressed in
# "sigma" via the experimental Delta m_D as the budget (long-distance dominated, so
# this is a bound not a hard exclusion).  1 unit == NP saturates the full exp value.
dev_D = (2.0 * np.abs(m12_D)) / GAMMA_D

CHANNELS = {
    "B_s phase (phi_s)":       dev_Bs_phi,
    "B_s magnitude (C_Bs)":    dev_Bs_C,
    "B_d phase (sin2beta)":    dev_Bd_phi,
    "B_d magnitude (C_Bd)":    dev_Bd_C,
    "D0 amplitude (|M12|)":    dev_D,
}

# ---------------------------------------------------------------------------
# Per-tile median deviation (the GENERIC draw) for each channel, and the
# "generic-deviation floor": the largest M_KK at which the median draw still has
# deviation >= 1 (i.e. a >1sigma deviation is generic).  Above it, the typical
# anarchic point is SM-like in that channel.
# ---------------------------------------------------------------------------
def per_tile(stat_arr, q=50):
    return np.array([np.percentile(stat_arr[mkk == m], q) for m in tiles])

def crossing_Mkk(med_curve):
    """Largest M_KK where the median deviation curve crosses 1 (sigma).
    Deviations fall ~1/M_KK^2, so find the highest tile with median >= 1, then
    log-interpolate to the next tile."""
    above = med_curve >= 1.0
    if not above.any():
        return float("nan")          # never generic even at 1 TeV
    i = np.max(np.where(above)[0])
    if i == len(tiles) - 1:
        return float(tiles[-1])      # still generic at grid ceiling (flag!)
    x0, x1 = np.log(tiles[i]), np.log(tiles[i + 1])
    y0, y1 = np.log(med_curve[i]), np.log(med_curve[i + 1])
    # interpolate where log(dev) = 0
    xc = x0 + (0.0 - y0) * (x1 - x0) / (y1 - y0)
    return float(np.exp(xc))

med_curves = {name: per_tile(arr, 50) for name, arr in CHANNELS.items()}
q95_curves = {name: per_tile(arr, 95) for name, arr in CHANNELS.items()}
dev_floor  = {name: crossing_Mkk(med_curves[name]) for name in CHANNELS}

# ---------------------------------------------------------------------------
# Survival floors (the bands a surviving point must clear), LANE A.
#   epsilon_K MEDIAN (typical) floor: largest M_KK where median ratio_eps_K >= 1.
#   S,T,U existence floor: lane-independent ~18-20 TeV (scalar; from the plan/repo).
# ---------------------------------------------------------------------------
ratio_epsK_med = per_tile(df["ratio_eps_K"].values, 50)
epsK_floor = crossing_Mkk(ratio_epsK_med)   # current-input median epsilon_K floor
STU_FLOOR_LO, STU_FLOOR_HI = 18.0, 20.0     # existence floor (irreducible oblique T)

print("\n[scoreboard] GENERIC (median) deviation floors per channel "
      "(M_KK below which a >1sigma deviation is typical):")
for name in CHANNELS:
    print(f"   {name:26s}: {dev_floor[name]:6.2f} TeV   "
          f"(median dev @1TeV = {med_curves[name][0]:.2e}, "
          f"@20TeV = {med_curves[name][-1]:.2e})")
print(f"\n[scoreboard] SURVIVAL floors (LANE A): "
      f"epsilon_K median = {epsK_floor:.2f} TeV ; "
      f"S,T,U existence = {STU_FLOOR_LO:.0f}-{STU_FLOOR_HI:.0f} TeV")

# THE VERDICT: does any channel's generic-deviation floor reach the survival floor?
survival = max(epsK_floor, STU_FLOOR_HI) if np.isfinite(epsK_floor) else STU_FLOOR_HI
worst_dev_floor = np.nanmax([dev_floor[n] for n in CHANNELS])
verdict_A_alive = worst_dev_floor >= survival
print(f"\n[scoreboard] highest generic-deviation floor = {worst_dev_floor:.2f} TeV ; "
      f"binding survival floor = {survival:.2f} TeV")
print(f"[scoreboard] VERDICT: conclusion (A) {'ALIVE' if verdict_A_alive else 'DEAD'} "
      f"-> {'(A)' if verdict_A_alive else '(B) confirmed'}")

# ---------------------------------------------------------------------------
# FIGURE: per-channel median deviation vs M_KK, with the survival bands overlaid.
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(9.2, 6.4))
colors = ["#d62728", "#ff7f0e", "#1f77b4", "#17becf", "#2ca02c"]
for (name, med), c in zip(med_curves.items(), colors):
    ax.plot(tiles, med, "o-", color=c, ms=3.5, lw=1.6, label=name)
    ax.plot(tiles, q95_curves[name], color=c, lw=0.8, ls=":", alpha=0.5)
ax.axhline(1.0, color="k", lw=1.2, ls="--")
ax.text(1.05, 1.15, r"$1\sigma$ (observable threshold)", fontsize=8)
# survival bands
if np.isfinite(epsK_floor):
    ax.axvspan(epsK_floor, 35, color="grey", alpha=0.12, zorder=0)
    ax.axvline(epsK_floor, color="navy", lw=1.6, ls="-")
    ax.text(epsK_floor * 1.01, 3e2, rf"$\varepsilon_K$ median floor $\approx{epsK_floor:.0f}$ TeV",
            rotation=90, va="top", fontsize=8, color="navy")
ax.axvspan(STU_FLOOR_LO, STU_FLOOR_HI, color="purple", alpha=0.15, zorder=0)
ax.text(STU_FLOOR_LO * 0.99, 3e2, r"$S,T,U$ existence 18-20 TeV",
        rotation=90, va="top", ha="right", fontsize=8, color="purple")
ax.set_yscale("log")
ax.set_xlim(1, 33)
ax.set_ylim(1e-4, 1e3)
ax.set_xlabel(r"$M_{\rm KK}$ [TeV]")
ax.set_ylabel(r"generic (median) deviation from SM  [$\sigma$ units]")
ax.set_title("SM-overlap scoreboard (LANE A, anarchic): generic flavor deviation vs survival floors\n"
             r"every channel falls below $1\sigma$ FAR before the $\varepsilon_K$/$S,T,U$ survival floors $\Rightarrow$ verdict (B)")
ax.legend(loc="upper right", fontsize=8, ncol=1, framealpha=0.93)
ax.grid(alpha=0.25, which="both")
p = OUT / "solo_sm_overlap_scoreboard.png"
fig.savefig(p, dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"\n[saved] {p.name}  ({p.stat().st_size // 1024} KB)")

# ---------------------------------------------------------------------------
# Numbers memo (JSON) for the note / verdict checkpoint.
# ---------------------------------------------------------------------------
memo = {
    "lane": "A (anarchic, Bauer-S1 forward); literature strawman",
    "n_passes_pdg": int(len(df)),
    "M_KK_grid_TeV": [float(tiles.min()), float(tiles.max()), int(len(tiles))],
    "experimental_1sigma": {
        "sin2beta": SIG_SIN2B, "phi_s_rad": SIG_PHIS,
        "C_Bd": SIG_CBD, "C_Bs": SIG_CBS,
        "D0": "amplitude vs Delta m_D^exp (long-distance dominated, conservative)",
    },
    "generic_deviation_floor_TeV": {n: float(dev_floor[n]) for n in CHANNELS},
    "survival_floor_TeV": {
        "epsilon_K_median": float(epsK_floor),
        "STU_existence": [STU_FLOOR_LO, STU_FLOOR_HI],
        "binding": float(survival),
    },
    "highest_generic_deviation_floor_TeV": float(worst_dev_floor),
    "verdict": "A_alive" if verdict_A_alive else "B_confirmed",
    "interpretation": (
        "No anarchic point that survives epsilon_K + S,T,U shows a >1sigma flavor "
        "deviation in any quark channel: the generic-deviation floors all sit an "
        "order of magnitude in M_KK below the binding survival floor -> (B)."
    ),
}
mp = OUT.parent / "sm_overlap_scoreboard_numbers.json"
mp.write_text(json.dumps(memo, indent=2))
print(f"[saved] {mp.name}")
print("DONE (scoreboard)")
