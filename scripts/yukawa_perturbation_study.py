#!/usr/bin/env python
"""Yukawa perturbation study for RS eps_K (RS-FLAVOR-ALIGNMENT-2026-07).

Treats the eps_K constraint as a function on Yukawa space and probes its local
structure around surviving points of four classes:
  - flat_typical : flat-anarchic survivor, magnitude-suppressed (small |C4|)
  - flat_tuned   : flat-anarchic survivor, phase-aligned (large |C4|, theta_K~0)
  - nelson_barr  : real Y_d (CP in up sector)
  - u2           : rank-one/U(2) light-family-symmetric ensemble

Three analyses:
  A) NAIVE multiplicative noise: Y_ij -> Y_ij*(1+sigma*z_ij), z~N(0,1) per entry;
     scan sigma, record eps_K pass-fraction AND median |eps_K|/bound (how badly it
     fails).
  B) LINEAR RESPONSE: finite-difference gradient g = d(Im C_4)/d(Y_d components),
     the local linear transformation (a covector on the 18-real-dim Y_d tangent
     space). Gives the per-direction tuning radius delta_crit(Zhat)=S_bound/|g.Zhat|.
  C) GRADIENT FIELD: collect g over an ensemble per class, PCA the unit gradient
     directions -> effective dimension of the "dangerous subspace" = codimension of
     the protected flat. Anarchy fills many dims; U(2) confines to a few.

Im(C_4)_12 is recovered as  S = C4abs_12 * sin(Phi_12)  (signed, smooth) from the
instrumented evaluator, so gradients are well defined even at theta_K~0.
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.anarchic_bauer_s1 import DEFAULT_K_GEV, SCENARIOS, _draw_bauer_matrix, _fn_c_values
from scripts.instrument_epsK_phase import _instrument_draw, draw_nelson_barr_yukawas
from scripts.rankone_u2_lane import _assign_u2_bulk_masses, _draw_rankone_u2_yukawas
from scripts.reproducible_seeds import stable_seed_offset
from scripts.run_rs_anarchy import _load_pdg_targets
from warpConfig.wavefuncs import f_IR

TARGETS = _load_pdg_targets()
MKK = 3000.0
XI = 1.0
EPSG = MKK / XI / DEFAULT_K_GEV
MF, CF, JF = 3.0, 3.0, 10.0
YMAX = SCENARIOS["S1"]["y_max"]


def evaluate(Yu, Yd, fQ, fu, fd):
    """Return (ratio_eps_K, ImC4, passes_pdg). ImC4 = C4abs_12*sin(Phi_12) (signed)."""
    r = _instrument_draw(Yu, Yd, fQ, fu, fd, MKK, XI, TARGETS, MF, CF, JF)
    imc4 = r["C4abs_12"] * math.sin(r["Phi_12"])
    return r["ratio_eps_K"], imc4, r["passes_pdg"]


def _flat_fit(Yu, Yd, rng):
    cQ, cu, cd = _fn_c_values(float(rng.uniform(-2.0, 0.5)), EPSG, TARGETS,
                              Y_u=Yu, Y_d=Yd, common_cd=False, rng=rng, c_jitter=0.0)
    return f_IR(cQ, EPSG), f_IR(cu, EPSG), f_IR(cd, EPSG)


def get_base_point(kind, rng, budget=40000):
    """Return (Yu, Yd, fQ, fu, fd) for a surviving point of the requested class.
    flat_tuned uses best-of-budget (largest |C4| survivor = most phase-aligned)."""
    best = None; best_c4 = -1.0
    for _ in range(budget):
        if kind == "nelson_barr":
            Yu, Yd = draw_nelson_barr_yukawas(rng, y_max=YMAX, rho_cp=1.0, eta_leak=0.0)
            try: fQ, fu, fd = _flat_fit(Yu, Yd, rng)
            except Exception: continue
        elif kind == "u2":
            d = _draw_rankone_u2_yukawas(rng, TARGETS, rh_down_leak=1.0, rh_up_leak=1.0,
                                         left_23_boost=4.0, left_13_boost=4.0,
                                         left_cp_boost=20.0, perturbativity_mode="reject")
            if d is None: continue
            try: b = _assign_u2_bulk_masses(d, EPSG, TARGETS, c_Q12=0.5, c_Q3=0.25)
            except Exception: continue
            Yu, Yd, fQ, fu, fd = d.Y_u, d.Y_d, b.f_Q, b.f_u, b.f_d
        else:
            Yu = _draw_bauer_matrix(rng, 0.1, YMAX); Yd = _draw_bauer_matrix(rng, 0.1, YMAX)
            try: fQ, fu, fd = _flat_fit(Yu, Yd, rng)
            except Exception: continue
        try:
            re, imc4, ok = evaluate(Yu, Yd, fQ, fu, fd)
        except Exception:
            continue
        if not (ok and re <= 1.0):
            continue
        c4 = _instrument_draw(Yu, Yd, fQ, fu, fd, MKK, XI, TARGETS, MF, CF, JF)["C4abs_12"]
        if kind == "flat_typical":
            if c4 < 3e-17:                      # first small-|C4| survivor
                return Yu, Yd, fQ, fu, fd
        elif kind == "flat_tuned":
            if c4 > best_c4:                    # keep the most phase-aligned so far
                best_c4 = c4; best = (Yu, Yd, fQ, fu, fd)
            if c4 > 5e-16:                      # good enough, stop early
                return Yu, Yd, fQ, fu, fd
        else:
            return Yu, Yd, fQ, fu, fd
    if best is not None:
        return best
    raise RuntimeError(f"no base point found for {kind}")


# ---------------------------------------------------------------- A: noise
def multiplicative_noise(Yu, Yd, fQ, fu, fd, sigmas, ndir, rng, perturb="down"):
    out = []
    for s in sigmas:
        passes, ratios = 0, []
        for _ in range(ndir):
            Zd = 1.0 + s * rng.normal(size=(3, 3))
            Yd2 = Yd * Zd
            Yu2 = Yu * (1.0 + s * rng.normal(size=(3, 3))) if perturb == "both" else Yu
            try:
                re, _, ok = evaluate(Yu2, Yd2, fQ, fu, fd)
            except Exception:
                continue
            ratios.append(re)
            if re <= 1.0:
                passes += 1
        n = max(1, len(ratios))
        out.append((s, passes / n, float(np.median(ratios)) if ratios else float("nan")))
    return out


# ---------------------------------------------------------------- B: gradient
def imc4_gradient(Yu, Yd, fQ, fu, fd, h=1e-4):
    """Finite-difference gradient of S=Im(C4)_12 wrt the 18 real components of Y_d.
    Returns (g[18], S0, labels[18]).  Components ordered (Re11,Im11,Re12,...)."""
    def S_of(Ydx):
        _, imc4, _ = evaluate(Yu, Ydx, fQ, fu, fd)
        return imc4
    S0 = S_of(Yd)
    g = np.zeros(18); labels = []
    k = 0
    for i in range(3):
        for j in range(3):
            for part in ("re", "im"):
                d = np.zeros((3, 3), complex)
                d[i, j] = (1.0 if part == "re" else 1j) * h
                g[k] = (S_of(Yd + d) - S_of(Yd - d)) / (2 * h)
                labels.append(f"{part}Y{i+1}{j+1}")
                k += 1
    return g, S0, labels


def s_bound(Yu, Yd, fQ, fu, fd):
    """Value of |S|=|Im C4| at which ratio_eps_K=1 (calibrate from current point)."""
    re, imc4, _ = evaluate(Yu, Yd, fQ, fu, fd)
    if re <= 0 or abs(imc4) == 0:
        return float("nan")
    return abs(imc4) / re  # since ratio_eps_K is linear in |Im C4|


# ---------------------------------------------------------------- C: gradient field
def draw_ensemble_point(kind, rng):
    """A PDG-passing point of the ensemble (NOT required to pass eps_K) -- fast."""
    for _ in range(4000):
        if kind == "nelson_barr":
            Yu, Yd = draw_nelson_barr_yukawas(rng, y_max=YMAX, rho_cp=1.0, eta_leak=0.0)
            try: fQ, fu, fd = _flat_fit(Yu, Yd, rng)
            except Exception: continue
        elif kind == "u2":
            d = _draw_rankone_u2_yukawas(rng, TARGETS, rh_down_leak=1.0, rh_up_leak=1.0,
                                         left_23_boost=4.0, left_13_boost=4.0,
                                         left_cp_boost=20.0, perturbativity_mode="reject")
            if d is None: continue
            try: b = _assign_u2_bulk_masses(d, EPSG, TARGETS, c_Q12=0.5, c_Q3=0.25)
            except Exception: continue
            Yu, Yd, fQ, fu, fd = d.Y_u, d.Y_d, b.f_Q, b.f_u, b.f_d
        else:
            Yu = _draw_bauer_matrix(rng, 0.1, YMAX); Yd = _draw_bauer_matrix(rng, 0.1, YMAX)
            try: fQ, fu, fd = _flat_fit(Yu, Yd, rng)
            except Exception: continue
        try:
            _, _, ok = evaluate(Yu, Yd, fQ, fu, fd)
        except Exception:
            continue
        if ok:
            return Yu, Yd, fQ, fu, fd
    return None


def gradient_field(kind, n, rng):
    """Collect unit gradient directions of Im(C4) over an ensemble; return (N,18)."""
    G = []
    for _ in range(n):
        p = draw_ensemble_point(kind, rng)
        if p is None: continue
        try:
            g, S0, _ = imc4_gradient(*p)
        except Exception:
            continue
        nrm = np.linalg.norm(g)
        if nrm > 0 and np.isfinite(nrm):
            G.append(g / nrm)
    return np.array(G)


def participation_dim(sv):
    """Effective dimension = (sum s^2)^2 / sum s^4 (participation ratio of variances)."""
    v = sv**2
    return float((v.sum()**2) / (v @ v)) if v.sum() > 0 else 0.0


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--part", default="A", choices=["A", "B", "C", "all"])
    ap.add_argument("--ndir", type=int, default=200)
    ap.add_argument("--nens", type=int, default=60)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--out", default=".orchestration/runs/RS-FLAVOR-ALIGNMENT-2026-07/perturb_study.npz")
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    classes = ["flat_typical", "flat_tuned", "nelson_barr", "u2"]
    save = {}

    if a.part in ("A", "B", "all"):
        print("finding base points...", flush=True)
        pts = {}
        for c in classes:
            class_seed = a.seed + stable_seed_offset(
                c,
                modulus=1000,
                namespace="yukawa_perturbation.base_point.v1",
            )
            pts[c] = get_base_point(c, np.random.default_rng(class_seed))
            re, imc4, ok = evaluate(*pts[c])
            print(f"  {c}: ratio_eps_K={re:.3f}  Im(C4)={imc4:.2e}", flush=True)

    if a.part in ("A", "all"):
        sigmas = np.array([3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1])
        print("\n=== A: multiplicative noise  Y_ij -> Y_ij*(1+sigma*N(0,1)) ===", flush=True)
        for c in classes:
            res = multiplicative_noise(*pts[c], sigmas, a.ndir, rng)
            save[f"noise_{c}"] = np.array(res)
            print(f"\n{c}:", flush=True)
            for s, pf, mr in res:
                print(f"  sigma={s:.1e}: pass={pf*100:5.1f}%  median ratio_eps_K={mr:.2e}", flush=True)

    if a.part in ("B", "all"):
        print("\n=== B: linear response (gradient of Im C4 wrt Y_d, tuning radius) ===", flush=True)
        for c in classes:
            g, S0, labels = imc4_gradient(*pts[c])
            Sb = s_bound(*pts[c])
            gnorm = np.linalg.norm(g)
            dmin = Sb / gnorm if gnorm > 0 else float("inf")  # tuning radius along steepest dir
            save[f"grad_{c}"] = g; save["gradlabels"] = np.array(labels)
            top = np.argsort(-np.abs(g))[:4]
            print(f"\n{c}: |grad|={gnorm:.2e}  S_bound={Sb:.2e}  tuning radius delta_min={dmin:.2e}", flush=True)
            print("   most dangerous entries: " + ", ".join(f"{labels[k]}({g[k]:+.1e})" for k in top), flush=True)

    if a.part in ("C", "all"):
        print("\n=== C: gradient-direction field, PCA -> dangerous-subspace dimension ===", flush=True)
        for c in ["flat_typical", "nelson_barr", "u2"]:
            class_seed = a.seed + 31 * stable_seed_offset(
                c,
                modulus=97,
                namespace="yukawa_perturbation.gradient_field.v1",
            )
            G = gradient_field(c, a.nens, np.random.default_rng(class_seed))
            if len(G) < 3:
                print(f"  {c}: too few gradients ({len(G)})"); continue
            sv = np.linalg.svd(G, compute_uv=False)
            save[f"pca_{c}"] = sv; save[f"nG_{c}"] = np.array([len(G)])
            print(f"  {c}: N={len(G)} gradients; PCA singular values (18) top: "
                  f"{np.round(sv[:6],3)}", flush=True)
            print(f"       effective dangerous-subspace dim (participation) = {participation_dim(sv):.1f}", flush=True)

    outp = Path(a.out)
    if not outp.is_absolute(): outp = REPO / outp
    np.savez(outp, **save)
    print(f"\nsaved -> {outp}", flush=True)
    print("DONE", flush=True)
