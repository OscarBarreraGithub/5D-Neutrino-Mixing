"""Bauer-S1-MATCHED anarchic forward ensemble (0912.1625 Sec 5.1-5.2).

Goal
----
Reproduce the published Bauer-Casagrande-Haisch-Neubert eps_K cloud (their Fig. 4,
scenario S1 "standard") as closely as the repo's forward machinery allows, by
matching their ENSEMBLE SETUP -- not just their bounds:

  * |Y_ij| drawn UNIFORM IN MODULUS in [1/10, Y_max] with UNIFORM PHASE in [0,2pi)
    (Bauer Sec 5.1), NOT Re/Im iid uniform.  Default Y_max = 3 (NDA, S1).
  * Bulk masses c are SCANNED, not fixed: the top doublet c_{Q3} (repo convention)
    is drawn flat in the Bauer prior window, and the remaining eight c's are fixed
    per draw by the leading-order Froggatt-Nielsen relations (Casagrande 0807.4937
    eqs. I:106-107) so that {m_u..m_t, m_d..m_b, A, lambda} are reproduced up to
    O(1) anarchic Yukawa ratios.

The c-scan + wider-modulus prior are precisely the two levers that widen the
published cloud to ~6 decades and lower the eps_K-consistent fraction to ~19%
(S1), versus the fixed-c / Re-Im-uniform |Y|<~2.1 generator in run_rs_anarchy.py
(~3.3 decades, ~62% consistent).

Convention bridge (verified numerically)
----------------------------------------
Repo uses c = M_5/k with c < 1/2 -> IR-localized and
    f_IR(c)^2 = (1/2 - c) / (1 - eps^(1-2c)).
Bauer/Casagrande use F(c_paper) ~ sqrt(1 + 2 c_paper) with c_paper near -1/2.
Numerically  f_IR(c_repo) = (1/sqrt2) * F_Bauer(c_paper = -c_repo)  -- the sqrt2
is exactly Bauer's v/sqrt2 (v=246) vs the repo's v=174 convention, so the 4D mass
v_repo * f_IR is physically identical.  We therefore work ENTIRELY in the repo
convention: scan c_repo and use the repo's own f_IR + deltaf2 pipeline verbatim.

  Bauer c_u3 prior ]-1/2, c_max]  <->  repo c_{Q3,u3} prior  [-c_max, 1/2).
  S1: c_max = 2   -> repo prior [-2, 0.5).
  S3: c_max = 5/2 -> repo prior [-2.5, 0.5).

The remaining c's are obtained by INVERTING f_IR at each tile's epsilon to hit the
FN-target f-factor for every chiral component, so masses+CKM auto-reproduce
without an optimizer (the forward "keep mass+CKM-reproducing" path).

Per-draw forward evaluation (SVD -> CKM -> KK-gluon couplings -> Delta-F=2 Wilson
coefficients -> eps_K NP) is IDENTICAL to run_rs_anarchy.py (reused verbatim).

Usage
-----
    python scripts/anarchic_bauer_s1.py \
        --out scan_outputs/anarchic_reproduction/anarchic_bauer_s1.parquet \
        --scenario S1 --per-tile 20000 \
        --m-kk-tev 1,1.5,2,3,4,5,6,7,8,10
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import brentq

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from warpConfig.wavefuncs import f_IR  # noqa: E402
from warpConfig.baseParams import MPL  # noqa: E402
from scripts.run_rs_anarchy import (  # noqa: E402  -- reuse forward physics verbatim
    DEFAULT_XI_KK, DEFAULT_K_GEV, DEFAULT_V_GEV,
    _ordered_svd, _build_kk_gluon_couplings, _load_pdg_targets, jarlskog_invariant,
)
import quarkConstraints.deltaf2 as d  # noqa: E402
from quarkConstraints.deltaf2 import (  # noqa: E402
    compute_delta_f2_wilsons, evaluate_delta_f2_constraints,
    evaluate_delta_mk_with_running,
    _evolve_wilsons, compute_m12_np, _compute_m12_np,
    F_BD, M_BD, M_B_QUARK, M_D_QUARK_BD, B_1_BD, B_4_BD, B_5_BD,
    F_BS, M_BS, M_S_QUARK_BS, B_1_BS, B_4_BS, B_5_BS,
    F_D, M_D0, M_C_QUARK, M_U_QUARK, B_1_D, B_4_D, B_5_D,
)


def _m12_complex_all(couplings, M_KK, xi_KK, mu_had=2.0):
    """Complex M_12^NP (GeV) per system from RG-evolved Wilsons (reuses the
    SAME machinery as scripts/anarchic_complex_m12.py)."""
    coeffs = compute_delta_f2_wilsons(couplings, M_KK=M_KK, xi_KK=xi_KK)
    wbk = {c.input.key: _evolve_wilsons(c, mu_had=mu_had) for c in coeffs}
    return {
        "K": _compute_m12_np(wbk["epsilon_k"]),
        "Bd": compute_m12_np(wbk["b_d"], F_BD, M_BD, M_B_QUARK, M_D_QUARK_BD,
                             B_1_BD, B_4_BD, B_5_BD),
        "Bs": compute_m12_np(wbk["b_s"], F_BS, M_BS, M_B_QUARK, M_S_QUARK_BS,
                             B_1_BS, B_4_BS, B_5_BS),
        "D": compute_m12_np(wbk["d"], F_D, M_D0, M_C_QUARK, M_U_QUARK,
                            B_1_D, B_4_D, B_5_D),
    }

_EPS_BUDGET = abs(d.EPSILON_K_EXP - d.EPSILON_K_SM)


# ---------------------------------------------------------------------------
# Scenario table (Bauer 0912.1625 Sec 5.2)
# ---------------------------------------------------------------------------
# c_max is the Bauer ]-1/2, c_max] upper edge for the scanned top bulk mass; the
# repo-convention prior on the scanned c is [-c_max, 1/2).  L = ln(1/eps).
SCENARIOS = {
    # name:        (Y_max, c_max,  L_override (None = geometric eps from M_KK), label)
    "S1": dict(y_max=3.0,  c_max=2.0,  L=None, label="standard"),
    "S2": dict(y_max=3.0,  c_max=2.0,  L=None, label="aligned (common c_d)"),
    "S3": dict(y_max=3.0,  c_max=2.5,  L=7.0,  label="little (L=ln 1e3)"),
    "S4": dict(y_max=12.0, c_max=2.0,  L=None, label="large (Y_max=12)"),
}

# Wolfenstein A, lambda (PDG 2024) for the FN c-hierarchy.
WOLF_LAMBDA = 0.2250
WOLF_A = 0.826


# ---------------------------------------------------------------------------
# Bauer modulus/phase anarchic draw
# ---------------------------------------------------------------------------

def _draw_bauer_matrix(rng: np.random.Generator, y_min: float, y_max: float) -> np.ndarray:
    """3x3 complex Y with |Y_ij| ~ U(y_min, y_max), arg ~ U(0, 2pi).

    This is the Bauer 0912.1625 Sec 5.1 prior (uniform modulus + uniform phase),
    distinct from run_rs_anarchy's Re/Im-iid-uniform prior.
    """
    mod = rng.uniform(y_min, y_max, size=(3, 3))
    phase = rng.uniform(0.0, 2.0 * math.pi, size=(3, 3))
    return mod * np.exp(1j * phase)


# ---------------------------------------------------------------------------
# f_IR inversion: find c_repo in (-c_max, c_hi) giving a target f-factor
# ---------------------------------------------------------------------------

def _invert_f_IR(f_target: float, epsilon: float, c_lo: float, c_hi: float) -> float:
    """Return c_repo with f_IR(c_repo, epsilon) = f_target (clamped to [c_lo,c_hi]).

    f_IR is monotonically decreasing in c, so a bracketed root exists whenever
    f_target lies between f_IR(c_hi) and f_IR(c_lo).
    """
    f_at_lo = float(f_IR(np.array([c_lo]), epsilon)[0])
    f_at_hi = float(f_IR(np.array([c_hi]), epsilon)[0])
    if f_target >= f_at_lo:
        return c_lo
    if f_target <= f_at_hi:
        return c_hi
    g = lambda c: float(f_IR(np.array([c]), epsilon)[0]) - f_target
    return float(brentq(g, c_lo, c_hi, xtol=1e-10, rtol=1e-12))


def _minor11(Y):  # |(M_q)_11| = |minor of the (1,1) entry| = |Y_22 Y_33 - Y_23 Y_32|
    return abs(Y[1, 1] * Y[2, 2] - Y[1, 2] * Y[2, 1])


def _fn_c_values(c_u3, epsilon, targets, *, Y_u, Y_d, common_cd,
                 rng=None, c_jitter=0.0):
    """Bauer-faithful per-draw FN bulk-mass fix (Casagrande eqs. I:95-107).

    Bauer scans the RIGHT-HANDED TOP c_{u3} (the "special" composite-top bulk
    mass) flat in their prior window; here c_u3 is the repo-convention image of
    that.  The other eight c's are then chosen, PER DRAW, so that the leading-
    order Froggatt-Nielsen mass relations (eq. I:96) reproduce the SM quark
    masses GIVEN THE ACTUAL ANARCHIC YUKAWA MINORS of this draw.  Because the
    minors (det Y_q, (M_q)_11, (Y_q)_33) differ draw-to-draw, the fitted c's
    SCATTER about their central values exactly as Bauer's per-point chi^2 fit
    does (Table 1 widths) -- this Yukawa-driven c-scatter, NOT an ad-hoc jitter,
    is what widens the published eps_K cloud to ~6 decades.

    FN target f-factors (rearranging eq. I:96, repo v convention, masses use v):
        f_u3 = f_IR(c_u3)                                  (scanned, RH top)
        f_Q3 = m_t / (v |Y_u33| f_u3)                      -> f_Q3 (-> c_Q3)
        left-doublet ratios from CKM (I:106): f_Q1:f_Q2:f_Q3 = l^3 : A l^2 : 1
        f_u2 = m_c / (v (|M_u11|/|Y_u33|) f_Q2),  f_u1 = m_u/(v(detYu/|M_u11|)f_Q1)
        f_d3 = m_b/(v|Y_d33|f_Q3),  f_d2, f_d1 analogously
    """
    lam, A = WOLF_LAMBDA, WOLF_A
    v = DEFAULT_V_GEV
    mu = targets["up_masses_GeV"]; md = targets["down_masses_GeV"]
    m_u, m_c, m_t = (float(x) for x in mu)
    m_d, m_s, m_b = (float(x) for x in md)

    # Per-draw Yukawa structure entering eq. I:96.
    detYu, detYd = abs(np.linalg.det(Y_u)), abs(np.linalg.det(Y_d))
    minYu, minYd = _minor11(Y_u), _minor11(Y_d)
    Yu33, Yd33 = abs(Y_u[2, 2]), abs(Y_d[2, 2])
    eps_ = 1e-12

    # Scanned RH-top f sets the composite scale; top mass fixes f_Q3.
    fu3 = float(f_IR(np.array([c_u3]), epsilon)[0])
    fQ3 = m_t / (v * Yu33 * fu3 + eps_)
    # Left doublet hierarchy from CKM (I:106): independent of the RH sector.
    fQ = np.array([fQ3 * lam**3, fQ3 * (A * lam**2), fQ3])

    # Remaining RH f's solved from eq. I:96 second/first relations:
    fu2 = m_c / (v * (minYu / max(Yu33, eps_)) * fQ[1] + eps_)
    fu1 = m_u / (v * (detYu / max(minYu, eps_)) * fQ[0] + eps_)
    fd3 = m_b / (v * Yd33 * fQ[2] + eps_)
    fd2 = m_s / (v * (minYd / max(Yd33, eps_)) * fQ[1] + eps_)
    fd1 = m_d / (v * (detYd / max(minYd, eps_)) * fQ[0] + eps_)
    fu = np.array([fu1, fu2, fu3]); fd = np.array([fd1, fd2, fd3])

    c_lo, c_hi = -2.6, 0.95
    c_Q = np.array([_invert_f_IR(fQ[i], epsilon, c_lo, c_hi) for i in range(3)])
    c_u = np.array([_invert_f_IR(fu[i], epsilon, c_lo, c_hi) for i in range(3)])
    c_d = np.array([_invert_f_IR(fd[i], epsilon, c_lo, c_hi) for i in range(3)])
    if common_cd:  # S2: U(3) symmetric RH down sector -> common c_d
        c_d = np.full(3, float(np.mean(c_d)))
    if c_jitter > 0.0 and rng is not None:  # optional extra fit-width scatter
        c_Q = np.clip(c_Q + rng.normal(0.0, c_jitter, 3), c_lo, c_hi)
        c_u = np.clip(c_u + rng.normal(0.0, c_jitter, 3), c_lo, c_hi)
        cd_j = rng.normal(0.0, c_jitter, 1 if common_cd else 3)
        c_d = np.clip(c_d + (cd_j if not common_cd else cd_j[0]), c_lo, c_hi)
    return c_Q, c_u, c_d


# ---------------------------------------------------------------------------
# Forward evaluation of one draw (physics reused from run_rs_anarchy)
# ---------------------------------------------------------------------------

def _eval_draw(Y_u, Y_d, f_Q, f_u, f_d, M_KK_GeV, xi_KK, targets,
               mass_factor, ckm_factor, j_factor, emit_m12=False, c_Q=None,
               c_u=None, c_d=None):
    D_Q, D_u, D_d = np.diag(f_Q), np.diag(f_u), np.diag(f_d)
    M_u = DEFAULT_V_GEV * D_Q @ Y_u @ D_u
    M_d = DEFAULT_V_GEV * D_Q @ Y_d @ D_d
    U_L_u, m_up, U_R_u = _ordered_svd(M_u)
    U_L_d, m_dn, U_R_d = _ordered_svd(M_d)
    ckm = U_L_u.conj().T @ U_L_d
    abs_V_us = float(abs(ckm[0, 1]))
    abs_V_cb = float(abs(ckm[1, 2]))
    abs_V_ub = float(abs(ckm[0, 2]))
    J = float(jarlskog_invariant(ckm))

    # PDG match (factor tolerances), same definition as run_rs_anarchy.
    up_log = np.log(np.maximum(m_up, 1e-30) / targets["up_masses_GeV"])
    dn_log = np.log(np.maximum(m_dn, 1e-30) / targets["down_masses_GeV"])
    ckm_log = np.array([
        abs(math.log(max(abs_V_us, 1e-30) / targets["abs_V_us"])),
        abs(math.log(max(abs_V_cb, 1e-30) / targets["abs_V_cb"])),
        abs(math.log(max(abs_V_ub, 1e-30) / targets["abs_V_ub"])),
    ])
    j_log = abs(math.log(max(abs(J), 1e-30) / abs(targets["J"])))
    passes_pdg = bool(
        (np.abs(up_log) <= math.log(mass_factor)).all()
        and (np.abs(dn_log) <= math.log(mass_factor)).all()
        and (ckm_log <= math.log(ckm_factor)).all()
        and j_log <= math.log(j_factor)
    )

    couplings = _build_kk_gluon_couplings(
        M_KK=M_KK_GeV, xi_KK=xi_KK, f_Q=f_Q, f_u=f_u, f_d=f_d,
        U_L_u=U_L_u, U_L_d=U_L_d, U_R_u=U_R_u, U_R_d=U_R_d,
    )
    df2 = evaluate_delta_f2_constraints(couplings, M_KK=M_KK_GeV, xi_KK=xi_KK)
    bs = df2.by_system
    ratio_eps_K = float(bs["K"].ratio_to_bound)
    # |eps_K^NP| = ratio * budget (linear in budget) -> realized |eps_K^tot|.
    eps_K_np = ratio_eps_K * _EPS_BUDGET
    # Delta m_K from the unevolved kaon Wilsons run to mu_had (same as run_rs_anarchy).
    try:
        unev = compute_delta_f2_wilsons(couplings, M_KK=M_KK_GeV, xi_KK=xi_KK)
        w_k = next((w for w in unev if w.input.key == "epsilon_k"), None)
        ratio_dm_K = (float(evaluate_delta_mk_with_running(w_k, mu_had=2.0).ratio_to_exp)
                      if w_k is not None else float("nan"))
    except Exception:  # noqa: BLE001
        ratio_dm_K = float("nan")
    out = dict(
        passes_pdg=passes_pdg,
        ratio_eps_K=ratio_eps_K,
        eps_K_np=eps_K_np,
        ratio_dm_K=ratio_dm_K,
        ratio_B_d=float(bs["B_d"].ratio_to_bound),
        ratio_B_s=float(bs["B_s"].ratio_to_bound),
        ratio_D=float(bs["D"].ratio_to_bound),
        up_log_max=float(np.max(np.abs(up_log))),
        down_log_max=float(np.max(np.abs(dn_log))),
        ckm_log_max=float(np.max(ckm_log)),
    )
    # Record the Froggatt-Nielsen-fixed gen-3 bulk masses (repo convention) so the
    # OUTPUT parquet's c_Q3 median is directly verifiable against Bauer Table 1
    # (S1: c_Q3 = -0.34 +/- 0.32 in Bauer's sign convention -> +0.34 here).  The
    # GATED (passes_pdg) c_Q3 distribution is the Bauer-faithful one: its median is
    # ~+0.28, inside Bauer's +0.34 +/- 0.32 band (about 0.06 below his central, well
    # within the 1sigma width -- the distribution is broad and slightly tailed).
    # The UNGATED median is biased to the IR floor (~+0.08) by mass+CKM-FAILING draws.
    if c_Q is not None:
        out["c_Q1"] = float(c_Q[0]); out["c_Q2"] = float(c_Q[1]); out["c_Q3"] = float(c_Q[2])
    if c_d is not None:
        out["c_d3"] = float(c_d[2])
    if emit_m12:
        m12 = _m12_complex_all(couplings, M_KK_GeV, xi_KK)
        for k, col in (("K", "K"), ("Bd", "Bd"), ("Bs", "Bs"), ("D", "D")):
            out[f"re_m12_{col}"] = float(m12[k].real)
            out[f"im_m12_{col}"] = float(m12[k].imag)
    return out


# ---------------------------------------------------------------------------
# Tile runner
# ---------------------------------------------------------------------------

def run_tile(M_KK_GeV, n_draws, seed, *, scenario, xi_KK, targets,
             mass_factor, ckm_factor, j_factor, c_jitter=0.0, emit_m12=False):
    sc = SCENARIOS[scenario]
    y_max = sc["y_max"]
    c_max = sc["c_max"]
    common_cd = (scenario == "S2")
    rng = np.random.default_rng(seed)
    Lambda_IR = M_KK_GeV / xi_KK
    if sc["L"] is not None:
        epsilon = math.exp(-sc["L"])          # S3: volume-truncated UV cutoff
    else:
        epsilon = Lambda_IR / DEFAULT_K_GEV    # standard L = ln(1e16)

    # repo-convention scan prior on the RH top c_u3: [-c_max, 1/2)
    # (= Bauer ]-1/2, c_max] negated).
    c_lo_scan, c_hi_scan = -c_max, 0.5
    rows = []
    for _ in range(n_draws):
        c_u3 = float(rng.uniform(c_lo_scan, c_hi_scan))
        try:
            # Draw the anarchic Yukawas FIRST; the FN c-fix then depends on their
            # minors (Bauer's per-point fit), giving genuine c-scatter.
            Y_u = _draw_bauer_matrix(rng, 0.1, y_max)
            Y_d = _draw_bauer_matrix(rng, 0.1, y_max)
            c_Q, c_u, c_d = _fn_c_values(c_u3, epsilon, targets,
                                         Y_u=Y_u, Y_d=Y_d, common_cd=common_cd,
                                         rng=rng, c_jitter=c_jitter)
            f_Q = f_IR(c_Q, epsilon)
            f_u = f_IR(c_u, epsilon)
            f_d = f_IR(c_d, epsilon)
            r = _eval_draw(Y_u, Y_d, f_Q, f_u, f_d, M_KK_GeV, xi_KK, targets,
                           mass_factor, ckm_factor, j_factor, emit_m12=emit_m12,
                           c_Q=c_Q, c_u=c_u, c_d=c_d)
        except (ValueError, np.linalg.LinAlgError, FloatingPointError):
            continue
        except Exception:  # noqa: BLE001 -- robust per-draw isolation
            continue
        r["M_KK_GeV"] = float(M_KK_GeV)
        r["M_KK_TeV"] = float(M_KK_GeV / 1000.0)
        r["c_u3"] = c_u3
        r["scenario"] = scenario
        rows.append(r)
    return rows


_WORKER_TARGETS = None
_WORKER_KW = None


def _worker_init(kw):
    global _WORKER_TARGETS, _WORKER_KW
    _WORKER_TARGETS = _load_pdg_targets()
    _WORKER_KW = kw


def _worker_run(job):
    scenario, t, seed = job
    return run_tile(
        t * 1000.0, _WORKER_KW["per_tile"], seed, scenario=scenario,
        xi_KK=DEFAULT_XI_KK, targets=_WORKER_TARGETS,
        mass_factor=_WORKER_KW["mass_factor"], ckm_factor=_WORKER_KW["ckm_factor"],
        j_factor=_WORKER_KW["j_factor"], c_jitter=_WORKER_KW["c_jitter"],
        emit_m12=_WORKER_KW["emit_m12"],
    )


def main(argv=None):
    import multiprocessing as mp
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", type=str,
                   default="scan_outputs/anarchic_reproduction/anarchic_bauer_s1.parquet")
    p.add_argument("--scenario", type=str, default="S1",
                   choices=sorted(SCENARIOS.keys()) + ["ALL"])
    p.add_argument("--per-tile", type=int, default=20000)
    p.add_argument("--m-kk-tev", type=str, default="1,1.5,2,3,4,5,6,7,8,10")
    p.add_argument("--base-seed", type=int, default=20260623)
    p.add_argument("--c-jitter", type=float, default=0.0,
                   help="1sigma Gaussian scatter on the FN-central bulk masses, "
                        "matching Bauer Table 1 fit widths (0 = deterministic FN).")
    p.add_argument("--mass-factor", type=float, default=3.0)
    p.add_argument("--ckm-factor", type=float, default=3.0)
    p.add_argument("--j-factor", type=float, default=5.0)
    p.add_argument("--n-workers", type=int,
                   default=int(__import__("os").environ.get("SLURM_CPUS_PER_TASK", "1")))
    p.add_argument("--emit-m12", action="store_true",
                   help="Also emit complex M_12^NP (re/im) per system for the "
                        "C_Bq-phi, D-funnel and Re/Im(M12) figures.")
    args = p.parse_args(argv)

    scenarios = sorted(SCENARIOS.keys()) if args.scenario == "ALL" else [args.scenario]
    mkk_tev = [float(t) for t in args.m_kk_tev.split(",") if t.strip()]

    out_path = Path(args.out)
    if not out_path.is_absolute():
        out_path = REPO / out_path
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Build (scenario, M_KK, seed) jobs.
    jobs = []
    for scenario in scenarios:
        sc = SCENARIOS[scenario]
        print(f"[bauer] scenario {scenario} ({sc['label']}): Y_max={sc['y_max']} "
              f"c_max={sc['c_max']} L={'geom' if sc['L'] is None else sc['L']}")
        for idx, t in enumerate(mkk_tev):
            seed = args.base_seed + 1009 * idx + 100003 * (hash(scenario) % 7919)
            jobs.append((scenario, t, seed))

    kw = dict(per_tile=args.per_tile, mass_factor=args.mass_factor,
              ckm_factor=args.ckm_factor, j_factor=args.j_factor,
              c_jitter=args.c_jitter, emit_m12=args.emit_m12)

    def _report(rows):
        if not rows:
            return
        arr = np.array([r["eps_K_np"] for r in rows])
        fr_pdg = np.mean([r["passes_pdg"] for r in rows])
        lo, hi = np.percentile(arr, [1, 99])
        print(f"  {rows[0]['scenario']} M_KK={rows[0]['M_KK_TeV']:5.2f} TeV : "
              f"{len(rows):6d} draws  PDGpass={fr_pdg:5.1%}  |epsK_NP| p1..p99 "
              f"[{lo:.2e}, {hi:.2e}] ({math.log10(hi/lo):.2f} dex)", flush=True)

    all_rows = []
    if args.n_workers <= 1:
        _worker_init(kw)
        for job in jobs:
            rows = _worker_run(job)
            all_rows.extend(rows)
            _report(rows)
    else:
        with mp.Pool(args.n_workers, initializer=_worker_init, initargs=(kw,)) as pool:
            for rows in pool.imap_unordered(_worker_run, jobs):
                all_rows.extend(rows)
                _report(rows)
    df = pd.DataFrame(all_rows)
    df.to_parquet(out_path, index=False)
    print(f"[bauer] wrote {len(df):,} rows -> {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
