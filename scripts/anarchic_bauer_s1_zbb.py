"""Per-draw Z->bb consistency for the Bauer-S1 anarchic ensemble (0912.1625 Fig. 4).

Companion to ``scripts/anarchic_bauer_s1.py``.  That script produces the eps_K
"fat cloud" but stores NO Z->bb information (only the Delta-F=2 ratios + c_u3).
Bauer Fig. 4 colours each anarchic point by Z->bb consistency, so to reproduce
the figure we must compute, per draw, the TOTAL b-coupling shift and the
resulting (R_b, A_b) pseudo-observables, then flag whether the point is
consistent with the Z-pole Z->bb measurements at 99% CL (Bauer's cut).

What this script does
---------------------
It RE-DRAWS the EXACT same anarchic ensemble as ``anarchic_bauer_s1.py``
(identical Bauer modulus/phase Yukawa prior + per-draw Froggatt-Nielsen bulk-mass
fix, reusing ``_draw_bauer_matrix`` and ``_fn_c_values`` verbatim from that
module), but on each draw it ALSO:

  1. builds the minimal-RS EW couplings from the drawn (c_Q, c_u, c_d, F's, Y_d,
     SVD rotations) via ``quarkConstraints.rs_ew_couplings.build_rs_ew_couplings``
     with ``include_fermion_kk_mixing=True`` -- the SAME builder
     ``point_builder.build_from_rs_ew_inputs`` uses (it is the byte-identical
     code path; ``extract_plot_quantities.py`` / ``gen_sidebyside_points.py``
     read the result the same way);
  2. takes the TOTAL b-coupling shift as the mass-basis [2,2] entry of
     ``z_delta_g_L_d`` / ``z_delta_g_R_d`` = gauge-KK + fermion-KK admixture
     (exactly what the T010/T011 Z-pole adapter consumes,
     ``_coupling_entry(rsc, "z_delta_g_L_d", 2, 2)``);
  3. evaluates the Z-pole R_b / A_b pseudo-observables with the SM-fit Zbb
     reference T010/T011 use (bottom radiator calibrated to the LEP/SLC SM-fit
     R_b), shifted by (delta_g_L_b, delta_g_R_b);
  4. flags ``passes_Zbb`` = both R_b and A_b lie within 99% CL of their LEP/SLC
     measurements (|pred - exp| <= z99 * combined_sigma, z99 = 2.5758, one
     d.o.f. two-sided) -- Bauer's "consistent with Z->bb at 99% CL", which is a
     LOOSER cut than our default M1 NP-shift budget on purpose so the figure
     matches the paper.

It ALSO recomputes ``ratio_eps_K`` (and hence ``eps_K_np``) on the SAME draw via
the verbatim run_rs_anarchy Delta-F=2 path, so every output row is internally
self-consistent (same Yukawas -> same eps_K AND same Zbb flag).  The plot then
colours points with that one parquet.

Because the per-draw RS-EW spectrum overlap a(c) is ~1.4 s/draw for the wide
anarchic c-range, this runs a REPRESENTATIVE SUBSET per M_KK tile (default 2000
draws/tile), parallelized with multiprocessing.  A subset is sufficient for a
clean scatter; the eps_K cloud shape itself still comes from the full
``anarchic_bauer_s1.py`` ensemble in the plot (grey/orange), and the blue
Z->bb-consistent overlay uses this subset.

Approximation note
------------------
The Z->bb shift here is the FULL minimal-RS tree result (gauge-KK light-Z shift
PLUS the Casagrande fermion-KK admixture, ``include_fermion_kk_mixing=True``) --
NOT a gauge-only approximation.  It is byte-identical to the builder path the
fitted-point plots use.  The only difference from the fitted scan is that the
bulk c's are the anarchic FN-fixed draws (not PDG-refit values) and the SVD
rotations come from the drawn mass matrices, which is the whole point of an
anarchic forward scan.

Usage
-----
    python scripts/anarchic_bauer_s1_zbb.py \
        --out scan_outputs/anarchic_reproduction/anarchic_bauer_s1_zbb.parquet \
        --per-tile 2000 --m-kk-tev 1,1.5,2,3,4,5,6,7,8,10 --workers 8
"""
from __future__ import annotations

import argparse
import math
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from warpConfig.wavefuncs import f_IR  # noqa: E402
from scripts.anarchic_bauer_s1 import (  # noqa: E402  -- reuse the EXACT draw machinery
    BAUER_REPO_MASS_PREFAC_GEV, DEFAULT_K_GEV, DEFAULT_XI_KK, SCENARIOS,
    _deltaf2_helpers, _draw_bauer_matrix, _fn_c_values, _run_rs_anarchy_helpers,
)
# RS-EW build knobs -- MUST match scripts/run_full_catalog_scan.py defaults so
# the Z->bb shift is the same object the fitted-point plots use.
K_GEV = 1.2209e19
N_GAUGE_MODES = 512
QUADRATURE_ORDER = 4096
OVERLAP_RTOL = 1.0e-3
MIN_OVERLAP_MODES = 16
MAX_OVERLAP_MODES = 512
EW_MODEL = "minimal_rs"

# 99% CL two-sided critical value (1 d.o.f.) -- Bauer's Z->bb cut.
Z99 = 2.5758293035489004

# Mass / CKM / Jarlskog factor tolerances for the per-draw Bauer-style accept
# (identical to scripts/run_rs_anarchy.py + anarchic_bauer_s1.py scan defaults).
MASS_FACTOR = 3.0
CKM_FACTOR = 3.0
J_FACTOR = 5.0


# ---------------------------------------------------------------------------
# Duck-typed quark_fit_result for build_rs_ew_couplings (no PDG refit needed).
# build_rs_ew_couplings reads: bulk_state.{c_Q,c_u,c_d,F_Q,F_u,F_d,Y_d_bulk_basis},
# and quark_fit_result.{U_L_u,U_R_u,U_L_d,U_R_d,masses_down}.
# ---------------------------------------------------------------------------
class _DuckBulkState:
    __slots__ = ("c_Q", "c_u", "c_d", "F_Q", "F_u", "F_d", "Y_d_bulk_basis")

    def __init__(self, c_Q, c_u, c_d, F_Q, F_u, F_d, Y_d_bulk_basis):
        self.c_Q = np.asarray(c_Q, float)
        self.c_u = np.asarray(c_u, float)
        self.c_d = np.asarray(c_d, float)
        self.F_Q = np.asarray(F_Q, float)
        self.F_u = np.asarray(F_u, float)
        self.F_d = np.asarray(F_d, float)
        self.Y_d_bulk_basis = np.asarray(Y_d_bulk_basis, np.complex128)


class _DuckFitResult:
    __slots__ = ("bulk_state", "U_L_u", "U_R_u", "U_L_d", "U_R_d", "masses_down")

    def __init__(self, bulk_state, U_L_u, U_R_u, U_L_d, U_R_d, masses_down):
        self.bulk_state = bulk_state
        self.U_L_u = np.asarray(U_L_u, np.complex128)
        self.U_R_u = np.asarray(U_R_u, np.complex128)
        self.U_L_d = np.asarray(U_L_d, np.complex128)
        self.U_R_d = np.asarray(U_R_d, np.complex128)
        self.masses_down = np.asarray(masses_down, float)


# ---------------------------------------------------------------------------
# Worker globals (built once per process / per tile)
# ---------------------------------------------------------------------------
_W: dict = {}


def _worker_init():
    """Per-process one-time setup: PDG targets + Z-pole SM-fit reference."""
    from flavor_catalog_constraints import registry
    registry.discover()
    from flavor_catalog_constraints.primary.top_higgs_ew.T010 import (
        _load_t010_anchor, _RB_OBSERVABLE, _AB_OBSERVABLE,
    )
    from flavor_catalog_constraints.physics_adapters.zpole import (
        zpole_default_sm_inputs, zpole_inputs_with_bottom_radiator,
        zpole_sm_couplings,
    )
    _, _, _load_pdg_targets, _ = _run_rs_anarchy_helpers()
    anc = _load_t010_anchor("T010")
    sm_inputs = zpole_inputs_with_bottom_radiator(
        anc.sm_fit_values[_RB_OBSERVABLE].value, zpole_default_sm_inputs()
    )
    _W["targets"] = _load_pdg_targets()
    _W["sm_inputs"] = sm_inputs
    _W["sm_b"] = zpole_sm_couplings("b", sm_inputs)
    _W["RB_EXP"] = float(anc.r_b.value)
    _W["RB_SIG"] = float(anc.budgets[_RB_OBSERVABLE].combined_sigma)
    _W["AB_EXP"] = float(anc.a_b.value)
    _W["AB_SIG"] = float(anc.budgets[_AB_OBSERVABLE].combined_sigma)


def _run_tile(job):
    """One (scenario, M_KK_TeV, seed) tile -> list of row dicts."""
    from quarkConstraints.rs_ew_couplings import build_rs_ew_couplings
    from quarkConstraints.rs_ew_spectrum import RSEWSpectrum
    from flavor_catalog_constraints.physics_adapters.zpole import (
        zpole_shifted_couplings, zpole_evaluate_quark,
    )
    _ordered_svd, _build_kk_gluon_couplings, _, jarlskog_invariant = (
        _run_rs_anarchy_helpers()
    )
    deltaf2 = _deltaf2_helpers()
    scenario, M_KK_TeV, seed, n_draws = job
    sc = SCENARIOS[scenario]
    y_max = sc["y_max"]
    c_max = sc["c_max"]
    common_cd = (scenario == "S2")
    M_KK_GeV = M_KK_TeV * 1000.0
    xi_KK = DEFAULT_XI_KK
    Lambda_IR = M_KK_GeV / xi_KK
    epsilon = math.exp(-sc["L"]) if sc["L"] is not None else Lambda_IR / DEFAULT_K_GEV

    targets = _W["targets"]
    sm_inputs = _W["sm_inputs"]
    sm_b = _W["sm_b"]
    RB_EXP, RB_SIG = _W["RB_EXP"], _W["RB_SIG"]
    AB_EXP, AB_SIG = _W["AB_EXP"], _W["AB_SIG"]

    spectrum = RSEWSpectrum.build(
        lambda_ir_gev=Lambda_IR, k_gev=K_GEV, n_gauge_modes=N_GAUGE_MODES,
        quadrature_order=QUADRATURE_ORDER, model_label=EW_MODEL,
    )

    rng = np.random.default_rng(seed)
    c_lo_scan, c_hi_scan = -c_max, 0.5
    rows = []
    for _ in range(n_draws):
        c_u3 = float(rng.uniform(c_lo_scan, c_hi_scan))
        try:
            Y_u = _draw_bauer_matrix(rng, 0.1, y_max)
            Y_d = _draw_bauer_matrix(rng, 0.1, y_max)
            c_Q, c_u, c_d = _fn_c_values(c_u3, epsilon, targets, Y_u=Y_u, Y_d=Y_d,
                                         common_cd=common_cd, rng=rng, c_jitter=0.0)
            f_Q = f_IR(c_Q, epsilon)
            f_u = f_IR(c_u, epsilon)
            f_d = f_IR(c_d, epsilon)
            D_Q, D_u, D_d = np.diag(f_Q), np.diag(f_u), np.diag(f_d)
            M_u = BAUER_REPO_MASS_PREFAC_GEV * D_Q @ Y_u @ D_u
            M_d = BAUER_REPO_MASS_PREFAC_GEV * D_Q @ Y_d @ D_d
            U_L_u, m_up, U_R_u = _ordered_svd(M_u)
            U_L_d, m_dn, U_R_d = _ordered_svd(M_d)

            # --- mass + CKM consistency (Bauer's per-point chi^2 acceptance) ---
            # Bauer (0912.1625 Sec 5.2) REJECTS draws whose Froggatt-Nielsen fit
            # does not reproduce {m_q, CKM, J} (they cut chi^2/dof > 11.5/10), so
            # their plotted ensemble is mass+CKM-consistent.  Our forward scan
            # otherwise keeps every draw, including the mass+CKM-FAILING ones whose
            # anarchic Yukawa minors rail c_Q3 to the IR floor (UNgated median repo
            # c_Q3 ~ +0.08 with a huge IR tail; the GATED median rises to ~+0.28,
            # inside Bauer's +0.34 +/- 0.32 band).  We therefore record passes_pdg
            # with the SAME factor tolerances as run_rs_anarchy / anarchic_bauer_s1
            # so the Z->bb overlay can be restricted to Bauer's accepted ensemble.
            ckm = U_L_u.conj().T @ U_L_d
            up_log = np.log(np.maximum(m_up, 1e-30) / targets["up_masses_GeV"])
            dn_log = np.log(np.maximum(m_dn, 1e-30) / targets["down_masses_GeV"])
            ckm_log = np.array([
                abs(math.log(max(abs(ckm[0, 1]), 1e-30) / targets["abs_V_us"])),
                abs(math.log(max(abs(ckm[1, 2]), 1e-30) / targets["abs_V_cb"])),
                abs(math.log(max(abs(ckm[0, 2]), 1e-30) / targets["abs_V_ub"])),
            ])
            j_log = abs(math.log(max(abs(jarlskog_invariant(ckm)), 1e-30)
                                 / abs(targets["J"])))
            passes_pdg = bool(
                (np.abs(up_log) <= math.log(MASS_FACTOR)).all()
                and (np.abs(dn_log) <= math.log(MASS_FACTOR)).all()
                and (ckm_log <= math.log(CKM_FACTOR)).all()
                and j_log <= math.log(J_FACTOR)
            )

            # --- eps_K (verbatim Delta-F=2 path, same as anarchic_bauer_s1) ---
            couplings = _build_kk_gluon_couplings(
                M_KK=M_KK_GeV, xi_KK=xi_KK, f_Q=f_Q, f_u=f_u, f_d=f_d,
                U_L_u=U_L_u, U_L_d=U_L_d, U_R_u=U_R_u, U_R_d=U_R_d,
            )
            df2 = deltaf2["evaluate_delta_f2_constraints"](
                couplings, M_KK=M_KK_GeV, xi_KK=xi_KK
            )
            eps_k_summary = df2.by_system["K"]
            ratio_eps_K = float(eps_k_summary.ratio_to_bound)
            eps_K_np = float(abs(eps_k_summary.effective_amplitude))

            # --- Z->bb: TOTAL mass-basis [2,2] shift (gauge-KK + fermion-KK) ---
            bs = _DuckBulkState(c_Q, c_u, c_d, f_Q, f_u, f_d, Y_d)
            fr = _DuckFitResult(bs, U_L_u, U_R_u, U_L_d, U_R_d, m_dn)
            rsc = build_rs_ew_couplings(
                fr, spectrum=spectrum, include_fermion_kk_mixing=True,
                overlap_rel_tol=OVERLAP_RTOL, min_overlap_modes=MIN_OVERLAP_MODES,
                max_overlap_modes=MAX_OVERLAP_MODES, model_label=EW_MODEL,
            )
            dgL = float(np.real(rsc.z_delta_g_L_d[2, 2]))
            dgR = float(np.real(rsc.z_delta_g_R_d[2, 2]))
            shifted = zpole_shifted_couplings(sm_b, delta_g_left=dgL, delta_g_right=dgR)
            pred = zpole_evaluate_quark("b", {"b": shifted}, inputs=sm_inputs)
            R_b = float(pred.r_q)
            A_b = float(pred.a_q)
            pull_Rb = abs(R_b - RB_EXP) / RB_SIG
            pull_Ab = abs(A_b - AB_EXP) / AB_SIG
            passes_Zbb = bool(pull_Rb <= Z99 and pull_Ab <= Z99)
        except (ValueError, np.linalg.LinAlgError, FloatingPointError):
            continue
        except Exception:  # noqa: BLE001 -- robust per-draw isolation
            continue
        rows.append(dict(
            scenario=scenario,
            M_KK_GeV=M_KK_GeV,
            M_KK_TeV=float(M_KK_TeV),
            c_u3=c_u3,
            c_Q3=float(c_Q[2]),
            ratio_eps_K=ratio_eps_K,
            eps_K_np=eps_K_np,
            delta_g_L_b=dgL,
            delta_g_R_b=dgR,
            R_b_pred=R_b,
            A_b_pred=A_b,
            pull_Rb=float(pull_Rb),
            pull_Ab=float(pull_Ab),
            passes_Zbb=passes_Zbb,
            passes_pdg=passes_pdg,
        ))
    return rows


def main(argv=None):
    import multiprocessing as mp
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", type=str,
                   default="scan_outputs/anarchic_reproduction/anarchic_bauer_s1_zbb.parquet")
    p.add_argument("--scenario", type=str, default="S1", choices=sorted(SCENARIOS.keys()))
    p.add_argument("--per-tile", type=int, default=2000)
    p.add_argument("--m-kk-tev", type=str, default="1,1.5,2,3,4,5,6,7,8,10")
    p.add_argument("--base-seed", type=int, default=20260623)
    p.add_argument("--workers", type=int, default=max(1, len(os.sched_getaffinity(0))))
    args = p.parse_args(argv)

    mkk_tev = [float(t) for t in args.m_kk_tev.split(",") if t.strip()]
    out_path = Path(args.out)
    if not out_path.is_absolute():
        out_path = REPO / out_path
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # One job per tile; each worker builds its spectrum once and runs all draws.
    jobs = []
    for idx, t in enumerate(mkk_tev):
        seed = args.base_seed + 1009 * idx  # deterministic, scenario-stable
        jobs.append((args.scenario, t, seed, args.per_tile))

    print(f"[zbb] scenario {args.scenario}: {len(jobs)} tiles x {args.per_tile} draws "
          f"over {args.workers} workers (RS-EW spectrum per tile)", flush=True)
    t0 = time.time()
    all_rows = []

    def _report(rows):
        if not rows:
            return
        m = rows[0]["M_KK_TeV"]
        fr = np.mean([r["passes_Zbb"] for r in rows])
        print(f"  M_KK={m:5.2f} TeV : {len(rows):5d} draws  Zbb-pass(99%)={fr:5.1%}  "
              f"({time.time()-t0:.0f}s elapsed)", flush=True)

    if args.workers <= 1:
        _worker_init()
        for job in jobs:
            rows = _run_tile(job)
            all_rows.extend(rows)
            _report(rows)
    else:
        ctx = mp.get_context("spawn")
        with ctx.Pool(args.workers, initializer=_worker_init) as pool:
            for rows in pool.imap_unordered(_run_tile, jobs):
                all_rows.extend(rows)
                _report(rows)

    df = pd.DataFrame(all_rows)
    df.to_parquet(out_path, index=False)
    print(f"\n[zbb] wrote {len(df):,} rows -> {out_path} in {time.time()-t0:.0f}s")

    # --- sanity: Zbb pass-fraction vs M_KK (should RISE with M_KK) ---
    print("\n==== Z->bb pass fraction (99% CL) vs M_KK ====")
    for m in sorted(df.M_KK_TeV.unique()):
        g = df[df.M_KK_TeV == m]
        print(f"  M_KK={m:5.2f} TeV : n={len(g):5d}  passes_Zbb={g.passes_Zbb.mean():6.1%}  "
              f"med|dgL|={g.delta_g_L_b.abs().median():.2e} med|dgR|={g.delta_g_R_b.abs().median():.2e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
