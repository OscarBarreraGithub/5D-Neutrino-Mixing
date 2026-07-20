#!/usr/bin/env python
"""Rank-one/U(2) "wrong ensemble" lane for RS epsilon_K.

This is a research diagnostic for RS-FLAVOR-ALIGNMENT-2026-07.  It compares
the usual flat Bauer-anarchic ensemble against a sequential low-rank Yukawa
ensemble with exact light-family bulk degeneracies

    c_Q1 = c_Q2,     c_d1 = c_d2

so light-family U(2) rotations cannot generate a down-sector 1-2 KK-gluon
coupling.  The light mass ratios are placed in Yukawa singular values, not in
profile splittings.
"""
from __future__ import annotations

import argparse
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
from scipy.optimize import brentq

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from quarkConstraints.deltaf2 import (  # noqa: E402
    compute_delta_f2_wilsons,
    evaluate_delta_f2_constraints,
    evaluate_delta_mk_with_running,
)
from scripts.anarchic_bauer_s1 import (  # noqa: E402
    SCENARIOS,
    _draw_bauer_matrix,
    _fn_c_values,
)
from scripts.run_rs_anarchy import (  # noqa: E402
    DEFAULT_K_GEV,
    DEFAULT_V_GEV,
    _build_kk_gluon_couplings,
    _load_pdg_targets,
    _ordered_svd,
    jarlskog_invariant,
)
from warpConfig.wavefuncs import f_IR  # noqa: E402

LAM = 0.2250
WOLF_A = 0.826
PERTURBATIVE_CAP = 3.0


@dataclass(frozen=True)
class RankOneDraw:
    Y_u: np.ndarray
    Y_d: np.ndarray
    s_u: np.ndarray
    s_d: np.ndarray
    q_u: np.ndarray
    q_d: np.ndarray
    r_u: np.ndarray
    r_d: np.ndarray
    max_abs_y: float


@dataclass(frozen=True)
class BulkAssignment:
    c_Q: np.ndarray
    c_u: np.ndarray
    c_d: np.ndarray
    f_Q: np.ndarray
    f_u: np.ndarray
    f_d: np.ndarray


def _rot(i: int, j: int, theta: float, phi: float = 0.0) -> np.ndarray:
    r = np.eye(3, dtype=np.complex128)
    c = math.cos(theta)
    s = math.sin(theta) * np.exp(-1j * phi)
    r[i, i] = c
    r[j, j] = c
    r[i, j] = s
    r[j, i] = -s.conjugate()
    return r


def _haar_u2(rng: np.random.Generator) -> np.ndarray:
    z = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
    q, rr = np.linalg.qr(z)
    q *= np.exp(-1j * np.angle(np.diag(rr)))
    u = np.eye(3, dtype=np.complex128)
    u[:2, :2] = q
    return u


def _ckm_target(
    delta_scale: float = 1.0,
    *,
    boost_23: float = 1.0,
    boost_13: float = 1.0,
) -> np.ndarray:
    """CKM-shaped unitary spurion.

    ``boost_23`` and ``boost_13`` compensate for the profile dressing in
    ``F_Q Y_f``: the input left frame needs somewhat larger third-family angles
    than the physical SVD CKM angles when the light doublet is UV-localized.
    """
    s12 = 0.2243
    s23 = min(0.95, 0.0408 * boost_23)
    s13 = min(0.95, 0.00382 * boost_13)
    c12 = math.sqrt(1.0 - s12 * s12)
    c23 = math.sqrt(1.0 - s23 * s23)
    c13 = math.sqrt(1.0 - s13 * s13)
    j_target = 3.0e-5
    denom = c12 * c23 * c13 * c13 * s12 * s23 * s13
    sin_delta = max(-0.999, min(0.999, delta_scale * j_target / denom))
    delta = math.asin(sin_delta)
    e_pos = np.exp(1j * delta)
    e_neg = np.exp(-1j * delta)
    return np.array(
        [
            [c12 * c13, s12 * c13, s13 * e_neg],
            [
                -s12 * c23 - c12 * s23 * s13 * e_pos,
                c12 * c23 - s12 * s23 * s13 * e_pos,
                s23 * c13,
            ],
            [
                s12 * s23 - c12 * c23 * s13 * e_pos,
                -c12 * s23 - s12 * c23 * s13 * e_pos,
                c23 * c13,
            ],
        ],
        dtype=np.complex128,
    )


def _third_family_leak(
    rng: np.random.Generator,
    *,
    leak: float,
    theta23_power: float = 2.0,
    theta13_power: float = 3.0,
) -> np.ndarray:
    if leak == 0.0:
        return np.eye(3, dtype=np.complex128)
    th23 = leak * (LAM**theta23_power) * rng.lognormal(0.0, 0.25)
    th13 = leak * (LAM**theta13_power) * rng.lognormal(0.0, 0.35)
    ph23 = rng.uniform(0.0, 2.0 * math.pi)
    ph13 = rng.uniform(0.0, 2.0 * math.pi)
    return _rot(1, 2, th23, ph23) @ _rot(0, 2, th13, ph13)


def _rankone_yukawa(left_frame: np.ndarray, right_frame: np.ndarray, s: np.ndarray) -> np.ndarray:
    y = np.zeros((3, 3), dtype=np.complex128)
    for a in range(3):
        y += s[a] * np.outer(left_frame[:, a], right_frame[:, a].conjugate())
    return y


def _draw_rankone_u2_yukawas(
    rng: np.random.Generator,
    targets: dict,
    *,
    rh_down_leak: float,
    rh_up_leak: float,
    left_23_boost: float,
    left_13_boost: float,
    left_cp_boost: float,
    perturbativity_mode: str,
) -> RankOneDraw | None:
    """Sequential low-rank draw with shared light-family U(2)_Q plane."""
    mu = targets["up_masses_GeV"]
    md = targets["down_masses_GeV"]

    su3 = rng.uniform(1.15, 3.0)
    sd3 = rng.uniform(0.80, 3.0)
    eps_u = float(np.clip(rng.lognormal(math.log(0.13), 0.25), 0.085, 0.28))
    eps_d = float(np.clip(rng.lognormal(math.log(0.12), 0.35), 0.045, 0.32))
    su2 = su3 * eps_u
    sd2 = sd3 * eps_d
    # The light 1-2 mass hierarchy is deliberately placed in the singular
    # values.  With degenerate light profiles this makes m1/m2 = s1/s2.
    su1 = su2 * float(mu[0] / mu[1])
    sd1 = sd2 * float(md[0] / md[1])
    s_u = np.array([su1, su2, su3], dtype=float)
    s_d = np.array([sd1, sd2, sd3], dtype=float)

    q_u = _haar_u2(rng)
    q_d = q_u @ _ckm_target(
        delta_scale=left_cp_boost * rng.lognormal(0.0, 0.08),
        boost_23=left_23_boost,
        boost_13=left_13_boost,
    )
    r_u = _haar_u2(rng) @ _third_family_leak(rng, leak=rh_up_leak)
    r_d = _haar_u2(rng) @ _third_family_leak(rng, leak=rh_down_leak)

    Y_u = _rankone_yukawa(q_u, r_u, s_u)
    Y_d = _rankone_yukawa(q_d, r_d, s_d)
    max_abs_y = float(max(np.abs(Y_u).max(), np.abs(Y_d).max()))
    if max_abs_y > PERTURBATIVE_CAP:
        if perturbativity_mode == "reject":
            return None
        if perturbativity_mode != "rescale":
            raise ValueError(f"unknown perturbativity mode: {perturbativity_mode}")
        scale = PERTURBATIVE_CAP / max_abs_y
        Y_u = Y_u * scale
        Y_d = Y_d * scale
        s_u = s_u * scale
        s_d = s_d * scale
        max_abs_y = PERTURBATIVE_CAP
    return RankOneDraw(
        Y_u=Y_u,
        Y_d=Y_d,
        s_u=s_u,
        s_d=s_d,
        q_u=q_u,
        q_d=q_d,
        r_u=r_u,
        r_d=r_d,
        max_abs_y=max_abs_y,
    )


def _invert_f_IR(f_target: float, epsilon: float, c_lo: float = -2.6, c_hi: float = 0.95) -> float:
    f_at_lo = float(f_IR(np.array([c_lo]), epsilon)[0])
    f_at_hi = float(f_IR(np.array([c_hi]), epsilon)[0])
    if not (f_at_hi <= f_target <= f_at_lo):
        raise ValueError(f"target f={f_target:.4g} outside [{f_at_hi:.4g}, {f_at_lo:.4g}]")
    g = lambda c: float(f_IR(np.array([c]), epsilon)[0]) - f_target
    return float(brentq(g, c_lo, c_hi, xtol=1e-10, rtol=1e-12))


def _assign_u2_bulk_masses(
    draw: RankOneDraw,
    epsilon: float,
    targets: dict,
    *,
    c_Q12: float,
    c_Q3: float,
) -> BulkAssignment:
    """Assign profiles with exact c_Q1=c_Q2 and c_d1=c_d2.

    The right-up light profiles are also taken degenerate by default here.  This
    is stronger than needed for epsilon_K, but keeps the "light mass hierarchy
    lives in singular values" logic explicit in both up and down sectors.
    """
    mu = targets["up_masses_GeV"]
    md = targets["down_masses_GeV"]
    f_Q = f_IR(np.array([c_Q12, c_Q12, c_Q3], dtype=float), epsilon)
    fQ_l = float(f_Q[0])
    fQ_3 = float(f_Q[2])
    su1, su2, su3 = (float(x) for x in draw.s_u)
    sd1, sd2, sd3 = (float(x) for x in draw.s_d)
    del su1, sd1  # ratios were used when drawing s1/s2

    f_u12 = float(mu[1] / (DEFAULT_V_GEV * fQ_l * su2))
    f_u3 = float(mu[2] / (DEFAULT_V_GEV * fQ_3 * su3))
    f_d12 = float(md[1] / (DEFAULT_V_GEV * fQ_l * sd2))
    f_d3 = float(md[2] / (DEFAULT_V_GEV * fQ_3 * sd3))

    c_u12 = _invert_f_IR(f_u12, epsilon)
    c_u3 = _invert_f_IR(f_u3, epsilon)
    c_d12 = _invert_f_IR(f_d12, epsilon)
    c_d3 = _invert_f_IR(f_d3, epsilon)
    c_Q = np.array([c_Q12, c_Q12, c_Q3], dtype=float)
    c_u = np.array([c_u12, c_u12, c_u3], dtype=float)
    c_d = np.array([c_d12, c_d12, c_d3], dtype=float)
    return BulkAssignment(
        c_Q=c_Q,
        c_u=c_u,
        c_d=c_d,
        f_Q=f_IR(c_Q, epsilon),
        f_u=f_IR(c_u, epsilon),
        f_d=f_IR(c_d, epsilon),
    )


def _pdg_residuals(
    masses_up: np.ndarray,
    masses_down: np.ndarray,
    ckm: np.ndarray,
    J: float,
    targets: dict,
    *,
    mass_factor: float,
    ckm_factor: float,
    j_factor: float,
) -> tuple[bool, dict]:
    up_log = np.log(np.maximum(masses_up, 1e-30) / targets["up_masses_GeV"])
    down_log = np.log(np.maximum(masses_down, 1e-30) / targets["down_masses_GeV"])
    ckm_log = np.array(
        [
            abs(math.log(max(abs(ckm[0, 1]), 1e-30) / targets["abs_V_us"])),
            abs(math.log(max(abs(ckm[1, 2]), 1e-30) / targets["abs_V_cb"])),
            abs(math.log(max(abs(ckm[0, 2]), 1e-30) / targets["abs_V_ub"])),
        ],
        dtype=float,
    )
    j_log = abs(math.log(max(abs(J), 1e-30) / abs(targets["J"])))
    passes = bool(
        (np.abs(up_log) <= math.log(mass_factor)).all()
        and (np.abs(down_log) <= math.log(mass_factor)).all()
        and (ckm_log <= math.log(ckm_factor)).all()
        and j_log <= math.log(j_factor)
    )
    return passes, {
        "up_log_max": float(np.max(np.abs(up_log))),
        "down_log_max": float(np.max(np.abs(down_log))),
        "ckm_log_max": float(np.max(ckm_log)),
        "j_log": float(j_log),
    }


def _evaluate_point(
    *,
    lane: str,
    Y_u: np.ndarray,
    Y_d: np.ndarray,
    bulk: BulkAssignment,
    M_KK_GeV: float,
    xi_KK: float,
    targets: dict,
    mass_factor: float,
    ckm_factor: float,
    j_factor: float,
    meta: dict | None = None,
) -> dict:
    D_Q = np.diag(bulk.f_Q)
    D_u = np.diag(bulk.f_u)
    D_d = np.diag(bulk.f_d)
    M_u = DEFAULT_V_GEV * D_Q @ Y_u @ D_u
    M_d = DEFAULT_V_GEV * D_Q @ Y_d @ D_d
    U_L_u, m_up, U_R_u = _ordered_svd(M_u)
    U_L_d, m_down, U_R_d = _ordered_svd(M_d)
    ckm = U_L_u.conjugate().T @ U_L_d
    J = float(jarlskog_invariant(ckm))
    passes_pdg, residuals = _pdg_residuals(
        m_up,
        m_down,
        ckm,
        J,
        targets,
        mass_factor=mass_factor,
        ckm_factor=ckm_factor,
        j_factor=j_factor,
    )
    couplings = _build_kk_gluon_couplings(
        M_KK=M_KK_GeV,
        xi_KK=xi_KK,
        f_Q=bulk.f_Q,
        f_u=bulk.f_u,
        f_d=bulk.f_d,
        U_L_u=U_L_u,
        U_L_d=U_L_d,
        U_R_u=U_R_u,
        U_R_d=U_R_d,
    )
    df2 = evaluate_delta_f2_constraints(couplings, M_KK=M_KK_GeV, xi_KK=xi_KK)
    bs = df2.by_system
    passes_dm_K = False
    try:
        unev = compute_delta_f2_wilsons(couplings, M_KK=M_KK_GeV, xi_KK=xi_KK)
        w_k = next((w for w in unev if w.input.key == "epsilon_k"), None)
        if w_k is not None:
            dmk = evaluate_delta_mk_with_running(w_k, mu_had=2.0)
            ratio_dm_K = float(dmk.ratio_to_exp)
            passes_dm_K = bool(dmk.passes)
        else:
            ratio_dm_K = float("nan")
    except Exception:
        ratio_dm_K = float("nan")

    GL = np.asarray(couplings.left_down)
    GR = np.asarray(couplings.right_down)
    prod = complex(GL[0, 1] * GR[0, 1])
    c4_abs = float(abs(prod) / (M_KK_GeV**2))
    phi = float(np.angle(prod)) if c4_abs > 0.0 else 0.0
    row = {
        "lane": lane,
        "passes_pdg": passes_pdg,
        "passes_eps_K": bool(bs["K"].passes),
        "passes_deltaf2_all": bool(df2.passes_all and passes_dm_K),
        "ratio_eps_K": float(bs["K"].ratio_to_bound),
        "ratio_dm_K": ratio_dm_K,
        "passes_dm_K": passes_dm_K,
        "ratio_B_d": float(bs["B_d"].ratio_to_bound),
        "ratio_B_s": float(bs["B_s"].ratio_to_bound),
        "ratio_D": float(bs["D"].ratio_to_bound),
        "GL12_abs": float(abs(GL[0, 1])),
        "GR12_abs": float(abs(GR[0, 1])),
        "C4abs_12": c4_abs,
        "Phi12": phi,
        "sinPhi12_abs": float(abs(math.sin(phi))),
        "abs_V_us": float(abs(ckm[0, 1])),
        "abs_V_cb": float(abs(ckm[1, 2])),
        "abs_V_ub": float(abs(ckm[0, 2])),
        "J": J,
        "m_u": float(m_up[0]),
        "m_c": float(m_up[1]),
        "m_t": float(m_up[2]),
        "m_d": float(m_down[0]),
        "m_s": float(m_down[1]),
        "m_b": float(m_down[2]),
        "c_Q1": float(bulk.c_Q[0]),
        "c_Q2": float(bulk.c_Q[1]),
        "c_Q3": float(bulk.c_Q[2]),
        "c_u1": float(bulk.c_u[0]),
        "c_u2": float(bulk.c_u[1]),
        "c_u3": float(bulk.c_u[2]),
        "c_d1": float(bulk.c_d[0]),
        "c_d2": float(bulk.c_d[1]),
        "c_d3": float(bulk.c_d[2]),
        "f_Q1": float(bulk.f_Q[0]),
        "f_Q2": float(bulk.f_Q[1]),
        "f_Q3": float(bulk.f_Q[2]),
        "f_d1": float(bulk.f_d[0]),
        "f_d2": float(bulk.f_d[1]),
        "f_d3": float(bulk.f_d[2]),
        **residuals,
    }
    if meta:
        row.update(meta)
    return row


def _run_u2_tile(
    *,
    M_KK_GeV: float,
    xi_KK: float,
    n_draws: int,
    seed: int,
    targets: dict,
    mass_factor: float,
    ckm_factor: float,
    j_factor: float,
    c_Q12: float,
    c_Q3: float,
    rh_down_leak: float,
    rh_up_leak: float,
    left_23_boost: float,
    left_13_boost: float,
    left_cp_boost: float,
    perturbativity_mode: str,
) -> tuple[list[dict], dict]:
    rng = np.random.default_rng(seed)
    epsilon = (M_KK_GeV / xi_KK) / DEFAULT_K_GEV
    rows: list[dict] = []
    perturb_rejects = 0
    fit_rejects = 0
    eval_rejects = 0
    for sample_idx in range(n_draws):
        draw = _draw_rankone_u2_yukawas(
            rng,
            targets,
            rh_down_leak=rh_down_leak,
            rh_up_leak=rh_up_leak,
            left_23_boost=left_23_boost,
            left_13_boost=left_13_boost,
            left_cp_boost=left_cp_boost,
            perturbativity_mode=perturbativity_mode,
        )
        if draw is None:
            perturb_rejects += 1
            continue
        try:
            bulk = _assign_u2_bulk_masses(
                draw,
                epsilon,
                targets,
                c_Q12=c_Q12,
                c_Q3=c_Q3,
            )
        except (ValueError, FloatingPointError):
            fit_rejects += 1
            continue
        try:
            row = _evaluate_point(
                lane="rankone_u2",
                Y_u=draw.Y_u,
                Y_d=draw.Y_d,
                bulk=bulk,
                M_KK_GeV=M_KK_GeV,
                xi_KK=xi_KK,
                targets=targets,
                mass_factor=mass_factor,
                ckm_factor=ckm_factor,
                j_factor=j_factor,
                meta={
                    "sample_idx": sample_idx,
                    "M_KK_GeV": float(M_KK_GeV),
                    "M_KK_TeV": float(M_KK_GeV / 1000.0),
                    "max_abs_y": draw.max_abs_y,
                    "s_u1": float(draw.s_u[0]),
                    "s_u2": float(draw.s_u[1]),
                    "s_u3": float(draw.s_u[2]),
                    "s_d1": float(draw.s_d[0]),
                    "s_d2": float(draw.s_d[1]),
                    "s_d3": float(draw.s_d[2]),
                },
            )
        except (ValueError, np.linalg.LinAlgError, FloatingPointError):
            eval_rejects += 1
            continue
        except Exception:
            eval_rejects += 1
            continue
        rows.append(row)
    summary = {
        "attempts": n_draws,
        "perturb_rejects": perturb_rejects,
        "fit_rejects": fit_rejects,
        "eval_rejects": eval_rejects,
        "evaluated": len(rows),
    }
    return rows, summary


def _run_flat_tile(
    *,
    M_KK_GeV: float,
    xi_KK: float,
    n_draws: int,
    seed: int,
    targets: dict,
    mass_factor: float,
    ckm_factor: float,
    j_factor: float,
) -> tuple[list[dict], dict]:
    sc = SCENARIOS["S1"]
    rng = np.random.default_rng(seed)
    epsilon = (M_KK_GeV / xi_KK) / DEFAULT_K_GEV
    c_lo_scan, c_hi_scan = -float(sc["c_max"]), 0.5
    rows: list[dict] = []
    fit_rejects = 0
    eval_rejects = 0
    for sample_idx in range(n_draws):
        try:
            c_u3 = float(rng.uniform(c_lo_scan, c_hi_scan))
            Y_u = _draw_bauer_matrix(rng, 0.1, float(sc["y_max"]))
            Y_d = _draw_bauer_matrix(rng, 0.1, float(sc["y_max"]))
            c_Q, c_u, c_d = _fn_c_values(
                c_u3,
                epsilon,
                targets,
                Y_u=Y_u,
                Y_d=Y_d,
                common_cd=False,
                rng=rng,
                c_jitter=0.0,
            )
            bulk = BulkAssignment(
                c_Q=c_Q,
                c_u=c_u,
                c_d=c_d,
                f_Q=f_IR(c_Q, epsilon),
                f_u=f_IR(c_u, epsilon),
                f_d=f_IR(c_d, epsilon),
            )
        except (ValueError, FloatingPointError):
            fit_rejects += 1
            continue
        try:
            rows.append(
                _evaluate_point(
                    lane="flat_bauer_s1",
                    Y_u=Y_u,
                    Y_d=Y_d,
                    bulk=bulk,
                    M_KK_GeV=M_KK_GeV,
                    xi_KK=xi_KK,
                    targets=targets,
                    mass_factor=mass_factor,
                    ckm_factor=ckm_factor,
                    j_factor=j_factor,
                    meta={
                        "sample_idx": sample_idx,
                        "M_KK_GeV": float(M_KK_GeV),
                        "M_KK_TeV": float(M_KK_GeV / 1000.0),
                        "max_abs_y": float(max(np.abs(Y_u).max(), np.abs(Y_d).max())),
                    },
                )
            )
        except (ValueError, np.linalg.LinAlgError, FloatingPointError):
            eval_rejects += 1
            continue
        except Exception:
            eval_rejects += 1
            continue
    return rows, {
        "attempts": n_draws,
        "fit_rejects": fit_rejects,
        "eval_rejects": eval_rejects,
        "evaluated": len(rows),
    }


def _finite(values: Iterable[float]) -> np.ndarray:
    arr = np.asarray(list(values), dtype=float)
    return arr[np.isfinite(arr)]


def _median(rows: list[dict], key: str, *, pdg_only: bool = True) -> float:
    vals = _finite(r[key] for r in rows if (r.get("passes_pdg") or not pdg_only))
    return float(np.median(vals)) if vals.size else float("nan")


def _frac(rows: list[dict], pred, *, pdg_only: bool = True) -> float:
    denom = [r for r in rows if (r.get("passes_pdg") or not pdg_only)]
    if not denom:
        return float("nan")
    return float(sum(1 for r in denom if pred(r)) / len(denom))


def _summarize_tile(
    *,
    mkk_tev: float,
    u2_rows: list[dict],
    flat_rows: list[dict],
    u2_counts: dict,
    flat_counts: dict,
) -> dict:
    flat_c4 = _median(flat_rows, "C4abs_12")
    flat_gl = _median(flat_rows, "GL12_abs")
    flat_gr = _median(flat_rows, "GR12_abs")
    u2_c4 = _median(u2_rows, "C4abs_12")
    u2_gl = _median(u2_rows, "GL12_abs")
    u2_gr = _median(u2_rows, "GR12_abs")
    u2_pdg = [r for r in u2_rows if r.get("passes_pdg")]
    flat_pdg = [r for r in flat_rows if r.get("passes_pdg")]
    cQ_split = _finite(abs(r["c_Q1"] - r["c_Q2"]) for r in u2_rows)
    cd_split = _finite(abs(r["c_d1"] - r["c_d2"]) for r in u2_rows)
    return {
        "M_KK_TeV": mkk_tev,
        "u2_evaluated": len(u2_rows),
        "flat_evaluated": len(flat_rows),
        "u2_pdg_pass": len(u2_pdg),
        "flat_pdg_pass": len(flat_pdg),
        "u2_pdg_pass_fraction": len(u2_pdg) / max(1, len(u2_rows)),
        "flat_pdg_pass_fraction": len(flat_pdg) / max(1, len(flat_rows)),
        "u2_epsK_pass_fraction_pdg": _frac(u2_rows, lambda r: r["passes_eps_K"]),
        "flat_epsK_pass_fraction_pdg": _frac(flat_rows, lambda r: r["passes_eps_K"]),
        "u2_all_deltaf2_pass_fraction_pdg": _frac(u2_rows, lambda r: r["passes_deltaf2_all"]),
        "flat_all_deltaf2_pass_fraction_pdg": _frac(flat_rows, lambda r: r["passes_deltaf2_all"]),
        "u2_median_C4abs_12": u2_c4,
        "flat_median_C4abs_12": flat_c4,
        "C4_suppression_u2_over_flat": u2_c4 / flat_c4 if flat_c4 > 0 else float("nan"),
        "u2_median_GL12_abs": u2_gl,
        "flat_median_GL12_abs": flat_gl,
        "GL12_suppression_u2_over_flat": u2_gl / flat_gl if flat_gl > 0 else float("nan"),
        "u2_median_GR12_abs": u2_gr,
        "flat_median_GR12_abs": flat_gr,
        "GR12_suppression_u2_over_flat": u2_gr / flat_gr if flat_gr > 0 else float("nan"),
        "u2_median_abs_sinPhi12": _median(u2_rows, "sinPhi12_abs"),
        "u2_epsK_pass_median_abs_sinPhi12": float(
            np.median(
                _finite(
                    r["sinPhi12_abs"]
                    for r in u2_rows
                    if r.get("passes_pdg") and r.get("passes_eps_K")
                )
            )
        )
        if any(r.get("passes_pdg") and r.get("passes_eps_K") for r in u2_rows)
        else float("nan"),
        "u2_epsK_pass_frac_sinPhi_lt_0p1": _frac(
            [r for r in u2_rows if r.get("passes_eps_K")],
            lambda r: r["sinPhi12_abs"] < 0.1,
        ),
        "u2_max_cQ12_split": float(np.max(cQ_split)) if cQ_split.size else float("nan"),
        "u2_max_cd12_split": float(np.max(cd_split)) if cd_split.size else float("nan"),
        "u2_counts": u2_counts,
        "flat_counts": flat_counts,
    }


def _fmt(x: float, kind: str = "g") -> str:
    if not np.isfinite(x):
        return "nan"
    if kind == "pct":
        return f"{100.0 * x:.1f}%"
    if kind == "sci":
        return f"{x:.3e}"
    if kind == "ratio":
        return f"{x:.3g}"
    return f"{x:.4g}"


def _build_report(args, summaries: list[dict]) -> str:
    suppress_ok = all(
        s["C4_suppression_u2_over_flat"] <= 1.0e-2
        for s in summaries
        if np.isfinite(s["C4_suppression_u2_over_flat"])
    )
    low_mkk = [s for s in summaries if s["M_KK_TeV"] in (2.0, 3.0)]
    eps_ok = any(s["u2_epsK_pass_fraction_pdg"] > 0.5 for s in low_mkk)
    phase_tuned = any(
        s["u2_epsK_pass_median_abs_sinPhi12"] < 0.1
        for s in low_mkk
        if np.isfinite(s["u2_epsK_pass_median_abs_sinPhi12"])
    )
    split_ok = all(
        s["u2_max_cQ12_split"] < 1e-12 and s["u2_max_cd12_split"] < 1e-12
        for s in summaries
    )
    verdict = (
        "PASS"
        if suppress_ok and eps_ok and not phase_tuned and split_ok
        else "INCONCLUSIVE/FAIL"
    )

    lines = [
        "# Rank-one/U(2) RS epsilon_K lane",
        "",
        f"Verdict: **{verdict}**",
        "",
        "Configuration:",
        f"- draws per lane per M_KK: `{args.n_draws}`",
        f"- M_KK grid [TeV]: `{','.join(str(x) for x in args.mkk_tev)}`",
        f"- xi_KK: `{args.xi_kk}`",
        f"- fixed U(2) doublet profiles: `c_Q1=c_Q2={args.c_Q12}`, `c_Q3={args.c_Q3}`",
        f"- right-down third-family leakage: `{args.rh_down_leak}` times CKM-sized angles",
        f"- profile-dressed left spurion boosts: `V23 x {args.left_23_boost}`, `V13 x {args.left_13_boost}`, `CP x {args.left_cp_boost}`",
        f"- perturbativity: max |Y_ij| <= {PERTURBATIVE_CAP} with `{args.perturbativity_mode}` mode",
        "",
        "Summary table, using PDG-passing draws for medians and pass fractions:",
        "",
        "| M_KK [TeV] | U2 PDG pass | U2 eps_K pass | flat eps_K pass | med \\|C4\\| U2/flat | med \\|GL12\\| U2/flat | med \\|GR12\\| U2/flat | med \\|sin Phi12\\| U2 | c splits |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for s in summaries:
        lines.append(
            "| "
            + " | ".join(
                [
                    _fmt(s["M_KK_TeV"]),
                    f"{s['u2_pdg_pass']}/{s['u2_evaluated']} ({_fmt(s['u2_pdg_pass_fraction'], 'pct')})",
                    _fmt(s["u2_epsK_pass_fraction_pdg"], "pct"),
                    _fmt(s["flat_epsK_pass_fraction_pdg"], "pct"),
                    _fmt(s["C4_suppression_u2_over_flat"], "ratio"),
                    _fmt(s["GL12_suppression_u2_over_flat"], "ratio"),
                    _fmt(s["GR12_suppression_u2_over_flat"], "ratio"),
                    _fmt(s["u2_median_abs_sinPhi12"], "ratio"),
                    f"Q={_fmt(s['u2_max_cQ12_split'], 'sci')}, d={_fmt(s['u2_max_cd12_split'], 'sci')}",
                ]
            )
            + " |"
        )
    lines.extend(
        [
            "",
            "Detailed medians:",
            "",
            "| M_KK [TeV] | med \\|C4\\| U2 | med \\|C4\\| flat | med \\|GL12\\| U2 | med \\|GL12\\| flat | med \\|GR12\\| U2 | med \\|GR12\\| flat |",
            "|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for s in summaries:
        lines.append(
            "| "
            + " | ".join(
                [
                    _fmt(s["M_KK_TeV"]),
                    _fmt(s["u2_median_C4abs_12"], "sci"),
                    _fmt(s["flat_median_C4abs_12"], "sci"),
                    _fmt(s["u2_median_GL12_abs"], "sci"),
                    _fmt(s["flat_median_GL12_abs"], "sci"),
                    _fmt(s["u2_median_GR12_abs"], "sci"),
                    _fmt(s["flat_median_GR12_abs"], "sci"),
                ]
            )
            + " |"
        )
    lines.extend(["", "Diagnostics:"])
    for s in summaries:
        u2_counts = s["u2_counts"]
        flat_counts = s["flat_counts"]
        reject_frac = u2_counts["perturb_rejects"] / max(1, u2_counts["attempts"])
        lines.append(
            f"- M_KK={s['M_KK_TeV']:.1f} TeV: rank-one perturbative reject fraction "
            f"{_fmt(reject_frac, 'pct')}; U2 all-DeltaF2 pass fraction "
            f"{_fmt(s['u2_all_deltaf2_pass_fraction_pdg'], 'pct')} vs flat "
            f"{_fmt(s['flat_all_deltaf2_pass_fraction_pdg'], 'pct')}; "
            f"U2 fit/eval rejects {u2_counts['fit_rejects']}/{u2_counts['eval_rejects']}, "
            f"flat fit/eval rejects {flat_counts['fit_rejects']}/{flat_counts['eval_rejects']}."
        )
    lines.extend(
        [
            "",
            "Interpretation:",
            f"- Median |C4| suppression >= 1e-2 at all measured M_KK: **{suppress_ok}**.",
            f"- eps_K pass fraction exceeds 50% at 2-3 TeV in the U(2) lane: **{eps_ok}**.",
            f"- Phase tuning diagnostic, median |sin Phi12| among low-M_KK eps_K-pass draws < 0.1: **{phase_tuned}**.",
            f"- Exact U(2) c-degeneracy retained in the fit: **{split_ok}**.",
            "",
            "The structural check is the last line: if `c_Q1-c_Q2` or `c_d1-c_d2` is nonzero, the lane has failed by construction because the c-fit reintroduced the anarchic source of the 1-2 coupling.",
        ]
    )
    return "\n".join(lines) + "\n"


def _parse_mkk(text: str) -> list[float]:
    return [float(x) for x in text.split(",") if x.strip()]


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--out",
        default=".orchestration/runs/RS-FLAVOR-ALIGNMENT-2026-07/rankone_u2_lane_report.md",
        help="Markdown report path.",
    )
    p.add_argument("--n-draws", type=int, default=1000)
    p.add_argument("--mkk-tev", type=_parse_mkk, default=_parse_mkk("2,3,5"))
    p.add_argument("--seed", type=int, default=20260702)
    p.add_argument("--xi-kk", type=float, default=1.0)
    p.add_argument("--mass-factor", type=float, default=3.0)
    p.add_argument("--ckm-factor", type=float, default=3.0)
    p.add_argument("--j-factor", type=float, default=5.0)
    p.add_argument("--c-Q12", type=float, default=0.50)
    p.add_argument("--c-Q3", type=float, default=0.25)
    p.add_argument("--left-23-boost", type=float, default=4.0)
    p.add_argument("--left-13-boost", type=float, default=4.0)
    p.add_argument("--left-cp-boost", type=float, default=20.0)
    p.add_argument("--rh-down-leak", type=float, default=1.0)
    p.add_argument("--rh-up-leak", type=float, default=1.0)
    p.add_argument("--perturbativity-mode", choices=["reject", "rescale"], default="reject")
    args = p.parse_args(argv)

    targets = _load_pdg_targets()
    summaries: list[dict] = []
    print("[rankone_u2] starting")
    print(f"[rankone_u2] M_KK grid TeV = {args.mkk_tev}")
    print(f"[rankone_u2] n_draws per lane/tile = {args.n_draws}")
    print(f"[rankone_u2] xi_KK = {args.xi_kk}")
    print(f"[rankone_u2] exact degeneracies: c_Q1=c_Q2={args.c_Q12}, c_d1=c_d2 per draw")
    for idx, mkk_tev in enumerate(args.mkk_tev):
        M_KK_GeV = 1000.0 * mkk_tev
        print(f"[rankone_u2] tile M_KK={mkk_tev:.1f} TeV: rank-one/U(2)")
        u2_rows, u2_counts = _run_u2_tile(
            M_KK_GeV=M_KK_GeV,
            xi_KK=args.xi_kk,
            n_draws=args.n_draws,
            seed=args.seed + 1009 * idx,
            targets=targets,
            mass_factor=args.mass_factor,
            ckm_factor=args.ckm_factor,
            j_factor=args.j_factor,
            c_Q12=args.c_Q12,
            c_Q3=args.c_Q3,
            rh_down_leak=args.rh_down_leak,
            rh_up_leak=args.rh_up_leak,
            left_23_boost=args.left_23_boost,
            left_13_boost=args.left_13_boost,
            left_cp_boost=args.left_cp_boost,
            perturbativity_mode=args.perturbativity_mode,
        )
        print(f"[rankone_u2] tile M_KK={mkk_tev:.1f} TeV: flat Bauer S1 baseline")
        flat_rows, flat_counts = _run_flat_tile(
            M_KK_GeV=M_KK_GeV,
            xi_KK=args.xi_kk,
            n_draws=args.n_draws,
            seed=args.seed + 500_000 + 1009 * idx,
            targets=targets,
            mass_factor=args.mass_factor,
            ckm_factor=args.ckm_factor,
            j_factor=args.j_factor,
        )
        summary = _summarize_tile(
            mkk_tev=mkk_tev,
            u2_rows=u2_rows,
            flat_rows=flat_rows,
            u2_counts=u2_counts,
            flat_counts=flat_counts,
        )
        summaries.append(summary)
        print(
            "[rankone_u2] "
            f"M_KK={mkk_tev:.1f} TeV "
            f"U2 epsK pass={_fmt(summary['u2_epsK_pass_fraction_pdg'], 'pct')} "
            f"flat epsK pass={_fmt(summary['flat_epsK_pass_fraction_pdg'], 'pct')} "
            f"C4 ratio={_fmt(summary['C4_suppression_u2_over_flat'], 'ratio')} "
            f"GL ratio={_fmt(summary['GL12_suppression_u2_over_flat'], 'ratio')} "
            f"GR ratio={_fmt(summary['GR12_suppression_u2_over_flat'], 'ratio')}"
        )

    report = _build_report(args, summaries)
    out = Path(args.out)
    if not out.is_absolute():
        out = REPO / out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(report)
    print("")
    print(report)
    print(f"[rankone_u2] wrote report -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
