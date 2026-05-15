"""RS flavor anarchy (ACPS) ensemble for the RS5D quark sector.

Literature setup (Agashe-Contino-Pomarol-Sundrum 2004 hep-ph/0408134;
Csaki-Falkowski-Weiler 0804.1954):

    1. **Fix** bulk-mass parameters c_Q, c_u, c_d to a hierarchical pattern.
    2. Draw anarchic O(1) Yukawas Y_u, Y_d.
    3. Build the 4-d quark mass matrices

           M_u = v . diag(f_Q) . Y_u . diag(f_u),
           M_d = v . diag(f_Q) . Y_d . diag(f_d),

       where the f's are the IR overlaps from warpConfig.wavefuncs.f_IR
       evaluated at the **fixed** c-values.
    4. SVD gives quark masses + L/R rotations -> CKM.
    5. Evaluate all five Delta-F=2 systems via the repo's mass-basis
       KK-gluon coupling -> Wilson-coefficient pipeline.

This is a *forward* ensemble: NO optimizer is invoked; NO call to
fit_quark_sector. The ensemble characterizes the natural anarchic
distribution at fixed geometric inputs.

Per-draw row written to a JSONL file; per-tile percentile / pass-rate
summaries written to a tile_summary.json file.

Usage
-----
    python scripts/run_rs_anarchy.py \
        --output-dir scan_outputs/rs_anarchy_<timestamp>/ \
        --n-draws 100000 \
        --m-kk-tev 3,5,7,10,15,20,30,50

    # Smoke test (one tile, light)
    python scripts/run_rs_anarchy.py \
        --output-dir /tmp/rs_anarchy_smoke_dir \
        --n-draws 500 \
        --m-kk-tev 10 \
        --smoke-json /tmp/rs_anarchy_smoke.json
"""

from __future__ import annotations

import argparse
import json
import math
import multiprocessing as mp
import os
import sys
import time
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

import numpy as np

# Make the repo importable.
_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from qcd import alpha_s
from quarkConstraints.deltaf2 import (
    compute_delta_f2_wilsons,
    evaluate_delta_f2_constraints,
    evaluate_delta_mk_with_running,
)
from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.fit import jarlskog_invariant
from warpConfig.baseParams import MPL
from warpConfig.wavefuncs import f_IR


# ---------------------------------------------------------------------------
# Canonical ACPS-like c-value pattern
# ---------------------------------------------------------------------------

# c-values fixed (independent of Y).  These are tuned (see "Verification" notes
# below) so f_IR(c) * v gives roughly SM-quark masses up to O(1) Y factors at
# Lambda_IR ~ 1 TeV.  See the smoke output for explicit f_q values.
DEFAULT_C_Q: Tuple[float, float, float] = (0.63, 0.57, 0.20)
DEFAULT_C_U: Tuple[float, float, float] = (0.66, 0.50, -0.50)
DEFAULT_C_D: Tuple[float, float, float] = (0.66, 0.61, 0.55)

# ACPS / first-KK convention M_KK = xi_KK * Lambda_IR with the gauge KK root.
DEFAULT_XI_KK: float = 2.4487
DEFAULT_K_GEV: float = MPL  # 1.2209e19 GeV
DEFAULT_V_GEV: float = 174.0

# Anarchic Y prior: real and imag parts iid Uniform(-Y_HALF_RANGE, +Y_HALF_RANGE).
DEFAULT_Y_HALF_RANGE: float = 1.5
DEFAULT_Y_FLOOR: float = 0.1  # reject draws with |Y_ij| < Y_FLOOR (literature convention)

# PDG-tolerance defaults for "anarchy works": every mass within factor of 3,
# every CKM element within factor of 3, J within factor of 5. These are the
# fall-back factors applied when the corresponding CLI flag is not given.
PDG_MASS_FACTOR_TOL: float = 3.0
PDG_CKM_FACTOR_TOL: float = 3.0
PDG_J_FACTOR_TOL: float = 5.0
PDG_LOG_MASS_TOL: float = math.log(PDG_MASS_FACTOR_TOL)  # ~1.0986

# Y-prior defaults (Gaussian variant)
DEFAULT_Y_PRIOR: str = "uniform"
DEFAULT_Y_SIGMA: float = 1.0
DEFAULT_Y_TRUNC_SIGMA: float = 3.0


# ---------------------------------------------------------------------------
# PDG targets at mu = m_t(m_t) = 163.5 GeV (matches repo convention)
# ---------------------------------------------------------------------------

def _load_pdg_targets() -> dict:
    """Return PDG target arrays (masses at mu = 163.5 GeV) and CKM observables.

    Lazy import of the repo's pdg_quark_masses helper to keep this script
    independent of the modern-inputs frozen schema.
    """
    from quarkConstraints.pdg_quark_masses import pdg_up_down_arrays_at_scale

    up, down = pdg_up_down_arrays_at_scale(163.5)
    # PDG 2024 |V_ij| central values (Workman et al.)
    abs_V_us = 0.2243
    abs_V_cb = 0.0408
    abs_V_ub = 0.00382
    # Jarlskog invariant J ~ 3.0e-5 (PDG 2024)
    J = 3.0e-5
    return {
        "up_masses_GeV": np.asarray(up, dtype=float),
        "down_masses_GeV": np.asarray(down, dtype=float),
        "abs_V_us": abs_V_us,
        "abs_V_cb": abs_V_cb,
        "abs_V_ub": abs_V_ub,
        "J": J,
    }


# ---------------------------------------------------------------------------
# Anarchic Y sampler
# ---------------------------------------------------------------------------

def _draw_anarchic_matrix(
    rng: np.random.Generator,
    *,
    half_range: float,
    floor: float,
    max_tries: int = 256,
    prior: str = DEFAULT_Y_PRIOR,
    sigma: float = DEFAULT_Y_SIGMA,
    trunc_sigma: float = DEFAULT_Y_TRUNC_SIGMA,
) -> np.ndarray:
    """Return a complex 3x3 Y. Two prior families are supported:

    - ``prior="uniform"`` (default): Re,Im iid in U(-half_range, +half_range)
      with entrywise |Y_ij| >= floor enforced by rejection sampling. This is
      the historical behaviour and is bit-for-bit reproducible.
    - ``prior="gaussian"``: Re,Im iid in N(0, sigma) truncated component-wise
      to |Re|,|Im| <= trunc_sigma * sigma, again with |Y_ij| >= floor.
    """
    if floor < 0.0:
        raise ValueError("require floor >= 0")
    prior = prior.lower()
    if prior == "uniform":
        if half_range <= 0.0:
            raise ValueError("require half_range > 0 for uniform prior")
        for _ in range(max_tries):
            re = rng.uniform(-half_range, half_range, size=(3, 3))
            im = rng.uniform(-half_range, half_range, size=(3, 3))
            Y = re + 1j * im
            if (np.abs(Y) >= floor).all():
                return Y
        # Fall back: keep direction, push tiny entries up to the floor.
        re = rng.uniform(-half_range, half_range, size=(3, 3))
        im = rng.uniform(-half_range, half_range, size=(3, 3))
        Y = re + 1j * im
        m = np.abs(Y)
        phase = Y / np.where(m > 0.0, m, 1.0)
        return phase * np.maximum(m, floor)
    elif prior == "gaussian":
        if sigma <= 0.0:
            raise ValueError("require sigma > 0 for gaussian prior")
        if trunc_sigma <= 0.0:
            raise ValueError("require trunc_sigma > 0 for gaussian prior")
        cap = trunc_sigma * sigma

        def _sample_truncated_normal() -> np.ndarray:
            # Draw a (3,3) array with each entry resampled until within cap.
            arr = rng.normal(loc=0.0, scale=sigma, size=(3, 3))
            for attempt in range(64):
                bad = np.abs(arr) > cap
                if not bad.any():
                    return arr
                arr = np.where(bad, rng.normal(loc=0.0, scale=sigma, size=(3, 3)), arr)
            return np.clip(arr, -cap, cap)

        for _ in range(max_tries):
            re = _sample_truncated_normal()
            im = _sample_truncated_normal()
            Y = re + 1j * im
            if (np.abs(Y) >= floor).all():
                return Y
        # Fall-back: bump magnitudes up to floor while preserving phase.
        re = _sample_truncated_normal()
        im = _sample_truncated_normal()
        Y = re + 1j * im
        m = np.abs(Y)
        phase = Y / np.where(m > 0.0, m, 1.0)
        return phase * np.maximum(m, floor)
    else:
        raise ValueError(f"unknown Y prior: {prior!r}")


# ---------------------------------------------------------------------------
# Tile / config dataclasses
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class TileSpec:
    M_KK_GeV: float
    Lambda_IR_GeV: float
    n_draws: int
    seed: int


@dataclass
class EnsembleConfig:
    mkk_values_GeV: Tuple[float, ...]
    n_draws_per_tile: int
    c_Q: Tuple[float, float, float] = DEFAULT_C_Q
    c_u: Tuple[float, float, float] = DEFAULT_C_U
    c_d: Tuple[float, float, float] = DEFAULT_C_D
    xi_KK: float = DEFAULT_XI_KK
    k_GeV: float = DEFAULT_K_GEV
    v_GeV: float = DEFAULT_V_GEV
    y_half_range: float = DEFAULT_Y_HALF_RANGE
    y_floor: float = DEFAULT_Y_FLOOR
    y_prior: str = DEFAULT_Y_PRIOR
    y_sigma: float = DEFAULT_Y_SIGMA
    y_trunc_sigma: float = DEFAULT_Y_TRUNC_SIGMA
    pdg_mass_factor: float = PDG_MASS_FACTOR_TOL
    pdg_ckm_factor: float = PDG_CKM_FACTOR_TOL
    pdg_j_factor: float = PDG_J_FACTOR_TOL
    base_seed: int = 20260506


# ---------------------------------------------------------------------------
# Forward evaluation
# ---------------------------------------------------------------------------

def _ordered_svd(matrix: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """SVD with singular values ascending (light -> heavy)."""
    U, s, Vh = np.linalg.svd(matrix)
    idx = np.argsort(s)
    return U[:, idx], s[idx], Vh.conj().T[:, idx]


def _hermitian(matrix: np.ndarray) -> np.ndarray:
    arr = np.asarray(matrix, dtype=np.complex128)
    return 0.5 * (arr + arr.conjugate().T)


def _mass_basis_overlap(rotation: np.ndarray, profile_values: np.ndarray) -> np.ndarray:
    profile = np.diag(np.asarray(profile_values, dtype=float) ** 2)
    return _hermitian(rotation.conjugate().T @ profile @ rotation)


def _build_kk_gluon_couplings(
    *,
    M_KK: float,
    xi_KK: float,
    f_Q: np.ndarray,
    f_u: np.ndarray,
    f_d: np.ndarray,
    U_L_u: np.ndarray,
    U_L_d: np.ndarray,
    U_R_u: np.ndarray,
    U_R_d: np.ndarray,
) -> QuarkMassBasisCouplings:
    """Construct mass-basis KK-gluon couplings using the SAME formulas as
    quarkConstraints.couplings.compute_quark_kk_gluon_couplings, but for an
    arbitrary fixed (c_Q, c_u, c_d) state — no QuarkFitResult needed."""
    running_alpha_s = float(alpha_s(M_KK, precision="high"))
    g_s = float(math.sqrt(4.0 * math.pi * running_alpha_s))
    left_overlap = _mass_basis_overlap(U_L_d, f_Q)
    right_down_overlap = _mass_basis_overlap(U_R_d, f_d)
    right_up_overlap = _mass_basis_overlap(U_R_u, f_u)
    left_up_overlap = _mass_basis_overlap(U_L_u, f_Q)
    return QuarkMassBasisCouplings(
        M_KK=float(M_KK),
        xi_KK=float(xi_KK),
        alpha_s=running_alpha_s,
        g_s=g_s,
        left_overlap=left_overlap,
        right_up_overlap=right_up_overlap,
        right_down_overlap=right_down_overlap,
        left_up=_hermitian(g_s * left_up_overlap),
        left_down=_hermitian(g_s * left_overlap),
        right_up=_hermitian(g_s * right_up_overlap),
        right_down=_hermitian(g_s * right_down_overlap),
    )


def _check_pdg_match(
    masses_up: np.ndarray,
    masses_down: np.ndarray,
    abs_V_us: float,
    abs_V_cb: float,
    abs_V_ub: float,
    J: float,
    targets: dict,
    *,
    mass_factor: float = PDG_MASS_FACTOR_TOL,
    ckm_factor: float = PDG_CKM_FACTOR_TOL,
    j_factor: float = PDG_J_FACTOR_TOL,
) -> Tuple[bool, dict]:
    """Apply per-observable PDG-tightness gates.

    Defaults reproduce the historical factor-of-3 mass + CKM tolerance and
    factor-of-5 J tolerance. Each tolerance is independently configurable.
    """
    up_log = np.log(np.maximum(masses_up, 1e-30) / targets["up_masses_GeV"])
    dn_log = np.log(np.maximum(masses_down, 1e-30) / targets["down_masses_GeV"])
    ckm_obs = np.array([
        abs(math.log(max(abs_V_us, 1e-30) / targets["abs_V_us"])),
        abs(math.log(max(abs_V_cb, 1e-30) / targets["abs_V_cb"])),
        abs(math.log(max(abs_V_ub, 1e-30) / targets["abs_V_ub"])),
    ])
    j_log = abs(math.log(max(abs(J), 1e-30) / abs(targets["J"])))
    log_mass_tol = math.log(mass_factor)
    log_ckm_tol = math.log(ckm_factor)
    log_j_tol = math.log(j_factor)
    passes_masses = (np.abs(up_log) <= log_mass_tol).all() and (
        np.abs(dn_log) <= log_mass_tol
    ).all()
    passes_ckm = (ckm_obs <= log_ckm_tol).all()
    passes_j = j_log <= log_j_tol
    return (
        bool(passes_masses and passes_ckm and passes_j),
        {
            "up_log": up_log,
            "down_log": dn_log,
            "ckm_log": ckm_obs,
            "j_log": j_log,
        },
    )


def _evaluate_one_draw(
    rng: np.random.Generator,
    *,
    M_KK_GeV: float,
    Lambda_IR_GeV: float,
    epsilon: float,
    f_Q: np.ndarray,
    f_u: np.ndarray,
    f_d: np.ndarray,
    cfg: EnsembleConfig,
    targets: dict,
) -> dict:
    Y_u = _draw_anarchic_matrix(
        rng,
        half_range=cfg.y_half_range,
        floor=cfg.y_floor,
        prior=cfg.y_prior,
        sigma=cfg.y_sigma,
        trunc_sigma=cfg.y_trunc_sigma,
    )
    Y_d = _draw_anarchic_matrix(
        rng,
        half_range=cfg.y_half_range,
        floor=cfg.y_floor,
        prior=cfg.y_prior,
        sigma=cfg.y_sigma,
        trunc_sigma=cfg.y_trunc_sigma,
    )

    out: dict = {
        "M_KK_GeV": float(M_KK_GeV),
        "Lambda_IR_GeV": float(Lambda_IR_GeV),
        "epsilon": float(epsilon),
        "xi_KK": float(cfg.xi_KK),
        "y_half_range": float(cfg.y_half_range),
        "y_floor": float(cfg.y_floor),
        "y_prior": str(cfg.y_prior),
        "ok": False,
    }

    try:
        D_Q = np.diag(f_Q)
        D_u = np.diag(f_u)
        D_d = np.diag(f_d)
        M_u = cfg.v_GeV * D_Q @ Y_u @ D_u
        M_d = cfg.v_GeV * D_Q @ Y_d @ D_d
        U_L_u, masses_up, U_R_u = _ordered_svd(M_u)
        U_L_d, masses_down, U_R_d = _ordered_svd(M_d)
        ckm = U_L_u.conj().T @ U_L_d
    except (ValueError, np.linalg.LinAlgError) as err:
        out["error"] = f"{type(err).__name__}: {err}"
        return out

    abs_V_us = float(abs(ckm[0, 1]))
    abs_V_cb = float(abs(ckm[1, 2]))
    abs_V_ub = float(abs(ckm[0, 2]))
    J = float(jarlskog_invariant(ckm))

    passes_pdg, residuals = _check_pdg_match(
        masses_up, masses_down, abs_V_us, abs_V_cb, abs_V_ub, J, targets,
        mass_factor=cfg.pdg_mass_factor,
        ckm_factor=cfg.pdg_ckm_factor,
        j_factor=cfg.pdg_j_factor,
    )

    out.update(
        {
            "ok": True,
            "passes_pdg": passes_pdg,
            "masses_up_GeV": [float(x) for x in masses_up],
            "masses_down_GeV": [float(x) for x in masses_down],
            "abs_V_us": abs_V_us,
            "abs_V_cb": abs_V_cb,
            "abs_V_ub": abs_V_ub,
            "J": J,
            "up_log_max": float(np.max(np.abs(residuals["up_log"]))),
            "down_log_max": float(np.max(np.abs(residuals["down_log"]))),
            "ckm_log_max": float(np.max(residuals["ckm_log"])),
            "j_log": float(residuals["j_log"]),
        }
    )

    # Forward-eval Delta-F=2 for *every* draw (PDG-pass independent).
    try:
        couplings = _build_kk_gluon_couplings(
            M_KK=M_KK_GeV,
            xi_KK=cfg.xi_KK,
            f_Q=f_Q,
            f_u=f_u,
            f_d=f_d,
            U_L_u=U_L_u,
            U_L_d=U_L_d,
            U_R_u=U_R_u,
            U_R_d=U_R_d,
        )
        df2 = evaluate_delta_f2_constraints(
            couplings, M_KK=M_KK_GeV, xi_KK=cfg.xi_KK
        )
        # Delta_m_K from the *unevolved* kaon Wilsons run down to mu_had.
        unevolved_wilsons = compute_delta_f2_wilsons(
            couplings, M_KK=M_KK_GeV, xi_KK=cfg.xi_KK
        )
        wilsons_eps_k = next(
            (w for w in unevolved_wilsons if w.input.key == "epsilon_k"), None
        )
        if wilsons_eps_k is not None:
            dmk = evaluate_delta_mk_with_running(wilsons_eps_k, mu_had=2.0)
            ratio_dm_K = float(dmk.ratio_to_exp)
            passes_dm_K = bool(dmk.passes)
        else:
            ratio_dm_K = float("nan")
            passes_dm_K = False
        by_system = df2.by_system
        ratios = {
            "epsilon_K": float(by_system["K"].ratio_to_bound),
            "Delta_m_K": ratio_dm_K,
            "Delta_m_Bd": float(by_system["B_d"].ratio_to_bound),
            "Delta_m_Bs": float(by_system["B_s"].ratio_to_bound),
            "Delta_m_D0": float(by_system["D"].ratio_to_bound),
        }
        passes = {
            "epsilon_K": bool(by_system["K"].passes),
            "Delta_m_K": passes_dm_K,
            "Delta_m_Bd": bool(by_system["B_d"].passes),
            "Delta_m_Bs": bool(by_system["B_s"].passes),
            "Delta_m_D0": bool(by_system["D"].passes),
        }
        max_ratio_5sys = float(max(v for v in ratios.values() if not math.isnan(v)))
        out.update(
            {
                "deltaf2_ratios": ratios,
                "deltaf2_passes": passes,
                "deltaf2_pass_all": bool(all(passes.values())),
                "max_ratio": max_ratio_5sys,
            }
        )
    except Exception as err:  # noqa: BLE001 — robust per-draw isolation
        out["deltaf2_error"] = f"{type(err).__name__}: {err}"

    return out


# ---------------------------------------------------------------------------
# Tile orchestration
# ---------------------------------------------------------------------------

def _run_tile(
    tile: TileSpec,
    cfg: EnsembleConfig,
    targets: dict,
) -> Tuple[List[dict], dict]:
    rng = np.random.default_rng(tile.seed)
    epsilon = tile.Lambda_IR_GeV / cfg.k_GeV
    f_Q = f_IR(np.asarray(cfg.c_Q), epsilon)
    f_u = f_IR(np.asarray(cfg.c_u), epsilon)
    f_d = f_IR(np.asarray(cfg.c_d), epsilon)

    rows: List[dict] = []
    for k in range(tile.n_draws):
        row = _evaluate_one_draw(
            rng,
            M_KK_GeV=tile.M_KK_GeV,
            Lambda_IR_GeV=tile.Lambda_IR_GeV,
            epsilon=epsilon,
            f_Q=f_Q,
            f_u=f_u,
            f_d=f_d,
            cfg=cfg,
            targets=targets,
        )
        row["sample_idx"] = k
        rows.append(row)

    # Tile summary: percentiles + pass rates + per-system fractions.
    ok_rows = [r for r in rows if r.get("ok") and "deltaf2_ratios" in r]
    pdg_rows = [r for r in ok_rows if r.get("passes_pdg")]
    max_ratios = np.array([r["max_ratio"] for r in ok_rows], dtype=float)
    pdg_max_ratios = np.array(
        [r["max_ratio"] for r in pdg_rows], dtype=float
    )
    summary = {
        "M_KK_GeV": tile.M_KK_GeV,
        "Lambda_IR_GeV": tile.Lambda_IR_GeV,
        "epsilon": epsilon,
        "n_draws": tile.n_draws,
        "n_ok": len(ok_rows),
        "n_pdg_pass": len(pdg_rows),
        "pdg_pass_fraction": (len(pdg_rows) / max(1, len(ok_rows))),
        "f_Q": [float(x) for x in f_Q],
        "f_u": [float(x) for x in f_u],
        "f_d": [float(x) for x in f_d],
    }
    if max_ratios.size > 0:
        summary["max_ratio_percentiles_all"] = {
            "p05": float(np.percentile(max_ratios, 5)),
            "p25": float(np.percentile(max_ratios, 25)),
            "p50": float(np.percentile(max_ratios, 50)),
            "p75": float(np.percentile(max_ratios, 75)),
            "p95": float(np.percentile(max_ratios, 95)),
        }
    if pdg_max_ratios.size > 0:
        summary["max_ratio_percentiles_pdg"] = {
            "p05": float(np.percentile(pdg_max_ratios, 5)),
            "p25": float(np.percentile(pdg_max_ratios, 25)),
            "p50": float(np.percentile(pdg_max_ratios, 50)),
            "p75": float(np.percentile(pdg_max_ratios, 75)),
            "p95": float(np.percentile(pdg_max_ratios, 95)),
        }
    # Per-system pass fractions (within PDG-passing subset)
    sys_keys = ("epsilon_K", "Delta_m_K", "Delta_m_Bd", "Delta_m_Bs", "Delta_m_D0")
    if pdg_rows:
        per_sys_pass = {
            k: float(
                sum(1 for r in pdg_rows if r["deltaf2_passes"].get(k, False))
                / len(pdg_rows)
            )
            for k in sys_keys
        }
        summary["per_system_pass_fraction_pdg"] = per_sys_pass
    if ok_rows:
        per_sys_pass_all = {
            k: float(
                sum(1 for r in ok_rows if r["deltaf2_passes"].get(k, False))
                / len(ok_rows)
            )
            for k in sys_keys
        }
        summary["per_system_pass_fraction_all"] = per_sys_pass_all

    return rows, summary


# ---------------------------------------------------------------------------
# Multiprocessing worker
# ---------------------------------------------------------------------------

_GLOBAL_CFG: Optional[EnsembleConfig] = None
_GLOBAL_TARGETS: Optional[dict] = None


def _worker_init(cfg_dict: dict) -> None:
    global _GLOBAL_CFG, _GLOBAL_TARGETS
    cfg = EnsembleConfig(
        mkk_values_GeV=tuple(cfg_dict["mkk_values_GeV"]),
        n_draws_per_tile=int(cfg_dict["n_draws_per_tile"]),
        c_Q=tuple(cfg_dict["c_Q"]),
        c_u=tuple(cfg_dict["c_u"]),
        c_d=tuple(cfg_dict["c_d"]),
        xi_KK=float(cfg_dict["xi_KK"]),
        k_GeV=float(cfg_dict["k_GeV"]),
        v_GeV=float(cfg_dict["v_GeV"]),
        y_half_range=float(cfg_dict["y_half_range"]),
        y_floor=float(cfg_dict["y_floor"]),
        y_prior=str(cfg_dict.get("y_prior", DEFAULT_Y_PRIOR)),
        y_sigma=float(cfg_dict.get("y_sigma", DEFAULT_Y_SIGMA)),
        y_trunc_sigma=float(cfg_dict.get("y_trunc_sigma", DEFAULT_Y_TRUNC_SIGMA)),
        pdg_mass_factor=float(cfg_dict.get("pdg_mass_factor", PDG_MASS_FACTOR_TOL)),
        pdg_ckm_factor=float(cfg_dict.get("pdg_ckm_factor", PDG_CKM_FACTOR_TOL)),
        pdg_j_factor=float(cfg_dict.get("pdg_j_factor", PDG_J_FACTOR_TOL)),
        base_seed=int(cfg_dict["base_seed"]),
    )
    _GLOBAL_CFG = cfg
    _GLOBAL_TARGETS = _load_pdg_targets()


def _worker_run_tile(tile_payload: dict) -> Tuple[List[dict], dict]:
    assert _GLOBAL_CFG is not None and _GLOBAL_TARGETS is not None
    tile = TileSpec(
        M_KK_GeV=float(tile_payload["M_KK_GeV"]),
        Lambda_IR_GeV=float(tile_payload["Lambda_IR_GeV"]),
        n_draws=int(tile_payload["n_draws"]),
        seed=int(tile_payload["seed"]),
    )
    return _run_tile(tile, _GLOBAL_CFG, _GLOBAL_TARGETS)


# ---------------------------------------------------------------------------
# CLI / driver
# ---------------------------------------------------------------------------

def _parse_csv_floats(text: str) -> Tuple[float, ...]:
    return tuple(float(t) for t in text.split(",") if t.strip())


def _build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="RS flavor anarchy (ACPS) ensemble.")
    p.add_argument("--output-dir", type=str, required=True)
    p.add_argument("--n-draws", type=int, default=100_000)
    p.add_argument(
        "--m-kk-tev",
        type=str,
        default="3,5,7,10,15,20,30,50",
        help="Comma-separated M_KK values in TeV.",
    )
    p.add_argument("--xi-kk", type=float, default=DEFAULT_XI_KK)
    p.add_argument("--k-gev", type=float, default=DEFAULT_K_GEV)
    p.add_argument("--v-gev", type=float, default=DEFAULT_V_GEV)
    p.add_argument("--y-half-range", type=float, default=DEFAULT_Y_HALF_RANGE)
    p.add_argument("--y-floor", type=float, default=DEFAULT_Y_FLOOR)
    p.add_argument("--base-seed", type=int, default=20260506)
    p.add_argument("--n-workers", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", "1")))
    p.add_argument("--smoke-json", type=str, default=None,
                   help="If set, write a smoke summary including sample rows.")
    p.add_argument("--c-Q", type=str, default=",".join(str(x) for x in DEFAULT_C_Q))
    p.add_argument("--c-u", type=str, default=",".join(str(x) for x in DEFAULT_C_U))
    p.add_argument("--c-d", type=str, default=",".join(str(x) for x in DEFAULT_C_D))
    # Y prior
    p.add_argument("--y-prior", type=str, default=DEFAULT_Y_PRIOR,
                   choices=["uniform", "gaussian"],
                   help="Y prior family. 'uniform' is the historical default and "
                        "is bit-for-bit reproducible.")
    p.add_argument("--y-sigma", type=float, default=DEFAULT_Y_SIGMA,
                   help="Standard deviation of Re/Im for the gaussian Y prior.")
    p.add_argument("--y-trunc-sigma", type=float, default=DEFAULT_Y_TRUNC_SIGMA,
                   help="Truncation in units of sigma for the gaussian Y prior.")
    # PDG-tightness gates
    p.add_argument("--pdg-mass-factor", type=float, default=PDG_MASS_FACTOR_TOL,
                   help="Allowed multiplicative deviation in quark masses.")
    p.add_argument("--pdg-ckm-factor", type=float, default=PDG_CKM_FACTOR_TOL,
                   help="Allowed multiplicative deviation in |V_us|, |V_cb|, |V_ub|.")
    p.add_argument("--pdg-j-factor", type=float, default=PDG_J_FACTOR_TOL,
                   help="Allowed multiplicative deviation in the Jarlskog invariant J.")
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _build_argparser().parse_args(argv)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    mkk_values_GeV = tuple(t * 1000.0 for t in _parse_csv_floats(args.m_kk_tev))
    cfg = EnsembleConfig(
        mkk_values_GeV=mkk_values_GeV,
        n_draws_per_tile=int(args.n_draws),
        c_Q=_parse_csv_floats(args.c_Q),
        c_u=_parse_csv_floats(args.c_u),
        c_d=_parse_csv_floats(args.c_d),
        xi_KK=args.xi_kk,
        k_GeV=args.k_gev,
        v_GeV=args.v_gev,
        y_half_range=args.y_half_range,
        y_floor=args.y_floor,
        y_prior=args.y_prior,
        y_sigma=args.y_sigma,
        y_trunc_sigma=args.y_trunc_sigma,
        pdg_mass_factor=args.pdg_mass_factor,
        pdg_ckm_factor=args.pdg_ckm_factor,
        pdg_j_factor=args.pdg_j_factor,
        base_seed=args.base_seed,
    )

    print(f"[rs_anarchy] output_dir = {output_dir}")
    print(f"[rs_anarchy] M_KK (TeV)  = {[m / 1000.0 for m in mkk_values_GeV]}")
    print(f"[rs_anarchy] n_draws/tile = {cfg.n_draws_per_tile}")
    print(f"[rs_anarchy] xi_KK       = {cfg.xi_KK}")
    print(f"[rs_anarchy] c_Q         = {cfg.c_Q}")
    print(f"[rs_anarchy] c_u         = {cfg.c_u}")
    print(f"[rs_anarchy] c_d         = {cfg.c_d}")
    if cfg.y_prior == "uniform":
        print(f"[rs_anarchy] Y prior: uniform; |Y| in [{cfg.y_floor}, sqrt(2)*{cfg.y_half_range}]; iid Re/Im")
    else:
        print(
            f"[rs_anarchy] Y prior: gaussian; sigma={cfg.y_sigma}, "
            f"trunc={cfg.y_trunc_sigma} sigma, |Y_ij| floor={cfg.y_floor}"
        )
    print(
        f"[rs_anarchy] PDG gates: mass<x{cfg.pdg_mass_factor}, "
        f"CKM<x{cfg.pdg_ckm_factor}, J<x{cfg.pdg_j_factor}"
    )
    print(f"[rs_anarchy] n_workers   = {args.n_workers}")

    # Build tile payloads.
    tiles: List[dict] = []
    for idx, M_KK in enumerate(mkk_values_GeV):
        Lambda_IR = M_KK / cfg.xi_KK
        tiles.append(
            {
                "M_KK_GeV": float(M_KK),
                "Lambda_IR_GeV": float(Lambda_IR),
                "n_draws": cfg.n_draws_per_tile,
                "seed": int(cfg.base_seed + 1009 * idx),
            }
        )

    rows_path = output_dir / "draws.jsonl"
    summary_path = output_dir / "tile_summary.json"
    targets = _load_pdg_targets()

    cfg_dict = {
        "mkk_values_GeV": list(cfg.mkk_values_GeV),
        "n_draws_per_tile": cfg.n_draws_per_tile,
        "c_Q": list(cfg.c_Q),
        "c_u": list(cfg.c_u),
        "c_d": list(cfg.c_d),
        "xi_KK": cfg.xi_KK,
        "k_GeV": cfg.k_GeV,
        "v_GeV": cfg.v_GeV,
        "y_half_range": cfg.y_half_range,
        "y_floor": cfg.y_floor,
        "y_prior": cfg.y_prior,
        "y_sigma": cfg.y_sigma,
        "y_trunc_sigma": cfg.y_trunc_sigma,
        "pdg_mass_factor": cfg.pdg_mass_factor,
        "pdg_ckm_factor": cfg.pdg_ckm_factor,
        "pdg_j_factor": cfg.pdg_j_factor,
        "base_seed": cfg.base_seed,
    }

    summaries: List[dict] = []
    t_start = time.time()
    if args.n_workers <= 1:
        # Inline single-process path.
        _worker_init(cfg_dict)
        with rows_path.open("w") as fh:
            for tile_payload in tiles:
                rows, tile_summary = _worker_run_tile(tile_payload)
                for r in rows:
                    fh.write(json.dumps(r) + "\n")
                summaries.append(tile_summary)
                print(
                    f"[rs_anarchy] tile M_KK={tile_payload['M_KK_GeV']/1000:.2f} TeV  "
                    f"n_pdg_pass={tile_summary['n_pdg_pass']}/{tile_summary['n_ok']}  "
                    f"PDG_frac={tile_summary['pdg_pass_fraction']:.3%}"
                )
    else:
        with mp.Pool(
            args.n_workers, initializer=_worker_init, initargs=(cfg_dict,)
        ) as pool, rows_path.open("w") as fh:
            for rows, tile_summary in pool.imap_unordered(
                _worker_run_tile, tiles
            ):
                for r in rows:
                    fh.write(json.dumps(r) + "\n")
                summaries.append(tile_summary)
                print(
                    f"[rs_anarchy] tile M_KK={tile_summary['M_KK_GeV']/1000:.2f} TeV  "
                    f"n_pdg_pass={tile_summary['n_pdg_pass']}/{tile_summary['n_ok']}  "
                    f"PDG_frac={tile_summary['pdg_pass_fraction']:.3%}"
                )
    elapsed = time.time() - t_start

    if cfg.y_prior == "uniform":
        y_prior_doc = (
            "Re,Im iid Uniform(-y_half_range, +y_half_range); "
            "rejection sampler enforces |Y_ij| >= y_floor"
        )
    else:
        y_prior_doc = (
            f"Re,Im iid Normal(0, sigma={cfg.y_sigma}) truncated component-wise "
            f"to |Re|,|Im| <= {cfg.y_trunc_sigma}*sigma; rejection sampler "
            "enforces |Y_ij| >= y_floor"
        )
    summary_payload = {
        "config": cfg_dict,
        "elapsed_seconds": elapsed,
        "tiles": sorted(summaries, key=lambda s: s["M_KK_GeV"]),
        "schema": "rs_anarchy_acps_v1",
        "convention": {
            "f_factor": "warpConfig.wavefuncs.f_IR (IR-localized formula)",
            "M_KK_to_Lambda_IR": "Lambda_IR = M_KK / xi_KK",
            "default_xi_KK": cfg.xi_KK,
            "Y_prior": y_prior_doc,
            "y_prior_kind": cfg.y_prior,
            "y_half_range": cfg.y_half_range,
            "y_floor": cfg.y_floor,
            "y_sigma": cfg.y_sigma,
            "y_trunc_sigma": cfg.y_trunc_sigma,
            "pdg_tolerance": {
                "mass_factor": cfg.pdg_mass_factor,
                "ckm_factor": cfg.pdg_ckm_factor,
                "j_factor": cfg.pdg_j_factor,
            },
        },
    }
    with summary_path.open("w") as fh:
        json.dump(summary_payload, fh, indent=2, default=float)

    print(f"[rs_anarchy] Done in {elapsed:.1f}s")
    print(f"[rs_anarchy] rows -> {rows_path}")
    print(f"[rs_anarchy] summary -> {summary_path}")

    # Optional smoke output
    if args.smoke_json:
        # Stream rows back from the JSONL to assemble sample_pass / sample_fail.
        sample_pass: Optional[dict] = None
        sample_fail: Optional[dict] = None
        with rows_path.open() as fh:
            for line in fh:
                try:
                    r = json.loads(line)
                except json.JSONDecodeError:
                    continue
                if not r.get("ok"):
                    continue
                if sample_pass is None and r.get("passes_pdg"):
                    sample_pass = r
                if sample_fail is None and not r.get("passes_pdg"):
                    sample_fail = r
                if sample_pass is not None and sample_fail is not None:
                    break
        smoke_payload = {
            "config": cfg_dict,
            "tile_summary": summary_payload["tiles"],
            "sample_passing_draw": sample_pass,
            "sample_failing_draw": sample_fail,
        }
        with open(args.smoke_json, "w") as fh:
            json.dump(smoke_payload, fh, indent=2, default=float)
        print(f"[rs_anarchy] smoke -> {args.smoke_json}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
