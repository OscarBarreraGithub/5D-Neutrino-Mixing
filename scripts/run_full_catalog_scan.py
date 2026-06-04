#!/usr/bin/env python3
"""Full flavor-catalog RS scan harness.

This driver is intentionally outside ``flavor_catalog_constraints``: it is a
scan harness, not constraint physics.  It builds one RS-EW spectrum/spline cache
per tile, injects that cache into the existing point builder, evaluates the full
constraint registry, and streams one compact JSON row per draw.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import multiprocessing as mp
import os
import subprocess
import sys
import time
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from flavor_catalog_constraints import point_builder, registry
from flavor_catalog_constraints.base import ConstraintResult, Severity
from quarkConstraints.benchmarks import default_quark_targets
from quarkConstraints.couplings import compute_quark_kk_gluon_couplings
from quarkConstraints.fit import QuarkFitResult, QuarkFitSeed, fit_quark_sector
from quarkConstraints.model import RotationParameters
from quarkConstraints.rs_ew_spectrum import (
    DEFAULT_MAX_TRUNCATION_MODES,
    DEFAULT_MIN_TRUNCATION_MODES,
    DEFAULT_N_GAUGE_MODES,
    DEFAULT_OVERLAP_RTOL,
    DEFAULT_QUADRATURE_ORDER,
    RSEWOverlapSplineCache,
    RSEWSpectrum,
)
from warpConfig.baseParams import MPL
from warpConfig.wavefuncs import f_IR
from yukawa.compute_yukawas import YukawaResult, compute_all_yukawas


EXPECTED_REGISTRY_COUNT = 103
QUARK_ONLY_CANDIDATE_IDS = (
    "K001", "K002", "B001", "B002", "B003", "B004", "C001", "C002",
    "B011", "B012", "B013", "B014",
    "B032", "B033", "B034", "C003", "K003", "K013",
    "T001", "T002", "T003", "T004", "T005", "T006", "T007", "T008",
    "T010", "T011", "T012", "T014",
    "EW001", "EW002", "EW003",
    "E004", "E006", "E007", "E008", "E009",
)
QUARK_ONLY_DEFERRED_LEPTON_FOLLOWUP = ("EW002",)
QUARK_ONLY_OPTIONAL_LEPTON_DIAGNOSTIC_IDS = ("EW003",)
QUARK_ONLY_ALLOWLIST_IDS = tuple(
    pid for pid in QUARK_ONLY_CANDIDATE_IDS if pid not in QUARK_ONLY_DEFERRED_LEPTON_FOLLOWUP
)
QUARK_ONLY_ALLOWLIST_SET = frozenset(QUARK_ONLY_ALLOWLIST_IDS)
QUARK_ONLY_LEPTON_SECTOR_LABEL = "dropped (not rigorous)"
QUARK_ONLY_DEFERRED_SCOPE_TAG = "deferred_lepton_followup"
QUARK_ONLY_FORBIDDEN_EXTRAS = frozenset(
    {
        "lepton_mass_basis_couplings",
        "rs_charged_current",
        "rs_higgs_yukawas",
    }
)
QUARK_ONLY_BUILD_INCLUDE_FLAGS = {
    "include_charged_current": False,
    "include_fermion_kk_mixing": True,
    "include_higgs_yukawas": False,
}
QUARK_ONLY_ALLOWLIST_EXTRAS: dict[str, tuple[str, ...]] = {
    "K001": ("quark_mass_basis_couplings",),
    "K002": ("quark_mass_basis_couplings",),
    "B001": ("quark_mass_basis_couplings",),
    "B002": ("quark_mass_basis_couplings",),
    "B003": ("quark_mass_basis_couplings",),
    "B004": ("quark_mass_basis_couplings",),
    "C001": ("quark_mass_basis_couplings",),
    "C002": ("quark_mass_basis_couplings",),
    "B011": ("quark_mass_basis_couplings", "kk_ew_mass_gev"),
    "B012": ("quark_mass_basis_couplings", "kk_ew_mass_gev"),
    "B013": ("quark_mass_basis_couplings", "kk_ew_mass_gev"),
    "B014": ("quark_mass_basis_couplings", "kk_ew_mass_gev"),
    "B032": (),
    "B033": (),
    "B034": (),
    "C003": (),
    "K003": (),
    "K013": (),
    "T001": ("quark_mass_basis_couplings", "kk_ew_mass_gev"),
    "T002": ("quark_mass_basis_couplings", "kk_ew_mass_gev"),
    "T003": ("quark_mass_basis_couplings", "kk_ew_mass_gev"),
    "T004": ("quark_mass_basis_couplings", "kk_ew_mass_gev"),
    "T005": ("quark_mass_basis_couplings", "kk_gluon_mass_gev"),
    "T006": ("quark_mass_basis_couplings", "kk_gluon_mass_gev"),
    "T007": ("quark_mass_basis_couplings", "kk_ew_mass_gev"),
    "T008": ("quark_mass_basis_couplings", "kk_ew_mass_gev"),
    "T010": ("rs_ew_couplings", "quark_mass_basis_couplings"),
    "T011": ("rs_ew_couplings", "quark_mass_basis_couplings"),
    "T012": ("rs_ew_couplings", "quark_mass_basis_couplings"),
    "T014": ("rs_ew_couplings", "quark_mass_basis_couplings"),
    "EW001": ("kk_ew_mass_gev", "kk_gluon_mass_gev", "quark_mass_basis_couplings"),
    "EW003": ("rs_charged_current",),
    "E004": (),
    "E006": (),
    "E007": (),
    "E008": (),
    "E009": (),
}
QUARK_ONLY_DEFERRED_EXTRAS: dict[str, tuple[str, ...]] = {
    "EW002": ("rs_charged_current",),
}
DEFAULT_XI_KK = 2.4487
DEFAULT_V_GEV = 174.0
DEFAULT_BASE_SEED = 20260604
DEFAULT_TILE_SEED_STRIDE = 1_000_003
DEFAULT_Y_HALF_RANGE = 1.5
DEFAULT_Y_FLOOR = 0.1
DEFAULT_Y_PRIOR = "uniform"
DEFAULT_Y_SIGMA = 1.0
DEFAULT_Y_TRUNC_SIGMA = 3.0
DEFAULT_C_MIN = 0.3
DEFAULT_C_MAX = 0.9
DEFAULT_M_N_MIN_GEV = 1.0e10
DEFAULT_M_N_MAX_GEV = 1.0e15
DEFAULT_LIGHTEST_NU_MIN_EV = 0.0
DEFAULT_LIGHTEST_NU_MAX_EV = 1.0e-2
DEFAULT_PERTURBATIVE_YBAR_MAX = 4.0
DEFAULT_QUARK_Y_ABS_MAX = 4.0
DEFAULT_QUARK_FIT_R = 0.25
DEFAULT_QUARK_FIT_MAX_NFEV = 120
DEFAULT_MAX_FIT_SCORE = 1.0e-2
DEFAULT_C_HALF_ATOL = 1.0e-10
DEFAULT_SPLINE_GRID_SIZE = 241
DEFAULT_SPLINE_VERIFY_POINTS = 41


@dataclass(frozen=True)
class ScanConfig:
    mkk_values_gev: tuple[float, ...]
    n_draws_per_tile: int
    quark_only: bool = False
    xi_kk: float = DEFAULT_XI_KK
    k_gev: float = MPL
    v_gev: float = DEFAULT_V_GEV
    base_seed: int = DEFAULT_BASE_SEED
    tile_seed_stride: int = DEFAULT_TILE_SEED_STRIDE
    y_half_range: float = DEFAULT_Y_HALF_RANGE
    y_floor: float = DEFAULT_Y_FLOOR
    y_prior: str = DEFAULT_Y_PRIOR
    y_sigma: float = DEFAULT_Y_SIGMA
    y_trunc_sigma: float = DEFAULT_Y_TRUNC_SIGMA
    c_min: float = DEFAULT_C_MIN
    c_max: float = DEFAULT_C_MAX
    m_n_min_gev: float = DEFAULT_M_N_MIN_GEV
    m_n_max_gev: float = DEFAULT_M_N_MAX_GEV
    lightest_nu_min_ev: float = DEFAULT_LIGHTEST_NU_MIN_EV
    lightest_nu_max_ev: float = DEFAULT_LIGHTEST_NU_MAX_EV
    normal_ordering_probability: float = 0.5
    perturbative_ybar_max: float = DEFAULT_PERTURBATIVE_YBAR_MAX
    quark_y_abs_max: float = DEFAULT_QUARK_Y_ABS_MAX
    quark_fit_r: float = DEFAULT_QUARK_FIT_R
    quark_fit_max_nfev: int = DEFAULT_QUARK_FIT_MAX_NFEV
    max_fit_score: float = DEFAULT_MAX_FIT_SCORE
    require_fit_success: bool = True
    n_gauge_modes: int = DEFAULT_N_GAUGE_MODES
    quadrature_order: int = DEFAULT_QUADRATURE_ORDER
    min_overlap_modes: int = DEFAULT_MIN_TRUNCATION_MODES
    max_overlap_modes: int = DEFAULT_MAX_TRUNCATION_MODES
    overlap_rel_tol: float = DEFAULT_OVERLAP_RTOL
    spline_grid_size: int = DEFAULT_SPLINE_GRID_SIZE
    spline_verify_points: int = DEFAULT_SPLINE_VERIFY_POINTS
    c_half_atol: float = DEFAULT_C_HALF_ATOL
    expected_registry_count: int = EXPECTED_REGISTRY_COUNT
    sanity_mkk_gev: float = 50_000.0
    sanity_c: float = 0.4


@dataclass(frozen=True)
class TileSpec:
    tile_id: int
    mkk_gev: float
    lambda_ir_gev: float
    n_draws: int
    seed: int


_GLOBAL_CFG: ScanConfig | None = None
_GLOBAL_PROVENANCE: dict[str, Any] | None = None
_GLOBAL_REGISTRY_COUNT: int | None = None


def _worker_init(cfg_payload: Mapping[str, Any], provenance: Mapping[str, Any]) -> None:
    global _GLOBAL_CFG, _GLOBAL_PROVENANCE, _GLOBAL_REGISTRY_COUNT
    _GLOBAL_CFG = _config_from_payload(cfg_payload)
    _GLOBAL_PROVENANCE = dict(provenance)
    registry.discover()
    failures = registry.import_failures()
    count = len(registry.all_constraints())
    if count != _GLOBAL_CFG.expected_registry_count or failures:
        failure_keys = sorted(failures)
        raise RuntimeError(
            "constraint registry discovery failed: "
            f"count={count}, expected={_GLOBAL_CFG.expected_registry_count}, "
            f"import_failures={failure_keys}"
        )
    _GLOBAL_REGISTRY_COUNT = count


def _worker_run_tile(payload: Mapping[str, Any]) -> dict[str, Any]:
    if _GLOBAL_CFG is None or _GLOBAL_PROVENANCE is None or _GLOBAL_REGISTRY_COUNT is None:
        raise RuntimeError("worker was not initialized")
    tile = TileSpec(
        tile_id=int(payload["tile_id"]),
        mkk_gev=float(payload["mkk_gev"]),
        lambda_ir_gev=float(payload["lambda_ir_gev"]),
        n_draws=int(payload["n_draws"]),
        seed=int(payload["seed"]),
    )
    output_dir = Path(str(payload["output_dir"]))
    config_hash = str(payload["config_hash"])
    return _run_tile(
        tile,
        cfg=_GLOBAL_CFG,
        provenance=_GLOBAL_PROVENANCE,
        registry_count=_GLOBAL_REGISTRY_COUNT,
        output_dir=output_dir,
        config_hash=config_hash,
    )


def _run_tile(
    tile: TileSpec,
    *,
    cfg: ScanConfig,
    provenance: Mapping[str, Any],
    registry_count: int,
    output_dir: Path,
    config_hash: str,
) -> dict[str, Any]:
    final_jsonl = output_dir / f"tile-{tile.tile_id:05d}.jsonl"
    final_summary = output_dir / f"tile-{tile.tile_id:05d}.summary.json"
    tmp_jsonl = final_jsonl.with_suffix(final_jsonl.suffix + f".tmp.{os.getpid()}")
    tmp_summary = final_summary.with_suffix(final_summary.suffix + f".tmp.{os.getpid()}")

    tile_started = time.perf_counter()
    cache_started = time.perf_counter()
    spectrum = RSEWSpectrum.build(
        lambda_ir_gev=tile.lambda_ir_gev,
        k_gev=cfg.k_gev,
        n_gauge_modes=cfg.n_gauge_modes,
        quadrature_order=cfg.quadrature_order,
        model_label="minimal_rs",
    )
    overlap_cache = RSEWOverlapSplineCache.build(
        spectrum,
        c_min=cfg.c_min,
        c_max=cfg.c_max,
        grid_size=cfg.spline_grid_size,
        include_omega=True,
        verify_points=cfg.spline_verify_points,
        rel_tol=cfg.overlap_rel_tol,
        min_modes=cfg.min_overlap_modes,
        max_modes=cfg.max_overlap_modes,
    )
    cache_seconds = time.perf_counter() - cache_started

    counters: Counter[str] = Counter()
    hard_veto_rigorous: Counter[str] = Counter()
    hard_veto_proxy: Counter[str] = Counter()
    hard_not_evaluated: Counter[str] = Counter()
    tag_counts: Counter[str] = Counter()
    exception_ids: Counter[str] = Counter()
    constraint_tallies = (
        _empty_constraint_tallies(QUARK_ONLY_ALLOWLIST_IDS) if cfg.quark_only else None
    )
    draw_loop_started = time.perf_counter()

    with tmp_jsonl.open("w", encoding="utf-8") as fh:
        for draw_idx in range(tile.n_draws):
            draw_seed = int(tile.seed + draw_idx)
            draw_rng = np.random.default_rng(draw_seed)
            row = _evaluate_draw(
                draw_rng,
                tile=tile,
                draw_idx=draw_idx,
                draw_seed=draw_seed,
                cfg=cfg,
                spectrum=spectrum,
                overlap_cache=overlap_cache,
                provenance=provenance,
                registry_count=registry_count,
                config_hash=config_hash,
            )
            fh.write(json.dumps(row, separators=(",", ":"), sort_keys=True) + "\n")
            if (draw_idx + 1) % 100 == 0:
                fh.flush()
            _accumulate_row(
                row,
                counters=counters,
                hard_veto_rigorous=hard_veto_rigorous,
                hard_veto_proxy=hard_veto_proxy,
                hard_not_evaluated=hard_not_evaluated,
                tag_counts=tag_counts,
                exception_ids=exception_ids,
                constraint_tallies=constraint_tallies,
            )

    draw_loop_seconds = time.perf_counter() - draw_loop_started
    elapsed = time.perf_counter() - tile_started
    summary = {
        "tile_id": tile.tile_id,
        "complete": True,
        "config_hash": config_hash,
        "seed": tile.seed,
        "n_requested": tile.n_draws,
        "n_rows": int(counters["rows"]),
        "n_evaluated_points": int(counters["evaluated_points"]),
        "n_skipped": int(counters["skipped"]),
        "skip_reasons": _counter_dict_with_prefix(counters, "skip:"),
        "survives_all_HARD_strict": int(counters["survives_strict"]),
        "survives_all_HARD_inclusive": int(counters["survives_inclusive"]),
        "constraint_evaluations": int(counters["constraint_results"]),
        "constraint_evaluated": int(counters["constraint_evaluated"]),
        "constraint_active": int(counters["constraint_active"]),
        "constraint_exceptions": int(counters["constraint_exceptions"]),
        "constraint_exception_rate": _safe_div(
            counters["constraint_exceptions"], counters["constraint_results"]
        ),
        "tag_counts": dict(sorted(tag_counts.items())),
        "hard_vetoes_rigorous": dict(sorted(hard_veto_rigorous.items())),
        "hard_vetoes_proxy": dict(sorted(hard_veto_proxy.items())),
        "hard_not_evaluated": dict(sorted(hard_not_evaluated.items())),
        "exception_ids": dict(sorted(exception_ids.items())),
        "timing": {
            "elapsed_seconds": float(elapsed),
            "cache_build_seconds": float(cache_seconds),
            "draw_loop_seconds": float(draw_loop_seconds),
            "post_cache_seconds_per_draw": _safe_div(draw_loop_seconds, tile.n_draws),
            "post_cache_seconds_per_evaluated_point": _safe_div(
                draw_loop_seconds, counters["evaluated_points"]
            ),
        },
        "cache_metrics": {
            "lambda_ir_gev": float(tile.lambda_ir_gev),
            "mkk_gev": float(tile.mkk_gev),
            "n_gauge_modes": int(spectrum.n_gauge_modes),
            "quadrature_order": int(spectrum.quadrature_order),
            "spline_grid_size": int(overlap_cache.c_grid.size),
            "spline_max_modes": int(overlap_cache.max_modes),
            "spline_max_a_rel_err": float(overlap_cache.max_a_rel_err),
            "spline_verify_points": int(overlap_cache.verification_points),
            "omega_precomputed": overlap_cache.omega_grid is not None,
        },
        "provenance": {**dict(provenance), "registry_count": int(registry_count)},
    }
    if cfg.quark_only:
        summary.update(
            {
                "mode": "quark_only",
                "allowlist": list(QUARK_ONLY_ALLOWLIST_IDS),
                "lepton_sector": QUARK_ONLY_LEPTON_SECTOR_LABEL,
                "deferred_scope_tag": QUARK_ONLY_DEFERRED_SCOPE_TAG,
                "deferred_lepton_followup": list(QUARK_ONLY_DEFERRED_LEPTON_FOLLOWUP),
                "build_include_flags": dict(QUARK_ONLY_BUILD_INCLUDE_FLAGS),
                "constraint_tallies": _finalize_constraint_tallies(
                    constraint_tallies or {}
                ),
            }
        )

    with tmp_summary.open("w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, sort_keys=True)
        fh.write("\n")
    os.replace(tmp_jsonl, final_jsonl)
    os.replace(tmp_summary, final_summary)
    return summary


def _evaluate_draw(
    rng: np.random.Generator,
    *,
    tile: TileSpec,
    draw_idx: int,
    draw_seed: int,
    cfg: ScanConfig,
    spectrum: RSEWSpectrum,
    overlap_cache: RSEWOverlapSplineCache,
    provenance: Mapping[str, Any],
    registry_count: int,
    config_hash: str,
) -> dict[str, Any]:
    params: dict[str, Any] = {
        "Lambda_IR": float(tile.lambda_ir_gev),
        "M_KK": float(tile.mkk_gev),
        "xi_KK": float(cfg.xi_kk),
        "k": float(cfg.k_gev),
    }
    fit_diagnostics: dict[str, Any] = {}
    try:
        y_u_seed, y_d_seed, quark_seed = _draw_quark_seed(rng, cfg)
        params["quark_yukawa_seed"] = {
            "Y_u_re": _matrix_real(y_u_seed.real),
            "Y_u_im": _matrix_real(y_u_seed.imag),
            "Y_d_re": _matrix_real(y_d_seed.real),
            "Y_d_im": _matrix_real(y_d_seed.imag),
        }
        if cfg.quark_only:
            quark_solution = fit_quark_sector(
                default_quark_targets(),
                r=cfg.quark_fit_r,
                seed=quark_seed,
                Lambda_IR=tile.lambda_ir_gev,
                k=cfg.k_gev,
                max_nfev=cfg.quark_fit_max_nfev,
                fit_orientation=True,
            )
            fit_result = quark_solution.result
            fit_diagnostics = _quark_fit_diagnostics(quark_solution)
            _require_valid_quark_fit(fit_result, fit_diagnostics, cfg)
            _require_no_c_half_singularity(
                [
                    *np.asarray(fit_result.bulk_state.c_Q, dtype=float),
                    *np.asarray(fit_result.bulk_state.c_u, dtype=float),
                    *np.asarray(fit_result.bulk_state.c_d, dtype=float),
                ],
                cfg=cfg,
            )
            rs_point = point_builder.build_from_rs_ew_inputs(
                fit_result,
                Lambda_IR=tile.lambda_ir_gev,
                k=cfg.k_gev,
                n_gauge_modes=cfg.n_gauge_modes,
                quadrature_order=cfg.quadrature_order,
                min_overlap_modes=cfg.min_overlap_modes,
                max_overlap_modes=cfg.max_overlap_modes,
                overlap_rel_tol=cfg.overlap_rel_tol,
                spectrum=spectrum,
                rs_ew_cache=overlap_cache,
                **QUARK_ONLY_BUILD_INCLUDE_FLAGS,
                lepton_yukawa_result=None,
                raw={"tile_id": tile.tile_id, "draw_id": draw_idx, "params": params},
            )
            quark_couplings = compute_quark_kk_gluon_couplings(
                fit_result,
                M_KK=tile.mkk_gev,
                xi_KK=cfg.xi_kk,
                g_s_star=None,
            )
            point = point_builder.make_point(
                raw=rs_point.raw,
                **dict(rs_point.extras),
                quark_mass_basis_couplings=quark_couplings,
                kk_gluon_mass_gev=float(tile.mkk_gev),
            )
            results = _evaluate_constraint_ids(point, QUARK_ONLY_ALLOWLIST_IDS)
            payload = _classify_results(results)
            return {
                "draw_id": int(draw_idx),
                "tile_id": int(tile.tile_id),
                "seed": int(draw_seed),
                "skipped": False,
                "mode": "quark_only",
                "params": params,
                "fit_diagnostics": fit_diagnostics,
                "constraints": payload["constraints"],
                "survives_all_HARD_strict": payload["survives_all_HARD_strict"],
                "survives_all_HARD_inclusive": payload["survives_all_HARD_inclusive"],
                "excluded_by_rigorous": payload["excluded_by_rigorous"],
                "excluded_by_proxy": payload["excluded_by_proxy"],
                "hard_not_evaluated": payload["hard_not_evaluated"],
                "coverage_complete": payload["coverage_complete"],
                "advisory_flags": [
                    QUARK_ONLY_DEFERRED_SCOPE_TAG,
                    *payload["advisory_flags"],
                ],
                "allowlist": list(QUARK_ONLY_ALLOWLIST_IDS),
                "lepton_sector": QUARK_ONLY_LEPTON_SECTOR_LABEL,
                "deferred_lepton_followup": list(QUARK_ONLY_DEFERRED_LEPTON_FOLLOWUP),
                "cache_metrics": {
                    "spectrum_injected": point.extras.get("rs_ew_spectrum") is spectrum,
                    "rs_ew_cache_injected": True,
                    "spline_max_a_rel_err": float(overlap_cache.max_a_rel_err),
                },
                "provenance": {
                    **dict(provenance),
                    "registry_count": int(registry_count),
                    "config_hash": config_hash,
                },
            }
        lepton_inputs = _draw_lepton_inputs(rng, cfg)
        params["lepton_inputs"] = dict(lepton_inputs)

        quark_solution = fit_quark_sector(
            default_quark_targets(),
            r=cfg.quark_fit_r,
            seed=quark_seed,
            Lambda_IR=tile.lambda_ir_gev,
            k=cfg.k_gev,
            max_nfev=cfg.quark_fit_max_nfev,
            fit_orientation=True,
        )
        fit_result = quark_solution.result
        fit_diagnostics = _quark_fit_diagnostics(quark_solution)
        _require_valid_quark_fit(fit_result, fit_diagnostics, cfg)
        _require_no_c_half_singularity(
            [
                *np.asarray(fit_result.bulk_state.c_Q, dtype=float),
                *np.asarray(fit_result.bulk_state.c_u, dtype=float),
                *np.asarray(fit_result.bulk_state.c_d, dtype=float),
                lepton_inputs["c_L"],
                *lepton_inputs["c_E"],
                lepton_inputs["c_N"],
            ],
            cfg=cfg,
        )

        lepton_yukawas = compute_all_yukawas(
            Lambda_IR=tile.lambda_ir_gev,
            c_L=float(lepton_inputs["c_L"]),
            c_E=tuple(float(x) for x in lepton_inputs["c_E"]),
            c_N=float(lepton_inputs["c_N"]),
            M_N=float(lepton_inputs["M_N"]),
            lightest_nu_mass=float(lepton_inputs["lightest_nu_mass"]),
            ordering=str(lepton_inputs["ordering"]),
            majorana_alpha=float(lepton_inputs["majorana_alpha"]),
            majorana_beta=float(lepton_inputs["majorana_beta"]),
            k=cfg.k_gev,
            v=cfg.v_gev,
        )
        _require_perturbative_leptons(lepton_yukawas, cfg)

        rs_point = point_builder.build_from_rs_ew_inputs(
            fit_result,
            Lambda_IR=tile.lambda_ir_gev,
            k=cfg.k_gev,
            n_gauge_modes=cfg.n_gauge_modes,
            quadrature_order=cfg.quadrature_order,
            min_overlap_modes=cfg.min_overlap_modes,
            max_overlap_modes=cfg.max_overlap_modes,
            overlap_rel_tol=cfg.overlap_rel_tol,
            spectrum=spectrum,
            rs_ew_cache=overlap_cache,
            include_charged_current=True,
            include_fermion_kk_mixing=True,
            include_higgs_yukawas=True,
            lepton_yukawa_result=lepton_yukawas,
            raw={"tile_id": tile.tile_id, "draw_id": draw_idx, "params": params},
        )
        quark_couplings = compute_quark_kk_gluon_couplings(
            fit_result,
            M_KK=tile.mkk_gev,
            xi_KK=cfg.xi_kk,
            g_s_star=None,
        )
        point = point_builder.make_point(
            raw=rs_point.raw,
            **dict(rs_point.extras),
            quark_mass_basis_couplings=quark_couplings,
            kk_gluon_mass_gev=float(tile.mkk_gev),
        )
        results = registry.evaluate_all(point)
        payload = _classify_results(results)
        return {
            "draw_id": int(draw_idx),
            "tile_id": int(tile.tile_id),
            "seed": int(draw_seed),
            "skipped": False,
            "params": params,
            "fit_diagnostics": fit_diagnostics,
            "constraints": payload["constraints"],
            "survives_all_HARD_strict": payload["survives_all_HARD_strict"],
            "survives_all_HARD_inclusive": payload["survives_all_HARD_inclusive"],
            "excluded_by_rigorous": payload["excluded_by_rigorous"],
            "excluded_by_proxy": payload["excluded_by_proxy"],
            "hard_not_evaluated": payload["hard_not_evaluated"],
            "coverage_complete": payload["coverage_complete"],
            "advisory_flags": payload["advisory_flags"],
            "cache_metrics": {
                "spectrum_injected": point.extras.get("rs_ew_spectrum") is spectrum,
                "rs_ew_cache_injected": True,
                "spline_max_a_rel_err": float(overlap_cache.max_a_rel_err),
            },
            "provenance": {
                **dict(provenance),
                "registry_count": int(registry_count),
                "config_hash": config_hash,
            },
        }
    except Exception as exc:  # noqa: BLE001 - a failed draw must not abort a tile
        reason = _skip_reason(exc)
        row = {
            "draw_id": int(draw_idx),
            "tile_id": int(tile.tile_id),
            "seed": int(draw_seed),
            "skipped": True,
            "skip_reason": reason,
            "skip_error": f"{type(exc).__name__}: {exc}",
            "params": params,
            "fit_diagnostics": fit_diagnostics,
            "constraints": {},
            "survives_all_HARD_strict": False,
            "survives_all_HARD_inclusive": False,
            "excluded_by_rigorous": [],
            "excluded_by_proxy": [],
            "hard_not_evaluated": [],
            "coverage_complete": False,
            "advisory_flags": [reason],
            "cache_metrics": {
                "spectrum_injected": True,
                "rs_ew_cache_injected": True,
                "spline_max_a_rel_err": float(overlap_cache.max_a_rel_err),
            },
            "provenance": {
                **dict(provenance),
                "registry_count": int(registry_count),
                "config_hash": config_hash,
            },
        }
        if cfg.quark_only:
            row.update(
                {
                    "mode": "quark_only",
                    "allowlist": list(QUARK_ONLY_ALLOWLIST_IDS),
                    "lepton_sector": QUARK_ONLY_LEPTON_SECTOR_LABEL,
                    "deferred_lepton_followup": list(QUARK_ONLY_DEFERRED_LEPTON_FOLLOWUP),
                    "advisory_flags": [
                        QUARK_ONLY_DEFERRED_SCOPE_TAG,
                        *row["advisory_flags"],
                    ],
                }
            )
        return row


def _draw_quark_seed(
    rng: np.random.Generator,
    cfg: ScanConfig,
) -> tuple[np.ndarray, np.ndarray, QuarkFitSeed]:
    y_u = _draw_anarchic_matrix(rng, cfg)
    y_d = _draw_anarchic_matrix(rng, cfg)
    up_s, up_left, up_right = _svd_seed_parts(y_u)
    down_s, down_left, down_right = _svd_seed_parts(y_d)
    return (
        y_u,
        y_d,
        QuarkFitSeed(
            up_singular_values=up_s,
            down_singular_values=down_s,
            overall_scale=1.0,
            up_left=up_left,
            up_right=up_right,
            down_left=down_left,
            down_right=down_right,
        ),
    )


def _draw_anarchic_matrix(rng: np.random.Generator, cfg: ScanConfig) -> np.ndarray:
    prior = cfg.y_prior.lower()
    if cfg.y_floor < 0.0:
        raise ValueError("y_floor must be non-negative")
    if prior == "uniform":
        if cfg.y_half_range <= 0.0:
            raise ValueError("y_half_range must be positive")
        for _ in range(256):
            y = rng.uniform(-cfg.y_half_range, cfg.y_half_range, size=(3, 3)) + 1j * rng.uniform(
                -cfg.y_half_range, cfg.y_half_range, size=(3, 3)
            )
            if (np.abs(y) >= cfg.y_floor).all():
                return y
        y = rng.uniform(-cfg.y_half_range, cfg.y_half_range, size=(3, 3)) + 1j * rng.uniform(
            -cfg.y_half_range, cfg.y_half_range, size=(3, 3)
        )
    elif prior == "gaussian":
        if cfg.y_sigma <= 0.0 or cfg.y_trunc_sigma <= 0.0:
            raise ValueError("gaussian prior requires positive sigma and trunc_sigma")
        cap = cfg.y_sigma * cfg.y_trunc_sigma
        for _ in range(256):
            re = np.clip(rng.normal(0.0, cfg.y_sigma, size=(3, 3)), -cap, cap)
            im = np.clip(rng.normal(0.0, cfg.y_sigma, size=(3, 3)), -cap, cap)
            y = re + 1j * im
            if (np.abs(y) >= cfg.y_floor).all():
                return y
        y = re + 1j * im
    else:
        raise ValueError(f"unknown Yukawa prior: {cfg.y_prior!r}")
    mag = np.abs(y)
    phase = y / np.where(mag > 0.0, mag, 1.0)
    return phase * np.maximum(mag, cfg.y_floor)


def _draw_lepton_inputs(rng: np.random.Generator, cfg: ScanConfig) -> dict[str, Any]:
    c_l = float(rng.uniform(cfg.c_min, cfg.c_max))
    c_e = [float(x) for x in rng.uniform(cfg.c_min, cfg.c_max, size=3)]
    c_n = float(rng.uniform(cfg.c_min, cfg.c_max))
    log_mn = rng.uniform(math.log10(cfg.m_n_min_gev), math.log10(cfg.m_n_max_gev))
    ordering = "normal" if rng.random() < cfg.normal_ordering_probability else "inverted"
    return {
        "c_L": c_l,
        "c_E": c_e,
        "c_N": c_n,
        "M_N": float(10.0**log_mn),
        "lightest_nu_mass": float(
            rng.uniform(cfg.lightest_nu_min_ev, cfg.lightest_nu_max_ev)
        ),
        "ordering": ordering,
        "majorana_alpha": float(rng.uniform(0.0, 2.0 * math.pi)),
        "majorana_beta": float(rng.uniform(0.0, 2.0 * math.pi)),
    }


def _svd_seed_parts(matrix: np.ndarray) -> tuple[np.ndarray, RotationParameters, RotationParameters]:
    u, s, vh = np.linalg.svd(np.asarray(matrix, dtype=np.complex128), full_matrices=True)
    order = np.argsort(s, kind="mergesort")
    s_sorted = np.asarray(s[order], dtype=float)
    u_sorted = u[:, order]
    v_sorted = vh.conjugate().T[:, order]
    return s_sorted, _rotation_from_unitary(u_sorted), _rotation_from_unitary(v_sorted)


def _rotation_from_unitary(matrix: np.ndarray) -> RotationParameters:
    u = np.asarray(matrix, dtype=np.complex128)
    if u.shape != (3, 3):
        raise ValueError("rotation unitary must have shape (3, 3)")
    theta13 = float(np.arcsin(np.clip(np.abs(u[0, 2]), 0.0, 1.0)))
    theta12 = float(np.arctan2(np.abs(u[0, 1]), max(np.abs(u[0, 0]), 1.0e-15)))
    theta23 = float(np.arctan2(np.abs(u[1, 2]), max(np.abs(u[2, 2]), 1.0e-15)))
    quartet = u[0, 0] * u[1, 2] * np.conjugate(u[0, 2]) * np.conjugate(u[2, 2])
    delta = float(((np.angle(quartet) + math.pi) % (2.0 * math.pi)) - math.pi)
    return RotationParameters(theta12=theta12, theta13=theta13, theta23=theta23, delta=delta)


def _quark_fit_diagnostics(solution: Any) -> dict[str, Any]:
    result = solution.result
    return {
        "success": bool(solution.success),
        "message": str(solution.message),
        "nfev": int(solution.nfev),
        "initial_score": _finite_or_none(solution.initial_score),
        "score": _finite_or_none(result.score),
        "residual_norm": _finite_or_none(result.residual_norm),
        "max_up_log_residual": _max_abs_or_none(result.mass_residuals_up),
        "max_down_log_residual": _max_abs_or_none(result.mass_residuals_down),
        "max_ckm_residual": _max_abs_or_none(result.ckm_residuals),
        "max_abs_quark_yukawa": _finite_or_none(
            max(
                float(np.max(np.abs(result.point.Y_u))),
                float(np.max(np.abs(result.point.Y_d))),
            )
        ),
        "bulk_c_Q": _array_real(result.bulk_state.c_Q),
        "bulk_c_u": _array_real(result.bulk_state.c_u),
        "bulk_c_d": _array_real(result.bulk_state.c_d),
    }


def _require_valid_quark_fit(
    fit_result: QuarkFitResult,
    fit_diagnostics: Mapping[str, Any],
    cfg: ScanConfig,
) -> None:
    if cfg.require_fit_success and not bool(fit_diagnostics.get("success")):
        raise RuntimeError("quark_fit_failed")
    score = float(fit_diagnostics.get("score", math.inf))
    if not math.isfinite(score) or score > cfg.max_fit_score:
        raise RuntimeError(f"quark_fit_score_exceeded:{score:.6g}")
    max_y = max(float(np.max(np.abs(fit_result.point.Y_u))), float(np.max(np.abs(fit_result.point.Y_d))))
    if not math.isfinite(max_y) or max_y > cfg.quark_y_abs_max:
        raise RuntimeError(f"nonperturbative_quark_yukawa:{max_y:.6g}")


def _require_perturbative_leptons(lepton_yukawas: YukawaResult, cfg: ScanConfig) -> None:
    max_ybar = max(
        float(np.max(np.abs(lepton_yukawas.Y_E_bar))),
        float(np.max(np.abs(lepton_yukawas.Y_N_bar))),
    )
    if not math.isfinite(max_ybar) or not lepton_yukawas.is_perturbative(cfg.perturbative_ybar_max):
        raise RuntimeError(f"nonperturbative_lepton_yukawa:{max_ybar:.6g}")


def _force_degenerate_neutrino_yukawas(lepton_yukawas: YukawaResult, *, k: float) -> None:
    """Mutate a sanity-only lepton result to remove PMNS-induced LFV spurions."""

    from neutrinos.neutrinoValues import get_pmns

    ybar = float(np.min(np.abs(lepton_yukawas.Y_N_bar)))
    lepton_yukawas.Y_N_bar = np.full(3, ybar, dtype=np.complex128)
    lepton_yukawas.Y_N = lepton_yukawas.Y_N_bar / (2.0 * float(k))
    params = dict(lepton_yukawas.params)
    pmns = get_pmns(
        str(params["ordering"]),
        float(params["majorana_alpha"]),
        float(params["majorana_beta"]),
    )
    lepton_yukawas.Y_N_matrix = pmns @ np.diag(lepton_yukawas.Y_N)


def _require_no_c_half_singularity(values: Sequence[float], *, cfg: ScanConfig) -> None:
    for value in values:
        if abs(float(value) - 0.5) <= cfg.c_half_atol:
            raise RuntimeError(f"c_half_singularity:{float(value):.17g}")


def _classify_results(results: Mapping[str, ConstraintResult]) -> dict[str, Any]:
    constraints: dict[str, Any] = {}
    excluded_by_rigorous: list[str] = []
    excluded_by_proxy: list[str] = []
    hard_not_evaluated: list[str] = []
    advisory_flags: list[str] = []

    for pid in sorted(results):
        result = results[pid]
        diag = dict(result.diagnostics)
        evaluated = _result_evaluated(result)
        active = _result_active(result, evaluated=evaluated)
        tag, matching_status, needs_human, proxy_flags = tag_result(result)
        severity = result.severity.value
        constraints[pid] = {
            "passes": bool(result.passes),
            "severity": severity,
            "ratio": _finite_or_none(result.ratio),
            "active": active,
            "evaluated": evaluated,
            "match_status": matching_status,
            "tag": tag,
            "needs_human_physics": needs_human,
            "proxy_flags": proxy_flags,
        }
        if "exception_type" in diag:
            constraints[pid]["exception_type"] = str(diag["exception_type"])
            constraints[pid]["exception"] = str(
                diag.get("exception", diag.get("exception_message", ""))
            )

        if result.severity is Severity.HARD:
            if evaluated and tag == "rigorous":
                if not result.passes:
                    excluded_by_rigorous.append(pid)
            elif evaluated and tag == "proxy":
                if not result.passes:
                    excluded_by_proxy.append(pid)
            else:
                hard_not_evaluated.append(pid)
        elif not result.passes:
            advisory_flags.append(f"{severity}:{pid}")

    exception_ids = [
        pid
        for pid, item in constraints.items()
        if item.get("exception_type")
    ]
    if exception_ids:
        advisory_flags.append("constraint_exceptions:" + ",".join(exception_ids))
    if hard_not_evaluated:
        advisory_flags.append("hard_coverage_gap:" + ",".join(hard_not_evaluated))

    return {
        "constraints": constraints,
        "survives_all_HARD_strict": not excluded_by_rigorous,
        "survives_all_HARD_inclusive": not excluded_by_rigorous and not excluded_by_proxy,
        "excluded_by_rigorous": excluded_by_rigorous,
        "excluded_by_proxy": excluded_by_proxy,
        "hard_not_evaluated": hard_not_evaluated,
        "coverage_complete": not hard_not_evaluated,
        "advisory_flags": advisory_flags,
    }


def _evaluate_constraint_ids(
    point: point_builder.ParameterPoint | Any,
    process_ids: Sequence[str],
) -> dict[str, ConstraintResult]:
    """Evaluate a selected process-id list with registry-style isolation."""

    available = registry.all_constraints()
    missing = sorted(set(process_ids) - set(available))
    if missing:
        raise RuntimeError(f"quark-only allowlist ids are not registered: {missing}")
    out: dict[str, ConstraintResult] = {}
    for pid in process_ids:
        constraint = available[str(pid)]
        try:
            out[str(pid)] = constraint.evaluate(point)
        except Exception as exc:  # noqa: BLE001 - mirror registry.evaluate_all isolation
            out[str(pid)] = ConstraintResult(
                process_id=str(pid),
                severity=getattr(constraint, "severity", Severity.HARD),
                passes=False,
                notes=f"evaluate() raised {type(exc).__name__}: {exc}",
                diagnostics={"exception_type": type(exc).__name__, "exception": str(exc)},
            )
    return out


def tag_result(result: ConstraintResult) -> tuple[str, str | None, str | None, dict[str, Any]]:
    """Return ``(tag, matching_status, needs_human, proxy_flags)`` for a result."""

    diag = dict(result.diagnostics)
    evaluated = _result_evaluated(result)
    matching_status = _matching_status(diag)
    needs_human = _string_or_none(diag.get("needs_human_physics"))
    proxy_flags = _proxy_flags(diag)
    status_text = "" if matching_status is None else matching_status.lower()
    needs_text = "" if needs_human is None else needs_human.lower()
    if not evaluated or "missing_extra" in diag or "exception_type" in diag:
        return "stub", matching_status, needs_human, proxy_flags
    if "stub" in status_text or "quarantine" in status_text:
        return "stub", matching_status, needs_human, proxy_flags
    if proxy_flags:
        return "proxy", matching_status, needs_human, proxy_flags
    if needs_human:
        if "proxy" in needs_text or "recast" in needs_text:
            return "proxy", matching_status, needs_human, proxy_flags
        return "partial", matching_status, needs_human, proxy_flags
    if "partial" in status_text or "needs-human" in status_text or "deferred" in status_text:
        return "partial", matching_status, needs_human, proxy_flags
    if "proxy" in status_text or "recast" in status_text:
        resolved_proxy_phrases = (
            "no proxy",
            "proxy resolved",
            "proxy is not used",
            "proxy not used",
            "not a proxy",
        )
        if "rigorous" not in status_text and not any(phrase in status_text for phrase in resolved_proxy_phrases):
            return "proxy", matching_status, needs_human, proxy_flags
    if "recast" in status_text:
        return "proxy", matching_status, needs_human, proxy_flags
    return "rigorous", matching_status, needs_human, proxy_flags


def _result_evaluated(result: ConstraintResult) -> bool:
    diag = dict(result.diagnostics)
    if "exception_type" in diag:
        return False
    if "evaluated" in diag:
        return bool(diag["evaluated"])
    if "missing_extra" in diag or "unevaluated_reason" in diag:
        return False
    return True


def _result_active(result: ConstraintResult, *, evaluated: bool) -> bool:
    diag = dict(result.diagnostics)
    if "active" in diag:
        return bool(diag["active"])
    if not evaluated:
        return False
    return result.ratio is not None or result.predicted is not None or result.budget is not None


def _is_sm_tension_only(result: ConstraintResult) -> bool:
    if result.passes or result.severity is not Severity.HARD or not _result_evaluated(result):
        return False
    predicted = _finite_or_none(result.predicted)
    sm_prediction = _finite_or_none(result.sm_prediction)
    if predicted is None or sm_prediction is None:
        return False
    scale = max(abs(predicted), abs(sm_prediction), 1.0e-30)
    if abs(predicted - sm_prediction) > 1.0e-9 * scale:
        return False
    diag = dict(result.diagnostics)
    for key, value in diag.items():
        key_l = str(key).lower()
        if "np_" not in key_l and not key_l.startswith("delta_"):
            continue
        numeric = _finite_or_none(value)
        if numeric is None:
            continue
        if abs(numeric) > 1.0e-12 * max(abs(sm_prediction), 1.0e-30):
            return False
    return True


def _matching_status(diag: Mapping[str, Any]) -> str | None:
    status_keys = [key for key in diag if key == "matching_status" or key.endswith("_matching_status")]
    if not status_keys:
        return None
    if len(status_keys) == 1:
        return _string_or_none(diag[status_keys[0]])
    payload = {key: _json_sanitize(diag[key]) for key in sorted(status_keys)}
    return json.dumps(payload, sort_keys=True, separators=(",", ":"))


def _proxy_flags(diag: Mapping[str, Any]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for key, value in diag.items():
        key_l = str(key).lower()
        if "proxy" not in key_l and "recast" not in key_l:
            continue
        if value is False or value is None:
            continue
        if isinstance(value, (int, float, str, bool)):
            out[str(key)] = _json_sanitize(value)
        else:
            out[str(key)] = _string_or_none(value)
    return out


def _accumulate_row(
    row: Mapping[str, Any],
    *,
    counters: Counter[str],
    hard_veto_rigorous: Counter[str],
    hard_veto_proxy: Counter[str],
    hard_not_evaluated: Counter[str],
    tag_counts: Counter[str],
    exception_ids: Counter[str],
    constraint_tallies: dict[str, dict[str, Any]] | None = None,
) -> None:
    counters["rows"] += 1
    if row.get("skipped"):
        counters["skipped"] += 1
        counters[f"skip:{row.get('skip_reason', 'unknown')}"] += 1
        return
    counters["evaluated_points"] += 1
    if row.get("survives_all_HARD_strict"):
        counters["survives_strict"] += 1
    if row.get("survives_all_HARD_inclusive"):
        counters["survives_inclusive"] += 1
    for pid, item in dict(row.get("constraints", {})).items():
        counters["constraint_results"] += 1
        if item.get("evaluated"):
            counters["constraint_evaluated"] += 1
        if item.get("active"):
            counters["constraint_active"] += 1
        tag_counts[str(item.get("tag", "unknown"))] += 1
        if item.get("exception_type"):
            counters["constraint_exceptions"] += 1
            exception_ids[str(pid)] += 1
        if constraint_tallies is not None:
            _accumulate_constraint_tally(constraint_tallies, str(pid), item)
    hard_veto_rigorous.update(str(x) for x in row.get("excluded_by_rigorous", []))
    hard_veto_proxy.update(str(x) for x in row.get("excluded_by_proxy", []))
    hard_not_evaluated.update(str(x) for x in row.get("hard_not_evaluated", []))


def _empty_constraint_tallies(process_ids: Sequence[str]) -> dict[str, dict[str, Any]]:
    return {
        str(pid): {
            "points": 0,
            "evaluated": 0,
            "active": 0,
            "failed": 0,
            "vetoed": 0,
            "severity": None,
            "tag": None,
            "tag_counts": {},
        }
        for pid in process_ids
    }


def _accumulate_constraint_tally(
    tallies: dict[str, dict[str, Any]],
    pid: str,
    item: Mapping[str, Any],
) -> None:
    tally = tallies.setdefault(pid, _empty_constraint_tallies((pid,))[pid])
    tally["points"] = int(tally.get("points", 0)) + 1
    if item.get("evaluated"):
        tally["evaluated"] = int(tally.get("evaluated", 0)) + 1
    if item.get("active"):
        tally["active"] = int(tally.get("active", 0)) + 1
    if not bool(item.get("passes", True)):
        tally["failed"] = int(tally.get("failed", 0)) + 1
        if (
            str(item.get("severity")) == Severity.HARD.value
            and bool(item.get("evaluated"))
            and str(item.get("tag")) in {"rigorous", "proxy"}
        ):
            tally["vetoed"] = int(tally.get("vetoed", 0)) + 1
    severity = item.get("severity")
    if severity is not None:
        tally["severity"] = str(severity)
    tag = str(item.get("tag", "unknown"))
    tag_counts = dict(tally.get("tag_counts", {}))
    tag_counts[tag] = int(tag_counts.get(tag, 0)) + 1
    tally["tag_counts"] = tag_counts
    if len(tag_counts) == 1:
        tally["tag"] = tag
    else:
        tally["tag"] = "mixed"


def _finalize_constraint_tallies(
    tallies: Mapping[str, Mapping[str, Any]]
) -> dict[str, dict[str, Any]]:
    return {
        str(pid): {
            "points": int(item.get("points", 0)),
            "evaluated": int(item.get("evaluated", 0)),
            "active": int(item.get("active", 0)),
            "failed": int(item.get("failed", 0)),
            "vetoed": int(item.get("vetoed", 0)),
            "severity": item.get("severity"),
            "tag": item.get("tag"),
            "tag_counts": dict(sorted(dict(item.get("tag_counts", {})).items())),
        }
        for pid, item in sorted(tallies.items())
    }


def _merge_constraint_tallies(
    target: dict[str, dict[str, Any]],
    source: Mapping[str, Any],
) -> None:
    for pid, raw_item in dict(source).items():
        item = dict(raw_item)
        tally = target.setdefault(str(pid), _empty_constraint_tallies((str(pid),))[str(pid)])
        for key in ("points", "evaluated", "active", "failed", "vetoed"):
            tally[key] = int(tally.get(key, 0)) + int(item.get(key, 0))
        severity = item.get("severity")
        if severity is not None:
            tally["severity"] = str(severity)
        tag_counts = dict(tally.get("tag_counts", {}))
        for tag, count in dict(item.get("tag_counts", {})).items():
            tag_counts[str(tag)] = int(tag_counts.get(str(tag), 0)) + int(count)
        tally["tag_counts"] = tag_counts
        if len(tag_counts) == 1:
            tally["tag"] = next(iter(tag_counts))
        elif tag_counts:
            tally["tag"] = "mixed"


def run_universal_c_sanity(
    cfg: ScanConfig,
    *,
    provenance: Mapping[str, Any],
    config_hash: str,
) -> dict[str, Any]:
    """Evaluate a decoupled universal-c point and require no HARD veto lists."""

    registry.discover()
    failures = registry.import_failures()
    count = len(registry.all_constraints())
    if count != cfg.expected_registry_count or failures:
        raise RuntimeError(f"registry sanity failed: count={count}, failures={sorted(failures)}")
    lambda_ir = cfg.sanity_mkk_gev / cfg.xi_kk
    spectrum = RSEWSpectrum.build(
        lambda_ir_gev=lambda_ir,
        k_gev=cfg.k_gev,
        n_gauge_modes=cfg.n_gauge_modes,
        quadrature_order=cfg.quadrature_order,
        model_label="minimal_rs",
    )
    cache = RSEWOverlapSplineCache.build(
        spectrum,
        c_min=cfg.c_min,
        c_max=cfg.c_max,
        grid_size=cfg.spline_grid_size,
        include_omega=True,
        verify_points=cfg.spline_verify_points,
        rel_tol=cfg.overlap_rel_tol,
        min_modes=cfg.min_overlap_modes,
        max_modes=cfg.max_overlap_modes,
    )
    fit = fit_quark_sector(
        default_quark_targets(),
        r=cfg.quark_fit_r,
        Lambda_IR=lambda_ir,
        k=cfg.k_gev,
        seed=None,
        max_nfev=cfg.quark_fit_max_nfev,
    ).result
    c = float(cfg.sanity_c)
    epsilon = lambda_ir / cfg.k_gev
    f = float(f_IR(c, epsilon))
    from dataclasses import replace

    universal_state = replace(
        fit.bulk_state,
        c_Q=np.full(3, c, dtype=float),
        c_u=np.full(3, c, dtype=float),
        c_d=np.full(3, c, dtype=float),
        F_Q=np.full(3, f, dtype=float),
        F_u=np.full(3, f, dtype=float),
        F_d=np.full(3, f, dtype=float),
    )
    fit = replace(
        fit,
        bulk_state=universal_state,
        U_L_u=np.eye(3, dtype=np.complex128),
        U_L_d=np.eye(3, dtype=np.complex128),
        U_R_u=np.eye(3, dtype=np.complex128),
        U_R_d=np.eye(3, dtype=np.complex128),
        ckm=np.eye(3, dtype=np.complex128),
    )
    lepton_yukawas = compute_all_yukawas(
        Lambda_IR=lambda_ir,
        c_L=c,
        c_E=(c, c, c),
        c_N=c,
        M_N=1.0e15,
        lightest_nu_mass=0.001,
        ordering="normal",
        k=cfg.k_gev,
        v=cfg.v_gev,
    )
    _force_degenerate_neutrino_yukawas(lepton_yukawas, k=cfg.k_gev)
    rs_point = point_builder.build_from_rs_ew_inputs(
        fit,
        Lambda_IR=lambda_ir,
        k=cfg.k_gev,
        n_gauge_modes=cfg.n_gauge_modes,
        quadrature_order=cfg.quadrature_order,
        min_overlap_modes=cfg.min_overlap_modes,
        max_overlap_modes=cfg.max_overlap_modes,
        overlap_rel_tol=cfg.overlap_rel_tol,
        spectrum=spectrum,
        rs_ew_cache=cache,
        include_charged_current=True,
        include_fermion_kk_mixing=True,
        include_higgs_yukawas=True,
        lepton_yukawa_result=lepton_yukawas,
        raw={"sanity": "universal_c_diagonal_leptons"},
    )
    point = point_builder.make_point(
        raw=rs_point.raw,
        **dict(rs_point.extras),
        quark_mass_basis_couplings=compute_quark_kk_gluon_couplings(
            fit,
            M_KK=cfg.sanity_mkk_gev,
            xi_KK=cfg.xi_kk,
            g_s_star=None,
        ),
        kk_gluon_mass_gev=float(cfg.sanity_mkk_gev),
    )
    results = registry.evaluate_all(point)
    payload = _classify_results(results)
    sm_tension_ids = [
        pid
        for pid in [*payload["excluded_by_rigorous"], *payload["excluded_by_proxy"]]
        if _is_sm_tension_only(results[pid])
    ]
    spurious_rigorous = [
        pid for pid in payload["excluded_by_rigorous"] if pid not in sm_tension_ids
    ]
    spurious_proxy = [
        pid for pid in payload["excluded_by_proxy"] if pid not in sm_tension_ids
    ]
    return {
        "name": "universal_c_diagonal_leptons",
        "mkk_gev": float(cfg.sanity_mkk_gev),
        "lambda_ir_gev": float(lambda_ir),
        "c": c,
        "registry_count": count,
        "excluded_by_rigorous": spurious_rigorous,
        "excluded_by_proxy": spurious_proxy,
        "raw_excluded_by_rigorous": payload["excluded_by_rigorous"],
        "raw_excluded_by_proxy": payload["excluded_by_proxy"],
        "sm_tension_hard_exclusions": sm_tension_ids,
        "hard_not_evaluated": payload["hard_not_evaluated"],
        "passes_no_spurious_hard_exclusions": (
            not spurious_rigorous and not spurious_proxy
        ),
        "coverage_complete": payload["coverage_complete"],
        "provenance": {**dict(provenance), "config_hash": config_hash},
    }


def _build_run_summary(
    summaries: Sequence[Mapping[str, Any]],
    *,
    cfg: ScanConfig,
    config_hash: str,
    provenance: Mapping[str, Any],
    sanity: Mapping[str, Any] | None,
) -> dict[str, Any]:
    totals: Counter[str] = Counter()
    rigorous = Counter()
    proxy = Counter()
    hard_gap = Counter()
    tags = Counter()
    exceptions = Counter()
    constraint_tallies = (
        _empty_constraint_tallies(QUARK_ONLY_ALLOWLIST_IDS) if cfg.quark_only else None
    )
    cache_seconds = 0.0
    draw_loop_seconds = 0.0
    elapsed_seconds = 0.0
    for summary in summaries:
        totals["n_requested"] += int(summary.get("n_requested", 0))
        totals["n_rows"] += int(summary.get("n_rows", 0))
        totals["n_evaluated_points"] += int(summary.get("n_evaluated_points", 0))
        totals["n_skipped"] += int(summary.get("n_skipped", 0))
        totals["constraint_results"] += int(summary.get("constraint_evaluations", 0))
        totals["constraint_evaluated"] += int(summary.get("constraint_evaluated", 0))
        totals["constraint_active"] += int(summary.get("constraint_active", 0))
        totals["constraint_exceptions"] += int(summary.get("constraint_exceptions", 0))
        totals["survives_strict"] += int(summary.get("survives_all_HARD_strict", 0))
        totals["survives_inclusive"] += int(summary.get("survives_all_HARD_inclusive", 0))
        rigorous.update(summary.get("hard_vetoes_rigorous", {}))
        proxy.update(summary.get("hard_vetoes_proxy", {}))
        hard_gap.update(summary.get("hard_not_evaluated", {}))
        tags.update(summary.get("tag_counts", {}))
        exceptions.update(summary.get("exception_ids", {}))
        if constraint_tallies is not None:
            _merge_constraint_tallies(
                constraint_tallies,
                summary.get("constraint_tallies", {}),
            )
        timing = dict(summary.get("timing", {}))
        cache_seconds += float(timing.get("cache_build_seconds", 0.0))
        draw_loop_seconds += float(timing.get("draw_loop_seconds", 0.0))
        elapsed_seconds += float(timing.get("elapsed_seconds", 0.0))

    post_cache_per_draw = _safe_div(draw_loop_seconds, totals["n_requested"])
    post_cache_per_evaluated = _safe_div(draw_loop_seconds, totals["n_evaluated_points"])
    out = {
        "schema": "full_catalog_scan_w6b_summary_v1",
        "config": _config_payload(cfg),
        "config_hash": config_hash,
        "provenance": dict(provenance),
        "n_tiles": len(summaries),
        "totals": dict(totals),
        "timing": {
            "tile_elapsed_seconds_sum": float(elapsed_seconds),
            "cache_build_seconds_sum": float(cache_seconds),
            "draw_loop_seconds_sum": float(draw_loop_seconds),
            "post_cache_seconds_per_draw": post_cache_per_draw,
            "post_cache_seconds_per_evaluated_point": post_cache_per_evaluated,
            "extrapolated_1e8_core_hours_per_draw": (
                None if post_cache_per_draw is None else float(post_cache_per_draw * 1.0e8 / 3600.0)
            ),
            "extrapolated_1e8_core_hours_per_evaluated_point": (
                None
                if post_cache_per_evaluated is None
                else float(post_cache_per_evaluated * 1.0e8 / 3600.0)
            ),
        },
        "constraint_counts": {
            "evaluated_per_evaluated_point": _safe_div(
                totals["constraint_evaluated"], totals["n_evaluated_points"]
            ),
            "active_per_evaluated_point": _safe_div(
                totals["constraint_active"], totals["n_evaluated_points"]
            ),
            "exception_rate": _safe_div(totals["constraint_exceptions"], totals["constraint_results"]),
            "tag_counts": dict(sorted(tags.items())),
        },
        "top_hard_vetoes_rigorous": rigorous.most_common(20),
        "top_hard_vetoes_proxy": proxy.most_common(20),
        "top_hard_not_evaluated": hard_gap.most_common(20),
        "exception_ids": exceptions.most_common(20),
        "universal_c_sanity": sanity,
        "tiles": list(summaries),
    }
    if cfg.quark_only:
        out.update(
            {
                "schema": "quark_only_bucket1_scan_wq_summary_v1",
                "mode": "quark_only",
                "allowlist": list(QUARK_ONLY_ALLOWLIST_IDS),
                "candidate_allowlist": list(QUARK_ONLY_CANDIDATE_IDS),
                "lepton_sector": QUARK_ONLY_LEPTON_SECTOR_LABEL,
                "deferred_scope_tag": QUARK_ONLY_DEFERRED_SCOPE_TAG,
                "deferred_lepton_followup": list(QUARK_ONLY_DEFERRED_LEPTON_FOLLOWUP),
                "build_include_flags": dict(QUARK_ONLY_BUILD_INCLUDE_FLAGS),
                "constraint_tallies": _finalize_constraint_tallies(
                    constraint_tallies or {}
                ),
            }
        )
    return out


def _write_markdown_report(path: Path, summary: Mapping[str, Any]) -> None:
    timing = dict(summary.get("timing", {}))
    totals = dict(summary.get("totals", {}))
    counts = dict(summary.get("constraint_counts", {}))
    sanity = dict(summary.get("universal_c_sanity") or {})
    lines = [
        "# W6b Full-Catalog Smoke Report",
        "",
        f"- draws requested: {totals.get('n_requested', 0)}",
        f"- evaluated points: {totals.get('n_evaluated_points', 0)}",
        f"- skipped draws: {totals.get('n_skipped', 0)}",
        f"- post-cache seconds/draw: {_fmt_optional(timing.get('post_cache_seconds_per_draw'), 6)}",
        (
            "- post-cache seconds/evaluated point: "
            f"{_fmt_optional(timing.get('post_cache_seconds_per_evaluated_point'), 6)}"
        ),
        (
            "- extrapolated 1e8 core-hours per draw: "
            f"{_fmt_optional(timing.get('extrapolated_1e8_core_hours_per_draw'), 2)}"
        ),
        (
            "- extrapolated 1e8 core-hours per evaluated point: "
            f"{_fmt_optional(timing.get('extrapolated_1e8_core_hours_per_evaluated_point'), 2)}"
        ),
        f"- evaluated constraints per point: {_fmt_optional(counts.get('evaluated_per_evaluated_point'), 3)}",
        f"- active constraints per point: {_fmt_optional(counts.get('active_per_evaluated_point'), 3)}",
        f"- constraint exception rate: {_fmt_optional(counts.get('exception_rate'), 6)}",
        f"- top rigorous HARD vetoes: {summary.get('top_hard_vetoes_rigorous', [])[:10]}",
        f"- top proxy HARD vetoes: {summary.get('top_hard_vetoes_proxy', [])[:10]}",
        (
            "- universal-c sanity: "
            f"{bool(sanity.get('passes_no_spurious_hard_exclusions'))}; "
            f"rigorous={sanity.get('excluded_by_rigorous', [])}; "
            f"proxy={sanity.get('excluded_by_proxy', [])}"
        ),
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def _completed_summary(path: Path, *, expected_hash: str, expected_draws: int) -> dict[str, Any] | None:
    if not path.is_file():
        return None
    try:
        summary = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    if not summary.get("complete"):
        return None
    if summary.get("config_hash") != expected_hash:
        return None
    if int(summary.get("n_requested", -1)) != int(expected_draws):
        return None
    if int(summary.get("n_rows", -1)) != int(expected_draws):
        return None
    return summary


def _config_payload(cfg: ScanConfig) -> dict[str, Any]:
    payload = asdict(cfg)
    payload["mkk_values_gev"] = list(cfg.mkk_values_gev)
    if not cfg.quark_only:
        payload.pop("quark_only", None)
    return payload


def _config_from_payload(payload: Mapping[str, Any]) -> ScanConfig:
    values = dict(payload)
    values["mkk_values_gev"] = tuple(float(x) for x in values["mkk_values_gev"])
    return ScanConfig(**values)


def _config_hash(cfg: ScanConfig) -> str:
    blob = json.dumps(_config_payload(cfg), sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(blob).hexdigest()[:16]


def _resolve_git_sha() -> str:
    try:
        out = subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            cwd=str(_REPO_ROOT),
            stderr=subprocess.DEVNULL,
        )
    except (subprocess.CalledProcessError, FileNotFoundError, OSError):
        return "unknown"
    return out.decode().strip() or "unknown"


def _git_dirty() -> bool:
    try:
        out = subprocess.check_output(
            ["git", "status", "--porcelain"],
            cwd=str(_REPO_ROOT),
            stderr=subprocess.DEVNULL,
        )
    except (subprocess.CalledProcessError, FileNotFoundError, OSError):
        return True
    return bool(out.decode().strip())


def _parse_csv_floats(text: str) -> tuple[float, ...]:
    values = tuple(float(part.strip()) for part in text.split(",") if part.strip())
    if not values:
        raise argparse.ArgumentTypeError("expected at least one comma-separated value")
    return values


def _build_tiles(cfg: ScanConfig) -> list[TileSpec]:
    out: list[TileSpec] = []
    for tile_id, mkk in enumerate(cfg.mkk_values_gev):
        if mkk <= 0.0:
            raise ValueError("M_KK values must be positive")
        lambda_ir = float(mkk / cfg.xi_kk)
        out.append(
            TileSpec(
                tile_id=tile_id,
                mkk_gev=float(mkk),
                lambda_ir_gev=lambda_ir,
                n_draws=cfg.n_draws_per_tile,
                seed=int(cfg.base_seed + cfg.tile_seed_stride * tile_id),
            )
        )
    return out


def _build_argparser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run a full flavor-catalog RS scan.")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--n-draws", type=int, default=100)
    parser.add_argument("--m-kk-tev", type=str, default="1,3,5,10")
    parser.add_argument("--quark-only", action="store_true")
    parser.add_argument("--xi-kk", type=float, default=DEFAULT_XI_KK)
    parser.add_argument("--k-gev", type=float, default=MPL)
    parser.add_argument("--v-gev", type=float, default=DEFAULT_V_GEV)
    parser.add_argument("--base-seed", type=int, default=DEFAULT_BASE_SEED)
    parser.add_argument("--tile-seed-stride", type=int, default=DEFAULT_TILE_SEED_STRIDE)
    parser.add_argument("--n-workers", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", "1")))
    parser.add_argument("--no-resume", action="store_true")
    parser.add_argument("--skip-sanity", action="store_true")
    parser.add_argument("--smoke-report", type=str, default=None)
    parser.add_argument("--smoke-report-md", type=str, default=None)
    parser.add_argument("--y-half-range", type=float, default=DEFAULT_Y_HALF_RANGE)
    parser.add_argument("--y-floor", type=float, default=DEFAULT_Y_FLOOR)
    parser.add_argument("--y-prior", choices=("uniform", "gaussian"), default=DEFAULT_Y_PRIOR)
    parser.add_argument("--y-sigma", type=float, default=DEFAULT_Y_SIGMA)
    parser.add_argument("--y-trunc-sigma", type=float, default=DEFAULT_Y_TRUNC_SIGMA)
    parser.add_argument("--c-min", type=float, default=DEFAULT_C_MIN)
    parser.add_argument("--c-max", type=float, default=DEFAULT_C_MAX)
    parser.add_argument("--m-n-min-gev", type=float, default=DEFAULT_M_N_MIN_GEV)
    parser.add_argument("--m-n-max-gev", type=float, default=DEFAULT_M_N_MAX_GEV)
    parser.add_argument("--lightest-nu-min-ev", type=float, default=DEFAULT_LIGHTEST_NU_MIN_EV)
    parser.add_argument("--lightest-nu-max-ev", type=float, default=DEFAULT_LIGHTEST_NU_MAX_EV)
    parser.add_argument("--normal-ordering-probability", type=float, default=0.5)
    parser.add_argument("--perturbative-ybar-max", type=float, default=DEFAULT_PERTURBATIVE_YBAR_MAX)
    parser.add_argument("--quark-y-abs-max", type=float, default=DEFAULT_QUARK_Y_ABS_MAX)
    parser.add_argument("--quark-fit-r", type=float, default=DEFAULT_QUARK_FIT_R)
    parser.add_argument("--quark-fit-max-nfev", type=int, default=DEFAULT_QUARK_FIT_MAX_NFEV)
    parser.add_argument("--max-fit-score", type=float, default=DEFAULT_MAX_FIT_SCORE)
    parser.add_argument("--allow-unsuccessful-fit", action="store_true")
    parser.add_argument("--n-gauge-modes", type=int, default=DEFAULT_N_GAUGE_MODES)
    parser.add_argument("--quadrature-order", type=int, default=DEFAULT_QUADRATURE_ORDER)
    parser.add_argument("--min-overlap-modes", type=int, default=DEFAULT_MIN_TRUNCATION_MODES)
    parser.add_argument("--max-overlap-modes", type=int, default=DEFAULT_MAX_TRUNCATION_MODES)
    parser.add_argument("--overlap-rel-tol", type=float, default=DEFAULT_OVERLAP_RTOL)
    parser.add_argument("--spline-grid-size", type=int, default=DEFAULT_SPLINE_GRID_SIZE)
    parser.add_argument("--spline-verify-points", type=int, default=DEFAULT_SPLINE_VERIFY_POINTS)
    parser.add_argument("--c-half-atol", type=float, default=DEFAULT_C_HALF_ATOL)
    parser.add_argument("--sanity-mkk-tev", type=float, default=50.0)
    parser.add_argument("--sanity-c", type=float, default=0.4)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _build_argparser().parse_args(argv)
    mkk_values_gev = tuple(1000.0 * x for x in _parse_csv_floats(args.m_kk_tev))
    cfg = ScanConfig(
        mkk_values_gev=mkk_values_gev,
        n_draws_per_tile=int(args.n_draws),
        quark_only=bool(args.quark_only),
        xi_kk=float(args.xi_kk),
        k_gev=float(args.k_gev),
        v_gev=float(args.v_gev),
        base_seed=int(args.base_seed),
        tile_seed_stride=int(args.tile_seed_stride),
        y_half_range=float(args.y_half_range),
        y_floor=float(args.y_floor),
        y_prior=str(args.y_prior),
        y_sigma=float(args.y_sigma),
        y_trunc_sigma=float(args.y_trunc_sigma),
        c_min=float(args.c_min),
        c_max=float(args.c_max),
        m_n_min_gev=float(args.m_n_min_gev),
        m_n_max_gev=float(args.m_n_max_gev),
        lightest_nu_min_ev=float(args.lightest_nu_min_ev),
        lightest_nu_max_ev=float(args.lightest_nu_max_ev),
        normal_ordering_probability=float(args.normal_ordering_probability),
        perturbative_ybar_max=float(args.perturbative_ybar_max),
        quark_y_abs_max=float(args.quark_y_abs_max),
        quark_fit_r=float(args.quark_fit_r),
        quark_fit_max_nfev=int(args.quark_fit_max_nfev),
        max_fit_score=float(args.max_fit_score),
        require_fit_success=not bool(args.allow_unsuccessful_fit),
        n_gauge_modes=int(args.n_gauge_modes),
        quadrature_order=int(args.quadrature_order),
        min_overlap_modes=int(args.min_overlap_modes),
        max_overlap_modes=int(args.max_overlap_modes),
        overlap_rel_tol=float(args.overlap_rel_tol),
        spline_grid_size=int(args.spline_grid_size),
        spline_verify_points=int(args.spline_verify_points),
        c_half_atol=float(args.c_half_atol),
        sanity_mkk_gev=float(args.sanity_mkk_tev) * 1000.0,
        sanity_c=float(args.sanity_c),
    )
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    config_hash = _config_hash(cfg)
    provenance = {
        "git_sha": _resolve_git_sha(),
        "dirty": _git_dirty(),
        "schema": "full_catalog_scan_w6b_row_v1",
    }
    if cfg.quark_only:
        provenance.update(
            {
                "schema": "quark_only_bucket1_scan_wq_row_v1",
                "mode": "quark_only",
                "allowlist": list(QUARK_ONLY_ALLOWLIST_IDS),
                "candidate_allowlist": list(QUARK_ONLY_CANDIDATE_IDS),
                "lepton_sector": QUARK_ONLY_LEPTON_SECTOR_LABEL,
                "deferred_scope_tag": QUARK_ONLY_DEFERRED_SCOPE_TAG,
                "deferred_lepton_followup": list(QUARK_ONLY_DEFERRED_LEPTON_FOLLOWUP),
                "build_include_flags": dict(QUARK_ONLY_BUILD_INCLUDE_FLAGS),
            }
        )
    log_prefix = "[quark-only]" if cfg.quark_only else "[full-catalog]"
    tiles = _build_tiles(cfg)
    pending: list[TileSpec] = []
    summaries: list[dict[str, Any]] = []
    for tile in tiles:
        summary_path = output_dir / f"tile-{tile.tile_id:05d}.summary.json"
        completed = None
        if not args.no_resume:
            completed = _completed_summary(
                summary_path,
                expected_hash=config_hash,
                expected_draws=tile.n_draws,
            )
        if completed is None:
            pending.append(tile)
        else:
            summaries.append(completed)
            print(f"{log_prefix} resume skip tile {tile.tile_id:05d}")

    cfg_payload = _config_payload(cfg)
    worker_payloads = [
        {
            "tile_id": tile.tile_id,
            "mkk_gev": tile.mkk_gev,
            "lambda_ir_gev": tile.lambda_ir_gev,
            "n_draws": tile.n_draws,
            "seed": tile.seed,
            "output_dir": str(output_dir),
            "config_hash": config_hash,
        }
        for tile in pending
    ]
    print(f"{log_prefix} output_dir={output_dir}")
    print(f"{log_prefix} tiles={len(tiles)} pending={len(pending)} n_draws/tile={cfg.n_draws_per_tile}")
    print(f"{log_prefix} registry_expected={cfg.expected_registry_count} workers={args.n_workers}")
    start = time.perf_counter()
    if worker_payloads:
        if int(args.n_workers) <= 1:
            _worker_init(cfg_payload, provenance)
            for payload in worker_payloads:
                summary = _worker_run_tile(payload)
                summaries.append(summary)
                print(_tile_status_line(summary, prefix=log_prefix))
        else:
            with mp.Pool(
                int(args.n_workers),
                initializer=_worker_init,
                initargs=(cfg_payload, provenance),
            ) as pool:
                for summary in pool.imap_unordered(_worker_run_tile, worker_payloads):
                    summaries.append(summary)
                    print(_tile_status_line(summary, prefix=log_prefix))

    sanity = None
    if cfg.quark_only:
        sanity = {
            "name": "universal_c_diagonal_leptons",
            "skipped": True,
            "reason": "quark_only mode drops the swept lepton sector",
        }
        print(f"{log_prefix} universal-c sanity skipped for quark-only mode")
    elif not args.skip_sanity:
        sanity = run_universal_c_sanity(cfg, provenance=provenance, config_hash=config_hash)
        print(
            f"{log_prefix} universal-c sanity "
            f"passes={sanity['passes_no_spurious_hard_exclusions']}"
        )
    run_summary = _build_run_summary(
        sorted(summaries, key=lambda item: int(item["tile_id"])),
        cfg=cfg,
        config_hash=config_hash,
        provenance={**provenance, "elapsed_wall_seconds": float(time.perf_counter() - start)},
        sanity=sanity,
    )
    summary_path = output_dir / "run_summary.json"
    summary_path.write_text(json.dumps(run_summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if args.smoke_report:
        Path(args.smoke_report).write_text(
            json.dumps(run_summary, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    if args.smoke_report_md:
        _write_markdown_report(Path(args.smoke_report_md), run_summary)
    print(f"{log_prefix} run_summary={summary_path}")
    return 0


def _tile_status_line(summary: Mapping[str, Any], *, prefix: str = "[full-catalog]") -> str:
    timing = dict(summary.get("timing", {}))
    return (
        f"{prefix} tile {int(summary['tile_id']):05d} "
        f"eval={summary.get('n_evaluated_points', 0)}/{summary.get('n_requested', 0)} "
        f"skip={summary.get('n_skipped', 0)} "
        f"s/eval={_fmt_optional(timing.get('post_cache_seconds_per_evaluated_point'), 4)}"
    )


def _skip_reason(exc: Exception) -> str:
    text = str(exc)
    if text.startswith("quark_fit_failed") or text.startswith("quark_fit_score_exceeded"):
        return "quark_fit_failed"
    if text.startswith("nonperturbative_quark_yukawa"):
        return "nonperturbative_quark_yukawa"
    if text.startswith("nonperturbative_lepton_yukawa"):
        return "nonperturbative_lepton_yukawa"
    if text.startswith("c_half_singularity"):
        return "c_half_singularity"
    return "builder_or_fit_error"


def _json_sanitize(value: Any) -> Any:
    if value is None or isinstance(value, (str, bool)):
        return value
    if isinstance(value, (int, np.integer)):
        return int(value)
    if isinstance(value, (float, np.floating)):
        return _finite_or_none(float(value))
    if isinstance(value, complex):
        return {"re": float(value.real), "im": float(value.imag)}
    if isinstance(value, np.ndarray):
        return _json_sanitize(value.tolist())
    if isinstance(value, Mapping):
        return {str(k): _json_sanitize(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_sanitize(v) for v in value]
    return str(value)


def _string_or_none(value: Any) -> str | None:
    if value is None:
        return None
    if isinstance(value, str):
        return value
    try:
        return json.dumps(_json_sanitize(value), sort_keys=True, separators=(",", ":"))
    except TypeError:
        return str(value)


def _finite_or_none(value: Any) -> float | None:
    if value is None:
        return None
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if math.isfinite(out) else None


def _safe_div(numerator: float, denominator: float) -> float | None:
    denominator = float(denominator)
    if denominator == 0.0:
        return None
    return float(numerator) / denominator


def _array_real(values: Any) -> list[float]:
    return [float(x) for x in np.asarray(values, dtype=float).reshape(-1)]


def _matrix_real(values: Any) -> list[list[float]]:
    return [[float(x) for x in row] for row in np.asarray(values, dtype=float)]


def _max_abs_or_none(values: Any) -> float | None:
    arr = np.asarray(values, dtype=float)
    if arr.size == 0 or not np.all(np.isfinite(arr)):
        return None
    return float(np.max(np.abs(arr)))


def _counter_dict_with_prefix(counter: Counter[str], prefix: str) -> dict[str, int]:
    return {
        key[len(prefix) :]: int(value)
        for key, value in sorted(counter.items())
        if key.startswith(prefix)
    }


def _fmt_optional(value: Any, digits: int) -> str:
    if value is None:
        return "n/a"
    return f"{float(value):.{digits}g}"


if __name__ == "__main__":
    raise SystemExit(main())
