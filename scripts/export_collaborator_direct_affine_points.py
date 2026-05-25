#!/usr/bin/env python3
"""Export collaborator points with direct affine Yukawa-derived bulk masses.

This exporter is intentionally separate from ``export_collaborator_5tev_points``.
It does not use ``BulkMassMap``.  The bulk masses are eigenvalues of explicit
FPR-like affine spurion matrices,

    C_Q = alpha_Q I + beta_Q Y_d Y_d^dagger + gamma_Q Y_u Y_u^dagger
    C_u = alpha_u I + beta_u Y_u^dagger Y_u
    C_d = alpha_d I + beta_d Y_d^dagger Y_d

and the coefficients are included in every CSV row.
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from types import SimpleNamespace
from typing import Iterable

import numpy as np
from scipy.optimize import least_squares

REPO_ROOT = Path(__file__).resolve().parents[1]

import sys

if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from warpConfig.baseParams import MPL, get_warp_params  # noqa: E402
from warpConfig.wavefuncs import f_IR  # noqa: E402

from quarkConstraints.benchmarks import default_quark_targets, default_spurion_seed  # noqa: E402
from quarkConstraints.fit import (  # noqa: E402
    _canonicalize_fit_seed,
    _decode_seed,
    _encode_seed,
    build_mass_matrices,
    ckm_observables,
    fit_residuals,
    mass_matrix_observables,
)
from quarkConstraints.model import (  # noqa: E402
    QuarkSpurionPoint,
    build_mfv_point_from_singular_values,
)
from quarkConstraints.modern.artifacts import build_modern_point_bridge_artifact  # noqa: E402
from quarkConstraints.modern.couplings import build_modern_point_couplings  # noqa: E402
from quarkConstraints.modern.evaluation import MODERN_POINT_EVALUATION_SCHEMA_ID  # noqa: E402
from quarkConstraints.modern.matching import build_modern_point_matching  # noqa: E402
from quarkConstraints.modern.phenomenology import (  # noqa: E402
    build_modern_point_phenomenology_artifact,
)
from quarkConstraints.scales import (  # noqa: E402
    DEFAULT_QUARK_FIT_SCALE_GEV,
    DEFAULT_QUARK_TARGET_SCALE_GEV,
    GAUGE_KK_ROOT_NN,
)


DEFAULT_OUTPUT = REPO_ROOT / "artifacts" / "collaborator_direct_affine_5_10tev_points.csv"
DEFAULT_PROVENANCE = (
    REPO_ROOT / "artifacts" / "collaborator_direct_affine_5_10tev_points.provenance.json"
)
SYSTEM_ORDER = ("epsilon_K", "K", "B_d", "B_s", "D0")
GEN_LABELS = ("1", "2", "3")
QUARK_UP_LABELS = ("u", "c", "t")
QUARK_DOWN_LABELS = ("d", "s", "b")


@dataclass(frozen=True)
class DirectAffinePolicy:
    name: str
    alpha_Q: float
    beta_Q: float
    gamma_Q: float
    alpha_u: float
    beta_u: float
    alpha_d: float
    beta_d: float
    note: str

    @property
    def r(self) -> float:
        return float(self.gamma_Q / self.beta_Q)


@dataclass(frozen=True)
class DirectAffineBulkState:
    point: QuarkSpurionPoint
    policy: DirectAffinePolicy
    epsilon: float
    C_Q: np.ndarray
    C_u: np.ndarray
    C_d: np.ndarray
    c_Q: np.ndarray
    c_u: np.ndarray
    c_d: np.ndarray
    F_Q: np.ndarray
    F_u: np.ndarray
    F_d: np.ndarray
    rotation_Q: np.ndarray
    rotation_u: np.ndarray
    rotation_d: np.ndarray
    eig_Q: np.ndarray
    eig_u: np.ndarray
    eig_d: np.ndarray
    Y_u_bulk_basis: np.ndarray
    Y_d_bulk_basis: np.ndarray


@dataclass(frozen=True)
class DirectAffineFitResult:
    bulk_state: DirectAffineBulkState
    M_u: np.ndarray
    M_d: np.ndarray
    U_L_u: np.ndarray
    U_R_u: np.ndarray
    U_L_d: np.ndarray
    U_R_d: np.ndarray
    masses_up: np.ndarray
    masses_down: np.ndarray
    ckm: np.ndarray
    up_log_residuals: np.ndarray
    down_log_residuals: np.ndarray
    ckm_observable_residuals: np.ndarray
    score: float

    @property
    def state(self) -> DirectAffineBulkState:
        return self.bulk_state

    @property
    def point(self) -> QuarkSpurionPoint:
        return self.bulk_state.point

    @property
    def mass_residuals_up(self) -> np.ndarray:
        return self.up_log_residuals

    @property
    def mass_residuals_down(self) -> np.ndarray:
        return self.down_log_residuals

    @property
    def ckm_residuals(self) -> np.ndarray:
        return self.ckm_observable_residuals

    @property
    def ckm_matrix(self) -> np.ndarray:
        return self.ckm

    @property
    def ckm_observables(self) -> np.ndarray:
        return ckm_observables(self.ckm)

    @property
    def residual_norm(self) -> float:
        return float(
            np.linalg.norm(
                np.concatenate(
                    [
                        self.up_log_residuals,
                        self.down_log_residuals,
                        self.ckm_observable_residuals,
                    ]
                )
            )
        )


@dataclass(frozen=True)
class CandidateResult:
    policy: DirectAffinePolicy
    m_gkk_GeV: float
    seed: object
    fit_result: DirectAffineFitResult
    phenomenology: object
    matching: object
    fit_residual_norm: float
    optimizer_success: bool
    optimizer_message: str
    optimizer_nfev: int

    @property
    def max_abs_Y_bulk(self) -> float:
        state = self.fit_result.state
        return float(
            max(
                np.max(np.abs(state.Y_u_bulk_basis)),
                np.max(np.abs(state.Y_d_bulk_basis)),
            )
        )

    @property
    def max_abs_Y_flavor(self) -> float:
        point = self.fit_result.point
        return float(max(np.max(np.abs(point.Y_u)), np.max(np.abs(point.Y_d))))

    @property
    def ratios(self) -> dict[str, float]:
        return {
            result.system_id: float(result.ratio_to_bound)
            for result in self.phenomenology.system_results
        }

    @property
    def bounds(self) -> dict[str, float]:
        return {
            result.system_id: float(result.bound)
            for result in self.phenomenology.system_results
        }

    @property
    def max_ratio_to_bound(self) -> float:
        return float(max(self.ratios.values()))

    @property
    def dominant_constraint(self) -> str:
        ratios = self.ratios
        return max(ratios, key=ratios.get)


def default_policies() -> tuple[DirectAffinePolicy, ...]:
    return (
        DirectAffinePolicy(
            name="down_aligned_flat_d",
            alpha_Q=0.64,
            beta_Q=-0.10,
            gamma_Q=-0.003,
            alpha_u=0.68,
            beta_u=-0.12,
            alpha_d=0.64,
            beta_d=-0.04,
            note="Small up-spurion admixture in C_Q; flatter down singlet sector.",
        ),
        DirectAffinePolicy(
            name="down_aligned_flatter_u",
            alpha_Q=0.64,
            beta_Q=-0.10,
            gamma_Q=-0.003,
            alpha_u=0.70,
            beta_u=-0.09,
            alpha_d=0.64,
            beta_d=-0.04,
            note="Same C_Q policy with a milder direct up-singlet slope.",
        ),
        DirectAffinePolicy(
            name="balanced_flat",
            alpha_Q=0.66,
            beta_Q=-0.08,
            gamma_Q=-0.004,
            alpha_u=0.70,
            beta_u=-0.08,
            alpha_d=0.64,
            beta_d=-0.04,
            note="Balanced direct-affine slopes with max |Y_bulk| close to 3.",
        ),
        DirectAffinePolicy(
            name="stronger_down_alignment",
            alpha_Q=0.64,
            beta_Q=-0.20,
            gamma_Q=-0.006,
            alpha_u=0.68,
            beta_u=-0.12,
            alpha_d=0.65,
            beta_d=-0.10,
            note="Larger direct down-spurion coefficient with the same gamma/beta ratio.",
        ),
        DirectAffinePolicy(
            name="moderate_up_admixture",
            alpha_Q=0.64,
            beta_Q=-0.10,
            gamma_Q=-0.010,
            alpha_u=0.68,
            beta_u=-0.12,
            alpha_d=0.65,
            beta_d=-0.10,
            note="Moderate C_Q up-spurion admixture; included only if flavor-safe.",
        ),
        DirectAffinePolicy(
            name="larger_Q_slope_high_scale",
            alpha_Q=0.66,
            beta_Q=-0.50,
            gamma_Q=-0.015,
            alpha_u=0.70,
            beta_u=-0.08,
            alpha_d=0.64,
            beta_d=-0.04,
            note="Larger Q-sector down slope; usually viable only at the high end of M_gkk.",
        ),
    )


def _format_float(value: float) -> str:
    return f"{float(value):.17g}"


def _format_bool(value: bool) -> str:
    return "true" if value else "false"


def _ordered_hermitian_spectrum(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    hermitian = 0.5 * (matrix + matrix.conjugate().T)
    values, vectors = np.linalg.eigh(hermitian)
    order = np.argsort(values.real, kind="stable")
    return values[order].real.astype(float), vectors[:, order].astype(np.complex128)


def _build_point(seed: object, policy: DirectAffinePolicy, m_gkk_GeV: float) -> QuarkSpurionPoint:
    lambda_ir = float(m_gkk_GeV / GAUGE_KK_ROOT_NN)
    return build_mfv_point_from_singular_values(
        up_singular_values=seed.up_singular_values,
        down_singular_values=seed.down_singular_values,
        overall_scale=seed.overall_scale,
        r=policy.r,
        up_left=seed.up_left,
        up_right=seed.up_right,
        down_left=seed.down_left,
        down_right=seed.down_right,
        Lambda_IR=lambda_ir,
        k=MPL,
        label=f"{policy.name}_{m_gkk_GeV / 1000.0:.1f}tev",
        metadata={
            "bulk_mass_derivation": "explicit_direct_affine_yukawa_spurion",
            "hidden_BulkMassMap_used": False,
            "nonlinear_eigenvalue_map_used": False,
            "physical_m_gkk_GeV": float(m_gkk_GeV),
            "xi_kk": GAUGE_KK_ROOT_NN,
        },
    )


def _derive_direct_affine_state(
    point: QuarkSpurionPoint,
    policy: DirectAffinePolicy,
) -> DirectAffineBulkState:
    identity = np.eye(3, dtype=np.complex128)
    Y_u = point.Y_u
    Y_d = point.Y_d
    C_Q = (
        policy.alpha_Q * identity
        + policy.beta_Q * (Y_d @ Y_d.conjugate().T)
        + policy.gamma_Q * (Y_u @ Y_u.conjugate().T)
    )
    C_u = policy.alpha_u * identity + policy.beta_u * (Y_u.conjugate().T @ Y_u)
    C_d = policy.alpha_d * identity + policy.beta_d * (Y_d.conjugate().T @ Y_d)

    eig_Q, rotation_Q = _ordered_hermitian_spectrum(C_Q)
    eig_u, rotation_u = _ordered_hermitian_spectrum(C_u)
    eig_d, rotation_d = _ordered_hermitian_spectrum(C_d)

    epsilon = float(get_warp_params(k=point.k, Lambda_IR=point.Lambda_IR)["epsilon"])
    F_Q = np.asarray(f_IR(eig_Q, epsilon), dtype=float)
    F_u = np.asarray(f_IR(eig_u, epsilon), dtype=float)
    F_d = np.asarray(f_IR(eig_d, epsilon), dtype=float)
    for name, arr in (("F_Q", F_Q), ("F_u", F_u), ("F_d", F_d)):
        if not np.all(np.isfinite(arr)) or np.any(arr <= 1.0e-14):
            raise ValueError(f"{name} contains invalid direct-affine profile values")

    return DirectAffineBulkState(
        point=point,
        policy=policy,
        epsilon=epsilon,
        C_Q=C_Q,
        C_u=C_u,
        C_d=C_d,
        c_Q=eig_Q,
        c_u=eig_u,
        c_d=eig_d,
        F_Q=F_Q,
        F_u=F_u,
        F_d=F_d,
        rotation_Q=rotation_Q,
        rotation_u=rotation_u,
        rotation_d=rotation_d,
        eig_Q=eig_Q,
        eig_u=eig_u,
        eig_d=eig_d,
        Y_u_bulk_basis=rotation_Q.conjugate().T @ Y_u @ rotation_u,
        Y_d_bulk_basis=rotation_Q.conjugate().T @ Y_d @ rotation_d,
    )


def _evaluate_direct_affine_fit(
    seed: object,
    policy: DirectAffinePolicy,
    m_gkk_GeV: float,
) -> DirectAffineFitResult:
    point = _build_point(seed, policy, m_gkk_GeV)
    state = _derive_direct_affine_state(point, policy)
    M_u, M_d = build_mass_matrices(state)
    observables = mass_matrix_observables(M_u, M_d)
    targets = default_quark_targets()
    residuals = fit_residuals(
        observables["masses_up"],
        observables["masses_down"],
        observables["ckm"],
        targets,
    )
    return DirectAffineFitResult(
        bulk_state=state,
        M_u=M_u,
        M_d=M_d,
        U_L_u=observables["U_L_u"],
        U_R_u=observables["U_R_u"],
        U_L_d=observables["U_L_d"],
        U_R_d=observables["U_R_d"],
        masses_up=observables["masses_up"],
        masses_down=observables["masses_down"],
        ckm=observables["ckm"],
        up_log_residuals=np.asarray(residuals["up_log_residuals"], dtype=float),
        down_log_residuals=np.asarray(residuals["down_log_residuals"], dtype=float),
        ckm_observable_residuals=np.asarray(
            residuals["ckm_observable_residuals"],
            dtype=float,
        ),
        score=float(residuals["total_score"]),
    )


def _fit_candidate(
    policy: DirectAffinePolicy,
    m_gkk_GeV: float,
    *,
    x0: np.ndarray,
    max_nfev: int,
) -> tuple[CandidateResult, np.ndarray]:
    seed_template = _canonicalize_fit_seed(default_spurion_seed(), overall_scale=None)

    def residual_vector(vector: np.ndarray) -> np.ndarray:
        try:
            candidate_seed = _decode_seed(vector, seed_template, True)
            result = _evaluate_direct_affine_fit(candidate_seed, policy, m_gkk_GeV)
            vector_out = np.concatenate(
                [
                    result.mass_residuals_up,
                    result.mass_residuals_down,
                    result.ckm_residuals,
                ]
            )
            if not np.all(np.isfinite(vector_out)):
                return np.full(10, 1.0e4)
            return vector_out
        except (ValueError, np.linalg.LinAlgError, FloatingPointError):
            return np.full(10, 1.0e4)

    optimum = least_squares(
        residual_vector,
        x0=x0,
        max_nfev=max_nfev,
        ftol=1.0e-10,
        xtol=1.0e-10,
        gtol=1.0e-10,
    )
    best_seed = _decode_seed(optimum.x, seed_template, True)
    fit_result = _evaluate_direct_affine_fit(best_seed, policy, m_gkk_GeV)
    phenomenology, matching = _physical_mgkk_phenomenology(
        fit_result,
        point_label=fit_result.point.label,
        m_gkk_GeV=m_gkk_GeV,
    )
    candidate = CandidateResult(
        policy=policy,
        m_gkk_GeV=float(m_gkk_GeV),
        seed=best_seed,
        fit_result=fit_result,
        phenomenology=phenomenology,
        matching=matching,
        fit_residual_norm=float(np.linalg.norm(residual_vector(optimum.x))),
        optimizer_success=bool(optimum.success),
        optimizer_message=str(optimum.message),
        optimizer_nfev=int(optimum.nfev),
    )
    return candidate, optimum.x


def _physical_mgkk_phenomenology(
    result: DirectAffineFitResult,
    *,
    point_label: str,
    m_gkk_GeV: float,
):
    couplings = build_modern_point_couplings(
        result,
        point_id=point_label,
        point_label=point_label,
        M_KK=m_gkk_GeV,
        g_s_star=3.0,
    )
    matching = build_modern_point_matching(couplings)
    source = SimpleNamespace(
        schema_id=MODERN_POINT_EVALUATION_SCHEMA_ID,
        point_id=point_label,
        point_label=point_label,
        input_bundle_schema_id=couplings.input_bundle_schema_id,
        input_bundle_id=couplings.input_bundle_id,
        input_provenance_id=couplings.input_provenance_id,
        input_resolution_policy_id=couplings.input_resolution_policy_id,
        qcd_metadata_id=couplings.qcd_metadata_id,
        alpha_s_policy_id=couplings.alpha_s_policy_id,
        coupling_schema_id=couplings.schema_id,
        matching_schema_id=matching.schema_id,
        couplings=couplings.as_dict(),
        matching=matching.as_dict(),
    )
    bridge = build_modern_point_bridge_artifact(source)
    return build_modern_point_phenomenology_artifact(bridge), matching


def scan_candidates(
    *,
    masses_GeV: Iterable[float],
    max_nfev: int,
) -> list[CandidateResult]:
    seed_template = _canonicalize_fit_seed(default_spurion_seed(), overall_scale=None)
    candidates: list[CandidateResult] = []
    for policy in default_policies():
        x0 = _encode_seed(seed_template, True)
        for m_gkk_GeV in masses_GeV:
            try:
                candidate, x0 = _fit_candidate(
                    policy,
                    float(m_gkk_GeV),
                    x0=x0,
                    max_nfev=max_nfev,
                )
                candidates.append(candidate)
            except (ValueError, np.linalg.LinAlgError, FloatingPointError):
                x0 = _encode_seed(seed_template, True)
    return candidates


def select_points(
    candidates: list[CandidateResult],
    *,
    count: int,
    max_y_bulk: float,
    max_fit_residual_norm: float,
) -> list[CandidateResult]:
    viable = [
        item
        for item in candidates
        if item.fit_residual_norm <= max_fit_residual_norm
        and item.max_abs_Y_bulk <= max_y_bulk
        and item.phenomenology.non_cp_passes
        and item.max_ratio_to_bound <= 1.0
        and 5000.0 <= item.m_gkk_GeV <= 10000.0
    ]

    preferred_slots = (
        ("down_aligned_flat_d", 5000.0),
        ("down_aligned_flatter_u", 6000.0),
        ("stronger_down_alignment", 7000.0),
        ("balanced_flat", 8000.0),
        ("larger_Q_slope_high_scale", 10000.0),
    )
    selected: list[CandidateResult] = []
    used_keys: set[tuple[str, float]] = set()
    by_key = {
        (item.policy.name, round(item.m_gkk_GeV, 6)): item
        for item in viable
    }
    for policy_name, mass in preferred_slots:
        key = (policy_name, round(mass, 6))
        if key in by_key:
            selected.append(by_key[key])
            used_keys.add(key)
    if len(selected) < count:
        for item in sorted(
            viable,
            key=lambda x: (
                x.max_ratio_to_bound,
                x.max_abs_Y_bulk,
                x.m_gkk_GeV,
                x.policy.name,
            ),
        ):
            key = (item.policy.name, round(item.m_gkk_GeV, 6))
            if key in used_keys:
                continue
            selected.append(item)
            used_keys.add(key)
            if len(selected) == count:
                break
    if len(selected) < count:
        raise RuntimeError(
            f"only found {len(selected)} viable direct-affine points; requested {count}"
        )
    return selected[:count]


def _add_vector(row: dict[str, str], prefix: str, values: np.ndarray, labels: Iterable[str]) -> None:
    for label, value in zip(labels, np.asarray(values), strict=True):
        row[f"{prefix}_{label}"] = _format_float(float(value))


def _add_complex_matrix(
    row: dict[str, str],
    prefix: str,
    matrix: np.ndarray,
    row_labels: Iterable[str],
    col_labels: Iterable[str],
    *,
    include_abs: bool = False,
) -> None:
    arr = np.asarray(matrix, dtype=np.complex128)
    for i, row_label in enumerate(row_labels):
        for j, col_label in enumerate(col_labels):
            value = arr[i, j]
            stem = f"{prefix}_{row_label}{col_label}"
            row[f"{stem}_re"] = _format_float(float(value.real))
            row[f"{stem}_im"] = _format_float(float(value.imag))
            if include_abs:
                row[f"{stem}_abs"] = _format_float(float(abs(value)))


def _row_for(candidate: CandidateResult) -> dict[str, str]:
    result = candidate.fit_result
    state = result.state
    policy = candidate.policy
    targets = default_quark_targets()
    ratios = candidate.ratios
    bounds = candidate.bounds
    row: dict[str, str] = {
        "coefficient_policy": policy.name,
        "m_gkk_TeV": _format_float(candidate.m_gkk_GeV / 1000.0),
        "m_gkk_GeV": _format_float(candidate.m_gkk_GeV),
        "Lambda_IR_GeV": _format_float(result.point.Lambda_IR),
        "xi_kk": _format_float(GAUGE_KK_ROOT_NN),
        "g_s_star": _format_float(3.0),
        "mass_target_scale_GeV": _format_float(DEFAULT_QUARK_FIT_SCALE_GEV),
        "wilson_reference_scale_GeV": _format_float(DEFAULT_QUARK_TARGET_SCALE_GEV),
        "constraint_propagator_mass_GeV": _format_float(candidate.m_gkk_GeV),
        "bulk_mass_derivation": "explicit_direct_affine_yukawa_spurion",
        "hidden_BulkMassMap_used": _format_bool(False),
        "nonlinear_eigenvalue_map_used": _format_bool(False),
        "bulk_mass_clipping_used": _format_bool(False),
        "C_Q_formula": "alpha_Q*I + beta_Q*Y_d*Y_d^dagger + gamma_Q*Y_u*Y_u^dagger",
        "C_u_formula": "alpha_u*I + beta_u*Y_u^dagger*Y_u",
        "C_d_formula": "alpha_d*I + beta_d*Y_d^dagger*Y_d",
        "alpha_Q": _format_float(policy.alpha_Q),
        "beta_Q": _format_float(policy.beta_Q),
        "gamma_Q": _format_float(policy.gamma_Q),
        "r_gamma_over_beta": _format_float(policy.r),
        "alpha_u": _format_float(policy.alpha_u),
        "beta_u": _format_float(policy.beta_u),
        "alpha_d": _format_float(policy.alpha_d),
        "beta_d": _format_float(policy.beta_d),
        "coefficient_note": policy.note,
        "fit_residual_norm": _format_float(candidate.fit_residual_norm),
        "fit_score_rms": _format_float(result.score),
        "max_abs_mass_log_residual": _format_float(
            max(
                float(np.max(np.abs(result.mass_residuals_up))),
                float(np.max(np.abs(result.mass_residuals_down))),
            )
        ),
        "max_abs_ckm_observable_residual": _format_float(
            float(np.max(np.abs(result.ckm_residuals)))
        ),
        "max_abs_Y_flavor": _format_float(candidate.max_abs_Y_flavor),
        "max_abs_Y_bulk": _format_float(candidate.max_abs_Y_bulk),
        "accepted_current_constraints": _format_bool(candidate.phenomenology.non_cp_passes),
        "max_ratio_to_bound": _format_float(candidate.max_ratio_to_bound),
        "dominant_constraint": candidate.dominant_constraint,
        "optimizer_success": _format_bool(candidate.optimizer_success),
        "optimizer_nfev": str(candidate.optimizer_nfev),
        "optimizer_message": candidate.optimizer_message,
    }
    for system in SYSTEM_ORDER:
        row[f"ratio_{system}"] = _format_float(ratios[system])
        row[f"passes_{system}"] = _format_bool(ratios[system] <= 1.0)
        row[f"bound_{system}"] = _format_float(bounds[system])

    _add_vector(row, "target_m_up_GeV", targets.up_masses, QUARK_UP_LABELS)
    _add_vector(row, "target_m_down_GeV", targets.down_masses, QUARK_DOWN_LABELS)
    _add_vector(row, "fit_m_up_GeV", result.masses_up, QUARK_UP_LABELS)
    _add_vector(row, "fit_m_down_GeV", result.masses_down, QUARK_DOWN_LABELS)
    _add_vector(row, "mass_log_residual_up", result.mass_residuals_up, QUARK_UP_LABELS)
    _add_vector(row, "mass_log_residual_down", result.mass_residuals_down, QUARK_DOWN_LABELS)
    _add_vector(
        row,
        "ckm_observable",
        result.ckm_observables,
        ("Vus_abs", "Vcb_abs", "Vub_abs", "J"),
    )
    _add_vector(
        row,
        "target_ckm_observable",
        targets.ckm_observables,
        ("Vus_abs", "Vcb_abs", "Vub_abs", "J"),
    )
    _add_vector(
        row,
        "ckm_observable_residual",
        result.ckm_residuals,
        ("Vus_abs", "Vcb_abs", "Vub_abs", "J"),
    )
    _add_vector(row, "c_Q", state.c_Q, GEN_LABELS)
    _add_vector(row, "c_u", state.c_u, GEN_LABELS)
    _add_vector(row, "c_d", state.c_d, GEN_LABELS)
    _add_vector(row, "F_Q", state.F_Q, GEN_LABELS)
    _add_vector(row, "F_u", state.F_u, GEN_LABELS)
    _add_vector(row, "F_d", state.F_d, GEN_LABELS)
    _add_vector(row, "eig_Q", state.eig_Q, GEN_LABELS)
    _add_vector(row, "eig_u", state.eig_u, GEN_LABELS)
    _add_vector(row, "eig_d", state.eig_d, GEN_LABELS)
    _add_vector(
        row,
        "Y_u_singular_values",
        np.sort(np.linalg.svd(result.point.Y_u, compute_uv=False)),
        GEN_LABELS,
    )
    _add_vector(
        row,
        "Y_d_singular_values",
        np.sort(np.linalg.svd(result.point.Y_d, compute_uv=False)),
        GEN_LABELS,
    )

    _add_complex_matrix(row, "CKM", result.ckm, QUARK_UP_LABELS, QUARK_DOWN_LABELS, include_abs=True)
    _add_complex_matrix(
        row,
        "target_CKM",
        targets.ckm,
        QUARK_UP_LABELS,
        QUARK_DOWN_LABELS,
        include_abs=True,
    )
    _add_complex_matrix(row, "Y_u_flavor", result.point.Y_u, GEN_LABELS, GEN_LABELS, include_abs=True)
    _add_complex_matrix(row, "Y_d_flavor", result.point.Y_d, GEN_LABELS, GEN_LABELS, include_abs=True)
    _add_complex_matrix(row, "Y_u_bulk", state.Y_u_bulk_basis, GEN_LABELS, GEN_LABELS, include_abs=True)
    _add_complex_matrix(row, "Y_d_bulk", state.Y_d_bulk_basis, GEN_LABELS, GEN_LABELS, include_abs=True)
    _add_complex_matrix(row, "C_Q", state.C_Q, GEN_LABELS, GEN_LABELS)
    _add_complex_matrix(row, "C_u", state.C_u, GEN_LABELS, GEN_LABELS)
    _add_complex_matrix(row, "C_d", state.C_d, GEN_LABELS, GEN_LABELS)
    _add_complex_matrix(row, "rotation_Q", state.rotation_Q, GEN_LABELS, GEN_LABELS)
    _add_complex_matrix(row, "rotation_u", state.rotation_u, GEN_LABELS, GEN_LABELS)
    _add_complex_matrix(row, "rotation_d", state.rotation_d, GEN_LABELS, GEN_LABELS)
    _add_complex_matrix(row, "M_u_GeV", result.M_u, GEN_LABELS, GEN_LABELS)
    _add_complex_matrix(row, "M_d_GeV", result.M_d, GEN_LABELS, GEN_LABELS)
    _add_complex_matrix(row, "U_L_u", result.U_L_u, GEN_LABELS, GEN_LABELS)
    _add_complex_matrix(row, "U_R_u", result.U_R_u, GEN_LABELS, GEN_LABELS)
    _add_complex_matrix(row, "U_L_d", result.U_L_d, GEN_LABELS, GEN_LABELS)
    _add_complex_matrix(row, "U_R_d", result.U_R_d, GEN_LABELS, GEN_LABELS)
    return row


def _git_value(args: list[str]) -> str | None:
    try:
        return subprocess.check_output(
            ["git", *args],
            cwd=REPO_ROOT,
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except Exception:
        return None


def write_outputs(
    selected: list[CandidateResult],
    all_candidates: list[CandidateResult],
    *,
    output: Path,
    provenance: Path,
    args: argparse.Namespace,
) -> None:
    rows = [_row_for(item) for item in selected]
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    provenance_payload = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "output_csv": str(output.relative_to(REPO_ROOT)),
        "git_commit": _git_value(["rev-parse", "HEAD"]),
        "git_dirty_status": _git_value(["status", "--short"]),
        "selection": {
            "requested_count": int(args.count),
            "max_abs_Y_bulk": float(args.max_y_bulk),
            "max_fit_residual_norm": float(args.max_fit_residual_norm),
            "mass_grid_GeV": [float(x) for x in args.mass_grid_GeV],
        },
        "construction": {
            "bulk_mass_derivation": "explicit_direct_affine_yukawa_spurion",
            "uses_hidden_BulkMassMap": False,
            "uses_nonlinear_eigenvalue_map": False,
            "C_Q_formula": "alpha_Q I + beta_Q Y_d Y_d^dagger + gamma_Q Y_u Y_u^dagger",
            "C_u_formula": "alpha_u I + beta_u Y_u^dagger Y_u",
            "C_d_formula": "alpha_d I + beta_d Y_d^dagger Y_d",
            "xi_kk": GAUGE_KK_ROOT_NN,
            "mass_target_scale_GeV": DEFAULT_QUARK_FIT_SCALE_GEV,
            "wilson_reference_scale_GeV": DEFAULT_QUARK_TARGET_SCALE_GEV,
            "constraint_propagator_mass": "physical m_gkk for each row",
        },
        "policies": [policy.__dict__ | {"r_gamma_over_beta": policy.r} for policy in default_policies()],
        "selected": [
            {
                "coefficient_policy": item.policy.name,
                "m_gkk_GeV": item.m_gkk_GeV,
                "max_ratio_to_bound": item.max_ratio_to_bound,
                "dominant_constraint": item.dominant_constraint,
                "max_abs_Y_bulk": item.max_abs_Y_bulk,
                "fit_residual_norm": item.fit_residual_norm,
            }
            for item in selected
        ],
        "all_viable_count": sum(
            1
            for item in all_candidates
            if item.phenomenology.non_cp_passes
            and item.max_abs_Y_bulk <= args.max_y_bulk
            and item.fit_residual_norm <= args.max_fit_residual_norm
        ),
        "all_candidate_count": len(all_candidates),
    }
    provenance.parent.mkdir(parents=True, exist_ok=True)
    provenance.write_text(json.dumps(provenance_payload, indent=2, sort_keys=True) + "\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--provenance", type=Path, default=DEFAULT_PROVENANCE)
    parser.add_argument("--count", type=int, default=5)
    parser.add_argument("--max-y-bulk", type=float, default=3.0)
    parser.add_argument("--max-fit-residual-norm", type=float, default=1.0e-6)
    parser.add_argument("--max-nfev", type=int, default=900)
    parser.add_argument(
        "--mass-grid-gev",
        dest="mass_grid_GeV",
        type=float,
        nargs="+",
        default=[5000.0, 6000.0, 7000.0, 8000.0, 9000.0, 10000.0],
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    candidates = scan_candidates(masses_GeV=args.mass_grid_GeV, max_nfev=args.max_nfev)
    selected = select_points(
        candidates,
        count=args.count,
        max_y_bulk=args.max_y_bulk,
        max_fit_residual_norm=args.max_fit_residual_norm,
    )
    write_outputs(
        selected,
        candidates,
        output=args.output,
        provenance=args.provenance,
        args=args,
    )
    print(f"wrote {len(selected)} rows to {args.output}")
    for item in selected:
        print(
            f"{item.policy.name:28s} m_gkk={item.m_gkk_GeV / 1000.0:.1f} TeV "
            f"max_ratio={item.max_ratio_to_bound:.3g} "
            f"dominant={item.dominant_constraint} "
            f"max|Y_bulk|={item.max_abs_Y_bulk:.3g} "
            f"fit_norm={item.fit_residual_norm:.3g}"
        )


if __name__ == "__main__":
    main()
