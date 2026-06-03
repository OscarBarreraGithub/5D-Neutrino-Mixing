"""Shared Phase-3b RS-EW test points."""

from __future__ import annotations

import math
from dataclasses import dataclass
from functools import lru_cache

import numpy as np

from flavor_catalog_constraints import point_builder
from quarkConstraints.rs_ew_couplings import DEFAULT_A_REF_C

GAUGE_ROOT_EPS_1E_MINUS_15 = 2.450509663813736
EPSILON_RS = 1.0e-15
N_GAUGE_MODES = 64
QUADRATURE_ORDER = 512
MIN_OVERLAP_MODES = 16
MAX_OVERLAP_MODES = 64
OVERLAP_REL_TOL = 1.0e-3


@dataclass(frozen=True)
class BulkState:
    c_Q: np.ndarray
    c_u: np.ndarray
    c_d: np.ndarray


@dataclass(frozen=True)
class QuarkFit:
    bulk_state: BulkState
    U_L_u: np.ndarray
    U_L_d: np.ndarray
    U_R_u: np.ndarray
    U_R_d: np.ndarray


def _scales_for_mkk(mkk_gev: float) -> tuple[float, float]:
    lambda_ir = float(mkk_gev) / GAUGE_ROOT_EPS_1E_MINUS_15
    return lambda_ir, lambda_ir / EPSILON_RS


def _rot12(theta: float) -> np.ndarray:
    c = math.cos(theta)
    s = math.sin(theta)
    return np.array([[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]], dtype=np.complex128)


def _rot23(theta: float) -> np.ndarray:
    c = math.cos(theta)
    s = math.sin(theta)
    return np.array([[1.0, 0.0, 0.0], [0.0, c, s], [0.0, -s, c]], dtype=np.complex128)


def _rot13(theta: float, phase: float) -> np.ndarray:
    c = math.cos(theta)
    s = math.sin(theta)
    e_ip = complex(math.cos(phase), math.sin(phase))
    return np.array(
        [[c, 0.0, s * np.conjugate(e_ip)], [0.0, 1.0, 0.0], [-s * e_ip, 0.0, c]],
        dtype=np.complex128,
    )


def _sample_fit() -> QuarkFit:
    return QuarkFit(
        bulk_state=BulkState(
            c_Q=np.array([0.64, 0.56, 0.43], dtype=float),
            c_u=np.array([0.62, 0.34, 0.18], dtype=float),
            c_d=np.array([0.66, 0.57, 0.20], dtype=float),
        ),
        U_L_u=_rot12(0.11) @ _rot23(-0.07) @ _rot13(0.025, 0.4),
        U_L_d=_rot12(-0.19) @ _rot23(0.08) @ _rot13(-0.015, -0.2),
        U_R_u=_rot12(-0.09) @ _rot23(0.045),
        U_R_d=_rot12(0.05) @ _rot23(-0.12),
    )


def _identity_fit(c_q: np.ndarray, c_u: np.ndarray, c_d: np.ndarray) -> QuarkFit:
    identity = np.eye(3, dtype=np.complex128)
    return QuarkFit(
        bulk_state=BulkState(
            c_Q=np.asarray(c_q, dtype=float),
            c_u=np.asarray(c_u, dtype=float),
            c_d=np.asarray(c_d, dtype=float),
        ),
        U_L_u=identity,
        U_L_d=identity,
        U_R_u=identity,
        U_R_d=identity,
    )


def _build_point(fit: QuarkFit, *, mkk_gev: float = 3000.0):
    lambda_ir, k = _scales_for_mkk(mkk_gev)
    return point_builder.build_from_rs_ew_inputs(
        fit,
        Lambda_IR=lambda_ir,
        k=k,
        n_gauge_modes=N_GAUGE_MODES,
        quadrature_order=QUADRATURE_ORDER,
        min_overlap_modes=MIN_OVERLAP_MODES,
        max_overlap_modes=MAX_OVERLAP_MODES,
        overlap_rel_tol=OVERLAP_REL_TOL,
    )


@lru_cache(maxsize=None)
def sample_rs_ew_point():
    return _build_point(_sample_fit())


@lru_cache(maxsize=None)
def sm_limit_rs_ew_point():
    ref = np.array([DEFAULT_A_REF_C, DEFAULT_A_REF_C, DEFAULT_A_REF_C], dtype=float)
    return _build_point(_identity_fit(ref, ref, ref))
