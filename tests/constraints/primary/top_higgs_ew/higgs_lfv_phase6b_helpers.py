"""Shared Phase-6b Higgs-LFV test helpers."""

from __future__ import annotations

import math
from dataclasses import replace
from functools import lru_cache

import numpy as np

from flavor_catalog_constraints import point_builder
from quarkConstraints.rs_higgs_yukawas import build_rs_higgs_yukawas
from tests.rs_ew_phase3b_helpers import (
    MAX_OVERLAP_MODES,
    MIN_OVERLAP_MODES,
    N_GAUGE_MODES,
    OVERLAP_REL_TOL,
    QUADRATURE_ORDER,
    _sample_fit,
    _scales_for_mkk,
)

_LEPTON_INDEX = {"e": 0, "mu": 1, "tau": 2}


def lepton_sweep_inputs(c_l: float = 0.58):
    return {
        "c_L": float(c_l),
        "c_E": [float(c_l), float(c_l), float(c_l)],
        "c_N": 0.27,
        "M_N": 1.22e18,
        "lightest_nu_mass": 0.002,
        "ordering": "normal",
        "majorana_alpha": 0.0,
        "majorana_beta": 0.0,
    }


@lru_cache(maxsize=None)
def diagonal_higgs_point(mkk_gev: float = 3000.0):
    lambda_ir, k = _scales_for_mkk(mkk_gev)
    return point_builder.build_from_rs_ew_inputs(
        _sample_fit(),
        Lambda_IR=lambda_ir,
        k=k,
        n_gauge_modes=N_GAUGE_MODES,
        quadrature_order=QUADRATURE_ORDER,
        min_overlap_modes=MIN_OVERLAP_MODES,
        max_overlap_modes=MAX_OVERLAP_MODES,
        overlap_rel_tol=OVERLAP_REL_TOL,
        lepton_sweep_inputs=lepton_sweep_inputs(),
        include_higgs_yukawas=True,
    )


def diagonal_higgs_yukawas():
    return diagonal_higgs_point().extras["rs_higgs_yukawas"]


def diagonal_constraint_point():
    return point_builder.make_point(rs_higgs_yukawas=diagonal_higgs_yukawas())


@lru_cache(maxsize=None)
def live_lepton_couplings(mkk_gev: float = 3000.0):
    point = diagonal_higgs_point(mkk_gev)
    lepton = point.extras["lepton_mass_basis_couplings"]
    anarchic_y_e = np.array(
        [
            [0.7, 0.2 + 0.1j, -0.3j],
            [0.1 - 0.2j, 1.1, 0.4 + 0.2j],
            [-0.15j, -0.25 + 0.05j, 1.4],
        ],
        dtype=np.complex128,
    )
    rotated = replace(
        lepton,
        c_L=np.array([0.51, 0.65, 0.74], dtype=float),
        c_E=np.array([0.68, 0.57, 0.48], dtype=float),
        U_e_L=_rot12(0.28) @ _rot23(-0.22) @ _rot13(0.18, 0.35),
        U_e_R=_rot12(-0.19) @ _rot23(0.24) @ _rot13(-0.16, -0.45),
        Y_E_bar_matrix=anarchic_y_e,
    )
    return rotated


@lru_cache(maxsize=None)
def live_higgs_yukawas(mkk_gev: float = 3000.0):
    point = diagonal_higgs_point(mkk_gev)
    return build_rs_higgs_yukawas(
        live_lepton_couplings(mkk_gev),
        spectrum=point.extras["rs_ew_spectrum"],
    )


def live_constraint_point(mkk_gev: float = 3000.0):
    return point_builder.make_point(rs_higgs_yukawas=live_higgs_yukawas(mkk_gev))


def constraint_point_with_pair(
    initial: str,
    final: str,
    yukawa_ij: complex,
    yukawa_ji: complex = 0.0j,
):
    source = live_higgs_yukawas()
    matrix = np.array(source.higgs_yukawa_matrix, dtype=np.complex128, copy=True)
    i = _LEPTON_INDEX[initial]
    j = _LEPTON_INDEX[final]
    matrix[i, j] = complex(yukawa_ij)
    matrix[j, i] = complex(yukawa_ji)
    shifted = replace(
        source,
        higgs_yukawa_matrix=matrix,
        diagnostics={
            **dict(source.diagnostics),
            "test_pair_matrix_override": f"{initial}_{final}",
        },
    )
    return point_builder.make_point(rs_higgs_yukawas=shifted)


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
