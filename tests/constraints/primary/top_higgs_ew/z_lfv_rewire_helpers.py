"""Shared helpers for T015-T017 Phase-4b Z-LFV rewire tests."""

from __future__ import annotations

import math
from dataclasses import replace
from functools import lru_cache

import numpy as np

from flavor_catalog_constraints import point_builder
from quarkConstraints.rs_ew_couplings import build_rs_ew_couplings
from quarkConstraints.zpole import default_sm_inputs, partial_width_weight, sm_couplings
from tests.rs_ew_phase3b_helpers import (
    MAX_OVERLAP_MODES,
    MIN_OVERLAP_MODES,
    N_GAUGE_MODES,
    OVERLAP_REL_TOL,
    QUADRATURE_ORDER,
    _sample_fit,
    _scales_for_mkk,
)

LFV_CHANNELS = {
    "T015": {
        "indices": (0, 1),
        "initial_flavor": "e",
        "final_flavor": "mu",
        "left_key": "delta_g_left_emu",
        "right_key": "delta_g_right_emu",
        "observable": "BR(Z -> e mu)",
    },
    "T016": {
        "indices": (0, 2),
        "initial_flavor": "e",
        "final_flavor": "tau",
        "left_key": "delta_g_left_etau",
        "right_key": "delta_g_right_etau",
        "observable": "BR(Z -> e tau)",
    },
    "T017": {
        "indices": (1, 2),
        "initial_flavor": "mu",
        "final_flavor": "tau",
        "left_key": "delta_g_left_mutau",
        "right_key": "delta_g_right_mutau",
        "observable": "BR(Z -> mu tau)",
    },
}


def lepton_sweep_inputs() -> dict[str, object]:
    return {
        "c_L": 0.58,
        "c_E": [0.61, 0.59, 0.52],
        "c_N": 0.27,
        "M_N": 1.22e18,
        "lightest_nu_mass": 0.002,
        "ordering": "normal",
        "majorana_alpha": 0.0,
        "majorana_beta": 0.0,
    }


@lru_cache(maxsize=None)
def diagonal_rs_ew_point(mkk_gev: float = 3000.0):
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
    )


@lru_cache(maxsize=None)
def lfv_live_rs_ew_point(mkk_gev: float = 3000.0):
    base_point = diagonal_rs_ew_point(mkk_gev)
    lepton = replace(
        base_point.extras["lepton_mass_basis_couplings"],
        c_L=np.array([0.51, 0.65, 0.74], dtype=float),
        c_E=np.array([0.68, 0.57, 0.48], dtype=float),
        U_e_L=_rot12(0.28) @ _rot23(-0.22) @ _rot13(0.18, 0.35),
        U_e_R=_rot12(-0.19) @ _rot23(0.24) @ _rot13(-0.16, -0.45),
    )
    couplings = build_rs_ew_couplings(
        _sample_fit(),
        spectrum=base_point.extras["rs_ew_spectrum"],
        lepton_mass_basis_couplings=lepton,
        overlap_rel_tol=OVERLAP_REL_TOL,
        min_overlap_modes=MIN_OVERLAP_MODES,
        max_overlap_modes=MAX_OVERLAP_MODES,
    )
    return point_builder.make_point(rs_ew_couplings=couplings)


def old_style_lepton_only_point():
    return point_builder.make_point(
        lepton_mass_basis_couplings=diagonal_rs_ew_point().extras[
            "lepton_mass_basis_couplings"
        ]
    )


def amplified_tree_point(process_id: str, delta_g: complex = 0.01):
    couplings = lfv_live_rs_ew_point().extras["rs_ew_couplings"]
    row, column = LFV_CHANNELS[process_id]["indices"]
    z_left = np.array(couplings.z_delta_g_L_e, dtype=np.complex128, copy=True)
    z_right = np.array(couplings.z_delta_g_R_e, dtype=np.complex128, copy=True)
    z_left[row, column] = complex(delta_g)
    z_left[column, row] = complex(delta_g).conjugate()
    z_right[row, column] = 0.0j
    z_right[column, row] = 0.0j
    amplified = replace(
        couplings,
        z_delta_g_L_e=z_left,
        z_delta_g_R_e=z_right,
    )
    return point_builder.make_point(rs_ew_couplings=amplified)


def manual_sm_width_weights() -> dict[str, float]:
    inputs = default_sm_inputs()
    weights = {}
    for flavor in ("u", "d", "s", "c", "b"):
        weights[flavor] = partial_width_weight(
            sm_couplings(flavor, inputs),
            radiator=inputs.radiator_for(flavor),
        )
    for flavor in ("e", "mu", "tau"):
        weights[flavor] = partial_width_weight(sm_couplings(flavor, inputs))
    for flavor in ("nu_e", "nu_mu", "nu_tau"):
        weights[flavor] = partial_width_weight(sm_couplings("nu", inputs))
    return weights


def manual_lfv_br(delta_g_left: complex, delta_g_right: complex) -> tuple[float, float, float]:
    norm = float(abs(delta_g_left) ** 2 + abs(delta_g_right) ** 2)
    sm_total = float(sum(manual_sm_width_weights().values()))
    lfv_weight = 2.0 * norm
    return float(lfv_weight / (sm_total + lfv_weight)), norm, sm_total


def channel_deltas(couplings, process_id: str) -> tuple[complex, complex]:
    row, column = LFV_CHANNELS[process_id]["indices"]
    return (
        complex(couplings.z_delta_g_L_e[row, column]),
        complex(couplings.z_delta_g_R_e[row, column]),
    )


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
