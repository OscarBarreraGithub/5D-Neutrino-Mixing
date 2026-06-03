"""Shared Phase-4d helpers for active-neutrino rare FCNC tests."""

from __future__ import annotations

from functools import lru_cache

import numpy as np

from flavor_catalog_constraints import point_builder
from quarkConstraints.rs_ew_couplings import DEFAULT_A_REF_C
from tests.rs_ew_phase3b_helpers import (
    MAX_OVERLAP_MODES,
    MIN_OVERLAP_MODES,
    N_GAUGE_MODES,
    OVERLAP_REL_TOL,
    QUADRATURE_ORDER,
    _identity_fit,
    _sample_fit,
    _scales_for_mkk,
)


def lepton_inputs(
    c_l: float = 0.58,
    *,
    alpha: float = 0.0,
    beta: float = 0.0,
) -> dict[str, object]:
    return {
        "c_L": float(c_l),
        "c_E": [float(c_l), float(c_l), float(c_l)],
        "c_N": 0.27,
        "M_N": 1.22e18,
        "lightest_nu_mass": 0.002,
        "ordering": "normal",
        "majorana_alpha": float(alpha),
        "majorana_beta": float(beta),
    }


def build_rs_ew_point(fit, *, alpha: float = 0.0, beta: float = 0.0):
    lambda_ir, k = _scales_for_mkk(3000.0)
    return point_builder.build_from_rs_ew_inputs(
        fit,
        Lambda_IR=lambda_ir,
        k=k,
        n_gauge_modes=N_GAUGE_MODES,
        quadrature_order=QUADRATURE_ORDER,
        min_overlap_modes=MIN_OVERLAP_MODES,
        max_overlap_modes=MAX_OVERLAP_MODES,
        overlap_rel_tol=OVERLAP_REL_TOL,
        lepton_sweep_inputs=lepton_inputs(alpha=alpha, beta=beta),
    )


@lru_cache(maxsize=None)
def rigorous_point(*, alpha: float = 0.0, beta: float = 0.0):
    return build_rs_ew_point(_sample_fit(), alpha=alpha, beta=beta)


@lru_cache(maxsize=None)
def sm_limit_point():
    ref = np.array([DEFAULT_A_REF_C, DEFAULT_A_REF_C, DEFAULT_A_REF_C], dtype=float)
    return build_rs_ew_point(_identity_fit(ref, ref, ref))


def nunu_block(point, transition: str):
    bundle = point.extras["rs_semileptonic_wilsons"]
    if transition == "b_to_s":
        return bundle.b_to_s_nunu
    if transition == "s_to_d":
        return bundle.s_to_d_nunu
    raise ValueError(f"unsupported transition {transition!r}")


def scalar_x_np(block) -> tuple[complex, complex, complex]:
    x_left = complex(np.trace(block.x_np_left) / 3.0)
    x_right = complex(np.trace(block.x_np_right) / 3.0)
    return x_left, x_right, complex(x_left + x_right)


def direct_contact_x_np(point, transition: str) -> tuple[complex, complex, complex]:
    if transition == "b_to_s":
        i, j = 1, 2
    elif transition == "s_to_d":
        i, j = 0, 1
    else:
        raise ValueError(f"unsupported transition {transition!r}")
    couplings = point.extras["rs_ew_couplings"]
    block = nunu_block(point, transition)
    left_contact = sum(
        couplings.nunu_contact("d", "L", i, j, a, a) for a in range(3)
    ) / 3.0
    right_contact = sum(
        couplings.nunu_contact("d", "R", i, j, a, a) for a in range(3)
    ) / 3.0
    x_left = complex(left_contact / block.g_sm_squared_gev_minus2)
    x_right = complex(right_contact / block.g_sm_squared_gev_minus2)
    return x_left, x_right, complex(x_left + x_right)
