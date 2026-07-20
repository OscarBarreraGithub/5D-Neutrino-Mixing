"""Shared Phase-5b charged-current test points."""

from __future__ import annotations

from dataclasses import replace
from functools import lru_cache

import numpy as np

from flavor_catalog_constraints import point_builder
from quarkConstraints.rs_charged_current import RSChargedCurrentCouplings
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


def lepton_inputs(c_l: float = 0.58) -> dict[str, object]:
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


def build_charged_point(fit, *, c_l: float = 0.58, mkk_gev: float = 3000.0):
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
        lepton_sweep_inputs=lepton_inputs(c_l),
        include_charged_current=True,
    )


@lru_cache(maxsize=None)
def sample_charged_point():
    return build_charged_point(_sample_fit(), c_l=0.58)


@lru_cache(maxsize=None)
def universal_charged_point():
    ref = np.array([DEFAULT_A_REF_C, DEFAULT_A_REF_C, DEFAULT_A_REF_C], dtype=float)
    return build_charged_point(_identity_fit(ref, ref, ref), c_l=DEFAULT_A_REF_C)


def charged_with_epsilon(
    base: RSChargedCurrentCouplings,
    updates: dict[tuple[int, int, int], complex],
) -> RSChargedCurrentCouplings:
    epsilon = np.array(base.epsilon, dtype=np.complex128, copy=True)
    for index, value in updates.items():
        epsilon[index] = complex(value)
    return replace(
        base,
        epsilon=epsilon,
        delta_abs_vij_over_vij=np.real(epsilon),
    )
