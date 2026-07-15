"""Regression tests for the Phase 2 Wilson-RG audit."""

from __future__ import annotations

import numpy as np
import pytest

from qcd.constants import ALPHA_S_MZ, M_TOP_MS
from qcd.running import alpha_s as reference_alpha_s
from quarkConstraints.qcd_running import (
    _nf_for_scale,
    evolve_deltaf2_wilsons,
    run_alpha_s,
)


MU_HIGH_GEV = 3000.0
MU_LOW_GEV = 2.0
M_T_GEV = M_TOP_MS
M_B_GEV = 4.18
M_C_GEV = 1.27


def _as_real_tuple(values: tuple[complex, complex, complex, complex]) -> tuple[float, ...]:
    assert all(abs(value.imag) < 1.0e-14 for value in values)
    return tuple(float(value.real) for value in values)


def test_unit_wilson_evolution_matches_audited_lo_reference() -> None:
    expected = {
        "C1_VLL": (0.7288737958332887, 0.0, 0.0, 0.0),
        "C1_VRR": (0.0, 0.7288737958332887, 0.0, 0.0),
        "C4_LR": (0.0, 0.0, 3.5431590851191723, 0.0),
        "C5_LR": (0.0, 0.0, 0.8964726753709399, 0.8537410590063527),
    }
    unit_vectors = {
        "C1_VLL": (1.0, 0.0, 0.0, 0.0),
        "C1_VRR": (0.0, 1.0, 0.0, 0.0),
        "C4_LR": (0.0, 0.0, 1.0, 0.0),
        "C5_LR": (0.0, 0.0, 0.0, 1.0),
    }

    for name, vector in unit_vectors.items():
        observed = _as_real_tuple(
            evolve_deltaf2_wilsons(*vector, MU_HIGH_GEV, MU_LOW_GEV)
        )
        assert observed == pytest.approx(expected[name], rel=2.0e-12, abs=1.0e-14)


def test_alpha_s_uses_top_bottom_thresholds_for_3tev_to_2gev_path() -> None:
    assert _nf_for_scale(3000.0, m_t=M_T_GEV, m_b=M_B_GEV, m_c=M_C_GEV) == 6
    assert _nf_for_scale(M_T_GEV, m_t=M_T_GEV, m_b=M_B_GEV, m_c=M_C_GEV) == 5
    assert _nf_for_scale(M_B_GEV, m_t=M_T_GEV, m_b=M_B_GEV, m_c=M_C_GEV) == 4
    assert _nf_for_scale(MU_LOW_GEV, m_t=M_T_GEV, m_b=M_B_GEV, m_c=M_C_GEV) == 4

    for mu in (MU_HIGH_GEV, M_T_GEV, M_B_GEV, MU_LOW_GEV):
        observed = run_alpha_s(mu)
        expected = reference_alpha_s(
            mu,
            n_loops=1,
            matching_loops=0,
            alpha_s_ref=ALPHA_S_MZ,
        )
        assert observed == pytest.approx(expected, rel=1.0e-10, abs=1.0e-14)

    no_top_alpha_3tev = 0.0784656638154406
    assert not np.isclose(run_alpha_s(MU_HIGH_GEV), no_top_alpha_3tev)
