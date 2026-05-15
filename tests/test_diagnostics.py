"""Tests for ``quarkConstraints.diagnostics``.

The synthetic round-trip uses a Yukawa configuration that *exactly*
reproduces a chosen MS-bar mass spectrum at the fitter scale and asserts
that ``extract_msbar_masses_from_yukawa_row`` returns those masses (to
running-numerical precision) when called with that scale.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from qcd.constants import M_TOP_MS
from quarkConstraints.diagnostics import extract_msbar_masses_from_yukawa_row
from warpConfig.baseParams import V_EWSB, get_warp_params
from warpConfig.wavefuncs import f_IR


def _build_synthetic_row(
    *,
    up_masses: tuple[float, float, float],
    down_masses: tuple[float, float, float],
    c_Q: tuple[float, float, float],
    c_u: tuple[float, float, float],
    c_d: tuple[float, float, float],
    Lambda_IR: float = 3000.0,
    k: float = 1.2209e19,
) -> pd.Series:
    """Build a CSV-shaped row whose Yukawas reproduce (up_masses, down_masses).

    We pick diagonal Y_u, Y_d so that
        m_q[i] = 2 v F_Q[i] Y[i,i] F_q[i].
    """
    params = get_warp_params(k=k, Lambda_IR=Lambda_IR)
    eps = float(params["epsilon"])
    F_Q = np.asarray(f_IR(np.array(c_Q), eps), dtype=float)
    F_u = np.asarray(f_IR(np.array(c_u), eps), dtype=float)
    F_d = np.asarray(f_IR(np.array(c_d), eps), dtype=float)

    Y_u = np.zeros((3, 3), dtype=np.complex128)
    Y_d = np.zeros((3, 3), dtype=np.complex128)
    for i in range(3):
        Y_u[i, i] = up_masses[i] / (2.0 * V_EWSB * F_Q[i] * F_u[i])
        Y_d[i, i] = down_masses[i] / (2.0 * V_EWSB * F_Q[i] * F_d[i])

    data: dict[str, float] = {}
    for i in range(3):
        for j in range(3):
            data[f"Y_u_{i+1}{j+1}_re"] = float(np.real(Y_u[i, j]))
            data[f"Y_u_{i+1}{j+1}_im"] = float(np.imag(Y_u[i, j]))
            data[f"Y_d_{i+1}{j+1}_re"] = float(np.real(Y_d[i, j]))
            data[f"Y_d_{i+1}{j+1}_im"] = float(np.imag(Y_d[i, j]))
    for i, val in enumerate(c_Q):
        data[f"c_Q{i+1}"] = float(val)
    for i, val in enumerate(c_u):
        data[f"c_u{i+1}"] = float(val)
    for i, val in enumerate(c_d):
        data[f"c_d{i+1}"] = float(val)
    return pd.Series(data)


def test_synthetic_row_round_trip_at_fitter_scale():
    # Pick masses already at m_t(m_t); request the same scale; expect equality.
    up = (1.0e-3, 0.6, 162.5)
    down = (3.0e-3, 0.05, 2.7)
    row = _build_synthetic_row(
        up_masses=up,
        down_masses=down,
        c_Q=(0.55, 0.50, 0.45),
        c_u=(0.65, 0.55, 0.40),
        c_d=(0.66, 0.58, 0.55),
    )
    out = extract_msbar_masses_from_yukawa_row(row, scale_GeV=M_TOP_MS)
    expected = {
        "u": up[0], "c": up[1], "t": up[2],
        "d": down[0], "s": down[1], "b": down[2],
    }
    for f, val in expected.items():
        assert out[f] == pytest.approx(val, rel=1e-6), f"{f}: got {out[f]} expected {val}"


def test_diagnostic_runs_and_returns_six_positive_finite_values():
    row = _build_synthetic_row(
        up_masses=(1e-3, 0.5, 160.0),
        down_masses=(2.5e-3, 0.05, 2.7),
        c_Q=(0.55, 0.50, 0.45),
        c_u=(0.65, 0.55, 0.40),
        c_d=(0.66, 0.58, 0.55),
    )
    out = extract_msbar_masses_from_yukawa_row(row, scale_GeV=M_TOP_MS)
    assert set(out.keys()) == {"u", "c", "t", "d", "s", "b"}
    for f, m in out.items():
        assert np.isfinite(m), f"{f} produced non-finite mass {m}"
        assert m > 0.0, f"{f} produced non-positive mass {m}"


def test_diagnostic_can_run_to_a_lower_scale():
    # Run output to mu = 2 GeV and confirm light masses move only modestly
    # relative to m_t scale (qualitative direction: heavier at lower mu).
    row = _build_synthetic_row(
        up_masses=(1e-3, 0.5, 160.0),
        down_masses=(2.5e-3, 0.05, 2.7),
        c_Q=(0.55, 0.50, 0.45),
        c_u=(0.65, 0.55, 0.40),
        c_d=(0.66, 0.58, 0.55),
    )
    out_mt = extract_msbar_masses_from_yukawa_row(row, scale_GeV=M_TOP_MS)
    out_2 = extract_msbar_masses_from_yukawa_row(row, scale_GeV=2.0)
    # Light quark masses are larger at lower mu (multiplicative factor > 1).
    assert out_2["s"] > out_mt["s"]
    assert out_2["d"] > out_mt["d"]
