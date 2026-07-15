import numpy as np
import pytest

from scripts.anarchic_bauer_s1 import (
    BAUER_REPO_MASS_PREFAC_GEV,
    DEFAULT_V_GEV,
    SCENARIO_SEED_STRIDE,
    _fn_c_values,
    _job_seed,
    _repo_wolfenstein_lam_A,
)
from warpConfig.wavefuncs import f_IR


def _toy_targets():
    return {
        "up_masses_GeV": np.array([1.0e-3, 0.60, 160.0], dtype=float),
        "down_masses_GeV": np.array([2.0e-3, 5.0e-2, 2.8], dtype=float),
        "abs_V_us": 0.2243,
        "abs_V_cb": 0.0408,
        "abs_V_ub": 0.00382,
        "J": 3.0e-5,
    }


def test_bauer_repo_bridge_uses_two_sqrt2_factors_in_fn_mass_inversion():
    targets = _toy_targets()
    epsilon = 1.0e-15
    y_u = np.eye(3, dtype=np.complex128)
    y_d = np.eye(3, dtype=np.complex128)

    c_Q, c_u, c_d = _fn_c_values(
        -0.45,
        epsilon,
        targets,
        Y_u=y_u,
        Y_d=y_d,
        common_cd=False,
    )
    f_Q = f_IR(c_Q, epsilon)
    f_u = f_IR(c_u, epsilon)
    f_d = f_IR(c_d, epsilon)
    d_Q, d_u, d_d = np.diag(f_Q), np.diag(f_u), np.diag(f_d)

    m_u = BAUER_REPO_MASS_PREFAC_GEV * d_Q @ y_u @ d_u
    m_d = BAUER_REPO_MASS_PREFAC_GEV * d_Q @ y_d @ d_d
    old_one_sqrt2_m_u = DEFAULT_V_GEV * d_Q @ y_u @ d_u
    old_one_sqrt2_m_d = DEFAULT_V_GEV * d_Q @ y_d @ d_d

    np.testing.assert_allclose(np.diag(m_u), targets["up_masses_GeV"], rtol=2.0e-9)
    np.testing.assert_allclose(np.diag(m_d), targets["down_masses_GeV"], rtol=2.0e-9)
    np.testing.assert_allclose(
        np.diag(old_one_sqrt2_m_u),
        0.5 * targets["up_masses_GeV"],
        rtol=2.0e-9,
    )
    np.testing.assert_allclose(
        np.diag(old_one_sqrt2_m_d),
        0.5 * targets["down_masses_GeV"],
        rtol=2.0e-9,
    )


def test_fn_left_doublet_hierarchy_includes_repo_wolfenstein_A():
    targets = _toy_targets()
    lam, A = _repo_wolfenstein_lam_A(targets)
    c_Q, _, _ = _fn_c_values(
        -0.45,
        1.0e-15,
        targets,
        Y_u=np.eye(3, dtype=np.complex128),
        Y_d=np.eye(3, dtype=np.complex128),
        common_cd=False,
    )
    f_Q = f_IR(c_Q, 1.0e-15)

    assert A == pytest.approx(targets["abs_V_cb"] / targets["abs_V_us"] ** 2)
    assert f_Q[0] / f_Q[2] == pytest.approx(A * lam**3, rel=2.0e-9)
    assert f_Q[1] / f_Q[2] == pytest.approx(A * lam**2, rel=2.0e-9)


def test_bauer_tile_seed_is_deterministic_and_scenario_indexed():
    base = 20260623

    assert _job_seed(base, 2, "S1") == base + 2 * 1009
    assert _job_seed(base, 2, "S4") == base + 2 * 1009 + 3 * SCENARIO_SEED_STRIDE
    assert _job_seed(base, 2, "S4") == _job_seed(base, 2, "S4")
    with pytest.raises(ValueError, match="unknown Bauer scenario"):
        _job_seed(base, 0, "not-a-scenario")
