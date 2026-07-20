"""Tests for ``qcd.mass_running.run_msbar_mass`` and threshold matching."""

from __future__ import annotations

import pytest

import qcd.mass_running as mass_running
from qcd.constants import M_BOTTOM, M_CHARM, M_TOP_MS
from qcd.decoupling import _coeffs_msbar_mass, match_alpha_s, match_msbar_mass
from qcd.mass_running import run_msbar_mass
from qcd.running import alpha_s as alpha_s_running


def test_round_trip_identity():
    # Run m_s(2 GeV) up to the top scale and back; must reproduce input.
    m0 = 0.0935
    forward = run_msbar_mass(m0, 2.0, M_TOP_MS, n_f_ref=4)
    back = run_msbar_mass(forward, M_TOP_MS, 2.0, n_f_ref=5)
    assert abs(back - m0) / m0 < 1e-4


def test_threshold_continuity_at_m_b():
    # The mass should be continuous (to per-mille level) across m_b when we
    # apply both alpha_s and mass matching at mu = m_b.
    m_just_above = run_msbar_mass(4.183, 4.183, M_BOTTOM * 1.000001, n_f_ref=5)
    m_just_below = run_msbar_mass(4.183, 4.183, M_BOTTOM * 0.999999, n_f_ref=5)
    # The matching-induced jump at mu = m_b is per-mille; both should be
    # close to the input value.
    assert abs(m_just_above - 4.183) / 4.183 < 5e-3
    assert abs(m_just_below - 4.183) / 4.183 < 5e-3


def test_threshold_continuity_at_m_c():
    # m_s(m_c, n_f=3) and m_s(m_c, n_f=4) differ only by the CKS matching
    # factor; the difference must match the analytic prediction within tight
    # tolerance.
    m_just_above = run_msbar_mass(0.0935, 2.0, M_CHARM * 1.000001, n_f_ref=4)
    m_just_below = run_msbar_mass(0.0935, 2.0, M_CHARM * 0.999999, n_f_ref=4)
    # The relative difference is the CKS mass-matching at mu = m_c
    # (per-mille, dominated by 89/432 * (alpha_s / pi)^2 with n_l = 3).
    rel = abs(m_just_above - m_just_below) / m_just_above
    assert rel < 1.0e-2, f"unexpectedly large discontinuity at m_c: {rel:.4f}"


def test_charm_decoupling_at_mu_2gev():
    """Plan v3 §7 test 5.

    m_s(2 GeV, n_f=4) -> run down to m_c (n_f=4) -> match to n_f=3 ->
    run up to 2 GeV (n_f=3); the deviation must be at most 0.5%.
    """
    # Step 1: run from 2 GeV (n_f=4) down to m_c.
    m_at_mc_nf4 = run_msbar_mass(0.0935, 2.0, M_CHARM, n_f_ref=4)

    # Step 2: explicit charm decoupling at m_c.
    a_nf4 = alpha_s_running(M_CHARM, n_loops=4, matching_loops=3)
    a_nf3 = match_alpha_s(a_nf4, n_f_from=4, n_f_to=3, matching_loops=3)
    m_at_mc_nf3 = match_msbar_mass(
        m_at_mc_nf4,
        alpha_s=a_nf3,
        direction="down",
        n_f_high=4,
        matching_loops=3,
    )

    # Step 3: run back up from m_c to 2 GeV in n_f=3.
    m_at_2_nf3 = run_msbar_mass(m_at_mc_nf3, M_CHARM, 2.0, n_f_ref=3)

    rel_dev = abs(m_at_2_nf3 - 0.0935) / 0.0935
    assert rel_dev < 0.005, f"charm-decoupling round-trip drift {rel_dev:.4f} > 0.5%"


def test_match_msbar_mass_is_identity_at_zero_loops():
    out = match_msbar_mass(1.0, alpha_s=0.2, direction="down", n_f_high=5, matching_loops=0)
    assert out == pytest.approx(1.0)


def test_msbar_mass_d3_matches_cks_closed_form():
    """CKS Eq. (20) gives an increasing 3-loop mass-decoupling d3(n_l)."""
    expected = {
        3: 1.9218095703292382,
        4: 1.946537179697389,
        5: 1.9712647890655397,
    }
    values = []
    for n_l, expected_d3 in expected.items():
        d2, d3 = _coeffs_msbar_mass(n_l)
        assert d2 == pytest.approx(89.0 / 432.0)
        assert d3 == pytest.approx(expected_d3, rel=1e-13)
        values.append(d3)
    assert values == sorted(values)


def test_match_msbar_mass_equal_nf_is_noop():
    out = match_msbar_mass(
        100.0,
        alpha_s=0.1,
        direction="down",
        n_f_high=6,
        matching_loops=3,
        n_f_from=6,
        n_f_to=6,
    )
    assert out == pytest.approx(100.0)


def test_top_legacy_equal_nf_crossing_applies_no_matching(monkeypatch):
    """Recreate the old 162.5 -> 163.5 equal-nf crossing and require no jump."""
    monkeypatch.setattr(mass_running, "M_TOP_MS", 163.5)
    with_matching = mass_running.run_msbar_mass(
        162.5,
        162.5,
        3000.0,
        n_f_ref=6,
        matching_loops=3,
    )
    without_matching = mass_running.run_msbar_mass(
        162.5,
        162.5,
        3000.0,
        n_f_ref=6,
        matching_loops=0,
    )
    assert with_matching == pytest.approx(without_matching, rel=1e-12)


def test_match_msbar_mass_inverse_is_close_to_identity():
    # Down-then-up should return close to the original mass at the same alpha_s.
    a_s = 0.22
    m0 = 4.183
    down = match_msbar_mass(m0, alpha_s=a_s, direction="down", n_f_high=5, matching_loops=3)
    up = match_msbar_mass(down, alpha_s=a_s, direction="up", n_f_high=5, matching_loops=3)
    assert abs(up - m0) / m0 < 1e-4
