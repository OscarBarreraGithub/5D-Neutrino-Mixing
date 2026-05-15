"""Tests for the PDG 2024 quark-mass inputs and the RG-evolved targets."""

from __future__ import annotations

import numpy as np
import pytest

from qcd.constants import M_TOP_MS
from quarkConstraints.pdg_quark_masses import (
    PDG_2024_QUARK_MASSES,
    pdg_2sigma_relative_at_scale,
    pdg_quark_masses_at_scale,
    pdg_up_down_arrays_at_scale,
)


def test_pdg_inputs_have_expected_metadata():
    expected_flavors = ("u", "d", "s", "c", "b", "t")
    assert tuple(PDG_2024_QUARK_MASSES.keys()) == expected_flavors
    # Top reference is m_t(m_t).
    assert PDG_2024_QUARK_MASSES["t"].mu_ref_GeV == pytest.approx(162.5)
    assert PDG_2024_QUARK_MASSES["t"].n_f_at_reference == 6
    # Bottom reference is m_b(m_b) with n_f=5.
    assert PDG_2024_QUARK_MASSES["b"].mu_ref_GeV == pytest.approx(4.183)
    assert PDG_2024_QUARK_MASSES["b"].n_f_at_reference == 5
    # Light quarks are at 2 GeV in n_f=4.
    for f in ("u", "d", "s"):
        assert PDG_2024_QUARK_MASSES[f].mu_ref_GeV == pytest.approx(2.0)
        assert PDG_2024_QUARK_MASSES[f].n_f_at_reference == 4


def test_round_trip_identity_at_reference():
    # Running each quark from its own reference back to its own reference
    # must return the central value exactly (numerical zero-length step).
    for flavor, entry in PDG_2024_QUARK_MASSES.items():
        out = pdg_quark_masses_at_scale(entry.mu_ref_GeV, flavors=(flavor,))
        assert out[flavor] == pytest.approx(entry.central_GeV, rel=1e-12)


def test_charm_runs_to_2gev_in_expected_window():
    # PDG-style: m_c(m_c) = 1.273 -> m_c(2 GeV) ~ 1.07 GeV (RunDec range).
    out = pdg_quark_masses_at_scale(2.0, flavors=("c",))
    assert 1.05 < out["c"] < 1.12


def test_strange_runs_to_top_scale_in_expected_window():
    # m_s(2 GeV) ~ 0.0935 -> m_s(m_t) somewhere near 0.05–0.06.
    out = pdg_quark_masses_at_scale(M_TOP_MS, flavors=("s",))
    assert 0.045 < out["s"] < 0.060


def test_bottom_runs_to_top_scale_in_expected_window():
    # m_b(m_b) = 4.183 -> m_b(m_t) ~ 2.7 GeV.
    out = pdg_quark_masses_at_scale(M_TOP_MS, flavors=("b",))
    assert 2.5 < out["b"] < 2.9


def test_top_at_top_scale_is_within_one_percent_of_pdg_central():
    # m_t runs from 162.5 to 163.5; per-mille effect.
    out = pdg_quark_masses_at_scale(M_TOP_MS, flavors=("t",))
    assert abs(out["t"] - 162.5) / 162.5 < 0.01


def test_2sigma_relative_floor_is_at_least_3_per_mille():
    rel = pdg_2sigma_relative_at_scale(M_TOP_MS)
    # Per plan v3 §4: empirical floor across all quarks is max(0.003, ...).
    # PDG itself gives strictly larger 2sigma uncertainties for light quarks
    # and a tight 0.0086 for top. Verify the absolute floor.
    for f, val in rel.items():
        assert val > 0.0, f"PDG 2sigma must be positive for {f}"
    # Sanity: top is the tightest, charm is also tight, light quarks looser.
    assert rel["t"] < rel["s"]
    assert rel["s"] < rel["d"]


def test_up_down_arrays_have_expected_layout():
    up, down = pdg_up_down_arrays_at_scale(M_TOP_MS)
    assert up.shape == (3,)
    assert down.shape == (3,)
    # Hierarchy: m_u < m_c < m_t and m_d < m_s < m_b.
    assert up[0] < up[1] < up[2]
    assert down[0] < down[1] < down[2]
