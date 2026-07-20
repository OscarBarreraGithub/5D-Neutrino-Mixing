"""Tests for CKM phase extraction helpers."""

from __future__ import annotations

import cmath
import math

import numpy as np
import pytest

from quarkConstraints.ckm_extraction import (
    ckm_phases_from_matrix,
    neutral_b_mixing_sm_amplitude,
    repo_default_ckm_matrix,
    repo_default_ckm_phases,
)


def test_neutral_b_sm_amplitude_restores_complex_top_box_phase():
    matrix = repo_default_ckm_matrix()

    for light_index in (0, 1):
        result = neutral_b_mixing_sm_amplitude(
            delta_m_sm_gev=2.5e-13,
            light_down_index=light_index,
            ckm=matrix,
            ckm_source="independent-test-matrix",
        )
        expected_ckm_factor = (
            np.conjugate(matrix[2, light_index]) * matrix[2, 2]
        ) ** 2
        expected = 1.25e-13 * expected_ckm_factor / abs(expected_ckm_factor)

        assert result.magnitude_gev == pytest.approx(1.25e-13)
        assert result.ckm_factor == pytest.approx(expected_ckm_factor)
        assert result.amplitude_gev == pytest.approx(expected)
        assert result.ckm_source == "independent-test-matrix"


def test_neutral_b_sm_amplitude_and_qb_squared_rephase_together():
    matrix = repo_default_ckm_matrix()
    rng = np.random.default_rng(7262026)
    qb_coupling = 0.003 - 0.007j

    for light_index in (0, 1):
        base = neutral_b_mixing_sm_amplitude(
            delta_m_sm_gev=2.5e-13,
            light_down_index=light_index,
            ckm=matrix,
        )
        base_ratio = qb_coupling**2 / base.amplitude_gev
        for _ in range(16):
            up_phases = np.exp(1j * rng.uniform(-math.pi, math.pi, size=3))
            down_phases = np.exp(1j * rng.uniform(-math.pi, math.pi, size=3))
            rephased_ckm = (
                np.diag(np.conjugate(up_phases))
                @ matrix
                @ np.diag(down_phases)
            )
            rephased_qb = (
                np.conjugate(down_phases[light_index])
                * qb_coupling
                * down_phases[2]
            )
            rephased = neutral_b_mixing_sm_amplitude(
                delta_m_sm_gev=2.5e-13,
                light_down_index=light_index,
                ckm=rephased_ckm,
            )

            assert rephased_qb**2 / rephased.amplitude_gev == pytest.approx(
                base_ratio,
                rel=1e-12,
                abs=1e-12,
            )


def test_repo_default_ckm_phases_match_pdg_windows_and_committed_values():
    phases = repo_default_ckm_phases()

    assert phases.two_beta_degrees == pytest.approx(45.5, abs=0.6)
    assert phases.sin_2beta == pytest.approx(0.70, abs=0.02)
    assert phases.beta_s == pytest.approx(0.019, abs=0.004)
    assert phases.phi_s == pytest.approx(-0.0368, abs=0.003)

    assert phases.beta_degrees == pytest.approx(22.549961677305074)
    assert phases.two_beta_degrees == pytest.approx(45.09992335461015)
    assert phases.sin_2beta == pytest.approx(0.7083388934693238)
    assert phases.phi_s == pytest.approx(-0.037945102831894784)


def test_ckm_phase_helper_matches_independent_matrix_ratio_recompute():
    matrix = repo_default_ckm_matrix()
    phases = ckm_phases_from_matrix(matrix)

    beta = cmath.phase(
        -(matrix[1, 0] * np.conjugate(matrix[1, 2]))
        / (matrix[2, 0] * np.conjugate(matrix[2, 2]))
    )
    beta_s = cmath.phase(
        -(matrix[2, 1] * np.conjugate(matrix[2, 2]))
        / (matrix[1, 1] * np.conjugate(matrix[1, 2]))
    )

    assert phases.beta == pytest.approx(beta)
    assert phases.two_beta == pytest.approx(2.0 * beta)
    assert phases.sin_2beta == pytest.approx(math.sin(2.0 * beta))
    assert phases.beta_s == pytest.approx(beta_s)
    assert phases.phi_s == pytest.approx(-2.0 * beta_s)


def test_ckm_phase_helper_rejects_bad_or_singular_matrices():
    with pytest.raises(ValueError, match="shape"):
        ckm_phases_from_matrix(np.eye(2))

    with pytest.raises(ValueError, match="non-zero"):
        ckm_phases_from_matrix(np.eye(3))
