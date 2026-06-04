"""Tests for CKM phase extraction helpers."""

from __future__ import annotations

import cmath
import math

import numpy as np
import pytest

from quarkConstraints.ckm_extraction import (
    ckm_phases_from_matrix,
    repo_default_ckm_matrix,
    repo_default_ckm_phases,
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
