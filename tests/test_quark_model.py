"""Tests for the quark-sector MFV model layer."""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from quarkConstraints.benchmarks import benchmark_spurion_input
from quarkConstraints.model import QuarkSpurionPoint, derive_bulk_state


def test_bulk_state_matrices_are_hermitian_and_overlaps_are_finite():
    """The derived spurion invariants should be Hermitian with finite overlaps."""
    state = derive_bulk_state(benchmark_spurion_input())

    assert np.allclose(state.C_Q, state.C_Q.conjugate().T)
    assert np.allclose(state.C_u, state.C_u.conjugate().T)
    assert np.allclose(state.C_d, state.C_d.conjugate().T)
    assert np.all(state.eig_Q >= -1e-12)
    assert np.all(state.eig_u >= -1e-12)
    assert np.all(state.eig_d >= -1e-12)
    assert np.all(np.isfinite(state.F_Q))
    assert np.all(np.isfinite(state.F_u))
    assert np.all(np.isfinite(state.F_d))
    assert np.all(state.F_Q > 0.0)
    assert np.all(state.F_u > 0.0)
    assert np.all(state.F_d > 0.0)


def test_bulk_basis_rotations_diagonalize_spurion_invariants():
    """The stored rotations should diagonalize the three C matrices."""
    state = derive_bulk_state(benchmark_spurion_input())

    diag_q = state.rotation_Q.conjugate().T @ state.C_Q @ state.rotation_Q
    diag_u = state.rotation_u.conjugate().T @ state.C_u @ state.rotation_u
    diag_d = state.rotation_d.conjugate().T @ state.C_d @ state.rotation_d

    assert np.allclose(diag_q, np.diag(state.eig_Q), atol=1e-10)
    assert np.allclose(diag_u, np.diag(state.eig_u), atol=1e-10)
    assert np.allclose(diag_d, np.diag(state.eig_d), atol=1e-10)


def test_negative_r_is_rejected_at_model_construction():
    point = benchmark_spurion_input()

    with pytest.raises(ValueError, match="non-negative"):
        QuarkSpurionPoint(Y_u=point.Y_u, Y_d=point.Y_d, r=-0.1)
