"""Tests for the quark-sector exact fit layer."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from quarkConstraints.benchmarks import default_quark_targets, default_spurion_seed
from quarkConstraints.fit import QuarkTargets, fit_quark_sector, fit_residuals


def test_fit_quark_sector_improves_seed_score_and_matches_targets():
    """The fit should improve sharply over the seed and reproduce the targets."""
    targets = default_quark_targets()
    solution = fit_quark_sector(
        targets,
        seed=default_spurion_seed(),
        overall_scale=3.0,
        max_nfev=120,
    )

    assert solution.success
    assert solution.result.score < solution.initial_score
    assert np.max(np.abs(solution.result.mass_residuals_up)) < 1.0e-2
    assert np.max(np.abs(solution.result.mass_residuals_down)) < 1.0e-2
    assert np.max(np.abs(solution.result.ckm_residuals)) < 1.0e-2


def test_fitted_ckm_matrix_is_unitary():
    """The extracted CKM matrix should remain unitary after exact diagonalization."""
    solution = fit_quark_sector(default_quark_targets(), overall_scale=3.0, max_nfev=120)
    ckm = solution.result.ckm_matrix
    assert np.allclose(ckm.conjugate().T @ ckm, np.eye(3), atol=1e-10)


def test_ckm_observable_residuals_are_invariant_under_jarlskog_sign_convention():
    solution = fit_quark_sector(default_quark_targets(), overall_scale=3.0, max_nfev=120)
    ckm = solution.result.ckm_matrix
    targets = default_quark_targets()
    flipped_targets = QuarkTargets(
        up_masses=targets.up_masses,
        down_masses=targets.down_masses,
        ckm=np.conjugate(targets.ckm),
        label="conjugated-targets",
    )

    residuals = fit_residuals(
        solution.result.masses_up,
        solution.result.masses_down,
        ckm,
        targets,
    )
    flipped = fit_residuals(
        solution.result.masses_up,
        solution.result.masses_down,
        np.conjugate(ckm),
        flipped_targets,
    )

    assert np.allclose(
        residuals["ckm_observable_residuals"][:3],
        flipped["ckm_observable_residuals"][:3],
        atol=1e-12,
    )
    assert np.isclose(
        residuals["ckm_observable_residuals"][3],
        -flipped["ckm_observable_residuals"][3],
        atol=1e-12,
    )
    assert np.isclose(residuals["total_score"], flipped["total_score"], atol=1e-12)
