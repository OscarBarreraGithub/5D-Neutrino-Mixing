"""Tests for the quark-sector exact fit layer.

R02-I2 (C04) — coverage scope note for the spurion-seed assertions
------------------------------------------------------------------

The C04 dispatch prompt asked whether these tests "only cover k=0" in
the same sense as the Wilson upper-limit test (R07-I2).  They do not.
There is no "k" parameter here: the spurion seed is a deterministic
benchmark point (``default_spurion_seed()`` in
``quarkConstraints.benchmarks``) used by both ``fit_quark_sector`` and
the quotient-equivalent regression family below.  The non-default
seeds exercised in this file include:

  * ``default_spurion_seed()`` itself (test_fit_quark_sector_*),
  * ``_quotient_equivalent_seed(default_spurion_seed())`` — a
    quotient-direction shift covering 2pi-translations of left/right
    rotation angles and a multiplicative rescale of the overall scale
    (test_fit_quark_sector_is_invariant_under_quotient_directions,
    test_fit_orientation_false_remains_deterministic_and_restricted),
  * ad-hoc seeds with non-trivial right-rotation angles
    (test_canonical_vector_quotients_out_*).

That spans the canonical fundamental domain (the "k=0" representative)
*and* multiple quotient-equivalent orbits, so the suite is not
analogous to the k=0-only Wilson-UL gap that R07-I2 closed.  The
deeper provenance for the seed *literals* themselves (where the
numerical values came from) is a separate issue routed to C13.
"""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from quarkConstraints.benchmarks import default_quark_targets, default_spurion_seed
from quarkConstraints.fit import (
    QUARK_FIT_CANONICAL_VECTOR_FIELDS,
    QUARK_FIT_CANONICAL_VECTOR_SIZE,
    QuarkFitSeed,
    QuarkTargets,
    _ordered_dirac_svd,
    _rephase_to_pdg_convention,
    ckm_observables,
    evaluate_quark_fit,
    decode_quark_fit_canonical_vector,
    encode_quark_fit_canonical_vector,
    fit_quark_sector,
    fit_residuals,
    jarlskog_invariant,
    mass_matrix_observables,
)
from quarkConstraints.model import RotationParameters, build_mfv_point_from_singular_values


def _fit_seed_from_benchmark() -> QuarkFitSeed:
    seed = default_spurion_seed()
    return QuarkFitSeed(
        up_singular_values=seed.up_singular_values,
        down_singular_values=seed.down_singular_values,
        overall_scale=seed.overall_scale,
        up_left=seed.up_left,
        up_right=seed.up_right,
        down_left=seed.down_left,
        down_right=seed.down_right,
    )


def _assert_rotation_close(actual: RotationParameters, expected: RotationParameters) -> None:
    assert np.isclose(actual.theta12, expected.theta12)
    assert np.isclose(actual.theta13, expected.theta13)
    assert np.isclose(actual.theta23, expected.theta23)
    assert np.isclose(actual.delta, expected.delta)


def _canonical_angle(angle: float) -> float:
    return float(((float(angle) + np.pi) % (2.0 * np.pi)) - np.pi)


def _identity_rotation() -> RotationParameters:
    return RotationParameters()


def _rotation_from_ckm_observables(observables: np.ndarray) -> RotationParameters:
    vus, vcb, vub, jarlskog = np.round(np.asarray(observables, dtype=float), 5)
    s13 = float(np.clip(abs(vub), 0.0, 1.0))
    c13 = float(np.sqrt(max(1.0 - s13 * s13, 0.0)))
    s12 = float(np.clip(abs(vus) / max(c13, 1e-15), 0.0, 1.0))
    s23 = float(np.clip(abs(vcb) / max(c13, 1e-15), 0.0, 1.0))
    c12 = float(np.sqrt(max(1.0 - s12 * s12, 0.0)))
    c23 = float(np.sqrt(max(1.0 - s23 * s23, 0.0)))
    denominator = max(c12 * c13 * c13 * c23 * s12 * s13 * s23, 1e-30)
    delta = float(np.arcsin(np.clip(jarlskog / denominator, -1.0, 1.0)))
    return RotationParameters(
        theta12=float(np.round(np.arcsin(s12), 5)),
        theta13=float(np.round(np.arcsin(s13), 5)),
        theta23=float(np.round(np.arcsin(s23), 5)),
        delta=float(np.round(((delta + np.pi) % (2.0 * np.pi)) - np.pi, 5)),
    )


def _quotient_equivalent_seed(seed: QuarkFitSeed) -> QuarkFitSeed:
    return QuarkFitSeed(
        up_singular_values=seed.up_singular_values / 7.0,
        down_singular_values=seed.down_singular_values / 7.0,
        overall_scale=seed.overall_scale * 7.0,
        up_left=RotationParameters(
            theta12=seed.up_left.theta12 + 2.0 * np.pi,
            theta13=seed.up_left.theta13 - 2.0 * np.pi,
            theta23=seed.up_left.theta23 + 2.0 * np.pi,
            delta=seed.up_left.delta - 2.0 * np.pi,
        ),
        up_right=RotationParameters(
            theta12=seed.up_right.theta12 + 0.31,
            theta13=seed.up_right.theta13 - 0.27,
            theta23=seed.up_right.theta23 + 0.19,
            delta=seed.up_right.delta - 0.11,
        ),
        down_left=RotationParameters(
            theta12=seed.down_left.theta12 - 2.0 * np.pi,
            theta13=seed.down_left.theta13 + 2.0 * np.pi,
            theta23=seed.down_left.theta23 - 2.0 * np.pi,
            delta=seed.down_left.delta + 2.0 * np.pi,
        ),
        down_right=RotationParameters(
            theta12=seed.down_right.theta12 - 0.17,
            theta13=seed.down_right.theta13 + 0.23,
            theta23=seed.down_right.theta23 - 0.29,
            delta=seed.down_right.delta + 0.41,
        ),
    )


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


def test_fit_quark_sector_returns_the_canonical_quotient_representative():
    solution = fit_quark_sector(default_quark_targets(), seed=default_spurion_seed(), overall_scale=3.0, max_nfev=120)

    assert np.isclose(solution.seed.overall_scale, 1.0)
    assert np.all(np.diff(solution.seed.up_singular_values) >= -1e-12)
    assert np.all(np.diff(solution.seed.down_singular_values) >= -1e-12)
    _assert_rotation_close(solution.seed.up_left, _identity_rotation())
    _assert_rotation_close(solution.seed.up_right, _identity_rotation())
    _assert_rotation_close(solution.seed.down_right, _identity_rotation())
    _assert_rotation_close(solution.seed.down_left, _rotation_from_ckm_observables(solution.result.ckm_observables))
    np.testing.assert_allclose(
        encode_quark_fit_canonical_vector(solution.seed),
        encode_quark_fit_canonical_vector(decode_quark_fit_canonical_vector(encode_quark_fit_canonical_vector(solution.seed))),
        rtol=0.0,
        atol=1e-12,
    )


def test_fitted_ckm_matrix_is_unitary():
    """The extracted CKM matrix should remain unitary after exact diagonalization."""
    solution = fit_quark_sector(default_quark_targets(), overall_scale=3.0, max_nfev=120)
    ckm = solution.result.ckm_matrix
    assert np.allclose(ckm.conjugate().T @ ckm, np.eye(3), atol=1e-10)


def test_fit_quark_sector_is_invariant_under_quotient_directions():
    targets = default_quark_targets()
    base_seed = _fit_seed_from_benchmark()
    shifted_seed = _quotient_equivalent_seed(base_seed)

    base_solution = fit_quark_sector(targets, seed=base_seed, max_nfev=2000)
    shifted_solution = fit_quark_sector(targets, seed=shifted_seed, max_nfev=2000)

    assert np.isclose(base_solution.initial_score, shifted_solution.initial_score, rtol=0.0, atol=1e-9)
    assert np.isclose(base_solution.result.score, shifted_solution.result.score, rtol=0.0, atol=1e-8)
    np.testing.assert_allclose(
        base_solution.result.mass_residuals_up,
        shifted_solution.result.mass_residuals_up,
        rtol=0.0,
        atol=2e-4,
    )
    np.testing.assert_allclose(
        base_solution.result.mass_residuals_down,
        shifted_solution.result.mass_residuals_down,
        rtol=0.0,
        atol=2e-4,
    )
    np.testing.assert_allclose(
        base_solution.result.ckm_residuals,
        shifted_solution.result.ckm_residuals,
        rtol=0.0,
        atol=3e-7,
    )
    np.testing.assert_allclose(
        encode_quark_fit_canonical_vector(base_solution.seed),
        encode_quark_fit_canonical_vector(shifted_solution.seed),
        rtol=0.0,
        atol=5e-6,
    )


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


def test_full_vector_shape_and_field_order_are_deterministic():
    seed = _fit_seed_from_benchmark()
    vector = encode_quark_fit_canonical_vector(seed)

    assert QUARK_FIT_CANONICAL_VECTOR_FIELDS == (
        "log_up_physical_singular_values[0]",
        "log_up_physical_singular_values[1]",
        "log_up_physical_singular_values[2]",
        "log_down_physical_singular_values[0]",
        "log_down_physical_singular_values[1]",
        "log_down_physical_singular_values[2]",
        "up_left.theta12",
        "up_left.theta13",
        "up_left.theta23",
        "up_left.delta",
        "down_left.theta12",
        "down_left.theta13",
        "down_left.theta23",
        "down_left.delta",
    )
    assert vector.shape == (QUARK_FIT_CANONICAL_VECTOR_SIZE,)
    np.testing.assert_allclose(vector[:3], np.log(seed.overall_scale * seed.up_singular_values))
    np.testing.assert_allclose(
        vector[3:6], np.log(seed.overall_scale * seed.down_singular_values)
    )
    np.testing.assert_allclose(
        vector[6:10],
        np.array(
            [
                _canonical_angle(seed.up_left.theta12),
                _canonical_angle(seed.up_left.theta13),
                _canonical_angle(seed.up_left.theta23),
                _canonical_angle(seed.up_left.delta),
            ],
            dtype=float,
        ),
    )
    np.testing.assert_allclose(
        vector[10:14],
        np.array(
            [
                _canonical_angle(seed.down_left.theta12),
                _canonical_angle(seed.down_left.theta13),
                _canonical_angle(seed.down_left.theta23),
                _canonical_angle(seed.down_left.delta),
            ],
            dtype=float,
        ),
    )


def test_canonical_vector_encode_decode_is_idempotent_on_normalized_representative():
    seed = _fit_seed_from_benchmark()
    vector = encode_quark_fit_canonical_vector(seed)
    recovered = decode_quark_fit_canonical_vector(vector)
    round_tripped = encode_quark_fit_canonical_vector(recovered)

    np.testing.assert_allclose(
        recovered.up_singular_values,
        seed.overall_scale * seed.up_singular_values,
    )
    np.testing.assert_allclose(
        recovered.down_singular_values,
        seed.overall_scale * seed.down_singular_values,
    )
    assert np.isclose(recovered.overall_scale, 1.0)
    _assert_rotation_close(recovered.up_left, seed.up_left)
    _assert_rotation_close(recovered.down_left, seed.down_left)
    assert np.isclose(recovered.up_right.theta12, 0.0)
    assert np.isclose(recovered.up_right.theta13, 0.0)
    assert np.isclose(recovered.up_right.theta23, 0.0)
    assert np.isclose(recovered.up_right.delta, 0.0)
    assert np.isclose(recovered.down_right.theta12, 0.0)
    assert np.isclose(recovered.down_right.theta13, 0.0)
    assert np.isclose(recovered.down_right.theta23, 0.0)
    assert np.isclose(recovered.down_right.delta, 0.0)
    np.testing.assert_allclose(
        round_tripped,
        vector,
        rtol=0.0,
        atol=np.finfo(float).eps,
    )


def test_canonical_vector_quotients_out_common_scale():
    seed = QuarkFitSeed(
        up_singular_values=np.array([0.25, 0.5, 1.0], dtype=float),
        down_singular_values=np.array([0.125, 0.25, 0.75], dtype=float),
        overall_scale=1.0,
        up_left=RotationParameters(theta12=0.1, theta13=0.2, theta23=0.3, delta=0.4),
        up_right=RotationParameters(theta12=0.05, theta13=0.06, theta23=0.07, delta=0.08),
        down_left=RotationParameters(theta12=0.11, theta13=0.12, theta23=0.13, delta=0.14),
        down_right=RotationParameters(theta12=0.09, theta13=0.04, theta23=0.03, delta=0.02),
    )
    altered = QuarkFitSeed(
        up_singular_values=seed.up_singular_values / 8.0,
        down_singular_values=seed.down_singular_values / 8.0,
        overall_scale=seed.overall_scale * 8.0,
        up_left=seed.up_left,
        up_right=seed.up_right,
        down_left=seed.down_left,
        down_right=seed.down_right,
    )

    base_vector = encode_quark_fit_canonical_vector(seed)
    altered_vector = encode_quark_fit_canonical_vector(altered)
    np.testing.assert_allclose(base_vector, altered_vector, rtol=0.0, atol=0.0)

    base_point = build_mfv_point_from_singular_values(
        up_singular_values=seed.up_singular_values,
        down_singular_values=seed.down_singular_values,
        overall_scale=seed.overall_scale,
        r=0.25,
        up_left=seed.up_left,
        up_right=seed.up_right,
        down_left=seed.down_left,
        down_right=seed.down_right,
        label="base-right-rotation-regression",
    )
    altered_point = build_mfv_point_from_singular_values(
        up_singular_values=altered.up_singular_values,
        down_singular_values=altered.down_singular_values,
        overall_scale=altered.overall_scale,
        r=0.25,
        up_left=altered.up_left,
        up_right=altered.up_right,
        down_left=altered.down_left,
        down_right=altered.down_right,
        label="altered-right-rotation-regression",
    )

    base_result = evaluate_quark_fit(base_point, default_quark_targets())
    altered_result = evaluate_quark_fit(altered_point, default_quark_targets())

    assert np.isclose(base_result.score, altered_result.score, rtol=0.0, atol=1e-9)
    np.testing.assert_allclose(
        base_result.masses_up,
        altered_result.masses_up,
        rtol=0.0,
        atol=1e-9,
    )
    np.testing.assert_allclose(
        base_result.masses_down,
        altered_result.masses_down,
        rtol=0.0,
        atol=1e-9,
    )
    np.testing.assert_allclose(
        base_result.ckm_observables,
        altered_result.ckm_observables,
        rtol=0.0,
        atol=2e-9,
    )
    np.testing.assert_allclose(base_result.state.c_u, altered_result.state.c_u, rtol=0.0, atol=1e-9)
    np.testing.assert_allclose(base_result.state.c_d, altered_result.state.c_d, rtol=0.0, atol=1e-9)


def test_canonical_vector_quotients_out_right_rotations_without_changing_observables():
    seed = _fit_seed_from_benchmark()
    altered = QuarkFitSeed(
        up_singular_values=seed.up_singular_values,
        down_singular_values=seed.down_singular_values,
        overall_scale=seed.overall_scale,
        up_left=seed.up_left,
        up_right=RotationParameters(
            theta12=seed.up_right.theta12 + 0.11,
            theta13=seed.up_right.theta13 + 0.07,
            theta23=seed.up_right.theta23 + 0.05,
            delta=seed.up_right.delta + 0.13,
        ),
        down_left=seed.down_left,
        down_right=RotationParameters(
            theta12=seed.down_right.theta12 - 0.09,
            theta13=seed.down_right.theta13 + 0.04,
            theta23=seed.down_right.theta23 + 0.02,
            delta=seed.down_right.delta - 0.15,
        ),
    )

    base_vector = encode_quark_fit_canonical_vector(seed)
    altered_vector = encode_quark_fit_canonical_vector(altered)
    np.testing.assert_allclose(base_vector, altered_vector, rtol=0.0, atol=0.0)

    base_point = build_mfv_point_from_singular_values(
        up_singular_values=seed.up_singular_values,
        down_singular_values=seed.down_singular_values,
        overall_scale=seed.overall_scale,
        r=0.25,
        up_left=seed.up_left,
        up_right=seed.up_right,
        down_left=seed.down_left,
        down_right=seed.down_right,
        label="base-right-rotation-regression",
    )
    altered_point = build_mfv_point_from_singular_values(
        up_singular_values=altered.up_singular_values,
        down_singular_values=altered.down_singular_values,
        overall_scale=altered.overall_scale,
        r=0.25,
        up_left=altered.up_left,
        up_right=altered.up_right,
        down_left=altered.down_left,
        down_right=altered.down_right,
        label="altered-right-rotation-regression",
    )

    base_result = evaluate_quark_fit(base_point, default_quark_targets())
    altered_result = evaluate_quark_fit(altered_point, default_quark_targets())

    assert np.isclose(base_result.score, altered_result.score, rtol=0.0, atol=1e-9)
    np.testing.assert_allclose(
        base_result.masses_up,
        altered_result.masses_up,
        rtol=0.0,
        atol=1e-9,
    )
    np.testing.assert_allclose(
        base_result.masses_down,
        altered_result.masses_down,
        rtol=0.0,
        atol=1e-9,
    )
    np.testing.assert_allclose(
        base_result.ckm_observables,
        altered_result.ckm_observables,
        rtol=0.0,
        atol=2e-9,
    )
    np.testing.assert_allclose(base_result.state.c_u, altered_result.state.c_u, rtol=0.0, atol=1e-9)
    np.testing.assert_allclose(base_result.state.c_d, altered_result.state.c_d, rtol=0.0, atol=1e-9)


def test_canonical_vector_wraps_left_angles_into_fundamental_domain():
    seed = _fit_seed_from_benchmark()
    shifted = QuarkFitSeed(
        up_singular_values=seed.up_singular_values,
        down_singular_values=seed.down_singular_values,
        overall_scale=seed.overall_scale,
        up_left=RotationParameters(
            theta12=seed.up_left.theta12 + 2.0 * np.pi,
            theta13=seed.up_left.theta13,
            theta23=seed.up_left.theta23,
            delta=seed.up_left.delta,
        ),
        up_right=seed.up_right,
        down_left=RotationParameters(
            theta12=seed.down_left.theta12,
            theta13=seed.down_left.theta13,
            theta23=seed.down_left.theta23,
            delta=seed.down_left.delta - 4.0 * np.pi,
        ),
        down_right=seed.down_right,
    )

    base_vector = encode_quark_fit_canonical_vector(seed)
    shifted_vector = encode_quark_fit_canonical_vector(shifted)
    np.testing.assert_allclose(base_vector, shifted_vector, rtol=0.0, atol=1e-15)

    recovered = decode_quark_fit_canonical_vector(shifted_vector)
    assert -np.pi <= recovered.up_left.theta12 < np.pi
    assert -np.pi <= recovered.down_left.theta12 < np.pi
    assert np.isclose(recovered.up_left.theta12, _canonical_angle(seed.up_left.theta12), rtol=0.0, atol=1e-12)
    assert np.isclose(
        recovered.down_left.theta12,
        _canonical_angle(seed.down_left.theta12),
        rtol=0.0,
        atol=1e-12,
    )


def test_fit_orientation_false_remains_deterministic_and_restricted():
    targets = default_quark_targets()
    base_seed = _fit_seed_from_benchmark()
    shifted_seed = _quotient_equivalent_seed(base_seed)

    base_solution = fit_quark_sector(
        targets,
        seed=base_seed,
        max_nfev=2000,
        fit_orientation=False,
    )
    shifted_solution = fit_quark_sector(
        targets,
        seed=shifted_seed,
        max_nfev=2000,
        fit_orientation=False,
    )

    np.testing.assert_allclose(
        base_solution.initial_score,
        shifted_solution.initial_score,
        rtol=0.0,
        atol=1e-9,
    )
    np.testing.assert_allclose(
        base_solution.result.masses_up,
        shifted_solution.result.masses_up,
        rtol=0.0,
        atol=2e-4,
    )
    np.testing.assert_allclose(
        base_solution.result.masses_down,
        shifted_solution.result.masses_down,
        rtol=0.0,
        atol=2e-4,
    )
    np.testing.assert_allclose(
        base_solution.result.ckm_residuals,
        shifted_solution.result.ckm_residuals,
        rtol=0.0,
        atol=3e-7,
    )
    np.testing.assert_allclose(
        encode_quark_fit_canonical_vector(base_solution.seed),
        encode_quark_fit_canonical_vector(shifted_solution.seed),
        rtol=0.0,
        atol=5e-6,
    )
    assert np.isclose(base_solution.seed.overall_scale, 1.0)
    assert np.all(np.diff(base_solution.seed.up_singular_values) >= -1e-12)
    assert np.all(np.diff(base_solution.seed.down_singular_values) >= -1e-12)
    _assert_rotation_close(base_solution.seed.up_left, _identity_rotation())
    _assert_rotation_close(base_solution.seed.up_right, _identity_rotation())
    _assert_rotation_close(base_solution.seed.down_right, _identity_rotation())
    _assert_rotation_close(
        base_solution.seed.down_left,
        _rotation_from_ckm_observables(base_solution.result.ckm_observables),
    )


def test_fit_residuals_reject_zero_ckm_observable_scales():
    targets = QuarkTargets(
        up_masses=np.array([1.0, 2.0, 3.0], dtype=float),
        down_masses=np.array([4.0, 5.0, 6.0], dtype=float),
        ckm=np.eye(3, dtype=np.complex128),
        label="zero-observable-scale-target",
    )
    masses_up = np.array([1.0, 2.0, 3.0], dtype=float)
    masses_down = np.array([4.0, 5.0, 6.0], dtype=float)

    with pytest.raises(ValueError, match="non-zero CKM observable targets"):
        fit_residuals(masses_up, masses_down, np.eye(3, dtype=np.complex128), targets)


def test_canonical_encoder_rejects_non_positive_physical_singular_values():
    seed = QuarkFitSeed(
        up_singular_values=np.array([0.0, 0.5, 1.0], dtype=float),
        down_singular_values=np.array([0.2, 0.3, 0.4], dtype=float),
        overall_scale=1.0,
        up_left=RotationParameters(),
        up_right=RotationParameters(),
        down_left=RotationParameters(),
        down_right=RotationParameters(),
    )

    with pytest.raises(ValueError, match="up physical singular values must be strictly positive"):
        encode_quark_fit_canonical_vector(seed)


def test_target_ckm_observables_match_pdg_2024():
    """Regression test: anchor the target CKM convention against PDG 2024.

    Catches future row/column permutation regressions in the SVD ordering
    or the V = U_uL^dagger U_dL convention. The expected values come from
    numerically evaluating ``ckm_like_unitary(theta12=0.2274, theta13=0.00368,
    theta23=0.0415, delta=1.196)`` and matching against PDG 2024 §12 to ~1sigma.
    """
    targets = default_quark_targets()
    observables = ckm_observables(targets.ckm)
    vus, vcb, vub, jarlskog = observables

    # PDG 2024 §12 averages with ~1sigma tolerance:
    #   |V_us| = 0.22501(50), |V_cb| = 0.04183(150),
    #   |V_ub| = 0.00382(20), J = 3.08e-5(13e-5).
    # The repo target unitary lands within 1sigma of each.
    assert vus == pytest.approx(0.2253, abs=1e-3), f"|V_us| = {vus}"
    assert vcb == pytest.approx(0.0415, abs=2e-3), f"|V_cb| = {vcb}"
    assert vub == pytest.approx(0.00368, abs=3e-4), f"|V_ub| = {vub}"
    assert jarlskog == pytest.approx(3.1e-5, abs=2e-6), f"J = {jarlskog}"
    # Sign of J must be positive (PDG convention; SM CP phase delta > 0).
    assert jarlskog > 0


# ===================================================================
# B2 — PDG SVD rephasing: convention-stability, determinism, invariants
# ===================================================================

def _b2_mass_matrix_pair(seed: int):
    """A hierarchical complex Dirac mass-matrix pair (up, down)."""
    rng = np.random.default_rng(seed)

    def _rand(scales):
        a = rng.standard_normal((3, 3)) + 1j * rng.standard_normal((3, 3))
        return a * np.array(scales)[None, :]

    return _rand([1.0e-3, 0.5, 170.0]), _rand([5.0e-3, 0.1, 4.2])


def test_b2_rephasing_makes_epsilon_k_inputs_convention_stable():
    """The audit's x6 demonstration, inverted into a regression guard: an
    arbitrary PHYSICAL column rephasing of the SVD factors (the raw-SVD gauge
    freedom) must leave the rephased CKM -- and hence Im(M12)/epsilon_K -- bit-
    stable to machine precision."""
    M_u, M_d = _b2_mass_matrix_pair(11)
    U_Lu, _, U_Ru = _ordered_dirac_svd(M_u)
    U_Ld, _, U_Rd = _ordered_dirac_svd(M_d)

    a, _, c, _ = _rephase_to_pdg_convention(U_Lu, U_Ru, U_Ld, U_Rd)
    V_base = a.conj().T @ c

    rng = np.random.default_rng(99)
    P_u = np.diag(np.exp(1j * rng.uniform(0.0, 2.0 * np.pi, 3)))
    P_d = np.diag(np.exp(1j * rng.uniform(0.0, 2.0 * np.pi, 3)))
    a2, _, c2, _ = _rephase_to_pdg_convention(
        U_Lu @ P_u, U_Ru @ P_u, U_Ld @ P_d, U_Rd @ P_d
    )
    V_phased = a2.conj().T @ c2

    # Im(M12) ~ Im of products of CKM entries: invariant CKM => invariant ImM12.
    np.testing.assert_allclose(V_phased, V_base, atol=1.0e-12, rtol=0.0)


def test_b2_rephasing_is_deterministic_bit_identical():
    """Paired-scan guard (the (r, mkk, draw_seed) join): the SAME input through
    the rephasing must give a BIT-identical CKM twice (==, not just rtol), and
    a random physical pre-rephasing must reproduce it bit-identically -- catching
    any np.angle/signed-zero branch nondeterminism."""
    M_u, M_d = _b2_mass_matrix_pair(7)
    U_Lu, _, U_Ru = _ordered_dirac_svd(M_u)
    U_Ld, _, U_Rd = _ordered_dirac_svd(M_d)

    a1, _, c1, _ = _rephase_to_pdg_convention(U_Lu, U_Ru, U_Ld, U_Rd)
    a2, _, c2, _ = _rephase_to_pdg_convention(U_Lu, U_Ru, U_Ld, U_Rd)
    assert np.array_equal(a1.conj().T @ c1, a2.conj().T @ c2)


def test_b2_pdg_convention_entries_real_nonnegative():
    """After the fix, V_ud, V_us, V_cb, V_tb are real and >= 0; the single
    physical phase is carried by V_ub (PLAN §1.3)."""
    M_u, M_d = _b2_mass_matrix_pair(3)
    obs = mass_matrix_observables(M_u, M_d)
    V = obs["ckm"]
    for entry in (V[0, 0], V[0, 1], V[1, 2], V[2, 2]):
        assert abs(entry.imag) < 1.0e-12
        assert entry.real >= 0.0
    # V_ub carries a nonzero physical phase (not pinned real).
    assert abs(V[0, 2]) > 0.0


def test_b2_rephasing_preserves_masses_ckm_magnitudes_and_jarlskog():
    """Invariants pin: masses, |CKM| element-by-element, and the Jarlskog
    invariant are unchanged by the rephasing (rtol <= 1e-12)."""
    M_u, M_d = _b2_mass_matrix_pair(5)
    U_Lu, s_u, U_Ru = _ordered_dirac_svd(M_u)
    U_Ld, s_d, U_Rd = _ordered_dirac_svd(M_d)
    V_raw = U_Lu.conj().T @ U_Ld

    obs = mass_matrix_observables(M_u, M_d)
    np.testing.assert_allclose(obs["masses_up"], s_u, rtol=1.0e-12)
    np.testing.assert_allclose(obs["masses_down"], s_d, rtol=1.0e-12)
    np.testing.assert_allclose(np.abs(obs["ckm"]), np.abs(V_raw), rtol=1.0e-12)
    assert jarlskog_invariant(obs["ckm"]) == pytest.approx(
        jarlskog_invariant(V_raw), rel=1.0e-10
    )


def test_b2_rephasing_is_idempotent():
    """Applying the rephasing to an already-canonical set is a no-op."""
    M_u, M_d = _b2_mass_matrix_pair(13)
    obs = mass_matrix_observables(M_u, M_d)
    a, b, c, d = _rephase_to_pdg_convention(
        obs["U_L_u"], obs["U_R_u"], obs["U_L_d"], obs["U_R_d"]
    )
    np.testing.assert_allclose(a.conj().T @ c, obs["ckm"], atol=1.0e-13, rtol=0.0)
