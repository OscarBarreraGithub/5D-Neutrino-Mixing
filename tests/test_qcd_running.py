"""Tests for quarkConstraints.qcd_running: alpha_s running and Wilson coefficient evolution.

Tests cover:
  - alpha_s running: reference value, asymptotic freedom, known values, threshold
    continuity, positivity
  - Wilson coefficient evolution: near-identity at close scales, VLL/VRR
    enhancement, LR mixing and eigenvalue structure, composition, linearity,
    phase preservation, operator sector isolation, known factors
  - Integration with epsilon_K: running modifies Wilson coefficients consistently

The module implements one-loop alpha_s running and leading-log RG evolution with
the anomalous dimensions gamma_VLL = +4 and the 2x2 LR matrix from Buras,
Misiak, and Urban (NPB 2000).  In the standard BMU convention, VLL/VRR are
*enhanced* when running to lower scales, and the LR sector has one enhanced
eigenvalue and one suppressed eigenvalue.
"""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from quarkConstraints.qcd_running import run_alpha_s, evolve_deltaf2_wilsons


# ---------------------------------------------------------------------------
# Physical constants used in tests
# ---------------------------------------------------------------------------

_ALPHA_S_MZ_MODULE = 0.1179  # internal default of qcd_running
_M_Z = 91.1876
_M_B = 4.18    # GeV, b-quark MS-bar mass
_M_C = 1.27    # GeV, c-quark MS-bar mass
_M_KK = 3000.0  # GeV, typical KK scale for this repo


# ===================================================================
# Section 1: alpha_s running tests
# ===================================================================


class TestRunAlphaS:
    """Tests for the run_alpha_s(mu) one-loop wrapper."""

    def test_alpha_s_at_mz_equals_input_value(self):
        """alpha_s(M_Z) should return the input value (0.1179)."""
        result = run_alpha_s(_M_Z)
        assert abs(result - _ALPHA_S_MZ_MODULE) < 1e-6, (
            f"alpha_s(M_Z) = {result}, expected ~{_ALPHA_S_MZ_MODULE}"
        )

    def test_asymptotic_freedom(self):
        """alpha_s should decrease with increasing mu (asymptotic freedom)."""
        scales = [2.0, 10.0, 91.0, 500.0, 3000.0, 10000.0]
        values = [run_alpha_s(mu) for mu in scales]
        for i in range(len(values) - 1):
            assert values[i] > values[i + 1], (
                f"alpha_s({scales[i]}) = {values[i]} should be > "
                f"alpha_s({scales[i+1]}) = {values[i+1]} (asymptotic freedom)"
            )

    def test_alpha_s_at_1_gev(self):
        """alpha_s(1 GeV) should be roughly 0.3-0.6 (known ballpark)."""
        result = run_alpha_s(1.0)
        assert 0.3 < result < 0.6, (
            f"alpha_s(1 GeV) = {result}, expected in range [0.3, 0.6]"
        )

    def test_alpha_s_at_mb(self):
        """alpha_s(m_b = 4.18 GeV) should be approximately 0.22."""
        result = run_alpha_s(_M_B)
        assert 0.19 < result < 0.25, (
            f"alpha_s(m_b) = {result}, expected in range [0.19, 0.25]"
        )

    def test_alpha_s_continuous_across_bottom_threshold(self):
        """alpha_s should be continuous across the b-quark threshold."""
        delta = 0.01  # GeV
        below = run_alpha_s(_M_B - delta)
        above = run_alpha_s(_M_B + delta)
        assert abs(below - above) < 0.005, (
            f"alpha_s discontinuity at m_b: "
            f"alpha_s({_M_B - delta}) = {below}, "
            f"alpha_s({_M_B + delta}) = {above}, "
            f"jump = {abs(below - above)}"
        )

    def test_alpha_s_continuous_across_charm_threshold(self):
        """alpha_s should be continuous across the c-quark threshold."""
        delta = 0.01  # GeV
        below = run_alpha_s(_M_C - delta)
        above = run_alpha_s(_M_C + delta)
        assert abs(below - above) < 0.01, (
            f"alpha_s discontinuity at m_c: "
            f"alpha_s({_M_C - delta}) = {below}, "
            f"alpha_s({_M_C + delta}) = {above}, "
            f"jump = {abs(below - above)}"
        )

    def test_alpha_s_positive_above_lambda_qcd(self):
        """alpha_s should be positive for all mu > Lambda_QCD (~0.3 GeV)."""
        for mu in [0.5, 1.0, 2.0, 5.0, 10.0, 100.0, 1000.0, 10000.0]:
            result = run_alpha_s(mu)
            assert result > 0, f"alpha_s({mu}) = {result} is not positive"

    def test_alpha_s_at_high_scale(self):
        """alpha_s at a high scale (3 TeV) should be small and perturbative."""
        result = run_alpha_s(_M_KK)
        assert 0.05 < result < 0.12, (
            f"alpha_s(3 TeV) = {result}, expected in range [0.05, 0.12]"
        )

    def test_alpha_s_rejects_non_positive_mu(self):
        """run_alpha_s should raise ValueError for non-positive mu."""
        with pytest.raises(ValueError):
            run_alpha_s(0.0)
        with pytest.raises(ValueError):
            run_alpha_s(-1.0)

    def test_alpha_s_monotone_across_full_range(self):
        """alpha_s should be strictly monotonically decreasing over a dense
        grid spanning all threshold regions.
        """
        scales = np.logspace(np.log10(0.5), np.log10(1e4), 50)
        values = [run_alpha_s(mu) for mu in scales]
        for i in range(len(values) - 1):
            assert values[i] > values[i + 1], (
                f"Monotonicity violated: alpha_s({scales[i]:.3f}) = {values[i]:.6f} "
                f"<= alpha_s({scales[i+1]:.3f}) = {values[i+1]:.6f}"
            )


# ===================================================================
# Section 2: Wilson coefficient evolution tests
# ===================================================================


class TestEvolveDeltaF2Wilsons:
    """Tests for the evolve_deltaf2_wilsons RG evolution.

    In the BMU convention, gamma_VLL = +4. This produces a VLL
    *enhancement* factor when running down (mu_high -> mu_low < mu_high).
    The LR anomalous dimension matrix has one enhanced eigenvalue and one
    suppressed eigenvalue.
    """

    def test_near_identity_at_close_scales(self):
        """Evolving between two very close scales should barely change the
        coefficients.  (The module requires mu_high > mu_low.)
        """
        c_vll = 1.0 + 0.5j
        c_vrr = 0.3 - 0.2j
        c4_lr = -2.0 + 1.0j
        c5_lr = 0.7 - 0.3j
        mu_high = 1000.0
        mu_low = 999.9

        result = evolve_deltaf2_wilsons(c_vll, c_vrr, c4_lr, c5_lr, mu_high, mu_low)

        assert np.isclose(result[0], c_vll, rtol=1e-3), (
            f"C_VLL changed too much: {result[0]} vs {c_vll}"
        )
        assert np.isclose(result[1], c_vrr, rtol=1e-3), (
            f"C_VRR changed too much: {result[1]} vs {c_vrr}"
        )
        assert np.isclose(result[2], c4_lr, rtol=1e-3), (
            f"C4_LR changed too much: {result[2]} vs {c4_lr}"
        )
        assert np.isclose(result[3], c5_lr, rtol=1e-3), (
            f"C5_LR changed too much: {result[3]} vs {c5_lr}"
        )

    def test_mu_high_must_exceed_mu_low(self):
        """The module should reject mu_high <= mu_low."""
        with pytest.raises(ValueError):
            evolve_deltaf2_wilsons(1.0, 0.0, 0.0, 0.0, 2.0, 2.0)
        with pytest.raises(ValueError):
            evolve_deltaf2_wilsons(1.0, 0.0, 0.0, 0.0, 1.0, 3.0)

    def test_vll_enhancement_from_high_to_low_scale(self):
        """C_VLL should be enhanced when running from 3 TeV to 2 GeV.

        With gamma_VLL = +4 in the BMU convention, the leading-log factor is
            (alpha_s(low)/alpha_s(high))^(gamma_VLL/(2*beta_0))
        Since alpha_s(low) > alpha_s(high), the ratio > 1. With gamma_VLL > 0,
        the exponent > 0, so the factor > 1 (enhancement).
        The expected enhancement factor is ~1.3-1.5 from 3 TeV to 2 GeV.
        """
        c_vll_in = 1.0
        result = evolve_deltaf2_wilsons(c_vll_in, 0.0, 0.0, 0.0, _M_KK, 2.0)
        c_vll_out = abs(result[0])

        enhancement = c_vll_out / abs(c_vll_in)
        assert 1.1 < enhancement < 2.0, (
            f"VLL enhancement factor = {enhancement}, expected in [1.1, 2.0]"
        )

    def test_vrr_same_factor_as_vll(self):
        """VRR should receive the same evolution factor as VLL since they
        share the same anomalous dimension gamma = -4.
        """
        result = evolve_deltaf2_wilsons(1.0, 1.0, 0.0, 0.0, _M_KK, 2.0)
        assert np.isclose(abs(result[0]), abs(result[1]), rtol=1e-10), (
            f"VLL factor = {abs(result[0])}, VRR factor = {abs(result[1])}"
        )

    def test_lr_eigenvalue_structure(self):
        """The LR evolution matrix should have one enhanced and one suppressed
        eigenvalue.

        From 3 TeV to 2 GeV, the dominant LR eigenvalue should be ~2.0-3.0
        (enhancement) and the subdominant should be < 1.0 (suppression).
        """
        # Build the 2x2 evolution matrix from two basis runs
        r_c4 = evolve_deltaf2_wilsons(0.0, 0.0, 1.0, 0.0, _M_KK, 2.0)
        r_c5 = evolve_deltaf2_wilsons(0.0, 0.0, 0.0, 1.0, _M_KK, 2.0)
        U = np.array([
            [r_c4[2], r_c5[2]],
            [r_c4[3], r_c5[3]],
        ])

        eigenvalues = np.linalg.eigvals(U)
        eig_abs = sorted([abs(e) for e in eigenvalues])

        # Subdominant: suppressed
        assert eig_abs[0] < 1.0, (
            f"Subdominant LR eigenvalue = {eig_abs[0]}, expected < 1.0"
        )
        # Dominant: enhanced
        assert eig_abs[1] > 1.5, (
            f"Dominant LR eigenvalue = {eig_abs[1]}, expected > 1.5"
        )
        assert eig_abs[1] < 4.0, (
            f"Dominant LR eigenvalue = {eig_abs[1]}, expected < 4.0"
        )

    def test_lr_mixing_c4_generates_c5(self):
        """C4 and C5 mix under running. Starting with C4=1, C5=0 should give
        nonzero C5 at the low scale.
        """
        result = evolve_deltaf2_wilsons(0.0, 0.0, 1.0, 0.0, _M_KK, 2.0)
        assert abs(result[3]) > 1e-3, (
            f"C5_LR should be nonzero after running with C4 input, got {result[3]}"
        )

    def test_lr_mixing_c5_generates_c4(self):
        """Starting with C5=1, C4=0 should give nonzero C4 at the low scale."""
        result = evolve_deltaf2_wilsons(0.0, 0.0, 0.0, 1.0, _M_KK, 2.0)
        assert abs(result[2]) > 1e-3, (
            f"C4_LR should be nonzero after running with C5 input, got {result[2]}"
        )

    def test_vll_does_not_mix_into_lr(self):
        """VLL should not generate LR coefficients."""
        result = evolve_deltaf2_wilsons(1.0, 0.0, 0.0, 0.0, _M_KK, 2.0)
        assert abs(result[2]) < 1e-12, (
            f"VLL should not generate C4_LR, got {result[2]}"
        )
        assert abs(result[3]) < 1e-12, (
            f"VLL should not generate C5_LR, got {result[3]}"
        )

    def test_lr_does_not_mix_into_vll_vrr(self):
        """LR operators should not generate VLL or VRR."""
        result = evolve_deltaf2_wilsons(0.0, 0.0, 1.0, 1.0, _M_KK, 2.0)
        assert abs(result[0]) < 1e-12, (
            f"LR should not generate C_VLL, got {result[0]}"
        )
        assert abs(result[1]) < 1e-12, (
            f"LR should not generate C_VRR, got {result[1]}"
        )

    def test_composition_property(self):
        """Evolving high->mid->low should equal evolving high->low directly.

        We pick an intermediate scale above m_b so that both evolution steps
        remain valid (mu_high > mu_mid > mu_low).
        """
        c_vll = 1.0 + 0.5j
        c_vrr = 0.3 - 0.2j
        c4_lr = -2.0 + 1.0j
        c5_lr = 0.7 - 0.3j
        mu_high = _M_KK
        mu_mid = 10.0
        mu_low = 2.0

        # Two-step evolution
        step1 = evolve_deltaf2_wilsons(c_vll, c_vrr, c4_lr, c5_lr, mu_high, mu_mid)
        step2 = evolve_deltaf2_wilsons(step1[0], step1[1], step1[2], step1[3], mu_mid, mu_low)

        # One-step evolution
        direct = evolve_deltaf2_wilsons(c_vll, c_vrr, c4_lr, c5_lr, mu_high, mu_low)

        for i, name in enumerate(["C_VLL", "C_VRR", "C4_LR", "C5_LR"]):
            assert np.isclose(step2[i], direct[i], rtol=1e-3), (
                f"{name}: two-step = {step2[i]}, direct = {direct[i]}, "
                f"relative diff = {abs(step2[i] - direct[i]) / max(abs(direct[i]), 1e-30)}"
            )

    def test_linearity(self):
        """Doubling the input coefficients should double the output."""
        c_vll = 1.0 + 0.5j
        c_vrr = 0.3 - 0.2j
        c4_lr = -2.0 + 1.0j
        c5_lr = 0.7 - 0.3j

        result_1x = evolve_deltaf2_wilsons(
            c_vll, c_vrr, c4_lr, c5_lr, _M_KK, 2.0
        )
        result_2x = evolve_deltaf2_wilsons(
            2 * c_vll, 2 * c_vrr, 2 * c4_lr, 2 * c5_lr, _M_KK, 2.0
        )

        for i, name in enumerate(["C_VLL", "C_VRR", "C4_LR", "C5_LR"]):
            assert np.isclose(result_2x[i], 2 * result_1x[i], rtol=1e-10), (
                f"{name}: 2x result = {result_2x[i]}, "
                f"2 * 1x result = {2 * result_1x[i]}"
            )

    def test_phase_preservation_for_vll(self):
        """VLL doesn't mix with other operators, so a purely imaginary C_VLL
        should remain purely imaginary after running (just scaled by a real
        factor).
        """
        c_vll = 1.0j  # purely imaginary
        result = evolve_deltaf2_wilsons(c_vll, 0.0, 0.0, 0.0, _M_KK, 2.0)

        assert abs(result[0].real) < 1e-10 * abs(result[0].imag), (
            f"C_VLL phase not preserved: real = {result[0].real}, "
            f"imag = {result[0].imag}"
        )
        # The imaginary part should remain positive (suppressed but same sign)
        assert result[0].imag > 0, (
            f"C_VLL imaginary part should be positive, got {result[0].imag}"
        )

    def test_phase_preservation_for_vrr(self):
        """VRR doesn't mix with other operators, so a purely imaginary C_VRR
        should remain purely imaginary after running.
        """
        c_vrr = 0.5j
        result = evolve_deltaf2_wilsons(0.0, c_vrr, 0.0, 0.0, _M_KK, 2.0)

        assert abs(result[1].real) < 1e-10 * abs(result[1].imag), (
            f"C_VRR phase not preserved: real = {result[1].real}, "
            f"imag = {result[1].imag}"
        )

    def test_known_vll_factor_from_leading_log_formula(self):
        """Verify VLL enhancement factor against an independent leading-log
        calculation using run_alpha_s.

        The factor for each nf segment is:
            (alpha_s(mu_low)/alpha_s(mu_high))^(gamma_VLL/(2*beta_0))
        with gamma_VLL = +4.

        We test the single-segment case (mu > m_b, nf=5) to avoid threshold
        subtleties.
        """
        mu_high = 1000.0
        mu_low = 10.0  # both above m_b, so single nf=5 segment
        gamma_vll = 4.0
        nf = 5
        beta_0 = (33.0 - 2.0 * nf) / 3.0

        alpha_high = run_alpha_s(mu_high)
        alpha_low = run_alpha_s(mu_low)
        expected_factor = (alpha_low / alpha_high) ** (gamma_vll / (2.0 * beta_0))

        result = evolve_deltaf2_wilsons(1.0, 0.0, 0.0, 0.0, mu_high, mu_low)
        actual_factor = abs(result[0])

        assert np.isclose(actual_factor, expected_factor, rtol=1e-6), (
            f"VLL factor: actual = {actual_factor}, expected = {expected_factor}"
        )

    def test_lr_evolution_matrix_determinant(self):
        """The determinant of the 2x2 LR evolution matrix should equal
        the product of eigenvalues from the anomalous dimension matrix.

        For a single segment, det(U_LR) = (alpha_low/alpha_high)^(Tr(gamma_LR)/(2*beta_0)).
        """
        mu_high = 1000.0
        mu_low = 10.0  # single nf=5 segment
        nf = 5
        beta_0 = (33.0 - 2.0 * nf) / 3.0

        # gamma_LR trace: 8 - 28/3 = (24-28)/3 = -4/3
        gamma_lr_trace = 8.0 - 28.0 / 3.0

        alpha_high = run_alpha_s(mu_high)
        alpha_low = run_alpha_s(mu_low)
        expected_det = (alpha_low / alpha_high) ** (gamma_lr_trace / (2.0 * beta_0))

        r_c4 = evolve_deltaf2_wilsons(0.0, 0.0, 1.0, 0.0, mu_high, mu_low)
        r_c5 = evolve_deltaf2_wilsons(0.0, 0.0, 0.0, 1.0, mu_high, mu_low)
        U = np.array([
            [r_c4[2], r_c5[2]],
            [r_c4[3], r_c5[3]],
        ])
        actual_det = abs(np.linalg.det(U))

        assert np.isclose(actual_det, expected_det, rtol=1e-6), (
            f"LR matrix det: actual = {actual_det}, expected = {expected_det}"
        )


# ===================================================================
# Section 3: Integration tests with epsilon_K
# ===================================================================


class TestEpsilonKWithRunning:
    """Integration tests: epsilon_K with and without QCD running.

    In the BMU convention, VLL is enhanced and the LR sector has a mixed
    enhancement/suppression eigenvalue structure. The net effect on
    epsilon_K depends on the relative sizes and phases of the tree-level
    Wilson coefficients. These tests verify that the evolution is applied
    consistently and produces a nontrivial change.
    """

    @pytest.fixture(scope="class")
    def fitted_point(self):
        """Fit the quark sector at r=0.2 and return the fit result."""
        from quarkConstraints.benchmarks import default_quark_targets
        from quarkConstraints.fit import fit_quark_sector

        targets = default_quark_targets()
        solution = fit_quark_sector(targets, r=0.2, max_nfev=200)
        return solution.result

    @pytest.fixture(scope="class")
    def tree_level_wilsons(self, fitted_point):
        """Compute tree-level Wilson coefficients (no running)."""
        from quarkConstraints.couplings import compute_quark_kk_gluon_couplings
        from quarkConstraints.deltaf2 import compute_delta_f2_wilsons

        couplings = compute_quark_kk_gluon_couplings(fitted_point)
        wilsons_tuple = compute_delta_f2_wilsons(couplings)
        for w in wilsons_tuple:
            if w.input.key == "epsilon_k":
                return w
        pytest.skip("No epsilon_k entry found in Wilson coefficients")

    def test_running_modifies_epsilon_k_nontrivially(self, tree_level_wilsons):
        """QCD running should produce a nontrivial change in epsilon_K.

        The ratio should be significantly different from 1 since the
        evolution factors are not unity for the multi-TeV to 2 GeV interval.
        """
        from quarkConstraints.deltaf2 import (
            DeltaF2WilsonCoefficients,
            evaluate_epsilon_k,
        )

        w = tree_level_wilsons
        eps_no_run = evaluate_epsilon_k(w)

        evolved = evolve_deltaf2_wilsons(
            w.c1_vll, w.c1_vrr, w.c4_lr, w.c5_lr,
            w.matching_scale, 2.0,
        )
        w_evolved = DeltaF2WilsonCoefficients(
            input=w.input,
            M_KK=w.M_KK,
            matching_scale=2.0,
            left_coupling=w.left_coupling,
            right_coupling=w.right_coupling,
            c1_vll=evolved[0],
            c1_vrr=evolved[1],
            c4_lr=evolved[2],
            c5_lr=evolved[3],
        )
        eps_with_run = evaluate_epsilon_k(w_evolved)

        if eps_no_run.epsilon_k_np > 0 and eps_with_run.epsilon_k_np > 0:
            ratio = eps_with_run.epsilon_k_np / eps_no_run.epsilon_k_np
            # The running should produce a ratio significantly different from 1
            assert abs(ratio - 1.0) > 0.1, (
                f"epsilon_K ratio = {ratio:.4f}, expected to be significantly "
                f"different from 1.0"
            )
        else:
            # At minimum, the evolved epsilon_K should be finite and non-negative
            assert eps_with_run.epsilon_k_np >= 0, (
                f"epsilon_K with running should be non-negative, "
                f"got {eps_with_run.epsilon_k_np}"
            )

    def test_evolved_wilsons_are_finite_and_well_behaved(self, tree_level_wilsons):
        """All evolved Wilson coefficients should be finite and have magnitudes
        that are within a reasonable range of the tree-level values (no
        runaway or numerical instability).
        """
        w = tree_level_wilsons
        evolved = evolve_deltaf2_wilsons(
            w.c1_vll, w.c1_vrr, w.c4_lr, w.c5_lr,
            w.matching_scale, 2.0,
        )

        for i, name in enumerate(["C_VLL", "C_VRR", "C4_LR", "C5_LR"]):
            assert np.isfinite(evolved[i]), (
                f"{name} is not finite after evolution: {evolved[i]}"
            )

        # The evolved VLL/VRR should be larger in magnitude (enhanced)
        if abs(w.c1_vll) > 1e-30:
            assert abs(evolved[0]) > abs(w.c1_vll), (
                f"VLL should be enhanced: |evolved| = {abs(evolved[0])}, "
                f"|tree| = {abs(w.c1_vll)}"
            )

    def test_running_preserves_linearity_at_physical_point(self, tree_level_wilsons):
        """Scaling the tree-level Wilsons by a constant factor and then
        evolving should give the same result as evolving first and then
        scaling (linearity check at a physical point).
        """
        from quarkConstraints.deltaf2 import (
            DeltaF2WilsonCoefficients,
            evaluate_epsilon_k,
        )

        w = tree_level_wilsons
        scale = 2.5

        # Evolve original, then scale
        evolved_orig = evolve_deltaf2_wilsons(
            w.c1_vll, w.c1_vrr, w.c4_lr, w.c5_lr,
            w.matching_scale, 2.0,
        )

        # Scale, then evolve
        evolved_scaled = evolve_deltaf2_wilsons(
            scale * w.c1_vll, scale * w.c1_vrr,
            scale * w.c4_lr, scale * w.c5_lr,
            w.matching_scale, 2.0,
        )

        for i, name in enumerate(["C_VLL", "C_VRR", "C4_LR", "C5_LR"]):
            assert np.isclose(evolved_scaled[i], scale * evolved_orig[i], rtol=1e-10), (
                f"{name}: scale*evolved = {scale * evolved_orig[i]}, "
                f"evolved(scale*input) = {evolved_scaled[i]}"
            )
