"""Physics sanity tests for the epsilon_K and Delta m_K implementation.

These tests verify dimensional analysis, chiral enhancement, CP structure,
scaling behavior, NP budget, and operator dominance of the kaon mixing
observables implemented in ``quarkConstraints.deltaf2``.
"""

import math
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.deltaf2 import (
    B_1_K,
    B_4_K,
    B_5_K,
    DELTA_M_K,
    EPSILON_K_EXP,
    EPSILON_K_SM,
    F_K,
    KAPPA_EPSILON,
    M_D_2GEV,
    M_K,
    M_S_2GEV,
    DeltaF2Input,
    DeltaF2WilsonCoefficients,
    _compute_m12_np,
    _kaon_matrix_elements,
    compute_delta_f2_wilsons,
    delta_f2_epsilon_k_budget_policy,
    evaluate_delta_mk,
    evaluate_epsilon_k,
)


# ---------------------------------------------------------------------------
# Helpers: build synthetic couplings that only populate the (s,d) = (0,1) slot
# ---------------------------------------------------------------------------

def _kaon_input() -> DeltaF2Input:
    """Return a minimal kaon DeltaF2Input for standalone Wilson tests."""
    return DeltaF2Input(
        key="epsilon_k",
        display_name="epsilon_K",
        column_name="epsilon_k_ratio",
        reject_reason="epsilon_k",
        sector="down",
        generations=(0, 1),
        bound=2.0e-8,
    )


def _make_wilsons(
    *,
    left: complex,
    right: complex,
    M_KK: float = 3000.0,
) -> DeltaF2WilsonCoefficients:
    """Build a DeltaF2WilsonCoefficients for the kaon system from explicit
    off-diagonal left and right couplings."""
    prefactor = 1.0 / M_KK**2
    return DeltaF2WilsonCoefficients(
        input=_kaon_input(),
        M_KK=M_KK,
        matching_scale=M_KK,
        left_coupling=left,
        right_coupling=right,
        c1_vll=left * left * prefactor / 6.0,
        c1_vrr=right * right * prefactor / 6.0,
        c4_lr=-(left * right) * prefactor,
        c5_lr=(left * right) * prefactor / 3.0,
    )


def _sd_couplings(
    left: complex,
    right: complex,
    M_KK: float = 3000.0,
) -> QuarkMassBasisCouplings:
    """Return a QuarkMassBasisCouplings with only the (s,d) off-diagonal
    entries populated."""
    zeros = np.zeros((3, 3), dtype=np.complex128)
    left_down = zeros.copy()
    right_down = zeros.copy()
    left_down[0, 1] = left
    left_down[1, 0] = np.conj(left)
    right_down[0, 1] = right
    right_down[1, 0] = np.conj(right)
    return QuarkMassBasisCouplings(
        M_KK=M_KK,
        xi_KK=1.0,
        alpha_s=0.09,
        g_s=1.0,
        left_overlap=zeros,
        right_up_overlap=zeros,
        right_down_overlap=zeros,
        left_up=zeros,
        left_down=left_down,
        right_up=zeros,
        right_down=right_down,
    )


# ===================================================================
# 1. Dimensional analysis
# ===================================================================

class TestDimensionalAnalysis:
    """Matrix elements should have dimensions GeV^3, M_12^NP in GeV,
    and epsilon_K should be dimensionless."""

    def test_kaon_matrix_elements_are_gev_cubed_ballpark(self):
        me = _kaon_matrix_elements()

        # f_K^2 * m_K ~ (0.156)^2 * 0.498 ~ 0.012 GeV^3
        fk2_mk = F_K**2 * M_K
        assert 0.005 < fk2_mk < 0.020, (
            f"f_K^2 * m_K = {fk2_mk:.4e} GeV^3 should be O(0.01)"
        )

        # O_VLL ~ (2/3) * fk2_mk * B_1 ~ few x 10^-3 GeV^3
        assert 1e-4 < me["O1_VLL"] < 1e-1, (
            f"O1_VLL = {me['O1_VLL']:.4e} should be a few x 10^-3 GeV^3"
        )

        # LR matrix elements should also be O(10^-2 -- 10^-1) GeV^3
        # because of the chiral enhancement factor
        for key in ("O4_LR", "O5_LR"):
            assert 1e-3 < me[key] < 1.0, (
                f"{key} = {me[key]:.4e} should be O(10^-2) GeV^3"
            )

    def test_chiral_enhancement_factor_is_about_25(self):
        """The ratio (m_K / (m_s + m_d))^2 should be roughly 25."""
        m_ratio_sq = (M_K / (M_S_2GEV + M_D_2GEV)) ** 2
        assert 20.0 < m_ratio_sq < 30.0, (
            f"Chiral enhancement factor = {m_ratio_sq:.1f}, expected ~25"
        )

    def test_m12_np_has_dimensions_of_gev(self):
        """M_12^NP = C_i [1/GeV^2] * <O_i> [GeV^3] should have dimensions GeV."""
        wilsons = _make_wilsons(left=0.01, right=0.05, M_KK=3000.0)
        m12 = _compute_m12_np(wilsons)

        # Wilsons are O(coupling^2 / M_KK^2) ~ 10^-3 / (9e6) ~ 10^-10 GeV^-2
        # matrix elements are O(10^-2) GeV^3
        # so M_12 ~ 10^-12 GeV, well below M_K ~ 0.5 GeV
        assert abs(m12) < M_K, (
            f"|M_12^NP| = {abs(m12):.3e} GeV should be << M_K = {M_K} GeV"
        )

    def test_epsilon_k_is_dimensionless_and_finite(self):
        """epsilon_K should be a dimensionless finite number.

        Note: even modest off-diagonal couplings at M_KK = 3 TeV can
        produce epsilon_K^NP >> epsilon_K^exp. This is the core physics
        of why RS models are so tightly constrained by kaon mixing. We
        use very small couplings here to verify the dimensionless nature.
        """
        # Use very small couplings so that epsilon_K^NP < 1
        wilsons = _make_wilsons(left=0.001 + 0.0005j, right=0.001 + 0.0002j)
        result = evaluate_epsilon_k(wilsons)
        assert np.isfinite(result.epsilon_k_np)
        # With these tiny couplings, epsilon_K should be a small number
        assert result.epsilon_k_np < 1.0, (
            f"epsilon_K = {result.epsilon_k_np:.3e} should be < 1 for tiny couplings"
        )


# ===================================================================
# 2. Chiral enhancement
# ===================================================================

class TestChiralEnhancement:
    """LR matrix elements should be ~25x larger than VLL due to the
    chiral enhancement factor (m_K / (m_s + m_d))^2."""

    def test_lr_matrix_elements_are_chirally_enhanced_over_vll(self):
        me = _kaon_matrix_elements()

        # O4_LR and O5_LR contain the factor (m_K/(m_s+m_d))^2 ~ 25
        # while O1_VLL does not. The ratio should be O(10-50).
        ratio_o4 = me["O4_LR"] / me["O1_VLL"]
        ratio_o5 = me["O5_LR"] / me["O1_VLL"]

        # The chiral enhancement factor is ~25, but the exact ratio depends
        # on bag parameters and numerical prefactors. Still, LR should be
        # clearly larger.
        assert ratio_o4 > 5.0, (
            f"O4_LR / O1_VLL = {ratio_o4:.1f}, expected large chiral enhancement"
        )
        assert ratio_o5 > 5.0, (
            f"O5_LR / O1_VLL = {ratio_o5:.1f}, expected large chiral enhancement"
        )


# ===================================================================
# 2b. B3 literature-anchored ABSOLUTE pins (GGMS hep-ph/9604387 Eq. 8)
# ===================================================================

class TestB3GGMSMatrixElementCoefficients:
    """Literature-anchored absolute pins for the corrected Delta-F=2 matrix
    elements (M7).  Expected values are built from the GGMS Eq. (8) closed-form
    RATIONALS in-test -- NEVER read off the production ``_kaon_matrix_elements``
    output -- so these are an independent oracle, not a self-pin.  In the
    M12-ready normalization the colour-SINGLET O4 carries the LARGE coefficient
    (R/4 + 1/24), the colour-CROSSED O5 the SMALL (R/12 + 1/8), and O1 is
    (1/3).  The previous code had O4/O5 swapped and each x2 too large.
    """

    def _R(self) -> float:
        return (M_K / (M_S_2GEV + M_D_2GEV)) ** 2

    def test_o4_o5_coefficient_ratio_matches_ggms_closed_form_rational(self):
        R = self._R()
        # GGMS Eq. (8) closed-form rationals, written literally (the oracle):
        coeff_o4 = R / 4.0 + 1.0 / 24.0   # colour-SINGLET, LARGE
        coeff_o5 = R / 12.0 + 1.0 / 8.0   # colour-CROSSED, SMALL
        expected_ratio = coeff_o4 / coeff_o5
        # ~2.85 at R ~ 25.7 (PLAN §2.3 / §6.1).
        assert expected_ratio == pytest.approx(2.853, rel=2.0e-3)
        # Production matrix elements, with the bag factors divided out, must
        # reproduce the same coefficient ratio.
        me = _kaon_matrix_elements()
        produced_ratio = (me["O4_LR"] / B_4_K) / (me["O5_LR"] / B_5_K)
        assert produced_ratio == pytest.approx(expected_ratio, rel=1.0e-12)

    def test_singlet_o4_carries_the_large_coefficient_at_large_R(self):
        # Guards against a re-swap: at R >> 1 the singlet coefficient dominates.
        R = self._R()
        assert R > 1.0
        coeff_o4 = R / 4.0 + 1.0 / 24.0
        coeff_o5 = R / 12.0 + 1.0 / 8.0
        assert coeff_o4 > coeff_o5
        me = _kaon_matrix_elements()
        assert (me["O4_LR"] / B_4_K) > (me["O5_LR"] / B_5_K)

    def test_o1_vll_is_one_third_normalization(self):
        # The M12-ready O1 is (1/3) f^2 m B1 -- half the legacy (2/3).  Required
        # for the SM box to reproduce the textbook epsilon_K master formula
        # (PLAN §2.2 SM-box anchor); this is the pin that catches the missing
        # 1/(2 m_M).
        expected_o1 = (1.0 / 3.0) * F_K**2 * M_K * B_1_K
        me = _kaon_matrix_elements()
        assert me["O1_VLL"] == pytest.approx(expected_o1, rel=1.0e-12)

    def test_absolute_o1_matrix_element_matches_literature_number(self):
        # Absolute M12-ready O1 ME ~ 2.21e-3 GeV^3 (PLAN §6.1, slice-4 M-2).
        me = _kaon_matrix_elements()
        assert me["O1_VLL"] == pytest.approx(2.213e-3, rel=5.0e-3)


# ===================================================================
# 3. CP structure
# ===================================================================

class TestCPStructure:
    """epsilon_K measures CP violation via Im(M_12). Real couplings give
    zero epsilon_K; purely imaginary couplings maximize it."""

    def test_real_couplings_give_zero_epsilon_k(self):
        """If all Wilson coefficients are real, Im(M_12) = 0 so epsilon_K = 0."""
        wilsons = _make_wilsons(left=0.01, right=0.05)
        result = evaluate_epsilon_k(wilsons)

        assert abs(result.im_m12_np) < 1e-30, (
            f"Im(M_12^NP) = {result.im_m12_np:.3e} should be 0 for real couplings"
        )
        assert result.epsilon_k_np < 1e-25, (
            f"epsilon_K = {result.epsilon_k_np:.3e} should be 0 for real couplings"
        )

    def test_purely_imaginary_couplings_maximize_epsilon_k(self):
        """If left and right couplings are purely imaginary, the product
        left*right is real and negative, so C4_LR and C5_LR are real.
        But C1_VLL = left^2 is real (and negative), C1_VRR = right^2 is real
        (and negative). All Wilsons are real, so Im(M_12) = 0 again.

        For maximal epsilon_K we need couplings that make the Wilsons complex.
        A left real + right imaginary coupling makes the LR Wilsons purely
        imaginary, maximizing Im(M_12) relative to Re(M_12).
        """
        # left real, right purely imaginary -> LR products are purely imaginary
        left_real = 0.01
        right_imag = 0.05j

        wilsons_max = _make_wilsons(left=left_real, right=right_imag)
        result_max = evaluate_epsilon_k(wilsons_max)

        # Same magnitudes but both real -> epsilon_K = 0
        wilsons_zero = _make_wilsons(left=left_real, right=0.05)
        result_zero = evaluate_epsilon_k(wilsons_zero)

        assert result_max.epsilon_k_np > result_zero.epsilon_k_np * 1e10, (
            f"Imaginary phase should produce large epsilon_K: "
            f"got {result_max.epsilon_k_np:.3e} vs {result_zero.epsilon_k_np:.3e}"
        )

    def test_complex_couplings_give_nonzero_epsilon_k(self):
        """Generic complex couplings should produce nonzero epsilon_K."""
        wilsons = _make_wilsons(left=0.01 + 0.003j, right=0.05 + 0.02j)
        result = evaluate_epsilon_k(wilsons)

        assert result.epsilon_k_np > 0.0
        assert result.im_m12_np != 0.0


# ===================================================================
# 4. Scaling with M_KK
# ===================================================================

class TestMKKScaling:
    """Wilson coefficients scale as 1/M_KK^2, so epsilon_K^NP ~ 1/M_KK^2."""

    def test_doubling_mkk_reduces_epsilon_k_by_factor_4(self):
        left = 0.01 + 0.005j
        right = 0.05 + 0.02j

        w1 = _make_wilsons(left=left, right=right, M_KK=3000.0)
        w2 = _make_wilsons(left=left, right=right, M_KK=6000.0)

        eps1 = evaluate_epsilon_k(w1)
        eps2 = evaluate_epsilon_k(w2)

        ratio = eps1.epsilon_k_np / eps2.epsilon_k_np
        assert np.isclose(ratio, 4.0, rtol=1e-6), (
            f"Doubling M_KK should reduce epsilon_K by 4x, got ratio = {ratio:.4f}"
        )

    def test_doubling_mkk_reduces_delta_mk_by_factor_4(self):
        left = 0.01 + 0.005j
        right = 0.05 + 0.02j

        w1 = _make_wilsons(left=left, right=right, M_KK=3000.0)
        w2 = _make_wilsons(left=left, right=right, M_KK=6000.0)

        dm1 = evaluate_delta_mk(w1)
        dm2 = evaluate_delta_mk(w2)

        ratio = dm1.abs_m12_np / dm2.abs_m12_np
        assert np.isclose(ratio, 4.0, rtol=1e-6), (
            f"Doubling M_KK should reduce |M_12| by 4x, got ratio = {ratio:.4f}"
        )


# ===================================================================
# 5. NP budget
# ===================================================================

class TestNPBudget:
    """The NP budget should use the shared signed BGS+exp one-sigma policy."""

    def test_np_budget_is_positive_and_correct_order(self):
        policy = delta_f2_epsilon_k_budget_policy()
        central_budget = abs(EPSILON_K_EXP - EPSILON_K_SM)
        combined_sigma = math.sqrt((0.18e-3) ** 2 + (0.011e-3) ** 2)

        assert policy.central_budget == pytest.approx(central_budget)
        assert policy.primary_combined_sigma == pytest.approx(combined_sigma)
        assert policy.budget_lowers_epsilon_k == pytest.approx(
            abs(combined_sigma - central_budget)
        )
        assert policy.budget_raises_epsilon_k == pytest.approx(
            central_budget + combined_sigma
        )
        assert policy.sm_choice_sensitivity == pytest.approx(0.15e-3)
        # exp > SM, so positive
        assert EPSILON_K_EXP > EPSILON_K_SM, (
            "epsilon_K^exp should exceed epsilon_K^SM"
        )

    def test_budget_stored_in_result(self):
        """The EpsilonKResult should carry the direction-selected budget."""
        wilsons = _make_wilsons(left=0.01 + 0.005j, right=0.05 + 0.02j)
        result = evaluate_epsilon_k(wilsons)
        policy = delta_f2_epsilon_k_budget_policy()
        expected_budget, expected_direction = policy.selected_signed_budget(
            result.epsilon_k_np_signed
        )
        assert np.isclose(result.epsilon_k_np_budget, expected_budget, rtol=1e-10)
        assert result.selected_budget_direction == expected_direction
        assert result.central_diagnostic_budget == pytest.approx(
            abs(EPSILON_K_EXP - EPSILON_K_SM)
        )
        assert result.budget_policy_id == policy.policy_id
        assert result.confidence_level == "68.27% one_sigma_sensitivity"


# ===================================================================
# 6. Comparison between Delta m_K and epsilon_K
# ===================================================================

class TestDeltaMKVsEpsilonK:
    """For the same Wilson coefficients, epsilon_K uses Im(M_12)
    while Delta m_K uses |M_12|. So |Im(M_12)| <= |M_12| always."""

    def test_im_m12_leq_abs_m12(self):
        """The imaginary part is bounded by the absolute value."""
        wilsons = _make_wilsons(left=0.01 + 0.005j, right=0.05 + 0.02j)

        eps_result = evaluate_epsilon_k(wilsons)
        dm_result = evaluate_delta_mk(wilsons)

        assert abs(eps_result.im_m12_np) <= dm_result.abs_m12_np * (1.0 + 1e-12), (
            f"|Im(M_12)| = {abs(eps_result.im_m12_np):.3e} should be <= "
            f"|M_12| = {dm_result.abs_m12_np:.3e}"
        )

    def test_purely_real_m12_gives_zero_epsilon_k_but_nonzero_delta_mk(self):
        """Real couplings -> real M_12 -> epsilon_K = 0 but Delta m_K != 0."""
        wilsons = _make_wilsons(left=0.01, right=0.05)

        eps_result = evaluate_epsilon_k(wilsons)
        dm_result = evaluate_delta_mk(wilsons)

        assert eps_result.epsilon_k_np < 1e-25
        assert dm_result.abs_m12_np > 0.0, (
            "Real couplings should still give nonzero |M_12|"
        )


# ===================================================================
# 7. Operator dominance
# ===================================================================

class TestOperatorDominance:
    """The C4_LR operator should dominate over VLL by roughly the chiral
    enhancement factor (~25) when left and right couplings are comparable."""

    def test_c4_lr_contribution_dominates_over_vll(self):
        """When left ~ right, the C4_LR contribution to M_12 should dominate
        because of chiral enhancement of the LR matrix elements."""
        left = 0.03
        right = 0.03
        M_KK = 3000.0
        prefactor = 1.0 / M_KK**2
        me = _kaon_matrix_elements()

        # VLL contribution: |C1_VLL * O1_VLL| = |left^2/(6*M_KK^2) * O1_VLL|
        vll_contrib = abs(left * left * prefactor / 6.0 * me["O1_VLL"])

        # C4_LR contribution: |C4_LR * O4_LR| = |left*right/M_KK^2 * O4_LR|
        c4_contrib = abs(left * right * prefactor * me["O4_LR"])

        ratio = c4_contrib / vll_contrib
        # Should be roughly 6 * O4_LR/O1_VLL ~ 6 * chiral_enhancement
        # The factor 6 comes from the 1/6 in C1_VLL
        assert ratio > 10.0, (
            f"C4_LR/C1_VLL contribution ratio = {ratio:.1f}, "
            f"expected >> 1 from chiral enhancement"
        )


# ===================================================================
# 8. Smoke test with a real fitted point
# ===================================================================

class TestFittedPointSmokeTest:
    """Use the existing fitter to get a benchmark point, compute couplings,
    and evaluate epsilon_K. Results should be physically sensible."""

    @pytest.fixture(scope="class")
    def fitted_result(self):
        """Run the quark sector fit once for all tests in this class."""
        from quarkConstraints.benchmarks import default_quark_targets, default_spurion_seed
        from quarkConstraints.fit import fit_quark_sector

        seed = default_spurion_seed()
        solution = fit_quark_sector(
            default_quark_targets(),
            r=0.2,
            seed=seed,
            overall_scale=seed.overall_scale,
            max_nfev=120,
        )
        return solution.result

    @pytest.fixture(scope="class")
    def kaon_wilsons(self, fitted_result):
        """Compute kaon Wilson coefficients from the fitted point."""
        wilsons_all = compute_delta_f2_wilsons(fitted_result, M_KK=3000.0)
        # Find the kaon system
        for w in wilsons_all:
            if w.input.key == "epsilon_k":
                return w
        pytest.fail("No epsilon_k Wilson coefficients found")

    def test_epsilon_k_result_is_finite(self, kaon_wilsons):
        result = evaluate_epsilon_k(kaon_wilsons)
        assert np.isfinite(result.epsilon_k_np)
        assert np.isfinite(result.im_m12_np)
        assert np.isfinite(result.ratio_to_budget)

    def test_epsilon_k_ratio_to_budget_is_positive(self, kaon_wilsons):
        result = evaluate_epsilon_k(kaon_wilsons)
        assert result.ratio_to_budget >= 0.0, (
            f"ratio_to_budget = {result.ratio_to_budget:.3e} should be >= 0"
        )

    def test_delta_mk_result_is_finite(self, kaon_wilsons):
        result = evaluate_delta_mk(kaon_wilsons)
        assert np.isfinite(result.abs_m12_np)
        assert np.isfinite(result.ratio_to_exp)

    def test_m12_np_is_much_smaller_than_meson_mass(self, kaon_wilsons):
        """M_12^NP should be many orders of magnitude below m_K."""
        m12 = _compute_m12_np(kaon_wilsons)
        assert abs(m12) < M_K, (
            f"|M_12^NP| = {abs(m12):.3e} GeV should be << m_K = {M_K} GeV"
        )

    def test_epsilon_k_budget_matches_constants(self, kaon_wilsons):
        result = evaluate_epsilon_k(kaon_wilsons)
        expected, _ = delta_f2_epsilon_k_budget_policy().selected_signed_budget(
            result.epsilon_k_np_signed
        )
        assert np.isclose(result.epsilon_k_np_budget, expected)
