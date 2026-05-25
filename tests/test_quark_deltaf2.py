"""Tests for the quark-sector ``Delta F = 2`` observable layer."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import quarkConstraints.deltaf2 as deltaf2
from quarkConstraints.benchmarks import default_quark_targets, default_spurion_seed
from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.deltaf2 import compute_delta_f2_wilsons, evaluate_delta_f2_constraints
from quarkConstraints.fit import fit_quark_sector


def _fit_result(r_value: float):
    seed = default_spurion_seed()
    return fit_quark_sector(
        default_quark_targets(),
        r=r_value,
        seed=seed,
        overall_scale=seed.overall_scale,
        max_nfev=120,
    ).result


def _zero_couplings(M_KK: float = 3000.0) -> QuarkMassBasisCouplings:
    zeros = np.zeros((3, 3), dtype=np.complex128)
    return QuarkMassBasisCouplings(
        M_KK=M_KK,
        xi_KK=1.0,
        alpha_s=0.09,
        g_s=1.0,
        left_overlap=zeros,
        right_up_overlap=zeros,
        right_down_overlap=zeros,
        left_up=zeros,
        left_down=zeros,
        right_up=zeros,
        right_down=zeros,
    )


def _sd_couplings(left: float, right: float, M_KK: float = 3000.0) -> QuarkMassBasisCouplings:
    zeros = np.zeros((3, 3), dtype=np.complex128)
    left_down = zeros.copy()
    right_down = zeros.copy()
    left_down[0, 1] = left_down[1, 0] = left
    right_down[0, 1] = right_down[1, 0] = right
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


def test_deltaf2_wilsons_vanish_for_trivial_aligned_couplings():
    summary = evaluate_delta_f2_constraints(_zero_couplings())

    assert summary.passes_all
    assert np.isclose(summary.worst_ratio, 0.0)
    for observable in summary.observables:
        coeffs = observable.wilsons
        assert coeffs.c1_vll == 0.0j
        assert coeffs.c1_vrr == 0.0j
        assert coeffs.c4_lr == 0.0j
        assert coeffs.c5_lr == 0.0j


def test_deltaf2_wilsons_scale_down_with_larger_mkk():
    result = _fit_result(0.25)
    nominal = compute_delta_f2_wilsons(result, M_KK=3000.0)
    heavier = compute_delta_f2_wilsons(result, M_KK=6000.0)

    for nominal_item, heavier_item in zip(nominal, heavier):
        assert abs(heavier_item.c1_vll) < abs(nominal_item.c1_vll)
        assert abs(heavier_item.c4_lr) < abs(nominal_item.c4_lr)

    nominal_summary = evaluate_delta_f2_constraints(result, M_KK=3000.0)
    heavier_summary = evaluate_delta_f2_constraints(result, M_KK=6000.0)
    assert heavier_summary.worst_ratio < nominal_summary.worst_ratio


def test_deltaf2_pass_fail_uses_dominant_operator_not_coherent_cancellation():
    # With the hadronic evaluation, epsilon_K uses Im(M_12) so real couplings
    # give epsilon_K = 0.  Use the operator-weight fallback to test the old
    # pass/fail invariant that the dominant operator controls acceptance.
    summary = evaluate_delta_f2_constraints(
        _sd_couplings(left=0.0012662771285475794, right=0.04808397328881469),
        use_hadronic=False,
    )
    kaon = summary.get("epsilon_k")

    # The important invariant is that pass/fail is decided by the dominant
    # single operator, not the coherent sum.
    assert np.isclose(kaon.effective_amplitude, kaon.dominant_operator_size)
    assert kaon.dominant_operator == max(
        kaon.weighted_operator_sizes,
        key=kaon.weighted_operator_sizes.get,
    )
    assert not kaon.passes


def test_audited_deltaf2_hadronic_constants_match_selected_sources():
    assert np.isclose(deltaf2.B_1_K, 0.5503, rtol=1e-3)
    assert np.isclose(deltaf2.B_4_K, 0.903, rtol=1e-3)
    assert np.isclose(deltaf2.B_5_K, 0.691, rtol=1e-3)
    assert np.isclose(deltaf2.EPSILON_K_SM, 2.161e-3, rtol=1e-3)


def test_modern_phenomenology_kaon_constants_match_deltaf2_canonical():
    """
    R03-I2: pin the modern-lane vendored kaon constants to the canonical
    deltaf2.py values. The modern lane is intentionally firewalled from
    importing deltaf2 (see tests/test_modern_phenomenology.py), so the
    constants are duplicated; this test enforces they stay in sync.

    The firewall guards against the modern *module* importing deltaf2 at
    source / sys.modules level. Test-side imports of both modules are
    fine because tests are not the firewall surface.
    """
    import quarkConstraints.modern.phenomenology as modern_phen

    # Use np.isclose (not ==) so float literals can differ in trailing digits
    # if a future update introduces more precision on one side.
    pairs = [
        ("_KAON_B_1", "B_1_K"),
        ("_KAON_B_4", "B_4_K"),
        ("_KAON_B_5", "B_5_K"),
        ("_KAON_EPSILON_K_SM", "EPSILON_K_SM"),
        ("_KAON_EPSILON_K_EXP", "EPSILON_K_EXP"),
        ("_KAON_F_K", "F_K"),
        ("_KAON_M_K", "M_K"),
        ("_KAON_DELTA_M_K", "DELTA_M_K"),
        ("_KAON_M_S_2GEV", "M_S_2GEV"),
        ("_KAON_M_D_2GEV", "M_D_2GEV"),
        ("_KAON_KAPPA_EPSILON", "KAPPA_EPSILON"),
    ]
    for modern_name, canonical_name in pairs:
        modern_val = getattr(modern_phen, modern_name)
        canonical_val = getattr(deltaf2, canonical_name)
        assert np.isclose(modern_val, canonical_val), (
            f"{modern_name}={modern_val} (vendored in modern/phenomenology.py) "
            f"out of sync with {canonical_name}={canonical_val} (canonical in deltaf2.py). "
            "Update both literals together."
        )


def test_default_benchmark_point_has_stable_deltaf2_outputs():
    summary = evaluate_delta_f2_constraints(_fit_result(0.25), M_KK=3000.0)

    # The Wilson-RG audit tightened the epsilon_K path; this benchmark is now
    # a stable epsilon_K fail while the other hadronic systems still pass.
    epsilon_k_ratio = summary.get("epsilon_k").ratio_to_bound
    assert np.isclose(epsilon_k_ratio, 1.9286313761001348)
    assert 1.92 <= epsilon_k_ratio < 2.37
    assert summary.get("epsilon_k").ratio_to_bound > 1.0
    assert summary.get("b_d").ratio_to_bound < 1.0
    assert summary.get("b_s").ratio_to_bound < 1.0
    assert summary.get("d").ratio_to_bound < 1.0
    assert not summary.passes_all
