"""Tests for the quark-sector ``Delta F = 2`` observable layer."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

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
    summary = evaluate_delta_f2_constraints(
        _sd_couplings(left=0.0012662771285475794, right=0.04808397328881469)
    )
    kaon = summary.get("epsilon_k")

    assert kaon.coherent_amplitude < kaon.bound
    assert np.isclose(kaon.effective_amplitude, kaon.dominant_operator_size)
    assert kaon.dominant_operator == "C4_LR"
    assert not kaon.passes


def test_default_benchmark_point_has_stable_deltaf2_outputs():
    summary = evaluate_delta_f2_constraints(_fit_result(0.25), M_KK=3000.0)

    assert summary.passes_all
    assert summary.get("epsilon_k").ratio_to_bound < 1.0
    assert summary.get("b_d").ratio_to_bound < 1.0
    assert summary.get("b_s").ratio_to_bound < 1.0
    assert summary.get("d").ratio_to_bound < 1.0
