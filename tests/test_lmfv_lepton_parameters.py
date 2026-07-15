"""Tests for the LMFV lepton carrier used by L001."""

from __future__ import annotations

from dataclasses import replace
import math

import numpy as np
import pytest

from flavorConstraints import (
    GAUGE_KK_ROOT_NN,
    PEREZ_RANDALL_LFV_M_KK_CONVENTION,
    perez_randall_lfv_m_kk_from_lambda_ir,
)
import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import point_builder
from flavor_catalog_constraints.physics_adapters.lepton import (
    LMFVLeptonParameters,
    lmfv_lepton_parameters_from_yukawa_result,
)
from scanParams import ScanConfig, run_scan
from tests.rs_ew_phase3b_helpers import (
    MAX_OVERLAP_MODES,
    MIN_OVERLAP_MODES,
    N_GAUGE_MODES,
    OVERLAP_REL_TOL,
    QUADRATURE_ORDER,
    _sample_fit,
    _scales_for_mkk,
)
from yukawa import compute_all_yukawas


def _benchmark_yukawa(*, lambda_ir: float = 3000.0, k: float | None = None):
    kwargs = {} if k is None else {"k": float(k)}
    return compute_all_yukawas(
        Lambda_IR=float(lambda_ir),
        c_L=0.58,
        c_E=[0.75, 0.60, 0.50],
        c_N=0.27,
        M_N=1.22e18,
        lightest_nu_mass=0.002,
        ordering="normal",
        **kwargs,
    )


def _carrier() -> LMFVLeptonParameters:
    return lmfv_lepton_parameters_from_yukawa_result(
        _benchmark_yukawa(),
        m_kk_gev=3000.0,
    )


def test_carrier_from_yukawa_result_is_readonly_and_records_lmfv_spurion():
    carrier = _carrier()

    for array in (
        carrier.Y_N,
        carrier.Y_N_bar,
        carrier.Y_N_matrix,
        carrier.Y_N_bar_matrix,
        carrier.pmns,
        carrier.lmfv_spurion,
        carrier.c_L,
        carrier.c_E,
        carrier.c_N,
    ):
        assert array.flags.writeable is False
    assert carrier.M_KK_gev == pytest.approx(3000.0)
    assert math.isfinite(carrier.epsilon)
    assert carrier.matching_status == "locked_lmfv_nda_carrier"
    np.testing.assert_allclose(
        carrier.lmfv_spurion,
        carrier.Y_N_bar_matrix @ carrier.Y_N_bar_matrix.conjugate().T,
        rtol=1.0e-12,
        atol=1.0e-14,
    )


@pytest.mark.parametrize(
    "updates, match",
    [
        (
            lambda c: {"Y_N_bar": np.array([np.nan, c.Y_N_bar[1], c.Y_N_bar[2]])},
            "finite",
        ),
        (lambda _c: {"pmns": np.eye(3, dtype=np.complex128) * 2.0}, "unitary"),
        (lambda _c: {"M_KK_gev": -1.0}, "positive finite"),
        (
            lambda c: {"Y_N_bar_matrix": c.Y_N_bar_matrix + np.eye(3) * 1.0e-3},
            "Y_N_bar_matrix",
        ),
    ],
)
def test_carrier_rejects_nonfinite_nonunitary_negative_mass_and_inconsistent_matrix(
    updates,
    match,
):
    carrier = _carrier()

    with pytest.raises(ValueError, match=match):
        replace(carrier, **updates(carrier))


def test_rs_ew_builder_attaches_lmfv_carrier_with_perez_randall_lfv_mkk_convention():
    target_mkk = 4500.0
    lambda_ir, k = _scales_for_mkk(target_mkk)
    point = point_builder.build_from_rs_ew_inputs(
        _sample_fit(),
        Lambda_IR=lambda_ir,
        k=k,
        n_gauge_modes=N_GAUGE_MODES,
        quadrature_order=QUADRATURE_ORDER,
        min_overlap_modes=MIN_OVERLAP_MODES,
        max_overlap_modes=MAX_OVERLAP_MODES,
        overlap_rel_tol=OVERLAP_REL_TOL,
        include_higgs_yukawas=False,
        lepton_sweep_inputs={
            "c_L": 0.58,
            "c_E": [0.75, 0.60, 0.50],
            "c_N": 0.27,
            "M_N": 1.22e18,
            "lightest_nu_mass": 0.002,
            "ordering": "normal",
        },
    )

    assert "lepton_mass_basis_couplings" in point.extras
    assert "lepton_lmfv_parameters" in point.extras
    carrier = point.extras["lepton_lmfv_parameters"]
    spectrum = point.extras["rs_ew_spectrum"]
    # Self-pin for M-14: the L001 Perez-Randall prefactor is calibrated in the
    # geometric LFV convention, while kk_ew_mass_gev remains the physical EW mass.
    assert PEREZ_RANDALL_LFV_M_KK_CONVENTION == "perez_randall_geometric_lambda_ir_v1"
    assert carrier.M_KK_gev == pytest.approx(
        perez_randall_lfv_m_kk_from_lambda_ir(lambda_ir)
    )
    assert carrier.M_KK_gev == pytest.approx(lambda_ir)
    assert carrier.M_KK_gev != pytest.approx(point.extras["kk_ew_mass_gev"])
    assert carrier.M_KK_gev != pytest.approx(spectrum.kk_ew_mass_gev)
    assert point.extras["kk_ew_mass_gev"] == pytest.approx(spectrum.kk_ew_mass_gev)
    assert point.extras["kk_ew_mass_gev"] == pytest.approx(target_mkk)
    assert point.extras["kk_ew_mass_gev"] / carrier.M_KK_gev == pytest.approx(
        GAUGE_KK_ROOT_NN,
        rel=1.0e-3,
    )
    assert math.isfinite(carrier.epsilon)
    assert carrier.Y_N_bar.flags.writeable is False
    np.testing.assert_allclose(
        carrier.lmfv_spurion,
        carrier.Y_N_bar_matrix @ carrier.Y_N_bar_matrix.conjugate().T,
        rtol=1.0e-12,
        atol=1.0e-14,
    )


def test_scan_and_catalog_l001_agree_on_perez_randall_lfv_mkk_convention():
    target_mkk = 4500.0
    lambda_ir, k = _scales_for_mkk(target_mkk)
    lepton_inputs = {
        "c_L": 0.58,
        "c_E": [0.75, 0.60, 0.50],
        "c_N": 0.27,
        "M_N": 1.22e18,
        "lightest_nu_mass": 0.002,
        "ordering": "normal",
    }

    scan_config = ScanConfig(
        record_git_metadata=False,
        k=k,
        Lambda_IR_values=np.array([lambda_ir]),
        c_L_values=np.array([lepton_inputs["c_L"]]),
        c_N_values=np.array([lepton_inputs["c_N"]]),
        c_E_fixed=lepton_inputs["c_E"],
        MN_mode="fixed_ratio",
        MN_over_k=lepton_inputs["M_N"] / k,
        lightest_nu_mass_values=np.array([lepton_inputs["lightest_nu_mass"]]),
    )
    scan_row = run_scan(scan_config, progress_every=0)[0]
    point = point_builder.build_from_rs_ew_inputs(
        _sample_fit(),
        Lambda_IR=lambda_ir,
        k=k,
        n_gauge_modes=N_GAUGE_MODES,
        quadrature_order=QUADRATURE_ORDER,
        min_overlap_modes=MIN_OVERLAP_MODES,
        max_overlap_modes=MAX_OVERLAP_MODES,
        overlap_rel_tol=OVERLAP_REL_TOL,
        include_higgs_yukawas=False,
        lepton_sweep_inputs=lepton_inputs,
    )

    result = fcc.get("L001").evaluate(point)

    assert scan_row["M_KK"] == pytest.approx(lambda_ir)
    assert result.diagnostics["m_kk_convention"] == PEREZ_RANDALL_LFV_M_KK_CONVENTION
    assert result.diagnostics["m_kk_gev"] == pytest.approx(scan_row["M_KK"])
    assert result.diagnostics["dipole_lhs"] == pytest.approx(scan_row["lfv_lhs"])
    assert result.ratio == pytest.approx(scan_row["lfv_ratio"] ** 2)
