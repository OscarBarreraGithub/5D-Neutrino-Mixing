import math
import importlib.util
from dataclasses import dataclass
from pathlib import Path
import sys

import numpy as np
import pytest

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import point_builder
from quarkConstraints.rs_ew_couplings import (
    DEFAULT_A_REF_C,
    build_rs_zbb_fermion_kk_mixing,
)
from quarkConstraints.rs_ew_spectrum import RSEWSpectrum
from warpConfig.wavefuncs import f_IR


GAUGE_ROOT_EPS_1E_MINUS_15 = 2.450509663813736
EPSILON_RS = 1.0e-15
N_GAUGE_MODES = 64
QUADRATURE_ORDER = 512
MIN_OVERLAP_MODES = 16
MAX_OVERLAP_MODES = 64
OVERLAP_REL_TOL = 1.0e-3
REPO_ROOT = Path(__file__).resolve().parents[1]
SCAN_SCRIPT_PATH = REPO_ROOT / "scripts" / "run_full_catalog_scan.py"


def _load_scan_harness():
    spec = importlib.util.spec_from_file_location(
        "run_full_catalog_scan_zbb_retag",
        SCAN_SCRIPT_PATH,
    )
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@dataclass(frozen=True)
class _BulkState:
    c_Q: np.ndarray
    c_u: np.ndarray
    c_d: np.ndarray
    F_Q: np.ndarray
    F_d: np.ndarray
    Y_d_bulk_basis: np.ndarray


@dataclass(frozen=True)
class _QuarkFit:
    bulk_state: _BulkState
    U_L_u: np.ndarray
    U_L_d: np.ndarray
    U_R_u: np.ndarray
    U_R_d: np.ndarray
    masses_down: np.ndarray


def _scales_for_mkk(mkk_gev: float) -> tuple[float, float]:
    lambda_ir = float(mkk_gev) / GAUGE_ROOT_EPS_1E_MINUS_15
    return lambda_ir, lambda_ir / EPSILON_RS


def _spectrum(mkk_gev: float = 3000.0) -> RSEWSpectrum:
    lambda_ir, k = _scales_for_mkk(mkk_gev)
    return RSEWSpectrum.build(
        lambda_ir_gev=lambda_ir,
        k_gev=k,
        n_gauge_modes=N_GAUGE_MODES,
        quadrature_order=QUADRATURE_ORDER,
        model_label="minimal_rs",
    )


def _fit(*, m_b_gev: float = 4.18, universal_sm: bool = False) -> _QuarkFit:
    if universal_sm:
        c_q = np.full(3, DEFAULT_A_REF_C, dtype=float)
        c_u = np.full(3, DEFAULT_A_REF_C, dtype=float)
        c_d = np.full(3, DEFAULT_A_REF_C, dtype=float)
        y_d = np.diag([0.1, 0.2, 1.0]).astype(np.complex128)
    else:
        c_q = np.array([0.64, 0.56, 0.43], dtype=float)
        c_u = np.array([0.62, 0.34, 0.18], dtype=float)
        c_d = np.array([0.66, 0.57, 0.20], dtype=float)
        y_d = np.array(
            [
                [1.0, 0.0, 1.0e-5],
                [0.0, 1.0, 2.0e-5],
                [3.0e-5, 4.0e-5, 1.0],
            ],
            dtype=np.complex128,
        )
    identity = np.eye(3, dtype=np.complex128)
    return _QuarkFit(
        bulk_state=_BulkState(
            c_Q=c_q,
            c_u=c_u,
            c_d=c_d,
            F_Q=np.asarray(f_IR(c_q, EPSILON_RS), dtype=float),
            F_d=np.asarray(f_IR(c_d, EPSILON_RS), dtype=float),
            Y_d_bulk_basis=y_d,
        ),
        U_L_u=identity,
        U_L_d=identity,
        U_R_u=identity,
        U_R_d=identity,
        masses_down=np.array([0.004, 0.09, float(m_b_gev)], dtype=float),
    )


def _point(
    fit: _QuarkFit,
    *,
    mkk_gev: float = 3000.0,
    include_fermion_kk_mixing: bool,
):
    lambda_ir, k = _scales_for_mkk(mkk_gev)
    return point_builder.build_from_rs_ew_inputs(
        fit,
        Lambda_IR=lambda_ir,
        k=k,
        n_gauge_modes=N_GAUGE_MODES,
        quadrature_order=QUADRATURE_ORDER,
        min_overlap_modes=MIN_OVERLAP_MODES,
        max_overlap_modes=MAX_OVERLAP_MODES,
        overlap_rel_tol=OVERLAP_REL_TOL,
        include_fermion_kk_mixing=include_fermion_kk_mixing,
    )


def _manual_B(c: float, F: float) -> float:
    return (1.0 / (1.0 - 2.0 * c)) * (
        1.0 / (F * F) - 1.0 + (F * F) / (3.0 + 2.0 * c)
    )


def _manual_fermion_shift(fit: _QuarkFit, lambda_ir_gev: float) -> tuple[float, float]:
    y_d = fit.bulk_state.Y_d_bulk_basis
    denom = abs(y_d[2, 2]) ** 2
    row_ratio = np.abs(y_d[2, :]) ** 2 / denom
    column_ratio = np.abs(y_d[:, 2]) ** 2 / denom
    profile_d = np.array(
        [_manual_B(c, f) for c, f in zip(fit.bulk_state.c_d, fit.bulk_state.F_d)]
    )
    profile_q = np.array(
        [_manual_B(c, f) for c, f in zip(fit.bulk_state.c_Q, fit.bulk_state.F_Q)]
    )
    prefactor = fit.masses_down[2] ** 2 / (2.0 * lambda_ir_gev**2)
    return (
        float(prefactor * np.dot(row_ratio, profile_d)),
        float(-prefactor * np.dot(column_ratio, profile_q)),
    )


def _zbb_metadata(point):
    return point.extras["rs_ew_couplings"].metadata["zbb_fermion_kk_mixing"]


def test_casagrande_zbb_fermion_shift_signs_and_independent_formula():
    fit = _fit()
    point = _point(fit, include_fermion_kk_mixing=True)
    couplings = point.extras["rs_ew_couplings"]
    zbb = _zbb_metadata(point)
    manual_left, manual_right = _manual_fermion_shift(
        fit,
        point.extras["rs_ew_spectrum"].lambda_ir_gev,
    )

    assert zbb["delta_g_L_b"] == pytest.approx(manual_left, rel=1.0e-14, abs=1.0e-18)
    assert zbb["delta_g_R_b"] == pytest.approx(manual_right, rel=1.0e-14, abs=1.0e-18)
    assert zbb["delta_g_L_b"] == pytest.approx(2.3518892841117975e-05)
    assert zbb["delta_g_R_b"] == pytest.approx(-5.49194949716315e-04)
    assert zbb["delta_g_L_b"] > 0.0
    assert zbb["delta_g_R_b"] < 0.0
    assert 1.0e-5 < abs(zbb["delta_g_L_b"]) < 1.0e-4
    assert 1.0e-4 < abs(zbb["delta_g_R_b"]) < 1.0e-3
    assert couplings.metadata["fermion_kk_mixing_included"] is True
    assert couplings.metadata["minimal_rs_tree_zbb_complete"] is True
    assert couplings.metadata["custodial_toppartner_zbL_needs_human"] is True
    assert "delta_g_L_b_top_partner" not in zbb


def test_m_b_zero_alignment_leaves_only_gauge_piece_and_no_double_count():
    fit = _fit(m_b_gev=0.0)
    gauge_only = _point(fit, include_fermion_kk_mixing=False).extras["rs_ew_couplings"]
    mixed = _point(fit, include_fermion_kk_mixing=True).extras["rs_ew_couplings"]
    zbb = mixed.metadata["zbb_fermion_kk_mixing"]

    assert zbb["delta_g_L_b"] == pytest.approx(0.0, rel=0.0, abs=0.0)
    assert zbb["delta_g_R_b"] == pytest.approx(0.0, rel=0.0, abs=0.0)
    np.testing.assert_allclose(mixed.z_delta_g_L_d, gauge_only.z_delta_g_L_d)
    np.testing.assert_allclose(mixed.z_delta_g_R_d, gauge_only.z_delta_g_R_d)

    fit = _fit()
    gauge_only = _point(fit, include_fermion_kk_mixing=False).extras["rs_ew_couplings"]
    mixed = _point(fit, include_fermion_kk_mixing=True).extras["rs_ew_couplings"]
    zbb = mixed.metadata["zbb_fermion_kk_mixing"]
    expected_left = np.zeros((3, 3), dtype=np.complex128)
    expected_right = np.zeros((3, 3), dtype=np.complex128)
    expected_left[2, 2] = zbb["delta_g_L_b"]
    expected_right[2, 2] = zbb["delta_g_R_b"]

    np.testing.assert_allclose(mixed.z_delta_g_L_d - gauge_only.z_delta_g_L_d, expected_left)
    np.testing.assert_allclose(mixed.z_delta_g_R_d - gauge_only.z_delta_g_R_d, expected_right)
    np.testing.assert_allclose(mixed.z_total_g_L_d - gauge_only.z_total_g_L_d, expected_left)
    np.testing.assert_allclose(mixed.z_total_g_R_d - gauge_only.z_total_g_R_d, expected_right)


def test_fermion_piece_scales_as_inverse_lambda_ir_squared_and_is_deterministic():
    fit = _fit()
    point_3tev = _point(fit, mkk_gev=3000.0, include_fermion_kk_mixing=True)
    point_6tev = _point(fit, mkk_gev=6000.0, include_fermion_kk_mixing=True)
    repeat = _point(fit, mkk_gev=3000.0, include_fermion_kk_mixing=True)
    zbb_3 = _zbb_metadata(point_3tev)
    zbb_6 = _zbb_metadata(point_6tev)
    couplings = point_3tev.extras["rs_ew_couplings"]
    repeat_couplings = repeat.extras["rs_ew_couplings"]

    assert zbb_6["delta_g_L_b"] / zbb_3["delta_g_L_b"] == pytest.approx(0.25)
    assert zbb_6["delta_g_R_b"] / zbb_3["delta_g_R_b"] == pytest.approx(0.25)
    for name in ("z_delta_g_L_d", "z_delta_g_R_d", "z_total_g_L_d", "z_total_g_R_d"):
        assert np.array_equal(getattr(couplings, name), getattr(repeat_couplings, name))
        assert np.all(np.isfinite(getattr(couplings, name)))


def test_universal_c_m_b_zero_sm_limit_keeps_t010_t011_at_committed_sm_values():
    point = _point(
        _fit(m_b_gev=0.0, universal_sm=True),
        include_fermion_kk_mixing=True,
    )
    couplings = point.extras["rs_ew_couplings"]
    zbb = couplings.metadata["zbb_fermion_kk_mixing"]

    assert np.max(np.abs(couplings.z_delta_g_L_d)) == pytest.approx(0.0, abs=1.0e-18)
    assert np.max(np.abs(couplings.z_delta_g_R_d)) == pytest.approx(0.0, abs=1.0e-18)
    assert zbb["delta_g_L_b"] == pytest.approx(0.0, abs=0.0)
    assert zbb["delta_g_R_b"] == pytest.approx(0.0, abs=0.0)

    t010 = fcc.get("T010").evaluate(point)
    t011 = fcc.get("T011").evaluate(point)
    for result in (t010, t011):
        for value in (
            result.predicted,
            result.sm_prediction,
            result.experimental,
            result.ratio,
            result.budget,
        ):
            assert isinstance(value, float)
            assert math.isfinite(value)
        assert result.diagnostics["minimal_rs_tree_complete"] is True
        assert result.diagnostics["minimal_rs_tree_veto_ready"] is True
        assert result.diagnostics["fermion_kk_mixing_included"] is True
        assert result.diagnostics["custodial_variant_deferred"] is True
        assert result.diagnostics["custodial_toppartner_zbL_deferred"] is True
        assert result.diagnostics["brane_kinetic_terms_deferred"] is True
        assert "needs_human_physics" not in result.diagnostics
        assert result.predicted == pytest.approx(result.sm_prediction)

    assert t010.diagnostics["delta_g_left_b"] == pytest.approx(0.0j)
    assert t010.diagnostics["delta_g_right_b"] == pytest.approx(0.0j)
    assert t011.diagnostics["delta_g_left_b"] == pytest.approx(0.0j)
    assert t011.diagnostics["delta_g_right_b"] == pytest.approx(0.0j)
    assert t011.ratio == pytest.approx(0.0)


def test_minimal_complete_zbb_is_rigorous_vetoing_and_numeric_values_are_unchanged():
    harness = _load_scan_harness()
    failing_point = _point(_fit(), include_fermion_kk_mixing=True)
    t010 = fcc.get("T010").evaluate(failing_point)
    t011 = fcc.get("T011").evaluate(failing_point)

    assert t010.predicted == 0.21328283177386856
    assert t010.ratio == 4.470421078527965
    assert t011.predicted == 0.9383244719075848
    assert t011.ratio == 0.08826504462857736

    for result in (t010, t011):
        assert result.diagnostics["minimal_rs_tree_complete"] is True
        assert result.diagnostics["minimal_rs_tree_veto_ready"] is True
        assert result.diagnostics["custodial_variant_deferred"] is True
        assert "Deferred refinement" in result.diagnostics["custodial_variant_deferred_note"]
        assert "needs_human_physics" not in result.diagnostics
        tag, _, needs_human, proxy_flags = harness.tag_result(result)
        assert tag == "rigorous"
        assert needs_human is None
        assert proxy_flags == {}

    failing_payload = harness._classify_results({"T010": t010})
    assert failing_payload["excluded_by_rigorous"] == ["T010"]
    assert failing_payload["hard_not_evaluated"] == []
    assert failing_payload["constraints"]["T010"]["tag"] == "rigorous"
    assert failing_payload["survives_all_HARD_strict"] is False

    passing_point = _point(
        _fit(m_b_gev=0.0, universal_sm=True),
        include_fermion_kk_mixing=True,
    )
    passing_payload = harness._classify_results(
        {
            "T010": fcc.get("T010").evaluate(passing_point),
            "T011": fcc.get("T011").evaluate(passing_point),
        }
    )
    assert passing_payload["excluded_by_rigorous"] == []
    assert passing_payload["hard_not_evaluated"] == []
    assert passing_payload["survives_all_HARD_strict"] is True
    assert passing_payload["constraints"]["T010"]["tag"] == "rigorous"
    assert passing_payload["constraints"]["T011"]["tag"] == "rigorous"


def test_zbb_fermion_mixing_fails_loud_on_missing_required_inputs():
    spectrum = _spectrum()
    fit = _fit()

    @dataclass(frozen=True)
    class _MissingBulk:
        c_Q: np.ndarray
        c_d: np.ndarray
        F_Q: np.ndarray
        F_d: np.ndarray

    @dataclass(frozen=True)
    class _MissingMassFit:
        bulk_state: _BulkState

    missing_y = _QuarkFit(
        bulk_state=_MissingBulk(
            c_Q=fit.bulk_state.c_Q,
            c_d=fit.bulk_state.c_d,
            F_Q=fit.bulk_state.F_Q,
            F_d=fit.bulk_state.F_d,
        ),
        U_L_u=fit.U_L_u,
        U_L_d=fit.U_L_d,
        U_R_u=fit.U_R_u,
        U_R_d=fit.U_R_d,
        masses_down=fit.masses_down,
    )

    with pytest.raises(ValueError, match="Y_d_bulk_basis"):
        build_rs_zbb_fermion_kk_mixing(missing_y, spectrum=spectrum)

    with pytest.raises(ValueError, match="masses_down"):
        build_rs_zbb_fermion_kk_mixing(_MissingMassFit(fit.bulk_state), spectrum=spectrum)
