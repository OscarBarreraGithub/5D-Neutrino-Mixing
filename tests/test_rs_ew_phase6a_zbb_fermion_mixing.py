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


def _cghnp_diagonal_bracket(c_repo: float, f_repo: float) -> float:
    """CGHNP (0807.4937) Z->bb ZMA DIAGONAL bracket, computed independently in
    the test from the convention dictionary (c_CGHNP = -c_repo,
    F^2_CGHNP = 2 f_IR,repo^2) -- NOT a byte-copy of the production line.

    Starting from the CGHNP bracket in its OWN variables,
        B(c_C, F_C) = 1/(1 - 2 c_C) * ( 1/F_C^2 - 1 + F_C^2/(3 + 2 c_C) ),
    substitute c_C = -c_repo and F_C^2 = 2 f_repo^2:
        1/(1 - 2(-c)) = 1/(1 + 2c),  1/F_C^2 = 1/(2 f^2),  F_C^2 = 2 f^2,
        3 + 2 c_C = 3 - 2c.
    => B(c, f) = 1/(1 + 2c) * ( 1/(2 f^2) - 1 + 2 f^2/(3 - 2c) ).

    This is the literature-anchored check (M7 refinement 4): the OLD _manual_B
    was byte-identical to the buggy production bracket, which is exactly why the
    suite passed with the sign-flipped term live.  We derive from the dictionary
    instead so the test is an independent oracle, not a re-pin against the code.
    """
    c_cghnp = -c_repo
    f_cghnp_sq = 2.0 * f_repo * f_repo
    return (1.0 / (1.0 - 2.0 * c_cghnp)) * (
        1.0 / f_cghnp_sq - 1.0 + f_cghnp_sq / (3.0 + 2.0 * c_cghnp)
    )


def _manual_fermion_shift(fit: _QuarkFit, lambda_ir_gev: float) -> tuple[float, float]:
    """Independent CGHNP fermion-KK shift (PLAN §4.2 structure).

    delta g_L^b = +(m_b^2/2 M_KK^2) * B_d,
    delta g_R^b = -(m_b^2/2 M_KK^2) * B_Q, with
      B_d = B_diag(c_d3, f_d3)
            + (1/(2 f_d3^2)) * sum_{i=1,2} |Y_d,3i|^2/|Y_d,33|^2 * 1/(1 + 2 c_d_i),
    and symmetric B_Q (column sum, c_Q, f_q3).  The diagonal (b_R/b_L) term
    uses the CGHNP bracket; the light-generation sum uses the COMMON 1/(2 f_3^2)
    factor, NOT a per-light-gen full bracket (PLAN §4.2 error 1a).

    The flavour-sum denominator is 1/(1 - 2 c_di) in CGHNP variables
    (Eq. 170); through the convention dictionary c_CGHNP = -c_repo it becomes
    1/(1 + 2 c_di,repo) -- the SAME conversion already applied to the diagonal
    bracket (_cghnp_diagonal_bracket).  The previous oracle used the
    un-converted 1/(1 - 2 c_repo), which shared the production bug and so could
    not catch it; it is now the independent, sign-correct check.
    """
    bs = fit.bulk_state
    y_d = bs.Y_d_bulk_basis
    denom = abs(y_d[2, 2]) ** 2
    row_ratio = np.abs(y_d[2, :]) ** 2 / denom
    column_ratio = np.abs(y_d[:, 2]) ** 2 / denom
    f_d3_sq = float(bs.F_d[2]) ** 2
    f_q3_sq = float(bs.F_Q[2]) ** 2

    light_sum_d = sum(
        float(row_ratio[i]) / (1.0 + 2.0 * float(bs.c_d[i])) for i in (0, 1)
    )
    light_sum_q = sum(
        float(column_ratio[i]) / (1.0 + 2.0 * float(bs.c_Q[i])) for i in (0, 1)
    )
    B_d = _cghnp_diagonal_bracket(float(bs.c_d[2]), float(bs.F_d[2])) + (
        1.0 / (2.0 * f_d3_sq)
    ) * light_sum_d
    B_Q = _cghnp_diagonal_bracket(float(bs.c_Q[2]), float(bs.F_Q[2])) + (
        1.0 / (2.0 * f_q3_sq)
    ) * light_sum_q
    prefactor = fit.masses_down[2] ** 2 / (2.0 * lambda_ir_gev**2)
    return float(prefactor * B_d), float(-prefactor * B_Q)


def _zbb_metadata(point):
    return point.extras["rs_ew_couplings"].metadata["zbb_fermion_kk_mixing"]


def test_cghnp_diagonal_bracket_matches_convention_dictionary_and_production():
    """Literature-anchored ABSOLUTE pin for the B1 diagonal bracket (M7).

    The expected value is built from the CGHNP convention dictionary
    (c_CGHNP = -c_repo, F^2_CGHNP = 2 f_repo^2) in-test, NOT read off the
    production function -- the production ``_casagrande_zbb_B_profile`` must
    AGREE with that independent oracle.  This is the pin that would have caught
    the sign-flipped / F^2-missing bracket the suite previously hid.
    """
    from quarkConstraints.rs_ew_couplings import _casagrande_zbb_B_profile

    fit = _fit()
    c_d3 = float(fit.bulk_state.c_d[2])  # 0.20: IR-localized
    f_d3 = float(fit.bulk_state.F_d[2])
    expected = _cghnp_diagonal_bracket(c_d3, f_d3)
    produced = _casagrande_zbb_B_profile(c_d3, f_d3, name="B_d[2]")
    assert produced == pytest.approx(expected, rel=1.0e-14)

    # Closed-form spot value from the dictionary: with c_C = -c, F_C^2 = 2 f^2,
    # B = 1/(1 + 2c) * (1/(2 f^2) - 1 + 2 f^2/(3 - 2c)).  For a UV-localized
    # b_R (c_r > 1/2), 1 + 2c > 0 so the bracket is the correct SIGN (the old
    # 1/(1 - 2c) was negative there, flipping delta g_L^b).
    c_uv = 0.557
    f_uv = float(_fit().bulk_state.F_d[1])  # arbitrary positive overlap
    manual = (1.0 / (1.0 + 2.0 * c_uv)) * (
        1.0 / (2.0 * f_uv * f_uv) - 1.0 + (2.0 * f_uv * f_uv) / (3.0 - 2.0 * c_uv)
    )
    assert _casagrande_zbb_B_profile(c_uv, f_uv, name="uv") == pytest.approx(manual)


def test_casagrande_zbb_fermion_piece_is_smaller_than_gauge_piece():
    """Magnitude guard (M7): the corrected CGHNP fermion piece must NOT exceed
    the independently-verified gauge piece (the old bug inflated it ~30x)."""
    fit = _fit()
    gauge_only = _point(fit, include_fermion_kk_mixing=False).extras["rs_ew_couplings"]
    mixed = _point(fit, include_fermion_kk_mixing=True).extras["rs_ew_couplings"]
    gauge_L = complex(gauge_only.z_delta_g_L_d[2, 2])
    fermion_L = complex(mixed.z_delta_g_L_d[2, 2]) - gauge_L
    gauge_R = complex(gauge_only.z_delta_g_R_d[2, 2])
    fermion_R = complex(mixed.z_delta_g_R_d[2, 2]) - gauge_R
    assert abs(fermion_L) < abs(gauge_L)
    assert abs(fermion_R) < abs(gauge_R)


def test_casagrande_zbb_fermion_shift_signs_and_independent_formula():
    fit = _fit()
    point = _point(fit, include_fermion_kk_mixing=True)
    couplings = point.extras["rs_ew_couplings"]
    zbb = _zbb_metadata(point)
    manual_left, manual_right = _manual_fermion_shift(
        fit,
        point.extras["rs_ew_spectrum"].lambda_ir_gev,
    )

    # Production matches the INDEPENDENT convention-dictionary recomputation
    # (CGHNP c=-c_repo, F^2=2 f^2), the literature-anchored oracle (M7).
    assert zbb["delta_g_L_b"] == pytest.approx(manual_left, rel=1.0e-14, abs=1.0e-18)
    assert zbb["delta_g_R_b"] == pytest.approx(manual_right, rel=1.0e-14, abs=1.0e-18)
    # Re-pinned after B1 (CGHNP retranslation: c-sign + F^2=2f^2 + flavour-sum
    # index).  The OLD pins (2.35e-5, -5.49e-4) were pins-of-the-bug; the
    # corrected fermion piece is sign-correct and much smaller.
    assert zbb["delta_g_L_b"] == pytest.approx(3.7365285051847189e-06)
    assert zbb["delta_g_R_b"] == pytest.approx(-1.9279722705643624e-05)
    # SM-Zbb sign pins (UV-localized b_R, c_r > 1/2): delta g_L^b > 0 (reduces
    # R_b), delta g_R^b < 0.  These FAILED under the old sign-flipped bracket.
    assert zbb["delta_g_L_b"] > 0.0
    assert zbb["delta_g_R_b"] < 0.0
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
    # M1 literature-anchored pin: at the SM-limit prediction the NP shift is 0,
    # so BOTH gates' veto ratios are 0 -> the SM passes by construction (no
    # 0.004sigma R_b knife-edge).  Pre-M1, T010's max|pull|-vs-experiment gate
    # gave a NON-zero SM ratio (~0.996) here -- the gate artifact this fixes.
    assert t011.ratio == pytest.approx(0.0)
    assert t010.ratio == pytest.approx(0.0)
    assert t010.passes is True


def test_minimal_complete_zbb_is_rigorous_vetoing_and_numeric_values_are_unchanged():
    harness = _load_scan_harness()
    failing_point = _point(_fit(), include_fermion_kk_mixing=True)
    t010 = fcc.get("T010").evaluate(failing_point)
    t011 = fcc.get("T011").evaluate(failing_point)

    # Re-pinned after B1 (corrected fermion-KK couplings) AND M1 (T010 R_b veto
    # now uses the T011-style loose-edge NP-shift budget, not max|pull| vs exp).
    # t010.ratio is now |predicted - SM_limit| / hard_veto_budget.
    assert t010.predicted == pytest.approx(0.21337155324655482)
    assert t010.ratio == pytest.approx(1.6745946465205102)
    assert t011.predicted == pytest.approx(0.9374824088474414)
    assert t011.ratio == pytest.approx(0.0616176902505238)

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


def test_zbb_flavour_sum_sign_convention_uv_light_singlets():
    """Regression pin for the CGHNP flavour-sum convention bug (2026-06-23).

    The light-generation flavour sum in B_d/B_Q uses 1/(1 - 2 c_di) in the
    CGHNP convention; with c_CGHNP = -c_repo it MUST become 1/(1 + 2 c_di,repo).
    Using the un-converted 1/(1 - 2 c_repo) made the denominator NEGATIVE for the
    scan-typical UV light singlets (c_repo ~ +0.6) and inflated |delta g_L^b| by
    ~5-8x with the WRONG SIGN.  The default ``_fit`` fixture hid this because its
    off-diagonal Yukawas are ~1e-5 (flavour sum ~1e-10, negligible); this test
    uses O(1) off-diagonal Yukawas + UV light down singlets so the flavour sum is
    O(1) and the sign/magnitude is exercised.
    """
    eps = 1.0e-15
    c_q = np.array([0.63, 0.57, 0.34], dtype=float)   # IR b-doublet (Bauer-S1)
    c_d = np.array([0.66, 0.62, 0.58], dtype=float)   # UV light down singlets
    y_d = np.array(
        [[1.0, 0.3, 0.5], [0.4, 1.0, 0.6], [1.5, 0.9, 1.0]],
        dtype=np.complex128,
    )

    @dataclass(frozen=True)
    class _BS:
        c_Q: np.ndarray
        c_u: np.ndarray
        c_d: np.ndarray
        F_Q: np.ndarray
        F_d: np.ndarray
        Y_d_bulk_basis: np.ndarray

    @dataclass(frozen=True)
    class _FR:
        bulk_state: object
        masses_down: np.ndarray

    bs = _BS(
        c_Q=c_q, c_u=c_q, c_d=c_d,
        F_Q=np.asarray(f_IR(c_q, eps), dtype=float),
        F_d=np.asarray(f_IR(c_d, eps), dtype=float),
        Y_d_bulk_basis=y_d,
    )
    fr = _FR(bulk_state=bs, masses_down=np.array([0.004, 0.09, 2.5]))
    spectrum = RSEWSpectrum.build(
        lambda_ir_gev=3000.0 / 2.45, k_gev=1.2209e19,
        n_gauge_modes=128, quadrature_order=1024, model_label="minimal_rs",
    )
    fk = build_rs_zbb_fermion_kk_mixing(fr, spectrum=spectrum)

    # The fix: positive, bounded delta g_L^b (the SM/Casagrande sign), magnitude
    # comparable to the gauge piece -- NOT a large negative O(0.01) blow-up.
    assert fk.delta_g_L_b == pytest.approx(5.847017e-03, rel=1e-4)
    assert fk.delta_g_L_b > 0.0
    assert 0.0 < fk.delta_g_L_b < 0.02   # bounded, not inflated 5-8x
