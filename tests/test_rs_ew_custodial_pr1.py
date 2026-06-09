import importlib.util
import math
import sys
from functools import lru_cache
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import point_builder
from quarkConstraints import oblique_stu as oblique_core
from quarkConstraints.rs_ew_couplings import CUSTODIAL_RS_PLR_EW_MODEL
from tests.rs_ew_phase3b_helpers import (
    MAX_OVERLAP_MODES,
    MIN_OVERLAP_MODES,
    N_GAUGE_MODES,
    OVERLAP_REL_TOL,
    QUADRATURE_ORDER,
    _sample_fit,
    _scales_for_mkk,
)
from tests.test_rs_ew_phase6a_zbb_fermion_mixing import _fit as _zbb_fit

REPO_ROOT = Path(__file__).resolve().parents[1]


_MATRIX_NAMES = (
    "z_delta_g_L_u",
    "z_delta_g_R_u",
    "z_delta_g_L_d",
    "z_delta_g_R_d",
    "z_delta_g_L_e",
    "z_delta_g_R_e",
    "z_delta_g_L_nu",
    "z_total_g_L_u",
    "z_total_g_R_u",
    "z_total_g_L_d",
    "z_total_g_R_d",
    "z_total_g_L_e",
    "z_total_g_R_e",
    "z_total_g_L_nu",
)
_DOWN_FCNC_PAIRS = {"bs": (1, 2), "bd": (0, 2), "sd": (0, 1)}


def _load_scan_harness():
    spec = importlib.util.spec_from_file_location(
        "run_full_catalog_scan_custodial_pr1",
        REPO_ROOT / "scripts" / "run_full_catalog_scan.py",
    )
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _build_point(fit, *, mkk_gev: float, ew_model: str = "minimal_rs", **kwargs):
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
        ew_model=ew_model,
        **kwargs,
    )


@lru_cache(maxsize=None)
def _zbb_point(
    *,
    mkk_gev: float = 3000.0,
    ew_model: str = "minimal_rs",
    custodial_residual: bool = False,
    kappa_b: float = 0.0,
):
    return _build_point(
        _zbb_fit(),
        mkk_gev=mkk_gev,
        ew_model=ew_model,
        include_fermion_kk_mixing=True,
        custodial_PLR_breaking_residual=custodial_residual,
        kappa_b=kappa_b,
    )


@lru_cache(maxsize=None)
def _fcnc_point(*, ew_model: str):
    return _build_point(
        _sample_fit(),
        mkk_gev=3000.0,
        ew_model=ew_model,
        include_fermion_kk_mixing=False,
    )


def _assert_all_coupling_arrays_equal(left, right) -> None:
    for name in _MATRIX_NAMES:
        np.testing.assert_array_equal(getattr(left, name), getattr(right, name))
    assert set(left.neutral_contacts) == set(right.neutral_contacts)
    for key in left.neutral_contacts:
        np.testing.assert_array_equal(left.neutral_contacts[key], right.neutral_contacts[key])
    assert set(left.a_profile_values) == set(right.a_profile_values)
    for key in left.a_profile_values:
        np.testing.assert_array_equal(left.a_profile_values[key], right.a_profile_values[key])
    assert set(left.a_mass_basis) == set(right.a_mass_basis)
    for key in left.a_mass_basis:
        np.testing.assert_array_equal(left.a_mass_basis[key], right.a_mass_basis[key])


def test_default_minimal_rs_is_byte_identical_to_explicit_minimal_and_keeps_scan_hash():
    default = _build_point(_sample_fit(), mkk_gev=3000.0).extras["rs_ew_couplings"]
    explicit = _build_point(
        _sample_fit(),
        mkk_gev=3000.0,
        ew_model="minimal_rs",
    ).extras["rs_ew_couplings"]
    harness = _load_scan_harness()
    cfg = harness.ScanConfig(
        mkk_values_gev=(1000.0, 3000.0),
        n_draws_per_tile=5,
        xi_kk=2.0,
        base_seed=11,
        tile_seed_stride=17,
    )

    _assert_all_coupling_arrays_equal(default, explicit)
    assert dict(default.metadata) == dict(explicit.metadata)
    assert default.metadata["ew_model"] == "minimal_rs"
    assert "ew_model" not in harness._config_payload(cfg)
    assert harness._config_hash(cfg) == "45e21a07585f7489"


def test_custodial_zeroes_only_down_left_diagonal_applies_residual_and_zeroes_b_right():
    minimal = _zbb_point().extras["rs_ew_couplings"]
    custodial = _zbb_point(
        ew_model=CUSTODIAL_RS_PLR_EW_MODEL,
        custodial_residual=True,
        kappa_b=0.7,
    ).extras["rs_ew_couplings"]
    metadata = custodial.metadata
    expected_b_residual = (
        0.7 / float(metadata["rs_volume_log"]) * minimal.z_delta_g_L_d[2, 2]
    )

    assert metadata["custodial_protection_included"] is True
    assert metadata["protected_down_left_diagonal_indices"] == [[0, 0], [1, 1], [2, 2]]
    assert custodial.z_delta_g_L_d[0, 0] == pytest.approx(0.0j, abs=0.0)
    assert custodial.z_delta_g_L_d[1, 1] == pytest.approx(0.0j, abs=0.0)
    assert custodial.z_delta_g_L_d[2, 2] == pytest.approx(expected_b_residual)
    assert metadata["custodial_residual_value"] == pytest.approx(expected_b_residual)
    assert metadata["custodial_residual_applied"] is True
    assert custodial.z_delta_g_R_d[2, 2] == pytest.approx(0.0j, abs=0.0)

    offdiag = ~np.eye(3, dtype=bool)
    np.testing.assert_array_equal(
        custodial.z_delta_g_L_d[offdiag],
        minimal.z_delta_g_L_d[offdiag],
    )
    right_mask = np.ones((3, 3), dtype=bool)
    right_mask[2, 2] = False
    np.testing.assert_array_equal(
        custodial.z_delta_g_R_d[right_mask],
        minimal.z_delta_g_R_d[right_mask],
    )


def test_t014_down_fcnc_offdiagonals_and_result_are_unchanged_in_custodial_branch():
    minimal_point = _fcnc_point(ew_model="minimal_rs")
    custodial_point = _fcnc_point(ew_model=CUSTODIAL_RS_PLR_EW_MODEL)
    minimal = minimal_point.extras["rs_ew_couplings"]
    custodial = custodial_point.extras["rs_ew_couplings"]

    for label, (i, j) in _DOWN_FCNC_PAIRS.items():
        assert abs(minimal.z_delta_g_L_d[i, j]) > 0.0, label
        assert custodial.z_delta_g_L_d[i, j] == minimal.z_delta_g_L_d[i, j]
        assert custodial.z_delta_g_R_d[i, j] == minimal.z_delta_g_R_d[i, j]

    minimal_result = fcc.get("T014").evaluate(minimal_point)
    custodial_result = fcc.get("T014").evaluate(custodial_point)

    assert custodial_result.predicted == pytest.approx(minimal_result.predicted)
    assert custodial_result.ratio == pytest.approx(minimal_result.ratio)
    assert custodial_result.passes is minimal_result.passes
    for label in _DOWN_FCNC_PAIRS:
        assert custodial_result.diagnostics["channels"][label]["delta_g_left"] == (
            pytest.approx(minimal_result.diagnostics["channels"][label]["delta_g_left"])
        )
        assert custodial_result.diagnostics["channels"][label]["ratio"] == (
            pytest.approx(minimal_result.diagnostics["channels"][label]["ratio"])
        )


def test_custodial_oblique_switches_only_t_coefficient_keeps_s_and_u_zero():
    constraint = fcc.get("EW001")
    mkk_gev = 6000.0
    minimal = constraint.evaluate(point_builder.make_point(kk_ew_mass_gev=mkk_gev))
    custodial = constraint.evaluate(
        point_builder.make_point(
            kk_ew_mass_gev=mkk_gev,
            rs_ew_couplings=SimpleNamespace(
                metadata={"ew_model": CUSTODIAL_RS_PLR_EW_MODEL}
            ),
        )
    )
    raw_custodial = constraint.evaluate(
        point_builder.make_point(
            raw={"ew_model": CUSTODIAL_RS_PLR_EW_MODEL},
            kk_ew_mass_gev=mkk_gev,
        )
    )
    scale = (oblique_core.DEFAULT_HIGGS_VEV_GEV / mkk_gev) ** 2

    assert custodial.diagnostics["ew_model"] == CUSTODIAL_RS_PLR_EW_MODEL
    assert custodial.diagnostics["s_prediction"] == pytest.approx(
        minimal.diagnostics["s_prediction"]
    )
    assert custodial.diagnostics["u_prediction"] == pytest.approx(0.0)
    assert custodial.diagnostics["t_prediction"] == pytest.approx(
        oblique_core.custodial_rs_plr_t_coefficient() * scale
    )
    assert custodial.diagnostics["t_coefficient"] == pytest.approx(
        oblique_core.custodial_rs_plr_t_coefficient()
    )
    assert custodial.diagnostics["t_coefficient"] / minimal.diagnostics[
        "t_coefficient"
    ] == pytest.approx(-1.0 / (2.0 * oblique_core.DEFAULT_RS_VOLUME_LOG**2))
    assert abs(custodial.diagnostics["t_prediction"]) < (
        abs(minimal.diagnostics["t_prediction"]) / 1000.0
    )
    assert raw_custodial.diagnostics["t_prediction"] == pytest.approx(
        custodial.diagnostics["t_prediction"]
    )


def test_t010_t011_custodial_tree_is_active_not_deferred_and_stays_rigorous():
    point = _zbb_point(ew_model=CUSTODIAL_RS_PLR_EW_MODEL)
    harness = _load_scan_harness()

    for pid in ("T010", "T011"):
        result = fcc.get(pid).evaluate(point)
        diagnostics = result.diagnostics
        tag, _, needs_human, proxy_flags = harness.tag_result(result)

        assert result.passes is True
        assert diagnostics["ew_model"] == CUSTODIAL_RS_PLR_EW_MODEL
        assert diagnostics["custodial_protection_included"] is True
        assert diagnostics["custodial_variant_deferred"] is False
        assert "custodial_variant_deferred_note" not in diagnostics
        assert diagnostics["custodial_toppartner_zbL_deferred"] is True
        assert diagnostics["top_partner_loop_numerics_included"] is False
        assert diagnostics["qL_rep"] == "all_gen_Q_L_bidoublet_(2,2)_{2/3}"
        assert diagnostics["bR_strategy"] == "elementary_zero"
        assert diagnostics["custodial_omissions"] == {
            "SU2_R_tower": False,
            "custodian_spectrum": False,
            "exact_NC_mixing": False,
            "BKT": False,
        }
        assert "needs_human_physics" not in diagnostics
        assert tag == "rigorous"
        assert needs_human is None
        assert proxy_flags == {}


def test_15tev_minimal_zbb_veto_survives_in_custodial_branch():
    minimal_point = _zbb_point(mkk_gev=15_000.0)
    custodial_point = _zbb_point(
        mkk_gev=15_000.0,
        ew_model=CUSTODIAL_RS_PLR_EW_MODEL,
    )
    minimal_t010 = fcc.get("T010").evaluate(minimal_point)
    custodial_t010 = fcc.get("T010").evaluate(custodial_point)
    custodial_t011 = fcc.get("T011").evaluate(custodial_point)

    assert minimal_t010.passes is False
    assert minimal_t010.ratio == pytest.approx(1.135208212413654)
    assert custodial_t010.passes is True
    assert custodial_t011.passes is True
    assert custodial_t010.ratio < 1.0
    assert custodial_t011.ratio < 1.0
    assert custodial_t010.diagnostics["delta_g_left_b"] == pytest.approx(0.0j)
    assert custodial_t010.diagnostics["delta_g_right_b"] == pytest.approx(0.0j)


def test_custodial_build_is_deterministic_for_arrays_and_zbb_results():
    first = _build_point(
        _zbb_fit(),
        mkk_gev=3000.0,
        ew_model=CUSTODIAL_RS_PLR_EW_MODEL,
        include_fermion_kk_mixing=True,
        custodial_PLR_breaking_residual=True,
        kappa_b=0.7,
    )
    second = _build_point(
        _zbb_fit(),
        mkk_gev=3000.0,
        ew_model=CUSTODIAL_RS_PLR_EW_MODEL,
        include_fermion_kk_mixing=True,
        custodial_PLR_breaking_residual=True,
        kappa_b=0.7,
    )

    _assert_all_coupling_arrays_equal(
        first.extras["rs_ew_couplings"],
        second.extras["rs_ew_couplings"],
    )
    assert fcc.get("T010").evaluate(first) == fcc.get("T010").evaluate(second)
    assert fcc.get("T011").evaluate(first) == fcc.get("T011").evaluate(second)
    assert math.isfinite(
        float(first.extras["rs_ew_couplings"].metadata["custodial_residual_value"].real)
    )
