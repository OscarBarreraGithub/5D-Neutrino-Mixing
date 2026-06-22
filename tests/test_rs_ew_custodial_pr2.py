from __future__ import annotations

import importlib.util
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import point_builder
from quarkConstraints import oblique_stu as oblique_core
from quarkConstraints.rs_ew_couplings import CUSTODIAL_RS_PLR_EW_MODEL
from tests.rs_ew_phase3b_helpers import (
    GAUGE_ROOT_EPS_1E_MINUS_15,
    MAX_OVERLAP_MODES,
    MIN_OVERLAP_MODES,
    N_GAUGE_MODES,
    OVERLAP_REL_TOL,
    QUADRATURE_ORDER,
    _sample_fit,
    _scales_for_mkk,
)
from warpConfig.wavefuncs import f_IR

REPO_ROOT = Path(__file__).resolve().parents[1]
EPSILON_RS = 1.0e-15


@dataclass(frozen=True)
class _LoopBulkState:
    c_Q: np.ndarray
    c_u: np.ndarray
    c_d: np.ndarray
    F_Q: np.ndarray
    F_u: np.ndarray


@dataclass(frozen=True)
class _LoopFit:
    bulk_state: _LoopBulkState
    U_L_u: np.ndarray
    U_L_d: np.ndarray
    U_R_u: np.ndarray
    U_R_d: np.ndarray
    masses_up: np.ndarray


def _load_scan_harness():
    spec = importlib.util.spec_from_file_location(
        "run_full_catalog_scan_custodial_pr2",
        REPO_ROOT / "scripts" / "run_full_catalog_scan.py",
    )
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _loop_fit() -> _LoopFit:
    c_q = np.array([0.64, 0.56, 0.43], dtype=float)
    c_u = np.array([0.62, 0.34, 0.18], dtype=float)
    c_d = np.array([0.66, 0.57, 0.20], dtype=float)
    identity = np.eye(3, dtype=np.complex128)
    return _LoopFit(
        bulk_state=_LoopBulkState(
            c_Q=c_q,
            c_u=c_u,
            c_d=c_d,
            F_Q=np.asarray(f_IR(c_q, EPSILON_RS), dtype=float),
            F_u=np.asarray(f_IR(c_u, EPSILON_RS), dtype=float),
        ),
        U_L_u=identity,
        U_L_d=identity,
        U_R_u=identity,
        U_R_d=identity,
        masses_up=np.array([0.002, 1.27, 173.0], dtype=float),
    )


def _rs_kwargs(mkk_gev: float = 3000.0) -> dict[str, object]:
    lambda_ir, k = _scales_for_mkk(mkk_gev)
    return {
        "Lambda_IR": lambda_ir,
        "k": k,
        "n_gauge_modes": N_GAUGE_MODES,
        "quadrature_order": QUADRATURE_ORDER,
        "min_overlap_modes": MIN_OVERLAP_MODES,
        "max_overlap_modes": MAX_OVERLAP_MODES,
        "overlap_rel_tol": OVERLAP_REL_TOL,
    }


def _loop_point(**kwargs):
    return point_builder.build_from_rs_ew_inputs(
        _loop_fit(),
        ew_model=CUSTODIAL_RS_PLR_EW_MODEL,
        include_top_partner_loops=True,
        **_rs_kwargs(),
        **kwargs,
    )


def _fcnc_point(**kwargs):
    return point_builder.build_from_rs_ew_inputs(
        _sample_fit(),
        ew_model=CUSTODIAL_RS_PLR_EW_MODEL,
        **_rs_kwargs(),
        **kwargs,
    )


def test_top_partner_loop_proxy_singlet_numeric_oracle_3tev():
    point = _loop_point(
        top_partner_loop_components="singlet",
        top_partner_loop_t_sign=+1,
    )
    couplings = point.extras["rs_ew_couplings"]
    metadata = couplings.metadata
    inputs = metadata["top_partner_loop_inputs"]

    assert inputs["F_Q3"] == pytest.approx(0.2656322304046161)
    assert inputs["F_u3"] == pytest.approx(0.5656854250202850)
    assert inputs["Y_t_eff"] == pytest.approx(3.308347352802294)
    assert metadata["top_partner_delta_t_singlet_magnitude"] == pytest.approx(
        0.15745209098112428
    )
    assert metadata["top_partner_delta_t_loop_applied"] == pytest.approx(
        0.15745209098112428
    )
    assert metadata["top_partner_delta_g_L_b_singlet"] == pytest.approx(
        0.00041018530641991833
    )
    assert metadata["top_partner_delta_g_L_b_loop_applied"] == pytest.approx(
        0.00041018530641991833
    )
    assert metadata["top_partner_delta_g_R_b_loop"] == pytest.approx(0.0)
    assert metadata["top_partner_zbb_loop_numerics_included"] is True
    assert metadata["top_partner_t_loop_numerics_included"] is True
    assert metadata["top_partner_loop_numerics_included"] is True
    assert couplings.z_delta_g_L_d[2, 2] == pytest.approx(
        metadata["z_delta_g_L_d_tree_pr1_b_before_top_partner"]
        + metadata["top_partner_delta_g_L_b_loop_applied"]
    )


def test_t010_singlet_zbb_pseudo_observable_oracle():
    delta_g = 0.00041018530641991833
    z_delta_l = np.zeros((3, 3), dtype=np.complex128)
    z_delta_r = np.zeros((3, 3), dtype=np.complex128)
    z_delta_l[2, 2] = delta_g
    couplings = SimpleNamespace(
        z_delta_g_L_d=z_delta_l,
        z_delta_g_R_d=z_delta_r,
        metadata={
            "ew_model": CUSTODIAL_RS_PLR_EW_MODEL,
            "custodial_protection_included": True,
            "minimal_rs_tree_zbb_complete": True,
            "top_partner_zbb_loop_numerics_included": True,
            "top_partner_loop_numerics_included": True,
            "custodial_toppartner_zbL_needs_human": False,
        },
        matching_assumption="synthetic additive Zbb oracle",
        model_label="synthetic",
        kk_ew_mass_gev=3000.0,
    )

    result = fcc.get("T010").evaluate(point_builder.make_point(rs_ew_couplings=couplings))

    shifted = result.diagnostics["shifted_zbb_couplings"]
    assert shifted["g_left"] == pytest.approx(-0.42242314802691344 + 0.0j)
    assert shifted["g_right"] == pytest.approx(0.07716666666666666 + 0.0j)
    assert result.predicted == pytest.approx(0.21530246427000885)
    assert result.diagnostics["observables"]["A_b"]["predicted"] == pytest.approx(
        0.9354140642148572
    )
    assert result.diagnostics["selected_observable"] == "R_b^0"
    # ratio re-pinned after M1: it is now |predicted - SM_limit| / hard_veto_budget
    # (T011-style loose-edge), not max|pull| vs experiment.  The shifted-coupling
    # / predicted / A_b oracle values above are unchanged (they are the Zbb
    # pseudo-observables, not the gate scalar).  This fixed delta_g now PASSES.
    assert result.ratio == pytest.approx(0.2364937629532263)
    assert result.passes is True


def test_bidoublet_vertex_numeric_oracle_no_fake_negative_t():
    point = _loop_point(
        top_partner_loop_components="bidoublet_vertex",
        top_partner_loop_mixing_scales={"xi_q": 0.5, "xi_chi": 1.0},
    )
    couplings = point.extras["rs_ew_couplings"]
    metadata = couplings.metadata

    assert metadata["top_partner_delta_g_L_b_bidoublet_vertex"] == pytest.approx(
        -0.0003174965326198927
    )
    assert metadata["top_partner_delta_g_L_b_loop_applied"] == pytest.approx(
        -0.0003174965326198927
    )
    assert couplings.z_delta_g_L_d[2, 2] == pytest.approx(
        metadata["z_delta_g_L_d_tree_pr1_b_before_top_partner"]
        + metadata["top_partner_delta_g_L_b_loop_applied"]
    )
    assert metadata["top_partner_delta_t_loop_applied"] == pytest.approx(0.0)
    assert metadata["top_partner_zbb_loop_numerics_included"] is True
    assert metadata["top_partner_t_loop_numerics_included"] is False
    assert metadata["top_partner_loop_numerics_included"] is False
    assert metadata["top_partner_loop_exact_Teq_Zbbeq_included"] is False
    assert metadata["exact_top_partner_mass_matrix_included"] is False
    assert "bidoublet" not in metadata["top_partner_loop_t_source"]


def test_negative_t_requires_numeric_override_and_threads_to_ew001():
    with pytest.raises(ValueError, match="requires a finite"):
        _loop_point(
            top_partner_loop_components="bidoublet_vertex",
            top_partner_loop_t_sign=-1,
        )

    point = _loop_point(
        top_partner_loop_components="bidoublet_vertex",
        top_partner_loop_t_sign=-1,
        top_partner_loop_delta_t_override=-0.05,
    )
    metadata = point.extras["rs_ew_couplings"].metadata
    ew001 = fcc.get("EW001").evaluate(point)

    assert metadata["top_partner_loop_t_source"] == "explicit_numeric_override"
    assert metadata["top_partner_loop_t_sign"] == "-1"
    assert metadata["top_partner_delta_t_loop_applied"] == pytest.approx(-0.05)
    assert metadata["top_partner_t_loop_numerics_included"] is True
    assert ew001.diagnostics["t_loop_prediction"] == pytest.approx(-0.05)
    assert ew001.diagnostics["top_partner_loop_t_source"] == "explicit_numeric_override"
    assert ew001.diagnostics["full_Teq_Zbbeq_loop_matching_included"] is False


def test_top_partner_loop_missing_sign_computes_magnitudes_but_does_not_apply():
    tree_point = point_builder.build_from_rs_ew_inputs(
        _loop_fit(),
        ew_model=CUSTODIAL_RS_PLR_EW_MODEL,
        **_rs_kwargs(),
    )
    missing = _loop_point(top_partner_loop_components="singlet")
    tree = tree_point.extras["rs_ew_couplings"]
    couplings = missing.extras["rs_ew_couplings"]
    metadata = couplings.metadata

    np.testing.assert_array_equal(couplings.z_delta_g_L_d, tree.z_delta_g_L_d)
    np.testing.assert_array_equal(couplings.z_delta_g_R_d, tree.z_delta_g_R_d)
    assert metadata["top_partner_loop_mode"] == "computed_not_applied_missing_t_input"
    assert metadata["top_partner_loop_magnitudes_computed"] is True
    assert metadata["top_partner_zbb_loop_numerics_included"] is False
    assert metadata["top_partner_loop_numerics_included"] is False
    assert metadata["top_partner_delta_t_loop_applied"] == pytest.approx(0.0)
    assert fcc.get("EW001").evaluate(missing).diagnostics["t_prediction"] == pytest.approx(
        fcc.get("EW001").evaluate(tree_point).diagnostics["t_prediction"]
    )


def test_minimal_rs_rejects_top_partner_loops_and_scan_hash_stays_pinned():
    with pytest.raises(ValueError, match="requires ew_model='custodial_rs_plr'"):
        point_builder.build_from_rs_ew_inputs(
            _loop_fit(),
            ew_model="minimal_rs",
            include_top_partner_loops=True,
            top_partner_loop_components="singlet",
            top_partner_loop_t_sign=+1,
            **_rs_kwargs(),
        )

    harness = _load_scan_harness()
    cfg = harness.ScanConfig(
        mkk_values_gev=(1000.0, 3000.0),
        n_draws_per_tile=5,
        xi_kk=2.0,
        base_seed=11,
        tile_seed_stride=17,
    )
    assert harness._config_hash(cfg) == "45e21a07585f7489"


def test_ew001_adds_loop_t_as_separate_term_and_deferred_is_tree_only():
    loop_point = _loop_point(
        top_partner_loop_components="singlet",
        top_partner_loop_t_sign=+1,
    )
    loop_result = fcc.get("EW001").evaluate(loop_point)
    diagnostics = loop_result.diagnostics

    assert diagnostics["t_prediction"] == pytest.approx(
        diagnostics["t_tree_prediction"] + diagnostics["t_loop_prediction"]
    )
    assert diagnostics["t_loop_prediction"] == pytest.approx(0.15745209098112428)
    assert diagnostics["s_prediction"] == pytest.approx(
        fcc.get("EW001")
        .evaluate(point_builder.make_point(kk_ew_mass_gev=3000.0))
        .diagnostics["s_prediction"]
    )
    assert diagnostics["u_prediction"] == pytest.approx(0.0)
    assert diagnostics["top_partner_t_loop_numerics_included"] is True

    missing = _loop_point(top_partner_loop_components="singlet")
    missing_result = fcc.get("EW001").evaluate(missing)
    scale = (oblique_core.DEFAULT_HIGGS_VEV_GEV / 3000.0) ** 2
    assert "t_loop_prediction" not in missing_result.diagnostics
    assert missing_result.diagnostics["t_prediction"] == pytest.approx(
        oblique_core.custodial_rs_plr_t_coefficient() * scale
    )


def test_t010_t011_consume_looped_zbb_matrix_and_report_components():
    point = _loop_point(
        top_partner_loop_components="singlet",
        top_partner_loop_t_sign=+1,
    )
    couplings = point.extras["rs_ew_couplings"]

    for pid in ("T010", "T011"):
        result = fcc.get(pid).evaluate(point)
        diagnostics = result.diagnostics
        assert diagnostics["delta_g_left_b"] == pytest.approx(
            couplings.z_delta_g_L_d[2, 2]
        )
        assert diagnostics["top_partner_zbb_loop_numerics_included"] is True
        assert diagnostics["top_partner_t_loop_numerics_included"] is True
        assert diagnostics["custodial_toppartner_zbL_deferred"] is False
        assert diagnostics["z_delta_g_L_d_tree_pr1_b_before_top_partner"] == pytest.approx(
            0.0j
        )
        assert diagnostics["top_partner_delta_g_L_b_loop_applied"] == pytest.approx(
            0.00041018530641991833
        )


def test_custodial_fcnc_all_gen_bidoublet_zero_lh_mode_t014_rh_minimal_oracle():
    minimal = point_builder.build_from_rs_ew_inputs(
        _sample_fit(),
        ew_model="minimal_rs",
        **_rs_kwargs(),
    ).extras["rs_ew_couplings"]
    point = _fcnc_point(
        custodial_fcnc_mode="all_gen_bidoublet_mass_basis_proxy",
        kappa_fcnc=0.0,
    )
    couplings = point.extras["rs_ew_couplings"]

    offdiag = ~np.eye(3, dtype=bool)
    np.testing.assert_array_equal(couplings.z_delta_g_L_d[offdiag], np.zeros(6))
    np.testing.assert_array_equal(
        couplings.z_delta_g_R_d[offdiag],
        minimal.z_delta_g_R_d[offdiag],
    )
    result = fcc.get("T014").evaluate(point)
    assert result.diagnostics["selected_channel"] == "bs"
    assert result.predicted == pytest.approx(6.369734635064622e-08)
    assert result.ratio == pytest.approx(2.196460218987801e-05)
    assert result.diagnostics["custodial_fcnc_leading_PLR_zeroed"] is True
    assert result.diagnostics["custodial_fcnc_rh_status"] == (
        "right_handed_down_offdiagonal_kept_minimal"
    )


def test_custodial_fcnc_residual_kappa_over_l_keeps_rh_minimal():
    minimal = point_builder.build_from_rs_ew_inputs(
        _sample_fit(),
        ew_model="minimal_rs",
        **_rs_kwargs(),
    ).extras["rs_ew_couplings"]
    zero_point = _fcnc_point(
        custodial_fcnc_mode="all_gen_bidoublet_mass_basis_proxy",
        kappa_fcnc=0.0,
    )
    residual_point = _fcnc_point(
        custodial_fcnc_mode="all_gen_bidoublet_mass_basis_proxy",
        kappa_fcnc=1.0,
    )
    residual = residual_point.extras["rs_ew_couplings"]
    volume = float(residual.metadata["rs_volume_log"])
    offdiag = ~np.eye(3, dtype=bool)

    np.testing.assert_allclose(
        residual.z_delta_g_L_d[offdiag],
        minimal.z_delta_g_L_d[offdiag] / volume,
    )
    np.testing.assert_array_equal(
        residual.z_delta_g_R_d[offdiag],
        minimal.z_delta_g_R_d[offdiag],
    )
    assert fcc.get("T014").evaluate(residual_point).predicted > fcc.get("T014").evaluate(
        zero_point
    ).predicted


def test_honest_omission_flags_present_for_deferred_and_loop_computed_branches():
    deferred = _fcnc_point().extras["rs_ew_couplings"].metadata
    computed = _loop_point(
        top_partner_loop_components="bidoublet_vertex",
        top_partner_loop_mixing_scales={"xi_q": 0.5, "xi_chi": 1.0},
    ).extras["rs_ew_couplings"].metadata

    for metadata in (deferred, computed):
        assert "custodial_omission_details" in metadata
        assert metadata["custodian_spectrum_inferred"] is False
        assert metadata["exact_top_partner_mass_matrix_included"] is False
        assert metadata["brane_kinetic_terms_included"] is False
        assert metadata["full_Teq_Zbbeq_loop_matching_included"] is False
        assert metadata["top_partner_loop_exact_Teq_Zbbeq_included"] is False


def test_top_partner_mass_ratio_metadata_uses_physical_mkk_denominator():
    rho_t = 1.0 / GAUGE_ROOT_EPS_1E_MINUS_15
    point = _loop_point(
        top_partner_loop_components="singlet",
        top_partner_loop_t_sign=+1,
        top_partner_loop_mass_ratios={"rho_t": rho_t},
    )
    metadata = point.extras["rs_ew_couplings"].metadata

    assert metadata["top_partner_loop_mass_ratios"]["rho_t"] == pytest.approx(rho_t)
    assert metadata["top_partner_loop_mass_ratios"]["ratio_denominator"] == (
        "physical_M_KK"
    )
    assert metadata["top_partner_loop_inputs"]["M_KK_convention"] == (
        "physical first electroweak gauge KK mass"
    )
