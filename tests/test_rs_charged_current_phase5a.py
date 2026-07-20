import math
from dataclasses import FrozenInstanceError

import numpy as np
import pytest

from flavor_catalog_constraints import point_builder
from quarkConstraints.rs_ew_couplings import DEFAULT_A_REF_C, RSEWNeutralCurrentInputs
from tests.rs_ew_phase3b_helpers import (
    GAUGE_ROOT_EPS_1E_MINUS_15,
    MAX_OVERLAP_MODES,
    MIN_OVERLAP_MODES,
    N_GAUGE_MODES,
    OVERLAP_REL_TOL,
    QUADRATURE_ORDER,
    _identity_fit,
    _sample_fit,
    _scales_for_mkk,
)


def _lepton_inputs(c_l: float):
    return {
        "c_L": float(c_l),
        "c_E": [float(c_l), float(c_l), float(c_l)],
        "c_N": 0.27,
        "M_N": 1.22e18,
        "lightest_nu_mass": 0.002,
        "ordering": "normal",
        "majorana_alpha": 0.0,
        "majorana_beta": 0.0,
    }


def _build_charged_point(fit, *, c_l: float = 0.58, mkk_gev: float = 3000.0):
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
        lepton_sweep_inputs=_lepton_inputs(c_l),
        include_charged_current=True,
    )


@pytest.fixture(scope="module")
def sample_point():
    return _build_charged_point(_sample_fit(), c_l=0.58)


@pytest.fixture(scope="module")
def universal_point():
    ref = np.array([DEFAULT_A_REF_C, DEFAULT_A_REF_C, DEFAULT_A_REF_C], dtype=float)
    return _build_charged_point(_identity_fit(ref, ref, ref), c_l=DEFAULT_A_REF_C)


def test_charged_current_extra_key_and_requested_leptons_fail_loudly():
    marker = object()
    point = point_builder.make_point(rs_charged_current=marker)
    assert point.get_extra("rs_charged_current") is marker

    lambda_ir, k = _scales_for_mkk(3000.0)
    with pytest.raises(ValueError, match="include_charged_current=True requires"):
        point_builder.build_from_rs_ew_inputs(
            _sample_fit(),
            Lambda_IR=lambda_ir,
            k=k,
            n_gauge_modes=N_GAUGE_MODES,
            quadrature_order=QUADRATURE_ORDER,
            min_overlap_modes=MIN_OVERLAP_MODES,
            max_overlap_modes=MAX_OVERLAP_MODES,
            overlap_rel_tol=OVERLAP_REL_TOL,
            include_charged_current=True,
        )


def test_charged_current_inherits_nondefault_neutral_a_ref_c():
    lambda_ir, k = _scales_for_mkk(3000.0)
    point = point_builder.build_from_rs_ew_inputs(
        _sample_fit(),
        Lambda_IR=lambda_ir,
        k=k,
        n_gauge_modes=N_GAUGE_MODES,
        quadrature_order=QUADRATURE_ORDER,
        min_overlap_modes=MIN_OVERLAP_MODES,
        max_overlap_modes=MAX_OVERLAP_MODES,
        overlap_rel_tol=OVERLAP_REL_TOL,
        lepton_sweep_inputs=_lepton_inputs(0.58),
        neutral_current_inputs=RSEWNeutralCurrentInputs(a_ref_c=0.60),
        include_charged_current=True,
    )
    charged = point.extras["rs_charged_current"]
    neutral = point.extras["rs_ew_couplings"]

    assert charged.a_ref_c == pytest.approx(0.60)
    assert charged.a_ref == neutral.a_ref


def test_exact_w_root_mass_matrix_and_eta_w_sign_are_pinned(sample_point):
    spectrum = sample_point.extras["rs_ew_spectrum"]
    charged = sample_point.extras["rs_charged_current"]

    assert spectrum.gauge_roots_x[0] == pytest.approx(
        GAUGE_ROOT_EPS_1E_MINUS_15,
        rel=0.0,
        abs=2.0e-14,
    )
    assert charged.m_w_gev == pytest.approx(77.50672575578955, rel=1.0e-12, abs=1.0e-12)
    assert charged.m_wprime_gev == pytest.approx(3071.5336376392315, rel=1.0e-12, abs=1.0e-10)
    assert 70.0 < charged.m_w_gev < 90.0
    assert spectrum.kk_ew_mass_gev < charged.m_wprime_gev < 1.05 * spectrum.kk_ew_mass_gev
    assert charged.w_diagonalization.zero_plus_gauge_masses_gev[0] == pytest.approx(0.0)
    assert charged.w_diagonalization.zero_plus_chi_ir[0] == pytest.approx(1.0)
    assert charged.eta_W == -1.0
    assert charged.w_diagonalization.metadata["eta_W_projection_value"] < 0.0
    assert charged.diagnostics["w_diagonalization_location"].endswith(
        "RSEWSpectrum.charged_w_diagonalization"
    )


def test_shared_a_ref_and_universal_c_recover_sm_after_one_gf_subtraction(universal_point):
    charged = universal_point.extras["rs_charged_current"]
    neutral = universal_point.extras["rs_ew_couplings"]

    assert charged.a_ref == neutral.a_ref
    assert charged.diagnostics["a_ref_source"] == "shared_a_ref_from_rs_ew_couplings"
    assert np.max(np.abs(charged.delta_g_W_ud_L)) == pytest.approx(0.0, abs=1.0e-18)
    assert np.max(np.abs(charged.delta_g_W_ud_R)) == pytest.approx(0.0, abs=0.0)
    assert np.max(np.abs(charged.delta_g_W_lnu_L)) == pytest.approx(0.0, abs=1.0e-18)
    assert np.max(np.abs(charged.delta_g_W_lnu_R)) == pytest.approx(0.0, abs=0.0)
    assert np.max(np.abs(charged.epsilon)) == pytest.approx(0.0, abs=1.0e-18)
    assert charged.delta_G_F_over_G_F == pytest.approx(2.898593923605594e-05, rel=1.0e-12)
    assert charged.charged_contact_LL[0, 0, 0, 0] / charged.C_SM_gev_minus2 == pytest.approx(
        charged.delta_G_F_over_G_F + 0.0j,
        rel=1.0e-12,
        abs=1.0e-18,
    )


def test_delta_gf_scale_sign_and_no_double_count_or_second_mkk(sample_point):
    spectrum = sample_point.extras["rs_ew_spectrum"]
    charged = sample_point.extras["rs_charged_current"]
    scale = (charged.m_w_gev / spectrum.kk_ew_mass_gev) ** 2

    assert scale == pytest.approx(0.0006674769485981234, rel=1.0e-12)
    assert charged.delta_G_F_over_G_F == pytest.approx(
        -1.8334382815196566e-05,
        rel=1.0e-12,
        abs=1.0e-16,
    )
    assert 1.0e-6 < abs(charged.delta_G_F_over_G_F) < scale
    assert charged.diagnostics["delta_GF_charged_contact"].real > 0.0
    assert charged.diagnostics["delta_GF_light_w_vertices"].real < 0.0

    i, j, a = 0, 1, 0
    unsubtracted = (
        charged.delta_g_W_ud_L[i, j] / charged.ckm[i, j]
        + charged.delta_g_W_lnu_L[a, a]
        + charged.charged_contact_LL[i, j, a, a] / charged.C_SM_gev_minus2
    )
    assert charged.epsilon[i, j, a] == pytest.approx(
        unsubtracted - charged.delta_G_F_over_G_F,
        rel=1.0e-12,
        abs=1.0e-16,
    )
    assert charged.epsilon[i, j, a] != pytest.approx(
        unsubtracted - 2.0 * charged.delta_G_F_over_G_F,
        rel=1.0e-6,
        abs=1.0e-16,
    )
    wrong_second_mkk = (
        charged.delta_g_W_ud_L[i, j] / charged.ckm[i, j]
        + charged.delta_g_W_lnu_L[a, a]
        + charged.charged_contact_LL[i, j, a, a]
        / (charged.C_SM_gev_minus2 * spectrum.kk_ew_mass_gev**2)
        - charged.delta_G_F_over_G_F
    )
    assert charged.epsilon[i, j, a] != pytest.approx(
        wrong_second_mkk,
        rel=1.0e-6,
        abs=1.0e-16,
    )


def test_shapes_finiteness_immutability_and_determinism(sample_point):
    repeat = _build_charged_point(_sample_fit(), c_l=0.58)
    charged = sample_point.extras["rs_charged_current"]
    repeat_charged = repeat.extras["rs_charged_current"]

    assert charged.delta_g_W_ud_L.shape == (3, 3)
    assert charged.delta_g_W_lnu_L.shape == (3, 3)
    assert charged.charged_contact_LL.shape == (3, 3, 3, 3)
    assert charged.epsilon.shape == (3, 3, 3)
    assert charged.delta_abs_vij_over_vij.shape == (3, 3, 3)
    for arr in (
        charged.delta_g_W_ud_L,
        charged.delta_g_W_lnu_L,
        charged.charged_contact_LL,
        charged.epsilon,
        charged.delta_abs_vij_over_vij,
    ):
        assert np.all(np.isfinite(arr))
        assert arr.flags.writeable is False
    assert charged.delta_g_W_ud_R_status == "minimal_rs_no_right_handed_charged_current"
    assert charged.delta_g_W_lnu_R_status == (
        "minimal_rs_no_right_handed_charged_lepton_neutrino_current"
    )
    with pytest.raises(ValueError):
        charged.epsilon[0, 0, 0] = 0.0
    assert charged.diagnostics["charged_contact_LL_actual"].flags.writeable is False
    with pytest.raises(ValueError):
        charged.diagnostics["charged_contact_LL_actual"][0, 0, 0, 0] = 0.0
    with pytest.raises(TypeError):
        charged.metadata["new"] = "blocked"
    with pytest.raises(FrozenInstanceError):
        charged.m_w_gev = 0.0

    assert np.array_equal(charged.delta_g_W_ud_L, repeat_charged.delta_g_W_ud_L)
    assert np.array_equal(charged.charged_contact_LL, repeat_charged.charged_contact_LL)
    assert np.array_equal(charged.epsilon, repeat_charged.epsilon)


def test_independent_recompute_one_delta_g_and_epsilon(sample_point):
    spectrum = sample_point.extras["rs_ew_spectrum"]
    lepton = sample_point.extras["lepton_mass_basis_couplings"]
    charged = sample_point.extras["rs_charged_current"]
    fit = _sample_fit()
    i, j, a = 0, 1, 0

    a_q = np.array(
        [
            spectrum.a(
                float(c),
                a_ref=charged.a_ref,
                rel_tol=OVERLAP_REL_TOL,
                min_modes=MIN_OVERLAP_MODES,
                max_modes=MAX_OVERLAP_MODES,
            )
            for c in fit.bulk_state.c_Q
        ],
        dtype=float,
    )
    a_l = np.array(
        [
            spectrum.a(
                float(c),
                a_ref=charged.a_ref,
                rel_tol=OVERLAP_REL_TOL,
                min_modes=MIN_OVERLAP_MODES,
                max_modes=MAX_OVERLAP_MODES,
            )
            for c in lepton.c_L
        ],
        dtype=float,
    )
    scale = (charged.m_w_gev / spectrum.kk_ew_mass_gev) ** 2
    manual_delta_ud = (
        charged.eta_W
        * scale
        * (fit.U_L_u.conjugate().T @ np.diag(a_q) @ fit.U_L_d)[i, j]
    )
    manual_delta_l = _hermitian(lepton.U_e_L.conjugate().T @ np.diag(a_l) @ lepton.U_e_L)
    manual_delta_l *= charged.eta_W * scale
    manual_contact, manual_mu_contact = _manual_charged_contact(
        spectrum,
        fit,
        lepton,
        charged.ckm,
    )
    manual_delta_gf = (
        manual_delta_l[0, 0]
        + manual_delta_l[1, 1]
        + manual_mu_contact / charged.C_SM_gev_minus2
    ).real
    manual_epsilon = (
        manual_delta_ud / charged.ckm[i, j]
        + manual_delta_l[a, a]
        + manual_contact[i, j, a, a] / charged.C_SM_gev_minus2
        - manual_delta_gf
    )

    assert charged.delta_g_W_ud_L[i, j] == pytest.approx(
        -8.169915813726022e-07 + 3.4880867617547286e-06j,
        rel=1.0e-12,
        abs=1.0e-16,
    )
    assert charged.epsilon[i, j, a] == pytest.approx(
        2.6427715369627687e-05 - 1.2332357024682375e-05j,
        rel=1.0e-12,
        abs=1.0e-16,
    )
    assert charged.delta_g_W_ud_L[i, j] == pytest.approx(
        manual_delta_ud,
        rel=1.0e-12,
        abs=1.0e-16,
    )
    assert charged.delta_G_F_over_G_F == pytest.approx(manual_delta_gf, rel=1.0e-12)
    assert charged.epsilon[i, j, a] == pytest.approx(
        manual_epsilon,
        rel=1.0e-12,
        abs=1.0e-16,
    )


def _manual_charged_contact(spectrum, fit, lepton, ckm):
    omega_q = np.stack(
        [spectrum.omega(float(c), max_modes=MAX_OVERLAP_MODES) for c in fit.bulk_state.c_Q],
        axis=1,
    )
    omega_l = np.stack(
        [spectrum.omega(float(c), max_modes=MAX_OVERLAP_MODES) for c in lepton.c_L],
        axis=1,
    )
    contact = np.zeros((3, 3, 3, 3), dtype=np.complex128)
    mu_contact = 0.0j
    g2 = math.sqrt(4.0 * math.pi * (1.0 / 127.952) / 0.23122)
    for mode_index in range(MAX_OVERLAP_MODES):
        prefactor = 0.5 * g2**2 / spectrum.gauge_masses_gev[mode_index] ** 2
        q_mode = fit.U_L_u.conjugate().T @ np.diag(omega_q[mode_index]) @ fit.U_L_d
        q_ratio = q_mode / ckm
        l_mode = _hermitian(
            lepton.U_e_L.conjugate().T @ np.diag(omega_l[mode_index]) @ lepton.U_e_L
        )
        contact += prefactor * q_ratio[:, :, None, None] * l_mode[None, None, :, :]
        mu_contact += prefactor * l_mode[0, 0] * l_mode[1, 1]
    return contact, mu_contact


def _hermitian(matrix):
    arr = np.asarray(matrix, dtype=np.complex128)
    return 0.5 * (arr + arr.conjugate().T)
