import math
from dataclasses import replace

import numpy as np
import pytest

from flavor_catalog_constraints import point_builder
from quarkConstraints import rare_kaon_dilepton, rare_kaon_snd
from quarkConstraints.rs_ew_couplings import (
    DEFAULT_A_REF_C,
    DEFAULT_S_Z,
    build_rs_ew_couplings,
)
from quarkConstraints.rs_semileptonic_wilsons import build_rs_semileptonic_wilsons
from tests.rs_ew_phase3b_helpers import (
    MAX_OVERLAP_MODES,
    MIN_OVERLAP_MODES,
    N_GAUGE_MODES,
    OVERLAP_REL_TOL,
    QUADRATURE_ORDER,
    _sample_fit,
    _scales_for_mkk,
)

SIN2_THETA_W = 0.23122
M_Z_GEV = 91.1876
ALPHA_EM_MZ = 1.0 / 127.952
GF_GEV_MINUS2 = 1.1663787e-5


def _lepton_inputs(c_l: float, *, c_e=None, alpha: float = 0.0, beta: float = 0.0):
    return {
        "c_L": float(c_l),
        "c_E": [float(c_l), float(c_l), float(c_l)] if c_e is None else c_e,
        "c_N": 0.27,
        "M_N": 1.22e18,
        "lightest_nu_mass": 0.002,
        "ordering": "normal",
        "majorana_alpha": float(alpha),
        "majorana_beta": float(beta),
    }


def _build_point(*, lepton_inputs=None, mkk_gev: float = 3000.0):
    lambda_ir, k = _scales_for_mkk(mkk_gev)
    return point_builder.build_from_rs_ew_inputs(
        _sample_fit(),
        Lambda_IR=lambda_ir,
        k=k,
        n_gauge_modes=N_GAUGE_MODES,
        quadrature_order=QUADRATURE_ORDER,
        min_overlap_modes=MIN_OVERLAP_MODES,
        max_overlap_modes=MAX_OVERLAP_MODES,
        overlap_rel_tol=OVERLAP_REL_TOL,
        lepton_sweep_inputs=lepton_inputs,
    )


def _rot12(theta: float) -> np.ndarray:
    c = math.cos(theta)
    s = math.sin(theta)
    return np.array([[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]], dtype=np.complex128)


def _rotated_nonuniversal_lfv_toy(*, mkk_gev: float):
    point = _build_point(lepton_inputs=_lepton_inputs(0.58), mkk_gev=mkk_gev)
    lepton = replace(
        point.extras["lepton_mass_basis_couplings"],
        c_L=np.array([0.55, 0.70, 0.58], dtype=float),
        U_e_L=_rot12(0.3),
    )
    couplings = build_rs_ew_couplings(
        _sample_fit(),
        spectrum=point.extras["rs_ew_spectrum"],
        lepton_mass_basis_couplings=lepton,
        overlap_rel_tol=OVERLAP_REL_TOL,
        min_overlap_modes=MIN_OVERLAP_MODES,
        max_overlap_modes=MAX_OVERLAP_MODES,
    )
    return couplings, build_rs_semileptonic_wilsons(couplings)


def _offdiag(matrix: np.ndarray) -> np.ndarray:
    arr = np.array(matrix, dtype=np.complex128, copy=True)
    np.fill_diagonal(arr, 0.0)
    return arr


def _max_charged_lfv_contact(couplings) -> float:
    max_value = 0.0
    for key, tensor in couplings.neutral_contacts.items():
        if key.endswith("_nunu"):
            continue
        for a in range(3):
            for b in range(3):
                if a != b:
                    max_value = max(max_value, float(np.max(np.abs(tensor[:, :, a, b]))))
    return max_value


def _max_lfv_wilson(bundle) -> float:
    return max(
        abs(value)
        for block in bundle.lfv_llqq.values()
        for coeff in block.values()
        for value in coeff.wilsons.values()
    )


def test_lepton_reference_c_zeroes_z_lfv_contacts_and_lfv_wilsons():
    point = _build_point(lepton_inputs=_lepton_inputs(DEFAULT_A_REF_C))
    couplings = point.extras["rs_ew_couplings"]
    wilsons = point.extras["rs_semileptonic_wilsons"]

    assert np.max(np.abs(couplings.z_delta_g_L_e)) == pytest.approx(0.0, abs=1.0e-18)
    assert np.max(np.abs(couplings.z_delta_g_R_e)) == pytest.approx(0.0, abs=1.0e-18)
    assert np.max(np.abs(couplings.z_delta_g_L_nu)) == pytest.approx(0.0, abs=1.0e-18)
    assert _max_charged_lfv_contact(couplings) == pytest.approx(0.0, abs=1.0e-22)
    assert _max_lfv_wilson(wilsons) == pytest.approx(0.0, abs=1.0e-14)


def test_family_universal_leptons_are_diagonal_nonzero_but_lfv_zero():
    point = _build_point(lepton_inputs=_lepton_inputs(0.58))
    couplings = point.extras["rs_ew_couplings"]
    wilsons = point.extras["rs_semileptonic_wilsons"]

    assert couplings.z_delta_g_L_e[0, 0] == pytest.approx(
        8.5400871097628e-06 + 0.0j,
        rel=1.0e-11,
        abs=1.0e-16,
    )
    assert couplings.z_delta_g_L_nu[0, 0] == pytest.approx(
        -1.588676075184686e-05 + 0.0j,
        rel=1.0e-11,
        abs=1.0e-16,
    )
    assert np.max(np.abs(_offdiag(couplings.z_delta_g_L_e))) == pytest.approx(0.0, abs=1.0e-20)
    assert np.max(np.abs(_offdiag(couplings.z_delta_g_R_e))) == pytest.approx(0.0, abs=1.0e-20)
    assert np.max(np.abs(_offdiag(couplings.z_delta_g_L_nu))) < 2.0e-20
    assert _max_charged_lfv_contact(couplings) == pytest.approx(0.0, abs=1.0e-22)
    assert _max_lfv_wilson(wilsons) == pytest.approx(0.0, abs=1.0e-14)


def test_rotated_nonuniversal_leptons_generate_lfv_wilsons_and_mkk_scaling():
    couplings, wilsons = _rotated_nonuniversal_lfv_toy(mkk_gev=3000.0)
    high_mkk_couplings, _ = _rotated_nonuniversal_lfv_toy(mkk_gev=6000.0)

    z_le_emu = couplings.z_delta_g_L_e[0, 1]
    lfv_sd_emu = wilsons.lfv_llqq["s_to_d"]["e_mu"]
    assert abs(z_le_emu) > 1.0e-8
    assert z_le_emu == pytest.approx(1.212164466423143e-05 + 0.0j, rel=1.0e-10, abs=1.0e-14)
    assert abs(lfv_sd_emu.c9_lfv_np) > 1.0e-7
    assert abs(lfv_sd_emu.c10_lfv_np) > 1.0e-7
    assert lfv_sd_emu.c9_lfv_np == pytest.approx(
        0.00023656625045856276 + 6.666935783109433e-05j,
        rel=1.0e-10,
        abs=1.0e-14,
    )
    assert lfv_sd_emu.c10_lfv_np == pytest.approx(
        -0.00023656625045856276 - 6.666935783109433e-05j,
        rel=1.0e-10,
        abs=1.0e-14,
    )

    profile_ratio = (
        abs(high_mkk_couplings.a_mass_basis["L_e"][0, 1]) ** 2
        / abs(couplings.a_mass_basis["L_e"][0, 1]) ** 2
    )
    rate_ratio = (
        abs(high_mkk_couplings.z_delta_g_L_e[0, 1]) ** 2
        / abs(couplings.z_delta_g_L_e[0, 1]) ** 2
    )
    expected_ratio = profile_ratio * (3000.0 / 6000.0) ** 4
    assert rate_ratio == pytest.approx(expected_ratio, rel=1.0e-12, abs=1.0e-15)
    assert rate_ratio == pytest.approx(0.0625, rel=1.0e-12, abs=1.0e-15)


def test_contacts_units_no_double_count_no_extra_gz_or_second_mkk_and_manual_wilsons():
    point = _build_point(lepton_inputs=_lepton_inputs(0.58))
    spectrum = point.extras["rs_ew_spectrum"]
    couplings = point.extras["rs_ew_couplings"]
    wilsons = point.extras["rs_semileptonic_wilsons"]
    fit = _sample_fit()

    assert couplings.contact_units == "GeV^-2"
    assert wilsons.contact_units == "GeV^-2"
    for matrix in (couplings.z_delta_g_L_e, couplings.z_delta_g_R_e, couplings.z_delta_g_L_nu):
        assert np.allclose(matrix, matrix.conjugate().T, rtol=0.0, atol=5.0e-13)

    manual_z_le = _manual_z_delta(
        spectrum,
        np.eye(3, dtype=np.complex128),
        np.array([0.58, 0.58, 0.58]),
        "e",
        "L",
    )
    assert couplings.z_delta_g_L_e[0, 0] == pytest.approx(manual_z_le[0, 0], rel=1.0e-12)
    assert couplings.z_delta_g_L_e[0, 0] != pytest.approx(
        couplings.g_z * manual_z_le[0, 0],
        rel=1.0e-6,
        abs=1.0e-16,
    )

    z_l_d = _manual_z_delta(spectrum, fit.U_L_d, fit.bulk_state.c_Q, "d", "L")
    z_l_e = manual_z_le
    z_l_nu = _manual_z_delta(
        spectrum,
        point.extras["lepton_mass_basis_couplings"].U_nu_L,
        np.array([0.58, 0.58, 0.58]),
        "nu",
        "L",
    )
    manual_lfv_contact = _manual_contact(
        z_l_d,
        z_l_e,
        quark_sector="d",
        quark_chirality="L",
        lepton_species="e",
        lepton_chirality="L",
        i=0,
        j=1,
        a=0,
        b=1,
    )
    manual_lfv_c9 = _manual_c9(manual_lfv_contact, 0.0j, transition="s_d")
    assert manual_lfv_c9 == pytest.approx(0.0j, abs=1.0e-18)
    assert wilsons.lfv_llqq["s_to_d"]["e_mu"].c9_lfv_np == pytest.approx(
        manual_lfv_c9,
        abs=1.0e-18,
    )

    manual_nunu_contact = _manual_contact(
        z_l_d,
        z_l_nu,
        quark_sector="d",
        quark_chirality="L",
        lepton_species="nu",
        lepton_chirality="L",
        i=0,
        j=1,
        a=0,
        b=0,
    )
    manual_x_np = manual_nunu_contact / rare_kaon_snd.g_sm_squared()
    assert wilsons.s_to_d_nunu.x_np_left[0, 0] == pytest.approx(
        0.000820388142621831 - 0.00011326866594393457j,
        rel=1.0e-12,
        abs=1.0e-15,
    )
    assert wilsons.s_to_d_nunu.x_np_left[0, 0] == pytest.approx(
        manual_x_np,
        rel=1.0e-12,
        abs=1.0e-15,
    )
    wrong_second_mkk = manual_nunu_contact / (
        rare_kaon_snd.g_sm_squared() * couplings.kk_ew_mass_gev**2
    )
    assert wilsons.s_to_d_nunu.x_np_left[0, 0] != pytest.approx(
        wrong_second_mkk,
        rel=1.0e-6,
        abs=1.0e-18,
    )


def test_completed_contact_reduces_to_phase3a_when_lepton_delta_g_is_zero():
    quark_only = _build_point()
    lepton_reference = _build_point(lepton_inputs=_lepton_inputs(DEFAULT_A_REF_C))
    old_wilson = quark_only.extras["rs_semileptonic_wilsons"].b_to_s_ll["mu"]
    new_wilson = lepton_reference.extras["rs_semileptonic_wilsons"].b_to_s_ll["mu"]
    old_contact = quark_only.extras["rs_ew_couplings"].contact("d", "L", "L", 1, 2, 1, 1)
    new_contact = lepton_reference.extras["rs_ew_couplings"].contact(
        "d",
        "L",
        "L",
        1,
        2,
        1,
        1,
    )

    assert new_contact == pytest.approx(old_contact, rel=0.0, abs=1.0e-22)
    assert new_wilson.c9_np == pytest.approx(
        0.3046272275031843 + 0.005803160503768459j,
        rel=1.0e-12,
        abs=1.0e-14,
    )
    assert new_wilson.c9_np == pytest.approx(old_wilson.c9_np, rel=0.0, abs=1.0e-14)
    assert new_wilson.c10_np == pytest.approx(old_wilson.c10_np, rel=0.0, abs=1.0e-14)


def test_nunu_majorana_phases_cancel_for_universal_active_cL():
    dirac = _build_point(lepton_inputs=_lepton_inputs(0.58, alpha=0.0, beta=0.0))
    majorana = _build_point(lepton_inputs=_lepton_inputs(0.58, alpha=1.1, beta=-0.7))
    c_dirac = dirac.extras["rs_ew_couplings"]
    c_majorana = majorana.extras["rs_ew_couplings"]
    w_dirac = dirac.extras["rs_semileptonic_wilsons"]
    w_majorana = majorana.extras["rs_semileptonic_wilsons"]

    assert np.max(np.abs(c_dirac.z_delta_g_L_nu - c_majorana.z_delta_g_L_nu)) < 1.0e-18
    assert np.max(np.abs(w_dirac.b_to_s_nunu.x_np_left - w_majorana.b_to_s_nunu.x_np_left)) < 1.0e-18
    assert np.max(np.abs(w_dirac.s_to_d_nunu.x_np_left - w_majorana.s_to_d_nunu.x_np_left)) < 1.0e-18


def test_determinism_finiteness_and_immutability_of_4a_extras():
    point = _build_point(lepton_inputs=_lepton_inputs(0.58))
    repeat = _build_point(lepton_inputs=_lepton_inputs(0.58))
    lepton = point.extras["lepton_mass_basis_couplings"]
    couplings = point.extras["rs_ew_couplings"]
    wilsons = point.extras["rs_semileptonic_wilsons"]

    assert np.array_equal(lepton.Y_N_bar_matrix, repeat.extras["lepton_mass_basis_couplings"].Y_N_bar_matrix)
    assert np.array_equal(couplings.z_delta_g_L_nu, repeat.extras["rs_ew_couplings"].z_delta_g_L_nu)
    assert np.array_equal(wilsons.b_to_s_nunu.x_np_left, repeat.extras["rs_semileptonic_wilsons"].b_to_s_nunu.x_np_left)
    for arr in (
        lepton.Y_E_bar_matrix,
        lepton.Y_N_bar_matrix,
        lepton.lfv_dipole_spurion,
        couplings.z_delta_g_L_e,
        couplings.z_delta_g_R_e,
        couplings.z_delta_g_L_nu,
        wilsons.b_to_s_nunu.x_np_left,
    ):
        assert np.all(np.isfinite(arr))
    with pytest.raises(ValueError):
        lepton.Y_N_bar_matrix[0, 0] = 0.0
    with pytest.raises(TypeError):
        lepton.params["c_L"] = 0.1
    for tensor in couplings.neutral_contacts.values():
        assert tensor.shape == (3, 3, 3, 3)
        assert tensor.flags.writeable is False
    with pytest.raises(TypeError):
        wilsons.lfv_llqq["s_to_d"] = {}


def _manual_z_delta(spectrum, rotation, c_values, species: str, chirality: str):
    a_ref = spectrum.a(DEFAULT_A_REF_C, rel_tol=OVERLAP_REL_TOL, min_modes=MIN_OVERLAP_MODES, max_modes=MAX_OVERLAP_MODES)
    a_values = np.array(
        [
            spectrum.a(
                float(c),
                a_ref=a_ref,
                rel_tol=OVERLAP_REL_TOL,
                min_modes=MIN_OVERLAP_MODES,
                max_modes=MAX_OVERLAP_MODES,
            )
            for c in c_values
        ],
        dtype=float,
    )
    a_mass = rotation.conjugate().T @ np.diag(a_values) @ rotation
    scale = (M_Z_GEV / spectrum.kk_ew_mass_gev) ** 2
    return DEFAULT_S_Z * _sm_z(species, chirality) * scale * a_mass


def _manual_contact(
    z_delta_q,
    z_delta_l,
    *,
    quark_sector: str,
    quark_chirality: str,
    lepton_species: str,
    lepton_chirality: str,
    i: int,
    j: int,
    a: int,
    b: int,
) -> complex:
    g_z = math.sqrt(4.0 * math.pi * ALPHA_EM_MZ / (SIN2_THETA_W * (1.0 - SIN2_THETA_W)))
    g_q = _sm_z(quark_sector, quark_chirality)
    g_l = _sm_z(lepton_species, lepton_chirality)
    delta_ij = 1.0 if i == j else 0.0
    delta_ab = 1.0 if a == b else 0.0
    return complex(
        g_z**2
        / M_Z_GEV**2
        * (
            (g_q * delta_ij + z_delta_q[i, j])
            * (g_l * delta_ab + z_delta_l[a, b])
            - g_q * g_l * delta_ij * delta_ab
        )
    )


def _manual_c9(contact_ll: complex, contact_lr: complex, *, transition: str) -> complex:
    lambda_ckm = rare_kaon_dilepton.ckm_factors().lambda_t
    if transition != "s_d":
        raise AssertionError("test helper only supports s_d")
    prefactor = -math.pi / (math.sqrt(2.0) * GF_GEV_MINUS2 * ALPHA_EM_MZ * lambda_ckm)
    return prefactor * (contact_ll + contact_lr)


def _sm_z(species: str, chirality: str) -> float:
    if species == "u":
        charge = 2.0 / 3.0
        t3 = 0.5 if chirality == "L" else 0.0
    elif species == "d":
        charge = -1.0 / 3.0
        t3 = -0.5 if chirality == "L" else 0.0
    elif species == "e":
        charge = -1.0
        t3 = -0.5 if chirality == "L" else 0.0
    elif species == "nu":
        charge = 0.0
        t3 = 0.5
    else:
        raise AssertionError("unsupported species")
    return float(t3 - charge * SIN2_THETA_W)
