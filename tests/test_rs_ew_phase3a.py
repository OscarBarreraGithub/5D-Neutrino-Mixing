import math
from dataclasses import dataclass

import numpy as np
import pytest

from flavor_catalog_constraints import point_builder
from quarkConstraints import rare_b_dilepton, rare_charm_dilepton, rare_kaon_dilepton
from quarkConstraints.rs_ew_couplings import DEFAULT_A_REF_C, DEFAULT_S_Z


GAUGE_ROOT_EPS_1E_MINUS_15 = 2.450509663813736
EPSILON_RS = 1.0e-15
TARGET_MKK_GEV = 3000.0
N_GAUGE_MODES = 64
QUADRATURE_ORDER = 512
MIN_OVERLAP_MODES = 16
MAX_OVERLAP_MODES = 64
OVERLAP_REL_TOL = 1.0e-3
ALPHA_EM_MZ = 1.0 / 127.952
SIN2_THETA_W = 0.23122
M_Z_GEV = 91.1876
GF_GEV_MINUS2 = 1.1663787e-5


@dataclass(frozen=True)
class _BulkState:
    c_Q: np.ndarray
    c_u: np.ndarray
    c_d: np.ndarray


@dataclass(frozen=True)
class _QuarkFit:
    bulk_state: _BulkState
    U_L_u: np.ndarray
    U_L_d: np.ndarray
    U_R_u: np.ndarray
    U_R_d: np.ndarray


def _scales_for_mkk(mkk_gev: float) -> tuple[float, float]:
    lambda_ir = float(mkk_gev) / GAUGE_ROOT_EPS_1E_MINUS_15
    return lambda_ir, lambda_ir / EPSILON_RS


def _rot12(theta: float) -> np.ndarray:
    c = math.cos(theta)
    s = math.sin(theta)
    return np.array([[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]], dtype=np.complex128)


def _rot23(theta: float) -> np.ndarray:
    c = math.cos(theta)
    s = math.sin(theta)
    return np.array([[1.0, 0.0, 0.0], [0.0, c, s], [0.0, -s, c]], dtype=np.complex128)


def _rot13(theta: float, phase: float) -> np.ndarray:
    c = math.cos(theta)
    s = math.sin(theta)
    e_ip = complex(math.cos(phase), math.sin(phase))
    return np.array(
        [[c, 0.0, s * np.conjugate(e_ip)], [0.0, 1.0, 0.0], [-s * e_ip, 0.0, c]],
        dtype=np.complex128,
    )


def _sample_fit() -> _QuarkFit:
    return _QuarkFit(
        bulk_state=_BulkState(
            c_Q=np.array([0.64, 0.56, 0.43], dtype=float),
            c_u=np.array([0.62, 0.34, 0.18], dtype=float),
            c_d=np.array([0.66, 0.57, 0.20], dtype=float),
        ),
        U_L_u=_rot12(0.11) @ _rot23(-0.07) @ _rot13(0.025, 0.4),
        U_L_d=_rot12(-0.19) @ _rot23(0.08) @ _rot13(-0.015, -0.2),
        U_R_u=_rot12(-0.09) @ _rot23(0.045),
        U_R_d=_rot12(0.05) @ _rot23(-0.12),
    )


def _identity_fit(c_q: np.ndarray, c_u: np.ndarray, c_d: np.ndarray) -> _QuarkFit:
    identity = np.eye(3, dtype=np.complex128)
    return _QuarkFit(
        bulk_state=_BulkState(c_Q=np.asarray(c_q), c_u=np.asarray(c_u), c_d=np.asarray(c_d)),
        U_L_u=identity,
        U_L_d=identity,
        U_R_u=identity,
        U_R_d=identity,
    )


def _build_point(fit: _QuarkFit, *, mkk_gev: float = TARGET_MKK_GEV):
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
    )


@pytest.fixture(scope="module")
def sample_fit() -> _QuarkFit:
    return _sample_fit()


@pytest.fixture(scope="module")
def sample_point(sample_fit):
    return _build_point(sample_fit)


def _max_abs(values) -> float:
    return float(max(np.max(np.abs(np.asarray(value))) for value in values))


def _all_quark_z_matrices(couplings):
    return [
        couplings.z_delta_g_L_u,
        couplings.z_delta_g_R_u,
        couplings.z_delta_g_L_d,
        couplings.z_delta_g_R_d,
        couplings.z_total_g_L_u,
        couplings.z_total_g_R_u,
        couplings.z_total_g_L_d,
        couplings.z_total_g_R_d,
    ]


def test_schema_fail_loud_and_declared_phase3a_keys_are_accepted():
    with pytest.raises(KeyError):
        point_builder.make_point(rs_ew_couplingz=object())

    marker = object()
    point = point_builder.make_point(
        rs_ew_spectrum=marker,
        rs_ew_couplings=marker,
        rs_semileptonic_wilsons=marker,
    )
    assert point.get_extra("rs_ew_spectrum") is marker
    assert point.get_extra("rs_ew_couplings") is marker
    assert point.get_extra("rs_semileptonic_wilsons") is marker


def test_builder_outputs_physical_ew_mass_and_hermitian_z_matrices(sample_point):
    spectrum = sample_point.extras["rs_ew_spectrum"]
    couplings = sample_point.extras["rs_ew_couplings"]

    assert sample_point.extras["kk_ew_mass_gev"] == pytest.approx(
        spectrum.kk_ew_mass_gev, rel=0.0, abs=1.0e-12
    )
    assert spectrum.kk_ew_mass_gev == pytest.approx(3000.000000000015, rel=0.0, abs=1.0e-9)
    for matrix in _all_quark_z_matrices(couplings):
        assert np.allclose(matrix, matrix.conjugate().T, rtol=0.0, atol=5.0e-13)


def test_full_subtraction_sm_limit_has_zero_z_contacts_and_wilsons():
    ref = np.array([DEFAULT_A_REF_C, DEFAULT_A_REF_C, DEFAULT_A_REF_C], dtype=float)
    point = _build_point(_identity_fit(ref, ref, ref))
    couplings = point.extras["rs_ew_couplings"]
    wilsons = point.extras["rs_semileptonic_wilsons"]

    max_z_delta = _max_abs(
        [
            couplings.z_delta_g_L_u,
            couplings.z_delta_g_R_u,
            couplings.z_delta_g_L_d,
            couplings.z_delta_g_R_d,
        ]
    )
    max_contact = _max_abs(couplings.neutral_contacts.values())
    max_wilson = max(
        abs(value)
        for block in (
            wilsons.b_to_s_ll,
            wilsons.b_to_d_ll,
            wilsons.s_to_d_ll,
            wilsons.c_to_u_ll,
        )
        for coeffs in block.values()
        for value in coeffs.wilsons.values()
    )

    assert max_z_delta == pytest.approx(0.0, rel=0.0, abs=1.0e-18)
    assert max_contact == pytest.approx(0.0, rel=0.0, abs=1.0e-22)
    assert max_wilson == pytest.approx(0.0, rel=0.0, abs=1.0e-14)


def test_species_family_universal_profiles_remove_fcnc_but_allow_diagonal_shifts():
    fit = _QuarkFit(
        bulk_state=_BulkState(
            c_Q=np.array([0.60, 0.60, 0.60], dtype=float),
            c_u=np.array([0.35, 0.35, 0.35], dtype=float),
            c_d=np.array([0.20, 0.20, 0.20], dtype=float),
        ),
        U_L_u=_rot12(0.31) @ _rot23(-0.17) @ _rot13(0.08, 0.7),
        U_L_d=_rot12(-0.22) @ _rot23(0.16) @ _rot13(-0.04, -0.5),
        U_R_u=_rot12(-0.18) @ _rot23(0.09),
        U_R_d=_rot12(0.14) @ _rot23(-0.21),
    )
    point = _build_point(fit)
    couplings = point.extras["rs_ew_couplings"]
    wilsons = point.extras["rs_semileptonic_wilsons"]

    offdiag = []
    for matrix in (
        couplings.z_delta_g_L_u,
        couplings.z_delta_g_R_u,
        couplings.z_delta_g_L_d,
        couplings.z_delta_g_R_d,
    ):
        masked = np.array(matrix, copy=True)
        np.fill_diagonal(masked, 0.0)
        offdiag.append(masked)

    assert _max_abs(offdiag) < 2.0e-18
    assert abs(couplings.z_delta_g_R_d[2, 2]) > 1.0e-3
    assert abs(couplings.contact("d", "R", "L", 2, 2, 1, 1)) > 1.0e-8
    assert max(abs(coeffs.c9_np) for coeffs in wilsons.b_to_s_ll.values()) < 1.0e-12
    assert max(abs(coeffs.c9_np) for coeffs in wilsons.s_to_d_ll.values()) < 1.0e-12
    assert max(abs(coeffs.c9_np) for coeffs in wilsons.c_to_u_ll.values()) < 1.0e-12


def test_ir_localized_b_right_sign_magnitude_and_mkk_scaling():
    ref = np.array([DEFAULT_A_REF_C, DEFAULT_A_REF_C, DEFAULT_A_REF_C], dtype=float)
    c_d = np.array([DEFAULT_A_REF_C, DEFAULT_A_REF_C, 0.20], dtype=float)
    fit = _identity_fit(ref, ref, c_d)
    point_3tev = _build_point(fit, mkk_gev=3000.0)
    point_6tev = _build_point(fit, mkk_gev=6000.0)

    z_3tev = float(point_3tev.extras["rs_ew_couplings"].z_delta_g_R_d[2, 2].real)
    z_6tev = float(point_6tev.extras["rs_ew_couplings"].z_delta_g_R_d[2, 2].real)

    assert z_3tev == pytest.approx(-0.0016604227573826689, rel=1.0e-10, abs=1.0e-15)
    assert z_6tev == pytest.approx(-0.0004151056893456672, rel=1.0e-10, abs=1.0e-15)
    assert z_3tev < 0.0
    assert 5.0e-4 < abs(z_3tev) < 3.0e-3
    assert z_6tev / z_3tev == pytest.approx(0.25, rel=1.0e-12, abs=1.0e-15)


def test_contact_units_immutability_determinism_and_finiteness(sample_fit, sample_point):
    repeat_point = _build_point(sample_fit)
    couplings = sample_point.extras["rs_ew_couplings"]
    repeat_couplings = repeat_point.extras["rs_ew_couplings"]
    wilsons = sample_point.extras["rs_semileptonic_wilsons"]
    repeat_wilsons = repeat_point.extras["rs_semileptonic_wilsons"]

    assert couplings.contact_units == "GeV^-2"
    assert wilsons.contact_units == "GeV^-2"
    assert couplings.includes_heavy_neutral_lepton is False
    assert wilsons.includes_heavy_neutral_lepton is False
    with pytest.raises(TypeError):
        sample_point.extras["rs_ew_couplings"] = None

    for matrix in _all_quark_z_matrices(couplings):
        assert np.all(np.isfinite(matrix))
    for contact in couplings.neutral_contacts.values():
        assert np.all(np.isfinite(contact))
    for block in (
        wilsons.b_to_s_ll,
        wilsons.b_to_d_ll,
        wilsons.s_to_d_ll,
        wilsons.c_to_u_ll,
    ):
        for coeffs in block.values():
            assert all(math.isfinite(value.real) and math.isfinite(value.imag) for value in coeffs.wilsons.values())

    for name in ("z_delta_g_L_u", "z_delta_g_R_u", "z_delta_g_L_d", "z_delta_g_R_d"):
        assert np.array_equal(getattr(couplings, name), getattr(repeat_couplings, name))
    assert wilsons.b_to_s_ll["mu"].c9_np == repeat_wilsons.b_to_s_ll["mu"].c9_np
    assert wilsons.c_to_u_ll["mu"].c10_np == repeat_wilsons.c_to_u_ll["mu"].c10_np


def test_independent_manual_recompute_matches_representative_wilsons(sample_fit, sample_point):
    spectrum = sample_point.extras["rs_ew_spectrum"]
    wilsons = sample_point.extras["rs_semileptonic_wilsons"]

    manual_bsmu = _manual_wilsons(
        spectrum,
        sample_fit,
        quark_sector="d",
        final_index=1,
        initial_index=2,
        lepton_index=1,
        transition="b_s",
    )
    manual_sdee = _manual_wilsons(
        spectrum,
        sample_fit,
        quark_sector="d",
        final_index=0,
        initial_index=1,
        lepton_index=0,
        transition="s_d",
    )
    manual_cumu = _manual_wilsons(
        spectrum,
        sample_fit,
        quark_sector="u",
        final_index=0,
        initial_index=1,
        lepton_index=1,
        transition="c_u",
    )

    assert manual_bsmu["C9_NP"] == pytest.approx(
        -0.15231361375159214 - 0.0029015802518842296j,
        rel=1.0e-11,
        abs=1.0e-13,
    )
    assert manual_bsmu["C10_NP"] == pytest.approx(
        2.027604016927477 + 0.03862593519547694j,
        rel=1.0e-11,
        abs=1.0e-13,
    )
    assert manual_sdee["C9_NP"] == pytest.approx(
        0.3665108412822379 + 0.10329048365544855j,
        rel=1.0e-11,
        abs=1.0e-13,
    )
    assert manual_cumu["C10_NP"] == pytest.approx(
        10.220387251251232 + 15.692970083295808j,
        rel=1.0e-11,
        abs=1.0e-12,
    )

    _assert_wilson_match(manual_bsmu, wilsons.b_to_s_ll["mu"])
    _assert_wilson_match(manual_sdee, wilsons.s_to_d_ll["e"])
    _assert_wilson_match(manual_cumu, wilsons.c_to_u_ll["mu"])


def _assert_wilson_match(manual: dict[str, complex], coeffs) -> None:
    for key, value in manual.items():
        assert coeffs.wilsons[key] == pytest.approx(value, rel=1.0e-12, abs=1.0e-14)


def _manual_wilsons(
    spectrum,
    fit: _QuarkFit,
    *,
    quark_sector: str,
    final_index: int,
    initial_index: int,
    lepton_index: int,
    transition: str,
) -> dict[str, complex]:
    if quark_sector == "d":
        z_l = _manual_z_delta(spectrum, fit.U_L_d, fit.bulk_state.c_Q, "d", "L")
        z_r = _manual_z_delta(spectrum, fit.U_R_d, fit.bulk_state.c_d, "d", "R")
    elif quark_sector == "u":
        z_l = _manual_z_delta(spectrum, fit.U_L_u, fit.bulk_state.c_Q, "u", "L")
        z_r = _manual_z_delta(spectrum, fit.U_R_u, fit.bulk_state.c_u, "u", "R")
    else:
        raise AssertionError("unsupported quark sector")

    g_l_e = -0.5 + SIN2_THETA_W
    g_r_e = SIN2_THETA_W
    c_ll = _manual_contact(z_l, quark_sector, "L", g_l_e, final_index, initial_index, lepton_index)
    c_lr = _manual_contact(z_l, quark_sector, "L", g_r_e, final_index, initial_index, lepton_index)
    c_rl = _manual_contact(z_r, quark_sector, "R", g_l_e, final_index, initial_index, lepton_index)
    c_rr = _manual_contact(z_r, quark_sector, "R", g_r_e, final_index, initial_index, lepton_index)
    prefactor = math.pi / (
        2.0
        * math.sqrt(2.0)
        * GF_GEV_MINUS2
        * ALPHA_EM_MZ
        * _manual_lambda(transition)
    )
    return {
        "C9_NP": prefactor * (c_ll + c_lr),
        "C10_NP": prefactor * (c_lr - c_ll),
        "C9p_NP": prefactor * (c_rl + c_rr),
        "C10p_NP": prefactor * (c_rr - c_rl),
    }


def _manual_z_delta(
    spectrum,
    rotation: np.ndarray,
    c_values: np.ndarray,
    species: str,
    chirality: str,
) -> np.ndarray:
    a_ref = spectrum.a(
        DEFAULT_A_REF_C,
        rel_tol=OVERLAP_REL_TOL,
        min_modes=MIN_OVERLAP_MODES,
        max_modes=MAX_OVERLAP_MODES,
    )
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
    a_mass = 0.5 * (a_mass + a_mass.conjugate().T)
    return DEFAULT_S_Z * _manual_sm_chiral_z(species, chirality) * (M_Z_GEV / spectrum.kk_ew_mass_gev) ** 2 * a_mass


def _manual_contact(
    z_delta_q: np.ndarray,
    quark_sector: str,
    quark_chirality: str,
    lepton_coupling: float,
    i: int,
    j: int,
    a: int,
) -> complex:
    g_z = math.sqrt(4.0 * math.pi * ALPHA_EM_MZ / (SIN2_THETA_W * (1.0 - SIN2_THETA_W)))
    g_q = _manual_sm_chiral_z(quark_sector, quark_chirality)
    delta_ij = 1.0 if i == j else 0.0
    return complex(
        g_z**2
        / M_Z_GEV**2
        * (
            (g_q * delta_ij + z_delta_q[i, j]) * lepton_coupling
            - g_q * lepton_coupling * delta_ij
        )
    )


def _manual_sm_chiral_z(species: str, chirality: str) -> float:
    if species == "u":
        charge = 2.0 / 3.0
        t3 = 0.5 if chirality == "L" else 0.0
    elif species == "d":
        charge = -1.0 / 3.0
        t3 = -0.5 if chirality == "L" else 0.0
    else:
        raise AssertionError("unsupported species")
    return float(t3 - charge * SIN2_THETA_W)


def _manual_lambda(transition: str) -> complex:
    if transition == "b_s":
        return rare_b_dilepton.ckm_factors("b_s").lambda_t
    if transition == "s_d":
        return rare_kaon_dilepton.ckm_factors().lambda_t
    if transition == "c_u":
        return rare_charm_dilepton.ckm_factors("c_u").lambda_b
    raise AssertionError("unsupported transition")
