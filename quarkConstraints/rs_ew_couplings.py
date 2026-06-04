"""RS electroweak neutral-current couplings.

This module is the builder bridge from the numerical RS-EW spectrum to
mass-basis light-Z coupling shifts and semileptonic contacts.  It intentionally
does not call any rare-decay proxy helper.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Any, Mapping

import numpy as np

from .rs_ew_spectrum import (
    DEFAULT_MAX_TRUNCATION_MODES,
    DEFAULT_MIN_TRUNCATION_MODES,
    DEFAULT_OVERLAP_RTOL,
    RSEWSpectrum,
)


RS_EW_COUPLINGS_INPUT_BUNDLE_V1 = "quarkConstraints.rs_ew_couplings.inputs.v1"
RS_EW_COUPLINGS_MODEL_V1 = "RS_EW_NEUTRAL_CURRENT_PHASE4A_V1"
RS_EW_COUPLINGS_MATCHING_ASSUMPTION_V1 = (
    "minimal-RS light-Z quark, charged-lepton, and active-neutrino neutral "
    "currents with full light-Z product-minus-SM semileptonic contacts; "
    "heavy neutral Z'/gamma' exchange remains deferred"
)
RS_ZBB_FERMION_KK_MIXING_MODEL_V1 = "RS_ZBB_FERMION_KK_MIXING_PHASE6A_V1"
RS_ZBB_FERMION_KK_MIXING_MATCHING_ASSUMPTION_V1 = (
    "minimal non-custodial Casagrande ZMA fermion-KK admixture for diagonal "
    "Z b bbar only; no exact tower diagonalization, custodial top partner, "
    "exotic fermion, or BKT contribution is inferred"
)
RS_LEPTON_MASS_BASIS_INPUT_BUNDLE_V1 = (
    "quarkConstraints.rs_ew_couplings.lepton_inputs.v1"
)
RS_LEPTON_MASS_BASIS_MODEL_V1 = "RS_LEPTON_MASS_BASIS_PHASE4A_V1"
RS_LEPTON_MASS_BASIS_MATCHING_ASSUMPTION_V1 = (
    "charged-lepton Yukawa fit is diagonal, so U_e_L=U_e_R=I in the "
    "charged-lepton mass basis; PMNS is stored once for active-neutrino "
    "basis metadata and dipole spurions"
)

DEFAULT_A_REF_C = 0.65
DEFAULT_S_Z = -1.0
LEPTON_FLAVORS: tuple[str, str, str] = ("e", "mu", "tau")
NEUTRINO_FLAVORS: tuple[str, str, str] = ("nu1", "nu2", "nu3")


@dataclass(frozen=True)
class RSEWNeutralCurrentInputs:
    """Shared numerical inputs for Phase-4a neutral-current matching."""

    input_bundle: str = RS_EW_COUPLINGS_INPUT_BUNDLE_V1
    alpha_em_mz: float = 1.0 / 127.952
    sin2_theta_w: float = 0.23122
    m_z_gev: float = 91.1876
    s_z: float = DEFAULT_S_Z
    a_ref_c: float = DEFAULT_A_REF_C

    def __post_init__(self) -> None:
        for name in ("alpha_em_mz", "sin2_theta_w", "m_z_gev"):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be positive and finite")
        if not 0.0 < float(self.sin2_theta_w) < 1.0:
            raise ValueError("sin2_theta_w must lie between zero and one")
        if not math.isfinite(float(self.s_z)):
            raise ValueError("s_z must be finite")

    @property
    def g_z(self) -> float:
        return float(
            math.sqrt(
                4.0
                * math.pi
                * float(self.alpha_em_mz)
                / (float(self.sin2_theta_w) * (1.0 - float(self.sin2_theta_w)))
            )
        )


@dataclass(frozen=True)
class RSLeptonMassBasisCouplings:
    """Charged-lepton mass-basis inputs for Phase-4a light-Z matching."""

    model_label: str
    input_bundle: str
    matching_assumption: str
    kk_ew_mass_gev: float
    c_L: np.ndarray
    c_E: np.ndarray
    c_N: np.ndarray
    M_N: float
    f_L: np.ndarray
    f_E: np.ndarray
    f_N: np.ndarray
    f_N_UV: np.ndarray
    Y_E_bar_vector: np.ndarray
    Y_E_bar_matrix: np.ndarray
    Y_N_bar_vector: np.ndarray
    Y_N_matrix: np.ndarray
    Y_N_bar_matrix: np.ndarray
    pmns: np.ndarray
    U_e_L: np.ndarray
    U_e_R: np.ndarray
    U_nu_L: np.ndarray
    lfv_dipole_spurion: np.ndarray
    params: Mapping[str, Any]
    basis_metadata: Mapping[str, Any] = field(default_factory=dict)
    matching_status: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not math.isfinite(float(self.kk_ew_mass_gev)) or float(self.kk_ew_mass_gev) <= 0.0:
            raise ValueError("kk_ew_mass_gev must be positive and finite")
        if not math.isfinite(float(self.M_N)) or float(self.M_N) <= 0.0:
            raise ValueError("M_N must be positive and finite")
        for name in ("c_L", "c_E", "c_N", "f_L", "f_E", "f_N", "f_N_UV"):
            object.__setattr__(self, name, _readonly_real_triplet(getattr(self, name), name))
        for name in ("Y_E_bar_vector", "Y_N_bar_vector"):
            object.__setattr__(
                self,
                name,
                _readonly_complex_triplet(getattr(self, name), name),
            )
        for name in (
            "Y_E_bar_matrix",
            "Y_N_matrix",
            "Y_N_bar_matrix",
            "pmns",
            "U_e_L",
            "U_e_R",
            "U_nu_L",
            "lfv_dipole_spurion",
        ):
            object.__setattr__(self, name, _readonly_complex_matrix(getattr(self, name), name))
        for name in ("pmns", "U_e_L", "U_e_R", "U_nu_L"):
            matrix = getattr(self, name)
            identity = np.eye(3, dtype=np.complex128)
            if not np.allclose(matrix.conjugate().T @ matrix, identity, rtol=0.0, atol=1.0e-8):
                raise ValueError(f"{name} must be unitary")
        if not np.allclose(
            self.lfv_dipole_spurion,
            self.lfv_dipole_spurion.conjugate().T,
            rtol=0.0,
            atol=5.0e-13,
        ):
            raise ValueError("lfv_dipole_spurion must be Hermitian")
        object.__setattr__(self, "params", MappingProxyType(dict(self.params)))
        object.__setattr__(
            self,
            "basis_metadata",
            MappingProxyType(dict(self.basis_metadata)),
        )
        object.__setattr__(
            self,
            "matching_status",
            MappingProxyType(dict(self.matching_status)),
        )

    @property
    def Y_N_bar(self) -> np.ndarray:
        return self.Y_N_bar_vector

    @property
    def Y_E_bar(self) -> np.ndarray:
        return self.Y_E_bar_vector

    @property
    def y_n_bar(self) -> np.ndarray:
        return self.Y_N_bar_vector

    @property
    def pmns_matrix(self) -> np.ndarray:
        return self.pmns


@dataclass(frozen=True)
class RSZbbFermionKKMixing:
    """Minimal-RS Casagrande Zbb fermion-KK admixture result."""

    model_label: str
    matching_assumption: str
    lambda_ir_gev: float
    m_b_gev: float
    c_Q: np.ndarray
    c_d: np.ndarray
    F_Q: np.ndarray
    F_d: np.ndarray
    Y_d_bulk_basis: np.ndarray
    profile_B_Q: np.ndarray
    profile_B_d: np.ndarray
    yukawa_ratio_column_b: np.ndarray
    yukawa_ratio_row_b: np.ndarray
    B_Q: float
    B_d: float
    delta_g_L_b: float
    delta_g_R_b: float
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        lambda_ir = float(self.lambda_ir_gev)
        if not math.isfinite(lambda_ir) or lambda_ir <= 0.0:
            raise ValueError("lambda_ir_gev must be positive and finite")
        m_b = float(self.m_b_gev)
        if not math.isfinite(m_b) or m_b < 0.0:
            raise ValueError("m_b_gev must be non-negative and finite")
        object.__setattr__(self, "lambda_ir_gev", lambda_ir)
        object.__setattr__(self, "m_b_gev", m_b)

        for name in ("c_Q", "c_d", "F_Q", "F_d", "profile_B_Q", "profile_B_d"):
            object.__setattr__(self, name, _readonly_real_triplet(getattr(self, name), name))
        for name in ("yukawa_ratio_column_b", "yukawa_ratio_row_b"):
            arr = _readonly_real_triplet(getattr(self, name), name)
            if np.any(arr < 0.0):
                raise ValueError(f"{name} must be non-negative")
            object.__setattr__(self, name, arr)
        object.__setattr__(
            self,
            "Y_d_bulk_basis",
            _readonly_complex_matrix(self.Y_d_bulk_basis, "Y_d_bulk_basis"),
        )
        for name in ("B_Q", "B_d", "delta_g_L_b", "delta_g_R_b"):
            value = float(getattr(self, name))
            if not math.isfinite(value):
                raise ValueError(f"{name} must be finite")
            object.__setattr__(self, name, value)
        object.__setattr__(self, "metadata", MappingProxyType(dict(self.metadata)))


@dataclass(frozen=True)
class RSEWMassBasisCouplings:
    """Mass-basis light-Z shifts and neutral semileptonic contacts."""

    model_label: str
    input_bundle: str
    matching_assumption: str
    kk_ew_mass_gev: float
    m_z_gev: float
    alpha_em_mz: float
    sin2_theta_w: float
    g_z: float
    s_z: float
    a_ref: float
    a_ref_c: float
    contact_units: str
    includes_heavy_neutral_exchange: bool
    includes_heavy_neutral_lepton: bool
    z_delta_g_L_u: np.ndarray
    z_delta_g_R_u: np.ndarray
    z_delta_g_L_d: np.ndarray
    z_delta_g_R_d: np.ndarray
    z_delta_g_L_e: np.ndarray
    z_delta_g_R_e: np.ndarray
    z_delta_g_L_nu: np.ndarray
    z_total_g_L_u: np.ndarray
    z_total_g_R_u: np.ndarray
    z_total_g_L_d: np.ndarray
    z_total_g_R_d: np.ndarray
    z_total_g_L_e: np.ndarray
    z_total_g_R_e: np.ndarray
    z_total_g_L_nu: np.ndarray
    neutral_contacts: Mapping[str, np.ndarray]
    a_profile_values: Mapping[str, np.ndarray]
    a_mass_basis: Mapping[str, np.ndarray]
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.contact_units != "GeV^-2":
            raise ValueError("contact_units must be 'GeV^-2'")
        for name in ("kk_ew_mass_gev", "m_z_gev", "alpha_em_mz", "g_z"):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be positive and finite")
        if not 0.0 < float(self.sin2_theta_w) < 1.0:
            raise ValueError("sin2_theta_w must lie between zero and one")

        matrix_names = (
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
        for name in matrix_names:
            arr = _readonly_complex_matrix(getattr(self, name), name)
            if not np.allclose(arr, arr.conjugate().T, rtol=0.0, atol=5.0e-13):
                raise ValueError(f"{name} must be Hermitian")
            object.__setattr__(self, name, arr)

        contacts: dict[str, np.ndarray] = {}
        for key, value in self.neutral_contacts.items():
            arr = np.array(value, dtype=np.complex128, copy=True)
            if arr.shape != (3, 3, 3, 3):
                raise ValueError(f"neutral contact {key!r} must have shape (3, 3, 3, 3)")
            if not np.all(np.isfinite(arr)):
                raise ValueError(f"neutral contact {key!r} contains non-finite values")
            arr.setflags(write=False)
            contacts[str(key)] = arr
        object.__setattr__(self, "neutral_contacts", MappingProxyType(contacts))

        object.__setattr__(
            self,
            "a_profile_values",
            MappingProxyType(
                {
                    str(key): _readonly_real_triplet(value, f"a_profile_values[{key!r}]")
                    for key, value in self.a_profile_values.items()
                }
            ),
        )
        object.__setattr__(
            self,
            "a_mass_basis",
            MappingProxyType(
                {
                    str(key): _readonly_complex_matrix(value, f"a_mass_basis[{key!r}]")
                    for key, value in self.a_mass_basis.items()
                }
            ),
        )
        object.__setattr__(self, "metadata", MappingProxyType(dict(self.metadata)))

    def contact(
        self,
        quark_sector: str,
        quark_chirality: str,
        lepton_chirality: str,
        i: int,
        j: int,
        a: int,
        b: int,
    ) -> complex:
        """Return ``C_AB(q_i q_j l_a l_b)`` in ``GeV^-2``."""

        key = _contact_key(quark_sector, quark_chirality, lepton_chirality)
        return complex(self.neutral_contacts[key][int(i), int(j), int(a), int(b)])

    def nunu_contact(
        self,
        quark_sector: str,
        quark_chirality: str,
        i: int,
        j: int,
        a: int,
        b: int,
    ) -> complex:
        """Return ``C_AB(q_i q_j nu_a nu_b)`` in ``GeV^-2``."""

        key = _nunu_contact_key(quark_sector, quark_chirality)
        return complex(self.neutral_contacts[key][int(i), int(j), int(a), int(b)])


def build_rs_lepton_mass_basis_couplings(
    yukawa_result: Any,
    *,
    kk_ew_mass_gev: float,
    rotation_unitarity_tolerance: float = 1.0e-8,
) -> RSLeptonMassBasisCouplings:
    """Build the Phase-4a charged-lepton mass-basis typed extra."""

    from neutrinos.neutrinoValues import get_pmns

    params = dict(getattr(yukawa_result, "params"))
    k = _positive_float(params["k"], "k")
    c_l = _broadcast_real_triplet(params["c_L"], "c_L")
    c_e = _readonly_real_triplet(params["c_E"], "c_E")
    c_n = _broadcast_real_triplet(params["c_N"], "c_N")
    m_n = _positive_float(params["M_N"], "M_N")
    y_e_bar = _readonly_complex_triplet(getattr(yukawa_result, "Y_E_bar"), "Y_E_bar")
    y_n_bar = _readonly_complex_triplet(getattr(yukawa_result, "Y_N_bar"), "Y_N_bar")
    y_n_matrix = _readonly_complex_matrix(
        getattr(yukawa_result, "Y_N_matrix"),
        "Y_N_matrix",
    )
    ordering = str(params["ordering"])
    pmns = _readonly_complex_matrix(
        get_pmns(
            ordering,
            float(params["majorana_alpha"]),
            float(params["majorana_beta"]),
        ),
        "pmns",
    )
    y_n_bar_matrix_from_pmns = pmns @ np.diag(y_n_bar)
    y_n_bar_matrix = np.array(2.0 * k * y_n_matrix, dtype=np.complex128, copy=True)
    if not np.allclose(
        y_n_bar_matrix,
        y_n_bar_matrix_from_pmns,
        rtol=1.0e-9,
        atol=1.0e-12,
    ):
        raise ValueError("Y_N_bar_matrix must equal both 2*k*Y_N_matrix and PMNS@diag(Y_N_bar)")

    identity = np.eye(3, dtype=np.complex128)
    if not np.allclose(identity.conjugate().T @ identity, identity, rtol=0.0, atol=rotation_unitarity_tolerance):
        raise ValueError("identity charged-lepton rotations failed unitarity validation")
    params_with_mkk = dict(params)
    params_with_mkk["M_KK"] = float(kk_ew_mass_gev)
    lfv_spurion = y_n_bar_matrix @ y_n_bar_matrix.conjugate().T
    return RSLeptonMassBasisCouplings(
        model_label=RS_LEPTON_MASS_BASIS_MODEL_V1,
        input_bundle=RS_LEPTON_MASS_BASIS_INPUT_BUNDLE_V1,
        matching_assumption=RS_LEPTON_MASS_BASIS_MATCHING_ASSUMPTION_V1,
        kk_ew_mass_gev=float(kk_ew_mass_gev),
        c_L=c_l,
        c_E=c_e,
        c_N=c_n,
        M_N=m_n,
        f_L=_broadcast_real_triplet(getattr(yukawa_result, "f_L"), "f_L"),
        f_E=_readonly_real_triplet(getattr(yukawa_result, "f_E"), "f_E"),
        f_N=_broadcast_real_triplet(getattr(yukawa_result, "f_N"), "f_N"),
        f_N_UV=_broadcast_real_triplet(getattr(yukawa_result, "f_N_UV"), "f_N_UV"),
        Y_E_bar_vector=y_e_bar,
        Y_E_bar_matrix=np.diag(y_e_bar),
        Y_N_bar_vector=y_n_bar,
        Y_N_matrix=y_n_matrix,
        Y_N_bar_matrix=y_n_bar_matrix,
        pmns=pmns,
        U_e_L=identity,
        U_e_R=identity,
        U_nu_L=pmns,
        lfv_dipole_spurion=lfv_spurion,
        params=params_with_mkk,
        basis_metadata={
            "charged_lepton_basis": "charged-lepton mass basis",
            "charged_lepton_fit": "diagonal",
            "U_e_L_source": "identity_from_diagonal_charged_lepton_fit",
            "U_e_R_source": "identity_from_diagonal_charged_lepton_fit",
            "U_nu_L_source": "PMNS from neutrinos.neutrinoValues.get_pmns",
            "pmns_second_rotation_applied": False,
            "active_neutrino_current_basis": "neutrino mass basis for nu_a nu_b blocks",
        },
        matching_status={
            "tree_level_light_z_lepton_currents": "available",
            "loop_dipole_matching": "spurion_only",
            "heavy_neutral_exchange": "deferred",
            "uses_yukawa_result_params": True,
        },
    )


def build_rs_ew_couplings(
    quark_fit_result: Any,
    *,
    spectrum: RSEWSpectrum,
    lepton_mass_basis_couplings: RSLeptonMassBasisCouplings | None = None,
    inputs: RSEWNeutralCurrentInputs | None = None,
    include_fermion_kk_mixing: bool = False,
    overlap_rel_tol: float = DEFAULT_OVERLAP_RTOL,
    min_overlap_modes: int = DEFAULT_MIN_TRUNCATION_MODES,
    max_overlap_modes: int | None = None,
    model_label: str = "minimal_rs",
    rotation_unitarity_tolerance: float = 1.0e-8,
) -> RSEWMassBasisCouplings:
    """Build Phase-4a RS-EW mass-basis couplings from a quark fit result."""

    p = RSEWNeutralCurrentInputs() if inputs is None else inputs
    max_modes = _resolve_max_overlap_modes(spectrum, max_overlap_modes)
    min_modes = int(min_overlap_modes)
    if min_modes >= max_modes:
        raise ValueError("min_overlap_modes must be smaller than max_overlap_modes")

    bulk_state = getattr(quark_fit_result, "bulk_state")
    c_q = _real_triplet_from_attr(bulk_state, "c_Q")
    c_u = _real_triplet_from_attr(bulk_state, "c_u")
    c_d = _real_triplet_from_attr(bulk_state, "c_d")

    rotations = {
        "U_L_u": _unitary_matrix_from_attr(
            quark_fit_result, "U_L_u", tolerance=rotation_unitarity_tolerance
        ),
        "U_L_d": _unitary_matrix_from_attr(
            quark_fit_result, "U_L_d", tolerance=rotation_unitarity_tolerance
        ),
        "U_R_u": _unitary_matrix_from_attr(
            quark_fit_result, "U_R_u", tolerance=rotation_unitarity_tolerance
        ),
        "U_R_d": _unitary_matrix_from_attr(
            quark_fit_result, "U_R_d", tolerance=rotation_unitarity_tolerance
        ),
    }

    a_ref = float(
        spectrum.a(
            p.a_ref_c,
            rel_tol=overlap_rel_tol,
            min_modes=min_modes,
            max_modes=max_modes,
        )
    )
    a_profiles = {
        "c_Q": _a_triplet(
            spectrum,
            c_q,
            a_ref=a_ref,
            rel_tol=overlap_rel_tol,
            min_modes=min_modes,
            max_modes=max_modes,
        ),
        "c_u": _a_triplet(
            spectrum,
            c_u,
            a_ref=a_ref,
            rel_tol=overlap_rel_tol,
            min_modes=min_modes,
            max_modes=max_modes,
        ),
        "c_d": _a_triplet(
            spectrum,
            c_d,
            a_ref=a_ref,
            rel_tol=overlap_rel_tol,
            min_modes=min_modes,
            max_modes=max_modes,
        ),
    }

    a_mass = {
        "L_u": _mass_basis_profile(rotations["U_L_u"], a_profiles["c_Q"]),
        "L_d": _mass_basis_profile(rotations["U_L_d"], a_profiles["c_Q"]),
        "R_u": _mass_basis_profile(rotations["U_R_u"], a_profiles["c_u"]),
        "R_d": _mass_basis_profile(rotations["U_R_d"], a_profiles["c_d"]),
    }

    scale = (float(p.m_z_gev) / float(spectrum.kk_ew_mass_gev)) ** 2
    g_l_u = _sm_chiral_z_coupling("u", "L", p.sin2_theta_w)
    g_r_u = _sm_chiral_z_coupling("u", "R", p.sin2_theta_w)
    g_l_d = _sm_chiral_z_coupling("d", "L", p.sin2_theta_w)
    g_r_d = _sm_chiral_z_coupling("d", "R", p.sin2_theta_w)
    g_l_e = _sm_chiral_z_coupling("e", "L", p.sin2_theta_w)
    g_r_e = _sm_chiral_z_coupling("e", "R", p.sin2_theta_w)
    g_l_nu = _sm_chiral_z_coupling("nu", "L", p.sin2_theta_w)

    z_delta_l_u = _z_delta(g_l_u, a_mass["L_u"], scale=scale, s_z=p.s_z)
    z_delta_r_u = _z_delta(g_r_u, a_mass["R_u"], scale=scale, s_z=p.s_z)
    z_delta_l_d = _z_delta(g_l_d, a_mass["L_d"], scale=scale, s_z=p.s_z)
    z_delta_r_d = _z_delta(g_r_d, a_mass["R_d"], scale=scale, s_z=p.s_z)
    zbb_fermion_mixing = None
    if include_fermion_kk_mixing:
        zbb_fermion_mixing = build_rs_zbb_fermion_kk_mixing(
            quark_fit_result,
            spectrum=spectrum,
        )
        z_delta_l_d = np.array(z_delta_l_d, dtype=np.complex128, copy=True)
        z_delta_r_d = np.array(z_delta_r_d, dtype=np.complex128, copy=True)
        z_delta_l_d[2, 2] += complex(zbb_fermion_mixing.delta_g_L_b, 0.0)
        z_delta_r_d[2, 2] += complex(zbb_fermion_mixing.delta_g_R_b, 0.0)
        z_delta_l_d = _hermitian(z_delta_l_d)
        z_delta_r_d = _hermitian(z_delta_r_d)
    lepton_matching_status: Mapping[str, Any]
    if lepton_mass_basis_couplings is None:
        z_delta_l_e = np.zeros((3, 3), dtype=np.complex128)
        z_delta_r_e = np.zeros((3, 3), dtype=np.complex128)
        z_delta_l_nu = np.zeros((3, 3), dtype=np.complex128)
        lepton_matching_status = {
            "lepton_mass_basis_couplings_present": False,
            "charged_lepton_delta_g": "zero_deferred_or_not_supplied",
            "active_neutrino_delta_g": "zero_deferred_or_not_supplied",
        }
    else:
        lepton_profiles = {
            "c_L": _a_triplet(
                spectrum,
                lepton_mass_basis_couplings.c_L,
                a_ref=a_ref,
                rel_tol=overlap_rel_tol,
                min_modes=min_modes,
                max_modes=max_modes,
            ),
            "c_E": _a_triplet(
                spectrum,
                lepton_mass_basis_couplings.c_E,
                a_ref=a_ref,
                rel_tol=overlap_rel_tol,
                min_modes=min_modes,
                max_modes=max_modes,
            ),
        }
        a_profiles.update(lepton_profiles)
        a_mass.update(
            {
                "L_e": _mass_basis_profile(
                    lepton_mass_basis_couplings.U_e_L,
                    lepton_profiles["c_L"],
                ),
                "R_e": _mass_basis_profile(
                    lepton_mass_basis_couplings.U_e_R,
                    lepton_profiles["c_E"],
                ),
                "L_nu": _mass_basis_profile(
                    lepton_mass_basis_couplings.U_nu_L,
                    lepton_profiles["c_L"],
                ),
            }
        )
        z_delta_l_e = _z_delta(g_l_e, a_mass["L_e"], scale=scale, s_z=p.s_z)
        z_delta_r_e = _z_delta(g_r_e, a_mass["R_e"], scale=scale, s_z=p.s_z)
        z_delta_l_nu = _z_delta(g_l_nu, a_mass["L_nu"], scale=scale, s_z=p.s_z)
        lepton_matching_status = {
            "lepton_mass_basis_couplings_present": True,
            "charged_lepton_delta_g": "included_with_U_e_identity",
            "active_neutrino_delta_g": "included_with_PMNS_basis_metadata",
            "singlet_c_N_used_for_z_current": False,
            "pmns_second_rotation_applied": False,
        }

    contacts = _neutral_contacts(
        g_z=p.g_z,
        m_z_gev=p.m_z_gev,
        sin2_theta_w=p.sin2_theta_w,
        z_delta_g_L_u=z_delta_l_u,
        z_delta_g_R_u=z_delta_r_u,
        z_delta_g_L_d=z_delta_l_d,
        z_delta_g_R_d=z_delta_r_d,
        z_delta_g_L_e=z_delta_l_e,
        z_delta_g_R_e=z_delta_r_e,
        z_delta_g_L_nu=z_delta_l_nu,
    )

    identity = np.eye(3, dtype=np.complex128)
    matching_assumption = RS_EW_COUPLINGS_MATCHING_ASSUMPTION_V1
    if zbb_fermion_mixing is not None:
        matching_assumption = (
            f"{matching_assumption}; "
            "diagonal Zbb minimal non-custodial Casagrande fermion-KK "
            "admixture included using Lambda_IR=spectrum.lambda_ir_gev"
        )
    return RSEWMassBasisCouplings(
        model_label=RS_EW_COUPLINGS_MODEL_V1,
        input_bundle=p.input_bundle,
        matching_assumption=matching_assumption,
        kk_ew_mass_gev=float(spectrum.kk_ew_mass_gev),
        m_z_gev=float(p.m_z_gev),
        alpha_em_mz=float(p.alpha_em_mz),
        sin2_theta_w=float(p.sin2_theta_w),
        g_z=float(p.g_z),
        s_z=float(p.s_z),
        a_ref=a_ref,
        a_ref_c=float(p.a_ref_c),
        contact_units="GeV^-2",
        includes_heavy_neutral_exchange=False,
        includes_heavy_neutral_lepton=False,
        z_delta_g_L_u=z_delta_l_u,
        z_delta_g_R_u=z_delta_r_u,
        z_delta_g_L_d=z_delta_l_d,
        z_delta_g_R_d=z_delta_r_d,
        z_delta_g_L_e=z_delta_l_e,
        z_delta_g_R_e=z_delta_r_e,
        z_delta_g_L_nu=z_delta_l_nu,
        z_total_g_L_u=g_l_u * identity + z_delta_l_u,
        z_total_g_R_u=g_r_u * identity + z_delta_r_u,
        z_total_g_L_d=g_l_d * identity + z_delta_l_d,
        z_total_g_R_d=g_r_d * identity + z_delta_r_d,
        z_total_g_L_e=g_l_e * identity + z_delta_l_e,
        z_total_g_R_e=g_r_e * identity + z_delta_r_e,
        z_total_g_L_nu=g_l_nu * identity + z_delta_l_nu,
        neutral_contacts=contacts,
        a_profile_values=a_profiles,
        a_mass_basis=a_mass,
        metadata={
            "requested_model_label": str(model_label),
            "ew_model": "minimal_rs",
            "custodial_protection_included": False,
            "brane_kinetic_terms_included": False,
            "fermion_kk_mixing_included": bool(zbb_fermion_mixing is not None),
            "minimal_rs_tree_zbb_complete": bool(zbb_fermion_mixing is not None),
            "custodial_variant_needs_human": True,
            "custodial_toppartner_zbL_needs_human": True,
            "overlap_rel_tol": float(overlap_rel_tol),
            "min_overlap_modes": int(min_modes),
            "max_overlap_modes": int(max_modes),
            "a_ref_interpretation": (
                "EW-universal subtraction spectrum.a(DEFAULT_A_REF_C=0.65)"
            ),
            "lepton_matching_status": dict(lepton_matching_status),
            "full_light_z_product_minus_sm_contacts": True,
            "z_delta_formula": (
                "s_Z * g_A_SM * (m_Z^2/M_KK^2) * U^dagger "
                "diag(a(c)-a_ref) U"
            ),
            "zbb_fermion_kk_mixing": (
                None
                if zbb_fermion_mixing is None
                else _zbb_fermion_kk_mixing_metadata(zbb_fermion_mixing)
            ),
        },
    )


def build_rs_zbb_fermion_kk_mixing(
    quark_fit_result: Any,
    *,
    spectrum: RSEWSpectrum,
) -> RSZbbFermionKKMixing:
    """Build the Phase-6a minimal non-custodial Zbb fermion-KK shift."""

    bulk_state = _required_attr(quark_fit_result, "bulk_state")
    c_q = _real_triplet_required_attr(bulk_state, "c_Q")
    c_d = _real_triplet_required_attr(bulk_state, "c_d")
    f_q = _positive_real_triplet_required_attr(bulk_state, "F_Q")
    f_d = _positive_real_triplet_required_attr(bulk_state, "F_d")
    y_d = _complex_matrix_required_attr(bulk_state, "Y_d_bulk_basis")
    masses_down = _real_triplet_required_attr(quark_fit_result, "masses_down")
    m_b = float(masses_down[2])
    if m_b < 0.0:
        raise ValueError("masses_down[2] must be non-negative")

    lambda_ir = _positive_float(
        _required_attr(spectrum, "lambda_ir_gev"),
        "spectrum.lambda_ir_gev",
    )
    y33_abs_sq = float(abs(y_d[2, 2]) ** 2)
    if not math.isfinite(y33_abs_sq) or y33_abs_sq <= 0.0:
        raise ValueError("Y_d_bulk_basis[2,2] must be non-zero and finite")

    profile_b_q = _casagrande_zbb_B_profile_triplet(c_q, f_q, name="B_Q")
    profile_b_d = _casagrande_zbb_B_profile_triplet(c_d, f_d, name="B_d")
    row_ratio = np.array(np.abs(y_d[2, :]) ** 2 / y33_abs_sq, dtype=float)
    column_ratio = np.array(np.abs(y_d[:, 2]) ** 2 / y33_abs_sq, dtype=float)
    if not np.all(np.isfinite(row_ratio)) or not np.all(np.isfinite(column_ratio)):
        raise ValueError("Y_d_bulk_basis Yukawa-ratio sums contain non-finite entries")

    # Casagrande's Zbb convention cross-assigns the singlet tower sum to
    # delta g_L^b and the doublet tower sum to delta g_R^b; keep B_d/B_Q
    # names explicit to avoid swapping the chiral shifts in future edits.
    B_d = float(np.dot(row_ratio, profile_b_d))
    B_Q = float(np.dot(column_ratio, profile_b_q))
    prefactor = float(m_b * m_b / (2.0 * lambda_ir * lambda_ir))
    delta_g_L_b = float(prefactor * B_d)
    delta_g_R_b = float(-prefactor * B_Q)

    return RSZbbFermionKKMixing(
        model_label=RS_ZBB_FERMION_KK_MIXING_MODEL_V1,
        matching_assumption=RS_ZBB_FERMION_KK_MIXING_MATCHING_ASSUMPTION_V1,
        lambda_ir_gev=lambda_ir,
        m_b_gev=m_b,
        c_Q=c_q,
        c_d=c_d,
        F_Q=f_q,
        F_d=f_d,
        Y_d_bulk_basis=y_d,
        profile_B_Q=profile_b_q,
        profile_B_d=profile_b_d,
        yukawa_ratio_column_b=column_ratio,
        yukawa_ratio_row_b=row_ratio,
        B_Q=B_Q,
        B_d=B_d,
        delta_g_L_b=delta_g_L_b,
        delta_g_R_b=delta_g_R_b,
        metadata={
            "source": "Casagrande et al. arXiv:0807.4937 minimal-RS Zbb ZMA",
            "M_KK_convention": "geometric Lambda_IR = spectrum.lambda_ir_gev",
            "physical_first_gauge_mass_gev": float(spectrum.kk_ew_mass_gev),
            "includes_exact_fermion_tower_diagonalization": False,
            "fermion_kk_mixing_included": True,
            "custodial_protection_included": False,
            "custodial_variant_needs_human": True,
            "custodial_toppartner_zbL_needs_human": True,
            "brane_kinetic_terms_included": False,
            "row_column_indexing": "zero-based b index 2",
        },
    )


def _resolve_max_overlap_modes(
    spectrum: RSEWSpectrum, requested: int | None
) -> int:
    if requested is not None:
        max_modes = int(requested)
    else:
        max_modes = min(int(spectrum.n_gauge_modes), DEFAULT_MAX_TRUNCATION_MODES)
        max_modes = 1 << (max_modes.bit_length() - 1)
    if max_modes < 32:
        raise ValueError("max_overlap_modes must be at least 32")
    if max_modes > int(spectrum.n_gauge_modes):
        raise ValueError("max_overlap_modes exceeds spectrum.n_gauge_modes")
    if max_modes & (max_modes - 1):
        raise ValueError("max_overlap_modes must be a power of two")
    return max_modes


def _real_triplet_from_attr(source: Any, name: str) -> np.ndarray:
    return _readonly_real_triplet(getattr(source, name), name).copy()


def _required_attr(source: Any, name: str) -> Any:
    try:
        return getattr(source, name)
    except AttributeError as exc:
        raise ValueError(f"{name} is required for RS Zbb fermion-KK mixing") from exc


def _real_triplet_required_attr(source: Any, name: str) -> np.ndarray:
    return _readonly_real_triplet(_required_attr(source, name), name).copy()


def _positive_real_triplet_required_attr(source: Any, name: str) -> np.ndarray:
    arr = _real_triplet_required_attr(source, name)
    if np.any(arr <= 0.0):
        raise ValueError(f"{name} must be strictly positive")
    return arr


def _complex_matrix_required_attr(source: Any, name: str) -> np.ndarray:
    return _readonly_complex_matrix(_required_attr(source, name), name).copy()


def _positive_float(value: Any, name: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0.0:
        raise ValueError(f"{name} must be positive and finite")
    return number


def _broadcast_real_triplet(values: Any, name: str) -> np.ndarray:
    arr = np.array(values, dtype=float, copy=True)
    if arr.shape == ():
        arr = np.full(3, float(arr), dtype=float)
    return _readonly_real_triplet(arr, name)


def _readonly_real_triplet(values: Any, name: str) -> np.ndarray:
    arr = np.array(values, dtype=float, copy=True)
    if arr.shape != (3,):
        raise ValueError(f"{name} must have shape (3,)")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} contains non-finite values")
    arr.setflags(write=False)
    return arr


def _readonly_complex_triplet(values: Any, name: str) -> np.ndarray:
    arr = np.array(values, dtype=np.complex128, copy=True)
    if arr.shape != (3,):
        raise ValueError(f"{name} must have shape (3,)")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} contains non-finite values")
    arr.setflags(write=False)
    return arr


def _unitary_matrix_from_attr(
    source: Any, name: str, *, tolerance: float
) -> np.ndarray:
    arr = np.array(getattr(source, name), dtype=np.complex128, copy=True)
    if arr.shape != (3, 3):
        raise ValueError(f"{name} must have shape (3, 3)")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} contains non-finite values")
    identity = np.eye(3, dtype=np.complex128)
    if not np.allclose(arr.conjugate().T @ arr, identity, rtol=0.0, atol=tolerance):
        raise ValueError(f"{name} must be unitary within {tolerance:g}")
    return arr


def _readonly_complex_matrix(values: Any, name: str) -> np.ndarray:
    arr = np.array(values, dtype=np.complex128, copy=True)
    if arr.shape != (3, 3):
        raise ValueError(f"{name} must have shape (3, 3)")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} contains non-finite values")
    arr.setflags(write=False)
    return arr


def _a_triplet(
    spectrum: RSEWSpectrum,
    c_values: np.ndarray,
    *,
    a_ref: float,
    rel_tol: float,
    min_modes: int,
    max_modes: int,
) -> np.ndarray:
    values = np.array(
        [
            spectrum.a(
                float(c),
                a_ref=float(a_ref),
                rel_tol=float(rel_tol),
                min_modes=int(min_modes),
                max_modes=int(max_modes),
            )
            for c in c_values
        ],
        dtype=float,
    )
    if not np.all(np.isfinite(values)):
        raise ValueError("a(c)-a_ref values contain non-finite entries")
    return values


def _mass_basis_profile(rotation: np.ndarray, a_values: np.ndarray) -> np.ndarray:
    return _hermitian(rotation.conjugate().T @ np.diag(a_values) @ rotation)


def _hermitian(matrix: np.ndarray) -> np.ndarray:
    arr = np.asarray(matrix, dtype=np.complex128)
    return 0.5 * (arr + arr.conjugate().T)


def _z_delta(g_sm: float, a_mass_basis: np.ndarray, *, scale: float, s_z: float) -> np.ndarray:
    return _hermitian(float(s_z) * float(g_sm) * float(scale) * a_mass_basis)


def _casagrande_zbb_B_profile_triplet(
    c_values: np.ndarray,
    F_values: np.ndarray,
    *,
    name: str,
) -> np.ndarray:
    values = np.array(
        [
            _casagrande_zbb_B_profile(float(c), float(f), name=f"{name}[{idx}]")
            for idx, (c, f) in enumerate(zip(c_values, F_values, strict=True))
        ],
        dtype=float,
    )
    values.setflags(write=False)
    return values


def _casagrande_zbb_B_profile(c: float, F: float, *, name: str) -> float:
    if not math.isfinite(c):
        raise ValueError(f"{name} c must be finite")
    if not math.isfinite(F) or F <= 0.0:
        raise ValueError(f"{name} F must be positive and finite")
    denom_left = 1.0 - 2.0 * c
    denom_right = 3.0 + 2.0 * c
    if denom_left == 0.0 or denom_right == 0.0:
        raise ValueError(f"{name} Casagrande B(c) denominator is singular")
    value = (1.0 / denom_left) * (1.0 / (F * F) - 1.0 + (F * F) / denom_right)
    if not math.isfinite(value):
        raise ValueError(f"{name} Casagrande B(c) is non-finite")
    return float(value)


def _zbb_fermion_kk_mixing_metadata(
    mixing: RSZbbFermionKKMixing,
) -> dict[str, Any]:
    return {
        "model_label": mixing.model_label,
        "matching_assumption": mixing.matching_assumption,
        "lambda_ir_gev": float(mixing.lambda_ir_gev),
        "m_b_gev": float(mixing.m_b_gev),
        "B_Q": float(mixing.B_Q),
        "B_d": float(mixing.B_d),
        "delta_g_L_b": float(mixing.delta_g_L_b),
        "delta_g_R_b": float(mixing.delta_g_R_b),
        "profile_B_Q": np.array(mixing.profile_B_Q, dtype=float, copy=True),
        "profile_B_d": np.array(mixing.profile_B_d, dtype=float, copy=True),
        "yukawa_ratio_column_b": np.array(
            mixing.yukawa_ratio_column_b,
            dtype=float,
            copy=True,
        ),
        "yukawa_ratio_row_b": np.array(
            mixing.yukawa_ratio_row_b,
            dtype=float,
            copy=True,
        ),
        **dict(mixing.metadata),
    }


def _sm_chiral_z_coupling(species: str, chirality: str, sin2_theta_w: float) -> float:
    s2 = float(sin2_theta_w)
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
        if chirality != "L":
            raise ValueError("active neutrino Z coupling is LH only")
        charge = 0.0
        t3 = 0.5
    else:
        raise ValueError(f"unsupported species {species!r}")
    return float(t3 - charge * s2)


def _contact_key(quark_sector: str, quark_chirality: str, lepton_chirality: str) -> str:
    sector = str(quark_sector)
    if sector not in {"u", "d"}:
        raise ValueError("quark_sector must be 'u' or 'd'")
    q_ch = str(quark_chirality).upper()
    l_ch = str(lepton_chirality).upper()
    if q_ch not in {"L", "R"} or l_ch not in {"L", "R"}:
        raise ValueError("chiralities must be 'L' or 'R'")
    return f"{sector}_{q_ch}{l_ch}"


def _nunu_contact_key(quark_sector: str, quark_chirality: str) -> str:
    sector = str(quark_sector)
    if sector not in {"u", "d"}:
        raise ValueError("quark_sector must be 'u' or 'd'")
    q_ch = str(quark_chirality).upper()
    if q_ch not in {"L", "R"}:
        raise ValueError("quark_chirality must be 'L' or 'R'")
    return f"{sector}_{q_ch}L_nunu"


def _neutral_contact_tensor(
    *,
    g_z: float,
    m_z_gev: float,
    g_q_sm: float,
    g_l_sm: float,
    z_delta_q: np.ndarray,
    z_delta_l: np.ndarray,
) -> np.ndarray:
    prefactor = float(g_z) ** 2 / float(m_z_gev) ** 2
    tensor = np.zeros((3, 3, 3, 3), dtype=np.complex128)
    for i in range(3):
        for j in range(3):
            delta_ij = 1.0 if i == j else 0.0
            for a in range(3):
                for b in range(3):
                    delta_ab = 1.0 if a == b else 0.0
                    tensor[i, j, a, b] = prefactor * (
                        (g_q_sm * delta_ij + z_delta_q[i, j])
                        * (g_l_sm * delta_ab + z_delta_l[a, b])
                        - g_q_sm * g_l_sm * delta_ij * delta_ab
                    )
    return tensor


def _neutral_contacts(
    *,
    g_z: float,
    m_z_gev: float,
    sin2_theta_w: float,
    z_delta_g_L_u: np.ndarray,
    z_delta_g_R_u: np.ndarray,
    z_delta_g_L_d: np.ndarray,
    z_delta_g_R_d: np.ndarray,
    z_delta_g_L_e: np.ndarray,
    z_delta_g_R_e: np.ndarray,
    z_delta_g_L_nu: np.ndarray,
) -> dict[str, np.ndarray]:
    g_l_e = _sm_chiral_z_coupling("e", "L", sin2_theta_w)
    g_r_e = _sm_chiral_z_coupling("e", "R", sin2_theta_w)
    g_l_nu = _sm_chiral_z_coupling("nu", "L", sin2_theta_w)
    return {
        "u_LL": _neutral_contact_tensor(
            g_z=g_z,
            m_z_gev=m_z_gev,
            g_q_sm=_sm_chiral_z_coupling("u", "L", sin2_theta_w),
            g_l_sm=g_l_e,
            z_delta_q=z_delta_g_L_u,
            z_delta_l=z_delta_g_L_e,
        ),
        "u_LR": _neutral_contact_tensor(
            g_z=g_z,
            m_z_gev=m_z_gev,
            g_q_sm=_sm_chiral_z_coupling("u", "L", sin2_theta_w),
            g_l_sm=g_r_e,
            z_delta_q=z_delta_g_L_u,
            z_delta_l=z_delta_g_R_e,
        ),
        "u_RL": _neutral_contact_tensor(
            g_z=g_z,
            m_z_gev=m_z_gev,
            g_q_sm=_sm_chiral_z_coupling("u", "R", sin2_theta_w),
            g_l_sm=g_l_e,
            z_delta_q=z_delta_g_R_u,
            z_delta_l=z_delta_g_L_e,
        ),
        "u_RR": _neutral_contact_tensor(
            g_z=g_z,
            m_z_gev=m_z_gev,
            g_q_sm=_sm_chiral_z_coupling("u", "R", sin2_theta_w),
            g_l_sm=g_r_e,
            z_delta_q=z_delta_g_R_u,
            z_delta_l=z_delta_g_R_e,
        ),
        "d_LL": _neutral_contact_tensor(
            g_z=g_z,
            m_z_gev=m_z_gev,
            g_q_sm=_sm_chiral_z_coupling("d", "L", sin2_theta_w),
            g_l_sm=g_l_e,
            z_delta_q=z_delta_g_L_d,
            z_delta_l=z_delta_g_L_e,
        ),
        "d_LR": _neutral_contact_tensor(
            g_z=g_z,
            m_z_gev=m_z_gev,
            g_q_sm=_sm_chiral_z_coupling("d", "L", sin2_theta_w),
            g_l_sm=g_r_e,
            z_delta_q=z_delta_g_L_d,
            z_delta_l=z_delta_g_R_e,
        ),
        "d_RL": _neutral_contact_tensor(
            g_z=g_z,
            m_z_gev=m_z_gev,
            g_q_sm=_sm_chiral_z_coupling("d", "R", sin2_theta_w),
            g_l_sm=g_l_e,
            z_delta_q=z_delta_g_R_d,
            z_delta_l=z_delta_g_L_e,
        ),
        "d_RR": _neutral_contact_tensor(
            g_z=g_z,
            m_z_gev=m_z_gev,
            g_q_sm=_sm_chiral_z_coupling("d", "R", sin2_theta_w),
            g_l_sm=g_r_e,
            z_delta_q=z_delta_g_R_d,
            z_delta_l=z_delta_g_R_e,
        ),
        "u_LL_nunu": _neutral_contact_tensor(
            g_z=g_z,
            m_z_gev=m_z_gev,
            g_q_sm=_sm_chiral_z_coupling("u", "L", sin2_theta_w),
            g_l_sm=g_l_nu,
            z_delta_q=z_delta_g_L_u,
            z_delta_l=z_delta_g_L_nu,
        ),
        "u_RL_nunu": _neutral_contact_tensor(
            g_z=g_z,
            m_z_gev=m_z_gev,
            g_q_sm=_sm_chiral_z_coupling("u", "R", sin2_theta_w),
            g_l_sm=g_l_nu,
            z_delta_q=z_delta_g_R_u,
            z_delta_l=z_delta_g_L_nu,
        ),
        "d_LL_nunu": _neutral_contact_tensor(
            g_z=g_z,
            m_z_gev=m_z_gev,
            g_q_sm=_sm_chiral_z_coupling("d", "L", sin2_theta_w),
            g_l_sm=g_l_nu,
            z_delta_q=z_delta_g_L_d,
            z_delta_l=z_delta_g_L_nu,
        ),
        "d_RL_nunu": _neutral_contact_tensor(
            g_z=g_z,
            m_z_gev=m_z_gev,
            g_q_sm=_sm_chiral_z_coupling("d", "R", sin2_theta_w),
            g_l_sm=g_l_nu,
            z_delta_q=z_delta_g_R_d,
            z_delta_l=z_delta_g_L_nu,
        ),
    }


__all__ = [
    "DEFAULT_A_REF_C",
    "DEFAULT_S_Z",
    "LEPTON_FLAVORS",
    "NEUTRINO_FLAVORS",
    "RSLeptonMassBasisCouplings",
    "RSEWMassBasisCouplings",
    "RSEWNeutralCurrentInputs",
    "RSZbbFermionKKMixing",
    "RS_EW_COUPLINGS_INPUT_BUNDLE_V1",
    "RS_EW_COUPLINGS_MATCHING_ASSUMPTION_V1",
    "RS_EW_COUPLINGS_MODEL_V1",
    "RS_LEPTON_MASS_BASIS_INPUT_BUNDLE_V1",
    "RS_LEPTON_MASS_BASIS_MATCHING_ASSUMPTION_V1",
    "RS_LEPTON_MASS_BASIS_MODEL_V1",
    "RS_ZBB_FERMION_KK_MIXING_MATCHING_ASSUMPTION_V1",
    "RS_ZBB_FERMION_KK_MIXING_MODEL_V1",
    "build_rs_ew_couplings",
    "build_rs_lepton_mass_basis_couplings",
    "build_rs_zbb_fermion_kk_mixing",
]
