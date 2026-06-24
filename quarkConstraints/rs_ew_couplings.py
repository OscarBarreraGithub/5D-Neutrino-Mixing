"""RS electroweak neutral-current couplings.

This module is the builder bridge from the numerical RS-EW spectrum to
mass-basis light-Z coupling shifts and semileptonic contacts.  It intentionally
does not call any rare-decay proxy helper.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Any, Mapping, Sequence

import numpy as np

from .rs_ew_spectrum import (
    DEFAULT_MAX_TRUNCATION_MODES,
    DEFAULT_MIN_TRUNCATION_MODES,
    DEFAULT_OVERLAP_RTOL,
    RSEWSpectrum,
)
from warpConfig.baseParams import V_EWSB
from warpConfig.wavefuncs import f_IR


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
MINIMAL_RS_EW_MODEL = "minimal_rs"
CUSTODIAL_RS_PLR_EW_MODEL = "custodial_rs_plr"
SUPPORTED_RS_EW_MODELS: tuple[str, str] = (
    MINIMAL_RS_EW_MODEL,
    CUSTODIAL_RS_PLR_EW_MODEL,
)
DEFAULT_CUSTODIAL_PROTECT_SCOPE = "all_gen_down_left_diagonal"
DEFAULT_CUSTODIAL_Q_L_REP = "all_gen_Q_L_bidoublet_(2,2)_{2/3}"
DEFAULT_CUSTODIAL_T_R_REP = "t_R_singlet_(1,1)_{2/3}"
DEFAULT_CUSTODIAL_B_R_REP = "b_R_elementary"
DEFAULT_CUSTODIAL_B_R_STRATEGY = "elementary_zero"
CUSTODIAL_FCNC_PR1_MINIMAL_OFFDIAG = "pr1_minimal_offdiag"
CUSTODIAL_FCNC_ALL_GEN_BIDOUBLET_PROXY = "all_gen_bidoublet_mass_basis_proxy"
TOP_PARTNER_LOOP_SOURCE = (
    "Carena-Ponton-Santiago-Wagner PhysRevD76 035006, arXiv:hep-ph/0701055"
)
TOP_PARTNER_LOOP_FORMULA_SET = (
    "leading_singlet_main_text_plus_bidoublet_vertex; "
    "full_Teq_Zbbeq_not_reconstructed"
)


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
class RSCustodialTopPartnerLoopProxy:
    """Tagged leading top-partner loop proxy for custodial RS."""

    top_partner_loop_mode: str
    include_top_partner_loops: bool
    top_partner_loop_magnitudes_computed: bool
    top_partner_zbb_loop_numerics_included: bool
    top_partner_t_loop_numerics_included: bool
    top_partner_loop_components_requested: str
    top_partner_loop_components_applied: tuple[str, ...]
    top_partner_loop_components_deferred: tuple[str, ...]
    top_partner_loop_t_sign: str
    top_partner_loop_t_source: str
    top_partner_loop_sign_convention: str
    top_partner_delta_t_singlet_magnitude: float
    top_partner_delta_t_override: float | None
    top_partner_delta_t_loop_applied: float
    top_partner_delta_g_L_b_singlet: float
    top_partner_delta_g_L_b_bidoublet_vertex: float
    top_partner_delta_g_L_b_loop_applied: float
    top_partner_delta_g_R_b_loop: float
    top_partner_delta_g_R_b_loop_reason: str
    top_partner_loop_inputs: Mapping[str, Any]
    top_partner_loop_mass_ratios: Mapping[str, Any]
    top_partner_loop_mixing_scales: Mapping[str, Any]

    def __post_init__(self) -> None:
        for name in (
            "top_partner_delta_t_singlet_magnitude",
            "top_partner_delta_t_loop_applied",
            "top_partner_delta_g_L_b_singlet",
            "top_partner_delta_g_L_b_bidoublet_vertex",
            "top_partner_delta_g_L_b_loop_applied",
            "top_partner_delta_g_R_b_loop",
        ):
            value = float(getattr(self, name))
            if not math.isfinite(value):
                raise ValueError(f"{name} must be finite")
            object.__setattr__(self, name, value)
        if self.top_partner_delta_t_override is not None:
            override = float(self.top_partner_delta_t_override)
            if not math.isfinite(override):
                raise ValueError("top_partner_delta_t_override must be finite")
            object.__setattr__(self, "top_partner_delta_t_override", override)
        object.__setattr__(
            self,
            "top_partner_loop_components_applied",
            tuple(str(item) for item in self.top_partner_loop_components_applied),
        )
        object.__setattr__(
            self,
            "top_partner_loop_components_deferred",
            tuple(str(item) for item in self.top_partner_loop_components_deferred),
        )
        object.__setattr__(
            self,
            "top_partner_loop_inputs",
            MappingProxyType(dict(self.top_partner_loop_inputs)),
        )
        object.__setattr__(
            self,
            "top_partner_loop_mass_ratios",
            MappingProxyType(dict(self.top_partner_loop_mass_ratios)),
        )
        object.__setattr__(
            self,
            "top_partner_loop_mixing_scales",
            MappingProxyType(dict(self.top_partner_loop_mixing_scales)),
        )

    @property
    def top_partner_loop_numerics_included(self) -> bool:
        """EW001-facing alias: true only when a ``Delta T`` loop is applied."""
        return bool(self.top_partner_t_loop_numerics_included)

    def metadata(self) -> dict[str, Any]:
        """Return serializable diagnostic metadata for the coupling bundle."""
        return {
            "top_partner_loop_mode": str(self.top_partner_loop_mode),
            "include_top_partner_loops": bool(self.include_top_partner_loops),
            "top_partner_loop_magnitudes_computed": bool(
                self.top_partner_loop_magnitudes_computed
            ),
            "top_partner_zbb_loop_numerics_included": bool(
                self.top_partner_zbb_loop_numerics_included
            ),
            "top_partner_t_loop_numerics_included": bool(
                self.top_partner_t_loop_numerics_included
            ),
            "top_partner_loop_numerics_included": bool(
                self.top_partner_loop_numerics_included
            ),
            "top_partner_loop_source": TOP_PARTNER_LOOP_SOURCE,
            "top_partner_loop_formula_set": TOP_PARTNER_LOOP_FORMULA_SET,
            "top_partner_loop_components_requested": str(
                self.top_partner_loop_components_requested
            ),
            "top_partner_loop_components_applied": list(
                self.top_partner_loop_components_applied
            ),
            "top_partner_loop_components_deferred": list(
                self.top_partner_loop_components_deferred
            ),
            "top_partner_loop_t_sign": str(self.top_partner_loop_t_sign),
            "top_partner_loop_t_source": str(self.top_partner_loop_t_source),
            "top_partner_loop_sign_convention": str(
                self.top_partner_loop_sign_convention
            ),
            "top_partner_loop_exact_Teq_Zbbeq_included": False,
            "top_partner_loop_proxy_status": "proxy_not_full_custodian_spectrum",
            "top_partner_delta_t_singlet_magnitude": float(
                self.top_partner_delta_t_singlet_magnitude
            ),
            "top_partner_delta_t_override": self.top_partner_delta_t_override,
            "top_partner_delta_t_loop_applied": float(
                self.top_partner_delta_t_loop_applied
            ),
            "top_partner_delta_g_L_b_singlet": float(
                self.top_partner_delta_g_L_b_singlet
            ),
            "top_partner_delta_g_L_b_bidoublet_vertex": float(
                self.top_partner_delta_g_L_b_bidoublet_vertex
            ),
            "top_partner_delta_g_L_b_loop_applied": float(
                self.top_partner_delta_g_L_b_loop_applied
            ),
            "top_partner_delta_g_R_b_loop": float(self.top_partner_delta_g_R_b_loop),
            "top_partner_delta_g_R_b_loop_reason": str(
                self.top_partner_delta_g_R_b_loop_reason
            ),
            "top_partner_loop_applied_to_z_delta_g_L_d_22": bool(
                self.top_partner_zbb_loop_numerics_included
            ),
            "top_partner_loop_inputs": dict(self.top_partner_loop_inputs),
            "top_partner_loop_mass_ratios": dict(self.top_partner_loop_mass_ratios),
            "top_partner_loop_mixing_scales": dict(self.top_partner_loop_mixing_scales),
        }


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
    qL_rep: str = DEFAULT_CUSTODIAL_Q_L_REP,
    tR_rep: str = DEFAULT_CUSTODIAL_T_R_REP,
    bR_rep: str = DEFAULT_CUSTODIAL_B_R_REP,
    protect_scope: Any = DEFAULT_CUSTODIAL_PROTECT_SCOPE,
    bR_strategy: str = DEFAULT_CUSTODIAL_B_R_STRATEGY,
    kappa_b: float = 0.0,
    custodial_PLR_breaking_residual: bool = False,
    include_top_partner_loops: bool = False,
    top_partner_loop_t_sign: Any | None = None,
    top_partner_loop_delta_t_override: float | None = None,
    top_partner_loop_components: str = "defer",
    top_partner_loop_mass_ratios: Mapping[str, Any] | None = None,
    top_partner_loop_mixing_scales: Mapping[str, Any] | None = None,
    custodial_fcnc_mode: str = CUSTODIAL_FCNC_PR1_MINIMAL_OFFDIAG,
    kappa_fcnc: float = 0.0,
    rotation_unitarity_tolerance: float = 1.0e-8,
) -> RSEWMassBasisCouplings:
    """Build Phase-4a RS-EW mass-basis couplings from a quark fit result."""

    ew_model = _validate_ew_model(model_label)
    if include_top_partner_loops and ew_model == MINIMAL_RS_EW_MODEL:
        raise ValueError("include_top_partner_loops=True requires ew_model='custodial_rs_plr'")
    if ew_model == CUSTODIAL_RS_PLR_EW_MODEL:
        _validate_custodial_options(
            qL_rep=qL_rep,
            tR_rep=tR_rep,
            bR_rep=bR_rep,
            bR_strategy=bR_strategy,
            kappa_b=kappa_b,
        )

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
    minimal_z_delta_l_d_full = np.array(z_delta_l_d, dtype=np.complex128, copy=True)
    minimal_z_delta_r_d_full = np.array(z_delta_r_d, dtype=np.complex128, copy=True)
    custodial_metadata: dict[str, Any] | None = None
    if ew_model == CUSTODIAL_RS_PLR_EW_MODEL:
        top_partner_loop_proxy = _build_custodial_top_partner_loop_proxy(
            quark_fit_result,
            spectrum=spectrum,
            inputs=p,
            include_top_partner_loops=bool(include_top_partner_loops),
            top_partner_loop_t_sign=top_partner_loop_t_sign,
            top_partner_loop_delta_t_override=top_partner_loop_delta_t_override,
            top_partner_loop_components=str(top_partner_loop_components),
            top_partner_loop_mass_ratios=top_partner_loop_mass_ratios,
            top_partner_loop_mixing_scales=top_partner_loop_mixing_scales,
        )
        z_delta_l_d, z_delta_r_d, custodial_metadata = _apply_custodial_rs_plr_proxy(
            z_delta_l_d=z_delta_l_d,
            z_delta_r_d=z_delta_r_d,
            minimal_z_delta_l_d_full=minimal_z_delta_l_d_full,
            minimal_z_delta_r_d_full=minimal_z_delta_r_d_full,
            spectrum=spectrum,
            qL_rep=qL_rep,
            tR_rep=tR_rep,
            bR_rep=bR_rep,
            protect_scope=protect_scope,
            bR_strategy=bR_strategy,
            kappa_b=float(kappa_b),
            custodial_PLR_breaking_residual=bool(custodial_PLR_breaking_residual),
            include_top_partner_loops=bool(include_top_partner_loops),
            top_partner_loop_proxy=top_partner_loop_proxy,
            custodial_fcnc_mode=str(custodial_fcnc_mode),
            kappa_fcnc=float(kappa_fcnc),
        )
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
            "ew_model": (
                MINIMAL_RS_EW_MODEL
                if ew_model == MINIMAL_RS_EW_MODEL
                else CUSTODIAL_RS_PLR_EW_MODEL
            ),
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
            **({} if custodial_metadata is None else custodial_metadata),
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

    # Per-generation diagonal CGHNP brackets (metadata; the b_R/b_L diagonal
    # term used in the sum is the 3rd-generation entry, profile_b_*[2]).
    profile_b_q = _casagrande_zbb_B_profile_triplet(c_q, f_q, name="B_Q")
    profile_b_d = _casagrande_zbb_B_profile_triplet(c_d, f_d, name="B_d")
    row_ratio = np.array(np.abs(y_d[2, :]) ** 2 / y33_abs_sq, dtype=float)
    column_ratio = np.array(np.abs(y_d[:, 2]) ** 2 / y33_abs_sq, dtype=float)
    if not np.all(np.isfinite(row_ratio)) or not np.all(np.isfinite(column_ratio)):
        raise ValueError("Y_d_bulk_basis Yukawa-ratio sums contain non-finite entries")

    # CGHNP (0807.4937) Z->bb ZMA, retranslated (PLAN §4.2):
    #   B_d = B_correct(c_d3, f_d3)
    #         + (1/(2 f_d3^2)) * sum_{i=1,2} |Y_d,3i|^2/|Y_d,33|^2 * 1/(1 + 2 c_d_i)
    #   (the 1/(1 - 2 c_di) of Eq. (170) is in the CGHNP c-convention; with
    #    c_CGHNP = -c_repo this becomes 1/(1 + 2 c_di,repo), matching the
    #    converted diagonal bracket below.)
    # and the symmetric expression for B_Q with c_Q, f_q3, column_ratio.  The
    # diagonal (b_R/b_L) term carries the 3rd-gen bracket; the light-generation
    # flavour sum carries the COMMON 1/(2 f_3^2) factor (the 3rd-gen singlet
    # overlap), NOT a per-light-gen full bracket B(c_d_i, F_d_i) -- the former
    # is the CGHNP structure, the latter inflated i=1,2 by F^2(c_d3)/F^2(c_d_i)
    # ~ 10^2..10^4 (PLAN §4.2 error 1a).
    f_d3_sq = float(f_d[2]) ** 2
    f_q3_sq = float(f_q[2]) ** 2
    light_sum_d = 0.0
    light_sum_q = 0.0
    for i in (0, 1):
        # CGHNP (0807.4937) Eq. (170) flavour-sum denominator is 1/(1 - 2 c_di)
        # in THEIR convention c_di,CGHNP.  The convention dictionary is
        # c_CGHNP = -c_repo (proved in audit slice 3), so 1 - 2 c_di,CGHNP =
        # 1 + 2 c_di,repo.  The diagonal bracket _casagrande_zbb_B_profile
        # already applies this conversion (it uses 1 + 2c); the light-generation
        # flavour sum MUST use the same converted denominator.  Using the
        # un-converted 1 - 2 c_repo here was a bug: for the scan-typical UV light
        # singlets (c_repo ~ +0.6, i.e. c_CGHNP ~ -0.6) it gives a NEGATIVE,
        # large denominator factor 1/(1 - 2*0.6) = -3.3, flipping the sign and
        # inflating |delta g_L^b| by ~5-8x (e.g. -0.044 vs the correct +0.007 at
        # the Bauer-S1 central point, M_KK = 3 TeV).
        denom_d = 1.0 + 2.0 * float(c_d[i])
        denom_q = 1.0 + 2.0 * float(c_q[i])
        if denom_d == 0.0 or denom_q == 0.0:
            raise ValueError("Zbb light-generation flavour-sum denominator (1 + 2c) is singular")
        light_sum_d += float(row_ratio[i]) / denom_d
        light_sum_q += float(column_ratio[i]) / denom_q
    B_d = float(profile_b_d[2] + (1.0 / (2.0 * f_d3_sq)) * light_sum_d)
    B_Q = float(profile_b_q[2] + (1.0 / (2.0 * f_q3_sq)) * light_sum_q)
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


def _validate_ew_model(model_label: str) -> str:
    model = str(model_label)
    if model not in SUPPORTED_RS_EW_MODELS:
        raise ValueError(
            f"unsupported ew_model/model_label {model!r}; "
            f"supported models are {SUPPORTED_RS_EW_MODELS}"
        )
    return model


def _validate_custodial_options(
    *,
    qL_rep: str,
    tR_rep: str,
    bR_rep: str,
    bR_strategy: str,
    kappa_b: float,
) -> None:
    if not str(qL_rep):
        raise ValueError("qL_rep must be a non-empty string")
    if not str(tR_rep):
        raise ValueError("tR_rep must be a non-empty string")
    if not str(bR_rep):
        raise ValueError("bR_rep must be a non-empty string")
    if str(bR_strategy) != DEFAULT_CUSTODIAL_B_R_STRATEGY:
        raise ValueError(
            "PR1 supports only bR_strategy='elementary_zero'; "
            "explicit custodial b_R representations are deferred to PR2"
        )
    if not math.isfinite(float(kappa_b)):
        raise ValueError("kappa_b must be finite")


def _build_custodial_top_partner_loop_proxy(
    quark_fit_result: Any,
    *,
    spectrum: RSEWSpectrum,
    inputs: RSEWNeutralCurrentInputs,
    include_top_partner_loops: bool,
    top_partner_loop_t_sign: Any | None,
    top_partner_loop_delta_t_override: float | None,
    top_partner_loop_components: str,
    top_partner_loop_mass_ratios: Mapping[str, Any] | None,
    top_partner_loop_mixing_scales: Mapping[str, Any] | None,
) -> RSCustodialTopPartnerLoopProxy:
    requested = str(top_partner_loop_components).strip()
    allowed = {"defer", "singlet", "bidoublet_vertex", "singlet_plus_bidoublet_vertex"}
    if requested not in allowed:
        raise ValueError(f"unsupported top_partner_loop_components {requested!r}")

    ratios = _top_partner_loop_knobs(
        top_partner_loop_mass_ratios,
        defaults={"rho_t": 1.0, "rho_q": 1.0, "rho_chi": 1.0},
        name="top_partner_loop_mass_ratios",
    )
    mixings = _top_partner_loop_knobs(
        top_partner_loop_mixing_scales,
        defaults={"xi_s": 1.0, "xi_q": 1.0, "xi_chi": 1.0},
        name="top_partner_loop_mixing_scales",
    )
    override = _finite_or_none(
        top_partner_loop_delta_t_override,
        "top_partner_loop_delta_t_override",
    )
    sign = _normalize_top_partner_t_sign(top_partner_loop_t_sign)
    if bool(include_top_partner_loops) and sign == -1 and override is None:
        raise ValueError(
            "top_partner_loop_t_sign=-1 requires a finite "
            "top_partner_loop_delta_t_override"
        )

    if not bool(include_top_partner_loops) or requested == "defer":
        return RSCustodialTopPartnerLoopProxy(
            top_partner_loop_mode="deferred",
            include_top_partner_loops=bool(include_top_partner_loops),
            top_partner_loop_magnitudes_computed=False,
            top_partner_zbb_loop_numerics_included=False,
            top_partner_t_loop_numerics_included=False,
            top_partner_loop_components_requested=requested,
            top_partner_loop_components_applied=(),
            top_partner_loop_components_deferred=(),
            top_partner_loop_t_sign=_top_partner_sign_metadata(sign, override),
            top_partner_loop_t_source="not_applicable_vertex_only",
            top_partner_loop_sign_convention="",
            top_partner_delta_t_singlet_magnitude=0.0,
            top_partner_delta_t_override=override,
            top_partner_delta_t_loop_applied=0.0,
            top_partner_delta_g_L_b_singlet=0.0,
            top_partner_delta_g_L_b_bidoublet_vertex=0.0,
            top_partner_delta_g_L_b_loop_applied=0.0,
            top_partner_delta_g_R_b_loop=0.0,
            top_partner_delta_g_R_b_loop_reason=(
                "not present in leading Carena Zb_L proxy"
            ),
            top_partner_loop_inputs={},
            top_partner_loop_mass_ratios={
                **ratios,
                "ratio_denominator": "physical_M_KK",
            },
            top_partner_loop_mixing_scales=mixings,
        )

    bulk_state = getattr(quark_fit_result, "bulk_state")
    c_q = _real_triplet_from_attr(bulk_state, "c_Q")
    c_u = _real_triplet_from_attr(bulk_state, "c_u")
    f_q = _profile_triplet_from_attr_or_fallback(
        bulk_state,
        "F_Q",
        c_q,
        spectrum=spectrum,
    )
    f_u = _profile_triplet_from_attr_or_fallback(
        bulk_state,
        "F_u",
        c_u,
        spectrum=spectrum,
    )
    masses_up = _real_triplet_for_top_partner_loop(quark_fit_result, "masses_up")
    m_t = float(masses_up[2])
    if m_t <= 0.0:
        raise ValueError("masses_up[2] must be positive for custodial top-partner loops")

    m_kk = _positive_float(getattr(spectrum, "kk_ew_mass_gev"), "spectrum.kk_ew_mass_gev")
    lambda_ir = _positive_float(getattr(spectrum, "lambda_ir_gev"), "spectrum.lambda_ir_gev")
    warp_log = _positive_float(getattr(spectrum, "warp_log"), "spectrum.warp_log")
    epsilon = _positive_float(getattr(spectrum, "epsilon"), "spectrum.epsilon")
    sin2 = float(inputs.sin2_theta_w)
    c_w2 = 1.0 - sin2
    if c_w2 <= 0.0:
        raise ValueError("inputs.sin2_theta_w must be less than one")
    m_z = _positive_float(inputs.m_z_gev, "inputs.m_z_gev")
    m_w = float(m_z * math.sqrt(c_w2))
    alpha = _positive_float(inputs.alpha_em_mz, "inputs.alpha_em_mz")
    f_q3 = _positive_float(f_q[2], "F_Q3")
    f_u3 = _positive_float(f_u[2], "F_u3")
    y_t_eff = float(m_t / (2.0 * float(V_EWSB) * f_q3 * f_u3))
    if not math.isfinite(y_t_eff) or y_t_eff <= 0.0:
        raise ValueError("Y_t_eff must be positive and finite")

    m_t_partner = float(ratios["rho_t"] * m_kk)
    m_q_partner = float(ratios["rho_q"] * m_kk)
    m_chi_partner = float(ratios["rho_chi"] * m_kk)
    m_q0t_t = float(mixings["xi_s"] * m_t / f_u3)
    m_qt_t = float(mixings["xi_q"] * m_t / f_q3)
    m_chid_t = float(mixings["xi_chi"] * m_t / f_q3)
    for name, value in (
        ("M_t", m_t_partner),
        ("M_q", m_q_partner),
        ("M_chi", m_chi_partner),
        ("m_q0t_t", m_q0t_t),
        ("m_qt_t", m_qt_t),
        ("m_chid_t", m_chid_t),
    ):
        _positive_float(value, name)

    t_top = float(3.0 * m_t * m_t / (16.0 * math.pi * sin2 * c_w2 * m_z * m_z))
    log_t = math.log((m_t_partner * m_t_partner) / (m_t * m_t))
    delta_t_singlet = float(
        t_top
        * (2.0 * m_q0t_t * m_q0t_t / (m_t_partner * m_t_partner))
        * (log_t - 1.0 + (m_q0t_t * m_q0t_t) / (2.0 * m_t * m_t))
    )
    delta_g_singlet = float(
        alpha
        / (16.0 * math.pi * sin2 * m_w * m_w)
        * (m_q0t_t**4 / (m_t_partner * m_t_partner))
        * (
            1.0
            + 2.0
            * m_t
            * m_t
            / (m_q0t_t * m_q0t_t)
            * (log_t - 1.0)
        )
    )
    delta_g_bidoublet = float(
        alpha
        / (32.0 * math.pi * sin2 * m_w * m_w)
        * m_t
        * m_t
        * (
            (m_qt_t * m_qt_t / (m_q_partner * m_q_partner))
            * math.log((m_q_partner * m_q_partner) / (m_t * m_t))
            - (m_chid_t * m_chid_t / (m_chi_partner * m_chi_partner))
            * math.log((m_chi_partner * m_chi_partner) / (m_t * m_t))
        )
    )
    for name, value in (
        ("top_partner_delta_t_singlet_magnitude", delta_t_singlet),
        ("top_partner_delta_g_L_b_singlet", delta_g_singlet),
        ("top_partner_delta_g_L_b_bidoublet_vertex", delta_g_bidoublet),
    ):
        if not math.isfinite(value):
            raise ValueError(f"{name} must be finite")

    requested_has_singlet = requested in {"singlet", "singlet_plus_bidoublet_vertex"}
    requested_has_bidoublet = requested in {
        "bidoublet_vertex",
        "singlet_plus_bidoublet_vertex",
    }
    t_input_valid = bool(override is not None or sign == +1)
    if requested_has_singlet and not t_input_valid:
        deferred = ("singlet",)
        if requested_has_bidoublet:
            deferred = ("singlet", "bidoublet_vertex")
        return RSCustodialTopPartnerLoopProxy(
            top_partner_loop_mode="computed_not_applied_missing_t_input",
            include_top_partner_loops=True,
            top_partner_loop_magnitudes_computed=True,
            top_partner_zbb_loop_numerics_included=False,
            top_partner_t_loop_numerics_included=False,
            top_partner_loop_components_requested=requested,
            top_partner_loop_components_applied=(),
            top_partner_loop_components_deferred=deferred,
            top_partner_loop_t_sign="None",
            top_partner_loop_t_source="missing_t_input_not_applied",
            top_partner_loop_sign_convention="",
            top_partner_delta_t_singlet_magnitude=delta_t_singlet,
            top_partner_delta_t_override=override,
            top_partner_delta_t_loop_applied=0.0,
            top_partner_delta_g_L_b_singlet=delta_g_singlet,
            top_partner_delta_g_L_b_bidoublet_vertex=delta_g_bidoublet,
            top_partner_delta_g_L_b_loop_applied=0.0,
            top_partner_delta_g_R_b_loop=0.0,
            top_partner_delta_g_R_b_loop_reason=(
                "not present in leading Carena Zb_L proxy"
            ),
            top_partner_loop_inputs=_top_partner_loop_inputs_metadata(
                m_t=m_t,
                c_q3=float(c_q[2]),
                f_q3=f_q3,
                c_u3=float(c_u[2]),
                f_u3=f_u3,
                y_t_eff=y_t_eff,
                m_z=m_z,
                m_w=m_w,
                alpha=alpha,
                sin2=sin2,
                c_w2=c_w2,
                m_kk=m_kk,
                lambda_ir=lambda_ir,
                warp_log=warp_log,
            ),
            top_partner_loop_mass_ratios={
                **ratios,
                "ratio_denominator": "physical_M_KK",
            },
            top_partner_loop_mixing_scales=mixings,
        )

    applied: list[str] = []
    deferred_list: list[str] = []
    delta_g_applied = 0.0
    if requested_has_singlet:
        delta_g_applied += delta_g_singlet
        applied.append("singlet")
    if requested_has_bidoublet:
        delta_g_applied += delta_g_bidoublet
        applied.append("bidoublet_vertex")

    if override is not None:
        delta_t_applied = float(override)
        t_included = True
        t_source = "explicit_numeric_override"
    elif requested_has_singlet and sign == +1:
        delta_t_applied = delta_t_singlet
        t_included = True
        t_source = "explicit_singlet_positive"
    else:
        delta_t_applied = 0.0
        t_included = False
        t_source = "not_applicable_vertex_only"

    zbb_included = bool(applied)
    sign_meta = _top_partner_sign_metadata(sign, override)
    sign_convention = (
        "Carena leading Zbb shifts are additive in repo g_L; Delta T uses "
        "positive singlet sign or explicit numeric override only"
        if (zbb_included or t_included)
        else ""
    )
    return RSCustodialTopPartnerLoopProxy(
        top_partner_loop_mode="carena_leading_proxy_applied",
        include_top_partner_loops=True,
        top_partner_loop_magnitudes_computed=True,
        top_partner_zbb_loop_numerics_included=zbb_included,
        top_partner_t_loop_numerics_included=t_included,
        top_partner_loop_components_requested=requested,
        top_partner_loop_components_applied=tuple(applied),
        top_partner_loop_components_deferred=tuple(deferred_list),
        top_partner_loop_t_sign=sign_meta,
        top_partner_loop_t_source=t_source,
        top_partner_loop_sign_convention=sign_convention,
        top_partner_delta_t_singlet_magnitude=delta_t_singlet,
        top_partner_delta_t_override=override,
        top_partner_delta_t_loop_applied=delta_t_applied,
        top_partner_delta_g_L_b_singlet=delta_g_singlet,
        top_partner_delta_g_L_b_bidoublet_vertex=delta_g_bidoublet,
        top_partner_delta_g_L_b_loop_applied=delta_g_applied,
        top_partner_delta_g_R_b_loop=0.0,
        top_partner_delta_g_R_b_loop_reason=(
            "not present in leading Carena Zb_L proxy"
        ),
        top_partner_loop_inputs=_top_partner_loop_inputs_metadata(
            m_t=m_t,
            c_q3=float(c_q[2]),
            f_q3=f_q3,
            c_u3=float(c_u[2]),
            f_u3=f_u3,
            y_t_eff=y_t_eff,
            m_z=m_z,
            m_w=m_w,
            alpha=alpha,
            sin2=sin2,
            c_w2=c_w2,
            m_kk=m_kk,
            lambda_ir=lambda_ir,
            warp_log=warp_log,
        ),
        top_partner_loop_mass_ratios={
            **ratios,
            "ratio_denominator": "physical_M_KK",
        },
        top_partner_loop_mixing_scales=mixings,
    )


def _finite_or_none(value: Any, name: str) -> float | None:
    if value is None:
        return None
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"{name} must be finite")
    return number


def _top_partner_loop_knobs(
    values: Mapping[str, Any] | None,
    *,
    defaults: Mapping[str, float],
    name: str,
) -> dict[str, float]:
    resolved = {str(key): float(value) for key, value in defaults.items()}
    if values is not None:
        unknown = set(values) - set(defaults)
        if unknown:
            raise ValueError(f"{name} contains unsupported keys: {sorted(unknown)}")
        for key, value in values.items():
            resolved[str(key)] = float(value)
    for key, value in resolved.items():
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError(f"{name}.{key} must be positive and finite")
    return resolved


def _normalize_top_partner_t_sign(value: Any | None) -> int | None:
    if value is None:
        return None
    if isinstance(value, str):
        cleaned = value.strip()
        if cleaned in {"+1", "1"}:
            return +1
        if cleaned == "-1":
            return -1
    number = float(value)
    if number == 1.0:
        return +1
    if number == -1.0:
        return -1
    raise ValueError("top_partner_loop_t_sign must be +1, -1, or None")


def _top_partner_sign_metadata(sign: int | None, override: float | None) -> str:
    if override is not None:
        if sign == -1:
            return "-1"
        return "override"
    if sign == +1:
        return "+1"
    if sign == -1:
        return "-1"
    return "None"


def _profile_triplet_from_attr_or_fallback(
    bulk_state: Any,
    attr_name: str,
    c_values: np.ndarray,
    *,
    spectrum: RSEWSpectrum,
) -> np.ndarray:
    if hasattr(bulk_state, attr_name):
        values = _readonly_real_triplet(getattr(bulk_state, attr_name), attr_name)
    else:
        values = _readonly_real_triplet(
            f_IR(c_values, _positive_float(getattr(spectrum, "epsilon"), "spectrum.epsilon")),
            attr_name,
        )
    if np.any(values <= 0.0):
        raise ValueError(f"{attr_name} must be strictly positive")
    return values


def _real_triplet_for_top_partner_loop(source: Any, name: str) -> np.ndarray:
    try:
        values = getattr(source, name)
    except AttributeError as exc:
        raise ValueError(f"{name} is required for custodial top-partner loops") from exc
    return _readonly_real_triplet(values, name).copy()


def _top_partner_loop_inputs_metadata(
    *,
    m_t: float,
    c_q3: float,
    f_q3: float,
    c_u3: float,
    f_u3: float,
    y_t_eff: float,
    m_z: float,
    m_w: float,
    alpha: float,
    sin2: float,
    c_w2: float,
    m_kk: float,
    lambda_ir: float,
    warp_log: float,
) -> dict[str, float | str]:
    return {
        "m_top_gev": float(m_t),
        "c_Q3": float(c_q3),
        "F_Q3": float(f_q3),
        "c_u3": float(c_u3),
        "F_u3": float(f_u3),
        "Y_t_eff": float(y_t_eff),
        "V_EWSB_gev": float(V_EWSB),
        "m_Z_gev": float(m_z),
        "m_W_gev": float(m_w),
        "alpha_em_mz": float(alpha),
        "sin2_theta_w": float(sin2),
        "cW2": float(c_w2),
        "M_KK_gev": float(m_kk),
        "lambda_ir_gev": float(lambda_ir),
        "warp_log": float(warp_log),
        "M_KK_convention": "physical first electroweak gauge KK mass",
    }


def _custodial_protected_down_left_mask(protect_scope: Any) -> np.ndarray:
    mask = np.zeros((3, 3), dtype=bool)
    if isinstance(protect_scope, str):
        scope = protect_scope.strip().lower().replace("-", "_").replace(" ", "_")
        if scope in {
            DEFAULT_CUSTODIAL_PROTECT_SCOPE,
            "all_gen_down_left_diag",
            "all_generation_down_left_diagonal",
        }:
            np.fill_diagonal(mask, True)
            return mask
        if scope in {"third_gen_down_left_diagonal", "b_left_diagonal", "zbb_only"}:
            mask[2, 2] = True
            return mask
        raise ValueError(f"unsupported custodial protect_scope {protect_scope!r}")

    arr = np.asarray(protect_scope)
    if arr.shape == (3, 3):
        mask = np.asarray(arr, dtype=bool)
    else:
        try:
            pairs: Sequence[Any] = list(protect_scope)
        except TypeError as exc:
            raise ValueError("protect_scope must be a supported string, mask, or index pairs") from exc
        for pair in pairs:
            if len(pair) != 2:
                raise ValueError("protect_scope index pairs must have length two")
            i = int(pair[0])
            j = int(pair[1])
            if not (0 <= i < 3 and 0 <= j < 3):
                raise ValueError("protect_scope indices must be in {0,1,2}")
            mask[i, j] = True

    offdiag = mask.copy()
    np.fill_diagonal(offdiag, False)
    if np.any(offdiag):
        raise ValueError("PR1 custodial protection can only target diagonal down-left entries")
    if not np.any(mask):
        raise ValueError("protect_scope must protect at least one diagonal entry")
    return mask


def _apply_custodial_rs_plr_proxy(
    *,
    z_delta_l_d: np.ndarray,
    z_delta_r_d: np.ndarray,
    minimal_z_delta_l_d_full: np.ndarray,
    minimal_z_delta_r_d_full: np.ndarray,
    spectrum: RSEWSpectrum,
    qL_rep: str,
    tR_rep: str,
    bR_rep: str,
    protect_scope: Any,
    bR_strategy: str,
    kappa_b: float,
    custodial_PLR_breaking_residual: bool,
    include_top_partner_loops: bool,
    top_partner_loop_proxy: RSCustodialTopPartnerLoopProxy,
    custodial_fcnc_mode: str,
    kappa_fcnc: float,
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    mask = _custodial_protected_down_left_mask(protect_scope)
    protected_indices = [
        [int(i), int(j)]
        for i in range(3)
        for j in range(3)
        if bool(mask[i, j])
    ]
    left = np.array(z_delta_l_d, dtype=np.complex128, copy=True)
    right = np.array(z_delta_r_d, dtype=np.complex128, copy=True)

    for i, j in protected_indices:
        left[i, j] = 0.0j

    volume_log = _positive_float(getattr(spectrum, "warp_log"), "spectrum.warp_log")
    residual_value = 0.0j
    residual_applied = False
    if bool(custodial_PLR_breaking_residual) and bool(mask[2, 2]):
        residual_value = complex(
            float(kappa_b) * (1.0 / volume_log) * minimal_z_delta_l_d_full[2, 2]
        )
        left[2, 2] = residual_value
        residual_applied = True

    if str(bR_strategy) == DEFAULT_CUSTODIAL_B_R_STRATEGY:
        right[2, 2] = 0.0j

    fcnc_mode = str(custodial_fcnc_mode)
    if fcnc_mode not in {
        CUSTODIAL_FCNC_PR1_MINIMAL_OFFDIAG,
        CUSTODIAL_FCNC_ALL_GEN_BIDOUBLET_PROXY,
    }:
        raise ValueError(f"unsupported custodial_fcnc_mode {custodial_fcnc_mode!r}")
    kappa_fcnc_value = float(kappa_fcnc)
    if not math.isfinite(kappa_fcnc_value):
        raise ValueError("kappa_fcnc must be finite")
    fcnc_metadata: dict[str, Any] = {}
    if fcnc_mode == CUSTODIAL_FCNC_ALL_GEN_BIDOUBLET_PROXY:
        before_left = np.array(left, dtype=np.complex128, copy=True)
        offdiag_mask = ~np.eye(3, dtype=bool)
        left[offdiag_mask] = 0.0j
        residual_matrix = np.zeros((3, 3), dtype=np.complex128)
        residual_matrix[offdiag_mask] = (
            kappa_fcnc_value / volume_log * minimal_z_delta_l_d_full[offdiag_mask]
        )
        left[offdiag_mask] = residual_matrix[offdiag_mask]
        fcnc_metadata = {
            "custodial_fcnc_modeling": "all_gen_bidoublet_mass_basis_proxy",
            "custodial_fcnc_mode": CUSTODIAL_FCNC_ALL_GEN_BIDOUBLET_PROXY,
            "custodial_fcnc_basis": "all_gen_bidoublet_mass_basis",
            "custodial_fcnc_leading_PLR_zeroed": True,
            "custodial_fcnc_residual_source": (
                "kappa_fcnc*(1/L)*minimal_z_delta_l_d_full[i,j]"
            ),
            "custodial_fcnc_residual_applied": bool(kappa_fcnc_value != 0.0),
            "custodial_fcnc_rh_status": (
                "right_handed_down_offdiagonal_kept_minimal"
            ),
            "kappa_fcnc": float(kappa_fcnc_value),
            "minimal_z_delta_l_d_full_offdiag": _offdiag_metadata(
                minimal_z_delta_l_d_full
            ),
            "minimal_z_delta_r_d_full_offdiag": _offdiag_metadata(
                minimal_z_delta_r_d_full
            ),
            "custodial_z_delta_l_d_offdiag_before_after": {
                "before": _offdiag_metadata(before_left),
                "after": _offdiag_metadata(left),
            },
        }

    z_delta_g_L_d_tree_pr1_b_before_top_partner = complex(left[2, 2])
    if bool(top_partner_loop_proxy.top_partner_zbb_loop_numerics_included):
        left[2, 2] += complex(
            top_partner_loop_proxy.top_partner_delta_g_L_b_loop_applied,
            0.0,
        )

    left = _hermitian(left)
    right = _hermitian(right)
    omissions = {
        "SU2_R_tower": False,
        "custodian_spectrum": False,
        "exact_NC_mixing": False,
        "BKT": False,
    }
    omission_details = {
        "SU2_R_tower": "not_included_in_PR2_proxy",
        "custodian_spectrum": "not_inferred_full_spectrum_absent",
        "exact_NC_mixing": "not_reconstructed",
        "BKT": "not_included",
        "full_Teq_Zbbeq_loop_matching": "not_reconstructed",
    }
    loop_metadata = top_partner_loop_proxy.metadata()
    metadata = {
        "ew_model": CUSTODIAL_RS_PLR_EW_MODEL,
        "custodial_protection_included": True,
        "custodial_proxy_scope": "tree_level_P_LR_Zbb_diagonal_only",
        "qL_rep": str(qL_rep),
        "tR_rep": str(tR_rep),
        "bR_rep": str(bR_rep),
        "protect_scope": protect_scope if isinstance(protect_scope, str) else protected_indices,
        "protected_down_left_diagonal_mask": mask.astype(bool).tolist(),
        "protected_down_left_diagonal_indices": protected_indices,
        "minimal_z_delta_l_d_full_source": (
            "post minimal gauge-profile shift plus optional Casagrande Zbb admixture"
        ),
        "minimal_z_delta_l_d_full_b": complex(minimal_z_delta_l_d_full[2, 2]),
        "minimal_z_delta_r_d_full_b": complex(minimal_z_delta_r_d_full[2, 2]),
        "custodial_PLR_breaking_residual": bool(custodial_PLR_breaking_residual),
        "custodial_residual_source": "kappa_b*(1/L)*minimal_z_delta_l_d_full[2,2]",
        "custodial_residual_value": residual_value,
        "custodial_residual_applied": residual_applied,
        "kappa_b": float(kappa_b),
        "rs_volume_log": float(volume_log),
        "bR_strategy": str(bR_strategy),
        "bR_elementary_zero_applied": True,
        "include_top_partner_loops": bool(include_top_partner_loops),
        "top_partner_loops": (
            "carena_leading_proxy"
            if (
                top_partner_loop_proxy.top_partner_zbb_loop_numerics_included
                or top_partner_loop_proxy.top_partner_t_loop_numerics_included
            )
            else "deferred"
        ),
        **loop_metadata,
        "z_delta_g_L_d_tree_pr1_b_before_top_partner": (
            z_delta_g_L_d_tree_pr1_b_before_top_partner
        ),
        "custodial_toppartner_zbL_needs_human": not bool(
            top_partner_loop_proxy.top_partner_zbb_loop_numerics_included
        ),
        "custodial_variant_needs_human": False,
        "custodial_fcnc_modeling": "deferred_PR2_off_diagonal_kept_minimal",
        "custodial_omissions": omissions,
        "custodial_omission_details": omission_details,
        "SU2_R_tower": False,
        "custodian_spectrum": False,
        "exact_NC_mixing": False,
        "BKT": False,
        "custodian_spectrum_inferred": False,
        "exact_top_partner_mass_matrix_included": False,
        "brane_kinetic_terms_included": False,
        "full_Teq_Zbbeq_loop_matching_included": False,
        "M_KK_convention": "physical first electroweak gauge KK mass",
        "physical_first_gauge_mass_gev": float(spectrum.kk_ew_mass_gev),
        "lambda_ir_gev": float(spectrum.lambda_ir_gev),
        **fcnc_metadata,
    }
    return left, right, metadata


def _offdiag_metadata(matrix: np.ndarray) -> dict[str, complex]:
    arr = np.asarray(matrix, dtype=np.complex128)
    return {
        f"{i}{j}": complex(arr[i, j])
        for i in range(3)
        for j in range(3)
        if i != j
    }


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
    """Diagonal (3rd-generation) CGHNP Z->bb bracket in repo variables.

    CGHNP (0807.4937) Z->bb ZMA, translated through the convention dictionary
    ``c_CGHNP = -c_repo``, ``F^2_CGHNP = 2 f_IR,repo^2`` (proved exactly in audit
    slice 3).  In repo variables the diagonal bracket is

        B(c, F) = 1/(1 + 2c) * ( 1/(2 F^2) - 1 + 2 F^2/(3 - 2c) ).

    The previous code used ``1/(1 - 2c) * (1/F^2 - 1 + F^2/(3 + 2c))`` -- the
    c-sign was wrong in BOTH denominators and the F^2 = 2 f^2 factor was
    missing.  For UV-localized b_R (c > 1/2, scan-typical) the old ``1/(1 - 2c)``
    is negative, giving the wrong SIGN of delta g_L^b (PLAN §4.2).
    """
    if not math.isfinite(c):
        raise ValueError(f"{name} c must be finite")
    if not math.isfinite(F) or F <= 0.0:
        raise ValueError(f"{name} F must be positive and finite")
    denom_left = 1.0 + 2.0 * c
    denom_right = 3.0 - 2.0 * c
    if denom_left == 0.0 or denom_right == 0.0:
        raise ValueError(f"{name} Casagrande B(c) denominator is singular")
    f_sq = F * F
    value = (1.0 / denom_left) * (1.0 / (2.0 * f_sq) - 1.0 + (2.0 * f_sq) / denom_right)
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
    "RSCustodialTopPartnerLoopProxy",
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
