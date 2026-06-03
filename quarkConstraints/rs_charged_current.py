"""RS electroweak charged-current couplings.

This module implements the Phase-5a data builder only.  It exposes the
charged W zero-plus-KK diagonalization, light-W vertex shifts, W' LL contacts,
and the muon-decay ``delta_G_F/G_F`` subtraction needed by later adapters.
No constraint is rewired here.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Any, Mapping

import numpy as np

from .rs_ew_couplings import DEFAULT_A_REF_C, RSLeptonMassBasisCouplings
from .rs_ew_spectrum import (
    DEFAULT_MAX_TRUNCATION_MODES,
    DEFAULT_MIN_TRUNCATION_MODES,
    DEFAULT_OVERLAP_RTOL,
    RSChargedGaugeDiagonalization,
    RSEWSpectrum,
)


RS_CHARGED_CURRENT_INPUT_BUNDLE_V1 = "quarkConstraints.rs_charged_current.inputs.v1"
RS_CHARGED_CURRENT_MODEL_V1 = "RS_EW_CHARGED_CURRENT_PHASE5A_V1"
RS_CHARGED_CURRENT_MATCHING_ASSUMPTION_V1 = (
    "minimal-RS LH charged currents with W zero-plus-KK diagonalization, "
    "eta_W pinned from the light-W eigenvector, CKM-normalized W' LL "
    "contacts, and a single muon-decay delta_G_F/G_F subtraction"
)

DEFAULT_ALPHA_EM_MZ = 1.0 / 127.952
DEFAULT_SIN2_THETA_W = 0.23122
DEFAULT_GF_GEV_MINUS2 = 1.1663787e-5
DEFAULT_V_HIGGS_GEV = 1.0 / math.sqrt(math.sqrt(2.0) * DEFAULT_GF_GEV_MINUS2)
_CKM_DENOM_FLOOR = 1.0e-14
_ZERO_ATOL = 1.0e-16


@dataclass(frozen=True)
class RSChargedCurrentInputs:
    """Numerical inputs for Phase-5a charged-current matching."""

    input_bundle: str = RS_CHARGED_CURRENT_INPUT_BUNDLE_V1
    alpha_em_mz: float = DEFAULT_ALPHA_EM_MZ
    sin2_theta_w: float = DEFAULT_SIN2_THETA_W
    gf_gev_minus2: float = DEFAULT_GF_GEV_MINUS2
    v_higgs_gev: float = DEFAULT_V_HIGGS_GEV
    a_ref_c: float = DEFAULT_A_REF_C

    def __post_init__(self) -> None:
        for name in ("alpha_em_mz", "sin2_theta_w", "gf_gev_minus2", "v_higgs_gev"):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be positive and finite")
            object.__setattr__(self, name, value)
        if not 0.0 < float(self.sin2_theta_w) < 1.0:
            raise ValueError("sin2_theta_w must lie between zero and one")
        object.__setattr__(self, "a_ref_c", float(self.a_ref_c))

    @property
    def g2(self) -> float:
        return float(math.sqrt(4.0 * math.pi * self.alpha_em_mz / self.sin2_theta_w))


@dataclass(frozen=True)
class RSChargedCurrentCouplings:
    """Mass-basis charged-current shifts and W' contact amplitudes."""

    model_label: str
    input_bundle: str
    matching_assumption: str
    kk_ew_mass_gev: float
    m_w_gev: float
    m_wprime_gev: float
    g2: float
    v_higgs_gev: float
    a_ref: float
    a_ref_c: float
    eta_W: float
    C_SM_gev_minus2: float
    delta_G_F_over_G_F: float
    contact_units: str
    ckm: np.ndarray
    w_diagonalization: RSChargedGaugeDiagonalization
    delta_g_W_ud_L: np.ndarray
    delta_g_W_ud_R: np.ndarray
    delta_g_W_ud_R_status: str
    delta_g_W_lnu_L: np.ndarray
    delta_g_W_lnu_R: np.ndarray
    delta_g_W_lnu_R_status: str
    charged_contact_LL: np.ndarray
    epsilon: np.ndarray
    delta_abs_vij_over_vij: np.ndarray
    a_profile_values: Mapping[str, np.ndarray]
    a_mass_basis: Mapping[str, np.ndarray]
    units: Mapping[str, str]
    input_parameters: Mapping[str, Any]
    diagnostics: Mapping[str, Any] = field(default_factory=dict)
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.contact_units != "GeV^-2":
            raise ValueError("contact_units must be 'GeV^-2'")
        for name in (
            "kk_ew_mass_gev",
            "m_w_gev",
            "m_wprime_gev",
            "g2",
            "v_higgs_gev",
            "C_SM_gev_minus2",
        ):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be positive and finite")
            object.__setattr__(self, name, value)
        for name in ("a_ref", "a_ref_c", "delta_G_F_over_G_F"):
            value = float(getattr(self, name))
            if not math.isfinite(value):
                raise ValueError(f"{name} must be finite")
            object.__setattr__(self, name, value)
        eta = float(self.eta_W)
        if eta not in {-1.0, 1.0}:
            raise ValueError("eta_W must be pinned to +/-1")
        object.__setattr__(self, "eta_W", eta)

        object.__setattr__(self, "ckm", _readonly_complex_matrix(self.ckm, "ckm"))
        for name in (
            "delta_g_W_ud_L",
            "delta_g_W_ud_R",
            "delta_g_W_lnu_L",
            "delta_g_W_lnu_R",
        ):
            object.__setattr__(
                self,
                name,
                _readonly_complex_matrix(getattr(self, name), name),
            )
        if not np.allclose(self.delta_g_W_ud_R, 0.0, rtol=0.0, atol=0.0):
            raise ValueError("delta_g_W_ud_R must be exactly zero in minimal RS")
        if not np.allclose(self.delta_g_W_lnu_R, 0.0, rtol=0.0, atol=0.0):
            raise ValueError("delta_g_W_lnu_R must be exactly zero in minimal RS")

        object.__setattr__(
            self,
            "charged_contact_LL",
            _readonly_complex_array(self.charged_contact_LL, "charged_contact_LL", (3, 3, 3, 3)),
        )
        object.__setattr__(
            self,
            "epsilon",
            _readonly_complex_array(self.epsilon, "epsilon", (3, 3, 3)),
        )
        object.__setattr__(
            self,
            "delta_abs_vij_over_vij",
            _readonly_real_array(
                self.delta_abs_vij_over_vij,
                "delta_abs_vij_over_vij",
                (3, 3, 3),
            ),
        )
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
        object.__setattr__(self, "units", MappingProxyType(dict(self.units)))
        object.__setattr__(self, "input_parameters", MappingProxyType(dict(self.input_parameters)))
        object.__setattr__(self, "diagnostics", _readonly_mapping(self.diagnostics))
        object.__setattr__(self, "metadata", MappingProxyType(dict(self.metadata)))

    def epsilon_ij_a(self, i: int, j: int, a: int) -> complex:
        """Return the charged-current amplitude shift epsilon_ij^a."""

        return complex(self.epsilon[int(i), int(j), int(a)])


def build_rs_charged_current(
    quark_fit_result: Any,
    *,
    spectrum: RSEWSpectrum,
    lepton_mass_basis_couplings: RSLeptonMassBasisCouplings,
    inputs: RSChargedCurrentInputs | None = None,
    shared_a_ref: float | None = None,
    overlap_rel_tol: float = DEFAULT_OVERLAP_RTOL,
    min_overlap_modes: int = DEFAULT_MIN_TRUNCATION_MODES,
    max_overlap_modes: int | None = None,
    model_label: str = "minimal_rs",
    rotation_unitarity_tolerance: float = 1.0e-8,
) -> RSChargedCurrentCouplings:
    """Build the Phase-5a charged-current typed extra."""

    if lepton_mass_basis_couplings is None:
        raise ValueError("charged-current matching requires lepton_mass_basis_couplings with c_L")

    p = RSChargedCurrentInputs() if inputs is None else inputs
    max_modes = _resolve_max_overlap_modes(spectrum, max_overlap_modes)
    min_modes = int(min_overlap_modes)
    if min_modes >= max_modes:
        raise ValueError("min_overlap_modes must be smaller than max_overlap_modes")

    bulk_state = getattr(quark_fit_result, "bulk_state")
    c_q = _readonly_real_triplet(getattr(bulk_state, "c_Q"), "c_Q")
    c_l = _readonly_real_triplet(lepton_mass_basis_couplings.c_L, "c_L")
    u_l_u = _unitary_matrix_from_attr(
        quark_fit_result,
        "U_L_u",
        tolerance=rotation_unitarity_tolerance,
    )
    u_l_d = _unitary_matrix_from_attr(
        quark_fit_result,
        "U_L_d",
        tolerance=rotation_unitarity_tolerance,
    )
    u_e_l = _unitary_matrix(
        lepton_mass_basis_couplings.U_e_L,
        "U_e_L",
        tolerance=rotation_unitarity_tolerance,
    )
    ckm = u_l_u.conjugate().T @ u_l_d
    if not np.allclose(ckm.conjugate().T @ ckm, np.eye(3), rtol=0.0, atol=rotation_unitarity_tolerance):
        raise ValueError("V_CKM=U_L_u^dagger U_L_d must be unitary")

    a_ref = (
        float(shared_a_ref)
        if shared_a_ref is not None
        else float(
            spectrum.a(
                p.a_ref_c,
                rel_tol=overlap_rel_tol,
                min_modes=min_modes,
                max_modes=max_modes,
            )
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
        "c_L": _a_triplet(
            spectrum,
            c_l,
            a_ref=a_ref,
            rel_tol=overlap_rel_tol,
            min_modes=min_modes,
            max_modes=max_modes,
        ),
    }
    a_q_ud = u_l_u.conjugate().T @ np.diag(a_profiles["c_Q"]) @ u_l_d
    a_l = _hermitian(u_e_l.conjugate().T @ np.diag(a_profiles["c_L"]) @ u_e_l)

    w_diag = spectrum.charged_w_diagonalization(
        g2=p.g2,
        v_higgs_gev=p.v_higgs_gev,
        max_modes=max_modes,
    )
    scale = (float(w_diag.m_w_gev) / float(spectrum.kk_ew_mass_gev)) ** 2
    delta_ud_l = complex(w_diag.eta_W * scale) * a_q_ud
    delta_lnu_l = _hermitian(complex(w_diag.eta_W * scale) * a_l)
    zero = np.zeros((3, 3), dtype=np.complex128)

    contact_norm, contact_actual, mu_decay_contact = _charged_contact_tensors(
        spectrum,
        c_q=c_q,
        c_l=c_l,
        U_L_u=u_l_u,
        U_L_d=u_l_d,
        U_e_L=u_e_l,
        ckm=ckm,
        g2=p.g2,
        max_modes=max_modes,
    )
    c_sm = float(p.g2 * p.g2 / (2.0 * w_diag.m_w_gev * w_diag.m_w_gev))
    delta_gf_light = complex(delta_lnu_l[0, 0] + delta_lnu_l[1, 1])
    delta_gf_contact = complex(mu_decay_contact / c_sm)
    delta_gf = _real_close(
        delta_gf_light + delta_gf_contact,
        "delta_G_F_over_G_F",
    )

    delta_ud_ratio = _safe_ratio_matrix(delta_ud_l, ckm, "delta_g_W_ud_L/V_CKM")
    epsilon = np.zeros((3, 3, 3), dtype=np.complex128)
    for i in range(3):
        for j in range(3):
            if abs(ckm[i, j]) <= _CKM_DENOM_FLOOR:
                if abs(delta_ud_l[i, j]) <= _ZERO_ATOL and np.max(np.abs(contact_actual[i, j])) <= _ZERO_ATOL:
                    continue
                raise ValueError(
                    "epsilon is undefined for a numerically zero V_CKM entry "
                    "with nonzero charged-current numerator"
                )
            for a in range(3):
                epsilon[i, j, a] = (
                    delta_ud_ratio[i, j]
                    + delta_lnu_l[a, a]
                    + contact_norm[i, j, a, a] / c_sm
                    - delta_gf
                )

    return RSChargedCurrentCouplings(
        model_label=RS_CHARGED_CURRENT_MODEL_V1,
        input_bundle=p.input_bundle,
        matching_assumption=RS_CHARGED_CURRENT_MATCHING_ASSUMPTION_V1,
        kk_ew_mass_gev=float(spectrum.kk_ew_mass_gev),
        m_w_gev=float(w_diag.m_w_gev),
        m_wprime_gev=float(w_diag.m_wprime_gev),
        g2=float(p.g2),
        v_higgs_gev=float(p.v_higgs_gev),
        a_ref=a_ref,
        a_ref_c=float(p.a_ref_c),
        eta_W=float(w_diag.eta_W),
        C_SM_gev_minus2=c_sm,
        delta_G_F_over_G_F=delta_gf,
        contact_units="GeV^-2",
        ckm=ckm,
        w_diagonalization=w_diag,
        delta_g_W_ud_L=delta_ud_l,
        delta_g_W_ud_R=zero,
        delta_g_W_ud_R_status="minimal_rs_no_right_handed_charged_current",
        delta_g_W_lnu_L=delta_lnu_l,
        delta_g_W_lnu_R=zero,
        delta_g_W_lnu_R_status="minimal_rs_no_right_handed_charged_lepton_neutrino_current",
        charged_contact_LL=contact_norm,
        epsilon=epsilon,
        delta_abs_vij_over_vij=np.real(epsilon),
        a_profile_values=a_profiles,
        a_mass_basis={
            "W_ud_L": a_q_ud,
            "W_lnu_L": a_l,
        },
        units={
            "m_w_gev": "GeV",
            "m_wprime_gev": "GeV",
            "C_SM_gev_minus2": "GeV^-2",
            "charged_contact_LL": "GeV^-2 (CKM-normalized C_Wprime/V_ij)",
            "delta_G_F_over_G_F": "dimensionless",
            "delta_g_W": "dimensionless relative coupling shift",
            "epsilon": "dimensionless",
        },
        input_parameters={
            "requested_model_label": str(model_label),
            "alpha_em_mz": float(p.alpha_em_mz),
            "sin2_theta_w": float(p.sin2_theta_w),
            "gf_gev_minus2": float(p.gf_gev_minus2),
            "v_higgs_gev": float(p.v_higgs_gev),
            "a_ref_c": float(p.a_ref_c),
            "overlap_rel_tol": float(overlap_rel_tol),
            "min_overlap_modes": int(min_modes),
            "max_overlap_modes": int(max_modes),
        },
        diagnostics={
            "w_diagonalization_location": "quarkConstraints.rs_ew_spectrum.RSEWSpectrum.charged_w_diagonalization",
            "a_ref_source": (
                "shared_a_ref_from_rs_ew_couplings"
                if shared_a_ref is not None
                else "spectrum.a(DEFAULT_A_REF_C)"
            ),
            "delta_GF_light_w_vertices": delta_gf_light,
            "delta_GF_charged_contact": delta_gf_contact,
            "delta_GF_charged_contact_gev_minus2": complex(mu_decay_contact),
            "delta_GF_C_SM_gev_minus2": c_sm,
            "charged_contact_LL_actual": contact_actual,
            "charged_contact_LL_storage": "CKM-normalized C_Wprime/V_ij",
            "epsilon_formula": (
                "delta_g_W_ud_L[i,j]/V_ij + delta_g_W_lnu_L[a,a] "
                "+ charged_contact_LL[i,j,a,a]/C_SM - delta_G_F_over_G_F"
            ),
            "pmns_handling_for_epsilon": (
                "charged-lepton flavor-basis W l nu shift after PMNS unitarity reduction"
            ),
        },
        metadata={
            "ew_model": "minimal_rs",
            "custodial_protection_included": False,
            "brane_kinetic_terms_included": False,
            "fermion_kk_mixing_included": False,
            "right_handed_W_couplings_included": False,
            "charged_contact_current_normalization": "L_W = g2/sqrt(2) W J",
            "C_SM_convention": "g2^2/(2*m_W^2) in the same L_W convention",
            "delta_G_F_subtraction": "subtract exactly once in consumed epsilon",
            "a_ref_interpretation": (
                "EW-universal subtraction shared with rs_ew_couplings.a_ref"
            ),
        },
    )


def _resolve_max_overlap_modes(spectrum: RSEWSpectrum, requested: int | None) -> int:
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


def _charged_contact_tensors(
    spectrum: RSEWSpectrum,
    *,
    c_q: np.ndarray,
    c_l: np.ndarray,
    U_L_u: np.ndarray,
    U_L_d: np.ndarray,
    U_e_L: np.ndarray,
    ckm: np.ndarray,
    g2: float,
    max_modes: int,
) -> tuple[np.ndarray, np.ndarray, complex]:
    omega_q = np.stack(
        [spectrum.omega(float(c), max_modes=max_modes) for c in c_q],
        axis=1,
    )
    omega_l = np.stack(
        [spectrum.omega(float(c), max_modes=max_modes) for c in c_l],
        axis=1,
    )
    normalized = np.zeros((3, 3, 3, 3), dtype=np.complex128)
    actual = np.zeros((3, 3, 3, 3), dtype=np.complex128)
    mu_decay_contact = 0.0j
    for mode_index in range(max_modes):
        prefactor = float(g2) * float(g2) / (
            2.0 * float(spectrum.gauge_masses_gev[mode_index]) ** 2
        )
        q_mode = U_L_u.conjugate().T @ np.diag(omega_q[mode_index]) @ U_L_d
        q_ratio = _safe_ratio_matrix(q_mode, ckm, "charged_contact_q_mode/V_CKM")
        l_mode = _hermitian(U_e_L.conjugate().T @ np.diag(omega_l[mode_index]) @ U_e_L)
        normalized += prefactor * q_ratio[:, :, None, None] * l_mode[None, None, :, :]
        actual += prefactor * q_mode[:, :, None, None] * l_mode[None, None, :, :]
        mu_decay_contact += prefactor * l_mode[0, 0] * l_mode[1, 1]
    return normalized, actual, complex(mu_decay_contact)


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


def _safe_ratio_matrix(numerator: np.ndarray, denominator: np.ndarray, label: str) -> np.ndarray:
    num = np.asarray(numerator, dtype=np.complex128)
    den = np.asarray(denominator, dtype=np.complex128)
    if num.shape != (3, 3) or den.shape != (3, 3):
        raise ValueError(f"{label} requires two 3x3 matrices")
    out = np.zeros((3, 3), dtype=np.complex128)
    mask = np.abs(den) > _CKM_DENOM_FLOOR
    out[mask] = num[mask] / den[mask]
    bad = (~mask) & (np.abs(num) > _ZERO_ATOL)
    if np.any(bad):
        raise ValueError(f"{label} has nonzero entries where V_CKM is numerically zero")
    return out


def _real_close(value: complex, name: str, *, atol: float = 1.0e-13) -> float:
    number = complex(value)
    if not math.isfinite(number.real) or not math.isfinite(number.imag):
        raise ValueError(f"{name} must be finite")
    if abs(number.imag) > atol:
        raise ValueError(f"{name} must be real within {atol:g}")
    return float(number.real)


def _unitary_matrix_from_attr(source: Any, name: str, *, tolerance: float) -> np.ndarray:
    return _unitary_matrix(getattr(source, name), name, tolerance=tolerance)


def _unitary_matrix(values: Any, name: str, *, tolerance: float) -> np.ndarray:
    arr = np.array(values, dtype=np.complex128, copy=True)
    if arr.shape != (3, 3):
        raise ValueError(f"{name} must have shape (3, 3)")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} contains non-finite values")
    identity = np.eye(3, dtype=np.complex128)
    if not np.allclose(arr.conjugate().T @ arr, identity, rtol=0.0, atol=tolerance):
        raise ValueError(f"{name} must be unitary within {tolerance:g}")
    return arr


def _hermitian(matrix: np.ndarray) -> np.ndarray:
    arr = np.asarray(matrix, dtype=np.complex128)
    return 0.5 * (arr + arr.conjugate().T)


def _readonly_real_triplet(values: Any, name: str) -> np.ndarray:
    arr = np.array(values, dtype=float, copy=True)
    if arr.shape != (3,):
        raise ValueError(f"{name} must have shape (3,)")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} contains non-finite values")
    arr.setflags(write=False)
    return arr


def _readonly_real_array(values: Any, name: str, shape: tuple[int, ...]) -> np.ndarray:
    arr = np.array(values, dtype=float, copy=True)
    if arr.shape != shape:
        raise ValueError(f"{name} must have shape {shape}")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} contains non-finite values")
    arr.setflags(write=False)
    return arr


def _readonly_complex_matrix(values: Any, name: str) -> np.ndarray:
    return _readonly_complex_array(values, name, (3, 3))


def _readonly_complex_array(values: Any, name: str, shape: tuple[int, ...]) -> np.ndarray:
    arr = np.array(values, dtype=np.complex128, copy=True)
    if arr.shape != shape:
        raise ValueError(f"{name} must have shape {shape}")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} contains non-finite values")
    arr.setflags(write=False)
    return arr


def _readonly_mapping(mapping: Mapping[str, Any]) -> Mapping[str, Any]:
    frozen: dict[str, Any] = {}
    for key, value in mapping.items():
        if isinstance(value, np.ndarray):
            arr = np.array(value, copy=True)
            arr.setflags(write=False)
            frozen[str(key)] = arr
        elif isinstance(value, Mapping):
            frozen[str(key)] = _readonly_mapping(value)
        else:
            frozen[str(key)] = value
    return MappingProxyType(frozen)


__all__ = [
    "DEFAULT_ALPHA_EM_MZ",
    "DEFAULT_GF_GEV_MINUS2",
    "DEFAULT_SIN2_THETA_W",
    "DEFAULT_V_HIGGS_GEV",
    "RSChargedCurrentCouplings",
    "RSChargedCurrentInputs",
    "RS_CHARGED_CURRENT_INPUT_BUNDLE_V1",
    "RS_CHARGED_CURRENT_MATCHING_ASSUMPTION_V1",
    "RS_CHARGED_CURRENT_MODEL_V1",
    "build_rs_charged_current",
]
