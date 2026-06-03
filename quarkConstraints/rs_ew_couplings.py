"""RS electroweak quark neutral-current couplings.

This module is the Phase-3a bridge from the numerical RS-EW spectrum to
mass-basis light-Z coupling shifts and same-flavor charged-lepton contacts.
It intentionally does not call any rare-decay proxy helper.
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
RS_EW_COUPLINGS_MODEL_V1 = "RS_EW_NEUTRAL_CURRENT_PHASE3A_V1"
RS_EW_COUPLINGS_MATCHING_ASSUMPTION_V1 = (
    "minimal-RS light-Z quark neutral currents with charged-lepton SM Z "
    "couplings; charged-lepton delta_g and heavy-neutral lepton exchange "
    "are deferred to Phase 4"
)

DEFAULT_A_REF_C = 0.65
DEFAULT_S_Z = -1.0
LEPTON_FLAVORS: tuple[str, str, str] = ("e", "mu", "tau")


@dataclass(frozen=True)
class RSEWNeutralCurrentInputs:
    """Shared numerical inputs for Phase-3a neutral-current matching."""

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
    z_total_g_L_u: np.ndarray
    z_total_g_R_u: np.ndarray
    z_total_g_L_d: np.ndarray
    z_total_g_R_d: np.ndarray
    z_total_g_L_e: np.ndarray
    z_total_g_R_e: np.ndarray
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
            "z_total_g_L_u",
            "z_total_g_R_u",
            "z_total_g_L_d",
            "z_total_g_R_d",
            "z_total_g_L_e",
            "z_total_g_R_e",
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


def build_rs_ew_couplings(
    quark_fit_result: Any,
    *,
    spectrum: RSEWSpectrum,
    inputs: RSEWNeutralCurrentInputs | None = None,
    overlap_rel_tol: float = DEFAULT_OVERLAP_RTOL,
    min_overlap_modes: int = DEFAULT_MIN_TRUNCATION_MODES,
    max_overlap_modes: int | None = None,
    model_label: str = "minimal_rs",
    rotation_unitarity_tolerance: float = 1.0e-8,
) -> RSEWMassBasisCouplings:
    """Build Phase-3a RS-EW mass-basis couplings from a quark fit result."""

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

    z_delta_l_u = _z_delta(g_l_u, a_mass["L_u"], scale=scale, s_z=p.s_z)
    z_delta_r_u = _z_delta(g_r_u, a_mass["R_u"], scale=scale, s_z=p.s_z)
    z_delta_l_d = _z_delta(g_l_d, a_mass["L_d"], scale=scale, s_z=p.s_z)
    z_delta_r_d = _z_delta(g_r_d, a_mass["R_d"], scale=scale, s_z=p.s_z)
    z_delta_l_e = np.zeros((3, 3), dtype=np.complex128)
    z_delta_r_e = np.zeros((3, 3), dtype=np.complex128)

    contacts = _neutral_contacts(
        g_z=p.g_z,
        m_z_gev=p.m_z_gev,
        sin2_theta_w=p.sin2_theta_w,
        z_delta_g_L_u=z_delta_l_u,
        z_delta_g_R_u=z_delta_r_u,
        z_delta_g_L_d=z_delta_l_d,
        z_delta_g_R_d=z_delta_r_d,
    )

    identity = np.eye(3, dtype=np.complex128)
    return RSEWMassBasisCouplings(
        model_label=RS_EW_COUPLINGS_MODEL_V1,
        input_bundle=p.input_bundle,
        matching_assumption=RS_EW_COUPLINGS_MATCHING_ASSUMPTION_V1,
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
        z_total_g_L_u=g_l_u * identity + z_delta_l_u,
        z_total_g_R_u=g_r_u * identity + z_delta_r_u,
        z_total_g_L_d=g_l_d * identity + z_delta_l_d,
        z_total_g_R_d=g_r_d * identity + z_delta_r_d,
        z_total_g_L_e=g_l_e * identity + z_delta_l_e,
        z_total_g_R_e=g_r_e * identity + z_delta_r_e,
        neutral_contacts=contacts,
        a_profile_values=a_profiles,
        a_mass_basis=a_mass,
        metadata={
            "requested_model_label": str(model_label),
            "ew_model": "minimal_rs",
            "custodial_protection_included": False,
            "brane_kinetic_terms_included": False,
            "fermion_kk_mixing_included": False,
            "overlap_rel_tol": float(overlap_rel_tol),
            "min_overlap_modes": int(min_modes),
            "max_overlap_modes": int(max_modes),
            "z_delta_formula": (
                "s_Z * g_A_SM * (m_Z^2/M_KK^2) * U^dagger "
                "diag(a(c)-a_ref) U"
            ),
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


def _readonly_real_triplet(values: Any, name: str) -> np.ndarray:
    arr = np.array(values, dtype=float, copy=True)
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


def _neutral_contact_tensor(
    *,
    g_z: float,
    m_z_gev: float,
    g_q_sm: float,
    g_l_sm: float,
    z_delta_q: np.ndarray,
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
                        * (g_l_sm * delta_ab)
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
) -> dict[str, np.ndarray]:
    g_l_e = _sm_chiral_z_coupling("e", "L", sin2_theta_w)
    g_r_e = _sm_chiral_z_coupling("e", "R", sin2_theta_w)
    return {
        "u_LL": _neutral_contact_tensor(
            g_z=g_z,
            m_z_gev=m_z_gev,
            g_q_sm=_sm_chiral_z_coupling("u", "L", sin2_theta_w),
            g_l_sm=g_l_e,
            z_delta_q=z_delta_g_L_u,
        ),
        "u_LR": _neutral_contact_tensor(
            g_z=g_z,
            m_z_gev=m_z_gev,
            g_q_sm=_sm_chiral_z_coupling("u", "L", sin2_theta_w),
            g_l_sm=g_r_e,
            z_delta_q=z_delta_g_L_u,
        ),
        "u_RL": _neutral_contact_tensor(
            g_z=g_z,
            m_z_gev=m_z_gev,
            g_q_sm=_sm_chiral_z_coupling("u", "R", sin2_theta_w),
            g_l_sm=g_l_e,
            z_delta_q=z_delta_g_R_u,
        ),
        "u_RR": _neutral_contact_tensor(
            g_z=g_z,
            m_z_gev=m_z_gev,
            g_q_sm=_sm_chiral_z_coupling("u", "R", sin2_theta_w),
            g_l_sm=g_r_e,
            z_delta_q=z_delta_g_R_u,
        ),
        "d_LL": _neutral_contact_tensor(
            g_z=g_z,
            m_z_gev=m_z_gev,
            g_q_sm=_sm_chiral_z_coupling("d", "L", sin2_theta_w),
            g_l_sm=g_l_e,
            z_delta_q=z_delta_g_L_d,
        ),
        "d_LR": _neutral_contact_tensor(
            g_z=g_z,
            m_z_gev=m_z_gev,
            g_q_sm=_sm_chiral_z_coupling("d", "L", sin2_theta_w),
            g_l_sm=g_r_e,
            z_delta_q=z_delta_g_L_d,
        ),
        "d_RL": _neutral_contact_tensor(
            g_z=g_z,
            m_z_gev=m_z_gev,
            g_q_sm=_sm_chiral_z_coupling("d", "R", sin2_theta_w),
            g_l_sm=g_l_e,
            z_delta_q=z_delta_g_R_d,
        ),
        "d_RR": _neutral_contact_tensor(
            g_z=g_z,
            m_z_gev=m_z_gev,
            g_q_sm=_sm_chiral_z_coupling("d", "R", sin2_theta_w),
            g_l_sm=g_r_e,
            z_delta_q=z_delta_g_R_d,
        ),
    }


__all__ = [
    "DEFAULT_A_REF_C",
    "DEFAULT_S_Z",
    "LEPTON_FLAVORS",
    "RSEWMassBasisCouplings",
    "RSEWNeutralCurrentInputs",
    "RS_EW_COUPLINGS_INPUT_BUNDLE_V1",
    "RS_EW_COUPLINGS_MATCHING_ASSUMPTION_V1",
    "RS_EW_COUPLINGS_MODEL_V1",
    "build_rs_ew_couplings",
]
