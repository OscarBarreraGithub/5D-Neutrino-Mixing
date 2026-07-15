"""RS charged-lepton Higgs-Yukawa couplings.

This module implements the Phase-6b leading minimal-RS Higgs-LFV matching
object.  The current repo lepton fit is diagonal, so the production builder
returns exactly zero off-diagonal Higgs Yukawas unless callers supply a
non-diagonal charged-lepton structure through the mass-basis lepton object.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from types import MappingProxyType
from typing import Any, Mapping

import numpy as np

from quarkConstraints.casagrande_profiles import casagrande_cghnp_B_profile
from yukawa.constants import LEPTON_MASSES


RS_HIGGS_YUKAWA_INPUT_BUNDLE_V1 = "quarkConstraints.rs_higgs_yukawas.inputs.v1"
RS_HIGGS_YUKAWA_MODEL_V1 = "RS_HIGGS_YUKAWA_PHASE6B_V1"
RS_HIGGS_YUKAWA_SOURCE_V1 = (
    "Casagrande-style minimal-RS brane-Higgs Delta g_h lepton analog"
)
RS_HIGGS_YUKAWA_MATCHING_ASSUMPTION_V1 = (
    "minimal-RS leading tree/ZMA charged-lepton Higgs-Yukawa correction from "
    "fermion-KK mixing using the supplied charged-lepton mass-basis rotations; "
    "no exact tower diagonalization, brane kinetic term, custodial variant, or "
    "new anarchic charged-lepton fit is inferred"
)
RS_HIGGS_YUKAWA_V1_HONESTY_NOTE = (
    "With the repo v1 diagonal charged-lepton fit (Y_E_bar_matrix diagonal and "
    "U_e_L=U_e_R=I), delta_L and delta_E are diagonal, so all tree Higgs-LFV "
    "off-diagonal Yukawas are exactly zero. Nonzero generic RS Higgs-LFV "
    "requires a non-diagonal charged-lepton Yukawa structure or rotations."
)


@dataclass(frozen=True)
class RSHiggsYukawaCouplings:
    """Minimal-RS charged-lepton Higgs-Yukawa matrix in the mass basis."""

    model_label: str
    input_bundle: str
    matching_assumption: str
    lambda_ir_gev: float
    kk_ew_mass_gev: float
    v_gev: float
    units: str
    source: str
    includes_fermion_kk_mixing: bool
    charged_lepton_masses_gev: np.ndarray
    c_L: np.ndarray
    c_E: np.ndarray
    F_L: np.ndarray
    F_E: np.ndarray
    x_e: np.ndarray
    profile_B_L: np.ndarray
    profile_B_E: np.ndarray
    delta_L: np.ndarray
    delta_E: np.ndarray
    higgs_yukawa_matrix: np.ndarray
    Y_E_bar_matrix: np.ndarray
    U_e_L: np.ndarray
    U_e_R: np.ndarray
    diagnostics: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        for name in ("lambda_ir_gev", "kk_ew_mass_gev", "v_gev"):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be positive and finite")
            object.__setattr__(self, name, value)
        if self.units != "dimensionless":
            raise ValueError("units must be 'dimensionless'")
        object.__setattr__(
            self,
            "includes_fermion_kk_mixing",
            bool(self.includes_fermion_kk_mixing),
        )

        for name in (
            "charged_lepton_masses_gev",
            "c_L",
            "c_E",
            "F_L",
            "F_E",
            "profile_B_L",
            "profile_B_E",
        ):
            object.__setattr__(self, name, _readonly_real_triplet(getattr(self, name), name))
        if np.any(self.charged_lepton_masses_gev < 0.0):
            raise ValueError("charged_lepton_masses_gev must be non-negative")
        if np.any(self.F_L <= 0.0) or np.any(self.F_E <= 0.0):
            raise ValueError("F_L and F_E must be strictly positive")

        for name in (
            "x_e",
            "delta_L",
            "delta_E",
            "higgs_yukawa_matrix",
            "Y_E_bar_matrix",
            "U_e_L",
            "U_e_R",
        ):
            object.__setattr__(
                self,
                name,
                _readonly_complex_matrix(getattr(self, name), name),
            )
        for name in ("U_e_L", "U_e_R"):
            _validate_unitary(getattr(self, name), name)
        for name in ("delta_L", "delta_E"):
            arr = getattr(self, name)
            if not np.allclose(arr, arr.conjugate().T, rtol=0.0, atol=5.0e-13):
                raise ValueError(f"{name} must be Hermitian")

        object.__setattr__(
            self,
            "diagnostics",
            MappingProxyType(dict(self.diagnostics)),
        )

    @property
    def Y_h_mass(self) -> np.ndarray:
        """Alias for callers using the explicit mass-basis label."""

        return self.higgs_yukawa_matrix


def build_rs_higgs_yukawas(
    lepton_mass_basis_couplings: Any,
    *,
    spectrum: Any,
    charged_lepton_masses_gev: tuple[float, float, float] = LEPTON_MASSES,
) -> RSHiggsYukawaCouplings:
    """Build the Phase-6b minimal-RS charged-lepton Higgs-Yukawa extra."""

    params = _required_mapping_attr(lepton_mass_basis_couplings, "params")
    v_gev = _positive_float(params.get("v"), "lepton_mass_basis_couplings.params['v']")
    lambda_ir = _positive_float(
        _required_attr(spectrum, "lambda_ir_gev"),
        "spectrum.lambda_ir_gev",
    )
    kk_ew_mass = _positive_float(
        _required_attr(spectrum, "kk_ew_mass_gev"),
        "spectrum.kk_ew_mass_gev",
    )

    masses = _readonly_real_triplet(charged_lepton_masses_gev, "charged_lepton_masses_gev")
    if np.any(masses < 0.0):
        raise ValueError("charged_lepton_masses_gev must be non-negative")
    c_l = _readonly_real_triplet(
        _required_attr(lepton_mass_basis_couplings, "c_L"),
        "c_L",
    )
    c_e = _readonly_real_triplet(
        _required_attr(lepton_mass_basis_couplings, "c_E"),
        "c_E",
    )
    f_l = _readonly_real_triplet(
        _required_attr(lepton_mass_basis_couplings, "f_L"),
        "f_L",
    )
    f_e = _readonly_real_triplet(
        _required_attr(lepton_mass_basis_couplings, "f_E"),
        "f_E",
    )
    u_e_l = _readonly_complex_matrix(
        _required_attr(lepton_mass_basis_couplings, "U_e_L"),
        "U_e_L",
    )
    u_e_r = _readonly_complex_matrix(
        _required_attr(lepton_mass_basis_couplings, "U_e_R"),
        "U_e_R",
    )
    _validate_unitary(u_e_l, "U_e_L")
    _validate_unitary(u_e_r, "U_e_R")
    y_e_bar_matrix = _readonly_complex_matrix(
        _required_attr(lepton_mass_basis_couplings, "Y_E_bar_matrix"),
        "Y_E_bar_matrix",
    )

    profile_b_l = _casagrande_higgs_B_profile_triplet(c_l, f_l, name="B(c_L)")
    profile_b_e = _casagrande_higgs_B_profile_triplet(c_e, f_e, name="B(c_E)")
    x_e = np.diag(masses / lambda_ir).astype(np.complex128)
    delta_l = _hermitian(
        x_e @ u_e_l.conjugate().T @ np.diag(profile_b_e) @ u_e_l @ x_e
    )
    delta_e = _hermitian(
        x_e @ u_e_r.conjugate().T @ np.diag(profile_b_l) @ u_e_r @ x_e
    )
    higgs_yukawa_matrix = _higgs_yukawa_matrix(
        charged_lepton_masses_gev=masses,
        v_gev=v_gev,
        delta_L=delta_l,
        delta_E=delta_e,
    )
    offdiag = _offdiag(higgs_yukawa_matrix)
    diagonal_yukawa = _is_diagonal(y_e_bar_matrix)
    identity_l = _is_identity(u_e_l)
    identity_r = _is_identity(u_e_r)
    offdiag_exact_zero = bool(np.count_nonzero(offdiag) == 0)

    return RSHiggsYukawaCouplings(
        model_label=RS_HIGGS_YUKAWA_MODEL_V1,
        input_bundle=RS_HIGGS_YUKAWA_INPUT_BUNDLE_V1,
        matching_assumption=RS_HIGGS_YUKAWA_MATCHING_ASSUMPTION_V1,
        lambda_ir_gev=lambda_ir,
        kk_ew_mass_gev=kk_ew_mass,
        v_gev=v_gev,
        units="dimensionless",
        source=RS_HIGGS_YUKAWA_SOURCE_V1,
        includes_fermion_kk_mixing=True,
        charged_lepton_masses_gev=masses,
        c_L=c_l,
        c_E=c_e,
        F_L=f_l,
        F_E=f_e,
        x_e=x_e,
        profile_B_L=profile_b_l,
        profile_B_E=profile_b_e,
        delta_L=delta_l,
        delta_E=delta_e,
        higgs_yukawa_matrix=higgs_yukawa_matrix,
        Y_E_bar_matrix=y_e_bar_matrix,
        U_e_L=u_e_l,
        U_e_R=u_e_r,
        diagnostics={
            "matching_status": "minimal_rs_tree_higgs_lfv_available",
            "v1_honesty_note": RS_HIGGS_YUKAWA_V1_HONESTY_NOTE,
            "charged_lepton_yukawa_is_diagonal": diagonal_yukawa,
            "U_e_L_is_identity": identity_l,
            "U_e_R_is_identity": identity_r,
            "diagonal_v1_tree_lfv_zero": bool(
                diagonal_yukawa and identity_l and identity_r and offdiag_exact_zero
            ),
            "higgs_lfv_offdiag_exact_zero": offdiag_exact_zero,
            "higgs_lfv_offdiag_max_abs": float(np.max(np.abs(offdiag))),
            "generic_rs_higgs_lfv_requires_non_diagonal_charged_lepton_structure": True,
            "includes_exact_fermion_tower_diagonalization": False,
            "brane_kinetic_terms_included": False,
            "custodial_variant_included": False,
            "M_KK_convention": "geometric Lambda_IR = spectrum.lambda_ir_gev",
            "physical_first_gauge_mass_gev": kk_ew_mass,
            "formula": (
                "Y_ij(offdiag)=-[(m_i/v)*(delta_E)_ij+(delta_L)_ij*(m_j/v)], "
                "delta_L=x_e U_e_L^dag diag[B(c_E)] U_e_L x_e, "
                "delta_E=x_e U_e_R^dag diag[B(c_L)] U_e_R x_e"
            ),
        },
    )


def _higgs_yukawa_matrix(
    *,
    charged_lepton_masses_gev: np.ndarray,
    v_gev: float,
    delta_L: np.ndarray,
    delta_E: np.ndarray,
) -> np.ndarray:
    masses = np.asarray(charged_lepton_masses_gev, dtype=float)
    m_over_v = masses / float(v_gev)
    matrix = np.diag(m_over_v).astype(np.complex128)
    for i in range(3):
        for j in range(3):
            matrix[i, j] += -(
                m_over_v[i] * complex(delta_E[i, j])
                + complex(delta_L[i, j]) * m_over_v[j]
            )
    return matrix


def _casagrande_higgs_B_profile_triplet(
    c_values: np.ndarray,
    F_values: np.ndarray,
    *,
    name: str,
) -> np.ndarray:
    values = np.array(
        [
            _casagrande_higgs_B_profile(float(c), float(f), name=f"{name}[{idx}]")
            for idx, (c, f) in enumerate(zip(c_values, F_values, strict=True))
        ],
        dtype=float,
    )
    values.setflags(write=False)
    return values


def _casagrande_higgs_B_profile(c: float, F: float, *, name: str) -> float:
    """CGHNP Higgs-Yukawa fermion-KK bracket in repo variables.

    This is the same diagonal CGHNP bracket used for the audited B1 Zbb fix:
    ``c_CGHNP = -c_repo`` and ``F_CGHNP^2 = 2 f_IR,repo^2``.  Keeping the
    Higgs and Zbb paths on the same helper prevents the old raw-CGHNP
    ``1/(1 - 2c)`` form from reappearing with repo-convention profiles.
    """

    return casagrande_cghnp_B_profile(c, F, name=name)


def _required_attr(source: Any, name: str) -> Any:
    try:
        return getattr(source, name)
    except AttributeError as exc:
        raise ValueError(f"{name} is required for RS Higgs-Yukawa matching") from exc


def _required_mapping_attr(source: Any, name: str) -> Mapping[str, Any]:
    value = _required_attr(source, name)
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must be a mapping for RS Higgs-Yukawa matching")
    return value


def _positive_float(value: Any, name: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be positive and finite") from exc
    if not math.isfinite(number) or number <= 0.0:
        raise ValueError(f"{name} must be positive and finite")
    return number


def _readonly_real_triplet(values: Any, name: str) -> np.ndarray:
    arr = np.array(values, dtype=float, copy=True)
    if arr.shape != (3,):
        raise ValueError(f"{name} must have shape (3,)")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} contains non-finite values")
    arr.setflags(write=False)
    return arr


def _readonly_complex_matrix(values: Any, name: str) -> np.ndarray:
    arr = np.array(values, dtype=np.complex128, copy=True)
    if arr.shape != (3, 3):
        raise ValueError(f"{name} must have shape (3, 3)")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} contains non-finite values")
    arr.setflags(write=False)
    return arr


def _validate_unitary(matrix: np.ndarray, name: str) -> None:
    identity = np.eye(3, dtype=np.complex128)
    if not np.allclose(matrix.conjugate().T @ matrix, identity, rtol=0.0, atol=1.0e-8):
        raise ValueError(f"{name} must be unitary")


def _hermitian(matrix: np.ndarray) -> np.ndarray:
    arr = np.asarray(matrix, dtype=np.complex128)
    return 0.5 * (arr + arr.conjugate().T)


def _offdiag(matrix: np.ndarray) -> np.ndarray:
    arr = np.array(matrix, dtype=np.complex128, copy=True)
    np.fill_diagonal(arr, 0.0)
    return arr


def _is_diagonal(matrix: np.ndarray) -> bool:
    return bool(np.count_nonzero(_offdiag(matrix)) == 0)


def _is_identity(matrix: np.ndarray) -> bool:
    return bool(np.array_equal(matrix, np.eye(3, dtype=np.complex128)))


__all__ = [
    "RS_HIGGS_YUKAWA_INPUT_BUNDLE_V1",
    "RS_HIGGS_YUKAWA_MATCHING_ASSUMPTION_V1",
    "RS_HIGGS_YUKAWA_MODEL_V1",
    "RS_HIGGS_YUKAWA_SOURCE_V1",
    "RS_HIGGS_YUKAWA_V1_HONESTY_NOTE",
    "RSHiggsYukawaCouplings",
    "build_rs_higgs_yukawas",
]
