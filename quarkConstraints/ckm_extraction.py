"""CKM extractions from source-level semileptonic inputs.

This module owns small, model-independent arithmetic used by catalog
constraints that need to turn measured products into CKM magnitudes.  It does
not implement RS charged-current matching.  For K_l3,

    |V_us| = (|V_us| f_+(0)) / f_+(0),

with uncorrelated product, lattice, and optional radiative/isospin uncertainty
components propagated in quadrature.
"""

from __future__ import annotations

import cmath
import math
from dataclasses import dataclass
from typing import Mapping

import numpy as np

__all__ = [
    "NeutralBMixingSMAmplitude",
    "CKMPhaseAngles",
    "KL3VusExtraction",
    "VusConsistencyResult",
    "ckm_phases_from_matrix",
    "neutral_b_mixing_sm_amplitude",
    "extract_vus_from_kl3",
    "repo_default_ckm_matrix",
    "repo_default_ckm_phases",
    "vus_consistency_pull",
]


@dataclass(frozen=True)
class KL3VusExtraction:
    """Extracted ``|V_us|`` and propagated uncertainty components."""

    fplus_vus_product: float
    fplus_vus_uncertainty: float
    fplus_zero: float
    fplus_zero_uncertainty: float
    vus: float
    experimental_uncertainty: float
    lattice_uncertainty: float
    radiative_isospin_uncertainty: float
    total_uncertainty: float

    @property
    def uncertainty_components(self) -> Mapping[str, float]:
        return {
            "experimental": float(self.experimental_uncertainty),
            "lattice": float(self.lattice_uncertainty),
            "radiative_isospin": float(self.radiative_isospin_uncertainty),
        }


@dataclass(frozen=True)
class VusConsistencyResult:
    """One-sigma consistency pull for two ``|V_us|`` determinations."""

    extracted_vus: float
    reference_vus: float
    delta_vus: float
    budget: float
    pull_sigma: float
    passes: bool


@dataclass(frozen=True)
class CKMPhaseAngles:
    """Rephasing-invariant CKM angles relevant for neutral B mixing."""

    source: str
    beta: float
    beta_degrees: float
    two_beta: float
    two_beta_degrees: float
    sin_2beta: float
    beta_s: float
    beta_s_degrees: float
    phi_s: float
    phi_s_degrees: float


@dataclass(frozen=True)
class NeutralBMixingSMAmplitude:
    """Complex SM dispersive amplitude for a neutral ``B_q`` system.

    The magnitude is fixed by ``Delta m_q^SM / 2`` and the weak phase by the
    top-box CKM factor ``(V_tq* V_tb)^2``.  Keeping this phase is essential:
    only the ratio of a new-physics amplitude to this complex SM amplitude is
    invariant under quark mass-eigenstate rephasings.
    """

    system: str
    light_down_index: int
    delta_m_sm_gev: float
    magnitude_gev: float
    ckm_factor: complex
    phase_factor: complex
    amplitude_gev: complex
    ckm_source: str


def _finite_float(name: str, value: object) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name}={value!r} is not numeric") from exc
    if not math.isfinite(number):
        raise ValueError(f"{name} must be finite")
    return number


def _positive_float(name: str, value: object) -> float:
    number = _finite_float(name, value)
    if number <= 0.0:
        raise ValueError(f"{name} must be positive")
    return number


def _nonnegative_float(name: str, value: object) -> float:
    number = _finite_float(name, value)
    if number < 0.0:
        raise ValueError(f"{name} must be nonnegative")
    return number


def _ckm_matrix(name: str, values: object) -> np.ndarray:
    matrix = np.asarray(values, dtype=np.complex128)
    if matrix.shape != (3, 3):
        raise ValueError(f"{name} must have shape (3, 3)")
    if not np.all(np.isfinite(matrix.real)) or not np.all(np.isfinite(matrix.imag)):
        raise ValueError(f"{name} must contain finite entries")
    return matrix


def _nonzero_complex(name: str, value: complex) -> complex:
    number = complex(value)
    if not math.isfinite(number.real) or not math.isfinite(number.imag):
        raise ValueError(f"{name} must be finite")
    if abs(number) <= 0.0:
        raise ValueError(f"{name} must be non-zero")
    return number


def ckm_phases_from_matrix(
    ckm: object,
    *,
    source: str = "ckm_matrix_input",
) -> CKMPhaseAngles:
    """Compute ``beta``, ``sin(2 beta)``, and ``phi_s = -2 beta_s``.

    The phase conventions are the standard rephasing-invariant definitions

    ``beta = arg(-(V_cd V_cb*) / (V_td V_tb*))``

    and

    ``beta_s = arg(-(V_ts V_tb*) / (V_cs V_cb*))``.
    """

    matrix = _ckm_matrix("ckm", ckm)
    beta_denominator = _nonzero_complex(
        "V_td V_tb*",
        matrix[2, 0] * np.conjugate(matrix[2, 2]),
    )
    beta_s_denominator = _nonzero_complex(
        "V_cs V_cb*",
        matrix[1, 1] * np.conjugate(matrix[1, 2]),
    )
    beta_ratio = -(
        matrix[1, 0] * np.conjugate(matrix[1, 2])
    ) / beta_denominator
    beta_s_ratio = -(
        matrix[2, 1] * np.conjugate(matrix[2, 2])
    ) / beta_s_denominator
    beta = float(cmath.phase(beta_ratio))
    beta_s = float(cmath.phase(beta_s_ratio))
    two_beta = 2.0 * beta
    phi_s = -2.0 * beta_s
    return CKMPhaseAngles(
        source=str(source),
        beta=beta,
        beta_degrees=math.degrees(beta),
        two_beta=two_beta,
        two_beta_degrees=math.degrees(two_beta),
        sin_2beta=float(math.sin(two_beta)),
        beta_s=beta_s,
        beta_s_degrees=math.degrees(beta_s),
        phi_s=phi_s,
        phi_s_degrees=math.degrees(phi_s),
    )


def repo_default_ckm_matrix() -> np.ndarray:
    """Return the repo-owned PDG-2024 CKM target unitary."""

    from quarkConstraints.model import RotationParameters, ckm_like_unitary
    from quarkConstraints.modern.inputs import ModernDefaultCKMTarget

    target = ModernDefaultCKMTarget()
    return ckm_like_unitary(
        RotationParameters(
            theta12=float(target.theta12),
            theta13=float(target.theta13),
            theta23=float(target.theta23),
            delta=float(target.delta),
        )
    )


def repo_default_ckm_phases() -> CKMPhaseAngles:
    """Return B-mixing CKM phases from the repo-owned PDG-2024 target."""

    from quarkConstraints.modern.inputs import ModernDefaultCKMTarget

    target = ModernDefaultCKMTarget()
    return ckm_phases_from_matrix(repo_default_ckm_matrix(), source=target.target_id)


def neutral_b_mixing_sm_amplitude(
    *,
    delta_m_sm_gev: float,
    light_down_index: int,
    ckm: object | None = None,
    ckm_source: str | None = None,
) -> NeutralBMixingSMAmplitude:
    """Return the complex ``M12_SM`` for ``B_d`` or ``B_s`` mixing.

    The short-distance and hadronic inputs already encoded in
    ``delta_m_sm_gev`` determine ``|M12_SM| = Delta m_SM / 2``.  The phase is
    restored from the SM top-box factor

    ``lambda_tq^2 = (conj(V_tq) * V_tb)^2``.

    This convention transforms exactly like Wilson coefficients built from
    the mass-basis ``q-b`` coupling squared.  Consequently
    ``M12_NP / M12_SM`` is unchanged by simultaneous quark-field rephasings.
    See arXiv:0809.1073, Eqs. (4.39)--(4.45), for the neutral-B amplitude and
    CP-phase convention used by the RS analysis.
    """

    delta_m = _positive_float("delta_m_sm_gev", delta_m_sm_gev)
    if light_down_index not in (0, 1):
        raise ValueError("light_down_index must be 0 (d) or 1 (s)")

    if ckm is None:
        matrix = repo_default_ckm_matrix()
        source = repo_default_ckm_phases().source
    else:
        matrix = _ckm_matrix("ckm", ckm)
        source = (
            "quark_mass_basis_couplings.ckm_matrix"
            if ckm_source is None
            else str(ckm_source)
        )
    if not source.strip():
        raise ValueError("ckm_source must be non-empty")

    lambda_tq = _nonzero_complex(
        "conj(V_tq) V_tb",
        np.conjugate(matrix[2, light_down_index]) * matrix[2, 2],
    )
    ckm_factor = lambda_tq * lambda_tq
    phase_factor = ckm_factor / abs(ckm_factor)
    magnitude = delta_m / 2.0
    system = "B_d" if light_down_index == 0 else "B_s"
    return NeutralBMixingSMAmplitude(
        system=system,
        light_down_index=int(light_down_index),
        delta_m_sm_gev=float(delta_m),
        magnitude_gev=float(magnitude),
        ckm_factor=complex(ckm_factor),
        phase_factor=complex(phase_factor),
        amplitude_gev=complex(magnitude * phase_factor),
        ckm_source=source,
    )


def extract_vus_from_kl3(
    *,
    fplus_vus_product: float,
    fplus_vus_uncertainty: float,
    fplus_zero: float,
    fplus_zero_uncertainty: float,
    radiative_isospin_uncertainty: float = 0.0,
) -> KL3VusExtraction:
    """Return ``|V_us|`` from the K_l3 product and lattice ``f_+(0)``."""

    product = _positive_float("fplus_vus_product", fplus_vus_product)
    product_sigma = _positive_float("fplus_vus_uncertainty", fplus_vus_uncertainty)
    fplus = _positive_float("fplus_zero", fplus_zero)
    fplus_sigma = _positive_float("fplus_zero_uncertainty", fplus_zero_uncertainty)
    radiative_sigma = _nonnegative_float(
        "radiative_isospin_uncertainty",
        radiative_isospin_uncertainty,
    )

    vus = product / fplus
    experimental_sigma = product_sigma / fplus
    lattice_sigma = product * fplus_sigma / (fplus * fplus)
    total_sigma = math.sqrt(
        experimental_sigma * experimental_sigma
        + lattice_sigma * lattice_sigma
        + radiative_sigma * radiative_sigma
    )
    return KL3VusExtraction(
        fplus_vus_product=float(product),
        fplus_vus_uncertainty=float(product_sigma),
        fplus_zero=float(fplus),
        fplus_zero_uncertainty=float(fplus_sigma),
        vus=float(vus),
        experimental_uncertainty=float(experimental_sigma),
        lattice_uncertainty=float(lattice_sigma),
        radiative_isospin_uncertainty=float(radiative_sigma),
        total_uncertainty=float(total_sigma),
    )


def vus_consistency_pull(
    *,
    extracted_vus: float,
    reference_vus: float,
    budget: float,
) -> VusConsistencyResult:
    """Compare an extracted ``|V_us|`` with a reference value."""

    extracted = _finite_float("extracted_vus", extracted_vus)
    reference = _finite_float("reference_vus", reference_vus)
    sigma = _positive_float("budget", budget)
    delta = extracted - reference
    pull = abs(delta) / sigma
    return VusConsistencyResult(
        extracted_vus=float(extracted),
        reference_vus=float(reference),
        delta_vus=float(delta),
        budget=float(sigma),
        pull_sigma=float(pull),
        passes=bool(pull <= 1.0),
    )
