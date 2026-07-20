"""Mass-basis KK-gluon couplings for the quark-sector MFV module.

This layer rotates the non-universal zero-mode overlap profiles into the exact
quark mass basis returned by :mod:`quarkConstraints.fit`. The resulting
Hermitian matrices are the building blocks for downstream ``Delta F = 2``
matching.

The repo default remains the bookkeeping convention ``M_KK = Lambda_IR``.
Callers that want a different physical KK-mass convention should pass an
explicit ``xi_KK`` or ``M_KK``.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import log, pi, sqrt

import numpy as np

from qcd import alpha_s

from .fit import QuarkFitResult
from .scales import (
    DEFAULT_QUARK_XI_KK,
    KKGluonMassConventionError,
    assert_kk_gluon_mass_convention,
)


COUPLING_POLICY_PERTURBATIVE_4D_LEGACY = "perturbative_4d_legacy"
COUPLING_POLICY_FIXED_GSSTAR_3 = "fixed_gsstar_3"
COUPLING_POLICY_RS_VOLUME_SQRT2L_PHYSICAL = "rs_volume_sqrt2L_physical"

OPERATOR_CONVENTION_PERTURBATIVE_4D_LEGACY = (
    "kk_gluon_tree.perturbative_4d_legacy.v1"
)
OPERATOR_CONVENTION_FIXED_GSSTAR_3 = "kk_gluon_tree.fixed_gsstar_3.v1"
OPERATOR_CONVENTION_RS_VOLUME_SQRT2L_PHYSICAL = (
    "kk_gluon_tree.rs_volume_sqrt2L_physical.v1"
)

_OPERATOR_CONVENTION_BY_POLICY = {
    COUPLING_POLICY_PERTURBATIVE_4D_LEGACY: OPERATOR_CONVENTION_PERTURBATIVE_4D_LEGACY,
    COUPLING_POLICY_FIXED_GSSTAR_3: OPERATOR_CONVENTION_FIXED_GSSTAR_3,
    COUPLING_POLICY_RS_VOLUME_SQRT2L_PHYSICAL: (
        OPERATOR_CONVENTION_RS_VOLUME_SQRT2L_PHYSICAL
    ),
}


@dataclass(frozen=True)
class KKGluonCouplingPolicy:
    """Resolved KK-gluon coupling normalization."""

    coupling_policy_id: str
    operator_convention_id: str
    g_s_4d: float
    g_eff: float
    g_s_multiplier: float


def operator_convention_id_for_coupling_policy(coupling_policy_id: str) -> str:
    """Return the Delta-F=2 operator convention for a coupling policy."""

    try:
        return _OPERATOR_CONVENTION_BY_POLICY[coupling_policy_id]
    except KeyError as exc:
        raise ValueError(f"unknown coupling_policy_id: {coupling_policy_id!r}") from exc


def _resolve_coupling_policy_id(
    *,
    coupling_policy_id: str | None,
    g_s_star: float | None,
) -> str:
    if coupling_policy_id is not None:
        operator_convention_id_for_coupling_policy(coupling_policy_id)
        return coupling_policy_id
    if g_s_star is None:
        return COUPLING_POLICY_PERTURBATIVE_4D_LEGACY
    if abs(float(g_s_star) - 3.0) <= 1e-12:
        return COUPLING_POLICY_FIXED_GSSTAR_3
    raise ValueError(
        "g_s_star values other than 3.0 must be represented by an explicit "
        "named coupling policy; supported policies are "
        f"{sorted(_OPERATOR_CONVENTION_BY_POLICY)}"
    )


def resolve_kk_gluon_coupling_policy(
    *,
    g_s_4d: float,
    lambda_ir_gev: float,
    k_gev: float,
    g_s_star: float | None = None,
    coupling_policy_id: str | None = None,
) -> KKGluonCouplingPolicy:
    """Resolve one of the explicit KK-gluon coupling policies."""

    resolved_policy_id = _resolve_coupling_policy_id(
        coupling_policy_id=coupling_policy_id,
        g_s_star=g_s_star,
    )
    g_s_4d = float(g_s_4d)
    if g_s_4d <= 0.0:
        raise ValueError("g_s_4d must be positive")
    if resolved_policy_id == COUPLING_POLICY_PERTURBATIVE_4D_LEGACY:
        if g_s_star is not None:
            raise ValueError("perturbative_4d_legacy requires g_s_star=None")
        g_eff = g_s_4d
        multiplier = 1.0
    elif resolved_policy_id == COUPLING_POLICY_FIXED_GSSTAR_3:
        if g_s_star is not None and abs(float(g_s_star) - 3.0) > 1e-12:
            raise ValueError("fixed_gsstar_3 requires g_s_star=3.0")
        g_eff = 3.0
        multiplier = float(g_eff / g_s_4d)
    elif resolved_policy_id == COUPLING_POLICY_RS_VOLUME_SQRT2L_PHYSICAL:
        if g_s_star is not None:
            raise ValueError(
                "rs_volume_sqrt2L_physical derives g_eff from g_s_4d*sqrt(2L); "
                "do not also pass g_s_star"
            )
        if lambda_ir_gev <= 0.0 or k_gev <= lambda_ir_gev:
            raise ValueError("physical RS coupling requires k_gev > lambda_ir_gev > 0")
        multiplier = float(sqrt(2.0 * log(float(k_gev) / float(lambda_ir_gev))))
        g_eff = float(g_s_4d * multiplier)
    else:  # pragma: no cover - guarded by _resolve_coupling_policy_id
        raise ValueError(f"unknown coupling_policy_id: {resolved_policy_id!r}")
    return KKGluonCouplingPolicy(
        coupling_policy_id=resolved_policy_id,
        operator_convention_id=operator_convention_id_for_coupling_policy(resolved_policy_id),
        g_s_4d=float(g_s_4d),
        g_eff=float(g_eff),
        g_s_multiplier=float(multiplier),
    )


def _default_quark_m_kk_from_lambda_ir(Lambda_IR: float, xi_KK: float = 1.0) -> float:
    """Resolve the KK scale using the public helper when it is available."""
    try:
        from .scales import default_quark_m_kk_from_lambda_ir

        return default_quark_m_kk_from_lambda_ir(Lambda_IR, xi_KK=xi_KK)
    except Exception as err:
        if Lambda_IR <= 0.0:
            raise ValueError("Lambda_IR must be positive") from err
        if xi_KK <= 0.0:
            raise ValueError("xi_KK must be positive") from err
        return float(xi_KK * Lambda_IR)


def _hermitian(matrix: np.ndarray) -> np.ndarray:
    """Project small numerical noise back onto the Hermitian subspace."""
    arr = np.asarray(matrix, dtype=np.complex128)
    return 0.5 * (arr + arr.conjugate().T)


def _mass_basis_overlap(rotation: np.ndarray, profile_values: np.ndarray) -> np.ndarray:
    """Rotate a diagonal overlap profile into the corresponding mass basis."""
    profile = np.diag(np.asarray(profile_values, dtype=float) ** 2)
    return _hermitian(rotation.conjugate().T @ profile @ rotation)


def _off_diagonal_norm(matrix: np.ndarray) -> float:
    """Return the Frobenius norm of off-diagonal entries only."""
    arr = np.asarray(matrix, dtype=np.complex128).copy()
    np.fill_diagonal(arr, 0.0)
    return float(np.linalg.norm(arr, ord="fro"))


@dataclass(frozen=True)
class QuarkMassBasisCouplings:
    """KK-gluon couplings in the exact quark mass basis.

    The ``*_overlap`` matrices are the dimensionless overlap structures. The
    ``left_*`` / ``right_*`` matrices multiply them by the running QCD gauge
    coupling ``g_s(M_KK)``.
    """

    M_KK: float
    xi_KK: float
    alpha_s: float
    g_s: float
    left_overlap: np.ndarray
    right_up_overlap: np.ndarray
    right_down_overlap: np.ndarray
    left_up: np.ndarray
    left_down: np.ndarray
    right_up: np.ndarray
    right_down: np.ndarray
    lambda_ir_gev: float | None = None
    m_kk_physical_gev: float | None = None
    mass_convention_id: str | None = None
    g_s_4d: float | None = None
    g_eff: float | None = None
    g_s_multiplier: float | None = None
    coupling_policy_id: str = COUPLING_POLICY_PERTURBATIVE_4D_LEGACY
    operator_convention_id: str | None = None
    ckm_matrix: np.ndarray | None = None
    ckm_source: str | None = None

    def __post_init__(self) -> None:
        lambda_ir = (
            float(self.M_KK) / float(self.xi_KK)
            if self.lambda_ir_gev is None
            else float(self.lambda_ir_gev)
        )
        m_kk_physical = (
            float(self.M_KK)
            if self.m_kk_physical_gev is None
            else float(self.m_kk_physical_gev)
        )
        scale = assert_kk_gluon_mass_convention(
            lambda_ir_gev=lambda_ir,
            m_kk_physical_gev=m_kk_physical,
            xi_kk=float(self.xi_KK),
            mass_convention_id=self.mass_convention_id,
            context="QuarkMassBasisCouplings",
        )
        g_s_4d = float(self.g_s if self.g_s_4d is None else self.g_s_4d)
        g_eff = float(self.g_s if self.g_eff is None else self.g_eff)
        multiplier = float(
            g_eff / g_s_4d if self.g_s_multiplier is None else self.g_s_multiplier
        )
        operator_convention_id = (
            operator_convention_id_for_coupling_policy(self.coupling_policy_id)
            if self.operator_convention_id is None
            else str(self.operator_convention_id)
        )
        expected_operator = operator_convention_id_for_coupling_policy(
            self.coupling_policy_id
        )
        if operator_convention_id != expected_operator:
            raise ValueError(
                "operator_convention_id must match coupling_policy_id: "
                f"{operator_convention_id!r} != {expected_operator!r}"
            )
        ckm_matrix = self.ckm_matrix
        ckm_source = self.ckm_source
        if ckm_matrix is None:
            if ckm_source is not None:
                raise ValueError("ckm_source requires ckm_matrix")
        else:
            ckm_matrix = np.asarray(ckm_matrix, dtype=np.complex128).copy()
            if ckm_matrix.shape != (3, 3):
                raise ValueError("ckm_matrix must have shape (3, 3)")
            if not np.all(np.isfinite(ckm_matrix)):
                raise ValueError("ckm_matrix must contain finite entries")
            if not np.allclose(
                ckm_matrix.conjugate().T @ ckm_matrix,
                np.eye(3, dtype=np.complex128),
                rtol=0.0,
                atol=1e-10,
            ):
                raise ValueError("ckm_matrix must be unitary")
            if ckm_source is None:
                ckm_source = "quark_mass_basis_couplings.ckm_matrix"
            elif not str(ckm_source).strip():
                raise ValueError("ckm_source must be non-empty")
            ckm_matrix.setflags(write=False)
        object.__setattr__(self, "M_KK", float(scale.m_kk_physical_gev))
        object.__setattr__(self, "xi_KK", float(scale.xi_kk))
        object.__setattr__(self, "lambda_ir_gev", float(scale.lambda_ir_gev))
        object.__setattr__(self, "m_kk_physical_gev", float(scale.m_kk_physical_gev))
        object.__setattr__(self, "mass_convention_id", scale.mass_convention_id)
        object.__setattr__(self, "g_s_4d", g_s_4d)
        object.__setattr__(self, "g_eff", g_eff)
        object.__setattr__(self, "g_s_multiplier", multiplier)
        object.__setattr__(self, "operator_convention_id", operator_convention_id)
        object.__setattr__(self, "ckm_matrix", ckm_matrix)
        object.__setattr__(self, "ckm_source", ckm_source)

    @property
    def left_down_offdiag_norm(self) -> float:
        return _off_diagonal_norm(self.left_down)

    @property
    def left_up_offdiag_norm(self) -> float:
        return _off_diagonal_norm(self.left_up)

    @property
    def right_down_offdiag_norm(self) -> float:
        return _off_diagonal_norm(self.right_down)

    @property
    def right_up_offdiag_norm(self) -> float:
        return _off_diagonal_norm(self.right_up)


def compute_quark_kk_gluon_couplings(
    fit_result: QuarkFitResult,
    *,
    M_KK: float | None = None,
    m_kk_physical_gev: float | None = None,
    lambda_ir_gev: float | None = None,
    xi_KK: float | None = None,
    g_s_star: float | None = None,
    coupling_policy_id: str | None = None,
) -> QuarkMassBasisCouplings:
    """Return mass-basis KK-gluon couplings for a fitted MFV point.

    Parameters
    ----------
    fit_result
        Exact quark fit result from :func:`quarkConstraints.fit.fit_quark_sector`.
    M_KK
        Explicit KK scale in GeV. When omitted, the helper uses the repo's
        default quark convention ``M_KK = xi_KK * Lambda_IR``.
    m_kk_physical_gev
        Explicit physical first KK-gluon mass in GeV. This is the typed
        replacement for the legacy ``M_KK`` alias; if both are supplied they
        must agree.
    lambda_ir_gev
        Explicit geometric IR scale in GeV. Defaults to the fitted point's
        ``Lambda_IR``.
    xi_KK
        Conversion factor between the geometric IR scale and the KK scale when
        ``M_KK`` is not supplied. The repo default is ``1.0``.
    g_s_star
        Enhanced 5D KK-gluon coupling. In RS models the first KK-gluon mode
        has a wavefunction peaked near the IR brane, enhancing its coupling to
        zero-mode quarks by a factor ``sqrt(2 k pi r_c)`` relative to the SM
        ``g_s``. Typical values are ``~6`` at tree level (Csaki, Falkowski,
        Weiler 2008) or ``~3`` at one-loop level (Gedalia, Grossman, Nir, Perez
        2009). When provided, this value is used directly as the gauge coupling
        multiplying the overlap matrices instead of the perturbative
        ``g_s(M_KK)``. When ``None`` (default), falls back to the 4D
        perturbative coupling ``g_s = sqrt(4 pi alpha_s(M_KK))``.
    coupling_policy_id
        Named coupling policy. Supported values are
        ``perturbative_4d_legacy``, ``fixed_gsstar_3``, and
        ``rs_volume_sqrt2L_physical``.
    """
    if xi_KK is not None and xi_KK <= 0.0:
        raise ValueError("xi_KK must be positive")
    if g_s_star is not None and float(g_s_star) <= 0.0:
        raise ValueError("g_s_star must be positive when provided")

    state = fit_result.state
    resolved_lambda_ir = (
        float(state.point.Lambda_IR) if lambda_ir_gev is None else float(lambda_ir_gev)
    )
    if M_KK is not None and m_kk_physical_gev is not None:
        if abs(float(M_KK) - float(m_kk_physical_gev)) > 1e-9:
            raise KKGluonMassConventionError(
                "compute_quark_kk_gluon_couplings: M_KK and "
                "m_kk_physical_gev must agree"
            )
    supplied_mkk = m_kk_physical_gev if m_kk_physical_gev is not None else M_KK
    resolved_xi_kk = (
        float(supplied_mkk) / resolved_lambda_ir
        if xi_KK is None and supplied_mkk is not None
        else float(DEFAULT_QUARK_XI_KK if xi_KK is None else xi_KK)
    )
    resolved_mkk = (
        float(supplied_mkk)
        if supplied_mkk is not None
        else _default_quark_m_kk_from_lambda_ir(resolved_lambda_ir, xi_KK=resolved_xi_kk)
    )
    scale = assert_kk_gluon_mass_convention(
        lambda_ir_gev=resolved_lambda_ir,
        m_kk_physical_gev=resolved_mkk,
        xi_kk=resolved_xi_kk,
        context="compute_quark_kk_gluon_couplings",
    )

    running_alpha_s = float(alpha_s(scale.m_kk_physical_gev, precision="high"))
    perturbative_g_s = float(sqrt(4.0 * pi * running_alpha_s))
    coupling_policy = resolve_kk_gluon_coupling_policy(
        g_s_4d=perturbative_g_s,
        lambda_ir_gev=scale.lambda_ir_gev,
        k_gev=state.point.k,
        g_s_star=g_s_star,
        coupling_policy_id=coupling_policy_id,
    )
    effective_g_s = coupling_policy.g_eff

    left_overlap = _mass_basis_overlap(fit_result.U_L_d, state.F_Q)
    right_down_overlap = _mass_basis_overlap(fit_result.U_R_d, state.F_d)
    right_up_overlap = _mass_basis_overlap(fit_result.U_R_u, state.F_u)
    left_up_overlap = _mass_basis_overlap(fit_result.U_L_u, state.F_Q)

    return QuarkMassBasisCouplings(
        M_KK=scale.m_kk_physical_gev,
        xi_KK=scale.xi_kk,
        alpha_s=running_alpha_s,
        g_s=effective_g_s,
        left_overlap=left_overlap,
        right_up_overlap=right_up_overlap,
        right_down_overlap=right_down_overlap,
        left_up=_hermitian(effective_g_s * left_up_overlap),
        left_down=_hermitian(effective_g_s * left_overlap),
        right_up=_hermitian(effective_g_s * right_up_overlap),
        right_down=_hermitian(effective_g_s * right_down_overlap),
        lambda_ir_gev=scale.lambda_ir_gev,
        m_kk_physical_gev=scale.m_kk_physical_gev,
        mass_convention_id=scale.mass_convention_id,
        g_s_4d=coupling_policy.g_s_4d,
        g_eff=coupling_policy.g_eff,
        g_s_multiplier=coupling_policy.g_s_multiplier,
        coupling_policy_id=coupling_policy.coupling_policy_id,
        operator_convention_id=coupling_policy.operator_convention_id,
        ckm_matrix=fit_result.ckm_matrix,
        ckm_source="QuarkFitResult.ckm_matrix (U_L_u^dagger U_L_d)",
    )
