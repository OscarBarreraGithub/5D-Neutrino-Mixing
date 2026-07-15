"""Paper-owned coupling contracts for the dedicated ``paper_0710_1869`` mode."""

from __future__ import annotations

import math
from dataclasses import asdict, dataclass

import numpy as np

from qcd import alpha_s
from warpConfig.baseParams import MPL, get_warp_params

from .conventions import PAPER_0710_1869_MODE_ID, PAPER_0710_1869_PAPER_ID
from .scales import Paper07101869ScalePoint, default_paper_0710_1869_scales
from .validation import (
    normalize_optional_positive_finite,
    require_known_schema_id,
    require_member,
    require_nonempty_identifier,
    require_positive_finite,
)

PAPER_0710_1869_COUPLINGS_SCHEMA_ID = "quarkConstraints.paper_0710_1869.couplings.v1"
PAPER_0710_1869_GAUGE_NORMALIZATION_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.gauge_normalization.v1"
)
PAPER_0710_1869_DIMENSIONLESS_MATRIX_POLICY_ID = (
    "dimensionless_flavor_matrices.separate_from_gs.normalization.v1"
)
PAPER_0710_1869_GS_NORMALIZATION_ID = (
    "explicit_mu_gs.g_s_sqrt_4pi_alpha_s.with_rs_volume_sqrt_2L.v2"
)
PAPER_0710_1869_MU_GS_SEMANTICS_ID = "alpha_s.high_precision.msbar.at_mu_gs_GeV.v1"
PAPER_0710_1869_RS_VOLUME_POLICY_ID = (
    "rs_volume.L_equals_log_k_over_lambda_ir.sqrt_2L_kk_gluon.v1"
)
PAPER_0710_1869_PROPAGATOR_MASS_RULE_ID = (
    "propagator_mass.scale_point_propagator_mass_GeV.v1"
)
PAPER_0710_1869_UNIVERSAL_SUBTRACTION_POLICY_ID = (
    "subtract_trace_average_identity.from_dimensionless_flavor_matrix.v1"
)
PAPER_0710_1869_SINGLE_MODE_SCALE_POLICY_ID = "single_mode.first_kk_gluon_only.v1"
PAPER_0710_1869_EFFECTIVE_SCALE_POLICY_ID = "effective_scale.explicit_m_KK_eff_GeV.v1"


@dataclass(frozen=True)
class Paper07101869CouplingContract:
    """Frozen paper-mode coupling contract for KK-gluon normalization semantics."""

    schema_id: str = PAPER_0710_1869_COUPLINGS_SCHEMA_ID
    mode_id: str = PAPER_0710_1869_MODE_ID
    paper_id: str = PAPER_0710_1869_PAPER_ID
    dimensionless_matrix_policy_id: str = PAPER_0710_1869_DIMENSIONLESS_MATRIX_POLICY_ID
    kk_gluon_normalization_id: str = PAPER_0710_1869_GS_NORMALIZATION_ID
    mu_gs_semantics_id: str = PAPER_0710_1869_MU_GS_SEMANTICS_ID
    alpha_s_precision: str = "high"
    rs_volume_policy_id: str = PAPER_0710_1869_RS_VOLUME_POLICY_ID
    propagator_mass_rule_id: str = PAPER_0710_1869_PROPAGATOR_MASS_RULE_ID
    universal_subtraction_policy_id: str = PAPER_0710_1869_UNIVERSAL_SUBTRACTION_POLICY_ID
    single_mode_scale_policy_id: str = PAPER_0710_1869_SINGLE_MODE_SCALE_POLICY_ID
    effective_scale_policy_id: str = PAPER_0710_1869_EFFECTIVE_SCALE_POLICY_ID
    mu_gs_semantics_note: str = (
        "Evaluate alpha_s only at scale_point.mu_gs_GeV with precision='high'. "
        "mu_match_GeV and propagator masses do not change the 4D g_s normalization."
    )
    rs_volume_note: str = (
        "The KK-gluon flavor matrices use the RS volume L = log(k/Lambda_IR) = "
        "pi*k*r_c from the repository warp geometry. Full matrices are normalized "
        "as g_s*(sqrt(2L)*raw - I/sqrt(2L)); FCNC-subtracted matrices are "
        "g_s*sqrt(2L)*subtracted."
    )
    propagator_mass_rule_note: str = (
        "Heavy propagator denominators use scale_point.propagator_mass_GeV. "
        "That equals m_g1_GeV unless an explicit m_KK_eff_GeV override is present."
    )
    universal_subtraction_note: str = (
        "The universal piece is the trace-average identity contribution in each "
        "dimensionless flavor matrix. Both raw and subtracted matrices are retained."
    )
    scale_policy_note: str = (
        "Without m_KK_eff_GeV the contract is first-mode only. "
        "Supplying m_KK_eff_GeV opts into an explicit effective-scale policy."
    )

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "schema_id",
            require_known_schema_id(
                "schema_id",
                self.schema_id,
                expected=PAPER_0710_1869_COUPLINGS_SCHEMA_ID,
            ),
        )
        object.__setattr__(self, "mode_id", require_nonempty_identifier("mode_id", self.mode_id))
        object.__setattr__(self, "paper_id", require_nonempty_identifier("paper_id", self.paper_id))
        require_member("mode_id", self.mode_id, (PAPER_0710_1869_MODE_ID,))
        require_member("paper_id", self.paper_id, (PAPER_0710_1869_PAPER_ID,))
        require_member("alpha_s_precision", self.alpha_s_precision, ("high",))
        for field_name in (
            "dimensionless_matrix_policy_id",
            "kk_gluon_normalization_id",
            "mu_gs_semantics_id",
            "rs_volume_policy_id",
            "propagator_mass_rule_id",
            "universal_subtraction_policy_id",
            "single_mode_scale_policy_id",
            "effective_scale_policy_id",
            "mu_gs_semantics_note",
            "rs_volume_note",
            "propagator_mass_rule_note",
            "universal_subtraction_note",
            "scale_policy_note",
        ):
            object.__setattr__(
                self,
                field_name,
                require_nonempty_identifier(field_name, getattr(self, field_name)),
            )

    def scale_policy_id_for(self, scale_point: Paper07101869ScalePoint) -> str:
        """Return the frozen scale-policy identifier for one explicit scale point."""
        if not isinstance(scale_point, Paper07101869ScalePoint):
            raise ValueError("scale_point must be a Paper07101869ScalePoint")
        if scale_point.has_explicit_effective_kk_scale:
            return self.effective_scale_policy_id
        return self.single_mode_scale_policy_id

    def as_dict(self) -> dict[str, str]:
        """Return a stable mapping representation."""
        return asdict(self)


@dataclass(frozen=True)
class Paper07101869GaugeCouplingNormalization:
    """Resolved QCD normalization bundle for one paper-mode scale point."""

    scale_label: str
    Lambda_IR_GeV: float
    mu_gs_GeV: float
    mu_match_GeV: float
    m_g1_GeV: float
    propagator_mass_GeV: float
    alpha_s_mu_gs: float
    g_s_mu_gs: float
    rs_volume_L: float
    rs_volume_sqrt_2L: float
    g_s_star_mu_gs: float
    schema_id: str = PAPER_0710_1869_GAUGE_NORMALIZATION_SCHEMA_ID
    contract_schema_id: str = PAPER_0710_1869_COUPLINGS_SCHEMA_ID
    mode_id: str = PAPER_0710_1869_MODE_ID
    paper_id: str = PAPER_0710_1869_PAPER_ID
    kk_gluon_normalization_id: str = PAPER_0710_1869_GS_NORMALIZATION_ID
    mu_gs_semantics_id: str = PAPER_0710_1869_MU_GS_SEMANTICS_ID
    alpha_s_precision: str = "high"
    rs_volume_policy_id: str = PAPER_0710_1869_RS_VOLUME_POLICY_ID
    scale_policy_id: str = PAPER_0710_1869_SINGLE_MODE_SCALE_POLICY_ID
    propagator_mass_rule_id: str = PAPER_0710_1869_PROPAGATOR_MASS_RULE_ID
    m_KK_eff_GeV: float | None = None

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "schema_id",
            require_known_schema_id(
                "schema_id",
                self.schema_id,
                expected=PAPER_0710_1869_GAUGE_NORMALIZATION_SCHEMA_ID,
            ),
        )
        object.__setattr__(
            self,
            "contract_schema_id",
            require_known_schema_id(
                "contract_schema_id",
                self.contract_schema_id,
                expected=PAPER_0710_1869_COUPLINGS_SCHEMA_ID,
            ),
        )
        object.__setattr__(
            self,
            "scale_label",
            require_nonempty_identifier("scale_label", self.scale_label),
        )
        object.__setattr__(self, "mode_id", require_nonempty_identifier("mode_id", self.mode_id))
        object.__setattr__(self, "paper_id", require_nonempty_identifier("paper_id", self.paper_id))
        require_member("mode_id", self.mode_id, (PAPER_0710_1869_MODE_ID,))
        require_member("paper_id", self.paper_id, (PAPER_0710_1869_PAPER_ID,))
        require_member("alpha_s_precision", self.alpha_s_precision, ("high",))
        require_member(
            "scale_policy_id",
            self.scale_policy_id,
            (
                PAPER_0710_1869_SINGLE_MODE_SCALE_POLICY_ID,
                PAPER_0710_1869_EFFECTIVE_SCALE_POLICY_ID,
            ),
        )
        for field_name in (
            "kk_gluon_normalization_id",
            "mu_gs_semantics_id",
            "rs_volume_policy_id",
            "propagator_mass_rule_id",
        ):
            object.__setattr__(
                self,
                field_name,
                require_nonempty_identifier(field_name, getattr(self, field_name)),
            )

        object.__setattr__(
            self,
            "Lambda_IR_GeV",
            require_positive_finite("Lambda_IR_GeV", self.Lambda_IR_GeV),
        )
        object.__setattr__(self, "mu_gs_GeV", require_positive_finite("mu_gs_GeV", self.mu_gs_GeV))
        object.__setattr__(
            self,
            "mu_match_GeV",
            require_positive_finite("mu_match_GeV", self.mu_match_GeV),
        )
        object.__setattr__(self, "m_g1_GeV", require_positive_finite("m_g1_GeV", self.m_g1_GeV))
        object.__setattr__(
            self,
            "propagator_mass_GeV",
            require_positive_finite("propagator_mass_GeV", self.propagator_mass_GeV),
        )
        object.__setattr__(
            self,
            "m_KK_eff_GeV",
            normalize_optional_positive_finite("m_KK_eff_GeV", self.m_KK_eff_GeV),
        )
        object.__setattr__(
            self,
            "alpha_s_mu_gs",
            require_positive_finite("alpha_s_mu_gs", self.alpha_s_mu_gs),
        )
        object.__setattr__(
            self,
            "g_s_mu_gs",
            require_positive_finite("g_s_mu_gs", self.g_s_mu_gs),
        )
        object.__setattr__(
            self,
            "rs_volume_L",
            require_positive_finite("rs_volume_L", self.rs_volume_L),
        )
        object.__setattr__(
            self,
            "rs_volume_sqrt_2L",
            require_positive_finite("rs_volume_sqrt_2L", self.rs_volume_sqrt_2L),
        )
        object.__setattr__(
            self,
            "g_s_star_mu_gs",
            require_positive_finite("g_s_star_mu_gs", self.g_s_star_mu_gs),
        )

        expected_scale_policy = (
            PAPER_0710_1869_EFFECTIVE_SCALE_POLICY_ID
            if self.m_KK_eff_GeV is not None
            else PAPER_0710_1869_SINGLE_MODE_SCALE_POLICY_ID
        )
        if self.scale_policy_id != expected_scale_policy:
            raise ValueError(
                "scale_policy_id is inconsistent with the presence of m_KK_eff_GeV"
            )

        expected_propagator_mass = (
            self.m_g1_GeV if self.m_KK_eff_GeV is None else self.m_KK_eff_GeV
        )
        if not math.isclose(
            self.propagator_mass_GeV,
            expected_propagator_mass,
            rel_tol=1.0e-12,
            abs_tol=0.0,
        ):
            raise ValueError(
                "propagator_mass_GeV must equal m_g1_GeV unless m_KK_eff_GeV is set"
            )

        expected_gs = math.sqrt(4.0 * math.pi * self.alpha_s_mu_gs)
        if not math.isclose(self.g_s_mu_gs, expected_gs, rel_tol=1.0e-12, abs_tol=0.0):
            raise ValueError("g_s_mu_gs must equal sqrt(4*pi*alpha_s_mu_gs)")
        expected_sqrt_2L = math.sqrt(2.0 * self.rs_volume_L)
        if not math.isclose(
            self.rs_volume_sqrt_2L,
            expected_sqrt_2L,
            rel_tol=1.0e-12,
            abs_tol=0.0,
        ):
            raise ValueError("rs_volume_sqrt_2L must equal sqrt(2*rs_volume_L)")
        expected_gs_star = self.g_s_mu_gs * self.rs_volume_sqrt_2L
        if not math.isclose(
            self.g_s_star_mu_gs,
            expected_gs_star,
            rel_tol=1.0e-12,
            abs_tol=0.0,
        ):
            raise ValueError("g_s_star_mu_gs must equal g_s_mu_gs * rs_volume_sqrt_2L")

    def as_dict(self) -> dict[str, float | str | None]:
        """Return a stable mapping representation."""
        return asdict(self)


def default_paper_0710_1869_coupling_contract() -> Paper07101869CouplingContract:
    """Return the frozen default paper-mode coupling contract."""
    return Paper07101869CouplingContract()


def evaluate_paper_0710_1869_gauge_coupling(
    scale_point: Paper07101869ScalePoint | None = None,
    *,
    contract: Paper07101869CouplingContract | None = None,
) -> Paper07101869GaugeCouplingNormalization:
    """Resolve the QCD KK-gluon normalization at ``mu_gs_GeV`` only."""
    resolved_scale = default_paper_0710_1869_scales() if scale_point is None else scale_point
    if not isinstance(resolved_scale, Paper07101869ScalePoint):
        raise ValueError("scale_point must be a Paper07101869ScalePoint")

    resolved_contract = (
        default_paper_0710_1869_coupling_contract() if contract is None else contract
    )
    if not isinstance(resolved_contract, Paper07101869CouplingContract):
        raise ValueError("contract must be a Paper07101869CouplingContract")

    resolved_alpha_s = require_positive_finite(
        "alpha_s(mu_gs_GeV)",
        alpha_s(resolved_scale.mu_gs_GeV, precision=resolved_contract.alpha_s_precision),
    )
    resolved_g_s = require_positive_finite(
        "g_s(mu_gs_GeV)",
        math.sqrt(4.0 * math.pi * resolved_alpha_s),
    )
    rs_volume_L = require_positive_finite(
        "rs_volume_L",
        get_warp_params(k=MPL, Lambda_IR=resolved_scale.Lambda_IR_GeV)["warp_log"],
    )
    rs_volume_sqrt_2L = require_positive_finite(
        "rs_volume_sqrt_2L",
        math.sqrt(2.0 * rs_volume_L),
    )
    return Paper07101869GaugeCouplingNormalization(
        scale_label=resolved_scale.label,
        Lambda_IR_GeV=resolved_scale.Lambda_IR_GeV,
        mu_gs_GeV=resolved_scale.mu_gs_GeV,
        mu_match_GeV=resolved_scale.mu_match_GeV,
        m_g1_GeV=resolved_scale.m_g1_GeV,
        propagator_mass_GeV=resolved_scale.propagator_mass_GeV,
        alpha_s_mu_gs=resolved_alpha_s,
        g_s_mu_gs=resolved_g_s,
        rs_volume_L=rs_volume_L,
        rs_volume_sqrt_2L=rs_volume_sqrt_2L,
        g_s_star_mu_gs=resolved_g_s * rs_volume_sqrt_2L,
        kk_gluon_normalization_id=resolved_contract.kk_gluon_normalization_id,
        mu_gs_semantics_id=resolved_contract.mu_gs_semantics_id,
        alpha_s_precision=resolved_contract.alpha_s_precision,
        rs_volume_policy_id=resolved_contract.rs_volume_policy_id,
        scale_policy_id=resolved_contract.scale_policy_id_for(resolved_scale),
        propagator_mass_rule_id=resolved_contract.propagator_mass_rule_id,
        m_KK_eff_GeV=resolved_scale.m_KK_eff_GeV,
    )


def build_coupling_contract_summary() -> dict[str, object]:
    """Return a deterministic summary of the frozen coupling contract."""
    contract = default_paper_0710_1869_coupling_contract()
    scale_point = default_paper_0710_1869_scales()
    shifted_match_scale = Paper07101869ScalePoint(
        label=f"{scale_point.label}_mu_match_shifted",
        Lambda_IR_GeV=scale_point.Lambda_IR_GeV,
        m_g1_GeV=scale_point.m_g1_GeV,
        xi_g=scale_point.xi_g,
        mu_match_GeV=scale_point.mu_match_GeV + 137.0,
        mu_gs_GeV=scale_point.mu_gs_GeV,
        m_KK_eff_GeV=scale_point.m_KK_eff_GeV,
    )
    normalization = evaluate_paper_0710_1869_gauge_coupling(
        scale_point,
        contract=contract,
    )
    shifted_match_normalization = evaluate_paper_0710_1869_gauge_coupling(
        shifted_match_scale,
        contract=contract,
    )
    return {
        "contract": contract.as_dict(),
        "default_scale_point": scale_point.as_dict(),
        "default_normalization": normalization.as_dict(),
        "checks": {
            "default_scale_point_exposes_mu_gs_and_mu_match": {
                "mu_gs_GeV",
                "mu_match_GeV",
            }.issubset(scale_point.as_dict()),
            "gs_is_independent_of_mu_match_variation": math.isclose(
                normalization.g_s_mu_gs,
                shifted_match_normalization.g_s_mu_gs,
                rel_tol=1.0e-12,
                abs_tol=0.0,
            )
            and math.isclose(
                normalization.alpha_s_mu_gs,
                shifted_match_normalization.alpha_s_mu_gs,
                rel_tol=1.0e-12,
                abs_tol=0.0,
            ),
            "propagator_mass_matches_scale_contract": (
                normalization.propagator_mass_GeV == scale_point.propagator_mass_GeV
            ),
            "single_vs_effective_policy_matches_scale_point": (
                normalization.scale_policy_id == contract.scale_policy_id_for(scale_point)
            ),
            "gs_is_derived_only_from_mu_gs": math.isclose(
                normalization.g_s_mu_gs,
                math.sqrt(4.0 * math.pi * normalization.alpha_s_mu_gs),
                rel_tol=1.0e-12,
                abs_tol=0.0,
            ),
            "rs_volume_sqrt_2L_matches_volume": math.isclose(
                normalization.rs_volume_sqrt_2L,
                math.sqrt(2.0 * normalization.rs_volume_L),
                rel_tol=1.0e-12,
                abs_tol=0.0,
            ),
            "gs_star_carries_rs_volume": math.isclose(
                normalization.g_s_star_mu_gs,
                normalization.g_s_mu_gs * normalization.rs_volume_sqrt_2L,
                rel_tol=1.0e-12,
                abs_tol=0.0,
            ),
        },
    }


def kk_gluon_normalization(
    scales: Paper07101869ScalePoint | None = None,
    *,
    contract: Paper07101869CouplingContract | None = None,
) -> dict[str, float | str | None]:
    """Return only the normalization data controlled by ``mu_gs_GeV``."""
    normalization = evaluate_paper_0710_1869_gauge_coupling(scales, contract=contract)
    return {
        "schema_id": normalization.schema_id,
        "contract_schema_id": normalization.contract_schema_id,
        "mode_id": normalization.mode_id,
        "paper_id": normalization.paper_id,
        "scale_label": normalization.scale_label,
        "Lambda_IR_GeV": normalization.Lambda_IR_GeV,
        "mu_gs_GeV": normalization.mu_gs_GeV,
        "m_g1_GeV": normalization.m_g1_GeV,
        "m_KK_eff_GeV": normalization.m_KK_eff_GeV,
        "propagator_mass_GeV": normalization.propagator_mass_GeV,
        "alpha_s_mu_gs": normalization.alpha_s_mu_gs,
        "g_s_mu_gs": normalization.g_s_mu_gs,
        "rs_volume_L": normalization.rs_volume_L,
        "rs_volume_sqrt_2L": normalization.rs_volume_sqrt_2L,
        "g_s_star_mu_gs": normalization.g_s_star_mu_gs,
        "kk_gluon_normalization_id": normalization.kk_gluon_normalization_id,
        "mu_gs_semantics_id": normalization.mu_gs_semantics_id,
        "alpha_s_precision": normalization.alpha_s_precision,
        "rs_volume_policy_id": normalization.rs_volume_policy_id,
        "scale_policy_id": normalization.scale_policy_id,
        "propagator_mass_rule_id": normalization.propagator_mass_rule_id,
    }


def kk_gluon_coupling_matrix(
    *,
    profiles: tuple[float, float, float] | list[float] | np.ndarray,
    u_left: list[list[float]] | list[list[complex]] | np.ndarray,
    u_right: list[list[float]] | list[list[complex]] | np.ndarray,
    scales: Paper07101869ScalePoint | None = None,
    contract: Paper07101869CouplingContract | None = None,
    subtract_universal: bool = False,
) -> list[list[float]]:
    """Return one real-valued compatibility matrix for contract tests."""
    normalization = evaluate_paper_0710_1869_gauge_coupling(scales, contract=contract)
    profile_vector = np.asarray(profiles, dtype=float)
    if profile_vector.shape != (3,):
        raise ValueError("profiles must have shape (3,)")
    if not np.all(np.isfinite(profile_vector)):
        raise ValueError("profiles must contain only finite values")
    if np.any(profile_vector <= 0.0):
        raise ValueError("profiles must contain only positive values")
    left = np.asarray(u_left, dtype=np.complex128)
    right = np.asarray(u_right, dtype=np.complex128)
    if left.shape != (3, 3) or right.shape != (3, 3):
        raise ValueError("u_left and u_right must both have shape (3, 3)")
    if not np.all(np.isfinite(left.real)) or not np.all(np.isfinite(left.imag)):
        raise ValueError("u_left must contain only finite values")
    if not np.all(np.isfinite(right.real)) or not np.all(np.isfinite(right.imag)):
        raise ValueError("u_right must contain only finite values")

    kernel = np.diag(profile_vector**2).astype(np.complex128)
    raw = 0.5 * (left.conjugate().T @ kernel @ right + right.conjugate().T @ kernel @ left)
    if subtract_universal:
        universal = float(np.trace(raw).real / 3.0)
        raw = raw - universal * np.eye(3, dtype=np.complex128)
        normalized = (
            normalization.g_s_mu_gs
            * normalization.rs_volume_sqrt_2L
            * 0.5
            * (raw + raw.conjugate().T)
        )
    else:
        hermitian_raw = 0.5 * (raw + raw.conjugate().T)
        normalized = normalization.g_s_mu_gs * (
            normalization.rs_volume_sqrt_2L * hermitian_raw
            - (np.eye(3, dtype=np.complex128) / normalization.rs_volume_sqrt_2L)
        )
    if not np.allclose(normalized.imag, 0.0, atol=1.0e-12):
        raise ValueError("compatibility matrix carries non-negligible imaginary parts")
    return [[float(value) for value in row] for row in normalized.real.tolist()]


def coupling_contract_summary() -> dict[str, object]:
    """Compatibility alias for acceptance tooling."""
    return build_coupling_contract_summary()


def paper_0710_1869_coupling_contract_summary() -> dict[str, object]:
    """Explicitly named compatibility alias for acceptance tooling."""
    return build_coupling_contract_summary()


__all__ = [
    "PAPER_0710_1869_COUPLINGS_SCHEMA_ID",
    "PAPER_0710_1869_DIMENSIONLESS_MATRIX_POLICY_ID",
    "PAPER_0710_1869_EFFECTIVE_SCALE_POLICY_ID",
    "PAPER_0710_1869_GAUGE_NORMALIZATION_SCHEMA_ID",
    "PAPER_0710_1869_GS_NORMALIZATION_ID",
    "PAPER_0710_1869_MU_GS_SEMANTICS_ID",
    "PAPER_0710_1869_PROPAGATOR_MASS_RULE_ID",
    "PAPER_0710_1869_RS_VOLUME_POLICY_ID",
    "PAPER_0710_1869_SINGLE_MODE_SCALE_POLICY_ID",
    "PAPER_0710_1869_UNIVERSAL_SUBTRACTION_POLICY_ID",
    "Paper07101869CouplingContract",
    "Paper07101869GaugeCouplingNormalization",
    "build_coupling_contract_summary",
    "coupling_contract_summary",
    "default_paper_0710_1869_coupling_contract",
    "evaluate_paper_0710_1869_gauge_coupling",
    "kk_gluon_coupling_matrix",
    "kk_gluon_normalization",
    "paper_0710_1869_coupling_contract_summary",
]
