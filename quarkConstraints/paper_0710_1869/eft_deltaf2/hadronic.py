"""Paper-owned hadronic inputs for the 0710.1869 ``Delta F = 2`` slice."""

from __future__ import annotations

import math
from dataclasses import dataclass

from qcd import alpha_s
from qcd.beta_function import beta_0

from ..conventions import PAPER_0710_1869_MODE_ID, PAPER_0710_1869_PAPER_ID
from ..validation import (
    require_known_schema_id,
    require_member,
    require_nonempty_identifier,
    require_positive_finite,
)
from .operators import (
    PAPER_0710_1869_DELTAF2_OPERATOR_BASIS_ID,
    PAPER_0710_1869_DELTAF2_OPERATOR_NORMALIZATION_ID,
    PAPER_0710_1869_DELTAF2_Q1_VLL,
    PAPER_0710_1869_DELTAF2_Q1_VRR,
    PAPER_0710_1869_DELTAF2_Q4_LR,
    PAPER_0710_1869_DELTAF2_Q5_LR,
)
from .rg_inputs import (
    PAPER_0710_1869_DELTAF2_DEFAULT_MU_HAD_GEV,
    PAPER_0710_1869_DELTAF2_RG_ALPHA_S_POLICY_ID,
    PAPER_0710_1869_DELTAF2_RG_SCHEME_ID,
    default_paper_0710_1869_rg_thresholds,
)

PAPER_0710_1869_DELTAF2_HADRONIC_SOURCE_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.eft_deltaf2.hadronic_source.v1"
)
PAPER_0710_1869_DELTAF2_KAON_HADRONIC_BUNDLE_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.eft_deltaf2.kaon_hadronic_bundle.v1"
)
PAPER_0710_1869_DELTAF2_KAON_HADRONIC_CONTRACT_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.eft_deltaf2.kaon_hadronic_contract.v1"
)
PAPER_0710_1869_DELTAF2_KAON_HADRONIC_SUMMARY_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.eft_deltaf2.kaon_hadronic_summary.v1"
)
PAPER_0710_1869_DELTAF2_BD_HADRONIC_BUNDLE_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.eft_deltaf2.bd_hadronic_bundle.v1"
)
PAPER_0710_1869_DELTAF2_BD_HADRONIC_CONTRACT_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.eft_deltaf2.bd_hadronic_contract.v1"
)
PAPER_0710_1869_DELTAF2_BD_HADRONIC_SUMMARY_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.eft_deltaf2.bd_hadronic_summary.v1"
)
PAPER_0710_1869_DELTAF2_BS_HADRONIC_BUNDLE_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.eft_deltaf2.bs_hadronic_bundle.v1"
)
PAPER_0710_1869_DELTAF2_BS_HADRONIC_CONTRACT_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.eft_deltaf2.bs_hadronic_contract.v1"
)
PAPER_0710_1869_DELTAF2_BS_HADRONIC_SUMMARY_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.eft_deltaf2.bs_hadronic_summary.v1"
)
PAPER_0710_1869_DELTAF2_KAON_LR_HADRONIC_INPUTS_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.eft_deltaf2.kaon_lr_hadronic_inputs.v1"
)
PAPER_0710_1869_DELTAF2_KAON_LR_HADRONIC_CONTRACT_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.eft_deltaf2.kaon_lr_hadronic_contract.v1"
)
PAPER_0710_1869_DELTAF2_KAON_LR_HADRONIC_SUMMARY_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.eft_deltaf2.kaon_lr_hadronic_summary.v1"
)
PAPER_0710_1869_DELTAF2_KAON_HADRONIC_SYSTEM_ID = "kaon"
PAPER_0710_1869_DELTAF2_BD_HADRONIC_SYSTEM_ID = "B_d"
PAPER_0710_1869_DELTAF2_BS_HADRONIC_SYSTEM_ID = "B_s"
PAPER_0710_1869_DELTAF2_HAMILTONIAN_CONVENTION_ID = (
    "heff.sum_ci_qi.no_hc_factor.v1"
)
PAPER_0710_1869_DELTAF2_KAON_MATRIX_ELEMENT_FORMULA_ID = (
    "kaon.q1_vll_vrr.8over3_fk2_mk2_bk_mu.v1"
)
PAPER_0710_1869_DELTAF2_BD_MATRIX_ELEMENT_FORMULA_ID = (
    "bd.q1_vll_vrr.8over3_fbd2_mbd2_bbd_mu.v1"
)
PAPER_0710_1869_DELTAF2_BS_MATRIX_ELEMENT_FORMULA_ID = (
    "bs.q1_vll_vrr.8over3_fbs2_mbs2_bbs_mu.v1"
)
PAPER_0710_1869_DELTAF2_KAON_LR_Q4_MATRIX_ELEMENT_FORMULA_ID = (
    "kaon.q4_lr.o4_scalar_lr.bv2004.eq5.matrix_element.mu_had.v1"
)
PAPER_0710_1869_DELTAF2_KAON_LR_Q5_MATRIX_ELEMENT_FORMULA_ID = (
    "kaon.q5_lr.o5_scalar_lr.bv2004.eq5.matrix_element.mu_had.v1"
)
PAPER_0710_1869_DELTAF2_KAON_LR_CHIRAL_RATIO_DEFINITION_ID = (
    "kaon.lr.chiral_ratio.r_chi.mu_had.custom_input.v1"
)
PAPER_0710_1869_DELTAF2_KAON_LR_FORMULA_SOURCE_ID = (
    "becirevic-villadoro.hep-lat-0408029.eq5.o4_o5.matrix_elements.v1"
)
PAPER_0710_1869_DELTAF2_KAON_LR_SCHEME_CONVERSION_POLICY_ID = (
    "none.caller_supplies_lr_inputs_in_declared_scheme_at_mu_had.v2"
)
PAPER_0710_1869_DELTAF2_KAON_LR_INPUT_POLICY_ID = (
    "custom_input_only.defaults_not_frozen.v1"
)
PAPER_0710_1869_DELTAF2_KAON_PARITY_RELATION_ID = (
    "kaon.q1_vll_equals_q1_vrr.by_parity.v1"
)
PAPER_0710_1869_DELTAF2_BD_PARITY_RELATION_ID = (
    "bd.q1_vll_equals_q1_vrr.by_parity.v1"
)
PAPER_0710_1869_DELTAF2_BS_PARITY_RELATION_ID = (
    "bs.q1_vll_equals_q1_vrr.by_parity.v1"
)
PAPER_0710_1869_DELTAF2_BAG_PARAMETER_SOURCE_SCHEME_ID = "rgi.flag21.average.v1"
PAPER_0710_1869_DELTAF2_BAG_PARAMETER_TRANSFORMATION_ID = (
    "hat_bk_to_bk_mu_had.lo_inverse_q1_running.v1"
)
PAPER_0710_1869_DELTAF2_HADRONIC_DEFAULT_PROVENANCE_MODE_ID = (
    "paper_sourced_defaults.v1"
)
PAPER_0710_1869_DELTAF2_HADRONIC_CUSTOM_PROVENANCE_MODE_ID = (
    "custom_inputs.validated_sources.v1"
)
PAPER_0710_1869_DELTAF2_LR_HADRONIC_CUSTOM_PROVENANCE_MODE_ID = (
    "custom_lr_inputs.validated_sources.v1"
)

_SUPPORTED_OPERATORS = (
    PAPER_0710_1869_DELTAF2_Q1_VLL,
    PAPER_0710_1869_DELTAF2_Q1_VRR,
)
_UNSUPPORTED_OPERATORS = (
    PAPER_0710_1869_DELTAF2_Q4_LR,
    PAPER_0710_1869_DELTAF2_Q5_LR,
)
_SUPPORTED_LR_OPERATORS = (
    PAPER_0710_1869_DELTAF2_Q4_LR,
    PAPER_0710_1869_DELTAF2_Q5_LR,
)
_UNSUPPORTED_LR_OPERATORS = (
    PAPER_0710_1869_DELTAF2_Q1_VLL,
    PAPER_0710_1869_DELTAF2_Q1_VRR,
)
_DEFAULT_M_K0_GEV = 0.497611
_DEFAULT_F_K_GEV = 0.1557
_DEFAULT_HAT_B_K_RGI = 0.7625
_DEFAULT_HADRONIC_BUNDLE_ID = "hadronic.kaon.pr5a.v1"
_DEFAULT_HADRONIC_SOURCE_ID = "hadronic.kaon.pdg2024_flag2024.v1"
_CUSTOM_HADRONIC_BUNDLE_ID = "hadronic.kaon.user_supplied.v1"
_CUSTOM_HADRONIC_SOURCE_ID = "hadronic.kaon.user_supplied.v1"
_CUSTOM_BD_HADRONIC_BUNDLE_ID = "hadronic.bd.user_supplied.v1"
_CUSTOM_BD_HADRONIC_SOURCE_ID = "hadronic.bd.user_supplied.v1"
_CUSTOM_BS_HADRONIC_BUNDLE_ID = "hadronic.bs.user_supplied.v1"
_CUSTOM_BS_HADRONIC_SOURCE_ID = "hadronic.bs.user_supplied.v1"
_CUSTOM_LR_HADRONIC_INPUTS_ID = "hadronic.kaon.lr.user_supplied.v1"
_CUSTOM_LR_HADRONIC_SOURCE_ID = "hadronic.kaon.lr.user_supplied.aggregate.v1"
_DEFAULT_MASS_SOURCE_ID = "pdg.2024.k0.mass.v1"
_DEFAULT_DECAY_CONSTANT_SOURCE_ID = "pdg.2024.fkplus.eq72.14.v1"
_DEFAULT_BAG_SOURCE_ID = "pdg.2024.flag.hat_bk.sec17.3.2.v1"
_CUSTOM_B_SYSTEMS = (
    PAPER_0710_1869_DELTAF2_BD_HADRONIC_SYSTEM_ID,
    PAPER_0710_1869_DELTAF2_BS_HADRONIC_SYSTEM_ID,
)


def _require_finite_float(name: str, value: float) -> float:
    numeric = float(value)
    if not math.isfinite(numeric):
        raise ValueError(f"{name} must be finite")
    return numeric


def _float_matches(lhs: float, rhs: float, *, atol: float = 1e-15) -> bool:
    return math.isclose(float(lhs), float(rhs), rel_tol=0.0, abs_tol=atol)


def _require_optional_scale(name: str, value: float | None) -> float | None:
    if value is None:
        return None
    return require_positive_finite(name, value)


def _n_f_at_scale(
    mu_GeV: float,
    thresholds: tuple[tuple[float, int, int], ...],
) -> int:
    if not thresholds:
        return 5
    n_f = thresholds[0][1]
    for mass, _, n_f_above in thresholds:
        if mu_GeV >= mass:
            n_f = n_f_above
        else:
            break
    return n_f


def _derive_b_k_mu_had(
    *,
    mu_had_GeV: float,
    hat_B_K_rgi: float,
) -> float:
    thresholds = default_paper_0710_1869_rg_thresholds()
    alpha_s_mu = float(
        alpha_s(
            mu_had_GeV,
            n_loops=1,
            matching_loops=0,
            thresholds=list(thresholds),
        )
    )
    n_f = _n_f_at_scale(mu_had_GeV, thresholds)
    exponent = 4.0 / (2.0 * float(beta_0(n_f)))
    return require_positive_finite(
        "B_K_mu_had",
        hat_B_K_rgi * (alpha_s_mu**exponent),
    )


@dataclass(frozen=True)
class Paper07101869HadronicSourceRef:
    """One externally sourced hadronic input reference."""

    source_id: str
    source_kind: str
    citation: str
    locator_label: str
    year: int
    renormalization_scheme_id: str | None = None
    scale_GeV: float | None = None
    transformation_id: str = "none"
    notes: str | None = None
    schema_id: str = PAPER_0710_1869_DELTAF2_HADRONIC_SOURCE_SCHEMA_ID

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "schema_id",
            require_known_schema_id(
                "schema_id",
                self.schema_id,
                expected=PAPER_0710_1869_DELTAF2_HADRONIC_SOURCE_SCHEMA_ID,
            ),
        )
        for field_name in (
            "source_id",
            "citation",
            "locator_label",
            "transformation_id",
        ):
            object.__setattr__(
                self,
                field_name,
                require_nonempty_identifier(field_name, getattr(self, field_name)),
            )
        object.__setattr__(
            self,
            "source_kind",
            require_member(
                "source_kind",
                self.source_kind,
                (
                    "particle-property",
                    "review-average",
                    "derived-from-rgi",
                    "lattice-result",
                    "derived-from-quark-masses",
                ),
            ),
        )
        if not isinstance(self.year, int) or isinstance(self.year, bool) or self.year < 1900:
            raise ValueError("year must be a positive integer year")
        if self.renormalization_scheme_id is not None:
            object.__setattr__(
                self,
                "renormalization_scheme_id",
                require_nonempty_identifier(
                    "renormalization_scheme_id", self.renormalization_scheme_id
                ),
            )
        object.__setattr__(
            self,
            "scale_GeV",
            _require_optional_scale("scale_GeV", self.scale_GeV),
        )
        if self.notes is not None:
            object.__setattr__(self, "notes", require_nonempty_identifier("notes", self.notes))

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_id": self.schema_id,
            "source_id": self.source_id,
            "source_kind": self.source_kind,
            "citation": self.citation,
            "locator_label": self.locator_label,
            "year": self.year,
            "renormalization_scheme_id": self.renormalization_scheme_id,
            "scale_GeV": self.scale_GeV,
            "transformation_id": self.transformation_id,
            "notes": self.notes,
        }


@dataclass(frozen=True)
class Paper07101869KaonHadronicContract:
    """Public hadronic contract view used by the kaon NP-only observable surface."""

    system_id: str = PAPER_0710_1869_DELTAF2_KAON_HADRONIC_SYSTEM_ID
    operator_basis_id: str = PAPER_0710_1869_DELTAF2_OPERATOR_BASIS_ID
    operator_normalization_id: str = PAPER_0710_1869_DELTAF2_OPERATOR_NORMALIZATION_ID
    renormalization_scheme_id: str = PAPER_0710_1869_DELTAF2_RG_SCHEME_ID
    mu_had_GeV: float = PAPER_0710_1869_DELTAF2_DEFAULT_MU_HAD_GEV
    evaluation_scale_GeV: float = PAPER_0710_1869_DELTAF2_DEFAULT_MU_HAD_GEV
    parity_relation_id: str = PAPER_0710_1869_DELTAF2_KAON_PARITY_RELATION_ID
    matrix_element_formula_id: str = PAPER_0710_1869_DELTAF2_KAON_MATRIX_ELEMENT_FORMULA_ID
    hamiltonian_convention_id: str = PAPER_0710_1869_DELTAF2_HAMILTONIAN_CONVENTION_ID
    schema_id: str = PAPER_0710_1869_DELTAF2_KAON_HADRONIC_CONTRACT_SCHEMA_ID

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "schema_id",
            require_known_schema_id(
                "schema_id",
                self.schema_id,
                expected=PAPER_0710_1869_DELTAF2_KAON_HADRONIC_CONTRACT_SCHEMA_ID,
            ),
        )
        for field_name in (
            "system_id",
            "operator_basis_id",
            "operator_normalization_id",
            "renormalization_scheme_id",
            "parity_relation_id",
            "matrix_element_formula_id",
            "hamiltonian_convention_id",
        ):
            object.__setattr__(
                self,
                field_name,
                require_nonempty_identifier(field_name, getattr(self, field_name)),
            )
        for field_name in ("mu_had_GeV", "evaluation_scale_GeV"):
            object.__setattr__(
                self,
                field_name,
                require_positive_finite(field_name, getattr(self, field_name)),
            )

    @property
    def scheme_id(self) -> str:
        return self.renormalization_scheme_id

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_id": self.schema_id,
            "system_id": self.system_id,
            "operator_basis_id": self.operator_basis_id,
            "operator_normalization_id": self.operator_normalization_id,
            "renormalization_scheme_id": self.renormalization_scheme_id,
            "scheme_id": self.scheme_id,
            "mu_had_GeV": self.mu_had_GeV,
            "evaluation_scale_GeV": self.evaluation_scale_GeV,
            "parity_relation_id": self.parity_relation_id,
            "matrix_element_formula_id": self.matrix_element_formula_id,
            "hamiltonian_convention_id": self.hamiltonian_convention_id,
        }


@dataclass(frozen=True)
class Paper07101869KaonHadronicBundle:
    """Frozen kaon hadronic inputs for the supported paper-mode observable subset."""

    m_K0_GeV: float
    f_K_GeV: float
    B_K_mu_had: float
    mass_source: Paper07101869HadronicSourceRef
    decay_constant_source: Paper07101869HadronicSourceRef
    bag_parameter_source: Paper07101869HadronicSourceRef
    contract: Paper07101869KaonHadronicContract | None = None
    system_id: str = PAPER_0710_1869_DELTAF2_KAON_HADRONIC_SYSTEM_ID
    operator_basis_id: str = PAPER_0710_1869_DELTAF2_OPERATOR_BASIS_ID
    operator_normalization_id: str = PAPER_0710_1869_DELTAF2_OPERATOR_NORMALIZATION_ID
    renormalization_scheme_id: str = PAPER_0710_1869_DELTAF2_RG_SCHEME_ID
    mu_had_GeV: float = PAPER_0710_1869_DELTAF2_DEFAULT_MU_HAD_GEV
    bundle_id: str = _DEFAULT_HADRONIC_BUNDLE_ID
    source_id: str = _DEFAULT_HADRONIC_SOURCE_ID
    provenance_ids: tuple[str, ...] = (_DEFAULT_HADRONIC_SOURCE_ID,)
    input_provenance_mode_id: str = PAPER_0710_1869_DELTAF2_HADRONIC_DEFAULT_PROVENANCE_MODE_ID
    hat_B_K_rgi_source_value: float = _DEFAULT_HAT_B_K_RGI
    bag_parameter_source_scheme_id: str = PAPER_0710_1869_DELTAF2_BAG_PARAMETER_SOURCE_SCHEME_ID
    bag_parameter_transformation_id: str = (
        PAPER_0710_1869_DELTAF2_BAG_PARAMETER_TRANSFORMATION_ID
    )
    alpha_s_policy_id: str = PAPER_0710_1869_DELTAF2_RG_ALPHA_S_POLICY_ID
    hamiltonian_convention_id: str = PAPER_0710_1869_DELTAF2_HAMILTONIAN_CONVENTION_ID
    matrix_element_formula_id: str = PAPER_0710_1869_DELTAF2_KAON_MATRIX_ELEMENT_FORMULA_ID
    parity_relation_id: str = PAPER_0710_1869_DELTAF2_KAON_PARITY_RELATION_ID
    supported_operator_names: tuple[str, ...] = _SUPPORTED_OPERATORS
    unsupported_operator_names: tuple[str, ...] = _UNSUPPORTED_OPERATORS
    schema_id: str = PAPER_0710_1869_DELTAF2_KAON_HADRONIC_BUNDLE_SCHEMA_ID
    mode_id: str = PAPER_0710_1869_MODE_ID
    paper_id: str = PAPER_0710_1869_PAPER_ID
    notes: str = (
        "PR5a kaon hadronic bundle: NP-only M12 support for Q1_VLL and Q1_VRR only. "
        "Non-zero LR operator contributions remain unsupported end-to-end."
    )

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "schema_id",
            require_known_schema_id(
                "schema_id",
                self.schema_id,
                expected=PAPER_0710_1869_DELTAF2_KAON_HADRONIC_BUNDLE_SCHEMA_ID,
            ),
        )
        require_member("mode_id", self.mode_id, (PAPER_0710_1869_MODE_ID,))
        require_member("paper_id", self.paper_id, (PAPER_0710_1869_PAPER_ID,))
        require_member(
            "system_id",
            self.system_id,
            (PAPER_0710_1869_DELTAF2_KAON_HADRONIC_SYSTEM_ID,),
        )
        object.__setattr__(
            self,
            "operator_basis_id",
            require_member(
                "operator_basis_id",
                self.operator_basis_id,
                (PAPER_0710_1869_DELTAF2_OPERATOR_BASIS_ID,),
            ),
        )
        object.__setattr__(
            self,
            "operator_normalization_id",
            require_member(
                "operator_normalization_id",
                self.operator_normalization_id,
                (PAPER_0710_1869_DELTAF2_OPERATOR_NORMALIZATION_ID,),
            ),
        )
        object.__setattr__(
            self,
            "renormalization_scheme_id",
            require_nonempty_identifier(
                "renormalization_scheme_id",
                self.renormalization_scheme_id,
            ),
        )
        if not isinstance(self.mass_source, Paper07101869HadronicSourceRef):
            raise ValueError("mass_source must be a Paper07101869HadronicSourceRef")
        if not isinstance(self.decay_constant_source, Paper07101869HadronicSourceRef):
            raise ValueError("decay_constant_source must be a Paper07101869HadronicSourceRef")
        if not isinstance(self.bag_parameter_source, Paper07101869HadronicSourceRef):
            raise ValueError("bag_parameter_source must be a Paper07101869HadronicSourceRef")
        resolved_contract = self.contract
        if resolved_contract is None:
            resolved_contract = Paper07101869KaonHadronicContract(
                system_id=self.system_id,
                operator_basis_id=self.operator_basis_id,
                operator_normalization_id=self.operator_normalization_id,
                renormalization_scheme_id=self.renormalization_scheme_id,
                mu_had_GeV=self.mu_had_GeV,
                evaluation_scale_GeV=self.mu_had_GeV,
                parity_relation_id=self.parity_relation_id,
                matrix_element_formula_id=self.matrix_element_formula_id,
                hamiltonian_convention_id=self.hamiltonian_convention_id,
            )
            object.__setattr__(self, "contract", resolved_contract)
        if not isinstance(resolved_contract, Paper07101869KaonHadronicContract):
            raise ValueError("contract must be a Paper07101869KaonHadronicContract")
        if resolved_contract.system_id != self.system_id:
            raise ValueError("contract.system_id must match the frozen bundle system_id")
        if resolved_contract.operator_basis_id != self.operator_basis_id:
            raise ValueError(
                "contract.operator_basis_id must match the frozen bundle operator_basis_id"
            )
        if (
            resolved_contract.operator_normalization_id
            != self.operator_normalization_id
        ):
            raise ValueError(
                "contract.operator_normalization_id must match the frozen bundle "
                "operator_normalization_id"
            )
        if (
            resolved_contract.renormalization_scheme_id
            != self.renormalization_scheme_id
        ):
            raise ValueError(
                "contract.renormalization_scheme_id must match the frozen bundle "
                "renormalization_scheme_id"
            )
        if resolved_contract.parity_relation_id != self.parity_relation_id:
            raise ValueError(
                "contract.parity_relation_id must match the bundle parity_relation_id"
            )
        if not _float_matches(resolved_contract.mu_had_GeV, self.mu_had_GeV, atol=1e-12):
            raise ValueError("contract.mu_had_GeV must match the frozen bundle mu_had_GeV")
        if not _float_matches(
            resolved_contract.evaluation_scale_GeV,
            self.mu_had_GeV,
            atol=1e-12,
        ):
            raise ValueError(
                "contract.evaluation_scale_GeV must match the frozen bundle mu_had_GeV"
            )
        for field_name in (
            "bundle_id",
            "source_id",
            "input_provenance_mode_id",
            "bag_parameter_source_scheme_id",
            "bag_parameter_transformation_id",
            "alpha_s_policy_id",
            "hamiltonian_convention_id",
            "matrix_element_formula_id",
            "parity_relation_id",
            "notes",
        ):
            object.__setattr__(
                self,
                field_name,
                require_nonempty_identifier(field_name, getattr(self, field_name)),
            )
        if not self.provenance_ids:
            raise ValueError("provenance_ids must contain at least one source reference")
        object.__setattr__(
            self,
            "provenance_ids",
            tuple(
                require_nonempty_identifier("provenance_id", provenance_id)
                for provenance_id in self.provenance_ids
            ),
        )
        if self.source_id not in self.provenance_ids:
            raise ValueError("source_id must be present in provenance_ids for every bundle mode")
        require_member(
            "input_provenance_mode_id",
            self.input_provenance_mode_id,
            (
                PAPER_0710_1869_DELTAF2_HADRONIC_DEFAULT_PROVENANCE_MODE_ID,
                PAPER_0710_1869_DELTAF2_HADRONIC_CUSTOM_PROVENANCE_MODE_ID,
            ),
        )
        for field_name in (
            "m_K0_GeV",
            "f_K_GeV",
            "B_K_mu_had",
            "mu_had_GeV",
            "hat_B_K_rgi_source_value",
        ):
            object.__setattr__(
                self,
                field_name,
                require_positive_finite(field_name, getattr(self, field_name)),
            )
        if self.bag_parameter_source.renormalization_scheme_id is None:
            raise ValueError(
                "bag_parameter_source must carry renormalization_scheme_id matching the "
                "bundle bag-parameter source scheme"
            )
        if (
            self.bag_parameter_source.renormalization_scheme_id
            != self.bag_parameter_source_scheme_id
        ):
            raise ValueError(
                "bag_parameter_source.renormalization_scheme_id must match "
                "bag_parameter_source_scheme_id"
            )
        if self.bag_parameter_source.transformation_id != self.bag_parameter_transformation_id:
            raise ValueError(
                "bag_parameter_source.transformation_id must match "
                "bag_parameter_transformation_id"
            )
        if (
            self.input_provenance_mode_id
            == PAPER_0710_1869_DELTAF2_HADRONIC_DEFAULT_PROVENANCE_MODE_ID
        ):
            require_member("bundle_id", self.bundle_id, (_DEFAULT_HADRONIC_BUNDLE_ID,))
            require_member("source_id", self.source_id, (_DEFAULT_HADRONIC_SOURCE_ID,))
            if self.provenance_ids != (_DEFAULT_HADRONIC_SOURCE_ID,):
                raise ValueError(
                    "default hadronic provenance mode must keep the frozen default provenance_ids"
                )
            if self.mass_source.source_id != _DEFAULT_MASS_SOURCE_ID:
                raise ValueError(
                    "default hadronic provenance mode must keep the frozen PDG mass source"
                )
            if self.decay_constant_source.source_id != _DEFAULT_DECAY_CONSTANT_SOURCE_ID:
                raise ValueError(
                    "default hadronic provenance mode must keep the frozen "
                    "PDG/FLAG decay-constant source"
                )
            if self.bag_parameter_source.source_id != _DEFAULT_BAG_SOURCE_ID:
                raise ValueError(
                    "default hadronic provenance mode must keep the frozen FLAG "
                    "bag-parameter source"
                )
            if not _float_matches(self.m_K0_GeV, _DEFAULT_M_K0_GEV):
                raise ValueError(
                    "default hadronic provenance mode requires the frozen PDG m_K0 value"
                )
            if not _float_matches(self.f_K_GeV, _DEFAULT_F_K_GEV):
                raise ValueError(
                    "default hadronic provenance mode requires the frozen PDG/FLAG f_K value"
                )
            if not _float_matches(self.hat_B_K_rgi_source_value, _DEFAULT_HAT_B_K_RGI):
                raise ValueError(
                    "default hadronic provenance mode requires the frozen FLAG hat_B_K source value"
                )
            expected_b_k = _derive_b_k_mu_had(
                mu_had_GeV=self.mu_had_GeV,
                hat_B_K_rgi=_DEFAULT_HAT_B_K_RGI,
            )
            if not _float_matches(self.B_K_mu_had, expected_b_k):
                raise ValueError(
                    "default hadronic provenance mode requires B_K_mu_had derived from the "
                    "frozen FLAG hat_B_K source value"
                )
        if (
            self.input_provenance_mode_id
            == PAPER_0710_1869_DELTAF2_HADRONIC_CUSTOM_PROVENANCE_MODE_ID
        ):
            if self.bundle_id == _DEFAULT_HADRONIC_BUNDLE_ID:
                raise ValueError(
                    "custom hadronic provenance mode must not reuse the default bundle_id"
                )
            if self.source_id == _DEFAULT_HADRONIC_SOURCE_ID:
                raise ValueError(
                    "custom hadronic provenance mode must not reuse the default source_id"
                )
            if self.provenance_ids == (_DEFAULT_HADRONIC_SOURCE_ID,):
                raise ValueError(
                    "custom hadronic provenance mode must not reuse the default provenance_ids"
                )
        if tuple(self.supported_operator_names) != _SUPPORTED_OPERATORS:
            raise ValueError("supported_operator_names must match the supported kaon subset")
        if tuple(self.unsupported_operator_names) != _UNSUPPORTED_OPERATORS:
            raise ValueError("unsupported_operator_names must match the guarded LR subset")

    @property
    def q1_matrix_element_GeV4(self) -> float:
        return _require_finite_float(
            "q1_matrix_element_GeV4",
            (8.0 / 3.0)
            * (self.f_K_GeV**2)
            * (self.m_K0_GeV**2)
            * self.B_K_mu_had,
        )

    @property
    def tags(self) -> dict[str, object]:
        return {
            "mode_id": self.mode_id,
            "paper_id": self.paper_id,
            "system_id": self.system_id,
            "operator_basis_id": self.operator_basis_id,
            "operator_normalization_id": self.operator_normalization_id,
            "renormalization_scheme_id": self.renormalization_scheme_id,
            "contract": self.contract.as_dict(),
            "bundle_id": self.bundle_id,
            "source_id": self.source_id,
            "provenance_ids": list(self.provenance_ids),
            "input_provenance_mode_id": self.input_provenance_mode_id,
            "mu_had_GeV": self.mu_had_GeV,
            "evaluation_scale_GeV": self.contract.evaluation_scale_GeV,
            "supported_operator_names": list(self.supported_operator_names),
            "unsupported_operator_names": list(self.unsupported_operator_names),
            "matrix_element_formula_id": self.matrix_element_formula_id,
            "hamiltonian_convention_id": self.hamiltonian_convention_id,
            "parity_relation_id": self.parity_relation_id,
        }

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_id": self.schema_id,
            "mode_id": self.mode_id,
            "paper_id": self.paper_id,
            "system_id": self.system_id,
            "operator_basis_id": self.operator_basis_id,
            "operator_normalization_id": self.operator_normalization_id,
            "renormalization_scheme_id": self.renormalization_scheme_id,
            "scheme_id": self.contract.scheme_id,
            "contract": self.contract.as_dict(),
            "bundle_id": self.bundle_id,
            "source_id": self.source_id,
            "provenance_ids": list(self.provenance_ids),
            "input_provenance_mode_id": self.input_provenance_mode_id,
            "mu_had_GeV": self.mu_had_GeV,
            "evaluation_scale_GeV": self.contract.evaluation_scale_GeV,
            "m_K0_GeV": self.m_K0_GeV,
            "f_K_GeV": self.f_K_GeV,
            "B_K_mu_had": self.B_K_mu_had,
            "hat_B_K_rgi_source_value": self.hat_B_K_rgi_source_value,
            "bag_parameter_source_scheme_id": self.bag_parameter_source_scheme_id,
            "bag_parameter_transformation_id": self.bag_parameter_transformation_id,
            "alpha_s_policy_id": self.alpha_s_policy_id,
            "matrix_element_formula_id": self.matrix_element_formula_id,
            "hamiltonian_convention_id": self.hamiltonian_convention_id,
            "parity_relation_id": self.parity_relation_id,
            "supported_operator_names": list(self.supported_operator_names),
            "unsupported_operator_names": list(self.unsupported_operator_names),
            "mass_source": self.mass_source.as_dict(),
            "decay_constant_source": self.decay_constant_source.as_dict(),
            "bag_parameter_source": self.bag_parameter_source.as_dict(),
            "q1_matrix_element_GeV4": self.q1_matrix_element_GeV4,
            "tags": self.tags,
            "notes": self.notes,
        }

    def summary(self) -> dict[str, object]:
        return {
            **self.as_dict(),
            "schema_id": PAPER_0710_1869_DELTAF2_KAON_HADRONIC_SUMMARY_SCHEMA_ID,
        }


_B_MESON_CONTRACT_SCHEMA_BY_SYSTEM = {
    PAPER_0710_1869_DELTAF2_BD_HADRONIC_SYSTEM_ID: (
        PAPER_0710_1869_DELTAF2_BD_HADRONIC_CONTRACT_SCHEMA_ID
    ),
    PAPER_0710_1869_DELTAF2_BS_HADRONIC_SYSTEM_ID: (
        PAPER_0710_1869_DELTAF2_BS_HADRONIC_CONTRACT_SCHEMA_ID
    ),
}
_B_MESON_BUNDLE_SCHEMA_BY_SYSTEM = {
    PAPER_0710_1869_DELTAF2_BD_HADRONIC_SYSTEM_ID: (
        PAPER_0710_1869_DELTAF2_BD_HADRONIC_BUNDLE_SCHEMA_ID
    ),
    PAPER_0710_1869_DELTAF2_BS_HADRONIC_SYSTEM_ID: (
        PAPER_0710_1869_DELTAF2_BS_HADRONIC_BUNDLE_SCHEMA_ID
    ),
}
_B_MESON_SUMMARY_SCHEMA_BY_SYSTEM = {
    PAPER_0710_1869_DELTAF2_BD_HADRONIC_SYSTEM_ID: (
        PAPER_0710_1869_DELTAF2_BD_HADRONIC_SUMMARY_SCHEMA_ID
    ),
    PAPER_0710_1869_DELTAF2_BS_HADRONIC_SYSTEM_ID: (
        PAPER_0710_1869_DELTAF2_BS_HADRONIC_SUMMARY_SCHEMA_ID
    ),
}
_B_MESON_MATRIX_ELEMENT_FORMULA_BY_SYSTEM = {
    PAPER_0710_1869_DELTAF2_BD_HADRONIC_SYSTEM_ID: (
        PAPER_0710_1869_DELTAF2_BD_MATRIX_ELEMENT_FORMULA_ID
    ),
    PAPER_0710_1869_DELTAF2_BS_HADRONIC_SYSTEM_ID: (
        PAPER_0710_1869_DELTAF2_BS_MATRIX_ELEMENT_FORMULA_ID
    ),
}
_B_MESON_PARITY_RELATION_BY_SYSTEM = {
    PAPER_0710_1869_DELTAF2_BD_HADRONIC_SYSTEM_ID: (
        PAPER_0710_1869_DELTAF2_BD_PARITY_RELATION_ID
    ),
    PAPER_0710_1869_DELTAF2_BS_HADRONIC_SYSTEM_ID: (
        PAPER_0710_1869_DELTAF2_BS_PARITY_RELATION_ID
    ),
}
_B_MESON_BUNDLE_ID_BY_SYSTEM = {
    PAPER_0710_1869_DELTAF2_BD_HADRONIC_SYSTEM_ID: _CUSTOM_BD_HADRONIC_BUNDLE_ID,
    PAPER_0710_1869_DELTAF2_BS_HADRONIC_SYSTEM_ID: _CUSTOM_BS_HADRONIC_BUNDLE_ID,
}
_B_MESON_SOURCE_ID_BY_SYSTEM = {
    PAPER_0710_1869_DELTAF2_BD_HADRONIC_SYSTEM_ID: _CUSTOM_BD_HADRONIC_SOURCE_ID,
    PAPER_0710_1869_DELTAF2_BS_HADRONIC_SYSTEM_ID: _CUSTOM_BS_HADRONIC_SOURCE_ID,
}


def _require_custom_b_system_id(name: str, value: str) -> str:
    return require_member(name, value, _CUSTOM_B_SYSTEMS)


@dataclass(frozen=True)
class Paper07101869BMesonHadronicContract:
    """Custom-input-only Q1 hadronic contract for B_d/B_s observables."""

    system_id: str
    renormalization_scheme_id: str
    mu_had_GeV: float
    evaluation_scale_GeV: float
    operator_basis_id: str = PAPER_0710_1869_DELTAF2_OPERATOR_BASIS_ID
    operator_normalization_id: str = PAPER_0710_1869_DELTAF2_OPERATOR_NORMALIZATION_ID
    parity_relation_id: str | None = None
    matrix_element_formula_id: str | None = None
    hamiltonian_convention_id: str = PAPER_0710_1869_DELTAF2_HAMILTONIAN_CONVENTION_ID
    schema_id: str | None = None

    def __post_init__(self) -> None:
        resolved_system_id = _require_custom_b_system_id("system_id", self.system_id)
        object.__setattr__(self, "system_id", resolved_system_id)
        expected_schema_id = _B_MESON_CONTRACT_SCHEMA_BY_SYSTEM[resolved_system_id]
        object.__setattr__(
            self,
            "schema_id",
            require_known_schema_id(
                "schema_id",
                expected_schema_id if self.schema_id is None else self.schema_id,
                expected=expected_schema_id,
            ),
        )
        object.__setattr__(
            self,
            "operator_basis_id",
            require_member(
                "operator_basis_id",
                self.operator_basis_id,
                (PAPER_0710_1869_DELTAF2_OPERATOR_BASIS_ID,),
            ),
        )
        object.__setattr__(
            self,
            "operator_normalization_id",
            require_member(
                "operator_normalization_id",
                self.operator_normalization_id,
                (PAPER_0710_1869_DELTAF2_OPERATOR_NORMALIZATION_ID,),
            ),
        )
        object.__setattr__(
            self,
            "renormalization_scheme_id",
            require_nonempty_identifier(
                "renormalization_scheme_id", self.renormalization_scheme_id
            ),
        )
        for field_name in ("mu_had_GeV", "evaluation_scale_GeV"):
            object.__setattr__(
                self,
                field_name,
                require_positive_finite(field_name, getattr(self, field_name)),
            )
        resolved_parity_relation_id = (
            _B_MESON_PARITY_RELATION_BY_SYSTEM[resolved_system_id]
            if self.parity_relation_id is None
            else self.parity_relation_id
        )
        object.__setattr__(
            self,
            "parity_relation_id",
            require_member(
                "parity_relation_id",
                resolved_parity_relation_id,
                (_B_MESON_PARITY_RELATION_BY_SYSTEM[resolved_system_id],),
            ),
        )
        resolved_formula_id = (
            _B_MESON_MATRIX_ELEMENT_FORMULA_BY_SYSTEM[resolved_system_id]
            if self.matrix_element_formula_id is None
            else self.matrix_element_formula_id
        )
        object.__setattr__(
            self,
            "matrix_element_formula_id",
            require_member(
                "matrix_element_formula_id",
                resolved_formula_id,
                (_B_MESON_MATRIX_ELEMENT_FORMULA_BY_SYSTEM[resolved_system_id],),
            ),
        )
        object.__setattr__(
            self,
            "hamiltonian_convention_id",
            require_member(
                "hamiltonian_convention_id",
                self.hamiltonian_convention_id,
                (PAPER_0710_1869_DELTAF2_HAMILTONIAN_CONVENTION_ID,),
            ),
        )

    @property
    def scheme_id(self) -> str:
        return self.renormalization_scheme_id

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_id": self.schema_id,
            "system_id": self.system_id,
            "operator_basis_id": self.operator_basis_id,
            "operator_normalization_id": self.operator_normalization_id,
            "renormalization_scheme_id": self.renormalization_scheme_id,
            "scheme_id": self.scheme_id,
            "mu_had_GeV": self.mu_had_GeV,
            "evaluation_scale_GeV": self.evaluation_scale_GeV,
            "parity_relation_id": self.parity_relation_id,
            "matrix_element_formula_id": self.matrix_element_formula_id,
            "hamiltonian_convention_id": self.hamiltonian_convention_id,
        }


@dataclass(frozen=True)
class Paper07101869BMesonHadronicBundle:
    """Custom-input-only Q1 hadronic bundle for B_d/B_s observables."""

    system_id: str
    meson_mass_GeV: float
    meson_decay_constant_GeV: float
    bag_parameter_mu_had: float
    mass_source: Paper07101869HadronicSourceRef
    decay_constant_source: Paper07101869HadronicSourceRef
    bag_parameter_source: Paper07101869HadronicSourceRef
    renormalization_scheme_id: str
    mu_had_GeV: float
    contract: Paper07101869BMesonHadronicContract | None = None
    operator_basis_id: str = PAPER_0710_1869_DELTAF2_OPERATOR_BASIS_ID
    operator_normalization_id: str = PAPER_0710_1869_DELTAF2_OPERATOR_NORMALIZATION_ID
    bundle_id: str | None = None
    source_id: str | None = None
    provenance_ids: tuple[str, ...] | None = None
    input_provenance_mode_id: str = PAPER_0710_1869_DELTAF2_HADRONIC_CUSTOM_PROVENANCE_MODE_ID
    hamiltonian_convention_id: str = PAPER_0710_1869_DELTAF2_HAMILTONIAN_CONVENTION_ID
    matrix_element_formula_id: str | None = None
    parity_relation_id: str | None = None
    supported_operator_names: tuple[str, ...] = _SUPPORTED_OPERATORS
    unsupported_operator_names: tuple[str, ...] = _UNSUPPORTED_OPERATORS
    schema_id: str | None = None
    mode_id: str = PAPER_0710_1869_MODE_ID
    paper_id: str = PAPER_0710_1869_PAPER_ID
    notes: str | None = None

    def __post_init__(self) -> None:
        resolved_system_id = _require_custom_b_system_id("system_id", self.system_id)
        object.__setattr__(self, "system_id", resolved_system_id)
        expected_schema_id = _B_MESON_BUNDLE_SCHEMA_BY_SYSTEM[resolved_system_id]
        object.__setattr__(
            self,
            "schema_id",
            require_known_schema_id(
                "schema_id",
                expected_schema_id if self.schema_id is None else self.schema_id,
                expected=expected_schema_id,
            ),
        )
        require_member("mode_id", self.mode_id, (PAPER_0710_1869_MODE_ID,))
        require_member("paper_id", self.paper_id, (PAPER_0710_1869_PAPER_ID,))
        object.__setattr__(
            self,
            "operator_basis_id",
            require_member(
                "operator_basis_id",
                self.operator_basis_id,
                (PAPER_0710_1869_DELTAF2_OPERATOR_BASIS_ID,),
            ),
        )
        object.__setattr__(
            self,
            "operator_normalization_id",
            require_member(
                "operator_normalization_id",
                self.operator_normalization_id,
                (PAPER_0710_1869_DELTAF2_OPERATOR_NORMALIZATION_ID,),
            ),
        )
        object.__setattr__(
            self,
            "renormalization_scheme_id",
            require_nonempty_identifier(
                "renormalization_scheme_id", self.renormalization_scheme_id
            ),
        )
        for field_name in (
            "meson_mass_GeV",
            "meson_decay_constant_GeV",
            "bag_parameter_mu_had",
            "mu_had_GeV",
        ):
            object.__setattr__(
                self,
                field_name,
                require_positive_finite(field_name, getattr(self, field_name)),
            )
        for field_name in ("mass_source", "decay_constant_source", "bag_parameter_source"):
            if not isinstance(getattr(self, field_name), Paper07101869HadronicSourceRef):
                raise ValueError(f"{field_name} must be a Paper07101869HadronicSourceRef")
        resolved_formula_id = (
            _B_MESON_MATRIX_ELEMENT_FORMULA_BY_SYSTEM[resolved_system_id]
            if self.matrix_element_formula_id is None
            else self.matrix_element_formula_id
        )
        object.__setattr__(
            self,
            "matrix_element_formula_id",
            require_member(
                "matrix_element_formula_id",
                resolved_formula_id,
                (_B_MESON_MATRIX_ELEMENT_FORMULA_BY_SYSTEM[resolved_system_id],),
            ),
        )
        resolved_parity_relation_id = (
            _B_MESON_PARITY_RELATION_BY_SYSTEM[resolved_system_id]
            if self.parity_relation_id is None
            else self.parity_relation_id
        )
        object.__setattr__(
            self,
            "parity_relation_id",
            require_member(
                "parity_relation_id",
                resolved_parity_relation_id,
                (_B_MESON_PARITY_RELATION_BY_SYSTEM[resolved_system_id],),
            ),
        )
        object.__setattr__(
            self,
            "hamiltonian_convention_id",
            require_member(
                "hamiltonian_convention_id",
                self.hamiltonian_convention_id,
                (PAPER_0710_1869_DELTAF2_HAMILTONIAN_CONVENTION_ID,),
            ),
        )
        object.__setattr__(
            self,
            "input_provenance_mode_id",
            require_member(
                "input_provenance_mode_id",
                self.input_provenance_mode_id,
                (PAPER_0710_1869_DELTAF2_HADRONIC_CUSTOM_PROVENANCE_MODE_ID,),
            ),
        )
        resolved_contract = self.contract
        if resolved_contract is None:
            resolved_contract = Paper07101869BMesonHadronicContract(
                system_id=resolved_system_id,
                operator_basis_id=self.operator_basis_id,
                operator_normalization_id=self.operator_normalization_id,
                renormalization_scheme_id=self.renormalization_scheme_id,
                mu_had_GeV=self.mu_had_GeV,
                evaluation_scale_GeV=self.mu_had_GeV,
                parity_relation_id=self.parity_relation_id,
                matrix_element_formula_id=self.matrix_element_formula_id,
                hamiltonian_convention_id=self.hamiltonian_convention_id,
            )
            object.__setattr__(self, "contract", resolved_contract)
        if not isinstance(resolved_contract, Paper07101869BMesonHadronicContract):
            raise ValueError("contract must be a Paper07101869BMesonHadronicContract")
        if resolved_contract.system_id != resolved_system_id:
            raise ValueError("contract.system_id must match the B-system bundle system_id")
        if resolved_contract.operator_basis_id != self.operator_basis_id:
            raise ValueError(
                "contract.operator_basis_id must match the B-system bundle operator_basis_id"
            )
        if resolved_contract.operator_normalization_id != self.operator_normalization_id:
            raise ValueError(
                "contract.operator_normalization_id must match the B-system bundle "
                "operator_normalization_id"
            )
        if resolved_contract.renormalization_scheme_id != self.renormalization_scheme_id:
            raise ValueError(
                "contract.renormalization_scheme_id must match the B-system bundle "
                "renormalization_scheme_id"
            )
        if resolved_contract.parity_relation_id != self.parity_relation_id:
            raise ValueError("contract.parity_relation_id must match the bundle parity_relation_id")
        if resolved_contract.matrix_element_formula_id != self.matrix_element_formula_id:
            raise ValueError(
                "contract.matrix_element_formula_id must match the bundle matrix_element_formula_id"
            )
        if resolved_contract.hamiltonian_convention_id != self.hamiltonian_convention_id:
            raise ValueError(
                "contract.hamiltonian_convention_id must match the bundle "
                "hamiltonian_convention_id"
            )
        if not _float_matches(resolved_contract.mu_had_GeV, self.mu_had_GeV, atol=1e-12):
            raise ValueError("contract.mu_had_GeV must match the bundle mu_had_GeV")
        if not _float_matches(
            resolved_contract.evaluation_scale_GeV,
            self.mu_had_GeV,
            atol=1e-12,
        ):
            raise ValueError("contract.evaluation_scale_GeV must match the bundle mu_had_GeV")
        resolved_bundle_id = (
            _B_MESON_BUNDLE_ID_BY_SYSTEM[resolved_system_id]
            if self.bundle_id is None
            else require_nonempty_identifier("bundle_id", self.bundle_id)
        )
        resolved_source_id = (
            _B_MESON_SOURCE_ID_BY_SYSTEM[resolved_system_id]
            if self.source_id is None
            else require_nonempty_identifier("source_id", self.source_id)
        )
        resolved_provenance_ids = (
            (resolved_source_id,)
            if self.provenance_ids is None
            else tuple(
                require_nonempty_identifier("provenance_id", provenance_id)
                for provenance_id in self.provenance_ids
            )
        )
        if resolved_bundle_id == _DEFAULT_HADRONIC_BUNDLE_ID:
            raise ValueError("custom B-system hadronic inputs must not reuse the default bundle_id")
        if resolved_source_id == _DEFAULT_HADRONIC_SOURCE_ID:
            raise ValueError("custom B-system hadronic inputs must not reuse the default source_id")
        if resolved_provenance_ids == (_DEFAULT_HADRONIC_SOURCE_ID,):
            raise ValueError(
                "custom B-system hadronic inputs must not reuse the default provenance_ids"
            )
        if resolved_source_id not in resolved_provenance_ids:
            raise ValueError("source_id must be present in provenance_ids for custom B bundles")
        object.__setattr__(self, "bundle_id", resolved_bundle_id)
        object.__setattr__(self, "source_id", resolved_source_id)
        object.__setattr__(self, "provenance_ids", resolved_provenance_ids)
        if self.bag_parameter_source.renormalization_scheme_id is None:
            raise ValueError(
                "bag_parameter_source must carry renormalization_scheme_id matching the "
                "B-system bundle renormalization scheme"
            )
        if self.bag_parameter_source.renormalization_scheme_id != self.renormalization_scheme_id:
            raise ValueError(
                "bag_parameter_source.renormalization_scheme_id must match "
                "renormalization_scheme_id"
            )
        if self.bag_parameter_source.scale_GeV is None:
            raise ValueError(
                "bag_parameter_source must carry scale_GeV matching the bundle mu_had_GeV"
            )
        if not _float_matches(self.bag_parameter_source.scale_GeV, self.mu_had_GeV, atol=1e-12):
            raise ValueError("bag_parameter_source.scale_GeV must match mu_had_GeV")
        if tuple(self.supported_operator_names) != _SUPPORTED_OPERATORS:
            raise ValueError("supported_operator_names must match the Q1-only supported subset")
        if tuple(self.unsupported_operator_names) != _UNSUPPORTED_OPERATORS:
            raise ValueError("unsupported_operator_names must match the guarded LR subset")
        resolved_notes = (
            f"Custom {resolved_system_id} hadronic bundle: NP-only Q1_VLL/Q1_VRR support at "
            "declared mu_had and renormalization scheme. No frozen default lattice dataset "
            "is assumed in this slice."
            if self.notes is None
            else self.notes
        )
        object.__setattr__(self, "notes", require_nonempty_identifier("notes", resolved_notes))

    @property
    def q1_matrix_element_GeV4(self) -> float:
        return _require_finite_float(
            "q1_matrix_element_GeV4",
            (8.0 / 3.0)
            * (self.meson_decay_constant_GeV**2)
            * (self.meson_mass_GeV**2)
            * self.bag_parameter_mu_had,
        )

    @property
    def tags(self) -> dict[str, object]:
        return {
            "mode_id": self.mode_id,
            "paper_id": self.paper_id,
            "system_id": self.system_id,
            "operator_basis_id": self.operator_basis_id,
            "operator_normalization_id": self.operator_normalization_id,
            "renormalization_scheme_id": self.renormalization_scheme_id,
            "contract": self.contract.as_dict(),
            "bundle_id": self.bundle_id,
            "source_id": self.source_id,
            "provenance_ids": list(self.provenance_ids),
            "input_provenance_mode_id": self.input_provenance_mode_id,
            "mu_had_GeV": self.mu_had_GeV,
            "evaluation_scale_GeV": self.contract.evaluation_scale_GeV,
            "supported_operator_names": list(self.supported_operator_names),
            "unsupported_operator_names": list(self.unsupported_operator_names),
            "matrix_element_formula_id": self.matrix_element_formula_id,
            "hamiltonian_convention_id": self.hamiltonian_convention_id,
            "parity_relation_id": self.parity_relation_id,
        }

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_id": self.schema_id,
            "mode_id": self.mode_id,
            "paper_id": self.paper_id,
            "system_id": self.system_id,
            "operator_basis_id": self.operator_basis_id,
            "operator_normalization_id": self.operator_normalization_id,
            "renormalization_scheme_id": self.renormalization_scheme_id,
            "scheme_id": self.contract.scheme_id,
            "contract": self.contract.as_dict(),
            "bundle_id": self.bundle_id,
            "source_id": self.source_id,
            "provenance_ids": list(self.provenance_ids),
            "input_provenance_mode_id": self.input_provenance_mode_id,
            "mu_had_GeV": self.mu_had_GeV,
            "evaluation_scale_GeV": self.contract.evaluation_scale_GeV,
            "meson_mass_GeV": self.meson_mass_GeV,
            "meson_decay_constant_GeV": self.meson_decay_constant_GeV,
            "bag_parameter_mu_had": self.bag_parameter_mu_had,
            "matrix_element_formula_id": self.matrix_element_formula_id,
            "hamiltonian_convention_id": self.hamiltonian_convention_id,
            "parity_relation_id": self.parity_relation_id,
            "supported_operator_names": list(self.supported_operator_names),
            "unsupported_operator_names": list(self.unsupported_operator_names),
            "mass_source": self.mass_source.as_dict(),
            "decay_constant_source": self.decay_constant_source.as_dict(),
            "bag_parameter_source": self.bag_parameter_source.as_dict(),
            "q1_matrix_element_GeV4": self.q1_matrix_element_GeV4,
            "tags": self.tags,
            "notes": self.notes,
        }

    def summary(self) -> dict[str, object]:
        return {
            **self.as_dict(),
            "schema_id": _B_MESON_SUMMARY_SCHEMA_BY_SYSTEM[self.system_id],
        }


@dataclass(frozen=True)
class Paper07101869KaonLRHadronicContract:
    """Frozen contract for custom LR kaon matrix elements at ``mu_had``."""

    system_id: str = PAPER_0710_1869_DELTAF2_KAON_HADRONIC_SYSTEM_ID
    operator_basis_id: str = PAPER_0710_1869_DELTAF2_OPERATOR_BASIS_ID
    operator_normalization_id: str = PAPER_0710_1869_DELTAF2_OPERATOR_NORMALIZATION_ID
    renormalization_scheme_id: str = PAPER_0710_1869_DELTAF2_RG_SCHEME_ID
    mu_had_GeV: float = PAPER_0710_1869_DELTAF2_DEFAULT_MU_HAD_GEV
    evaluation_scale_GeV: float = PAPER_0710_1869_DELTAF2_DEFAULT_MU_HAD_GEV
    hamiltonian_convention_id: str = PAPER_0710_1869_DELTAF2_HAMILTONIAN_CONVENTION_ID
    q4_matrix_element_formula_id: str = (
        PAPER_0710_1869_DELTAF2_KAON_LR_Q4_MATRIX_ELEMENT_FORMULA_ID
    )
    q5_matrix_element_formula_id: str = (
        PAPER_0710_1869_DELTAF2_KAON_LR_Q5_MATRIX_ELEMENT_FORMULA_ID
    )
    chiral_ratio_definition_id: str = (
        PAPER_0710_1869_DELTAF2_KAON_LR_CHIRAL_RATIO_DEFINITION_ID
    )
    scheme_conversion_policy_id: str = (
        PAPER_0710_1869_DELTAF2_KAON_LR_SCHEME_CONVERSION_POLICY_ID
    )
    input_policy_id: str = PAPER_0710_1869_DELTAF2_KAON_LR_INPUT_POLICY_ID
    formula_source_id: str = PAPER_0710_1869_DELTAF2_KAON_LR_FORMULA_SOURCE_ID
    formula_citation: str = "Becirevic and Villadoro, hep-lat/0408029, Eq. (5)"
    formula_locator_label: str = "Eq. (5)"
    schema_id: str = PAPER_0710_1869_DELTAF2_KAON_LR_HADRONIC_CONTRACT_SCHEMA_ID
    notes: str = (
        "LR-HAD-1 freezes only the O4/O5 scalar LR matrix-element formulas at mu_had. "
        "B4, B5, and R_chi must be supplied explicitly in the bundle's declared "
        "renormalization scheme and at mu_had; LR-HAD-1 performs no hidden scheme "
        "conversion. Defaults for LR hadronic inputs are not frozen yet."
    )

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "schema_id",
            require_known_schema_id(
                "schema_id",
                self.schema_id,
                expected=PAPER_0710_1869_DELTAF2_KAON_LR_HADRONIC_CONTRACT_SCHEMA_ID,
            ),
        )
        for field_name in (
            "system_id",
            "operator_basis_id",
            "operator_normalization_id",
            "renormalization_scheme_id",
            "hamiltonian_convention_id",
            "q4_matrix_element_formula_id",
            "q5_matrix_element_formula_id",
            "chiral_ratio_definition_id",
            "scheme_conversion_policy_id",
            "input_policy_id",
            "formula_source_id",
            "formula_citation",
            "formula_locator_label",
            "notes",
        ):
            object.__setattr__(
                self,
                field_name,
                require_nonempty_identifier(field_name, getattr(self, field_name)),
            )
        for field_name in ("mu_had_GeV", "evaluation_scale_GeV"):
            object.__setattr__(
                self,
                field_name,
                require_positive_finite(field_name, getattr(self, field_name)),
            )

    @property
    def scheme_id(self) -> str:
        return self.renormalization_scheme_id

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_id": self.schema_id,
            "system_id": self.system_id,
            "operator_basis_id": self.operator_basis_id,
            "operator_normalization_id": self.operator_normalization_id,
            "renormalization_scheme_id": self.renormalization_scheme_id,
            "scheme_id": self.scheme_id,
            "mu_had_GeV": self.mu_had_GeV,
            "evaluation_scale_GeV": self.evaluation_scale_GeV,
            "hamiltonian_convention_id": self.hamiltonian_convention_id,
            "q4_matrix_element_formula_id": self.q4_matrix_element_formula_id,
            "q5_matrix_element_formula_id": self.q5_matrix_element_formula_id,
            "chiral_ratio_definition_id": self.chiral_ratio_definition_id,
            "scheme_conversion_policy_id": self.scheme_conversion_policy_id,
            "input_policy_id": self.input_policy_id,
            "formula_source_id": self.formula_source_id,
            "formula_citation": self.formula_citation,
            "formula_locator_label": self.formula_locator_label,
            "notes": self.notes,
        }


@dataclass(frozen=True)
class Paper07101869KaonLRHadronicInputs:
    """Custom LR hadronic inputs and matrix elements in the paper O4/O5 basis."""

    B4_mu_had: float
    B5_mu_had: float
    R_chi_mu_had: float
    b4_source: Paper07101869HadronicSourceRef
    b5_source: Paper07101869HadronicSourceRef
    r_chi_source: Paper07101869HadronicSourceRef
    mass_source: Paper07101869HadronicSourceRef
    decay_constant_source: Paper07101869HadronicSourceRef
    m_K0_GeV: float = _DEFAULT_M_K0_GEV
    f_K_GeV: float = _DEFAULT_F_K_GEV
    contract: Paper07101869KaonLRHadronicContract | None = None
    system_id: str = PAPER_0710_1869_DELTAF2_KAON_HADRONIC_SYSTEM_ID
    operator_basis_id: str = PAPER_0710_1869_DELTAF2_OPERATOR_BASIS_ID
    operator_normalization_id: str = PAPER_0710_1869_DELTAF2_OPERATOR_NORMALIZATION_ID
    renormalization_scheme_id: str = PAPER_0710_1869_DELTAF2_RG_SCHEME_ID
    mu_had_GeV: float = PAPER_0710_1869_DELTAF2_DEFAULT_MU_HAD_GEV
    bundle_id: str = _CUSTOM_LR_HADRONIC_INPUTS_ID
    source_id: str = _CUSTOM_LR_HADRONIC_SOURCE_ID
    provenance_ids: tuple[str, ...] = (_CUSTOM_LR_HADRONIC_SOURCE_ID,)
    input_provenance_mode_id: str = (
        PAPER_0710_1869_DELTAF2_LR_HADRONIC_CUSTOM_PROVENANCE_MODE_ID
    )
    hamiltonian_convention_id: str = PAPER_0710_1869_DELTAF2_HAMILTONIAN_CONVENTION_ID
    q4_matrix_element_formula_id: str = (
        PAPER_0710_1869_DELTAF2_KAON_LR_Q4_MATRIX_ELEMENT_FORMULA_ID
    )
    q5_matrix_element_formula_id: str = (
        PAPER_0710_1869_DELTAF2_KAON_LR_Q5_MATRIX_ELEMENT_FORMULA_ID
    )
    chiral_ratio_definition_id: str = (
        PAPER_0710_1869_DELTAF2_KAON_LR_CHIRAL_RATIO_DEFINITION_ID
    )
    input_policy_id: str = PAPER_0710_1869_DELTAF2_KAON_LR_INPUT_POLICY_ID
    supported_operator_names: tuple[str, ...] = _SUPPORTED_LR_OPERATORS
    unsupported_operator_names: tuple[str, ...] = _UNSUPPORTED_LR_OPERATORS
    schema_id: str = PAPER_0710_1869_DELTAF2_KAON_LR_HADRONIC_INPUTS_SCHEMA_ID
    mode_id: str = PAPER_0710_1869_MODE_ID
    paper_id: str = PAPER_0710_1869_PAPER_ID
    notes: str = (
        "LR-HAD-1 custom-input-only bundle: O4/O5 scalar LR matrix elements from BV 2004 "
        "Eq. (5). The caller supplies B4, B5, and R_chi in the bundle's declared "
        "renormalization scheme and at mu_had; LR-HAD-1 performs no hidden scheme "
        "conversion. Defaults for LR hadronic inputs are not frozen yet, and LR "
        "observables remain blocked."
    )

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "schema_id",
            require_known_schema_id(
                "schema_id",
                self.schema_id,
                expected=PAPER_0710_1869_DELTAF2_KAON_LR_HADRONIC_INPUTS_SCHEMA_ID,
            ),
        )
        require_member("mode_id", self.mode_id, (PAPER_0710_1869_MODE_ID,))
        require_member("paper_id", self.paper_id, (PAPER_0710_1869_PAPER_ID,))
        require_member(
            "system_id",
            self.system_id,
            (PAPER_0710_1869_DELTAF2_KAON_HADRONIC_SYSTEM_ID,),
        )
        object.__setattr__(
            self,
            "operator_basis_id",
            require_member(
                "operator_basis_id",
                self.operator_basis_id,
                (PAPER_0710_1869_DELTAF2_OPERATOR_BASIS_ID,),
            ),
        )
        object.__setattr__(
            self,
            "operator_normalization_id",
            require_member(
                "operator_normalization_id",
                self.operator_normalization_id,
                (PAPER_0710_1869_DELTAF2_OPERATOR_NORMALIZATION_ID,),
            ),
        )
        object.__setattr__(
            self,
            "renormalization_scheme_id",
            require_nonempty_identifier(
                "renormalization_scheme_id",
                self.renormalization_scheme_id,
            ),
        )
        for source_name in (
            "b4_source",
            "b5_source",
            "r_chi_source",
            "mass_source",
            "decay_constant_source",
        ):
            if not isinstance(getattr(self, source_name), Paper07101869HadronicSourceRef):
                raise ValueError(f"{source_name} must be a Paper07101869HadronicSourceRef")
        resolved_contract = self.contract
        if resolved_contract is None:
            resolved_contract = Paper07101869KaonLRHadronicContract(
                system_id=self.system_id,
                operator_basis_id=self.operator_basis_id,
                operator_normalization_id=self.operator_normalization_id,
                renormalization_scheme_id=self.renormalization_scheme_id,
                mu_had_GeV=self.mu_had_GeV,
                evaluation_scale_GeV=self.mu_had_GeV,
                hamiltonian_convention_id=self.hamiltonian_convention_id,
                q4_matrix_element_formula_id=self.q4_matrix_element_formula_id,
                q5_matrix_element_formula_id=self.q5_matrix_element_formula_id,
                chiral_ratio_definition_id=self.chiral_ratio_definition_id,
                input_policy_id=self.input_policy_id,
            )
            object.__setattr__(self, "contract", resolved_contract)
        if not isinstance(resolved_contract, Paper07101869KaonLRHadronicContract):
            raise ValueError("contract must be a Paper07101869KaonLRHadronicContract")
        if resolved_contract.system_id != self.system_id:
            raise ValueError("contract.system_id must match the LR hadronic system_id")
        if resolved_contract.operator_basis_id != self.operator_basis_id:
            raise ValueError(
                "contract.operator_basis_id must match the LR hadronic operator_basis_id"
            )
        if resolved_contract.operator_normalization_id != self.operator_normalization_id:
            raise ValueError(
                "contract.operator_normalization_id must match the LR hadronic "
                "operator_normalization_id"
            )
        if resolved_contract.renormalization_scheme_id != self.renormalization_scheme_id:
            raise ValueError(
                "contract.renormalization_scheme_id must match the LR hadronic "
                "renormalization_scheme_id"
            )
        if not _float_matches(resolved_contract.mu_had_GeV, self.mu_had_GeV, atol=1e-12):
            raise ValueError("contract.mu_had_GeV must match the LR hadronic mu_had_GeV")
        if not _float_matches(
            resolved_contract.evaluation_scale_GeV,
            self.mu_had_GeV,
            atol=1e-12,
        ):
            raise ValueError(
                "contract.evaluation_scale_GeV must match the LR hadronic mu_had_GeV"
            )
        for field_name in (
            "bundle_id",
            "source_id",
            "input_provenance_mode_id",
            "hamiltonian_convention_id",
            "q4_matrix_element_formula_id",
            "q5_matrix_element_formula_id",
            "chiral_ratio_definition_id",
            "input_policy_id",
            "notes",
        ):
            object.__setattr__(
                self,
                field_name,
                require_nonempty_identifier(field_name, getattr(self, field_name)),
            )
        if self.bundle_id == _DEFAULT_HADRONIC_BUNDLE_ID:
            raise ValueError("LR hadronic inputs must not reuse the default kaon hadronic bundle")
        if self.source_id == _DEFAULT_HADRONIC_SOURCE_ID:
            raise ValueError("LR hadronic inputs must not reuse the default kaon hadronic source")
        if not self.provenance_ids:
            raise ValueError("provenance_ids must contain at least one source reference")
        object.__setattr__(
            self,
            "provenance_ids",
            tuple(
                require_nonempty_identifier("provenance_id", provenance_id)
                for provenance_id in self.provenance_ids
            ),
        )
        if self.source_id not in self.provenance_ids:
            raise ValueError("source_id must be present in provenance_ids for LR hadronic inputs")
        for field_name, source in (
            ("mass_source", self.mass_source),
            ("decay_constant_source", self.decay_constant_source),
            ("b4_source", self.b4_source),
            ("b5_source", self.b5_source),
            ("r_chi_source", self.r_chi_source),
        ):
            if source.source_id not in self.provenance_ids:
                raise ValueError(
                    f"{field_name}.source_id must be present in provenance_ids for LR "
                    "hadronic inputs"
                )
        require_member(
            "input_provenance_mode_id",
            self.input_provenance_mode_id,
            (PAPER_0710_1869_DELTAF2_LR_HADRONIC_CUSTOM_PROVENANCE_MODE_ID,),
        )
        for field_name in (
            "B4_mu_had",
            "B5_mu_had",
            "R_chi_mu_had",
            "m_K0_GeV",
            "f_K_GeV",
            "mu_had_GeV",
        ):
            object.__setattr__(
                self,
                field_name,
                require_positive_finite(field_name, getattr(self, field_name)),
            )
        for field_name, source in (
            ("b4_source", self.b4_source),
            ("b5_source", self.b5_source),
            ("r_chi_source", self.r_chi_source),
        ):
            if source.renormalization_scheme_id is None:
                raise ValueError(
                    f"{field_name} must carry renormalization_scheme_id matching the bundle "
                    "scheme"
                )
            if source.renormalization_scheme_id != self.renormalization_scheme_id:
                raise ValueError(
                    f"{field_name}.renormalization_scheme_id must match the LR hadronic "
                    "renormalization_scheme_id"
                )
            if source.scale_GeV is None:
                raise ValueError(f"{field_name} must carry scale_GeV matching mu_had_GeV")
            if not _float_matches(source.scale_GeV, self.mu_had_GeV, atol=1e-12):
                raise ValueError(f"{field_name}.scale_GeV must match mu_had_GeV")
        if tuple(self.supported_operator_names) != _SUPPORTED_LR_OPERATORS:
            raise ValueError("supported_operator_names must match the LR hadronic subset")
        if tuple(self.unsupported_operator_names) != _UNSUPPORTED_LR_OPERATORS:
            raise ValueError("unsupported_operator_names must match the non-LR held-back subset")

    @property
    def q4_matrix_element_GeV4(self) -> float:
        return _require_finite_float(
            "q4_matrix_element_GeV4",
            2.0 * (self.m_K0_GeV**2) * (self.f_K_GeV**2) * self.R_chi_mu_had * self.B4_mu_had,
        )

    @property
    def q5_matrix_element_GeV4(self) -> float:
        return _require_finite_float(
            "q5_matrix_element_GeV4",
            (2.0 / 3.0)
            * (self.m_K0_GeV**2)
            * (self.f_K_GeV**2)
            * self.R_chi_mu_had
            * self.B5_mu_had,
        )

    @property
    def tags(self) -> dict[str, object]:
        return {
            "mode_id": self.mode_id,
            "paper_id": self.paper_id,
            "system_id": self.system_id,
            "operator_basis_id": self.operator_basis_id,
            "operator_normalization_id": self.operator_normalization_id,
            "renormalization_scheme_id": self.renormalization_scheme_id,
            "contract": self.contract.as_dict(),
            "bundle_id": self.bundle_id,
            "source_id": self.source_id,
            "provenance_ids": list(self.provenance_ids),
            "input_provenance_mode_id": self.input_provenance_mode_id,
            "mu_had_GeV": self.mu_had_GeV,
            "evaluation_scale_GeV": self.contract.evaluation_scale_GeV,
            "supported_operator_names": list(self.supported_operator_names),
            "unsupported_operator_names": list(self.unsupported_operator_names),
            "q4_matrix_element_formula_id": self.q4_matrix_element_formula_id,
            "q5_matrix_element_formula_id": self.q5_matrix_element_formula_id,
            "chiral_ratio_definition_id": self.chiral_ratio_definition_id,
            "hamiltonian_convention_id": self.hamiltonian_convention_id,
            "input_policy_id": self.input_policy_id,
        }

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_id": self.schema_id,
            "mode_id": self.mode_id,
            "paper_id": self.paper_id,
            "system_id": self.system_id,
            "operator_basis_id": self.operator_basis_id,
            "operator_normalization_id": self.operator_normalization_id,
            "renormalization_scheme_id": self.renormalization_scheme_id,
            "scheme_id": self.contract.scheme_id,
            "contract": self.contract.as_dict(),
            "bundle_id": self.bundle_id,
            "source_id": self.source_id,
            "provenance_ids": list(self.provenance_ids),
            "input_provenance_mode_id": self.input_provenance_mode_id,
            "mu_had_GeV": self.mu_had_GeV,
            "evaluation_scale_GeV": self.contract.evaluation_scale_GeV,
            "m_K0_GeV": self.m_K0_GeV,
            "f_K_GeV": self.f_K_GeV,
            "B4_mu_had": self.B4_mu_had,
            "B5_mu_had": self.B5_mu_had,
            "R_chi_mu_had": self.R_chi_mu_had,
            "q4_matrix_element_formula_id": self.q4_matrix_element_formula_id,
            "q5_matrix_element_formula_id": self.q5_matrix_element_formula_id,
            "chiral_ratio_definition_id": self.chiral_ratio_definition_id,
            "input_policy_id": self.input_policy_id,
            "hamiltonian_convention_id": self.hamiltonian_convention_id,
            "supported_operator_names": list(self.supported_operator_names),
            "unsupported_operator_names": list(self.unsupported_operator_names),
            "mass_source": self.mass_source.as_dict(),
            "decay_constant_source": self.decay_constant_source.as_dict(),
            "b4_source": self.b4_source.as_dict(),
            "b5_source": self.b5_source.as_dict(),
            "r_chi_source": self.r_chi_source.as_dict(),
            "q4_matrix_element_GeV4": self.q4_matrix_element_GeV4,
            "q5_matrix_element_GeV4": self.q5_matrix_element_GeV4,
            "tags": self.tags,
            "notes": self.notes,
        }

    def summary(self) -> dict[str, object]:
        return {
            **self.as_dict(),
            "schema_id": PAPER_0710_1869_DELTAF2_KAON_LR_HADRONIC_SUMMARY_SCHEMA_ID,
        }


def _default_mass_source() -> Paper07101869HadronicSourceRef:
    return Paper07101869HadronicSourceRef(
        source_id=_DEFAULT_MASS_SOURCE_ID,
        source_kind="particle-property",
        citation="PDG 2024 neutral kaon listing",
        locator_label="K0 mass listing",
        year=2024,
        transformation_id="none",
        notes="Uses the PDG neutral-kaon mass m_K0 = 497.611 MeV.",
    )


def _default_decay_constant_source() -> Paper07101869HadronicSourceRef:
    return Paper07101869HadronicSourceRef(
        source_id=_DEFAULT_DECAY_CONSTANT_SOURCE_ID,
        source_kind="review-average",
        citation="PDG 2024 Leptonic Decays of Charged Pseudoscalar Mesons, Eq. (72.14)",
        locator_label="Eq. (72.14)",
        year=2024,
        transformation_id="charged_kaon_decay_constant_reused_as_fK.v1",
        notes=(
            "PR5a reuses the PDG/FLAG charged-kaon decay constant fK+ = 155.7 MeV "
            "as f_K in the kaon NP-only slice. No separate isospin correction is applied."
        ),
    )


def _default_bag_parameter_source() -> Paper07101869HadronicSourceRef:
    return Paper07101869HadronicSourceRef(
        source_id=_DEFAULT_BAG_SOURCE_ID,
        source_kind="derived-from-rgi",
        citation="PDG 2024 Lattice QCD review, Sec. 17.3.2, FLAG average",
        locator_label="Sec. 17.3.2, kaon bag parameter discussion",
        year=2024,
        renormalization_scheme_id=PAPER_0710_1869_DELTAF2_BAG_PARAMETER_SOURCE_SCHEME_ID,
        transformation_id=PAPER_0710_1869_DELTAF2_BAG_PARAMETER_TRANSFORMATION_ID,
        notes=(
            "Uses the FLAG average hat{B}_K = 0.7625 and converts it to B_K(mu_had) "
            "with the paper-mode LO Q1 inverse-running convention."
        ),
    )


def build_paper_0710_1869_kaon_lr_hadronic_inputs(
    *,
    B4_mu_had: float,
    B5_mu_had: float,
    R_chi_mu_had: float,
    b4_source: Paper07101869HadronicSourceRef,
    b5_source: Paper07101869HadronicSourceRef,
    r_chi_source: Paper07101869HadronicSourceRef,
    renormalization_scheme_id: str = PAPER_0710_1869_DELTAF2_RG_SCHEME_ID,
    mu_had_GeV: float = PAPER_0710_1869_DELTAF2_DEFAULT_MU_HAD_GEV,
    m_K0_GeV: float = _DEFAULT_M_K0_GEV,
    f_K_GeV: float = _DEFAULT_F_K_GEV,
    mass_source: Paper07101869HadronicSourceRef | None = None,
    decay_constant_source: Paper07101869HadronicSourceRef | None = None,
    bundle_id: str | None = None,
    source_id: str | None = None,
    provenance_ids: tuple[str, ...] | list[str] | None = None,
) -> Paper07101869KaonLRHadronicInputs:
    """Return the custom LR hadronic inputs in the paper O4/O5 scalar LR basis."""

    resolved_mu_had = require_positive_finite("mu_had_GeV", mu_had_GeV)
    resolved_renormalization_scheme_id = require_nonempty_identifier(
        "renormalization_scheme_id",
        renormalization_scheme_id,
    )
    custom_mass_value = not _float_matches(m_K0_GeV, _DEFAULT_M_K0_GEV)
    custom_decay_constant_value = not _float_matches(f_K_GeV, _DEFAULT_F_K_GEV)
    if custom_mass_value and mass_source is None:
        raise ValueError(
            "m_K0_GeV override requires a matching mass_source override so the LR hadronic "
            "builder does not emit stale default PDG provenance"
        )
    if custom_decay_constant_value and decay_constant_source is None:
        raise ValueError(
            "f_K_GeV override requires a matching decay_constant_source override so the LR "
            "hadronic builder does not emit stale default PDG/FLAG provenance"
        )
    resolved_mass_source = _default_mass_source() if mass_source is None else mass_source
    resolved_decay_constant_source = (
        _default_decay_constant_source()
        if decay_constant_source is None
        else decay_constant_source
    )
    resolved_bundle_id = (
        _CUSTOM_LR_HADRONIC_INPUTS_ID
        if bundle_id is None
        else require_nonempty_identifier("bundle_id", bundle_id)
    )
    resolved_source_id = (
        _CUSTOM_LR_HADRONIC_SOURCE_ID
        if source_id is None
        else require_nonempty_identifier("source_id", source_id)
    )
    if provenance_ids is None:
        resolved_provenance_ids = tuple(
            dict.fromkeys(
                (
                    resolved_source_id,
                    resolved_mass_source.source_id,
                    resolved_decay_constant_source.source_id,
                    b4_source.source_id,
                    b5_source.source_id,
                    r_chi_source.source_id,
                )
            )
        )
    else:
        resolved_provenance_ids = tuple(
            require_nonempty_identifier("provenance_id", provenance_id)
            for provenance_id in tuple(provenance_ids)
        )
    return Paper07101869KaonLRHadronicInputs(
        B4_mu_had=B4_mu_had,
        B5_mu_had=B5_mu_had,
        R_chi_mu_had=R_chi_mu_had,
        b4_source=b4_source,
        b5_source=b5_source,
        r_chi_source=r_chi_source,
        renormalization_scheme_id=resolved_renormalization_scheme_id,
        mu_had_GeV=resolved_mu_had,
        m_K0_GeV=m_K0_GeV,
        f_K_GeV=f_K_GeV,
        mass_source=resolved_mass_source,
        decay_constant_source=resolved_decay_constant_source,
        bundle_id=resolved_bundle_id,
        source_id=resolved_source_id,
        provenance_ids=resolved_provenance_ids,
    )


def build_paper_0710_1869_kaon_lr_hadronic(
    *,
    B4_mu_had: float,
    B5_mu_had: float,
    R_chi_mu_had: float,
    b4_source: Paper07101869HadronicSourceRef,
    b5_source: Paper07101869HadronicSourceRef,
    r_chi_source: Paper07101869HadronicSourceRef,
    renormalization_scheme_id: str = PAPER_0710_1869_DELTAF2_RG_SCHEME_ID,
    mu_had_GeV: float = PAPER_0710_1869_DELTAF2_DEFAULT_MU_HAD_GEV,
    m_K0_GeV: float = _DEFAULT_M_K0_GEV,
    f_K_GeV: float = _DEFAULT_F_K_GEV,
    mass_source: Paper07101869HadronicSourceRef | None = None,
    decay_constant_source: Paper07101869HadronicSourceRef | None = None,
    bundle_id: str | None = None,
    source_id: str | None = None,
    provenance_ids: tuple[str, ...] | list[str] | None = None,
) -> Paper07101869KaonLRHadronicInputs:
    """Compatibility alias for the custom LR hadronic input bundle."""

    return build_paper_0710_1869_kaon_lr_hadronic_inputs(
        B4_mu_had=B4_mu_had,
        B5_mu_had=B5_mu_had,
        R_chi_mu_had=R_chi_mu_had,
        b4_source=b4_source,
        b5_source=b5_source,
        r_chi_source=r_chi_source,
        renormalization_scheme_id=renormalization_scheme_id,
        mu_had_GeV=mu_had_GeV,
        m_K0_GeV=m_K0_GeV,
        f_K_GeV=f_K_GeV,
        mass_source=mass_source,
        decay_constant_source=decay_constant_source,
        bundle_id=bundle_id,
        source_id=source_id,
        provenance_ids=provenance_ids,
    )


def build_paper_0710_1869_kaon_lr_hadronic_summary(
    *,
    B4_mu_had: float,
    B5_mu_had: float,
    R_chi_mu_had: float,
    b4_source: Paper07101869HadronicSourceRef,
    b5_source: Paper07101869HadronicSourceRef,
    r_chi_source: Paper07101869HadronicSourceRef,
    renormalization_scheme_id: str = PAPER_0710_1869_DELTAF2_RG_SCHEME_ID,
    mu_had_GeV: float = PAPER_0710_1869_DELTAF2_DEFAULT_MU_HAD_GEV,
    m_K0_GeV: float = _DEFAULT_M_K0_GEV,
    f_K_GeV: float = _DEFAULT_F_K_GEV,
    mass_source: Paper07101869HadronicSourceRef | None = None,
    decay_constant_source: Paper07101869HadronicSourceRef | None = None,
    bundle_id: str | None = None,
    source_id: str | None = None,
    provenance_ids: tuple[str, ...] | list[str] | None = None,
) -> dict[str, object]:
    """Return the deterministic summary mapping for the custom LR hadronic inputs."""

    return build_paper_0710_1869_kaon_lr_hadronic_inputs(
        B4_mu_had=B4_mu_had,
        B5_mu_had=B5_mu_had,
        R_chi_mu_had=R_chi_mu_had,
        b4_source=b4_source,
        b5_source=b5_source,
        r_chi_source=r_chi_source,
        renormalization_scheme_id=renormalization_scheme_id,
        mu_had_GeV=mu_had_GeV,
        m_K0_GeV=m_K0_GeV,
        f_K_GeV=f_K_GeV,
        mass_source=mass_source,
        decay_constant_source=decay_constant_source,
        bundle_id=bundle_id,
        source_id=source_id,
        provenance_ids=provenance_ids,
    ).summary()


def build_paper_0710_1869_kaon_hadronic_bundle(
    *,
    mu_had_GeV: float = PAPER_0710_1869_DELTAF2_DEFAULT_MU_HAD_GEV,
    m_K0_GeV: float = _DEFAULT_M_K0_GEV,
    f_K_GeV: float = _DEFAULT_F_K_GEV,
    B_K_mu_had: float | None = None,
    hat_B_K_rgi_source_value: float = _DEFAULT_HAT_B_K_RGI,
    mass_source: Paper07101869HadronicSourceRef | None = None,
    decay_constant_source: Paper07101869HadronicSourceRef | None = None,
    bag_parameter_source: Paper07101869HadronicSourceRef | None = None,
    bundle_id: str | None = None,
    source_id: str | None = None,
    provenance_ids: tuple[str, ...] | list[str] | None = None,
    bag_parameter_source_scheme_id: str | None = None,
    bag_parameter_transformation_id: str | None = None,
) -> Paper07101869KaonHadronicBundle:
    """Return the deterministic kaon hadronic bundle for the supported PR5a subset."""

    resolved_mu_had = require_positive_finite("mu_had_GeV", mu_had_GeV)
    resolved_hat_bk = require_positive_finite(
        "hat_B_K_rgi_source_value", hat_B_K_rgi_source_value
    )
    default_b_k_from_defaults = _derive_b_k_mu_had(
        mu_had_GeV=resolved_mu_had,
        hat_B_K_rgi=_DEFAULT_HAT_B_K_RGI,
    )
    resolved_b_k = (
        _derive_b_k_mu_had(
            mu_had_GeV=resolved_mu_had,
            hat_B_K_rgi=resolved_hat_bk,
        )
        if B_K_mu_had is None
        else require_positive_finite("B_K_mu_had", B_K_mu_had)
    )
    uses_derived_b_k_from_hat = B_K_mu_had is None
    custom_mass_value = not _float_matches(m_K0_GeV, _DEFAULT_M_K0_GEV)
    custom_decay_constant_value = not _float_matches(f_K_GeV, _DEFAULT_F_K_GEV)
    custom_hat_bk_value = not _float_matches(resolved_hat_bk, _DEFAULT_HAT_B_K_RGI)
    custom_b_k_value = B_K_mu_had is not None and not _float_matches(
        resolved_b_k,
        default_b_k_from_defaults,
    )
    if custom_mass_value and mass_source is None:
        raise ValueError(
            "m_K0_GeV override requires a matching mass_source override so the builder "
            "does not emit stale default PDG provenance"
        )
    if custom_decay_constant_value and decay_constant_source is None:
        raise ValueError(
            "f_K_GeV override requires a matching decay_constant_source override so the "
            "builder does not emit stale default PDG/FLAG provenance"
        )
    if (custom_hat_bk_value or custom_b_k_value) and bag_parameter_source is None:
        raise ValueError(
            "B_K_mu_had or hat_B_K_rgi_source_value override requires a matching "
            "bag_parameter_source override so the builder does not emit stale "
            "default FLAG provenance"
        )
    custom_provenance_requested = any(
        item is not None
        for item in (
            mass_source,
            decay_constant_source,
            bag_parameter_source,
            bundle_id,
            source_id,
            provenance_ids,
            bag_parameter_source_scheme_id,
            bag_parameter_transformation_id,
        )
    ) or any(
        (
            custom_mass_value,
            custom_decay_constant_value,
            custom_hat_bk_value,
            custom_b_k_value,
        )
    )
    if provenance_ids is None:
        normalized_provenance_ids: tuple[str, ...] | None = None
    else:
        normalized_provenance_ids = tuple(
            require_nonempty_identifier("provenance_id", provenance_id)
            for provenance_id in tuple(provenance_ids)
        )
    resolved_mass_source = _default_mass_source() if mass_source is None else mass_source
    resolved_decay_constant_source = (
        _default_decay_constant_source()
        if decay_constant_source is None
        else decay_constant_source
    )
    resolved_bag_parameter_source = (
        _default_bag_parameter_source()
        if bag_parameter_source is None
        else bag_parameter_source
    )
    resolved_bag_parameter_source_scheme_id = (
        PAPER_0710_1869_DELTAF2_BAG_PARAMETER_SOURCE_SCHEME_ID
        if bag_parameter_source_scheme_id is None and bag_parameter_source is None
        else (
            resolved_bag_parameter_source.renormalization_scheme_id
            if bag_parameter_source_scheme_id is None
            else require_nonempty_identifier(
                "bag_parameter_source_scheme_id", bag_parameter_source_scheme_id
            )
        )
    )
    if resolved_bag_parameter_source_scheme_id is None:
        raise ValueError(
            "bag_parameter_source must carry renormalization_scheme_id or the builder must "
            "receive bag_parameter_source_scheme_id explicitly"
        )
    resolved_bag_parameter_transformation_id = (
        PAPER_0710_1869_DELTAF2_BAG_PARAMETER_TRANSFORMATION_ID
        if bag_parameter_transformation_id is None and bag_parameter_source is None
        else (
            resolved_bag_parameter_source.transformation_id
            if bag_parameter_transformation_id is None
            else require_nonempty_identifier(
                "bag_parameter_transformation_id", bag_parameter_transformation_id
            )
        )
    )
    if uses_derived_b_k_from_hat:
        require_member(
            "bag_parameter_source_scheme_id",
            resolved_bag_parameter_source_scheme_id,
            (PAPER_0710_1869_DELTAF2_BAG_PARAMETER_SOURCE_SCHEME_ID,),
        )
        require_member(
            "bag_parameter_transformation_id",
            resolved_bag_parameter_transformation_id,
            (PAPER_0710_1869_DELTAF2_BAG_PARAMETER_TRANSFORMATION_ID,),
        )
    if custom_provenance_requested:
        if bundle_id == _DEFAULT_HADRONIC_BUNDLE_ID:
            raise ValueError("custom hadronic inputs must not reuse the default bundle_id")
        if source_id == _DEFAULT_HADRONIC_SOURCE_ID:
            raise ValueError("custom hadronic inputs must not reuse the default source_id")
        if normalized_provenance_ids == (_DEFAULT_HADRONIC_SOURCE_ID,):
            raise ValueError(
                "custom hadronic inputs must not reuse the default provenance_ids"
            )
        resolved_bundle_id = (
            _CUSTOM_HADRONIC_BUNDLE_ID
            if bundle_id is None
            else require_nonempty_identifier("bundle_id", bundle_id)
        )
        resolved_source_id = (
            _CUSTOM_HADRONIC_SOURCE_ID
            if source_id is None
            else require_nonempty_identifier("source_id", source_id)
        )
        resolved_provenance_ids = (
            (resolved_source_id,)
            if normalized_provenance_ids is None
            else normalized_provenance_ids
        )
        resolved_input_provenance_mode_id = (
            PAPER_0710_1869_DELTAF2_HADRONIC_CUSTOM_PROVENANCE_MODE_ID
        )
    else:
        resolved_bundle_id = _DEFAULT_HADRONIC_BUNDLE_ID
        resolved_source_id = _DEFAULT_HADRONIC_SOURCE_ID
        resolved_provenance_ids = (_DEFAULT_HADRONIC_SOURCE_ID,)
        resolved_input_provenance_mode_id = (
            PAPER_0710_1869_DELTAF2_HADRONIC_DEFAULT_PROVENANCE_MODE_ID
        )
    return Paper07101869KaonHadronicBundle(
        mu_had_GeV=resolved_mu_had,
        m_K0_GeV=m_K0_GeV,
        f_K_GeV=f_K_GeV,
        B_K_mu_had=resolved_b_k,
        bundle_id=resolved_bundle_id,
        source_id=resolved_source_id,
        provenance_ids=resolved_provenance_ids,
        input_provenance_mode_id=resolved_input_provenance_mode_id,
        hat_B_K_rgi_source_value=resolved_hat_bk,
        bag_parameter_source_scheme_id=resolved_bag_parameter_source_scheme_id,
        bag_parameter_transformation_id=resolved_bag_parameter_transformation_id,
        mass_source=resolved_mass_source,
        decay_constant_source=resolved_decay_constant_source,
        bag_parameter_source=resolved_bag_parameter_source,
    )


def default_paper_0710_1869_kaon_hadronic_bundle() -> Paper07101869KaonHadronicBundle:
    """Return the default kaon hadronic bundle at the frozen paper ``mu_had``."""

    return build_paper_0710_1869_kaon_hadronic_bundle()


def build_paper_0710_1869_kaon_hadronic(
    *,
    mu_had_GeV: float = PAPER_0710_1869_DELTAF2_DEFAULT_MU_HAD_GEV,
    m_K0_GeV: float = _DEFAULT_M_K0_GEV,
    f_K_GeV: float = _DEFAULT_F_K_GEV,
    B_K_mu_had: float | None = None,
    hat_B_K_rgi_source_value: float = _DEFAULT_HAT_B_K_RGI,
    mass_source: Paper07101869HadronicSourceRef | None = None,
    decay_constant_source: Paper07101869HadronicSourceRef | None = None,
    bag_parameter_source: Paper07101869HadronicSourceRef | None = None,
    bundle_id: str | None = None,
    source_id: str | None = None,
    provenance_ids: tuple[str, ...] | list[str] | None = None,
    bag_parameter_source_scheme_id: str | None = None,
    bag_parameter_transformation_id: str | None = None,
) -> Paper07101869KaonHadronicBundle:
    """Compatibility alias for the deterministic kaon hadronic bundle."""

    return build_paper_0710_1869_kaon_hadronic_bundle(
        mu_had_GeV=mu_had_GeV,
        m_K0_GeV=m_K0_GeV,
        f_K_GeV=f_K_GeV,
        B_K_mu_had=B_K_mu_had,
        hat_B_K_rgi_source_value=hat_B_K_rgi_source_value,
        mass_source=mass_source,
        decay_constant_source=decay_constant_source,
        bag_parameter_source=bag_parameter_source,
        bundle_id=bundle_id,
        source_id=source_id,
        provenance_ids=provenance_ids,
        bag_parameter_source_scheme_id=bag_parameter_source_scheme_id,
        bag_parameter_transformation_id=bag_parameter_transformation_id,
    )


def default_paper_0710_1869_kaon_hadronic() -> Paper07101869KaonHadronicBundle:
    """Compatibility alias for the default kaon hadronic bundle."""

    return default_paper_0710_1869_kaon_hadronic_bundle()


def build_paper_0710_1869_kaon_hadronic_summary(
    *,
    mu_had_GeV: float = PAPER_0710_1869_DELTAF2_DEFAULT_MU_HAD_GEV,
    m_K0_GeV: float = _DEFAULT_M_K0_GEV,
    f_K_GeV: float = _DEFAULT_F_K_GEV,
    B_K_mu_had: float | None = None,
    hat_B_K_rgi_source_value: float = _DEFAULT_HAT_B_K_RGI,
    mass_source: Paper07101869HadronicSourceRef | None = None,
    decay_constant_source: Paper07101869HadronicSourceRef | None = None,
    bag_parameter_source: Paper07101869HadronicSourceRef | None = None,
    bundle_id: str | None = None,
    source_id: str | None = None,
    provenance_ids: tuple[str, ...] | list[str] | None = None,
    bag_parameter_source_scheme_id: str | None = None,
    bag_parameter_transformation_id: str | None = None,
) -> dict[str, object]:
    """Return the deterministic summary mapping for the default kaon hadronic bundle."""

    return build_paper_0710_1869_kaon_hadronic_bundle(
        mu_had_GeV=mu_had_GeV,
        m_K0_GeV=m_K0_GeV,
        f_K_GeV=f_K_GeV,
        B_K_mu_had=B_K_mu_had,
        hat_B_K_rgi_source_value=hat_B_K_rgi_source_value,
        mass_source=mass_source,
        decay_constant_source=decay_constant_source,
        bag_parameter_source=bag_parameter_source,
        bundle_id=bundle_id,
        source_id=source_id,
        provenance_ids=provenance_ids,
        bag_parameter_source_scheme_id=bag_parameter_source_scheme_id,
        bag_parameter_transformation_id=bag_parameter_transformation_id,
    ).summary()


def default_paper_0710_1869_kaon_hadronic_summary() -> dict[str, object]:
    """Return the summary mapping for the default kaon hadronic bundle."""

    return default_paper_0710_1869_kaon_hadronic_bundle().summary()


def build_paper_0710_1869_kaon_hadronic_summary_alias(
    *,
    mu_had_GeV: float = PAPER_0710_1869_DELTAF2_DEFAULT_MU_HAD_GEV,
    m_K0_GeV: float = _DEFAULT_M_K0_GEV,
    f_K_GeV: float = _DEFAULT_F_K_GEV,
    B_K_mu_had: float | None = None,
    hat_B_K_rgi_source_value: float = _DEFAULT_HAT_B_K_RGI,
    mass_source: Paper07101869HadronicSourceRef | None = None,
    decay_constant_source: Paper07101869HadronicSourceRef | None = None,
    bag_parameter_source: Paper07101869HadronicSourceRef | None = None,
    bundle_id: str | None = None,
    source_id: str | None = None,
    provenance_ids: tuple[str, ...] | list[str] | None = None,
    bag_parameter_source_scheme_id: str | None = None,
    bag_parameter_transformation_id: str | None = None,
) -> dict[str, object]:
    """Compatibility alias for the deterministic kaon hadronic summary."""

    return build_paper_0710_1869_kaon_hadronic_summary(
        mu_had_GeV=mu_had_GeV,
        m_K0_GeV=m_K0_GeV,
        f_K_GeV=f_K_GeV,
        B_K_mu_had=B_K_mu_had,
        hat_B_K_rgi_source_value=hat_B_K_rgi_source_value,
        mass_source=mass_source,
        decay_constant_source=decay_constant_source,
        bag_parameter_source=bag_parameter_source,
        bundle_id=bundle_id,
        source_id=source_id,
        provenance_ids=provenance_ids,
        bag_parameter_source_scheme_id=bag_parameter_source_scheme_id,
        bag_parameter_transformation_id=bag_parameter_transformation_id,
    )


def _build_paper_0710_1869_b_hadronic_bundle(
    *,
    system_id: str,
    mu_had_GeV: float,
    meson_mass_GeV: float,
    meson_decay_constant_GeV: float,
    bag_parameter_mu_had: float,
    renormalization_scheme_id: str,
    mass_source: Paper07101869HadronicSourceRef,
    decay_constant_source: Paper07101869HadronicSourceRef,
    bag_parameter_source: Paper07101869HadronicSourceRef,
    bundle_id: str | None = None,
    source_id: str | None = None,
    provenance_ids: tuple[str, ...] | list[str] | None = None,
) -> Paper07101869BMesonHadronicBundle:
    resolved_system_id = _require_custom_b_system_id("system_id", system_id)
    resolved_provenance_ids = None if provenance_ids is None else tuple(provenance_ids)
    return Paper07101869BMesonHadronicBundle(
        system_id=resolved_system_id,
        meson_mass_GeV=meson_mass_GeV,
        meson_decay_constant_GeV=meson_decay_constant_GeV,
        bag_parameter_mu_had=bag_parameter_mu_had,
        renormalization_scheme_id=renormalization_scheme_id,
        mu_had_GeV=mu_had_GeV,
        mass_source=mass_source,
        decay_constant_source=decay_constant_source,
        bag_parameter_source=bag_parameter_source,
        bundle_id=bundle_id,
        source_id=source_id,
        provenance_ids=resolved_provenance_ids,
    )


def build_paper_0710_1869_bd_hadronic_bundle(
    *,
    mu_had_GeV: float,
    renormalization_scheme_id: str,
    m_Bd_GeV: float,
    f_Bd_GeV: float,
    B_Bd_mu_had: float,
    mass_source: Paper07101869HadronicSourceRef,
    decay_constant_source: Paper07101869HadronicSourceRef,
    bag_parameter_source: Paper07101869HadronicSourceRef,
    bundle_id: str | None = None,
    source_id: str | None = None,
    provenance_ids: tuple[str, ...] | list[str] | None = None,
) -> Paper07101869BMesonHadronicBundle:
    """Build the custom-input-only B_d Q1 hadronic bundle."""

    return _build_paper_0710_1869_b_hadronic_bundle(
        system_id=PAPER_0710_1869_DELTAF2_BD_HADRONIC_SYSTEM_ID,
        mu_had_GeV=mu_had_GeV,
        renormalization_scheme_id=renormalization_scheme_id,
        meson_mass_GeV=m_Bd_GeV,
        meson_decay_constant_GeV=f_Bd_GeV,
        bag_parameter_mu_had=B_Bd_mu_had,
        mass_source=mass_source,
        decay_constant_source=decay_constant_source,
        bag_parameter_source=bag_parameter_source,
        bundle_id=bundle_id,
        source_id=source_id,
        provenance_ids=provenance_ids,
    )


def build_paper_0710_1869_bd_hadronic(
    *,
    mu_had_GeV: float,
    renormalization_scheme_id: str,
    m_Bd_GeV: float,
    f_Bd_GeV: float,
    B_Bd_mu_had: float,
    mass_source: Paper07101869HadronicSourceRef,
    decay_constant_source: Paper07101869HadronicSourceRef,
    bag_parameter_source: Paper07101869HadronicSourceRef,
    bundle_id: str | None = None,
    source_id: str | None = None,
    provenance_ids: tuple[str, ...] | list[str] | None = None,
) -> Paper07101869BMesonHadronicBundle:
    """Compatibility alias for the custom-input-only B_d Q1 hadronic bundle."""

    return build_paper_0710_1869_bd_hadronic_bundle(
        mu_had_GeV=mu_had_GeV,
        renormalization_scheme_id=renormalization_scheme_id,
        m_Bd_GeV=m_Bd_GeV,
        f_Bd_GeV=f_Bd_GeV,
        B_Bd_mu_had=B_Bd_mu_had,
        mass_source=mass_source,
        decay_constant_source=decay_constant_source,
        bag_parameter_source=bag_parameter_source,
        bundle_id=bundle_id,
        source_id=source_id,
        provenance_ids=provenance_ids,
    )


def build_paper_0710_1869_bd_hadronic_summary(
    *,
    mu_had_GeV: float,
    renormalization_scheme_id: str,
    m_Bd_GeV: float,
    f_Bd_GeV: float,
    B_Bd_mu_had: float,
    mass_source: Paper07101869HadronicSourceRef,
    decay_constant_source: Paper07101869HadronicSourceRef,
    bag_parameter_source: Paper07101869HadronicSourceRef,
    bundle_id: str | None = None,
    source_id: str | None = None,
    provenance_ids: tuple[str, ...] | list[str] | None = None,
) -> dict[str, object]:
    """Return the summary mapping for the custom-input-only B_d Q1 hadronic bundle."""

    return build_paper_0710_1869_bd_hadronic_bundle(
        mu_had_GeV=mu_had_GeV,
        renormalization_scheme_id=renormalization_scheme_id,
        m_Bd_GeV=m_Bd_GeV,
        f_Bd_GeV=f_Bd_GeV,
        B_Bd_mu_had=B_Bd_mu_had,
        mass_source=mass_source,
        decay_constant_source=decay_constant_source,
        bag_parameter_source=bag_parameter_source,
        bundle_id=bundle_id,
        source_id=source_id,
        provenance_ids=provenance_ids,
    ).summary()


def build_paper_0710_1869_bs_hadronic_bundle(
    *,
    mu_had_GeV: float,
    renormalization_scheme_id: str,
    m_Bs_GeV: float,
    f_Bs_GeV: float,
    B_Bs_mu_had: float,
    mass_source: Paper07101869HadronicSourceRef,
    decay_constant_source: Paper07101869HadronicSourceRef,
    bag_parameter_source: Paper07101869HadronicSourceRef,
    bundle_id: str | None = None,
    source_id: str | None = None,
    provenance_ids: tuple[str, ...] | list[str] | None = None,
) -> Paper07101869BMesonHadronicBundle:
    """Build the custom-input-only B_s Q1 hadronic bundle."""

    return _build_paper_0710_1869_b_hadronic_bundle(
        system_id=PAPER_0710_1869_DELTAF2_BS_HADRONIC_SYSTEM_ID,
        mu_had_GeV=mu_had_GeV,
        renormalization_scheme_id=renormalization_scheme_id,
        meson_mass_GeV=m_Bs_GeV,
        meson_decay_constant_GeV=f_Bs_GeV,
        bag_parameter_mu_had=B_Bs_mu_had,
        mass_source=mass_source,
        decay_constant_source=decay_constant_source,
        bag_parameter_source=bag_parameter_source,
        bundle_id=bundle_id,
        source_id=source_id,
        provenance_ids=provenance_ids,
    )


def build_paper_0710_1869_bs_hadronic(
    *,
    mu_had_GeV: float,
    renormalization_scheme_id: str,
    m_Bs_GeV: float,
    f_Bs_GeV: float,
    B_Bs_mu_had: float,
    mass_source: Paper07101869HadronicSourceRef,
    decay_constant_source: Paper07101869HadronicSourceRef,
    bag_parameter_source: Paper07101869HadronicSourceRef,
    bundle_id: str | None = None,
    source_id: str | None = None,
    provenance_ids: tuple[str, ...] | list[str] | None = None,
) -> Paper07101869BMesonHadronicBundle:
    """Compatibility alias for the custom-input-only B_s Q1 hadronic bundle."""

    return build_paper_0710_1869_bs_hadronic_bundle(
        mu_had_GeV=mu_had_GeV,
        renormalization_scheme_id=renormalization_scheme_id,
        m_Bs_GeV=m_Bs_GeV,
        f_Bs_GeV=f_Bs_GeV,
        B_Bs_mu_had=B_Bs_mu_had,
        mass_source=mass_source,
        decay_constant_source=decay_constant_source,
        bag_parameter_source=bag_parameter_source,
        bundle_id=bundle_id,
        source_id=source_id,
        provenance_ids=provenance_ids,
    )


def build_paper_0710_1869_bs_hadronic_summary(
    *,
    mu_had_GeV: float,
    renormalization_scheme_id: str,
    m_Bs_GeV: float,
    f_Bs_GeV: float,
    B_Bs_mu_had: float,
    mass_source: Paper07101869HadronicSourceRef,
    decay_constant_source: Paper07101869HadronicSourceRef,
    bag_parameter_source: Paper07101869HadronicSourceRef,
    bundle_id: str | None = None,
    source_id: str | None = None,
    provenance_ids: tuple[str, ...] | list[str] | None = None,
) -> dict[str, object]:
    """Return the summary mapping for the custom-input-only B_s Q1 hadronic bundle."""

    return build_paper_0710_1869_bs_hadronic_bundle(
        mu_had_GeV=mu_had_GeV,
        renormalization_scheme_id=renormalization_scheme_id,
        m_Bs_GeV=m_Bs_GeV,
        f_Bs_GeV=f_Bs_GeV,
        B_Bs_mu_had=B_Bs_mu_had,
        mass_source=mass_source,
        decay_constant_source=decay_constant_source,
        bag_parameter_source=bag_parameter_source,
        bundle_id=bundle_id,
        source_id=source_id,
        provenance_ids=provenance_ids,
    ).summary()


__all__ = [
    "PAPER_0710_1869_DELTAF2_BAG_PARAMETER_SOURCE_SCHEME_ID",
    "PAPER_0710_1869_DELTAF2_BAG_PARAMETER_TRANSFORMATION_ID",
    "PAPER_0710_1869_DELTAF2_BD_HADRONIC_BUNDLE_SCHEMA_ID",
    "PAPER_0710_1869_DELTAF2_BD_HADRONIC_CONTRACT_SCHEMA_ID",
    "PAPER_0710_1869_DELTAF2_BD_HADRONIC_SUMMARY_SCHEMA_ID",
    "PAPER_0710_1869_DELTAF2_BD_HADRONIC_SYSTEM_ID",
    "PAPER_0710_1869_DELTAF2_BD_MATRIX_ELEMENT_FORMULA_ID",
    "PAPER_0710_1869_DELTAF2_BD_PARITY_RELATION_ID",
    "PAPER_0710_1869_DELTAF2_BS_HADRONIC_BUNDLE_SCHEMA_ID",
    "PAPER_0710_1869_DELTAF2_BS_HADRONIC_CONTRACT_SCHEMA_ID",
    "PAPER_0710_1869_DELTAF2_BS_HADRONIC_SUMMARY_SCHEMA_ID",
    "PAPER_0710_1869_DELTAF2_BS_HADRONIC_SYSTEM_ID",
    "PAPER_0710_1869_DELTAF2_BS_MATRIX_ELEMENT_FORMULA_ID",
    "PAPER_0710_1869_DELTAF2_BS_PARITY_RELATION_ID",
    "PAPER_0710_1869_DELTAF2_HADRONIC_SOURCE_SCHEMA_ID",
    "PAPER_0710_1869_DELTAF2_HAMILTONIAN_CONVENTION_ID",
    "PAPER_0710_1869_DELTAF2_KAON_HADRONIC_BUNDLE_SCHEMA_ID",
    "PAPER_0710_1869_DELTAF2_KAON_HADRONIC_CONTRACT_SCHEMA_ID",
    "PAPER_0710_1869_DELTAF2_KAON_HADRONIC_SUMMARY_SCHEMA_ID",
    "PAPER_0710_1869_DELTAF2_KAON_HADRONIC_SYSTEM_ID",
    "PAPER_0710_1869_DELTAF2_KAON_LR_CHIRAL_RATIO_DEFINITION_ID",
    "PAPER_0710_1869_DELTAF2_KAON_LR_FORMULA_SOURCE_ID",
    "PAPER_0710_1869_DELTAF2_KAON_LR_HADRONIC_CONTRACT_SCHEMA_ID",
    "PAPER_0710_1869_DELTAF2_KAON_LR_HADRONIC_INPUTS_SCHEMA_ID",
    "PAPER_0710_1869_DELTAF2_KAON_LR_HADRONIC_SUMMARY_SCHEMA_ID",
    "PAPER_0710_1869_DELTAF2_KAON_LR_INPUT_POLICY_ID",
    "PAPER_0710_1869_DELTAF2_KAON_LR_Q4_MATRIX_ELEMENT_FORMULA_ID",
    "PAPER_0710_1869_DELTAF2_KAON_LR_Q5_MATRIX_ELEMENT_FORMULA_ID",
    "PAPER_0710_1869_DELTAF2_KAON_LR_SCHEME_CONVERSION_POLICY_ID",
    "PAPER_0710_1869_DELTAF2_KAON_MATRIX_ELEMENT_FORMULA_ID",
    "PAPER_0710_1869_DELTAF2_KAON_PARITY_RELATION_ID",
    "PAPER_0710_1869_DELTAF2_LR_HADRONIC_CUSTOM_PROVENANCE_MODE_ID",
    "Paper07101869BMesonHadronicBundle",
    "Paper07101869BMesonHadronicContract",
    "Paper07101869HadronicSourceRef",
    "Paper07101869KaonHadronicContract",
    "Paper07101869KaonHadronicBundle",
    "Paper07101869KaonLRHadronicContract",
    "Paper07101869KaonLRHadronicInputs",
    "build_paper_0710_1869_bd_hadronic",
    "build_paper_0710_1869_bd_hadronic_bundle",
    "build_paper_0710_1869_bd_hadronic_summary",
    "build_paper_0710_1869_bs_hadronic",
    "build_paper_0710_1869_bs_hadronic_bundle",
    "build_paper_0710_1869_bs_hadronic_summary",
    "build_paper_0710_1869_kaon_hadronic",
    "build_paper_0710_1869_kaon_hadronic_bundle",
    "build_paper_0710_1869_kaon_hadronic_summary",
    "build_paper_0710_1869_kaon_hadronic_summary_alias",
    "build_paper_0710_1869_kaon_lr_hadronic",
    "build_paper_0710_1869_kaon_lr_hadronic_inputs",
    "build_paper_0710_1869_kaon_lr_hadronic_summary",
    "default_paper_0710_1869_kaon_hadronic",
    "default_paper_0710_1869_kaon_hadronic_bundle",
    "default_paper_0710_1869_kaon_hadronic_summary",
]
