"""Minimal paper-owned ``Delta F = 2`` operator and Wilson contracts."""

from __future__ import annotations

import math
from dataclasses import dataclass

from ..conventions import PAPER_0710_1869_MODE_ID, PAPER_0710_1869_PAPER_ID
from ..couplings import (
    PAPER_0710_1869_DIMENSIONLESS_MATRIX_POLICY_ID,
    PAPER_0710_1869_GS_NORMALIZATION_ID,
    PAPER_0710_1869_PROPAGATOR_MASS_RULE_ID,
    PAPER_0710_1869_UNIVERSAL_SUBTRACTION_POLICY_ID,
)
from ..validation import (
    require_known_schema_id,
    require_member,
    require_nonempty_identifier,
    require_positive_finite,
)

PAPER_0710_1869_DELTAF2_OPERATOR_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.eft_deltaf2.operator.v1"
)
PAPER_0710_1869_DELTAF2_OPERATOR_BASIS_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.eft_deltaf2.operator_basis.v1"
)
PAPER_0710_1869_DELTAF2_WILSON_CONTRACT_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.eft_deltaf2.wilson_contract.v1"
)
PAPER_0710_1869_DELTAF2_WILSON_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.eft_deltaf2.wilson_coefficients.v1"
)
PAPER_0710_1869_DELTAF2_OPERATOR_DEFINITION_Q1_VLL_ID = (
    "paper_0710_1869.deltaf2.q1_vll.current_current.color_singlet.v1"
)
PAPER_0710_1869_DELTAF2_OPERATOR_DEFINITION_Q1_VRR_ID = (
    "paper_0710_1869.deltaf2.q1_vrr.current_current.color_singlet.v1"
)
PAPER_0710_1869_DELTAF2_OPERATOR_DEFINITION_Q4_LR_ID = (
    "paper_0710_1869.deltaf2.q4_lr.susy_o4.scalar_lr.color_singlet.v1"
)
PAPER_0710_1869_DELTAF2_OPERATOR_DEFINITION_Q5_LR_ID = (
    "paper_0710_1869.deltaf2.q5_lr.susy_o5.scalar_lr.color_mixed.v1"
)
PAPER_0710_1869_DELTAF2_PROJECTOR_NORMALIZATION_ID = (
    "paper_0710_1869.deltaf2.projectors.pl_pr_equals_1mp_gamma5_over_2.v1"
)
PAPER_0710_1869_DELTAF2_PROJECTOR_NORMALIZATION_NOTE = (
    "All chiral bilinears use P_L = (1-gamma5)/2 and P_R = (1+gamma5)/2. "
    "The legacy susy_o4/susy_o5 identifiers are retained for compatibility only. "
    "They label the O4/O5 scalar LR basis used in the DeltaF=2 literature and do "
    "not imply a bare (1-gamma5) or (1+gamma5) normalization without the "
    "factor-of-two projector convention."
)
PAPER_0710_1869_DELTAF2_OPERATOR_METADATA_COMPATIBILITY_NOTE = (
    "This LR-freeze slice is a metadata clarification only. It does not change the "
    "frozen operator ordering, projector normalization convention already used by "
    "the formulas, supported Q1-only observable subset, or any matching/RG/"
    "observable numeric behavior. Legacy susy_o4/susy_o5 identifiers are preserved "
    "as stable machine ids and are not statements about the underlying mediator "
    "sector."
)

# This matches the frozen paper-mode conventions bundle while the explicit
# four-operator content lives in the basis object below.
PAPER_0710_1869_DELTAF2_OPERATOR_BASIS_ID = "kk_gluon_tree_np_only.v1"
PAPER_0710_1869_DELTAF2_RENORMALIZATION_SCHEME_ID = (
    "paper_0710_1869.deltaf2.tree_level_matching_scheme.v1"
)
PAPER_0710_1869_DELTAF2_OPERATOR_NORMALIZATION_ID = (
    "paper_0710_1869.deltaf2.kk_gluon_tree_color_normalization.v1"
)
PAPER_0710_1869_DELTAF2_MATCHING_ID = "eft_matching.lo.v1"

PAPER_0710_1869_DELTAF2_Q1_VLL = "Q1_VLL"
PAPER_0710_1869_DELTAF2_Q1_VRR = "Q1_VRR"
PAPER_0710_1869_DELTAF2_Q4_LR = "Q4_LR"
PAPER_0710_1869_DELTAF2_Q5_LR = "Q5_LR"
PAPER_0710_1869_DELTAF2_OPERATOR_ORDER = (
    PAPER_0710_1869_DELTAF2_Q1_VLL,
    PAPER_0710_1869_DELTAF2_Q1_VRR,
    PAPER_0710_1869_DELTAF2_Q4_LR,
    PAPER_0710_1869_DELTAF2_Q5_LR,
)


def _require_finite_complex(name: str, value: complex) -> complex:
    complex_value = complex(value)
    if not math.isfinite(complex_value.real) or not math.isfinite(complex_value.imag):
        raise ValueError(f"{name} must be finite")
    real_part = 0.0 if complex_value.real == 0.0 else float(complex_value.real)
    imag_part = 0.0 if complex_value.imag == 0.0 else float(complex_value.imag)
    return complex(real_part, imag_part)


def _complex_as_dict(value: complex) -> dict[str, float]:
    complex_value = _require_finite_complex("value", value)
    return {
        "real": float(complex_value.real),
        "imag": float(complex_value.imag),
    }


def _require_generations(
    generations: tuple[int, int] | list[int],
) -> tuple[int, int]:
    values = tuple(generations)
    if len(values) != 2:
        raise ValueError("generations must contain exactly two flavor indices")
    if any(not isinstance(item, int) or isinstance(item, bool) for item in values):
        raise ValueError("generations must contain integer flavor indices")
    i, j = values
    if i == j or i not in (0, 1, 2) or j not in (0, 1, 2):
        raise ValueError("generations must be distinct indices in {0, 1, 2}")
    return (i, j)


@dataclass(frozen=True)
class Paper07101869DeltaF2Operator:
    """One frozen operator entry in the minimal paper-owned basis."""

    name: str
    definition_id: str
    chirality: str
    lorentz_structure: str
    color_structure: str
    operator_formula: str
    schema_id: str = PAPER_0710_1869_DELTAF2_OPERATOR_SCHEMA_ID

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "schema_id",
            require_known_schema_id(
                "schema_id",
                self.schema_id,
                expected=PAPER_0710_1869_DELTAF2_OPERATOR_SCHEMA_ID,
            ),
        )
        require_member("name", self.name, PAPER_0710_1869_DELTAF2_OPERATOR_ORDER)
        require_member("chirality", self.chirality, ("VLL", "VRR", "LR"))
        for field_name in (
            "definition_id",
            "lorentz_structure",
            "color_structure",
            "operator_formula",
        ):
            object.__setattr__(
                self,
                field_name,
                require_nonempty_identifier(field_name, getattr(self, field_name)),
            )

    def as_dict(self) -> dict[str, str]:
        return {
            "schema_id": self.schema_id,
            "name": self.name,
            "definition_id": self.definition_id,
            "chirality": self.chirality,
            "lorentz_structure": self.lorentz_structure,
            "color_structure": self.color_structure,
            "operator_formula": self.operator_formula,
        }


@dataclass(frozen=True)
class Paper07101869DeltaF2OperatorBasis:
    """Explicit minimal four-operator basis used by the paper-owned matching slice."""

    operators: tuple[Paper07101869DeltaF2Operator, ...]
    basis_id: str = PAPER_0710_1869_DELTAF2_OPERATOR_BASIS_ID
    renormalization_scheme_id: str = PAPER_0710_1869_DELTAF2_RENORMALIZATION_SCHEME_ID
    operator_normalization_id: str = PAPER_0710_1869_DELTAF2_OPERATOR_NORMALIZATION_ID
    projector_normalization_id: str = PAPER_0710_1869_DELTAF2_PROJECTOR_NORMALIZATION_ID
    projector_normalization_note: str = PAPER_0710_1869_DELTAF2_PROJECTOR_NORMALIZATION_NOTE
    matching_id: str = PAPER_0710_1869_DELTAF2_MATCHING_ID
    schema_id: str = PAPER_0710_1869_DELTAF2_OPERATOR_BASIS_SCHEMA_ID
    mode_id: str = PAPER_0710_1869_MODE_ID
    paper_id: str = PAPER_0710_1869_PAPER_ID
    metadata_compatibility_note: str = (
        PAPER_0710_1869_DELTAF2_OPERATOR_METADATA_COMPATIBILITY_NOTE
    )
    notes: str = (
        "Freeze the tree-level basis to Q1_VLL, Q1_VRR, and explicit O4/O5 scalar "
        "LR operators Q4_LR and Q5_LR in the basis used throughout the DeltaF=2 "
        "literature."
    )

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "schema_id",
            require_known_schema_id(
                "schema_id",
                self.schema_id,
                expected=PAPER_0710_1869_DELTAF2_OPERATOR_BASIS_SCHEMA_ID,
            ),
        )
        require_member("mode_id", self.mode_id, (PAPER_0710_1869_MODE_ID,))
        require_member("paper_id", self.paper_id, (PAPER_0710_1869_PAPER_ID,))
        for field_name in (
            "basis_id",
            "renormalization_scheme_id",
            "operator_normalization_id",
            "projector_normalization_id",
            "projector_normalization_note",
            "matching_id",
            "metadata_compatibility_note",
            "notes",
        ):
            object.__setattr__(
                self,
                field_name,
                require_nonempty_identifier(field_name, getattr(self, field_name)),
            )
        if len(self.operators) != len(PAPER_0710_1869_DELTAF2_OPERATOR_ORDER):
            raise ValueError(
                "operators must contain exactly the frozen four-operator Delta F = 2 basis"
            )
        if not all(isinstance(item, Paper07101869DeltaF2Operator) for item in self.operators):
            raise ValueError(
                "operators must contain only Paper07101869DeltaF2Operator entries"
            )
        operator_names = tuple(item.name for item in self.operators)
        if operator_names != PAPER_0710_1869_DELTAF2_OPERATOR_ORDER:
            raise ValueError(
                "operators must follow the frozen order "
                f"{PAPER_0710_1869_DELTAF2_OPERATOR_ORDER!r}"
            )

    @property
    def operator_names(self) -> tuple[str, ...]:
        return tuple(item.name for item in self.operators)

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_id": self.schema_id,
            "mode_id": self.mode_id,
            "paper_id": self.paper_id,
            "basis_id": self.basis_id,
            "renormalization_scheme_id": self.renormalization_scheme_id,
            "operator_normalization_id": self.operator_normalization_id,
            "projector_normalization_id": self.projector_normalization_id,
            "projector_normalization_note": self.projector_normalization_note,
            "matching_id": self.matching_id,
            "operators": [item.as_dict() for item in self.operators],
            "metadata_compatibility_note": self.metadata_compatibility_note,
            "notes": self.notes,
        }


@dataclass(frozen=True)
class Paper07101869DeltaF2WilsonContract:
    """Self-describing Wilson contract with stable tag identifiers."""

    operator_basis: Paper07101869DeltaF2OperatorBasis
    kk_gluon_normalization_id: str = PAPER_0710_1869_GS_NORMALIZATION_ID
    dimensionless_matrix_policy_id: str = PAPER_0710_1869_DIMENSIONLESS_MATRIX_POLICY_ID
    universal_subtraction_policy_id: str = PAPER_0710_1869_UNIVERSAL_SUBTRACTION_POLICY_ID
    propagator_mass_rule_id: str = PAPER_0710_1869_PROPAGATOR_MASS_RULE_ID
    matching_scale_name: str = "mu_match_GeV"
    propagator_mass_name: str = "propagator_mass_GeV"
    schema_id: str = PAPER_0710_1869_DELTAF2_WILSON_CONTRACT_SCHEMA_ID
    mode_id: str = PAPER_0710_1869_MODE_ID
    paper_id: str = PAPER_0710_1869_PAPER_ID
    notes: str = (
        "Wilson objects carry basis, scheme, normalization, matching, and scale tags "
        "explicitly so downstream artifact export can remain deterministic while the "
        "frozen LR operator definitions stay audit-ready."
    )

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "schema_id",
            require_known_schema_id(
                "schema_id",
                self.schema_id,
                expected=PAPER_0710_1869_DELTAF2_WILSON_CONTRACT_SCHEMA_ID,
            ),
        )
        require_member("mode_id", self.mode_id, (PAPER_0710_1869_MODE_ID,))
        require_member("paper_id", self.paper_id, (PAPER_0710_1869_PAPER_ID,))
        if not isinstance(self.operator_basis, Paper07101869DeltaF2OperatorBasis):
            raise ValueError("operator_basis must be a Paper07101869DeltaF2OperatorBasis")
        for field_name in (
            "kk_gluon_normalization_id",
            "dimensionless_matrix_policy_id",
            "universal_subtraction_policy_id",
            "propagator_mass_rule_id",
            "matching_scale_name",
            "propagator_mass_name",
            "notes",
        ):
            object.__setattr__(
                self,
                field_name,
                require_nonempty_identifier(field_name, getattr(self, field_name)),
            )

    @property
    def operator_basis_id(self) -> str:
        return self.operator_basis.basis_id

    @property
    def renormalization_scheme_id(self) -> str:
        return self.operator_basis.renormalization_scheme_id

    @property
    def operator_normalization_id(self) -> str:
        return self.operator_basis.operator_normalization_id

    @property
    def projector_normalization_id(self) -> str:
        return self.operator_basis.projector_normalization_id

    @property
    def matching_id(self) -> str:
        return self.operator_basis.matching_id

    @property
    def operator_order(self) -> tuple[str, ...]:
        return self.operator_basis.operator_names

    @property
    def tags(self) -> dict[str, str]:
        return {
            "mode_id": self.mode_id,
            "paper_id": self.paper_id,
            "operator_basis_id": self.operator_basis_id,
            "renormalization_scheme_id": self.renormalization_scheme_id,
            "operator_normalization_id": self.operator_normalization_id,
            "projector_normalization_id": self.projector_normalization_id,
            "matching_id": self.matching_id,
            "kk_gluon_normalization_id": self.kk_gluon_normalization_id,
            "dimensionless_matrix_policy_id": self.dimensionless_matrix_policy_id,
            "universal_subtraction_policy_id": self.universal_subtraction_policy_id,
            "propagator_mass_rule_id": self.propagator_mass_rule_id,
            "matching_scale_name": self.matching_scale_name,
            "propagator_mass_name": self.propagator_mass_name,
        }

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_id": self.schema_id,
            "mode_id": self.mode_id,
            "paper_id": self.paper_id,
            "operator_basis": self.operator_basis.as_dict(),
            "kk_gluon_normalization_id": self.kk_gluon_normalization_id,
            "dimensionless_matrix_policy_id": self.dimensionless_matrix_policy_id,
            "universal_subtraction_policy_id": self.universal_subtraction_policy_id,
            "propagator_mass_rule_id": self.propagator_mass_rule_id,
            "matching_scale_name": self.matching_scale_name,
            "propagator_mass_name": self.propagator_mass_name,
            "projector_normalization_id": self.projector_normalization_id,
            "tags": self.tags,
            "notes": self.notes,
        }


@dataclass(frozen=True)
class Paper07101869DeltaF2WilsonCoefficients:
    """Wilson coefficients for one paper-owned ``Delta F = 2`` benchmark system."""

    contract: Paper07101869DeltaF2WilsonContract
    benchmark_id: str
    scale_label: str
    system_id: str
    sector_id: str
    generations: tuple[int, int]
    matching_scale_GeV: float
    propagator_mass_GeV: float
    left_coupling: complex
    right_coupling: complex
    q1_vll: complex
    q1_vrr: complex
    q4_lr: complex
    q5_lr: complex
    schema_id: str = PAPER_0710_1869_DELTAF2_WILSON_SCHEMA_ID

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "schema_id",
            require_known_schema_id(
                "schema_id",
                self.schema_id,
                expected=PAPER_0710_1869_DELTAF2_WILSON_SCHEMA_ID,
            ),
        )
        if not isinstance(self.contract, Paper07101869DeltaF2WilsonContract):
            raise ValueError("contract must be a Paper07101869DeltaF2WilsonContract")
        for field_name in ("benchmark_id", "scale_label", "system_id"):
            object.__setattr__(
                self,
                field_name,
                require_nonempty_identifier(field_name, getattr(self, field_name)),
            )
        require_member("sector_id", self.sector_id, ("down", "up"))
        object.__setattr__(self, "generations", _require_generations(self.generations))
        object.__setattr__(
            self,
            "matching_scale_GeV",
            require_positive_finite("matching_scale_GeV", self.matching_scale_GeV),
        )
        object.__setattr__(
            self,
            "propagator_mass_GeV",
            require_positive_finite("propagator_mass_GeV", self.propagator_mass_GeV),
        )
        for field_name in (
            "left_coupling",
            "right_coupling",
            "q1_vll",
            "q1_vrr",
            "q4_lr",
            "q5_lr",
        ):
            object.__setattr__(
                self,
                field_name,
                _require_finite_complex(field_name, getattr(self, field_name)),
            )

    @property
    def coefficients(self) -> dict[str, complex]:
        return {
            PAPER_0710_1869_DELTAF2_Q1_VLL: self.q1_vll,
            PAPER_0710_1869_DELTAF2_Q1_VRR: self.q1_vrr,
            PAPER_0710_1869_DELTAF2_Q4_LR: self.q4_lr,
            PAPER_0710_1869_DELTAF2_Q5_LR: self.q5_lr,
        }

    @property
    def tags(self) -> dict[str, object]:
        return {
            **self.contract.tags,
            "benchmark_id": self.benchmark_id,
            "scale_label": self.scale_label,
            "system_id": self.system_id,
            "sector_id": self.sector_id,
            "generations": list(self.generations),
            "matching_scale_GeV": self.matching_scale_GeV,
            "propagator_mass_GeV": self.propagator_mass_GeV,
        }

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_id": self.schema_id,
            "benchmark_id": self.benchmark_id,
            "scale_label": self.scale_label,
            "system_id": self.system_id,
            "sector_id": self.sector_id,
            "generations": list(self.generations),
            "matching_scale_GeV": self.matching_scale_GeV,
            "propagator_mass_GeV": self.propagator_mass_GeV,
            "left_coupling": _complex_as_dict(self.left_coupling),
            "right_coupling": _complex_as_dict(self.right_coupling),
            "coefficients": {
                name: _complex_as_dict(value) for name, value in self.coefficients.items()
            },
            "tags": self.tags,
        }


def default_paper_0710_1869_deltaf2_operator_basis() -> Paper07101869DeltaF2OperatorBasis:
    """Return the frozen minimal four-operator basis."""
    return Paper07101869DeltaF2OperatorBasis(
        operators=(
            Paper07101869DeltaF2Operator(
                name=PAPER_0710_1869_DELTAF2_Q1_VLL,
                definition_id=PAPER_0710_1869_DELTAF2_OPERATOR_DEFINITION_Q1_VLL_ID,
                chirality="VLL",
                lorentz_structure="(V-A)x(V-A)",
                color_structure="color_singlet_current_current",
                operator_formula=(
                    "(bar q_hi^alpha gamma_mu P_L q_lo^alpha)"
                    "(bar q_hi^beta gamma^mu P_L q_lo^beta)"
                ),
            ),
            Paper07101869DeltaF2Operator(
                name=PAPER_0710_1869_DELTAF2_Q1_VRR,
                definition_id=PAPER_0710_1869_DELTAF2_OPERATOR_DEFINITION_Q1_VRR_ID,
                chirality="VRR",
                lorentz_structure="(V+A)x(V+A)",
                color_structure="color_singlet_current_current",
                operator_formula=(
                    "(bar q_hi^alpha gamma_mu P_R q_lo^alpha)"
                    "(bar q_hi^beta gamma^mu P_R q_lo^beta)"
                ),
            ),
            Paper07101869DeltaF2Operator(
                name=PAPER_0710_1869_DELTAF2_Q4_LR,
                definition_id=PAPER_0710_1869_DELTAF2_OPERATOR_DEFINITION_Q4_LR_ID,
                chirality="LR",
                lorentz_structure="(S-P)x(S+P)",
                color_structure="color_singlet_density_density",
                operator_formula=(
                    "(bar q_hi^alpha P_L q_lo^alpha)"
                    "(bar q_hi^beta P_R q_lo^beta)"
                ),
            ),
            Paper07101869DeltaF2Operator(
                name=PAPER_0710_1869_DELTAF2_Q5_LR,
                definition_id=PAPER_0710_1869_DELTAF2_OPERATOR_DEFINITION_Q5_LR_ID,
                chirality="LR",
                lorentz_structure="(S-P)x(S+P)",
                color_structure="color_mixed_density_density",
                operator_formula=(
                    "(bar q_hi^alpha P_L q_lo^beta)"
                    "(bar q_hi^beta P_R q_lo^alpha)"
                ),
            ),
        )
    )


def default_paper_0710_1869_deltaf2_wilson_contract() -> Paper07101869DeltaF2WilsonContract:
    """Return the default self-describing Wilson contract."""
    return Paper07101869DeltaF2WilsonContract(
        operator_basis=default_paper_0710_1869_deltaf2_operator_basis()
    )


__all__ = [
    "PAPER_0710_1869_DELTAF2_MATCHING_ID",
    "PAPER_0710_1869_DELTAF2_OPERATOR_BASIS_ID",
    "PAPER_0710_1869_DELTAF2_OPERATOR_BASIS_SCHEMA_ID",
    "PAPER_0710_1869_DELTAF2_OPERATOR_DEFINITION_Q1_VLL_ID",
    "PAPER_0710_1869_DELTAF2_OPERATOR_DEFINITION_Q1_VRR_ID",
    "PAPER_0710_1869_DELTAF2_OPERATOR_DEFINITION_Q4_LR_ID",
    "PAPER_0710_1869_DELTAF2_OPERATOR_DEFINITION_Q5_LR_ID",
    "PAPER_0710_1869_DELTAF2_OPERATOR_METADATA_COMPATIBILITY_NOTE",
    "PAPER_0710_1869_DELTAF2_OPERATOR_NORMALIZATION_ID",
    "PAPER_0710_1869_DELTAF2_OPERATOR_ORDER",
    "PAPER_0710_1869_DELTAF2_PROJECTOR_NORMALIZATION_ID",
    "PAPER_0710_1869_DELTAF2_PROJECTOR_NORMALIZATION_NOTE",
    "PAPER_0710_1869_DELTAF2_OPERATOR_SCHEMA_ID",
    "PAPER_0710_1869_DELTAF2_Q1_VLL",
    "PAPER_0710_1869_DELTAF2_Q1_VRR",
    "PAPER_0710_1869_DELTAF2_Q4_LR",
    "PAPER_0710_1869_DELTAF2_Q5_LR",
    "PAPER_0710_1869_DELTAF2_RENORMALIZATION_SCHEME_ID",
    "PAPER_0710_1869_DELTAF2_WILSON_CONTRACT_SCHEMA_ID",
    "PAPER_0710_1869_DELTAF2_WILSON_SCHEMA_ID",
    "Paper07101869DeltaF2Operator",
    "Paper07101869DeltaF2OperatorBasis",
    "Paper07101869DeltaF2WilsonCoefficients",
    "Paper07101869DeltaF2WilsonContract",
    "default_paper_0710_1869_deltaf2_operator_basis",
    "default_paper_0710_1869_deltaf2_wilson_contract",
]
