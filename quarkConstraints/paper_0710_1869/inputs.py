"""Sourced benchmark inputs owned by the dedicated ``paper_0710_1869`` mode."""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from .conventions import (
    PAPER_0710_1869_MODE_ID,
    PAPER_0710_1869_PAPER_ID,
    PAPER_0710_1869_PROFILE_DERIVATION_POLICY_ID,
    PAPER_0710_1869_SEED_TO_PROFILE_MAPPING_POLICY_ID,
    PAPER_0710_1869_UNIVERSAL_TERM_COEFFICIENT_POLICY_ID,
)
from .validation import (
    require_known_schema_id,
    require_member,
    require_nonempty_identifier,
    require_positive_finite,
)

PAPER_0710_1869_SOURCE_SCHEMA_ID = "quarkConstraints.paper_0710_1869.source_ref.v1"
PAPER_0710_1869_TABLE_I_SCHEMA_ID = "quarkConstraints.paper_0710_1869.table_i_inputs.v1"
PAPER_0710_1869_EQ3_SCHEMA_ID = "quarkConstraints.paper_0710_1869.eq3_example.v1"
PAPER_0710_1869_AFFINE_BULK_MASS_SECTOR_POLICY_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.affine_bulk_mass_sector_policy.v1"
)
PAPER_0710_1869_UNIVERSAL_TERM_POLICY_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.universal_term_policy.v1"
)
PAPER_0710_1869_SEED_TO_PROFILE_MAPPING_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.seed_to_profile_mapping_policy.v1"
)
PAPER_0710_1869_PHYSICAL_SEED_TO_PROFILE_CONTRACT_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.seed_to_profile_contract.v1"
)
PAPER_0710_1869_AFFINE_BULK_MASS_LEADING_TERM_COEFFICIENT = -1.0
PAPER_0710_1869_AFFINE_BULK_MASS_UNIVERSAL_OFFSET = 0.6
PAPER_0710_1869_BULK_MASS_PARAMETERIZATION_ID = (
    "quarkConstraints.paper_0710_1869.bulk_mass_parameterization."
    "affine_per_sector_eigenvalues.provisional_negative_slope.v2"
)
PAPER_0710_1869_EIGENVALUE_ORDERING_ID = (
    "quarkConstraints.paper_0710_1869.bulk_mass_ordering.ascending_eigenvalues.v1"
)


def _require_finite_float(name: str, value: float) -> float:
    numeric = float(value)
    if not math.isfinite(numeric):
        raise ValueError(f"{name} must be a finite float")
    return numeric


def _require_exact_finite_float(name: str, value: float, *, expected: float) -> float:
    numeric = _require_finite_float(name, value)
    if numeric != expected:
        raise ValueError(f"{name} must be exactly {expected!r}")
    return numeric


def _require_negative_finite_float(name: str, value: float) -> float:
    numeric = _require_finite_float(name, value)
    if numeric >= 0.0:
        raise ValueError(f"{name} must be negative")
    return numeric


def _require_bool(name: str, value: bool) -> bool:
    if not isinstance(value, bool):
        raise ValueError(f"{name} must be a bool")
    return value


@dataclass(frozen=True)
class Paper07101869SourceRef:
    """Frozen citation metadata for one sourced paper value block."""

    schema_id: str = PAPER_0710_1869_SOURCE_SCHEMA_ID
    paper_id: str = PAPER_0710_1869_PAPER_ID
    source_kind: str = "table"
    locator_label: str = "Table I"
    detail: str = "c_Q/f_Q eigenvalues"
    citation: str = "arXiv:0710.1869, Table I"
    notes: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "schema_id",
            require_known_schema_id(
                "schema_id",
                self.schema_id,
                expected=PAPER_0710_1869_SOURCE_SCHEMA_ID,
            ),
        )
        object.__setattr__(self, "paper_id", require_nonempty_identifier("paper_id", self.paper_id))
        require_member("paper_id", self.paper_id, (PAPER_0710_1869_PAPER_ID,))
        object.__setattr__(
            self,
            "source_kind",
            require_member("source_kind", self.source_kind, ("table", "equation")),
        )
        object.__setattr__(
            self,
            "locator_label",
            require_nonempty_identifier("locator_label", self.locator_label),
        )
        object.__setattr__(self, "detail", require_nonempty_identifier("detail", self.detail))
        object.__setattr__(self, "citation", require_nonempty_identifier("citation", self.citation))
        if self.notes is not None:
            object.__setattr__(self, "notes", require_nonempty_identifier("notes", self.notes))

    def as_dict(self) -> dict[str, str | None]:
        return {
            "schema_id": self.schema_id,
            "paper_id": self.paper_id,
            "source_kind": self.source_kind,
            "locator_label": self.locator_label,
            "detail": self.detail,
            "citation": self.citation,
            "notes": self.notes,
        }


@dataclass(frozen=True)
class Paper07101869EigenvalueEntry:
    """One generation entry for the sourced Table I profile eigenvalues."""

    generation: int
    c_value: float
    f_value: float

    def __post_init__(self) -> None:
        if not isinstance(self.generation, int) or isinstance(self.generation, bool):
            raise ValueError("generation must be an integer")
        require_member("generation", str(self.generation), ("1", "2", "3"))
        object.__setattr__(self, "c_value", _require_finite_float("c_value", self.c_value))
        object.__setattr__(self, "f_value", require_positive_finite("f_value", self.f_value))

    def as_dict(self) -> dict[str, int | float]:
        return {
            "generation": self.generation,
            "c_value": self.c_value,
            "f_value": self.f_value,
        }


@dataclass(frozen=True)
class Paper07101869SectorTableInputs:
    """One sourced Table I sector block."""

    sector_id: str
    source: Paper07101869SourceRef
    entries: tuple[Paper07101869EigenvalueEntry, ...]

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "sector_id",
            require_member("sector_id", self.sector_id, ("Q", "u", "d")),
        )
        if not isinstance(self.source, Paper07101869SourceRef):
            raise ValueError("source must be a Paper07101869SourceRef")
        if len(self.entries) != 3:
            raise ValueError("entries must contain exactly three generations")
        if not all(isinstance(entry, Paper07101869EigenvalueEntry) for entry in self.entries):
            raise ValueError("entries must contain only Paper07101869EigenvalueEntry items")
        generations = tuple(entry.generation for entry in self.entries)
        if generations != (1, 2, 3):
            raise ValueError("entries must be ordered by generation 1, 2, 3")

    def as_dict(self) -> dict[str, object]:
        return {
            "sector_id": self.sector_id,
            "source": self.source.as_dict(),
            "entries": [entry.as_dict() for entry in self.entries],
        }

    @classmethod
    def from_pairs(
        cls,
        *,
        sector_id: str,
        source: Paper07101869SourceRef,
        pairs: tuple[tuple[float, float], tuple[float, float], tuple[float, float]],
    ) -> "Paper07101869SectorTableInputs":
        entries = tuple(
            Paper07101869EigenvalueEntry(
                generation=index,
                c_value=values[0],
                f_value=values[1],
            )
            for index, values in enumerate(pairs, start=1)
        )
        return cls(sector_id=sector_id, source=source, entries=entries)


@dataclass(frozen=True)
class Paper07101869TableIInputs:
    """Closed, sourced copy of the Table I eigenvalue inputs."""

    schema_id: str = PAPER_0710_1869_TABLE_I_SCHEMA_ID
    mode_id: str = PAPER_0710_1869_MODE_ID
    sectors: tuple[Paper07101869SectorTableInputs, ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "schema_id",
            require_known_schema_id(
                "schema_id",
                self.schema_id,
                expected=PAPER_0710_1869_TABLE_I_SCHEMA_ID,
            ),
        )
        object.__setattr__(self, "mode_id", require_nonempty_identifier("mode_id", self.mode_id))
        require_member("mode_id", self.mode_id, (PAPER_0710_1869_MODE_ID,))
        if len(self.sectors) != 3:
            raise ValueError("sectors must contain exactly Q, u, and d")
        if not all(isinstance(sector, Paper07101869SectorTableInputs) for sector in self.sectors):
            raise ValueError("sectors must contain only Paper07101869SectorTableInputs")
        sector_ids = tuple(sector.sector_id for sector in self.sectors)
        if sector_ids != ("Q", "u", "d"):
            raise ValueError("sectors must be ordered as Q, u, d")

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_id": self.schema_id,
            "mode_id": self.mode_id,
            "sectors": [sector.as_dict() for sector in self.sectors],
        }


@dataclass(frozen=True)
class Paper07101869Eq3Example:
    """Sourced structural copy of the Eq. (3) example parameters."""

    schema_id: str = PAPER_0710_1869_EQ3_SCHEMA_ID
    mode_id: str = PAPER_0710_1869_MODE_ID
    source: Paper07101869SourceRef = field(
        default_factory=lambda: Paper07101869SourceRef(
            source_kind="equation",
            locator_label="Eq. (3)",
            detail="example parameters (a, r, theta12, theta23, theta13, delta)",
            citation="arXiv:0710.1869, Eq. (3)",
            notes=(
                "Angles are carried in degrees exactly as extracted in the local paper notes. "
                "delta is carried as the quoted paper parameter without any extra conversion."
            ),
        )
    )
    a: float = 0.8
    r: float = 0.3
    theta12_deg: float = 115.0
    theta23_deg: float = 65.0
    theta13_deg: float = 70.0
    delta: float = 0.6

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "schema_id",
            require_known_schema_id(
                "schema_id",
                self.schema_id,
                expected=PAPER_0710_1869_EQ3_SCHEMA_ID,
            ),
        )
        object.__setattr__(self, "mode_id", require_nonempty_identifier("mode_id", self.mode_id))
        require_member("mode_id", self.mode_id, (PAPER_0710_1869_MODE_ID,))
        if not isinstance(self.source, Paper07101869SourceRef):
            raise ValueError("source must be a Paper07101869SourceRef")
        object.__setattr__(self, "a", require_positive_finite("a", self.a))
        object.__setattr__(self, "r", require_positive_finite("r", self.r))
        object.__setattr__(
            self,
            "theta12_deg",
            _require_finite_float("theta12_deg", self.theta12_deg),
        )
        object.__setattr__(
            self,
            "theta23_deg",
            _require_finite_float("theta23_deg", self.theta23_deg),
        )
        object.__setattr__(
            self,
            "theta13_deg",
            _require_finite_float("theta13_deg", self.theta13_deg),
        )
        object.__setattr__(self, "delta", _require_finite_float("delta", self.delta))

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_id": self.schema_id,
            "mode_id": self.mode_id,
            "source": self.source.as_dict(),
            "a": self.a,
            "r": self.r,
            "theta12_deg": self.theta12_deg,
            "theta23_deg": self.theta23_deg,
            "theta13_deg": self.theta13_deg,
            "delta": self.delta,
        }


@dataclass(frozen=True)
class Paper07101869AffineBulkMassSectorPolicy:
    """Affine bulk-mass map for one MFV sector."""

    schema_id: str = PAPER_0710_1869_AFFINE_BULK_MASS_SECTOR_POLICY_SCHEMA_ID
    sector_id: str = "Q"
    leading_term_coefficient: float = PAPER_0710_1869_AFFINE_BULK_MASS_LEADING_TERM_COEFFICIENT
    universal_offset: float = PAPER_0710_1869_AFFINE_BULK_MASS_UNIVERSAL_OFFSET
    notes: str = (
        "RESIDUAL(C-6): exact Table-I affine coefficients pending paper 0710.1869. "
        "Provisional negative-slope affine eigenvalue map for one sector; larger "
        "Yukawa eigenvalues map to smaller c values so the heaviest quark is the "
        "most IR-localized."
    )

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "schema_id",
            require_known_schema_id(
                "schema_id",
                self.schema_id,
                expected=PAPER_0710_1869_AFFINE_BULK_MASS_SECTOR_POLICY_SCHEMA_ID,
            ),
        )
        object.__setattr__(
            self,
            "sector_id",
            require_member("sector_id", self.sector_id, ("Q", "u", "d")),
        )
        object.__setattr__(
            self,
            "leading_term_coefficient",
            _require_negative_finite_float(
                "leading_term_coefficient",
                self.leading_term_coefficient,
            ),
        )
        object.__setattr__(
            self,
            "universal_offset",
            _require_finite_float(
                "universal_offset",
                self.universal_offset,
            ),
        )
        object.__setattr__(self, "notes", require_nonempty_identifier("notes", self.notes))

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_id": self.schema_id,
            "sector_id": self.sector_id,
            "leading_term_coefficient": self.leading_term_coefficient,
            "universal_offset": self.universal_offset,
            "notes": self.notes,
        }


@dataclass(frozen=True)
class Paper07101869UniversalTermPolicy:
    """Explicit defaults for omitted universal terms and leading coefficients."""

    schema_id: str = PAPER_0710_1869_UNIVERSAL_TERM_POLICY_SCHEMA_ID
    mode_id: str = PAPER_0710_1869_MODE_ID
    policy_id: str = PAPER_0710_1869_UNIVERSAL_TERM_COEFFICIENT_POLICY_ID
    source: Paper07101869SourceRef = field(
        default_factory=lambda: Paper07101869SourceRef(
            source_kind="equation",
            locator_label="Eq. (1)",
            detail="MFV bulk-mass construction with omitted universal terms and coefficients",
            citation="arXiv:0710.1869, Eq. (1)",
            notes=(
                "The paper text suppresses or omits additive universal terms and "
                "order-one coefficients. This package freezes them as explicit "
                "policy data instead of inheriting a hidden BulkMassMap surrogate."
            ),
        )
    )
    sector_policies: tuple[
        Paper07101869AffineBulkMassSectorPolicy,
        Paper07101869AffineBulkMassSectorPolicy,
        Paper07101869AffineBulkMassSectorPolicy,
    ] = field(
        default_factory=lambda: (
            Paper07101869AffineBulkMassSectorPolicy(sector_id="Q"),
            Paper07101869AffineBulkMassSectorPolicy(sector_id="u"),
            Paper07101869AffineBulkMassSectorPolicy(sector_id="d"),
        )
    )
    omitted_terms_treatment_id: str = (
        "quarkConstraints.paper_0710_1869.omitted_terms.provisional_affine_offsets.v2"
    )
    notes: str = (
        "RESIDUAL(C-6): exact Table-I affine coefficients pending paper 0710.1869. "
        "Closed default policy for the missing paper-honest physical bridge: "
        "sector policies are explicit affine maps with negative leading "
        "coefficients until a later reviewed version installs the sourced values."
    )

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "schema_id",
            require_known_schema_id(
                "schema_id",
                self.schema_id,
                expected=PAPER_0710_1869_UNIVERSAL_TERM_POLICY_SCHEMA_ID,
            ),
        )
        object.__setattr__(self, "mode_id", require_nonempty_identifier("mode_id", self.mode_id))
        require_member("mode_id", self.mode_id, (PAPER_0710_1869_MODE_ID,))
        object.__setattr__(
            self,
            "policy_id",
            require_known_schema_id(
                "policy_id",
                self.policy_id,
                expected=PAPER_0710_1869_UNIVERSAL_TERM_COEFFICIENT_POLICY_ID,
            ),
        )
        if not isinstance(self.source, Paper07101869SourceRef):
            raise ValueError("source must be a Paper07101869SourceRef")
        policies = tuple(self.sector_policies)
        if len(policies) != 3:
            raise ValueError("sector_policies must contain exactly Q, u, and d")
        if not all(isinstance(policy, Paper07101869AffineBulkMassSectorPolicy) for policy in policies):
            raise ValueError(
                "sector_policies must contain only Paper07101869AffineBulkMassSectorPolicy items"
            )
        sector_ids = tuple(policy.sector_id for policy in policies)
        if sector_ids != ("Q", "u", "d"):
            raise ValueError("sector_policies must be ordered as Q, u, d")
        for policy in policies:
            if policy.leading_term_coefficient >= 0.0:
                raise ValueError("sector_policies must use negative leading_term_coefficient")
            _require_finite_float("sector_policy.universal_offset", policy.universal_offset)
        object.__setattr__(self, "sector_policies", policies)
        object.__setattr__(
            self,
            "omitted_terms_treatment_id",
            require_nonempty_identifier(
                "omitted_terms_treatment_id", self.omitted_terms_treatment_id
            ),
        )
        object.__setattr__(self, "notes", require_nonempty_identifier("notes", self.notes))

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_id": self.schema_id,
            "mode_id": self.mode_id,
            "policy_id": self.policy_id,
            "source": self.source.as_dict(),
            "sector_policies": [policy.as_dict() for policy in self.sector_policies],
            "omitted_terms_treatment_id": self.omitted_terms_treatment_id,
            "notes": self.notes,
        }


@dataclass(frozen=True)
class Paper07101869SeedToProfileMappingPolicy:
    """Closed schema for the future physical seed-to-profile mapping path."""

    schema_id: str = PAPER_0710_1869_SEED_TO_PROFILE_MAPPING_SCHEMA_ID
    mode_id: str = PAPER_0710_1869_MODE_ID
    policy_id: str = PAPER_0710_1869_SEED_TO_PROFILE_MAPPING_POLICY_ID
    source: Paper07101869SourceRef = field(
        default_factory=lambda: Paper07101869SourceRef(
            source_kind="equation",
            locator_label="Eq. (1)",
            detail="seeded MFV bulk-mass matrices to affine bulk masses and geometry-derived profiles",
            citation="arXiv:0710.1869, Eq. (1)",
            notes=(
                "This is a schema freeze only. The later model slice must consume "
                "this contract directly rather than recover physical profiles from "
                "Table I aliases or any hidden BulkMassMap path."
            ),
        )
    )
    bulk_mass_parameterization_id: str = (
        PAPER_0710_1869_BULK_MASS_PARAMETERIZATION_ID
    )
    eigenvalue_ordering_id: str = (
        PAPER_0710_1869_EIGENVALUE_ORDERING_ID
    )
    profile_derivation_policy_id: str = PAPER_0710_1869_PROFILE_DERIVATION_POLICY_ID
    mapping_formula: str = "c_i^(x) = alpha_x * lambda_i(C_x) + beta_x"
    profile_formula_note: str = (
        "After the affine bulk masses are fixed, the physical profile factors F_i^(x) "
        "must be derived from the geometry using those c_i^(x) values. This schema "
        "does not permit quoted Table I profile attachment on the physical path."
    )
    uses_hidden_bulk_mass_map_surrogate: bool = False
    notes: str = (
        "Closed mapping contract for the missing QS1 physical bridge only. "
        "Masses and CKM must still be obtained later by full matrix "
        "diagonalization, not by overlap-ratio surrogates."
    )

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "schema_id",
            require_known_schema_id(
                "schema_id",
                self.schema_id,
                expected=PAPER_0710_1869_SEED_TO_PROFILE_MAPPING_SCHEMA_ID,
            ),
        )
        object.__setattr__(self, "mode_id", require_nonempty_identifier("mode_id", self.mode_id))
        require_member("mode_id", self.mode_id, (PAPER_0710_1869_MODE_ID,))
        object.__setattr__(
            self,
            "policy_id",
            require_known_schema_id(
                "policy_id",
                self.policy_id,
                expected=PAPER_0710_1869_SEED_TO_PROFILE_MAPPING_POLICY_ID,
            ),
        )
        if not isinstance(self.source, Paper07101869SourceRef):
            raise ValueError("source must be a Paper07101869SourceRef")
        object.__setattr__(
            self,
            "bulk_mass_parameterization_id",
            require_known_schema_id(
                "bulk_mass_parameterization_id",
                self.bulk_mass_parameterization_id,
                expected=PAPER_0710_1869_BULK_MASS_PARAMETERIZATION_ID,
            ),
        )
        object.__setattr__(
            self,
            "eigenvalue_ordering_id",
            require_known_schema_id(
                "eigenvalue_ordering_id",
                self.eigenvalue_ordering_id,
                expected=PAPER_0710_1869_EIGENVALUE_ORDERING_ID,
            ),
        )
        object.__setattr__(
            self,
            "profile_derivation_policy_id",
            require_known_schema_id(
                "profile_derivation_policy_id",
                self.profile_derivation_policy_id,
                expected=PAPER_0710_1869_PROFILE_DERIVATION_POLICY_ID,
            ),
        )
        object.__setattr__(
            self,
            "mapping_formula",
            require_nonempty_identifier("mapping_formula", self.mapping_formula),
        )
        object.__setattr__(
            self,
            "profile_formula_note",
            require_nonempty_identifier("profile_formula_note", self.profile_formula_note),
        )
        object.__setattr__(
            self,
            "uses_hidden_bulk_mass_map_surrogate",
            _require_bool(
                "uses_hidden_bulk_mass_map_surrogate",
                self.uses_hidden_bulk_mass_map_surrogate,
            ),
        )
        if self.uses_hidden_bulk_mass_map_surrogate:
            raise ValueError("uses_hidden_bulk_mass_map_surrogate must remain False")
        object.__setattr__(self, "notes", require_nonempty_identifier("notes", self.notes))

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_id": self.schema_id,
            "mode_id": self.mode_id,
            "policy_id": self.policy_id,
            "source": self.source.as_dict(),
            "bulk_mass_parameterization_id": self.bulk_mass_parameterization_id,
            "eigenvalue_ordering_id": self.eigenvalue_ordering_id,
            "profile_derivation_policy_id": self.profile_derivation_policy_id,
            "mapping_formula": self.mapping_formula,
            "profile_formula_note": self.profile_formula_note,
            "uses_hidden_bulk_mass_map_surrogate": self.uses_hidden_bulk_mass_map_surrogate,
            "notes": self.notes,
        }


@dataclass(frozen=True)
class Paper07101869PhysicalSeedToProfileContract:
    """Typed bundle for the frozen QS1 physical bridge schema only."""

    schema_id: str = PAPER_0710_1869_PHYSICAL_SEED_TO_PROFILE_CONTRACT_SCHEMA_ID
    mode_id: str = PAPER_0710_1869_MODE_ID
    mapping_policy: Paper07101869SeedToProfileMappingPolicy = field(
        default_factory=lambda: Paper07101869SeedToProfileMappingPolicy()
    )
    universal_term_policy: Paper07101869UniversalTermPolicy = field(
        default_factory=lambda: Paper07101869UniversalTermPolicy()
    )
    notes: str = (
        "Schema freeze for the future paper-honest physical bridge. "
        "This bundle is policy/config only and does not itself implement "
        "bulk states, fits, scans, couplings, or observables."
    )

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "schema_id",
            require_known_schema_id(
                "schema_id",
                self.schema_id,
                expected=PAPER_0710_1869_PHYSICAL_SEED_TO_PROFILE_CONTRACT_SCHEMA_ID,
            ),
        )
        object.__setattr__(self, "mode_id", require_nonempty_identifier("mode_id", self.mode_id))
        require_member("mode_id", self.mode_id, (PAPER_0710_1869_MODE_ID,))
        if not isinstance(self.mapping_policy, Paper07101869SeedToProfileMappingPolicy):
            raise ValueError("mapping_policy must be a Paper07101869SeedToProfileMappingPolicy")
        if not isinstance(self.universal_term_policy, Paper07101869UniversalTermPolicy):
            raise ValueError("universal_term_policy must be a Paper07101869UniversalTermPolicy")
        self.mapping_policy.__post_init__()
        self.universal_term_policy.__post_init__()
        if self.mapping_policy.mode_id != self.mode_id:
            raise ValueError("mapping_policy.mode_id must match contract mode_id")
        if self.universal_term_policy.mode_id != self.mode_id:
            raise ValueError("universal_term_policy.mode_id must match contract mode_id")
        if self.mapping_policy.policy_id != PAPER_0710_1869_SEED_TO_PROFILE_MAPPING_POLICY_ID:
            raise ValueError("mapping_policy.policy_id must remain the frozen exact ID")
        if self.universal_term_policy.policy_id != PAPER_0710_1869_UNIVERSAL_TERM_COEFFICIENT_POLICY_ID:
            raise ValueError("universal_term_policy.policy_id must remain the frozen exact ID")
        object.__setattr__(self, "notes", require_nonempty_identifier("notes", self.notes))

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_id": self.schema_id,
            "mode_id": self.mode_id,
            "mapping_policy": self.mapping_policy.as_dict(),
            "universal_term_policy": self.universal_term_policy.as_dict(),
            "notes": self.notes,
        }


def default_paper_0710_1869_table_i_inputs() -> Paper07101869TableIInputs:
    """Return the frozen Table I eigenvalue inputs for PR1."""
    q_source = Paper07101869SourceRef(
        source_kind="table",
        locator_label="Table I",
        detail="c_Q/f_Q eigenvalues",
        citation="arXiv:0710.1869, Table I",
    )
    u_source = Paper07101869SourceRef(
        source_kind="table",
        locator_label="Table I",
        detail="c_u/f_u eigenvalues",
        citation="arXiv:0710.1869, Table I",
    )
    d_source = Paper07101869SourceRef(
        source_kind="table",
        locator_label="Table I",
        detail="c_d/f_d eigenvalues",
        citation="arXiv:0710.1869, Table I",
    )
    return Paper07101869TableIInputs(
        sectors=(
            Paper07101869SectorTableInputs.from_pairs(
                sector_id="Q",
                source=q_source,
                pairs=((0.64, 0.002), (0.59, 0.01), (0.46, 0.2)),
            ),
            Paper07101869SectorTableInputs.from_pairs(
                sector_id="u",
                source=u_source,
                pairs=((0.68, 0.0007), (0.53, 0.06), (-0.06, 0.8)),
            ),
            Paper07101869SectorTableInputs.from_pairs(
                sector_id="d",
                source=d_source,
                pairs=((0.65, 0.002), (0.60, 0.008), (0.58, 0.02)),
            ),
        )
    )


def default_paper_0710_1869_eq3_example() -> Paper07101869Eq3Example:
    """Return the frozen Eq. (3) example parameters for PR1."""
    return Paper07101869Eq3Example()


def default_paper_0710_1869_universal_term_policy() -> Paper07101869UniversalTermPolicy:
    """Return the frozen default universal-term/coefficient policy for QS1."""
    return Paper07101869UniversalTermPolicy()


def default_paper_0710_1869_seed_to_profile_mapping_policy() -> (
    Paper07101869SeedToProfileMappingPolicy
):
    """Return the frozen default seed-to-profile mapping policy for QS1."""
    return Paper07101869SeedToProfileMappingPolicy()


def default_paper_0710_1869_physical_seed_to_profile_contract() -> (
    Paper07101869PhysicalSeedToProfileContract
):
    """Return the closed policy bundle for the future physical QS1 bridge."""
    return Paper07101869PhysicalSeedToProfileContract()
