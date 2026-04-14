"""Canonical conventions for the paper-facing 0710.1869 mode."""

from __future__ import annotations

from dataclasses import asdict, dataclass

from .validation import require_known_schema_id, require_member, require_nonempty_identifier

PAPER_0710_1869_MODE_ID = "paper_0710_1869"
PAPER_0710_1869_PAPER_ID = "arXiv:0710.1869"
PAPER_0710_1869_CONVENTIONS_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.conventions.v1"
)
PAPER_0710_1869_OPERATOR_BASIS_ID = "kk_gluon_tree_np_only.v1"
PAPER_0710_1869_MATCHING_ID = "eft_matching.lo.v1"
PAPER_0710_1869_KK_GLUON_NORMALIZATION_ID = "explicit_mu_gs.normalization.v1"
PAPER_0710_1869_PROVENANCE_POLICY_ID = "deterministic_bundle.required.v1"
PAPER_0710_1869_VERIFIER_POLICY_ID = "independent_verifier.required.v1"
PAPER_0710_1869_SEED_TO_PROFILE_MAPPING_POLICY_ID = (
    "quarkConstraints.paper_0710_1869.seed_to_profile_mapping."
    "affine_per_sector_eigenvalue.v1"
)
PAPER_0710_1869_UNIVERSAL_TERM_COEFFICIENT_POLICY_ID = (
    "quarkConstraints.paper_0710_1869.universal_terms."
    "explicit_zero_offset_unit_leading_coefficients.v1"
)
PAPER_0710_1869_PROFILE_DERIVATION_POLICY_ID = (
    "quarkConstraints.paper_0710_1869.profile_derivation."
    "geometry_from_affine_bulk_masses.v1"
)


@dataclass(frozen=True)
class Paper07101869Conventions:
    """Closed convention bundle for the canonical paper-mode path."""

    schema_id: str = PAPER_0710_1869_CONVENTIONS_SCHEMA_ID
    mode_id: str = PAPER_0710_1869_MODE_ID
    paper_id: str = PAPER_0710_1869_PAPER_ID
    operator_basis_id: str = PAPER_0710_1869_OPERATOR_BASIS_ID
    matching_id: str = PAPER_0710_1869_MATCHING_ID
    rg_order: str = "lo"
    observable_scope: str = "np_only"
    kk_gluon_normalization_id: str = PAPER_0710_1869_KK_GLUON_NORMALIZATION_ID
    provenance_policy_id: str = PAPER_0710_1869_PROVENANCE_POLICY_ID
    verifier_policy_id: str = PAPER_0710_1869_VERIFIER_POLICY_ID
    seed_to_profile_mapping_policy_id: str = PAPER_0710_1869_SEED_TO_PROFILE_MAPPING_POLICY_ID
    universal_term_coefficient_policy_id: str = (
        PAPER_0710_1869_UNIVERSAL_TERM_COEFFICIENT_POLICY_ID
    )
    profile_derivation_policy_id: str = PAPER_0710_1869_PROFILE_DERIVATION_POLICY_ID
    notes: str = (
        "Structural paper-mode conventions only. "
        "This package does not expose the repo-v1 deltaf2 observable engine. "
        "The future physical seed-to-profile bridge is frozen only as an explicit "
        "affine-per-sector schema with surfaced universal-term/coefficient defaults; "
        "it must not route through a hidden BulkMassMap surrogate."
    )

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "schema_id",
            require_known_schema_id(
                "schema_id",
                self.schema_id,
                expected=PAPER_0710_1869_CONVENTIONS_SCHEMA_ID,
            ),
        )
        object.__setattr__(self, "mode_id", require_nonempty_identifier("mode_id", self.mode_id))
        object.__setattr__(self, "paper_id", require_nonempty_identifier("paper_id", self.paper_id))
        object.__setattr__(
            self,
            "operator_basis_id",
            require_known_schema_id(
                "operator_basis_id",
                self.operator_basis_id,
                expected=PAPER_0710_1869_OPERATOR_BASIS_ID,
            ),
        )
        object.__setattr__(
            self,
            "matching_id",
            require_known_schema_id(
                "matching_id",
                self.matching_id,
                expected=PAPER_0710_1869_MATCHING_ID,
            ),
        )
        object.__setattr__(
            self,
            "kk_gluon_normalization_id",
            require_known_schema_id(
                "kk_gluon_normalization_id",
                self.kk_gluon_normalization_id,
                expected=PAPER_0710_1869_KK_GLUON_NORMALIZATION_ID,
            ),
        )
        object.__setattr__(
            self,
            "provenance_policy_id",
            require_known_schema_id(
                "provenance_policy_id",
                self.provenance_policy_id,
                expected=PAPER_0710_1869_PROVENANCE_POLICY_ID,
            ),
        )
        object.__setattr__(
            self,
            "verifier_policy_id",
            require_known_schema_id(
                "verifier_policy_id",
                self.verifier_policy_id,
                expected=PAPER_0710_1869_VERIFIER_POLICY_ID,
            ),
        )
        object.__setattr__(
            self,
            "seed_to_profile_mapping_policy_id",
            require_known_schema_id(
                "seed_to_profile_mapping_policy_id",
                self.seed_to_profile_mapping_policy_id,
                expected=PAPER_0710_1869_SEED_TO_PROFILE_MAPPING_POLICY_ID,
            ),
        )
        object.__setattr__(
            self,
            "universal_term_coefficient_policy_id",
            require_known_schema_id(
                "universal_term_coefficient_policy_id",
                self.universal_term_coefficient_policy_id,
                expected=PAPER_0710_1869_UNIVERSAL_TERM_COEFFICIENT_POLICY_ID,
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
        object.__setattr__(self, "notes", require_nonempty_identifier("notes", self.notes))
        require_member("mode_id", self.mode_id, (PAPER_0710_1869_MODE_ID,))
        require_member("rg_order", self.rg_order, ("lo", "nlo"))
        require_member("observable_scope", self.observable_scope, ("np_only",))

    def as_dict(self) -> dict[str, str]:
        """Return a stable mapping representation."""
        return asdict(self)


def default_paper_0710_1869_conventions() -> Paper07101869Conventions:
    """Return the default closed convention bundle."""
    return Paper07101869Conventions()
