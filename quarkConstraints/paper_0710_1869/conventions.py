"""Canonical conventions for the paper-facing 0710.1869 mode."""

from __future__ import annotations

from dataclasses import asdict, dataclass

from .validation import (
    require_known_schema_id,
    require_member,
    require_nonempty_identifier,
)

PAPER_0710_1869_MODE_ID = "paper_0710_1869"
PAPER_0710_1869_PAPER_ID = "arXiv:0710.1869"
PAPER_0710_1869_CONVENTIONS_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.conventions.v1"
)


@dataclass(frozen=True)
class Paper07101869Conventions:
    """Closed convention bundle for the canonical paper-mode path."""

    schema_id: str = PAPER_0710_1869_CONVENTIONS_SCHEMA_ID
    mode_id: str = PAPER_0710_1869_MODE_ID
    paper_id: str = PAPER_0710_1869_PAPER_ID
    operator_basis_id: str = "kk_gluon_tree_np_only.v1"
    matching_id: str = "eft_matching.lo.v1"
    rg_order: str = "lo"
    observable_scope: str = "np_only"
    kk_gluon_normalization_id: str = "explicit_mu_gs.normalization.v1"
    provenance_policy_id: str = "deterministic_bundle.required.v1"
    verifier_policy_id: str = "independent_verifier.required.v1"
    notes: str = (
        "Structural paper-mode conventions only. "
        "This package does not expose the repo-v1 deltaf2 observable engine."
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
            require_nonempty_identifier("operator_basis_id", self.operator_basis_id),
        )
        object.__setattr__(
            self, "matching_id", require_nonempty_identifier("matching_id", self.matching_id)
        )
        object.__setattr__(
            self,
            "kk_gluon_normalization_id",
            require_nonempty_identifier(
                "kk_gluon_normalization_id", self.kk_gluon_normalization_id
            ),
        )
        object.__setattr__(
            self,
            "provenance_policy_id",
            require_nonempty_identifier("provenance_policy_id", self.provenance_policy_id),
        )
        object.__setattr__(
            self,
            "verifier_policy_id",
            require_nonempty_identifier("verifier_policy_id", self.verifier_policy_id),
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
