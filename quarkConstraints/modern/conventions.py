"""Lane-level conventions for the modern quark registry."""

from __future__ import annotations

from dataclasses import asdict, dataclass

MODERN_LANE_ID = "modern"
MODERN_CONVENTIONS_SCHEMA_ID = "quarkConstraints.modern.conventions.v1"
MODERN_INPUT_REGISTRY_SCHEMA_ID = "quarkConstraints.modern.inputs.registry.v1"
MODERN_INPUT_NAMESPACE = "quarkConstraints.modern.inputs"
MODERN_PROVENANCE_NAMESPACE = "quarkConstraints.modern.provenance"
MODERN_BUNDLE_ID_TEMPLATE = f"{MODERN_INPUT_NAMESPACE}.{{family_id}}.{{bundle_slug}}.v1"
MODERN_PROVENANCE_ID_TEMPLATE = (
    f"{MODERN_PROVENANCE_NAMESPACE}.{{family_id}}.{{bundle_slug}}.v1"
)
MODERN_STRICT_PAPER_BUNDLE_FAMILY_ID = "strict_paper_inputs"
MODERN_DEFAULT_BUNDLE_FAMILY_ID = "modern_default_inputs"
MODERN_STRICT_PAPER_BUNDLE_SLUG = "paper_sourced"
MODERN_DEFAULT_BUNDLE_SLUG = "default"
MODERN_STRICT_PAPER_SOURCE_SNAPSHOT_ID = "quarkConstraints.modern.inputs.strict_paper_snapshot.v1"
MODERN_STRICT_PAPER_SOURCE_SNAPSHOT_NAMES = (
    "table_i_snapshot",
    "eq3_example_snapshot",
)
MODERN_STRICT_PAPER_SOURCE_SCHEMA_IDS = (
    "quarkConstraints.paper_0710_1869.table_i_inputs.v1",
    "quarkConstraints.paper_0710_1869.eq3_example.v1",
)
MODERN_STRICT_PAPER_RESOLUTION_POLICY_ID = (
    "quarkConstraints.modern.inputs.strict_paper_inputs.lazy_resolver.v1"
)
MODERN_STRICT_PAPER_INPUT_BUNDLE_ID = MODERN_BUNDLE_ID_TEMPLATE.format(
    family_id=MODERN_STRICT_PAPER_BUNDLE_FAMILY_ID,
    bundle_slug=MODERN_STRICT_PAPER_BUNDLE_SLUG,
)
MODERN_STRICT_PAPER_INPUT_PROVENANCE_ID = MODERN_PROVENANCE_ID_TEMPLATE.format(
    family_id=MODERN_STRICT_PAPER_BUNDLE_FAMILY_ID,
    bundle_slug=MODERN_STRICT_PAPER_BUNDLE_SLUG,
)
MODERN_DEFAULT_INPUT_BUNDLE_ID = MODERN_BUNDLE_ID_TEMPLATE.format(
    family_id=MODERN_DEFAULT_BUNDLE_FAMILY_ID,
    bundle_slug=MODERN_DEFAULT_BUNDLE_SLUG,
)
MODERN_DEFAULT_INPUT_PROVENANCE_ID = MODERN_PROVENANCE_ID_TEMPLATE.format(
    family_id=MODERN_DEFAULT_BUNDLE_FAMILY_ID,
    bundle_slug=MODERN_DEFAULT_BUNDLE_SLUG,
)


def _require_text(name: str, value: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a non-empty string")
    return value.strip()


def _require_exact(name: str, value: str, *, expected: str) -> str:
    if value != expected:
        raise ValueError(f"{name} must be exactly {expected!r}")
    return value


def _require_member(name: str, value: str, allowed: tuple[str, ...]) -> str:
    if value not in allowed:
        raise ValueError(f"{name} must be one of {allowed!r}")
    return value


def build_modern_bundle_id(family_id: str, bundle_slug: str) -> str:
    """Return the canonical modern bundle identifier for one family."""
    family = _require_text("family_id", family_id)
    slug = _require_text("bundle_slug", bundle_slug)
    return MODERN_BUNDLE_ID_TEMPLATE.format(family_id=family, bundle_slug=slug)


def build_modern_provenance_id(family_id: str, bundle_slug: str) -> str:
    """Return the canonical modern provenance identifier for one family."""
    family = _require_text("family_id", family_id)
    slug = _require_text("bundle_slug", bundle_slug)
    return MODERN_PROVENANCE_ID_TEMPLATE.format(family_id=family, bundle_slug=slug)


@dataclass(frozen=True)
class ModernLaneConventions:
    """Frozen conventions for the modern production lane."""

    schema_id: str = MODERN_CONVENTIONS_SCHEMA_ID
    lane_id: str = MODERN_LANE_ID
    input_registry_schema_id: str = MODERN_INPUT_REGISTRY_SCHEMA_ID
    input_namespace: str = MODERN_INPUT_NAMESPACE
    provenance_namespace: str = MODERN_PROVENANCE_NAMESPACE
    bundle_id_template: str = MODERN_BUNDLE_ID_TEMPLATE
    provenance_id_template: str = MODERN_PROVENANCE_ID_TEMPLATE
    strict_paper_bundle_family_id: str = MODERN_STRICT_PAPER_BUNDLE_FAMILY_ID
    modern_default_bundle_family_id: str = MODERN_DEFAULT_BUNDLE_FAMILY_ID
    notes: str = (
        "Modern-lane registry conventions only. strict_paper_inputs may adapt "
        "paper-lane sourced inputs, while modern_default_inputs must remain a "
        "separate modern namespace."
    )

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "schema_id",
            _require_exact(
                "schema_id",
                self.schema_id,
                expected=MODERN_CONVENTIONS_SCHEMA_ID,
            ),
        )
        object.__setattr__(self, "lane_id", _require_exact("lane_id", self.lane_id, expected=MODERN_LANE_ID))
        object.__setattr__(
            self,
            "input_registry_schema_id",
            _require_exact(
                "input_registry_schema_id",
                self.input_registry_schema_id,
                expected=MODERN_INPUT_REGISTRY_SCHEMA_ID,
            ),
        )
        object.__setattr__(
            self,
            "input_namespace",
            _require_exact("input_namespace", self.input_namespace, expected=MODERN_INPUT_NAMESPACE),
        )
        object.__setattr__(
            self,
            "provenance_namespace",
            _require_exact(
                "provenance_namespace",
                self.provenance_namespace,
                expected=MODERN_PROVENANCE_NAMESPACE,
            ),
        )
        object.__setattr__(
            self,
            "bundle_id_template",
            _require_exact(
                "bundle_id_template",
                self.bundle_id_template,
                expected=MODERN_BUNDLE_ID_TEMPLATE,
            ),
        )
        object.__setattr__(
            self,
            "provenance_id_template",
            _require_exact(
                "provenance_id_template",
                self.provenance_id_template,
                expected=MODERN_PROVENANCE_ID_TEMPLATE,
            ),
        )
        object.__setattr__(
            self,
            "strict_paper_bundle_family_id",
            _require_exact(
                "strict_paper_bundle_family_id",
                self.strict_paper_bundle_family_id,
                expected=MODERN_STRICT_PAPER_BUNDLE_FAMILY_ID,
            ),
        )
        object.__setattr__(
            self,
            "modern_default_bundle_family_id",
            _require_exact(
                "modern_default_bundle_family_id",
                self.modern_default_bundle_family_id,
                expected=MODERN_DEFAULT_BUNDLE_FAMILY_ID,
            ),
        )
        object.__setattr__(self, "notes", _require_text("notes", self.notes))
        _require_member("lane_id", self.lane_id, (MODERN_LANE_ID,))
        _require_member(
            "strict_paper_bundle_family_id",
            self.strict_paper_bundle_family_id,
            (MODERN_STRICT_PAPER_BUNDLE_FAMILY_ID,),
        )
        _require_member(
            "modern_default_bundle_family_id",
            self.modern_default_bundle_family_id,
            (MODERN_DEFAULT_BUNDLE_FAMILY_ID,),
        )

    def as_dict(self) -> dict[str, str]:
        """Return a stable mapping representation."""
        return asdict(self)


def default_modern_lane_conventions() -> ModernLaneConventions:
    """Return the frozen default modern-lane conventions."""
    return ModernLaneConventions()
