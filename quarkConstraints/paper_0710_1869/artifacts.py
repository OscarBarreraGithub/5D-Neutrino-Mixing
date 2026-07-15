"""Versioned artifact schemas and canonical export builders for ``paper_0710_1869``."""

from __future__ import annotations

import json
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace
from math import isclose, isfinite
from pathlib import Path
from typing import Any, TypeAlias

_REPO_ROOT = Path(__file__).resolve().parents[2]

PAPER_MODE = "paper_0710_1869"
ARTIFACT_SCHEMA_VERSION = 1
WILSON_BUNDLE_SCHEMA = "paper_0710_1869.wilson_bundle"
HADRONIC_BUNDLE_SCHEMA = "paper_0710_1869.hadronic_bundle"
OBSERVABLE_BUNDLE_SCHEMA = "paper_0710_1869.observable_bundle"
PROVENANCE_BUNDLE_SCHEMA = "paper_0710_1869.provenance_bundle"

DEFAULT_KAON_ARTIFACT_POINT_ID = "pr1.table_i_eq3_example.v1.kaon.central.np_only.v1"
DEFAULT_KAON_WILSON_BUNDLE_ID = "wilson.kaon.pr1.table_i_eq3_example.central.v1"
DEFAULT_KAON_HADRONIC_BUNDLE_ID = "hadronic.kaon.pr5a.v1"
DEFAULT_KAON_OBSERVABLE_BUNDLE_ID = "observable.kaon.pr1.table_i_eq3_example.central.v1"
DEFAULT_KAON_PROVENANCE_BUNDLE_ID = "provenance.kaon.pr1.table_i_eq3_example.central.v1"
DEFAULT_KAON_ARTIFACT_DIR = _REPO_ROOT / "results" / "paper_0710_1869" / "default_kaon"
STRICT_PAPER_ARTIFACT_POINT_ID = "pr1.table_i_eq3_example.strict_paper.v1"
STRICT_PAPER_WILSON_BUNDLE_ID = "wilson.kaon.pr1.table_i_eq3_example.strict_paper.v1"
STRICT_PAPER_HADRONIC_BUNDLE_ID = "hadronic.kaon.pr1.table_i_eq3_example.strict_paper.v1"
STRICT_PAPER_OBSERVABLE_BUNDLE_ID = "observable.kaon.pr1.table_i_eq3_example.strict_paper.v1"
STRICT_PAPER_PROVENANCE_BUNDLE_ID = "provenance.kaon.pr1.table_i_eq3_example.strict_paper.v1"
STRICT_PAPER_ARTIFACT_DIR = _REPO_ROOT / "results" / "paper_0710_1869" / "strict_paper_kaon"
DEFAULT_KAON_SUPPORTED_OPERATORS = ("Q1_VLL", "Q1_VRR")
DEFAULT_KAON_UNSUPPORTED_OPERATORS = ("Q4_LR", "Q5_LR")
DELTA_M_RELATION_ID = "delta_m_k_np.equals.2_re_m12_k_np.v1"
DEFAULT_BAG_PARAMETER_CONVERSION_FORMULA_ID = (
    "b_k_mu_had.equals.hat_b_k_rgi_times_alpha_s_mu_had_pow_2_over_beta0_lo.v1"
)
DEFAULT_KAON_FROZEN_MATCHING_Q1_VLL = complex(
    -1.1004498898491606e-12,
    -7.601544537269597e-13,
)
DEFAULT_KAON_FROZEN_RG_Q1_VLL = complex(
    -1.509786848232995e-12,
    -1.0429109107548845e-12,
)
DEFAULT_KAON_FROZEN_BAG_PARAMETER_CONVERSION_ALPHA_S_MU_HAD = 0.26770024052060604
DEFAULT_KAON_FROZEN_M12_K_NP = complex(
    -1.3495753042583394e-14,
    -9.322420653906486e-15,
)
DEFAULT_KAON_FROZEN_DELTA_M_K_NP_GEV = -2.6991506085166787e-14
DEFAULT_KAON_ARTIFACT_FILENAMES = {
    "wilsons": "wilsons.json",
    "observables": "observables.json",
    "wilson": "wilsons.json",
    "hadronic": "hadronic.json",
    "observable": "observables.json",
    "provenance": "provenance.json",
}

DEFAULT_PAPER_PROVENANCE_ID = "paper.0710.1869.v1"
DEFAULT_MATCHING_PROVENANCE_ID = "matching.kaon.default.v1"
DEFAULT_RG_PROVENANCE_ID = "rg.kaon.default_lo.v1"
STRICT_PAPER_CONVENTIONS_PROVENANCE_ID = "strict.paper.conventions.v1"
STRICT_PAPER_MATCHING_PROVENANCE_ID = "strict.paper.matching.v1"
STRICT_PAPER_RG_PROVENANCE_ID = "strict.paper.rg.v1"
STRICT_PAPER_HADRONIC_SOURCE_ID = "strict.paper.hadronic.source.v1"
STRICT_PAPER_MASS_SOURCE_ID = "strict.paper.mass.source.v1"
STRICT_PAPER_DECAY_CONSTANT_SOURCE_ID = "strict.paper.decay.constant.source.v1"
STRICT_PAPER_BAG_SOURCE_ID = "strict.paper.bag.source.v1"

JSONDict: TypeAlias = dict[str, Any]


class ArtifactSchemaError(ValueError):
    """Raised when a paper artifact violates the frozen JSON contract."""


def _require_mapping(value: Any, *, context: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ArtifactSchemaError(f"{context} must be a JSON object")
    return value


def _require_sequence(value: Any, *, context: str) -> Sequence[Any]:
    if isinstance(value, (str, bytes, bytearray)) or not isinstance(value, Sequence):
        raise ArtifactSchemaError(f"{context} must be a JSON array")
    return value


def _require_text(value: Any, *, context: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ArtifactSchemaError(f"{context} must be a non-empty string")
    return value


def _require_optional_text(value: Any, *, context: str) -> str | None:
    if value is None:
        return None
    return _require_text(value, context=context)


def _require_int(value: Any, *, context: str) -> int:
    if not isinstance(value, int) or isinstance(value, bool):
        raise ArtifactSchemaError(f"{context} must be an integer")
    return value


def _require_float(value: Any, *, context: str) -> float:
    if not isinstance(value, (int, float)) or isinstance(value, bool):
        raise ArtifactSchemaError(f"{context} must be a finite float")
    coerced = float(value)
    if not isfinite(coerced):
        raise ArtifactSchemaError(f"{context} must be a finite float")
    return coerced


def _require_optional_float(value: Any, *, context: str) -> float | None:
    if value is None:
        return None
    return _require_float(value, context=context)


def _dump_json(payload: JSONDict) -> str:
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"


def _ensure_unique_strings(values: Sequence[str], *, context: str) -> None:
    seen: set[str] = set()
    duplicates = {value for value in values if value in seen or seen.add(value)}
    if duplicates:
        dup_list = ", ".join(sorted(duplicates))
        raise ArtifactSchemaError(f"{context} contains duplicate entries: {dup_list}")


def _scales_to_dict(scales: Sequence["ArtifactScale"]) -> list[JSONDict]:
    return [item.to_dict() for item in scales]


def _ids_to_list(values: Sequence[str]) -> list[str]:
    return list(values)


def _parse_scales(value: Any, *, context: str) -> tuple["ArtifactScale", ...]:
    entries = _require_sequence(value, context=context)
    return tuple(
        ArtifactScale.from_dict(_require_mapping(item, context=context)) for item in entries
    )


def _parse_string_ids(value: Any, *, context: str) -> tuple[str, ...]:
    entries = _require_sequence(value, context=context)
    resolved = tuple(_require_text(item, context=context) for item in entries)
    _ensure_unique_strings(resolved, context=context)
    return resolved


def _parse_optional_string_ids(value: Any, *, context: str) -> tuple[str, ...]:
    if value is None:
        return ()
    return _parse_string_ids(value, context=context)


def _qcd_beta0_from_n_f(n_f: int) -> float:
    return 11.0 - (2.0 / 3.0) * float(n_f)


def _n_f_at_scale(
    mu_gev: float,
    thresholds: tuple[tuple[float, int, int], ...],
) -> int:
    if not thresholds:
        return 5
    n_f = thresholds[0][1]
    for mass, _, n_f_above in thresholds:
        if mu_gev >= mass:
            n_f = n_f_above
        else:
            break
    return n_f


@dataclass(frozen=True)
class ArtifactMetadata:
    """Frozen bundle header shared by all paper-facing exports."""

    schema_name: str
    schema_version: int
    bundle_id: str
    point_id: str
    mode: str = PAPER_MODE

    def __post_init__(self) -> None:
        _require_text(self.schema_name, context="metadata.schema_name")
        _require_int(self.schema_version, context="metadata.schema_version")
        if self.schema_version < 1:
            raise ArtifactSchemaError("metadata.schema_version must be positive")
        _require_text(self.bundle_id, context="metadata.bundle_id")
        _require_text(self.point_id, context="metadata.point_id")
        if self.mode != PAPER_MODE:
            raise ArtifactSchemaError(f"metadata.mode must be {PAPER_MODE!r}")

    @classmethod
    def create(
        cls,
        *,
        schema_name: str,
        bundle_id: str,
        point_id: str,
        schema_version: int = ARTIFACT_SCHEMA_VERSION,
    ) -> "ArtifactMetadata":
        return cls(
            schema_name=schema_name,
            schema_version=schema_version,
            bundle_id=bundle_id,
            point_id=point_id,
        )

    def to_dict(self) -> JSONDict:
        return {
            "schema_name": self.schema_name,
            "schema_version": self.schema_version,
            "bundle_id": self.bundle_id,
            "point_id": self.point_id,
            "mode": self.mode,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "ArtifactMetadata":
        mapping = _require_mapping(data, context="metadata")
        return cls(
            schema_name=_require_text(mapping.get("schema_name"), context="metadata.schema_name"),
            schema_version=_require_int(
                mapping.get("schema_version"), context="metadata.schema_version"
            ),
            bundle_id=_require_text(mapping.get("bundle_id"), context="metadata.bundle_id"),
            point_id=_require_text(mapping.get("point_id"), context="metadata.point_id"),
            mode=_require_text(mapping.get("mode"), context="metadata.mode"),
        )


@dataclass(frozen=True)
class ArtifactScale:
    """Named scale object carried across paper exports."""

    name: str
    role: str
    value_gev: float

    def __post_init__(self) -> None:
        object.__setattr__(self, "name", _require_text(self.name, context="scale.name"))
        object.__setattr__(self, "role", _require_text(self.role, context="scale.role"))
        object.__setattr__(
            self,
            "value_gev",
            _require_float(self.value_gev, context="scale.value_gev"),
        )
        if self.value_gev <= 0.0:
            raise ArtifactSchemaError("scale.value_gev must be positive")

    def to_dict(self) -> JSONDict:
        return {
            "name": self.name,
            "role": self.role,
            "value_gev": self.value_gev,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "ArtifactScale":
        mapping = _require_mapping(data, context="scale")
        return cls(
            name=_require_text(mapping.get("name"), context="scale.name"),
            role=_require_text(mapping.get("role"), context="scale.role"),
            value_gev=_require_float(mapping.get("value_gev"), context="scale.value_gev"),
        )


@dataclass(frozen=True)
class ComplexValue:
    """JSON-safe complex scalar container."""

    real: float
    imag: float = 0.0

    def __post_init__(self) -> None:
        _require_float(self.real, context="complex.real")
        _require_float(self.imag, context="complex.imag")

    @classmethod
    def from_complex(cls, value: complex) -> "ComplexValue":
        return cls(real=float(value.real), imag=float(value.imag))

    def to_complex(self) -> complex:
        return complex(self.real, self.imag)

    def to_dict(self) -> JSONDict:
        return {"real": self.real, "imag": self.imag}

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "ComplexValue":
        mapping = _require_mapping(data, context="complex")
        return cls(
            real=_require_float(mapping.get("real"), context="complex.real"),
            imag=_require_float(mapping.get("imag", 0.0), context="complex.imag"),
        )


@dataclass(frozen=True)
class WilsonCoefficientRecord:
    """One exported Wilson coefficient entry."""

    sector: str
    system: str
    operator: str
    value: ComplexValue

    def __post_init__(self) -> None:
        _require_text(self.sector, context="coefficient.sector")
        _require_text(self.system, context="coefficient.system")
        _require_text(self.operator, context="coefficient.operator")
        if not isinstance(self.value, ComplexValue):
            raise ArtifactSchemaError("coefficient.value must be a ComplexValue")

    def to_dict(self) -> JSONDict:
        return {
            "sector": self.sector,
            "system": self.system,
            "operator": self.operator,
            "value": self.value.to_dict(),
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "WilsonCoefficientRecord":
        mapping = _require_mapping(data, context="coefficient")
        return cls(
            sector=_require_text(mapping.get("sector"), context="coefficient.sector"),
            system=_require_text(mapping.get("system"), context="coefficient.system"),
            operator=_require_text(mapping.get("operator"), context="coefficient.operator"),
            value=ComplexValue.from_dict(
                _require_mapping(mapping.get("value"), context="coefficient.value")
            ),
        )


@dataclass(frozen=True)
class ObservableRecord:
    """One exported observable entry."""

    name: str
    system: str
    value: float
    units: str
    reference_bound: float | None = None
    ratio_to_bound: float | None = None
    passes: bool | None = None

    def __post_init__(self) -> None:
        _require_text(self.name, context="observable.name")
        _require_text(self.system, context="observable.system")
        _require_text(self.units, context="observable.units")
        _require_float(self.value, context="observable.value")
        if self.reference_bound is not None:
            _require_float(self.reference_bound, context="observable.reference_bound")
        if self.ratio_to_bound is not None:
            _require_float(self.ratio_to_bound, context="observable.ratio_to_bound")
        if self.passes is not None and not isinstance(self.passes, bool):
            raise ArtifactSchemaError("observable.passes must be a boolean when provided")

    def to_dict(self) -> JSONDict:
        payload: JSONDict = {
            "name": self.name,
            "system": self.system,
            "value": self.value,
            "units": self.units,
        }
        if self.reference_bound is not None:
            payload["reference_bound"] = self.reference_bound
        if self.ratio_to_bound is not None:
            payload["ratio_to_bound"] = self.ratio_to_bound
        if self.passes is not None or self.passes is False:
            payload["passes"] = self.passes
        return payload

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "ObservableRecord":
        mapping = _require_mapping(data, context="observable")
        passes = mapping.get("passes")
        if passes is not None and not isinstance(passes, bool):
            raise ArtifactSchemaError("observable.passes must be a boolean when provided")
        return cls(
            name=_require_text(mapping.get("name"), context="observable.name"),
            system=_require_text(mapping.get("system"), context="observable.system"),
            value=_require_float(mapping.get("value"), context="observable.value"),
            units=_require_text(mapping.get("units"), context="observable.units"),
            reference_bound=(
                None
                if "reference_bound" not in mapping
                else _require_float(
                    mapping.get("reference_bound"), context="observable.reference_bound"
                )
            ),
            ratio_to_bound=(
                None
                if "ratio_to_bound" not in mapping
                else _require_float(
                    mapping.get("ratio_to_bound"), context="observable.ratio_to_bound"
                )
            ),
            passes=passes,
        )


@dataclass(frozen=True)
class ProvenanceRecord:
    """One provenance entry referenced by artifact bundles."""

    record_id: str
    category: str
    label: str
    version: str
    source: str
    citation: str | None = None
    digest: str | None = None

    def __post_init__(self) -> None:
        _require_text(self.record_id, context="provenance.record_id")
        _require_text(self.category, context="provenance.category")
        _require_text(self.label, context="provenance.label")
        _require_text(self.version, context="provenance.version")
        _require_text(self.source, context="provenance.source")
        if self.citation is not None:
            _require_text(self.citation, context="provenance.citation")
        if self.digest is not None:
            _require_text(self.digest, context="provenance.digest")

    def to_dict(self) -> JSONDict:
        payload: JSONDict = {
            "record_id": self.record_id,
            "category": self.category,
            "label": self.label,
            "version": self.version,
            "source": self.source,
        }
        if self.citation is not None:
            payload["citation"] = self.citation
        if self.digest is not None:
            payload["digest"] = self.digest
        return payload

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "ProvenanceRecord":
        mapping = _require_mapping(data, context="provenance")
        return cls(
            record_id=_require_text(mapping.get("record_id"), context="provenance.record_id"),
            category=_require_text(mapping.get("category"), context="provenance.category"),
            label=_require_text(mapping.get("label"), context="provenance.label"),
            version=_require_text(mapping.get("version"), context="provenance.version"),
            source=_require_text(mapping.get("source"), context="provenance.source"),
            citation=(
                None
                if "citation" not in mapping
                else _require_text(mapping.get("citation"), context="provenance.citation")
            ),
            digest=(
                None
                if "digest" not in mapping
                else _require_text(mapping.get("digest"), context="provenance.digest")
            ),
        )


@dataclass(frozen=True)
class ArtifactSourceRecord:
    """One cited source entry embedded in a hadronic artifact bundle."""

    source_id: str
    source_kind: str
    citation: str
    locator_label: str
    year: int
    renormalization_scheme_id: str | None = None
    scale_GeV: float | None = None
    transformation_id: str = "none"
    notes: str | None = None

    def __post_init__(self) -> None:
        _require_text(self.source_id, context="source.source_id")
        _require_text(self.source_kind, context="source.source_kind")
        _require_text(self.citation, context="source.citation")
        _require_text(self.locator_label, context="source.locator_label")
        _require_int(self.year, context="source.year")
        if self.year < 1900:
            raise ArtifactSchemaError("source.year must be >= 1900")
        if self.renormalization_scheme_id is not None:
            _require_text(
                self.renormalization_scheme_id,
                context="source.renormalization_scheme_id",
            )
        if self.scale_GeV is not None:
            value = _require_float(self.scale_GeV, context="source.scale_GeV")
            if value <= 0.0:
                raise ArtifactSchemaError("source.scale_GeV must be positive")
        _require_text(self.transformation_id, context="source.transformation_id")
        if self.notes is not None:
            _require_text(self.notes, context="source.notes")

    def to_dict(self) -> JSONDict:
        payload: JSONDict = {
            "source_id": self.source_id,
            "source_kind": self.source_kind,
            "citation": self.citation,
            "locator_label": self.locator_label,
            "year": self.year,
            "transformation_id": self.transformation_id,
        }
        if self.renormalization_scheme_id is not None:
            payload["renormalization_scheme_id"] = self.renormalization_scheme_id
        if self.scale_GeV is not None:
            payload["scale_GeV"] = self.scale_GeV
        if self.notes is not None:
            payload["notes"] = self.notes
        return payload

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "ArtifactSourceRecord":
        mapping = _require_mapping(data, context="source")
        return cls(
            source_id=_require_text(mapping.get("source_id"), context="source.source_id"),
            source_kind=_require_text(mapping.get("source_kind"), context="source.source_kind"),
            citation=_require_text(mapping.get("citation"), context="source.citation"),
            locator_label=_require_text(
                mapping.get("locator_label"), context="source.locator_label"
            ),
            year=_require_int(mapping.get("year"), context="source.year"),
            renormalization_scheme_id=_require_optional_text(
                mapping.get("renormalization_scheme_id"),
                context="source.renormalization_scheme_id",
            ),
            scale_GeV=_require_optional_float(mapping.get("scale_GeV"), context="source.scale_GeV"),
            transformation_id=_require_text(
                mapping.get("transformation_id", "none"),
                context="source.transformation_id",
            ),
            notes=_require_optional_text(mapping.get("notes"), context="source.notes"),
        )


@dataclass(frozen=True)
class WilsonArtifactBundleV1:
    """Exported Wilson-coefficient bundle for paper mode."""

    metadata: ArtifactMetadata
    operator_basis: str
    renormalization_scheme: str
    scales: tuple[ArtifactScale, ...]
    coefficients: tuple[WilsonCoefficientRecord, ...]
    operator_normalization: str | None = None
    coefficient_scale_name: str = "mu_had"
    matching_scale_name: str = "mu_match"
    matching_coefficients: tuple[WilsonCoefficientRecord, ...] = ()
    supported_operator_names: tuple[str, ...] = ()
    unsupported_operator_names: tuple[str, ...] = ()
    provenance_ids: tuple[str, ...] = ()
    notes: str | None = None

    def __post_init__(self) -> None:
        if self.metadata.schema_name != WILSON_BUNDLE_SCHEMA:
            raise ArtifactSchemaError("wilson bundle metadata.schema_name is not supported")
        if self.metadata.schema_version != ARTIFACT_SCHEMA_VERSION:
            raise ArtifactSchemaError("wilson bundle metadata.schema_version is not supported")
        _require_text(self.operator_basis, context="wilson.operator_basis")
        _require_text(self.renormalization_scheme, context="wilson.renormalization_scheme")
        if self.operator_normalization is not None:
            _require_text(
                self.operator_normalization,
                context="wilson.operator_normalization",
            )
        _require_text(self.coefficient_scale_name, context="wilson.coefficient_scale_name")
        _require_text(self.matching_scale_name, context="wilson.matching_scale_name")
        if not self.scales:
            raise ArtifactSchemaError("wilson.scales must not be empty")
        if not self.coefficients:
            raise ArtifactSchemaError("wilson.coefficients must not be empty")
        _ensure_unique_strings([item.name for item in self.scales], context="wilson.scales")
        if self.supported_operator_names:
            _ensure_unique_strings(
                self.supported_operator_names,
                context="wilson.supported_operator_names",
            )
        if self.unsupported_operator_names:
            _ensure_unique_strings(
                self.unsupported_operator_names,
                context="wilson.unsupported_operator_names",
            )
        if set(self.supported_operator_names) & set(self.unsupported_operator_names):
            raise ArtifactSchemaError(
                "wilson supported_operator_names and unsupported_operator_names must be disjoint"
            )
        supported_set = set(self.supported_operator_names)
        unsupported_set = set(self.unsupported_operator_names)
        allowed_names = supported_set | unsupported_set
        if allowed_names:
            coefficient_names = {item.operator for item in self.coefficients}
            matching_names = {item.operator for item in self.matching_coefficients}
            if not coefficient_names <= allowed_names:
                unsupported = ", ".join(sorted(coefficient_names - allowed_names))
                raise ArtifactSchemaError(
                    f"wilson.coefficients contain undeclared operators: {unsupported}"
                )
            if not matching_names <= allowed_names:
                unsupported = ", ".join(sorted(matching_names - allowed_names))
                raise ArtifactSchemaError(
                    f"wilson.matching_coefficients contain undeclared operators: {unsupported}"
                )
        _ensure_unique_strings(self.provenance_ids, context="wilson.provenance_ids")
        if self.notes is not None:
            _require_text(self.notes, context="wilson.notes")

    def to_dict(self) -> JSONDict:
        payload: JSONDict = {
            "metadata": self.metadata.to_dict(),
            "operator_basis": self.operator_basis,
            "renormalization_scheme": self.renormalization_scheme,
            "scales": _scales_to_dict(self.scales),
            "coefficients": [item.to_dict() for item in self.coefficients],
            "coefficient_scale_name": self.coefficient_scale_name,
            "matching_scale_name": self.matching_scale_name,
            "matching_coefficients": [item.to_dict() for item in self.matching_coefficients],
            "provenance_ids": _ids_to_list(self.provenance_ids),
        }
        if self.operator_normalization is not None:
            payload["operator_normalization"] = self.operator_normalization
        if self.supported_operator_names:
            payload["supported_operator_names"] = list(self.supported_operator_names)
        if self.unsupported_operator_names:
            payload["unsupported_operator_names"] = list(self.unsupported_operator_names)
        if self.notes is not None:
            payload["notes"] = self.notes
        return payload

    def to_json(self) -> str:
        return _dump_json(self.to_dict())

    def write_json(self, path: str | Path) -> Path:
        destination = Path(path)
        destination.write_text(self.to_json(), encoding="utf-8")
        return destination

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "WilsonArtifactBundleV1":
        mapping = _require_mapping(data, context="wilson_bundle")
        coefficient_entries = _require_sequence(
            mapping.get("coefficients"), context="wilson_bundle.coefficients"
        )
        matching_entries = _require_sequence(
            mapping.get("matching_coefficients", ()),
            context="wilson_bundle.matching_coefficients",
        )
        return cls(
            metadata=ArtifactMetadata.from_dict(
                _require_mapping(mapping.get("metadata"), context="wilson_bundle.metadata")
            ),
            operator_basis=_require_text(
                mapping.get("operator_basis"), context="wilson_bundle.operator_basis"
            ),
            renormalization_scheme=_require_text(
                mapping.get("renormalization_scheme"),
                context="wilson_bundle.renormalization_scheme",
            ),
            scales=_parse_scales(mapping.get("scales"), context="wilson_bundle.scales"),
            coefficients=tuple(
                WilsonCoefficientRecord.from_dict(
                    _require_mapping(item, context="wilson_bundle.coefficients")
                )
                for item in coefficient_entries
            ),
            operator_normalization=_require_optional_text(
                mapping.get("operator_normalization"),
                context="wilson_bundle.operator_normalization",
            ),
            coefficient_scale_name=_require_text(
                mapping.get("coefficient_scale_name", "mu_had"),
                context="wilson_bundle.coefficient_scale_name",
            ),
            matching_scale_name=_require_text(
                mapping.get("matching_scale_name", "mu_match"),
                context="wilson_bundle.matching_scale_name",
            ),
            matching_coefficients=tuple(
                WilsonCoefficientRecord.from_dict(
                    _require_mapping(item, context="wilson_bundle.matching_coefficients")
                )
                for item in matching_entries
            ),
            supported_operator_names=_parse_optional_string_ids(
                mapping.get("supported_operator_names"),
                context="wilson_bundle.supported_operator_names",
            ),
            unsupported_operator_names=_parse_optional_string_ids(
                mapping.get("unsupported_operator_names"),
                context="wilson_bundle.unsupported_operator_names",
            ),
            provenance_ids=(
                ()
                if "provenance_ids" not in mapping
                else _parse_string_ids(
                    mapping.get("provenance_ids"), context="wilson_bundle.provenance_ids"
                )
            ),
            notes=_require_optional_text(mapping.get("notes"), context="wilson_bundle.notes"),
        )

    @classmethod
    def from_json(cls, payload: str | bytes) -> "WilsonArtifactBundleV1":
        return cls.from_dict(json.loads(payload))

    @classmethod
    def read_json(cls, path: str | Path) -> "WilsonArtifactBundleV1":
        return cls.from_json(Path(path).read_text(encoding="utf-8"))


@dataclass(frozen=True)
class ObservableArtifactBundleV1:
    """Exported observable bundle for paper mode."""

    metadata: ArtifactMetadata
    source_wilson_bundle_id: str
    interpretation: str
    operator_basis: str
    renormalization_scheme: str
    scales: tuple[ArtifactScale, ...]
    observables: tuple[ObservableRecord, ...]
    operator_normalization: str | None = None
    source_hadronic_bundle_id: str | None = None
    observable_scope_id: str | None = None
    m12_observable_id: str | None = None
    delta_m_observable_id: str | None = None
    delta_m_relation_id: str | None = None
    m12_value: ComplexValue | None = None
    supported_operator_names: tuple[str, ...] = ()
    unsupported_operator_names: tuple[str, ...] = ()
    provenance_ids: tuple[str, ...] = ()
    notes: str | None = None

    def __post_init__(self) -> None:
        if self.metadata.schema_name != OBSERVABLE_BUNDLE_SCHEMA:
            raise ArtifactSchemaError("observable bundle metadata.schema_name is not supported")
        if self.metadata.schema_version != ARTIFACT_SCHEMA_VERSION:
            raise ArtifactSchemaError("observable bundle metadata.schema_version is not supported")
        _require_text(
            self.source_wilson_bundle_id,
            context="observable.source_wilson_bundle_id",
        )
        _require_text(self.interpretation, context="observable.interpretation")
        _require_text(self.operator_basis, context="observable.operator_basis")
        _require_text(
            self.renormalization_scheme,
            context="observable.renormalization_scheme",
        )
        if self.operator_normalization is not None:
            _require_text(
                self.operator_normalization,
                context="observable.operator_normalization",
            )
        if self.source_hadronic_bundle_id is not None:
            _require_text(
                self.source_hadronic_bundle_id,
                context="observable.source_hadronic_bundle_id",
            )
        if self.observable_scope_id is not None:
            _require_text(
                self.observable_scope_id,
                context="observable.observable_scope_id",
            )
        if self.m12_observable_id is not None:
            _require_text(self.m12_observable_id, context="observable.m12_observable_id")
        if self.delta_m_observable_id is not None:
            _require_text(
                self.delta_m_observable_id,
                context="observable.delta_m_observable_id",
            )
        if self.delta_m_relation_id is not None:
            _require_text(
                self.delta_m_relation_id,
                context="observable.delta_m_relation_id",
            )
        if self.m12_value is not None and not isinstance(self.m12_value, ComplexValue):
            raise ArtifactSchemaError("observable.m12_value must be a ComplexValue when provided")
        if not self.scales:
            raise ArtifactSchemaError("observable.scales must not be empty")
        if not self.observables:
            raise ArtifactSchemaError("observable.observables must not be empty")
        _ensure_unique_strings([item.name for item in self.scales], context="observable.scales")
        if self.supported_operator_names:
            _ensure_unique_strings(
                self.supported_operator_names,
                context="observable.supported_operator_names",
            )
        if self.unsupported_operator_names:
            _ensure_unique_strings(
                self.unsupported_operator_names,
                context="observable.unsupported_operator_names",
            )
        _ensure_unique_strings(self.provenance_ids, context="observable.provenance_ids")
        if self.notes is not None:
            _require_text(self.notes, context="observable.notes")
        if self.m12_value is not None and self.delta_m_observable_id is not None:
            delta_m_records = [
                record.value
                for record in self.observables
                if record.name == self.delta_m_observable_id
            ]
            if delta_m_records:
                expected_delta_m = 2.0 * self.m12_value.real
                if not isclose(delta_m_records[0], expected_delta_m, rel_tol=0.0, abs_tol=1e-18):
                    raise ArtifactSchemaError(
                        "observable delta_m claim must equal 2 * Re(m12_value)"
                    )

    def to_dict(self) -> JSONDict:
        payload: JSONDict = {
            "metadata": self.metadata.to_dict(),
            "source_wilson_bundle_id": self.source_wilson_bundle_id,
            "interpretation": self.interpretation,
            "operator_basis": self.operator_basis,
            "renormalization_scheme": self.renormalization_scheme,
            "scales": _scales_to_dict(self.scales),
            "observables": [item.to_dict() for item in self.observables],
            "provenance_ids": _ids_to_list(self.provenance_ids),
        }
        if self.operator_normalization is not None:
            payload["operator_normalization"] = self.operator_normalization
        if self.source_hadronic_bundle_id is not None:
            payload["source_hadronic_bundle_id"] = self.source_hadronic_bundle_id
        if self.observable_scope_id is not None:
            payload["observable_scope_id"] = self.observable_scope_id
        if self.m12_observable_id is not None:
            payload["m12_observable_id"] = self.m12_observable_id
        if self.delta_m_observable_id is not None:
            payload["delta_m_observable_id"] = self.delta_m_observable_id
        if self.delta_m_relation_id is not None:
            payload["delta_m_relation_id"] = self.delta_m_relation_id
        if self.m12_value is not None:
            payload["m12_value"] = self.m12_value.to_dict()
        if self.supported_operator_names:
            payload["supported_operator_names"] = list(self.supported_operator_names)
        if self.unsupported_operator_names:
            payload["unsupported_operator_names"] = list(self.unsupported_operator_names)
        if self.notes is not None:
            payload["notes"] = self.notes
        return payload

    def to_json(self) -> str:
        return _dump_json(self.to_dict())

    def write_json(self, path: str | Path) -> Path:
        destination = Path(path)
        destination.write_text(self.to_json(), encoding="utf-8")
        return destination

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "ObservableArtifactBundleV1":
        mapping = _require_mapping(data, context="observable_bundle")
        observable_entries = _require_sequence(
            mapping.get("observables"), context="observable_bundle.observables"
        )
        return cls(
            metadata=ArtifactMetadata.from_dict(
                _require_mapping(mapping.get("metadata"), context="observable_bundle.metadata")
            ),
            source_wilson_bundle_id=_require_text(
                mapping.get("source_wilson_bundle_id"),
                context="observable_bundle.source_wilson_bundle_id",
            ),
            interpretation=_require_text(
                mapping.get("interpretation"), context="observable_bundle.interpretation"
            ),
            operator_basis=_require_text(
                mapping.get("operator_basis"), context="observable_bundle.operator_basis"
            ),
            renormalization_scheme=_require_text(
                mapping.get("renormalization_scheme"),
                context="observable_bundle.renormalization_scheme",
            ),
            scales=_parse_scales(mapping.get("scales"), context="observable_bundle.scales"),
            observables=tuple(
                ObservableRecord.from_dict(
                    _require_mapping(item, context="observable_bundle.observables")
                )
                for item in observable_entries
            ),
            operator_normalization=_require_optional_text(
                mapping.get("operator_normalization"),
                context="observable_bundle.operator_normalization",
            ),
            source_hadronic_bundle_id=_require_optional_text(
                mapping.get("source_hadronic_bundle_id"),
                context="observable_bundle.source_hadronic_bundle_id",
            ),
            observable_scope_id=_require_optional_text(
                mapping.get("observable_scope_id"),
                context="observable_bundle.observable_scope_id",
            ),
            m12_observable_id=_require_optional_text(
                mapping.get("m12_observable_id"),
                context="observable_bundle.m12_observable_id",
            ),
            delta_m_observable_id=_require_optional_text(
                mapping.get("delta_m_observable_id"),
                context="observable_bundle.delta_m_observable_id",
            ),
            delta_m_relation_id=_require_optional_text(
                mapping.get("delta_m_relation_id"),
                context="observable_bundle.delta_m_relation_id",
            ),
            m12_value=(
                None
                if "m12_value" not in mapping
                else ComplexValue.from_dict(
                    _require_mapping(
                        mapping.get("m12_value"),
                        context="observable_bundle.m12_value",
                    )
                )
            ),
            supported_operator_names=_parse_optional_string_ids(
                mapping.get("supported_operator_names"),
                context="observable_bundle.supported_operator_names",
            ),
            unsupported_operator_names=_parse_optional_string_ids(
                mapping.get("unsupported_operator_names"),
                context="observable_bundle.unsupported_operator_names",
            ),
            provenance_ids=(
                ()
                if "provenance_ids" not in mapping
                else _parse_string_ids(
                    mapping.get("provenance_ids"),
                    context="observable_bundle.provenance_ids",
                )
            ),
            notes=_require_optional_text(
                mapping.get("notes"),
                context="observable_bundle.notes",
            ),
        )

    @classmethod
    def from_json(cls, payload: str | bytes) -> "ObservableArtifactBundleV1":
        return cls.from_dict(json.loads(payload))

    @classmethod
    def read_json(cls, path: str | Path) -> "ObservableArtifactBundleV1":
        return cls.from_json(Path(path).read_text(encoding="utf-8"))


@dataclass(frozen=True)
class HadronicArtifactBundleV1:
    """Exported hadronic-input bundle for the kaon NP-only paper scope."""

    metadata: ArtifactMetadata
    operator_basis: str
    renormalization_scheme: str
    matrix_element_formula_id: str
    hamiltonian_convention_id: str
    parity_relation_id: str
    m_K0_GeV: float
    f_K_GeV: float
    B_K_mu_had: float
    q1_matrix_element_GeV4: float
    benchmark_id: str | None = None
    system: str | None = None
    system_id: str | None = None
    scales: tuple[ArtifactScale, ...] = ()
    mu_had_GeV: float | None = None
    source_id: str | None = None
    hadronic_source_id: str | None = None
    supported_operator_names: tuple[str, ...] = DEFAULT_KAON_SUPPORTED_OPERATORS
    unsupported_operator_names: tuple[str, ...] = DEFAULT_KAON_UNSUPPORTED_OPERATORS
    operator_normalization: str | None = None
    hat_B_K_rgi_source_value: float | None = None
    bag_parameter_source_scheme_id: str | None = None
    bag_parameter_transformation_id: str | None = None
    input_provenance_mode_id: str | None = None
    alpha_s_policy_id: str | None = None
    bag_parameter_conversion_formula_id: str | None = None
    bag_parameter_conversion_alpha_s_mu_had: float | None = None
    bag_parameter_conversion_n_f: int | None = None
    bag_parameter_conversion_beta0: float | None = None
    bag_parameter_conversion_exponent: float | None = None
    mass_source_id: str | None = None
    decay_constant_source_id: str | None = None
    bag_parameter_source_id: str | None = None
    source_wilson_bundle_id: str | None = None
    mass_source: ArtifactSourceRecord | None = None
    decay_constant_source: ArtifactSourceRecord | None = None
    bag_parameter_source: ArtifactSourceRecord | None = None
    provenance_ids: tuple[str, ...] = ()
    notes: str | None = None

    def __post_init__(self) -> None:
        if self.metadata.schema_name != HADRONIC_BUNDLE_SCHEMA:
            raise ArtifactSchemaError("hadronic bundle metadata.schema_name is not supported")
        if self.metadata.schema_version != ARTIFACT_SCHEMA_VERSION:
            raise ArtifactSchemaError("hadronic bundle metadata.schema_version is not supported")
        if self.benchmark_id is not None:
            _require_text(self.benchmark_id, context="hadronic.benchmark_id")
        resolved_system = self.system if self.system is not None else self.system_id
        _require_text(resolved_system, context="hadronic.system")
        if self.system is not None and self.system_id is not None and self.system != self.system_id:
            raise ArtifactSchemaError("hadronic.system and hadronic.system_id must match")
        object.__setattr__(self, "system", resolved_system)
        object.__setattr__(self, "system_id", resolved_system)
        resolved_source_id = (
            self.source_id if self.source_id is not None else self.hadronic_source_id
        )
        _require_text(resolved_source_id, context="hadronic.source_id")
        if (
            self.source_id is not None
            and self.hadronic_source_id is not None
            and self.source_id != self.hadronic_source_id
        ):
            raise ArtifactSchemaError(
                "hadronic.source_id and hadronic.hadronic_source_id must match"
            )
        object.__setattr__(self, "source_id", resolved_source_id)
        object.__setattr__(self, "hadronic_source_id", resolved_source_id)
        _require_text(self.operator_basis, context="hadronic.operator_basis")
        _require_text(self.renormalization_scheme, context="hadronic.renormalization_scheme")
        if self.operator_normalization is not None:
            _require_text(
                self.operator_normalization,
                context="hadronic.operator_normalization",
            )
        _require_text(
            self.matrix_element_formula_id,
            context="hadronic.matrix_element_formula_id",
        )
        _require_text(
            self.hamiltonian_convention_id,
            context="hadronic.hamiltonian_convention_id",
        )
        _require_text(self.parity_relation_id, context="hadronic.parity_relation_id")
        resolved_mu_had = self.mu_had_GeV
        if resolved_mu_had is None and self.scales:
            scale_map = {item.name: item.value_gev for item in self.scales}
            resolved_mu_had = scale_map.get("mu_had")
        if resolved_mu_had is None:
            raise ArtifactSchemaError(
                "hadronic.mu_had_GeV must be provided explicitly or via scales"
            )
        resolved_mu_had = _require_float(resolved_mu_had, context="hadronic.mu_had_GeV")
        if resolved_mu_had <= 0.0:
            raise ArtifactSchemaError("hadronic.mu_had_GeV must be positive")
        object.__setattr__(self, "mu_had_GeV", resolved_mu_had)
        if not self.scales:
            object.__setattr__(
                self,
                "scales",
                (
                    ArtifactScale(
                        name="mu_had",
                        role="hadronic evaluation scale",
                        value_gev=resolved_mu_had,
                    ),
                ),
            )
        else:
            scale_map = {item.name: item.value_gev for item in self.scales}
            if "mu_had" in scale_map and not isclose(
                scale_map["mu_had"],
                resolved_mu_had,
                rel_tol=0.0,
                abs_tol=1e-18,
            ):
                raise ArtifactSchemaError("hadronic.mu_had_GeV must match scales['mu_had']")
        if not self.supported_operator_names:
            raise ArtifactSchemaError("hadronic.supported_operator_names must not be empty")
        _ensure_unique_strings([item.name for item in self.scales], context="hadronic.scales")
        _ensure_unique_strings(
            self.supported_operator_names,
            context="hadronic.supported_operator_names",
        )
        _ensure_unique_strings(
            self.unsupported_operator_names,
            context="hadronic.unsupported_operator_names",
        )
        _ensure_unique_strings(self.provenance_ids, context="hadronic.provenance_ids")
        if tuple(self.supported_operator_names) != DEFAULT_KAON_SUPPORTED_OPERATORS:
            raise ArtifactSchemaError(
                "hadronic.supported_operator_names must match the kaon NP-only Q1 subset"
            )
        if tuple(self.unsupported_operator_names) != DEFAULT_KAON_UNSUPPORTED_OPERATORS:
            raise ArtifactSchemaError(
                "hadronic.unsupported_operator_names must match the guarded LR subset"
            )
        for field_name in ("m_K0_GeV", "f_K_GeV", "B_K_mu_had", "q1_matrix_element_GeV4"):
            value = _require_float(getattr(self, field_name), context=f"hadronic.{field_name}")
            if value <= 0.0:
                raise ArtifactSchemaError(f"hadronic.{field_name} must be positive")
        expected_q1_matrix_element = (
            (2.0 / 3.0) * (self.f_K_GeV**2) * (self.m_K0_GeV**2) * self.B_K_mu_had
        )
        if not isclose(
            self.q1_matrix_element_GeV4,
            expected_q1_matrix_element,
            rel_tol=0.0,
            abs_tol=1e-18,
        ):
            raise ArtifactSchemaError(
                "hadronic.q1_matrix_element_GeV4 must match the exported kaon Q1 formula"
            )
        if self.hat_B_K_rgi_source_value is not None:
            value = _require_float(
                self.hat_B_K_rgi_source_value,
                context="hadronic.hat_B_K_rgi_source_value",
            )
            if value <= 0.0:
                raise ArtifactSchemaError(
                    "hadronic.hat_B_K_rgi_source_value must be positive"
                )
        if self.bag_parameter_source_scheme_id is not None:
            _require_text(
                self.bag_parameter_source_scheme_id,
                context="hadronic.bag_parameter_source_scheme_id",
            )
        if self.bag_parameter_transformation_id is not None:
            _require_text(
                self.bag_parameter_transformation_id,
                context="hadronic.bag_parameter_transformation_id",
            )
        if self.input_provenance_mode_id is not None:
            _require_text(
                self.input_provenance_mode_id,
                context="hadronic.input_provenance_mode_id",
            )
        if self.alpha_s_policy_id is not None:
            _require_text(
                self.alpha_s_policy_id,
                context="hadronic.alpha_s_policy_id",
            )
        if self.bag_parameter_conversion_formula_id is not None:
            _require_text(
                self.bag_parameter_conversion_formula_id,
                context="hadronic.bag_parameter_conversion_formula_id",
            )
        if self.bag_parameter_conversion_alpha_s_mu_had is not None:
            alpha_s_value = _require_float(
                self.bag_parameter_conversion_alpha_s_mu_had,
                context="hadronic.bag_parameter_conversion_alpha_s_mu_had",
            )
            if alpha_s_value <= 0.0:
                raise ArtifactSchemaError(
                    "hadronic.bag_parameter_conversion_alpha_s_mu_had must be positive"
                )
        if self.bag_parameter_conversion_n_f is not None:
            n_f_value = _require_int(
                self.bag_parameter_conversion_n_f,
                context="hadronic.bag_parameter_conversion_n_f",
            )
            if n_f_value <= 0:
                raise ArtifactSchemaError(
                    "hadronic.bag_parameter_conversion_n_f must be positive"
                )
        if self.bag_parameter_conversion_beta0 is not None:
            beta0_value = _require_float(
                self.bag_parameter_conversion_beta0,
                context="hadronic.bag_parameter_conversion_beta0",
            )
            if beta0_value <= 0.0:
                raise ArtifactSchemaError(
                    "hadronic.bag_parameter_conversion_beta0 must be positive"
                )
        if self.bag_parameter_conversion_exponent is not None:
            _require_float(
                self.bag_parameter_conversion_exponent,
                context="hadronic.bag_parameter_conversion_exponent",
            )
        if self.source_wilson_bundle_id is not None:
            _require_text(
                self.source_wilson_bundle_id,
                context="hadronic.source_wilson_bundle_id",
            )
        for field_name in ("mass_source", "decay_constant_source", "bag_parameter_source"):
            source_value = getattr(self, field_name)
            if source_value is not None and not isinstance(source_value, ArtifactSourceRecord):
                raise ArtifactSchemaError(
                    f"hadronic.{field_name} must be an ArtifactSourceRecord when provided"
                )
        if self.mass_source_id is not None:
            _require_text(self.mass_source_id, context="hadronic.mass_source_id")
        if self.decay_constant_source_id is not None:
            _require_text(
                self.decay_constant_source_id,
                context="hadronic.decay_constant_source_id",
            )
        if self.bag_parameter_source_id is not None:
            _require_text(
                self.bag_parameter_source_id,
                context="hadronic.bag_parameter_source_id",
            )
        resolved_mass_source_id = (
            self.mass_source.source_id if self.mass_source is not None else self.mass_source_id
        )
        resolved_decay_source_id = (
            self.decay_constant_source.source_id
            if self.decay_constant_source is not None
            else self.decay_constant_source_id
        )
        resolved_bag_source_id = (
            self.bag_parameter_source.source_id
            if self.bag_parameter_source is not None
            else self.bag_parameter_source_id
        )
        if resolved_mass_source_id is not None:
            object.__setattr__(self, "mass_source_id", resolved_mass_source_id)
        if resolved_decay_source_id is not None:
            object.__setattr__(self, "decay_constant_source_id", resolved_decay_source_id)
        if resolved_bag_source_id is not None:
            object.__setattr__(self, "bag_parameter_source_id", resolved_bag_source_id)
        if (
            self.bag_parameter_source is not None
            and self.bag_parameter_source_scheme_id is not None
            and self.bag_parameter_source.renormalization_scheme_id
            != self.bag_parameter_source_scheme_id
        ):
            raise ArtifactSchemaError(
                "hadronic.bag_parameter_source.renormalization_scheme_id must match "
                "hadronic.bag_parameter_source_scheme_id"
            )
        if (
            self.bag_parameter_source is not None
            and self.bag_parameter_transformation_id is not None
            and self.bag_parameter_source.transformation_id
            != self.bag_parameter_transformation_id
        ):
            raise ArtifactSchemaError(
                "hadronic.bag_parameter_source.transformation_id must match "
                "hadronic.bag_parameter_transformation_id"
            )
        conversion_fields_present = any(
            value is not None
            for value in (
                self.bag_parameter_conversion_formula_id,
                self.bag_parameter_conversion_alpha_s_mu_had,
                self.bag_parameter_conversion_n_f,
                self.bag_parameter_conversion_beta0,
                self.bag_parameter_conversion_exponent,
            )
        )
        if conversion_fields_present:
            missing_fields = [
                field_name
                for field_name in (
                    "bag_parameter_conversion_formula_id",
                    "bag_parameter_conversion_alpha_s_mu_had",
                    "bag_parameter_conversion_n_f",
                    "bag_parameter_conversion_beta0",
                    "bag_parameter_conversion_exponent",
                    "hat_B_K_rgi_source_value",
                )
                if getattr(self, field_name) is None
            ]
            if missing_fields:
                missing = ", ".join(missing_fields)
                raise ArtifactSchemaError(
                    "hadronic bag-parameter conversion claim is incomplete; "
                    f"missing {missing}"
                )
            expected_beta0 = _qcd_beta0_from_n_f(self.bag_parameter_conversion_n_f)
            if not isclose(
                self.bag_parameter_conversion_beta0,
                expected_beta0,
                rel_tol=0.0,
                abs_tol=1e-15,
            ):
                raise ArtifactSchemaError(
                    "hadronic.bag_parameter_conversion_beta0 must match 11 - 2*n_f/3"
                )
            expected_b_k_from_hat = self.hat_B_K_rgi_source_value * (
                self.bag_parameter_conversion_alpha_s_mu_had
                ** self.bag_parameter_conversion_exponent
            )
            if not isclose(
                self.B_K_mu_had,
                expected_b_k_from_hat,
                rel_tol=0.0,
                abs_tol=1e-15,
            ):
                raise ArtifactSchemaError(
                    "hadronic B_K_mu_had must match the declared hat_B_K conversion audit data"
                )
        if self.notes is not None:
            _require_text(self.notes, context="hadronic.notes")

    def to_dict(self) -> JSONDict:
        payload: JSONDict = {
            "metadata": self.metadata.to_dict(),
            "source_id": self.source_id,
            "hadronic_source_id": self.hadronic_source_id,
            "system": self.system,
            "system_id": self.system_id,
            "operator_basis": self.operator_basis,
            "renormalization_scheme": self.renormalization_scheme,
            "scales": _scales_to_dict(self.scales),
            "matrix_element_formula_id": self.matrix_element_formula_id,
            "hamiltonian_convention_id": self.hamiltonian_convention_id,
            "parity_relation_id": self.parity_relation_id,
            "supported_operator_names": list(self.supported_operator_names),
            "unsupported_operator_names": list(self.unsupported_operator_names),
            "mu_had_GeV": self.mu_had_GeV,
            "m_K0_GeV": self.m_K0_GeV,
            "f_K_GeV": self.f_K_GeV,
            "B_K_mu_had": self.B_K_mu_had,
            "q1_matrix_element_GeV4": self.q1_matrix_element_GeV4,
            "provenance_ids": _ids_to_list(self.provenance_ids),
        }
        if self.benchmark_id is not None:
            payload["benchmark_id"] = self.benchmark_id
        if self.operator_normalization is not None:
            payload["operator_normalization"] = self.operator_normalization
        if self.hat_B_K_rgi_source_value is not None:
            payload["hat_B_K_rgi_source_value"] = self.hat_B_K_rgi_source_value
        if self.bag_parameter_source_scheme_id is not None:
            payload["bag_parameter_source_scheme_id"] = self.bag_parameter_source_scheme_id
        if self.bag_parameter_transformation_id is not None:
            payload["bag_parameter_transformation_id"] = self.bag_parameter_transformation_id
        if self.input_provenance_mode_id is not None:
            payload["input_provenance_mode_id"] = self.input_provenance_mode_id
        if self.alpha_s_policy_id is not None:
            payload["alpha_s_policy_id"] = self.alpha_s_policy_id
        if self.bag_parameter_conversion_formula_id is not None:
            payload["bag_parameter_conversion_formula_id"] = (
                self.bag_parameter_conversion_formula_id
            )
        if self.bag_parameter_conversion_alpha_s_mu_had is not None:
            payload["bag_parameter_conversion_alpha_s_mu_had"] = (
                self.bag_parameter_conversion_alpha_s_mu_had
            )
        if self.bag_parameter_conversion_n_f is not None:
            payload["bag_parameter_conversion_n_f"] = self.bag_parameter_conversion_n_f
        if self.bag_parameter_conversion_beta0 is not None:
            payload["bag_parameter_conversion_beta0"] = self.bag_parameter_conversion_beta0
        if self.bag_parameter_conversion_exponent is not None:
            payload["bag_parameter_conversion_exponent"] = (
                self.bag_parameter_conversion_exponent
            )
        if self.mass_source_id is not None:
            payload["mass_source_id"] = self.mass_source_id
        if self.decay_constant_source_id is not None:
            payload["decay_constant_source_id"] = self.decay_constant_source_id
        if self.bag_parameter_source_id is not None:
            payload["bag_parameter_source_id"] = self.bag_parameter_source_id
        if self.source_wilson_bundle_id is not None:
            payload["source_wilson_bundle_id"] = self.source_wilson_bundle_id
        if self.mass_source is not None:
            payload["mass_source"] = self.mass_source.to_dict()
        if self.decay_constant_source is not None:
            payload["decay_constant_source"] = self.decay_constant_source.to_dict()
        if self.bag_parameter_source is not None:
            payload["bag_parameter_source"] = self.bag_parameter_source.to_dict()
        if self.notes is not None:
            payload["notes"] = self.notes
        return payload

    def to_json(self) -> str:
        return _dump_json(self.to_dict())

    def write_json(self, path: str | Path) -> Path:
        destination = Path(path)
        destination.write_text(self.to_json(), encoding="utf-8")
        return destination

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "HadronicArtifactBundleV1":
        mapping = _require_mapping(data, context="hadronic_bundle")
        return cls(
            metadata=ArtifactMetadata.from_dict(
                _require_mapping(mapping.get("metadata"), context="hadronic_bundle.metadata")
            ),
            operator_basis=_require_text(
                mapping.get("operator_basis"), context="hadronic_bundle.operator_basis"
            ),
            renormalization_scheme=_require_text(
                mapping.get("renormalization_scheme"),
                context="hadronic_bundle.renormalization_scheme",
            ),
            scales=_parse_scales(
                mapping.get("scales", ()),
                context="hadronic_bundle.scales",
            ),
            matrix_element_formula_id=_require_text(
                mapping.get("matrix_element_formula_id"),
                context="hadronic_bundle.matrix_element_formula_id",
            ),
            hamiltonian_convention_id=_require_text(
                mapping.get("hamiltonian_convention_id"),
                context="hadronic_bundle.hamiltonian_convention_id",
            ),
            parity_relation_id=_require_text(
                mapping.get("parity_relation_id"),
                context="hadronic_bundle.parity_relation_id",
            ),
            supported_operator_names=(
                DEFAULT_KAON_SUPPORTED_OPERATORS
                if "supported_operator_names" not in mapping
                else _parse_string_ids(
                    mapping.get("supported_operator_names"),
                    context="hadronic_bundle.supported_operator_names",
                )
            ),
            unsupported_operator_names=(
                DEFAULT_KAON_UNSUPPORTED_OPERATORS
                if "unsupported_operator_names" not in mapping
                else _parse_string_ids(
                    mapping.get("unsupported_operator_names"),
                    context="hadronic_bundle.unsupported_operator_names",
                )
            ),
            m_K0_GeV=_require_float(mapping.get("m_K0_GeV"), context="hadronic_bundle.m_K0_GeV"),
            f_K_GeV=_require_float(mapping.get("f_K_GeV"), context="hadronic_bundle.f_K_GeV"),
            B_K_mu_had=_require_float(
                mapping.get("B_K_mu_had"),
                context="hadronic_bundle.B_K_mu_had",
            ),
            q1_matrix_element_GeV4=_require_float(
                mapping.get("q1_matrix_element_GeV4"),
                context="hadronic_bundle.q1_matrix_element_GeV4",
            ),
            benchmark_id=_require_optional_text(
                mapping.get("benchmark_id"),
                context="hadronic_bundle.benchmark_id",
            ),
            system=_require_optional_text(
                mapping.get("system"),
                context="hadronic_bundle.system",
            ),
            system_id=_require_optional_text(
                mapping.get("system_id"),
                context="hadronic_bundle.system_id",
            ),
            mu_had_GeV=_require_optional_float(
                mapping.get("mu_had_GeV"),
                context="hadronic_bundle.mu_had_GeV",
            ),
            source_id=_require_optional_text(
                mapping.get("source_id"),
                context="hadronic_bundle.source_id",
            ),
            hadronic_source_id=_require_optional_text(
                mapping.get("hadronic_source_id"),
                context="hadronic_bundle.hadronic_source_id",
            ),
            operator_normalization=_require_optional_text(
                mapping.get("operator_normalization"),
                context="hadronic_bundle.operator_normalization",
            ),
            hat_B_K_rgi_source_value=_require_optional_float(
                mapping.get("hat_B_K_rgi_source_value"),
                context="hadronic_bundle.hat_B_K_rgi_source_value",
            ),
            bag_parameter_source_scheme_id=_require_optional_text(
                mapping.get("bag_parameter_source_scheme_id"),
                context="hadronic_bundle.bag_parameter_source_scheme_id",
            ),
            bag_parameter_transformation_id=_require_optional_text(
                mapping.get("bag_parameter_transformation_id"),
                context="hadronic_bundle.bag_parameter_transformation_id",
            ),
            input_provenance_mode_id=_require_optional_text(
                mapping.get("input_provenance_mode_id"),
                context="hadronic_bundle.input_provenance_mode_id",
            ),
            alpha_s_policy_id=_require_optional_text(
                mapping.get("alpha_s_policy_id"),
                context="hadronic_bundle.alpha_s_policy_id",
            ),
            bag_parameter_conversion_formula_id=_require_optional_text(
                mapping.get("bag_parameter_conversion_formula_id"),
                context="hadronic_bundle.bag_parameter_conversion_formula_id",
            ),
            bag_parameter_conversion_alpha_s_mu_had=_require_optional_float(
                mapping.get("bag_parameter_conversion_alpha_s_mu_had"),
                context="hadronic_bundle.bag_parameter_conversion_alpha_s_mu_had",
            ),
            bag_parameter_conversion_n_f=(
                None
                if "bag_parameter_conversion_n_f" not in mapping
                else _require_int(
                    mapping.get("bag_parameter_conversion_n_f"),
                    context="hadronic_bundle.bag_parameter_conversion_n_f",
                )
            ),
            bag_parameter_conversion_beta0=_require_optional_float(
                mapping.get("bag_parameter_conversion_beta0"),
                context="hadronic_bundle.bag_parameter_conversion_beta0",
            ),
            bag_parameter_conversion_exponent=_require_optional_float(
                mapping.get("bag_parameter_conversion_exponent"),
                context="hadronic_bundle.bag_parameter_conversion_exponent",
            ),
            mass_source_id=_require_optional_text(
                mapping.get("mass_source_id"),
                context="hadronic_bundle.mass_source_id",
            ),
            decay_constant_source_id=_require_optional_text(
                mapping.get("decay_constant_source_id"),
                context="hadronic_bundle.decay_constant_source_id",
            ),
            bag_parameter_source_id=_require_optional_text(
                mapping.get("bag_parameter_source_id"),
                context="hadronic_bundle.bag_parameter_source_id",
            ),
            source_wilson_bundle_id=_require_optional_text(
                mapping.get("source_wilson_bundle_id"),
                context="hadronic_bundle.source_wilson_bundle_id",
            ),
            mass_source=(
                None
                if "mass_source" not in mapping
                else ArtifactSourceRecord.from_dict(
                    _require_mapping(
                        mapping.get("mass_source"),
                        context="hadronic_bundle.mass_source",
                    )
                )
            ),
            decay_constant_source=(
                None
                if "decay_constant_source" not in mapping
                else ArtifactSourceRecord.from_dict(
                    _require_mapping(
                        mapping.get("decay_constant_source"),
                        context="hadronic_bundle.decay_constant_source",
                    )
                )
            ),
            bag_parameter_source=(
                None
                if "bag_parameter_source" not in mapping
                else ArtifactSourceRecord.from_dict(
                    _require_mapping(
                        mapping.get("bag_parameter_source"),
                        context="hadronic_bundle.bag_parameter_source",
                    )
                )
            ),
            provenance_ids=(
                ()
                if "provenance_ids" not in mapping
                else _parse_string_ids(
                    mapping.get("provenance_ids"),
                    context="hadronic_bundle.provenance_ids",
                )
            ),
            notes=_require_optional_text(mapping.get("notes"), context="hadronic_bundle.notes"),
        )

    @classmethod
    def from_json(cls, payload: str | bytes) -> "HadronicArtifactBundleV1":
        return cls.from_dict(json.loads(payload))

    @classmethod
    def read_json(cls, path: str | Path) -> "HadronicArtifactBundleV1":
        return cls.from_json(Path(path).read_text(encoding="utf-8"))


@dataclass(frozen=True)
class ProvenanceBundleV1:
    """Frozen provenance manifest for paper exports."""

    metadata: ArtifactMetadata
    records: tuple[ProvenanceRecord, ...]

    def __post_init__(self) -> None:
        if self.metadata.schema_name != PROVENANCE_BUNDLE_SCHEMA:
            raise ArtifactSchemaError("provenance bundle metadata.schema_name is not supported")
        if self.metadata.schema_version != ARTIFACT_SCHEMA_VERSION:
            raise ArtifactSchemaError("provenance bundle metadata.schema_version is not supported")
        if not self.records:
            raise ArtifactSchemaError("provenance.records must not be empty")
        _ensure_unique_strings(
            [item.record_id for item in self.records],
            context="provenance.records",
        )

    def to_dict(self) -> JSONDict:
        return {
            "metadata": self.metadata.to_dict(),
            "records": [item.to_dict() for item in self.records],
        }

    def to_json(self) -> str:
        return _dump_json(self.to_dict())

    def write_json(self, path: str | Path) -> Path:
        destination = Path(path)
        destination.write_text(self.to_json(), encoding="utf-8")
        return destination

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "ProvenanceBundleV1":
        mapping = _require_mapping(data, context="provenance_bundle")
        record_entries = _require_sequence(
            mapping.get("records"), context="provenance_bundle.records"
        )
        return cls(
            metadata=ArtifactMetadata.from_dict(
                _require_mapping(mapping.get("metadata"), context="provenance_bundle.metadata")
            ),
            records=tuple(
                ProvenanceRecord.from_dict(
                    _require_mapping(item, context="provenance_bundle.records")
                )
                for item in record_entries
            ),
        )

    @classmethod
    def from_json(cls, payload: str | bytes) -> "ProvenanceBundleV1":
        return cls.from_dict(json.loads(payload))

    @classmethod
    def read_json(cls, path: str | Path) -> "ProvenanceBundleV1":
        return cls.from_json(Path(path).read_text(encoding="utf-8"))


PaperArtifactBundle: TypeAlias = (
    WilsonArtifactBundleV1
    | HadronicArtifactBundleV1
    | ObservableArtifactBundleV1
    | ProvenanceBundleV1
)


def artifact_from_dict(data: Mapping[str, Any]) -> PaperArtifactBundle:
    """Deserialize one supported paper artifact from a JSON dictionary."""

    mapping = _require_mapping(data, context="artifact")
    metadata = ArtifactMetadata.from_dict(
        _require_mapping(mapping.get("metadata"), context="artifact.metadata")
    )
    if metadata.schema_name == WILSON_BUNDLE_SCHEMA:
        return WilsonArtifactBundleV1.from_dict(mapping)
    if metadata.schema_name == HADRONIC_BUNDLE_SCHEMA:
        return HadronicArtifactBundleV1.from_dict(mapping)
    if metadata.schema_name == OBSERVABLE_BUNDLE_SCHEMA:
        return ObservableArtifactBundleV1.from_dict(mapping)
    if metadata.schema_name == PROVENANCE_BUNDLE_SCHEMA:
        return ProvenanceBundleV1.from_dict(mapping)
    raise ArtifactSchemaError(f"unsupported artifact schema: {metadata.schema_name}")


def read_artifact(path: str | Path) -> PaperArtifactBundle:
    """Read one supported paper artifact from disk."""

    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    return artifact_from_dict(payload)


@dataclass(frozen=True)
class Paper07101869ArtifactExportSet:
    """One deterministic export set for the default paper-facing kaon benchmark."""

    wilson_bundle: WilsonArtifactBundleV1
    hadronic_bundle: HadronicArtifactBundleV1
    observable_bundle: ObservableArtifactBundleV1
    provenance_bundle: ProvenanceBundleV1

    def __post_init__(self) -> None:
        point_ids = {
            self.wilson_bundle.metadata.point_id,
            self.hadronic_bundle.metadata.point_id,
            self.observable_bundle.metadata.point_id,
            self.provenance_bundle.metadata.point_id,
        }
        if len(point_ids) != 1:
            raise ArtifactSchemaError("artifact export set point_id values must match")
        if self.observable_bundle.source_wilson_bundle_id != self.wilson_bundle.metadata.bundle_id:
            raise ArtifactSchemaError(
                "observable bundle must reference the exported Wilson bundle"
            )
        if (
            self.observable_bundle.source_hadronic_bundle_id
            != self.hadronic_bundle.metadata.bundle_id
        ):
            raise ArtifactSchemaError(
                "observable bundle must reference the exported hadronic bundle"
            )
        if (
            self.wilson_bundle.supported_operator_names
            and tuple(self.wilson_bundle.supported_operator_names)
            != tuple(self.hadronic_bundle.supported_operator_names)
        ):
            raise ArtifactSchemaError(
                "Wilson and hadronic bundles must advertise the same supported operator subset"
            )
        if (
            self.wilson_bundle.unsupported_operator_names
            and tuple(self.wilson_bundle.unsupported_operator_names)
            != tuple(self.hadronic_bundle.unsupported_operator_names)
        ):
            raise ArtifactSchemaError(
                "Wilson and hadronic bundles must advertise the same unsupported operator subset"
            )
        if (
            self.observable_bundle.supported_operator_names
            and tuple(self.observable_bundle.supported_operator_names)
            != tuple(self.hadronic_bundle.supported_operator_names)
        ):
            raise ArtifactSchemaError(
                "Observable and hadronic bundles must advertise the same supported operator subset"
            )
        if (
            self.observable_bundle.unsupported_operator_names
            and tuple(self.observable_bundle.unsupported_operator_names)
            != tuple(self.hadronic_bundle.unsupported_operator_names)
        ):
            raise ArtifactSchemaError(
                "Observable and hadronic bundles must advertise the same "
                "unsupported operator subset"
            )
        provenance_ids = {record.record_id for record in self.provenance_bundle.records}
        referenced_ids = (
            set(self.wilson_bundle.provenance_ids)
            | set(self.hadronic_bundle.provenance_ids)
            | set(self.observable_bundle.provenance_ids)
        )
        missing_ids = referenced_ids - provenance_ids
        if missing_ids:
            missing = ", ".join(sorted(missing_ids))
            raise ArtifactSchemaError(f"export set is missing provenance records: {missing}")


@dataclass(frozen=True)
class Paper07101869ArtifactExportPaths:
    """Filesystem paths for one deterministic paper artifact export."""

    root_dir: Path
    wilson_path: Path
    hadronic_path: Path
    observable_path: Path
    provenance_path: Path

    def as_dict(self) -> JSONDict:
        return {
            "root_dir": str(self.root_dir),
            "wilson_path": str(self.wilson_path),
            "hadronic_path": str(self.hadronic_path),
            "observable_path": str(self.observable_path),
            "provenance_path": str(self.provenance_path),
        }


def _artifact_source_from_mapping(source_mapping: Mapping[str, Any]) -> ArtifactSourceRecord:
    return ArtifactSourceRecord.from_dict(source_mapping)


def _artifact_source_with_id(
    source: ArtifactSourceRecord,
    *,
    source_id: str,
) -> ArtifactSourceRecord:
    return replace(source, source_id=source_id)


def _build_default_wilson_scales(
    *,
    matching_scale_gev: float,
    hadronic_scale_gev: float,
    propagator_mass_gev: float,
) -> tuple[ArtifactScale, ...]:
    return (
        ArtifactScale(name="mu_match", role="matching scale", value_gev=matching_scale_gev),
        ArtifactScale(
            name="mu_had",
            role="hadronic evaluation scale",
            value_gev=hadronic_scale_gev,
        ),
        ArtifactScale(name="m_g1", role="KK-gluon propagator mass", value_gev=propagator_mass_gev),
    )


def _supported_wilson_records(
    *,
    coefficient_map: Mapping[str, complex],
    sector: str,
    system: str,
) -> tuple[WilsonCoefficientRecord, ...]:
    return tuple(
        WilsonCoefficientRecord(
            sector=sector,
            system=system,
            operator=operator_name,
            value=ComplexValue.from_complex(coefficient_map[operator_name]),
        )
        for operator_name in DEFAULT_KAON_SUPPORTED_OPERATORS
    )


def _provenance_record_from_source(source: ArtifactSourceRecord) -> ProvenanceRecord:
    return ProvenanceRecord(
        record_id=source.source_id,
        category=source.source_kind,
        label=source.locator_label,
        version=str(source.year),
        source=source.citation,
        citation=source.notes if source.notes is not None else source.citation,
    )


def _default_bag_parameter_conversion_audit(
    *,
    mu_had_gev: float,
) -> JSONDict:
    from qcd import alpha_s
    from qcd.beta_function import beta_0

    from .eft_deltaf2.rg_inputs import default_paper_0710_1869_rg_thresholds

    thresholds = default_paper_0710_1869_rg_thresholds()
    alpha_s_mu_had = float(
        alpha_s(
            mu_had_gev,
            n_loops=1,
            matching_loops=0,
            thresholds=list(thresholds),
        )
    )
    n_f = _n_f_at_scale(mu_had_gev, thresholds)
    beta0_value = float(beta_0(n_f))
    exponent = 4.0 / (2.0 * beta0_value)
    return {
        "formula_id": DEFAULT_BAG_PARAMETER_CONVERSION_FORMULA_ID,
        "alpha_s_mu_had": alpha_s_mu_had,
        "n_f": n_f,
        "beta0": beta0_value,
        "exponent": exponent,
    }


def _require_frozen_default_kaon_real(
    *,
    context: str,
    actual: float,
    expected: float,
    abs_tol: float,
) -> float:
    if not isclose(actual, expected, rel_tol=0.0, abs_tol=abs_tol):
        raise ArtifactSchemaError(
            f"{context} drifted beyond the frozen default-kaon export tolerance"
        )
    return expected


def _require_frozen_default_kaon_complex(
    *,
    context: str,
    actual: complex,
    expected: complex,
    abs_tol: float,
) -> complex:
    _require_frozen_default_kaon_real(
        context=f"{context}.real",
        actual=float(actual.real),
        expected=float(expected.real),
        abs_tol=abs_tol,
    )
    _require_frozen_default_kaon_real(
        context=f"{context}.imag",
        actual=float(actual.imag),
        expected=float(expected.imag),
        abs_tol=abs_tol,
    )
    return expected


def build_default_paper_0710_1869_kaon_artifact_export_set() -> Paper07101869ArtifactExportSet:
    """Build the canonical four-file artifact set for the default kaon benchmark."""

    from .eft_deltaf2.hadronic import default_paper_0710_1869_kaon_hadronic_bundle
    from .eft_deltaf2.matching_kkgluon import default_paper_0710_1869_kaon_matching
    from .eft_deltaf2.observables import compute_kaon_np_observables
    from .eft_deltaf2.rg import run_deltaf2_wilsons_lo

    matching = default_paper_0710_1869_kaon_matching()
    rg_result = run_deltaf2_wilsons_lo(matching.wilsons)
    hadronic = default_paper_0710_1869_kaon_hadronic_bundle()
    observable = compute_kaon_np_observables(rg_result.wilsons, hadronic_bundle=hadronic)
    bag_parameter_conversion = _default_bag_parameter_conversion_audit(
        mu_had_gev=hadronic.mu_had_GeV
    )
    frozen_matching_q1_vll = _require_frozen_default_kaon_complex(
        context="default kaon matching Q1_VLL",
        actual=matching.wilsons.q1_vll,
        expected=DEFAULT_KAON_FROZEN_MATCHING_Q1_VLL,
        abs_tol=1e-24,
    )
    frozen_rg_q1_vll = _require_frozen_default_kaon_complex(
        context="default kaon RG Q1_VLL",
        actual=rg_result.wilsons.q1_vll,
        expected=DEFAULT_KAON_FROZEN_RG_Q1_VLL,
        abs_tol=1e-24,
    )
    frozen_m12 = _require_frozen_default_kaon_complex(
        context="default kaon M12_K_NP",
        actual=observable.M12_K_NP_GeV,
        expected=DEFAULT_KAON_FROZEN_M12_K_NP,
        abs_tol=1e-24,
    )
    frozen_delta_m = _require_frozen_default_kaon_real(
        context="default kaon Delta_m_K_NP",
        actual=observable.delta_m_K_NP_GeV,
        expected=DEFAULT_KAON_FROZEN_DELTA_M_K_NP_GEV,
        abs_tol=1e-24,
    )
    frozen_bag_alpha_s = _require_frozen_default_kaon_real(
        context="default kaon bag-parameter conversion alpha_s(mu_had)",
        actual=float(bag_parameter_conversion["alpha_s_mu_had"]),
        expected=DEFAULT_KAON_FROZEN_BAG_PARAMETER_CONVERSION_ALPHA_S_MU_HAD,
        abs_tol=1e-15,
    )
    matching_coefficients = dict(matching.wilsons.coefficients)
    matching_coefficients["Q1_VLL"] = frozen_matching_q1_vll
    rg_coefficients = dict(rg_result.coefficients)
    rg_coefficients["Q1_VLL"] = frozen_rg_q1_vll
    hadronic_scale_gev = hadronic.mu_had_GeV
    if not isclose(
        rg_result.wilsons.matching_scale_GeV,
        hadronic_scale_gev,
        rel_tol=0.0,
        abs_tol=1e-18,
    ):
        raise ArtifactSchemaError(
            "RG Wilson evaluation scale must match hadronic.mu_had_GeV for the "
            "canonical paper artifact export"
        )

    shared_scales = _build_default_wilson_scales(
        matching_scale_gev=matching.matching_scale_GeV,
        hadronic_scale_gev=hadronic_scale_gev,
        propagator_mass_gev=matching.propagator_mass_GeV,
    )

    wilson_bundle = WilsonArtifactBundleV1(
        metadata=ArtifactMetadata.create(
            schema_name=WILSON_BUNDLE_SCHEMA,
            bundle_id=DEFAULT_KAON_WILSON_BUNDLE_ID,
            point_id=DEFAULT_KAON_ARTIFACT_POINT_ID,
        ),
        operator_basis=rg_result.wilsons.operator_basis_id,
        renormalization_scheme=rg_result.wilsons.renormalization_scheme_id,
        scales=shared_scales,
        coefficients=_supported_wilson_records(
            coefficient_map=rg_coefficients,
            sector=rg_result.wilsons.sector_id,
            system=rg_result.wilsons.system_id,
        ),
        operator_normalization=rg_result.wilsons.operator_normalization_id,
        coefficient_scale_name="mu_had",
        matching_scale_name="mu_match",
        matching_coefficients=_supported_wilson_records(
            coefficient_map=matching_coefficients,
            sector=matching.wilsons.sector_id,
            system=matching.wilsons.system_id,
        ),
        supported_operator_names=DEFAULT_KAON_SUPPORTED_OPERATORS,
        unsupported_operator_names=DEFAULT_KAON_UNSUPPORTED_OPERATORS,
        provenance_ids=(
            DEFAULT_PAPER_PROVENANCE_ID,
            DEFAULT_MATCHING_PROVENANCE_ID,
            DEFAULT_RG_PROVENANCE_ID,
        ),
        notes=(
            "Kaon NP-only Wilson surface at mu_had. Exported coefficients are restricted to "
            "Q1_VLL and Q1_VRR only; LR operators remain unsupported and are carried only as "
            "explicit out-of-scope labels."
        ),
    )

    mass_source = _artifact_source_from_mapping(hadronic.mass_source.as_dict())
    decay_constant_source = _artifact_source_from_mapping(hadronic.decay_constant_source.as_dict())
    bag_parameter_source = _artifact_source_from_mapping(hadronic.bag_parameter_source.as_dict())
    hadronic_provenance_ids = (
        hadronic.source_id,
        mass_source.source_id,
        decay_constant_source.source_id,
        bag_parameter_source.source_id,
    )
    hadronic_bundle = HadronicArtifactBundleV1(
        metadata=ArtifactMetadata.create(
            schema_name=HADRONIC_BUNDLE_SCHEMA,
            bundle_id=hadronic.bundle_id,
            point_id=DEFAULT_KAON_ARTIFACT_POINT_ID,
        ),
        benchmark_id="kaon.np_only.default.v1",
        source_id=hadronic.source_id,
        hadronic_source_id=hadronic.source_id,
        system=hadronic.system_id,
        system_id=hadronic.system_id,
        operator_basis=hadronic.operator_basis_id,
        renormalization_scheme=hadronic.renormalization_scheme_id,
        scales=shared_scales,
        mu_had_GeV=hadronic.mu_had_GeV,
        matrix_element_formula_id=hadronic.matrix_element_formula_id,
        hamiltonian_convention_id=hadronic.hamiltonian_convention_id,
        parity_relation_id=hadronic.parity_relation_id,
        supported_operator_names=hadronic.supported_operator_names,
        unsupported_operator_names=hadronic.unsupported_operator_names,
        m_K0_GeV=hadronic.m_K0_GeV,
        f_K_GeV=hadronic.f_K_GeV,
        B_K_mu_had=hadronic.B_K_mu_had,
        q1_matrix_element_GeV4=hadronic.q1_matrix_element_GeV4,
        operator_normalization=hadronic.operator_normalization_id,
        hat_B_K_rgi_source_value=hadronic.hat_B_K_rgi_source_value,
        bag_parameter_source_scheme_id=hadronic.bag_parameter_source_scheme_id,
        bag_parameter_transformation_id=hadronic.bag_parameter_transformation_id,
        input_provenance_mode_id=hadronic.input_provenance_mode_id,
        alpha_s_policy_id=hadronic.alpha_s_policy_id,
        bag_parameter_conversion_formula_id=bag_parameter_conversion["formula_id"],
        bag_parameter_conversion_alpha_s_mu_had=frozen_bag_alpha_s,
        bag_parameter_conversion_n_f=bag_parameter_conversion["n_f"],
        bag_parameter_conversion_beta0=bag_parameter_conversion["beta0"],
        bag_parameter_conversion_exponent=bag_parameter_conversion["exponent"],
        mass_source_id=mass_source.source_id,
        decay_constant_source_id=decay_constant_source.source_id,
        bag_parameter_source_id=bag_parameter_source.source_id,
        source_wilson_bundle_id=wilson_bundle.metadata.bundle_id,
        mass_source=mass_source,
        decay_constant_source=decay_constant_source,
        bag_parameter_source=bag_parameter_source,
        provenance_ids=hadronic_provenance_ids,
        notes=hadronic.notes,
    )

    observable_bundle = ObservableArtifactBundleV1(
        metadata=ArtifactMetadata.create(
            schema_name=OBSERVABLE_BUNDLE_SCHEMA,
            bundle_id=DEFAULT_KAON_OBSERVABLE_BUNDLE_ID,
            point_id=DEFAULT_KAON_ARTIFACT_POINT_ID,
        ),
        source_wilson_bundle_id=wilson_bundle.metadata.bundle_id,
        interpretation=observable.interpretation,
        operator_basis=observable.operator_basis_id,
        renormalization_scheme=observable.renormalization_scheme_id,
        scales=shared_scales,
        observables=(
            ObservableRecord(
                name=f"{observable.m12_observable_id}.re",
                system=observable.system_id,
                value=float(frozen_m12.real),
                units="GeV",
            ),
            ObservableRecord(
                name=f"{observable.m12_observable_id}.im",
                system=observable.system_id,
                value=float(frozen_m12.imag),
                units="GeV",
            ),
            ObservableRecord(
                name=observable.delta_m_observable_id,
                system=observable.system_id,
                value=frozen_delta_m,
                units="GeV",
            ),
        ),
        operator_normalization=observable.operator_normalization_id,
        source_hadronic_bundle_id=hadronic_bundle.metadata.bundle_id,
        observable_scope_id=observable.observable_scope_id,
        m12_observable_id=observable.m12_observable_id,
        delta_m_observable_id=observable.delta_m_observable_id,
        delta_m_relation_id=DELTA_M_RELATION_ID,
        m12_value=ComplexValue.from_complex(frozen_m12),
        supported_operator_names=DEFAULT_KAON_SUPPORTED_OPERATORS,
        unsupported_operator_names=DEFAULT_KAON_UNSUPPORTED_OPERATORS,
        provenance_ids=(
            DEFAULT_PAPER_PROVENANCE_ID,
            DEFAULT_MATCHING_PROVENANCE_ID,
            DEFAULT_RG_PROVENANCE_ID,
            *hadronic_provenance_ids,
        ),
        notes=observable.notes,
    )

    provenance_bundle = ProvenanceBundleV1(
        metadata=ArtifactMetadata.create(
            schema_name=PROVENANCE_BUNDLE_SCHEMA,
            bundle_id=DEFAULT_KAON_PROVENANCE_BUNDLE_ID,
            point_id=DEFAULT_KAON_ARTIFACT_POINT_ID,
        ),
        records=(
            ProvenanceRecord(
                record_id=DEFAULT_PAPER_PROVENANCE_ID,
                category="paper",
                label="0710.1869 paper-mode scope",
                version="arXiv:0710.1869",
                source="doi:10.48550/arXiv.0710.1869",
                citation="arXiv:0710.1869",
            ),
            ProvenanceRecord(
                record_id=DEFAULT_MATCHING_PROVENANCE_ID,
                category="code",
                label="Default kaon KK-gluon matching export",
                version="paper-mode:v1",
                source="quarkConstraints.paper_0710_1869.eft_deltaf2.matching_kkgluon",
                citation="arXiv:0710.1869 matching path",
            ),
            ProvenanceRecord(
                record_id=DEFAULT_RG_PROVENANCE_ID,
                category="code",
                label="Default kaon LO RG export",
                version="paper-mode:v1",
                source="quarkConstraints.paper_0710_1869.eft_deltaf2.rg",
                citation="BMU hep-ph/0005183 Q1 running path",
            ),
            ProvenanceRecord(
                record_id=hadronic.source_id,
                category="hadronic-input-bundle",
                label="Default kaon hadronic input bundle",
                version=hadronic.input_provenance_mode_id,
                source=hadronic.source_id,
                citation=hadronic.notes,
            ),
            _provenance_record_from_source(mass_source),
            _provenance_record_from_source(decay_constant_source),
            _provenance_record_from_source(bag_parameter_source),
        ),
    )

    return Paper07101869ArtifactExportSet(
        wilson_bundle=wilson_bundle,
        hadronic_bundle=hadronic_bundle,
        observable_bundle=observable_bundle,
        provenance_bundle=provenance_bundle,
    )


def build_strict_paper_0710_1869_kaon_artifact_export_set() -> Paper07101869ArtifactExportSet:
    """Build the strict-paper artifact quartet with default kaon numerics and strict IDs."""

    default_export_set = build_default_paper_0710_1869_kaon_artifact_export_set()
    default_wilson = default_export_set.wilson_bundle
    default_hadronic = default_export_set.hadronic_bundle
    default_observable = default_export_set.observable_bundle

    strict_mass_source = _artifact_source_with_id(
        default_hadronic.mass_source,
        source_id=STRICT_PAPER_MASS_SOURCE_ID,
    )
    strict_decay_constant_source = _artifact_source_with_id(
        default_hadronic.decay_constant_source,
        source_id=STRICT_PAPER_DECAY_CONSTANT_SOURCE_ID,
    )
    strict_bag_parameter_source = _artifact_source_with_id(
        default_hadronic.bag_parameter_source,
        source_id=STRICT_PAPER_BAG_SOURCE_ID,
    )

    strict_wilson_bundle = replace(
        default_wilson,
        metadata=ArtifactMetadata.create(
            schema_name=WILSON_BUNDLE_SCHEMA,
            bundle_id=STRICT_PAPER_WILSON_BUNDLE_ID,
            point_id=STRICT_PAPER_ARTIFACT_POINT_ID,
        ),
        provenance_ids=(
            STRICT_PAPER_CONVENTIONS_PROVENANCE_ID,
            STRICT_PAPER_MATCHING_PROVENANCE_ID,
            STRICT_PAPER_RG_PROVENANCE_ID,
        ),
    )
    strict_hadronic_bundle = replace(
        default_hadronic,
        metadata=ArtifactMetadata.create(
            schema_name=HADRONIC_BUNDLE_SCHEMA,
            bundle_id=STRICT_PAPER_HADRONIC_BUNDLE_ID,
            point_id=STRICT_PAPER_ARTIFACT_POINT_ID,
        ),
        benchmark_id="kaon.np_only.strict_paper.v1",
        source_id=STRICT_PAPER_HADRONIC_SOURCE_ID,
        hadronic_source_id=STRICT_PAPER_HADRONIC_SOURCE_ID,
        mass_source_id=STRICT_PAPER_MASS_SOURCE_ID,
        decay_constant_source_id=STRICT_PAPER_DECAY_CONSTANT_SOURCE_ID,
        bag_parameter_source_id=STRICT_PAPER_BAG_SOURCE_ID,
        source_wilson_bundle_id=STRICT_PAPER_WILSON_BUNDLE_ID,
        mass_source=strict_mass_source,
        decay_constant_source=strict_decay_constant_source,
        bag_parameter_source=strict_bag_parameter_source,
        provenance_ids=(
            STRICT_PAPER_HADRONIC_SOURCE_ID,
            STRICT_PAPER_MASS_SOURCE_ID,
            STRICT_PAPER_DECAY_CONSTANT_SOURCE_ID,
            STRICT_PAPER_BAG_SOURCE_ID,
        ),
    )
    strict_observable_bundle = replace(
        default_observable,
        metadata=ArtifactMetadata.create(
            schema_name=OBSERVABLE_BUNDLE_SCHEMA,
            bundle_id=STRICT_PAPER_OBSERVABLE_BUNDLE_ID,
            point_id=STRICT_PAPER_ARTIFACT_POINT_ID,
        ),
        source_wilson_bundle_id=STRICT_PAPER_WILSON_BUNDLE_ID,
        source_hadronic_bundle_id=STRICT_PAPER_HADRONIC_BUNDLE_ID,
        provenance_ids=(
            STRICT_PAPER_CONVENTIONS_PROVENANCE_ID,
            STRICT_PAPER_MATCHING_PROVENANCE_ID,
            STRICT_PAPER_RG_PROVENANCE_ID,
            STRICT_PAPER_HADRONIC_SOURCE_ID,
            STRICT_PAPER_MASS_SOURCE_ID,
            STRICT_PAPER_DECAY_CONSTANT_SOURCE_ID,
            STRICT_PAPER_BAG_SOURCE_ID,
        ),
    )
    strict_provenance_bundle = ProvenanceBundleV1(
        metadata=ArtifactMetadata.create(
            schema_name=PROVENANCE_BUNDLE_SCHEMA,
            bundle_id=STRICT_PAPER_PROVENANCE_BUNDLE_ID,
            point_id=STRICT_PAPER_ARTIFACT_POINT_ID,
        ),
        records=(
            ProvenanceRecord(
                record_id=STRICT_PAPER_CONVENTIONS_PROVENANCE_ID,
                category="paper",
                label="strict paper conventions",
                version="arXiv:0710.1869",
                source="doi:10.48550/arXiv.0710.1869",
                citation="arXiv:0710.1869",
            ),
            ProvenanceRecord(
                record_id=STRICT_PAPER_MATCHING_PROVENANCE_ID,
                category="code",
                label="strict paper matching export",
                version="paper-mode:v1",
                source="quarkConstraints.paper_0710_1869.eft_deltaf2.matching_kkgluon",
                citation="arXiv:0710.1869 matching path",
            ),
            ProvenanceRecord(
                record_id=STRICT_PAPER_RG_PROVENANCE_ID,
                category="code",
                label="strict paper LO RG export",
                version="paper-mode:v1",
                source="quarkConstraints.paper_0710_1869.eft_deltaf2.rg",
                citation="BMU hep-ph/0005183 Q1 running path",
            ),
            ProvenanceRecord(
                record_id=STRICT_PAPER_HADRONIC_SOURCE_ID,
                category="hadronic-input-bundle",
                label="strict paper hadronic input bundle",
                version=str(strict_hadronic_bundle.input_provenance_mode_id),
                source="strict paper hadronic source",
                citation=strict_hadronic_bundle.notes,
            ),
            _provenance_record_from_source(strict_mass_source),
            _provenance_record_from_source(strict_decay_constant_source),
            _provenance_record_from_source(strict_bag_parameter_source),
        ),
    )
    return Paper07101869ArtifactExportSet(
        wilson_bundle=strict_wilson_bundle,
        hadronic_bundle=strict_hadronic_bundle,
        observable_bundle=strict_observable_bundle,
        provenance_bundle=strict_provenance_bundle,
    )


def build_default_paper_0710_1869_kaon_hadronic_artifact_bundle() -> HadronicArtifactBundleV1:
    """Return the canonical hadronic artifact bundle for the default kaon benchmark."""

    return build_default_paper_0710_1869_kaon_artifact_export_set().hadronic_bundle


def default_paper_0710_1869_kaon_artifact_export_set() -> Paper07101869ArtifactExportSet:
    """Compatibility alias for the canonical default kaon artifact export set."""

    return build_default_paper_0710_1869_kaon_artifact_export_set()


def strict_paper_0710_1869_kaon_artifact_export_set() -> Paper07101869ArtifactExportSet:
    """Compatibility alias for the strict-paper kaon artifact export set."""

    return build_strict_paper_0710_1869_kaon_artifact_export_set()


def default_paper_0710_1869_kaon_artifact_export_paths(
    root_dir: str | Path | None = None,
) -> Paper07101869ArtifactExportPaths:
    """Return the deterministic filesystem layout for the default kaon artifact set."""

    resolved_root = DEFAULT_KAON_ARTIFACT_DIR if root_dir is None else Path(root_dir)
    return Paper07101869ArtifactExportPaths(
        root_dir=resolved_root,
        wilson_path=resolved_root / DEFAULT_KAON_ARTIFACT_FILENAMES["wilsons"],
        hadronic_path=resolved_root / DEFAULT_KAON_ARTIFACT_FILENAMES["hadronic"],
        observable_path=resolved_root / DEFAULT_KAON_ARTIFACT_FILENAMES["observables"],
        provenance_path=resolved_root / DEFAULT_KAON_ARTIFACT_FILENAMES["provenance"],
    )


def strict_paper_0710_1869_kaon_artifact_export_paths(
    root_dir: str | Path | None = None,
) -> Paper07101869ArtifactExportPaths:
    """Return the deterministic filesystem layout for the strict-paper kaon artifact set."""

    resolved_root = STRICT_PAPER_ARTIFACT_DIR if root_dir is None else Path(root_dir)
    return Paper07101869ArtifactExportPaths(
        root_dir=resolved_root,
        wilson_path=resolved_root / DEFAULT_KAON_ARTIFACT_FILENAMES["wilsons"],
        hadronic_path=resolved_root / DEFAULT_KAON_ARTIFACT_FILENAMES["hadronic"],
        observable_path=resolved_root / DEFAULT_KAON_ARTIFACT_FILENAMES["observables"],
        provenance_path=resolved_root / DEFAULT_KAON_ARTIFACT_FILENAMES["provenance"],
    )


def write_default_paper_0710_1869_kaon_artifact_exports(
    root_dir: str | Path | None = None,
) -> Paper07101869ArtifactExportPaths:
    """Write the canonical default kaon artifact set under a deterministic results path."""

    export_set = build_default_paper_0710_1869_kaon_artifact_export_set()
    paths = default_paper_0710_1869_kaon_artifact_export_paths(root_dir=root_dir)
    paths.root_dir.mkdir(parents=True, exist_ok=True)
    export_set.wilson_bundle.write_json(paths.wilson_path)
    export_set.hadronic_bundle.write_json(paths.hadronic_path)
    export_set.observable_bundle.write_json(paths.observable_path)
    export_set.provenance_bundle.write_json(paths.provenance_path)
    return paths


def write_strict_paper_0710_1869_kaon_artifact_exports(
    root_dir: str | Path | None = None,
) -> Paper07101869ArtifactExportPaths:
    """Write the strict-paper kaon artifact set under a deterministic results path."""

    export_set = build_strict_paper_0710_1869_kaon_artifact_export_set()
    paths = strict_paper_0710_1869_kaon_artifact_export_paths(root_dir=root_dir)
    paths.root_dir.mkdir(parents=True, exist_ok=True)
    export_set.wilson_bundle.write_json(paths.wilson_path)
    export_set.hadronic_bundle.write_json(paths.hadronic_path)
    export_set.observable_bundle.write_json(paths.observable_path)
    export_set.provenance_bundle.write_json(paths.provenance_path)
    return paths


__all__ = [
    "ARTIFACT_SCHEMA_VERSION",
    "ArtifactMetadata",
    "ArtifactScale",
    "ArtifactSchemaError",
    "ArtifactSourceRecord",
    "ComplexValue",
    "DEFAULT_BAG_PARAMETER_CONVERSION_FORMULA_ID",
    "DEFAULT_KAON_ARTIFACT_DIR",
    "DEFAULT_KAON_ARTIFACT_FILENAMES",
    "DEFAULT_KAON_ARTIFACT_POINT_ID",
    "DEFAULT_KAON_HADRONIC_BUNDLE_ID",
    "DEFAULT_KAON_OBSERVABLE_BUNDLE_ID",
    "DEFAULT_KAON_PROVENANCE_BUNDLE_ID",
    "DEFAULT_KAON_SUPPORTED_OPERATORS",
    "DEFAULT_KAON_UNSUPPORTED_OPERATORS",
    "DEFAULT_KAON_WILSON_BUNDLE_ID",
    "DEFAULT_MATCHING_PROVENANCE_ID",
    "DEFAULT_PAPER_PROVENANCE_ID",
    "DEFAULT_RG_PROVENANCE_ID",
    "DELTA_M_RELATION_ID",
    "HADRONIC_BUNDLE_SCHEMA",
    "HadronicArtifactBundleV1",
    "OBSERVABLE_BUNDLE_SCHEMA",
    "ObservableArtifactBundleV1",
    "ObservableRecord",
    "PAPER_MODE",
    "PROVENANCE_BUNDLE_SCHEMA",
    "Paper07101869ArtifactExportPaths",
    "Paper07101869ArtifactExportSet",
    "PaperArtifactBundle",
    "ProvenanceBundleV1",
    "ProvenanceRecord",
    "STRICT_PAPER_ARTIFACT_DIR",
    "STRICT_PAPER_ARTIFACT_POINT_ID",
    "STRICT_PAPER_BAG_SOURCE_ID",
    "STRICT_PAPER_CONVENTIONS_PROVENANCE_ID",
    "STRICT_PAPER_DECAY_CONSTANT_SOURCE_ID",
    "STRICT_PAPER_HADRONIC_BUNDLE_ID",
    "STRICT_PAPER_HADRONIC_SOURCE_ID",
    "STRICT_PAPER_MASS_SOURCE_ID",
    "STRICT_PAPER_MATCHING_PROVENANCE_ID",
    "STRICT_PAPER_OBSERVABLE_BUNDLE_ID",
    "STRICT_PAPER_PROVENANCE_BUNDLE_ID",
    "STRICT_PAPER_RG_PROVENANCE_ID",
    "STRICT_PAPER_WILSON_BUNDLE_ID",
    "WILSON_BUNDLE_SCHEMA",
    "WilsonArtifactBundleV1",
    "WilsonCoefficientRecord",
    "artifact_from_dict",
    "build_default_paper_0710_1869_kaon_artifact_export_set",
    "build_default_paper_0710_1869_kaon_hadronic_artifact_bundle",
    "build_strict_paper_0710_1869_kaon_artifact_export_set",
    "default_paper_0710_1869_kaon_artifact_export_paths",
    "default_paper_0710_1869_kaon_artifact_export_set",
    "read_artifact",
    "strict_paper_0710_1869_kaon_artifact_export_paths",
    "strict_paper_0710_1869_kaon_artifact_export_set",
    "write_default_paper_0710_1869_kaon_artifact_exports",
    "write_strict_paper_0710_1869_kaon_artifact_exports",
]
