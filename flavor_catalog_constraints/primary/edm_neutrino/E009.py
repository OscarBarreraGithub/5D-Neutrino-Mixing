"""E009 - Weinberg CP-odd three-gluon operator reference bounds.

Physics
-------
The catalogued E009 entries are neutron-EDM-derived benchmark translations
for the CP-odd Weinberg three-gluon operator.  The sidecar records both the
Pospelov-Ritz ``w(1 GeV)`` convention and the Haisch-Hala ``C_6``/``O_6``
convention:

    |w(1 GeV)| < 4.1e-11 GeV^-2
    |C_6|      < 1.2e-11 GeV^-2

These are central, single-source-at-a-time translations from the neutron EDM
limit in ``flavor_catalog/processes/edm_neutrino/E009.yaml``.  A rigorous RS
prediction is **not** available in this scaffold.  It would require

1. RS CP-odd gluonic matching onto the Weinberg operator in a chosen basis;
   and
2. non-perturbative hadronic matrix elements translating that gluonic CP-odd
   source into ``d_n``.

Both ingredients are flagged ``NEEDS-HUMAN-PHYSICS`` in the returned
diagnostics.  This file records the catalogued anchors only; it does not fake
a hadronic calculation and does not read colored CP-odd inputs from
``ParameterPoint``.

Severity
--------
INFO.  The neutron EDM is a stringent observed bound, but applying E009 as a
hard veto would require the missing RS gluonic matching and highly uncertain
hadronic matrix-element convention.  The returned ``passes`` value is
bookkeeping only and must not veto scan points.

Catalog sidecar
---------------
``flavor_catalog/processes/edm_neutrino/E009.yaml`` is the source of truth for
the neutron EDM anchor and Weinberg benchmark bounds.  Numeric values below
are loaded through the scaffold anchor helpers, including virtualized views
over the sidecar's nested ``measured_experimental_anchor`` mapping and
``values`` list.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Mapping

from flavor_catalog_constraints import anchors as anchor_scaffold
from flavor_catalog_constraints.anchors import (
    Anchor,
    AnchorError,
    find_block,
    load_anchor,
    load_full_yaml,
    load_pdg_block,
)
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.weinberg_operator import (
    WEINBERG_OPERATOR_HADRONIC_NEEDS_HUMAN_PHYSICS,
    WEINBERG_OPERATOR_RS_GLUONIC_MATCHING_NEEDS_HUMAN_PHYSICS,
    WEINBERG_OPERATOR_STUB_MODEL_V1,
    weinberg_operator_reference_bound,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "edm_neutrino"
_MEASURED_EXPERIMENTAL_ANCHOR_BLOCK = "measured_experimental_anchor"
_NEUTRON_EDM_BLOCK = "neutron_edm"
_PRIMARY_EXPERIMENT_BLOCK = "primary_experiment"
_VALUES_BLOCK = "values"
_AUXILIARY_THEORY_INPUTS_BLOCK = "auxiliary_theory_inputs"
_PAPER_ERA_REFERENCE_BLOCK = "paper_era_reference"
_PR_RESPONSE_VALUE_ID = "PospelovRitz2005:E009:weinberg_normalization"
_PR_BOUND_VALUE_ID = "PDG2026-PospelovRitz2005:E009:w_bound"
_HH_RESPONSE_VALUE_ID = "HaischHala2019:E009:o6_normalization"
_HH_BOUND_VALUE_ID = "PDG2026-HaischHala2019:E009:c6_bound"
_EXPECTED_LIMIT_OPERATOR = "<"
_EXPECTED_NEUTRON_EDM_UNITS = "e cm"
_EXPECTED_PRIMARY_MEASUREMENT_UNITS = "10^-26 e cm"
_EXPECTED_COEFFICIENT_UNITS = "MeV"
_EXPECTED_BOUND_UNITS = "GeV^-2"
_MEV_TO_GEV = 1.0e-3
_ROUNDING_REL_TOL = 0.05
_REFERENCE_BOUND_VALUE_ID = _HH_BOUND_VALUE_ID
_BUDGET_SOURCE = (
    "flavor_catalog/processes/edm_neutrino/E009.yaml values."
    "PDG2026-HaischHala2019:E009:c6_bound"
)
_PARAMETRIZATION_CITATION = (
    "Weinberg Phys. Rev. Lett. 63, 2333 (1989); "
    "Pospelov and Ritz Annals Phys. 318, 119 (2005), arXiv:hep-ph/0504231; "
    "Haisch and Hala JHEP 11, 154 (2019), arXiv:1909.08955"
)
_NEEDS_HUMAN_PHYSICS = (
    WEINBERG_OPERATOR_HADRONIC_NEEDS_HUMAN_PHYSICS,
    WEINBERG_OPERATOR_RS_GLUONIC_MATCHING_NEEDS_HUMAN_PHYSICS,
)


@dataclass(frozen=True)
class ExperimentalNeutronEDMAnchor:
    """Typed neutron EDM anchor used to derive the Weinberg benchmarks."""

    anchor: Anchor
    parent_block_key: str
    nested_block_key: str
    limit_operator: str | None
    confidence_level: str | None
    used_measurement: str | None
    table_value: float | None
    table_units: str | None
    measurement_value_e_cm: float
    statistical_uncertainty_e_cm: float
    systematic_uncertainty_e_cm: float
    total_uncertainty_e_cm: float
    measurement_units: str | None
    value_summary: str | None

    @property
    def block_key(self) -> str:
        """Dotted sidecar block path for diagnostics."""
        return f"{self.parent_block_key}.{self.nested_block_key}"

    @property
    def value(self) -> float:
        """Upper limit on ``|d_n|`` in e cm."""
        return self.anchor.value

    @property
    def source_url(self) -> str | None:
        """Primary source URL for the neutron EDM anchor."""
        return self.anchor.source_url


@dataclass(frozen=True)
class PrimaryNeutronEDMMeasurement:
    """Typed primary neutron EDM measurement, converted to e cm."""

    anchor: Anchor
    parent_block_key: str
    nested_block_key: str
    central_value_e_cm: float
    statistical_uncertainty_e_cm: float
    systematic_uncertainty_e_cm: float
    total_uncertainty_e_cm: float
    limit_value_e_cm: float
    limit_operator: str | None
    confidence_level: str | None
    units: str | None
    limit_units: str | None
    value_summary: str | None

    @property
    def block_key(self) -> str:
        """Dotted sidecar block path for diagnostics."""
        return f"{self.parent_block_key}.{self.nested_block_key}"


@dataclass(frozen=True)
class WeinbergResponseBenchmark:
    """Typed neutron-EDM response coefficient for one Weinberg convention."""

    anchor: Anchor
    value_id: str
    coefficient_gev: float
    coefficient_units: str | None
    coefficient_scale: str | None
    coefficient_relative_uncertainty: float | None
    operator_normalization: str | None
    neutron_edm_formula: str | None
    assumptions: str | None

    @property
    def block_key(self) -> str:
        """Virtual sidecar block key used by ``load_anchor``."""
        return self.anchor.block_key


@dataclass(frozen=True)
class WeinbergBenchmarkBound:
    """Typed single-source benchmark bound on a Weinberg coefficient."""

    anchor: Anchor
    value_id: str
    derived_bound_observable: str | None
    derived_bound_operator: str | None
    derived_bound_value: float
    derived_bound_units: str | None
    derivation: str | None
    assumptions: str | None
    snapshot_paths: tuple[str, ...]

    @property
    def block_key(self) -> str:
        """Virtual sidecar block key used by ``load_anchor``."""
        return self.anchor.block_key

    @property
    def value(self) -> float:
        """Upper bound on this Weinberg coefficient convention."""
        return self.anchor.value

    @property
    def source_url(self) -> str | None:
        """Primary source URL for this benchmark bound."""
        return self.anchor.source_url


@dataclass(frozen=True)
class NeutronLimitConversion:
    """Sidecar-recorded neutron EDM unit conversion used by E009."""

    source: str | None
    neutron_limit_e_cm: float
    neutron_limit_gev_inverse: float
    conversion_summary: str | None


@dataclass(frozen=True)
class E009Context:
    """Theory context carried for provenance, not calculation."""

    block_key: str
    source: str | None
    year: int | None
    source_url: str | None
    value_summary: str | None
    values: tuple[Mapping[str, Any], ...]
    snapshot_path: str | None


@dataclass(frozen=True)
class E009Anchor:
    """Typed E009 anchor bundle and non-vetoing budget convention."""

    neutron_edm_limit: ExperimentalNeutronEDMAnchor
    primary_experiment: PrimaryNeutronEDMMeasurement
    pospelov_ritz_response: WeinbergResponseBenchmark
    pospelov_ritz_w_bound: WeinbergBenchmarkBound
    haisch_hala_response: WeinbergResponseBenchmark
    haisch_hala_c6_bound: WeinbergBenchmarkBound
    neutron_limit_conversion: NeutronLimitConversion
    original_operator_context: E009Context
    cfw_context: E009Context
    global_analysis_context: E009Context
    composite_warped_context: E009Context
    lattice_gradient_flow_context: E009Context

    @property
    def reference_bound(self) -> WeinbergBenchmarkBound:
        """Default bookkeeping bound: the Haisch-Hala ``C_6`` benchmark."""
        return self.haisch_hala_c6_bound

    @property
    def value(self) -> float:
        """Default benchmark bound in GeV^-2."""
        return self.reference_bound.value

    @property
    def budget(self) -> float:
        """Non-vetoing bookkeeping budget: the default benchmark bound."""
        return self.reference_bound.value

    @property
    def source_url(self) -> str | None:
        """Source URL for the default benchmark bound."""
        return self.reference_bound.source_url


def _optional_str(value: Any) -> str | None:
    return None if value is None else str(value)


def _optional_int(value: Any) -> int | None:
    if value is None:
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _optional_float(value: Any, *, process_id: str, field_name: str) -> float | None:
    if value is None:
        return None
    return _required_float(value, process_id=process_id, field_name=field_name)


def _required_float(value: Any, *, process_id: str, field_name: str) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError) as exc:
        raise AnchorError(
            f"{process_id}: E009 anchor field {field_name!r}={value!r} "
            "is not numeric"
        ) from exc
    if not math.isfinite(out):
        raise AnchorError(f"{process_id}: E009 anchor field {field_name!r} is not finite")
    return out


def _positive_float(value: Any, *, process_id: str, field_name: str) -> float:
    out = _required_float(value, process_id=process_id, field_name=field_name)
    if out <= 0.0:
        raise AnchorError(f"{process_id}: E009 anchor field {field_name!r} must be positive")
    return out


def _scale_to_e_cm(units: str | None, *, process_id: str, field_name: str) -> float:
    if units == "e cm":
        return 1.0
    if units == _EXPECTED_PRIMARY_MEASUREMENT_UNITS:
        return 1.0e-26
    raise AnchorError(f"{process_id}: unsupported E009 EDM units for {field_name}: {units!r}")


def _snapshot_paths(value: Any) -> tuple[str, ...]:
    if value is None:
        return ()
    if isinstance(value, list):
        return tuple(str(item) for item in value)
    return (str(value),)


def _load_virtual_anchor(
    process_id: str,
    *,
    block_key: str,
    sub: Mapping[str, Any],
    value_key: str = "value",
) -> Anchor:
    virtual_block = {block_key: dict(sub)}
    original_load_pdg_block = anchor_scaffold.load_pdg_block

    def _load_virtual_pdg_block(
        request_process_id: str,
        **kwargs: Any,
    ) -> Mapping[str, Any]:
        if request_process_id == process_id and kwargs.get("family") == _FAMILY:
            return virtual_block
        return original_load_pdg_block(request_process_id, **kwargs)

    anchor_scaffold.load_pdg_block = _load_virtual_pdg_block
    try:
        anchor = load_anchor(
            process_id,
            family=_FAMILY,
            candidates=(block_key,),
            value_key=value_key,
        )
    finally:
        anchor_scaffold.load_pdg_block = original_load_pdg_block

    if anchor.block_key != block_key:
        raise AnchorError(
            f"{process_id}: load_anchor selected {anchor.block_key!r}, "
            f"expected {block_key!r} for E009 anchor"
        )
    return anchor


def _load_nested_mapping(
    process_id: str,
    *,
    parent_key: str,
    nested_key: str,
) -> Mapping[str, Any]:
    pdg = load_pdg_block(process_id, family=_FAMILY)
    parent = find_block(pdg, (parent_key,), process_id=process_id)
    return find_block(parent, (nested_key,), process_id=f"{process_id}.{parent_key}")


def _load_nested_anchor(
    process_id: str,
    *,
    parent_key: str,
    nested_key: str,
    value_key: str = "value",
) -> tuple[Anchor, Mapping[str, Any]]:
    sub = _load_nested_mapping(
        process_id,
        parent_key=parent_key,
        nested_key=nested_key,
    )
    block_key = f"{parent_key}.{nested_key}"
    return (
        _load_virtual_anchor(
            process_id,
            block_key=block_key,
            sub=sub,
            value_key=value_key,
        ),
        sub,
    )


def _load_value_entry(process_id: str, value_id: str) -> Mapping[str, Any]:
    pdg = load_pdg_block(process_id, family=_FAMILY)
    values = pdg.get(_VALUES_BLOCK)
    if not isinstance(values, list):
        raise AnchorError(f"{process_id}: pdg_or_equivalent.values is not a list")
    for item in values:
        if not isinstance(item, Mapping):
            raise AnchorError(f"{process_id}: pdg_or_equivalent.values contains non-mapping")
        if item.get("value_id") == value_id:
            return item
    raise AnchorError(f"{process_id}: missing E009 value_id {value_id!r}")


def _load_value_anchor(
    process_id: str,
    *,
    value_id: str,
    value_key: str = "value",
) -> tuple[Anchor, Mapping[str, Any]]:
    sub = _load_value_entry(process_id, value_id)
    block_key = f"{_VALUES_BLOCK}.{value_id}"
    return (
        _load_virtual_anchor(
            process_id,
            block_key=block_key,
            sub=sub,
            value_key=value_key,
        ),
        sub,
    )


def _load_neutron_edm_anchor(process_id: str) -> ExperimentalNeutronEDMAnchor:
    anchor, sub = _load_nested_anchor(
        process_id,
        parent_key=_MEASURED_EXPERIMENTAL_ANCHOR_BLOCK,
        nested_key=_NEUTRON_EDM_BLOCK,
    )
    if anchor.value <= 0.0:
        raise AnchorError(f"{process_id}: neutron EDM limit must be positive")
    if anchor.units != _EXPECTED_NEUTRON_EDM_UNITS:
        raise AnchorError(
            f"{process_id}: expected neutron EDM units {_EXPECTED_NEUTRON_EDM_UNITS!r}, "
            f"got {anchor.units!r}"
        )
    limit_operator = _optional_str(sub.get("limit_operator"))
    if limit_operator != _EXPECTED_LIMIT_OPERATOR:
        raise AnchorError(
            f"{process_id}: expected limit_operator {_EXPECTED_LIMIT_OPERATOR!r}, "
            f"got {limit_operator!r}"
        )
    measurement_units = _optional_str(sub.get("measurement_units"))
    measurement_scale = _scale_to_e_cm(
        measurement_units,
        process_id=process_id,
        field_name=f"{anchor.block_key}.measurement_units",
    )
    measurement = _required_float(
        sub.get("measurement_value"),
        process_id=process_id,
        field_name=f"{anchor.block_key}.measurement_value",
    )
    stat = _positive_float(
        sub.get("statistical_uncertainty"),
        process_id=process_id,
        field_name=f"{anchor.block_key}.statistical_uncertainty",
    )
    syst = _positive_float(
        sub.get("systematic_uncertainty"),
        process_id=process_id,
        field_name=f"{anchor.block_key}.systematic_uncertainty",
    )
    return ExperimentalNeutronEDMAnchor(
        anchor=anchor,
        parent_block_key=_MEASURED_EXPERIMENTAL_ANCHOR_BLOCK,
        nested_block_key=_NEUTRON_EDM_BLOCK,
        limit_operator=limit_operator,
        confidence_level=_optional_str(sub.get("confidence_level")),
        used_measurement=_optional_str(sub.get("used_measurement")),
        table_value=_optional_float(
            sub.get("table_value"),
            process_id=process_id,
            field_name=f"{anchor.block_key}.table_value",
        ),
        table_units=_optional_str(sub.get("table_units")),
        measurement_value_e_cm=float(measurement * measurement_scale),
        statistical_uncertainty_e_cm=float(stat * measurement_scale),
        systematic_uncertainty_e_cm=float(syst * measurement_scale),
        total_uncertainty_e_cm=float(
            math.sqrt(stat * stat + syst * syst) * measurement_scale
        ),
        measurement_units=measurement_units,
        value_summary=_optional_str(sub.get("value_summary")),
    )


def _load_primary_measurement(
    process_id: str,
    *,
    canonical_limit_e_cm: float,
    limit_operator: str | None,
    confidence_level: str | None,
) -> PrimaryNeutronEDMMeasurement:
    anchor, sub = _load_nested_anchor(
        process_id,
        parent_key=_MEASURED_EXPERIMENTAL_ANCHOR_BLOCK,
        nested_key=_PRIMARY_EXPERIMENT_BLOCK,
    )
    units = _optional_str(sub.get("units"))
    limit_units = _optional_str(sub.get("limit_units"))
    measurement_scale = _scale_to_e_cm(
        units,
        process_id=process_id,
        field_name=f"{anchor.block_key}.units",
    )
    limit_scale = _scale_to_e_cm(
        limit_units,
        process_id=process_id,
        field_name=f"{anchor.block_key}.limit_units",
    )
    stat = _positive_float(
        sub.get("statistical_uncertainty"),
        process_id=process_id,
        field_name=f"{anchor.block_key}.statistical_uncertainty",
    )
    syst = _positive_float(
        sub.get("systematic_uncertainty"),
        process_id=process_id,
        field_name=f"{anchor.block_key}.systematic_uncertainty",
    )
    central = _required_float(
        anchor.value,
        process_id=process_id,
        field_name=f"{anchor.block_key}.value",
    ) * measurement_scale
    limit_value = _positive_float(
        sub.get("limit_value"),
        process_id=process_id,
        field_name=f"{anchor.block_key}.limit_value",
    ) * limit_scale
    if not math.isclose(limit_value, canonical_limit_e_cm, rel_tol=1.0e-12):
        raise AnchorError(
            f"{process_id}: primary experiment limit {limit_value} e cm "
            f"does not match neutron anchor {canonical_limit_e_cm} e cm"
        )
    primary_operator = _optional_str(sub.get("limit_operator"))
    if primary_operator != limit_operator:
        raise AnchorError(
            f"{process_id}: primary experiment limit_operator {primary_operator!r} "
            f"does not match neutron anchor {limit_operator!r}"
        )
    primary_cl = _optional_str(sub.get("confidence_level"))
    if primary_cl != confidence_level:
        raise AnchorError(
            f"{process_id}: primary experiment confidence_level {primary_cl!r} "
            f"does not match neutron anchor {confidence_level!r}"
        )
    return PrimaryNeutronEDMMeasurement(
        anchor=anchor,
        parent_block_key=_MEASURED_EXPERIMENTAL_ANCHOR_BLOCK,
        nested_block_key=_PRIMARY_EXPERIMENT_BLOCK,
        central_value_e_cm=float(central),
        statistical_uncertainty_e_cm=float(stat * measurement_scale),
        systematic_uncertainty_e_cm=float(syst * measurement_scale),
        total_uncertainty_e_cm=float(math.sqrt(stat * stat + syst * syst) * measurement_scale),
        limit_value_e_cm=float(limit_value),
        limit_operator=primary_operator,
        confidence_level=primary_cl,
        units=units,
        limit_units=limit_units,
        value_summary=_optional_str(sub.get("value_summary")),
    )


def _load_response_benchmark(
    process_id: str,
    *,
    value_id: str,
) -> WeinbergResponseBenchmark:
    anchor, sub = _load_value_anchor(process_id, value_id=value_id)
    coefficient_units = _optional_str(sub.get("coefficient_units", anchor.units))
    if coefficient_units != _EXPECTED_COEFFICIENT_UNITS:
        raise AnchorError(
            f"{process_id}: expected coefficient units {_EXPECTED_COEFFICIENT_UNITS!r} "
            f"for {value_id}, got {coefficient_units!r}"
        )
    coefficient_mev = _positive_float(
        sub.get("coefficient_value", anchor.value),
        process_id=process_id,
        field_name=f"{anchor.block_key}.coefficient_value",
    )
    if not math.isclose(coefficient_mev, anchor.value, rel_tol=1.0e-12):
        raise AnchorError(
            f"{process_id}: response anchor value {anchor.value} does not match "
            f"coefficient_value {coefficient_mev} for {value_id}"
        )
    return WeinbergResponseBenchmark(
        anchor=anchor,
        value_id=value_id,
        coefficient_gev=float(coefficient_mev * _MEV_TO_GEV),
        coefficient_units=coefficient_units,
        coefficient_scale=_optional_str(sub.get("coefficient_scale")),
        coefficient_relative_uncertainty=_optional_float(
            sub.get("coefficient_relative_uncertainty"),
            process_id=process_id,
            field_name=f"{anchor.block_key}.coefficient_relative_uncertainty",
        ),
        operator_normalization=_optional_str(sub.get("operator_normalization")),
        neutron_edm_formula=_optional_str(sub.get("neutron_edm_formula")),
        assumptions=_optional_str(sub.get("assumptions")),
    )


def _load_bound_benchmark(
    process_id: str,
    *,
    value_id: str,
) -> WeinbergBenchmarkBound:
    anchor, sub = _load_value_anchor(
        process_id,
        value_id=value_id,
        value_key="derived_bound_value",
    )
    if anchor.value <= 0.0:
        raise AnchorError(f"{process_id}: Weinberg derived bound must be positive")
    operator = _optional_str(sub.get("derived_bound_operator"))
    if operator != _EXPECTED_LIMIT_OPERATOR:
        raise AnchorError(
            f"{process_id}: expected derived_bound_operator {_EXPECTED_LIMIT_OPERATOR!r} "
            f"for {value_id}, got {operator!r}"
        )
    units = _optional_str(sub.get("derived_bound_units", anchor.units))
    if units != _EXPECTED_BOUND_UNITS:
        raise AnchorError(
            f"{process_id}: expected Weinberg bound units {_EXPECTED_BOUND_UNITS!r} "
            f"for {value_id}, got {units!r}"
        )
    derived_bound = _positive_float(
        sub.get("derived_bound_value"),
        process_id=process_id,
        field_name=f"{anchor.block_key}.derived_bound_value",
    )
    if not math.isclose(derived_bound, anchor.value, rel_tol=1.0e-12):
        raise AnchorError(
            f"{process_id}: bound anchor value {anchor.value} does not match "
            f"derived_bound_value {derived_bound} for {value_id}"
        )
    return WeinbergBenchmarkBound(
        anchor=anchor,
        value_id=value_id,
        derived_bound_observable=_optional_str(sub.get("derived_bound_observable")),
        derived_bound_operator=operator,
        derived_bound_value=float(derived_bound),
        derived_bound_units=units,
        derivation=_optional_str(sub.get("derivation")),
        assumptions=_optional_str(sub.get("assumptions")),
        snapshot_paths=_snapshot_paths(sub.get("snapshot_paths", sub.get("snapshot_path"))),
    )


def _load_neutron_limit_conversion(process_id: str) -> NeutronLimitConversion:
    data = load_full_yaml(process_id, family=_FAMILY)
    auxiliary = data.get(_AUXILIARY_THEORY_INPUTS_BLOCK)
    if not isinstance(auxiliary, Mapping):
        raise AnchorError(f"{process_id}: missing auxiliary_theory_inputs mapping")
    sub = auxiliary.get("current_neutron_limit_conversion")
    if not isinstance(sub, Mapping):
        raise AnchorError(
            f"{process_id}: auxiliary_theory_inputs.current_neutron_limit_conversion "
            "is not a mapping"
        )
    return NeutronLimitConversion(
        source=_optional_str(sub.get("source")),
        neutron_limit_e_cm=_positive_float(
            sub.get("neutron_limit_e_cm"),
            process_id=process_id,
            field_name="current_neutron_limit_conversion.neutron_limit_e_cm",
        ),
        neutron_limit_gev_inverse=_positive_float(
            sub.get("neutron_limit_gev_inverse"),
            process_id=process_id,
            field_name="current_neutron_limit_conversion.neutron_limit_gev_inverse",
        ),
        conversion_summary=_optional_str(sub.get("conversion_summary")),
    )


def _load_context(process_id: str, *, parent_key: str, nested_key: str) -> E009Context:
    data = load_full_yaml(process_id, family=_FAMILY)
    parent = data.get(parent_key)
    if not isinstance(parent, Mapping):
        raise AnchorError(f"{process_id}: missing {parent_key} mapping")
    sub = parent.get(nested_key)
    if not isinstance(sub, Mapping):
        raise AnchorError(f"{process_id}: {parent_key}.{nested_key!r} is not a mapping")
    values_raw = sub.get("values")
    if values_raw is None:
        values: tuple[Mapping[str, Any], ...] = ()
    elif isinstance(values_raw, list):
        if not all(isinstance(item, Mapping) for item in values_raw):
            raise AnchorError(f"{process_id}: {parent_key}.{nested_key}.values invalid")
        values = tuple(dict(item) for item in values_raw)
    else:
        raise AnchorError(f"{process_id}: {parent_key}.{nested_key}.values is not a list")
    return E009Context(
        block_key=nested_key,
        source=_optional_str(sub.get("source")),
        year=_optional_int(sub.get("year")),
        source_url=_optional_str(sub.get("source_url")),
        value_summary=_optional_str(sub.get("value_summary")),
        values=values,
        snapshot_path=_optional_str(sub.get("snapshot_path")),
    )


def _assert_derived_bound_consistency(
    *,
    process_id: str,
    bound: WeinbergBenchmarkBound,
    response: WeinbergResponseBenchmark,
    neutron_limit_gev_inverse: float,
) -> float:
    unrounded = neutron_limit_gev_inverse / response.coefficient_gev
    if not math.isclose(bound.value, unrounded, rel_tol=_ROUNDING_REL_TOL, abs_tol=0.0):
        raise AnchorError(
            f"{process_id}: {bound.value_id} derived bound {bound.value} is "
            f"inconsistent with neutron limit/response coefficient ({unrounded})"
        )
    return float(unrounded)


def _load_e009_anchor(process_id: str) -> E009Anchor:
    neutron_edm_limit = _load_neutron_edm_anchor(process_id)
    primary_experiment = _load_primary_measurement(
        process_id,
        canonical_limit_e_cm=neutron_edm_limit.value,
        limit_operator=neutron_edm_limit.limit_operator,
        confidence_level=neutron_edm_limit.confidence_level,
    )
    neutron_limit_conversion = _load_neutron_limit_conversion(process_id)
    if not math.isclose(
        neutron_limit_conversion.neutron_limit_e_cm,
        neutron_edm_limit.value,
        rel_tol=1.0e-12,
    ):
        raise AnchorError(
            f"{process_id}: neutron limit conversion "
            f"{neutron_limit_conversion.neutron_limit_e_cm} e cm does not match "
            f"neutron EDM anchor {neutron_edm_limit.value} e cm"
        )

    pr_response = _load_response_benchmark(process_id, value_id=_PR_RESPONSE_VALUE_ID)
    pr_bound = _load_bound_benchmark(process_id, value_id=_PR_BOUND_VALUE_ID)
    hh_response = _load_response_benchmark(process_id, value_id=_HH_RESPONSE_VALUE_ID)
    hh_bound = _load_bound_benchmark(process_id, value_id=_HH_BOUND_VALUE_ID)
    _assert_derived_bound_consistency(
        process_id=process_id,
        bound=pr_bound,
        response=pr_response,
        neutron_limit_gev_inverse=neutron_limit_conversion.neutron_limit_gev_inverse,
    )
    _assert_derived_bound_consistency(
        process_id=process_id,
        bound=hh_bound,
        response=hh_response,
        neutron_limit_gev_inverse=neutron_limit_conversion.neutron_limit_gev_inverse,
    )
    return E009Anchor(
        neutron_edm_limit=neutron_edm_limit,
        primary_experiment=primary_experiment,
        pospelov_ritz_response=pr_response,
        pospelov_ritz_w_bound=pr_bound,
        haisch_hala_response=hh_response,
        haisch_hala_c6_bound=hh_bound,
        neutron_limit_conversion=neutron_limit_conversion,
        original_operator_context=_load_context(
            process_id,
            parent_key=_AUXILIARY_THEORY_INPUTS_BLOCK,
            nested_key="original_operator",
        ),
        cfw_context=_load_context(
            process_id,
            parent_key=_PAPER_ERA_REFERENCE_BLOCK,
            nested_key="cfw_baseline",
        ),
        global_analysis_context=_load_context(
            process_id,
            parent_key=_AUXILIARY_THEORY_INPUTS_BLOCK,
            nested_key="global_analysis_context",
        ),
        composite_warped_context=_load_context(
            process_id,
            parent_key=_AUXILIARY_THEORY_INPUTS_BLOCK,
            nested_key="composite_warped_context",
        ),
        lattice_gradient_flow_context=_load_context(
            process_id,
            parent_key=_AUXILIARY_THEORY_INPUTS_BLOCK,
            nested_key="lattice_gradient_flow_context",
        ),
    )


@register
class Constraint:
    """Catalogued non-vetoing Weinberg three-gluon operator stub."""

    process_id = "E009"
    severity = Severity.INFO
    observable = "Weinberg three-gluon operator reference bounds"

    def __init__(self) -> None:
        self.anchor = _load_e009_anchor(self.process_id)

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        del point
        reference = self.anchor.reference_bound
        reference_bound = weinberg_operator_reference_bound(
            bound_observable=(
                reference.derived_bound_observable
                or "Weinberg three-gluon coefficient"
            ),
            reference_bound_gev_minus2=self.anchor.budget,
            convention="Haisch-Hala O6/C6 benchmark convention",
        )

        return ConstraintResult(
            process_id=self.process_id,
            severity=self.severity,
            passes=True,
            predicted=None,
            sm_prediction=None,
            experimental=float(reference_bound.reference_bound_gev_minus2),
            ratio=None,
            budget=float(reference_bound.reference_bound_gev_minus2),
            notes=(
                "INFO-only E009 stub: loads the catalogued neutron-EDM-derived "
                "Weinberg three-gluon benchmark bounds and records the "
                "Haisch-Hala C6 convention as the default bookkeeping budget. "
                "No RS CP-odd gluonic matching and no hadronic neutron-EDM "
                "matrix-element calculation are performed; both required "
                "physics inputs are flagged NEEDS-HUMAN-PHYSICS, so this "
                "result is non-vetoing."
            ),
            diagnostics={
                "non_vetoing": True,
                "weinberg_prediction_evaluated": False,
                "hadronic_edm_calculation_evaluated": False,
                "no_hadronic_calculation": True,
                "no_rs_cp_odd_gluonic_matching": True,
                "stub_model": WEINBERG_OPERATOR_STUB_MODEL_V1,
                "budget_source": _BUDGET_SOURCE,
                "parametrization_citation": _PARAMETRIZATION_CITATION,
                "needs_human_physics": _NEEDS_HUMAN_PHYSICS,
                "needs_human_physics_hadronic_matrix_elements": (
                    WEINBERG_OPERATOR_HADRONIC_NEEDS_HUMAN_PHYSICS
                ),
                "needs_human_physics_rs_cp_odd_gluonic_matching": (
                    WEINBERG_OPERATOR_RS_GLUONIC_MATCHING_NEEDS_HUMAN_PHYSICS
                ),
                "hadronic_matrix_elements_available": False,
                "rs_cp_odd_gluonic_matching_available": False,
                "parameter_point_inputs_used": (),
                "required_low_energy_operators": (
                    "weinberg_three_gluon",
                    "quark_chromo_edm_thresholds",
                    "theta_term_policy",
                    "operator_mixing_convention",
                ),
                "reference_bound_value_id": reference.value_id,
                "reference_bound_block": reference.block_key,
                "reference_bound_observable": reference_bound.bound_observable,
                "reference_bound_gev_minus2": float(
                    reference_bound.reference_bound_gev_minus2
                ),
                "reference_bound_operator": reference.derived_bound_operator,
                "reference_bound_units": reference.derived_bound_units,
                "reference_bound_convention": reference_bound.convention,
                "pospelov_ritz_w_bound_block": self.anchor.pospelov_ritz_w_bound.block_key,
                "pospelov_ritz_w_bound_gev_minus2": float(
                    self.anchor.pospelov_ritz_w_bound.value
                ),
                "pospelov_ritz_response_coefficient_gev": float(
                    self.anchor.pospelov_ritz_response.coefficient_gev
                ),
                "pospelov_ritz_formula": (
                    self.anchor.pospelov_ritz_response.neutron_edm_formula
                ),
                "pospelov_ritz_assumptions": (
                    self.anchor.pospelov_ritz_w_bound.assumptions
                ),
                "haisch_hala_c6_bound_block": self.anchor.haisch_hala_c6_bound.block_key,
                "haisch_hala_c6_bound_gev_minus2": float(
                    self.anchor.haisch_hala_c6_bound.value
                ),
                "haisch_hala_response_coefficient_gev": float(
                    self.anchor.haisch_hala_response.coefficient_gev
                ),
                "haisch_hala_response_relative_uncertainty": (
                    self.anchor.haisch_hala_response.coefficient_relative_uncertainty
                ),
                "haisch_hala_formula": (
                    self.anchor.haisch_hala_response.neutron_edm_formula
                ),
                "haisch_hala_assumptions": self.anchor.haisch_hala_c6_bound.assumptions,
                "neutron_edm_anchor_block": self.anchor.neutron_edm_limit.block_key,
                "neutron_edm_limit_e_cm": float(self.anchor.neutron_edm_limit.value),
                "neutron_edm_limit_gev_inverse": float(
                    self.anchor.neutron_limit_conversion.neutron_limit_gev_inverse
                ),
                "neutron_edm_limit_source_url": (
                    self.anchor.neutron_edm_limit.source_url
                ),
                "primary_experiment_block": self.anchor.primary_experiment.block_key,
                "primary_experiment_central_e_cm": float(
                    self.anchor.primary_experiment.central_value_e_cm
                ),
                "primary_experiment_total_uncertainty_e_cm": float(
                    self.anchor.primary_experiment.total_uncertainty_e_cm
                ),
                "limit_operator": self.anchor.neutron_edm_limit.limit_operator,
                "confidence_level": self.anchor.neutron_edm_limit.confidence_level,
                "used_measurement": self.anchor.neutron_edm_limit.used_measurement,
                "current_neutron_limit_conversion": (
                    self.anchor.neutron_limit_conversion.conversion_summary
                ),
                "original_operator_context_summary": (
                    self.anchor.original_operator_context.value_summary
                ),
                "cfw_context_summary": self.anchor.cfw_context.value_summary,
                "global_analysis_context_summary": (
                    self.anchor.global_analysis_context.value_summary
                ),
                "composite_warped_context_summary": (
                    self.anchor.composite_warped_context.value_summary
                ),
                "lattice_gradient_flow_context_summary": (
                    self.anchor.lattice_gradient_flow_context.value_summary
                ),
            },
        )
