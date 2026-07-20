"""CR014 - four-top production with a top-philic vector mediator.

Physics
-------
The observable is the RS top-philic vector resonance mass proxy compared with
the catalogued CMS four-top two-lepton simplified-model exclusion,

    m(Z') > 850 GeV  for Gamma/m = 50%, observed 95% CL.

The reusable mass-vs-limit comparison machinery lives in
``quarkConstraints.collider_resonance`` and is reached only through the
``flavor_catalog_constraints.physics_adapters.collider_resonance`` boundary.

RS matching status
------------------
NEEDS-HUMAN-PHYSICS.  This is a width-dependent low-mass benchmark exclusion.
A faithful RS recast needs ``sigma(pp -> ttbar Z') * BR(Z' -> ttbar)``, the
top-philic coupling normalization, the total-width convention, interference,
four-top acceptance, and the SM four-top background likelihood.  The current
``ParameterPoint`` carries only a declared vector/Kaluza-Klein mass scale, so
CR014 uses the documented proxy ``m_Zprime_top_philic = kk_ew_mass_gev or
M_KK`` and flags the sigma*BR, width, acceptance, and background recasts in
the result diagnostics.

Severity
--------
HARD.  The active budget is the strongest applicable observed CMS
top-philic vector ``Gamma/m=50%`` mass-exclusion row loaded from
``CR014.yaml``.  The YAML value is 850 GeV and is converted to a 0.85 TeV
mass lower-bound budget for the shared collider-resonance comparison.  The
ratio is ``m_limit / m_Zprime_top_philic`` so ``ratio <= 1`` passes.

Catalog sidecar
---------------
``CR014.yaml`` stores numerical entries in a ``pdg_or_equivalent.values``
list.  This module adapts selected list entries into the scaffold
``load_anchor`` path and fails loudly if the expected value IDs, units, or
width-benchmark metadata are missing.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Mapping, Sequence

from flavor_catalog_constraints import anchors as anchor_scaffold
from flavor_catalog_constraints.anchors import (
    Anchor,
    AnchorError,
    load_anchor,
    load_full_yaml,
)
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.collider_resonance import (
    MASS_LOWER_BOUND,
    TOP_PHILIC_VECTOR_FOUR_TOP_PROXY_ASSUMPTION_V1,
    ColliderResonanceLimit,
    evaluate_collider_resonance_limit,
    resolve_top_philic_vector_mass_gev,
    top_philic_vector_four_top_prediction_from_mass_gev,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "collider_rs"
_MASS_EXTRA = "kk_ew_mass_gev"
_COUPLINGS_EXTRA = "quark_mass_basis_couplings"
_PARENT_KEY = "pdg_or_equivalent"
_SCAFFOLD_UNCERTAINTY_KEY = "__cr014_uncertainty_is_not_used__"
_MASS_UNITS_YAML = "GeV"
_MASS_UNITS_COMPARISON = "TeV"
_CROSS_SECTION_UNITS = "fb"
_RESONANCE = "top-philic vector mediator Z'"
_FINAL_STATE = "t tbar t tbar two-lepton"

_ACTIVE_CMS_VECTOR_50PCT = (
    "CMSB2G25005:CR014:top_philic_vector_50pct_width_excluded_up_to"
)
_ATLAS_4TOP_CROSS_SECTION = "PDG2025:CR014:ATLAS_4top_cross_section"
_CMS_4TOP_CROSS_SECTION = "PDG2025:CR014:CMS_4top_cross_section"
_ATLAS_TOP_PHILIC_XSEC_RANGE = (
    "ATLAS2024:CR014:top_philic_Zprime_cross_section_limit_range"
)
_REFERENCE_WIDTH_OVER_MASS = 0.50
_BUDGET_POLICY = (
    "active HARD budget is the strongest applicable observed CMS four-top "
    "two-lepton top-philic vector mass-exclusion row in CR014.yaml: "
    "Gamma/m=50%, converted from GeV to TeV. The low-mass, width-dependent "
    "mass edge is only a documented proxy for an RS four-top recast."
)


@dataclass(frozen=True)
class CR014MassLimitAnchor:
    """Typed CR014 mass-exclusion row routed through scaffold ``load_anchor``."""

    anchor: Anchor
    value_id: str
    entry_index: int
    display: str | None
    limit_type: str | None
    cl: str | None
    experiment: str | None
    source_key: str | None
    access_date: str | None
    sha256: str | None
    resonance_type: str | None
    width_over_mass: float
    coupling_pattern: str | None
    topology: str | None

    @property
    def value(self) -> float:
        return self.anchor.value

    @property
    def value_gev(self) -> float:
        return self.anchor.value

    @property
    def value_tev(self) -> float:
        return float(self.anchor.value / 1000.0)

    @property
    def budget(self) -> float:
        return self.value_tev

    @property
    def units(self) -> str | None:
        return self.anchor.units

    @property
    def source(self) -> str | None:
        return self.anchor.source

    @property
    def source_url(self) -> str | None:
        return self.anchor.source_url

    @property
    def snapshot_path(self) -> str | None:
        return self.anchor.snapshot_path


@dataclass(frozen=True)
class CR014CrossSectionAnchor:
    """Typed diagnostic four-top cross-section measurement."""

    anchor: Anchor
    value_id: str
    entry_index: int
    display: str | None
    limit_type: str | None
    experiment: str | None
    uncertainty_plus: float | None
    uncertainty_minus: float | None
    uncertainty_type: str | None
    source_key: str | None
    arxiv_url: str | None

    @property
    def value_fb(self) -> float:
        return self.anchor.value

    @property
    def units(self) -> str | None:
        return self.anchor.units

    @property
    def source(self) -> str | None:
        return self.anchor.source

    @property
    def source_url(self) -> str | None:
        return self.anchor.source_url

    @property
    def snapshot_path(self) -> str | None:
        return self.anchor.snapshot_path


@dataclass(frozen=True)
class CR014CrossSectionRange:
    """Diagnostic-only ATLAS top-philic cross-section upper-limit range."""

    value_id: str
    entry_index: int
    display: str | None
    limit_type: str | None
    cl: str | None
    experiment: str | None
    value_min_fb: float
    value_max_fb: float
    units: str | None
    source: str | None
    source_key: str | None
    source_url: str | None
    snapshot_path: str | None
    assumptions: Mapping[str, Any]


@dataclass(frozen=True)
class CR014Anchor:
    """YAML-loaded CR014 mass-limit bundle."""

    active_limit: CR014MassLimitAnchor
    all_applicable_mass_limits: tuple[CR014MassLimitAnchor, ...]
    atlas_four_top_cross_section: CR014CrossSectionAnchor
    cms_four_top_cross_section: CR014CrossSectionAnchor
    atlas_top_philic_cross_section_range: CR014CrossSectionRange
    budget_policy: str

    @property
    def value(self) -> float:
        return self.active_limit.value_tev

    @property
    def value_gev(self) -> float:
        return self.active_limit.value_gev

    @property
    def budget(self) -> float:
        return self.active_limit.budget

    @property
    def source_url(self) -> str | None:
        return self.active_limit.source_url

    @property
    def width_over_mass(self) -> float:
        return self.active_limit.width_over_mass


def _optional_str(value: Any) -> str | None:
    return None if value is None else str(value)


def _required_float(value: Any, *, process_id: str, field_name: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise AnchorError(
            f"{process_id}: CR014 field {field_name!r}={value!r} is not numeric"
        ) from exc
    if not math.isfinite(number):
        raise AnchorError(
            f"{process_id}: CR014 field {field_name!r}={value!r} is not finite"
        )
    return number


def _optional_float(value: Any, *, process_id: str, field_name: str) -> float | None:
    if value is None:
        return None
    return _required_float(value, process_id=process_id, field_name=field_name)


def _parent(process_id: str) -> Mapping[str, Any]:
    data = load_full_yaml(process_id, family=_FAMILY)
    parent = data.get(_PARENT_KEY)
    if not isinstance(parent, Mapping):
        raise AnchorError(f"{process_id}: expected mapping-shaped {_PARENT_KEY}")
    return parent


def _value_entries(process_id: str) -> Sequence[Mapping[str, Any]]:
    values = _parent(process_id).get("values")
    if not isinstance(values, list) or not values:
        raise AnchorError(f"{process_id}: expected non-empty {_PARENT_KEY}.values")
    for index, entry in enumerate(values):
        if not isinstance(entry, Mapping):
            raise AnchorError(
                f"{process_id}: {_PARENT_KEY}.values[{index}] is not a mapping"
            )
    return values


def _entry_by_value_id(
    process_id: str,
    value_id: str,
) -> tuple[int, Mapping[str, Any]]:
    matches: list[tuple[int, Mapping[str, Any]]] = []
    for index, entry in enumerate(_value_entries(process_id)):
        if entry.get("value_id") == value_id:
            matches.append((index, entry))
    if not matches:
        present = [str(entry.get("value_id")) for entry in _value_entries(process_id)]
        raise AnchorError(
            f"{process_id}: value_id {value_id!r} not found in "
            f"{_PARENT_KEY}.values (present: {present})"
        )
    if len(matches) > 1:
        raise AnchorError(f"{process_id}: duplicate value_id {value_id!r}")
    return matches[0]


def _assumptions(
    entry: Mapping[str, Any],
    *,
    process_id: str,
    value_id: str,
    required: bool,
) -> Mapping[str, Any]:
    assumptions = entry.get("assumptions")
    if assumptions is None and not required:
        return {}
    if not isinstance(assumptions, Mapping):
        raise AnchorError(f"{process_id}: {value_id}.assumptions must be a mapping")
    return assumptions


def _route_value_through_scaffold(
    entry: Mapping[str, Any],
    *,
    process_id: str,
    entry_index: int,
    value: float,
) -> Anchor:
    block_key = f"{_PARENT_KEY}.values[{entry_index}]"
    virtual_entry = dict(entry)
    virtual_entry["value"] = value
    virtual_block = {block_key: virtual_entry}
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
        scaffold_anchor = load_anchor(
            process_id,
            family=_FAMILY,
            candidates=(block_key,),
            uncertainty_key=_SCAFFOLD_UNCERTAINTY_KEY,
        )
    finally:
        anchor_scaffold.load_pdg_block = original_load_pdg_block

    if scaffold_anchor.block_key != block_key:
        raise AnchorError(
            f"{process_id}: load_anchor selected {scaffold_anchor.block_key!r}, "
            f"expected {block_key!r}"
        )
    return scaffold_anchor


def _load_mass_limit_anchor(
    value_id: str,
    *,
    process_id: str,
) -> CR014MassLimitAnchor:
    index, entry = _entry_by_value_id(process_id, value_id)
    value = _required_float(
        entry.get("value"),
        process_id=process_id,
        field_name=f"{value_id}.value",
    )
    if value <= 0.0:
        raise AnchorError(f"{process_id}: {value_id}.value must be positive")
    if entry.get("units") != _MASS_UNITS_YAML:
        raise AnchorError(
            f"{process_id}: expected units {_MASS_UNITS_YAML!r} for "
            f"{value_id}, got {entry.get('units')!r}"
        )
    if entry.get("limit_type") != "mass_excluded_up_to":
        raise AnchorError(f"{process_id}: {value_id} is not a mass-exclusion row")

    assumptions = _assumptions(entry, process_id=process_id, value_id=value_id, required=True)
    width_over_mass = _required_float(
        assumptions.get("width_over_mass"),
        process_id=process_id,
        field_name=f"{value_id}.assumptions.width_over_mass",
    )
    if width_over_mass <= 0.0:
        raise AnchorError(f"{process_id}: {value_id}.width_over_mass must be positive")

    scaffold_anchor = _route_value_through_scaffold(
        entry,
        process_id=process_id,
        entry_index=index,
        value=value,
    )

    return CR014MassLimitAnchor(
        anchor=scaffold_anchor,
        value_id=value_id,
        entry_index=index,
        display=_optional_str(entry.get("display")),
        limit_type=_optional_str(entry.get("limit_type")),
        cl=_optional_str(entry.get("cl")),
        experiment=_optional_str(entry.get("experiment")),
        source_key=_optional_str(entry.get("source_key")),
        access_date=_optional_str(entry.get("access_date")),
        sha256=_optional_str(entry.get("sha256")),
        resonance_type=_optional_str(assumptions.get("resonance_type")),
        width_over_mass=float(width_over_mass),
        coupling_pattern=_optional_str(assumptions.get("coupling_pattern")),
        topology=_optional_str(assumptions.get("topology")),
    )


def _load_cross_section_anchor(
    value_id: str,
    *,
    process_id: str,
) -> CR014CrossSectionAnchor:
    index, entry = _entry_by_value_id(process_id, value_id)
    value = _required_float(
        entry.get("value"),
        process_id=process_id,
        field_name=f"{value_id}.value",
    )
    if value <= 0.0:
        raise AnchorError(f"{process_id}: {value_id}.value must be positive")
    if entry.get("units") != _CROSS_SECTION_UNITS:
        raise AnchorError(
            f"{process_id}: expected units {_CROSS_SECTION_UNITS!r} for "
            f"{value_id}, got {entry.get('units')!r}"
        )
    scaffold_anchor = _route_value_through_scaffold(
        entry,
        process_id=process_id,
        entry_index=index,
        value=value,
    )
    return CR014CrossSectionAnchor(
        anchor=scaffold_anchor,
        value_id=value_id,
        entry_index=index,
        display=_optional_str(entry.get("display")),
        limit_type=_optional_str(entry.get("limit_type")),
        experiment=_optional_str(entry.get("experiment")),
        uncertainty_plus=_optional_float(
            entry.get("uncertainty_plus"),
            process_id=process_id,
            field_name=f"{value_id}.uncertainty_plus",
        ),
        uncertainty_minus=_optional_float(
            entry.get("uncertainty_minus"),
            process_id=process_id,
            field_name=f"{value_id}.uncertainty_minus",
        ),
        uncertainty_type=_optional_str(entry.get("uncertainty_type")),
        source_key=_optional_str(entry.get("source_key")),
        arxiv_url=_optional_str(entry.get("arxiv_url")),
    )


def _load_cross_section_range(
    value_id: str,
    *,
    process_id: str,
) -> CR014CrossSectionRange:
    index, entry = _entry_by_value_id(process_id, value_id)
    if entry.get("units") != _CROSS_SECTION_UNITS:
        raise AnchorError(
            f"{process_id}: expected units {_CROSS_SECTION_UNITS!r} for "
            f"{value_id}, got {entry.get('units')!r}"
        )
    value_min = _required_float(
        entry.get("value_min"),
        process_id=process_id,
        field_name=f"{value_id}.value_min",
    )
    value_max = _required_float(
        entry.get("value_max"),
        process_id=process_id,
        field_name=f"{value_id}.value_max",
    )
    if value_min <= 0.0 or value_max < value_min:
        raise AnchorError(
            f"{process_id}: {value_id} must satisfy 0 < value_min <= value_max"
        )
    return CR014CrossSectionRange(
        value_id=value_id,
        entry_index=index,
        display=_optional_str(entry.get("display")),
        limit_type=_optional_str(entry.get("limit_type")),
        cl=_optional_str(entry.get("cl")),
        experiment=_optional_str(entry.get("experiment")),
        value_min_fb=float(value_min),
        value_max_fb=float(value_max),
        units=_optional_str(entry.get("units")),
        source=_optional_str(entry.get("source")),
        source_key=_optional_str(entry.get("source_key")),
        source_url=_optional_str(entry.get("source_url")),
        snapshot_path=_optional_str(entry.get("snapshot_path")),
        assumptions=dict(
            _assumptions(
                entry,
                process_id=process_id,
                value_id=value_id,
                required=False,
            )
        ),
    )


def _applicable_mass_exclusion_value_ids(process_id: str) -> tuple[str, ...]:
    value_ids: list[str] = []
    for entry in _value_entries(process_id):
        if entry.get("limit_type") != "mass_excluded_up_to":
            continue
        value_id = _optional_str(entry.get("value_id"))
        if value_id is None:
            raise AnchorError(f"{process_id}: mass-exclusion row missing value_id")
        assumptions = _assumptions(
            entry,
            process_id=process_id,
            value_id=value_id,
            required=True,
        )
        resonance_type = _optional_str(assumptions.get("resonance_type"))
        width_over_mass = _required_float(
            assumptions.get("width_over_mass"),
            process_id=process_id,
            field_name=f"{value_id}.assumptions.width_over_mass",
        )
        if resonance_type == "vector mediator" and math.isclose(
            width_over_mass,
            _REFERENCE_WIDTH_OVER_MASS,
            rel_tol=0.0,
            abs_tol=1.0e-12,
        ):
            value_ids.append(value_id)
    if not value_ids:
        raise AnchorError(
            f"{process_id}: no vector mediator mass-exclusion row found at "
            f"width_over_mass={_REFERENCE_WIDTH_OVER_MASS}"
        )
    return tuple(value_ids)


def _load_cr014_anchor(process_id: str) -> CR014Anchor:
    applicable_limits = tuple(
        _load_mass_limit_anchor(value_id, process_id=process_id)
        for value_id in _applicable_mass_exclusion_value_ids(process_id)
    )
    strongest = max(applicable_limits, key=lambda limit: limit.value_gev)
    if strongest.value_id != _ACTIVE_CMS_VECTOR_50PCT:
        raise AnchorError(
            f"{process_id}: active CR014 limit {_ACTIVE_CMS_VECTOR_50PCT!r} is not "
            f"the strongest applicable row; strongest is {strongest.value_id!r}"
        )
    if strongest.value_tev <= 0.0:
        raise AnchorError(f"{process_id}: CR014 active budget must be positive")

    return CR014Anchor(
        active_limit=strongest,
        all_applicable_mass_limits=applicable_limits,
        atlas_four_top_cross_section=_load_cross_section_anchor(
            _ATLAS_4TOP_CROSS_SECTION,
            process_id=process_id,
        ),
        cms_four_top_cross_section=_load_cross_section_anchor(
            _CMS_4TOP_CROSS_SECTION,
            process_id=process_id,
        ),
        atlas_top_philic_cross_section_range=_load_cross_section_range(
            _ATLAS_TOP_PHILIC_XSEC_RANGE,
            process_id=process_id,
        ),
        budget_policy=_BUDGET_POLICY,
    )


def _limit_from_anchor(
    anchor: CR014MassLimitAnchor,
    *,
    process_id: str,
) -> ColliderResonanceLimit:
    benchmark = (
        f"{anchor.resonance_type}, Gamma/m={anchor.width_over_mass:g}, "
        "top-philic simplified interpretation"
    )
    return ColliderResonanceLimit(
        process_id=process_id,
        resonance=_RESONANCE,
        final_state=_FINAL_STATE,
        limit_kind=MASS_LOWER_BOUND,
        value=float(anchor.value_tev),
        units=_MASS_UNITS_COMPARISON,
        cl=anchor.cl,
        source=anchor.source,
        source_url=anchor.source_url,
        limit_type=anchor.limit_type,
        benchmark_model=benchmark,
        diagnostics={
            "value_id": anchor.value_id,
            "entry_index": anchor.entry_index,
            "display": anchor.display,
            "experiment": anchor.experiment,
            "source_key": anchor.source_key,
            "snapshot_path": anchor.snapshot_path,
            "access_date": anchor.access_date,
            "sha256": anchor.sha256,
            "yaml_units": anchor.units,
            "yaml_observable": anchor.anchor.observable,
            "yaml_value_gev": float(anchor.value_gev),
            "width_over_mass": float(anchor.width_over_mass),
            "resonance_type": anchor.resonance_type,
            "coupling_pattern": anchor.coupling_pattern,
            "topology": anchor.topology,
        },
    )


@register
class Constraint:
    """Catalogued CMS top-philic vector four-top mass-exclusion constraint."""

    process_id = "CR014"
    severity = Severity.HARD
    observable = "m(top-philic vector Z' -> t tbar t tbar)"

    def __init__(self) -> None:
        self.anchor = _load_cr014_anchor(self.process_id)
        self.limit = _limit_from_anchor(
            self.anchor.active_limit,
            process_id=self.process_id,
        )

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        try:
            mass_gev, mass_source = resolve_top_philic_vector_mass_gev(
                ew_mass_extra=point.get_extra(_MASS_EXTRA),
                couplings=point.get_extra(_COUPLINGS_EXTRA),
            )
        except (TypeError, ValueError, AttributeError) as exc:
            return ConstraintResult(
                process_id=self.process_id,
                severity=self.severity,
                passes=False,
                experimental=float(self.anchor.value),
                budget=float(self.anchor.budget),
                notes=f"invalid top-philic vector mass proxy input; CR014 not evaluated: {exc}",
                diagnostics={
                    "exception_type": type(exc).__name__,
                    "exception": str(exc),
                    "required_mass_sources": (
                        _MASS_EXTRA,
                        f"{_COUPLINGS_EXTRA}.M_KK",
                    ),
                    "active_value_id": self.anchor.active_limit.value_id,
                    "active_width_over_mass": float(self.anchor.width_over_mass),
                    "needs_human_physics": (
                        TOP_PHILIC_VECTOR_FOUR_TOP_PROXY_ASSUMPTION_V1
                    ),
                    "sigma_times_br_recast_status": "NEEDS-HUMAN-PHYSICS",
                    "width_dependence_recast_status": "NEEDS-HUMAN-PHYSICS",
                },
            )

        if mass_gev is None:
            return ConstraintResult(
                process_id=self.process_id,
                severity=self.severity,
                passes=True,
                experimental=float(self.anchor.value),
                budget=float(self.anchor.budget),
                notes=(
                    f"extras {_MASS_EXTRA!r} and {_COUPLINGS_EXTRA!r} absent; "
                    "top-philic vector four-top constraint was not evaluated."
                ),
                diagnostics={
                    "missing_extras": (_MASS_EXTRA, _COUPLINGS_EXTRA),
                    "required_mass_sources": (
                        _MASS_EXTRA,
                        f"{_COUPLINGS_EXTRA}.M_KK",
                    ),
                    "budget_policy": self.anchor.budget_policy,
                    "active_value_id": self.anchor.active_limit.value_id,
                    "active_limit_gev": float(self.anchor.value_gev),
                    "active_width_over_mass": float(self.anchor.width_over_mass),
                    "needs_human_physics": (
                        TOP_PHILIC_VECTOR_FOUR_TOP_PROXY_ASSUMPTION_V1
                    ),
                    "sigma_times_br_recast_status": "NEEDS-HUMAN-PHYSICS",
                    "width_dependence_recast_status": "NEEDS-HUMAN-PHYSICS",
                    "four_top_acceptance_recast_status": "NEEDS-HUMAN-PHYSICS",
                    "sm_four_top_background_recast_status": "NEEDS-HUMAN-PHYSICS",
                },
            )

        prediction = top_philic_vector_four_top_prediction_from_mass_gev(
            float(mass_gev),
            resonance=_RESONANCE,
            final_state=_FINAL_STATE,
            width_over_mass=float(self.anchor.width_over_mass),
        )
        comparison = evaluate_collider_resonance_limit(prediction, self.limit)
        diagnostics = dict(comparison.diagnostics)
        diagnostics.update(
            {
                "input_mass_gev": float(mass_gev),
                "m_top_philic_vector_proxy_gev": float(mass_gev),
                "m_top_philic_vector_proxy_tev": float(comparison.predicted_mass_tev),
                "mass_proxy": "m_Zprime_top_philic = kk_ew_mass_gev or M_KK",
                "mass_source": mass_source,
                "active_value_id": self.anchor.active_limit.value_id,
                "active_experiment": self.anchor.active_limit.experiment,
                "active_display": self.anchor.active_limit.display,
                "active_limit_gev": float(self.anchor.value_gev),
                "active_width_over_mass": float(self.anchor.width_over_mass),
                "active_coupling_pattern": self.anchor.active_limit.coupling_pattern,
                "active_topology": self.anchor.active_limit.topology,
                "budget_policy": self.anchor.budget_policy,
                "all_applicable_mass_limits_gev": {
                    limit.value_id: float(limit.value_gev)
                    for limit in self.anchor.all_applicable_mass_limits
                },
                "atlas_four_top_cross_section_fb": float(
                    self.anchor.atlas_four_top_cross_section.value_fb
                ),
                "atlas_four_top_cross_section_uncertainty_plus_fb": (
                    self.anchor.atlas_four_top_cross_section.uncertainty_plus
                ),
                "atlas_four_top_cross_section_uncertainty_minus_fb": (
                    self.anchor.atlas_four_top_cross_section.uncertainty_minus
                ),
                "cms_four_top_cross_section_fb": float(
                    self.anchor.cms_four_top_cross_section.value_fb
                ),
                "cms_four_top_cross_section_uncertainty_plus_fb": (
                    self.anchor.cms_four_top_cross_section.uncertainty_plus
                ),
                "cms_four_top_cross_section_uncertainty_minus_fb": (
                    self.anchor.cms_four_top_cross_section.uncertainty_minus
                ),
                "atlas_top_philic_cross_section_limit_range_fb": (
                    self.anchor.atlas_top_philic_cross_section_range.value_min_fb,
                    self.anchor.atlas_top_philic_cross_section_range.value_max_fb,
                ),
                "low_mass_benchmark_warning": (
                    "The active 0.85 TeV observed limit is a low-mass, "
                    "width-dependent CMS simplified-model benchmark."
                ),
                "needs_human_physics": (
                    TOP_PHILIC_VECTOR_FOUR_TOP_PROXY_ASSUMPTION_V1
                ),
                "sigma_times_br_recast_status": "NEEDS-HUMAN-PHYSICS",
                "width_dependence_recast_status": "NEEDS-HUMAN-PHYSICS",
                "top_philic_coupling_recast_status": "NEEDS-HUMAN-PHYSICS",
                "four_top_acceptance_recast_status": "NEEDS-HUMAN-PHYSICS",
                "sm_four_top_background_recast_status": "NEEDS-HUMAN-PHYSICS",
                "sigma_times_br_proxy_fb": comparison.predicted_sigma_times_br,
            }
        )

        return ConstraintResult(
            process_id=self.process_id,
            severity=self.severity,
            passes=bool(comparison.passes),
            predicted=float(comparison.predicted_mass_tev),
            experimental=float(comparison.experimental_limit),
            ratio=float(comparison.ratio_to_budget),
            budget=float(comparison.budget),
            notes=(
                "Four-top top-philic vector mediator uses the documented mass "
                "proxy m_Zprime_top_philic = kk_ew_mass_gev or M_KK and the "
                "strongest applicable observed CMS Gamma/m=50% mass-exclusion "
                "edge from CR014.yaml. HARD ratio is m_limit/m_proxy; "
                "sigma*BR, width, coupling, acceptance, and SM-four-top "
                "background recasts are marked NEEDS-HUMAN-PHYSICS."
            ),
            diagnostics=diagnostics,
        )
