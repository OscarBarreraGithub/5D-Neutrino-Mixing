"""CR012 - spin-1 diboson high-mass resonance limit.

Physics
-------
The observable is the first spin-1 electroweak/custodial KK-vector mass in

    pp -> V_KK^(1) -> WW, WZ, ZZ

compared with the catalogued HVT model-B diboson exclusions in
``flavor_catalog/processes/collider_rs/CR012.yaml``.  The reusable
collider-resonance mass-limit comparison machinery lives in
``quarkConstraints.collider_resonance`` and is reached only through the
``flavor_catalog_constraints.physics_adapters.collider_resonance`` boundary.

RS matching status
------------------
NEEDS-HUMAN-PHYSICS.  A faithful diboson recast needs
``sigma(pp -> W'/Z'/V_KK) * BR(V -> WW/WZ/ZZ)``, the VV branching surface,
total width, production-mode mixture, interference, acceptance, and the
experiment's mass-dependent limit curve.  The current ``ParameterPoint`` may
carry an EW KK mass but not that collider signal model, so CR012 uses the
documented mass-exclusion proxy ``m_spin1_diboson = kk_ew_mass_gev or M_KK``.

Severity
--------
HARD.  The active budget is the observed 95% CL HVT model-B ``W' -> WZ`` mass
lower limit loaded from ``CR012.yaml``.  The mass-degenerate ``V' -> VV``
diboson number, the ``VV+VH`` combined number, and the ``Z' -> WW`` excluded
intervals are retained as diagnostics because they require additional
benchmark assumptions beyond the single charged-vector mass proxy.  The ratio
is ``m_limit / m_spin1_diboson`` so ``ratio <= 1`` passes.

Catalog sidecar
---------------
``CR012.yaml`` stores numerical limits in a ``pdg_or_equivalent.values`` list.
This module adapts selected numeric list entries into the scaffold
``load_anchor`` path and fails loudly if the expected value IDs or units are
missing.
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
    KK_DIBOSON_SPIN1_MASS_PROXY_ASSUMPTION_V1,
    MASS_LOWER_BOUND,
    ColliderResonanceLimit,
    evaluate_collider_resonance_limit,
    kk_diboson_spin1_prediction_from_m_kk_gev,
    resolve_kk_diboson_spin1_mass_gev,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "collider_rs"
_MASS_EXTRA = "kk_ew_mass_gev"
_COUPLINGS_EXTRA = "quark_mass_basis_couplings"
_PARENT_KEY = "pdg_or_equivalent"
_SCAFFOLD_UNCERTAINTY_KEY = "__cr012_uncertainty_is_not_used__"
_MASS_UNITS = "TeV"
_RESONANCE = "V_KK^(1) spin-1"
_FINAL_STATE = "WW/WZ/ZZ"

_ACTIVE_WPRIME_WZ = "PDG2025:CR012:HVTB_Wprime_WZ_mass_lower"
_MASS_DEGENERATE_VV = "CMS2023:CR012:HVTB_Vprime_VV_mass_lower"
_VV_VH_COMBINED = "CMS2023:CR012:HVTB_Vprime_VV_VH_mass_lower"
_ZPRIME_WW_INTERVALS = "CMS2023:CR012:HVTB_Zprime_WW_excluded_intervals"
_NUMERIC_MASS_LIMIT_IDS = (
    _ACTIVE_WPRIME_WZ,
    _MASS_DEGENERATE_VV,
    _VV_VH_COMBINED,
)
_BUDGET_POLICY = (
    "active HARD budget is the observed 95% CL HVT model-B W'->WZ mass "
    "lower limit in CR012.yaml. The mass-degenerate V'->VV limit is kept as "
    "a stronger pure-diboson diagnostic because it assumes a degenerate HVT "
    "triplet, while the VV+VH combined number is not used as a pure-diboson "
    "threshold."
)


@dataclass(frozen=True)
class CR012ValueAnchor:
    """Typed CR012 numeric value entry routed through scaffold ``load_anchor``."""

    anchor: Anchor
    value_id: str
    entry_index: int
    display: str | None
    limit_type: str | None
    cl: str | None
    confidence_level: float | None
    experiment: str | None
    model: str | None
    production_mode: str | None
    decay_mode: str | None
    source_key: str | None
    access_date: str | None
    sha256: str | None

    @property
    def value(self) -> float:
        return self.anchor.value

    @property
    def value_tev(self) -> float:
        return self.anchor.value

    @property
    def value_gev(self) -> float:
        return float(self.anchor.value * 1000.0)

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
class CR012IntervalAnchor:
    """Typed diagnostic-only CR012 excluded-interval entry."""

    value_id: str
    entry_index: int
    display: str | None
    limit_type: str | None
    cl: str | None
    confidence_level: float | None
    experiment: str | None
    model: str | None
    production_mode: str | None
    decay_mode: str | None
    excluded_ranges_tev: tuple[tuple[float, float], ...]
    expected_limit_tev: float | None
    expected_limit_units: str | None
    source: str | None
    source_key: str | None
    source_url: str | None
    snapshot_path: str | None
    access_date: str | None
    sha256: str | None


@dataclass(frozen=True)
class CR012Anchor:
    """YAML-loaded CR012 diboson mass-limit bundle."""

    active_limit: CR012ValueAnchor
    mass_degenerate_vv: CR012ValueAnchor
    vv_vh_combined: CR012ValueAnchor
    zprime_ww_intervals: CR012IntervalAnchor
    all_mass_limits: tuple[CR012ValueAnchor, ...]
    parent_canonical_source: str | None
    parent_source_key: str | None
    parent_source_url: str | None
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
    def strongest_pure_diboson_limit(self) -> CR012ValueAnchor:
        return max(
            (self.active_limit, self.mass_degenerate_vv),
            key=lambda limit: limit.value_tev,
        )


def _optional_str(value: Any) -> str | None:
    return None if value is None else str(value)


def _optional_float(value: Any, *, process_id: str, field_name: str) -> float | None:
    if value is None:
        return None
    return _required_float(value, process_id=process_id, field_name=field_name)


def _required_float(value: Any, *, process_id: str, field_name: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise AnchorError(
            f"{process_id}: CR012 field {field_name!r}={value!r} is not numeric"
        ) from exc
    if not math.isfinite(number):
        raise AnchorError(
            f"{process_id}: CR012 field {field_name!r}={value!r} is not finite"
        )
    return number


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


def _load_mass_limit_anchor(value_id: str, *, process_id: str) -> CR012ValueAnchor:
    index, entry = _entry_by_value_id(process_id, value_id)
    value = _required_float(
        entry.get("value"),
        process_id=process_id,
        field_name=f"{value_id}.value",
    )
    if value <= 0.0:
        raise AnchorError(f"{process_id}: {value_id}.value must be positive")
    if entry.get("units") != _MASS_UNITS:
        raise AnchorError(
            f"{process_id}: expected units {_MASS_UNITS!r} for "
            f"{value_id}, got {entry.get('units')!r}"
        )

    block_key = f"{_PARENT_KEY}.values[{index}]"
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
            f"expected {block_key!r} for CR012 value_id {value_id!r}"
        )

    return CR012ValueAnchor(
        anchor=scaffold_anchor,
        value_id=value_id,
        entry_index=index,
        display=_optional_str(entry.get("display")),
        limit_type=_optional_str(entry.get("limit_type")),
        cl=_optional_str(entry.get("cl")),
        confidence_level=_optional_float(
            entry.get("confidence_level"),
            process_id=process_id,
            field_name=f"{value_id}.confidence_level",
        ),
        experiment=_optional_str(entry.get("experiment")),
        model=_optional_str(entry.get("model")),
        production_mode=_optional_str(entry.get("production_mode")),
        decay_mode=_optional_str(entry.get("decay_mode")),
        source_key=_optional_str(entry.get("source_key")),
        access_date=_optional_str(entry.get("access_date")),
        sha256=_optional_str(entry.get("sha256")),
    )


def _load_interval_anchor(value_id: str, *, process_id: str) -> CR012IntervalAnchor:
    index, entry = _entry_by_value_id(process_id, value_id)
    if entry.get("units") != _MASS_UNITS:
        raise AnchorError(
            f"{process_id}: expected units {_MASS_UNITS!r} for "
            f"{value_id}, got {entry.get('units')!r}"
        )
    ranges_raw = entry.get("excluded_ranges_TeV")
    if not isinstance(ranges_raw, list) or not ranges_raw:
        raise AnchorError(f"{process_id}: {value_id}.excluded_ranges_TeV is required")

    ranges: list[tuple[float, float]] = []
    for range_index, item in enumerate(ranges_raw):
        if not isinstance(item, (list, tuple)) or len(item) != 2:
            raise AnchorError(
                f"{process_id}: {value_id}.excluded_ranges_TeV[{range_index}] "
                "must be a two-element sequence"
            )
        low = _required_float(
            item[0],
            process_id=process_id,
            field_name=f"{value_id}.excluded_ranges_TeV[{range_index}][0]",
        )
        high = _required_float(
            item[1],
            process_id=process_id,
            field_name=f"{value_id}.excluded_ranges_TeV[{range_index}][1]",
        )
        if low <= 0.0 or high <= low:
            raise AnchorError(
                f"{process_id}: {value_id}.excluded_ranges_TeV[{range_index}] "
                f"must satisfy 0 < low < high, got {(low, high)!r}"
            )
        ranges.append((float(low), float(high)))

    expected_units = _optional_str(entry.get("expected_limit_units"))
    if expected_units is not None and expected_units != _MASS_UNITS:
        raise AnchorError(
            f"{process_id}: expected_limit_units for {value_id} must be "
            f"{_MASS_UNITS!r}, got {expected_units!r}"
        )

    return CR012IntervalAnchor(
        value_id=value_id,
        entry_index=index,
        display=_optional_str(entry.get("display")),
        limit_type=_optional_str(entry.get("limit_type")),
        cl=_optional_str(entry.get("cl")),
        confidence_level=_optional_float(
            entry.get("confidence_level"),
            process_id=process_id,
            field_name=f"{value_id}.confidence_level",
        ),
        experiment=_optional_str(entry.get("experiment")),
        model=_optional_str(entry.get("model")),
        production_mode=_optional_str(entry.get("production_mode")),
        decay_mode=_optional_str(entry.get("decay_mode")),
        excluded_ranges_tev=tuple(ranges),
        expected_limit_tev=_optional_float(
            entry.get("expected_limit"),
            process_id=process_id,
            field_name=f"{value_id}.expected_limit",
        ),
        expected_limit_units=expected_units,
        source=_optional_str(entry.get("source")),
        source_key=_optional_str(entry.get("source_key")),
        source_url=_optional_str(entry.get("source_url")),
        snapshot_path=_optional_str(entry.get("snapshot_path")),
        access_date=_optional_str(entry.get("access_date")),
        sha256=_optional_str(entry.get("sha256")),
    )


def _load_cr012_anchor(process_id: str) -> CR012Anchor:
    limits = tuple(
        _load_mass_limit_anchor(value_id, process_id=process_id)
        for value_id in _NUMERIC_MASS_LIMIT_IDS
    )
    by_id = {limit.value_id: limit for limit in limits}
    parent = _parent(process_id)
    anchor = CR012Anchor(
        active_limit=by_id[_ACTIVE_WPRIME_WZ],
        mass_degenerate_vv=by_id[_MASS_DEGENERATE_VV],
        vv_vh_combined=by_id[_VV_VH_COMBINED],
        zprime_ww_intervals=_load_interval_anchor(
            _ZPRIME_WW_INTERVALS,
            process_id=process_id,
        ),
        all_mass_limits=limits,
        parent_canonical_source=_optional_str(parent.get("canonical_source")),
        parent_source_key=_optional_str(parent.get("source_key")),
        parent_source_url=_optional_str(parent.get("source_url")),
        budget_policy=_BUDGET_POLICY,
    )
    if anchor.budget <= 0.0:
        raise AnchorError(f"{process_id}: CR012 active budget must be positive")
    return anchor


def _limit_from_anchor(
    anchor: CR012ValueAnchor,
    *,
    process_id: str,
) -> ColliderResonanceLimit:
    return ColliderResonanceLimit(
        process_id=process_id,
        resonance=_RESONANCE,
        final_state=_FINAL_STATE,
        limit_kind=MASS_LOWER_BOUND,
        value=float(anchor.value_tev),
        units=_MASS_UNITS,
        cl=anchor.cl,
        source=anchor.source,
        source_url=anchor.source_url,
        limit_type=anchor.limit_type,
        benchmark_model=anchor.model or anchor.anchor.observable,
        diagnostics={
            "value_id": anchor.value_id,
            "entry_index": anchor.entry_index,
            "display": anchor.display,
            "experiment": anchor.experiment,
            "confidence_level": anchor.confidence_level,
            "model": anchor.model,
            "production_mode": anchor.production_mode,
            "decay_mode": anchor.decay_mode,
            "source_key": anchor.source_key,
            "snapshot_path": anchor.snapshot_path,
            "access_date": anchor.access_date,
            "sha256": anchor.sha256,
            "yaml_units": anchor.units,
            "yaml_observable": anchor.anchor.observable,
        },
    )


@register
class Constraint:
    """Catalogued spin-1 diboson high-mass resonance constraint."""

    process_id = "CR012"
    severity = Severity.HARD
    observable = "m(V_KK^(1) spin-1 -> WW/WZ/ZZ)"

    def __init__(self) -> None:
        self.anchor = _load_cr012_anchor(self.process_id)
        self.limit = _limit_from_anchor(
            self.anchor.active_limit,
            process_id=self.process_id,
        )

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        try:
            m_kk_gev, mass_source = resolve_kk_diboson_spin1_mass_gev(
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
                notes=f"invalid spin-1 diboson mass proxy input; CR012 not evaluated: {exc}",
                diagnostics={
                    "exception_type": type(exc).__name__,
                    "exception": str(exc),
                    "required_mass_sources": (
                        _MASS_EXTRA,
                        f"{_COUPLINGS_EXTRA}.M_KK",
                    ),
                    "needs_human_physics": (
                        KK_DIBOSON_SPIN1_MASS_PROXY_ASSUMPTION_V1
                    ),
                    "sigma_times_br_recast_status": "NEEDS-HUMAN-PHYSICS",
                },
            )

        if m_kk_gev is None:
            return ConstraintResult(
                process_id=self.process_id,
                severity=self.severity,
                passes=True,
                experimental=float(self.anchor.value),
                budget=float(self.anchor.budget),
                notes=(
                    f"extras {_MASS_EXTRA!r} and {_COUPLINGS_EXTRA!r} absent; "
                    "spin-1 diboson resonance constraint was not evaluated."
                ),
                diagnostics={
                    "missing_extras": (_MASS_EXTRA, _COUPLINGS_EXTRA),
                    "required_mass_sources": (
                        _MASS_EXTRA,
                        f"{_COUPLINGS_EXTRA}.M_KK",
                    ),
                    "budget_policy": self.anchor.budget_policy,
                    "active_value_id": self.anchor.active_limit.value_id,
                    "needs_human_physics": (
                        KK_DIBOSON_SPIN1_MASS_PROXY_ASSUMPTION_V1
                    ),
                    "sigma_times_br_recast_status": "NEEDS-HUMAN-PHYSICS",
                    "branching_surface_recast_status": "NEEDS-HUMAN-PHYSICS",
                    "width_acceptance_recast_status": "NEEDS-HUMAN-PHYSICS",
                },
            )

        prediction = kk_diboson_spin1_prediction_from_m_kk_gev(
            float(m_kk_gev),
            resonance=_RESONANCE,
            final_state=_FINAL_STATE,
        )
        comparison = evaluate_collider_resonance_limit(prediction, self.limit)
        diagnostics = dict(comparison.diagnostics)
        diagnostics.update(
            {
                "m_kk_gev": float(m_kk_gev),
                "m_spin1_diboson_proxy_gev": float(m_kk_gev),
                "m_spin1_diboson_proxy_tev": float(comparison.predicted_mass_tev),
                "mass_proxy": "m_spin1_diboson = kk_ew_mass_gev or M_KK",
                "mass_source": mass_source,
                "active_value_id": self.anchor.active_limit.value_id,
                "active_experiment": self.anchor.active_limit.experiment,
                "active_display": self.anchor.active_limit.display,
                "active_model": self.anchor.active_limit.model,
                "active_decay_mode": self.anchor.active_limit.decay_mode,
                "active_confidence_level": self.anchor.active_limit.confidence_level,
                "parent_source_key": self.anchor.parent_source_key,
                "parent_canonical_source": self.anchor.parent_canonical_source,
                "budget_policy": self.anchor.budget_policy,
                "all_mass_limits_tev": {
                    limit.value_id: float(limit.value_tev)
                    for limit in self.anchor.all_mass_limits
                },
                "pure_diboson_mass_limits_tev": {
                    self.anchor.active_limit.value_id: float(
                        self.anchor.active_limit.value_tev
                    ),
                    self.anchor.mass_degenerate_vv.value_id: float(
                        self.anchor.mass_degenerate_vv.value_tev
                    ),
                },
                "strongest_pure_diboson_diagnostic_value_id": (
                    self.anchor.strongest_pure_diboson_limit.value_id
                ),
                "strongest_pure_diboson_diagnostic_tev": float(
                    self.anchor.strongest_pure_diboson_limit.value_tev
                ),
                "nonactive_vv_vh_combined_limit_tev": float(
                    self.anchor.vv_vh_combined.value_tev
                ),
                "zprime_ww_excluded_ranges_tev": (
                    self.anchor.zprime_ww_intervals.excluded_ranges_tev
                ),
                "zprime_ww_expected_limit_tev": (
                    self.anchor.zprime_ww_intervals.expected_limit_tev
                ),
                "zprime_ww_limit_type": self.anchor.zprime_ww_intervals.limit_type,
                "color_octet_diboson_applicability": (
                    "No direct KK-gluon constraint is inferred without a model "
                    "supplying appreciable diboson decays."
                ),
                "needs_human_physics": (
                    KK_DIBOSON_SPIN1_MASS_PROXY_ASSUMPTION_V1
                ),
                "sigma_times_br_recast_status": "NEEDS-HUMAN-PHYSICS",
                "branching_surface_recast_status": "NEEDS-HUMAN-PHYSICS",
                "width_acceptance_recast_status": "NEEDS-HUMAN-PHYSICS",
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
                "Spin-1 diboson resonance uses the documented mass proxy "
                "m_spin1_diboson = kk_ew_mass_gev or M_KK and the active "
                "observed HVT model-B W'->WZ mass lower bound. HARD ratio is "
                "m_limit/m_spin1_diboson; sigma*BR, branching-surface, width, "
                "and acceptance recasts are marked NEEDS-HUMAN-PHYSICS."
            ),
            diagnostics=diagnostics,
        )
