"""CR013 - diphoton high-mass resonance limit.

Physics
-------
The observable is the first spin-2 Randall-Sundrum KK-graviton mass in

    pp -> G_KK^(1) -> gamma gamma

compared with the strongest applicable catalogued RS-graviton diphoton
mass lower limit in ``flavor_catalog/processes/collider_rs/CR013.yaml``:
the CMS/PDG 2025 row for ``ktilde = k/M_Pl = 0.1``.

RS matching status
------------------
NEEDS-HUMAN-PHYSICS.  A faithful diphoton recast needs
``sigma(pp -> G_KK) * BR(G_KK -> gamma gamma)``, the branching surface,
``ktilde = k/M_Pl``, total width, acceptance, interference, and the
experiment's mass-dependent limit curve.  The current ``ParameterPoint`` has
no dedicated graviton collider signal model, so CR013 reuses the CR007 spin-2
spectrum map ``m_GKK = x_G1 * Lambda_IR`` and compares that mass to the
YAML-loaded benchmark exclusion.  The sigma*BR / branching-surface / coupling
/ width / acceptance recast is explicitly flagged ``NEEDS-HUMAN-PHYSICS`` in
the result diagnostics.

Severity
--------
HARD.  The active budget is the observed 95% CL CMS/PDG RS-graviton
diphoton mass lower limit at ``ktilde = 0.1`` loaded from ``CR013.yaml``.
The ratio is ``m_limit / m_GKK`` so ``ratio <= 1`` passes.  Spin-0 resonance
interpretations are not evaluated because the scan point has no scalar
``X`` spectrum or couplings.

Catalog sidecar
---------------
``CR013.yaml`` stores numerical limits in a ``pdg_or_equivalent.values`` list.
This module adapts selected numeric list entries into the scaffold
``load_anchor`` path and fails loudly if the expected value IDs, units, or
benchmark-coupling metadata are missing.
"""

from __future__ import annotations

import math
import re
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
from flavor_catalog_constraints.physics_adapters.kk_graviton_resonance import (
    GAUGE_KK_ROOT_NN,
    KK_GRAVITON_DIPHOTON_RECAST_PROXY_ASSUMPTION_V1,
    KK_GRAVITON_SPIN2_SPECTRUM_ASSUMPTION_V1,
    MASS_LOWER_BOUND,
    SPIN2_GRAVITON_KK_ROOT,
    ColliderResonanceLimit,
    evaluate_collider_resonance_limit,
    kk_graviton_diphoton_prediction_from_lambda_ir_gev,
    resolve_kk_graviton_mass_mapping,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "collider_rs"
_EW_MASS_EXTRA = "kk_ew_mass_gev"
_GLUON_MASS_EXTRA = "kk_gluon_mass_gev"
_COUPLINGS_EXTRA = "quark_mass_basis_couplings"
_PARENT_KEY = "pdg_or_equivalent"
_SCAFFOLD_UNCERTAINTY_KEY = "__cr013_uncertainty_is_not_used__"
_MASS_UNITS = "TeV"
_RESONANCE = "G_KK^(1)"
_FINAL_STATE = "gamma gamma"

_ACTIVE_CMS_KTILDE_0P1 = "PDG2025:CR013:CMS2024_RSG_diphoton_kMPl_0p1"
_CMS_2024_RANGE = "CMS2024:CR013:RSG_diphoton_exclusion_range_ktilde_0p01_0p2"
_ATLAS_2021_KTILDE_0P1 = "ATLAS2021:CR013:RSG_diphoton_kMPl_0p1"
_ATLAS_2021_KTILDE_0P05 = "ATLAS2021:CR013:RSG_diphoton_kMPl_0p05"
_ATLAS_2021_KTILDE_0P01 = "ATLAS2021:CR013:RSG_diphoton_kMPl_0p01"
_CMS_2018_RANGE = "CMS2018:CR013:RSG_diphoton_exclusion_range_ktilde_0p01_0p2"
_ATLAS_2015_KTILDE_0P1 = "ATLAS2015:CR013:RSG_diphoton_kMPl_0p1"
_KNOWN_VALUE_IDS = (
    _ACTIVE_CMS_KTILDE_0P1,
    _CMS_2024_RANGE,
    _ATLAS_2021_KTILDE_0P1,
    _ATLAS_2021_KTILDE_0P05,
    _ATLAS_2021_KTILDE_0P01,
    _CMS_2018_RANGE,
    _ATLAS_2015_KTILDE_0P1,
)
_BUDGET_POLICY = (
    "active HARD budget is the strongest exact ktilde=0.1 RS-graviton "
    "diphoton mass lower limit in CR013.yaml: the CMS/PDG 2025 observed "
    "95% CL row. Ktilde-range endpoints are retained as diagnostics because "
    "using them requires the experiment's ktilde-dependent exclusion curve."
)
_KTILDE_RE = re.compile(
    r"(?:k/M_P(?:l)?|ktilde)\s*=\s*(?P<value>[0-9]+(?:\.[0-9]+)?)"
)


@dataclass(frozen=True)
class CR013ValueAnchor:
    """Typed CR013 value entry routed through the scaffold ``load_anchor``."""

    anchor: Anchor
    value_id: str
    entry_index: int
    display: str | None
    limit_type: str | None
    cl: str | None
    experiment: str | None
    analysis_year: int | None
    model: str | None
    benchmark_ktilde: float | None
    range_min_tev: float | None
    range_max_tev: float | None
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
class CR013Anchor:
    """YAML-loaded CR013 diphoton mass-limit bundle."""

    active_limit: CR013ValueAnchor
    cms_2024_range: CR013ValueAnchor
    atlas_2021_ktilde_0p1: CR013ValueAnchor
    atlas_2021_ktilde_0p05: CR013ValueAnchor
    atlas_2021_ktilde_0p01: CR013ValueAnchor
    cms_2018_range: CR013ValueAnchor
    atlas_2015_ktilde_0p1: CR013ValueAnchor
    all_limits: tuple[CR013ValueAnchor, ...]
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
    def benchmark_ktilde(self) -> float:
        if self.active_limit.benchmark_ktilde is None:
            raise AnchorError("CR013: active limit is missing benchmark_ktilde")
        return float(self.active_limit.benchmark_ktilde)

    @property
    def strongest_reported_endpoint(self) -> CR013ValueAnchor:
        return max(self.all_limits, key=lambda limit: limit.value_tev)


def _optional_str(value: Any) -> str | None:
    return None if value is None else str(value)


def _optional_int(value: Any) -> int | None:
    if value is None:
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _required_float(value: Any, *, process_id: str, field_name: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise AnchorError(
            f"{process_id}: CR013 field {field_name!r}={value!r} is not numeric"
        ) from exc
    if not math.isfinite(number):
        raise AnchorError(
            f"{process_id}: CR013 field {field_name!r}={value!r} is not finite"
        )
    return number


def _optional_float(value: Any, *, process_id: str, field_name: str) -> float | None:
    if value is None:
        return None
    return _required_float(value, process_id=process_id, field_name=field_name)


def _parse_benchmark_ktilde(
    model: str | None,
    *,
    process_id: str,
    value_id: str,
) -> float | None:
    if model is None:
        return None
    match = _KTILDE_RE.search(model)
    if match is None:
        return None
    value = _required_float(
        match.group("value"),
        process_id=process_id,
        field_name=f"{value_id}.model.ktilde",
    )
    if value <= 0.0:
        raise AnchorError(f"{process_id}: {value_id}.model.ktilde must be positive")
    return float(value)


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


def _load_mass_limit_anchor(value_id: str, *, process_id: str) -> CR013ValueAnchor:
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

    range_min = _optional_float(
        entry.get("range_min"),
        process_id=process_id,
        field_name=f"{value_id}.range_min",
    )
    range_max = _optional_float(
        entry.get("range_max"),
        process_id=process_id,
        field_name=f"{value_id}.range_max",
    )
    if range_min is not None or range_max is not None:
        if range_min is None or range_max is None:
            raise AnchorError(
                f"{process_id}: {value_id} must provide both range_min and range_max"
            )
        if range_min <= 0.0 or range_max < range_min:
            raise AnchorError(
                f"{process_id}: {value_id} range must satisfy 0 < min <= max"
            )
        if not math.isclose(value, range_max, rel_tol=0.0, abs_tol=1.0e-12):
            raise AnchorError(
                f"{process_id}: {value_id}.value must equal range_max for "
                "lower_limit_range entries"
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
            f"expected {block_key!r} for CR013 value_id {value_id!r}"
        )

    model = _optional_str(entry.get("model"))
    return CR013ValueAnchor(
        anchor=scaffold_anchor,
        value_id=value_id,
        entry_index=index,
        display=_optional_str(entry.get("display")),
        limit_type=_optional_str(entry.get("limit_type")),
        cl=_optional_str(entry.get("cl")),
        experiment=_optional_str(entry.get("experiment")),
        analysis_year=_optional_int(entry.get("analysis_year")),
        model=model,
        benchmark_ktilde=_parse_benchmark_ktilde(
            model,
            process_id=process_id,
            value_id=value_id,
        ),
        range_min_tev=range_min,
        range_max_tev=range_max,
        source_key=_optional_str(entry.get("source_key")),
        access_date=_optional_str(entry.get("access_date")),
        sha256=_optional_str(entry.get("sha256")),
    )


def _load_cr013_anchor(process_id: str) -> CR013Anchor:
    limits = tuple(
        _load_mass_limit_anchor(value_id, process_id=process_id)
        for value_id in _KNOWN_VALUE_IDS
    )
    by_id = {limit.value_id: limit for limit in limits}
    parent = _parent(process_id)
    anchor = CR013Anchor(
        active_limit=by_id[_ACTIVE_CMS_KTILDE_0P1],
        cms_2024_range=by_id[_CMS_2024_RANGE],
        atlas_2021_ktilde_0p1=by_id[_ATLAS_2021_KTILDE_0P1],
        atlas_2021_ktilde_0p05=by_id[_ATLAS_2021_KTILDE_0P05],
        atlas_2021_ktilde_0p01=by_id[_ATLAS_2021_KTILDE_0P01],
        cms_2018_range=by_id[_CMS_2018_RANGE],
        atlas_2015_ktilde_0p1=by_id[_ATLAS_2015_KTILDE_0P1],
        all_limits=limits,
        parent_canonical_source=_optional_str(parent.get("canonical_source")),
        parent_source_key=_optional_str(parent.get("source_key")),
        parent_source_url=_optional_str(parent.get("source_url")),
        budget_policy=_BUDGET_POLICY,
    )
    if anchor.budget <= 0.0:
        raise AnchorError(f"{process_id}: CR013 active budget must be positive")
    active_ktilde = anchor.benchmark_ktilde
    exact_ktilde_limits = tuple(
        limit
        for limit in anchor.all_limits
        if limit.benchmark_ktilde is not None
        and math.isclose(limit.benchmark_ktilde, active_ktilde)
    )
    if not exact_ktilde_limits:
        raise AnchorError(f"{process_id}: no exact ktilde limits found")
    strongest_exact = max(exact_ktilde_limits, key=lambda limit: limit.value_tev)
    if strongest_exact.value_id != anchor.active_limit.value_id:
        raise AnchorError(
            f"{process_id}: active limit {anchor.active_limit.value_id!r} is not "
            f"the strongest exact ktilde={active_ktilde} row; strongest is "
            f"{strongest_exact.value_id!r}"
        )
    return anchor


def _limit_from_anchor(
    anchor: CR013ValueAnchor,
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
            "analysis_year": anchor.analysis_year,
            "model": anchor.model,
            "benchmark_ktilde": anchor.benchmark_ktilde,
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
    """Catalogued RS-graviton diphoton high-mass resonance constraint."""

    process_id = "CR013"
    severity = Severity.HARD
    observable = "m(G_KK^(1) -> gamma gamma)"

    def __init__(self) -> None:
        self.anchor = _load_cr013_anchor(self.process_id)
        self.limit = _limit_from_anchor(
            self.anchor.active_limit,
            process_id=self.process_id,
        )

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        try:
            mass_mapping = resolve_kk_graviton_mass_mapping(
                ew_mass_extra=point.get_extra(_EW_MASS_EXTRA),
                gluon_mass_extra=point.get_extra(_GLUON_MASS_EXTRA),
                couplings=point.get_extra(_COUPLINGS_EXTRA),
            )
        except (TypeError, ValueError, AttributeError) as exc:
            return ConstraintResult(
                process_id=self.process_id,
                severity=self.severity,
                passes=False,
                experimental=float(self.anchor.value),
                budget=float(self.anchor.budget),
                notes=f"invalid KK-graviton diphoton mass proxy input; CR013 not evaluated: {exc}",
                diagnostics={
                    "exception_type": type(exc).__name__,
                    "exception": str(exc),
                    "required_mass_sources": (
                        _EW_MASS_EXTRA,
                        _GLUON_MASS_EXTRA,
                        f"{_COUPLINGS_EXTRA}.M_KK",
                    ),
                    "active_value_id": self.anchor.active_limit.value_id,
                    "benchmark_ktilde": float(self.anchor.benchmark_ktilde),
                    "needs_human_physics": (
                        KK_GRAVITON_DIPHOTON_RECAST_PROXY_ASSUMPTION_V1
                    ),
                    "sigma_times_br_recast_status": "NEEDS-HUMAN-PHYSICS",
                },
            )

        if mass_mapping is None:
            return ConstraintResult(
                process_id=self.process_id,
                severity=self.severity,
                passes=True,
                experimental=float(self.anchor.value),
                budget=float(self.anchor.budget),
                notes=(
                    f"extras {_EW_MASS_EXTRA!r}, {_GLUON_MASS_EXTRA!r}, and "
                    f"{_COUPLINGS_EXTRA!r} absent; diphoton KK-graviton "
                    "resonance constraint was not evaluated."
                ),
                diagnostics={
                    "missing_extras": (
                        _EW_MASS_EXTRA,
                        _GLUON_MASS_EXTRA,
                        _COUPLINGS_EXTRA,
                    ),
                    "required_mass_sources": (
                        _EW_MASS_EXTRA,
                        _GLUON_MASS_EXTRA,
                        f"{_COUPLINGS_EXTRA}.M_KK",
                    ),
                    "budget_policy": self.anchor.budget_policy,
                    "active_value_id": self.anchor.active_limit.value_id,
                    "benchmark_ktilde": float(self.anchor.benchmark_ktilde),
                    "needs_human_physics": (
                        KK_GRAVITON_DIPHOTON_RECAST_PROXY_ASSUMPTION_V1
                    ),
                    "sigma_times_br_recast_status": "NEEDS-HUMAN-PHYSICS",
                    "branching_surface_recast_status": "NEEDS-HUMAN-PHYSICS",
                    "coupling_recast_status": "NEEDS-HUMAN-PHYSICS",
                    "width_acceptance_recast_status": "NEEDS-HUMAN-PHYSICS",
                },
            )

        prediction = kk_graviton_diphoton_prediction_from_lambda_ir_gev(
            float(mass_mapping.lambda_ir_gev),
            resonance=_RESONANCE,
            final_state=_FINAL_STATE,
            benchmark_ktilde=float(self.anchor.benchmark_ktilde),
        )
        comparison = evaluate_collider_resonance_limit(prediction, self.limit)
        diagnostics = dict(comparison.diagnostics)
        diagnostics.update(
            {
                "lambda_ir_gev": float(mass_mapping.lambda_ir_gev),
                "lambda_ir_source": mass_mapping.lambda_ir_source,
                "input_mass_gev": float(mass_mapping.input_mass_gev),
                "input_mass_role": mass_mapping.input_mass_role,
                "m_graviton_gev": float(mass_mapping.m_graviton_gev),
                "m_graviton_proxy_gev": float(mass_mapping.m_graviton_gev),
                "m_graviton_proxy_tev": float(comparison.predicted_mass_tev),
                "graviton_spin2_root": float(SPIN2_GRAVITON_KK_ROOT),
                "gauge_kk_root_nn": float(GAUGE_KK_ROOT_NN),
                "spectrum_assumption": KK_GRAVITON_SPIN2_SPECTRUM_ASSUMPTION_V1,
                "spectrum_mapping_status": "NEEDS-HUMAN-PHYSICS",
                "mass_proxy": (
                    "m_GKK = SPIN2_GRAVITON_KK_ROOT * Lambda_IR"
                ),
                "mass_source": mass_mapping.mass_source,
                "active_value_id": self.anchor.active_limit.value_id,
                "active_experiment": self.anchor.active_limit.experiment,
                "active_display": self.anchor.active_limit.display,
                "active_model": self.anchor.active_limit.model,
                "active_analysis_year": self.anchor.active_limit.analysis_year,
                "benchmark_ktilde": float(self.anchor.benchmark_ktilde),
                "parent_source_key": self.anchor.parent_source_key,
                "parent_canonical_source": self.anchor.parent_canonical_source,
                "budget_policy": self.anchor.budget_policy,
                "all_diphoton_mass_limits_tev": {
                    limit.value_id: float(limit.value_tev)
                    for limit in self.anchor.all_limits
                },
                "exact_ktilde_mass_limits_tev": {
                    limit.value_id: float(limit.value_tev)
                    for limit in self.anchor.all_limits
                    if limit.benchmark_ktilde is not None
                },
                "range_endpoint_limits_tev": {
                    limit.value_id: {
                        "range_min_tev": limit.range_min_tev,
                        "range_max_tev": limit.range_max_tev,
                        "value_tev": float(limit.value_tev),
                    }
                    for limit in self.anchor.all_limits
                    if limit.range_min_tev is not None
                },
                "strongest_reported_endpoint_value_id": (
                    self.anchor.strongest_reported_endpoint.value_id
                ),
                "strongest_reported_endpoint_tev": float(
                    self.anchor.strongest_reported_endpoint.value_tev
                ),
                "strongest_applicable_value_id": self.anchor.active_limit.value_id,
                "spin0_resonance_recast_status": (
                    "NEEDS-HUMAN-PHYSICS: not evaluated without a scalar X "
                    "spectrum and gamma gamma couplings on ParameterPoint."
                ),
                "needs_human_physics": (
                    KK_GRAVITON_DIPHOTON_RECAST_PROXY_ASSUMPTION_V1
                ),
                "sigma_times_br_recast_status": "NEEDS-HUMAN-PHYSICS",
                "branching_surface_recast_status": "NEEDS-HUMAN-PHYSICS",
                "coupling_recast_status": "NEEDS-HUMAN-PHYSICS",
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
                "Diphoton high-mass resonance uses the documented spin-2 "
                "m_GKK = x_G1 * Lambda_IR spectrum map and the strongest "
                "applicable CMS/PDG RS-graviton gamma gamma mass lower bound "
                "at the YAML benchmark ktilde. HARD ratio is m_limit/m_GKK; "
                "sigma*BR, branching-surface, coupling, width, and acceptance "
                "recasts are marked NEEDS-HUMAN-PHYSICS."
            ),
            diagnostics=diagnostics,
        )
