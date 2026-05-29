"""CR007 - spin-2 KK-graviton resonance mass limit.

Physics
-------
The observable is the first spin-2 Randall-Sundrum KK-graviton mass in

    pp -> G_KK^(1) -> ZZ, WW, gamma gamma, ell ell, t tbar

compared with the catalogued RS/bulk-RS benchmark exclusions in
``flavor_catalog/processes/collider_rs/CR007.yaml``.  The reusable
collider-resonance mass-limit comparison machinery is reached only through the
``flavor_catalog_constraints.physics_adapters.kk_graviton_resonance`` adapter,
which reuses the existing ``collider_resonance`` mass-vs-limit recast path.

RS matching status
------------------
NEEDS-HUMAN-PHYSICS.  A faithful spin-2 recast needs
``sigma(pp -> G_KK) * BR(G_KK -> WW/ZZ/gamma gamma/ll/ttbar)``, the
curvature coupling ``k/Mbar_Pl``, total width, branching fractions,
acceptance, interference, and the experiment's mass-dependent limit curve.
The current ``ParameterPoint`` has no dedicated graviton spectrum, so CR007
uses the documented spectrum map ``m_GKK = x_G1 * Lambda_IR`` with the first
spin-2 graviton root.  Gauge-mode mass extras are converted back to
``Lambda_IR`` with the repo gauge root, while quark-sector ``M_KK`` is treated
as the repo ``Lambda_IR`` bookkeeping scale unless the couplings object carries
``xi_KK``.  The precise RS graviton-spectrum mapping and the sigma*BR recast
are explicitly marked ``NEEDS-HUMAN-PHYSICS`` in the result diagnostics.

Severity
--------
HARD.  The active budget is the PDG 2025 / CMS 2023 bulk-RS all-jets
``WW/ZZ`` spin-2 graviton mass lower limit loaded from ``CR007.yaml``.  The
ratio is ``m_limit / m_GKK`` so ``ratio <= 1`` passes.  Stronger generic
diphoton/dilepton RS-graviton limits are retained as diagnostics because they
assume different light-field coupling benchmarks.

Catalog sidecar
---------------
``CR007.yaml`` stores its numerical limits in a ``pdg_or_equivalent.values``
list.  This module adapts selected list entries into the scaffold
``load_anchor`` path and fails loudly if the expected value IDs or units are
missing.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
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
    KK_GRAVITON_SPIN2_MASS_PROXY_ASSUMPTION_V1,
    KK_GRAVITON_SPIN2_SPECTRUM_ASSUMPTION_V1,
    MASS_LOWER_BOUND,
    SPIN2_GRAVITON_KK_ROOT,
    ColliderResonanceLimit,
    evaluate_collider_resonance_limit,
    kk_graviton_spin2_prediction_from_lambda_ir_gev,
    resolve_kk_graviton_mass_mapping,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "collider_rs"
_EW_MASS_EXTRA = "kk_ew_mass_gev"
_GLUON_MASS_EXTRA = "kk_gluon_mass_gev"
_COUPLINGS_EXTRA = "quark_mass_basis_couplings"
_PARENT_KEY = "pdg_or_equivalent"
_SCAFFOLD_UNCERTAINTY_KEY = "__cr007_uncertainty_is_not_used__"
_MASS_UNITS = "TeV"
_RESONANCE = "G_KK^(1)"
_FINAL_STATE = "WW, ZZ"

_ACTIVE_BULK_DIBOSON = "PDG2025:CR007:bulk_graviton_diboson_mass_kMpl_0p5"
_GENERIC_DIPHOTON = "PDG2025:CR007:RS_graviton_diphoton_mass_kMpl_0p1"
_GENERIC_DILEPTON = "PDG2025:CR007:RS_graviton_dilepton_mass_kMpl_0p1"
_ATLAS_BULK_COMBINED = "ATLAS2018:CR007:bulk_graviton_combined_mass_kMpl_1"
_KNOWN_VALUE_IDS = (
    _GENERIC_DIPHOTON,
    _GENERIC_DILEPTON,
    _ACTIVE_BULK_DIBOSON,
    _ATLAS_BULK_COMBINED,
)
_BUDGET_POLICY = (
    "active HARD budget is the PDG2025/CMS2023 bulk-RS spin-2 graviton "
    "WW/ZZ all-jets mass lower limit in CR007.yaml. Generic RS diphoton and "
    "dilepton limits are kept as diagnostics because their light-field "
    "coupling assumptions are not automatically valid for bulk-RS points."
)


@dataclass(frozen=True)
class CR007ValueAnchor:
    """Typed CR007 value entry routed through the scaffold ``load_anchor``."""

    anchor: Anchor
    value_id: str
    entry_index: int
    display: str | None
    limit_type: str | None
    cl: str | None
    experiment: str | None
    channel: str | None
    result_year: int | None
    source_key: str | None
    model_assumptions: str | None
    primary_paper_key: str | None
    primary_paper_snapshot_path: str | None
    primary_paper_sha256: str | None

    @property
    def value(self) -> float:
        return self.anchor.value

    @property
    def value_tev(self) -> float:
        return self.anchor.value

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
class CR007Anchor:
    """YAML-loaded CR007 graviton mass-limit bundle."""

    active_limit: CR007ValueAnchor
    generic_diphoton: CR007ValueAnchor
    generic_dilepton: CR007ValueAnchor
    atlas_bulk_combined: CR007ValueAnchor
    all_limits: tuple[CR007ValueAnchor, ...]
    parent_source: str | None
    parent_source_key: str | None
    parent_source_url: str | None
    summary: str | None
    budget_policy: str

    @property
    def value(self) -> float:
        return self.active_limit.value_tev

    @property
    def budget(self) -> float:
        return self.active_limit.budget

    @property
    def source_url(self) -> str | None:
        return self.active_limit.source_url


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
            f"{process_id}: CR007 field {field_name!r}={value!r} is not numeric"
        ) from exc
    if not math.isfinite(number):
        raise AnchorError(
            f"{process_id}: CR007 field {field_name!r}={value!r} is not finite"
        )
    return number


def _parent(process_id: str) -> Mapping[str, Any]:
    data = load_full_yaml(process_id, family=_FAMILY)
    parent = data.get(_PARENT_KEY)
    if not isinstance(parent, Mapping):
        raise AnchorError(f"{process_id}: expected mapping-shaped {_PARENT_KEY}")
    return parent


def _value_entries(process_id: str) -> Sequence[Mapping[str, Any]]:
    parent = _parent(process_id)
    values = parent.get("values")
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


def _load_value_anchor(value_id: str, *, process_id: str) -> CR007ValueAnchor:
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
            f"expected {block_key!r} for CR007 value_id {value_id!r}"
        )

    return CR007ValueAnchor(
        anchor=scaffold_anchor,
        value_id=value_id,
        entry_index=index,
        display=_optional_str(entry.get("display")),
        limit_type=_optional_str(entry.get("limit_type")),
        cl=_optional_str(entry.get("cl")),
        experiment=_optional_str(entry.get("experiment")),
        channel=_optional_str(entry.get("channel")),
        result_year=_optional_int(entry.get("result_year")),
        source_key=_optional_str(entry.get("source_key")),
        model_assumptions=_optional_str(entry.get("model_assumptions")),
        primary_paper_key=_optional_str(entry.get("primary_paper_key")),
        primary_paper_snapshot_path=_optional_str(
            entry.get("primary_paper_snapshot_path")
        ),
        primary_paper_sha256=_optional_str(entry.get("primary_paper_sha256")),
    )


def _load_cr007_anchor(process_id: str) -> CR007Anchor:
    limits = tuple(
        _load_value_anchor(value_id, process_id=process_id)
        for value_id in _KNOWN_VALUE_IDS
    )
    by_id = {limit.value_id: limit for limit in limits}
    parent = _parent(process_id)
    anchor = CR007Anchor(
        active_limit=by_id[_ACTIVE_BULK_DIBOSON],
        generic_diphoton=by_id[_GENERIC_DIPHOTON],
        generic_dilepton=by_id[_GENERIC_DILEPTON],
        atlas_bulk_combined=by_id[_ATLAS_BULK_COMBINED],
        all_limits=limits,
        parent_source=_optional_str(parent.get("source")),
        parent_source_key=_optional_str(parent.get("source_key")),
        parent_source_url=_optional_str(parent.get("source_url")),
        summary=_optional_str(parent.get("summary")),
        budget_policy=_BUDGET_POLICY,
    )
    if anchor.budget <= 0.0:
        raise AnchorError(f"{process_id}: CR007 active budget must be positive")
    return anchor


def _limit_from_anchor(
    anchor: CR007ValueAnchor,
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
        benchmark_model=anchor.model_assumptions or anchor.anchor.observable,
        diagnostics={
            "value_id": anchor.value_id,
            "entry_index": anchor.entry_index,
            "display": anchor.display,
            "experiment": anchor.experiment,
            "channel": anchor.channel,
            "result_year": anchor.result_year,
            "source_key": anchor.source_key,
            "snapshot_path": anchor.snapshot_path,
            "primary_paper_key": anchor.primary_paper_key,
            "primary_paper_snapshot_path": anchor.primary_paper_snapshot_path,
            "primary_paper_sha256": anchor.primary_paper_sha256,
            "yaml_units": anchor.units,
            "yaml_observable": anchor.anchor.observable,
        },
    )


@register
class Constraint:
    """Catalogued bulk-RS spin-2 KK-graviton mass constraint."""

    process_id = "CR007"
    severity = Severity.HARD
    observable = "m(G_KK^(1) -> WW, ZZ)"

    def __init__(self) -> None:
        self.anchor = _load_cr007_anchor(self.process_id)
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
                notes=f"invalid KK-graviton mass proxy input; CR007 not evaluated: {exc}",
                diagnostics={
                    "exception_type": type(exc).__name__,
                    "exception": str(exc),
                    "required_mass_sources": (
                        _EW_MASS_EXTRA,
                        _GLUON_MASS_EXTRA,
                        f"{_COUPLINGS_EXTRA}.M_KK",
                    ),
                    "needs_human_physics": (
                        KK_GRAVITON_SPIN2_MASS_PROXY_ASSUMPTION_V1
                    ),
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
                    f"{_COUPLINGS_EXTRA!r} absent; spin-2 KK-graviton "
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
                    "needs_human_physics": (
                        KK_GRAVITON_SPIN2_MASS_PROXY_ASSUMPTION_V1
                    ),
                },
            )

        prediction = kk_graviton_spin2_prediction_from_lambda_ir_gev(
            float(mass_mapping.lambda_ir_gev),
            resonance=_RESONANCE,
            final_state=_FINAL_STATE,
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
                "active_channel": self.anchor.active_limit.channel,
                "active_model_assumptions": (
                    self.anchor.active_limit.model_assumptions
                ),
                "parent_source_key": self.anchor.parent_source_key,
                "budget_policy": self.anchor.budget_policy,
                "all_mass_limits_tev": {
                    limit.value_id: float(limit.value_tev)
                    for limit in self.anchor.all_limits
                },
                "generic_rs_mass_limits_tev": {
                    self.anchor.generic_diphoton.value_id: float(
                        self.anchor.generic_diphoton.value_tev
                    ),
                    self.anchor.generic_dilepton.value_id: float(
                        self.anchor.generic_dilepton.value_tev
                    ),
                },
                "bulk_rs_mass_limits_tev": {
                    self.anchor.active_limit.value_id: float(
                        self.anchor.active_limit.value_tev
                    ),
                    self.anchor.atlas_bulk_combined.value_id: float(
                        self.anchor.atlas_bulk_combined.value_tev
                    ),
                },
                "needs_human_physics": (
                    KK_GRAVITON_SPIN2_MASS_PROXY_ASSUMPTION_V1
                ),
                "sigma_times_br_recast_status": "NEEDS-HUMAN-PHYSICS",
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
                "Spin-2 KK-graviton resonance uses the documented "
                "m_GKK = x_G1 * Lambda_IR spectrum map and the active PDG/CMS "
                "bulk-RS WW/ZZ mass lower bound. HARD ratio is "
                "m_limit/m_GKK; spectrum mapping and full sigma*BR recast are "
                "marked NEEDS-HUMAN-PHYSICS."
            ),
            diagnostics=diagnostics,
        )
