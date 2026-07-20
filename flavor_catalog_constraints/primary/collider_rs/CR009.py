"""CR009 - high-mass Drell-Yan contact-interaction scale limit.

Physics
-------
The observable is the effective ``llqq`` contact scale ``Lambda`` constrained
by high-mass Drell-Yan dilepton tails,

    pp -> l+ l-  (m_ll high),    Lambda_llqq > Lambda_limit.

The reusable scale-vs-lower-limit comparison is the ``collider_resonance``
mass-limit machinery in ``quarkConstraints.collider_resonance``, reached only
through ``flavor_catalog_constraints.physics_adapters.collider_resonance``.

RS matching status
------------------
NEEDS-HUMAN-PHYSICS.  A rigorous RS recast needs neutral electroweak KK/Z/Z'
exchange matched to semileptonic SMEFT operators, light-quark and charged-
lepton chiral couplings, the interference sign, EFT-validity treatment, PDFs,
EW corrections, and a binned high-mass dilepton likelihood.  ``ParameterPoint``
does not carry those collider/operator inputs, so CR009 uses the documented
proxy ``Lambda_RS = kk_ew_mass_gev`` with fallbacks to ``kk_gluon_mass_gev`` or
``quark_mass_basis_couplings.M_KK`` and compares that scale to the catalogued
contact-interaction lower limit.

Severity
--------
HARD.  The active budget is the conservative floor of the PDG 2025 full-Run-2
``llqq`` contact-scale lower limits loaded from ``CR009.yaml``.  The ratio is
``Lambda_limit / Lambda_RS_proxy`` so ``ratio <= 1`` passes.  Stronger
helicity and interference benchmarks remain in diagnostics because
``ParameterPoint`` does not carry an operator decomposition.

Catalog sidecar
---------------
``flavor_catalog/processes/collider_rs/CR009.yaml`` is the source of truth for
all contact-scale limits and provenance.  CR009 stores its limits in a
``pdg_or_equivalent.values`` list, so this module adapts selected list entries
into the scaffold ``load_anchor`` path and fails loudly if the expected value
IDs or units are missing.
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
    DY_CONTACT_OPERATOR_PROXY_ASSUMPTION_V1,
    MASS_LOWER_BOUND,
    ColliderResonanceLimit,
    dy_contact_prediction_from_scale_gev,
    evaluate_collider_resonance_limit,
    resolve_dy_contact_scale_gev,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "collider_rs"
_EW_MASS_EXTRA = "kk_ew_mass_gev"
_GLUON_MASS_EXTRA = "kk_gluon_mass_gev"
_COUPLINGS_EXTRA = "quark_mass_basis_couplings"
_PARENT_KEY = "pdg_or_equivalent"
_SCAFFOLD_UNCERTAINTY_KEY = "__cr009_uncertainty_is_not_used__"
_SCALE_UNITS = "TeV"
_RESONANCE = "llqq contact operator"
_FINAL_STATE = "ee + mumu high-mass tail"

_ATLAS_LL_CONSTRUCTIVE = "PDG2025:CR009:ATLAS_LL_constructive"
_ATLAS_LL_DESTRUCTIVE = "PDG2025:CR009:ATLAS_LL_destructive"
_ATLAS_RR_CONSTRUCTIVE = "PDG2025:CR009:ATLAS_RR_constructive"
_ATLAS_RR_DESTRUCTIVE = "PDG2025:CR009:ATLAS_RR_destructive"
_ATLAS_LR_CONSTRUCTIVE = "PDG2025:CR009:ATLAS_LR_constructive"
_ATLAS_LR_DESTRUCTIVE = "PDG2025:CR009:ATLAS_LR_destructive"
_CMS_LL_DESTRUCTIVE_ENDPOINT = "PDG2025:CR009:CMS_LL_destructive_range_endpoint"
_CMS_RR_CONSTRUCTIVE_ENDPOINT = "PDG2025:CR009:CMS_RR_constructive_range_endpoint"
_KNOWN_VALUE_IDS = (
    _ATLAS_LL_CONSTRUCTIVE,
    _ATLAS_LL_DESTRUCTIVE,
    _ATLAS_RR_CONSTRUCTIVE,
    _ATLAS_RR_DESTRUCTIVE,
    _ATLAS_LR_CONSTRUCTIVE,
    _ATLAS_LR_DESTRUCTIVE,
    _CMS_LL_DESTRUCTIVE_ENDPOINT,
    _CMS_RR_CONSTRUCTIVE_ENDPOINT,
)
_BUDGET_POLICY = (
    "conservative DY contact policy: active HARD budget intentionally uses "
    "the weakest PDG2025 full-Run-2 high-mass dilepton llqq "
    "contact-interaction lower limit in CR009.yaml because ParameterPoint "
    "does not provide RS llqq helicity/interference matching. Stronger "
    "helicity and interference benchmarks are retained as diagnostics."
)
_HELICITY_INTERFERENCE_MATCHING_STATUS = "NEEDS-HUMAN-PHYSICS"


@dataclass(frozen=True)
class CR009ValueAnchor:
    """Typed CR009 value entry routed through the scaffold ``load_anchor``."""

    anchor: Anchor
    value_id: str
    entry_index: int
    display: str | None
    limit_type: str | None
    cl: str | None
    experiment: str | None
    source_key: str | None
    model_assumptions: Mapping[str, str]

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

    @property
    def helicity(self) -> str | None:
        return self.model_assumptions.get("helicity")

    @property
    def interference(self) -> str | None:
        return self.model_assumptions.get("interference")

    @property
    def normalization(self) -> str | None:
        return self.model_assumptions.get("normalization")


@dataclass(frozen=True)
class CR009Anchor:
    """YAML-loaded CR009 contact-scale limit bundle."""

    active_limit: CR009ValueAnchor
    all_limits: tuple[CR009ValueAnchor, ...]
    parent_source: str | None
    parent_source_key: str | None
    parent_source_url: str | None
    parent_snapshot_path: str | None
    contact_matching_relation: str | None
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


def _required_float(value: Any, *, process_id: str, field_name: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise AnchorError(
            f"{process_id}: CR009 field {field_name!r}={value!r} is not numeric"
        ) from exc
    if not math.isfinite(number):
        raise AnchorError(
            f"{process_id}: CR009 field {field_name!r}={value!r} is not finite"
        )
    return number


def _string_mapping(
    value: Any,
    *,
    process_id: str,
    field_name: str,
) -> Mapping[str, str]:
    if not isinstance(value, Mapping):
        raise AnchorError(f"{process_id}: {field_name} must be a mapping")
    return {str(key): str(item) for key, item in value.items()}


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


def _contact_matching_relation(process_id: str) -> str | None:
    data = load_full_yaml(process_id, family=_FAMILY)
    paper_ref = data.get("paper_era_reference")
    if not isinstance(paper_ref, Mapping):
        return None
    matching = paper_ref.get("contact_scale_to_rs_vector_matching")
    if not isinstance(matching, Mapping):
        return None
    return _optional_str(matching.get("relation"))


def _load_contact_limit_anchor(
    value_id: str,
    *,
    process_id: str,
) -> CR009ValueAnchor:
    index, entry = _entry_by_value_id(process_id, value_id)
    value = _required_float(
        entry.get("value"),
        process_id=process_id,
        field_name=f"{value_id}.value",
    )
    if value <= 0.0:
        raise AnchorError(f"{process_id}: {value_id}.value must be positive")
    if entry.get("units") != _SCALE_UNITS:
        raise AnchorError(
            f"{process_id}: expected units {_SCALE_UNITS!r} for "
            f"{value_id}, got {entry.get('units')!r}"
        )
    if entry.get("limit_type") != "lower_limit":
        raise AnchorError(
            f"{process_id}: expected lower_limit for {value_id}, "
            f"got {entry.get('limit_type')!r}"
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
            f"expected {block_key!r} for CR009 value_id {value_id!r}"
        )

    return CR009ValueAnchor(
        anchor=scaffold_anchor,
        value_id=value_id,
        entry_index=index,
        display=_optional_str(entry.get("display")),
        limit_type=_optional_str(entry.get("limit_type")),
        cl=_optional_str(entry.get("cl")),
        experiment=_optional_str(entry.get("experiment")),
        source_key=_optional_str(entry.get("source_key")),
        model_assumptions=_string_mapping(
            entry.get("model_assumptions"),
            process_id=process_id,
            field_name=f"{value_id}.model_assumptions",
        ),
    )


def _load_cr009_anchor(process_id: str) -> CR009Anchor:
    limits = tuple(
        _load_contact_limit_anchor(value_id, process_id=process_id)
        for value_id in _KNOWN_VALUE_IDS
    )
    active_limit = min(limits, key=lambda limit: limit.value_tev)
    parent = _parent(process_id)
    anchor = CR009Anchor(
        active_limit=active_limit,
        all_limits=limits,
        parent_source=_optional_str(parent.get("canonical_source")),
        parent_source_key=_optional_str(parent.get("source_key")),
        parent_source_url=_optional_str(parent.get("source_url")),
        parent_snapshot_path=_optional_str(parent.get("snapshot_path")),
        contact_matching_relation=_contact_matching_relation(process_id),
        budget_policy=_BUDGET_POLICY,
    )
    if anchor.budget <= 0.0:
        raise AnchorError(f"{process_id}: CR009 active budget must be positive")
    return anchor


def _limit_from_anchor(
    anchor: CR009ValueAnchor,
    *,
    process_id: str,
) -> ColliderResonanceLimit:
    return ColliderResonanceLimit(
        process_id=process_id,
        resonance=_RESONANCE,
        final_state=_FINAL_STATE,
        limit_kind=MASS_LOWER_BOUND,
        value=float(anchor.value_tev),
        units=_SCALE_UNITS,
        cl=anchor.cl,
        source=anchor.source,
        source_url=anchor.source_url,
        limit_type=anchor.limit_type,
        benchmark_model=anchor.anchor.observable,
        diagnostics={
            "value_id": anchor.value_id,
            "entry_index": anchor.entry_index,
            "display": anchor.display,
            "experiment": anchor.experiment,
            "source_key": anchor.source_key,
            "snapshot_path": anchor.snapshot_path,
            "yaml_units": anchor.units,
            "yaml_observable": anchor.anchor.observable,
            "model_assumptions": dict(anchor.model_assumptions),
            "helicity": anchor.helicity,
            "interference": anchor.interference,
            "normalization": anchor.normalization,
        },
    )


@register
class Constraint:
    """Catalogued high-mass Drell-Yan contact-operator scale constraint."""

    process_id = "CR009"
    severity = Severity.HARD
    observable = "Lambda(llqq contact operator)"

    def __init__(self) -> None:
        self.anchor = _load_cr009_anchor(self.process_id)
        self.limit = _limit_from_anchor(
            self.anchor.active_limit,
            process_id=self.process_id,
        )

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        try:
            lambda_proxy_gev, scale_source = resolve_dy_contact_scale_gev(
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
                notes=f"invalid Drell-Yan contact-scale proxy input; CR009 not evaluated: {exc}",
                diagnostics={
                    "exception_type": type(exc).__name__,
                    "exception": str(exc),
                    "required_scale_sources": (
                        _EW_MASS_EXTRA,
                        _GLUON_MASS_EXTRA,
                        f"{_COUPLINGS_EXTRA}.M_KK",
                    ),
                    "budget_policy": self.anchor.budget_policy,
                    "active_value_id": self.anchor.active_limit.value_id,
                    "helicity_interference_matching_status": (
                        _HELICITY_INTERFERENCE_MATCHING_STATUS
                    ),
                    "needs_human_physics": DY_CONTACT_OPERATOR_PROXY_ASSUMPTION_V1,
                    "eft_recast_status": "NEEDS-HUMAN-PHYSICS",
                },
            )

        if lambda_proxy_gev is None:
            return ConstraintResult(
                process_id=self.process_id,
                severity=self.severity,
                passes=True,
                experimental=float(self.anchor.value),
                budget=float(self.anchor.budget),
                notes=(
                    f"extras {_EW_MASS_EXTRA!r}, {_GLUON_MASS_EXTRA!r}, and "
                    f"{_COUPLINGS_EXTRA!r} absent; Drell-Yan contact-operator "
                    "constraint was not evaluated."
                ),
                diagnostics={
                    "missing_extras": (
                        _EW_MASS_EXTRA,
                        _GLUON_MASS_EXTRA,
                        _COUPLINGS_EXTRA,
                    ),
                    "required_scale_sources": (
                        _EW_MASS_EXTRA,
                        _GLUON_MASS_EXTRA,
                        f"{_COUPLINGS_EXTRA}.M_KK",
                    ),
                    "budget_policy": self.anchor.budget_policy,
                    "active_value_id": self.anchor.active_limit.value_id,
                    "helicity_interference_matching_status": (
                        _HELICITY_INTERFERENCE_MATCHING_STATUS
                    ),
                    "needs_human_physics": DY_CONTACT_OPERATOR_PROXY_ASSUMPTION_V1,
                    "eft_recast_status": "NEEDS-HUMAN-PHYSICS",
                },
            )

        prediction = dy_contact_prediction_from_scale_gev(
            float(lambda_proxy_gev),
            resonance=_RESONANCE,
            final_state=_FINAL_STATE,
        )
        comparison = evaluate_collider_resonance_limit(prediction, self.limit)
        diagnostics = dict(comparison.diagnostics)
        diagnostics.update(
            {
                "lambda_rs_proxy_gev": float(lambda_proxy_gev),
                "lambda_rs_proxy_tev": float(comparison.predicted_mass_tev),
                "contact_scale_proxy": (
                    "Lambda_RS = kk_ew_mass_gev, kk_gluon_mass_gev, or M_KK"
                ),
                "scale_source": scale_source,
                "active_value_id": self.anchor.active_limit.value_id,
                "active_experiment": self.anchor.active_limit.experiment,
                "active_display": self.anchor.active_limit.display,
                "active_helicity": self.anchor.active_limit.helicity,
                "active_interference": self.anchor.active_limit.interference,
                "active_normalization": self.anchor.active_limit.normalization,
                "parent_source": self.anchor.parent_source,
                "parent_source_key": self.anchor.parent_source_key,
                "parent_source_url": self.anchor.parent_source_url,
                "parent_snapshot_path": self.anchor.parent_snapshot_path,
                "budget_policy": self.anchor.budget_policy,
                "helicity_interference_matching_status": (
                    _HELICITY_INTERFERENCE_MATCHING_STATUS
                ),
                "contact_scale_to_rs_vector_relation": (
                    self.anchor.contact_matching_relation
                ),
                "all_contact_limits_tev": {
                    limit.value_id: float(limit.value_tev)
                    for limit in self.anchor.all_limits
                },
                "stronger_contact_limits_tev": {
                    limit.value_id: float(limit.value_tev)
                    for limit in self.anchor.all_limits
                    if limit.value_tev > self.anchor.active_limit.value_tev
                },
                "contact_limits_by_benchmark": {
                    limit.value_id: {
                        "value_tev": float(limit.value_tev),
                        "experiment": limit.experiment,
                        "helicity": limit.helicity,
                        "interference": limit.interference,
                    }
                    for limit in self.anchor.all_limits
                },
                "needs_human_physics": DY_CONTACT_OPERATOR_PROXY_ASSUMPTION_V1,
                "eft_recast_status": "NEEDS-HUMAN-PHYSICS",
                "sigma_or_likelihood_recast_status": "NEEDS-HUMAN-PHYSICS",
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
                "High-mass Drell-Yan llqq contact scale uses the documented "
                "Lambda_RS proxy and the conservative weakest catalogued "
                "PDG2025 full-Run-2 contact-interaction lower limit. HARD "
                "ratio is Lambda_limit/Lambda_RS_proxy; helicity/interference "
                "matching and the EFT or binned-likelihood recast are marked "
                "NEEDS-HUMAN-PHYSICS."
            ),
            diagnostics=diagnostics,
        )
