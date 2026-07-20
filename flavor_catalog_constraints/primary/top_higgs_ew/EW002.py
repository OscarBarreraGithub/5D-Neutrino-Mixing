"""EW002 - first-row CKM unitarity.

Physics
-------
EW002 tests the data-level Standard Model relation

    |V_ud|^2 + |V_us|^2 + |V_ub|^2 = 1,

reported as ``Delta_CKM = sum - 1``.  The current implementation uses the
EW002.yaml PDG first-row-sum anchor as the rigorous observable, with the
catalogued PDG ``V_ud`` and kaon-sector ``V_us`` entries carried as provenance
diagnostics.  EW002.yaml does not provide a standalone ``|V_ub|`` value block;
the quoted PDG first-row sum already includes the numerically tiny ``|V_ub|^2``
term.

RS matching status
------------------
For points carrying the Phase-5a ``rs_charged_current`` extra, EW002 consumes
the stored minimal left-handed charged-current epsilons and evaluates

    Delta_CKM_NP ~= 2 |V_ud|^2 Re(epsilon_ud^e)
                  + 2 |V_us|^2 Re(epsilon_us^light).

No proxy mass scaling is used.  If the extra is absent, EW002 is explicitly
not evaluated and remains non-vetoing.

Severity
--------
SOFT.  This is an observed SM-vs-data tension, not a standalone RS veto.  The
reported ``ratio`` is the one-sigma pull ``|Delta_CKM| / sigma_unitarity``.

Catalog sidecar
---------------
``flavor_catalog/processes/top_higgs_ew/EW002.yaml`` is the source of truth.
Unlike the older mapping-shaped sidecars, EW002 stores ``pdg_or_equivalent`` as
a list of value blocks.  The loader below adapts those list entries into the
scaffold ``load_anchor`` path and fails loudly on missing or mismatched value
IDs.
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
from flavor_catalog_constraints.physics_adapters.charged_current import (
    CHARGED_CURRENT_NONMINIMAL_NEEDS_HUMAN,
    charged_current_source_diagnostics,
    ckm_first_row_delta_np,
)
from flavor_catalog_constraints.physics_adapters.ckm_unitarity import (
    evaluate_first_row_sum,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "top_higgs_ew"
_VUD_VALUE_ID = "PDG2025:EW002:Vud_superallowed"
_VUS_VALUE_ID = "PDG2025:EW002:Vus_kaon_average"
_FIRST_ROW_SUM_VALUE_ID = "PDG2025:EW002:first_row_sum"
_BUDGET_COMPONENT_KEY = "combined_quoted_elsewhere"
_TARGET_SUM = 1.0
_REQUIRED_EXTRA = "rs_charged_current"


@dataclass(frozen=True)
class FirstRowSumAnchor:
    """Typed PDG first-row-sum anchor with component uncertainty metadata."""

    anchor: Anchor
    uncertainty_components: Mapping[str, float]
    budget: float
    budget_source: str

    @property
    def value(self) -> float:
        return self.anchor.value

    @property
    def source_url(self) -> str | None:
        return self.anchor.source_url

    @property
    def block_key(self) -> str:
        return self.anchor.block_key


@dataclass(frozen=True)
class EW002Anchor:
    """Typed EW002 anchor bundle."""

    vud: Anchor
    vus: Anchor
    first_row_sum: FirstRowSumAnchor

    @property
    def value(self) -> float:
        return self.first_row_sum.value

    @property
    def budget(self) -> float:
        return self.first_row_sum.budget


def _required_float(value: Any, *, process_id: str, field_name: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise AnchorError(
            f"{process_id}: EW002 field {field_name!r}={value!r} is not numeric"
        ) from exc
    if not math.isfinite(number):
        raise AnchorError(
            f"{process_id}: EW002 field {field_name!r}={value!r} is not finite"
        )
    return number


def _positive_float(value: Any, *, process_id: str, field_name: str) -> float:
    number = _required_float(value, process_id=process_id, field_name=field_name)
    if number <= 0.0:
        raise AnchorError(f"{process_id}: EW002 field {field_name!r} must be positive")
    return number


def _pdg_entries(process_id: str) -> Sequence[Mapping[str, Any]]:
    data = load_full_yaml(process_id, family=_FAMILY)
    values = data.get("pdg_or_equivalent")
    if not isinstance(values, list):
        raise AnchorError(
            f"{process_id}: expected list-shaped 'pdg_or_equivalent' for EW002, "
            f"got {type(values).__name__}"
        )
    for index, entry in enumerate(values):
        if not isinstance(entry, Mapping):
            raise AnchorError(
                f"{process_id}: pdg_or_equivalent[{index}] is not a mapping "
                f"(got {type(entry).__name__})"
            )
    return values


def _entry_by_value_id(
    process_id: str,
    value_id: str,
) -> tuple[int, Mapping[str, Any]]:
    values = _pdg_entries(process_id)
    matches: list[Mapping[str, Any]] = []
    indexes: list[int] = []
    for index, entry in enumerate(values):
        if entry.get("value_id") == value_id:
            matches.append(entry)
            indexes.append(index)
    if not matches:
        present = [
            str(entry.get("value_id"))
            for entry in values
            if isinstance(entry, Mapping) and entry.get("value_id") is not None
        ]
        raise AnchorError(
            f"{process_id}: value_id {value_id!r} not found in pdg_or_equivalent "
            f"(present value_ids: {present})"
        )
    if len(matches) > 1:
        raise AnchorError(f"{process_id}: duplicate value_id {value_id!r}")
    return indexes[0], matches[0]


def _load_scaffold_list_anchor(
    value_id: str,
    *,
    process_id: str,
    uncertainty_key: str = "uncertainty",
) -> tuple[Anchor, Mapping[str, Any], int]:
    index, entry = _entry_by_value_id(process_id, value_id)
    block_key = f"pdg_or_equivalent[{index}]"
    virtual_block = {block_key: dict(entry)}
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
            uncertainty_key=uncertainty_key,
        )
    finally:
        anchor_scaffold.load_pdg_block = original_load_pdg_block

    if scaffold_anchor.block_key != block_key:
        raise AnchorError(
            f"{process_id}: load_anchor selected {scaffold_anchor.block_key!r}, "
            f"expected {block_key!r} for EW002 value_id {value_id!r}"
        )
    return scaffold_anchor, entry, index


def _load_value_anchor(process_id: str, value_id: str) -> Anchor:
    scaffold_anchor, _, _ = _load_scaffold_list_anchor(
        value_id,
        process_id=process_id,
    )
    return scaffold_anchor


def _load_first_row_sum_anchor(process_id: str) -> FirstRowSumAnchor:
    anchor, entry, _ = _load_scaffold_list_anchor(
        _FIRST_ROW_SUM_VALUE_ID,
        process_id=process_id,
        uncertainty_key="__ew002_no_scalar_uncertainty__",
    )
    raw_components = entry.get("uncertainties")
    if not isinstance(raw_components, Mapping):
        raise AnchorError(
            f"{process_id}: {_FIRST_ROW_SUM_VALUE_ID} must provide an "
            "'uncertainties' mapping"
        )
    required = ("Vud_squared", "Vus_squared", _BUDGET_COMPONENT_KEY)
    components = {
        key: _positive_float(
            raw_components.get(key),
            process_id=process_id,
            field_name=f"{_FIRST_ROW_SUM_VALUE_ID}.uncertainties.{key}",
        )
        for key in required
    }
    budget = components[_BUDGET_COMPONENT_KEY]
    return FirstRowSumAnchor(
        anchor=anchor,
        uncertainty_components=components,
        budget=float(budget),
        budget_source=(
            "flavor_catalog/processes/top_higgs_ew/EW002.yaml "
            f"{_FIRST_ROW_SUM_VALUE_ID}.uncertainties.{_BUDGET_COMPONENT_KEY}"
        ),
    )


def _load_ew002_anchor(process_id: str) -> EW002Anchor:
    vud = _load_value_anchor(process_id, _VUD_VALUE_ID)
    vus = _load_value_anchor(process_id, _VUS_VALUE_ID)
    first_row_sum = _load_first_row_sum_anchor(process_id)
    for name, anchor in (("Vud", vud), ("Vus", vus)):
        if anchor.uncertainty is None or anchor.uncertainty <= 0.0:
            raise AnchorError(
                f"{process_id}: {name} anchor must carry a positive uncertainty"
            )
    return EW002Anchor(vud=vud, vus=vus, first_row_sum=first_row_sum)


@register
class Constraint:
    """Catalogued first-row CKM unitarity consistency test (EW002)."""

    process_id = "EW002"
    severity = Severity.SOFT
    observable = "Delta_CKM first-row CKM unitarity"

    def __init__(self) -> None:
        self.anchor = _load_ew002_anchor(self.process_id)

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        data_result = evaluate_first_row_sum(
            self.anchor.first_row_sum.value,
            self.anchor.budget,
            target_sum=_TARGET_SUM,
        )
        vud_squared = float(self.anchor.vud.value * self.anchor.vud.value)
        vus_squared = float(self.anchor.vus.value * self.anchor.vus.value)
        sum_without_vub = float(vud_squared + vus_squared)
        common_diagnostics: dict[str, Any] = {
            "sm_vs_data_delta_ckm": float(data_result.delta_ckm),
            "sm_vs_data_abs_delta_ckm": float(abs(data_result.delta_ckm)),
            "sm_vs_data_pull_sigma": float(data_result.pull_sigma),
            "target_sum": float(data_result.target_sum),
            "vud": float(self.anchor.vud.value),
            "vud_uncertainty": float(self.anchor.vud.uncertainty),
            "vus": float(self.anchor.vus.value),
            "vus_uncertainty": float(self.anchor.vus.uncertainty),
            "vud_squared_from_yaml_value": vud_squared,
            "vus_squared_from_yaml_value": vus_squared,
            "sum_vud_vus_squared_without_vub": sum_without_vub,
            "pdg_first_row_sum_minus_vud_vus_squared": float(
                self.anchor.first_row_sum.value - sum_without_vub
            ),
            "vub_value_block_in_ew002_yaml": False,
            "vub_policy": (
                "EW002.yaml does not contain a standalone |Vub| value block; "
                "the PDG first-row-sum anchor is used as the central "
                "observable and already includes the negligible |Vub|^2 term. "
                "The optional epsilon_ub term is therefore not included in the "
                "NP shift."
            ),
            "first_row_sum_block": self.anchor.first_row_sum.block_key,
            "vud_block": self.anchor.vud.block_key,
            "vus_block": self.anchor.vus.block_key,
            "budget_source": self.anchor.first_row_sum.budget_source,
            "budget_uncertainty_components": dict(
                self.anchor.first_row_sum.uncertainty_components
            ),
            "nonminimal_charged_current_status": (
                CHARGED_CURRENT_NONMINIMAL_NEEDS_HUMAN
            ),
        }

        charged = point.get_extra(_REQUIRED_EXTRA)
        if charged is None:
            return ConstraintResult(
                process_id=self.process_id,
                severity=self.severity,
                passes=True,
                predicted=None,
                sm_prediction=float(_TARGET_SUM),
                experimental=float(self.anchor.first_row_sum.value),
                ratio=None,
                budget=float(self.anchor.budget),
                notes=(
                    f"extra {_REQUIRED_EXTRA!r} absent; EW002 charged-current "
                    "epsilon shift was not evaluated."
                ),
                diagnostics={
                    "evaluated": False,
                    "missing_extra": _REQUIRED_EXTRA,
                    "np_shift_delta_ckm": 0.0,
                    **common_diagnostics,
                },
            )

        try:
            shift = ckm_first_row_delta_np(
                charged,
                vud_abs=self.anchor.vud.value,
                vus_abs=self.anchor.vus.value,
                vub_abs=None,
            )
            source_diag = charged_current_source_diagnostics(charged)
        except Exception as exc:  # noqa: BLE001 - constraints degrade cleanly
            return ConstraintResult(
                process_id=self.process_id,
                severity=self.severity,
                passes=True,
                predicted=None,
                sm_prediction=float(_TARGET_SUM),
                experimental=float(self.anchor.first_row_sum.value),
                ratio=None,
                budget=float(self.anchor.budget),
                notes=(
                    "NOT EVALUATED - invalid rs_charged_current extra for "
                    "EW002 charged-current epsilon shift."
                ),
                diagnostics={
                    "evaluated": False,
                    "invalid_extra": _REQUIRED_EXTRA,
                    "exception_type": type(exc).__name__,
                    "exception_message": str(exc),
                    "np_shift_delta_ckm": 0.0,
                    **common_diagnostics,
                },
            )

        predicted_sum = float(_TARGET_SUM + shift.delta_ckm_np)
        delta_vs_data = float(predicted_sum - self.anchor.first_row_sum.value)
        pull = float(abs(delta_vs_data) / self.anchor.budget)

        return ConstraintResult(
            process_id=self.process_id,
            severity=self.severity,
            passes=bool(pull <= 1.0),
            predicted=predicted_sum,
            sm_prediction=float(_TARGET_SUM),
            experimental=float(self.anchor.first_row_sum.value),
            ratio=pull,
            budget=float(self.anchor.budget),
            notes=(
                "First-row CKM sum uses Delta_CKM_NP = 2 sum |V_ij|^2 "
                "Re(epsilon_ij) from rs_charged_current. The SOFT ratio "
                "compares 1 + Delta_CKM_NP with the YAML first-row sum."
            ),
            diagnostics={
                "evaluated": True,
                "delta_ckm": delta_vs_data,
                "abs_delta_ckm": float(abs(delta_vs_data)),
                "pull_sigma": pull,
                "np_shift_delta_ckm": float(shift.delta_ckm_np),
                "np_shift_terms": dict(shift.terms),
                "epsilon_ud_e": shift.epsilon_ud_e.epsilon,
                "epsilon_us_light_average_e_mu": shift.epsilon_us_light.epsilon,
                "epsilon_ub_included": False,
                **source_diag,
                **common_diagnostics,
            },
        )
