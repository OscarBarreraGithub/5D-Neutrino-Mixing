"""K018 - K_l3 extraction of ``|V_us|`` from ``|V_us| f_+(0)``.

Physics
-------
K_l3 determines the source-level semileptonic product
``|V_us| f_+(0)``.  Combining that product with the FLAG lattice
normalization ``f_+(0)`` gives

    |V_us| = (|V_us| f_+(0)) / f_+(0).

The arithmetic and uncertainty propagation live in
``quarkConstraints.ckm_extraction`` and are reached only through
``flavor_catalog_constraints.physics_adapters.ckm_extraction``.

RS matching
-----------
For points carrying ``rs_charged_current``, K018 evaluates the minimal
left-handed W/W' shift

    |V_us|_app = |V_us| |1 + epsilon_us^light|.

The light epsilon is the e/mu average because the K018 YAML extraction does
not provide semileptonic mode weights.  Radiative/isospin bookkeeping remains
diagnostic-only; no mass proxy is used.

Severity
--------
SOFT.  K018 is a source-level CKM input and data/lattice consistency check,
not a standalone FCNC veto.  The reported ratio is the pull between the
recomputed extraction and the YAML-derived ``|V_us|`` reference, using the
K018.yaml total ``|V_us|`` uncertainty as the budget.

Catalog sidecar
---------------
``flavor_catalog/processes/kaon/K018.yaml`` is the source of truth for
``|V_us| f_+(0)``, FLAG ``f_+(0)``, and derived ``|V_us|`` values.  K018 uses
a list-shaped ``pdg_or_equivalent`` block, so selected value IDs are adapted
through the scaffold ``load_anchor`` route and fail loudly on mismatch.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Mapping, Sequence

from flavor_catalog_constraints import anchors as anchor_scaffold
from flavor_catalog_constraints.anchors import Anchor, AnchorError, load_anchor, load_full_yaml
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.charged_current import (
    CHARGED_CURRENT_MINIMAL_LH_STATUS,
    CHARGED_CURRENT_NONMINIMAL_NEEDS_HUMAN,
    charged_current_light_epsilon,
    charged_current_source_diagnostics,
    shifted_abs_ckm,
)
from flavor_catalog_constraints.physics_adapters.ckm_extraction import (
    extract_vus_from_kl3,
    vus_consistency_pull,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "kaon"
_EXPECTED_UNITS = "dimensionless"
_FPLUS_VUS_VALUE_ID = "PDG2025:K018:average_fplus_vus"
_FPLUS_VALUE_ID = "FLAG2024:K018:fplus_Nf211"
_PDG_DERIVED_VUS_VALUE_ID = "PDG2025:K018:Vus_from_Kl3"
_FLAG_DERIVED_VUS_VALUE_ID = "FLAG2024:K018:Vus_from_fplus_Nf211"
_BUDGET_SOURCE = (
    "flavor_catalog/processes/kaon/K018.yaml "
    "PDG2025:K018:Vus_from_Kl3.uncertainty_total"
)
_REQUIRED_EXTRA = "rs_charged_current"
_RADIATIVE_ISOSPIN_MODE_STATUS = (
    "PARTIAL: K018 consumes the rigorous minimal-LH epsilon_us shift, but "
    "radiative/isospin conventions and e/mu mode weights are not exposed as "
    "point-dependent inputs."
)


@dataclass(frozen=True)
class K018ValueAnchor:
    """Typed view over one K018 list entry loaded through ``load_anchor``."""

    anchor: Anchor
    value_id: str
    entry_index: int
    uncertainty_total: float | None
    uncertainty_breakdown: Mapping[str, float]
    nf: str | None
    value_summary: str | None

    @property
    def block_key(self) -> str:
        return self.anchor.block_key

    @property
    def value(self) -> float:
        return self.anchor.value

    @property
    def uncertainty(self) -> float | None:
        return self.anchor.uncertainty

    @property
    def source_url(self) -> str | None:
        return self.anchor.source_url

    @property
    def source(self) -> str | None:
        return self.anchor.source

    @property
    def observable(self) -> str | None:
        return self.anchor.observable


@dataclass(frozen=True)
class K018Anchor:
    """All YAML-loaded values needed by the K018 constraint."""

    fplus_vus: K018ValueAnchor
    fplus_zero: K018ValueAnchor
    pdg_derived_vus: K018ValueAnchor
    flag_derived_vus: K018ValueAnchor
    budget: float
    budget_source: str

    @property
    def value(self) -> float:
        return self.pdg_derived_vus.value

    @property
    def source_url(self) -> str | None:
        return self.pdg_derived_vus.source_url


def _required_float(value: Any, *, process_id: str, field_name: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise AnchorError(
            f"{process_id}: K018 field {field_name!r}={value!r} is not numeric"
        ) from exc
    if not math.isfinite(number):
        raise AnchorError(
            f"{process_id}: K018 field {field_name!r}={value!r} is not finite"
        )
    return number


def _positive_float(value: Any, *, process_id: str, field_name: str) -> float:
    number = _required_float(value, process_id=process_id, field_name=field_name)
    if number <= 0.0:
        raise AnchorError(f"{process_id}: K018 field {field_name!r} must be positive")
    return number


def _optional_str(value: Any) -> str | None:
    return None if value is None else str(value)


def _pdg_entries(process_id: str) -> Sequence[Mapping[str, Any]]:
    data = load_full_yaml(process_id, family=_FAMILY)
    entries = data.get("pdg_or_equivalent")
    if not isinstance(entries, list) or not entries:
        raise AnchorError(f"{process_id}: expected non-empty pdg_or_equivalent list")
    for index, entry in enumerate(entries):
        if not isinstance(entry, Mapping):
            raise AnchorError(
                f"{process_id}: pdg_or_equivalent[{index}] is not a mapping"
            )
    return entries


def _entry_by_value_id(
    process_id: str,
    value_id: str,
) -> tuple[int, Mapping[str, Any]]:
    matches: list[tuple[int, Mapping[str, Any]]] = []
    for index, entry in enumerate(_pdg_entries(process_id)):
        if entry.get("value_id") == value_id:
            matches.append((index, entry))
    if not matches:
        present = [
            str(entry.get("value_id"))
            for entry in _pdg_entries(process_id)
            if entry.get("value_id") is not None
        ]
        raise AnchorError(
            f"{process_id}: value_id {value_id!r} not found in "
            f"pdg_or_equivalent list (present value_ids: {present})"
        )
    if len(matches) > 1:
        raise AnchorError(f"{process_id}: duplicate value_id {value_id!r}")
    return matches[0]


def _load_scaffold_list_anchor(
    value_id: str,
    *,
    process_id: str,
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
        )
    finally:
        anchor_scaffold.load_pdg_block = original_load_pdg_block

    if scaffold_anchor.block_key != block_key:
        raise AnchorError(
            f"{process_id}: load_anchor selected {scaffold_anchor.block_key!r}, "
            f"expected {block_key!r} for K018 value_id {value_id!r}"
        )
    return scaffold_anchor, entry, index


def _components_mapping(
    value: Any,
    *,
    process_id: str,
    field_name: str,
) -> Mapping[str, float]:
    if value is None:
        return {}
    if not isinstance(value, Mapping):
        raise AnchorError(f"{process_id}: {field_name} is not a mapping")
    return {
        str(key): _required_float(
            component,
            process_id=process_id,
            field_name=f"{field_name}.{key}",
        )
        for key, component in value.items()
    }


def _load_value_anchor(value_id: str, *, process_id: str) -> K018ValueAnchor:
    anchor, entry, index = _load_scaffold_list_anchor(value_id, process_id=process_id)
    if anchor.units != _EXPECTED_UNITS:
        raise AnchorError(
            f"{process_id}: expected units {_EXPECTED_UNITS!r} for {value_id}, "
            f"got {anchor.units!r}"
        )
    if anchor.uncertainty is not None and anchor.uncertainty <= 0.0:
        raise AnchorError(f"{process_id}: {value_id} uncertainty must be positive")
    total = entry.get("uncertainty_total")
    uncertainty_total = (
        None
        if total is None
        else _positive_float(
            total,
            process_id=process_id,
            field_name=f"pdg_or_equivalent[{index}].uncertainty_total",
        )
    )
    return K018ValueAnchor(
        anchor=anchor,
        value_id=value_id,
        entry_index=index,
        uncertainty_total=uncertainty_total,
        uncertainty_breakdown=_components_mapping(
            entry.get("uncertainty_breakdown"),
            process_id=process_id,
            field_name=f"pdg_or_equivalent[{index}].uncertainty_breakdown",
        ),
        nf=_optional_str(entry.get("nf")),
        value_summary=_optional_str(entry.get("value_summary")),
    )


def _load_k018_anchor(process_id: str) -> K018Anchor:
    fplus_vus = _load_value_anchor(_FPLUS_VUS_VALUE_ID, process_id=process_id)
    fplus_zero = _load_value_anchor(_FPLUS_VALUE_ID, process_id=process_id)
    pdg_derived_vus = _load_value_anchor(
        _PDG_DERIVED_VUS_VALUE_ID,
        process_id=process_id,
    )
    flag_derived_vus = _load_value_anchor(
        _FLAG_DERIVED_VUS_VALUE_ID,
        process_id=process_id,
    )
    for label, value_anchor in (
        ("fplus_vus", fplus_vus),
        ("fplus_zero", fplus_zero),
        ("pdg_derived_vus", pdg_derived_vus),
        ("flag_derived_vus", flag_derived_vus),
    ):
        if value_anchor.uncertainty is None or value_anchor.uncertainty <= 0.0:
            raise AnchorError(f"{process_id}: {label} must carry an uncertainty")
    budget = pdg_derived_vus.uncertainty_total or pdg_derived_vus.uncertainty
    if budget is None or budget <= 0.0:
        raise AnchorError(f"{process_id}: derived |V_us| budget must be positive")
    return K018Anchor(
        fplus_vus=fplus_vus,
        fplus_zero=fplus_zero,
        pdg_derived_vus=pdg_derived_vus,
        flag_derived_vus=flag_derived_vus,
        budget=float(budget),
        budget_source=_BUDGET_SOURCE,
    )


@register
class Constraint:
    """Catalogued K_l3 ``|V_us|`` extraction consistency test (K018)."""

    process_id = "K018"
    severity = Severity.SOFT
    observable = "K_l3 |V_us| f_+(0) extraction"

    def __init__(self) -> None:
        self.anchor = _load_k018_anchor(self.process_id)

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        radiative_isospin = self.anchor.pdg_derived_vus.uncertainty_breakdown.get(
            "radiative_isospin",
            0.0,
        )
        extraction = extract_vus_from_kl3(
            fplus_vus_product=self.anchor.fplus_vus.value,
            fplus_vus_uncertainty=float(self.anchor.fplus_vus.uncertainty),
            fplus_zero=self.anchor.fplus_zero.value,
            fplus_zero_uncertainty=float(self.anchor.fplus_zero.uncertainty),
            radiative_isospin_uncertainty=float(radiative_isospin),
        )
        consistency = vus_consistency_pull(
            extracted_vus=extraction.vus,
            reference_vus=self.anchor.pdg_derived_vus.value,
            budget=self.anchor.budget,
        )
        common_diagnostics: dict[str, Any] = {
            "qcd_running_applied": False,
            "fplus_vus_value_id": self.anchor.fplus_vus.value_id,
            "fplus_value_id": self.anchor.fplus_zero.value_id,
            "pdg_derived_vus_value_id": self.anchor.pdg_derived_vus.value_id,
            "flag_derived_vus_value_id": self.anchor.flag_derived_vus.value_id,
            "fplus_vus": float(extraction.fplus_vus_product),
            "fplus_vus_uncertainty": float(extraction.fplus_vus_uncertainty),
            "fplus_zero": float(extraction.fplus_zero),
            "fplus_zero_uncertainty": float(extraction.fplus_zero_uncertainty),
            "extracted_vus_uncertainty": float(extraction.total_uncertainty),
            "extracted_vus_uncertainty_components": dict(
                extraction.uncertainty_components
            ),
            "pdg_yaml_uncertainty_breakdown": dict(
                self.anchor.pdg_derived_vus.uncertainty_breakdown
            ),
            "sm_delta_vus_vs_pdg_derived": float(consistency.delta_vus),
            "sm_abs_delta_vus_vs_pdg_derived": float(abs(consistency.delta_vus)),
            "flag_derived_vus": float(self.anchor.flag_derived_vus.value),
            "sm_delta_vus_vs_flag_derived": float(
                extraction.vus - self.anchor.flag_derived_vus.value
            ),
            "budget_source": self.anchor.budget_source,
            "fplus_vus_block": self.anchor.fplus_vus.block_key,
            "fplus_block": self.anchor.fplus_zero.block_key,
            "pdg_derived_vus_block": self.anchor.pdg_derived_vus.block_key,
            "flag_derived_vus_block": self.anchor.flag_derived_vus.block_key,
            "minimal_lh_vector_matching_status": CHARGED_CURRENT_MINIMAL_LH_STATUS,
            "radiative_isospin_mode_weight_status": _RADIATIVE_ISOSPIN_MODE_STATUS,
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
                sm_prediction=None,
                experimental=float(self.anchor.pdg_derived_vus.value),
                ratio=None,
                budget=float(consistency.budget),
                notes=(
                    f"extra {_REQUIRED_EXTRA!r} absent; K_l3 epsilon_us shift "
                    "was not evaluated."
                ),
                diagnostics={
                    "evaluated": False,
                    "missing_extra": _REQUIRED_EXTRA,
                    "parameter_point_used": False,
                    "np_shift_vus": 0.0,
                    **common_diagnostics,
                },
            )

        try:
            epsilon_us = charged_current_light_epsilon(charged, up="u", down="s")
            shifted_vus = shifted_abs_ckm(extraction.vus, epsilon_us)
            shifted_consistency = vus_consistency_pull(
                extracted_vus=shifted_vus,
                reference_vus=self.anchor.pdg_derived_vus.value,
                budget=self.anchor.budget,
            )
            source_diag = charged_current_source_diagnostics(charged)
        except Exception as exc:  # noqa: BLE001 - constraints degrade cleanly
            return ConstraintResult(
                process_id=self.process_id,
                severity=self.severity,
                passes=True,
                predicted=None,
                sm_prediction=None,
                experimental=float(self.anchor.pdg_derived_vus.value),
                ratio=None,
                budget=float(consistency.budget),
                notes=(
                    "NOT EVALUATED - invalid rs_charged_current extra for K018 "
                    "epsilon_us shift."
                ),
                diagnostics={
                    "evaluated": False,
                    "invalid_extra": _REQUIRED_EXTRA,
                    "exception_type": type(exc).__name__,
                    "exception_message": str(exc),
                    "parameter_point_used": False,
                    "np_shift_vus": 0.0,
                    **common_diagnostics,
                },
            )

        return ConstraintResult(
            process_id=self.process_id,
            severity=self.severity,
            passes=bool(shifted_consistency.passes),
            predicted=float(shifted_vus),
            sm_prediction=None,
            experimental=float(self.anchor.pdg_derived_vus.value),
            ratio=float(shifted_consistency.pull_sigma),
            budget=float(shifted_consistency.budget),
            notes=(
                "K_l3 |V_us| is recomputed from YAML |V_us| f_+(0) and FLAG "
                "f_+(0), then shifted as |V_us|_app = |V_us| |1 + "
                "epsilon_us^light| from rs_charged_current. The SOFT ratio "
                "compares the shifted value with the YAML-derived |V_us|."
            ),
            diagnostics={
                "evaluated": True,
                "parameter_point_used": True,
                "epsilon_us_light_average_e_mu": epsilon_us.epsilon,
                "abs_1_plus_epsilon_us_light": float(epsilon_us.abs_multiplier),
                "np_shift_vus": float(shifted_vus - extraction.vus),
                "relative_shift_vus": float((shifted_vus / extraction.vus) - 1.0),
                "delta_vus_vs_pdg_derived": float(
                    shifted_consistency.delta_vus
                ),
                "abs_delta_vus_vs_pdg_derived": float(
                    abs(shifted_consistency.delta_vus)
                ),
                **source_diag,
                **common_diagnostics,
            },
        )
