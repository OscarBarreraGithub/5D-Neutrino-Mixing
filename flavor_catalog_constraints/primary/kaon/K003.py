"""K003 - direct CP violation in ``K -> pi pi``.

Physics
-------
``Re(epsilon'/epsilon)`` measures direct CP violation in kaon decays to two
pions.  This is a Delta S = 1 nonleptonic observable.  The Standard-Model
prediction involves a large cancellation between QCD and electroweak penguins,
Delta I = 1/2 enhancement, isospin-breaking corrections, and RBC/UKQCD
lattice matrix elements with large uncertainties.  The RS new-physics side
would require Delta S = 1 penguin and chromomagnetic matching, RG evolution,
and K -> pi pi matrix elements that are not present on ``ParameterPoint``.

This file is therefore an honest non-vetoing stub, like C003.  It loads the
PDG/NA48/KTeV measured ``Re(epsilon'/epsilon)`` anchor from
``flavor_catalog/processes/kaon/K003.yaml``, records the RBC/UKQCD 2020 SM
context value, and reports the full observed ``|Re(epsilon'/epsilon)|`` as
the documented non-vetoing NP room.  It does not fake a penguin calculation.

Severity
--------
INFO.  Both the SM side and the RS NP side are flagged
``NEEDS-HUMAN-PHYSICS``.  The returned ``passes`` value is advisory only and
must not veto scan points.

Catalog sidecar
---------------
``flavor_catalog/processes/kaon/K003.yaml`` is the source of truth for the
PDG average, KTeV/NA48 experimental inputs, RBC/UKQCD lattice SM context, and
Aebischer-Bobeth-Buras SM benchmarks.  The sidecar stores the primary
``pdg_or_equivalent`` anchor as a flat value block, so K003 routes that block
through the scaffold ``load_anchor`` path via a small virtual-block adapter
instead of hardcoding values.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from typing import Any, Mapping

from flavor_catalog_constraints import anchors as anchor_scaffold
from flavor_catalog_constraints.anchors import Anchor, AnchorError, load_anchor, load_full_yaml
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.kaon_direct_cp import (
    KAON_DIRECT_CP_RS_NEEDS_HUMAN_PHYSICS,
    KAON_DIRECT_CP_SM_NEEDS_HUMAN_PHYSICS,
    KAON_DIRECT_CP_STUB_MODEL_V1,
    compare_epsilon_prime_np_room_to_measurement,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "kaon"
_FLAT_PDG_ANCHOR_KEY = "pdg_or_equivalent"
_EXPECTED_UNITS = "dimensionless"
_SM_CONTEXT_KEY = "rbc_ukqcd_2020"
_REQUIRED_SUPPORTING_KEYS = (
    "pdg_datablock_fit",
    "pdg_datablock_average",
    "ktev_2011",
    "na48_2002",
    "rbc_ukqcd_2020",
    "aebischer_buras_2020_octet",
    "aebischer_buras_2020_nonet",
)
_BUDGET_SOURCE = (
    "flavor_catalog/processes/kaon/K003.yaml pdg_or_equivalent; full "
    "observed |Re(epsilon'/epsilon)| is reported as non-vetoing NP room "
    "because no reliable SM subtraction or RS Delta S=1 penguin matching "
    "exists in this scaffold"
)
_PARAMETRIZATION_CITATION = (
    "PDG 2025 K0 listing S013EPS; KTeV arXiv:1011.0127; NA48 "
    "arXiv:hep-ex/0208009; RBC/UKQCD arXiv:2004.09440; "
    "Aebischer-Bobeth-Buras arXiv:2005.05978"
)
_NEEDS_HUMAN_PHYSICS = (
    KAON_DIRECT_CP_SM_NEEDS_HUMAN_PHYSICS,
    KAON_DIRECT_CP_RS_NEEDS_HUMAN_PHYSICS,
)
_NUMBER_RE = re.compile(r"[0-9]+(?:\.[0-9]*)?(?:[eE][+-]?[0-9]+)?")


@dataclass(frozen=True)
class SupportingValueAnchor:
    """Typed supporting K003 value, scaled to dimensionless units."""

    key: str
    source: str | None
    raw_value: float
    raw_uncertainty: float | None
    raw_uncertainty_components: tuple[float, ...]
    scale: float
    value: float
    uncertainty: float | None
    uncertainty_components: tuple[float, ...]
    units: str | None
    display_value: str | None
    source_url: str | None
    snapshot_path: str | None
    sha256: str | None


@dataclass(frozen=True)
class K003Anchor:
    """Typed K003 measured anchor and non-vetoing room convention."""

    block_key: str
    source: str | None
    year: int | None
    value: float
    uncertainty: float
    units: str | None
    observable: str | None
    source_url: str | None
    snapshot_path: str | None
    display_value: str | None
    sha256: str | None
    supporting_values: tuple[SupportingValueAnchor, ...]

    @property
    def budget(self) -> float:
        """Non-vetoing NP room: the full observed absolute value."""
        return abs(self.value)

    def supporting(self, key: str) -> SupportingValueAnchor:
        """Return one parsed supporting value by sidecar key."""
        for item in self.supporting_values:
            if item.key == key:
                return item
        raise KeyError(key)

    @property
    def sm_context(self) -> SupportingValueAnchor:
        """RBC/UKQCD 2020 SM context value, not a veto-grade prediction."""
        return self.supporting(_SM_CONTEXT_KEY)


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
        out = float(value)
    except (TypeError, ValueError) as exc:
        raise AnchorError(
            f"{process_id}: K003 anchor field {field_name!r}={value!r} "
            "is not numeric"
        ) from exc
    if not math.isfinite(out):
        raise AnchorError(f"{process_id}: K003 anchor field {field_name!r} is not finite")
    return out


def _positive_float(value: Any, *, process_id: str, field_name: str) -> float:
    out = _required_float(value, process_id=process_id, field_name=field_name)
    if out <= 0.0:
        raise AnchorError(f"{process_id}: K003 anchor field {field_name!r} must be positive")
    return out


def _scale_for_units(units: str | None, *, process_id: str, field_name: str) -> float:
    if units == "dimensionless":
        return 1.0
    if units == "10^-3":
        return 1.0e-3
    if units == "10^-4":
        return 1.0e-4
    raise AnchorError(f"{process_id}: unsupported units {units!r} for {field_name}")


def _parse_uncertainty(
    value: Any,
    *,
    process_id: str,
    field_name: str,
) -> tuple[float | None, tuple[float, ...]]:
    if isinstance(value, str):
        components = tuple(float(match.group(0)) for match in _NUMBER_RE.finditer(value))
        if not components:
            raise AnchorError(
                f"{process_id}: K003 uncertainty field {field_name!r}={value!r} "
                "contains no numeric components"
            )
        for index, component in enumerate(components):
            _positive_float(
                component,
                process_id=process_id,
                field_name=f"{field_name}.component[{index}]",
            )
        return None, components

    raw = _positive_float(value, process_id=process_id, field_name=field_name)
    return raw, (raw,)


def _flat_pdg_block(process_id: str) -> Mapping[str, Any]:
    data = load_full_yaml(process_id, family=_FAMILY)
    block = data.get("pdg_or_equivalent")
    if not isinstance(block, Mapping):
        raise AnchorError(f"{process_id}: expected flat pdg_or_equivalent mapping")
    return block


def _load_flat_pdg_anchor(
    process_id: str,
    *,
    value_key: str = "value",
    uncertainty_key: str = "uncertainty",
) -> tuple[Anchor, Mapping[str, Any]]:
    """Route K003's flat ``pdg_or_equivalent`` block through ``load_anchor``."""

    pdg_block = _flat_pdg_block(process_id)
    virtual_block = {_FLAT_PDG_ANCHOR_KEY: dict(pdg_block)}
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
            candidates=(_FLAT_PDG_ANCHOR_KEY,),
            value_key=value_key,
            uncertainty_key=uncertainty_key,
        )
    finally:
        anchor_scaffold.load_pdg_block = original_load_pdg_block

    if scaffold_anchor.block_key != _FLAT_PDG_ANCHOR_KEY:
        raise AnchorError(
            f"{process_id}: load_anchor selected {scaffold_anchor.block_key!r}, "
            f"expected {_FLAT_PDG_ANCHOR_KEY!r}"
        )
    return scaffold_anchor, pdg_block


def _load_supporting_values(process_id: str) -> tuple[SupportingValueAnchor, ...]:
    data = load_full_yaml(process_id, family=_FAMILY)
    raw_values = data.get("supporting_values")
    if not isinstance(raw_values, list) or not raw_values:
        raise AnchorError(f"{process_id}: supporting_values must be a non-empty list")

    parsed: list[SupportingValueAnchor] = []
    for index, entry in enumerate(raw_values):
        if not isinstance(entry, Mapping):
            raise AnchorError(f"{process_id}: supporting_values[{index}] is not a mapping")
        key = _optional_str(entry.get("key"))
        if not key:
            raise AnchorError(f"{process_id}: supporting_values[{index}].key is required")
        units = _optional_str(entry.get("units"))
        scale = _scale_for_units(
            units,
            process_id=process_id,
            field_name=f"supporting_values[{index}].units",
        )
        raw_value = _required_float(
            entry.get("value"),
            process_id=process_id,
            field_name=f"supporting_values[{index}].value",
        )
        raw_uncertainty, raw_components = _parse_uncertainty(
            entry.get("uncertainty"),
            process_id=process_id,
            field_name=f"supporting_values[{index}].uncertainty",
        )
        components = tuple(float(component * scale) for component in raw_components)
        uncertainty = math.sqrt(sum(component * component for component in components))
        parsed.append(
            SupportingValueAnchor(
                key=key,
                source=_optional_str(entry.get("source")),
                raw_value=float(raw_value),
                raw_uncertainty=None if raw_uncertainty is None else float(raw_uncertainty),
                raw_uncertainty_components=tuple(float(c) for c in raw_components),
                scale=float(scale),
                value=float(raw_value * scale),
                uncertainty=float(uncertainty),
                uncertainty_components=components,
                units=units,
                display_value=_optional_str(entry.get("display_value")),
                source_url=_optional_str(entry.get("source_url")),
                snapshot_path=_optional_str(entry.get("snapshot_path")),
                sha256=_optional_str(entry.get("sha256_of_text_snapshot")),
            )
        )

    present = {item.key for item in parsed}
    missing = set(_REQUIRED_SUPPORTING_KEYS) - present
    if missing:
        raise AnchorError(
            f"{process_id}: missing required supporting_values keys {sorted(missing)}"
        )
    return tuple(parsed)


def _load_k003_anchor(process_id: str) -> K003Anchor:
    scaffold_anchor, pdg_block = _load_flat_pdg_anchor(process_id)
    if scaffold_anchor.units != _EXPECTED_UNITS:
        raise AnchorError(
            f"{process_id}: expected units {_EXPECTED_UNITS!r}, "
            f"got {scaffold_anchor.units!r}"
        )
    if scaffold_anchor.uncertainty is None:
        raise AnchorError(f"{process_id}: K003 PDG average must provide uncertainty")
    value = _required_float(
        scaffold_anchor.value,
        process_id=process_id,
        field_name="pdg_or_equivalent.value",
    )
    uncertainty = _positive_float(
        scaffold_anchor.uncertainty,
        process_id=process_id,
        field_name="pdg_or_equivalent.uncertainty",
    )
    anchor = K003Anchor(
        block_key=scaffold_anchor.block_key,
        source=_optional_str(scaffold_anchor.source),
        year=scaffold_anchor.year,
        value=float(value),
        uncertainty=float(uncertainty),
        units=scaffold_anchor.units,
        observable=_optional_str(scaffold_anchor.observable),
        source_url=_optional_str(scaffold_anchor.source_url),
        snapshot_path=_optional_str(scaffold_anchor.snapshot_path),
        display_value=_optional_str(pdg_block.get("display_value")),
        sha256=_optional_str(pdg_block.get("sha256_of_text_snapshot")),
        supporting_values=_load_supporting_values(process_id),
    )
    if anchor.budget <= 0.0:
        raise AnchorError(f"{process_id}: K003 non-vetoing NP room must be positive")
    return anchor


def _supporting_values_diagnostics(anchor: K003Anchor) -> tuple[Mapping[str, Any], ...]:
    return tuple(
        {
            "key": item.key,
            "source": item.source,
            "raw_value": float(item.raw_value),
            "raw_uncertainty": item.raw_uncertainty,
            "raw_uncertainty_components": tuple(
                float(value) for value in item.raw_uncertainty_components
            ),
            "scale": float(item.scale),
            "value": float(item.value),
            "uncertainty": None if item.uncertainty is None else float(item.uncertainty),
            "uncertainty_components": tuple(
                float(value) for value in item.uncertainty_components
            ),
            "units": item.units,
            "display_value": item.display_value,
            "source_url": item.source_url,
            "snapshot_path": item.snapshot_path,
            "sha256": item.sha256,
        }
        for item in anchor.supporting_values
    )


@register
class Constraint:
    """Catalogued non-vetoing epsilon-prime kaon-CP stub (process_id K003)."""

    process_id = "K003"
    severity = Severity.INFO
    observable = "Re(epsilon'/epsilon)"

    def __init__(self) -> None:
        self.anchor = _load_k003_anchor(self.process_id)

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        del point
        sm_context = self.anchor.sm_context
        comparison = compare_epsilon_prime_np_room_to_measurement(
            measured_re_epsilon_prime_over_epsilon=self.anchor.value,
            experimental_uncertainty=self.anchor.uncertainty,
            documented_np_room_abs=self.anchor.budget,
        )

        return ConstraintResult(
            process_id=self.process_id,
            severity=self.severity,
            passes=bool(comparison.passes),
            predicted=None,
            sm_prediction=float(sm_context.value),
            experimental=float(self.anchor.value),
            ratio=float(comparison.ratio_to_room),
            budget=float(comparison.documented_np_room_abs),
            notes=(
                "INFO-only K003 stub: loads the measured Re(epsilon'/epsilon) "
                "anchor and records the RBC/UKQCD 2020 SM context value. "
                "No SM penguin-cancellation calculation or RS Delta S=1 "
                "penguin/chromomagnetic matching is computed; both sides are "
                "flagged NEEDS-HUMAN-PHYSICS, so this result is non-vetoing."
            ),
            diagnostics={
                "non_vetoing": True,
                "passes_semantics": (
                    "non-vetoing only; no Re(epsilon'/epsilon) SM or NP "
                    "prediction was evaluated"
                ),
                "no_penguin_calculation": True,
                "no_rs_delta_s1_matching": True,
                "stub_model": KAON_DIRECT_CP_STUB_MODEL_V1,
                "budget_source": _BUDGET_SOURCE,
                "budget_is_full_observed_epsilon_prime_over_epsilon": True,
                "parametrization_citation": _PARAMETRIZATION_CITATION,
                "needs_human_physics": _NEEDS_HUMAN_PHYSICS,
                "needs_human_physics_sm": KAON_DIRECT_CP_SM_NEEDS_HUMAN_PHYSICS,
                "needs_human_physics_np": KAON_DIRECT_CP_RS_NEEDS_HUMAN_PHYSICS,
                "sm_prediction_field_semantics": (
                    "RBC/UKQCD 2020 context value only; not a veto-grade SM "
                    "prediction or SM subtraction."
                ),
                "sm_qcd_electroweak_penguin_cancellation": True,
                "sm_delta_i_half_lattice_matrix_elements": True,
                "rs_delta_s1_penguin_matching_available": False,
                "rs_delta_s1_chromomagnetic_matching_available": False,
                "parameter_point_inputs_used": (),
                "experimental_block": self.anchor.block_key,
                "experimental_display_value": self.anchor.display_value,
                "experimental_source": self.anchor.source,
                "experimental_source_url": self.anchor.source_url,
                "experimental_snapshot_path": self.anchor.snapshot_path,
                "experimental_sha256": self.anchor.sha256,
                "experimental_uncertainty": float(self.anchor.uncertainty),
                "measurement_abs": float(comparison.measurement_abs),
                "documented_np_room_abs": float(comparison.documented_np_room_abs),
                "measurement_to_np_room_ratio": float(comparison.ratio_to_room),
                "sm_context_key": sm_context.key,
                "sm_context_value": float(sm_context.value),
                "sm_context_uncertainty": float(sm_context.uncertainty or 0.0),
                "sm_context_raw_value_x1e_minus4": float(sm_context.raw_value),
                "sm_context_uncertainty_components": tuple(
                    float(value) for value in sm_context.uncertainty_components
                ),
                "ktev_2011_value": float(self.anchor.supporting("ktev_2011").value),
                "na48_2002_value": float(self.anchor.supporting("na48_2002").value),
                "aebischer_buras_2020_octet_value": float(
                    self.anchor.supporting("aebischer_buras_2020_octet").value
                ),
                "aebischer_buras_2020_nonet_value": float(
                    self.anchor.supporting("aebischer_buras_2020_nonet").value
                ),
                "supporting_values": _supporting_values_diagnostics(self.anchor),
            },
        )
