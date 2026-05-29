"""K013 - radiative long-lived neutral kaon decay ``K_L -> pi0 gamma gamma``.

Physics
-------
``BR(K_L -> pi0 gamma gamma)`` is dominated by long-distance ChPT hadronic
dynamics: O(p^4)/O(p^6) chiral loops plus vector-meson exchange.  The SM value
is a dedicated chiral calculation, not a first-principles short-distance
prediction, and the catalog has no RS Delta S = 1 radiative matching or ChPT
counterterm matching on ``ParameterPoint``.

This file is therefore an honest non-vetoing stub, like C003.  It loads the
PDG/KTeV/NA48 measured branching-fraction anchor from
``flavor_catalog/processes/kaon/K013.yaml``, records it, and reports the full
observed branching fraction as a non-vetoing bookkeeping budget.  It does not
fake a ChPT calculation or an RS amplitude.

Severity
--------
INFO.  Both the SM ChPT side and the RS NP side are flagged
``NEEDS-HUMAN-PHYSICS``.  The returned ``passes`` value is advisory only and
must not veto scan points.

Catalog sidecar
---------------
``flavor_catalog/processes/kaon/K013.yaml`` is the source of truth for the
PDG 2025 average and its KTeV/NA48 inputs.  The sidecar stores
``pdg_or_equivalent`` as a flat value block, so K013 routes that block through
the scaffold ``load_anchor`` path via a small virtual-block adapter instead of
hardcoding values.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Mapping

from flavor_catalog_constraints import anchors as anchor_scaffold
from flavor_catalog_constraints.anchors import Anchor, AnchorError, load_anchor, load_full_yaml
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.radiative_kaon import (
    RADIATIVE_KAON_RS_NEEDS_HUMAN_PHYSICS,
    RADIATIVE_KAON_SM_NEEDS_HUMAN_PHYSICS,
    RADIATIVE_KAON_STUB_MODEL_V1,
    compare_kl_pi0gammagamma_np_room_to_measurement,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "kaon"
_FLAT_PDG_ANCHOR_KEY = "pdg_or_equivalent"
_EXPECTED_UNITS = "dimensionless branching fraction"
_BUDGET_SOURCE = (
    "flavor_catalog/processes/kaon/K013.yaml pdg_or_equivalent; full "
    "observed BR is reported as non-vetoing NP room because no reliable "
    "SM ChPT subtraction or RS radiative matching exists"
)
_PARAMETRIZATION_CITATION = (
    "PDG 2025 K_L listing; KTeV Phys. Rev. D77 112004 (2008); "
    "NA48 Phys. Lett. B536 229 (2002); radiative-kaon ChPT/vector-exchange "
    "context from Cappiello-Cata-D'Ambrosio EPJC 2018 and "
    "Gabbiani-Valencia PRD64 094008 (2001)"
)
_NEEDS_HUMAN_PHYSICS = (
    RADIATIVE_KAON_SM_NEEDS_HUMAN_PHYSICS,
    RADIATIVE_KAON_RS_NEEDS_HUMAN_PHYSICS,
)


@dataclass(frozen=True)
class AverageInputAnchor:
    """Typed KTeV/NA48 input to the PDG K013 average."""

    label: str | None
    experiment: str | None
    value: float
    uncertainties: tuple[float, ...]
    original_reported_branching_fraction: str | None
    vector_exchange_parameter: str | None
    snapshot_path: str | None
    sha256: str | None


@dataclass(frozen=True)
class K013Anchor:
    """Typed K013 PDG average and non-vetoing room convention."""

    block_key: str
    source: str | None
    year: int | None
    value: float
    uncertainty: float
    units: str | None
    source_url: str | None
    snapshot_path: str | None
    display: str | None
    sha256: str | None
    average_inputs: tuple[AverageInputAnchor, ...]

    @property
    def budget(self) -> float:
        """Non-vetoing NP room: the full observed branching fraction."""
        return abs(self.value)


def _optional_str(value: Any) -> str | None:
    return None if value is None else str(value)


def _required_float(value: Any, *, process_id: str, field_name: str) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError) as exc:
        raise AnchorError(
            f"{process_id}: K013 anchor field {field_name!r}={value!r} "
            "is not numeric"
        ) from exc
    if not math.isfinite(out):
        raise AnchorError(f"{process_id}: K013 anchor field {field_name!r} is not finite")
    return out


def _positive_float(value: Any, *, process_id: str, field_name: str) -> float:
    out = _required_float(value, process_id=process_id, field_name=field_name)
    if out <= 0.0:
        raise AnchorError(f"{process_id}: K013 anchor field {field_name!r} must be positive")
    return out


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
    """Route K013's flat ``pdg_or_equivalent`` block through ``load_anchor``."""

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


def _load_average_inputs(
    raw_inputs: Any,
    *,
    process_id: str,
) -> tuple[AverageInputAnchor, ...]:
    if not isinstance(raw_inputs, list) or not raw_inputs:
        raise AnchorError(f"{process_id}: pdg_or_equivalent.average_inputs must be a list")

    parsed: list[AverageInputAnchor] = []
    for index, entry in enumerate(raw_inputs):
        if not isinstance(entry, Mapping):
            raise AnchorError(
                f"{process_id}: average_inputs[{index}] is not a mapping"
            )
        raw_uncertainties = entry.get("uncertainties")
        if not isinstance(raw_uncertainties, list) or not raw_uncertainties:
            raise AnchorError(
                f"{process_id}: average_inputs[{index}].uncertainties must be "
                "a non-empty list"
            )
        uncertainties = tuple(
            _positive_float(
                uncertainty,
                process_id=process_id,
                field_name=f"average_inputs[{index}].uncertainties[{u_index}]",
            )
            for u_index, uncertainty in enumerate(raw_uncertainties)
        )
        parsed.append(
            AverageInputAnchor(
                label=_optional_str(entry.get("label")),
                experiment=_optional_str(entry.get("experiment")),
                value=_positive_float(
                    entry.get("value"),
                    process_id=process_id,
                    field_name=f"average_inputs[{index}].value",
                ),
                uncertainties=uncertainties,
                original_reported_branching_fraction=_optional_str(
                    entry.get("original_reported_branching_fraction")
                ),
                vector_exchange_parameter=_optional_str(
                    entry.get("vector_exchange_parameter")
                ),
                snapshot_path=_optional_str(entry.get("snapshot_path")),
                sha256=_optional_str(entry.get("sha256")),
            )
        )
    return tuple(parsed)


def _load_k013_anchor(process_id: str) -> K013Anchor:
    scaffold_anchor, pdg_block = _load_flat_pdg_anchor(process_id)
    if scaffold_anchor.units != _EXPECTED_UNITS:
        raise AnchorError(
            f"{process_id}: expected units {_EXPECTED_UNITS!r}, "
            f"got {scaffold_anchor.units!r}"
        )
    if scaffold_anchor.uncertainty is None:
        raise AnchorError(f"{process_id}: K013 PDG average must provide uncertainty")
    value = _positive_float(
        scaffold_anchor.value,
        process_id=process_id,
        field_name="pdg_or_equivalent.value",
    )
    uncertainty = _positive_float(
        scaffold_anchor.uncertainty,
        process_id=process_id,
        field_name="pdg_or_equivalent.uncertainty",
    )
    anchor = K013Anchor(
        block_key=scaffold_anchor.block_key,
        source=_optional_str(scaffold_anchor.source),
        year=scaffold_anchor.year,
        value=float(value),
        uncertainty=float(uncertainty),
        units=scaffold_anchor.units,
        source_url=_optional_str(scaffold_anchor.source_url),
        snapshot_path=_optional_str(scaffold_anchor.snapshot_path),
        display=_optional_str(pdg_block.get("display")),
        sha256=_optional_str(pdg_block.get("sha256")),
        average_inputs=_load_average_inputs(
            pdg_block.get("average_inputs"),
            process_id=process_id,
        ),
    )
    if anchor.budget <= 0.0:
        raise AnchorError(f"{process_id}: K013 non-vetoing NP room must be positive")
    return anchor


def _average_inputs_diagnostics(anchor: K013Anchor) -> tuple[Mapping[str, Any], ...]:
    return tuple(
        {
            "label": item.label,
            "experiment": item.experiment,
            "value": float(item.value),
            "uncertainties": tuple(float(value) for value in item.uncertainties),
            "original_reported_branching_fraction": (
                item.original_reported_branching_fraction
            ),
            "vector_exchange_parameter": item.vector_exchange_parameter,
            "snapshot_path": item.snapshot_path,
            "sha256": item.sha256,
        }
        for item in anchor.average_inputs
    )


@register
class Constraint:
    """Catalogued non-vetoing radiative-kaon ChPT stub (process_id K013)."""

    process_id = "K013"
    severity = Severity.INFO
    observable = "BR(K_L -> pi0 gamma gamma)"

    def __init__(self) -> None:
        self.anchor = _load_k013_anchor(self.process_id)

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        del point
        comparison = compare_kl_pi0gammagamma_np_room_to_measurement(
            measured_branching_fraction=self.anchor.value,
            experimental_uncertainty=self.anchor.uncertainty,
            documented_np_room_abs=self.anchor.budget,
        )

        return ConstraintResult(
            process_id=self.process_id,
            severity=self.severity,
            passes=bool(comparison.passes),
            predicted=None,
            sm_prediction=None,
            experimental=float(self.anchor.value),
            ratio=float(comparison.ratio_to_room),
            budget=float(comparison.documented_np_room_abs),
            notes=(
                "INFO-only K013 stub: loads the measured BR(K_L -> pi0 gamma "
                "gamma) anchor and reports the full observed branching "
                "fraction as non-vetoing NP room. No SM ChPT chiral-loop "
                "amplitude or RS Delta S=1 radiative amplitude is computed; "
                "both sides are flagged NEEDS-HUMAN-PHYSICS."
            ),
            diagnostics={
                "non_vetoing": True,
                "passes_semantics": (
                    "non-vetoing only; no BR(K_L -> pi0 gamma gamma) SM or "
                    "NP prediction was evaluated"
                ),
                "no_chpt_calculation": True,
                "no_np_matching": True,
                "stub_model": RADIATIVE_KAON_STUB_MODEL_V1,
                "budget_source": _BUDGET_SOURCE,
                "budget_is_full_observed_branching_fraction": True,
                "parametrization_citation": _PARAMETRIZATION_CITATION,
                "needs_human_physics": _NEEDS_HUMAN_PHYSICS,
                "needs_human_physics_sm": RADIATIVE_KAON_SM_NEEDS_HUMAN_PHYSICS,
                "needs_human_physics_np": RADIATIVE_KAON_RS_NEEDS_HUMAN_PHYSICS,
                "sm_chpt_chiral_loop_dominated": True,
                "sm_dedicated_chpt_calculation_available": False,
                "rs_delta_s1_radiative_matching_available": False,
                "parameter_point_inputs_used": (),
                "experimental_block": self.anchor.block_key,
                "experimental_display": self.anchor.display,
                "experimental_source": self.anchor.source,
                "experimental_source_url": self.anchor.source_url,
                "experimental_snapshot_path": self.anchor.snapshot_path,
                "experimental_sha256": self.anchor.sha256,
                "experimental_uncertainty": float(self.anchor.uncertainty),
                "measurement_abs": float(comparison.measurement_abs),
                "documented_np_room_abs": float(comparison.documented_np_room_abs),
                "measurement_to_np_room_ratio": float(comparison.ratio_to_room),
                "pdg_average_input_count": float(len(self.anchor.average_inputs)),
                "pdg_average_inputs": _average_inputs_diagnostics(self.anchor),
                "pdg_average_input_experiments": tuple(
                    item.experiment for item in self.anchor.average_inputs
                ),
            },
        )
