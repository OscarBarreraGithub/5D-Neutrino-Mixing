"""C003 - direct CP violation in ``D0 -> K+K-`` and ``D0 -> pi+pi-``.

Physics
-------
``Delta A_CP = A_CP(D0 -> K+K-) - A_CP(D0 -> pi+pi-)`` is an observed direct
charm-CP asymmetry.  Unlike neutral-D mixing, this is a Delta C = 1
nonleptonic observable.  The Standard Model contribution is dominated by
long-distance penguin matrix elements and is not reliably computable from
first principles in this scaffold.  A rigorous RS contribution would require
Delta C = 1 penguin matching, RG evolution, and nonleptonic hadronic matrix
elements that are not present on ``ParameterPoint``.

This file is therefore an honest non-vetoing stub.  It loads the LHCb 2019
``Delta A_CP`` measurement from ``flavor_catalog/processes/charm/C003.yaml``,
records the observable, and reports the full measured ``|Delta A_CP|`` as the
documented NP room for bookkeeping.  It does not fake a penguin calculation.

Severity
--------
INFO.  Both the SM side and the RS NP side are flagged
``NEEDS-HUMAN-PHYSICS``.  The returned ``passes`` value is advisory only and
must not veto scan points.

Catalog sidecar
---------------
``flavor_catalog/processes/charm/C003.yaml`` is the source of truth for the
LHCb discovery measurement, HFLAV direct/indirect CPV companion values, and
CFW paper-era RS context.  Numeric values below are loaded through the
scaffold anchor loader and then scale-converted where the sidecar supplies a
``scale`` field.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Mapping

from flavor_catalog_constraints.anchors import (
    Anchor,
    AnchorError,
    load_anchor,
    load_full_yaml,
)
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.charm_direct_cp import (
    CHARM_DIRECT_CP_RS_NEEDS_HUMAN_PHYSICS,
    CHARM_DIRECT_CP_SM_NEEDS_HUMAN_PHYSICS,
    CHARM_DIRECT_CP_STUB_MODEL_V1,
    compare_delta_acp_np_room_to_measurement,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "charm"
_LHCb_DISCOVERY_CANDIDATES = ("lhcb2019_discovery_average",)
_HFLAV_DIRECT_CPV_CANDIDATES = ("hflav_direct_cpv_average",)
_HFLAV_INDIRECT_CPV_CANDIDATES = ("hflav_indirect_cpv_companion",)
_HFLAV_NO_CPV_KEY = "hflav_no_cpv_test"
_CFW_CONTEXT_KEY = "cfw2008_context"
_EXPECTED_LHCB_UNITS = "dimensionless asymmetry"
_EXPECTED_HFLAV_UNITS = "percent"
_BUDGET_SOURCE = (
    "flavor_catalog/processes/charm/C003.yaml lhcb2019_discovery_average; "
    "full observed |Delta A_CP| is reported as non-vetoing NP room because "
    "no reliable SM subtraction exists"
)
_PARAMETRIZATION_CITATION = (
    "LHCb Phys. Rev. Lett. 122, 211803 (2019), arXiv:1903.08726; "
    "HFLAV direct/indirect charm-CP combination updated 2025"
)
_NEEDS_HUMAN_PHYSICS = (
    CHARM_DIRECT_CP_SM_NEEDS_HUMAN_PHYSICS,
    CHARM_DIRECT_CP_RS_NEEDS_HUMAN_PHYSICS,
)


@dataclass(frozen=True)
class ScaledAsymmetryAnchor:
    """Typed LHCb ``Delta A_CP`` anchor in raw and dimensionless units."""

    block_key: str
    source: str | None
    year: int | None
    raw_value: float
    raw_uncertainty: float
    scale: float
    value: float
    uncertainty: float
    value_percent: float
    uncertainty_percent: float
    units: str | None
    source_url: str | None
    snapshot_path: str | None
    value_summary: str | None


@dataclass(frozen=True)
class PercentAsymmetryAnchor:
    """Typed HFLAV percent-level companion anchor."""

    block_key: str
    source: str | None
    year: int | None
    value_percent: float
    uncertainty_percent: float
    value: float
    uncertainty: float
    units: str | None
    source_url: str | None
    snapshot_path: str | None
    value_summary: str | None


@dataclass(frozen=True)
class NoCPVContext:
    """Typed HFLAV no-CPV significance context."""

    block_key: str
    source: str | None
    year: int | None
    delta_chi2: float
    significance_sigma: float
    confidence_level_table: str | None
    confidence_level_text: str | None
    source_url: str | None
    snapshot_path: str | None


@dataclass(frozen=True)
class CFWContext:
    """Typed paper-era RS context carried only for provenance."""

    block_key: str
    source: str | None
    year: int | None
    value_summary: str | None
    source_url: str | None
    snapshot_path: str | None


@dataclass(frozen=True)
class C003Anchor:
    """Typed C003 anchor bundle and non-vetoing room convention."""

    lhcb_discovery: ScaledAsymmetryAnchor
    hflav_direct_cpv: PercentAsymmetryAnchor
    hflav_indirect_cpv: PercentAsymmetryAnchor
    hflav_no_cpv: NoCPVContext
    cfw_context: CFWContext

    @property
    def value(self) -> float:
        """LHCb ``Delta A_CP`` central value as a dimensionless asymmetry."""
        return self.lhcb_discovery.value

    @property
    def uncertainty(self) -> float:
        """LHCb ``Delta A_CP`` uncertainty as a dimensionless asymmetry."""
        return self.lhcb_discovery.uncertainty

    @property
    def source_url(self) -> str | None:
        """LHCb source URL."""
        return self.lhcb_discovery.source_url

    @property
    def budget(self) -> float:
        """Non-vetoing NP room: the full observed ``|Delta A_CP|``."""
        return abs(self.lhcb_discovery.value)


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
            f"{process_id}: C003 anchor field {field_name!r}={value!r} "
            "is not numeric"
        ) from exc
    if not math.isfinite(out):
        raise AnchorError(f"{process_id}: C003 anchor field {field_name!r} is not finite")
    return out


def _positive_float(value: Any, *, process_id: str, field_name: str) -> float:
    out = _required_float(value, process_id=process_id, field_name=field_name)
    if out <= 0.0:
        raise AnchorError(f"{process_id}: C003 anchor field {field_name!r} must be positive")
    return out


def _pdg_block(process_id: str) -> Mapping[str, Any]:
    data = load_full_yaml(process_id, family=_FAMILY)
    pdg_block = data.get("pdg_or_equivalent")
    if not isinstance(pdg_block, Mapping):
        raise AnchorError(
            f"{process_id}: missing or invalid 'pdg_or_equivalent' mapping"
        )
    return pdg_block


def _pdg_subblock_for_anchor(anchor: Anchor, *, process_id: str) -> Mapping[str, Any]:
    pdg_block = _pdg_block(process_id)
    sub = pdg_block.get(anchor.block_key)
    if not isinstance(sub, Mapping):
        present = sorted(str(key) for key in pdg_block)
        raise AnchorError(
            f"{process_id}: load_anchor selected {anchor.block_key!r}, but that "
            f"block is not available as a mapping (present keys: {present})"
        )
    return sub


def _load_lhcb_discovery_anchor(process_id: str) -> ScaledAsymmetryAnchor:
    scaffold_anchor = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=_LHCb_DISCOVERY_CANDIDATES,
    )
    sub = _pdg_subblock_for_anchor(scaffold_anchor, process_id=process_id)
    if scaffold_anchor.units != _EXPECTED_LHCB_UNITS:
        raise AnchorError(
            f"{process_id}: {scaffold_anchor.block_key} must use units "
            f"{_EXPECTED_LHCB_UNITS!r}, got {scaffold_anchor.units!r}"
        )
    if scaffold_anchor.uncertainty is None:
        raise AnchorError(
            f"{process_id}: {scaffold_anchor.block_key} must provide uncertainty"
        )
    scale = _positive_float(
        sub.get("scale"),
        process_id=process_id,
        field_name=f"{scaffold_anchor.block_key}.scale",
    )
    raw_value = float(scaffold_anchor.value)
    raw_uncertainty = _positive_float(
        scaffold_anchor.uncertainty,
        process_id=process_id,
        field_name=f"{scaffold_anchor.block_key}.uncertainty",
    )
    value = raw_value * scale
    uncertainty = raw_uncertainty * scale
    value_percent = _required_float(
        sub.get("value_percent"),
        process_id=process_id,
        field_name=f"{scaffold_anchor.block_key}.value_percent",
    )
    uncertainty_percent = _positive_float(
        sub.get("uncertainty_percent"),
        process_id=process_id,
        field_name=f"{scaffold_anchor.block_key}.uncertainty_percent",
    )
    if not math.isclose(value_percent / 100.0, value, rel_tol=0.0, abs_tol=1.0e-15):
        raise AnchorError(
            f"{process_id}: value_percent is inconsistent with value*scale"
        )
    if not math.isclose(
        uncertainty_percent / 100.0,
        uncertainty,
        rel_tol=0.0,
        abs_tol=1.0e-15,
    ):
        raise AnchorError(
            f"{process_id}: uncertainty_percent is inconsistent with "
            "uncertainty*scale"
        )
    return ScaledAsymmetryAnchor(
        block_key=scaffold_anchor.block_key,
        source=_optional_str(scaffold_anchor.source),
        year=scaffold_anchor.year,
        raw_value=raw_value,
        raw_uncertainty=raw_uncertainty,
        scale=float(scale),
        value=float(value),
        uncertainty=float(uncertainty),
        value_percent=float(value_percent),
        uncertainty_percent=float(uncertainty_percent),
        units=scaffold_anchor.units,
        source_url=_optional_str(scaffold_anchor.source_url),
        snapshot_path=_optional_str(scaffold_anchor.snapshot_path),
        value_summary=_optional_str(sub.get("value_summary")),
    )


def _load_percent_anchor(
    candidates: tuple[str, ...],
    *,
    process_id: str,
) -> PercentAsymmetryAnchor:
    scaffold_anchor = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=candidates,
    )
    sub = _pdg_subblock_for_anchor(scaffold_anchor, process_id=process_id)
    if scaffold_anchor.units != _EXPECTED_HFLAV_UNITS:
        raise AnchorError(
            f"{process_id}: {scaffold_anchor.block_key} must use units "
            f"{_EXPECTED_HFLAV_UNITS!r}, got {scaffold_anchor.units!r}"
        )
    if scaffold_anchor.uncertainty is None:
        raise AnchorError(
            f"{process_id}: {scaffold_anchor.block_key} must provide uncertainty"
        )
    value_percent = float(scaffold_anchor.value)
    uncertainty_percent = _positive_float(
        scaffold_anchor.uncertainty,
        process_id=process_id,
        field_name=f"{scaffold_anchor.block_key}.uncertainty",
    )
    return PercentAsymmetryAnchor(
        block_key=scaffold_anchor.block_key,
        source=_optional_str(scaffold_anchor.source),
        year=scaffold_anchor.year,
        value_percent=value_percent,
        uncertainty_percent=float(uncertainty_percent),
        value=float(value_percent / 100.0),
        uncertainty=float(uncertainty_percent / 100.0),
        units=scaffold_anchor.units,
        source_url=_optional_str(scaffold_anchor.source_url),
        snapshot_path=_optional_str(scaffold_anchor.snapshot_path),
        value_summary=_optional_str(sub.get("value_summary")),
    )


def _load_no_cpv_context(process_id: str) -> NoCPVContext:
    pdg_block = _pdg_block(process_id)
    sub = pdg_block.get(_HFLAV_NO_CPV_KEY)
    if not isinstance(sub, Mapping):
        raise AnchorError(f"{process_id}: {_HFLAV_NO_CPV_KEY!r} is not a mapping")
    return NoCPVContext(
        block_key=_HFLAV_NO_CPV_KEY,
        source=_optional_str(sub.get("source")),
        year=_optional_int(sub.get("year")),
        delta_chi2=_required_float(
            sub.get("delta_chi2"),
            process_id=process_id,
            field_name=f"{_HFLAV_NO_CPV_KEY}.delta_chi2",
        ),
        significance_sigma=_required_float(
            sub.get("significance_sigma"),
            process_id=process_id,
            field_name=f"{_HFLAV_NO_CPV_KEY}.significance_sigma",
        ),
        confidence_level_table=_optional_str(sub.get("confidence_level_table")),
        confidence_level_text=_optional_str(sub.get("confidence_level_text")),
        source_url=_optional_str(sub.get("source_url")),
        snapshot_path=_optional_str(sub.get("snapshot_path")),
    )


def _load_cfw_context(process_id: str) -> CFWContext:
    pdg_block = _pdg_block(process_id)
    sub = pdg_block.get(_CFW_CONTEXT_KEY)
    if not isinstance(sub, Mapping):
        raise AnchorError(f"{process_id}: {_CFW_CONTEXT_KEY!r} is not a mapping")
    return CFWContext(
        block_key=_CFW_CONTEXT_KEY,
        source=_optional_str(sub.get("source")),
        year=_optional_int(sub.get("year")),
        value_summary=_optional_str(sub.get("value_summary")),
        source_url=_optional_str(sub.get("source_url")),
        snapshot_path=_optional_str(sub.get("snapshot_path")),
    )


def _load_c003_anchor(process_id: str) -> C003Anchor:
    anchor = C003Anchor(
        lhcb_discovery=_load_lhcb_discovery_anchor(process_id),
        hflav_direct_cpv=_load_percent_anchor(
            _HFLAV_DIRECT_CPV_CANDIDATES,
            process_id=process_id,
        ),
        hflav_indirect_cpv=_load_percent_anchor(
            _HFLAV_INDIRECT_CPV_CANDIDATES,
            process_id=process_id,
        ),
        hflav_no_cpv=_load_no_cpv_context(process_id),
        cfw_context=_load_cfw_context(process_id),
    )
    if anchor.budget <= 0.0:
        raise AnchorError(f"{process_id}: C003 non-vetoing NP room must be positive")
    return anchor


@register
class Constraint:
    """Catalogued non-vetoing direct charm-CP stub (process_id C003)."""

    process_id = "C003"
    severity = Severity.INFO
    observable = "Delta A_CP(D0 -> K+K-, pi+pi-)"

    def __init__(self) -> None:
        self.anchor = _load_c003_anchor(self.process_id)

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        del point
        comparison = compare_delta_acp_np_room_to_measurement(
            measured_delta_acp=self.anchor.value,
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
                "INFO-only C003 stub: loads the LHCb Delta A_CP measurement "
                "and reports the full observed |Delta A_CP| as documented "
                "NP room. No SM or RS penguin amplitude is computed; both "
                "sides are flagged NEEDS-HUMAN-PHYSICS, so this result is "
                "non-vetoing."
            ),
            diagnostics={
                "non_vetoing": True,
                "no_penguin_calculation": True,
                "stub_model": CHARM_DIRECT_CP_STUB_MODEL_V1,
                "budget_source": _BUDGET_SOURCE,
                "budget_is_full_observed_asymmetry": True,
                "parametrization_citation": _PARAMETRIZATION_CITATION,
                "needs_human_physics": _NEEDS_HUMAN_PHYSICS,
                "needs_human_physics_sm": CHARM_DIRECT_CP_SM_NEEDS_HUMAN_PHYSICS,
                "needs_human_physics_np": CHARM_DIRECT_CP_RS_NEEDS_HUMAN_PHYSICS,
                "sm_long_distance_penguin_dominated": True,
                "rs_penguin_matching_available": False,
                "parameter_point_inputs_used": (),
                "experimental_block": self.anchor.lhcb_discovery.block_key,
                "experimental_value_raw_x1e_minus4": float(
                    self.anchor.lhcb_discovery.raw_value
                ),
                "experimental_uncertainty_raw_x1e_minus4": float(
                    self.anchor.lhcb_discovery.raw_uncertainty
                ),
                "experimental_scale": float(self.anchor.lhcb_discovery.scale),
                "experimental_value_percent": float(
                    self.anchor.lhcb_discovery.value_percent
                ),
                "experimental_uncertainty_percent": float(
                    self.anchor.lhcb_discovery.uncertainty_percent
                ),
                "experimental_uncertainty": float(self.anchor.uncertainty),
                "measurement_abs": float(comparison.measurement_abs),
                "documented_np_room_abs": float(comparison.documented_np_room_abs),
                "measurement_to_np_room_ratio": float(comparison.ratio_to_room),
                "lhcb_discovery_significance_sigma_naive": float(
                    abs(self.anchor.lhcb_discovery.raw_value)
                    / self.anchor.lhcb_discovery.raw_uncertainty
                ),
                "hflav_direct_cpv_block": self.anchor.hflav_direct_cpv.block_key,
                "hflav_direct_cpv_value_percent": float(
                    self.anchor.hflav_direct_cpv.value_percent
                ),
                "hflav_direct_cpv_uncertainty_percent": float(
                    self.anchor.hflav_direct_cpv.uncertainty_percent
                ),
                "hflav_indirect_cpv_block": self.anchor.hflav_indirect_cpv.block_key,
                "hflav_indirect_cpv_value_percent": float(
                    self.anchor.hflav_indirect_cpv.value_percent
                ),
                "hflav_indirect_cpv_uncertainty_percent": float(
                    self.anchor.hflav_indirect_cpv.uncertainty_percent
                ),
                "hflav_no_cpv_delta_chi2": float(self.anchor.hflav_no_cpv.delta_chi2),
                "hflav_no_cpv_significance_sigma": float(
                    self.anchor.hflav_no_cpv.significance_sigma
                ),
                "cfw_context_summary": self.anchor.cfw_context.value_summary,
            },
        )
