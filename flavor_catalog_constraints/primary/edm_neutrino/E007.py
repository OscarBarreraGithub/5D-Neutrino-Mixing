"""E007 - Radium-225 and Xenon-129 atomic EDMs ``d_Ra, d_Xe``.

Physics
-------
The catalogued observables are direct atomic EDM limits for ``225Ra`` and
``129Xe``, anchored to the experiment-paper rows in
``flavor_catalog/processes/edm_neutrino/E007.yaml``:

    |d(225Ra)|     < 1.4e-23 e cm  (95% CL),
    |d_A(129Xe)|   < 1.4e-27 e cm  (95% CL).

A rigorous RS prediction is **not** available in this scaffold.  It would
require

1. separate nuclear Schiff-moment and relativistic atomic-structure inputs for
   Ra and Xe; and
2. RS matching onto CP-odd quark, gluon, and semileptonic operators.

Both ingredients are flagged ``NEEDS-HUMAN-PHYSICS`` in the returned
diagnostics.  This file does not fake a Schiff-moment, atomic-structure, or
RS CP-odd matching calculation and does not read model inputs from
``ParameterPoint``.

Severity
--------
INFO.  The Ra/Xe EDM limits are observed bounds, but applying them as a hard
veto here would require nuclear/atomic response inputs and RS CP-odd matching
that are absent from the current catalog contract.  The returned result is
non-vetoing and records the direct limits only.

Catalog sidecar
---------------
``flavor_catalog/processes/edm_neutrino/E007.yaml`` is the source of truth for
the Bishof ``225Ra`` direct limit, the Sachdeva ``129Xe`` direct limit, the
independent Heidelberg ``129Xe`` cross-check, and the contextual Ra/Xe program
metadata.  Numeric values below are loaded through the scaffold anchor loader.
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
    load_pdg_block,
)
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.atomic_edm import (
    ATOMIC_EDM_RS_MATCHING_NEEDS_HUMAN_PHYSICS,
    ATOMIC_EDM_SCHIFF_ATOMIC_NEEDS_HUMAN_PHYSICS,
    ATOMIC_EDM_STUB_MODEL_V1,
    atomic_edm_direct_limit,
    compare_atomic_edm_to_limit,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "edm_neutrino"
_RA225_CURRENT_BLOCK = "ra225_current_direct_limit"
_RA225_FIRST_BLOCK = "ra225_first_measurement"
_XE129_BEST_BLOCK = "xe129_best_direct_limit"
_XE129_INDEPENDENT_BLOCK = "xe129_independent_heidelberg_limit"
_EXPECTED_UNITS = "e cm"
_EXPECTED_LIMIT_OPERATOR = "<"
_BUDGET_SOURCE = (
    "flavor_catalog/processes/edm_neutrino/E007.yaml strongest current direct "
    "atomic-EDM limit among ra225_current_direct_limit and xe129_best_direct_limit"
)
_PARAMETRIZATION_CITATION = (
    "Bishof et al. Phys. Rev. C 94, 025501 (2016), arXiv:1606.04931; "
    "Sachdeva et al. Phys. Rev. Lett. 123, 143003 (2019), arXiv:1909.12800; "
    "Allmendinger et al. Phys. Rev. A 100, 022505 (2019), arXiv:1904.12295"
)
_NEEDS_HUMAN_PHYSICS = (
    ATOMIC_EDM_SCHIFF_ATOMIC_NEEDS_HUMAN_PHYSICS,
    ATOMIC_EDM_RS_MATCHING_NEEDS_HUMAN_PHYSICS,
)


@dataclass(frozen=True)
class AtomicEDMLimitAnchor:
    """Typed direct atomic-EDM limit in ``e cm``."""

    isotope: str
    experimental: Anchor
    limit_operator: str | None
    confidence_level: str | None
    central_value_e_cm: float | None
    statistical_uncertainty_e_cm: float | None
    systematic_uncertainty_e_cm: float | None
    total_uncertainty_e_cm: float | None
    improvement_factor: float | None
    value_summary: str | None

    @property
    def block_key(self) -> str:
        """Catalog block selected by the scaffold anchor loader."""
        return self.experimental.block_key

    @property
    def value(self) -> float:
        """Upper limit on the isotope's atomic EDM in ``e cm``."""
        return self.experimental.value

    @property
    def source_url(self) -> str | None:
        """Primary source URL for the direct limit."""
        return self.experimental.source_url

    @property
    def source(self) -> str | None:
        """Primary source label for the direct limit."""
        return self.experimental.source


@dataclass(frozen=True)
class E007Context:
    """Post-2008 context carried for provenance, not calculation."""

    key: str | None
    source_url: str | None
    year: int | None
    year_range: str | None
    value_summary: str | None
    snapshot_path: str | None
    numeric_values: Mapping[str, float]


@dataclass(frozen=True)
class E007Anchor:
    """Typed E007 anchor bundle and non-vetoing budget convention."""

    ra225_current_direct_limit: AtomicEDMLimitAnchor
    ra225_first_measurement: AtomicEDMLimitAnchor
    xe129_best_direct_limit: AtomicEDMLimitAnchor
    xe129_independent_heidelberg_limit: AtomicEDMLimitAnchor
    baseline_context: E007Context
    argonne_ra225_program_context: E007Context
    ra225_future_sensitivity_context: E007Context
    heidelberg_xe_upgrade_context: E007Context
    kuleuven_raf_acf_molecule_context: E007Context

    @property
    def reference_limit(self) -> AtomicEDMLimitAnchor:
        """Strongest current direct Ra/Xe atomic-EDM limit in ``e cm``."""
        return min(
            (self.ra225_current_direct_limit, self.xe129_best_direct_limit),
            key=lambda limit: limit.value,
        )

    @property
    def value(self) -> float:
        """Selected scalar reference limit in ``e cm`` for ConstraintResult."""
        return self.reference_limit.value

    @property
    def budget(self) -> float:
        """Non-vetoing bookkeeping budget: strongest current direct limit."""
        return self.reference_limit.value

    @property
    def source_url(self) -> str | None:
        """Source URL for the selected scalar reference limit."""
        return self.reference_limit.source_url


def _optional_str(value: Any) -> str | None:
    return None if value is None else str(value)


def _optional_int(value: Any) -> int | None:
    if value is None:
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _optional_float(value: Any, *, process_id: str, field_name: str) -> float | None:
    if value is None:
        return None
    return _required_float(value, process_id=process_id, field_name=field_name)


def _required_float(value: Any, *, process_id: str, field_name: str) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError) as exc:
        raise AnchorError(
            f"{process_id}: E007 anchor field {field_name!r}={value!r} "
            "is not numeric"
        ) from exc
    if not math.isfinite(out):
        raise AnchorError(f"{process_id}: E007 anchor field {field_name!r} is not finite")
    return out


def _positive_float(value: Any, *, process_id: str, field_name: str) -> float:
    out = _required_float(value, process_id=process_id, field_name=field_name)
    if out <= 0.0:
        raise AnchorError(f"{process_id}: E007 anchor field {field_name!r} must be positive")
    return out


def _pdg_block(process_id: str) -> Mapping[str, Any]:
    return load_pdg_block(process_id, family=_FAMILY)


def _subblock_for_anchor(
    process_id: str,
    anchor: Anchor,
) -> Mapping[str, Any]:
    pdg = _pdg_block(process_id)
    sub = pdg.get(anchor.block_key)
    if not isinstance(sub, Mapping):
        raise AnchorError(
            f"{process_id}: selected anchor block {anchor.block_key!r} "
            "is not available as a mapping"
        )
    return sub


def _load_atomic_limit_anchor(
    process_id: str,
    *,
    block_key: str,
    isotope: str,
    value_key: str,
) -> AtomicEDMLimitAnchor:
    experimental = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=(block_key,),
        value_key=value_key,
    )
    if experimental.block_key != block_key:
        raise AnchorError(
            f"{process_id}: load_anchor selected {experimental.block_key!r}, "
            f"expected {block_key!r} for an E007 direct limit"
        )
    if experimental.value <= 0.0:
        raise AnchorError(f"{process_id}: {block_key} limit must be positive")
    if experimental.units != _EXPECTED_UNITS:
        raise AnchorError(
            f"{process_id}: expected units {_EXPECTED_UNITS!r} for {block_key}, "
            f"got {experimental.units!r}"
        )

    sub = _subblock_for_anchor(process_id, experimental)
    limit_operator = _optional_str(sub.get("limit_operator"))
    if limit_operator != _EXPECTED_LIMIT_OPERATOR:
        raise AnchorError(
            f"{process_id}: expected limit_operator {_EXPECTED_LIMIT_OPERATOR!r} "
            f"for {block_key}, got {limit_operator!r}"
        )

    stat = _optional_float(
        sub.get("statistical_uncertainty"),
        process_id=process_id,
        field_name=f"{block_key}.statistical_uncertainty",
    )
    syst = _optional_float(
        sub.get("systematic_uncertainty"),
        process_id=process_id,
        field_name=f"{block_key}.systematic_uncertainty",
    )
    total = _optional_float(
        sub.get("total_uncertainty"),
        process_id=process_id,
        field_name=f"{block_key}.total_uncertainty",
    )
    if total is None and stat is not None and syst is not None:
        total = math.sqrt(stat * stat + syst * syst)

    return AtomicEDMLimitAnchor(
        isotope=isotope,
        experimental=experimental,
        limit_operator=limit_operator,
        confidence_level=_optional_str(sub.get("confidence_level")),
        central_value_e_cm=_optional_float(
            sub.get("central_value"),
            process_id=process_id,
            field_name=f"{block_key}.central_value",
        ),
        statistical_uncertainty_e_cm=stat,
        systematic_uncertainty_e_cm=syst,
        total_uncertainty_e_cm=total,
        improvement_factor=_optional_float(
            sub.get("improvement_factor"),
            process_id=process_id,
            field_name=f"{block_key}.improvement_factor",
        ),
        value_summary=_optional_str(sub.get("value_summary")),
    )


def _load_context(process_id: str, key: str) -> E007Context:
    data = load_full_yaml(process_id, family=_FAMILY)
    post_context = data.get("post_2008_context")
    if not isinstance(post_context, Mapping):
        raise AnchorError(f"{process_id}: missing post_2008_context mapping")
    sub = post_context.get(key)
    if not isinstance(sub, Mapping):
        raise AnchorError(f"{process_id}: post_2008_context.{key!r} is not a mapping")
    numeric_values = {
        str(name): float(value)
        for name, value in sub.items()
        if isinstance(value, (int, float)) and not isinstance(value, bool)
    }
    return E007Context(
        key=_optional_str(sub.get("key")),
        source_url=_optional_str(sub.get("source_url")),
        year=_optional_int(sub.get("year")),
        year_range=_optional_str(sub.get("year_range")),
        value_summary=_optional_str(sub.get("value_summary")),
        snapshot_path=_optional_str(sub.get("snapshot_path")),
        numeric_values=numeric_values,
    )


def _load_e007_anchor(process_id: str) -> E007Anchor:
    return E007Anchor(
        ra225_current_direct_limit=_load_atomic_limit_anchor(
            process_id,
            block_key=_RA225_CURRENT_BLOCK,
            isotope="225Ra",
            value_key="value",
        ),
        ra225_first_measurement=_load_atomic_limit_anchor(
            process_id,
            block_key=_RA225_FIRST_BLOCK,
            isotope="225Ra",
            value_key="value",
        ),
        xe129_best_direct_limit=_load_atomic_limit_anchor(
            process_id,
            block_key=_XE129_BEST_BLOCK,
            isotope="129Xe",
            value_key="limit_value",
        ),
        xe129_independent_heidelberg_limit=_load_atomic_limit_anchor(
            process_id,
            block_key=_XE129_INDEPENDENT_BLOCK,
            isotope="129Xe",
            value_key="limit_value",
        ),
        baseline_context=_load_context(process_id, "baseline"),
        argonne_ra225_program_context=_load_context(
            process_id,
            "argonne_ra225_program_context",
        ),
        ra225_future_sensitivity_context=_load_context(
            process_id,
            "ra225_future_sensitivity_from_bishof",
        ),
        heidelberg_xe_upgrade_context=_load_context(
            process_id,
            "heidelberg_xe_upgrade",
        ),
        kuleuven_raf_acf_molecule_context=_load_context(
            process_id,
            "kuleuven_raf_acf_molecule_context",
        ),
    )


@register
class Constraint:
    """Catalogued non-vetoing Ra-225/Xe-129 atomic EDM stub."""

    process_id = "E007"
    severity = Severity.INFO
    observable = "|d_Ra|, |d_Xe|"

    def __init__(self) -> None:
        self.anchor = _load_e007_anchor(self.process_id)

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        del point
        reference = self.anchor.reference_limit
        direct_limit = atomic_edm_direct_limit(
            isotope=reference.isotope,
            observable=reference.experimental.observable or self.observable,
            experimental_limit_e_cm=self.anchor.budget,
        )
        comparison = None
        if reference.central_value_e_cm is not None:
            comparison = compare_atomic_edm_to_limit(
                atomic_edm_e_cm=reference.central_value_e_cm,
                experimental_limit_e_cm=reference.value,
                isotope=reference.isotope,
            )
        ratio = 0.0 if comparison is None else float(comparison.ratio_to_limit)
        selected_abs_central = (
            None if comparison is None else float(comparison.abs_atomic_edm_e_cm)
        )

        return ConstraintResult(
            process_id=self.process_id,
            severity=self.severity,
            passes=True,
            predicted=None,
            sm_prediction=None,
            experimental=float(direct_limit.experimental_limit_e_cm),
            ratio=float(ratio),
            budget=float(direct_limit.experimental_limit_e_cm),
            notes=(
                "INFO-only E007 stub: loads the direct 225Ra and 129Xe atomic "
                "EDM limits and records the strongest current direct limit. "
                "No Schiff-moment, atomic-structure, or RS CP-odd matching "
                "calculation is performed; both required physics inputs are "
                "flagged NEEDS-HUMAN-PHYSICS, so this result is non-vetoing."
            ),
            diagnostics={
                "non_vetoing": True,
                "edm_prediction_evaluated": False,
                "atomic_edm_prediction_evaluated": False,
                "no_schiff_moment_calculation": True,
                "no_atomic_structure_calculation": True,
                "no_rs_cp_odd_matching": True,
                "stub_model": ATOMIC_EDM_STUB_MODEL_V1,
                "budget_source": _BUDGET_SOURCE,
                "budget_is_strongest_current_direct_ra_xe_limit": True,
                "parametrization_citation": _PARAMETRIZATION_CITATION,
                "needs_human_physics": _NEEDS_HUMAN_PHYSICS,
                "needs_human_physics_schiff_atomic_structure": (
                    ATOMIC_EDM_SCHIFF_ATOMIC_NEEDS_HUMAN_PHYSICS
                ),
                "needs_human_physics_rs_cp_odd_matching": (
                    ATOMIC_EDM_RS_MATCHING_NEEDS_HUMAN_PHYSICS
                ),
                "schiff_moment_atomic_structure_available": False,
                "rs_cp_odd_matching_available": False,
                "parameter_point_inputs_used": (),
                "required_low_energy_operators": (
                    "quark_edm",
                    "quark_chromo_edm",
                    "weinberg_three_gluon",
                    "cp_odd_semileptonic",
                    "cp_odd_four_quark",
                ),
                "selected_limit_block": reference.block_key,
                "selected_limit_isotope": reference.isotope,
                "selected_limit_e_cm": float(reference.value),
                "selected_limit_source_url": reference.source_url,
                "selected_measurement_central_e_cm": reference.central_value_e_cm,
                "selected_measurement_abs_e_cm": selected_abs_central,
                "selected_measurement_to_limit_ratio": float(ratio),
                "ra225_current_limit_block": (
                    self.anchor.ra225_current_direct_limit.block_key
                ),
                "ra225_current_limit_e_cm": float(
                    self.anchor.ra225_current_direct_limit.value
                ),
                "ra225_current_improvement_factor": float(
                    self.anchor.ra225_current_direct_limit.improvement_factor or 0.0
                ),
                "ra225_first_measurement_limit_e_cm": float(
                    self.anchor.ra225_first_measurement.value
                ),
                "xe129_best_limit_block": self.anchor.xe129_best_direct_limit.block_key,
                "xe129_best_limit_e_cm": float(
                    self.anchor.xe129_best_direct_limit.value
                ),
                "xe129_best_central_value_e_cm": float(
                    self.anchor.xe129_best_direct_limit.central_value_e_cm or 0.0
                ),
                "xe129_best_statistical_uncertainty_e_cm": float(
                    self.anchor.xe129_best_direct_limit.statistical_uncertainty_e_cm
                    or 0.0
                ),
                "xe129_best_systematic_uncertainty_e_cm": float(
                    self.anchor.xe129_best_direct_limit.systematic_uncertainty_e_cm
                    or 0.0
                ),
                "xe129_best_total_uncertainty_e_cm": float(
                    self.anchor.xe129_best_direct_limit.total_uncertainty_e_cm
                    or 0.0
                ),
                "xe129_best_measurement_to_limit_ratio": float(
                    abs(self.anchor.xe129_best_direct_limit.central_value_e_cm or 0.0)
                    / self.anchor.xe129_best_direct_limit.value
                ),
                "xe129_best_improvement_factor": float(
                    self.anchor.xe129_best_direct_limit.improvement_factor or 0.0
                ),
                "xe129_independent_limit_e_cm": float(
                    self.anchor.xe129_independent_heidelberg_limit.value
                ),
                "xe129_independent_central_value_e_cm": float(
                    self.anchor
                    .xe129_independent_heidelberg_limit
                    .central_value_e_cm
                    or 0.0
                ),
                "xe129_independent_total_uncertainty_e_cm": float(
                    self.anchor
                    .xe129_independent_heidelberg_limit
                    .total_uncertainty_e_cm
                    or 0.0
                ),
                "limit_operator": reference.limit_operator,
                "confidence_level": reference.confidence_level,
                "value_summary": reference.value_summary,
                "argonne_ra225_half_life_days": float(
                    self.anchor
                    .argonne_ra225_program_context
                    .numeric_values
                    .get("half_life_days", 0.0)
                ),
                "ra225_projected_statistical_sensitivity_e_cm": float(
                    self.anchor
                    .ra225_future_sensitivity_context
                    .numeric_values
                    .get("projected_statistical_sensitivity", 0.0)
                ),
                "ra225_projected_systematic_uncertainty_e_cm": float(
                    self.anchor
                    .ra225_future_sensitivity_context
                    .numeric_values
                    .get("projected_systematic_uncertainty", 0.0)
                ),
                "heidelberg_xe_upgrade_starting_limit_e_cm": float(
                    self.anchor
                    .heidelberg_xe_upgrade_context
                    .numeric_values
                    .get("starting_limit_value", 0.0)
                ),
                "baseline_context_key": self.anchor.baseline_context.key,
                "argonne_ra225_program_context_key": (
                    self.anchor.argonne_ra225_program_context.key
                ),
                "heidelberg_xe_upgrade_context_key": (
                    self.anchor.heidelberg_xe_upgrade_context.key
                ),
                "kuleuven_raf_acf_molecule_context_key": (
                    self.anchor.kuleuven_raf_acf_molecule_context.key
                ),
            },
        )
