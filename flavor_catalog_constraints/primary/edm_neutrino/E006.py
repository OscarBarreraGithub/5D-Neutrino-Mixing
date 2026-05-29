"""E006 - Mercury-199 atomic electric dipole moment ``|d_Hg|``.

Physics
-------
The catalogued observable is the direct ``199Hg`` atomic EDM upper limit

    |d_Hg| < 7.4e-30 e cm  (95% CL),

anchored to the Graner et al. measurement in
``flavor_catalog/processes/edm_neutrino/E006.yaml``.  A rigorous RS
prediction is **not** available in this scaffold.  It would require

1. nuclear Schiff-moment and relativistic atomic-structure inputs translating
   CP-odd low-energy operators into the atomic EDM; and
2. RS matching onto CP-odd quark, gluon, and semileptonic operators.

Both ingredients are flagged ``NEEDS-HUMAN-PHYSICS`` in the returned
diagnostics.  This file does not fake a Schiff-moment, atomic-structure, or
RS CP-odd matching calculation and does not read model inputs from
``ParameterPoint``.

Severity
--------
INFO.  The Hg EDM is a stringent observed bound, but applying it as a hard
veto here would require nuclear/atomic response inputs and RS CP-odd matching
that are absent from the current catalog contract.  The returned ``passes``
value is bookkeeping only and must not veto scan points.

Catalog sidecar
---------------
``flavor_catalog/processes/edm_neutrino/E006.yaml`` is the source of truth for
the Graner et al. canonical direct limit, improvement factor, PDG cross
reference, and Sahoo theory-translation context.  Numeric values below are
loaded through the scaffold anchor loader and then unit-converted where the
sidecar stores central values and uncertainties in ``10^-30 e cm``.
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
from flavor_catalog_constraints.physics_adapters.mercury_edm import (
    MERCURY_EDM_RS_MATCHING_NEEDS_HUMAN_PHYSICS,
    MERCURY_EDM_SCHIFF_ATOMIC_NEEDS_HUMAN_PHYSICS,
    MERCURY_EDM_STUB_MODEL_V1,
    compare_mercury_edm_measurement_to_limit,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "edm_neutrino"
_CURRENT_LIMIT_BLOCK = "canonical_direct_limit"
_IMPROVEMENT_BLOCK = "post_2008_improvement_factor"
_PDG_CROSS_REFERENCE_BLOCK = "pdg_cross_reference_neutron"
_THEORY_TRANSLATION_BLOCK = "theory_translation_context"
_LIMIT_ANCHOR_CANDIDATES = (_CURRENT_LIMIT_BLOCK,)
_EXPECTED_LIMIT_UNITS = "e cm"
_EXPECTED_MEASUREMENT_UNITS = "10^-30 e cm"
_EXPECTED_LIMIT_OPERATOR = "<"
_BUDGET_SOURCE = (
    "flavor_catalog/processes/edm_neutrino/E006.yaml canonical_direct_limit "
    "(Graner et al. 2016)"
)
_PARAMETRIZATION_CITATION = (
    "Graner et al. Phys. Rev. Lett. 116, 161601 (2016), arXiv:1601.04339; "
    "Sahoo Phys. Rev. D 95, 013002 (2017), arXiv:1612.09371"
)
_NEEDS_HUMAN_PHYSICS = (
    MERCURY_EDM_SCHIFF_ATOMIC_NEEDS_HUMAN_PHYSICS,
    MERCURY_EDM_RS_MATCHING_NEEDS_HUMAN_PHYSICS,
)


@dataclass(frozen=True)
class MercuryEDMLimitAnchor:
    """Typed current direct upper limit on ``|d_Hg|`` in e cm."""

    experimental: Anchor
    limit_operator: str | None
    confidence_level: str | None
    central_value_e_cm: float
    statistical_uncertainty_e_cm: float
    systematic_uncertainty_e_cm: float
    total_uncertainty_e_cm: float
    measurement_units: str | None
    value_summary: str | None

    @property
    def value(self) -> float:
        """Upper limit on ``|d_Hg|`` in e cm."""
        return self.experimental.value

    @property
    def budget(self) -> float:
        """Non-vetoing bookkeeping budget: the measured Hg EDM upper limit."""
        return self.experimental.value

    @property
    def source_url(self) -> str | None:
        """Primary source URL for the canonical direct limit."""
        return self.experimental.source_url


@dataclass(frozen=True)
class ScalarContext:
    """Typed scalar context block carried for provenance, not calculation."""

    block_key: str
    source: str | None
    year: int | None
    value: float
    units: str | None
    value_summary: str | None
    source_url: str | None
    snapshot_path: str | None


@dataclass(frozen=True)
class PDGCrossReference:
    """Typed PDG neutron-EDM cross-reference derived from Hg theory rows."""

    block_key: str
    source: str | None
    year: int | None
    graner_row_value: float
    sahoo_row_value: float
    table_units: str | None
    confidence_level: str | None
    value_summary: str | None
    source_url: str | None
    snapshot_path: str | None


@dataclass(frozen=True)
class TheoryTranslationContext:
    """Typed Sahoo interpretation-level context, not used as a constraint."""

    block_key: str
    source: str | None
    year: int | None
    neutron_edm_limit: float
    proton_edm_limit: float
    theta_bar_limit: float
    up_down_cedm_difference_limit: float
    edm_units: str | None
    value_summary: str | None
    source_url: str | None
    snapshot_path: str | None


@dataclass(frozen=True)
class Post2008Context:
    """Typed post-2008 context carried only for provenance."""

    key: str | None
    source_url: str | None
    year: int | None
    value_summary: str | None
    snapshot_path: str | None


@dataclass(frozen=True)
class E006Anchor:
    """Typed E006 anchor bundle and non-vetoing budget convention."""

    limit: MercuryEDMLimitAnchor
    improvement_factor: ScalarContext
    pdg_cross_reference: PDGCrossReference
    theory_translation_context: TheoryTranslationContext
    baseline_context: Post2008Context
    direct_hg_experiment_context: Post2008Context
    atomic_nuclear_translation_context: Post2008Context

    @property
    def value(self) -> float:
        """Canonical current upper limit on ``|d_Hg|`` in e cm."""
        return self.limit.value

    @property
    def budget(self) -> float:
        """Non-vetoing bookkeeping budget: the canonical direct limit."""
        return self.limit.budget

    @property
    def source_url(self) -> str | None:
        """Source URL for the canonical direct limit."""
        return self.limit.source_url


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
            f"{process_id}: E006 anchor field {field_name!r}={value!r} "
            "is not numeric"
        ) from exc
    if not math.isfinite(out):
        raise AnchorError(f"{process_id}: E006 anchor field {field_name!r} is not finite")
    return out


def _positive_float(value: Any, *, process_id: str, field_name: str) -> float:
    out = _required_float(value, process_id=process_id, field_name=field_name)
    if out <= 0.0:
        raise AnchorError(f"{process_id}: E006 anchor field {field_name!r} must be positive")
    return out


def _scale_to_e_cm(units: str | None, *, process_id: str, field_name: str) -> float:
    if units == _EXPECTED_LIMIT_UNITS:
        return 1.0
    if units == _EXPECTED_MEASUREMENT_UNITS:
        return 1.0e-30
    raise AnchorError(
        f"{process_id}: unsupported E006 units for {field_name}: {units!r}"
    )


def _pdg_block(process_id: str) -> Mapping[str, Any]:
    return load_pdg_block(process_id, family=_FAMILY)


def _load_limit_anchor(process_id: str) -> MercuryEDMLimitAnchor:
    experimental = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=_LIMIT_ANCHOR_CANDIDATES,
    )
    if experimental.block_key != _CURRENT_LIMIT_BLOCK:
        raise AnchorError(
            f"{process_id}: load_anchor selected {experimental.block_key!r}, "
            f"expected {_CURRENT_LIMIT_BLOCK!r} for the current Hg EDM limit"
        )
    if experimental.value <= 0.0:
        raise AnchorError(f"{process_id}: Hg EDM limit must be positive")
    if experimental.units != _EXPECTED_LIMIT_UNITS:
        raise AnchorError(
            f"{process_id}: expected units {_EXPECTED_LIMIT_UNITS!r}, "
            f"got {experimental.units!r}"
        )

    pdg = _pdg_block(process_id)
    sub = pdg.get(experimental.block_key)
    if not isinstance(sub, Mapping):
        raise AnchorError(
            f"{process_id}: selected anchor block {experimental.block_key!r} "
            "is not available as a mapping"
        )
    limit_operator = _optional_str(sub.get("limit_operator"))
    if limit_operator != _EXPECTED_LIMIT_OPERATOR:
        raise AnchorError(
            f"{process_id}: expected limit_operator {_EXPECTED_LIMIT_OPERATOR!r}, "
            f"got {limit_operator!r}"
        )
    measurement_units = _optional_str(sub.get("measurement_units"))
    measurement_scale = _scale_to_e_cm(
        measurement_units,
        process_id=process_id,
        field_name=f"{experimental.block_key}.measurement_units",
    )
    stat = _positive_float(
        sub.get("statistical_uncertainty"),
        process_id=process_id,
        field_name=f"{experimental.block_key}.statistical_uncertainty",
    )
    syst = _positive_float(
        sub.get("systematic_uncertainty"),
        process_id=process_id,
        field_name=f"{experimental.block_key}.systematic_uncertainty",
    )
    central_value_e_cm = _required_float(
        sub.get("central_value_from_arxiv_v4"),
        process_id=process_id,
        field_name=f"{experimental.block_key}.central_value_from_arxiv_v4",
    ) * measurement_scale

    return MercuryEDMLimitAnchor(
        experimental=experimental,
        limit_operator=limit_operator,
        confidence_level=_optional_str(sub.get("confidence_level")),
        central_value_e_cm=float(central_value_e_cm),
        statistical_uncertainty_e_cm=float(stat * measurement_scale),
        systematic_uncertainty_e_cm=float(syst * measurement_scale),
        total_uncertainty_e_cm=float(math.sqrt(stat * stat + syst * syst) * measurement_scale),
        measurement_units=measurement_units,
        value_summary=_optional_str(sub.get("value_summary")),
    )


def _load_improvement_factor(process_id: str) -> ScalarContext:
    anchor = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=(_IMPROVEMENT_BLOCK,),
    )
    if anchor.value <= 0.0:
        raise AnchorError(f"{process_id}: improvement factor must be positive")
    pdg = _pdg_block(process_id)
    sub = pdg.get(anchor.block_key)
    if not isinstance(sub, Mapping):
        raise AnchorError(f"{process_id}: {_IMPROVEMENT_BLOCK!r} is not a mapping")
    return ScalarContext(
        block_key=anchor.block_key,
        source=_optional_str(anchor.source),
        year=anchor.year,
        value=float(anchor.value),
        units=anchor.units,
        value_summary=_optional_str(sub.get("value_summary")),
        source_url=_optional_str(anchor.source_url),
        snapshot_path=_optional_str(anchor.snapshot_path),
    )


def _load_pdg_cross_reference(process_id: str) -> PDGCrossReference:
    pdg = _pdg_block(process_id)
    sub = pdg.get(_PDG_CROSS_REFERENCE_BLOCK)
    if not isinstance(sub, Mapping):
        raise AnchorError(f"{process_id}: {_PDG_CROSS_REFERENCE_BLOCK!r} is not a mapping")
    return PDGCrossReference(
        block_key=_PDG_CROSS_REFERENCE_BLOCK,
        source=_optional_str(sub.get("source")),
        year=_optional_int(sub.get("year")),
        graner_row_value=_required_float(
            sub.get("graner_row_value"),
            process_id=process_id,
            field_name=f"{_PDG_CROSS_REFERENCE_BLOCK}.graner_row_value",
        ),
        sahoo_row_value=_required_float(
            sub.get("sahoo_row_value"),
            process_id=process_id,
            field_name=f"{_PDG_CROSS_REFERENCE_BLOCK}.sahoo_row_value",
        ),
        table_units=_optional_str(sub.get("table_units")),
        confidence_level=_optional_str(sub.get("confidence_level")),
        value_summary=_optional_str(sub.get("value_summary")),
        source_url=_optional_str(sub.get("source_url")),
        snapshot_path=_optional_str(sub.get("snapshot_path")),
    )


def _load_theory_translation_context(process_id: str) -> TheoryTranslationContext:
    pdg = _pdg_block(process_id)
    sub = pdg.get(_THEORY_TRANSLATION_BLOCK)
    if not isinstance(sub, Mapping):
        raise AnchorError(f"{process_id}: {_THEORY_TRANSLATION_BLOCK!r} is not a mapping")
    return TheoryTranslationContext(
        block_key=_THEORY_TRANSLATION_BLOCK,
        source=_optional_str(sub.get("source")),
        year=_optional_int(sub.get("year")),
        neutron_edm_limit=_required_float(
            sub.get("neutron_edm_limit"),
            process_id=process_id,
            field_name=f"{_THEORY_TRANSLATION_BLOCK}.neutron_edm_limit",
        ),
        proton_edm_limit=_required_float(
            sub.get("proton_edm_limit"),
            process_id=process_id,
            field_name=f"{_THEORY_TRANSLATION_BLOCK}.proton_edm_limit",
        ),
        theta_bar_limit=_required_float(
            sub.get("theta_bar_limit"),
            process_id=process_id,
            field_name=f"{_THEORY_TRANSLATION_BLOCK}.theta_bar_limit",
        ),
        up_down_cedm_difference_limit=_required_float(
            sub.get("up_down_cedm_difference_limit"),
            process_id=process_id,
            field_name=f"{_THEORY_TRANSLATION_BLOCK}.up_down_cedm_difference_limit",
        ),
        edm_units=_optional_str(sub.get("edm_units")),
        value_summary=_optional_str(sub.get("value_summary")),
        source_url=_optional_str(sub.get("source_url")),
        snapshot_path=_optional_str(sub.get("snapshot_path")),
    )


def _load_post_2008_context(process_id: str, key: str) -> Post2008Context:
    data = load_full_yaml(process_id, family=_FAMILY)
    post_context = data.get("post_2008_context")
    if not isinstance(post_context, Mapping):
        raise AnchorError(f"{process_id}: missing post_2008_context mapping")
    sub = post_context.get(key)
    if not isinstance(sub, Mapping):
        raise AnchorError(f"{process_id}: post_2008_context.{key!r} is not a mapping")
    return Post2008Context(
        key=_optional_str(sub.get("key")),
        source_url=_optional_str(sub.get("source_url")),
        year=_optional_int(sub.get("year")),
        value_summary=_optional_str(sub.get("value_summary")),
        snapshot_path=_optional_str(sub.get("snapshot_path")),
    )


def _load_e006_anchor(process_id: str) -> E006Anchor:
    return E006Anchor(
        limit=_load_limit_anchor(process_id),
        improvement_factor=_load_improvement_factor(process_id),
        pdg_cross_reference=_load_pdg_cross_reference(process_id),
        theory_translation_context=_load_theory_translation_context(process_id),
        baseline_context=_load_post_2008_context(process_id, "baseline"),
        direct_hg_experiment_context=_load_post_2008_context(
            process_id,
            "direct_hg_experiment",
        ),
        atomic_nuclear_translation_context=_load_post_2008_context(
            process_id,
            "atomic_nuclear_translation",
        ),
    )


@register
class Constraint:
    """Catalogued non-vetoing Mercury-199 EDM stub."""

    process_id = "E006"
    severity = Severity.INFO
    observable = "|d_Hg|"

    def __init__(self) -> None:
        self.anchor = _load_e006_anchor(self.process_id)

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        del point
        comparison = compare_mercury_edm_measurement_to_limit(
            measured_mercury_edm_e_cm=self.anchor.limit.central_value_e_cm,
            experimental_limit_e_cm=self.anchor.budget,
        )

        return ConstraintResult(
            process_id=self.process_id,
            severity=self.severity,
            passes=bool(comparison.passes),
            predicted=None,
            sm_prediction=None,
            experimental=float(self.anchor.value),
            ratio=float(comparison.ratio_to_limit),
            budget=float(comparison.experimental_limit_e_cm),
            notes=(
                "INFO-only E006 stub: loads the direct 199Hg EDM Graner et al. "
                "limit and records the observable. No Schiff-moment, "
                "atomic-structure, or RS CP-odd matching calculation is "
                "performed; both required physics inputs are flagged "
                "NEEDS-HUMAN-PHYSICS, so this result is non-vetoing."
            ),
            diagnostics={
                "non_vetoing": True,
                "edm_prediction_evaluated": False,
                "no_schiff_moment_calculation": True,
                "no_atomic_structure_calculation": True,
                "no_rs_cp_odd_matching": True,
                "stub_model": MERCURY_EDM_STUB_MODEL_V1,
                "budget_source": _BUDGET_SOURCE,
                "budget_is_current_mercury_edm_limit": True,
                "parametrization_citation": _PARAMETRIZATION_CITATION,
                "needs_human_physics": _NEEDS_HUMAN_PHYSICS,
                "needs_human_physics_schiff_atomic_structure": (
                    MERCURY_EDM_SCHIFF_ATOMIC_NEEDS_HUMAN_PHYSICS
                ),
                "needs_human_physics_rs_cp_odd_matching": (
                    MERCURY_EDM_RS_MATCHING_NEEDS_HUMAN_PHYSICS
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
                "experimental_block": self.anchor.limit.experimental.block_key,
                "experimental_limit_e_cm": float(self.anchor.value),
                "measurement_central_e_cm": float(
                    self.anchor.limit.central_value_e_cm
                ),
                "measurement_abs_e_cm": float(comparison.measurement_abs_e_cm),
                "measurement_to_limit_ratio": float(comparison.ratio_to_limit),
                "statistical_uncertainty_e_cm": float(
                    self.anchor.limit.statistical_uncertainty_e_cm
                ),
                "systematic_uncertainty_e_cm": float(
                    self.anchor.limit.systematic_uncertainty_e_cm
                ),
                "total_uncertainty_e_cm": float(
                    self.anchor.limit.total_uncertainty_e_cm
                ),
                "measurement_units": self.anchor.limit.measurement_units,
                "limit_operator": self.anchor.limit.limit_operator,
                "confidence_level": self.anchor.limit.confidence_level,
                "value_summary": self.anchor.limit.value_summary,
                "improvement_factor": float(self.anchor.improvement_factor.value),
                "improvement_factor_units": self.anchor.improvement_factor.units,
                "pdg_cross_reference_graner_row_value": float(
                    self.anchor.pdg_cross_reference.graner_row_value
                ),
                "pdg_cross_reference_sahoo_row_value": float(
                    self.anchor.pdg_cross_reference.sahoo_row_value
                ),
                "pdg_cross_reference_table_units": (
                    self.anchor.pdg_cross_reference.table_units
                ),
                "theory_translation_neutron_edm_limit": float(
                    self.anchor.theory_translation_context.neutron_edm_limit
                ),
                "theory_translation_proton_edm_limit": float(
                    self.anchor.theory_translation_context.proton_edm_limit
                ),
                "theory_translation_theta_bar_limit": float(
                    self.anchor.theory_translation_context.theta_bar_limit
                ),
                "theory_translation_up_down_cedm_difference_limit": float(
                    self.anchor
                    .theory_translation_context
                    .up_down_cedm_difference_limit
                ),
                "baseline_context_key": self.anchor.baseline_context.key,
                "direct_hg_experiment_context_key": (
                    self.anchor.direct_hg_experiment_context.key
                ),
                "atomic_nuclear_translation_context_key": (
                    self.anchor.atomic_nuclear_translation_context.key
                ),
                "atomic_nuclear_translation_context_summary": (
                    self.anchor.atomic_nuclear_translation_context.value_summary
                ),
            },
        )
