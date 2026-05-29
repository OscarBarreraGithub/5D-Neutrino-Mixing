"""E008 - quark chromo-electric dipole moment reference bounds.

Physics
-------
The catalogued E008 entries are low-energy qCEDM translation benchmarks:

    |tilde d_d + 0.5 tilde d_u| < 1.6e-26 cm
    |tilde d_u - tilde d_d|     < 1.1e-27 cm

derived in ``flavor_catalog/processes/edm_neutrino/E008.yaml`` from neutron
and mercury EDM anchors.  A rigorous RS prediction is **not** available in
this scaffold.  It would require

1. RS CP-odd chromo-dipole matching onto quark cEDM Wilson coefficients; and
2. hadronic/nuclear matrix elements translating qCEDMs into neutron and
   diamagnetic-atom EDMs.

Both ingredients are flagged ``NEEDS-HUMAN-PHYSICS`` in the returned
diagnostics.  This file records the catalogued qCEDM anchors only; it does not
fake a hadronic calculation and does not read quark-sector inputs from
``ParameterPoint``.

Severity
--------
INFO.  These are stringent observed EDM-derived bounds, but applying them as
a scan veto would require the missing hadronic/nuclear inputs and RS qCEDM
matching.  The returned ``passes`` value is bookkeeping only and must not veto
scan points.

Catalog sidecar
---------------
``flavor_catalog/processes/edm_neutrino/E008.yaml`` is the source of truth for
the neutron/Hg experimental anchors and the qCEDM translation bounds.  Numeric
values below are loaded through the scaffold anchor helpers and converted
where the sidecar stores EDM values in scaled ``e cm`` units.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Mapping

from flavor_catalog_constraints import anchors as anchor_scaffold
from flavor_catalog_constraints.anchors import (
    Anchor,
    AnchorError,
    find_block,
    load_anchor,
    load_pdg_block,
)
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.quark_cedm import (
    QUARK_CEDM_HADRONIC_NEEDS_HUMAN_PHYSICS,
    QUARK_CEDM_RS_MATCHING_NEEDS_HUMAN_PHYSICS,
    QUARK_CEDM_STUB_MODEL_V1,
    quark_cedm_reference_bound,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "edm_neutrino"
_EXPERIMENTAL_ANCHORS_BLOCK = "experimental_anchors"
_QCEDM_TRANSLATIONS_BLOCK = "qcedm_translations"
_POST_2008_CONTEXT_BLOCK = "post_2008_context"
_NEUTRON_EDM_BLOCK = "neutron_edm"
_MERCURY_EDM_BLOCK = "mercury_199_edm"
_NEUTRON_QCEDM_BLOCK = "neutron_combination"
_MERCURY_QCEDM_BLOCK = "mercury_isovector_combination"
_NEUTRON_QCEDM_COEFFICIENT_KEY = "coefficient_for_qcedm_combination"
_MERCURY_QCEDM_COEFFICIENT_KEY = "coefficient_for_isovector_qcedm"
_EXPECTED_BOUND_OPERATOR = "<"
_EXPECTED_QCEDM_UNITS = "cm"
_ROUNDING_REL_TOL = 0.10
_BUDGET_SOURCE = (
    "flavor_catalog/processes/edm_neutrino/E008.yaml "
    "qcedm_translations.mercury_isovector_combination"
)
_PARAMETRIZATION_CITATION = (
    "Pospelov and Ritz Phys. Rev. D 63, 073015 (2001), arXiv:hep-ph/0010037; "
    "Olive, Pospelov, Ritz, and Santoso Phys. Rev. D 72, 075001 (2005), "
    "arXiv:hep-ph/0506106"
)
_NEEDS_HUMAN_PHYSICS = (
    QUARK_CEDM_HADRONIC_NEEDS_HUMAN_PHYSICS,
    QUARK_CEDM_RS_MATCHING_NEEDS_HUMAN_PHYSICS,
)


@dataclass(frozen=True)
class ExperimentalEDMLimitAnchor:
    """Typed neutron/mercury EDM limit anchor, converted to e cm."""

    anchor: Anchor
    parent_block_key: str
    nested_block_key: str
    limit_value_e_cm: float
    limit_value_raw: float
    limit_units: str | None
    limit_operator: str | None
    confidence_level: str | None
    value_summary: str | None
    table_value: float | None
    table_units: str | None

    @property
    def source_url(self) -> str | None:
        """Primary source URL for the experimental EDM anchor."""
        return self.anchor.source_url


@dataclass(frozen=True)
class QuarkCEDMTranslationBound:
    """Typed qCEDM translation bound from the E008 sidecar."""

    anchor: Anchor
    parent_block_key: str
    nested_block_key: str
    coefficient_key: str
    coefficient_value: float
    coefficient_uncertainty_summary: str | None
    derived_bound_operator: str | None
    derived_bound_observable: str | None
    derived_bound_units: str | None
    formula: str | None
    assumptions: str | None
    derivation: str | None

    @property
    def block_key(self) -> str:
        """Dotted sidecar block path for diagnostics."""
        return f"{self.parent_block_key}.{self.nested_block_key}"

    @property
    def value(self) -> float:
        """Upper bound on the qCEDM combination in cm."""
        return self.anchor.value

    @property
    def source_url(self) -> str | None:
        """Primary source URL for the translation formula."""
        return self.anchor.source_url


@dataclass(frozen=True)
class E008Context:
    """Post-2008 context carried for provenance, not calculation."""

    block_key: str
    source: str | None
    year: int | None
    source_url: str | None
    value_summary: str | None
    values: tuple[Mapping[str, Any], ...]
    snapshot_path: str | None


@dataclass(frozen=True)
class E008Anchor:
    """Typed E008 anchor bundle and non-vetoing budget convention."""

    neutron_edm_limit: ExperimentalEDMLimitAnchor
    mercury_edm_limit: ExperimentalEDMLimitAnchor
    neutron_combination: QuarkCEDMTranslationBound
    mercury_isovector_combination: QuarkCEDMTranslationBound
    neutron_unrounded_bound_cm: float
    mercury_unrounded_bound_cm: float
    cfw_context: E008Context
    composite_dipoles_context: E008Context
    lattice_qcedm_context: E008Context

    @property
    def reference_bound(self) -> QuarkCEDMTranslationBound:
        """Most stringent catalogued qCEDM translation bound."""
        return min(
            (self.neutron_combination, self.mercury_isovector_combination),
            key=lambda bound: bound.value,
        )

    @property
    def value(self) -> float:
        """Most stringent catalogued qCEDM reference bound in cm."""
        return self.reference_bound.value

    @property
    def budget(self) -> float:
        """Non-vetoing bookkeeping budget: the strongest qCEDM bound."""
        return self.reference_bound.value

    @property
    def source_url(self) -> str | None:
        """Source URL for the selected reference qCEDM bound."""
        return self.reference_bound.source_url


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
            f"{process_id}: E008 anchor field {field_name!r}={value!r} "
            "is not numeric"
        ) from exc
    if not math.isfinite(out):
        raise AnchorError(f"{process_id}: E008 anchor field {field_name!r} is not finite")
    return out


def _positive_float(value: Any, *, process_id: str, field_name: str) -> float:
    out = _required_float(value, process_id=process_id, field_name=field_name)
    if out <= 0.0:
        raise AnchorError(f"{process_id}: E008 anchor field {field_name!r} must be positive")
    return out


def _scale_to_e_cm(units: str | None, *, process_id: str, field_name: str) -> float:
    if units == "e cm":
        return 1.0
    if units == "10^-26 e cm":
        return 1.0e-26
    if units == "10^-30 e cm":
        return 1.0e-30
    raise AnchorError(f"{process_id}: unsupported E008 EDM units for {field_name}: {units!r}")


def _load_nested_mapping(
    process_id: str,
    *,
    parent_key: str,
    nested_key: str,
) -> Mapping[str, Any]:
    pdg = load_pdg_block(process_id, family=_FAMILY)
    parent = find_block(pdg, (parent_key,), process_id=process_id)
    return find_block(
        parent,
        (nested_key,),
        process_id=f"{process_id}.{parent_key}",
    )


def _load_nested_anchor(
    process_id: str,
    *,
    parent_key: str,
    nested_key: str,
    value_key: str,
) -> tuple[Anchor, Mapping[str, Any]]:
    sub = _load_nested_mapping(
        process_id,
        parent_key=parent_key,
        nested_key=nested_key,
    )
    block_key = f"{parent_key}.{nested_key}"
    virtual_block = {block_key: dict(sub)}
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
        anchor = load_anchor(
            process_id,
            family=_FAMILY,
            candidates=(block_key,),
            value_key=value_key,
        )
    finally:
        anchor_scaffold.load_pdg_block = original_load_pdg_block

    if anchor.block_key != block_key:
        raise AnchorError(
            f"{process_id}: load_anchor selected {anchor.block_key!r}, "
            f"expected {block_key!r} for E008 anchor"
        )
    return anchor, sub


def _load_experimental_limit(
    process_id: str,
    *,
    nested_key: str,
    value_key: str,
) -> ExperimentalEDMLimitAnchor:
    anchor, sub = _load_nested_anchor(
        process_id,
        parent_key=_EXPERIMENTAL_ANCHORS_BLOCK,
        nested_key=nested_key,
        value_key=value_key,
    )
    limit_units = _optional_str(sub.get("limit_units", sub.get("units")))
    limit_value_raw = _positive_float(
        anchor.value,
        process_id=process_id,
        field_name=f"{anchor.block_key}.{value_key}",
    )
    scale = _scale_to_e_cm(
        limit_units,
        process_id=process_id,
        field_name=f"{anchor.block_key}.limit_units",
    )
    limit_operator = _optional_str(sub.get("limit_operator"))
    if limit_operator != _EXPECTED_BOUND_OPERATOR:
        raise AnchorError(
            f"{process_id}: expected limit_operator {_EXPECTED_BOUND_OPERATOR!r} "
            f"for {anchor.block_key}, got {limit_operator!r}"
        )
    return ExperimentalEDMLimitAnchor(
        anchor=anchor,
        parent_block_key=_EXPERIMENTAL_ANCHORS_BLOCK,
        nested_block_key=nested_key,
        limit_value_e_cm=float(limit_value_raw * scale),
        limit_value_raw=float(limit_value_raw),
        limit_units=limit_units,
        limit_operator=limit_operator,
        confidence_level=_optional_str(sub.get("confidence_level")),
        value_summary=_optional_str(sub.get("value_summary")),
        table_value=_optional_float(
            sub.get("table_value"),
            process_id=process_id,
            field_name=f"{anchor.block_key}.table_value",
        ),
        table_units=_optional_str(sub.get("table_units")),
    )


def _load_qcedm_translation_bound(
    process_id: str,
    *,
    nested_key: str,
    coefficient_key: str,
) -> QuarkCEDMTranslationBound:
    anchor, sub = _load_nested_anchor(
        process_id,
        parent_key=_QCEDM_TRANSLATIONS_BLOCK,
        nested_key=nested_key,
        value_key="derived_bound_value",
    )
    if anchor.value <= 0.0:
        raise AnchorError(f"{process_id}: qCEDM derived bound must be positive")
    coefficient = _positive_float(
        sub.get(coefficient_key),
        process_id=process_id,
        field_name=f"{anchor.block_key}.{coefficient_key}",
    )
    operator = _optional_str(sub.get("derived_bound_operator"))
    if operator != _EXPECTED_BOUND_OPERATOR:
        raise AnchorError(
            f"{process_id}: expected derived_bound_operator "
            f"{_EXPECTED_BOUND_OPERATOR!r} for {anchor.block_key}, got {operator!r}"
        )
    units = _optional_str(sub.get("derived_bound_units"))
    if units != _EXPECTED_QCEDM_UNITS:
        raise AnchorError(
            f"{process_id}: expected qCEDM units {_EXPECTED_QCEDM_UNITS!r} "
            f"for {anchor.block_key}, got {units!r}"
        )

    return QuarkCEDMTranslationBound(
        anchor=anchor,
        parent_block_key=_QCEDM_TRANSLATIONS_BLOCK,
        nested_block_key=nested_key,
        coefficient_key=coefficient_key,
        coefficient_value=float(coefficient),
        coefficient_uncertainty_summary=_optional_str(
            sub.get("coefficient_uncertainty_summary")
        ),
        derived_bound_operator=operator,
        derived_bound_observable=_optional_str(sub.get("derived_bound_observable")),
        derived_bound_units=units,
        formula=_optional_str(sub.get("formula")),
        assumptions=_optional_str(sub.get("assumptions")),
        derivation=_optional_str(sub.get("derivation")),
    )


def _assert_rounded_bound_consistency(
    *,
    process_id: str,
    bound: QuarkCEDMTranslationBound,
    unrounded_bound_cm: float,
) -> None:
    if not math.isclose(
        bound.value,
        unrounded_bound_cm,
        rel_tol=_ROUNDING_REL_TOL,
        abs_tol=0.0,
    ):
        raise AnchorError(
            f"{process_id}: {bound.block_key} derived_bound_value={bound.value} "
            f"is inconsistent with the experimental anchor/coefficient "
            f"({unrounded_bound_cm}) beyond rounding tolerance"
        )


def _load_context(process_id: str, key: str) -> E008Context:
    pdg = load_pdg_block(process_id, family=_FAMILY)
    post_context = find_block(
        pdg,
        (_POST_2008_CONTEXT_BLOCK,),
        process_id=process_id,
    )
    sub = find_block(
        post_context,
        (key,),
        process_id=f"{process_id}.{_POST_2008_CONTEXT_BLOCK}",
    )
    values_raw = sub.get("values", ())
    if values_raw is None:
        values: tuple[Mapping[str, Any], ...] = ()
    elif isinstance(values_raw, list):
        if not all(isinstance(item, Mapping) for item in values_raw):
            raise AnchorError(
                f"{process_id}: post_2008_context.{key}.values must contain mappings"
            )
        values = tuple(dict(item) for item in values_raw)
    else:
        raise AnchorError(f"{process_id}: post_2008_context.{key}.values is not a list")
    return E008Context(
        block_key=key,
        source=_optional_str(sub.get("source")),
        year=_optional_int(sub.get("year")),
        source_url=_optional_str(sub.get("source_url")),
        value_summary=_optional_str(sub.get("value_summary")),
        values=values,
        snapshot_path=_optional_str(sub.get("snapshot_path")),
    )


def _load_e008_anchor(process_id: str) -> E008Anchor:
    neutron_edm_limit = _load_experimental_limit(
        process_id,
        nested_key=_NEUTRON_EDM_BLOCK,
        value_key="value",
    )
    mercury_edm_limit = _load_experimental_limit(
        process_id,
        nested_key=_MERCURY_EDM_BLOCK,
        value_key="limit_value",
    )
    neutron_combination = _load_qcedm_translation_bound(
        process_id,
        nested_key=_NEUTRON_QCEDM_BLOCK,
        coefficient_key=_NEUTRON_QCEDM_COEFFICIENT_KEY,
    )
    mercury_isovector_combination = _load_qcedm_translation_bound(
        process_id,
        nested_key=_MERCURY_QCEDM_BLOCK,
        coefficient_key=_MERCURY_QCEDM_COEFFICIENT_KEY,
    )

    neutron_unrounded = (
        neutron_edm_limit.limit_value_e_cm / neutron_combination.coefficient_value
    )
    mercury_unrounded = (
        mercury_edm_limit.limit_value_e_cm
        / mercury_isovector_combination.coefficient_value
    )
    _assert_rounded_bound_consistency(
        process_id=process_id,
        bound=neutron_combination,
        unrounded_bound_cm=neutron_unrounded,
    )
    _assert_rounded_bound_consistency(
        process_id=process_id,
        bound=mercury_isovector_combination,
        unrounded_bound_cm=mercury_unrounded,
    )

    return E008Anchor(
        neutron_edm_limit=neutron_edm_limit,
        mercury_edm_limit=mercury_edm_limit,
        neutron_combination=neutron_combination,
        mercury_isovector_combination=mercury_isovector_combination,
        neutron_unrounded_bound_cm=float(neutron_unrounded),
        mercury_unrounded_bound_cm=float(mercury_unrounded),
        cfw_context=_load_context(process_id, "cfw_baseline"),
        composite_dipoles_context=_load_context(process_id, "composite_dipoles"),
        lattice_qcedm_context=_load_context(process_id, "lattice_qcedm"),
    )


@register
class Constraint:
    """Catalogued non-vetoing quark-cEDM stub."""

    process_id = "E008"
    severity = Severity.INFO
    observable = "quark chromo-EDM reference bounds"

    def __init__(self) -> None:
        self.anchor = _load_e008_anchor(self.process_id)

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        del point
        reference = self.anchor.reference_bound
        reference_bound = quark_cedm_reference_bound(
            bound_observable=reference.derived_bound_observable
            or "quark cEDM combination",
            reference_bound_cm=self.anchor.budget,
        )

        return ConstraintResult(
            process_id=self.process_id,
            severity=self.severity,
            passes=True,
            predicted=None,
            sm_prediction=None,
            experimental=float(reference_bound.reference_bound_cm),
            ratio=None,
            budget=float(reference_bound.reference_bound_cm),
            notes=(
                "INFO-only E008 stub: loads the catalogued neutron/Hg-derived "
                "quark cEDM reference bounds and records the strongest bound. "
                "No RS CP-odd chromo-dipole matching and no hadronic/nuclear "
                "EDM matrix-element calculation are performed; both required "
                "physics inputs are flagged NEEDS-HUMAN-PHYSICS, so this "
                "result is non-vetoing."
            ),
            diagnostics={
                "non_vetoing": True,
                "qcedm_prediction_evaluated": False,
                "hadronic_edm_calculation_evaluated": False,
                "no_hadronic_calculation": True,
                "no_rs_cp_odd_cedm_matching": True,
                "stub_model": QUARK_CEDM_STUB_MODEL_V1,
                "budget_source": _BUDGET_SOURCE,
                "parametrization_citation": _PARAMETRIZATION_CITATION,
                "needs_human_physics": _NEEDS_HUMAN_PHYSICS,
                "needs_human_physics_hadronic_matrix_elements": (
                    QUARK_CEDM_HADRONIC_NEEDS_HUMAN_PHYSICS
                ),
                "needs_human_physics_rs_cp_odd_matching": (
                    QUARK_CEDM_RS_MATCHING_NEEDS_HUMAN_PHYSICS
                ),
                "hadronic_matrix_elements_available": False,
                "rs_cp_odd_quark_cedm_matching_available": False,
                "parameter_point_inputs_used": (),
                "required_low_energy_operators": (
                    "quark_chromo_edm",
                    "quark_edm",
                    "weinberg_three_gluon",
                    "theta_term_policy",
                ),
                "reference_bound_block": reference.block_key,
                "reference_bound_observable": reference_bound.bound_observable,
                "reference_bound_cm": float(reference_bound.reference_bound_cm),
                "reference_bound_operator": reference.derived_bound_operator,
                "reference_bound_units": reference.derived_bound_units,
                "neutron_qcedm_bound_block": self.anchor.neutron_combination.block_key,
                "neutron_qcedm_bound_cm": float(self.anchor.neutron_combination.value),
                "neutron_qcedm_unrounded_bound_cm": float(
                    self.anchor.neutron_unrounded_bound_cm
                ),
                "neutron_qcedm_coefficient": float(
                    self.anchor.neutron_combination.coefficient_value
                ),
                "neutron_qcedm_formula": self.anchor.neutron_combination.formula,
                "neutron_qcedm_assumptions": self.anchor.neutron_combination.assumptions,
                "neutron_edm_limit_e_cm": float(
                    self.anchor.neutron_edm_limit.limit_value_e_cm
                ),
                "neutron_edm_limit_source_url": self.anchor.neutron_edm_limit.source_url,
                "mercury_qcedm_bound_block": (
                    self.anchor.mercury_isovector_combination.block_key
                ),
                "mercury_qcedm_bound_cm": float(
                    self.anchor.mercury_isovector_combination.value
                ),
                "mercury_qcedm_unrounded_bound_cm": float(
                    self.anchor.mercury_unrounded_bound_cm
                ),
                "mercury_qcedm_coefficient": float(
                    self.anchor.mercury_isovector_combination.coefficient_value
                ),
                "mercury_qcedm_formula": (
                    self.anchor.mercury_isovector_combination.formula
                ),
                "mercury_qcedm_assumptions": (
                    self.anchor.mercury_isovector_combination.assumptions
                ),
                "mercury_edm_limit_e_cm": float(
                    self.anchor.mercury_edm_limit.limit_value_e_cm
                ),
                "mercury_edm_limit_source_url": self.anchor.mercury_edm_limit.source_url,
                "cfw_context_summary": self.anchor.cfw_context.value_summary,
                "composite_dipoles_context_summary": (
                    self.anchor.composite_dipoles_context.value_summary
                ),
                "lattice_qcedm_context_summary": (
                    self.anchor.lattice_qcedm_context.value_summary
                ),
            },
        )
