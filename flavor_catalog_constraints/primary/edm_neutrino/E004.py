"""E004 - neutron electric dipole moment ``|d_n|``.

Physics
-------
The catalogued observable is the neutron EDM upper limit

    |d_n| < 1.8e-26 e cm  (90% CL),

anchored to the PDG/nEDM@PSI Abel et al. measurement in
``flavor_catalog/processes/edm_neutrino/E004.yaml``.  A rigorous RS
prediction is **not** available in this scaffold.  It would require

1. RS matching onto CP-odd quark EDMs, quark chromo-EDMs, and gluonic
   operators; and
2. non-perturbative neutron matrix elements translating those operators into
   ``d_n``.

Both ingredients are flagged ``NEEDS-HUMAN-PHYSICS`` in the returned
diagnostics.  This file does not fake a hadronic calculation and does not
read quark-sector inputs from ``ParameterPoint``.

Severity
--------
INFO.  The neutron EDM is a stringent observed bound, but applying it as a
hard veto here would require non-perturbative hadronic input and RS CP-odd
quark-dipole matching that are absent from the current catalog contract.  The
returned ``passes`` value is bookkeeping only and must not veto scan points.

Catalog sidecar
---------------
``flavor_catalog/processes/edm_neutrino/E004.yaml`` is the source of truth for
the PDG Live 2026 canonical limit and the Abel et al. PSI primary measurement.
Numeric values below are loaded through the scaffold anchor loader and then
unit-converted where the sidecar stores values in ``10^-26 e cm``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Mapping

from flavor_catalog_constraints.anchors import (
    Anchor,
    AnchorError,
    load_anchor,
    load_full_yaml,
    load_pdg_block,
)
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.neutron_edm import (
    NEUTRON_EDM_HADRONIC_NEEDS_HUMAN_PHYSICS,
    NEUTRON_EDM_RS_MATCHING_NEEDS_HUMAN_PHYSICS,
    NEUTRON_EDM_STUB_MODEL_V1,
    compare_neutron_edm_measurement_to_limit,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "edm_neutrino"
_CURRENT_LIMIT_BLOCK = "canonical_limit"
_PRIMARY_MEASUREMENT_BLOCK = "primary_measurement"
_LIMIT_ANCHOR_CANDIDATES = (_CURRENT_LIMIT_BLOCK,)
_EXPECTED_LIMIT_UNITS = "e cm"
_EXPECTED_PRIMARY_MEASUREMENT_UNITS = "10^-26 e cm"
_EXPECTED_LIMIT_OPERATOR = "<"
_BUDGET_SOURCE = (
    "flavor_catalog/processes/edm_neutrino/E004.yaml canonical_limit "
    "(PDG Live 2026, Abel et al. 2020 nEDM@PSI)"
)
_PARAMETRIZATION_CITATION = (
    "Abel et al. Phys. Rev. Lett. 124, 081803 (2020), arXiv:2001.11966; "
    "PDG Live 2026 S017EDM neutron electric dipole moment datablock"
)
_NEEDS_HUMAN_PHYSICS = (
    NEUTRON_EDM_HADRONIC_NEEDS_HUMAN_PHYSICS,
    NEUTRON_EDM_RS_MATCHING_NEEDS_HUMAN_PHYSICS,
)


@dataclass(frozen=True)
class NeutronEDMLimitAnchor:
    """Typed current upper limit on ``|d_n|`` in e cm."""

    experimental: Anchor
    limit_operator: str | None
    confidence_level: str | None
    used_measurement: str | None
    table_value: float | None
    table_units: str | None
    value_summary: str | None

    @property
    def value(self) -> float:
        """Upper limit on ``|d_n|`` in e cm."""
        return self.experimental.value

    @property
    def budget(self) -> float:
        """Non-vetoing bookkeeping budget: the measured upper limit."""
        return self.experimental.value

    @property
    def source_url(self) -> str | None:
        """Primary source URL for the canonical limit."""
        return self.experimental.source_url


@dataclass(frozen=True)
class NeutronEDMPrimaryMeasurement:
    """Typed Abel et al. PSI measurement, converted to e cm."""

    block_key: str
    source: str | None
    year: int | None
    central_value_e_cm: float
    statistical_uncertainty_e_cm: float
    systematic_uncertainty_e_cm: float
    total_uncertainty_e_cm: float
    limit_value_e_cm: float
    limit_operator: str | None
    confidence_level: str | None
    units: str | None
    limit_units: str | None
    source_url: str | None
    snapshot_path: str | None
    measurement_summary: str | None


@dataclass(frozen=True)
class E004Context:
    """Post-2008 context carried for provenance, not for calculation."""

    key: str | None
    source_url: str | None
    year: int | None
    value_summary: str | None
    snapshot_path: str | None


@dataclass(frozen=True)
class E004Anchor:
    """Typed E004 anchor bundle and non-vetoing budget convention."""

    limit: NeutronEDMLimitAnchor
    primary_measurement: NeutronEDMPrimaryMeasurement
    composite_dipoles_context: E004Context
    lattice_hadronic_context: E004Context

    @property
    def value(self) -> float:
        """Canonical current upper limit on ``|d_n|`` in e cm."""
        return self.limit.value

    @property
    def budget(self) -> float:
        """Non-vetoing bookkeeping budget: the canonical current limit."""
        return self.limit.budget

    @property
    def source_url(self) -> str | None:
        """Source URL for the canonical current limit."""
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


def _optional_float(value: Any, *, process_id: str, field_name: str) -> float | None:
    if value is None:
        return None
    return _required_float(value, process_id=process_id, field_name=field_name)


def _required_float(value: Any, *, process_id: str, field_name: str) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError) as exc:
        raise AnchorError(
            f"{process_id}: E004 anchor field {field_name!r}={value!r} "
            "is not numeric"
        ) from exc
    if not math.isfinite(out):
        raise AnchorError(f"{process_id}: E004 anchor field {field_name!r} is not finite")
    return out


def _positive_float(value: Any, *, process_id: str, field_name: str) -> float:
    out = _required_float(value, process_id=process_id, field_name=field_name)
    if out <= 0.0:
        raise AnchorError(f"{process_id}: E004 anchor field {field_name!r} must be positive")
    return out


def _scale_to_e_cm(units: str | None, *, process_id: str, field_name: str) -> float:
    if units == "e cm":
        return 1.0
    if units == _EXPECTED_PRIMARY_MEASUREMENT_UNITS:
        return 1.0e-26
    raise AnchorError(
        f"{process_id}: unsupported E004 units for {field_name}: {units!r}"
    )


def _load_limit_anchor(process_id: str) -> NeutronEDMLimitAnchor:
    experimental = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=_LIMIT_ANCHOR_CANDIDATES,
    )
    if experimental.block_key != _CURRENT_LIMIT_BLOCK:
        raise AnchorError(
            f"{process_id}: load_anchor selected {experimental.block_key!r}, "
            f"expected {_CURRENT_LIMIT_BLOCK!r} for the current neutron EDM limit"
        )
    if experimental.value <= 0.0:
        raise AnchorError(f"{process_id}: neutron EDM limit must be positive")
    if experimental.units != _EXPECTED_LIMIT_UNITS:
        raise AnchorError(
            f"{process_id}: expected units {_EXPECTED_LIMIT_UNITS!r}, "
            f"got {experimental.units!r}"
        )

    pdg = load_pdg_block(process_id, family=_FAMILY)
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
    return NeutronEDMLimitAnchor(
        experimental=experimental,
        limit_operator=limit_operator,
        confidence_level=_optional_str(sub.get("confidence_level")),
        used_measurement=_optional_str(sub.get("used_measurement")),
        table_value=_optional_float(
            sub.get("table_value"),
            process_id=process_id,
            field_name=f"{experimental.block_key}.table_value",
        ),
        table_units=_optional_str(sub.get("table_units")),
        value_summary=_optional_str(sub.get("value_summary")),
    )


def _load_primary_measurement(
    process_id: str,
    *,
    canonical_limit_e_cm: float,
) -> NeutronEDMPrimaryMeasurement:
    pdg = load_pdg_block(process_id, family=_FAMILY)
    sub = pdg.get(_PRIMARY_MEASUREMENT_BLOCK)
    if not isinstance(sub, Mapping):
        raise AnchorError(
            f"{process_id}: {_PRIMARY_MEASUREMENT_BLOCK!r} is not a mapping"
        )

    units = _optional_str(sub.get("units"))
    limit_units = _optional_str(sub.get("limit_units"))
    measurement_scale = _scale_to_e_cm(
        units,
        process_id=process_id,
        field_name=f"{_PRIMARY_MEASUREMENT_BLOCK}.units",
    )
    limit_scale = _scale_to_e_cm(
        limit_units,
        process_id=process_id,
        field_name=f"{_PRIMARY_MEASUREMENT_BLOCK}.limit_units",
    )
    stat = _positive_float(
        sub.get("statistical_uncertainty"),
        process_id=process_id,
        field_name=f"{_PRIMARY_MEASUREMENT_BLOCK}.statistical_uncertainty",
    )
    syst = _positive_float(
        sub.get("systematic_uncertainty"),
        process_id=process_id,
        field_name=f"{_PRIMARY_MEASUREMENT_BLOCK}.systematic_uncertainty",
    )
    central_value_e_cm = _required_float(
        sub.get("value"),
        process_id=process_id,
        field_name=f"{_PRIMARY_MEASUREMENT_BLOCK}.value",
    ) * measurement_scale
    limit_value_e_cm = _positive_float(
        sub.get("limit_value"),
        process_id=process_id,
        field_name=f"{_PRIMARY_MEASUREMENT_BLOCK}.limit_value",
    ) * limit_scale
    if not math.isclose(
        limit_value_e_cm,
        canonical_limit_e_cm,
        rel_tol=1.0e-12,
        abs_tol=0.0,
    ):
        raise AnchorError(
            f"{process_id}: primary measurement limit {limit_value_e_cm} e cm "
            f"does not match canonical limit {canonical_limit_e_cm} e cm"
        )

    return NeutronEDMPrimaryMeasurement(
        block_key=_PRIMARY_MEASUREMENT_BLOCK,
        source=_optional_str(sub.get("source")),
        year=_optional_int(sub.get("year")),
        central_value_e_cm=float(central_value_e_cm),
        statistical_uncertainty_e_cm=float(stat * measurement_scale),
        systematic_uncertainty_e_cm=float(syst * measurement_scale),
        total_uncertainty_e_cm=float(math.sqrt(stat * stat + syst * syst) * measurement_scale),
        limit_value_e_cm=float(limit_value_e_cm),
        limit_operator=_optional_str(sub.get("limit_operator")),
        confidence_level=_optional_str(sub.get("confidence_level")),
        units=units,
        limit_units=limit_units,
        source_url=_optional_str(sub.get("source_url")),
        snapshot_path=_optional_str(sub.get("snapshot_path")),
        measurement_summary=_optional_str(sub.get("measurement_summary")),
    )


def _load_context(process_id: str, key: str) -> E004Context:
    data = load_full_yaml(process_id, family=_FAMILY)
    post_context = data.get("post_2008_context")
    if not isinstance(post_context, Mapping):
        raise AnchorError(f"{process_id}: missing post_2008_context mapping")
    sub = post_context.get(key)
    if not isinstance(sub, Mapping):
        raise AnchorError(f"{process_id}: post_2008_context.{key!r} is not a mapping")
    return E004Context(
        key=_optional_str(sub.get("key")),
        source_url=_optional_str(sub.get("source_url")),
        year=_optional_int(sub.get("year")),
        value_summary=_optional_str(sub.get("value_summary")),
        snapshot_path=_optional_str(sub.get("snapshot_path")),
    )


def _load_e004_anchor(process_id: str) -> E004Anchor:
    limit = _load_limit_anchor(process_id)
    primary_measurement = _load_primary_measurement(
        process_id,
        canonical_limit_e_cm=limit.value,
    )
    if primary_measurement.limit_operator != limit.limit_operator:
        raise AnchorError(
            f"{process_id}: primary measurement limit operator "
            f"{primary_measurement.limit_operator!r} does not match canonical "
            f"{limit.limit_operator!r}"
        )
    if primary_measurement.confidence_level != limit.confidence_level:
        raise AnchorError(
            f"{process_id}: primary measurement confidence level "
            f"{primary_measurement.confidence_level!r} does not match canonical "
            f"{limit.confidence_level!r}"
        )
    return E004Anchor(
        limit=limit,
        primary_measurement=primary_measurement,
        composite_dipoles_context=_load_context(process_id, "composite_dipoles"),
        lattice_hadronic_context=_load_context(process_id, "lattice_hadronic_translation"),
    )


@register
class Constraint:
    """Catalogued non-vetoing neutron EDM stub."""

    process_id = "E004"
    severity = Severity.INFO
    observable = "|d_n|"

    def __init__(self) -> None:
        self.anchor = _load_e004_anchor(self.process_id)

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        del point
        comparison = compare_neutron_edm_measurement_to_limit(
            measured_neutron_edm_e_cm=(
                self.anchor.primary_measurement.central_value_e_cm
            ),
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
                "INFO-only E004 stub: loads the neutron EDM nEDM@PSI/PDG "
                "limit and records the observable. No quark EDM/CEDM/"
                "Weinberg matching or neutron hadronic matrix-element "
                "calculation is performed; both required physics inputs are "
                "flagged NEEDS-HUMAN-PHYSICS, so this result is non-vetoing."
            ),
            diagnostics={
                "non_vetoing": True,
                "edm_prediction_evaluated": False,
                "no_hadronic_calculation": True,
                "no_quark_dipole_matching": True,
                "stub_model": NEUTRON_EDM_STUB_MODEL_V1,
                "budget_source": _BUDGET_SOURCE,
                "budget_is_current_neutron_edm_limit": True,
                "parametrization_citation": _PARAMETRIZATION_CITATION,
                "needs_human_physics": _NEEDS_HUMAN_PHYSICS,
                "needs_human_physics_hadronic_matrix_elements": (
                    NEUTRON_EDM_HADRONIC_NEEDS_HUMAN_PHYSICS
                ),
                "needs_human_physics_rs_cp_odd_matching": (
                    NEUTRON_EDM_RS_MATCHING_NEEDS_HUMAN_PHYSICS
                ),
                "hadronic_matrix_elements_available": False,
                "rs_cp_odd_quark_dipole_matching_available": False,
                "parameter_point_inputs_used": (),
                "required_low_energy_operators": (
                    "quark_edm",
                    "quark_chromo_edm",
                    "weinberg_three_gluon",
                ),
                "experimental_block": self.anchor.limit.experimental.block_key,
                "primary_measurement_block": (
                    self.anchor.primary_measurement.block_key
                ),
                "experimental_limit_e_cm": float(self.anchor.value),
                "measurement_central_e_cm": float(
                    self.anchor.primary_measurement.central_value_e_cm
                ),
                "measurement_abs_e_cm": float(comparison.measurement_abs_e_cm),
                "measurement_to_limit_ratio": float(comparison.ratio_to_limit),
                "statistical_uncertainty_e_cm": float(
                    self.anchor.primary_measurement.statistical_uncertainty_e_cm
                ),
                "systematic_uncertainty_e_cm": float(
                    self.anchor.primary_measurement.systematic_uncertainty_e_cm
                ),
                "total_uncertainty_e_cm": float(
                    self.anchor.primary_measurement.total_uncertainty_e_cm
                ),
                "limit_operator": self.anchor.limit.limit_operator,
                "confidence_level": self.anchor.limit.confidence_level,
                "used_measurement": self.anchor.limit.used_measurement,
                "value_summary": self.anchor.limit.value_summary,
                "table_value": self.anchor.limit.table_value,
                "table_units": self.anchor.limit.table_units,
                "primary_measurement_source": (
                    self.anchor.primary_measurement.source
                ),
                "primary_measurement_source_url": (
                    self.anchor.primary_measurement.source_url
                ),
                "primary_measurement_summary": (
                    self.anchor.primary_measurement.measurement_summary
                ),
                "composite_dipoles_context_key": (
                    self.anchor.composite_dipoles_context.key
                ),
                "composite_dipoles_context_summary": (
                    self.anchor.composite_dipoles_context.value_summary
                ),
                "lattice_hadronic_context_key": (
                    self.anchor.lattice_hadronic_context.key
                ),
                "lattice_hadronic_context_summary": (
                    self.anchor.lattice_hadronic_context.value_summary
                ),
            },
        )
