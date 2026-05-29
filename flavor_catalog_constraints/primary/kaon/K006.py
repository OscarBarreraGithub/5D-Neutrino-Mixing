"""K006 - short-distance component of ``K_L -> mu+ mu-``.

Physics
-------
The measured ``K_L -> mu+ mu-`` branching fraction is dominated by the
long-distance two-photon amplitude.  The constrained quantity here is only the
short-distance dispersive component,

    BR(K_L -> mu+ mu-)_SD = kappa_mu [
        Re(Y_eff) / lambda^5 + Re(lambda_c) P_c(Y) / lambda
    ]^2,

with ``Y_eff = lambda_t (Y_L - Y_R)``.  The low-level formula and the
documented RS matching assumption live in
``quarkConstraints.rare_kaon_dilepton`` and are reached only through the
``flavor_catalog_constraints.physics_adapters.rare_kaon_dilepton`` boundary.
This module intentionally does not touch ``quarkConstraints.rare_kaon_snd``,
which is reserved for ``K -> pi nu nubar``.

Severity
--------
HARD for the short-distance-bound verdict.  The active budget is the
Isidori-Unterdorfer conservative upper bound on the short-distance component
loaded from ``K006.yaml``.  The total PDG rate is long-distance dominated and
is kept diagnostic-only rather than being used as the veto observable.

Catalog sidecar
---------------
``flavor_catalog/processes/kaon/K006.yaml`` is the source of truth for the PDG
total branching fraction, the SM short-distance theory anchor, and the
conservative short-distance bound.  K006 stores the total value as a flat
``pdg_or_equivalent`` mapping and the SD bound in a supporting list, so this
module adapts both into the scaffold ``load_anchor`` path rather than reading
value-bearing anchors ad hoc.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
import re
from typing import Any, Mapping, Sequence

from flavor_catalog_constraints import anchors as anchor_scaffold
from flavor_catalog_constraints.anchors import (
    Anchor,
    AnchorError,
    load_anchor,
    load_full_yaml,
)
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.rare_kaon_dilepton import (
    RARE_KAON_DILEPTON_PARAMETRIZATION_CITATION,
    RARE_KAON_DILEPTON_RS_MATCHING_ASSUMPTION_V1,
    klong_mumu_short_distance_from_couplings,
    klong_mumu_short_distance_sm,
    rare_kaon_dilepton_default_sm_inputs,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "kaon"
_REQUIRED_EXTRA = "quark_mass_basis_couplings"
_OPTIONAL_EW_MASS_EXTRA = "kk_ew_mass_gev"
_PDG_TOTAL_BLOCK = "pdg_or_equivalent"
_E871_BLOCK = "supporting_numeric_values[0]"
_SM_SHORT_DISTANCE_BLOCK = "standard_model_short_distance_prediction"
_SHORT_DISTANCE_BOUND_BLOCK = "supporting_numeric_values[2]"
_EXPECTED_UNITS = ("branching fraction", "dimensionless branching fraction")
_SM_ANCHOR_SIGMA_TOLERANCE = 2.0
_SHORT_DISTANCE_BOUND_OBSERVABLE = (
    "Conservative upper bound on B(K_L -> mu+ mu-)_short"
)
_BUDGET_SOURCE = (
    "flavor_catalog/processes/kaon/K006.yaml "
    f"{_SHORT_DISTANCE_BOUND_BLOCK} Isidori-Unterdorfer conservative SD bound"
)
_UPPER_LIMIT_RE = re.compile(r"^\s*<\s*(?P<value>[0-9.eE+-]+)\s*$")
_LIST_BLOCK_RE = re.compile(
    r"^(?P<section>[A-Za-z_][A-Za-z0-9_]*)\[(?P<index>[0-9]+)\]$"
)


@dataclass(frozen=True)
class BranchingFractionAnchor:
    """Branching-fraction value or upper-limit entry from K006.yaml."""

    section_key: str
    entry_index: int
    block_key: str
    source: str | None
    year: int | None
    value: float
    uncertainty: float | None
    is_upper_limit: bool
    confidence_level: str | None
    units: str | None
    source_url: str | None
    snapshot_path: str | None
    value_raw: str | None
    observable: str | None


@dataclass(frozen=True)
class KLongMuMuShortDistanceBudget:
    """K006 short-distance upper-bound budget."""

    source: str
    construction: str
    hard_veto_budget: float
    confidence_level: str | None
    total_rate_value: float
    total_rate_uncertainty: float | None
    sm_anchor_value: float
    sm_anchor_uncertainty: float
    sm_formula_value: float
    sm_formula_minus_anchor: float
    sm_formula_anchor_pull: float
    bound_minus_sm_formula: float
    sm_subtracted: bool
    long_distance_dominated: bool


@dataclass(frozen=True)
class K006Anchor:
    """Typed K006 anchor: total, E871, SM-SD, SD bound, and budget."""

    total_rate: BranchingFractionAnchor
    e871_measurement: BranchingFractionAnchor
    standard_model_short_distance: Anchor
    short_distance_bound: BranchingFractionAnchor
    budget_band: KLongMuMuShortDistanceBudget

    @property
    def value(self) -> float:
        return self.short_distance_bound.value

    @property
    def uncertainty(self) -> float | None:
        return None

    @property
    def total_value(self) -> float:
        return self.total_rate.value

    @property
    def sm_value(self) -> float:
        return self.standard_model_short_distance.value

    @property
    def budget(self) -> float:
        return self.budget_band.hard_veto_budget


def _optional_str(value: Any) -> str | None:
    return None if value is None else str(value)


def _required_float(value: Any, *, process_id: str, field_name: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise AnchorError(
            f"{process_id}: K006 anchor field {field_name!r}={value!r} "
            "is not numeric"
        ) from exc
    if not math.isfinite(number):
        raise AnchorError(
            f"{process_id}: K006 anchor field {field_name!r}={value!r} "
            "is not finite"
        )
    return number


def _parse_value(
    value: Any,
    *,
    process_id: str,
    field_name: str,
    expect_upper_limit: bool,
) -> tuple[float, bool, str | None]:
    raw = _optional_str(value)
    if isinstance(value, str):
        match = _UPPER_LIMIT_RE.match(value)
        if match is not None:
            number = _required_float(
                match.group("value"),
                process_id=process_id,
                field_name=field_name,
            )
            return number, True, raw
    if expect_upper_limit:
        raise AnchorError(
            f"{process_id}: expected an upper-limit value in {field_name}, "
            f"got {value!r}"
        )
    return (
        _required_float(value, process_id=process_id, field_name=field_name),
        False,
        raw,
    )


def _validate_units(units: str | None, *, process_id: str, block_key: str) -> None:
    if units not in _EXPECTED_UNITS:
        raise AnchorError(
            f"{process_id}: expected branching-fraction units for {block_key}, "
            f"got {units!r}"
        )


def _parse_list_block_key(block_key: str, *, process_id: str) -> tuple[str, int]:
    match = _LIST_BLOCK_RE.match(block_key)
    if match is None:
        raise AnchorError(f"{process_id}: invalid K006 list anchor key {block_key!r}")
    return match.group("section"), int(match.group("index"))


def _section_entries(
    data: Mapping[str, Any],
    section_key: str,
    *,
    process_id: str,
) -> Sequence[Mapping[str, Any]]:
    entries = data.get(section_key)
    if not isinstance(entries, list) or not entries:
        raise AnchorError(
            f"{process_id}: expected non-empty list section {section_key!r}"
        )
    for index, entry in enumerate(entries):
        if not isinstance(entry, Mapping):
            raise AnchorError(
                f"{process_id}: {section_key}[{index}] is not a mapping"
            )
    return entries


def _list_entry_for_block(
    data: Mapping[str, Any],
    block_key: str,
    *,
    process_id: str,
) -> tuple[str, int, Mapping[str, Any]]:
    section_key, entry_index = _parse_list_block_key(
        block_key,
        process_id=process_id,
    )
    entries = _section_entries(data, section_key, process_id=process_id)
    try:
        entry = entries[entry_index]
    except IndexError as exc:
        raise AnchorError(
            f"{process_id}: K006 list anchor {block_key!r} is out of range"
        ) from exc
    return section_key, entry_index, entry


def _load_virtual_anchor(
    process_id: str,
    *,
    block_key: str,
    entry: Mapping[str, Any],
    expect_upper_limit: bool,
) -> tuple[Anchor, bool, str | None]:
    value, is_upper_limit, raw = _parse_value(
        entry.get("value"),
        process_id=process_id,
        field_name=f"{block_key}.value",
        expect_upper_limit=expect_upper_limit,
    )
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
        )
    finally:
        anchor_scaffold.load_pdg_block = original_load_pdg_block

    return scaffold_anchor, is_upper_limit, raw


def _branching_anchor_from_scaffold(
    scaffold_anchor: Anchor,
    *,
    entry: Mapping[str, Any],
    section_key: str,
    entry_index: int,
    is_upper_limit: bool,
    value_raw: str | None,
    process_id: str,
) -> BranchingFractionAnchor:
    _validate_units(
        scaffold_anchor.units,
        process_id=process_id,
        block_key=scaffold_anchor.block_key,
    )
    return BranchingFractionAnchor(
        section_key=section_key,
        entry_index=entry_index,
        block_key=scaffold_anchor.block_key,
        source=_optional_str(scaffold_anchor.source),
        year=scaffold_anchor.year,
        value=float(scaffold_anchor.value),
        uncertainty=(
            None
            if scaffold_anchor.uncertainty is None
            else float(scaffold_anchor.uncertainty)
        ),
        is_upper_limit=is_upper_limit,
        confidence_level=_optional_str(entry.get("confidence_level")),
        units=scaffold_anchor.units,
        source_url=_optional_str(scaffold_anchor.source_url),
        snapshot_path=_optional_str(scaffold_anchor.snapshot_path),
        value_raw=value_raw,
        observable=_optional_str(scaffold_anchor.observable),
    )


def _load_flat_pdg_anchor(process_id: str) -> BranchingFractionAnchor:
    data = load_full_yaml(process_id, family=_FAMILY)
    entry = data.get(_PDG_TOTAL_BLOCK)
    if not isinstance(entry, Mapping):
        raise AnchorError(
            f"{process_id}: expected flat mapping {_PDG_TOTAL_BLOCK!r}"
        )
    scaffold_anchor, is_upper_limit, raw = _load_virtual_anchor(
        process_id,
        block_key=_PDG_TOTAL_BLOCK,
        entry=entry,
        expect_upper_limit=False,
    )
    return _branching_anchor_from_scaffold(
        scaffold_anchor,
        entry=entry,
        section_key=_PDG_TOTAL_BLOCK,
        entry_index=-1,
        is_upper_limit=is_upper_limit,
        value_raw=raw,
        process_id=process_id,
    )


def _load_supporting_anchor(
    process_id: str,
    *,
    block_key: str,
    expect_upper_limit: bool,
    observable: str | None = None,
) -> BranchingFractionAnchor:
    data = load_full_yaml(process_id, family=_FAMILY)
    section_key, entry_index, entry = _list_entry_for_block(
        data,
        block_key,
        process_id=process_id,
    )
    if observable is not None and entry.get("observable") != observable:
        raise AnchorError(
            f"{process_id}: expected {block_key}.observable={observable!r}, "
            f"got {entry.get('observable')!r}"
        )
    scaffold_anchor, is_upper_limit, raw = _load_virtual_anchor(
        process_id,
        block_key=block_key,
        entry=entry,
        expect_upper_limit=expect_upper_limit,
    )
    return _branching_anchor_from_scaffold(
        scaffold_anchor,
        entry=entry,
        section_key=section_key,
        entry_index=entry_index,
        is_upper_limit=is_upper_limit,
        value_raw=raw,
        process_id=process_id,
    )


def _load_sm_short_distance_anchor(process_id: str) -> Anchor:
    anchor = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=(_SM_SHORT_DISTANCE_BLOCK,),
    )
    _validate_units(
        anchor.units,
        process_id=process_id,
        block_key=anchor.block_key,
    )
    if anchor.value <= 0.0 or not math.isfinite(anchor.value):
        raise AnchorError(
            f"{process_id}: SM short-distance anchor must be positive and finite"
        )
    if (
        anchor.uncertainty is None
        or anchor.uncertainty <= 0.0
        or not math.isfinite(anchor.uncertainty)
    ):
        raise AnchorError(
            f"{process_id}: SM short-distance anchor requires a positive uncertainty"
        )
    return anchor


def _validate_sm_formula_against_anchor(
    *,
    process_id: str,
    standard_model_short_distance: Anchor,
    sm_formula_value: float,
) -> None:
    delta = abs(float(sm_formula_value) - float(standard_model_short_distance.value))
    tolerance = (
        _SM_ANCHOR_SIGMA_TOLERANCE
        * float(standard_model_short_distance.uncertainty)
    )
    if delta > tolerance:
        raise AnchorError(
            f"{process_id}: SM short-distance formula {sm_formula_value:.6e} "
            f"differs from YAML anchor "
            f"{standard_model_short_distance.value:.6e} by more than "
            f"{_SM_ANCHOR_SIGMA_TOLERANCE:g} sigma"
        )


def _build_budget_band(
    *,
    total_rate: BranchingFractionAnchor,
    standard_model_short_distance: Anchor,
    short_distance_bound: BranchingFractionAnchor,
    sm_formula_value: float,
) -> KLongMuMuShortDistanceBudget:
    budget = float(short_distance_bound.value)
    if not short_distance_bound.is_upper_limit or budget <= 0.0:
        raise AnchorError("K006: short-distance budget must be a positive upper limit")
    sm_uncertainty = float(standard_model_short_distance.uncertainty)
    sm_delta = float(sm_formula_value - standard_model_short_distance.value)
    return KLongMuMuShortDistanceBudget(
        source=_BUDGET_SOURCE,
        construction=(
            "Compare predicted total short-distance branching fraction "
            "against the conservative Isidori-Unterdorfer bound; do not "
            "use the long-distance-dominated total PDG rate as a clean veto."
        ),
        hard_veto_budget=budget,
        confidence_level=short_distance_bound.confidence_level,
        total_rate_value=float(total_rate.value),
        total_rate_uncertainty=total_rate.uncertainty,
        sm_anchor_value=float(standard_model_short_distance.value),
        sm_anchor_uncertainty=sm_uncertainty,
        sm_formula_value=float(sm_formula_value),
        sm_formula_minus_anchor=sm_delta,
        sm_formula_anchor_pull=float(sm_delta / sm_uncertainty),
        bound_minus_sm_formula=float(budget - sm_formula_value),
        sm_subtracted=False,
        long_distance_dominated=True,
    )


def _load_k006_anchor(process_id: str, *, sm_formula_value: float) -> K006Anchor:
    total_rate = _load_flat_pdg_anchor(process_id)
    e871_measurement = _load_supporting_anchor(
        process_id,
        block_key=_E871_BLOCK,
        expect_upper_limit=False,
    )
    standard_model_short_distance = _load_sm_short_distance_anchor(process_id)
    _validate_sm_formula_against_anchor(
        process_id=process_id,
        standard_model_short_distance=standard_model_short_distance,
        sm_formula_value=sm_formula_value,
    )
    short_distance_bound = _load_supporting_anchor(
        process_id,
        block_key=_SHORT_DISTANCE_BOUND_BLOCK,
        expect_upper_limit=True,
        observable=_SHORT_DISTANCE_BOUND_OBSERVABLE,
    )
    return K006Anchor(
        total_rate=total_rate,
        e871_measurement=e871_measurement,
        standard_model_short_distance=standard_model_short_distance,
        short_distance_bound=short_distance_bound,
        budget_band=_build_budget_band(
            total_rate=total_rate,
            standard_model_short_distance=standard_model_short_distance,
            short_distance_bound=short_distance_bound,
            sm_formula_value=sm_formula_value,
        ),
    )


@register
class Constraint:
    """Catalogued short-distance ``K_L -> mu+ mu-`` constraint (K006)."""

    process_id = "K006"
    severity = Severity.HARD
    observable = "BR(K_L -> mu+ mu-)_SD"

    def __init__(self) -> None:
        self.sm_inputs = rare_kaon_dilepton_default_sm_inputs()
        self.sm_result = klong_mumu_short_distance_sm(self.sm_inputs)
        self.anchor = _load_k006_anchor(
            self.process_id,
            sm_formula_value=float(self.sm_result.branching_fraction),
        )

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        couplings = point.get_extra(_REQUIRED_EXTRA)
        if couplings is None:
            return ConstraintResult(
                process_id=self.process_id,
                severity=self.severity,
                passes=True,
                sm_prediction=float(self.sm_result.branching_fraction),
                experimental=float(self.anchor.value),
                budget=float(self.anchor.budget),
                notes=(
                    f"extra {_REQUIRED_EXTRA!r} absent; K_L -> mu+ mu- "
                    "short-distance constraint was not evaluated."
                ),
                diagnostics={
                    "missing_extra": _REQUIRED_EXTRA,
                    "pdg_total_branching_fraction": float(
                        self.anchor.total_rate.value
                    ),
                    "pdg_total_uncertainty": self.anchor.total_rate.uncertainty,
                    "short_distance_bound": float(
                        self.anchor.short_distance_bound.value
                    ),
                    "sm_anchor_branching_fraction": float(
                        self.anchor.standard_model_short_distance.value
                    ),
                    "sm_anchor_uncertainty": (
                        self.anchor.standard_model_short_distance.uncertainty
                    ),
                    "sm_formula_minus_anchor": float(
                        self.anchor.budget_band.sm_formula_minus_anchor
                    ),
                    "sm_formula_anchor_pull": float(
                        self.anchor.budget_band.sm_formula_anchor_pull
                    ),
                    "budget_source": self.anchor.budget_band.source,
                    "long_distance_dominated_total_rate": True,
                    "severity_justification": (
                        "HARD: the predicted SD component is compared directly "
                        "with the conservative long-distance-separated "
                        "dispersive upper bound."
                    ),
                },
            )

        kk_ew_mass = point.get_extra(_OPTIONAL_EW_MASS_EXTRA)
        result = klong_mumu_short_distance_from_couplings(
            couplings,
            m_kk_gev=None if kk_ew_mass is None else float(kk_ew_mass),
            inputs=self.sm_inputs,
        )
        predicted = float(result.branching_fraction)
        budget = float(self.anchor.budget)
        ratio = predicted / budget if budget > 0.0 else float("inf")
        passes = bool(ratio <= 1.0)

        diagnostics = dict(result.diagnostics)
        diagnostics.update(
            {
                "pdg_total_branching_fraction": float(
                    self.anchor.total_rate.value
                ),
                "pdg_total_uncertainty": self.anchor.total_rate.uncertainty,
                "pdg_total_block": self.anchor.total_rate.block_key,
                "e871_branching_fraction": float(
                    self.anchor.e871_measurement.value
                ),
                "short_distance_bound": float(
                    self.anchor.short_distance_bound.value
                ),
                "short_distance_bound_block": (
                    self.anchor.short_distance_bound.block_key
                ),
                "sm_anchor_branching_fraction": float(
                    self.anchor.standard_model_short_distance.value
                ),
                "sm_anchor_uncertainty": (
                    self.anchor.standard_model_short_distance.uncertainty
                ),
                "sm_anchor_block": (
                    self.anchor.standard_model_short_distance.block_key
                ),
                "sm_formula_minus_anchor": float(
                    self.anchor.budget_band.sm_formula_minus_anchor
                ),
                "sm_formula_anchor_pull": float(
                    self.anchor.budget_band.sm_formula_anchor_pull
                ),
                "budget_bound_minus_sm_formula": float(
                    self.anchor.budget_band.bound_minus_sm_formula
                ),
                "budget_sm_subtracted": self.anchor.budget_band.sm_subtracted,
                "budget_long_distance_dominated": (
                    self.anchor.budget_band.long_distance_dominated
                ),
                "budget_source": self.anchor.budget_band.source,
                "parametrization_citation": (
                    RARE_KAON_DILEPTON_PARAMETRIZATION_CITATION
                ),
                "rs_matching_assumption": (
                    RARE_KAON_DILEPTON_RS_MATCHING_ASSUMPTION_V1
                ),
                "needs_human_physics": (
                    "NEEDS-HUMAN-PHYSICS: full RS electroweak KK/Z/Z' tower "
                    "and muon axial-coupling matching are not available on "
                    "ParameterPoint; v1 uses the documented Z-like proxy."
                ),
                "severity_justification": (
                    "HARD: K_L -> mu+ mu- is long-distance dominated, but "
                    "the tested observable is the extracted conservative SD "
                    "upper bound."
                ),
                "kk_ew_mass_extra_used": kk_ew_mass is not None,
                "kappa_mu": float(result.kappa_mu),
                "p_c_y": float(result.p_c_y),
                "y_t": float(result.y_t),
                "lambda_wolfenstein": float(result.lambda_wolfenstein),
                "lambda_c": complex(result.lambda_c),
                "lambda_t": complex(result.lambda_t),
                "y_eff": complex(result.y_eff),
                "y_np_total": complex(result.y_np_total),
                "short_distance_amplitude": float(
                    result.short_distance_amplitude
                ),
                "np_shift_branching_fraction": float(
                    result.np_shift_branching_fraction
                ),
            }
        )

        return ConstraintResult(
            process_id=self.process_id,
            severity=self.severity,
            passes=passes,
            predicted=predicted,
            sm_prediction=float(result.sm_branching_fraction),
            experimental=float(self.anchor.value),
            ratio=float(ratio),
            budget=budget,
            notes=(
                "BR(K_L -> mu+ mu-)_SD uses the Buras/Isidori "
                "short-distance parametrization. RS contribution is a "
                "documented Z-like Y-function shift from mass-basis s-d "
                "couplings; full EW/lepton matching is marked "
                "NEEDS-HUMAN-PHYSICS. The HARD budget is the conservative "
                "short-distance bound from K006.yaml, not the total PDG rate."
            ),
            diagnostics=diagnostics,
        )
