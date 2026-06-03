"""K020 - LFV charged-kaon semileptonic decay ``K+ -> pi+ mu+ e-``.

Physics
-------
``K+ -> pi+ e mu`` is charged-lepton-flavor violating, so the Standard Model
rate is zero for catalog purposes.  K020 applies a pure-NP HARD upper bound,

    BR_NP(K+ -> pi+ mu+ e-) <= BR_exp^upper.

The short-distance tree-level prediction reuses the shared ``s -> d l l``
rare-kaon machinery and the Phase-4a LFV
``rs_semileptonic_wilsons.lfv_llqq`` e-mu block through
``flavor_catalog_constraints.physics_adapters.rare_kaon_lfv_semileptonic``.
That adapter adds the exclusive ``K -> pi`` vector form factor and full
three-body phase-space integration for one charge mode.

Severity
--------
HARD.  The K020 YAML anchor is a direct 90% CL upper limit on a charged-LFV
branching fraction from the PDG K+ listing / BNL E865+E777 combination.  A
point with an evaluated pure-NP prediction above that limit is excluded.

Catalog sidecar
---------------
``flavor_catalog/processes/secondary/kaon/K020.yaml`` is the source of truth
for the PDG/BNL branching-fraction limit and provenance.  Numeric values below
are loaded through the scaffold anchor loader, not hardcoded here.

Phase-4c status
---------------
The tree-level LFV lepton coupling is now read from the Phase-4a lepton-aware
semileptonic Wilson bundle.  For the current diagonal charged-lepton fit it is
rigorously zero, so the tree-level LFV rate is zero and non-vetoing.  Nonzero
tree-level rates require non-diagonal lepton structure; loop-induced LFV is
deferred.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Mapping

from flavor_catalog_constraints import anchors as anchor_scaffold
from flavor_catalog_constraints.anchors import (
    Anchor,
    AnchorError,
    load_anchor,
    load_full_yaml,
)
from flavor_catalog_constraints.base import (
    ConstraintLevel,
    ConstraintResult,
    ParameterPoint,
    Severity,
)
from flavor_catalog_constraints.physics_adapters.rare_kaon_lfv_semileptonic import (
    RARE_KAON_KTOPI_EMU_PARAMETRIZATION_CITATION,
    RARE_KAON_KTOPI_EMU_Q2_TREATMENT_V1,
    RARE_KAON_LFV_TREE_LEVEL_NOTE_V1,
    kplus_piplus_emu_from_rs_semileptonic_wilsons,
    kplus_piplus_emu_sm,
    rare_kaon_ktopi_emu_default_inputs,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "kaon"
_TIER = ConstraintLevel.SECONDARY
_REQUIRED_EXTRA = "rs_semileptonic_wilsons"
_OPTIONAL_EW_MASS_EXTRA = "kk_ew_mass_gev"
_PDG_VALUES_SECTION = "pdg_or_equivalent.values"
_CURRENT_LIMIT_VALUE_ID = "PDG2025:K020:Kplus_piplus_mup_em_limit"
_SHER_E865_VALUE_ID = "Sher2005:K020:Kplus_piplus_mup_em_e865_only_limit"
_OPPOSITE_CHARGE_VALUE_ID = "PDG2025:K020:Kplus_piplus_mum_ep_limit"
_CURRENT_OBSERVABLE = "BR(K+ -> pi+ mu+ e-)"
_EXPECTED_LIMIT_TYPE = "upper_limit"
_EXPECTED_CL = "90% CL"
_EXPECTED_UNITS = "dimensionless branching fraction"
_LOOSE_DIMENSIONLESS_UNITS = "(dimensionless)"
_CHARGE_MODE = "muplus_eminus"
_BUDGET_SOURCE = (
    "flavor_catalog/processes/secondary/kaon/K020.yaml "
    "PDG2025:K020:Kplus_piplus_mup_em_limit "
    "(PDG 2025 K+ listing / BNL E865+E777 90% CL limit)"
)
_UNEVALUATED_REASON = "missing rs_semileptonic_wilsons LFV llqq block"
_UNEVALUATED_NOTES = f"NOT EVALUATED - {_UNEVALUATED_REASON}"


@dataclass(frozen=True)
class BranchingLimitAnchor:
    """Typed upper-limit branching-fraction anchor."""

    block_key: str
    value_id: str
    observable: str | None
    source: str | None
    year: int | None
    value: float
    limit_type: str
    confidence_level: str | None
    units: str | None
    source_url: str | None
    snapshot_path: str | None
    display: str | None = None
    experimental_reference: str | None = None
    doi: str | None = None
    arxiv_id: str | None = None
    note: str | None = None


@dataclass(frozen=True)
class K020Anchor:
    """Typed K020 anchor bundle."""

    current_limit: BranchingLimitAnchor
    sher_e865_limit: BranchingLimitAnchor
    opposite_charge_limit: BranchingLimitAnchor

    @property
    def value(self) -> float:
        """Upper limit used as the experimental value."""

        return self.current_limit.value

    @property
    def budget(self) -> float:
        """HARD pure-NP branching-fraction budget."""

        return self.current_limit.value


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
            f"{process_id}: K020 anchor field {field_name!r}={value!r} "
            "is not numeric"
        ) from exc
    if not math.isfinite(out):
        raise AnchorError(f"{process_id}: K020 anchor field {field_name!r} is not finite")
    return out


def _positive_float(value: Any, *, process_id: str, field_name: str) -> float:
    out = _required_float(value, process_id=process_id, field_name=field_name)
    if out <= 0.0:
        raise AnchorError(f"{process_id}: K020 anchor field {field_name!r} must be positive")
    return out


def _pdg_value_entries(process_id: str) -> list[Mapping[str, Any]]:
    data = load_full_yaml(process_id, family=_FAMILY, tier=_TIER)
    pdg = data.get("pdg_or_equivalent")
    if not isinstance(pdg, Mapping):
        raise AnchorError(f"{process_id}: expected pdg_or_equivalent mapping")
    values = pdg.get("values")
    if not isinstance(values, list) or not values:
        raise AnchorError(f"{process_id}: expected non-empty pdg_or_equivalent.values")
    for index, entry in enumerate(values):
        if not isinstance(entry, Mapping):
            raise AnchorError(
                f"{process_id}: pdg_or_equivalent.values[{index}] is not a mapping"
            )
    return values


def _entry_by_value_id(
    process_id: str,
    value_id: str,
) -> tuple[int, Mapping[str, Any]]:
    matches = [
        (index, entry)
        for index, entry in enumerate(_pdg_value_entries(process_id))
        if entry.get("value_id") == value_id
    ]
    if not matches:
        present = [
            str(entry.get("value_id"))
            for entry in _pdg_value_entries(process_id)
            if entry.get("value_id") is not None
        ]
        raise AnchorError(
            f"{process_id}: value_id {value_id!r} not found in "
            f"pdg_or_equivalent.values (present value_ids: {present})"
        )
    if len(matches) > 1:
        raise AnchorError(f"{process_id}: duplicate value_id {value_id!r}")
    return matches[0]


def _load_scaffold_value_anchor(
    value_id: str,
    *,
    process_id: str,
) -> tuple[Anchor, Mapping[str, Any], int]:
    index, entry = _entry_by_value_id(process_id, value_id)
    block_key = f"pdg_or_equivalent.values[{index}]"
    virtual_block = {block_key: dict(entry)}
    original_load_pdg_block = anchor_scaffold.load_pdg_block

    def _load_virtual_pdg_block(
        request_process_id: str,
        **kwargs: Any,
    ) -> Mapping[str, Any]:
        if (
            request_process_id == process_id
            and kwargs.get("family") == _FAMILY
            and kwargs.get("tier", ConstraintLevel.PRIMARY) == _TIER
        ):
            return virtual_block
        return original_load_pdg_block(request_process_id, **kwargs)

    anchor_scaffold.load_pdg_block = _load_virtual_pdg_block
    try:
        scaffold_anchor = load_anchor(
            process_id,
            family=_FAMILY,
            tier=_TIER,
            candidates=(block_key,),
        )
    finally:
        anchor_scaffold.load_pdg_block = original_load_pdg_block

    if scaffold_anchor.block_key != block_key:
        raise AnchorError(
            f"{process_id}: load_anchor selected {scaffold_anchor.block_key!r}, "
            f"expected {block_key!r} for K020 value_id {value_id!r}"
        )
    return scaffold_anchor, entry, index


def _validate_limit_entry(
    entry: Mapping[str, Any],
    *,
    process_id: str,
    block_key: str,
    expected_units: tuple[str, ...],
) -> None:
    limit_type = _optional_str(entry.get("limit_type"))
    if limit_type != _EXPECTED_LIMIT_TYPE:
        raise AnchorError(
            f"{process_id}: {block_key}.limit_type must be "
            f"{_EXPECTED_LIMIT_TYPE!r}, got {limit_type!r}"
        )
    confidence_level = _optional_str(entry.get("cl"))
    if confidence_level != _EXPECTED_CL:
        raise AnchorError(
            f"{process_id}: {block_key}.cl must be {_EXPECTED_CL!r}, "
            f"got {confidence_level!r}"
        )
    units = _optional_str(entry.get("units"))
    if units not in expected_units:
        raise AnchorError(
            f"{process_id}: {block_key}.units must be one of {expected_units!r}, "
            f"got {units!r}"
        )


def _load_limit_anchor(
    value_id: str,
    *,
    process_id: str,
    expected_units: tuple[str, ...] = (_EXPECTED_UNITS,),
) -> BranchingLimitAnchor:
    anchor, entry, index = _load_scaffold_value_anchor(
        value_id,
        process_id=process_id,
    )
    block_key = f"pdg_or_equivalent.values[{index}]"
    if entry.get("value_id") != value_id:
        raise AnchorError(f"{process_id}: selected wrong K020 value_id")
    _validate_limit_entry(
        entry,
        process_id=process_id,
        block_key=block_key,
        expected_units=expected_units,
    )
    return BranchingLimitAnchor(
        block_key=anchor.block_key,
        value_id=value_id,
        observable=_optional_str(anchor.observable),
        source=_optional_str(anchor.source or entry.get("source")),
        year=anchor.year if anchor.year is not None else _optional_int(entry.get("year")),
        value=_positive_float(
            anchor.value,
            process_id=process_id,
            field_name=f"{block_key}.value",
        ),
        limit_type=str(entry.get("limit_type")),
        confidence_level=_optional_str(entry.get("cl")),
        units=_optional_str(anchor.units or entry.get("units")),
        source_url=_optional_str(anchor.source_url or entry.get("source_url")),
        snapshot_path=_optional_str(anchor.snapshot_path or entry.get("snapshot_path")),
        display=_optional_str(entry.get("display")),
        experimental_reference=_optional_str(entry.get("experimental_reference")),
        doi=_optional_str(entry.get("doi")),
        arxiv_id=_optional_str(entry.get("arxiv_id")),
        note=_optional_str(entry.get("note")),
    )


def _load_k020_anchor(process_id: str) -> K020Anchor:
    current = _load_limit_anchor(_CURRENT_LIMIT_VALUE_ID, process_id=process_id)
    sher = _load_limit_anchor(
        _SHER_E865_VALUE_ID,
        process_id=process_id,
        expected_units=(_EXPECTED_UNITS, _LOOSE_DIMENSIONLESS_UNITS),
    )
    opposite = _load_limit_anchor(_OPPOSITE_CHARGE_VALUE_ID, process_id=process_id)
    if current.observable != _CURRENT_OBSERVABLE:
        raise AnchorError(
            f"{process_id}: current K020 limit observable {current.observable!r} "
            f"does not match {_CURRENT_OBSERVABLE!r}"
        )
    if current.value > sher.value:
        raise AnchorError(
            f"{process_id}: PDG combined K020 limit should not be looser than "
            "the Sher/E865-only provenance row"
        )
    if current.value <= 0.0:
        raise AnchorError(f"{process_id}: K020 HARD budget must be positive")
    return K020Anchor(
        current_limit=current,
        sher_e865_limit=sher,
        opposite_charge_limit=opposite,
    )


def _complex_wilsons(result: Any) -> dict[str, complex]:
    if result.wilsons is None:
        return {}
    return {key: complex(value) for key, value in result.wilsons.wilsons.items()}


@register
class Constraint:
    """Catalogued pure-NP ``K+ -> pi+ mu+ e-`` LFV constraint."""

    process_id = "K020"
    severity = Severity.HARD
    observable = "BR(K+ -> pi+ mu+ e-) LFV"

    def __init__(self) -> None:
        self.anchor = _load_k020_anchor(self.process_id)
        self.sm_inputs = rare_kaon_ktopi_emu_default_inputs()
        self.sm_result = kplus_piplus_emu_sm(
            self.sm_inputs,
            charge_mode=_CHARGE_MODE,
        )

    def _base_diagnostics(self) -> dict[str, Any]:
        return {
            "experimental_block": self.anchor.current_limit.block_key,
            "experimental_value_id": self.anchor.current_limit.value_id,
            "experimental_confidence_level": self.anchor.current_limit.confidence_level,
            "experimental_source": self.anchor.current_limit.source,
            "experimental_source_url": self.anchor.current_limit.source_url,
            "sher_e865_value_id": self.anchor.sher_e865_limit.value_id,
            "sher_e865_limit": float(self.anchor.sher_e865_limit.value),
            "sher_e865_source": self.anchor.sher_e865_limit.source,
            "sher_e865_source_url": self.anchor.sher_e865_limit.source_url,
            "opposite_charge_value_id": self.anchor.opposite_charge_limit.value_id,
            "opposite_charge_limit": float(self.anchor.opposite_charge_limit.value),
            "opposite_charge_source": self.anchor.opposite_charge_limit.source,
            "budget_source": _BUDGET_SOURCE,
            "budget_is_upper_limit": True,
            "sm_branching_fraction": 0.0,
            "sm_lfv_policy": (
                "K+ -> pi+ e mu is charged-LFV and has zero SM rate for "
                "catalog purposes; the HARD budget is applied to the pure-NP "
                "short-distance tree-level prediction when the LFV llqq block "
                "is present."
            ),
            "charge_mode": _CHARGE_MODE,
            "charge_conjugate_modes_included": False,
            "parametrization_citation": (
                RARE_KAON_KTOPI_EMU_PARAMETRIZATION_CITATION
            ),
            "q2_treatment": RARE_KAON_KTOPI_EMU_Q2_TREATMENT_V1,
            "lfv_tree_level_note": RARE_KAON_LFV_TREE_LEVEL_NOTE_V1,
            "loop_lfv_status": "loop_induced_lfv_deferred",
        }

    def _unevaluated_result(self, diagnostics: Mapping[str, Any]) -> ConstraintResult:
        return ConstraintResult(
            process_id=self.process_id,
            severity=self.severity,
            passes=True,
            predicted=None,
            sm_prediction=0.0,
            experimental=float(self.anchor.value),
            ratio=None,
            budget=float(self.anchor.budget),
            notes=_UNEVALUATED_NOTES,
            diagnostics={
                "evaluated": False,
                "unevaluated_reason": _UNEVALUATED_REASON,
                "passes_semantics": (
                    "non-vetoing only; no K+ -> pi+ mu+ e- NP prediction was "
                    "evaluated"
                ),
                **self._base_diagnostics(),
                **dict(diagnostics),
            },
        )

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        wilson_input = point.get_extra(_REQUIRED_EXTRA)
        if wilson_input is None:
            return self._unevaluated_result({"missing_extra": _REQUIRED_EXTRA})

        kk_ew_mass = point.get_extra(_OPTIONAL_EW_MASS_EXTRA)
        try:
            result = kplus_piplus_emu_from_rs_semileptonic_wilsons(
                wilson_input,
                matching_scale_gev=None if kk_ew_mass is None else float(kk_ew_mass),
                inputs=self.sm_inputs,
                charge_mode=_CHARGE_MODE,
            )
        except (AttributeError, KeyError, TypeError, ValueError) as exc:
            return self._unevaluated_result(
                {
                    "invalid_extra": _REQUIRED_EXTRA,
                    "exception_type": type(exc).__name__,
                    "exception": str(exc),
                }
            )

        predicted = float(result.branching_fraction)
        budget = float(self.anchor.budget)
        ratio = predicted / budget if budget > 0.0 else float("inf")
        diagnostics = dict(result.diagnostics)
        diagnostics.update(
            {
                "evaluated": True,
                **self._base_diagnostics(),
                "np_shift_branching_fraction": float(
                    result.np_shift_branching_fraction
                ),
                "sm_formula_branching_fraction": float(
                    result.sm_branching_fraction
                ),
                "kk_ew_mass_extra_used": kk_ew_mass is not None,
                "wilson_coefficients": _complex_wilsons(result),
            }
        )

        return ConstraintResult(
            process_id=self.process_id,
            severity=self.severity,
            passes=bool(ratio <= 1.0),
            predicted=predicted,
            sm_prediction=0.0,
            experimental=float(self.anchor.value),
            ratio=float(ratio),
            budget=budget,
            notes=(
                "Pure-NP BR(K+ -> pi+ mu+ e-) bound using Phase-4a LFV llqq "
                "Wilsons additively mapped into the rare-kaon K->pi inputs. "
                "Tree-level LFV is rigorous and zero for the diagonal "
                "charged-lepton fit; loop-induced LFV is deferred."
            ),
            diagnostics=diagnostics,
        )
