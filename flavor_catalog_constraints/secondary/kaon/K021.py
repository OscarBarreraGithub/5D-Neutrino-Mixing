"""K021 - LFV neutral-kaon semileptonic decay ``K_L -> pi0 e+- mu-+``.

Physics
-------
``K_L -> pi0 e mu`` is charged-lepton-flavor violating, so the Standard Model
rate is zero for catalog purposes.  K021 applies a pure-NP HARD upper bound,

    BR_NP(K_L -> pi0 e+- mu-+) <= BR_exp^upper.

The short-distance tree-level prediction reuses the shared ``s -> d l l``
rare-kaon machinery and the Phase-4a LFV
``rs_semileptonic_wilsons.lfv_llqq`` e-mu block through
``flavor_catalog_constraints.physics_adapters.rare_kaon_lfv_semileptonic``.
The K021 neutral-mode adapter helper is append-only relative to K020: it keeps
the same K->pi vector-current form-factor and three-body q2/angular integral,
but supplies the K_L mass/lifetime, pi0 mass, neutral K_l3 f_+(0)
normalization, and the charge-summed e+- mu-+ final state.

Severity
--------
HARD.  The K021 YAML anchor is a direct 90% CL upper limit on the summed
charged-LFV branching fraction from the PDG K_L listing / KTeV result.  A
point with an evaluated pure-NP prediction above that limit is excluded.

Catalog sidecar
---------------
``flavor_catalog/processes/secondary/kaon/K021.yaml`` is the source of truth
for the PDG/KTeV branching-fraction limit and provenance.  Numeric values
below are loaded through the scaffold anchor loader, not hardcoded here.

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
    RARE_KAON_KL_PI0_EMU_PARAMETRIZATION_CITATION,
    RARE_KAON_KL_PI0_EMU_Q2_TREATMENT_V1,
    RARE_KAON_LFV_TREE_LEVEL_NOTE_V1,
    klong_pi0_emu_from_rs_semileptonic_wilsons,
    klong_pi0_emu_sm,
    rare_kaon_klong_pi0_emu_default_inputs,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "kaon"
_TIER = ConstraintLevel.SECONDARY
_REQUIRED_EXTRA = "rs_semileptonic_wilsons"
_OPTIONAL_EW_MASS_EXTRA = "kk_ew_mass_gev"
_CURRENT_LIMIT_VALUE_ID = "PDG2025:K021:KL_pi0emu_limit"
_KTEV_VALUE_ID = "KTeV2008:K021:KL_pi0mue_limit"
_CURRENT_OBSERVABLE = "BR(K_L -> pi0 e^+- mu^-+), summed charge states"
_EXPECTED_LIMIT_TYPE = "upper_limit"
_EXPECTED_CL = "90% CL"
_EXPECTED_UNITS = "dimensionless branching fraction"
_CHARGE_MODE = "muplus_eminus"
_BUDGET_SOURCE = (
    "flavor_catalog/processes/secondary/kaon/K021.yaml "
    "PDG2025:K021:KL_pi0emu_limit "
    "(PDG 2025 K_L listing / KTeV 90% CL limit, charge-summed)"
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
    note: str | None = None


@dataclass(frozen=True)
class K021Anchor:
    """Typed K021 anchor bundle."""

    current_limit: BranchingLimitAnchor
    ktev_limit: BranchingLimitAnchor

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
            f"{process_id}: K021 anchor field {field_name!r}={value!r} "
            "is not numeric"
        ) from exc
    if not math.isfinite(out):
        raise AnchorError(f"{process_id}: K021 anchor field {field_name!r} is not finite")
    return out


def _positive_float(value: Any, *, process_id: str, field_name: str) -> float:
    out = _required_float(value, process_id=process_id, field_name=field_name)
    if out <= 0.0:
        raise AnchorError(f"{process_id}: K021 anchor field {field_name!r} must be positive")
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
    entries = _pdg_value_entries(process_id)
    matches = [
        (index, entry)
        for index, entry in enumerate(entries)
        if entry.get("value_id") == value_id
    ]
    if not matches:
        present = [
            str(entry.get("value_id"))
            for entry in entries
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
            f"expected {block_key!r} for K021 value_id {value_id!r}"
        )
    return scaffold_anchor, entry, index


def _validate_limit_entry(
    entry: Mapping[str, Any],
    *,
    process_id: str,
    block_key: str,
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
    if units != _EXPECTED_UNITS:
        raise AnchorError(
            f"{process_id}: {block_key}.units must be {_EXPECTED_UNITS!r}, "
            f"got {units!r}"
        )


def _load_limit_anchor(
    value_id: str,
    *,
    process_id: str,
) -> BranchingLimitAnchor:
    anchor, entry, index = _load_scaffold_value_anchor(
        value_id,
        process_id=process_id,
    )
    block_key = f"pdg_or_equivalent.values[{index}]"
    if entry.get("value_id") != value_id:
        raise AnchorError(f"{process_id}: selected wrong K021 value_id")
    _validate_limit_entry(entry, process_id=process_id, block_key=block_key)
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
        note=_optional_str(entry.get("charge_state_note") or entry.get("note")),
    )


def _load_k021_anchor(process_id: str) -> K021Anchor:
    current = _load_limit_anchor(_CURRENT_LIMIT_VALUE_ID, process_id=process_id)
    ktev = _load_limit_anchor(_KTEV_VALUE_ID, process_id=process_id)
    if current.observable != _CURRENT_OBSERVABLE:
        raise AnchorError(
            f"{process_id}: current K021 limit observable {current.observable!r} "
            f"does not match {_CURRENT_OBSERVABLE!r}"
        )
    if current.value <= 0.0 or ktev.value <= 0.0:
        raise AnchorError(f"{process_id}: K021 HARD budget/provenance must be positive")
    return K021Anchor(current_limit=current, ktev_limit=ktev)


def _complex_wilsons(result: Any) -> dict[str, complex]:
    if result.wilsons is None:
        return {}
    return {key: complex(value) for key, value in result.wilsons.wilsons.items()}


@register
class Constraint:
    """Catalogued pure-NP ``K_L -> pi0 e+- mu-+`` LFV constraint."""

    process_id = "K021"
    severity = Severity.HARD
    observable = "BR(K_L -> pi0 e+- mu-+) LFV"

    def __init__(self) -> None:
        self.anchor = _load_k021_anchor(self.process_id)
        self.sm_inputs = rare_kaon_klong_pi0_emu_default_inputs()
        self.sm_result = klong_pi0_emu_sm(
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
            "experimental_charge_state_note": self.anchor.current_limit.note,
            "ktev_value_id": self.anchor.ktev_limit.value_id,
            "ktev_limit": float(self.anchor.ktev_limit.value),
            "ktev_source": self.anchor.ktev_limit.source,
            "ktev_source_url": self.anchor.ktev_limit.source_url,
            "budget_source": _BUDGET_SOURCE,
            "budget_is_upper_limit": True,
            "sm_branching_fraction": 0.0,
            "sm_lfv_policy": (
                "K_L -> pi0 e mu is charged-LFV and has zero SM rate for "
                "catalog purposes; the HARD budget is applied to the pure-NP "
                "short-distance tree-level prediction when the LFV llqq block "
                "is present."
            ),
            "charge_mode": _CHARGE_MODE,
            "charge_conjugate_modes_included": True,
            "summed_charge_states": "e+- mu-+",
            "parametrization_citation": (
                RARE_KAON_KL_PI0_EMU_PARAMETRIZATION_CITATION
            ),
            "q2_treatment": RARE_KAON_KL_PI0_EMU_Q2_TREATMENT_V1,
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
                    "non-vetoing only; no K_L -> pi0 e+- mu-+ NP prediction "
                    "was evaluated"
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
            result = klong_pi0_emu_from_rs_semileptonic_wilsons(
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
                "Pure-NP BR(K_L -> pi0 e+- mu-+) bound using Phase-4a LFV "
                "llqq Wilsons additively mapped into the neutral rare-kaon "
                "K->pi inputs and the YAML charge-summed limit. Tree-level "
                "LFV is rigorous and zero for the diagonal charged-lepton "
                "fit; loop-induced LFV is deferred."
            ),
            diagnostics=diagnostics,
        )
