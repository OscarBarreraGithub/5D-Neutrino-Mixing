"""C008 - LFV semileptonic decay ``D+ -> pi+ e+- mu-+``.

Physics
-------
The Standard Model rate is zero for catalog purposes, so C008 is a pure-NP
upper-bound constraint.  This module reuses the C007 D-to-pi form-factor
machinery and the Phase-4a LFV ``rs_semileptonic_wilsons.lfv_llqq`` block through
``flavor_catalog_constraints.physics_adapters.rare_charm_lfv_semileptonic``.

The constrained prediction is a one-charge-mode smooth full-q2 short-distance
tree-level rate with vector/axial ``C9_LFV + C9p_LFV`` and
``C10_LFV + C10p_LFV``.  The same tree-level e-mu prediction is compared to
both PDG/LHCb charge-mode limits from the C008 YAML, and the larger saturation
is used for the HARD veto.

Severity
--------
HARD.  Both active anchors are observed 90% CL upper limits on charged-LFV
branching fractions.  A point with an evaluated pure-NP prediction above either
charge-mode budget is excluded.

Catalog sidecar
---------------
``flavor_catalog/processes/charm/C008.yaml`` is the source of truth for the
two PDG/LHCb branching-fraction limits and provenance.  Numeric limit values
below are loaded from that sidecar through the scaffold anchor loader.

Phase-4c status
---------------
The tree-level LFV lepton coupling is now read from the Phase-4a lepton-aware
semileptonic Wilson bundle.  For the current diagonal charged-lepton fit it is
rigorously zero, so the tree-level LFV rate is zero and non-vetoing.  Nonzero
tree-level rates require non-diagonal lepton structure; loop-induced LFV is
deferred.  Scalar/tensor, resonance, and LHCb acceptance effects remain
outside this short-distance tree-level rewire.
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
)
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.rare_charm_lfv_semileptonic import (
    RARE_CHARM_DTOPI_EMU_PARAMETRIZATION_CITATION,
    RARE_CHARM_DTOPI_EMU_Q2_TREATMENT_V1,
    RARE_CHARM_LFV_TREE_LEVEL_NOTE_V1,
    dplus_piplus_emu_from_rs_semileptonic_wilsons,
    dplus_piplus_emu_sm,
    rare_charm_dtopi_emu_default_inputs,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "charm"
_REQUIRED_EXTRA = "rs_semileptonic_wilsons"
_OPTIONAL_EW_MASS_EXTRA = "kk_ew_mass_gev"
_EXPECTED_UNITS = "branching fraction"

_EPLUS_MUMINUS_CANDIDATES = ("dplus_piplus_eplus_muminus",)
_EMINUS_MUPLUS_CANDIDATES = ("dplus_piplus_eminus_muplus",)
_LHCB_SCOPE_CANDIDATES = ("lhcb_2021_search_scope",)
_BABAR_PREDECESSOR_CANDIDATES = ("babar_2011_predecessor",)
_PAPER_ERA_SECTION = "paper_era_reference"
_RS_BASELINE_KEY = "rs_baseline"
_BUDGET_SOURCE = (
    "flavor_catalog/processes/charm/C008.yaml "
    "dplus_piplus_eplus_muminus and dplus_piplus_eminus_muplus "
    "(PDG Live/API S031.110/S031.111, LHCb 2021, 90% CL)"
)
_UNEVALUATED_REASON = "missing rs_semileptonic_wilsons LFV llqq block"
_UNEVALUATED_NOTES = f"NOT EVALUATED - {_UNEVALUATED_REASON}"


@dataclass(frozen=True)
class BranchingLimitAnchor:
    """Typed upper-limit branching-fraction anchor."""

    block_key: str
    source: str | None
    year: int | None
    observable: str | None
    value: float
    confidence_level: float
    limit_type: str
    units: str | None
    source_url: str | None
    snapshot_path: str | None
    current_experimental_source: str | None = None
    companion_95cl_limit: float | None = None


@dataclass(frozen=True)
class SearchScopeContext:
    """Typed LHCb search-scope provenance from the C008 sidecar."""

    block_key: str
    source: str | None
    year: int | None
    decay_modes_investigated: float
    integrated_luminosity_fb_inv: float
    dataset_year: int | None
    collision_system: str | None
    no_significant_deviation: bool
    source_url: str | None
    snapshot_path: str | None


@dataclass(frozen=True)
class BaBarPredecessorContext:
    """Typed predecessor C008 limits from BaBar 2011."""

    block_key: str
    source: str | None
    year: int | None
    integrated_luminosity_fb_inv: float
    confidence_level: float
    eplus_muminus_limit: float
    eminus_muplus_limit: float
    units: str | None
    source_url: str | None
    snapshot_path: str | None


@dataclass(frozen=True)
class RSBaselineContext:
    """Typed paper-era RS provenance carried for proxy context."""

    section_key: str
    block_key: str
    source: str | None
    year: int | None
    use: str | None
    source_url: str | None
    snapshot_path: str | None


@dataclass(frozen=True)
class C008Anchor:
    """Typed C008 anchor: two charge-mode limits and provenance."""

    eplus_muminus_limit: BranchingLimitAnchor
    eminus_muplus_limit: BranchingLimitAnchor
    search_scope: SearchScopeContext
    babar_predecessor: BaBarPredecessorContext
    rs_baseline: RSBaselineContext

    @property
    def value(self) -> float:
        """Stricter active upper limit used as the displayed experimental value."""

        return self.active_limit.value

    @property
    def budget(self) -> float:
        """Stricter HARD branching-fraction budget."""

        return self.active_limit.value

    @property
    def active_limit(self) -> BranchingLimitAnchor:
        """Return the tighter of the two current charge-mode limits."""

        return min(
            (self.eplus_muminus_limit, self.eminus_muplus_limit),
            key=lambda item: item.value,
        )


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
            f"{process_id}: C008 anchor field {field_name!r}={value!r} is not numeric"
        ) from exc
    if not math.isfinite(out):
        raise AnchorError(f"{process_id}: C008 anchor field {field_name!r} is not finite")
    return out


def _positive_float(value: Any, *, process_id: str, field_name: str) -> float:
    out = _required_float(value, process_id=process_id, field_name=field_name)
    if out <= 0.0:
        raise AnchorError(f"{process_id}: C008 anchor field {field_name!r} must be positive")
    return out


def _pdg_subblock_for_anchor(anchor: Anchor, *, process_id: str) -> Mapping[str, Any]:
    data = load_full_yaml(process_id, family=_FAMILY)
    pdg_block = data.get("pdg_or_equivalent")
    if not isinstance(pdg_block, Mapping):
        raise AnchorError(
            f"{process_id}: 'pdg_or_equivalent' is not a mapping while loading "
            f"{anchor.block_key}"
        )
    sub = pdg_block.get(anchor.block_key)
    if not isinstance(sub, Mapping):
        present = sorted(str(key) for key in pdg_block)
        raise AnchorError(
            f"{process_id}: load_anchor selected {anchor.block_key!r}, but that "
            f"block is not available as a mapping (present keys: {present})"
        )
    return sub


def _load_limit_anchor(
    candidates: tuple[str, ...],
    *,
    process_id: str,
) -> BranchingLimitAnchor:
    scaffold_anchor = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=candidates,
    )
    sub = _pdg_subblock_for_anchor(scaffold_anchor, process_id=process_id)
    if scaffold_anchor.units != _EXPECTED_UNITS:
        raise AnchorError(
            f"{process_id}: {scaffold_anchor.block_key} must use units "
            f"{_EXPECTED_UNITS!r}, got {scaffold_anchor.units!r}"
        )
    value = _positive_float(
        scaffold_anchor.value,
        process_id=process_id,
        field_name=f"{scaffold_anchor.block_key}.value",
    )
    confidence_level = _positive_float(
        sub.get("confidence_level"),
        process_id=process_id,
        field_name=f"{scaffold_anchor.block_key}.confidence_level",
    )
    limit_type = str(sub.get("limit_type", "")).lower()
    if limit_type != "upper":
        raise AnchorError(f"{process_id}: {scaffold_anchor.block_key}.limit_type must be 'upper'")
    return BranchingLimitAnchor(
        block_key=scaffold_anchor.block_key,
        source=_optional_str(scaffold_anchor.source),
        year=scaffold_anchor.year,
        observable=_optional_str(scaffold_anchor.observable),
        value=float(value),
        confidence_level=float(confidence_level),
        limit_type=limit_type,
        units=scaffold_anchor.units,
        source_url=_optional_str(scaffold_anchor.source_url),
        snapshot_path=_optional_str(scaffold_anchor.snapshot_path),
        current_experimental_source=_optional_str(sub.get("current_experimental_source")),
        companion_95cl_limit=_optional_float(
            sub.get("companion_95cl_limit"),
            process_id=process_id,
            field_name=f"{scaffold_anchor.block_key}.companion_95cl_limit",
        ),
    )


def _load_search_scope(process_id: str) -> SearchScopeContext:
    scope = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=_LHCB_SCOPE_CANDIDATES,
        value_key="decay_modes_investigated",
    )
    sub = _pdg_subblock_for_anchor(scope, process_id=process_id)
    return SearchScopeContext(
        block_key=scope.block_key,
        source=_optional_str(scope.source),
        year=scope.year,
        decay_modes_investigated=_positive_float(
            scope.value,
            process_id=process_id,
            field_name=f"{scope.block_key}.decay_modes_investigated",
        ),
        integrated_luminosity_fb_inv=_positive_float(
            sub.get("integrated_luminosity_fb_inv"),
            process_id=process_id,
            field_name=f"{scope.block_key}.integrated_luminosity_fb_inv",
        ),
        dataset_year=_optional_int(sub.get("dataset_year")),
        collision_system=_optional_str(sub.get("collision_system")),
        no_significant_deviation=bool(sub.get("no_significant_deviation")),
        source_url=_optional_str(scope.source_url),
        snapshot_path=_optional_str(scope.snapshot_path),
    )


def _load_babar_predecessor(process_id: str) -> BaBarPredecessorContext:
    eplus = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=_BABAR_PREDECESSOR_CANDIDATES,
        value_key="dplus_piplus_eplus_muminus_limit",
    )
    eminus = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=_BABAR_PREDECESSOR_CANDIDATES,
        value_key="dplus_piplus_eminus_muplus_limit",
    )
    sub = _pdg_subblock_for_anchor(eplus, process_id=process_id)
    if eplus.units != _EXPECTED_UNITS or eminus.units != _EXPECTED_UNITS:
        raise AnchorError(f"{process_id}: BaBar predecessor limits must be branching fractions")
    return BaBarPredecessorContext(
        block_key=eplus.block_key,
        source=_optional_str(eplus.source),
        year=eplus.year,
        integrated_luminosity_fb_inv=_positive_float(
            sub.get("integrated_luminosity_fb_inv"),
            process_id=process_id,
            field_name=f"{eplus.block_key}.integrated_luminosity_fb_inv",
        ),
        confidence_level=_positive_float(
            sub.get("confidence_level"),
            process_id=process_id,
            field_name=f"{eplus.block_key}.confidence_level",
        ),
        eplus_muminus_limit=_positive_float(
            eplus.value,
            process_id=process_id,
            field_name=f"{eplus.block_key}.dplus_piplus_eplus_muminus_limit",
        ),
        eminus_muplus_limit=_positive_float(
            eminus.value,
            process_id=process_id,
            field_name=f"{eminus.block_key}.dplus_piplus_eminus_muplus_limit",
        ),
        units=eplus.units,
        source_url=_optional_str(eplus.source_url),
        snapshot_path=_optional_str(eplus.snapshot_path),
    )


def _top_level_subblock(
    process_id: str,
    *,
    section_key: str,
    block_key: str,
) -> Mapping[str, Any]:
    data = load_full_yaml(process_id, family=_FAMILY)
    section = data.get(section_key)
    if not isinstance(section, Mapping):
        raise AnchorError(
            f"{process_id}: top-level {section_key!r} is not a mapping while "
            f"loading {block_key!r}"
        )
    sub = section.get(block_key)
    if not isinstance(sub, Mapping):
        present = sorted(str(key) for key in section)
        raise AnchorError(
            f"{process_id}: {section_key}.{block_key} is missing or not a mapping "
            f"(present keys: {present})"
        )
    return sub


def _load_rs_baseline(process_id: str) -> RSBaselineContext:
    sub = _top_level_subblock(
        process_id,
        section_key=_PAPER_ERA_SECTION,
        block_key=_RS_BASELINE_KEY,
    )
    return RSBaselineContext(
        section_key=_PAPER_ERA_SECTION,
        block_key=_RS_BASELINE_KEY,
        source=_optional_str(sub.get("source")),
        year=_optional_int(sub.get("year")),
        use=_optional_str(sub.get("use")),
        source_url=_optional_str(sub.get("source_url")),
        snapshot_path=_optional_str(sub.get("snapshot_path")),
    )


def _load_c008_anchor(process_id: str) -> C008Anchor:
    anchor = C008Anchor(
        eplus_muminus_limit=_load_limit_anchor(
            _EPLUS_MUMINUS_CANDIDATES,
            process_id=process_id,
        ),
        eminus_muplus_limit=_load_limit_anchor(
            _EMINUS_MUPLUS_CANDIDATES,
            process_id=process_id,
        ),
        search_scope=_load_search_scope(process_id),
        babar_predecessor=_load_babar_predecessor(process_id),
        rs_baseline=_load_rs_baseline(process_id),
    )
    if anchor.budget <= 0.0:
        raise AnchorError(f"{process_id}: C008 HARD budget must be positive")
    return anchor


def _complex_wilsons(result: Any) -> dict[str, complex]:
    if result.wilsons is None:
        return {}
    return {key: complex(value) for key, value in result.wilsons.wilsons.items()}


@register
class Constraint:
    """Catalogued pure-NP ``D+ -> pi+ e+- mu-+`` LFV constraint."""

    process_id = "C008"
    severity = Severity.HARD
    observable = "BR(D+ -> pi+ e+- mu-+) LFV full-q2 tree-level"

    def __init__(self) -> None:
        self.anchor = _load_c008_anchor(self.process_id)
        self.sd_inputs = rare_charm_dtopi_emu_default_inputs()
        self.sm_result = dplus_piplus_emu_sm(self.sd_inputs)

    def _base_diagnostics(self) -> dict[str, Any]:
        scope = self.anchor.search_scope
        babar = self.anchor.babar_predecessor
        rs = self.anchor.rs_baseline
        return {
            "eplus_muminus_block": self.anchor.eplus_muminus_limit.block_key,
            "eminus_muplus_block": self.anchor.eminus_muplus_limit.block_key,
            "eplus_muminus_90cl_limit": float(self.anchor.eplus_muminus_limit.value),
            "eminus_muplus_90cl_limit": float(self.anchor.eminus_muplus_limit.value),
            "active_limit_block": self.anchor.active_limit.block_key,
            "active_limit_is_strictest_charge_mode": True,
            "experimental_confidence_level": float(
                self.anchor.active_limit.confidence_level
            ),
            "eplus_muminus_95cl_limit": float(
                self.anchor.eplus_muminus_limit.companion_95cl_limit or 0.0
            ),
            "eminus_muplus_95cl_limit": float(
                self.anchor.eminus_muplus_limit.companion_95cl_limit or 0.0
            ),
            "lhcb_current_experimental_source": (
                self.anchor.active_limit.current_experimental_source
            ),
            "lhcb_2021_search_scope_block": scope.block_key,
            "lhcb_2021_decay_modes_investigated": float(
                scope.decay_modes_investigated
            ),
            "lhcb_2021_integrated_luminosity_fb_inv": float(
                scope.integrated_luminosity_fb_inv
            ),
            "lhcb_2021_dataset_year": scope.dataset_year,
            "lhcb_2021_collision_system": scope.collision_system,
            "lhcb_2021_no_significant_deviation": scope.no_significant_deviation,
            "babar_2011_block": babar.block_key,
            "babar_2011_integrated_luminosity_fb_inv": float(
                babar.integrated_luminosity_fb_inv
            ),
            "babar_2011_eplus_muminus_90cl_limit": float(
                babar.eplus_muminus_limit
            ),
            "babar_2011_eminus_muplus_90cl_limit": float(
                babar.eminus_muplus_limit
            ),
            "budget_source": _BUDGET_SOURCE,
            "budget_is_upper_limit": True,
            "sm_branching_fraction": 0.0,
            "sm_lfv_policy": (
                "D+ -> pi+ e mu is charged-LFV and has zero SM rate for "
                "catalog purposes; the HARD budget is applied to the pure-NP "
                "tree-level prediction when the LFV llqq block is present."
            ),
            "charge_mode_prediction_policy": (
                "The tree-level lfv_llqq e-mu block is not orientation-"
                "specific; the same one-charge-mode prediction is compared "
                "to both C008 limits."
            ),
            "full_q2_tree_level_prediction_constrained": True,
            "q2_treatment": RARE_CHARM_DTOPI_EMU_Q2_TREATMENT_V1,
            "lhcb_window_acceptance_applied": False,
            "resonance_amplitudes_included": False,
            "scalar_tensor_operators_included": False,
            "rs_baseline_block": f"{rs.section_key}.{rs.block_key}",
            "rs_baseline_source": rs.source,
            "rs_baseline_year": rs.year,
            "rs_baseline_use": rs.use,
            "parametrization_citation": RARE_CHARM_DTOPI_EMU_PARAMETRIZATION_CITATION,
            "lfv_tree_level_note": RARE_CHARM_LFV_TREE_LEVEL_NOTE_V1,
            "loop_lfv_status": "loop_induced_lfv_deferred",
            "deferred_short_distance_effects": (
                "scalar/tensor operators, resonance amplitudes, and LHCb "
                "window/acceptance matching"
            ),
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
                    "non-vetoing only; no D+ -> pi+ e mu NP prediction was evaluated"
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
            result = dplus_piplus_emu_from_rs_semileptonic_wilsons(
                wilson_input,
                matching_scale_gev=None if kk_ew_mass is None else float(kk_ew_mass),
                inputs=self.sd_inputs,
                charge_mode="eplus_muminus",
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
        eplus_ratio = predicted / self.anchor.eplus_muminus_limit.value
        eminus_ratio = predicted / self.anchor.eminus_muplus_limit.value
        ratio = max(eplus_ratio, eminus_ratio)
        budget = float(self.anchor.budget)
        diagnostics = dict(result.diagnostics)
        diagnostics.update(
            {
                "evaluated": True,
                **self._base_diagnostics(),
                "prediction_per_charge_mode_branching_fraction": predicted,
                "eplus_muminus_ratio": float(eplus_ratio),
                "eminus_muplus_ratio": float(eminus_ratio),
                "active_ratio_is_max_charge_mode_saturation": True,
                "np_shift_branching_fraction": float(
                    result.np_shift_branching_fraction
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
                "Pure-NP BR(D+ -> pi+ e+- mu-+) bound using the shared C007 "
                "D->pi form-factor machinery and Phase-4a LFV llqq Wilsons. "
                "The same one-charge-mode full-q2 short-distance prediction is "
                "compared to both C008 YAML 90% CL charge-mode limits; the "
                "larger saturation sets the HARD ratio. Tree-level LFV is "
                "rigorous and zero for the diagonal charged-lepton fit; "
                "loop-induced LFV is deferred."
            ),
            diagnostics=diagnostics,
        )
