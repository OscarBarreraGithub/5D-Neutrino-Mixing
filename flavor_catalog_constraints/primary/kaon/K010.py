"""K010 - ``K_S -> pi0 e+ e-``.

Physics
-------
The K_S electron mode is a long-distance, CP-conserving rare-kaon input that
sets the indirect-CP normalization used in K_L -> pi0 e+e- analyses.  K010 uses
the electron-mode a_S parametrization

    BR(K_S -> pi0 e+e-) = 5.2e-9 |a_S^eff|^2

for the full kinematic region.  The central ``|a_S|`` is extracted from the
K010 full-region YAML anchor, while a possible RS contribution is a documented
shift of the CP-conserving vector amplitude from the Phase-3a s -> d e e
C9/C9p Wilsons mapped into the K008 y7V slot.

Severity
--------
HARD.  The predicted full-region branching fraction is compared with the
K010.yaml extrapolated-total NA48/PDG value using its asymmetric experimental
uncertainty.  The directly measured m_ee > 0.165 GeV partial rate is retained
as a diagnostic anchor because the full-rate number assumes a vector matrix
element and unit form factor.

Catalog sidecar
---------------
``flavor_catalog/processes/kaon/K010.yaml`` is the source of truth for the
partial-rate measurement, full-region extrapolation, event counts, and budget.

NEEDS-HUMAN-PHYSICS
-------------------
The semileptonic short-distance input is the Phase-3a light-Z Wilson path.
Complete K010 rigor still needs the K_S long-distance form factor and
EW KK/Z/Z'/photon-penguin/electron-sector inputs not available in Phase 3a.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Mapping

from flavor_catalog_constraints import anchors as anchor_scaffold
from flavor_catalog_constraints.anchors import Anchor, AnchorError, load_anchor, load_full_yaml
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.rare_kaon_dilepton_ks import (
    KSHORT_PI0EE_A_S_BRANCHING_COEFFICIENT,
    KShortPi0EEChPTInputs,
    RARE_KAON_KS_PI0EE_PARAMETRIZATION_CITATION,
    RARE_KAON_KS_PI0EE_RS_MATCHING_ASSUMPTION_V1,
    kshort_pi0ee_a_s_from_rs_semileptonic_wilsons,
    kshort_pi0ee_a_s_sm,
)
from flavor_catalog_constraints.physics_adapters.rare_kaon_dilepton import (
    RARE_KAON_RS_SEMILEPTONIC_VECTOR_MATCHING_STATUS_V1,
    RARE_KAON_PI0EE_SEMILEPTONIC_RUNNING_DIAGNOSTIC_V1,
    RS_SEMILEPTONIC_MATCHING_ASSUMPTION_V1,
    rare_kaon_dilepton_default_sm_inputs,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "kaon"
_REQUIRED_EXTRA = "rs_semileptonic_wilsons"
_OPTIONAL_EW_MASS_EXTRA = "kk_ew_mass_gev"
_SCAFFOLD_UNCERTAINTY_KEY = "__k010_uncertainty_is_parsed_below__"
_EXPECTED_BRANCHING_UNITS = "branching fraction"
_PARTIAL_RATE_BLOCK = "canonical_measurement"
_FULL_RATE_BLOCK = "extrapolated_total_rate"
_VALUES_SECTION = "pdg_or_equivalent.values"
_OBSERVED_CANDIDATES_OBS = "NA48/1 observed signal candidates"
_EXPECTED_BACKGROUND_OBS = "NA48/1 expected background"
_BUDGET_SOURCE = (
    "flavor_catalog/processes/kaon/K010.yaml extrapolated_total_rate "
    "PDG/NA48 full-region branching fraction"
)


@dataclass(frozen=True)
class BranchingFractionAnchor:
    """K010 branching-fraction value with asymmetric errors."""

    block_key: str
    source: str | None
    source_key: str | None
    year: int | None
    value: float
    uncertainty_upper: float
    uncertainty_lower: float
    units: str | None
    source_url: str | None
    snapshot_path: str | None
    assumptions: str | None = None

    @property
    def uncertainty(self) -> float:
        return float(0.5 * (self.uncertainty_upper + self.uncertainty_lower))


@dataclass(frozen=True)
class CatalogValueAnchor:
    """Value-bearing K010 YAML list entry routed through ``load_anchor``."""

    section_key: str
    entry_index: int
    block_key: str
    observable: str | None
    source: str | None
    source_key: str | None
    year: int | None
    value: float
    uncertainty: float | None
    units: str | None
    source_url: str | None
    snapshot_path: str | None
    value_raw: str | None


@dataclass(frozen=True)
class K010BudgetBand:
    """Direction-aware K010 full-rate budget."""

    source: str
    construction: str
    hard_veto_budget: float
    experimental_sigma_upper: float
    experimental_sigma_lower: float
    lower_edge: float
    upper_edge: float
    sm_subtracted: bool


@dataclass(frozen=True)
class K010Anchor:
    """Typed K010 anchor bundle."""

    partial_rate: BranchingFractionAnchor
    full_rate: BranchingFractionAnchor
    observed_candidates: CatalogValueAnchor
    expected_background: CatalogValueAnchor
    a_s_inputs: KShortPi0EEChPTInputs
    budget_band: K010BudgetBand

    @property
    def value(self) -> float:
        """Full-region branching fraction used by the HARD K010 verdict."""
        return self.full_rate.value

    @property
    def uncertainty(self) -> float:
        return self.full_rate.uncertainty

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
            f"{process_id}: K010 anchor field {field_name!r}={value!r} "
            "is not numeric"
        ) from exc
    if not math.isfinite(number):
        raise AnchorError(
            f"{process_id}: K010 anchor field {field_name!r}={value!r} is not finite"
        )
    return number


def _positive_float(value: Any, *, process_id: str, field_name: str) -> float:
    number = _required_float(value, process_id=process_id, field_name=field_name)
    if number <= 0.0:
        raise AnchorError(
            f"{process_id}: K010 anchor field {field_name!r} must be positive"
        )
    return number


def _optional_float(value: Any, *, process_id: str, field_name: str) -> float | None:
    if value is None:
        return None
    try:
        return _required_float(value, process_id=process_id, field_name=field_name)
    except AnchorError:
        return None


def _pdg_subblock_for_anchor(
    anchor: Anchor,
    *,
    process_id: str,
) -> Mapping[str, Any]:
    data = load_full_yaml(process_id, family=_FAMILY)
    pdg_block = data.get("pdg_or_equivalent")
    if not isinstance(pdg_block, Mapping):
        raise AnchorError(f"{process_id}: expected pdg_or_equivalent mapping")
    sub = pdg_block.get(anchor.block_key)
    if not isinstance(sub, Mapping):
        present = sorted(str(key) for key in pdg_block)
        raise AnchorError(
            f"{process_id}: load_anchor selected {anchor.block_key!r}, but that "
            f"block is not available as a mapping (present keys: {present})"
        )
    return sub


def _load_branching_subblock_anchor(
    process_id: str,
    block_key: str,
) -> BranchingFractionAnchor:
    scaffold_anchor = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=(block_key,),
        uncertainty_key=_SCAFFOLD_UNCERTAINTY_KEY,
    )
    sub = _pdg_subblock_for_anchor(scaffold_anchor, process_id=process_id)
    if scaffold_anchor.units != _EXPECTED_BRANCHING_UNITS:
        raise AnchorError(
            f"{process_id}: expected branching-fraction units for {block_key}, "
            f"got {scaffold_anchor.units!r}"
        )

    if block_key == _PARTIAL_RATE_BLOCK:
        stat_up = _positive_float(
            sub.get("uncertainty_stat_positive"),
            process_id=process_id,
            field_name=f"{block_key}.uncertainty_stat_positive",
        )
        stat_down = _positive_float(
            sub.get("uncertainty_stat_negative"),
            process_id=process_id,
            field_name=f"{block_key}.uncertainty_stat_negative",
        )
        syst = _positive_float(
            sub.get("uncertainty_syst"),
            process_id=process_id,
            field_name=f"{block_key}.uncertainty_syst",
        )
        uncertainty_upper = math.sqrt(stat_up * stat_up + syst * syst)
        uncertainty_lower = math.sqrt(stat_down * stat_down + syst * syst)
    else:
        uncertainty_upper = _positive_float(
            sub.get("uncertainty_positive"),
            process_id=process_id,
            field_name=f"{block_key}.uncertainty_positive",
        )
        uncertainty_lower = _positive_float(
            sub.get("uncertainty_negative"),
            process_id=process_id,
            field_name=f"{block_key}.uncertainty_negative",
        )

    return BranchingFractionAnchor(
        block_key=scaffold_anchor.block_key,
        source=_optional_str(scaffold_anchor.source),
        source_key=_optional_str(sub.get("source_key")),
        year=scaffold_anchor.year,
        value=float(scaffold_anchor.value),
        uncertainty_upper=float(uncertainty_upper),
        uncertainty_lower=float(uncertainty_lower),
        units=scaffold_anchor.units,
        source_url=_optional_str(scaffold_anchor.source_url),
        snapshot_path=_optional_str(scaffold_anchor.snapshot_path),
        assumptions=_optional_str(sub.get("assumptions")),
    )


def _values_entries(data: Mapping[str, Any], *, process_id: str) -> list[Mapping[str, Any]]:
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


def _value_entry_for_observable(
    process_id: str,
    observable: str,
) -> tuple[int, Mapping[str, Any]]:
    entries = _values_entries(load_full_yaml(process_id, family=_FAMILY), process_id=process_id)
    for index, entry in enumerate(entries):
        if entry.get("observable") == observable:
            return index, entry
    raise AnchorError(f"{process_id}: no pdg_or_equivalent.values entry {observable!r}")


def _load_virtual_anchor(
    process_id: str,
    *,
    block_key: str,
    entry: Mapping[str, Any],
) -> tuple[Anchor, str | None]:
    raw = _optional_str(entry.get("value"))
    value = _required_float(entry.get("value"), process_id=process_id, field_name=f"{block_key}.value")
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
            uncertainty_key=_SCAFFOLD_UNCERTAINTY_KEY,
        )
    finally:
        anchor_scaffold.load_pdg_block = original_load_pdg_block

    return scaffold_anchor, raw


def _load_value_anchor(
    process_id: str,
    observable: str,
    *,
    expected_units: str | None = None,
) -> CatalogValueAnchor:
    index, entry = _value_entry_for_observable(process_id, observable)
    block_key = f"pdg_or_equivalent.values[{index}]"
    scaffold_anchor, raw = _load_virtual_anchor(
        process_id,
        block_key=block_key,
        entry=entry,
    )
    if scaffold_anchor.observable != observable:
        raise AnchorError(
            f"{process_id}: selected {block_key} observable "
            f"{scaffold_anchor.observable!r}, expected {observable!r}"
        )
    if expected_units is not None and scaffold_anchor.units != expected_units:
        raise AnchorError(
            f"{process_id}: expected units {expected_units!r} for {observable!r}, "
            f"got {scaffold_anchor.units!r}"
        )
    return CatalogValueAnchor(
        section_key=_VALUES_SECTION,
        entry_index=index,
        block_key=scaffold_anchor.block_key,
        observable=_optional_str(scaffold_anchor.observable),
        source=_optional_str(scaffold_anchor.source),
        source_key=_optional_str(entry.get("source_key")),
        year=scaffold_anchor.year,
        value=float(scaffold_anchor.value),
        uncertainty=_optional_float(
            entry.get("uncertainty"),
            process_id=process_id,
            field_name=f"{scaffold_anchor.block_key}.uncertainty",
        ),
        units=scaffold_anchor.units,
        source_url=_optional_str(scaffold_anchor.source_url),
        snapshot_path=_optional_str(scaffold_anchor.snapshot_path),
        value_raw=raw,
    )


def _build_budget_band(full_rate: BranchingFractionAnchor) -> K010BudgetBand:
    lower_edge = float(full_rate.value - full_rate.uncertainty_lower)
    upper_edge = float(full_rate.value + full_rate.uncertainty_upper)
    hard_budget = max(full_rate.uncertainty_upper, full_rate.uncertainty_lower)
    if lower_edge <= 0.0 or hard_budget <= 0.0:
        raise AnchorError("K010: constructed full-rate budget is invalid")
    return K010BudgetBand(
        source=_BUDGET_SOURCE,
        construction=(
            "Observed full-region BR(K_S -> pi0 e+ e-) from K010.yaml; HARD "
            "ratio is |prediction - central| divided by the direction-aware "
            "asymmetric experimental uncertainty. No SM subtraction is applied."
        ),
        hard_veto_budget=float(hard_budget),
        experimental_sigma_upper=float(full_rate.uncertainty_upper),
        experimental_sigma_lower=float(full_rate.uncertainty_lower),
        lower_edge=lower_edge,
        upper_edge=upper_edge,
        sm_subtracted=False,
    )


def _a_s_inputs_from_full_rate(full_rate: BranchingFractionAnchor) -> KShortPi0EEChPTInputs:
    a_s_abs = math.sqrt(float(full_rate.value) / KSHORT_PI0EE_A_S_BRANCHING_COEFFICIENT)
    return KShortPi0EEChPTInputs(
        rate_coefficient=KSHORT_PI0EE_A_S_BRANCHING_COEFFICIENT,
        a_s_abs=float(a_s_abs),
    )


def _load_k010_anchor(process_id: str) -> K010Anchor:
    partial_rate = _load_branching_subblock_anchor(process_id, _PARTIAL_RATE_BLOCK)
    full_rate = _load_branching_subblock_anchor(process_id, _FULL_RATE_BLOCK)
    return K010Anchor(
        partial_rate=partial_rate,
        full_rate=full_rate,
        observed_candidates=_load_value_anchor(
            process_id,
            _OBSERVED_CANDIDATES_OBS,
            expected_units="events",
        ),
        expected_background=_load_value_anchor(
            process_id,
            _EXPECTED_BACKGROUND_OBS,
            expected_units="events",
        ),
        a_s_inputs=_a_s_inputs_from_full_rate(full_rate),
        budget_band=_build_budget_band(full_rate),
    )


def _selected_budget(
    predicted: float,
    anchor: K010Anchor,
) -> tuple[float, float, bool]:
    pull = float(predicted - anchor.value)
    budget = (
        anchor.budget_band.experimental_sigma_upper
        if pull >= 0.0
        else anchor.budget_band.experimental_sigma_lower
    )
    ratio = abs(pull) / budget if budget > 0.0 else float("inf")
    return float(budget), float(ratio), bool(ratio <= 1.0)


def _a_s_sign_branch_verdicts(
    result: Any,
    anchor: K010Anchor,
) -> tuple[
    dict[str, float | str | bool],
    dict[str, dict[str, float | str | bool]],
]:
    """Return (+|a_S|, -|a_S|) verdicts and the least-excluded branch."""

    branches: list[dict[str, float | str | bool]] = []
    for label, sign, predicted, a_s_effective in (
        (
            "positive",
            1.0,
            float(result.branching_fraction_positive_a_s),
            float(result.a_s_effective_positive),
        ),
        (
            "negative",
            -1.0,
            float(result.branching_fraction_negative_a_s),
            float(result.a_s_effective_negative),
        ),
    ):
        budget, ratio, passes = _selected_budget(predicted, anchor)
        branches.append(
            {
                "branch": label,
                "sign": float(sign),
                "predicted": float(predicted),
                "a_s_effective": float(a_s_effective),
                "budget": float(budget),
                "ratio": float(ratio),
                "passes": bool(passes),
            }
        )
    selected = min(branches, key=lambda item: float(item["ratio"]))
    return selected, {str(item["branch"]): item for item in branches}


@register
class Constraint:
    """Catalogued K_S -> pi0 e+ e- full-rate constraint (K010)."""

    process_id = "K010"
    severity = Severity.HARD
    observable = "BR(K_S -> pi0 e+ e-)"

    def __init__(self) -> None:
        self.anchor = _load_k010_anchor(self.process_id)
        self.chpt_inputs = self.anchor.a_s_inputs
        self.sm_inputs = rare_kaon_dilepton_default_sm_inputs()
        self.sm_result = kshort_pi0ee_a_s_sm(
            self.chpt_inputs,
            inputs=self.sm_inputs,
        )

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        rs_wilsons = point.get_extra(_REQUIRED_EXTRA)
        if rs_wilsons is None:
            return ConstraintResult(
                process_id=self.process_id,
                severity=self.severity,
                passes=True,
                sm_prediction=float(self.sm_result.sm_branching_fraction),
                experimental=float(self.anchor.value),
                budget=float(self.anchor.budget),
                notes=(
                    f"extra {_REQUIRED_EXTRA!r} absent; K_S -> pi0 e+ e- "
                    "constraint was not evaluated."
                ),
                diagnostics={
                    "evaluated": False,
                    "missing_extra": _REQUIRED_EXTRA,
                    "legacy_quark_mass_basis_couplings_present": (
                        point.get_extra("quark_mass_basis_couplings") is not None
                    ),
                    "partial_rate_branching_fraction": float(
                        self.anchor.partial_rate.value
                    ),
                    "full_rate_branching_fraction": float(self.anchor.full_rate.value),
                    "budget_source": self.anchor.budget_band.source,
                    "needs_human_physics": (
                        RS_SEMILEPTONIC_MATCHING_ASSUMPTION_V1
                        + " "
                        + RARE_KAON_KS_PI0EE_RS_MATCHING_ASSUMPTION_V1
                    ),
                    "semileptonic_qcd_running_applied": False,
                    "semileptonic_qcd_running_multiplicative_factor": 1.0,
                    "semileptonic_qcd_running_effect_fraction": 0.0,
                    "semileptonic_qcd_running_diagnostic": (
                        RARE_KAON_PI0EE_SEMILEPTONIC_RUNNING_DIAGNOSTIC_V1
                    ),
                },
            )

        kk_ew_mass = point.get_extra(_OPTIONAL_EW_MASS_EXTRA)
        try:
            result = kshort_pi0ee_a_s_from_rs_semileptonic_wilsons(
                rs_wilsons,
                chpt_inputs=self.chpt_inputs,
                lepton="e",
                matching_scale_gev=None if kk_ew_mass is None else float(kk_ew_mass),
                inputs=self.sm_inputs,
            )
        except Exception as exc:  # noqa: BLE001 - constraints degrade cleanly
            return ConstraintResult(
                process_id=self.process_id,
                severity=self.severity,
                passes=True,
                sm_prediction=float(self.sm_result.sm_branching_fraction),
                experimental=float(self.anchor.value),
                budget=float(self.anchor.budget),
                notes=(
                    "NOT EVALUATED - invalid rs_semileptonic_wilsons for "
                    "K_S -> pi0 e+ e-"
                ),
                diagnostics={
                    "evaluated": False,
                    "invalid_extra": _REQUIRED_EXTRA,
                    "exception_type": type(exc).__name__,
                    "exception_message": str(exc),
                    "budget_source": self.anchor.budget_band.source,
                },
            )
        selected_branch, sign_branches = _a_s_sign_branch_verdicts(
            result,
            self.anchor,
        )
        predicted = float(selected_branch["predicted"])
        budget = float(selected_branch["budget"])
        ratio = float(selected_branch["ratio"])
        passes = bool(selected_branch["passes"])
        positive_branch = sign_branches["positive"]
        negative_branch = sign_branches["negative"]

        diagnostics = dict(result.diagnostics)
        diagnostics.update(
            {
                "partial_rate_branching_fraction": float(
                    self.anchor.partial_rate.value
                ),
                "partial_rate_sigma_upper": float(
                    self.anchor.partial_rate.uncertainty_upper
                ),
                "partial_rate_sigma_lower": float(
                    self.anchor.partial_rate.uncertainty_lower
                ),
                "full_rate_branching_fraction": float(self.anchor.full_rate.value),
                "full_rate_sigma_upper": float(
                    self.anchor.full_rate.uncertainty_upper
                ),
                "full_rate_sigma_lower": float(
                    self.anchor.full_rate.uncertainty_lower
                ),
                "full_rate_assumptions": self.anchor.full_rate.assumptions,
                "observed_candidates": float(self.anchor.observed_candidates.value),
                "expected_background": float(self.anchor.expected_background.value),
                "budget": float(budget),
                "budget_source": self.anchor.budget_band.source,
                "budget_construction": self.anchor.budget_band.construction,
                "budget_lower_edge": float(self.anchor.budget_band.lower_edge),
                "budget_upper_edge": float(self.anchor.budget_band.upper_edge),
                "budget_sm_subtracted": self.anchor.budget_band.sm_subtracted,
                "parametrization_citation": RARE_KAON_KS_PI0EE_PARAMETRIZATION_CITATION,
                "rs_matching_assumption": RS_SEMILEPTONIC_MATCHING_ASSUMPTION_V1,
                "needs_human_physics": (
                    RS_SEMILEPTONIC_MATCHING_ASSUMPTION_V1
                    + " "
                    + RARE_KAON_KS_PI0EE_RS_MATCHING_ASSUMPTION_V1
                ),
                "rs_semileptonic_vector_matching_status": (
                    RARE_KAON_RS_SEMILEPTONIC_VECTOR_MATCHING_STATUS_V1
                ),
                "kk_ew_mass_extra_used": kk_ew_mass is not None,
                "sm_branching_fraction": float(result.sm_branching_fraction),
                "np_shift_branching_fraction": float(
                    predicted - result.sm_branching_fraction
                ),
                "a_s_sign_envelope_used_for_hard_verdict": True,
                "a_s_sign_envelope_selection_rule": (
                    "select the sign branch with the smaller direction-aware "
                    "HARD pull ratio"
                ),
                "selected_a_s_branch": str(selected_branch["branch"]),
                "selected_a_s_sign": float(selected_branch["sign"]),
                "a_s_effective": float(selected_branch["a_s_effective"]),
                "selected_a_s_effective": float(selected_branch["a_s_effective"]),
                "positive_a_s_effective": float(positive_branch["a_s_effective"]),
                "negative_a_s_effective": float(negative_branch["a_s_effective"]),
                "positive_a_s_branching_fraction": float(
                    positive_branch["predicted"]
                ),
                "negative_a_s_branching_fraction": float(
                    negative_branch["predicted"]
                ),
                "positive_a_s_np_shift_branching_fraction": float(
                    positive_branch["predicted"] - result.sm_branching_fraction
                ),
                "negative_a_s_np_shift_branching_fraction": float(
                    negative_branch["predicted"] - result.sm_branching_fraction
                ),
                "positive_a_s_budget": float(positive_branch["budget"]),
                "negative_a_s_budget": float(negative_branch["budget"]),
                "positive_a_s_ratio": float(positive_branch["ratio"]),
                "negative_a_s_ratio": float(negative_branch["ratio"]),
                "positive_a_s_passes": bool(positive_branch["passes"]),
                "negative_a_s_passes": bool(negative_branch["passes"]),
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
                "BR(K_S -> pi0 e+ e-) uses the a_S-driven full-region "
                "Buchalla-D'Ambrosio-Isidori / ISU parametrization. The HARD "
                "budget is the K010.yaml PDG/NA48 full-rate uncertainty; the "
                "direct partial-rate measurement is diagnostic. The RS shift "
                "is evaluated for both a_S signs and the HARD verdict uses "
                "the conservative least-excluded sign branch. Phase-3a RS "
                "C9 maps additively into the y7V a_S shift; this is still not "
                "a complete long-distance K_S matching."
            ),
            diagnostics=diagnostics,
        )
