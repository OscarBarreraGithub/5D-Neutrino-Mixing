"""B023 - invisible vector rare beauty decay ``B -> K* nu nubar``.

Physics
-------
The short-distance response is reused from the shared ``b -> s nu nubar``
adapter built for B022.  The core supplies the Wilson response and the vector
coefficient

    R_K* = (1 + 1.31 eta) epsilon^2,

where ``epsilon`` and ``eta`` are the standard left/right-current response
parameters.  The exclusive SM normalization is not hardcoded here: it is
loaded from ``B023.yaml`` and applied as

    BR(B -> K* nu nubar) = BR_SM(B -> K* nu nubar) R_K*.

The RS contribution is the same Phase-4d ``b_to_s_nunu`` Wilson block as B022.

Severity
--------
HARD.  The predicted total branching fraction is compared with the Belle 2017
combined vector 90% CL upper limit recorded in ``B023.yaml``.  Following the
upper-limit policy used by K005, the HARD ratio is ``BR_total / BR_limit``;
the SM-subtracted room is retained only as a diagnostic.

Catalog sidecar
---------------
``flavor_catalog/processes/beauty/B023.yaml`` is the source of truth for the
Belle combined vector limit and Buras et al. 2015 SM reference.  The sidecar
does not contain a Belle II vector-mode limit; this constraint preserves the
catalogued Belle 2017 provenance.

Phase-4d Rewire
---------------
The implemented NP term is ``rs_semileptonic_wilsons.b_to_s_nunu`` mapped
directly as ``X_NP=C/g_SM^2``; the old one-Z-like proxy is not used.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
import re
from typing import Any, Mapping

from flavor_catalog_constraints import anchors as anchor_scaffold
from flavor_catalog_constraints.anchors import Anchor, AnchorError, load_full_yaml
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.rare_b_nunu import (
    RARE_B_NUNU_KSTAR_ETA_COEFFICIENT,
    RARE_B_NUNU_RS_SEMILEPTONIC_MATCHING_STATUS_V1,
    RS_SEMILEPTONIC_MATCHING_ASSUMPTION_V1,
    b_to_kstar_nunu_branching_fraction,
    b_to_kstar_nunu_from_rs_semileptonic_wilsons,
    rare_b_nunu_default_sm_inputs,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "beauty"
_REQUIRED_EXTRA = "rs_semileptonic_wilsons"
_OPTIONAL_EW_MASS_EXTRA = "kk_ew_mass_gev"
_EXPECTED_UNITS = "branching fraction"
_BELLE_COMBINED_VALUE_ID = "Belle2017:B023:combined_vector_limit"
_BURAS_SM_VALUE_ID = "Buras2015:B023:sm_b0_kstar0_nunu"
_BUDGET_SOURCE = (
    "flavor_catalog/processes/beauty/B023.yaml "
    "pdg_or_equivalent.values[Belle2017:B023:combined_vector_limit] "
    "+ pdg_or_equivalent.values[Buras2015:B023:sm_b0_kstar0_nunu]"
)
_PARAMETRIZATION_CITATION = (
    "Buras-Girrbach-Noe-Niehoff-Straub arXiv:1409.4557 "
    "B -> K* nu nubar SM normalization and epsilon/eta/r_kstar response; "
    "shared Inami-Lim X_t b -> s nu nubar core"
)
_TREE_LEVEL_STATUS = (
    "Phase-4d rigorous light-Z active-neutrino Wilson from "
    "rs_semileptonic_wilsons.b_to_s_nunu; old one-Z-like proxy resolved."
)
_UNEVALUATED_REASON = "rs_semileptonic_wilsons.b_to_s_nunu not provided"
_UPPER_LIMIT_RE = re.compile(r"^\s*<\s*(?P<value>[0-9.eE+-]+)\s*$")


@dataclass(frozen=True)
class BranchingFractionAnchor:
    """Branching-fraction value or upper-limit entry from B023.yaml."""

    block_key: str
    entry_index: int
    entry_id: str
    observable: str | None
    source_key: str | None
    year: int | None
    value: float
    uncertainty: float | None
    units: str | None
    is_upper_limit: bool
    confidence_level: str | None
    source_url: str | None
    snapshot_path: str | None
    value_raw: str | None


@dataclass(frozen=True)
class BKstarNuNuBudgetBand:
    """B023 upper-limit budget for the HARD veto."""

    source: str
    construction: str
    hard_veto_budget: float
    confidence_level: str | None
    sm_theory_sigma: float | None
    limit_minus_sm_anchor: float
    sm_subtracted: bool


@dataclass(frozen=True)
class B023Anchor:
    """Typed B023 anchor: Belle limit, SM reference, and budget."""

    experimental: BranchingFractionAnchor
    standard_model: BranchingFractionAnchor
    budget_band: BKstarNuNuBudgetBand

    @property
    def value(self) -> float:
        return self.experimental.value

    @property
    def uncertainty(self) -> float | None:
        return None

    @property
    def sm_value(self) -> float:
        return self.standard_model.value

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
            f"{process_id}: B023 anchor field {field_name!r}={value!r} "
            "is not numeric"
        ) from exc
    if not math.isfinite(number):
        raise AnchorError(
            f"{process_id}: B023 anchor field {field_name!r}={value!r} "
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


def _pdg_block(process_id: str) -> Mapping[str, Any]:
    data = load_full_yaml(process_id, family=_FAMILY)
    block = data.get("pdg_or_equivalent")
    if not isinstance(block, Mapping):
        raise AnchorError(
            f"{process_id}: expected mapping-shaped 'pdg_or_equivalent' for B023, "
            f"got {type(block).__name__}"
        )
    return block


def _values_entries(
    pdg_block: Mapping[str, Any],
    *,
    process_id: str,
) -> tuple[Mapping[str, Any], ...]:
    entries = pdg_block.get("values")
    if not isinstance(entries, list) or not entries:
        raise AnchorError(f"{process_id}: pdg_or_equivalent.values must be a list")
    out: list[Mapping[str, Any]] = []
    for index, entry in enumerate(entries):
        if not isinstance(entry, Mapping):
            raise AnchorError(
                f"{process_id}: pdg_or_equivalent.values[{index}] is not a mapping"
            )
        out.append(entry)
    return tuple(out)


def _find_value_entry(
    entries: tuple[Mapping[str, Any], ...],
    *,
    value_id: str,
    process_id: str,
) -> tuple[int, Mapping[str, Any]]:
    for index, entry in enumerate(entries):
        if str(entry.get("value_id")) == value_id:
            return index, entry
    raise AnchorError(
        f"{process_id}: no pdg_or_equivalent.values entry with "
        f"value_id={value_id!r}"
    )


def _load_scaffold_value_anchor(
    value_id: str,
    *,
    process_id: str,
    expect_upper_limit: bool,
) -> tuple[Anchor, Mapping[str, Any], bool, str | None]:
    pdg_block = _pdg_block(process_id)
    index, entry = _find_value_entry(
        _values_entries(pdg_block, process_id=process_id),
        value_id=value_id,
        process_id=process_id,
    )
    block_key = f"pdg_or_equivalent.values[{index}]"
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
        scaffold_anchor = anchor_scaffold.load_anchor(
            process_id,
            family=_FAMILY,
            candidates=(block_key,),
        )
    finally:
        anchor_scaffold.load_pdg_block = original_load_pdg_block

    if scaffold_anchor.block_key != block_key:
        raise AnchorError(
            f"{process_id}: load_anchor selected {scaffold_anchor.block_key!r}, "
            f"expected {block_key!r} for B023 value entry {value_id!r}"
        )
    return scaffold_anchor, entry, is_upper_limit, raw


def _load_branching_anchor(
    value_id: str,
    *,
    process_id: str,
    expect_upper_limit: bool = False,
) -> BranchingFractionAnchor:
    scaffold_anchor, entry, is_upper_limit, raw = _load_scaffold_value_anchor(
        value_id,
        process_id=process_id,
        expect_upper_limit=expect_upper_limit,
    )
    if scaffold_anchor.units != _EXPECTED_UNITS:
        raise AnchorError(
            f"{process_id}: expected units {_EXPECTED_UNITS!r} in "
            f"{scaffold_anchor.block_key}, got {scaffold_anchor.units!r}"
        )
    value = _required_float(
        scaffold_anchor.value,
        process_id=process_id,
        field_name=f"{scaffold_anchor.block_key}.value",
    )
    uncertainty = (
        None
        if scaffold_anchor.uncertainty is None
        else _required_float(
            scaffold_anchor.uncertainty,
            process_id=process_id,
            field_name=f"{scaffold_anchor.block_key}.uncertainty",
        )
    )
    if value <= 0.0:
        raise AnchorError(f"{process_id}: B023 branching value must be positive")
    if uncertainty is not None and uncertainty <= 0.0:
        raise AnchorError(f"{process_id}: B023 uncertainty must be positive")
    if expect_upper_limit and not is_upper_limit:
        raise AnchorError(f"{process_id}: {value_id} must be an upper limit")
    return BranchingFractionAnchor(
        block_key=scaffold_anchor.block_key,
        entry_index=int(scaffold_anchor.block_key.rsplit("[", 1)[1].rstrip("]")),
        entry_id=value_id,
        observable=_optional_str(scaffold_anchor.observable),
        source_key=_optional_str(entry.get("source_key")),
        year=scaffold_anchor.year,
        value=float(value),
        uncertainty=None if uncertainty is None else float(uncertainty),
        units=scaffold_anchor.units,
        is_upper_limit=is_upper_limit,
        confidence_level=_optional_str(entry.get("confidence_level")),
        source_url=_optional_str(scaffold_anchor.source_url),
        snapshot_path=_optional_str(scaffold_anchor.snapshot_path),
        value_raw=raw,
    )


def _build_budget_band(
    *,
    experimental: BranchingFractionAnchor,
    standard_model: BranchingFractionAnchor,
) -> BKstarNuNuBudgetBand:
    if not experimental.is_upper_limit:
        raise AnchorError("B023: experimental anchor must be an upper limit")
    limit_minus_sm = experimental.value - standard_model.value
    if experimental.value <= 0.0 or limit_minus_sm <= 0.0:
        raise AnchorError("B023: constructed branching-ratio budget is invalid")
    return BKstarNuNuBudgetBand(
        source=_BUDGET_SOURCE,
        construction=(
            "Published Belle 2017 90% CL combined upper limit on the total "
            "B -> K* nu nubar branching fraction; ratio = BR_total / BR_limit "
            "with no central-residual subtraction."
        ),
        hard_veto_budget=float(experimental.value),
        confidence_level=experimental.confidence_level,
        sm_theory_sigma=standard_model.uncertainty,
        limit_minus_sm_anchor=float(limit_minus_sm),
        sm_subtracted=False,
    )


def _load_b023_anchor(process_id: str) -> B023Anchor:
    experimental = _load_branching_anchor(
        _BELLE_COMBINED_VALUE_ID,
        process_id=process_id,
        expect_upper_limit=True,
    )
    standard_model = _load_branching_anchor(
        _BURAS_SM_VALUE_ID,
        process_id=process_id,
    )
    return B023Anchor(
        experimental=experimental,
        standard_model=standard_model,
        budget_band=_build_budget_band(
            experimental=experimental,
            standard_model=standard_model,
        ),
    )


@register
class Constraint:
    """Catalogued ``B -> K* nu nubar`` branching-ratio constraint (B023)."""

    process_id = "B023"
    severity = Severity.HARD
    observable = "BR(B -> K* nu nubar)"

    def __init__(self) -> None:
        self.anchor = _load_b023_anchor(self.process_id)
        self.sm_inputs = rare_b_nunu_default_sm_inputs()
        self.sm_result = b_to_kstar_nunu_branching_fraction(
            sm_branching_fraction=self.anchor.sm_value,
            inputs=self.sm_inputs,
        )

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        couplings = point.get_extra(_REQUIRED_EXTRA)
        if couplings is None:
            return ConstraintResult(
                process_id=self.process_id,
                severity=self.severity,
                passes=True,
                predicted=None,
                sm_prediction=float(self.sm_result.branching_fraction),
                experimental=float(self.anchor.value),
                ratio=None,
                budget=float(self.anchor.budget),
                notes=(
                    "NOT EVALUATED - rs_semileptonic_wilsons.b_to_s_nunu "
                    "absent; B -> K* nu nubar constraint is non-vetoing."
                ),
                diagnostics={
                    "evaluated": False,
                    "unevaluated_reason": _UNEVALUATED_REASON,
                    "missing_extra": _REQUIRED_EXTRA,
                    "legacy_quark_mass_basis_couplings_present": (
                        point.get_extra("quark_mass_basis_couplings") is not None
                    ),
                    "sm_anchor_branching_fraction": float(self.anchor.sm_value),
                    "sm_formula_branching_fraction": float(
                        self.sm_result.branching_fraction
                    ),
                    "budget_source": self.anchor.budget_band.source,
                    "budget_construction": self.anchor.budget_band.construction,
                    "tree_level_status": _TREE_LEVEL_STATUS,
                },
            )

        kk_ew_mass = point.get_extra(_OPTIONAL_EW_MASS_EXTRA)
        try:
            result = b_to_kstar_nunu_from_rs_semileptonic_wilsons(
                couplings,
                sm_branching_fraction=self.anchor.sm_value,
                matching_scale_gev=None if kk_ew_mass is None else float(kk_ew_mass),
                inputs=self.sm_inputs,
            )
        except (AttributeError, KeyError, TypeError, ValueError) as exc:
            return ConstraintResult(
                process_id=self.process_id,
                severity=self.severity,
                passes=True,
                predicted=None,
                sm_prediction=float(self.sm_result.branching_fraction),
                experimental=float(self.anchor.value),
                ratio=None,
                budget=float(self.anchor.budget),
                notes=(
                    "NOT EVALUATED - invalid rs_semileptonic_wilsons for "
                    "B -> K* nu nubar"
                ),
                diagnostics={
                    "evaluated": False,
                    "unevaluated_reason": _UNEVALUATED_REASON,
                    "invalid_extra": _REQUIRED_EXTRA,
                    "exception_type": type(exc).__name__,
                    "exception_message": str(exc),
                    "tree_level_status": _TREE_LEVEL_STATUS,
                },
            )
        predicted = float(result.branching_fraction)
        budget = float(self.anchor.budget)
        ratio = predicted / budget if budget > 0.0 else float("inf")
        passes = bool(ratio <= 1.0)

        diagnostics = dict(result.diagnostics)
        diagnostics.update(
            {
                "evaluated": True,
                "sm_anchor_branching_fraction": float(self.anchor.sm_value),
                "sm_formula_branching_fraction": float(
                    result.sm_branching_fraction
                ),
                "sm_formula_minus_anchor": float(
                    result.sm_branching_fraction - self.anchor.sm_value
                ),
                "experimental_block": self.anchor.experimental.block_key,
                "experimental_entry_id": self.anchor.experimental.entry_id,
                "sm_block": self.anchor.standard_model.block_key,
                "sm_entry_id": self.anchor.standard_model.entry_id,
                "experimental_confidence_level": (
                    self.anchor.experimental.confidence_level
                ),
                "budget": float(budget),
                "budget_source": self.anchor.budget_band.source,
                "budget_construction": self.anchor.budget_band.construction,
                "budget_confidence_level": self.anchor.budget_band.confidence_level,
                "budget_limit_minus_sm_anchor": float(
                    self.anchor.budget_band.limit_minus_sm_anchor
                ),
                "budget_sm_theory_sigma": (
                    self.anchor.budget_band.sm_theory_sigma
                ),
                "budget_sm_subtracted": self.anchor.budget_band.sm_subtracted,
                "parametrization_citation": _PARAMETRIZATION_CITATION,
                "rs_matching_assumption": RS_SEMILEPTONIC_MATCHING_ASSUMPTION_V1,
                "tree_level_status": _TREE_LEVEL_STATUS,
                "rs_semileptonic_nunu_matching_status": (
                    RARE_B_NUNU_RS_SEMILEPTONIC_MATCHING_STATUS_V1
                ),
                "kk_ew_mass_extra_used": kk_ew_mass is not None,
                "lambda_wolfenstein": float(result.lambda_wolfenstein),
                "lambda_t_bs": complex(result.lambda_t_bs),
                "c_l_sm": complex(result.c_l_sm),
                "c_l_total": complex(result.c_l_total),
                "c_r_total": complex(result.c_r_total),
                "x_eff_left": complex(result.x_eff_left),
                "x_eff_right": complex(result.x_eff_right),
                "x_t": float(result.x_t),
                "epsilon": float(result.epsilon),
                "eta": float(result.eta),
                "r_k": float(result.r_k),
                "r_kstar": float(result.r_kstar),
                "kstar_eta_coefficient": float(RARE_B_NUNU_KSTAR_ETA_COEFFICIENT),
                "np_shift_branching_fraction": float(
                    result.np_shift_branching_fraction
                ),
                "down_sector_indices": (1, 2),
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
                "BR(B -> K* nu nubar) uses the shared b -> s nu nubar "
                "epsilon/eta response with the vector coefficient "
                "r_kstar=(1+1.31 eta) epsilon^2, normalized to the Buras "
                "2015 SM branching fraction from B023.yaml. RS contribution "
                "is the Phase-4d b_to_s_nunu Wilson block mapped as "
                "X_NP=C/g_SM^2; the old one-Z-like proxy is not used. "
                "Budget is the Belle 2017 90% CL combined upper limit."
            ),
            diagnostics=diagnostics,
        )
