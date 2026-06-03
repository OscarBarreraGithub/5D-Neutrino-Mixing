"""K005 - neutral rare kaon decay ``K_L -> pi0 nu nubar``.

Physics
-------
The SM short-distance prediction is the purely CP-violating sibling of K004,

    BR(K_L -> pi0 nu nubar) =
        kappa_L [ Im(lambda_t X_t) / lambda^5 ]^2.

The Phase-4d rewire reuses the K004 ``s -> d nu nubar`` Wilson machinery,
adding the same complex ``X_NP`` shift to ``lambda_t X_t`` but retaining only the
imaginary part,

    Im(lambda_t X_t + X_NP).

The low-level formula and the documented RS matching assumption live in
``quarkConstraints.rare_kaon_snd`` and are reached only through the
``flavor_catalog_constraints.physics_adapters.rare_kaon`` boundary.  K005's
extension of the K004 rare-kaon module is intentionally append-only and is the
designed location for this physics, not an isolation violation.

Severity
--------
HARD.  The predicted total branching fraction is compared with the KOTO 90% CL
upper limit loaded from ``K005.yaml``.  This is an observed upper bound, so the
budget is the published limit itself rather than a central residual.

Catalog sidecar
---------------
``flavor_catalog/processes/kaon/K005.yaml`` is the source of truth for the
KOTO/PDG limit and SM reference values.  K005 stores these as list entries, so
this module adapts the documented list entries into the scaffold ``load_anchor``
path rather than reading value-bearing anchors ad hoc.
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
from flavor_catalog_constraints.physics_adapters.rare_kaon import (
    RARE_KAON_KAPPA_L_CITATION,
    RARE_KAON_NUNU_RS_SEMILEPTONIC_MATCHING_STATUS_V1,
    RS_SEMILEPTONIC_MATCHING_ASSUMPTION_V1,
    klong_pi0_nunu_from_rs_semileptonic_wilsons,
    rare_kaon_default_sm_inputs,
    rare_kaon_neutral_sm_branching_fraction,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "kaon"
_REQUIRED_EXTRA = "rs_semileptonic_wilsons"
_OPTIONAL_EW_MASS_EXTRA = "kk_ew_mass_gev"
_EXPECTED_UNITS = "dimensionless branching fraction"
_KOTO_LIMIT_BLOCK = "experimental_inputs[0]"
_PDG_LIMIT_BLOCK = "pdg_or_equivalent[0]"
_SM_BURAS_VENTURINI_BLOCK = "theory_inputs[1]"
_SM_BGS_VALIDATION_BLOCK = "theory_inputs[0]"
_BUDGET_SOURCE = (
    f"flavor_catalog/processes/kaon/K005.yaml {_KOTO_LIMIT_BLOCK} "
    "KOTO 2025 90% CL upper limit"
)
_PARAMETRIZATION_CITATION = (
    "Buras et al. JHEP 11 (2015) 033, arXiv:1503.02693; "
    "Brod-Gorbahn-Stamou arXiv:2105.02868; "
    "Buras-Venturini arXiv:2203.10099"
)
_TREE_LEVEL_STATUS = (
    "Phase-4d rigorous light-Z active-neutrino Wilson from "
    "rs_semileptonic_wilsons.s_to_d_nunu; old one-Z-like proxy resolved."
)
_UNEVALUATED_REASON = "rs_semileptonic_wilsons.s_to_d_nunu not provided"

_UPPER_LIMIT_RE = re.compile(
    r"^\s*<\s*(?P<value>[0-9.eE+-]+)\s*$"
)
_LIST_BLOCK_RE = re.compile(
    r"^(?P<section>[A-Za-z_][A-Za-z0-9_]*)\[(?P<index>[0-9]+)\]$"
)


@dataclass(frozen=True)
class BranchingFractionAnchor:
    """Branching-fraction value or upper-limit entry from K005.yaml."""

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
class NeutralRareKaonBudgetBand:
    """K005 upper-limit budget for the HARD veto."""

    source: str
    construction: str
    hard_veto_budget: float
    confidence_level: str | None
    sm_theory_sigma: float | None
    limit_minus_sm_anchor: float
    sm_subtracted: bool


@dataclass(frozen=True)
class K005Anchor:
    """Typed K005 anchor: KOTO limit, PDG mirror, and SM references."""

    experimental: BranchingFractionAnchor
    pdg_limit: BranchingFractionAnchor
    standard_model: BranchingFractionAnchor
    validation_standard_model: BranchingFractionAnchor
    budget_band: NeutralRareKaonBudgetBand

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
            f"{process_id}: K005 anchor field {field_name!r}={value!r} "
            "is not numeric"
        ) from exc
    if not math.isfinite(number):
        raise AnchorError(
            f"{process_id}: K005 anchor field {field_name!r}={value!r} "
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


def _parse_list_block_key(block_key: str, *, process_id: str) -> tuple[str, int]:
    match = _LIST_BLOCK_RE.match(block_key)
    if match is None:
        raise AnchorError(f"{process_id}: invalid K005 list anchor key {block_key!r}")
    return match.group("section"), int(match.group("index"))


def _list_entry_for_block(
    data: Mapping[str, Any],
    block_key: str,
    *,
    process_id: str,
) -> tuple[str, int, Mapping[str, Any]]:
    section_key, entry_index = _parse_list_block_key(block_key, process_id=process_id)
    entries = _section_entries(data, section_key, process_id=process_id)
    try:
        entry = entries[entry_index]
    except IndexError as exc:
        raise AnchorError(
            f"{process_id}: K005 list anchor {block_key!r} is out of range"
        ) from exc
    return section_key, entry_index, entry


def _load_scaffold_list_anchor(
    block_key: str,
    *,
    process_id: str,
    expect_upper_limit: bool,
) -> tuple[Anchor, Mapping[str, Any], bool, str | None, str, int]:
    data = load_full_yaml(process_id, family=_FAMILY)
    section_key, entry_index, entry = _list_entry_for_block(
        data,
        block_key,
        process_id=process_id,
    )
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

    def _load_virtual_pdg_block(request_process_id: str, **kwargs: Any) -> Mapping[str, Any]:
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

    return scaffold_anchor, entry, is_upper_limit, raw, section_key, entry_index


def _load_branching_anchor(
    block_key: str,
    *,
    process_id: str,
    expect_upper_limit: bool = False,
) -> BranchingFractionAnchor:
    scaffold_anchor, entry, is_upper_limit, raw, section_key, index = (
        _load_scaffold_list_anchor(
            block_key,
            process_id=process_id,
            expect_upper_limit=expect_upper_limit,
        )
    )
    units = scaffold_anchor.units
    if units != _EXPECTED_UNITS:
        raise AnchorError(
            f"{process_id}: expected units {_EXPECTED_UNITS!r} in "
            f"{section_key}[{index}], got {units!r}"
        )
    value = _required_float(
        scaffold_anchor.value,
        process_id=process_id,
        field_name=f"{block_key}.value",
    )
    uncertainty = (
        None
        if scaffold_anchor.uncertainty is None
        else _required_float(
            scaffold_anchor.uncertainty,
            process_id=process_id,
            field_name=f"{block_key}.uncertainty",
        )
    )
    if value <= 0.0:
        raise AnchorError(f"{process_id}: K005 branching value must be positive")
    if uncertainty is not None and uncertainty <= 0.0:
        raise AnchorError(f"{process_id}: K005 uncertainty must be positive")
    return BranchingFractionAnchor(
        section_key=section_key,
        entry_index=index,
        block_key=block_key,
        source=_optional_str(scaffold_anchor.source),
        year=scaffold_anchor.year,
        value=float(value),
        uncertainty=None if uncertainty is None else float(uncertainty),
        is_upper_limit=is_upper_limit,
        confidence_level=_optional_str(entry.get("confidence_level")),
        units=units,
        source_url=_optional_str(scaffold_anchor.source_url),
        snapshot_path=_optional_str(scaffold_anchor.snapshot_path),
        value_raw=raw,
        observable=_optional_str(scaffold_anchor.observable),
    )


def _build_budget_band(
    *,
    experimental: BranchingFractionAnchor,
    standard_model: BranchingFractionAnchor,
) -> NeutralRareKaonBudgetBand:
    if not experimental.is_upper_limit:
        raise AnchorError("K005: experimental anchor must be an upper limit")
    limit_minus_sm = experimental.value - standard_model.value
    if experimental.value <= 0.0 or limit_minus_sm <= 0.0:
        raise AnchorError("K005: constructed branching-ratio budget is invalid")
    return NeutralRareKaonBudgetBand(
        source=_BUDGET_SOURCE,
        construction=(
            "Published KOTO 90% CL upper limit on the total branching fraction; "
            "ratio = BR_total / BR_limit with no central-residual subtraction."
        ),
        hard_veto_budget=float(experimental.value),
        confidence_level=experimental.confidence_level,
        sm_theory_sigma=standard_model.uncertainty,
        limit_minus_sm_anchor=float(limit_minus_sm),
        sm_subtracted=False,
    )


def _load_k005_anchor(process_id: str) -> K005Anchor:
    experimental = _load_branching_anchor(
        _KOTO_LIMIT_BLOCK,
        process_id=process_id,
        expect_upper_limit=True,
    )
    pdg_limit = _load_branching_anchor(
        _PDG_LIMIT_BLOCK,
        process_id=process_id,
        expect_upper_limit=True,
    )
    standard_model = _load_branching_anchor(
        _SM_BURAS_VENTURINI_BLOCK,
        process_id=process_id,
    )
    validation_standard_model = _load_branching_anchor(
        _SM_BGS_VALIDATION_BLOCK,
        process_id=process_id,
    )
    if not math.isclose(experimental.value, pdg_limit.value, rel_tol=0.0, abs_tol=1e-30):
        raise AnchorError(
            f"{process_id}: KOTO limit {experimental.value} does not match "
            f"PDG mirror {pdg_limit.value}"
        )
    return K005Anchor(
        experimental=experimental,
        pdg_limit=pdg_limit,
        standard_model=standard_model,
        validation_standard_model=validation_standard_model,
        budget_band=_build_budget_band(
            experimental=experimental,
            standard_model=standard_model,
        ),
    )


@register
class Constraint:
    """Catalogued neutral rare-kaon branching-ratio constraint (K005)."""

    process_id = "K005"
    severity = Severity.HARD
    observable = "BR(K_L -> pi0 nu nubar)"

    def __init__(self) -> None:
        self.anchor = _load_k005_anchor(self.process_id)
        self.sm_inputs = rare_kaon_default_sm_inputs()
        self.sm_result = rare_kaon_neutral_sm_branching_fraction(self.sm_inputs)

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
                    "NOT EVALUATED - rs_semileptonic_wilsons.s_to_d_nunu "
                    "absent; K_L -> pi0 nu nubar constraint is non-vetoing."
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
            result = klong_pi0_nunu_from_rs_semileptonic_wilsons(
                couplings,
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
                    "K_L -> pi0 nu nubar"
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
                "validation_sm_bgs2021_branching_fraction": float(
                    self.anchor.validation_standard_model.value
                ),
                "validation_sm_bgs2021_uncertainty": (
                    self.anchor.validation_standard_model.uncertainty
                ),
                "experimental_block": self.anchor.experimental.block_key,
                "pdg_limit_block": self.anchor.pdg_limit.block_key,
                "sm_block": self.anchor.standard_model.block_key,
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
                "kappa_l_citation": RARE_KAON_KAPPA_L_CITATION,
                "rs_matching_assumption": RS_SEMILEPTONIC_MATCHING_ASSUMPTION_V1,
                "tree_level_status": _TREE_LEVEL_STATUS,
                "rs_semileptonic_nunu_matching_status": (
                    RARE_KAON_NUNU_RS_SEMILEPTONIC_MATCHING_STATUS_V1
                ),
                "kk_ew_mass_extra_used": kk_ew_mass is not None,
                "qcd_running_applied": False,
                "qcd_running_note": (
                    "No Delta-F=2 four-quark matrix element is evaluated; "
                    "K_L -> pi0 nu nubar uses the clean semileptonic "
                    "kappa_L normalization and CP-odd short-distance "
                    "projection."
                ),
                "kappa_l": float(result.kappa_l),
                "x_t": float(result.x_t),
                "lambda_wolfenstein": float(result.lambda_wolfenstein),
                "lambda_t": complex(result.lambda_t),
                "x_eff_top": complex(result.x_eff_top),
                "x_np_total": complex(result.x_np_total),
                "im_x_eff_top": float(result.x_eff_top.imag),
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
                "BR(K_L -> pi0 nu nubar) uses the purely CP-violating "
                "Buras short-distance parametrization with Im(lambda_t X_t "
                "+ X_NP). RS contribution is the Phase-4d s_to_d_nunu Wilson "
                "block mapped as X_NP=C/g_SM^2 and restricted to the "
                "imaginary part; the old one-Z-like proxy is not used. "
                "Budget is the KOTO 90% CL upper limit from K005.yaml."
            ),
            diagnostics=diagnostics,
        )
