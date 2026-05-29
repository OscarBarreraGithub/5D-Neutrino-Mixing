"""B015 - inclusive rare decay ``B -> X_s ell+ ell-``.

Physics
-------
The active observable is the low-q2 inclusive bin in ``B015.yaml`` because the
sidecar supplies both the weighted experimental average and the Huber-Hurth-
Lunghi SM prediction for that bin:

    1 GeV^2 < q^2 < 6 GeV^2.

The prediction uses the inclusive extension of the shared ``b -> s l l`` core,
reached only through ``flavor_catalog_constraints.physics_adapters.rare_b_meson``.
It evaluates a leading partonic C7/C9/C10 shape,

    (1 - shat)^2 [(1 + 2 shat)(|C9|^2 + |C10|^2)
        + 4(1 + 2/shat)|C7|^2 + 12 Re(C7 C9*)] + primed terms,

normalized to the YAML SM bin branching fraction.  The same C9/C10 RS proxy
used by B016 is combined with the existing B011 C7 proxy.

Severity
--------
HARD.  The predicted low-q2 inclusive branching fraction is compared to the
YAML low-q2 experimental average with a B016-style proxy-theory budget:
``|BR_exp - BR_SM(formula)| + sqrt(sigma_exp^2 + sigma_SM^2) + sigma_proxy``.
The proxy term is the same documented 30% C9/C10 envelope used by B016.

Catalog sidecar
---------------
``flavor_catalog/processes/beauty/B015.yaml`` is the source of truth for the
HFLAV total inclusive average, low/high-q2 experimental averages, and low-q2
SM anchors.  Numeric values below are loaded from that sidecar, not hardcoded.

NEEDS-HUMAN-PHYSICS
-------------------
The RS contribution is a documented proxy.  ``ParameterPoint`` lacks the full
EW KK/Z/Z', lepton, dipole-loop, NNLO/QED, charm-veto, and covariance inputs
needed for a production inclusive ``b -> s l l`` likelihood.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
import re
from typing import Any, Mapping

from flavor_catalog_constraints import anchors as anchor_scaffold
from flavor_catalog_constraints.anchors import Anchor, AnchorError, load_full_yaml
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.rare_b_meson import (
    RARE_B_DILEPTON_INCLUSIVE_XS_LIMITATION_V1,
    RARE_B_DILEPTON_INCLUSIVE_XS_PROXY_THEORY_UNCERTAINTY_FRACTION,
    RARE_B_DILEPTON_INCLUSIVE_XS_PROXY_THEORY_UNCERTAINTY_RATIONALE,
    RARE_B_DILEPTON_RS_MATCHING_ASSUMPTION_V1,
    inclusive_b_to_xs_mumu_from_couplings,
    rare_b_inclusive_xs_dilepton_default_inputs,
    rare_b_inclusive_xs_mumu_sm_branching_fraction,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "beauty"
_REQUIRED_EXTRA = "quark_mass_basis_couplings"
_OPTIONAL_EW_MASS_EXTRA = "kk_ew_mass_gev"
_TOTAL_OBSERVABLE = "BR(B -> X_s ell+ ell-)"
_PDG_TOTAL_OBSERVABLE = "PDG-listed BR(B -> X_s ell+ ell-)"
_LOW_Q2_EXPERIMENTAL_OBSERVABLE = "Experimental low-q2 weighted average"
_HIGH_Q2_EXPERIMENTAL_OBSERVABLE = "Experimental high-q2 weighted average"
_LOW_Q2_SM_MUMU_OBSERVABLE = "SM BR[1,6]_mumu"
_LOW_Q2_SM_EE_OBSERVABLE = "SM BR[1,6]_ee"
_SCAFFOLD_NO_SCALAR_UNCERTAINTY_KEY = "__b015_uncertainty_parsed_below__"
_BUDGET_SOURCE = (
    "flavor_catalog/processes/beauty/B015.yaml "
    "Experimental low-q2 weighted average + SM BR[1,6]_mumu "
    "+ B016 proxy-theory envelope"
)
_PARAMETRIZATION_CITATION = (
    "Huber-Hurth-Lunghi JHEP 06 (2015) 176 for B -> X_s l l low-q2 SM "
    "validation; LO partonic C7/C9/C10 shape normalized to the YAML SM bin"
)
_NEEDS_HUMAN_PHYSICS = (
    "NEEDS-HUMAN-PHYSICS: full RS electroweak KK/Z/Z', lepton, C7 loop, "
    "NNLO/QED, charm-veto, power-correction and experimental-covariance "
    "matching is not available on ParameterPoint; v1 uses the documented "
    "C9/C10 proxy plus the existing C7 proxy."
)


@dataclass(frozen=True)
class BranchingFractionObservable:
    """Scaled B015 observable loaded through the scaffold anchor path."""

    block_key: str
    name: str
    source: str | None
    year: int | None
    value: float
    uncertainty: float
    raw_value: float
    raw_uncertainty: float
    units: str | None
    scale: float
    source_key: str | None
    source_url: str | None
    snapshot_path: str | None
    q2_region: str | None = None


@dataclass(frozen=True)
class InclusiveDileptonBudgetBand:
    """Low-q2 B015 total-BR budget for the HARD verdict."""

    source: str
    central_residual: float
    experimental_sigma: float
    sm_theory_sigma: float
    combined_exp_sm_sigma: float
    proxy_theory_sigma: float
    proxy_theory_fraction: float
    proxy_theory_rationale: str
    hard_veto_budget: float
    lower_edge: float
    upper_edge: float
    construction: str


@dataclass(frozen=True)
class B015Anchor:
    """Typed B015 anchor: active low-q2 bin plus total/high-q2 diagnostics."""

    total_hflav: BranchingFractionObservable
    total_pdg_listed: BranchingFractionObservable
    low_q2_experimental: BranchingFractionObservable
    high_q2_experimental: BranchingFractionObservable
    low_q2_sm_mumu: BranchingFractionObservable
    low_q2_sm_ee: BranchingFractionObservable
    q2_min_gev2: float
    q2_max_gev2: float
    budget_band: InclusiveDileptonBudgetBand

    @property
    def value(self) -> float:
        return self.low_q2_experimental.value

    @property
    def uncertainty(self) -> float:
        return self.low_q2_experimental.uncertainty

    @property
    def sm_value(self) -> float:
        return self.low_q2_sm_mumu.value

    @property
    def sm_uncertainty(self) -> float:
        return self.low_q2_sm_mumu.uncertainty

    @property
    def budget(self) -> float:
        return self.budget_band.hard_veto_budget


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
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise AnchorError(
            f"{process_id}: B015 anchor field {field_name!r}={value!r} is not numeric"
        ) from exc
    if not math.isfinite(number):
        raise AnchorError(
            f"{process_id}: B015 anchor field {field_name!r}={value!r} is not finite"
        )
    return number


def _positive_float(value: Any, *, process_id: str, field_name: str) -> float:
    number = _required_float(value, process_id=process_id, field_name=field_name)
    if number <= 0.0:
        raise AnchorError(f"{process_id}: B015 anchor field {field_name!r} <= 0")
    return number


_SYMMETRIC_UNCERTAINTY_RE = re.compile(r"^\s*(?:\+/-|±)\s*(?P<sigma>[0-9.eE+-]+)\s*$")
_ASYMMETRIC_UNCERTAINTY_RE = re.compile(
    r"^\s*\+?(?P<upper>[0-9.eE+-]+)\s*/\s*-(?P<lower>[0-9.eE+-]+)\s*$"
)
_LOW_Q2_REGION_RE = re.compile(
    r"(?P<lower>[0-9.]+)\s*<\s*q\^2\s*<\s*(?P<upper>[0-9.]+)\s*GeV\^2"
)


def _parse_symmetric_uncertainty(
    value: Any,
    *,
    process_id: str,
    field_name: str,
) -> float:
    if isinstance(value, str):
        match = _SYMMETRIC_UNCERTAINTY_RE.match(value)
        if match is not None:
            return _positive_float(
                match.group("sigma"),
                process_id=process_id,
                field_name=field_name,
            )
        asym = _ASYMMETRIC_UNCERTAINTY_RE.match(value)
        if asym is not None:
            upper = _positive_float(
                asym.group("upper"),
                process_id=process_id,
                field_name=f"{field_name}.upper",
            )
            lower = _positive_float(
                asym.group("lower"),
                process_id=process_id,
                field_name=f"{field_name}.lower",
            )
            return float(0.5 * (upper + lower))
    return _positive_float(value, process_id=process_id, field_name=field_name)


def _scale_from_units(units: str | None, *, process_id: str, field_name: str) -> float:
    if units == "10^-6":
        return 1.0e-6
    if units == "branching fraction":
        return 1.0
    raise AnchorError(
        f"{process_id}: {field_name} must use units '10^-6' or 'branching fraction', "
        f"got {units!r}"
    )


def _pdg_block(process_id: str) -> Mapping[str, Any]:
    data = load_full_yaml(process_id, family=_FAMILY)
    block = data.get("pdg_or_equivalent")
    if not isinstance(block, Mapping):
        raise AnchorError(
            f"{process_id}: expected mapping-shaped 'pdg_or_equivalent' for B015, "
            f"got {type(block).__name__}"
        )
    return block


def _entry_list(
    pdg_block: Mapping[str, Any],
    key: str,
    *,
    process_id: str,
) -> tuple[Mapping[str, Any], ...]:
    entries = pdg_block.get(key)
    if not isinstance(entries, list | tuple):
        raise AnchorError(f"{process_id}: pdg_or_equivalent.{key} must be a list")
    out: list[Mapping[str, Any]] = []
    for index, entry in enumerate(entries):
        if not isinstance(entry, Mapping):
            raise AnchorError(
                f"{process_id}: pdg_or_equivalent.{key}[{index}] is not a mapping"
            )
        out.append(entry)
    return tuple(out)


def _find_entry_index(
    entries: tuple[Mapping[str, Any], ...],
    *,
    key: str,
    expected: str,
    process_id: str,
    list_name: str,
) -> tuple[int, Mapping[str, Any]]:
    for index, entry in enumerate(entries):
        if str(entry.get(key)) == expected:
            return index, entry
    raise AnchorError(
        f"{process_id}: no pdg_or_equivalent.{list_name} entry with "
        f"{key}={expected!r}"
    )


def _load_scaffold_list_anchor(
    *,
    process_id: str,
    observable_name: str,
) -> tuple[Anchor, Mapping[str, Any], Mapping[str, Any]]:
    pdg_block = _pdg_block(process_id)
    index, entry = _find_entry_index(
        _entry_list(pdg_block, "observables", process_id=process_id),
        key="name",
        expected=observable_name,
        process_id=process_id,
        list_name="observables",
    )
    block_key = f"pdg_or_equivalent.observables[{index}]"
    virtual_block = {block_key: dict(entry)}
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
            uncertainty_key=_SCAFFOLD_NO_SCALAR_UNCERTAINTY_KEY,
        )
    finally:
        anchor_scaffold.load_pdg_block = original_load_pdg_block

    if scaffold_anchor.block_key != block_key:
        raise AnchorError(
            f"{process_id}: load_anchor selected {scaffold_anchor.block_key!r}, "
            f"expected {block_key!r} for B015 observable {observable_name!r}"
        )
    return scaffold_anchor, entry, pdg_block


def _load_branching_observable(
    *,
    process_id: str,
    observable_name: str,
) -> BranchingFractionObservable:
    scaffold_anchor, entry, pdg_block = _load_scaffold_list_anchor(
        process_id=process_id,
        observable_name=observable_name,
    )
    units = _optional_str(scaffold_anchor.units)
    scale = _scale_from_units(
        units,
        process_id=process_id,
        field_name=f"{observable_name}.units",
    )
    raw_value = _positive_float(
        scaffold_anchor.value,
        process_id=process_id,
        field_name=f"{observable_name}.value",
    )
    raw_uncertainty = _parse_symmetric_uncertainty(
        entry.get("uncertainty"),
        process_id=process_id,
        field_name=f"{observable_name}.uncertainty",
    )
    return BranchingFractionObservable(
        block_key=scaffold_anchor.block_key,
        name=observable_name,
        source=_optional_str(pdg_block.get("source")),
        year=_optional_int(entry.get("year", pdg_block.get("year"))),
        value=float(raw_value * scale),
        uncertainty=float(raw_uncertainty * scale),
        raw_value=float(raw_value),
        raw_uncertainty=float(raw_uncertainty),
        units=units,
        scale=float(scale),
        source_key=_optional_str(entry.get("source_key")),
        source_url=_optional_str(entry.get("source_url")),
        snapshot_path=_optional_str(entry.get("snapshot_path")),
        q2_region=_optional_str(entry.get("q2_region")),
    )


def _parse_low_q2_bounds(
    low_q2_experimental: BranchingFractionObservable,
    *,
    process_id: str,
) -> tuple[float, float]:
    region = low_q2_experimental.q2_region
    if region is None:
        raise AnchorError(f"{process_id}: low-q2 experimental anchor lacks q2_region")
    match = _LOW_Q2_REGION_RE.search(region)
    if match is None:
        raise AnchorError(
            f"{process_id}: cannot parse low-q2 q2_region {region!r}"
        )
    lower = _positive_float(
        match.group("lower"),
        process_id=process_id,
        field_name="low_q2.q2_min_gev2",
    )
    upper = _positive_float(
        match.group("upper"),
        process_id=process_id,
        field_name="low_q2.q2_max_gev2",
    )
    if upper <= lower:
        raise AnchorError(f"{process_id}: low-q2 upper edge must exceed lower edge")
    return float(lower), float(upper)


def _build_budget_band(
    *,
    experimental: BranchingFractionObservable,
    standard_model: BranchingFractionObservable,
) -> InclusiveDileptonBudgetBand:
    if experimental.uncertainty <= 0.0 or standard_model.uncertainty <= 0.0:
        raise AnchorError("B015: low-q2 experimental and SM uncertainties are required")
    central = abs(float(experimental.value) - float(standard_model.value))
    combined = math.sqrt(experimental.uncertainty**2 + standard_model.uncertainty**2)
    proxy_sigma = abs(float(standard_model.value)) * float(
        RARE_B_DILEPTON_INCLUSIVE_XS_PROXY_THEORY_UNCERTAINTY_FRACTION
    )
    budget = central + combined + proxy_sigma
    if budget <= 0.0 or not math.isfinite(budget):
        raise AnchorError("B015: constructed inclusive B -> X_s ll budget is invalid")
    return InclusiveDileptonBudgetBand(
        source=_BUDGET_SOURCE,
        central_residual=float(central),
        experimental_sigma=float(experimental.uncertainty),
        sm_theory_sigma=float(standard_model.uncertainty),
        combined_exp_sm_sigma=float(combined),
        proxy_theory_sigma=float(proxy_sigma),
        proxy_theory_fraction=float(
            RARE_B_DILEPTON_INCLUSIVE_XS_PROXY_THEORY_UNCERTAINTY_FRACTION
        ),
        proxy_theory_rationale=(
            RARE_B_DILEPTON_INCLUSIVE_XS_PROXY_THEORY_UNCERTAINTY_RATIONALE
        ),
        hard_veto_budget=float(budget),
        lower_edge=float(experimental.value - budget),
        upper_edge=float(experimental.value + budget),
        construction=(
            "|BR_exp(low-q2) - BR_SM_mumu[1,6]| + "
            "sqrt(sigma_exp^2 + sigma_SM^2) + "
            f"{RARE_B_DILEPTON_INCLUSIVE_XS_PROXY_THEORY_UNCERTAINTY_FRACTION:.0%}"
            "*BR_SM_mumu[1,6]; proxy term reuses the B016 C9/C10 "
            "NEEDS-HUMAN-PHYSICS treatment"
        ),
    )


def _load_b015_anchor(process_id: str) -> B015Anchor:
    total_hflav = _load_branching_observable(
        process_id=process_id,
        observable_name=_TOTAL_OBSERVABLE,
    )
    total_pdg = _load_branching_observable(
        process_id=process_id,
        observable_name=_PDG_TOTAL_OBSERVABLE,
    )
    low_exp = _load_branching_observable(
        process_id=process_id,
        observable_name=_LOW_Q2_EXPERIMENTAL_OBSERVABLE,
    )
    high_exp = _load_branching_observable(
        process_id=process_id,
        observable_name=_HIGH_Q2_EXPERIMENTAL_OBSERVABLE,
    )
    sm_mumu = _load_branching_observable(
        process_id=process_id,
        observable_name=_LOW_Q2_SM_MUMU_OBSERVABLE,
    )
    sm_ee = _load_branching_observable(
        process_id=process_id,
        observable_name=_LOW_Q2_SM_EE_OBSERVABLE,
    )
    q2_min, q2_max = _parse_low_q2_bounds(low_exp, process_id=process_id)
    return B015Anchor(
        total_hflav=total_hflav,
        total_pdg_listed=total_pdg,
        low_q2_experimental=low_exp,
        high_q2_experimental=high_exp,
        low_q2_sm_mumu=sm_mumu,
        low_q2_sm_ee=sm_ee,
        q2_min_gev2=q2_min,
        q2_max_gev2=q2_max,
        budget_band=_build_budget_band(
            experimental=low_exp,
            standard_model=sm_mumu,
        ),
    )


def _budget_result(predicted: float, anchor: B015Anchor) -> tuple[float, float, bool]:
    budget = float(anchor.budget)
    pull = float(predicted - anchor.value)
    ratio = abs(pull) / budget if budget > 0.0 else float("inf")
    return float(budget), float(ratio), bool(ratio <= 1.0)


@register
class Constraint:
    """Catalogued inclusive ``B -> X_s ell ell`` low-q2 branching-ratio constraint."""

    process_id = "B015"
    severity = Severity.HARD
    observable = "BR(B -> X_s ell+ ell-) low-q2 inclusive"

    def __init__(self) -> None:
        self.anchor = _load_b015_anchor(self.process_id)
        self.sm_inputs = rare_b_inclusive_xs_dilepton_default_inputs()
        self.sm_result = rare_b_inclusive_xs_mumu_sm_branching_fraction(
            sm_branching_fraction=self.anchor.sm_value,
            q2_min_gev2=self.anchor.q2_min_gev2,
            q2_max_gev2=self.anchor.q2_max_gev2,
            inputs=self.sm_inputs,
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
                    f"extra {_REQUIRED_EXTRA!r} absent; inclusive B -> X_s ell ell "
                    "constraint was not evaluated."
                ),
                diagnostics={
                    "evaluated": False,
                    "missing_extra": _REQUIRED_EXTRA,
                    "active_bin": "low_q2_1_to_6_gev2_mumu",
                    "q2_min_gev2": float(self.anchor.q2_min_gev2),
                    "q2_max_gev2": float(self.anchor.q2_max_gev2),
                    "sm_anchor_branching_fraction": float(self.anchor.sm_value),
                    "sm_formula_branching_fraction": float(
                        self.sm_result.branching_fraction
                    ),
                    "budget_source": self.anchor.budget_band.source,
                    "budget_proxy_theory_sigma": float(
                        self.anchor.budget_band.proxy_theory_sigma
                    ),
                    "budget_proxy_theory_fraction": float(
                        self.anchor.budget_band.proxy_theory_fraction
                    ),
                    "needs_human_physics": _NEEDS_HUMAN_PHYSICS,
                },
            )

        kk_ew_mass = point.get_extra(_OPTIONAL_EW_MASS_EXTRA)
        try:
            result = inclusive_b_to_xs_mumu_from_couplings(
                couplings,
                sm_branching_fraction=self.anchor.sm_value,
                q2_min_gev2=self.anchor.q2_min_gev2,
                q2_max_gev2=self.anchor.q2_max_gev2,
                m_kk_gev=None if kk_ew_mass is None else float(kk_ew_mass),
                inputs=self.sm_inputs,
            )
        except Exception as exc:  # noqa: BLE001 - constraints degrade cleanly
            return ConstraintResult(
                process_id=self.process_id,
                severity=self.severity,
                passes=True,
                sm_prediction=float(self.sm_result.branching_fraction),
                experimental=float(self.anchor.value),
                budget=float(self.anchor.budget),
                notes=(
                    "NOT EVALUATED - invalid quark_mass_basis_couplings for "
                    "the inclusive B -> X_s mu mu C7/C9/C10 proxy"
                ),
                diagnostics={
                    "evaluated": False,
                    "exception_type": type(exc).__name__,
                    "exception_message": str(exc),
                    "needs_human_physics": _NEEDS_HUMAN_PHYSICS,
                },
            )

        predicted = float(result.branching_fraction)
        budget, ratio, passes = _budget_result(predicted, self.anchor)

        diagnostics = dict(result.diagnostics)
        diagnostics.update(
            {
                "evaluated": True,
                "active_bin": "low_q2_1_to_6_gev2_mumu",
                "active_experimental_block": (
                    self.anchor.low_q2_experimental.block_key
                ),
                "active_sm_block": self.anchor.low_q2_sm_mumu.block_key,
                "active_experimental_source_key": (
                    self.anchor.low_q2_experimental.source_key
                ),
                "active_sm_source_key": self.anchor.low_q2_sm_mumu.source_key,
                "experimental_units": self.anchor.low_q2_experimental.units,
                "experimental_scale": float(self.anchor.low_q2_experimental.scale),
                "experimental_raw_value": float(
                    self.anchor.low_q2_experimental.raw_value
                ),
                "experimental_raw_uncertainty": float(
                    self.anchor.low_q2_experimental.raw_uncertainty
                ),
                "sm_anchor_branching_fraction": float(self.anchor.sm_value),
                "sm_formula_branching_fraction": float(
                    result.sm_branching_fraction
                ),
                "sm_formula_minus_anchor": float(
                    result.sm_branching_fraction - self.anchor.sm_value
                ),
                "sm_low_q2_ee_branching_fraction": float(
                    self.anchor.low_q2_sm_ee.value
                ),
                "total_hflav_branching_fraction": float(
                    self.anchor.total_hflav.value
                ),
                "total_hflav_uncertainty": float(
                    self.anchor.total_hflav.uncertainty
                ),
                "total_pdg_listed_branching_fraction": float(
                    self.anchor.total_pdg_listed.value
                ),
                "high_q2_experimental_branching_fraction": float(
                    self.anchor.high_q2_experimental.value
                ),
                "high_q2_experimental_uncertainty": float(
                    self.anchor.high_q2_experimental.uncertainty
                ),
                "high_q2_scope": (
                    "diagnostic_only_no_matching_SM_high_q2_anchor_in_B015_yaml"
                ),
                "budget_experimental_sigma": float(
                    self.anchor.budget_band.experimental_sigma
                ),
                "budget_sm_theory_sigma": float(
                    self.anchor.budget_band.sm_theory_sigma
                ),
                "budget_combined_exp_sm_sigma": float(
                    self.anchor.budget_band.combined_exp_sm_sigma
                ),
                "budget_central_residual": float(
                    self.anchor.budget_band.central_residual
                ),
                "budget_proxy_theory_sigma": float(
                    self.anchor.budget_band.proxy_theory_sigma
                ),
                "budget_proxy_theory_fraction": float(
                    self.anchor.budget_band.proxy_theory_fraction
                ),
                "budget_proxy_theory_rationale": (
                    self.anchor.budget_band.proxy_theory_rationale
                ),
                "budget_source": self.anchor.budget_band.source,
                "budget_construction": self.anchor.budget_band.construction,
                "parametrization_citation": _PARAMETRIZATION_CITATION,
                "rs_matching_assumption": RARE_B_DILEPTON_RS_MATCHING_ASSUMPTION_V1,
                "inclusive_limitations": RARE_B_DILEPTON_INCLUSIVE_XS_LIMITATION_V1,
                "needs_human_physics": _NEEDS_HUMAN_PHYSICS,
                "kk_ew_mass_extra_used": kk_ew_mass is not None,
                "wilson_coefficients": {
                    **dict(result.diagnostics.get("dilepton_wilson_coefficients", {})),
                    **dict(result.diagnostics.get("dipole_wilson_coefficients", {})),
                },
            }
        )

        return ConstraintResult(
            process_id=self.process_id,
            severity=self.severity,
            passes=passes,
            predicted=predicted,
            sm_prediction=float(result.sm_branching_fraction),
            experimental=float(self.anchor.value),
            ratio=ratio,
            budget=budget,
            notes=(
                "BR(B -> X_s ell+ ell-) uses the B015 low-q2 inclusive bin "
                "with a normalized partonic C7/C9/C10 shape. The RS "
                "contribution is the documented mass-basis b-s C9/C10 proxy "
                "plus the existing C7 proxy; full NNLO/QED/charm/covariance "
                "treatment is marked NEEDS-HUMAN-PHYSICS. The HARD ratio is "
                "the total low-q2 BR pull over the YAML exp/SM/proxy budget."
            ),
            diagnostics=diagnostics,
        )
