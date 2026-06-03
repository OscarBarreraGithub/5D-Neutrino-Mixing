"""B017 - composite ``B -> K(*) ell+ ell-`` / ``R_K(*)`` context.

Physics
-------
The active scalar implemented here is the charged-mode ``R_K central-q2`` row
in ``B017.yaml``:

    1.1 < q^2 < 6.0 GeV^2/c^4.

The prediction reuses the B016 exclusive ``B+ -> K+ mu+ mu-`` C9/C10
form-factor machinery through
``flavor_catalog_constraints.physics_adapters.rare_b_meson``:

    R_K^proxy = BR(B+ -> K+ mu+ mu-; C9/C10, bin)
                / BR(B+ -> K+ mu+ mu-; SM C9/C10, same bin).

This is a proxy for the LFU ratio, not a production ``R_K`` calculation,
because ``ParameterPoint`` does not carry electron/muon non-universal
electroweak/lepton-sector inputs.  The B017 inclusive, charged total, neutral
``K*`` total, and ``R_K*`` rows are loaded from the YAML and reported in
diagnostics.  The neutral ``K*`` angular/global-fit likelihood is not evaluated
until a dedicated K* form-factor, nonlocal-charm, and covariance backend exists.

Severity
--------
HARD.  The active ``R_K`` proxy is compared with the B017 central-q2 LHCb anchor
using the B016/B015 30% C9/C10 proxy-theory envelope:
``|R_exp - R_SM| + sigma_exp + 0.30*R_SM``.

Catalog sidecar
---------------
``flavor_catalog/processes/beauty/B017.yaml`` is the source of truth for all
observable values and provenance.  Numeric anchors below are routed through the
scaffold ``load_anchor`` path and then parsed into the list-entry shape used by
this sidecar.

NEEDS-HUMAN-PHYSICS
-------------------
The Phase-3a RS light-Z contribution supplies vector/axial
``C9/C10/C9'/C10'`` Wilsons for both the muon numerator and electron
denominator.  K* form-factor, angular-basis, nonlocal-charm, and covariance
inputs remain deferred.
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
    RARE_B_DILEPTON_EXCLUSIVE_BK_LIMITATION_V1,
    RARE_B_DILEPTON_EXCLUSIVE_BK_PROXY_THEORY_UNCERTAINTY_FRACTION,
    RARE_B_DILEPTON_EXCLUSIVE_BK_PROXY_THEORY_UNCERTAINTY_RATIONALE,
    RARE_B_DILEPTON_INCLUSIVE_XS_LIMITATION_V1,
    RARE_B_DILEPTON_RS_MATCHING_ASSUMPTION_V1,
    bplus_kplus_mumu_from_rs_semileptonic_wilsons,
    rare_b_to_k_dilepton_default_inputs,
    rare_b_to_k_mumu_sm_branching_fraction,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "beauty"
_REQUIRED_EXTRA = "rs_semileptonic_wilsons"
_OPTIONAL_EW_MASS_EXTRA = "kk_ew_mass_gev"
_INCLUSIVE_TOTAL_OBSERVABLE = "BR(B -> X_s ell+ ell-)"
_CHARGED_TOTAL_OBSERVABLE = "BR(B+ -> K+ ell+ ell-)"
_KSTAR_TOTAL_OBSERVABLE = "BR(B0 -> K*(892)0 ell+ ell-)"
_RK_LOW_OBSERVABLE = "R_K low-q2"
_RK_CENTRAL_OBSERVABLE = "R_K central-q2"
_RKSTAR_LOW_OBSERVABLE = "R_K* low-q2"
_RKSTAR_CENTRAL_OBSERVABLE = "R_K* central-q2"
_SCAFFOLD_NO_SCALAR_UNCERTAINTY_KEY = "__b017_uncertainty_parsed_below__"
_SM_RK_PROXY = 1.0
_BUDGET_SOURCE = (
    "flavor_catalog/processes/beauty/B017.yaml "
    "pdg_or_equivalent.observables[R_K central-q2] "
    "+ B016 C9/C10 form-factor proxy-theory envelope"
)
_PARAMETRIZATION_CITATION = (
    "LHCb 2023 central-q2 R_K row from B017.yaml; Buras b->sll Hamiltonian "
    "and BCL-like B -> K f_+(q^2) q^2 integration reused from B016"
)
_NEEDS_HUMAN_PHYSICS = (
    "NEEDS-HUMAN-PHYSICS: Phase-3a supplies light-Z vector/axial "
    "C9/C10/C9'/C10' terms for the muon numerator and electron denominator; "
    "C7/dipole, K* form factors, angular basis, nonlocal-charm and "
    "experimental-covariance inputs remain deferred."
)


@dataclass(frozen=True)
class B017ObservableAnchor:
    """One scaled B017 list observable loaded through the scaffold path."""

    block_key: str
    name: str
    source: str | None
    year: int | None
    value: float
    uncertainty_upper: float
    uncertainty_lower: float
    raw_value: float
    raw_uncertainty_upper: float
    raw_uncertainty_lower: float
    units: str | None
    scale: float
    source_url: str | None
    snapshot_path: str | None
    q2_region: str | None = None
    note: str | None = None

    @property
    def uncertainty(self) -> float:
        return float(0.5 * (self.uncertainty_upper + self.uncertainty_lower))


@dataclass(frozen=True)
class B017BudgetBand:
    """Active central-q2 ``R_K`` proxy budget."""

    source: str
    central_residual: float
    experimental_sigma: float
    experimental_sigma_upper: float
    experimental_sigma_lower: float
    proxy_theory_sigma: float
    proxy_theory_fraction: float
    proxy_theory_rationale: str
    hard_veto_budget: float
    lower_edge: float
    upper_edge: float
    construction: str


@dataclass(frozen=True)
class B017Anchor:
    """Typed B017 anchor: active ``R_K`` row plus umbrella diagnostics."""

    inclusive_total: B017ObservableAnchor
    charged_total: B017ObservableAnchor
    kstar_total: B017ObservableAnchor
    rk_low: B017ObservableAnchor
    rk_central: B017ObservableAnchor
    rkstar_low: B017ObservableAnchor
    rkstar_central: B017ObservableAnchor
    q2_min_gev2: float
    q2_max_gev2: float
    budget_band: B017BudgetBand

    @property
    def value(self) -> float:
        return self.rk_central.value

    @property
    def uncertainty(self) -> float:
        return self.rk_central.uncertainty

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
            f"{process_id}: B017 anchor field {field_name!r}={value!r} is not numeric"
        ) from exc
    if not math.isfinite(number):
        raise AnchorError(
            f"{process_id}: B017 anchor field {field_name!r}={value!r} is not finite"
        )
    return number


def _positive_float(value: Any, *, process_id: str, field_name: str) -> float:
    number = _required_float(value, process_id=process_id, field_name=field_name)
    if number <= 0.0:
        raise AnchorError(f"{process_id}: B017 anchor field {field_name!r} <= 0")
    return number


_SYMMETRIC_UNCERTAINTY_RE = re.compile(r"(?:\+/-|±)\s*(?P<sigma>[0-9.eE+-]+)")
_ASYMMETRIC_UNCERTAINTY_RE = re.compile(
    r"\+(?P<upper>[0-9.eE+-]+)\s*/\s*-(?P<lower>[0-9.eE+-]+)"
)
_Q2_REGION_RE = re.compile(
    r"(?P<lower>[0-9.]+)\s*<\s*q\^2\s*<\s*(?P<upper>[0-9.]+)\s*GeV\^2"
)


def _parse_uncertainty_pair(
    value: Any,
    *,
    process_id: str,
    field_name: str,
) -> tuple[float, float]:
    if isinstance(value, str):
        upper_sq = 0.0
        lower_sq = 0.0
        matched = False
        for component in value.split(","):
            asym = _ASYMMETRIC_UNCERTAINTY_RE.search(component)
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
                upper_sq += upper * upper
                lower_sq += lower * lower
                matched = True
                continue
            sym = _SYMMETRIC_UNCERTAINTY_RE.search(component)
            if sym is not None:
                sigma = _positive_float(
                    sym.group("sigma"),
                    process_id=process_id,
                    field_name=field_name,
                )
                upper_sq += sigma * sigma
                lower_sq += sigma * sigma
                matched = True
        if matched:
            return float(math.sqrt(upper_sq)), float(math.sqrt(lower_sq))
    sigma = _positive_float(value, process_id=process_id, field_name=field_name)
    return float(sigma), float(sigma)


def _scale_from_units(units: str | None, *, process_id: str, field_name: str) -> float:
    if units == "10^-6":
        return 1.0e-6
    if units == "10^-7":
        return 1.0e-7
    if units in {"dimensionless", "branching fraction"}:
        return 1.0
    raise AnchorError(
        f"{process_id}: {field_name} has unsupported B017 units {units!r}"
    )


def _pdg_block(process_id: str) -> Mapping[str, Any]:
    data = load_full_yaml(process_id, family=_FAMILY)
    block = data.get("pdg_or_equivalent")
    if not isinstance(block, Mapping):
        raise AnchorError(
            f"{process_id}: expected mapping-shaped 'pdg_or_equivalent' for B017, "
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
            f"expected {block_key!r} for B017 observable {observable_name!r}"
        )
    return scaffold_anchor, entry, pdg_block


def _load_observable(
    *,
    process_id: str,
    observable_name: str,
) -> B017ObservableAnchor:
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
    raw_value = _required_float(
        scaffold_anchor.value,
        process_id=process_id,
        field_name=f"{observable_name}.value",
    )
    raw_uncertainty_upper, raw_uncertainty_lower = _parse_uncertainty_pair(
        entry.get("uncertainty"),
        process_id=process_id,
        field_name=f"{observable_name}.uncertainty",
    )
    return B017ObservableAnchor(
        block_key=scaffold_anchor.block_key,
        name=observable_name,
        source=_optional_str(pdg_block.get("source")),
        year=_optional_int(entry.get("year", pdg_block.get("year"))),
        value=float(raw_value * scale),
        uncertainty_upper=float(raw_uncertainty_upper * scale),
        uncertainty_lower=float(raw_uncertainty_lower * scale),
        raw_value=float(raw_value),
        raw_uncertainty_upper=float(raw_uncertainty_upper),
        raw_uncertainty_lower=float(raw_uncertainty_lower),
        units=units,
        scale=float(scale),
        source_url=_optional_str(entry.get("source_url")),
        snapshot_path=_optional_str(entry.get("snapshot_path")),
        q2_region=_optional_str(entry.get("q2_region")),
        note=_optional_str(entry.get("note")),
    )


def _parse_q2_bounds(
    observable: B017ObservableAnchor,
    *,
    process_id: str,
) -> tuple[float, float]:
    if observable.q2_region is None:
        raise AnchorError(f"{process_id}: active B017 R_K row lacks q2_region")
    match = _Q2_REGION_RE.search(observable.q2_region)
    if match is None:
        raise AnchorError(
            f"{process_id}: cannot parse active B017 q2_region "
            f"{observable.q2_region!r}"
        )
    lower = _positive_float(
        match.group("lower"),
        process_id=process_id,
        field_name="rk_central.q2_min_gev2",
    )
    upper = _positive_float(
        match.group("upper"),
        process_id=process_id,
        field_name="rk_central.q2_max_gev2",
    )
    if upper <= lower:
        raise AnchorError(f"{process_id}: active B017 q2 upper edge <= lower edge")
    return float(lower), float(upper)


def _build_budget_band(active_rk: B017ObservableAnchor) -> B017BudgetBand:
    if active_rk.uncertainty <= 0.0:
        raise AnchorError("B017: active R_K experimental uncertainty is required")
    central = abs(float(active_rk.value) - _SM_RK_PROXY)
    proxy_sigma = abs(_SM_RK_PROXY) * float(
        RARE_B_DILEPTON_EXCLUSIVE_BK_PROXY_THEORY_UNCERTAINTY_FRACTION
    )
    budget = central + float(active_rk.uncertainty) + proxy_sigma
    if budget <= 0.0 or not math.isfinite(budget):
        raise AnchorError("B017: constructed R_K proxy budget is invalid")
    return B017BudgetBand(
        source=_BUDGET_SOURCE,
        central_residual=float(central),
        experimental_sigma=float(active_rk.uncertainty),
        experimental_sigma_upper=float(active_rk.uncertainty_upper),
        experimental_sigma_lower=float(active_rk.uncertainty_lower),
        proxy_theory_sigma=float(proxy_sigma),
        proxy_theory_fraction=float(
            RARE_B_DILEPTON_EXCLUSIVE_BK_PROXY_THEORY_UNCERTAINTY_FRACTION
        ),
        proxy_theory_rationale=(
            RARE_B_DILEPTON_EXCLUSIVE_BK_PROXY_THEORY_UNCERTAINTY_RATIONALE
        ),
        hard_veto_budget=float(budget),
        lower_edge=float(active_rk.value - budget),
        upper_edge=float(active_rk.value + budget),
        construction=(
            "|R_K^exp - R_K^SM| + sigma_exp + "
            f"{RARE_B_DILEPTON_EXCLUSIVE_BK_PROXY_THEORY_UNCERTAINTY_FRACTION:.0%}"
            "*R_K^SM; the proxy term reuses the B016 C9/C10 "
            "NEEDS-HUMAN-PHYSICS budget treatment"
        ),
    )


def _load_b017_anchor(process_id: str) -> B017Anchor:
    inclusive_total = _load_observable(
        process_id=process_id,
        observable_name=_INCLUSIVE_TOTAL_OBSERVABLE,
    )
    charged_total = _load_observable(
        process_id=process_id,
        observable_name=_CHARGED_TOTAL_OBSERVABLE,
    )
    kstar_total = _load_observable(
        process_id=process_id,
        observable_name=_KSTAR_TOTAL_OBSERVABLE,
    )
    rk_low = _load_observable(process_id=process_id, observable_name=_RK_LOW_OBSERVABLE)
    rk_central = _load_observable(
        process_id=process_id,
        observable_name=_RK_CENTRAL_OBSERVABLE,
    )
    rkstar_low = _load_observable(
        process_id=process_id,
        observable_name=_RKSTAR_LOW_OBSERVABLE,
    )
    rkstar_central = _load_observable(
        process_id=process_id,
        observable_name=_RKSTAR_CENTRAL_OBSERVABLE,
    )
    q2_min, q2_max = _parse_q2_bounds(rk_central, process_id=process_id)
    return B017Anchor(
        inclusive_total=inclusive_total,
        charged_total=charged_total,
        kstar_total=kstar_total,
        rk_low=rk_low,
        rk_central=rk_central,
        rkstar_low=rkstar_low,
        rkstar_central=rkstar_central,
        q2_min_gev2=q2_min,
        q2_max_gev2=q2_max,
        budget_band=_build_budget_band(rk_central),
    )


def _budget_result(predicted: float, anchor: B017Anchor) -> tuple[float, float, bool]:
    budget = float(anchor.budget)
    pull = float(predicted - anchor.value)
    ratio = abs(pull) / budget if budget > 0.0 else float("inf")
    return float(budget), float(ratio), bool(ratio <= 1.0)


def _anchor_diagnostics(anchor: B017Anchor) -> dict[str, Any]:
    return {
        "active_observable": anchor.rk_central.name,
        "active_experimental_block": anchor.rk_central.block_key,
        "active_q2_region": anchor.rk_central.q2_region,
        "q2_min_gev2": float(anchor.q2_min_gev2),
        "q2_max_gev2": float(anchor.q2_max_gev2),
        "active_experimental_uncertainty_upper": float(
            anchor.rk_central.uncertainty_upper
        ),
        "active_experimental_uncertainty_lower": float(
            anchor.rk_central.uncertainty_lower
        ),
        "active_experimental_raw_value": float(anchor.rk_central.raw_value),
        "active_experimental_units": anchor.rk_central.units,
        "charged_total_branching_fraction_anchor": float(anchor.charged_total.value),
        "charged_total_uncertainty": float(anchor.charged_total.uncertainty),
        "inclusive_total_branching_fraction_anchor": float(
            anchor.inclusive_total.value
        ),
        "inclusive_total_uncertainty": float(anchor.inclusive_total.uncertainty),
        "kstar_total_branching_fraction_anchor": float(anchor.kstar_total.value),
        "kstar_total_uncertainty": float(anchor.kstar_total.uncertainty),
        "rk_low_q2_value": float(anchor.rk_low.value),
        "rk_low_q2_uncertainty": float(anchor.rk_low.uncertainty),
        "rkstar_low_q2_value": float(anchor.rkstar_low.value),
        "rkstar_low_q2_uncertainty": float(anchor.rkstar_low.uncertainty),
        "rkstar_central_q2_value": float(anchor.rkstar_central.value),
        "rkstar_central_q2_uncertainty": float(
            anchor.rkstar_central.uncertainty
        ),
        "budget_experimental_sigma": float(anchor.budget_band.experimental_sigma),
        "budget_experimental_sigma_upper": float(
            anchor.budget_band.experimental_sigma_upper
        ),
        "budget_experimental_sigma_lower": float(
            anchor.budget_band.experimental_sigma_lower
        ),
        "budget_central_residual": float(anchor.budget_band.central_residual),
        "budget_proxy_theory_sigma": float(anchor.budget_band.proxy_theory_sigma),
        "budget_proxy_theory_fraction": float(
            anchor.budget_band.proxy_theory_fraction
        ),
        "budget_proxy_theory_rationale": anchor.budget_band.proxy_theory_rationale,
        "budget_source": anchor.budget_band.source,
        "budget_construction": anchor.budget_band.construction,
    }


@register
class Constraint:
    """Catalogued B017 central-q2 ``R_K`` proxy constraint."""

    process_id = "B017"
    severity = Severity.HARD
    observable = "R_K central-q2 in B -> K(*) ell+ ell- context"

    def __init__(self) -> None:
        self.anchor = _load_b017_anchor(self.process_id)
        self.sm_inputs = rare_b_to_k_dilepton_default_inputs()
        self.sm_result = rare_b_to_k_mumu_sm_branching_fraction(
            mode="bplus_kplus",
            q2_min_gev2=self.anchor.q2_min_gev2,
            q2_max_gev2=self.anchor.q2_max_gev2,
            inputs=self.sm_inputs,
        )
        self.charged_total_sm_result = rare_b_to_k_mumu_sm_branching_fraction(
            mode="bplus_kplus",
            inputs=self.sm_inputs,
        )

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        rs_wilsons = point.get_extra(_REQUIRED_EXTRA)
        if rs_wilsons is None:
            diagnostics = _anchor_diagnostics(self.anchor)
            diagnostics.update(
                {
                    "evaluated": False,
                    "missing_extra": _REQUIRED_EXTRA,
                    "legacy_quark_mass_basis_couplings_present": (
                        point.get_extra("quark_mass_basis_couplings") is not None
                    ),
                    "rk_proxy_definition": (
                        "BR(B+ -> K+ mu mu; C9/C10, bin) / "
                        "BR(B+ -> K+ e e; C9/C10, same bin)"
                    ),
                    "rk_proxy_denominator_branching_fraction": float(
                        self.sm_result.branching_fraction
                    ),
                    "charged_total_sm_formula_branching_fraction": float(
                        self.charged_total_sm_result.branching_fraction
                    ),
                    "inclusive_scope": (
                        "context_anchor_only; B015 is the dedicated inclusive "
                        "low-q2 implementation"
                    ),
                    "kstar_angular_scope": (
                        "diagnostic_only; K* angular/global-fit likelihood "
                        "needs a dedicated backend"
                    ),
                    "needs_human_physics": _NEEDS_HUMAN_PHYSICS,
                }
            )
            return ConstraintResult(
                process_id=self.process_id,
                severity=self.severity,
                passes=True,
                sm_prediction=float(_SM_RK_PROXY),
                experimental=float(self.anchor.value),
                budget=float(self.anchor.budget),
                notes=(
                    f"extra {_REQUIRED_EXTRA!r} absent; B017 central-q2 R_K "
                    "proxy was not evaluated."
                ),
                diagnostics=diagnostics,
            )

        kk_ew_mass = point.get_extra(_OPTIONAL_EW_MASS_EXTRA)
        try:
            result = bplus_kplus_mumu_from_rs_semileptonic_wilsons(
                rs_wilsons,
                lepton="mu",
                q2_min_gev2=self.anchor.q2_min_gev2,
                q2_max_gev2=self.anchor.q2_max_gev2,
                matching_scale_gev=None if kk_ew_mass is None else float(kk_ew_mass),
                inputs=self.sm_inputs,
            )
            denominator_result = bplus_kplus_mumu_from_rs_semileptonic_wilsons(
                rs_wilsons,
                lepton="e",
                q2_min_gev2=self.anchor.q2_min_gev2,
                q2_max_gev2=self.anchor.q2_max_gev2,
                matching_scale_gev=None if kk_ew_mass is None else float(kk_ew_mass),
                inputs=self.sm_inputs,
            )
            charged_total_result = bplus_kplus_mumu_from_rs_semileptonic_wilsons(
                rs_wilsons,
                lepton="mu",
                matching_scale_gev=None if kk_ew_mass is None else float(kk_ew_mass),
                inputs=self.sm_inputs,
            )
        except Exception as exc:  # noqa: BLE001 - constraints degrade cleanly
            diagnostics = _anchor_diagnostics(self.anchor)
            diagnostics.update(
                {
                    "evaluated": False,
                    "invalid_extra": _REQUIRED_EXTRA,
                    "exception_type": type(exc).__name__,
                    "exception_message": str(exc),
                    "needs_human_physics": _NEEDS_HUMAN_PHYSICS,
                }
            )
            return ConstraintResult(
                process_id=self.process_id,
                severity=self.severity,
                passes=True,
                sm_prediction=float(_SM_RK_PROXY),
                experimental=float(self.anchor.value),
                budget=float(self.anchor.budget),
                notes=(
                    "NOT EVALUATED - invalid rs_semileptonic_wilsons for "
                    "the B017 B -> K central-q2 C9/C10 ratio"
                ),
                diagnostics=diagnostics,
            )

        predicted = float(result.branching_fraction / denominator_result.branching_fraction)
        budget, ratio, passes = _budget_result(predicted, self.anchor)

        diagnostics = dict(result.diagnostics)
        diagnostics.update(_anchor_diagnostics(self.anchor))
        diagnostics.update(
            {
                "evaluated": True,
                "rk_proxy_definition": (
                    "BR(B+ -> K+ mu mu; C9/C10, bin) / "
                    "BR(B+ -> K+ e e; C9/C10, same bin)"
                ),
                "rk_proxy_numerator_branching_fraction": float(
                    result.branching_fraction
                ),
                "rk_proxy_denominator_branching_fraction": float(
                    denominator_result.branching_fraction
                ),
                "rk_proxy_sm_branching_fraction": float(
                    result.sm_branching_fraction
                ),
                "rk_proxy_ratio_to_sm": float(predicted / _SM_RK_PROXY),
                "rk_proxy_muon_ratio_to_sm": float(result.ratio_to_sm),
                "rk_proxy_electron_ratio_to_sm": float(denominator_result.ratio_to_sm),
                "charged_total_branching_fraction": float(
                    charged_total_result.branching_fraction
                ),
                "charged_total_sm_formula_branching_fraction": float(
                    charged_total_result.sm_branching_fraction
                ),
                "charged_total_ratio_to_anchor": float(
                    charged_total_result.branching_fraction
                    / self.anchor.charged_total.value
                ),
                "inclusive_scope": (
                    "context_anchor_only; B015 is the dedicated inclusive "
                    "low-q2 implementation"
                ),
                "inclusive_limitations": RARE_B_DILEPTON_INCLUSIVE_XS_LIMITATION_V1,
                "kstar_angular_scope": (
                    "diagnostic_only; no K* angular semileptonic form-factor, "
                    "nonlocal-charm or covariance backend is available"
                ),
                "parametrization_citation": _PARAMETRIZATION_CITATION,
                "rs_matching_assumption": result.diagnostics.get(
                    "rs_semileptonic_matching_assumption",
                    RARE_B_DILEPTON_RS_MATCHING_ASSUMPTION_V1,
                ),
                "exclusive_limitations": RARE_B_DILEPTON_EXCLUSIVE_BK_LIMITATION_V1,
                "needs_human_physics": _NEEDS_HUMAN_PHYSICS,
                "kk_ew_mass_extra_used": kk_ew_mass is not None,
                "wilson_coefficients": (
                    {}
                    if result.wilsons is None
                    else {
                        key: complex(value)
                        for key, value in result.wilsons.wilsons.items()
                    }
                ),
                "denominator_wilson_coefficients": (
                    {}
                    if denominator_result.wilsons is None
                    else {
                        key: complex(value)
                        for key, value in denominator_result.wilsons.wilsons.items()
                    }
                ),
            }
        )

        return ConstraintResult(
            process_id=self.process_id,
            severity=self.severity,
            passes=passes,
            predicted=predicted,
            sm_prediction=float(_SM_RK_PROXY),
            experimental=float(self.anchor.value),
            ratio=ratio,
            budget=budget,
            notes=(
                "B017 uses the B017.yaml central-q2 charged R_K row as a HARD "
                "single-number proxy. The prediction is the B016 B+ -> K+ mu "
                "mu C9/C10 form-factor bin divided by the matching electron "
                "bin so lepton-universal Wilsons cancel. K* "
                "angular/R_K* and inclusive rows are carried as diagnostics; "
                "the full global-fit treatment is marked NEEDS-HUMAN-PHYSICS."
            ),
            diagnostics=diagnostics,
        )
