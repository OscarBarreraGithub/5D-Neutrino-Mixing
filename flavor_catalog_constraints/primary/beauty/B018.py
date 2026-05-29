"""B018 - LFU ratio ``R_K`` in ``B+ -> K+ ell+ ell-``.

Physics
-------
The active scalar is the HFLAV central-q2 row in ``B018.yaml``,

    R_K = BR(B+ -> K+ mu+ mu-) / BR(B+ -> K+ e+ e-),
    1.1 < q^2 < 6.0 GeV^2/c^4.

The SM LFU reference is evaluated as a same-bin ratio using the shared B016
exclusive ``B+ -> K+ mu+ mu-`` C9/C10 form-factor machinery reached through
``flavor_catalog_constraints.physics_adapters.rare_b_meson``.  In the v1
catalog proxy the electron denominator is the same SM short-distance rate, so
``R_K^SM = 1`` up to QED/electron-mass effects not implemented in the core.

The NP proxy follows the B016/B017 machinery:

    R_K^proxy = BR(B+ -> K+ mu+ mu-; SM + NP_mu proxy, bin)
                / BR(B+ -> K+ e+ e-; SM proxy, bin).

This is a lepton-nonuniversal stress proxy, not a model-complete RS LFU
prediction.  A truly lepton-universal C9/C10 shift would largely cancel in
``R_K``; that cancellation and the missing electron/muon-specific matching are
reported explicitly in diagnostics.

Severity
--------
HARD.  The ratio is an observed LFU measurement.  The HARD budget is built from
the B018 central-q2 YAML row as ``|R_exp - R_SM(proxy)| + sigma_exp`` so the
current SM-compatible measurement does not veto the SM point while sizeable LFU
proxy shifts can still be excluded.

Catalog sidecar
---------------
``flavor_catalog/processes/beauty/B018.yaml`` is the source of truth for the
HFLAV R_K anchors and provenance.  Numeric experimental values are loaded
through the scaffold ``load_anchor`` path and then parsed into the list-entry
shape used by this sidecar.

NEEDS-HUMAN-PHYSICS
-------------------
Full RS ``R_K`` matching requires electroweak KK/Z/Z', electron-vs-muon lepton
couplings, QED/radiative treatment, C7/nonlocal charm, scalar/tensor operators,
and experimental covariance inputs that are not available on ``ParameterPoint``.
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
    RARE_B_DILEPTON_RS_MATCHING_ASSUMPTION_V1,
    bplus_kplus_mumu_from_couplings,
    rare_b_to_k_dilepton_default_inputs,
    rare_b_to_k_mumu_sm_branching_fraction,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "beauty"
_REQUIRED_EXTRA = "quark_mass_basis_couplings"
_OPTIONAL_EW_MASS_EXTRA = "kk_ew_mass_gev"
_RK_LOW_OBSERVABLE = "R_K low-q2"
_RK_CENTRAL_OBSERVABLE = "R_K central-q2"
_RK_FULL_BELLE_OBSERVABLE = "R_K full-q2 Belle-only"
_SUPERSEDED_RK_OBSERVABLE = "superseded LHCb 2021 central-q2 R_K"
_SCAFFOLD_NO_SCALAR_UNCERTAINTY_KEY = "__b018_uncertainty_parsed_below__"
_BUDGET_SOURCE = (
    "flavor_catalog/processes/beauty/B018.yaml "
    "pdg_or_equivalent.observables[R_K central-q2]"
)
_PARAMETRIZATION_CITATION = (
    "HFLAV Dec. 2025 central-q2 R_K row from B018.yaml; Buras b->sll "
    "Hamiltonian and BCL-like B -> K f_+(q^2) q^2 integration reused from "
    "B016/B017"
)
_NEEDS_HUMAN_PHYSICS = (
    "NEEDS-HUMAN-PHYSICS: full RS R_K matching needs EW KK/Z/Z', "
    "electron-vs-muon lepton couplings, QED/radiative corrections, C7/dipole, "
    "nonlocal-charm, scalar/tensor and covariance inputs not available on "
    "ParameterPoint; v1 uses a documented muon-only-over-SM-electron C9/C10 "
    "stress proxy. A truly lepton-universal C9/C10 shift would largely cancel "
    "in R_K."
)


@dataclass(frozen=True)
class B018ObservableAnchor:
    """One scaled B018 list observable loaded through the scaffold path."""

    block_key: str
    name: str
    definition: str | None
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
    source_key: str | None
    source_url: str | None
    snapshot_path: str | None
    q2_region: str | None = None
    note: str | None = None
    significance: str | None = None

    @property
    def uncertainty(self) -> float:
        return float(0.5 * (self.uncertainty_upper + self.uncertainty_lower))


@dataclass(frozen=True)
class B018ObservableSet:
    """B018 observable rows before the SM-ratio budget is attached."""

    rk_low: B018ObservableAnchor
    rk_central: B018ObservableAnchor
    rk_full_belle: B018ObservableAnchor
    superseded_rk_2021: B018ObservableAnchor
    q2_min_gev2: float
    q2_max_gev2: float


@dataclass(frozen=True)
class B018BudgetBand:
    """Active central-q2 ``R_K`` budget from the B018 YAML measurement."""

    source: str
    central_residual: float
    experimental_sigma: float
    experimental_sigma_upper: float
    experimental_sigma_lower: float
    sm_lfu_ratio: float
    hard_veto_budget: float
    lower_edge: float
    upper_edge: float
    construction: str


@dataclass(frozen=True)
class B018Anchor:
    """Typed B018 anchor: active central-q2 ``R_K`` row and context rows."""

    rk_low: B018ObservableAnchor
    rk_central: B018ObservableAnchor
    rk_full_belle: B018ObservableAnchor
    superseded_rk_2021: B018ObservableAnchor
    q2_min_gev2: float
    q2_max_gev2: float
    budget_band: B018BudgetBand

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
            f"{process_id}: B018 anchor field {field_name!r}={value!r} is not numeric"
        ) from exc
    if not math.isfinite(number):
        raise AnchorError(
            f"{process_id}: B018 anchor field {field_name!r}={value!r} is not finite"
        )
    return number


def _positive_float(value: Any, *, process_id: str, field_name: str) -> float:
    number = _required_float(value, process_id=process_id, field_name=field_name)
    if number <= 0.0:
        raise AnchorError(f"{process_id}: B018 anchor field {field_name!r} <= 0")
    return number


_SYMMETRIC_UNCERTAINTY_RE = re.compile(r"\+/-\s*(?P<sigma>[0-9.eE+-]+)")
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
    if units == "dimensionless":
        return 1.0
    raise AnchorError(f"{process_id}: {field_name} has unsupported B018 units {units!r}")


def _pdg_block(process_id: str) -> Mapping[str, Any]:
    data = load_full_yaml(process_id, family=_FAMILY)
    block = data.get("pdg_or_equivalent")
    if not isinstance(block, Mapping):
        raise AnchorError(
            f"{process_id}: expected mapping-shaped 'pdg_or_equivalent' for B018, "
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
            f"expected {block_key!r} for B018 observable {observable_name!r}"
        )
    return scaffold_anchor, entry, pdg_block


def _load_observable(
    *,
    process_id: str,
    observable_name: str,
) -> B018ObservableAnchor:
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
    return B018ObservableAnchor(
        block_key=scaffold_anchor.block_key,
        name=observable_name,
        definition=_optional_str(entry.get("definition")),
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
        source_key=_optional_str(entry.get("source_key")),
        source_url=_optional_str(entry.get("source_url")),
        snapshot_path=_optional_str(entry.get("snapshot_path")),
        q2_region=_optional_str(entry.get("q2_region")),
        note=_optional_str(entry.get("note")),
        significance=_optional_str(entry.get("significance")),
    )


def _parse_q2_bounds(
    observable: B018ObservableAnchor,
    *,
    process_id: str,
) -> tuple[float, float]:
    if observable.q2_region is None:
        raise AnchorError(f"{process_id}: active B018 R_K row lacks q2_region")
    match = _Q2_REGION_RE.search(observable.q2_region)
    if match is None:
        raise AnchorError(
            f"{process_id}: cannot parse active B018 q2_region "
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
        raise AnchorError(f"{process_id}: active B018 q2 upper edge <= lower edge")
    return float(lower), float(upper)


def _load_b018_observables(process_id: str) -> B018ObservableSet:
    rk_low = _load_observable(process_id=process_id, observable_name=_RK_LOW_OBSERVABLE)
    rk_central = _load_observable(
        process_id=process_id,
        observable_name=_RK_CENTRAL_OBSERVABLE,
    )
    rk_full_belle = _load_observable(
        process_id=process_id,
        observable_name=_RK_FULL_BELLE_OBSERVABLE,
    )
    superseded_rk = _load_observable(
        process_id=process_id,
        observable_name=_SUPERSEDED_RK_OBSERVABLE,
    )
    q2_min, q2_max = _parse_q2_bounds(rk_central, process_id=process_id)
    return B018ObservableSet(
        rk_low=rk_low,
        rk_central=rk_central,
        rk_full_belle=rk_full_belle,
        superseded_rk_2021=superseded_rk,
        q2_min_gev2=q2_min,
        q2_max_gev2=q2_max,
    )


def _build_budget_band(
    *,
    active_rk: B018ObservableAnchor,
    sm_lfu_ratio: float,
) -> B018BudgetBand:
    sm_ratio = float(sm_lfu_ratio)
    if not math.isfinite(sm_ratio) or sm_ratio <= 0.0:
        raise AnchorError("B018: SM LFU ratio proxy must be positive and finite")
    if active_rk.uncertainty <= 0.0:
        raise AnchorError("B018: active R_K experimental uncertainty is required")
    central = abs(float(active_rk.value) - sm_ratio)
    budget = central + float(active_rk.uncertainty)
    if budget <= 0.0 or not math.isfinite(budget):
        raise AnchorError("B018: constructed R_K budget is invalid")
    return B018BudgetBand(
        source=_BUDGET_SOURCE,
        central_residual=float(central),
        experimental_sigma=float(active_rk.uncertainty),
        experimental_sigma_upper=float(active_rk.uncertainty_upper),
        experimental_sigma_lower=float(active_rk.uncertainty_lower),
        sm_lfu_ratio=sm_ratio,
        hard_veto_budget=float(budget),
        lower_edge=float(active_rk.value - budget),
        upper_edge=float(active_rk.value + budget),
        construction=(
            "|R_K^exp - R_K^SM(proxy)| + sigma_exp from the active B018.yaml "
            "central-q2 row; no B016 proxy-theory envelope is added because "
            "the LFU ratio budget is taken from the direct R_K measurement"
        ),
    )


def _build_b018_anchor(
    observables: B018ObservableSet,
    *,
    sm_lfu_ratio: float,
) -> B018Anchor:
    return B018Anchor(
        rk_low=observables.rk_low,
        rk_central=observables.rk_central,
        rk_full_belle=observables.rk_full_belle,
        superseded_rk_2021=observables.superseded_rk_2021,
        q2_min_gev2=float(observables.q2_min_gev2),
        q2_max_gev2=float(observables.q2_max_gev2),
        budget_band=_build_budget_band(
            active_rk=observables.rk_central,
            sm_lfu_ratio=sm_lfu_ratio,
        ),
    )


def _budget_result(predicted: float, anchor: B018Anchor) -> tuple[float, float, bool]:
    budget = float(anchor.budget)
    pull = float(predicted - anchor.value)
    ratio = abs(pull) / budget if budget > 0.0 else float("inf")
    return float(budget), float(ratio), bool(ratio <= 1.0)


def _anchor_diagnostics(anchor: B018Anchor) -> dict[str, Any]:
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
        "active_source_key": anchor.rk_central.source_key,
        "rk_low_q2_value": float(anchor.rk_low.value),
        "rk_low_q2_uncertainty": float(anchor.rk_low.uncertainty),
        "rk_low_q2_region": anchor.rk_low.q2_region,
        "rk_full_q2_belle_value": float(anchor.rk_full_belle.value),
        "rk_full_q2_belle_uncertainty": float(anchor.rk_full_belle.uncertainty),
        "superseded_2021_rk_value": float(anchor.superseded_rk_2021.value),
        "superseded_2021_rk_uncertainty": float(
            anchor.superseded_rk_2021.uncertainty
        ),
        "superseded_2021_significance": anchor.superseded_rk_2021.significance,
        "budget_experimental_sigma": float(anchor.budget_band.experimental_sigma),
        "budget_experimental_sigma_upper": float(
            anchor.budget_band.experimental_sigma_upper
        ),
        "budget_experimental_sigma_lower": float(
            anchor.budget_band.experimental_sigma_lower
        ),
        "budget_central_residual": float(anchor.budget_band.central_residual),
        "budget_sm_lfu_ratio": float(anchor.budget_band.sm_lfu_ratio),
        "budget_source": anchor.budget_band.source,
        "budget_construction": anchor.budget_band.construction,
    }


@register
class Constraint:
    """Catalogued B018 central-q2 ``R_K`` LFU ratio proxy."""

    process_id = "B018"
    severity = Severity.HARD
    observable = "R_K central-q2"

    def __init__(self) -> None:
        observables = _load_b018_observables(self.process_id)
        self.sm_inputs = rare_b_to_k_dilepton_default_inputs()
        self.sm_muon_result = rare_b_to_k_mumu_sm_branching_fraction(
            mode="bplus_kplus",
            q2_min_gev2=observables.q2_min_gev2,
            q2_max_gev2=observables.q2_max_gev2,
            inputs=self.sm_inputs,
        )
        self.sm_electron_proxy_branching_fraction = float(
            self.sm_muon_result.branching_fraction
        )
        self.sm_lfu_ratio = float(
            self.sm_muon_result.branching_fraction
            / self.sm_electron_proxy_branching_fraction
        )
        self.anchor = _build_b018_anchor(
            observables,
            sm_lfu_ratio=self.sm_lfu_ratio,
        )

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        couplings = point.get_extra(_REQUIRED_EXTRA)
        if couplings is None:
            diagnostics = _anchor_diagnostics(self.anchor)
            diagnostics.update(
                {
                    "evaluated": False,
                    "missing_extra": _REQUIRED_EXTRA,
                    "rk_proxy_definition": (
                        "BR(B+ -> K+ mu mu; SM+NP_mu proxy, bin) / "
                        "BR(B+ -> K+ e e; SM proxy, same bin)"
                    ),
                    "sm_muon_branching_fraction": float(
                        self.sm_muon_result.branching_fraction
                    ),
                    "sm_electron_proxy_branching_fraction": float(
                        self.sm_electron_proxy_branching_fraction
                    ),
                    "sm_lfu_ratio": float(self.sm_lfu_ratio),
                    "electron_denominator_treatment": (
                        "SM B -> K ll rate reused as the electron denominator; "
                        "electron mass and QED corrections are not implemented"
                    ),
                    "lepton_universal_cancellation_note": (
                        "A truly lepton-universal C9/C10 shift would largely "
                        "cancel in R_K; the v1 result is a muon-only stress proxy"
                    ),
                    "lepton_universal_rk_proxy": 1.0,
                    "needs_human_physics": _NEEDS_HUMAN_PHYSICS,
                }
            )
            return ConstraintResult(
                process_id=self.process_id,
                severity=self.severity,
                passes=True,
                sm_prediction=float(self.sm_lfu_ratio),
                experimental=float(self.anchor.value),
                budget=float(self.anchor.budget),
                notes=(
                    f"extra {_REQUIRED_EXTRA!r} absent; B018 central-q2 R_K "
                    "proxy was not evaluated."
                ),
                diagnostics=diagnostics,
            )

        kk_ew_mass = point.get_extra(_OPTIONAL_EW_MASS_EXTRA)
        try:
            result = bplus_kplus_mumu_from_couplings(
                couplings,
                q2_min_gev2=self.anchor.q2_min_gev2,
                q2_max_gev2=self.anchor.q2_max_gev2,
                m_kk_gev=None if kk_ew_mass is None else float(kk_ew_mass),
                inputs=self.sm_inputs,
            )
        except Exception as exc:  # noqa: BLE001 - constraints degrade cleanly
            diagnostics = _anchor_diagnostics(self.anchor)
            diagnostics.update(
                {
                    "evaluated": False,
                    "exception_type": type(exc).__name__,
                    "exception_message": str(exc),
                    "needs_human_physics": _NEEDS_HUMAN_PHYSICS,
                }
            )
            return ConstraintResult(
                process_id=self.process_id,
                severity=self.severity,
                passes=True,
                sm_prediction=float(self.sm_lfu_ratio),
                experimental=float(self.anchor.value),
                budget=float(self.anchor.budget),
                notes=(
                    "NOT EVALUATED - invalid quark_mass_basis_couplings for "
                    "the B018 B -> K central-q2 C9/C10 LFU proxy"
                ),
                diagnostics=diagnostics,
            )

        predicted = float(
            result.branching_fraction / self.sm_electron_proxy_branching_fraction
        )
        budget, ratio, passes = _budget_result(predicted, self.anchor)

        diagnostics = dict(result.diagnostics)
        diagnostics.update(_anchor_diagnostics(self.anchor))
        diagnostics.update(
            {
                "evaluated": True,
                "rk_proxy_definition": (
                    "BR(B+ -> K+ mu mu; SM+NP_mu proxy, bin) / "
                    "BR(B+ -> K+ e e; SM proxy, same bin)"
                ),
                "rk_proxy_numerator_branching_fraction": float(
                    result.branching_fraction
                ),
                "rk_proxy_denominator_branching_fraction": float(
                    self.sm_electron_proxy_branching_fraction
                ),
                "sm_muon_branching_fraction": float(
                    self.sm_muon_result.branching_fraction
                ),
                "sm_electron_proxy_branching_fraction": float(
                    self.sm_electron_proxy_branching_fraction
                ),
                "sm_lfu_ratio": float(self.sm_lfu_ratio),
                "rk_proxy_ratio_to_sm": float(result.ratio_to_sm),
                "np_shift_branching_fraction": float(
                    result.np_shift_branching_fraction
                ),
                "electron_denominator_treatment": (
                    "SM B -> K ll rate reused as the electron denominator; "
                    "electron mass and QED corrections are not implemented"
                ),
                "lepton_universal_cancellation_note": (
                    "A truly lepton-universal C9/C10 shift would largely "
                    "cancel in R_K; the v1 result is a muon-only stress proxy"
                ),
                "lepton_universal_rk_proxy": 1.0,
                "parametrization_citation": _PARAMETRIZATION_CITATION,
                "rs_matching_assumption": RARE_B_DILEPTON_RS_MATCHING_ASSUMPTION_V1,
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
            }
        )

        return ConstraintResult(
            process_id=self.process_id,
            severity=self.severity,
            passes=passes,
            predicted=predicted,
            sm_prediction=float(self.sm_lfu_ratio),
            experimental=float(self.anchor.value),
            ratio=ratio,
            budget=budget,
            notes=(
                "B018 uses the B018.yaml central-q2 charged R_K row. The "
                "prediction is a B016/B017 B+ -> K+ mu mu C9/C10 q^2-bin "
                "integral divided by an SM electron-denominator proxy. Full "
                "electron-vs-muon, QED, C7, nonlocal charm and covariance "
                "matching is marked NEEDS-HUMAN-PHYSICS; a truly "
                "lepton-universal C9/C10 shift would largely cancel in R_K."
            ),
            diagnostics=diagnostics,
        )
