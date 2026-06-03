"""B019 - LFU ratio ``R_K*`` in ``B0 -> K*(892)0 ell+ ell-``.

Physics
-------
The active scalar is the HFLAV central-q2 row in ``B019.yaml``,

    R_K* = BR(B0 -> K*(892)0 mu+ mu-) / BR(B0 -> K*(892)0 e+ e-),
    1.1 < q^2 < 6.0 GeV^2/c^4.

The SM LFU reference is evaluated as a same-bin ratio using the dedicated
``B -> K*`` vector-mode adapter, which reuses the shared
``quarkConstraints.rare_b_dilepton`` C9/C10 Wilson matching and adds only a
lightweight K* form-factor normalization.  In the v1 catalog proxy the electron
denominator is the same SM short-distance K* rate, so ``R_K*^SM = 1`` up to
the QED/electron-mass caveat stored in ``B019.yaml``.

The NP ratio follows the B018 pattern with lepton-specific Wilson entries:

    R_K* = BR(B0 -> K*0 mu+ mu-; SM + NP_mu, bin)
           / BR(B0 -> K*0 e+ e-; SM + NP_e, same bin).

The electron denominator uses the same lightweight B -> K* core and the
``e`` entry of ``rs_semileptonic_wilsons.b_to_s_ll``; electron mass, QED,
and angular/form-factor covariance remain residual caveats.

Severity
--------
HARD.  The ratio is an observed LFU measurement.  The HARD budget is built from
the B019 central-q2 measurement and BIP SM/QED sidecar row as
``|R_exp - R_SM(sidecar)| + sigma_exp + sigma_SM_QED``.

Catalog sidecar
---------------
``flavor_catalog/processes/beauty/B019.yaml`` is the source of truth for the
HFLAV ``R_K*`` anchors, historical context rows, and BIP SM/QED reference.
Numeric values are loaded through the scaffold ``load_anchor`` path and then
parsed into the list-entry shape used by this sidecar.

NEEDS-HUMAN-PHYSICS
-------------------
Full RS ``R_K*`` matching requires electroweak KK/Z/Z', electron-vs-muon
lepton couplings, QED/radiative treatment, C7/nonlocal charm, scalar/tensor
operators, angular/form-factor covariance, S-wave treatment, and experimental
covariance inputs that are not available on ``ParameterPoint``.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
import re
from typing import Any, Mapping

from flavor_catalog_constraints import anchors as anchor_scaffold
from flavor_catalog_constraints.anchors import Anchor, AnchorError, load_full_yaml
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.rare_b_kstar_dilepton import (
    RARE_B_DILEPTON_EXCLUSIVE_BKSTAR_LIMITATION_V1,
    RARE_B_DILEPTON_RS_MATCHING_ASSUMPTION_V1,
    bzero_kstarzero_mumu_from_rs_semileptonic_wilsons,
    rare_b_to_kstar_dilepton_default_inputs,
    rare_b_to_kstar_mumu_sm_branching_fraction,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "beauty"
_REQUIRED_EXTRA = "rs_semileptonic_wilsons"
_OPTIONAL_EW_MASS_EXTRA = "kk_ew_mass_gev"
_RKSTAR_LOW_OBSERVABLE = "R_K* low-q2"
_RKSTAR_CENTRAL_OBSERVABLE = "R_K* central-q2"
_LHCB_2017_LOW_OBSERVABLE = "LHCb 2017 low-q2 R_K*"
_LHCB_2017_CENTRAL_OBSERVABLE = "LHCb 2017 central-q2 R_K*"
_SM_CENTRAL_THEORY = "BIP 2016 central-q2 SM R_K*"
_SM_LOW_THEORY = "BIP 2016 low-q2 SM R_K* for LHCb 2023 bin"
_SCAFFOLD_NO_SCALAR_UNCERTAINTY_KEY = "__b019_uncertainty_parsed_below__"
_BUDGET_SOURCE = (
    "flavor_catalog/processes/beauty/B019.yaml "
    "pdg_or_equivalent.observables[R_K* central-q2] "
    "+ pdg_or_equivalent.theory_inputs[BIP 2016 central-q2 SM R_K*]"
)
_PARAMETRIZATION_CITATION = (
    "HFLAV Dec. 2025 central-q2 R_K* row and BIP 2016 SM/QED row from "
    "B019.yaml; C9/C10 Wilson proxy reused from rare_b_dilepton with the "
    "append-only B -> K* vector form-factor normalization"
)
_NEEDS_HUMAN_PHYSICS = (
    "NEEDS-HUMAN-PHYSICS: Phase-3a supplies light-Z vector/axial "
    "C9/C10/C9'/C10' terms for the muon numerator and electron denominator; "
    "QED/radiative corrections, C7/dipole, nonlocal charm, scalar/tensor "
    "operators, K* angular/form-factor covariance, S-wave treatment and "
    "experimental covariance inputs remain deferred. Lepton-universal C9/C10 "
    "shifts largely cancel in R_K*."
)


@dataclass(frozen=True)
class B019ObservableAnchor:
    """One scaled B019 observable row loaded through the scaffold path."""

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
class B019TheoryAnchor:
    """One B019 SM/QED theory row loaded through the scaffold path."""

    block_key: str
    name: str
    value: float
    uncertainty: float
    raw_value: float
    raw_uncertainty: float
    source_key: str | None
    snapshot_path: str | None
    q2_region: str | None


@dataclass(frozen=True)
class B019ObservableSet:
    """B019 observable and theory rows before the budget is attached."""

    rkstar_low: B019ObservableAnchor
    rkstar_central: B019ObservableAnchor
    lhcb_2017_low: B019ObservableAnchor
    lhcb_2017_central: B019ObservableAnchor
    sm_central: B019TheoryAnchor
    sm_low: B019TheoryAnchor
    q2_min_gev2: float
    q2_max_gev2: float


@dataclass(frozen=True)
class B019BudgetBand:
    """Active central-q2 ``R_K*`` budget from B019 YAML rows."""

    source: str
    central_residual: float
    experimental_sigma: float
    experimental_sigma_upper: float
    experimental_sigma_lower: float
    sm_theory_sigma: float
    sm_lfu_ratio: float
    hard_veto_budget: float
    lower_edge: float
    upper_edge: float
    construction: str


@dataclass(frozen=True)
class B019Anchor:
    """Typed B019 anchor: active central-q2 ``R_K*`` row and context rows."""

    rkstar_low: B019ObservableAnchor
    rkstar_central: B019ObservableAnchor
    lhcb_2017_low: B019ObservableAnchor
    lhcb_2017_central: B019ObservableAnchor
    sm_central: B019TheoryAnchor
    sm_low: B019TheoryAnchor
    q2_min_gev2: float
    q2_max_gev2: float
    budget_band: B019BudgetBand

    @property
    def value(self) -> float:
        return self.rkstar_central.value

    @property
    def uncertainty(self) -> float:
        return self.rkstar_central.uncertainty

    @property
    def sm_value(self) -> float:
        return self.sm_central.value

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
            f"{process_id}: B019 anchor field {field_name!r}={value!r} is not numeric"
        ) from exc
    if not math.isfinite(number):
        raise AnchorError(
            f"{process_id}: B019 anchor field {field_name!r}={value!r} is not finite"
        )
    return number


def _positive_float(value: Any, *, process_id: str, field_name: str) -> float:
    number = _required_float(value, process_id=process_id, field_name=field_name)
    if number <= 0.0:
        raise AnchorError(f"{process_id}: B019 anchor field {field_name!r} <= 0")
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


def _parse_symmetric_uncertainty(
    value: Any,
    *,
    process_id: str,
    field_name: str,
) -> float:
    if isinstance(value, str):
        match = _SYMMETRIC_UNCERTAINTY_RE.search(value)
        if match is not None:
            return _positive_float(
                match.group("sigma"),
                process_id=process_id,
                field_name=field_name,
            )
    return _positive_float(value, process_id=process_id, field_name=field_name)


def _scale_from_units(units: str | None, *, process_id: str, field_name: str) -> float:
    if units in {None, "dimensionless"}:
        return 1.0
    raise AnchorError(f"{process_id}: {field_name} has unsupported B019 units {units!r}")


def _pdg_block(process_id: str) -> Mapping[str, Any]:
    data = load_full_yaml(process_id, family=_FAMILY)
    block = data.get("pdg_or_equivalent")
    if not isinstance(block, Mapping):
        raise AnchorError(
            f"{process_id}: expected mapping-shaped 'pdg_or_equivalent' for B019, "
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
    list_name: str,
    entry_name: str,
) -> tuple[Anchor, Mapping[str, Any], Mapping[str, Any]]:
    pdg_block = _pdg_block(process_id)
    index, entry = _find_entry_index(
        _entry_list(pdg_block, list_name, process_id=process_id),
        key="name",
        expected=entry_name,
        process_id=process_id,
        list_name=list_name,
    )
    block_key = f"pdg_or_equivalent.{list_name}[{index}]"
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
            f"expected {block_key!r} for B019 {list_name} entry {entry_name!r}"
        )
    return scaffold_anchor, entry, pdg_block


def _load_observable(
    *,
    process_id: str,
    observable_name: str,
) -> B019ObservableAnchor:
    scaffold_anchor, entry, pdg_block = _load_scaffold_list_anchor(
        process_id=process_id,
        list_name="observables",
        entry_name=observable_name,
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
    return B019ObservableAnchor(
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


def _load_theory(
    *,
    process_id: str,
    theory_name: str,
) -> B019TheoryAnchor:
    scaffold_anchor, entry, _ = _load_scaffold_list_anchor(
        process_id=process_id,
        list_name="theory_inputs",
        entry_name=theory_name,
    )
    raw_value = _required_float(
        scaffold_anchor.value,
        process_id=process_id,
        field_name=f"{theory_name}.value",
    )
    raw_uncertainty = _parse_symmetric_uncertainty(
        entry.get("uncertainty"),
        process_id=process_id,
        field_name=f"{theory_name}.uncertainty",
    )
    return B019TheoryAnchor(
        block_key=scaffold_anchor.block_key,
        name=theory_name,
        value=float(raw_value),
        uncertainty=float(raw_uncertainty),
        raw_value=float(raw_value),
        raw_uncertainty=float(raw_uncertainty),
        source_key=_optional_str(entry.get("source_key")),
        snapshot_path=_optional_str(entry.get("snapshot_path")),
        q2_region=_optional_str(entry.get("q2_region")),
    )


def _parse_q2_bounds(
    observable: B019ObservableAnchor,
    *,
    process_id: str,
) -> tuple[float, float]:
    if observable.q2_region is None:
        raise AnchorError(f"{process_id}: active B019 R_K* row lacks q2_region")
    match = _Q2_REGION_RE.search(observable.q2_region)
    if match is None:
        raise AnchorError(
            f"{process_id}: cannot parse active B019 q2_region "
            f"{observable.q2_region!r}"
        )
    lower = _positive_float(
        match.group("lower"),
        process_id=process_id,
        field_name="rkstar_central.q2_min_gev2",
    )
    upper = _positive_float(
        match.group("upper"),
        process_id=process_id,
        field_name="rkstar_central.q2_max_gev2",
    )
    if upper <= lower:
        raise AnchorError(f"{process_id}: active B019 q2 upper edge <= lower edge")
    return float(lower), float(upper)


def _load_b019_observables(process_id: str) -> B019ObservableSet:
    rkstar_low = _load_observable(
        process_id=process_id,
        observable_name=_RKSTAR_LOW_OBSERVABLE,
    )
    rkstar_central = _load_observable(
        process_id=process_id,
        observable_name=_RKSTAR_CENTRAL_OBSERVABLE,
    )
    lhcb_2017_low = _load_observable(
        process_id=process_id,
        observable_name=_LHCB_2017_LOW_OBSERVABLE,
    )
    lhcb_2017_central = _load_observable(
        process_id=process_id,
        observable_name=_LHCB_2017_CENTRAL_OBSERVABLE,
    )
    sm_central = _load_theory(
        process_id=process_id,
        theory_name=_SM_CENTRAL_THEORY,
    )
    sm_low = _load_theory(
        process_id=process_id,
        theory_name=_SM_LOW_THEORY,
    )
    q2_min, q2_max = _parse_q2_bounds(rkstar_central, process_id=process_id)
    return B019ObservableSet(
        rkstar_low=rkstar_low,
        rkstar_central=rkstar_central,
        lhcb_2017_low=lhcb_2017_low,
        lhcb_2017_central=lhcb_2017_central,
        sm_central=sm_central,
        sm_low=sm_low,
        q2_min_gev2=q2_min,
        q2_max_gev2=q2_max,
    )


def _build_budget_band(
    *,
    active_rkstar: B019ObservableAnchor,
    sm_central: B019TheoryAnchor,
) -> B019BudgetBand:
    if active_rkstar.uncertainty <= 0.0:
        raise AnchorError("B019: active R_K* experimental uncertainty is required")
    if sm_central.uncertainty <= 0.0:
        raise AnchorError("B019: active R_K* SM/QED uncertainty is required")
    sm_ratio = float(sm_central.value)
    if not math.isfinite(sm_ratio) or sm_ratio <= 0.0:
        raise AnchorError("B019: SM LFU ratio anchor must be positive and finite")
    central = abs(float(active_rkstar.value) - sm_ratio)
    budget = central + float(active_rkstar.uncertainty) + float(sm_central.uncertainty)
    if budget <= 0.0 or not math.isfinite(budget):
        raise AnchorError("B019: constructed R_K* budget is invalid")
    return B019BudgetBand(
        source=_BUDGET_SOURCE,
        central_residual=float(central),
        experimental_sigma=float(active_rkstar.uncertainty),
        experimental_sigma_upper=float(active_rkstar.uncertainty_upper),
        experimental_sigma_lower=float(active_rkstar.uncertainty_lower),
        sm_theory_sigma=float(sm_central.uncertainty),
        sm_lfu_ratio=sm_ratio,
        hard_veto_budget=float(budget),
        lower_edge=float(active_rkstar.value - budget),
        upper_edge=float(active_rkstar.value + budget),
        construction=(
            "|R_K*^exp - R_K*^SM(sidecar)| + sigma_exp + sigma_SM_QED from "
            "the active B019.yaml central-q2 row and BIP SM/QED theory row"
        ),
    )


def _build_b019_anchor(observables: B019ObservableSet) -> B019Anchor:
    return B019Anchor(
        rkstar_low=observables.rkstar_low,
        rkstar_central=observables.rkstar_central,
        lhcb_2017_low=observables.lhcb_2017_low,
        lhcb_2017_central=observables.lhcb_2017_central,
        sm_central=observables.sm_central,
        sm_low=observables.sm_low,
        q2_min_gev2=float(observables.q2_min_gev2),
        q2_max_gev2=float(observables.q2_max_gev2),
        budget_band=_build_budget_band(
            active_rkstar=observables.rkstar_central,
            sm_central=observables.sm_central,
        ),
    )


def _budget_result(predicted: float, anchor: B019Anchor) -> tuple[float, float, bool]:
    budget = float(anchor.budget)
    pull = float(predicted - anchor.value)
    ratio = abs(pull) / budget if budget > 0.0 else float("inf")
    return float(budget), float(ratio), bool(ratio <= 1.0)


def _anchor_diagnostics(anchor: B019Anchor) -> dict[str, Any]:
    return {
        "active_observable": anchor.rkstar_central.name,
        "active_experimental_block": anchor.rkstar_central.block_key,
        "active_q2_region": anchor.rkstar_central.q2_region,
        "q2_min_gev2": float(anchor.q2_min_gev2),
        "q2_max_gev2": float(anchor.q2_max_gev2),
        "active_experimental_uncertainty_upper": float(
            anchor.rkstar_central.uncertainty_upper
        ),
        "active_experimental_uncertainty_lower": float(
            anchor.rkstar_central.uncertainty_lower
        ),
        "active_experimental_raw_value": float(anchor.rkstar_central.raw_value),
        "active_experimental_units": anchor.rkstar_central.units,
        "active_source_key": anchor.rkstar_central.source_key,
        "rkstar_low_q2_value": float(anchor.rkstar_low.value),
        "rkstar_low_q2_uncertainty": float(anchor.rkstar_low.uncertainty),
        "rkstar_low_q2_region": anchor.rkstar_low.q2_region,
        "lhcb_2017_low_q2_value": float(anchor.lhcb_2017_low.value),
        "lhcb_2017_low_q2_uncertainty": float(anchor.lhcb_2017_low.uncertainty),
        "lhcb_2017_low_significance": anchor.lhcb_2017_low.significance,
        "lhcb_2017_central_q2_value": float(anchor.lhcb_2017_central.value),
        "lhcb_2017_central_q2_uncertainty": float(
            anchor.lhcb_2017_central.uncertainty
        ),
        "lhcb_2017_central_significance": anchor.lhcb_2017_central.significance,
        "sm_theory_block": anchor.sm_central.block_key,
        "sm_theory_q2_region": anchor.sm_central.q2_region,
        "sm_theory_lfu_ratio": float(anchor.sm_central.value),
        "sm_theory_qed_uncertainty": float(anchor.sm_central.uncertainty),
        "sm_low_q2_lfu_ratio": float(anchor.sm_low.value),
        "sm_low_q2_uncertainty": float(anchor.sm_low.uncertainty),
        "budget_experimental_sigma": float(anchor.budget_band.experimental_sigma),
        "budget_experimental_sigma_upper": float(
            anchor.budget_band.experimental_sigma_upper
        ),
        "budget_experimental_sigma_lower": float(
            anchor.budget_band.experimental_sigma_lower
        ),
        "budget_central_residual": float(anchor.budget_band.central_residual),
        "budget_sm_lfu_ratio": float(anchor.budget_band.sm_lfu_ratio),
        "budget_sm_theory_sigma": float(anchor.budget_band.sm_theory_sigma),
        "budget_source": anchor.budget_band.source,
        "budget_construction": anchor.budget_band.construction,
    }


@register
class Constraint:
    """Catalogued B019 central-q2 ``R_K*`` LFU ratio proxy."""

    process_id = "B019"
    severity = Severity.HARD
    observable = "R_K* central-q2"

    def __init__(self) -> None:
        observables = _load_b019_observables(self.process_id)
        self.sm_inputs = rare_b_to_kstar_dilepton_default_inputs()
        self.sm_muon_result = rare_b_to_kstar_mumu_sm_branching_fraction(
            mode="bzero_kstarzero",
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
        self.anchor = _build_b019_anchor(observables)

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
                    "rkstar_proxy_definition": (
                        "BR(B0 -> K*0 mu mu; SM+NP_mu proxy, bin) / "
                        "BR(B0 -> K*0 e e; SM+NP_e proxy, same bin)"
                    ),
                    "sm_muon_branching_fraction": float(
                        self.sm_muon_result.branching_fraction
                    ),
                    "sm_electron_proxy_branching_fraction": float(
                        self.sm_electron_proxy_branching_fraction
                    ),
                    "sm_lfu_ratio": float(self.sm_lfu_ratio),
                    "sm_formula_minus_sidecar": float(
                        self.sm_lfu_ratio - self.anchor.sm_value
                    ),
                    "electron_denominator_treatment": (
                        "B -> K* ll rate with the electron Wilson entry reused "
                        "as the electron denominator; electron mass and QED "
                        "corrections are represented only by the B019.yaml "
                        "SM/QED budget term"
                    ),
                    "lepton_universal_cancellation_note": (
                        "A lepton-universal C9/C10 shift largely cancels in R_K*"
                    ),
                    "lepton_universal_rkstar_proxy": 1.0,
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
                    f"extra {_REQUIRED_EXTRA!r} absent; B019 central-q2 R_K* "
                    "proxy was not evaluated."
                ),
                diagnostics=diagnostics,
            )

        kk_ew_mass = point.get_extra(_OPTIONAL_EW_MASS_EXTRA)
        try:
            result = bzero_kstarzero_mumu_from_rs_semileptonic_wilsons(
                rs_wilsons,
                lepton="mu",
                q2_min_gev2=self.anchor.q2_min_gev2,
                q2_max_gev2=self.anchor.q2_max_gev2,
                matching_scale_gev=None if kk_ew_mass is None else float(kk_ew_mass),
                inputs=self.sm_inputs,
            )
            denominator_result = bzero_kstarzero_mumu_from_rs_semileptonic_wilsons(
                rs_wilsons,
                lepton="e",
                q2_min_gev2=self.anchor.q2_min_gev2,
                q2_max_gev2=self.anchor.q2_max_gev2,
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
                sm_prediction=float(self.sm_lfu_ratio),
                experimental=float(self.anchor.value),
                budget=float(self.anchor.budget),
                notes=(
                    "NOT EVALUATED - invalid rs_semileptonic_wilsons for "
                    "the B019 B -> K* central-q2 C9/C10 LFU ratio"
                ),
                diagnostics=diagnostics,
            )

        predicted = float(
            result.branching_fraction / denominator_result.branching_fraction
        )
        budget, ratio, passes = _budget_result(predicted, self.anchor)

        diagnostics = dict(result.diagnostics)
        diagnostics.update(_anchor_diagnostics(self.anchor))
        diagnostics.update(
            {
                "evaluated": True,
                "rkstar_proxy_definition": (
                    "BR(B0 -> K*0 mu mu; SM+NP_mu proxy, bin) / "
                    "BR(B0 -> K*0 e e; SM+NP_e proxy, same bin)"
                ),
                "rkstar_proxy_numerator_branching_fraction": float(
                    result.branching_fraction
                ),
                "rkstar_proxy_denominator_branching_fraction": float(
                    denominator_result.branching_fraction
                ),
                "sm_muon_branching_fraction": float(
                    self.sm_muon_result.branching_fraction
                ),
                "sm_electron_proxy_branching_fraction": float(
                    self.sm_electron_proxy_branching_fraction
                ),
                "sm_lfu_ratio": float(self.sm_lfu_ratio),
                "sm_formula_minus_sidecar": float(
                    self.sm_lfu_ratio - self.anchor.sm_value
                ),
                "rkstar_proxy_ratio_to_sm": float(predicted / self.sm_lfu_ratio),
                "rkstar_proxy_muon_ratio_to_sm": float(result.ratio_to_sm),
                "rkstar_proxy_electron_ratio_to_sm": float(
                    denominator_result.ratio_to_sm
                ),
                "np_shift_branching_fraction": float(
                    result.np_shift_branching_fraction
                ),
                "electron_denominator_treatment": (
                    "B -> K* ll rate with the electron Wilson entry reused as "
                    "the electron denominator; electron mass and QED corrections "
                    "are represented only by the B019.yaml SM/QED budget term"
                ),
                "qed_caveat": (
                    "SM R_K* is unity in the implemented C9/C10 proxy; B019.yaml "
                    "carries the BIP +/-0.01 QED uncertainty."
                ),
                "lepton_universal_cancellation_note": (
                    "A lepton-universal C9/C10 shift largely cancels in R_K*"
                ),
                "lepton_universal_rkstar_proxy": 1.0,
                "parametrization_citation": _PARAMETRIZATION_CITATION,
                "rs_matching_assumption": result.diagnostics.get(
                    "rs_semileptonic_matching_assumption",
                    RARE_B_DILEPTON_RS_MATCHING_ASSUMPTION_V1,
                ),
                "exclusive_kstar_limitations": (
                    RARE_B_DILEPTON_EXCLUSIVE_BKSTAR_LIMITATION_V1
                ),
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
            sm_prediction=float(self.sm_lfu_ratio),
            experimental=float(self.anchor.value),
            ratio=ratio,
            budget=budget,
            notes=(
                "B019 uses the B019.yaml central-q2 R_K* row. The prediction "
                "is a B -> K* mu mu C9/C10 q^2-bin integral divided by the "
                "matching electron Wilson denominator. Full electron-vs-muon, QED, C7, "
                "nonlocal charm, angular/form-factor covariance and S-wave "
                "matching beyond Phase-3a remains deferred; lepton-universal "
                "C9/C10 shifts largely cancel in R_K*."
            ),
            diagnostics=diagnostics,
        )
