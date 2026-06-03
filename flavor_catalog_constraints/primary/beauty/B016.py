"""B016 - exclusive charged rare decay ``B+ -> K+ ell+ ell-``.

Physics
-------
The active observable is the charged HFLAV total exclusive ``B+ -> K+ ell ell``
average in ``B016.yaml``.  The companion neutral HFLAV average in the same YAML
sidecar is intentionally out of scope for B016 and must be handled by a
separate future entry.  The prediction uses the shared ``b -> s l l`` core
introduced for B005, reached through
``flavor_catalog_constraints.physics_adapters.rare_b_meson``.  The exclusive
rate is a q^2 integral of the C9/C10 short-distance amplitude with a BCL-like
``B -> K`` vector form factor ``f_+(q^2)``:

    dBR/dq^2 ~ tau_B G_F^2 alpha^2 |V_tb V_ts^*|^2
              lambda(B,K,q^2)^(3/2) (|C9 f_+|^2 + |C10 f_+|^2).

The stored low-q^2 number over 1.1 <= q^2 <= 6.0 GeV^2 is only an internal
formula benchmark of this C9/C10 integrator.  ``B016.yaml`` does not provide a
sourced partial-BR anchor for that bin, so the benchmark is not a data
validation and is not used in the HARD verdict.

Severity
--------
HARD.  The predicted charged total branching fraction is compared to the HFLAV
charged average with budget
``|BR_exp - BR_SM(formula)| + sigma_exp + sigma_proxy``.  ``sigma_proxy`` is a
documented 30% C9/C10-proxy theory envelope carried because B016 has no YAML SM
anchor, SM theory uncertainty, or covariance and the v1 formula omits C7,
nonlocal charm, and full form-factor/covariance inputs.

Catalog sidecar
---------------
``flavor_catalog/processes/beauty/B016.yaml`` is the source of truth for the
HFLAV charged average used by B016.  Numeric experimental values are loaded
from that sidecar, not hardcoded here.  The neutral average is catalogued there
but is not an anchor for this charged-mode constraint.

NEEDS-HUMAN-PHYSICS
-------------------
The Phase-3a RS light-Z contribution supplies rigorous vector/axial
``C9/C10/C9'/C10'`` Wilsons.  Dipole, nonlocal-charm, scalar/tensor, and
covariance inputs remain deferred for a global exclusive ``b -> s mu mu``
implementation.
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
    bplus_kplus_mumu_from_rs_semileptonic_wilsons,
    rare_b_to_k_dilepton_default_inputs,
    rare_b_to_k_mumu_sm_branching_fraction,
)
from flavor_catalog_constraints.registry import register
from quarkConstraints import rare_b_dilepton as rare_b_dilepton_core

_FAMILY = "beauty"
_REQUIRED_EXTRA = "rs_semileptonic_wilsons"
_OPTIONAL_EW_MASS_EXTRA = "kk_ew_mass_gev"
_CHARGED_OBSERVABLE = "BR(B+ -> K+ ell+ ell-)"
_LOW_Q2_BENCHMARK_Q2_MIN_GEV2 = 1.1
_LOW_Q2_BENCHMARK_Q2_MAX_GEV2 = 6.0
_LOW_Q2_BENCHMARK_LABEL = "internal_c9_c10_formula_benchmark"
_LOW_Q2_BENCHMARK_NOTE = (
    "Internal formula benchmark over 1.1 <= q^2 <= 6.0 GeV^2; B016.yaml has "
    "no sourced partial-BR anchor for this bin, so this is not a data validation."
)
_SCAFFOLD_NO_SCALAR_UNCERTAINTY_KEY = "__b016_uncertainty_parsed_below__"
_PROXY_THEORY_UNCERTAINTY_FRACTION = (
    rare_b_dilepton_core.RARE_B_DILEPTON_EXCLUSIVE_BK_PROXY_THEORY_UNCERTAINTY_FRACTION
)
_PROXY_THEORY_UNCERTAINTY_RATIONALE = (
    rare_b_dilepton_core.RARE_B_DILEPTON_EXCLUSIVE_BK_PROXY_THEORY_UNCERTAINTY_RATIONALE
)
_BUDGET_SOURCE = (
    "flavor_catalog/processes/beauty/B016.yaml "
    "pdg_or_equivalent.observables[BR(B+ -> K+ ell+ ell-)] "
    "+ C9/C10 form-factor SM residual + documented proxy theory envelope"
)
_PARAMETRIZATION_CITATION = (
    "Buras b->sll Hamiltonian; BCL-like B -> K f_+(q^2) q^2 integration; "
    "low-q^2 number is an internal formula benchmark, not a data validation"
)
_NEEDS_HUMAN_PHYSICS = (
    "NEEDS-HUMAN-PHYSICS: Phase-3a supplies the light-Z vector/axial "
    "C9/C10/C9'/C10' terms; C7/dipole, nonlocal-charm, scalar/tensor and "
    "experimental-covariance matching remain deferred."
)


@dataclass(frozen=True)
class BranchingFractionObservable:
    """Scaled B016 HFLAV observable loaded through the scaffold anchor path."""

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
    note: str | None = None


@dataclass(frozen=True)
class BToKBudgetBand:
    """Total-BR budget for the active charged B016 observable."""

    source: str
    central_residual: float
    experimental_sigma: float
    proxy_theory_sigma: float
    proxy_theory_fraction: float
    proxy_theory_rationale: str
    hard_veto_budget: float
    lower_edge: float
    upper_edge: float
    construction: str


@dataclass(frozen=True)
class LowQ2FormulaBenchmark:
    """Internal low-q^2 formula benchmark; not a YAML data validation."""

    label: str
    q2_min_gev2: float
    q2_max_gev2: float
    note: str


@dataclass(frozen=True)
class B016Anchor:
    """Typed B016 charged-mode anchor and active budget."""

    charged: BranchingFractionObservable
    low_q2_formula_benchmark: LowQ2FormulaBenchmark
    budget_band: BToKBudgetBand

    @property
    def value(self) -> float:
        return self.charged.value

    @property
    def uncertainty(self) -> float:
        return self.charged.uncertainty

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
            f"{process_id}: B016 anchor field {field_name!r}={value!r} is not numeric"
        ) from exc
    if not math.isfinite(number):
        raise AnchorError(
            f"{process_id}: B016 anchor field {field_name!r}={value!r} is not finite"
        )
    return number


def _positive_float(value: Any, *, process_id: str, field_name: str) -> float:
    number = _required_float(value, process_id=process_id, field_name=field_name)
    if number <= 0.0:
        raise AnchorError(f"{process_id}: B016 anchor field {field_name!r} <= 0")
    return number


_SYMMETRIC_UNCERTAINTY_RE = re.compile(r"^\s*(?:\+/-|±)\s*(?P<sigma>[0-9.eE+-]+)\s*$")


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
    return _positive_float(value, process_id=process_id, field_name=field_name)


def _scale_from_units(units: str | None, *, process_id: str, field_name: str) -> float:
    if units == "10^-7":
        return 1.0e-7
    if units == "branching fraction":
        return 1.0
    raise AnchorError(
        f"{process_id}: {field_name} must use units '10^-7' or 'branching fraction', "
        f"got {units!r}"
    )


def _pdg_block(process_id: str) -> Mapping[str, Any]:
    data = load_full_yaml(process_id, family=_FAMILY)
    block = data.get("pdg_or_equivalent")
    if not isinstance(block, Mapping):
        raise AnchorError(
            f"{process_id}: expected mapping-shaped 'pdg_or_equivalent' for B016, "
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
            f"expected {block_key!r} for B016 observable {observable_name!r}"
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
        note=_optional_str(entry.get("note")),
    )


def _build_low_q2_formula_benchmark() -> LowQ2FormulaBenchmark:
    return LowQ2FormulaBenchmark(
        label=_LOW_Q2_BENCHMARK_LABEL,
        q2_min_gev2=float(_LOW_Q2_BENCHMARK_Q2_MIN_GEV2),
        q2_max_gev2=float(_LOW_Q2_BENCHMARK_Q2_MAX_GEV2),
        note=_LOW_Q2_BENCHMARK_NOTE,
    )


def _build_budget_band(
    *,
    experimental: BranchingFractionObservable,
    formula_sm: float,
) -> BToKBudgetBand:
    if experimental.uncertainty <= 0.0:
        raise AnchorError("B016: charged branching-fraction uncertainty must be positive")
    central = abs(float(experimental.value) - float(formula_sm))
    proxy_theory_sigma = abs(float(formula_sm)) * float(
        _PROXY_THEORY_UNCERTAINTY_FRACTION
    )
    budget = central + float(experimental.uncertainty) + proxy_theory_sigma
    if budget <= 0.0 or not math.isfinite(budget):
        raise AnchorError("B016: constructed B -> K mu mu budget is invalid")
    return BToKBudgetBand(
        source=_BUDGET_SOURCE,
        central_residual=float(central),
        experimental_sigma=float(experimental.uncertainty),
        proxy_theory_sigma=float(proxy_theory_sigma),
        proxy_theory_fraction=float(_PROXY_THEORY_UNCERTAINTY_FRACTION),
        proxy_theory_rationale=_PROXY_THEORY_UNCERTAINTY_RATIONALE,
        hard_veto_budget=float(budget),
        lower_edge=float(experimental.value - budget),
        upper_edge=float(experimental.value + budget),
        construction=(
            "|BR_exp - BR_SM(formula)| + sigma_exp + "
            f"{_PROXY_THEORY_UNCERTAINTY_FRACTION:.0%}*BR_SM(formula); "
            "the proxy term covers the documented C9/C10-only, no-covariance "
            "NEEDS-HUMAN-PHYSICS limitation"
        ),
    )


def _load_b016_anchor(process_id: str, *, formula_sm: float) -> B016Anchor:
    charged = _load_branching_observable(
        process_id=process_id,
        observable_name=_CHARGED_OBSERVABLE,
    )
    return B016Anchor(
        charged=charged,
        low_q2_formula_benchmark=_build_low_q2_formula_benchmark(),
        budget_band=_build_budget_band(
            experimental=charged,
            formula_sm=float(formula_sm),
        ),
    )


def _budget_result(predicted: float, anchor: B016Anchor) -> tuple[float, float, bool]:
    budget = float(anchor.budget)
    pull = float(predicted - anchor.value)
    ratio = abs(pull) / budget if budget > 0.0 else float("inf")
    return float(budget), float(ratio), bool(ratio <= 1.0)


@register
class Constraint:
    """Catalogued charged ``B+ -> K+ ell ell`` branching-ratio constraint."""

    process_id = "B016"
    severity = Severity.HARD
    observable = "BR(B+ -> K+ ell+ ell-)"

    def __init__(self) -> None:
        self.sm_inputs = rare_b_to_k_dilepton_default_inputs()
        self.sm_result = rare_b_to_k_mumu_sm_branching_fraction(
            mode="bplus_kplus",
            inputs=self.sm_inputs,
        )
        self.low_q2_formula_benchmark_result = rare_b_to_k_mumu_sm_branching_fraction(
            mode="bplus_kplus",
            q2_min_gev2=_LOW_Q2_BENCHMARK_Q2_MIN_GEV2,
            q2_max_gev2=_LOW_Q2_BENCHMARK_Q2_MAX_GEV2,
            inputs=self.sm_inputs,
        )
        self.anchor = _load_b016_anchor(
            self.process_id,
            formula_sm=self.sm_result.branching_fraction,
        )

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        rs_wilsons = point.get_extra(_REQUIRED_EXTRA)
        if rs_wilsons is None:
            return ConstraintResult(
                process_id=self.process_id,
                severity=self.severity,
                passes=True,
                sm_prediction=float(self.sm_result.branching_fraction),
                experimental=float(self.anchor.value),
                budget=float(self.anchor.budget),
                notes=(
                    f"extra {_REQUIRED_EXTRA!r} absent; B+ -> K+ ell+ ell- "
                    "constraint was not evaluated."
                ),
                diagnostics={
                    "evaluated": False,
                    "missing_extra": _REQUIRED_EXTRA,
                    "legacy_quark_mass_basis_couplings_present": (
                        point.get_extra("quark_mass_basis_couplings") is not None
                    ),
                    "sm_formula_branching_fraction": float(
                        self.sm_result.branching_fraction
                    ),
                    "active_mode": "B+ -> K+ mu+ mu-",
                    "neutral_mode_scope": "out_of_scope_separate_entry_future",
                    "low_q2_formula_benchmark_branching_fraction": float(
                        self.low_q2_formula_benchmark_result.branching_fraction
                    ),
                    "low_q2_formula_benchmark_label": (
                        self.anchor.low_q2_formula_benchmark.label
                    ),
                    "low_q2_formula_benchmark_note": (
                        self.anchor.low_q2_formula_benchmark.note
                    ),
                    "budget_source": self.anchor.budget_band.source,
                    "budget_proxy_theory_sigma": float(
                        self.anchor.budget_band.proxy_theory_sigma
                    ),
                    "budget_proxy_theory_fraction": float(
                        self.anchor.budget_band.proxy_theory_fraction
                    ),
                    "budget_proxy_theory_rationale": (
                        self.anchor.budget_band.proxy_theory_rationale
                    ),
                    "needs_human_physics": _NEEDS_HUMAN_PHYSICS,
                },
            )

        kk_ew_mass = point.get_extra(_OPTIONAL_EW_MASS_EXTRA)
        try:
            result = bplus_kplus_mumu_from_rs_semileptonic_wilsons(
                rs_wilsons,
                lepton="mu",
                matching_scale_gev=None if kk_ew_mass is None else float(kk_ew_mass),
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
                    "NOT EVALUATED - invalid rs_semileptonic_wilsons for "
                    "the B -> K mu mu C9/C10 path"
                ),
                diagnostics={
                    "evaluated": False,
                    "invalid_extra": _REQUIRED_EXTRA,
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
                "experimental_block": self.anchor.charged.block_key,
                "experimental_source_key": self.anchor.charged.source_key,
                "experimental_units": self.anchor.charged.units,
                "experimental_scale": float(self.anchor.charged.scale),
                "experimental_raw_value": float(self.anchor.charged.raw_value),
                "experimental_raw_uncertainty": float(
                    self.anchor.charged.raw_uncertainty
                ),
                "active_mode": "B+ -> K+ mu+ mu-",
                "neutral_mode_scope": "out_of_scope_separate_entry_future",
                "sm_formula_branching_fraction": float(
                    result.sm_branching_fraction
                ),
                "sm_formula_minus_experiment": float(
                    result.sm_branching_fraction - self.anchor.value
                ),
                "low_q2_formula_benchmark_branching_fraction": float(
                    self.low_q2_formula_benchmark_result.branching_fraction
                ),
                "low_q2_formula_benchmark_q2_min_gev2": float(
                    self.anchor.low_q2_formula_benchmark.q2_min_gev2
                ),
                "low_q2_formula_benchmark_q2_max_gev2": float(
                    self.anchor.low_q2_formula_benchmark.q2_max_gev2
                ),
                "low_q2_formula_benchmark_label": (
                    self.anchor.low_q2_formula_benchmark.label
                ),
                "low_q2_formula_benchmark_note": (
                    self.anchor.low_q2_formula_benchmark.note
                ),
                "budget_experimental_sigma": float(
                    self.anchor.budget_band.experimental_sigma
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
                "rs_matching_assumption": (
                    result.wilsons.matching_assumption
                    if result.wilsons is not None
                    else RARE_B_DILEPTON_RS_MATCHING_ASSUMPTION_V1
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
                "BR(B+ -> K+ ell+ ell-) uses a C9/C10-only B -> K form-factor "
                "q^2 integral. Phase-3a RS semileptonic C9/C10/C9'/C10' "
                "Wilsons enter additively; full C7, nonlocal charm, scalar/tensor "
                "and covariance treatment remains deferred. The HARD "
                "ratio is the total BR pull from the HFLAV charged average over "
                "the charged-mode budget with an explicit proxy theory envelope; "
                "the neutral B0 mode is out of scope for B016."
            ),
            diagnostics=diagnostics,
        )
