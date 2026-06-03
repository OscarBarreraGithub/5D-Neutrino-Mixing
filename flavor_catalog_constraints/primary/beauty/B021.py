"""B021 - baryonic rare decay ``Lambda_b -> Lambda mu+ mu-``.

Physics
-------
The active observable is the LHCb high-q2 differential branching fraction in
``B021.yaml``,

    15 GeV^2 < q2 < 20 GeV^2,

converted to a bin-integrated branching fraction by multiplying by the catalog
bin width.  The prediction uses the baryonic extension of the shared
``b -> s l l`` C9/C10 machinery introduced for B005/B016, reached only through
``flavor_catalog_constraints.physics_adapters.rare_b_baryon``.  The
``Lambda_b -> Lambda`` form-factor normalization is the leading-HQET
Detmold-Lin-Meinel-Wingate ``F_+``/``F_-`` dipole fit.

Severity
--------
HARD.  The predicted high-q2 bin branching fraction is compared with the LHCb
bin anchor using the B016-style 30% proxy-theory treatment:
``|BR_exp - BR_SM(formula)| + sigma_exp + 0.30*BR_SM(formula)``.  For this
baryonic mode, that proxy envelope explicitly covers the leading-HQET
form-factor normalization/covariance limitation in addition to the C9/C10-only
RS proxy.

Catalog sidecar
---------------
``flavor_catalog/processes/beauty/B021.yaml`` is the source of truth for the
PDG total branching fraction diagnostic and the active high-q2 LHCb bin.  All
experimental numbers below are loaded through the scaffold anchor loader, not
hardcoded here.

NEEDS-HUMAN-PHYSICS
-------------------
The Phase-3a RS light-Z contribution supplies rigorous vector/axial
``C9/C10/C9'/C10'`` Wilsons.  A complete ``Lambda_b -> Lambda mu mu``
likelihood still needs C7/nonlocal charm, scalar/tensor terms, full baryonic
form-factor covariance, and experimental bin correlations.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
import re
from typing import Any, Mapping

from flavor_catalog_constraints import anchors as anchor_scaffold
from flavor_catalog_constraints.anchors import Anchor, AnchorError, load_full_yaml
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.rare_b_baryon import (
    RARE_B_BARYONIC_DILEPTON_FORM_FACTOR_CITATION,
    RARE_B_BARYONIC_DILEPTON_LIMITATION_V1,
    RARE_B_BARYONIC_DILEPTON_PROXY_THEORY_UNCERTAINTY_FRACTION,
    RARE_B_BARYONIC_DILEPTON_PROXY_THEORY_UNCERTAINTY_RATIONALE,
    lambdab_lambda_mumu_from_rs_semileptonic_wilsons,
    lambdab_lambda_mumu_sm_branching_fraction,
    rare_b_baryon_dilepton_default_inputs,
)
from flavor_catalog_constraints.physics_adapters.rare_b_meson import (
    RARE_B_DILEPTON_RS_MATCHING_ASSUMPTION_V1,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "beauty"
_REQUIRED_EXTRA = "rs_semileptonic_wilsons"
_OPTIONAL_EW_MASS_EXTRA = "kk_ew_mass_gev"
_TOTAL_OBSERVABLE = "BR(Lambda_b0 -> Lambda mu+ mu-)"
_HIGH_Q2_DIFFERENTIAL_OBSERVABLE = "dBR/dq2(Lambda_b0 -> Lambda mu+ mu-)"
_SCAFFOLD_NO_SCALAR_UNCERTAINTY_KEY = "__b021_uncertainty_parsed_below__"
_BUDGET_SOURCE = (
    "flavor_catalog/processes/beauty/B021.yaml "
    "dBR/dq2(Lambda_b0 -> Lambda mu+ mu-) high-q2 bin "
    "+ baryonic C9/C10 formula + 30% proxy-theory envelope"
)
_PARAMETRIZATION_CITATION = (
    "Detmold-Lin-Meinel-Wingate arXiv:1212.4827 leading-HQET "
    "Lambda_b -> Lambda F_+/- high-q2 form factors; Buras b->sll C9/C10 "
    "Hamiltonian convention reused from B005/B016"
)
_NEEDS_HUMAN_PHYSICS = (
    "NEEDS-HUMAN-PHYSICS: Phase-3a supplies the light-Z vector/axial "
    "C9/C10/C9'/C10' terms; C7/nonlocal charm, scalar/tensor, full baryonic "
    "form-factor covariance, and experimental-correlation matching remain "
    "deferred."
)


@dataclass(frozen=True)
class BaryonicBranchingObservable:
    """Scaled B021 observable loaded through the scaffold anchor path."""

    block_key: str
    name: str
    observable_type: str | None
    source: str | None
    year: int | None
    value: float
    uncertainty_upper: float
    uncertainty_lower: float
    raw_value: float
    raw_uncertainty_upper: float
    raw_uncertainty_lower: float
    average_differential_value: float | None
    average_differential_uncertainty_upper: float | None
    average_differential_uncertainty_lower: float | None
    units: str | None
    scale: float
    source_key: str | None
    source_url: str | None
    snapshot_path: str | None
    q2_region: str | None = None
    q2_min_gev2: float | None = None
    q2_max_gev2: float | None = None
    note: str | None = None

    @property
    def uncertainty(self) -> float:
        return float(0.5 * (self.uncertainty_upper + self.uncertainty_lower))

    @property
    def q2_bin_width_gev2(self) -> float | None:
        if self.q2_min_gev2 is None or self.q2_max_gev2 is None:
            return None
        return float(self.q2_max_gev2 - self.q2_min_gev2)


@dataclass(frozen=True)
class BaryonicDileptonBudgetBand:
    """High-q2 B021 bin budget for the HARD verdict."""

    source: str
    central_residual: float
    experimental_sigma_upper: float
    experimental_sigma_lower: float
    experimental_sigma_used: float
    proxy_theory_sigma: float
    proxy_theory_fraction: float
    proxy_theory_rationale: str
    hard_veto_budget: float
    lower_edge: float
    upper_edge: float
    construction: str


@dataclass(frozen=True)
class B021ObservableBundle:
    """B021 observables before the formula-dependent budget is constructed."""

    total_branching_fraction: BaryonicBranchingObservable
    high_q2_bin: BaryonicBranchingObservable

    @property
    def q2_min_gev2(self) -> float:
        if self.high_q2_bin.q2_min_gev2 is None:
            raise AnchorError("B021: high-q2 bin lacks q2_min_gev2")
        return float(self.high_q2_bin.q2_min_gev2)

    @property
    def q2_max_gev2(self) -> float:
        if self.high_q2_bin.q2_max_gev2 is None:
            raise AnchorError("B021: high-q2 bin lacks q2_max_gev2")
        return float(self.high_q2_bin.q2_max_gev2)


@dataclass(frozen=True)
class B021Anchor:
    """Typed B021 anchor: active high-q2 bin plus total-BR diagnostic."""

    total_branching_fraction: BaryonicBranchingObservable
    high_q2_bin: BaryonicBranchingObservable
    budget_band: BaryonicDileptonBudgetBand

    @property
    def value(self) -> float:
        return self.high_q2_bin.value

    @property
    def uncertainty(self) -> float:
        return self.high_q2_bin.uncertainty

    @property
    def q2_min_gev2(self) -> float:
        return float(self.high_q2_bin.q2_min_gev2)

    @property
    def q2_max_gev2(self) -> float:
        return float(self.high_q2_bin.q2_max_gev2)

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
            f"{process_id}: B021 anchor field {field_name!r}={value!r} is not numeric"
        ) from exc
    if not math.isfinite(number):
        raise AnchorError(
            f"{process_id}: B021 anchor field {field_name!r}={value!r} is not finite"
        )
    return number


def _positive_float(value: Any, *, process_id: str, field_name: str) -> float:
    number = _required_float(value, process_id=process_id, field_name=field_name)
    if number <= 0.0:
        raise AnchorError(f"{process_id}: B021 anchor field {field_name!r} <= 0")
    return number


_SYMMETRIC_UNCERTAINTY_RE = re.compile(r"\+/-\s*(?P<sigma>[0-9.eE+-]+)")
_ASYMMETRIC_UNCERTAINTY_RE = re.compile(
    r"\+(?P<upper>[0-9.eE+-]+)\s*/\s*-(?P<lower>[0-9.eE+-]+)"
)
_Q2_REGION_RE = re.compile(
    r"(?P<lower>[0-9.]+)\s*<\s*q\^?2\s*<\s*(?P<upper>[0-9.]+)\s*GeV\^2"
)


def _parse_uncertainty_pair(
    value: Any,
    *,
    process_id: str,
    field_name: str,
) -> tuple[float, float]:
    if isinstance(value, str):
        upper2 = 0.0
        lower2 = 0.0
        asym = _ASYMMETRIC_UNCERTAINTY_RE.search(value)
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
            upper2 += upper * upper
            lower2 += lower * lower
        for match in _SYMMETRIC_UNCERTAINTY_RE.finditer(value):
            sigma = _positive_float(
                match.group("sigma"),
                process_id=process_id,
                field_name=field_name,
            )
            upper2 += sigma * sigma
            lower2 += sigma * sigma
        if upper2 > 0.0 and lower2 > 0.0:
            return float(math.sqrt(upper2)), float(math.sqrt(lower2))
    sigma = _positive_float(value, process_id=process_id, field_name=field_name)
    return float(sigma), float(sigma)


def _scale_from_units(units: str | None, *, process_id: str, field_name: str) -> float:
    if units == "10^-6":
        return 1.0e-6
    if units == "10^-7 (GeV^2/c^4)^-1":
        return 1.0e-7
    if units == "branching fraction":
        return 1.0
    raise AnchorError(
        f"{process_id}: {field_name} has unsupported B021 units {units!r}"
    )


def _parse_q2_bounds(
    q2_region: str | None,
    *,
    process_id: str,
    field_name: str,
) -> tuple[float, float] | tuple[None, None]:
    if q2_region is None:
        return None, None
    match = _Q2_REGION_RE.search(q2_region)
    if match is None:
        raise AnchorError(f"{process_id}: cannot parse {field_name} {q2_region!r}")
    lower = _positive_float(
        match.group("lower"),
        process_id=process_id,
        field_name=f"{field_name}.q2_min_gev2",
    )
    upper = _positive_float(
        match.group("upper"),
        process_id=process_id,
        field_name=f"{field_name}.q2_max_gev2",
    )
    if upper <= lower:
        raise AnchorError(f"{process_id}: {field_name} upper edge must exceed lower")
    return float(lower), float(upper)


def _pdg_block(process_id: str) -> Mapping[str, Any]:
    data = load_full_yaml(process_id, family=_FAMILY)
    block = data.get("pdg_or_equivalent")
    if not isinstance(block, Mapping):
        raise AnchorError(
            f"{process_id}: expected mapping-shaped 'pdg_or_equivalent' for B021, "
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
            f"expected {block_key!r} for B021 observable {observable_name!r}"
        )
    return scaffold_anchor, entry, pdg_block


def _load_branching_observable(
    *,
    process_id: str,
    observable_name: str,
) -> BaryonicBranchingObservable:
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
    raw_upper, raw_lower = _parse_uncertainty_pair(
        entry.get("uncertainty"),
        process_id=process_id,
        field_name=f"{observable_name}.uncertainty",
    )
    q2_min, q2_max = _parse_q2_bounds(
        _optional_str(entry.get("q2_region")),
        process_id=process_id,
        field_name=f"{observable_name}.q2_region",
    )
    observable_type = _optional_str(entry.get("type"))
    is_differential = observable_type == "differential_branching_fraction"
    if is_differential:
        if q2_min is None or q2_max is None:
            raise AnchorError(f"{process_id}: differential B021 observable lacks q2_region")
        bin_width = float(q2_max - q2_min)
        value = raw_value * scale * bin_width
        upper = raw_upper * scale * bin_width
        lower = raw_lower * scale * bin_width
        avg_value = raw_value * scale
        avg_upper = raw_upper * scale
        avg_lower = raw_lower * scale
    else:
        value = raw_value * scale
        upper = raw_upper * scale
        lower = raw_lower * scale
        avg_value = avg_upper = avg_lower = None

    return BaryonicBranchingObservable(
        block_key=scaffold_anchor.block_key,
        name=observable_name,
        observable_type=observable_type,
        source=_optional_str(pdg_block.get("source")),
        year=_optional_int(entry.get("year", pdg_block.get("year"))),
        value=float(value),
        uncertainty_upper=float(upper),
        uncertainty_lower=float(lower),
        raw_value=float(raw_value),
        raw_uncertainty_upper=float(raw_upper),
        raw_uncertainty_lower=float(raw_lower),
        average_differential_value=avg_value,
        average_differential_uncertainty_upper=avg_upper,
        average_differential_uncertainty_lower=avg_lower,
        units=units,
        scale=float(scale),
        source_key=_optional_str(entry.get("source_key")),
        source_url=_optional_str(entry.get("source_url")),
        snapshot_path=_optional_str(entry.get("snapshot_path")),
        q2_region=_optional_str(entry.get("q2_region")),
        q2_min_gev2=q2_min,
        q2_max_gev2=q2_max,
        note=_optional_str(entry.get("note")),
    )


def _load_b021_observables(process_id: str) -> B021ObservableBundle:
    total = _load_branching_observable(
        process_id=process_id,
        observable_name=_TOTAL_OBSERVABLE,
    )
    high_q2 = _load_branching_observable(
        process_id=process_id,
        observable_name=_HIGH_Q2_DIFFERENTIAL_OBSERVABLE,
    )
    return B021ObservableBundle(
        total_branching_fraction=total,
        high_q2_bin=high_q2,
    )


def _build_budget_band(
    *,
    experimental: BaryonicBranchingObservable,
    formula_sm: float,
) -> BaryonicDileptonBudgetBand:
    if experimental.uncertainty_upper <= 0.0 or experimental.uncertainty_lower <= 0.0:
        raise AnchorError("B021: high-q2 experimental uncertainty must be positive")
    central = abs(float(experimental.value) - float(formula_sm))
    exp_sigma = max(
        float(experimental.uncertainty_upper),
        float(experimental.uncertainty_lower),
    )
    proxy_sigma = abs(float(formula_sm)) * float(
        RARE_B_BARYONIC_DILEPTON_PROXY_THEORY_UNCERTAINTY_FRACTION
    )
    budget = central + exp_sigma + proxy_sigma
    if budget <= 0.0 or not math.isfinite(budget):
        raise AnchorError("B021: constructed Lambda_b -> Lambda mu mu budget is invalid")
    return BaryonicDileptonBudgetBand(
        source=_BUDGET_SOURCE,
        central_residual=float(central),
        experimental_sigma_upper=float(experimental.uncertainty_upper),
        experimental_sigma_lower=float(experimental.uncertainty_lower),
        experimental_sigma_used=float(exp_sigma),
        proxy_theory_sigma=float(proxy_sigma),
        proxy_theory_fraction=float(
            RARE_B_BARYONIC_DILEPTON_PROXY_THEORY_UNCERTAINTY_FRACTION
        ),
        proxy_theory_rationale=(
            RARE_B_BARYONIC_DILEPTON_PROXY_THEORY_UNCERTAINTY_RATIONALE
        ),
        hard_veto_budget=float(budget),
        lower_edge=float(experimental.value - budget),
        upper_edge=float(experimental.value + budget),
        construction=(
            "|BR_exp(high-q2) - BR_SM(formula)| + max(sigma_exp_up, "
            f"sigma_exp_down) + "
            f"{RARE_B_BARYONIC_DILEPTON_PROXY_THEORY_UNCERTAINTY_FRACTION:.0%}"
            "*BR_SM(formula); proxy term covers C9/C10-only RS matching and "
            "baryonic form-factor normalization/covariance NEEDS-HUMAN-PHYSICS "
            "limitations"
        ),
    )


def _build_b021_anchor(
    observables: B021ObservableBundle,
    *,
    formula_sm: float,
) -> B021Anchor:
    return B021Anchor(
        total_branching_fraction=observables.total_branching_fraction,
        high_q2_bin=observables.high_q2_bin,
        budget_band=_build_budget_band(
            experimental=observables.high_q2_bin,
            formula_sm=float(formula_sm),
        ),
    )


def _load_b021_anchor(process_id: str, *, formula_sm: float) -> B021Anchor:
    return _build_b021_anchor(
        _load_b021_observables(process_id),
        formula_sm=float(formula_sm),
    )


def _budget_result(predicted: float, anchor: B021Anchor) -> tuple[float, float, bool]:
    budget = float(anchor.budget)
    pull = float(predicted - anchor.value)
    ratio = abs(pull) / budget if budget > 0.0 else float("inf")
    return float(budget), float(ratio), bool(ratio <= 1.0)


@register
class Constraint:
    """Catalogued baryonic ``Lambda_b -> Lambda mu mu`` high-q2 constraint."""

    process_id = "B021"
    severity = Severity.HARD
    observable = "BR(Lambda_b0 -> Lambda mu+ mu-) high-q2"

    def __init__(self) -> None:
        self.sm_inputs = rare_b_baryon_dilepton_default_inputs()
        observables = _load_b021_observables(self.process_id)
        self.sm_result = lambdab_lambda_mumu_sm_branching_fraction(
            q2_min_gev2=observables.q2_min_gev2,
            q2_max_gev2=observables.q2_max_gev2,
            inputs=self.sm_inputs,
        )
        self.anchor = _build_b021_anchor(
            observables,
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
                    f"extra {_REQUIRED_EXTRA!r} absent; Lambda_b -> Lambda mu mu "
                    "constraint was not evaluated."
                ),
                diagnostics={
                    "evaluated": False,
                    "missing_extra": _REQUIRED_EXTRA,
                    "legacy_quark_mass_basis_couplings_present": (
                        point.get_extra("quark_mass_basis_couplings") is not None
                    ),
                    "active_bin": "high_q2_15_to_20_gev2_mumu",
                    "q2_min_gev2": float(self.anchor.q2_min_gev2),
                    "q2_max_gev2": float(self.anchor.q2_max_gev2),
                    "sm_formula_branching_fraction": float(
                        self.sm_result.branching_fraction
                    ),
                    "sm_average_differential_branching_fraction": float(
                        self.sm_result.average_differential_branching_fraction
                    ),
                    "experimental_average_differential_branching_fraction": float(
                        self.anchor.high_q2_bin.average_differential_value
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
                    "form_factor_uncertainty": (
                        RARE_B_BARYONIC_DILEPTON_PROXY_THEORY_UNCERTAINTY_RATIONALE
                    ),
                    "needs_human_physics": _NEEDS_HUMAN_PHYSICS,
                },
            )

        kk_ew_mass = point.get_extra(_OPTIONAL_EW_MASS_EXTRA)
        try:
            result = lambdab_lambda_mumu_from_rs_semileptonic_wilsons(
                rs_wilsons,
                lepton="mu",
                q2_min_gev2=self.anchor.q2_min_gev2,
                q2_max_gev2=self.anchor.q2_max_gev2,
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
                    "the Lambda_b -> Lambda mu mu C9/C10 path"
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
                "active_bin": "high_q2_15_to_20_gev2_mumu",
                "experimental_block": self.anchor.high_q2_bin.block_key,
                "experimental_source_key": self.anchor.high_q2_bin.source_key,
                "experimental_units": self.anchor.high_q2_bin.units,
                "experimental_scale": float(self.anchor.high_q2_bin.scale),
                "experimental_raw_value": float(self.anchor.high_q2_bin.raw_value),
                "experimental_raw_uncertainty_upper": float(
                    self.anchor.high_q2_bin.raw_uncertainty_upper
                ),
                "experimental_raw_uncertainty_lower": float(
                    self.anchor.high_q2_bin.raw_uncertainty_lower
                ),
                "experimental_average_differential_branching_fraction": float(
                    self.anchor.high_q2_bin.average_differential_value
                ),
                "experimental_average_differential_uncertainty_upper": float(
                    self.anchor.high_q2_bin.average_differential_uncertainty_upper
                ),
                "experimental_average_differential_uncertainty_lower": float(
                    self.anchor.high_q2_bin.average_differential_uncertainty_lower
                ),
                "sm_formula_branching_fraction": float(result.sm_branching_fraction),
                "sm_formula_minus_experiment": float(
                    result.sm_branching_fraction - self.anchor.value
                ),
                "pdg_total_branching_fraction": float(
                    self.anchor.total_branching_fraction.value
                ),
                "pdg_total_branching_fraction_uncertainty": float(
                    self.anchor.total_branching_fraction.uncertainty
                ),
                "pdg_total_scope": "diagnostic_only_not_the_active_high_q2_bin",
                "budget_experimental_sigma_upper": float(
                    self.anchor.budget_band.experimental_sigma_upper
                ),
                "budget_experimental_sigma_lower": float(
                    self.anchor.budget_band.experimental_sigma_lower
                ),
                "budget_experimental_sigma_used": float(
                    self.anchor.budget_band.experimental_sigma_used
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
                "form_factor_citation": RARE_B_BARYONIC_DILEPTON_FORM_FACTOR_CITATION,
                "form_factor_uncertainty": (
                    RARE_B_BARYONIC_DILEPTON_PROXY_THEORY_UNCERTAINTY_RATIONALE
                ),
                "rs_matching_assumption": result.diagnostics.get(
                    "rs_semileptonic_matching_assumption",
                    RARE_B_DILEPTON_RS_MATCHING_ASSUMPTION_V1,
                ),
                "baryonic_limitations": RARE_B_BARYONIC_DILEPTON_LIMITATION_V1,
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
                "BR(Lambda_b -> Lambda mu mu) uses the B021 LHCb high-q2 "
                "bin and a C9/C10-only leading-HQET baryonic form-factor "
                "integral. Phase-3a RS semileptonic C9/C10/C9'/C10' Wilsons "
                "enter additively; C7/nonlocal charm and baryonic form-factor "
                "covariance remain deferred. The HARD ratio is the high-q2 bin BR "
                "pull over the 30% proxy-theory budget."
            ),
            diagnostics=diagnostics,
        )
