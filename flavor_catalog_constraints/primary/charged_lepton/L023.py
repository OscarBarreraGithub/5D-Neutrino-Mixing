"""L023 - muon-neutrino trident production.

Physics
-------
The measured observable is the ratio
``sigma(nu_mu N -> nu_mu N mu+ mu-) / sigma_SM``.  L023 evaluates the standard
heavy-mediator trident response

    sigma/sigma_SM = (C_V^2 + C_A^2) / (C_V_SM^2 + C_A_SM^2),
    C_V_SM = 1 + 4 sin^2(theta_W), C_A_SM = 1,

through ``flavor_catalog_constraints.physics_adapters.neutrino_trident``.

NEEDS-HUMAN-PHYSICS
-------------------
The current ``ParameterPoint`` does not carry the full RS EW KK/Z/Z' tower or
the lepton/neutrino neutral-current couplings needed for rigorous matching to
the ``nu_mu nu_mu mu mu`` vertex.  The NP path is therefore an explicit
effective-coupling or heavy-Z'-like proxy and is flagged in diagnostics.

Severity
--------
HARD.  CCFR is the strongest practical historical trident measurement in the
catalog sidecar; the active veto compares the predicted ratio to the CCFR
central value with the catalogued 95% CL reinterpretation converted to a
two-sided Gaussian budget.  CHARM-II and NuTeV are retained as diagnostics.

Catalog sidecar
---------------
``flavor_catalog/processes/charged_lepton/L023.yaml`` is the source of truth
for the measured ratios and the 95% CL CCFR reinterpretation.  List entries are
virtualized and routed through the scaffold ``load_anchor`` path before local
typed validation; no experimental number is hardcoded here.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from statistics import NormalDist
from typing import Any, Mapping

from flavor_catalog_constraints import anchors as anchor_scaffold
from flavor_catalog_constraints.anchors import (
    Anchor,
    AnchorError,
    load_anchor,
    load_full_yaml,
)
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.neutrino_trident import (
    TRIDENT_PARAMETRIZATION_CITATION,
    TRIDENT_PROXY_ASSUMPTION_V1,
    trident_sigma_ratio_from_lepton_input,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "charged_lepton"
_REQUIRED_EXTRA = "lepton_mass_basis_couplings"
_OPTIONAL_EW_MASS_EXTRA = "kk_ew_mass_gev"
_MEASURED_OBSERVABLES_LIST_KEY = "measured_observables"
_VALUES_LIST_KEY = "values"
_ACTIVE_EXPERIMENT = "CCFR"
_CCFR_CL_VALUE_ID = "Altmannshofer2014:L023:ccfr_exclusion_cl"
_EXPECTED_RATIO_UNITS = "dimensionless ratio"
_NO_SCALAR_UNCERTAINTY_KEY = "__l023_no_scalar_uncertainty__"
_UNEVALUATED_REASON = (
    "no neutrino-trident NP prediction available "
    "(RS nu_mu nu_mu mu mu neutral-current inputs not on ParameterPoint)"
)
_UNEVALUATED_NOTES = f"NOT EVALUATED - {_UNEVALUATED_REASON}"


@dataclass(frozen=True)
class TridentMeasurementAnchor:
    """Typed trident ``sigma_exp/sigma_SM`` measurement from L023.yaml."""

    block_key: str
    experiment: str
    observable: str
    value: float
    uncertainty_upper: float
    uncertainty_lower: float
    uncertainty_type: str
    units: str | None
    source: str | None
    primary_experiment_source: str | None
    source_url: str | None
    primary_source_url: str | None
    snapshot_path: str | None
    supporting_snapshot_path: str | None
    notes: str | None

    @property
    def uncertainty(self) -> float:
        return float(0.5 * (self.uncertainty_upper + self.uncertainty_lower))


@dataclass(frozen=True)
class TridentCLAnchor:
    """Typed confidence-level anchor for the CCFR reinterpretation."""

    block_key: str
    value_id: str
    value_percent: float
    source_url: str | None
    snapshot_path: str | None
    source_key: str | None


@dataclass(frozen=True)
class TridentBudgetBand:
    """Direction-aware L023 ratio budget around the active CCFR measurement."""

    source: str
    active_experiment: str
    confidence_level_percent: float
    gaussian_z: float
    sm_ratio: float
    sm_pull: float
    experimental_sigma_upper: float
    experimental_sigma_lower: float
    hard_veto_budget_upper: float
    hard_veto_budget_lower: float
    lower_edge: float
    upper_edge: float

    @property
    def hard_veto_budget(self) -> float:
        return float(max(self.hard_veto_budget_upper, self.hard_veto_budget_lower))


@dataclass(frozen=True)
class L023Anchor:
    """Typed L023 anchor bundle: measurements plus active CCFR budget."""

    measurements: tuple[TridentMeasurementAnchor, ...]
    active_measurement: TridentMeasurementAnchor
    ccfr_exclusion_cl: TridentCLAnchor
    budget_band: TridentBudgetBand

    @property
    def value(self) -> float:
        return self.active_measurement.value

    @property
    def uncertainty(self) -> float:
        return self.active_measurement.uncertainty

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
            f"{process_id}: L023 anchor field {field_name!r}={value!r} "
            "is not numeric"
        ) from exc
    if not math.isfinite(number):
        raise AnchorError(
            f"{process_id}: L023 anchor field {field_name!r}={value!r} "
            "is not finite"
        )
    return number


def _positive_float(value: Any, *, process_id: str, field_name: str) -> float:
    number = _required_float(value, process_id=process_id, field_name=field_name)
    if number <= 0.0:
        raise AnchorError(f"{process_id}: L023 anchor field {field_name!r} <= 0")
    return number


def _pdg_block(process_id: str) -> Mapping[str, Any]:
    data = load_full_yaml(process_id, family=_FAMILY)
    block = data.get("pdg_or_equivalent")
    if not isinstance(block, Mapping):
        raise AnchorError(
            f"{process_id}: 'pdg_or_equivalent' is not a mapping "
            f"(got {type(block).__name__})"
        )
    return block


def _list_entries(
    pdg_block: Mapping[str, Any],
    list_key: str,
    *,
    process_id: str,
) -> tuple[Mapping[str, Any], ...]:
    raw = pdg_block.get(list_key)
    if not isinstance(raw, list) or not raw:
        raise AnchorError(f"{process_id}: pdg_or_equivalent.{list_key} is not a list")
    entries: list[Mapping[str, Any]] = []
    for index, entry in enumerate(raw):
        if not isinstance(entry, Mapping):
            raise AnchorError(
                f"{process_id}: {list_key}[{index}] is not a mapping "
                f"(got {type(entry).__name__})"
            )
        entries.append(entry)
    return tuple(entries)


def _scaffold_anchor_for_entry(
    entry: Mapping[str, Any],
    *,
    process_id: str,
    list_key: str,
    index: int,
    uncertainty_key: str = "uncertainty",
) -> Anchor:
    block_key = f"{list_key}[{index}]"
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
        scaffold_anchor = load_anchor(
            process_id,
            family=_FAMILY,
            candidates=(block_key,),
            uncertainty_key=uncertainty_key,
        )
    finally:
        anchor_scaffold.load_pdg_block = original_load_pdg_block

    if scaffold_anchor.block_key != block_key:
        raise AnchorError(
            f"{process_id}: load_anchor selected {scaffold_anchor.block_key!r}, "
            f"expected {block_key!r} for L023 {list_key}"
        )
    return scaffold_anchor


def _measurement_uncertainties(
    entry: Mapping[str, Any],
    scaffold_anchor: Anchor,
    *,
    process_id: str,
    block_key: str,
) -> tuple[float, float, str]:
    uncertainty_type = str(entry.get("uncertainty_type", "symmetric"))
    if uncertainty_type == "symmetric":
        if scaffold_anchor.uncertainty is None:
            raise AnchorError(f"{process_id}: {block_key} has no uncertainty")
        sigma = _positive_float(
            scaffold_anchor.uncertainty,
            process_id=process_id,
            field_name=f"{block_key}.uncertainty",
        )
        return sigma, sigma, uncertainty_type
    if uncertainty_type == "asymmetric":
        upper = _positive_float(
            entry.get("uncertainty_plus"),
            process_id=process_id,
            field_name=f"{block_key}.uncertainty_plus",
        )
        lower = _positive_float(
            entry.get("uncertainty_minus"),
            process_id=process_id,
            field_name=f"{block_key}.uncertainty_minus",
        )
        return upper, lower, uncertainty_type
    raise AnchorError(
        f"{process_id}: {block_key} has unsupported uncertainty_type "
        f"{uncertainty_type!r}"
    )


def _measurement_anchor(
    entry: Mapping[str, Any],
    *,
    process_id: str,
    index: int,
) -> TridentMeasurementAnchor:
    scaffold_anchor = _scaffold_anchor_for_entry(
        entry,
        process_id=process_id,
        list_key=_MEASURED_OBSERVABLES_LIST_KEY,
        index=index,
    )
    block_key = scaffold_anchor.block_key
    units = scaffold_anchor.units
    if units != _EXPECTED_RATIO_UNITS:
        raise AnchorError(
            f"{process_id}: expected units {_EXPECTED_RATIO_UNITS!r} for "
            f"{block_key}, got {units!r}"
        )
    experiment = entry.get("experiment")
    observable = entry.get("observable")
    if not isinstance(experiment, str) or not experiment:
        raise AnchorError(f"{process_id}: {block_key} has no experiment")
    if not isinstance(observable, str) or not observable:
        raise AnchorError(f"{process_id}: {block_key} has no observable")
    upper, lower, uncertainty_type = _measurement_uncertainties(
        entry,
        scaffold_anchor,
        process_id=process_id,
        block_key=block_key,
    )
    return TridentMeasurementAnchor(
        block_key=block_key,
        experiment=experiment,
        observable=observable,
        value=_required_float(
            scaffold_anchor.value,
            process_id=process_id,
            field_name=f"{block_key}.value",
        ),
        uncertainty_upper=float(upper),
        uncertainty_lower=float(lower),
        uncertainty_type=uncertainty_type,
        units=units,
        source=_optional_str(scaffold_anchor.source),
        primary_experiment_source=_optional_str(
            entry.get("primary_experiment_source")
        ),
        source_url=_optional_str(scaffold_anchor.source_url),
        primary_source_url=_optional_str(entry.get("primary_source_url")),
        snapshot_path=_optional_str(scaffold_anchor.snapshot_path),
        supporting_snapshot_path=_optional_str(
            entry.get("supporting_snapshot_path")
        ),
        notes=_optional_str(entry.get("notes")),
    )


def _load_measurements(process_id: str) -> tuple[TridentMeasurementAnchor, ...]:
    pdg = _pdg_block(process_id)
    entries = _list_entries(
        pdg,
        _MEASURED_OBSERVABLES_LIST_KEY,
        process_id=process_id,
    )
    return tuple(
        _measurement_anchor(entry, process_id=process_id, index=index)
        for index, entry in enumerate(entries)
    )


def _find_measurement(
    measurements: tuple[TridentMeasurementAnchor, ...],
    experiment: str,
    *,
    process_id: str,
) -> TridentMeasurementAnchor:
    for measurement in measurements:
        if measurement.experiment == experiment:
            return measurement
    present = [measurement.experiment for measurement in measurements]
    raise AnchorError(
        f"{process_id}: no {experiment!r} trident measurement "
        f"(present experiments: {present})"
    )


def _load_ccfr_cl_anchor(process_id: str) -> TridentCLAnchor:
    pdg = _pdg_block(process_id)
    entries = _list_entries(pdg, _VALUES_LIST_KEY, process_id=process_id)
    for index, entry in enumerate(entries):
        if entry.get("value_id") != _CCFR_CL_VALUE_ID:
            continue
        scaffold_anchor = _scaffold_anchor_for_entry(
            entry,
            process_id=process_id,
            list_key=_VALUES_LIST_KEY,
            index=index,
            uncertainty_key=_NO_SCALAR_UNCERTAINTY_KEY,
        )
        return TridentCLAnchor(
            block_key=scaffold_anchor.block_key,
            value_id=str(entry["value_id"]),
            value_percent=_positive_float(
                scaffold_anchor.value,
                process_id=process_id,
                field_name=f"{scaffold_anchor.block_key}.value",
            ),
            source_url=_optional_str(scaffold_anchor.source_url),
            snapshot_path=_optional_str(scaffold_anchor.snapshot_path),
            source_key=_optional_str(entry.get("source_key")),
        )
    raise AnchorError(f"{process_id}: missing values entry {_CCFR_CL_VALUE_ID!r}")


def _build_budget_band(
    *,
    active: TridentMeasurementAnchor,
    cl_anchor: TridentCLAnchor,
    process_id: str,
) -> TridentBudgetBand:
    cl_fraction = cl_anchor.value_percent / 100.0
    if not 0.0 < cl_fraction < 1.0:
        raise AnchorError(
            f"{process_id}: confidence level must be between 0 and 100 percent"
        )
    gaussian_z = NormalDist().inv_cdf(0.5 + 0.5 * cl_fraction)
    upper_budget = gaussian_z * active.uncertainty_upper
    lower_budget = gaussian_z * active.uncertainty_lower
    if upper_budget <= 0.0 or lower_budget <= 0.0:
        raise AnchorError(f"{process_id}: L023 trident budget must be positive")
    return TridentBudgetBand(
        source=(
            "flavor_catalog/processes/charged_lepton/L023.yaml "
            "measured_observables[CCFR] + "
            f"{_CCFR_CL_VALUE_ID}"
        ),
        active_experiment=active.experiment,
        confidence_level_percent=float(cl_anchor.value_percent),
        gaussian_z=float(gaussian_z),
        sm_ratio=1.0,
        sm_pull=float((1.0 - active.value) / upper_budget),
        experimental_sigma_upper=float(active.uncertainty_upper),
        experimental_sigma_lower=float(active.uncertainty_lower),
        hard_veto_budget_upper=float(upper_budget),
        hard_veto_budget_lower=float(lower_budget),
        lower_edge=float(active.value - lower_budget),
        upper_edge=float(active.value + upper_budget),
    )


def _load_l023_anchor(process_id: str) -> L023Anchor:
    measurements = _load_measurements(process_id)
    active = _find_measurement(
        measurements,
        _ACTIVE_EXPERIMENT,
        process_id=process_id,
    )
    cl_anchor = _load_ccfr_cl_anchor(process_id)
    return L023Anchor(
        measurements=measurements,
        active_measurement=active,
        ccfr_exclusion_cl=cl_anchor,
        budget_band=_build_budget_band(
            active=active,
            cl_anchor=cl_anchor,
            process_id=process_id,
        ),
    )


def _selected_budget(
    predicted: float,
    anchor: L023Anchor,
) -> tuple[float, float, bool]:
    pull = float(predicted - anchor.value)
    budget = (
        anchor.budget_band.hard_veto_budget_upper
        if pull >= 0.0
        else anchor.budget_band.hard_veto_budget_lower
    )
    ratio = abs(pull) / budget if budget > 0.0 else float("inf")
    return float(budget), float(ratio), bool(ratio <= 1.0)


def _measurement_diagnostics(
    anchor: L023Anchor,
) -> tuple[dict[str, Any], ...]:
    return tuple(
        {
            "block_key": measurement.block_key,
            "experiment": measurement.experiment,
            "value": float(measurement.value),
            "uncertainty_upper": float(measurement.uncertainty_upper),
            "uncertainty_lower": float(measurement.uncertainty_lower),
            "uncertainty_type": measurement.uncertainty_type,
            "source_url": measurement.source_url,
            "primary_source_url": measurement.primary_source_url,
        }
        for measurement in anchor.measurements
    )


@register
class Constraint:
    """Catalogued neutrino-trident ratio constraint."""

    process_id = "L023"
    severity = Severity.HARD
    observable = "sigma(nu_mu N -> nu_mu N mu+ mu-) / sigma_SM"

    def __init__(self) -> None:
        self.anchor = _load_l023_anchor(self.process_id)

    def _base_diagnostics(self) -> dict[str, Any]:
        charmii = _find_measurement(
            self.anchor.measurements,
            "CHARM-II",
            process_id=self.process_id,
        )
        return {
            "active_experiment": self.anchor.active_measurement.experiment,
            "active_measurement_block": self.anchor.active_measurement.block_key,
            "active_measurement_source": self.anchor.active_measurement.source_url,
            "all_measurements": _measurement_diagnostics(self.anchor),
            "charmii_ratio": float(charmii.value),
            "charmii_uncertainty": float(charmii.uncertainty),
            "ccfr_exclusion_cl_block": self.anchor.ccfr_exclusion_cl.block_key,
            "ccfr_exclusion_cl_percent": float(
                self.anchor.ccfr_exclusion_cl.value_percent
            ),
            "ccfr_exclusion_cl_source": self.anchor.ccfr_exclusion_cl.source_url,
            "budget_source": self.anchor.budget_band.source,
            "budget_confidence_level_percent": float(
                self.anchor.budget_band.confidence_level_percent
            ),
            "budget_gaussian_z": float(self.anchor.budget_band.gaussian_z),
            "budget_upper_edge": float(self.anchor.budget_band.upper_edge),
            "budget_lower_edge": float(self.anchor.budget_band.lower_edge),
            "budget_hard_veto_upper": float(
                self.anchor.budget_band.hard_veto_budget_upper
            ),
            "budget_hard_veto_lower": float(
                self.anchor.budget_band.hard_veto_budget_lower
            ),
            "sm_ratio_by_definition": 1.0,
            "sm_pull_to_ccfr_budget": float(self.anchor.budget_band.sm_pull),
            "parametrization_citation": TRIDENT_PARAMETRIZATION_CITATION,
            "needs_human_physics": TRIDENT_PROXY_ASSUMPTION_V1,
        }

    def _unevaluated_result(
        self,
        *,
        diagnostics: Mapping[str, Any],
    ) -> ConstraintResult:
        merged = {
            "evaluated": False,
            "unevaluated_reason": _UNEVALUATED_REASON,
            "passes_semantics": (
                "non-vetoing only; no trident NP proxy was evaluated"
            ),
            **self._base_diagnostics(),
            **dict(diagnostics),
        }
        return ConstraintResult(
            process_id=self.process_id,
            severity=self.severity,
            passes=True,
            predicted=None,
            sm_prediction=1.0,
            experimental=float(self.anchor.value),
            ratio=None,
            budget=float(self.anchor.budget),
            notes=_UNEVALUATED_NOTES,
            diagnostics=merged,
        )

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        lepton_input = point.get_extra(_REQUIRED_EXTRA)
        if lepton_input is None:
            return self._unevaluated_result(
                diagnostics={"missing_extra": _REQUIRED_EXTRA},
            )

        kk_ew_mass = point.get_extra(_OPTIONAL_EW_MASS_EXTRA)
        try:
            result = trident_sigma_ratio_from_lepton_input(
                lepton_input,
                m_kk_gev=None if kk_ew_mass is None else float(kk_ew_mass),
            )
        except (AttributeError, KeyError, TypeError, ValueError) as exc:
            return self._unevaluated_result(
                diagnostics={
                    "invalid_extra": _REQUIRED_EXTRA,
                    "exception_type": type(exc).__name__,
                    "exception": str(exc),
                },
            )

        predicted = float(result.sigma_ratio)
        budget, ratio, passes = _selected_budget(predicted, self.anchor)
        diagnostics = dict(result.diagnostics)
        diagnostics.update(
            {
                "evaluated": True,
                **self._base_diagnostics(),
                "kk_ew_mass_extra_used": kk_ew_mass is not None,
                "sin2_theta_w": float(result.sin2_theta_w),
                "c_vector_sm": float(result.c_vector_sm),
                "c_axial_sm": float(result.c_axial_sm),
                "c_vector_total": float(result.c_vector_total),
                "c_axial_total": float(result.c_axial_total),
                "delta_c_vector": float(result.delta_c_vector),
                "delta_c_axial": float(result.delta_c_axial),
                "pull_to_active_measurement": float(predicted - self.anchor.value),
            }
        )

        return ConstraintResult(
            process_id=self.process_id,
            severity=self.severity,
            passes=passes,
            predicted=predicted,
            sm_prediction=float(result.sm_sigma_ratio),
            experimental=float(self.anchor.value),
            ratio=ratio,
            budget=budget,
            notes=(
                "Neutrino-trident sigma/sigma_SM from the effective C_V/C_A "
                "response, compared with the CCFR measured ratio from "
                "L023.yaml using the catalogued 95% CL reinterpretation as a "
                "Gaussian HARD budget. The RS nu_mu nu_mu mu mu matching is a "
                "documented proxy and is flagged NEEDS-HUMAN-PHYSICS."
            ),
            diagnostics=diagnostics,
        )
