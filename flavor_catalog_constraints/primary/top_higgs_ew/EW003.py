"""EW003 - inclusive-vs-exclusive ``|V_cb|`` and ``|V_ub|`` tensions.

Physics
-------
EW003 is an observable-level consistency diagnostic, not an RS new-physics
amplitude calculation.  It compares the inclusive and exclusive semileptonic
determinations in ``EW003.yaml`` with the direction-aware Gaussian pull

    pull = |Vxb_inc - Vxb_exc| / sqrt(sigma_inc^2 + sigma_exc^2),

using asymmetric uncertainty components when the YAML provides them.  The
statistical machinery lives in
``flavor_catalog_constraints.physics_adapters.semileptonic_ckm``.

For points carrying ``rs_charged_current``, EW003 also reports the minimal-LH
``epsilon_cb`` and ``epsilon_ub`` charged-current diagnostics from the shared
adapter.  These diagnostics do not alter the scalar pull: a common vector
rescaling ``|1 + epsilon|`` cancels in an inclusive/exclusive ratio for the
same CKM element.  A differential treatment is only possible with a supplied
covariance/scheme model.

NEEDS-HUMAN-PHYSICS
-------------------
Rigorous RS matching for this observable would require inclusive-scheme,
exclusive form-factor, and CKM-fit covariance inputs that are not available on
``ParameterPoint``.  No ungrounded RS proxy is applied; the result reports the
data tension and flags this matching gap in diagnostics.

Severity
--------
SOFT.  The current inclusive/exclusive spread is a theory/experimental
semileptonic tension, not an observed NP upper bound or veto.

Catalog sidecar
---------------
``flavor_catalog/processes/top_higgs_ew/EW003.yaml`` is the source of truth
for the PDG/HFLAV/FLAG values and provenance.  Its ``pdg_or_equivalent`` block
is a list, so this module adapts selected list entries into the scaffold
``load_anchor`` path rather than reading value-bearing anchors ad hoc.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Mapping, Sequence

from flavor_catalog_constraints import anchors as anchor_scaffold
from flavor_catalog_constraints.anchors import Anchor, AnchorError, load_anchor, load_full_yaml
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.charged_current import (
    EW003_CKM_TENSION_COVARIANCE_NEEDS_HUMAN,
    EW003_CKM_TENSION_DIAGNOSTIC_STATUS,
    charged_current_ckm_tension_diagnostics,
    charged_current_source_diagnostics,
)
from flavor_catalog_constraints.physics_adapters.semileptonic_ckm import (
    CKMDetermination,
    InclusiveExclusivePull,
    inclusive_exclusive_pull,
    summarize_ckm_tensions,
    uncertainty_from_components,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "top_higgs_ew"
_EXPECTED_CKM_UNITS = "10^-3"
_EXPECTED_BUDGET_UNITS = "sigma"
_PRIMARY_BUDGET_SOURCE = (
    "flavor_catalog/processes/top_higgs_ew/EW003.yaml: "
    "PDG 2024 |V_cb| inclusive-exclusive marginal consistency"
)
_OPTIONAL_CHARGED_CURRENT_EXTRA = "rs_charged_current"

_PDG_VCB_INCLUSIVE = "PDG 2024 |V_cb| inclusive"
_PDG_VCB_EXCLUSIVE = "PDG 2024 |V_cb| exclusive"
_PDG_VCB_AVERAGE = "PDG 2024 |V_cb| scaled average"
_PDG_VUB_INCLUSIVE = "PDG 2024 |V_ub| inclusive"
_PDG_VUB_EXCLUSIVE = "PDG 2024 |V_ub| exclusive"
_PDG_VUB_AVERAGE = "PDG 2024 |V_ub| scaled average"
_PDG_VCB_CONSISTENCY = "PDG 2024 |V_cb| inclusive-exclusive marginal consistency"
_HFLAV_VCB_INCLUSIVE = "HFLAV 2023 preferred inclusive |V_cb|"
_HFLAV_VCB_EXCLUSIVE = "HFLAV 2023 combined exclusive |V_cb|"
_HFLAV_VUB_EXCLUSIVE = "HFLAV 2023 combined exclusive |V_ub|"
_FLAG_VCB_EXCLUSIVE = "FLAG 2024 B -> D* l nu |V_cb|, Nf=2+1"
_FLAG_VUB_EXCLUSIVE = "FLAG 2024 B -> pi l nu |V_ub|, Nf=2+1"


@dataclass(frozen=True)
class EW003ValueAnchor:
    """Typed view over one EW003 value entry."""

    block_key: str
    entry_index: int
    observable: str
    value: float
    uncertainty_upper: float | None
    uncertainty_lower: float | None
    uncertainty_components: Mapping[str, float]
    units: str | None
    source: str | None
    source_url: str | None
    year: int | None
    snapshot_path: str | None

    @property
    def uncertainty(self) -> float | None:
        if self.uncertainty_upper is None or self.uncertainty_lower is None:
            return None
        return float(0.5 * (self.uncertainty_upper + self.uncertainty_lower))

    def as_determination(self) -> CKMDetermination:
        if self.uncertainty_upper is None or self.uncertainty_lower is None:
            raise AnchorError(
                f"EW003: anchor {self.observable!r} has no usable uncertainty"
            )
        return CKMDetermination(
            label=self.observable,
            value=float(self.value),
            uncertainty_upper=float(self.uncertainty_upper),
            uncertainty_lower=float(self.uncertainty_lower),
            units=self.units,
            source=self.source,
        )


@dataclass(frozen=True)
class EW003Anchor:
    """All YAML-loaded EW003 anchors used by the constraint."""

    pdg_vcb_inclusive: EW003ValueAnchor
    pdg_vcb_exclusive: EW003ValueAnchor
    pdg_vcb_scaled_average: EW003ValueAnchor
    pdg_vub_inclusive: EW003ValueAnchor
    pdg_vub_exclusive: EW003ValueAnchor
    pdg_vub_scaled_average: EW003ValueAnchor
    pdg_vcb_consistency_sigma: EW003ValueAnchor
    hflav_vcb_inclusive: EW003ValueAnchor
    hflav_vcb_exclusive: EW003ValueAnchor
    hflav_vub_exclusive: EW003ValueAnchor
    flag_vcb_exclusive: EW003ValueAnchor
    flag_vub_exclusive: EW003ValueAnchor

    @property
    def budget(self) -> float:
        return float(self.pdg_vcb_consistency_sigma.value)

    @property
    def value(self) -> float:
        return self.budget

    @property
    def source_url(self) -> str | None:
        return self.pdg_vcb_consistency_sigma.source_url


def _required_float(value: Any, *, process_id: str, field_name: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise AnchorError(
            f"{process_id}: EW003 field {field_name!r}={value!r} is not numeric"
        ) from exc
    if not math.isfinite(number):
        raise AnchorError(
            f"{process_id}: EW003 field {field_name!r}={value!r} is not finite"
        )
    return number


def _optional_str(value: Any) -> str | None:
    return None if value is None else str(value)


def _pdg_entries(process_id: str) -> Sequence[Mapping[str, Any]]:
    data = load_full_yaml(process_id, family=_FAMILY)
    entries = data.get("pdg_or_equivalent")
    if not isinstance(entries, list) or not entries:
        raise AnchorError(f"{process_id}: expected non-empty pdg_or_equivalent list")
    for index, entry in enumerate(entries):
        if not isinstance(entry, Mapping):
            raise AnchorError(
                f"{process_id}: pdg_or_equivalent[{index}] is not a mapping"
            )
    return entries


def _find_entry(
    process_id: str,
    observable: str,
) -> tuple[int, Mapping[str, Any]]:
    for index, entry in enumerate(_pdg_entries(process_id)):
        if entry.get("observable") == observable:
            return index, entry
    present = [str(entry.get("observable")) for entry in _pdg_entries(process_id)]
    raise AnchorError(
        f"{process_id}: observable {observable!r} not found in "
        f"pdg_or_equivalent list (present: {present})"
    )


def _load_scaffold_list_anchor(
    observable: str,
    *,
    process_id: str,
) -> tuple[Anchor, Mapping[str, Any], int]:
    index, entry = _find_entry(process_id, observable)
    block_key = f"pdg_or_equivalent[{index}]"
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
        )
    finally:
        anchor_scaffold.load_pdg_block = original_load_pdg_block

    return scaffold_anchor, entry, index


def _components_mapping(
    value: Any,
    *,
    process_id: str,
    field_name: str,
) -> Mapping[str, float]:
    if value is None:
        return {}
    if not isinstance(value, Mapping):
        raise AnchorError(f"{process_id}: {field_name} is not a mapping")
    return {
        str(key): _required_float(
            component,
            process_id=process_id,
            field_name=f"{field_name}.{key}",
        )
        for key, component in value.items()
    }


def _load_value_anchor(
    observable: str,
    *,
    process_id: str,
    expected_units: str,
    require_uncertainty: bool = True,
) -> EW003ValueAnchor:
    scaffold_anchor, entry, index = _load_scaffold_list_anchor(
        observable,
        process_id=process_id,
    )
    if scaffold_anchor.units != expected_units:
        raise AnchorError(
            f"{process_id}: expected units {expected_units!r} for {observable}, "
            f"got {scaffold_anchor.units!r}"
        )
    components = _components_mapping(
        entry.get("uncertainty_components"),
        process_id=process_id,
        field_name=f"pdg_or_equivalent[{index}].uncertainty_components",
    )
    uncertainty = None
    if scaffold_anchor.uncertainty is not None:
        uncertainty = uncertainty_from_components(symmetric=scaffold_anchor.uncertainty)
    elif components:
        uncertainty = uncertainty_from_components(components)
    elif require_uncertainty:
        raise AnchorError(f"{process_id}: {observable} has no uncertainty")

    return EW003ValueAnchor(
        block_key=scaffold_anchor.block_key,
        entry_index=index,
        observable=str(entry.get("observable")),
        value=_required_float(
            scaffold_anchor.value,
            process_id=process_id,
            field_name=f"pdg_or_equivalent[{index}].value",
        ),
        uncertainty_upper=None if uncertainty is None else float(uncertainty.upper),
        uncertainty_lower=None if uncertainty is None else float(uncertainty.lower),
        uncertainty_components={} if uncertainty is None else dict(uncertainty.components),
        units=scaffold_anchor.units,
        source=_optional_str(scaffold_anchor.source),
        source_url=_optional_str(scaffold_anchor.source_url),
        year=scaffold_anchor.year,
        snapshot_path=_optional_str(scaffold_anchor.snapshot_path),
    )


def _load_ew003_anchor(process_id: str) -> EW003Anchor:
    anchor = EW003Anchor(
        pdg_vcb_inclusive=_load_value_anchor(
            _PDG_VCB_INCLUSIVE,
            process_id=process_id,
            expected_units=_EXPECTED_CKM_UNITS,
        ),
        pdg_vcb_exclusive=_load_value_anchor(
            _PDG_VCB_EXCLUSIVE,
            process_id=process_id,
            expected_units=_EXPECTED_CKM_UNITS,
        ),
        pdg_vcb_scaled_average=_load_value_anchor(
            _PDG_VCB_AVERAGE,
            process_id=process_id,
            expected_units=_EXPECTED_CKM_UNITS,
        ),
        pdg_vub_inclusive=_load_value_anchor(
            _PDG_VUB_INCLUSIVE,
            process_id=process_id,
            expected_units=_EXPECTED_CKM_UNITS,
        ),
        pdg_vub_exclusive=_load_value_anchor(
            _PDG_VUB_EXCLUSIVE,
            process_id=process_id,
            expected_units=_EXPECTED_CKM_UNITS,
        ),
        pdg_vub_scaled_average=_load_value_anchor(
            _PDG_VUB_AVERAGE,
            process_id=process_id,
            expected_units=_EXPECTED_CKM_UNITS,
        ),
        pdg_vcb_consistency_sigma=_load_value_anchor(
            _PDG_VCB_CONSISTENCY,
            process_id=process_id,
            expected_units=_EXPECTED_BUDGET_UNITS,
            require_uncertainty=False,
        ),
        hflav_vcb_inclusive=_load_value_anchor(
            _HFLAV_VCB_INCLUSIVE,
            process_id=process_id,
            expected_units=_EXPECTED_CKM_UNITS,
        ),
        hflav_vcb_exclusive=_load_value_anchor(
            _HFLAV_VCB_EXCLUSIVE,
            process_id=process_id,
            expected_units=_EXPECTED_CKM_UNITS,
        ),
        hflav_vub_exclusive=_load_value_anchor(
            _HFLAV_VUB_EXCLUSIVE,
            process_id=process_id,
            expected_units=_EXPECTED_CKM_UNITS,
        ),
        flag_vcb_exclusive=_load_value_anchor(
            _FLAG_VCB_EXCLUSIVE,
            process_id=process_id,
            expected_units=_EXPECTED_CKM_UNITS,
        ),
        flag_vub_exclusive=_load_value_anchor(
            _FLAG_VUB_EXCLUSIVE,
            process_id=process_id,
            expected_units=_EXPECTED_CKM_UNITS,
        ),
    )
    if anchor.budget <= 0.0:
        raise AnchorError(f"{process_id}: EW003 sigma budget must be positive")
    return anchor


def _pull_diagnostics(pull: InclusiveExclusivePull) -> dict[str, Any]:
    return {
        "inclusive_value": float(pull.inclusive.value),
        "exclusive_value": float(pull.exclusive.value),
        "difference": float(pull.difference),
        "combined_sigma": float(pull.combined_sigma),
        "pull_sigma": float(pull.pull_sigma),
        "inclusive_sigma_used": float(pull.inclusive_sigma_used),
        "exclusive_sigma_used": float(pull.exclusive_sigma_used),
        "budget_sigma": pull.budget_sigma,
        "ratio_to_budget": pull.ratio_to_budget,
        "passes_budget": pull.passes_budget,
        "units": pull.inclusive.units,
        "inclusive_source": pull.inclusive.source,
        "exclusive_source": pull.exclusive.source,
    }


@register
class Constraint:
    """Catalogued semileptonic CKM inclusive/exclusive tension (EW003)."""

    process_id = "EW003"
    severity = Severity.SOFT
    observable = "|V_cb| and |V_ub| inclusive-exclusive tension"

    def __init__(self) -> None:
        self.anchor = _load_ew003_anchor(self.process_id)
        self.pdg_vcb_pull = inclusive_exclusive_pull(
            "PDG 2024 |V_cb|",
            self.anchor.pdg_vcb_inclusive.as_determination(),
            self.anchor.pdg_vcb_exclusive.as_determination(),
            budget_sigma=self.anchor.budget,
        )
        self.pdg_vub_pull = inclusive_exclusive_pull(
            "PDG 2024 |V_ub|",
            self.anchor.pdg_vub_inclusive.as_determination(),
            self.anchor.pdg_vub_exclusive.as_determination(),
            budget_sigma=self.anchor.budget,
        )
        self.summary = summarize_ckm_tensions(
            (self.pdg_vcb_pull, self.pdg_vub_pull),
            budget_sigma=self.anchor.budget,
        )
        self.hflav_vcb_pull = inclusive_exclusive_pull(
            "HFLAV 2023 |V_cb|",
            self.anchor.hflav_vcb_inclusive.as_determination(),
            self.anchor.hflav_vcb_exclusive.as_determination(),
            budget_sigma=self.anchor.budget,
        )

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        diagnostics: dict[str, Any] = {
            "primary_pulls": {
                self.pdg_vcb_pull.observable: _pull_diagnostics(self.pdg_vcb_pull),
                self.pdg_vub_pull.observable: _pull_diagnostics(self.pdg_vub_pull),
            },
            "validation_pulls": {
                self.hflav_vcb_pull.observable: _pull_diagnostics(
                    self.hflav_vcb_pull
                ),
            },
            "max_pull_observable": self.summary.max_pull.observable,
            "pdg_vcb_scaled_average": float(self.anchor.pdg_vcb_scaled_average.value),
            "pdg_vcb_scaled_average_uncertainty": (
                self.anchor.pdg_vcb_scaled_average.uncertainty
            ),
            "pdg_vub_scaled_average": float(self.anchor.pdg_vub_scaled_average.value),
            "pdg_vub_scaled_average_uncertainty": (
                self.anchor.pdg_vub_scaled_average.uncertainty
            ),
            "hflav_vub_exclusive": float(self.anchor.hflav_vub_exclusive.value),
            "flag_vcb_exclusive": float(self.anchor.flag_vcb_exclusive.value),
            "flag_vub_exclusive": float(self.anchor.flag_vub_exclusive.value),
            "budget_source": _PRIMARY_BUDGET_SOURCE,
            "budget_block": self.anchor.pdg_vcb_consistency_sigma.block_key,
            "matching_coverage": "PARTIAL",
            "ew003_charged_current_status": EW003_CKM_TENSION_DIAGNOSTIC_STATUS,
            "ew003_covariance_scheme_status": (
                EW003_CKM_TENSION_COVARIANCE_NEEDS_HUMAN
            ),
            "ew003_covariance_scheme_input_supplied": False,
            "needs_human_physics": EW003_CKM_TENSION_COVARIANCE_NEEDS_HUMAN,
            "required_parameter_point_extras": [],
            "optional_parameter_point_extras": [_OPTIONAL_CHARGED_CURRENT_EXTRA],
            "parameter_point_used": False,
            "parameter_point_used_for_scalar_pull": False,
            "qcd_running_applied": False,
        }
        charged = point.get_extra(_OPTIONAL_CHARGED_CURRENT_EXTRA)
        if charged is None:
            diagnostics.update(
                {
                    "charged_current_diagnostics_available": False,
                    "charged_current_missing_optional_extra": (
                        _OPTIONAL_CHARGED_CURRENT_EXTRA
                    ),
                }
            )
        else:
            try:
                diagnostics["charged_current_diagnostics"] = {
                    **charged_current_ckm_tension_diagnostics(charged),
                    **charged_current_source_diagnostics(charged),
                }
                diagnostics["charged_current_diagnostics_available"] = True
                diagnostics["charged_current_missing_optional_extra"] = None
            except Exception as exc:  # noqa: BLE001 - diagnostics degrade cleanly
                diagnostics.update(
                    {
                        "charged_current_diagnostics_available": False,
                        "charged_current_invalid_optional_extra": (
                            _OPTIONAL_CHARGED_CURRENT_EXTRA
                        ),
                        "charged_current_exception_type": type(exc).__name__,
                        "charged_current_exception_message": str(exc),
                    }
                )

        return ConstraintResult(
            process_id=self.process_id,
            severity=self.severity,
            passes=bool(self.summary.passes),
            predicted=float(self.summary.max_pull.pull_sigma),
            sm_prediction=0.0,
            experimental=float(self.anchor.pdg_vcb_consistency_sigma.value),
            ratio=float(self.summary.ratio_to_budget),
            budget=float(self.summary.budget_sigma),
            notes=(
                "Data-only inclusive/exclusive semileptonic CKM pull. "
                "Primary scalar is max(PDG |V_cb| pull, PDG |V_ub| pull) "
                "compared with the PDG 3.0 sigma marginal-consistency anchor; "
                "SOFT failure flags tension but is not an NP veto. "
                "Minimal-LH charged-current epsilons are diagnostic only "
                "because a universal rescaling cancels without a "
                "covariance/scheme input."
            ),
            diagnostics=diagnostics,
        )
