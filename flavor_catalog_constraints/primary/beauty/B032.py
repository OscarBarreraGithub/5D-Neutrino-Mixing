"""B032 - charmless nonleptonic ``Bbar -> pi Kbar`` observables.

Physics
-------
The ``B -> pi K`` branching fractions and direct CP asymmetries are classic
charmless nonleptonic observables.  Their SM interpretation requires a
hadronic amplitude framework such as QCDF/SCET or a phenomenological fit of
tree, QCD-penguin, electroweak-penguin, and strong-phase terms.  The "B -> pi
K puzzle" is therefore hadronic-uncertainty limited.  A rigorous RS
prediction would additionally require Delta B = 1 QCD/electroweak penguin
matching, RG evolution, and nonleptonic matrix elements.

This file is an honest non-vetoing stub, like C003.  It loads the HFLAV
End-of-December 2025 ``B -> K pi`` anchors from
``flavor_catalog/processes/beauty/B032.yaml`` and records the direct-CP puzzle
combination

    Delta A_CP(K pi) = A_CP(B+ -> K+ pi0) - A_CP(B0 -> K+ pi-).

No SM or RS hadronic amplitude is computed.

Severity
--------
INFO.  Both the SM side and the RS NP side are flagged
``NEEDS-HUMAN-PHYSICS``.  The returned result is non-vetoing and must not be
used as a hard exclusion.

Catalog sidecar
---------------
``B032.yaml`` stores the relevant HFLAV/PDG values in list-shaped blocks under
``pdg_or_equivalent``.  Selected list entries are virtualized as one-entry
mapping blocks and routed through the scaffold ``load_anchor`` path, then
checked with typed, loud list-entry validation.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Iterable, Mapping

from flavor_catalog_constraints import anchors as anchor_scaffold
from flavor_catalog_constraints.anchors import Anchor, AnchorError, load_anchor, load_full_yaml
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.charmless_b_pik import (
    B_TO_PI_K_NONLEPTONIC_STUB_MODEL_V1,
    B_TO_PI_K_RS_NEEDS_HUMAN_PHYSICS,
    B_TO_PI_K_SM_NEEDS_HUMAN_PHYSICS,
    compare_b_to_pi_k_np_room_to_measurement,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "beauty"
_BR_LIST_KEY = "branching_fractions"
_DIRECT_CP_LIST_KEY = "direct_cp_asymmetries"
_TIME_DEPENDENT_CP_LIST_KEY = "time_dependent_cp"
_POST_2008_LIST_KEY = "post_2008_measurements"

_EXPECTED_BR_UNITS = "branching fraction"
_EXPECTED_DIRECT_CP_UNITS = "percent"
_EXPECTED_DIMENSIONLESS_UNITS = "dimensionless"
_NO_SCALAR_UNCERTAINTY_KEY = "__b032_no_scalar_uncertainty__"

_ACP_CHARGED_PI0_OBSERVABLE = "A_CP(B+ -> K+ pi0)"
_ACP_NEUTRAL_KPLUS_PIMINUS_OBSERVABLE = "A_CP(B0 -> K+ pi-)"
_BELLE_II_SUM_RULE_OBSERVABLE = "Belle II Kpi sum-rule test"

_BUDGET_SOURCE = (
    "flavor_catalog/processes/beauty/B032.yaml direct_cp_asymmetries; "
    "full observed |Delta A_CP(K pi)| is reported as non-vetoing NP room "
    "because no reliable SM subtraction exists"
)
_PARAMETRIZATION_CITATION = (
    "HFLAV Rare B decays End-of-December 2025 B -> K pi branching fractions "
    "and direct CP averages; Delta A_CP(K pi) is a bookkeeping observable, "
    "not a QCDF/SCET amplitude calculation"
)
_NEEDS_HUMAN_PHYSICS = (
    B_TO_PI_K_SM_NEEDS_HUMAN_PHYSICS,
    B_TO_PI_K_RS_NEEDS_HUMAN_PHYSICS,
)


@dataclass(frozen=True)
class BranchingFractionAnchor:
    """Typed HFLAV branching-fraction entry from the B032 list block."""

    block_key: str
    observable: str
    year: int | None
    value: float
    uncertainty: float
    units: str | None
    source_key: str | None
    source_url: str | None
    access_date: str | None
    snapshot_path: str | None
    sha256: str | None


@dataclass(frozen=True)
class DirectCPAnchor:
    """Typed HFLAV direct-CP entry in both percent and dimensionless units."""

    block_key: str
    observable: str
    year: int | None
    value_percent: float
    uncertainty_percent: float
    value: float
    uncertainty: float
    units: str | None
    source_key: str | None
    source_url: str | None
    access_date: str | None
    snapshot_path: str | None
    sha256: str | None


@dataclass(frozen=True)
class TimeDependentCPAnchor:
    """Typed PDG time-dependent CP coefficient for ``B0 -> K0 pi0``."""

    block_key: str
    observable: str
    year: int | None
    value: float
    uncertainty: float
    units: str | None
    source_key: str | None
    source_url: str | None
    access_date: str | None
    snapshot_path: str | None
    sha256: str | None


@dataclass(frozen=True)
class Post2008Measurement:
    """Typed post-2008 companion measurement retained for diagnostics."""

    block_key: str
    observable: str
    year: int | None
    value: float
    uncertainty_text: str | None
    units: str | None
    source_key: str | None
    source_url: str | None
    access_date: str | None
    snapshot_path: str | None
    sha256: str | None


@dataclass(frozen=True)
class DeltaACPAnchor:
    """Direct-CP puzzle combination built from two HFLAV entries."""

    charged_pi0: DirectCPAnchor
    neutral_kplus_piminus: DirectCPAnchor
    value_percent: float
    uncertainty_percent: float
    value: float
    uncertainty: float

    @property
    def budget(self) -> float:
        """Non-vetoing NP room: the full observed ``|Delta A_CP(K pi)|``."""
        return abs(self.value)


@dataclass(frozen=True)
class B032Anchor:
    """Typed B032 anchor bundle and non-vetoing room convention."""

    branching_fractions: tuple[BranchingFractionAnchor, ...]
    direct_cp_asymmetries: tuple[DirectCPAnchor, ...]
    time_dependent_cp: tuple[TimeDependentCPAnchor, ...]
    post_2008_measurements: tuple[Post2008Measurement, ...]
    delta_acp_kpi: DeltaACPAnchor

    @property
    def value(self) -> float:
        """``Delta A_CP(K pi)`` central value as a dimensionless asymmetry."""
        return self.delta_acp_kpi.value

    @property
    def uncertainty(self) -> float:
        """``Delta A_CP(K pi)`` uncertainty as a dimensionless asymmetry."""
        return self.delta_acp_kpi.uncertainty

    @property
    def budget(self) -> float:
        """Non-vetoing NP room for ``Delta A_CP(K pi)``."""
        return self.delta_acp_kpi.budget


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
            f"{process_id}: B032 anchor field {field_name!r}={value!r} "
            "is not numeric"
        ) from exc
    if not math.isfinite(number):
        raise AnchorError(
            f"{process_id}: B032 anchor field {field_name!r}={value!r} "
            "is not finite"
        )
    return number


def _positive_float(value: Any, *, process_id: str, field_name: str) -> float:
    number = _required_float(value, process_id=process_id, field_name=field_name)
    if number <= 0.0:
        raise AnchorError(f"{process_id}: B032 anchor field {field_name!r} <= 0")
    return number


def _entry_observable(
    entry: Mapping[str, Any],
    *,
    process_id: str,
    list_key: str,
    index: int,
) -> str:
    observable = entry.get("observable")
    if not isinstance(observable, str) or not observable:
        raise AnchorError(
            f"{process_id}: {list_key}[{index}] has no non-empty 'observable'"
        )
    return observable


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
) -> Anchor:
    observable = _entry_observable(
        entry,
        process_id=process_id,
        list_key=list_key,
        index=index,
    )
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
            uncertainty_key=_NO_SCALAR_UNCERTAINTY_KEY
            if list_key == _POST_2008_LIST_KEY
            else "uncertainty",
        )
    finally:
        anchor_scaffold.load_pdg_block = original_load_pdg_block

    if scaffold_anchor.block_key != block_key:
        raise AnchorError(
            f"{process_id}: load_anchor selected {scaffold_anchor.block_key!r}, "
            f"expected {block_key!r} for B032 {list_key} observable {observable!r}"
        )
    return scaffold_anchor


def _load_branching_fraction_entries(
    pdg_block: Mapping[str, Any],
    *,
    process_id: str,
) -> tuple[BranchingFractionAnchor, ...]:
    out: list[BranchingFractionAnchor] = []
    for index, entry in enumerate(
        _list_entries(pdg_block, _BR_LIST_KEY, process_id=process_id)
    ):
        anchor = _scaffold_anchor_for_entry(
            entry,
            process_id=process_id,
            list_key=_BR_LIST_KEY,
            index=index,
        )
        if anchor.units != _EXPECTED_BR_UNITS:
            raise AnchorError(
                f"{process_id}: {anchor.block_key} expected units "
                f"{_EXPECTED_BR_UNITS!r}, got {anchor.units!r}"
            )
        if anchor.uncertainty is None:
            raise AnchorError(f"{process_id}: {anchor.block_key} lacks uncertainty")
        out.append(
            BranchingFractionAnchor(
                block_key=anchor.block_key,
                observable=str(anchor.observable),
                year=anchor.year,
                value=_positive_float(
                    anchor.value,
                    process_id=process_id,
                    field_name=f"{anchor.block_key}.value",
                ),
                uncertainty=_positive_float(
                    anchor.uncertainty,
                    process_id=process_id,
                    field_name=f"{anchor.block_key}.uncertainty",
                ),
                units=anchor.units,
                source_key=_optional_str(entry.get("source_key")),
                source_url=_optional_str(anchor.source_url),
                access_date=_optional_str(entry.get("access_date")),
                snapshot_path=_optional_str(anchor.snapshot_path),
                sha256=_optional_str(entry.get("sha256")),
            )
        )
    return tuple(out)


def _load_direct_cp_entries(
    pdg_block: Mapping[str, Any],
    *,
    process_id: str,
) -> tuple[DirectCPAnchor, ...]:
    out: list[DirectCPAnchor] = []
    for index, entry in enumerate(
        _list_entries(pdg_block, _DIRECT_CP_LIST_KEY, process_id=process_id)
    ):
        anchor = _scaffold_anchor_for_entry(
            entry,
            process_id=process_id,
            list_key=_DIRECT_CP_LIST_KEY,
            index=index,
        )
        if anchor.units != _EXPECTED_DIRECT_CP_UNITS:
            raise AnchorError(
                f"{process_id}: {anchor.block_key} expected units "
                f"{_EXPECTED_DIRECT_CP_UNITS!r}, got {anchor.units!r}"
            )
        if anchor.uncertainty is None:
            raise AnchorError(f"{process_id}: {anchor.block_key} lacks uncertainty")
        value_percent = _required_float(
            entry.get("value_percent"),
            process_id=process_id,
            field_name=f"{anchor.block_key}.value_percent",
        )
        uncertainty_percent = _positive_float(
            entry.get("uncertainty_percent"),
            process_id=process_id,
            field_name=f"{anchor.block_key}.uncertainty_percent",
        )
        if not math.isclose(value_percent, anchor.value, rel_tol=0.0, abs_tol=1e-15):
            raise AnchorError(f"{process_id}: {anchor.block_key} value_percent mismatch")
        if not math.isclose(
            uncertainty_percent,
            anchor.uncertainty,
            rel_tol=0.0,
            abs_tol=1e-15,
        ):
            raise AnchorError(
                f"{process_id}: {anchor.block_key} uncertainty_percent mismatch"
            )
        out.append(
            DirectCPAnchor(
                block_key=anchor.block_key,
                observable=str(anchor.observable),
                year=anchor.year,
                value_percent=float(value_percent),
                uncertainty_percent=float(uncertainty_percent),
                value=float(value_percent / 100.0),
                uncertainty=float(uncertainty_percent / 100.0),
                units=anchor.units,
                source_key=_optional_str(entry.get("source_key")),
                source_url=_optional_str(anchor.source_url),
                access_date=_optional_str(entry.get("access_date")),
                snapshot_path=_optional_str(anchor.snapshot_path),
                sha256=_optional_str(entry.get("sha256")),
            )
        )
    return tuple(out)


def _load_time_dependent_cp_entries(
    pdg_block: Mapping[str, Any],
    *,
    process_id: str,
) -> tuple[TimeDependentCPAnchor, ...]:
    out: list[TimeDependentCPAnchor] = []
    for index, entry in enumerate(
        _list_entries(pdg_block, _TIME_DEPENDENT_CP_LIST_KEY, process_id=process_id)
    ):
        anchor = _scaffold_anchor_for_entry(
            entry,
            process_id=process_id,
            list_key=_TIME_DEPENDENT_CP_LIST_KEY,
            index=index,
        )
        if anchor.units != _EXPECTED_DIMENSIONLESS_UNITS:
            raise AnchorError(
                f"{process_id}: {anchor.block_key} expected units "
                f"{_EXPECTED_DIMENSIONLESS_UNITS!r}, got {anchor.units!r}"
            )
        if anchor.uncertainty is None:
            raise AnchorError(f"{process_id}: {anchor.block_key} lacks uncertainty")
        out.append(
            TimeDependentCPAnchor(
                block_key=anchor.block_key,
                observable=str(anchor.observable),
                year=anchor.year,
                value=float(anchor.value),
                uncertainty=_positive_float(
                    anchor.uncertainty,
                    process_id=process_id,
                    field_name=f"{anchor.block_key}.uncertainty",
                ),
                units=anchor.units,
                source_key=_optional_str(entry.get("source_key")),
                source_url=_optional_str(anchor.source_url),
                access_date=_optional_str(entry.get("access_date")),
                snapshot_path=_optional_str(anchor.snapshot_path),
                sha256=_optional_str(entry.get("sha256")),
            )
        )
    return tuple(out)


def _load_post_2008_entries(
    pdg_block: Mapping[str, Any],
    *,
    process_id: str,
) -> tuple[Post2008Measurement, ...]:
    out: list[Post2008Measurement] = []
    for index, entry in enumerate(
        _list_entries(pdg_block, _POST_2008_LIST_KEY, process_id=process_id)
    ):
        anchor = _scaffold_anchor_for_entry(
            entry,
            process_id=process_id,
            list_key=_POST_2008_LIST_KEY,
            index=index,
        )
        if anchor.units != _EXPECTED_DIMENSIONLESS_UNITS:
            raise AnchorError(
                f"{process_id}: {anchor.block_key} expected units "
                f"{_EXPECTED_DIMENSIONLESS_UNITS!r}, got {anchor.units!r}"
            )
        out.append(
            Post2008Measurement(
                block_key=anchor.block_key,
                observable=str(anchor.observable),
                year=anchor.year,
                value=float(anchor.value),
                uncertainty_text=_optional_str(entry.get("uncertainty")),
                units=anchor.units,
                source_key=_optional_str(entry.get("source_key")),
                source_url=_optional_str(anchor.source_url),
                access_date=_optional_str(entry.get("access_date")),
                snapshot_path=_optional_str(anchor.snapshot_path),
                sha256=_optional_str(entry.get("sha256")),
            )
        )
    return tuple(out)


def _find_observable_anchor(
    anchors: Iterable[Any],
    observable: str,
    *,
    process_id: str,
    block_name: str,
) -> Any:
    for anchor in anchors:
        if getattr(anchor, "observable", None) == observable:
            return anchor
    raise AnchorError(
        f"{process_id}: {block_name} has no observable {observable!r}"
    )


def _build_delta_acp_anchor(
    direct_cp_asymmetries: tuple[DirectCPAnchor, ...],
    *,
    process_id: str,
) -> DeltaACPAnchor:
    charged = _find_observable_anchor(
        direct_cp_asymmetries,
        _ACP_CHARGED_PI0_OBSERVABLE,
        process_id=process_id,
        block_name=_DIRECT_CP_LIST_KEY,
    )
    neutral = _find_observable_anchor(
        direct_cp_asymmetries,
        _ACP_NEUTRAL_KPLUS_PIMINUS_OBSERVABLE,
        process_id=process_id,
        block_name=_DIRECT_CP_LIST_KEY,
    )
    value_percent = charged.value_percent - neutral.value_percent
    uncertainty_percent = math.sqrt(
        charged.uncertainty_percent**2 + neutral.uncertainty_percent**2
    )
    if uncertainty_percent <= 0.0:
        raise AnchorError(f"{process_id}: Delta A_CP(K pi) uncertainty <= 0")
    return DeltaACPAnchor(
        charged_pi0=charged,
        neutral_kplus_piminus=neutral,
        value_percent=float(value_percent),
        uncertainty_percent=float(uncertainty_percent),
        value=float(value_percent / 100.0),
        uncertainty=float(uncertainty_percent / 100.0),
    )


def _load_b032_anchor(process_id: str) -> B032Anchor:
    pdg_block = _pdg_block(process_id)
    direct_cp = _load_direct_cp_entries(pdg_block, process_id=process_id)
    anchor = B032Anchor(
        branching_fractions=_load_branching_fraction_entries(
            pdg_block,
            process_id=process_id,
        ),
        direct_cp_asymmetries=direct_cp,
        time_dependent_cp=_load_time_dependent_cp_entries(
            pdg_block,
            process_id=process_id,
        ),
        post_2008_measurements=_load_post_2008_entries(
            pdg_block,
            process_id=process_id,
        ),
        delta_acp_kpi=_build_delta_acp_anchor(
            direct_cp,
            process_id=process_id,
        ),
    )
    if anchor.budget <= 0.0:
        raise AnchorError(f"{process_id}: B032 non-vetoing NP room must be positive")
    return anchor


def _float_map(anchors: Iterable[Any]) -> dict[str, dict[str, float]]:
    return {
        anchor.observable: {
            "value": float(anchor.value),
            "uncertainty": float(anchor.uncertainty),
        }
        for anchor in anchors
    }


@register
class Constraint:
    """Catalogued non-vetoing charmless ``B -> pi K`` stub (process_id B032)."""

    process_id = "B032"
    severity = Severity.INFO
    observable = "Delta A_CP(B -> K pi)"

    def __init__(self) -> None:
        self.anchor = _load_b032_anchor(self.process_id)

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        del point
        comparison = compare_b_to_pi_k_np_room_to_measurement(
            measured_observable=self.anchor.value,
            experimental_uncertainty=self.anchor.uncertainty,
            documented_np_room_abs=self.anchor.budget,
        )
        sum_rule = _find_observable_anchor(
            self.anchor.post_2008_measurements,
            _BELLE_II_SUM_RULE_OBSERVABLE,
            process_id=self.process_id,
            block_name=_POST_2008_LIST_KEY,
        )

        return ConstraintResult(
            process_id=self.process_id,
            severity=self.severity,
            passes=bool(comparison.passes),
            predicted=None,
            sm_prediction=None,
            experimental=float(self.anchor.value),
            ratio=float(comparison.ratio_to_room),
            budget=float(comparison.documented_np_room_abs),
            notes=(
                "INFO-only B032 stub: loads HFLAV B -> K pi direct-CP anchors "
                "and reports the full observed |Delta A_CP(K pi)| as "
                "documented NP room. No QCDF/SCET SM amplitude or RS Delta B=1 "
                "penguin amplitude is computed; both sides are flagged "
                "NEEDS-HUMAN-PHYSICS, so this result is non-vetoing."
            ),
            diagnostics={
                "non_vetoing": True,
                "no_hadronic_amplitude": True,
                "no_penguin_calculation": True,
                "stub_model": B_TO_PI_K_NONLEPTONIC_STUB_MODEL_V1,
                "budget_source": _BUDGET_SOURCE,
                "budget_is_full_observed_delta_acp": True,
                "parametrization_citation": _PARAMETRIZATION_CITATION,
                "needs_human_physics": _NEEDS_HUMAN_PHYSICS,
                "needs_human_physics_sm": B_TO_PI_K_SM_NEEDS_HUMAN_PHYSICS,
                "needs_human_physics_np": B_TO_PI_K_RS_NEEDS_HUMAN_PHYSICS,
                "sm_qcdf_scet_hadronic_limited": True,
                "rs_penguin_matching_available": False,
                "parameter_point_inputs_used": (),
                "experimental_observable": "Delta A_CP(K pi)",
                "delta_acp_value_percent": float(
                    self.anchor.delta_acp_kpi.value_percent
                ),
                "delta_acp_uncertainty_percent": float(
                    self.anchor.delta_acp_kpi.uncertainty_percent
                ),
                "delta_acp_value": float(self.anchor.value),
                "delta_acp_uncertainty": float(self.anchor.uncertainty),
                "measurement_abs": float(comparison.measurement_abs),
                "documented_np_room_abs": float(comparison.documented_np_room_abs),
                "measurement_to_np_room_ratio": float(comparison.ratio_to_room),
                "charged_pi0_acp_percent": float(
                    self.anchor.delta_acp_kpi.charged_pi0.value_percent
                ),
                "charged_pi0_acp_uncertainty_percent": float(
                    self.anchor.delta_acp_kpi.charged_pi0.uncertainty_percent
                ),
                "neutral_kplus_piminus_acp_percent": float(
                    self.anchor.delta_acp_kpi.neutral_kplus_piminus.value_percent
                ),
                "neutral_kplus_piminus_acp_uncertainty_percent": float(
                    self.anchor.delta_acp_kpi.neutral_kplus_piminus
                    .uncertainty_percent
                ),
                "branching_fractions": _float_map(self.anchor.branching_fractions),
                "direct_cp_asymmetries": _float_map(self.anchor.direct_cp_asymmetries),
                "time_dependent_cp": _float_map(self.anchor.time_dependent_cp),
                "belle_ii_kpi_sum_rule_value": float(sum_rule.value),
                "belle_ii_kpi_sum_rule_uncertainty_text": sum_rule.uncertainty_text,
            },
        )
