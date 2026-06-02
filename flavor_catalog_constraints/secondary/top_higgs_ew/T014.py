"""T014 - flavor-changing down-sector Z decays.

Physics
-------
The Standard Model rates for ``Z -> b sbar + bbar s``,
``Z -> b dbar + bbar d``, and ``Z -> s dbar + sbar d`` are loop, GIM, and CKM
suppressed far below the direct non-standard hadronic-Z-width limit.  This
constraint therefore applies a pure-NP HARD upper bound to the largest
charge-summed FCNC branching-fraction proxy,

    BR = Gamma(Z -> q_i qbar_j + q_j qbar_i)
         / (Gamma_Z^{SM,total} + Gamma_FCNC).

The width arithmetic reuses the shared Z-pole effective-coupling convention
from ``quarkConstraints.zpole`` through the
``flavor_catalog_constraints.physics_adapters.zpole`` boundary.

RS matching status
------------------
NEEDS-HUMAN-PHYSICS.  A rigorous RS prediction needs the off-diagonal
down-sector ``Z d_i d_j`` effective coupling after EW KK/Z/Z' mixing,
custodial-embedding choices, brane kinetic terms, and down-quark mass-basis
rotation.  Those inputs are not present on ``ParameterPoint``.  This v1
implementation uses a documented proxy:
``delta g_ij = (m_Z/M_KK)^2 * down_overlap_ij`` from the available quark
mass-basis overlap matrices.

Severity
--------
HARD.  The catalogued rows are observed 95% CL direct upper limits on
non-standard hadronic Z widths.  A point with an evaluated pure-NP proxy above
the relevant limit is excluded.

Catalog sidecar
---------------
``flavor_catalog/processes/secondary/top_higgs_ew/T014.yaml`` is the source of
truth for the three 95% CL branching-fraction limits and provenance.  Numeric
limits are parsed from the YAML upper-limit strings and routed through the
scaffold ``load_anchor`` path; no experimental number is hardcoded here.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
import re
from typing import Any, Mapping, Sequence

from flavor_catalog_constraints import anchors as anchor_scaffold
from flavor_catalog_constraints.anchors import Anchor, AnchorError, load_anchor, load_full_yaml
from flavor_catalog_constraints.base import (
    ConstraintLevel,
    ConstraintResult,
    ParameterPoint,
    Severity,
)
from flavor_catalog_constraints.physics_adapters.zpole import (
    ZPOLE_DOWN_FCNC_PROXY_V1,
    ZPoleDownFCNCBranchingResult,
    ZPoleDownFCNCCouplingProxy,
    zpole_down_fcnc_branching_fraction_from_couplings,
    zpole_down_fcnc_branching_fraction_with_proxy,
    zpole_down_fcnc_effective_coupling_limit,
    zpole_down_fcnc_sm_hadronic_width_weight,
    zpole_down_fcnc_sm_total_width_weight,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "top_higgs_ew"
_TIER = ConstraintLevel.SECONDARY
_REQUIRED_EXTRA = "quark_mass_basis_couplings"
_OPTIONAL_EW_MASS_EXTRA = "kk_ew_mass_gev"

_BS_VALUE_ID = "ECFA2025:T014:zbs_direct_hadronic_width_limit"
_BD_VALUE_ID = "ECFA2025:T014:zbd_direct_hadronic_width_limit"
_SD_VALUE_ID = "AbuAjamieh2026:T014:zsd_direct_hadronic_width_limit"
_CHANNEL_SPECS: tuple[tuple[str, str, str, str], ...] = (
    ("bs", "b", "s", _BS_VALUE_ID),
    ("bd", "b", "d", _BD_VALUE_ID),
    ("sd", "s", "d", _SD_VALUE_ID),
)

_EXPECTED_UNITS = "dimensionless branching fraction"
_EXPECTED_CL = "95% CL"
_EXPECTED_LIMIT_TYPE = "upper_limit_from_inclusive_nonstandard_hadronic_Z_width"
_NUMBER_RE = re.compile(
    r"^\s*<?\s*(?P<number>[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)"
    r"(?:[eE][+-]?[0-9]+)?)\s*$"
)
_UNEVALUATED_REASON = (
    "no off-diagonal down-sector Z coupling prediction available "
    "(EW FCNC-Z matching not on ParameterPoint)"
)
_UNEVALUATED_NOTES = f"NOT EVALUATED - {_UNEVALUATED_REASON}"
_BUDGET_SOURCE = (
    "flavor_catalog/processes/secondary/top_higgs_ew/T014.yaml "
    "pdg_or_equivalent direct non-standard hadronic-Z-width rows"
)
_NEEDS_HUMAN_PHYSICS = (
    "NEEDS-HUMAN-PHYSICS: rigorous RS FCNC-Z matching requires EW KK/Z/Z' "
    "mixing, custodial representations, brane kinetic terms, and the "
    "down-sector Z neutral-current mass-basis coupling not present on "
    "ParameterPoint; this v1 uses the documented down-overlap proxy."
)


@dataclass(frozen=True)
class DownFCNCLimitEntry:
    """Typed view over one T014 direct branching-fraction upper-limit row."""

    anchor: Anchor
    channel_key: str
    flavor_i: str
    flavor_j: str
    value_id: str
    observable: str | None
    confidence_level: str | None
    limit_type: str | None
    effective_coupling_limit: float
    raw_row: Mapping[str, Any]

    @property
    def block_key(self) -> str:
        return self.anchor.block_key

    @property
    def value(self) -> float:
        return self.anchor.value

    @property
    def units(self) -> str | None:
        return self.anchor.units

    @property
    def source_url(self) -> str | None:
        return self.anchor.source_url

    @property
    def snapshot_path(self) -> str | None:
        return self.anchor.snapshot_path


@dataclass(frozen=True)
class T014Anchor:
    """Typed T014 anchor bundle for the three down-sector FCNC channels."""

    channels: Mapping[str, DownFCNCLimitEntry]
    sm_hadronic_width_weight: float
    sm_total_width_weight: float
    sm_hadronic_to_total_width_ratio: float

    @property
    def value(self) -> float:
        return min(entry.value for entry in self.channels.values())

    @property
    def budget(self) -> float:
        return self.value


def _optional_str(value: Any) -> str | None:
    return None if value is None else str(value)


def _optional_int(value: Any) -> int | None:
    if value is None:
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _parse_limit_value(
    value: Any,
    *,
    process_id: str,
    field_name: str,
) -> float:
    if isinstance(value, (int, float)):
        number = float(value)
    elif isinstance(value, str):
        match = _NUMBER_RE.match(value)
        if match is None:
            raise AnchorError(
                f"{process_id}: T014 field {field_name!r}={value!r} is not "
                "a numeric upper-limit string"
            )
        number = float(match.group("number"))
    else:
        raise AnchorError(
            f"{process_id}: T014 field {field_name!r}={value!r} is not numeric"
        )
    if not math.isfinite(number) or number <= 0.0:
        raise AnchorError(f"{process_id}: T014 field {field_name!r} must be positive")
    return number


def _pdg_values(process_id: str) -> Sequence[Mapping[str, Any]]:
    data = load_full_yaml(process_id, family=_FAMILY, tier=_TIER)
    entries = data.get("pdg_or_equivalent")
    if not isinstance(entries, list) or not entries:
        raise AnchorError(f"{process_id}: expected non-empty pdg_or_equivalent list")
    for index, entry in enumerate(entries):
        if not isinstance(entry, Mapping):
            raise AnchorError(
                f"{process_id}: pdg_or_equivalent[{index}] is not a mapping"
            )
    return entries


def _find_value_entry(
    process_id: str,
    value_id: str,
) -> tuple[int, Mapping[str, Any]]:
    for index, entry in enumerate(_pdg_values(process_id)):
        if entry.get("value_id") == value_id:
            return index, entry
    present = [str(entry.get("value_id")) for entry in _pdg_values(process_id)]
    raise AnchorError(
        f"{process_id}: value_id {value_id!r} not found in T014 "
        f"pdg_or_equivalent list (present: {present})"
    )


def _load_scaffold_limit_anchor(
    value_id: str,
    *,
    process_id: str,
) -> Anchor:
    index, entry = _find_value_entry(process_id, value_id)
    block_key = f"pdg_or_equivalent[{index}]"
    numeric_entry = dict(entry)
    numeric_entry["value"] = _parse_limit_value(
        entry.get("normalized_value", entry.get("value")),
        process_id=process_id,
        field_name=f"{value_id}.normalized_value",
    )
    virtual_block = {block_key: numeric_entry}
    original_load_pdg_block = anchor_scaffold.load_pdg_block

    def _load_virtual_pdg_block(
        request_process_id: str,
        **kwargs: Any,
    ) -> Mapping[str, Any]:
        if (
            request_process_id == process_id
            and kwargs.get("family") == _FAMILY
            and kwargs.get("tier", ConstraintLevel.PRIMARY) == _TIER
        ):
            return virtual_block
        return original_load_pdg_block(request_process_id, **kwargs)

    anchor_scaffold.load_pdg_block = _load_virtual_pdg_block
    try:
        scaffold_anchor = load_anchor(
            process_id,
            family=_FAMILY,
            tier=_TIER,
            candidates=(block_key,),
        )
    finally:
        anchor_scaffold.load_pdg_block = original_load_pdg_block

    if scaffold_anchor.block_key != block_key:
        raise AnchorError(
            f"{process_id}: load_anchor selected {scaffold_anchor.block_key!r}, "
            f"expected {block_key!r} for value_id {value_id!r}"
        )
    return scaffold_anchor


def _limit_entry(
    process_id: str,
    *,
    channel_key: str,
    flavor_i: str,
    flavor_j: str,
    value_id: str,
) -> DownFCNCLimitEntry:
    _, row = _find_value_entry(process_id, value_id)
    anchor = _load_scaffold_limit_anchor(value_id, process_id=process_id)
    if anchor.units != _EXPECTED_UNITS:
        raise AnchorError(
            f"{process_id}: T014 {value_id!r} units must be "
            f"{_EXPECTED_UNITS!r}, got {anchor.units!r}"
        )
    confidence_level = _optional_str(row.get("confidence_level", row.get("cl")))
    if confidence_level != _EXPECTED_CL:
        raise AnchorError(
            f"{process_id}: T014 {value_id!r} confidence level must be "
            f"{_EXPECTED_CL!r}, got {confidence_level!r}"
        )
    limit_type = _optional_str(row.get("limit_type"))
    if limit_type != _EXPECTED_LIMIT_TYPE:
        raise AnchorError(
            f"{process_id}: T014 {value_id!r} limit_type must be "
            f"{_EXPECTED_LIMIT_TYPE!r}, got {limit_type!r}"
        )
    if anchor.value <= 0.0 or not math.isfinite(anchor.value):
        raise AnchorError(f"{process_id}: T014 {value_id!r} limit must be positive")
    return DownFCNCLimitEntry(
        anchor=anchor,
        channel_key=channel_key,
        flavor_i=flavor_i,
        flavor_j=flavor_j,
        value_id=value_id,
        observable=_optional_str(row.get("observable")),
        confidence_level=confidence_level,
        limit_type=limit_type,
        effective_coupling_limit=zpole_down_fcnc_effective_coupling_limit(
            anchor.value,
            flavor_i=flavor_i,
            flavor_j=flavor_j,
        ),
        raw_row=row,
    )


def _load_t014_anchor(process_id: str) -> T014Anchor:
    channels = {
        channel_key: _limit_entry(
            process_id,
            channel_key=channel_key,
            flavor_i=flavor_i,
            flavor_j=flavor_j,
            value_id=value_id,
        )
        for channel_key, flavor_i, flavor_j, value_id in _CHANNEL_SPECS
    }
    sm_hadronic = float(sum(zpole_down_fcnc_sm_hadronic_width_weight().values()))
    sm_total = float(sum(zpole_down_fcnc_sm_total_width_weight().values()))
    return T014Anchor(
        channels=channels,
        sm_hadronic_width_weight=sm_hadronic,
        sm_total_width_weight=sm_total,
        sm_hadronic_to_total_width_ratio=float(sm_hadronic / sm_total),
    )


def _channel_diagnostics(
    limit: DownFCNCLimitEntry,
    result: ZPoleDownFCNCBranchingResult,
    proxy: ZPoleDownFCNCCouplingProxy,
) -> dict[str, Any]:
    return {
        "value_id": limit.value_id,
        "observable": limit.observable,
        "experimental_block": limit.block_key,
        "confidence_level": limit.confidence_level,
        "source_url": limit.source_url,
        "snapshot_path": limit.snapshot_path,
        "flavor_i": result.flavor_i,
        "flavor_j": result.flavor_j,
        "predicted": float(result.branching_fraction),
        "sm_prediction": 0.0,
        "budget": float(result.br_limit),
        "ratio": float(result.ratio_to_limit),
        "passes": bool(result.passes),
        "delta_g_left": complex(result.delta_g_left),
        "delta_g_right": complex(result.delta_g_right),
        "effective_coupling_norm": float(result.coupling_norm),
        "effective_coupling_limit": float(limit.effective_coupling_limit),
        "fcnc_width_weight": float(result.fcnc_width_weight),
        "sm_hadronic_width_weight": float(result.sm_hadronic_width_weight),
        "sm_total_width_weight": float(result.sm_total_width_weight),
        "total_hadronic_width_weight": float(result.total_hadronic_width_weight),
        "total_width_weight": float(result.total_width_weight),
        "sm_hadronic_to_total_width_ratio": float(
            result.sm_hadronic_to_total_width_ratio
        ),
        "charge_state_factor": float(result.charge_state_factor),
        "n_color": float(result.n_color),
        "radiator": float(result.radiator),
        "proxy": dict(proxy.diagnostics),
    }


@register
class Constraint:
    """Catalogued direct down-sector FCNC-Z hadronic-width constraint."""

    process_id = "T014"
    severity = Severity.HARD
    observable = "max BR(Z -> bs, bd, sd) down-sector FCNC"

    def __init__(self) -> None:
        self.anchor = _load_t014_anchor(self.process_id)
        self.sm_results = {
            channel_key: zpole_down_fcnc_branching_fraction_from_couplings(
                flavor_i=entry.flavor_i,
                flavor_j=entry.flavor_j,
                delta_g_left=0.0j,
                delta_g_right=0.0j,
                br_limit=entry.value,
            )
            for channel_key, entry in self.anchor.channels.items()
        }

    def _unevaluated_result(self, diagnostics: Mapping[str, Any]) -> ConstraintResult:
        return ConstraintResult(
            process_id=self.process_id,
            severity=self.severity,
            passes=True,
            predicted=None,
            sm_prediction=0.0,
            experimental=float(self.anchor.value),
            ratio=None,
            budget=float(self.anchor.budget),
            notes=_UNEVALUATED_NOTES,
            diagnostics={
                "evaluated": False,
                "unevaluated_reason": _UNEVALUATED_REASON,
                "passes_semantics": (
                    "non-vetoing only; no BR(Z -> down-sector FCNC) NP "
                    "prediction was evaluated"
                ),
                "required_parameter_point_extras": [_REQUIRED_EXTRA],
                "needs_human_physics": _NEEDS_HUMAN_PHYSICS,
                "rs_matching_assumption": ZPOLE_DOWN_FCNC_PROXY_V1,
                "budget_source": _BUDGET_SOURCE,
                "sm_prediction_policy": (
                    "SM FCNC-Z rates are loop/GIM/CKM suppressed and treated "
                    "as zero for this pure-NP direct-width bound."
                ),
                **dict(diagnostics),
            },
        )

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        couplings = point.get_extra(_REQUIRED_EXTRA)
        if couplings is None:
            return self._unevaluated_result({"missing_extra": _REQUIRED_EXTRA})

        kk_ew_mass = point.get_extra(_OPTIONAL_EW_MASS_EXTRA)
        channel_results: dict[str, ZPoleDownFCNCBranchingResult] = {}
        channel_proxies: dict[str, ZPoleDownFCNCCouplingProxy] = {}
        try:
            for channel_key, entry in self.anchor.channels.items():
                result, proxy = zpole_down_fcnc_branching_fraction_with_proxy(
                    couplings,
                    flavor_i=entry.flavor_i,
                    flavor_j=entry.flavor_j,
                    br_limit=entry.value,
                    m_kk_gev=None if kk_ew_mass is None else float(kk_ew_mass),
                )
                channel_results[channel_key] = result
                channel_proxies[channel_key] = proxy
        except (AttributeError, KeyError, TypeError, ValueError) as exc:
            return self._unevaluated_result(
                {
                    "invalid_extra": _REQUIRED_EXTRA,
                    "exception_type": type(exc).__name__,
                    "exception": str(exc),
                },
            )

        selected_key = max(
            channel_results,
            key=lambda key: float(channel_results[key].ratio_to_limit),
        )
        selected = channel_results[selected_key]
        selected_limit = self.anchor.channels[selected_key]
        passes = all(bool(result.passes) for result in channel_results.values())
        diagnostics: dict[str, Any] = {
            "evaluated": True,
            "selected_channel": selected_key,
            "channels": {
                key: _channel_diagnostics(
                    self.anchor.channels[key],
                    channel_results[key],
                    channel_proxies[key],
                )
                for key in channel_results
            },
            "sm_validation": {
                "zero_fcnc_coupling_branching_fractions": {
                    key: float(result.branching_fraction)
                    for key, result in self.sm_results.items()
                },
                "sm_hadronic_width_weight": float(
                    self.anchor.sm_hadronic_width_weight
                ),
                "sm_total_width_weight": float(self.anchor.sm_total_width_weight),
                "sm_hadronic_to_total_width_ratio": float(
                    self.anchor.sm_hadronic_to_total_width_ratio
                ),
                "sm_prediction_policy": (
                    "Pure-NP direct-width bound; catalog SM FCNC-Z rates are "
                    "negligible relative to 2.9e-3-scale direct limits."
                ),
            },
            "branching_formula": selected.diagnostics["branching_formula"],
            "sm_hadronic_width_weights": selected.diagnostics[
                "sm_hadronic_width_weights"
            ],
            "sm_total_width_weights": selected.diagnostics["sm_total_width_weights"],
            "budget_source": _BUDGET_SOURCE,
            "needs_human_physics": _NEEDS_HUMAN_PHYSICS,
            "rs_matching_assumption": ZPOLE_DOWN_FCNC_PROXY_V1,
            "required_parameter_point_extras": [_REQUIRED_EXTRA],
            "kk_ew_mass_extra_used": kk_ew_mass is not None,
        }

        return ConstraintResult(
            process_id=self.process_id,
            severity=self.severity,
            passes=passes,
            predicted=float(selected.branching_fraction),
            sm_prediction=0.0,
            experimental=float(selected_limit.value),
            ratio=float(selected.ratio_to_limit),
            budget=float(selected.br_limit),
            notes=(
                "Pure-NP direct hadronic-width bound on down-sector "
                "FCNC Z decays. The off-diagonal Z q_i q_j coupling is a "
                "documented down-overlap proxy and is flagged "
                "NEEDS-HUMAN-PHYSICS."
            ),
            diagnostics=diagnostics,
        )
