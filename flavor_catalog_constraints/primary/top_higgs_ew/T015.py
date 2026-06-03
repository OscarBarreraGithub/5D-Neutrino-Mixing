"""T015 - lepton-flavor-violating ``Z -> e mu`` branching fraction.

Physics
-------
The Standard Model contribution is zero for catalog purposes.  This constraint
therefore applies a pure-NP HARD upper bound to the charge-summed branching
fraction ``BR(Z -> e+- mu-+)``.  The rate uses the shared Z-pole
effective-coupling convention introduced for T010,

    L_Z = g_Z Z_mu lbar_i gamma^mu (delta g_L P_L + delta g_R P_R) l_j,

and the LFV extension computes

    BR = 2 (|delta g_L|^2 + |delta g_R|^2)
         / (SM total Z-width weight + LFV width weight),

where the factor of two implements the charge-summed ``e mu`` final state.

RS matching status
------------------
PARTIAL/NEEDS-HUMAN-PHYSICS.  The tree-level light-Z prediction is read from
``rs_ew_couplings.z_delta_g_L/R_e[0,1]`` and is rigorous within the Phase-4
minimal-RS neutral-current builder.  In the current diagonal charged-lepton
fit this matrix is diagonal, so the v1 tree-level LFV prediction is exactly
zero and non-vetoing.  Loop-induced Z-LFV from the lepton dipole spurion is
deferred to Phase 7 and is not faked here.

Catalog sidecar
---------------
``flavor_catalog/processes/top_higgs_ew/T015.yaml`` is the source of truth for
the CMS 2025 95% CL upper limit and provenance.  Its numeric limit is parsed
from the YAML entry and routed through the scaffold ``load_anchor`` path; no
experimental number is hardcoded here.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
import re
from typing import Any, Mapping, Sequence

from flavor_catalog_constraints import anchors as anchor_scaffold
from flavor_catalog_constraints.anchors import Anchor, AnchorError, load_anchor, load_full_yaml
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.zpole_lfv import (
    z_lfv_branching_fraction_from_couplings,
    zpole_lfv_effective_coupling_limit,
    zpole_lfv_sm_total_width_weight,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "top_higgs_ew"
_REQUIRED_EXTRA = "rs_ew_couplings"
_CMS_LIMIT_VALUE_ID = "CMS2025:T015:zemu_limit"
_PDG_LIMIT_VALUE_ID = "PDG2025:T015:zemu_limit"
_DATASET_VALUE_ID = "CMS2025:T015:dataset"
_NUMBER_RE = re.compile(
    r"^\s*<?\s*(?P<number>[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)"
    r"(?:[eE][+-]?[0-9]+)?)\s*$"
)
_UNEVALUATED_REASON = (
    "no off-diagonal lepton Z e-mu prediction available "
    "(rs_ew_couplings not provided on ParameterPoint)"
)
_UNEVALUATED_NOTES = f"NOT EVALUATED - {_UNEVALUATED_REASON}"
_NEEDS_HUMAN_PHYSICS = (
    "PARTIAL/NEEDS-HUMAN-PHYSICS: tree-level light-Z LFV uses rigorous "
    "rs_ew_couplings.z_delta_g_L/R_e[0,1]; with the v1 diagonal "
    "charged-lepton fit this off-diagonal entry is zero. Loop-induced Z-LFV "
    "from the lfv_dipole_spurion is deferred to Phase 7."
)
_TREE_LEVEL_STATUS = "rigorous_tree_light_z_from_rs_ew_couplings"
_LOOP_LFV_STATUS = "deferred_to_phase_7_lfv_dipole_spurion"


@dataclass(frozen=True)
class LFVLimitEntry:
    """Typed view over one value entry in the T015 YAML ``values`` list."""

    value_id: str
    observable: str | None
    year: int | None
    limit: float
    expected_limit: float | None
    confidence_level: str | None
    units: str | None
    source_url: str | None
    snapshot_path: str | None
    value_summary: str | None


@dataclass(frozen=True)
class T015Anchor:
    """Typed T015 anchor bundle: CMS limit plus PDG context."""

    experimental: Anchor
    cms_limit: LFVLimitEntry
    pdg_limit: LFVLimitEntry
    dataset: Mapping[str, Any]
    effective_coupling_limit: float
    sm_total_width_weight: float

    @property
    def value(self) -> float:
        return self.experimental.value

    @property
    def budget(self) -> float:
        return self.experimental.value

    @property
    def source_url(self) -> str | None:
        return self.experimental.source_url


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
                f"{process_id}: T015 field {field_name!r}={value!r} is not "
                "a numeric upper-limit string"
            )
        number = float(match.group("number"))
    else:
        raise AnchorError(
            f"{process_id}: T015 field {field_name!r}={value!r} is not numeric"
        )
    if not math.isfinite(number) or number <= 0.0:
        raise AnchorError(f"{process_id}: T015 field {field_name!r} must be positive")
    return number


def _pdg_values(process_id: str) -> Sequence[Mapping[str, Any]]:
    data = load_full_yaml(process_id, family=_FAMILY)
    pdg = data.get("pdg_or_equivalent")
    if not isinstance(pdg, Mapping):
        raise AnchorError(f"{process_id}: expected pdg_or_equivalent mapping")
    values = pdg.get("values")
    if not isinstance(values, list) or not values:
        raise AnchorError(f"{process_id}: expected non-empty pdg_or_equivalent.values")
    for index, entry in enumerate(values):
        if not isinstance(entry, Mapping):
            raise AnchorError(
                f"{process_id}: pdg_or_equivalent.values[{index}] is not a mapping"
            )
    return values


def _find_value_entry(
    process_id: str,
    value_id: str,
) -> tuple[int, Mapping[str, Any]]:
    for index, entry in enumerate(_pdg_values(process_id)):
        if entry.get("value_id") == value_id:
            return index, entry
    present = [str(entry.get("value_id")) for entry in _pdg_values(process_id)]
    raise AnchorError(
        f"{process_id}: value_id {value_id!r} not found in T015 values "
        f"(present: {present})"
    )


def _load_scaffold_value_anchor(
    value_id: str,
    *,
    process_id: str,
) -> Anchor:
    index, entry = _find_value_entry(process_id, value_id)
    block_key = f"pdg_or_equivalent.values[{index}]"
    numeric_entry = dict(entry)
    numeric_entry["value"] = _parse_limit_value(
        entry.get("normalized_value", entry.get("value")),
        process_id=process_id,
        field_name=f"{value_id}.normalized_value",
    )
    if "snapshot_path" not in numeric_entry:
        pdg = load_full_yaml(process_id, family=_FAMILY)["pdg_or_equivalent"]
        if isinstance(pdg, Mapping) and "snapshot_path" in pdg:
            numeric_entry["snapshot_path"] = pdg["snapshot_path"]
    virtual_block = {block_key: numeric_entry}
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

    if scaffold_anchor.block_key != block_key:
        raise AnchorError(
            f"{process_id}: load_anchor selected {scaffold_anchor.block_key!r}, "
            f"expected {block_key!r} for value_id {value_id!r}"
        )
    return scaffold_anchor


def _entry_view(
    process_id: str,
    value_id: str,
    *,
    value_summary: str | None,
) -> LFVLimitEntry:
    _, entry = _find_value_entry(process_id, value_id)
    return LFVLimitEntry(
        value_id=value_id,
        observable=_optional_str(entry.get("observable")),
        year=_optional_int(entry.get("year")),
        limit=_parse_limit_value(
            entry.get("normalized_value", entry.get("value")),
            process_id=process_id,
            field_name=f"{value_id}.normalized_value",
        ),
        expected_limit=(
            None
            if entry.get("expected_value") is None
            else _parse_limit_value(
                entry.get("expected_value"),
                process_id=process_id,
                field_name=f"{value_id}.expected_value",
            )
        ),
        confidence_level=_optional_str(entry.get("confidence_level")),
        units=_optional_str(entry.get("units")),
        source_url=_optional_str(entry.get("source_url")),
        snapshot_path=_optional_str(entry.get("snapshot_path")),
        value_summary=value_summary,
    )


def _dataset_entry(process_id: str) -> Mapping[str, Any]:
    _, entry = _find_value_entry(process_id, _DATASET_VALUE_ID)
    return dict(entry)


def _load_t015_anchor(process_id: str) -> T015Anchor:
    data = load_full_yaml(process_id, family=_FAMILY)
    pdg = data.get("pdg_or_equivalent")
    if not isinstance(pdg, Mapping):
        raise AnchorError(f"{process_id}: expected pdg_or_equivalent mapping")
    value_summary = _optional_str(pdg.get("value_summary"))
    experimental = _load_scaffold_value_anchor(
        _CMS_LIMIT_VALUE_ID,
        process_id=process_id,
    )
    if experimental.units != "branching fraction":
        raise AnchorError(
            f"{process_id}: CMS limit units must be 'branching fraction', "
            f"got {experimental.units!r}"
        )
    cms = _entry_view(process_id, _CMS_LIMIT_VALUE_ID, value_summary=value_summary)
    pdg_limit = _entry_view(process_id, _PDG_LIMIT_VALUE_ID, value_summary=value_summary)
    if not math.isclose(experimental.value, cms.limit, rel_tol=0.0, abs_tol=0.0):
        raise AnchorError(f"{process_id}: scaffold anchor and CMS entry disagree")
    weights = zpole_lfv_sm_total_width_weight()
    sm_total = float(sum(weights.values()))
    return T015Anchor(
        experimental=experimental,
        cms_limit=cms,
        pdg_limit=pdg_limit,
        dataset=_dataset_entry(process_id),
        effective_coupling_limit=zpole_lfv_effective_coupling_limit(
            experimental.value
        ),
        sm_total_width_weight=sm_total,
    )


def _coupling_entry(source: Any, matrix_name: str, row: int, column: int) -> complex:
    try:
        value = complex(getattr(source, matrix_name)[row, column])
    except (AttributeError, TypeError, KeyError, IndexError) as exc:
        raise ValueError(f"{matrix_name}[{row},{column}] is not available") from exc
    if not math.isfinite(value.real) or not math.isfinite(value.imag):
        raise ValueError(f"{matrix_name}[{row},{column}] must be finite")
    return value


@register
class Constraint:
    """Catalogued ``Z -> e mu`` LFV pure-NP branching-fraction constraint."""

    process_id = "T015"
    severity = Severity.HARD
    observable = "BR(Z -> e mu)"

    def __init__(self) -> None:
        self.anchor = _load_t015_anchor(self.process_id)
        self.sm_result = z_lfv_branching_fraction_from_couplings(
            delta_g_left=0.0j,
            delta_g_right=0.0j,
            initial_flavor="e",
            final_flavor="mu",
            br_limit=self.anchor.budget,
        )

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
                    "non-vetoing only; no BR(Z -> e mu) NP prediction was evaluated"
                ),
                "needs_human_physics": _NEEDS_HUMAN_PHYSICS,
                "required_parameter_point_extras": [_REQUIRED_EXTRA],
                "budget_source": self.anchor.source_url,
                "experimental_block": self.anchor.experimental.block_key,
                "sm_branching_fraction": 0.0,
                **dict(diagnostics),
            },
        )

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        rs_ew_couplings = point.get_extra(_REQUIRED_EXTRA)
        if rs_ew_couplings is None:
            return self._unevaluated_result(
                {
                    "missing_extra": _REQUIRED_EXTRA,
                    "legacy_lepton_mass_basis_couplings_present": (
                        point.get_extra("lepton_mass_basis_couplings") is not None
                    ),
                },
            )

        try:
            delta_g_left = _coupling_entry(rs_ew_couplings, "z_delta_g_L_e", 0, 1)
            delta_g_right = _coupling_entry(rs_ew_couplings, "z_delta_g_R_e", 0, 1)
            result = z_lfv_branching_fraction_from_couplings(
                delta_g_left=delta_g_left,
                delta_g_right=delta_g_right,
                initial_flavor="e",
                final_flavor="mu",
                br_limit=self.anchor.budget,
            )
        except (AttributeError, KeyError, TypeError, ValueError) as exc:
            return self._unevaluated_result(
                {
                    "invalid_extra": _REQUIRED_EXTRA,
                    "exception_type": type(exc).__name__,
                    "exception": str(exc),
                },
            )

        diagnostics = dict(result.diagnostics)
        diagnostics.update(
            {
                "evaluated": True,
                "sm_prediction_policy": (
                    "Charged-LFV Z decay has zero SM rate for catalog purposes; "
                    "the HARD budget is applied to the pure-NP branching fraction."
                ),
                "needs_human_physics": _NEEDS_HUMAN_PHYSICS,
                "tree_level_matching_status": _TREE_LEVEL_STATUS,
                "loop_lfv_status": _LOOP_LFV_STATUS,
                "rs_matching_assumption": getattr(
                    rs_ew_couplings, "matching_assumption", None
                ),
                "rs_ew_model_label": getattr(rs_ew_couplings, "model_label", None),
                "rs_ew_kk_mass_gev": float(
                    getattr(rs_ew_couplings, "kk_ew_mass_gev")
                ),
                "z_delta_g_indices": {
                    "left": "z_delta_g_L_e[0,1]",
                    "right": "z_delta_g_R_e[0,1]",
                },
                "initial_flavor": result.initial_flavor,
                "final_flavor": result.final_flavor,
                "delta_g_left_emu": complex(result.delta_g_left),
                "delta_g_right_emu": complex(result.delta_g_right),
                "effective_coupling_norm": float(result.coupling_norm),
                "effective_coupling_limit": float(
                    self.anchor.effective_coupling_limit
                ),
                "lfv_width_weight": float(result.lfv_width_weight),
                "sm_total_width_weight": float(result.sm_total_width_weight),
                "total_width_weight": float(result.total_width_weight),
                "charge_state_factor": float(result.charge_state_factor),
                "experimental_block": self.anchor.experimental.block_key,
                "cms_limit_value_id": self.anchor.cms_limit.value_id,
                "cms_expected_limit": self.anchor.cms_limit.expected_limit,
                "pdg_listing_limit": float(self.anchor.pdg_limit.limit),
                "pdg_listing_value_id": self.anchor.pdg_limit.value_id,
                "dataset": dict(self.anchor.dataset),
                "budget_source": self.anchor.source_url,
                "required_parameter_point_extras": [_REQUIRED_EXTRA],
            }
        )

        return ConstraintResult(
            process_id=self.process_id,
            severity=self.severity,
            passes=bool(result.passes),
            predicted=float(result.branching_fraction),
            sm_prediction=0.0,
            experimental=float(self.anchor.value),
            ratio=float(result.ratio_to_limit),
            budget=float(result.br_limit),
            notes=(
                "Pure-NP BR(Z -> e mu) bound using the shared zpole "
                "effective-coupling width weights. The off-diagonal Z e-mu "
                "coupling is the Phase-4 tree-level light-Z "
                "rs_ew_couplings.z_delta_g_e[0,1] entry; loop-induced Z-LFV "
                "from the dipole spurion is deferred to Phase 7."
            ),
            diagnostics=diagnostics,
        )
