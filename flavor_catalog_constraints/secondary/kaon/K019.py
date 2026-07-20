"""K019 - LFV long-lived neutral-kaon decay ``K_L -> e+- mu-+``.

Physics
-------
``K_L -> e mu`` is charged-lepton-flavor violating, so the Standard Model rate
is zero for catalog purposes.  K019 therefore applies a pure-NP upper bound,

    BR_NP(K_L -> e+- mu-+) <= BR_exp^upper.

This module reuses the shared ``s -> d l l`` rare-kaon dilepton machinery
through the
``flavor_catalog_constraints.physics_adapters.rare_kaon_lfv_dilepton``
boundary.  The K019 adapter keeps the K006 ``kappa_mu`` normalization and
feeds the Phase-4a rigorous ``rs_semileptonic_wilsons.lfv_llqq`` e-mu block
into the unequal-lepton two-body phase-space normalization.

Severity
--------
HARD.  The K019 YAML anchor is a direct 90% CL upper limit on a charged-LFV
branching fraction from BNL E871 / the PDG K_L listing.  A point with an
evaluated pure-NP prediction above that limit is excluded.

Catalog sidecar
---------------
``flavor_catalog/processes/secondary/kaon/K019.yaml`` is the source of truth
for the BNL E871/PDG branching-fraction limit and provenance.  Numeric values
below are loaded from that sidecar through the scaffold anchor loader, not
hardcoded here.

Phase-4c status
---------------
The tree-level LFV lepton coupling is now read from the Phase-4a lepton-aware
semileptonic Wilson bundle.  For the current diagonal charged-lepton fit it is
rigorously zero, so the tree-level LFV rate is zero and non-vetoing.  Nonzero
tree-level rates require non-diagonal lepton structure; loop-induced LFV is
deferred.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Mapping

from flavor_catalog_constraints import anchors as anchor_scaffold
from flavor_catalog_constraints.anchors import (
    Anchor,
    AnchorError,
    load_anchor,
    load_full_yaml,
)
from flavor_catalog_constraints.base import (
    ConstraintLevel,
    ConstraintResult,
    ParameterPoint,
    Severity,
)
from flavor_catalog_constraints.physics_adapters.rare_kaon_lfv_dilepton import (
    RARE_KAON_LFV_DILEPTON_PARAMETRIZATION_CITATION,
    RARE_KAON_LFV_TREE_LEVEL_NOTE_V1,
    klong_emu_from_rs_semileptonic_wilsons,
    rare_kaon_lfv_default_sm_inputs,
    rare_kaon_lfv_sm_branching_fraction,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "kaon"
_TIER = ConstraintLevel.SECONDARY
_REQUIRED_EXTRA = "rs_semileptonic_wilsons"
_OPTIONAL_EW_MASS_EXTRA = "kk_ew_mass_gev"
_PDG_BLOCK = "pdg_or_equivalent"
_EXPERIMENT_INPUT_KEY = "experiment_input"
_EXPECTED_UNITS = "dimensionless branching fraction"
_EXPECTED_LIMIT_TYPE = "upper_limit"
_EXPECTED_CL = "90% CL"
_BUDGET_SOURCE = (
    "flavor_catalog/processes/secondary/kaon/K019.yaml pdg_or_equivalent "
    "(PDG 2025 K_L listing / BNL E871 final 90% CL limit)"
)
_UNEVALUATED_REASON = "missing rs_semileptonic_wilsons LFV llqq block"
_UNEVALUATED_NOTES = f"NOT EVALUATED - {_UNEVALUATED_REASON}"


@dataclass(frozen=True)
class BranchingLimitAnchor:
    """Typed upper-limit branching-fraction anchor."""

    block_key: str
    source: str | None
    year: int | None
    value: float
    limit_type: str
    confidence_level: str | None
    units: str | None
    source_url: str | None
    snapshot_path: str | None
    display: str | None = None
    experiment: str | None = None


@dataclass(frozen=True)
class K019Anchor:
    """Typed K019 anchor: PDG current limit plus BNL E871 provenance."""

    current_limit: BranchingLimitAnchor
    bnl_e871_limit: BranchingLimitAnchor

    @property
    def value(self) -> float:
        """Upper limit used as the experimental value."""

        return self.current_limit.value

    @property
    def budget(self) -> float:
        """HARD pure-NP branching-fraction budget."""

        return self.current_limit.value


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
        out = float(value)
    except (TypeError, ValueError) as exc:
        raise AnchorError(
            f"{process_id}: K019 anchor field {field_name!r}={value!r} "
            "is not numeric"
        ) from exc
    if not math.isfinite(out):
        raise AnchorError(f"{process_id}: K019 anchor field {field_name!r} is not finite")
    return out


def _positive_float(value: Any, *, process_id: str, field_name: str) -> float:
    out = _required_float(value, process_id=process_id, field_name=field_name)
    if out <= 0.0:
        raise AnchorError(f"{process_id}: K019 anchor field {field_name!r} must be positive")
    return out


def _validate_limit_block(
    sub: Mapping[str, Any],
    *,
    process_id: str,
    block_key: str,
) -> None:
    units = _optional_str(sub.get("units"))
    if units != _EXPECTED_UNITS:
        raise AnchorError(
            f"{process_id}: {block_key} must use units {_EXPECTED_UNITS!r}, "
            f"got {units!r}"
        )
    limit_type = _optional_str(sub.get("limit_type"))
    if limit_type != _EXPECTED_LIMIT_TYPE:
        raise AnchorError(
            f"{process_id}: {block_key}.limit_type must be "
            f"{_EXPECTED_LIMIT_TYPE!r}, got {limit_type!r}"
        )
    confidence_level = _optional_str(sub.get("cl"))
    if confidence_level != _EXPECTED_CL:
        raise AnchorError(
            f"{process_id}: {block_key}.cl must be {_EXPECTED_CL!r}, "
            f"got {confidence_level!r}"
        )


def _limit_from_anchor(
    anchor: Anchor,
    sub: Mapping[str, Any],
    *,
    process_id: str,
    block_key: str,
) -> BranchingLimitAnchor:
    _validate_limit_block(sub, process_id=process_id, block_key=block_key)
    return BranchingLimitAnchor(
        block_key=block_key,
        source=_optional_str(anchor.source or sub.get("source")),
        year=anchor.year if anchor.year is not None else _optional_int(sub.get("year")),
        value=_positive_float(
            anchor.value,
            process_id=process_id,
            field_name=f"{block_key}.value",
        ),
        limit_type=str(sub.get("limit_type")),
        confidence_level=_optional_str(sub.get("cl")),
        units=_optional_str(anchor.units or sub.get("units")),
        source_url=_optional_str(anchor.source_url or sub.get("source_url")),
        snapshot_path=_optional_str(anchor.snapshot_path or sub.get("snapshot_path")),
        display=_optional_str(sub.get("display")),
        experiment=_optional_str(sub.get("experiment")),
    )


def _limit_from_subblock(
    sub: Mapping[str, Any],
    *,
    process_id: str,
    block_key: str,
) -> BranchingLimitAnchor:
    _validate_limit_block(sub, process_id=process_id, block_key=block_key)
    return BranchingLimitAnchor(
        block_key=block_key,
        source=_optional_str(sub.get("source")),
        year=_optional_int(sub.get("year")),
        value=_positive_float(
            sub.get("value"),
            process_id=process_id,
            field_name=f"{block_key}.value",
        ),
        limit_type=str(sub.get("limit_type")),
        confidence_level=_optional_str(sub.get("cl")),
        units=_optional_str(sub.get("units")),
        source_url=_optional_str(sub.get("source_url")),
        snapshot_path=_optional_str(sub.get("snapshot_path")),
        display=_optional_str(sub.get("display")),
        experiment=_optional_str(sub.get("experiment")),
    )


def _load_flat_pdg_anchor(process_id: str) -> tuple[Anchor, Mapping[str, Any]]:
    data = load_full_yaml(process_id, family=_FAMILY, tier=_TIER)
    pdg = data.get(_PDG_BLOCK)
    if not isinstance(pdg, Mapping):
        raise AnchorError(f"{process_id}: expected {_PDG_BLOCK} mapping")
    virtual_block = {_PDG_BLOCK: dict(pdg)}
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
            candidates=(_PDG_BLOCK,),
        )
    finally:
        anchor_scaffold.load_pdg_block = original_load_pdg_block

    if scaffold_anchor.block_key != _PDG_BLOCK:
        raise AnchorError(
            f"{process_id}: load_anchor selected {scaffold_anchor.block_key!r}, "
            f"expected {_PDG_BLOCK!r} for K019 current limit"
        )
    return scaffold_anchor, pdg


def _load_k019_anchor(process_id: str) -> K019Anchor:
    scaffold_anchor, pdg = _load_flat_pdg_anchor(process_id)
    current = _limit_from_anchor(
        scaffold_anchor,
        pdg,
        process_id=process_id,
        block_key=_PDG_BLOCK,
    )
    experiment_input = pdg.get(_EXPERIMENT_INPUT_KEY)
    if not isinstance(experiment_input, Mapping):
        raise AnchorError(
            f"{process_id}: expected {_PDG_BLOCK}.{_EXPERIMENT_INPUT_KEY} mapping"
        )
    bnl = _limit_from_subblock(
        experiment_input,
        process_id=process_id,
        block_key=f"{_PDG_BLOCK}.{_EXPERIMENT_INPUT_KEY}",
    )
    if not math.isclose(current.value, bnl.value, rel_tol=0.0, abs_tol=0.0):
        raise AnchorError(f"{process_id}: PDG and BNL E871 K019 limits disagree")
    if current.value <= 0.0:
        raise AnchorError(f"{process_id}: K019 HARD budget must be positive")
    return K019Anchor(current_limit=current, bnl_e871_limit=bnl)


def _complex_wilsons(result: Any) -> dict[str, complex]:
    if result.wilsons is None:
        return {}
    return {key: complex(value) for key, value in result.wilsons.wilsons.items()}


@register
class Constraint:
    """Catalogued pure-NP ``K_L -> e+- mu-+`` LFV constraint."""

    process_id = "K019"
    severity = Severity.HARD
    observable = "BR(K_L -> e+- mu-+) LFV"

    def __init__(self) -> None:
        self.anchor = _load_k019_anchor(self.process_id)
        self.sm_inputs = rare_kaon_lfv_default_sm_inputs()
        self.sm_result = rare_kaon_lfv_sm_branching_fraction(self.sm_inputs)

    def _base_diagnostics(self) -> dict[str, Any]:
        return {
            "experimental_block": self.anchor.current_limit.block_key,
            "experimental_confidence_level": self.anchor.current_limit.confidence_level,
            "bnl_e871_block": self.anchor.bnl_e871_limit.block_key,
            "bnl_e871_source": self.anchor.bnl_e871_limit.source,
            "bnl_e871_source_url": self.anchor.bnl_e871_limit.source_url,
            "budget_source": _BUDGET_SOURCE,
            "budget_is_upper_limit": True,
            "sm_branching_fraction": 0.0,
            "sm_lfv_policy": (
                "K_L -> e mu is charged-LFV and has zero SM rate for catalog "
                "purposes; the HARD budget is applied to the pure-NP "
                "tree-level prediction when the LFV llqq block is present."
            ),
            "charge_conjugate_modes_included": True,
            "parametrization_citation": (
                RARE_KAON_LFV_DILEPTON_PARAMETRIZATION_CITATION
            ),
            "lfv_tree_level_note": RARE_KAON_LFV_TREE_LEVEL_NOTE_V1,
            "loop_lfv_status": "loop_induced_lfv_deferred",
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
                    "non-vetoing only; no K_L -> e mu NP prediction was evaluated"
                ),
                **self._base_diagnostics(),
                **dict(diagnostics),
            },
        )

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        wilson_input = point.get_extra(_REQUIRED_EXTRA)
        if wilson_input is None:
            return self._unevaluated_result({"missing_extra": _REQUIRED_EXTRA})

        kk_ew_mass = point.get_extra(_OPTIONAL_EW_MASS_EXTRA)
        try:
            result = klong_emu_from_rs_semileptonic_wilsons(
                wilson_input,
                matching_scale_gev=None if kk_ew_mass is None else float(kk_ew_mass),
                inputs=self.sm_inputs,
            )
        except (AttributeError, KeyError, TypeError, ValueError) as exc:
            return self._unevaluated_result(
                {
                    "invalid_extra": _REQUIRED_EXTRA,
                    "exception_type": type(exc).__name__,
                    "exception": str(exc),
                }
            )

        predicted = float(result.branching_fraction)
        budget = float(self.anchor.budget)
        ratio = predicted / budget if budget > 0.0 else float("inf")
        diagnostics = dict(result.diagnostics)
        diagnostics.update(
            {
                "evaluated": True,
                **self._base_diagnostics(),
                "np_shift_branching_fraction": float(
                    result.np_shift_branching_fraction
                ),
                "kk_ew_mass_extra_used": kk_ew_mass is not None,
                "wilson_coefficients": _complex_wilsons(result),
            }
        )

        return ConstraintResult(
            process_id=self.process_id,
            severity=self.severity,
            passes=bool(ratio <= 1.0),
            predicted=predicted,
            sm_prediction=0.0,
            experimental=float(self.anchor.value),
            ratio=float(ratio),
            budget=budget,
            notes=(
                "Pure-NP BR(K_L -> e+- mu-+) bound using Phase-4a LFV llqq "
                "Wilsons additively mapped into the rare-kaon Y inputs. "
                "Tree-level LFV is rigorous and zero for the diagonal "
                "charged-lepton fit; loop-induced LFV is deferred."
            ),
            diagnostics=diagnostics,
        )
