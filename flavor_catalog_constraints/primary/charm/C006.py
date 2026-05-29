"""C006 - LFV neutral-charm decay ``D0 -> e+- mu-+``.

Physics
-------
The Standard Model rate is zero for catalog purposes, so C006 is a pure-NP
upper-bound constraint:

    BR_NP(D0 -> e+- mu-+) <= BR_exp^upper.

This module reuses the shared C004/C005 ``c -> u l l`` machinery through the
``flavor_catalog_constraints.physics_adapters.rare_charm_lfv_dilepton``
boundary.  The LFV extension keeps the same rare-charm Hamiltonian and uses
the unequal-lepton pseudoscalar two-body formula with ``C9-C9'`` and
``C10-C10'``.

Severity
--------
HARD.  The C006 YAML anchor is a direct observed 90% CL upper limit on a
charged-LFV branching fraction.  A point with an evaluated pure-NP prediction
above that limit is excluded.

Catalog sidecar
---------------
``flavor_catalog/processes/charm/C006.yaml`` is the source of truth for the
PDG/LHCb branching-fraction limit and provenance.  Numeric values below are
loaded from that sidecar, not hardcoded here.

NEEDS-HUMAN-PHYSICS
-------------------
A rigorous RS prediction needs the off-diagonal charged-lepton neutral-current
``e mu`` coupling after EW KK/Z/Z' mixing and charged-lepton mass-basis
rotation.  That coupling is not available as a standard ParameterPoint input.
This v1 constraint accepts an explicit ``lepton_mass_basis_couplings`` proxy
and marks the result NEEDS-HUMAN-PHYSICS.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Mapping

from flavor_catalog_constraints.anchors import (
    Anchor,
    AnchorError,
    load_anchor,
    load_full_yaml,
    load_pdg_block,
)
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.rare_charm_lfv_dilepton import (
    RARE_CHARM_LFV_DILEPTON_PROXY_V1,
    d0_emu_from_couplings,
    rare_charm_lfv_default_sm_inputs,
    rare_charm_lfv_sm_branching_fraction,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "charm"
_REQUIRED_QUARK_EXTRA = "quark_mass_basis_couplings"
_REQUIRED_LEPTON_EXTRA = "lepton_mass_basis_couplings"
_OPTIONAL_EW_MASS_EXTRA = "kk_ew_mass_gev"
_EXPECTED_UNITS = "branching fraction"

_CURRENT_LIMIT_CANDIDATES = ("canonical_current_limit",)
_LHCB_LIMIT_CANDIDATES = ("lhcb_primary_result",)
_PREVIOUS_LIMITS_KEY = "previous_limits"
_BELLE_PREVIOUS_KEY = "belle_2010"
_BABAR_PREVIOUS_KEY = "babar_2012"
_PAPER_ERA_SECTION = "paper_era_reference"
_RS_BASELINE_KEY = "rs_baseline"
_BUDGET_SOURCE = (
    "flavor_catalog/processes/charm/C006.yaml canonical_current_limit "
    "(PDG Live/API S032.40, 90% CL; LHCb 2016)"
)
_PARAMETRIZATION_CITATION = (
    "Rare-charm c->u l l C9/C10 Hamiltonian reused from C004/C005; "
    "unequal-lepton pseudoscalar D0 -> e mu two-body normalization"
)
_UNEVALUATED_REASON = (
    "missing c-u quark coupling or explicit e-mu lepton LFV proxy"
)
_UNEVALUATED_NOTES = f"NOT EVALUATED - {_UNEVALUATED_REASON}"
_NEEDS_HUMAN_PHYSICS = (
    "NEEDS-HUMAN-PHYSICS: the off-diagonal e-mu lepton neutral-current "
    "coupling is not a standard ParameterPoint input; C006 v1 requires an "
    "explicit lepton_mass_basis_couplings proxy and reuses the documented "
    "rare-charm LFV Z-like matching."
)


@dataclass(frozen=True)
class BranchingLimitAnchor:
    """Typed upper-limit branching-fraction anchor."""

    block_key: str
    source: str | None
    year: int | None
    value: float
    confidence_level: float
    limit_type: str
    units: str | None
    source_url: str | None
    snapshot_path: str | None
    integrated_luminosity_fb_inv: float | None = None
    collision_energy_tev: tuple[float, ...] = ()


@dataclass(frozen=True)
class RSBaselineContext:
    """Typed paper-era RS context carried for proxy provenance."""

    section_key: str
    block_key: str
    source: str | None
    year: int | None
    generic_rs_kk_gluon_scale_tev: float
    composite_pseudo_goldstone_kk_gluon_scale_tev: float
    source_url: str | None
    snapshot_path: str | None


@dataclass(frozen=True)
class C006Anchor:
    """Typed C006 anchor: current limit, previous limits, and RS context."""

    current_limit: BranchingLimitAnchor
    lhcb_primary_limit: BranchingLimitAnchor
    belle_previous_limit: BranchingLimitAnchor
    babar_previous_limit: BranchingLimitAnchor
    rs_baseline: RSBaselineContext

    @property
    def value(self) -> float:
        """Current PDG/LHCb upper limit used as the experimental value."""

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


def _optional_float(
    value: Any,
    *,
    process_id: str,
    field_name: str,
) -> float | None:
    if value is None:
        return None
    return _required_float(value, process_id=process_id, field_name=field_name)


def _required_float(value: Any, *, process_id: str, field_name: str) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError) as exc:
        raise AnchorError(
            f"{process_id}: C006 anchor field {field_name!r}={value!r} is not numeric"
        ) from exc
    if not math.isfinite(out):
        raise AnchorError(f"{process_id}: C006 anchor field {field_name!r} is not finite")
    return out


def _positive_float(value: Any, *, process_id: str, field_name: str) -> float:
    out = _required_float(value, process_id=process_id, field_name=field_name)
    if out <= 0.0:
        raise AnchorError(f"{process_id}: C006 anchor field {field_name!r} must be positive")
    return out


def _float_tuple(
    value: Any,
    *,
    process_id: str,
    field_name: str,
) -> tuple[float, ...]:
    if value is None:
        return ()
    if not isinstance(value, (list, tuple)):
        return (
            _required_float(value, process_id=process_id, field_name=field_name),
        )
    return tuple(
        _required_float(item, process_id=process_id, field_name=f"{field_name}[{idx}]")
        for idx, item in enumerate(value)
    )


def _pdg_subblock_for_anchor(anchor: Anchor, *, process_id: str) -> Mapping[str, Any]:
    data = load_full_yaml(process_id, family=_FAMILY)
    pdg_block = data.get("pdg_or_equivalent")
    if not isinstance(pdg_block, Mapping):
        raise AnchorError(
            f"{process_id}: 'pdg_or_equivalent' is not a mapping while loading "
            f"{anchor.block_key}"
        )
    sub = pdg_block.get(anchor.block_key)
    if not isinstance(sub, Mapping):
        present = sorted(str(key) for key in pdg_block)
        raise AnchorError(
            f"{process_id}: load_anchor selected {anchor.block_key!r}, but that "
            f"block is not available as a mapping (present keys: {present})"
        )
    return sub


def _limit_from_subblock(
    sub: Mapping[str, Any],
    *,
    process_id: str,
    block_key: str,
) -> BranchingLimitAnchor:
    units = _optional_str(sub.get("units"))
    if units != _EXPECTED_UNITS:
        raise AnchorError(
            f"{process_id}: {block_key} must use units {_EXPECTED_UNITS!r}, "
            f"got {units!r}"
        )
    value = _positive_float(
        sub.get("value"),
        process_id=process_id,
        field_name=f"{block_key}.value",
    )
    confidence_level = _positive_float(
        sub.get("confidence_level"),
        process_id=process_id,
        field_name=f"{block_key}.confidence_level",
    )
    limit_type = str(sub.get("limit_type", "")).lower()
    if limit_type != "upper":
        raise AnchorError(f"{process_id}: {block_key}.limit_type must be 'upper'")
    return BranchingLimitAnchor(
        block_key=block_key,
        source=_optional_str(sub.get("source")),
        year=_optional_int(sub.get("year")),
        value=float(value),
        confidence_level=float(confidence_level),
        limit_type=limit_type,
        units=units,
        source_url=_optional_str(sub.get("source_url")),
        snapshot_path=_optional_str(sub.get("snapshot_path")),
        integrated_luminosity_fb_inv=_optional_float(
            sub.get("integrated_luminosity_fb_inv"),
            process_id=process_id,
            field_name=f"{block_key}.integrated_luminosity_fb_inv",
        ),
        collision_energy_tev=_float_tuple(
            sub.get("collision_energy_TeV"),
            process_id=process_id,
            field_name=f"{block_key}.collision_energy_TeV",
        ),
    )


def _load_limit_anchor(
    candidates: tuple[str, ...],
    *,
    process_id: str,
) -> BranchingLimitAnchor:
    scaffold_anchor = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=candidates,
    )
    sub = _pdg_subblock_for_anchor(scaffold_anchor, process_id=process_id)
    return _limit_from_subblock(
        sub,
        process_id=process_id,
        block_key=scaffold_anchor.block_key,
    )


def _load_previous_limit(
    process_id: str,
    previous_key: str,
) -> BranchingLimitAnchor:
    pdg = load_pdg_block(process_id, family=_FAMILY)
    previous = pdg.get(_PREVIOUS_LIMITS_KEY)
    if not isinstance(previous, Mapping):
        raise AnchorError(
            f"{process_id}: pdg_or_equivalent.{_PREVIOUS_LIMITS_KEY} missing "
            "or not a mapping"
        )
    sub = previous.get(previous_key)
    if not isinstance(sub, Mapping):
        present = sorted(str(key) for key in previous)
        raise AnchorError(
            f"{process_id}: previous_limits.{previous_key} missing or not a "
            f"mapping (present keys: {present})"
        )
    return _limit_from_subblock(
        sub,
        process_id=process_id,
        block_key=f"{_PREVIOUS_LIMITS_KEY}.{previous_key}",
    )


def _top_level_subblock(
    process_id: str,
    *,
    section_key: str,
    block_key: str,
) -> Mapping[str, Any]:
    data = load_full_yaml(process_id, family=_FAMILY)
    section = data.get(section_key)
    if not isinstance(section, Mapping):
        raise AnchorError(
            f"{process_id}: top-level {section_key!r} is not a mapping while "
            f"loading {block_key!r}"
        )
    sub = section.get(block_key)
    if not isinstance(sub, Mapping):
        present = sorted(str(key) for key in section)
        raise AnchorError(
            f"{process_id}: {section_key}.{block_key} is missing or not a mapping "
            f"(present keys: {present})"
        )
    return sub


def _load_rs_baseline(process_id: str) -> RSBaselineContext:
    sub = _top_level_subblock(
        process_id,
        section_key=_PAPER_ERA_SECTION,
        block_key=_RS_BASELINE_KEY,
    )
    return RSBaselineContext(
        section_key=_PAPER_ERA_SECTION,
        block_key=_RS_BASELINE_KEY,
        source=_optional_str(sub.get("source")),
        year=_optional_int(sub.get("year")),
        generic_rs_kk_gluon_scale_tev=_positive_float(
            sub.get("generic_rs_kk_gluon_scale_TeV"),
            process_id=process_id,
            field_name="paper_era_reference.rs_baseline.generic_rs_kk_gluon_scale_TeV",
        ),
        composite_pseudo_goldstone_kk_gluon_scale_tev=_positive_float(
            sub.get("composite_pseudo_goldstone_kk_gluon_scale_TeV"),
            process_id=process_id,
            field_name=(
                "paper_era_reference.rs_baseline."
                "composite_pseudo_goldstone_kk_gluon_scale_TeV"
            ),
        ),
        source_url=_optional_str(sub.get("source_url")),
        snapshot_path=_optional_str(sub.get("snapshot_path")),
    )


def _load_c006_anchor(process_id: str) -> C006Anchor:
    anchor = C006Anchor(
        current_limit=_load_limit_anchor(
            _CURRENT_LIMIT_CANDIDATES,
            process_id=process_id,
        ),
        lhcb_primary_limit=_load_limit_anchor(
            _LHCB_LIMIT_CANDIDATES,
            process_id=process_id,
        ),
        belle_previous_limit=_load_previous_limit(process_id, _BELLE_PREVIOUS_KEY),
        babar_previous_limit=_load_previous_limit(process_id, _BABAR_PREVIOUS_KEY),
        rs_baseline=_load_rs_baseline(process_id),
    )
    if anchor.budget <= 0.0:
        raise AnchorError(f"{process_id}: C006 HARD budget must be positive")
    return anchor


def _complex_wilsons(result: Any) -> dict[str, complex]:
    if result.wilsons is None:
        return {}
    return {key: complex(value) for key, value in result.wilsons.wilsons.items()}


@register
class Constraint:
    """Catalogued pure-NP ``D0 -> e+- mu-+`` LFV constraint."""

    process_id = "C006"
    severity = Severity.HARD
    observable = "BR(D0 -> e+- mu-+) LFV"

    def __init__(self) -> None:
        self.anchor = _load_c006_anchor(self.process_id)
        self.sm_inputs = rare_charm_lfv_default_sm_inputs()
        self.sm_result = rare_charm_lfv_sm_branching_fraction(self.sm_inputs)

    def _base_diagnostics(self) -> dict[str, Any]:
        rs = self.anchor.rs_baseline
        return {
            "experimental_block": self.anchor.current_limit.block_key,
            "experimental_confidence_level": float(
                self.anchor.current_limit.confidence_level
            ),
            "lhcb_primary_90cl_limit": float(self.anchor.lhcb_primary_limit.value),
            "lhcb_primary_integrated_luminosity_fb_inv": float(
                self.anchor.lhcb_primary_limit.integrated_luminosity_fb_inv or 0.0
            ),
            "lhcb_primary_collision_energy_tev": (
                self.anchor.lhcb_primary_limit.collision_energy_tev
            ),
            "belle_2010_90cl_limit": float(self.anchor.belle_previous_limit.value),
            "babar_2012_90cl_limit": float(self.anchor.babar_previous_limit.value),
            "budget_source": _BUDGET_SOURCE,
            "budget_is_upper_limit": True,
            "sm_branching_fraction": 0.0,
            "sm_lfv_policy": (
                "D0 -> e mu is charged-LFV and has zero SM rate for catalog "
                "purposes; the HARD budget is applied to the pure-NP proxy."
            ),
            "charge_conjugate_modes_included": True,
            "rs_baseline_block": f"{rs.section_key}.{rs.block_key}",
            "rs_baseline_source": rs.source,
            "rs_baseline_year": rs.year,
            "rs_baseline_generic_kk_gluon_scale_tev": float(
                rs.generic_rs_kk_gluon_scale_tev
            ),
            "rs_baseline_composite_pseudo_goldstone_kk_gluon_scale_tev": float(
                rs.composite_pseudo_goldstone_kk_gluon_scale_tev
            ),
            "parametrization_citation": _PARAMETRIZATION_CITATION,
            "rs_matching_assumption": RARE_CHARM_LFV_DILEPTON_PROXY_V1,
            "needs_human_physics": _NEEDS_HUMAN_PHYSICS,
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
                    "non-vetoing only; no D0 -> e mu NP prediction was evaluated"
                ),
                **self._base_diagnostics(),
                **dict(diagnostics),
            },
        )

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        quark_input = point.get_extra(_REQUIRED_QUARK_EXTRA)
        lepton_input = point.get_extra(_REQUIRED_LEPTON_EXTRA)
        missing = [
            key
            for key, value in (
                (_REQUIRED_QUARK_EXTRA, quark_input),
                (_REQUIRED_LEPTON_EXTRA, lepton_input),
            )
            if value is None
        ]
        if missing:
            return self._unevaluated_result({"missing_extras": tuple(missing)})

        kk_ew_mass = point.get_extra(_OPTIONAL_EW_MASS_EXTRA)
        try:
            result = d0_emu_from_couplings(
                quark_input,
                lepton_input,
                m_kk_gev=None if kk_ew_mass is None else float(kk_ew_mass),
                inputs=self.sm_inputs,
            )
        except (AttributeError, KeyError, TypeError, ValueError) as exc:
            return self._unevaluated_result(
                {
                    "invalid_extra": (
                        _REQUIRED_QUARK_EXTRA,
                        _REQUIRED_LEPTON_EXTRA,
                    ),
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
                "Pure-NP BR(D0 -> e+- mu-+) bound using the shared rare-charm "
                "c->u l l C9/C10 Hamiltonian and an unequal-lepton D0 two-body "
                "rate. The e-mu lepton coupling is an explicit documented "
                "proxy and is flagged NEEDS-HUMAN-PHYSICS."
            ),
            diagnostics=diagnostics,
        )
