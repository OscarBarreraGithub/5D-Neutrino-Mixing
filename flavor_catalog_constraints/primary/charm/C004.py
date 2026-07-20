"""C004 - rare leptonic decay ``D0 -> mu+ mu-``.

Physics
-------
``D0 -> mu+ mu-`` is long-distance dominated in the Standard Model, mainly
through two-photon effects.  The experimentally clean catalog veto implemented
here is therefore the short-distance ``c -> u mu mu`` piece:

    BR_SD(D0 -> mu mu) <= BR_exp^upper.

The short-distance rate uses the shared ``c -> u l l`` core built for C004.
It follows the rare-charm Hamiltonian convention

    H_eff = -4 G_F/sqrt(2) lambda_b alpha/(4 pi) [C9 O9 + C10 O10 + ...],

with the leptonic D0 rate controlled by ``C10 - C10'``.  Long-distance SM
numbers from the sidecar are reported as diagnostics and are not subtracted
from the experimental upper limit.

Severity
--------
HARD.  The observed upper bound is a direct experimental limit on the total
branching fraction, and any RS short-distance contribution larger than the
catalogued PDG 90% CL bound would overfill it.

Catalog sidecar
---------------
``flavor_catalog/processes/charm/C004.yaml`` is the source of truth for the
PDG/CMS/LHCb branching-fraction limits and Burdman-Golowich-Hewett-Pakvasa
long-distance context.  Numeric values below are loaded through the scaffold
anchor loader, not hardcoded in this constraint.

NEEDS-HUMAN-PHYSICS
-------------------
The short-distance vector/axial RS contribution consumes Phase-3a
``rs_semileptonic_wilsons.c_to_u_ll`` C9/C10/C9'/C10' additively.  Charm
long-distance context and scalar/pseudoscalar matching remain diagnostic and
deferred.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Mapping

from flavor_catalog_constraints.anchors import Anchor, AnchorError, load_anchor, load_full_yaml
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.rare_charm_dilepton import (
    RARE_CHARM_RS_SEMILEPTONIC_VECTOR_MATCHING_STATUS_V1,
    RS_SEMILEPTONIC_MATCHING_ASSUMPTION_V1,
    d0_mumu_from_rs_semileptonic_wilsons,
    rare_charm_dilepton_default_sm_inputs,
    rare_charm_dilepton_sm_branching_fraction,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "charm"
_REQUIRED_EXTRA = "rs_semileptonic_wilsons"
_OPTIONAL_EW_MASS_EXTRA = "kk_ew_mass_gev"
_EXPECTED_UNITS = "branching fraction"

_CURRENT_LIMIT_CANDIDATES = ("canonical_current_limit",)
_CMS_LIMIT_CANDIDATES = ("cms_current_95cl_limit",)
_LHCB_LIMIT_CANDIDATES = ("lhcb_previous_limit",)
_LONG_DISTANCE_CANDIDATES = ("standard_model_long_distance_context",)
_BUDGET_SOURCE = (
    "flavor_catalog/processes/charm/C004.yaml canonical_current_limit "
    "(PDG Live/API S032.28, 90% CL)"
)
_PARAMETRIZATION_CITATION = (
    "Rare-charm c->u l l C9/C10 Hamiltonian; "
    "Burdman-Golowich-Hewett-Pakvasa hep-ph/0112235 long-distance caveat; "
    "Gisbert-Hiller-Suelmann JHEP 2024 rare-charm EFT context"
)
_NEEDS_HUMAN_PHYSICS = (
    "NEEDS-HUMAN-PHYSICS: C004 uses rigorous Phase-3a light-Z "
    "rs_semileptonic_wilsons.c_to_u_ll vector/axial C9/C10/C9'/C10' NP; "
    "long-distance charm, heavy-neutral/lepton-tower completion, Higgs/radion, "
    "scalar and pseudoscalar matching remain deferred."
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


@dataclass(frozen=True)
class LongDistanceContext:
    """Typed Burdman et al. long-distance context from the C004 sidecar."""

    block_key: str
    source: str | None
    year: int | None
    relation_to_d0_gammagamma_factor: float
    quoted_minimum_branching_fraction: float
    d0_gammagamma_vmd_estimate: float
    d0_gammagamma_vmd_uncertainty_plus: float
    d0_gammagamma_vmd_uncertainty_minus: float
    vmd_implied_branching_fraction: float
    vmd_implied_uncertainty_plus: float
    vmd_implied_uncertainty_minus: float
    units: str | None
    source_url: str | None
    snapshot_path: str | None


@dataclass(frozen=True)
class C004Anchor:
    """Typed C004 anchor: experimental limits, LD SM context, and budget."""

    current_limit: BranchingLimitAnchor
    cms_95_limit: BranchingLimitAnchor
    lhcb_previous_limit: BranchingLimitAnchor
    long_distance_context: LongDistanceContext

    @property
    def value(self) -> float:
        """Current PDG upper limit used as the experimental value."""
        return self.current_limit.value

    @property
    def budget(self) -> float:
        """HARD short-distance branching-fraction budget."""
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
            f"{process_id}: C004 anchor field {field_name!r}={value!r} is not numeric"
        ) from exc
    if not math.isfinite(out):
        raise AnchorError(f"{process_id}: C004 anchor field {field_name!r} is not finite")
    return out


def _positive_float(value: Any, *, process_id: str, field_name: str) -> float:
    out = _required_float(value, process_id=process_id, field_name=field_name)
    if out <= 0.0:
        raise AnchorError(f"{process_id}: C004 anchor field {field_name!r} must be positive")
    return out


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
    if scaffold_anchor.units != _EXPECTED_UNITS:
        raise AnchorError(
            f"{process_id}: {scaffold_anchor.block_key} must use units "
            f"{_EXPECTED_UNITS!r}, got {scaffold_anchor.units!r}"
        )
    value = _positive_float(
        scaffold_anchor.value,
        process_id=process_id,
        field_name=f"{scaffold_anchor.block_key}.value",
    )
    confidence_level = _positive_float(
        sub.get("confidence_level"),
        process_id=process_id,
        field_name=f"{scaffold_anchor.block_key}.confidence_level",
    )
    limit_type = str(sub.get("limit_type", "")).lower()
    if limit_type != "upper":
        raise AnchorError(
            f"{process_id}: {scaffold_anchor.block_key}.limit_type must be 'upper'"
        )
    return BranchingLimitAnchor(
        block_key=scaffold_anchor.block_key,
        source=_optional_str(scaffold_anchor.source),
        year=scaffold_anchor.year,
        value=float(value),
        confidence_level=float(confidence_level),
        limit_type=limit_type,
        units=scaffold_anchor.units,
        source_url=_optional_str(scaffold_anchor.source_url),
        snapshot_path=_optional_str(scaffold_anchor.snapshot_path),
    )


def _load_long_distance_context(process_id: str) -> LongDistanceContext:
    relation = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=_LONG_DISTANCE_CANDIDATES,
        value_key="relation_to_d0_gammagamma_factor",
    )
    quoted_minimum = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=_LONG_DISTANCE_CANDIDATES,
        value_key="quoted_minimum_branching_fraction",
    )
    vmd_gamma = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=_LONG_DISTANCE_CANDIDATES,
        value_key="d0_gammagamma_vmd_estimate",
        uncertainty_key="d0_gammagamma_vmd_uncertainty_plus",
    )
    sub = _pdg_subblock_for_anchor(relation, process_id=process_id)
    minus = _positive_float(
        sub.get("d0_gammagamma_vmd_uncertainty_minus"),
        process_id=process_id,
        field_name="standard_model_long_distance_context."
        "d0_gammagamma_vmd_uncertainty_minus",
    )
    factor = float(relation.value)
    gamma_vmd = float(vmd_gamma.value)
    plus = float(vmd_gamma.uncertainty or 0.0)
    if factor <= 0.0 or gamma_vmd <= 0.0 or quoted_minimum.value <= 0.0:
        raise AnchorError(f"{process_id}: C004 long-distance context must be positive")
    return LongDistanceContext(
        block_key=relation.block_key,
        source=_optional_str(relation.source),
        year=relation.year,
        relation_to_d0_gammagamma_factor=factor,
        quoted_minimum_branching_fraction=float(quoted_minimum.value),
        d0_gammagamma_vmd_estimate=gamma_vmd,
        d0_gammagamma_vmd_uncertainty_plus=plus,
        d0_gammagamma_vmd_uncertainty_minus=float(minus),
        vmd_implied_branching_fraction=float(factor * gamma_vmd),
        vmd_implied_uncertainty_plus=float(factor * plus),
        vmd_implied_uncertainty_minus=float(factor * minus),
        units=relation.units,
        source_url=_optional_str(relation.source_url),
        snapshot_path=_optional_str(relation.snapshot_path),
    )


def _load_c004_anchor(process_id: str) -> C004Anchor:
    anchor = C004Anchor(
        current_limit=_load_limit_anchor(
            _CURRENT_LIMIT_CANDIDATES,
            process_id=process_id,
        ),
        cms_95_limit=_load_limit_anchor(
            _CMS_LIMIT_CANDIDATES,
            process_id=process_id,
        ),
        lhcb_previous_limit=_load_limit_anchor(
            _LHCB_LIMIT_CANDIDATES,
            process_id=process_id,
        ),
        long_distance_context=_load_long_distance_context(process_id),
    )
    if anchor.budget <= 0.0:
        raise AnchorError(f"{process_id}: C004 HARD budget must be positive")
    return anchor


def _complex_wilsons(result: Any) -> dict[str, complex]:
    if result.wilsons is None:
        return {}
    return {key: complex(value) for key, value in result.wilsons.wilsons.items()}


@register
class Constraint:
    """Catalogued short-distance ``D0 -> mu+ mu-`` constraint (C004)."""

    process_id = "C004"
    severity = Severity.HARD
    observable = "BR(D0 -> mu+ mu-) short-distance"

    def __init__(self) -> None:
        self.anchor = _load_c004_anchor(self.process_id)
        self.sm_inputs = rare_charm_dilepton_default_sm_inputs()
        self.sm_result = rare_charm_dilepton_sm_branching_fraction("mu", self.sm_inputs)

    def _base_diagnostics(self) -> dict[str, Any]:
        ld = self.anchor.long_distance_context
        return {
            "experimental_block": self.anchor.current_limit.block_key,
            "experimental_confidence_level": float(
                self.anchor.current_limit.confidence_level
            ),
            "cms_95cl_limit": float(self.anchor.cms_95_limit.value),
            "lhcb_previous_90cl_limit": float(self.anchor.lhcb_previous_limit.value),
            "budget_source": _BUDGET_SOURCE,
            "budget_is_upper_limit": True,
            "sm_short_distance_branching_fraction": float(
                self.sm_result.branching_fraction
            ),
            "sm_long_distance_block": ld.block_key,
            "sm_long_distance_quoted_minimum": float(
                ld.quoted_minimum_branching_fraction
            ),
            "sm_long_distance_d0_gammagamma_factor": float(
                ld.relation_to_d0_gammagamma_factor
            ),
            "sm_long_distance_d0_gammagamma_vmd_estimate": float(
                ld.d0_gammagamma_vmd_estimate
            ),
            "sm_long_distance_vmd_implied_branching_fraction": float(
                ld.vmd_implied_branching_fraction
            ),
            "sm_long_distance_vmd_implied_uncertainty_plus": float(
                ld.vmd_implied_uncertainty_plus
            ),
            "sm_long_distance_vmd_implied_uncertainty_minus": float(
                ld.vmd_implied_uncertainty_minus
            ),
            "long_distance_not_subtracted": True,
            "charm_theory_caveat": (
                "D0 -> mu+ mu- is long-distance dominated; C004 constrains "
                "only the short-distance c->u mu mu component."
            ),
            "parametrization_citation": _PARAMETRIZATION_CITATION,
            "rs_matching_assumption": RS_SEMILEPTONIC_MATCHING_ASSUMPTION_V1,
            "rs_semileptonic_vector_matching_status": (
                RARE_CHARM_RS_SEMILEPTONIC_VECTOR_MATCHING_STATUS_V1
            ),
            "needs_human_physics": _NEEDS_HUMAN_PHYSICS,
        }

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
                    f"extra {_REQUIRED_EXTRA!r} absent; D0 -> mu+ mu- "
                    "short-distance constraint was not evaluated."
                ),
                diagnostics={
                    "evaluated": False,
                    "missing_extra": _REQUIRED_EXTRA,
                    "legacy_quark_mass_basis_couplings_present": (
                        point.get_extra("quark_mass_basis_couplings") is not None
                    ),
                    **self._base_diagnostics(),
                },
            )

        kk_ew_mass = point.get_extra(_OPTIONAL_EW_MASS_EXTRA)
        try:
            result = d0_mumu_from_rs_semileptonic_wilsons(
                rs_wilsons,
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
                    "D0 -> mu+ mu-"
                ),
                diagnostics={
                    "evaluated": False,
                    "invalid_extra": _REQUIRED_EXTRA,
                    "exception_type": type(exc).__name__,
                    "exception_message": str(exc),
                    **self._base_diagnostics(),
                },
            )
        predicted = float(result.branching_fraction)
        budget = float(self.anchor.budget)
        ratio = predicted / budget if budget > 0.0 else float("inf")

        diagnostics = dict(result.diagnostics)
        diagnostics.update(
            {
                **self._base_diagnostics(),
                "evaluated": True,
                "np_shift_branching_fraction": float(
                    result.np_shift_branching_fraction
                ),
                "short_distance_plus_vmd_long_distance": float(
                    predicted
                    + self.anchor.long_distance_context.vmd_implied_branching_fraction
                ),
                "short_distance_plus_vmd_ratio_to_limit": float(
                    (
                        predicted
                        + self.anchor.long_distance_context.vmd_implied_branching_fraction
                    )
                    / budget
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
            sm_prediction=float(result.sm_branching_fraction),
            experimental=float(self.anchor.value),
            ratio=float(ratio),
            budget=budget,
            notes=(
                "BR_SD(D0 -> mu+ mu-) uses the shared c->u l l short-distance "
                "C10-dominant formula. Phase-3a rs_semileptonic_wilsons.c_to_u_ll "
                "C9/C10/C9'/C10' enter additively with no old proxy prefactor "
                "or second M_KK suppression. The HARD ratio is the short-distance "
                "branching fraction over the C004 YAML PDG 90% CL upper limit; "
                "long-distance charm context is reported but not subtracted."
            ),
            diagnostics=diagnostics,
        )
