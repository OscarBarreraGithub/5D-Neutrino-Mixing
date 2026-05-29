"""C005 - rare leptonic decay ``D0 -> e+ e-``.

Physics
-------
``D0 -> e+ e-`` is the electron mode of the shared rare-charm
``c -> u l l`` machinery used by C004.  The Standard Model short-distance
piece is numerically zero in the current catalog implementation
(``C10_SM=0``), while the physical mode is long-distance dominated and
strongly helicity suppressed by ``m_e^2``:

    BR_SD(D0 -> e e) <= BR_exp^upper.

The short-distance rate is evaluated with the shared rare-charm Hamiltonian
convention

    H_eff = -4 G_F/sqrt(2) lambda_b alpha/(4 pi) [C9 O9 + C10 O10 + ...],

with the leptonic D0 rate controlled by ``C10 - C10'`` and the electron mass
selected from the common input bundle.

Severity
--------
HARD.  The observed upper bound is a direct experimental limit on the total
branching fraction, and any RS short-distance contribution larger than the
catalogued PDG/Belle 90% CL bound would overfill it.

Catalog sidecar
---------------
``flavor_catalog/processes/charm/C005.yaml`` is the source of truth for the
PDG/Belle/BABAR branching-fraction limits and rare-charm theory provenance.
Numeric values below are loaded through the scaffold anchor loader, not
hardcoded in this constraint.

NEEDS-HUMAN-PHYSICS
-------------------
The RS contribution reuses the shared documented Z/KK-penguin proxy because
the ``ParameterPoint`` does not provide the full electroweak KK/Z/Z',
charged-lepton, Higgs/radion, scalar, or pseudoscalar matching inputs needed
for a rigorous RS ``c -> u e e`` calculation.
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
)
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.rare_charm_dilepton import (
    RARE_CHARM_DILEPTON_RS_MATCHING_ASSUMPTION_V1,
    d0_ee_from_couplings,
    rare_charm_dilepton_default_sm_inputs,
    rare_charm_dilepton_sm_branching_fraction,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "charm"
_REQUIRED_EXTRA = "quark_mass_basis_couplings"
_OPTIONAL_EW_MASS_EXTRA = "kk_ew_mass_gev"
_EXPECTED_UNITS = "branching fraction"

_CURRENT_LIMIT_CANDIDATES = ("canonical_current_limit",)
_BELLE_LIMIT_CANDIDATES = ("belle_current_result",)
_BABAR_LIMIT_CANDIDATES = ("babar_independent_result",)
_THEORY_CONTEXT_SECTION = "theory_context"
_THEORY_CONTEXT_KEY = "rare_charm_model_dependence"
_PAPER_ERA_SECTION = "paper_era_reference"
_RS_BASELINE_KEY = "rs_baseline"
_BUDGET_SOURCE = (
    "flavor_catalog/processes/charm/C005.yaml canonical_current_limit "
    "(PDG Live/API S032.39, 90% CL; Belle 2010 current limit)"
)
_PARAMETRIZATION_CITATION = (
    "Rare-charm c->u l l C9/C10 Hamiltonian; "
    "Burdman-Golowich-Hewett-Pakvasa hep-ph/0112235 long-distance caveat"
)
_NEEDS_HUMAN_PHYSICS = (
    "NEEDS-HUMAN-PHYSICS: full RS electroweak KK/Z/Z', charged-lepton, "
    "Higgs/radion, scalar and pseudoscalar c->u e e matching is not "
    "available on ParameterPoint; v1 uses the documented Z/KK-penguin "
    "C9/C10 proxy."
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
class RareCharmTheoryContext:
    """Typed rare-charm theory provenance from the C005 sidecar."""

    section_key: str
    block_key: str
    source: str | None
    year: int | None
    topic: str | None
    source_url: str | None
    snapshot_path: str | None


@dataclass(frozen=True)
class RSBaselineContext:
    """Typed paper-era RS reference block carried for proxy provenance."""

    section_key: str
    block_key: str
    source: str | None
    year: int | None
    generic_rs_kk_gluon_scale_tev: float
    composite_pseudo_goldstone_kk_gluon_scale_tev: float
    source_url: str | None
    snapshot_path: str | None


@dataclass(frozen=True)
class C005Anchor:
    """Typed C005 anchor: experimental limits, theory context, and budget."""

    current_limit: BranchingLimitAnchor
    belle_current_limit: BranchingLimitAnchor
    babar_independent_limit: BranchingLimitAnchor
    theory_context: RareCharmTheoryContext
    rs_baseline: RSBaselineContext

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
            f"{process_id}: C005 anchor field {field_name!r}={value!r} is not numeric"
        ) from exc
    if not math.isfinite(out):
        raise AnchorError(f"{process_id}: C005 anchor field {field_name!r} is not finite")
    return out


def _positive_float(value: Any, *, process_id: str, field_name: str) -> float:
    out = _required_float(value, process_id=process_id, field_name=field_name)
    if out <= 0.0:
        raise AnchorError(f"{process_id}: C005 anchor field {field_name!r} must be positive")
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


def _load_theory_context(process_id: str) -> RareCharmTheoryContext:
    sub = _top_level_subblock(
        process_id,
        section_key=_THEORY_CONTEXT_SECTION,
        block_key=_THEORY_CONTEXT_KEY,
    )
    return RareCharmTheoryContext(
        section_key=_THEORY_CONTEXT_SECTION,
        block_key=_THEORY_CONTEXT_KEY,
        source=_optional_str(sub.get("source")),
        year=_optional_int(sub.get("year")),
        topic=_optional_str(sub.get("topic")),
        source_url=_optional_str(sub.get("source_url")),
        snapshot_path=_optional_str(sub.get("snapshot_path")),
    )


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


def _load_c005_anchor(process_id: str) -> C005Anchor:
    anchor = C005Anchor(
        current_limit=_load_limit_anchor(
            _CURRENT_LIMIT_CANDIDATES,
            process_id=process_id,
        ),
        belle_current_limit=_load_limit_anchor(
            _BELLE_LIMIT_CANDIDATES,
            process_id=process_id,
        ),
        babar_independent_limit=_load_limit_anchor(
            _BABAR_LIMIT_CANDIDATES,
            process_id=process_id,
        ),
        theory_context=_load_theory_context(process_id),
        rs_baseline=_load_rs_baseline(process_id),
    )
    if anchor.budget <= 0.0:
        raise AnchorError(f"{process_id}: C005 HARD budget must be positive")
    return anchor


def _complex_wilsons(result: Any) -> dict[str, complex]:
    if result.wilsons is None:
        return {}
    return {key: complex(value) for key, value in result.wilsons.wilsons.items()}


@register
class Constraint:
    """Catalogued short-distance ``D0 -> e+ e-`` constraint (C005)."""

    process_id = "C005"
    severity = Severity.HARD
    observable = "BR(D0 -> e+ e-) short-distance"

    def __init__(self) -> None:
        self.anchor = _load_c005_anchor(self.process_id)
        self.sm_inputs = rare_charm_dilepton_default_sm_inputs()
        self.sm_result = rare_charm_dilepton_sm_branching_fraction("e", self.sm_inputs)

    def _base_diagnostics(self) -> dict[str, Any]:
        theory = self.anchor.theory_context
        rs = self.anchor.rs_baseline
        electron_mass = float(self.sm_inputs.electron.mass_gev)
        muon_mass = float(self.sm_inputs.muon.mass_gev)
        return {
            "experimental_block": self.anchor.current_limit.block_key,
            "experimental_confidence_level": float(
                self.anchor.current_limit.confidence_level
            ),
            "belle_90cl_limit": float(self.anchor.belle_current_limit.value),
            "babar_90cl_limit": float(self.anchor.babar_independent_limit.value),
            "budget_source": _BUDGET_SOURCE,
            "budget_is_upper_limit": True,
            "sm_short_distance_branching_fraction": float(
                self.sm_result.branching_fraction
            ),
            "sm_short_distance_policy": (
                "No nonzero short-distance SM C10 anchor is provided in the "
                "C005 catalog sidecar; the clean constrained piece is evaluated "
                "with C10_SM=0 and long-distance rare-charm context is carried "
                "as provenance rather than subtracted."
            ),
            "electron_mass_gev": electron_mass,
            "muon_mass_gev": muon_mass,
            "electron_to_muon_helicity_suppression_m2": float(
                (electron_mass / muon_mass) ** 2
            ),
            "sm_long_distance_numeric_anchor_available": False,
            "long_distance_not_subtracted": True,
            "theory_context_block": f"{theory.section_key}.{theory.block_key}",
            "theory_context_source": theory.source,
            "theory_context_year": theory.year,
            "theory_context_topic": theory.topic,
            "rs_baseline_block": f"{rs.section_key}.{rs.block_key}",
            "rs_baseline_source": rs.source,
            "rs_baseline_year": rs.year,
            "rs_baseline_generic_kk_gluon_scale_tev": float(
                rs.generic_rs_kk_gluon_scale_tev
            ),
            "rs_baseline_composite_pseudo_goldstone_kk_gluon_scale_tev": float(
                rs.composite_pseudo_goldstone_kk_gluon_scale_tev
            ),
            "charm_theory_caveat": (
                "D0 -> e+ e- is long-distance dominated and additionally "
                "helicity suppressed; C005 constrains only the short-distance "
                "c->u e e component."
            ),
            "parametrization_citation": _PARAMETRIZATION_CITATION,
            "rs_matching_assumption": RARE_CHARM_DILEPTON_RS_MATCHING_ASSUMPTION_V1,
            "needs_human_physics": _NEEDS_HUMAN_PHYSICS,
        }

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        couplings = point.get_extra(_REQUIRED_EXTRA)
        if couplings is None:
            return ConstraintResult(
                process_id=self.process_id,
                severity=self.severity,
                passes=True,
                sm_prediction=float(self.sm_result.branching_fraction),
                experimental=float(self.anchor.value),
                budget=float(self.anchor.budget),
                notes=(
                    f"extra {_REQUIRED_EXTRA!r} absent; D0 -> e+ e- "
                    "short-distance constraint was not evaluated."
                ),
                diagnostics={
                    "missing_extra": _REQUIRED_EXTRA,
                    **self._base_diagnostics(),
                },
            )

        kk_ew_mass = point.get_extra(_OPTIONAL_EW_MASS_EXTRA)
        result = d0_ee_from_couplings(
            couplings,
            m_kk_gev=None if kk_ew_mass is None else float(kk_ew_mass),
            inputs=self.sm_inputs,
        )
        predicted = float(result.branching_fraction)
        budget = float(self.anchor.budget)
        ratio = predicted / budget if budget > 0.0 else float("inf")

        diagnostics = dict(result.diagnostics)
        diagnostics.update(
            {
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
            sm_prediction=float(result.sm_branching_fraction),
            experimental=float(self.anchor.value),
            ratio=float(ratio),
            budget=budget,
            notes=(
                "BR_SD(D0 -> e+ e-) uses the shared c->u l l short-distance "
                "C10-dominant formula with the electron mass, giving the "
                "expected m_e^2 helicity suppression. The RS contribution is "
                "the documented Z/KK-penguin C9/C10 proxy from mass-basis "
                "u-c couplings; full EW/lepton/scalar matching is marked "
                "NEEDS-HUMAN-PHYSICS. The HARD ratio is the short-distance "
                "branching fraction over the C005 YAML PDG/Belle 90% CL "
                "upper limit; long-distance charm context is reported but "
                "not subtracted."
            ),
            diagnostics=diagnostics,
        )
