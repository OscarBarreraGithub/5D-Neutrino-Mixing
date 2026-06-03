"""C007 - rare semileptonic decay ``D+ -> pi+ mu+ mu-``.

Physics
-------
This constraint reuses the shared C004 ``c -> u l l`` Wilson machinery through
``flavor_catalog_constraints.physics_adapters.rare_charm_semileptonic`` and
adds only the exclusive ``D -> pi`` short-distance form-factor normalization.
The smooth short-distance rate uses the vector/axial ``D -> pi mu mu``
distribution with ``C9 + C9'`` and ``C10 + C10'``.  Long-distance
``rho/omega/phi`` resonance amplitudes dominate the physical spectrum and are
not modeled.

Severity
--------
HARD.  The active budget is the C007 YAML PDG/LHCb 90% CL upper limit.  The
constrained object here is explicitly a smooth full-kinematic-q2 proxy: the
catalogued LHCb SM short-distance scale is added as a tiny incoherent SM floor,
while resonance long-distance contributions and experimental dimuon-window
acceptances are documented and not applied.

Catalog sidecar
---------------
``flavor_catalog/processes/charm/C007.yaml`` is the source of truth for the
PDG/LHCb limit, LHCb short-distance SM scale, search scope, theory context,
and RS-baseline provenance.  Numeric anchors are loaded through the scaffold
loader, not hardcoded in this constraint.

NEEDS-HUMAN-PHYSICS
-------------------
The short-distance vector/axial RS contribution consumes Phase-3a
``rs_semileptonic_wilsons.c_to_u_ll`` C9/C10/C9'/C10' additively.  The HARD
ratio remains a smooth full-q2 proxy comparison, not a LHCb windowed recast,
because long-distance resonances, scalar/tensor terms, and exact dimuon-window
acceptance remain diagnostic and deferred.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
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
    RARE_CHARM_RS_SEMILEPTONIC_VECTOR_MATCHING_STATUS_V1,
    RS_SEMILEPTONIC_MATCHING_ASSUMPTION_V1,
)
from flavor_catalog_constraints.physics_adapters.rare_charm_semileptonic import (
    RARE_CHARM_DTOPI_MUMU_PARAMETRIZATION_CITATION,
    RARE_CHARM_DTOPI_MUMU_RESONANCE_LIMITATION_V1,
    RareCharmDToPiFormFactorInputs,
    dplus_piplus_mumu_from_rs_semileptonic_wilsons,
    dplus_piplus_mumu_sm,
    rare_charm_dtopi_mumu_default_inputs,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "charm"
_REQUIRED_EXTRA = "rs_semileptonic_wilsons"
_OPTIONAL_EW_MASS_EXTRA = "kk_ew_mass_gev"
_EXPECTED_BRANCHING_UNITS = "branching fraction"

_CURRENT_LIMIT_CANDIDATES = ("canonical_current_limit",)
_LHCB_CURRENT_CANDIDATES = ("lhcb_current_result",)
_SEARCH_SCOPE_CANDIDATES = ("lhcb_2021_search_scope",)
_SM_SD_SCALE_CANDIDATES = ("lhcb_2021_short_distance_sm_scale",)
_PREVIOUS_LIMIT_CANDIDATES = ("lhcb_previous_result",)
_FF_FPLUS0_CANDIDATES = ("dtopi_form_factor_fplus0",)
_FF_POLE_MASS_CANDIDATES = ("dtopi_form_factor_pole_mass",)
_FF_FPLUS_SHAPE_A_CANDIDATES = ("dtopi_form_factor_fplus_shape_a",)
_FF_FZERO_SHAPE_B_CANDIDATES = ("dtopi_form_factor_fzero_shape_b",)
_THEORY_CONTEXT_KEY = "theory_context"
_RS_BASELINE_KEY = "rs_baseline"
_BUDGET_SOURCE = (
    "flavor_catalog/processes/charm/C007.yaml canonical_current_limit "
    "(PDG Live/API S031.42, LHCb 2021, 90% CL); compared to a smooth "
    "full-q2 short-distance proxy with no LHCb window/acceptance mask"
)
_FULL_Q2_PROXY_TREATMENT = (
    "full_kinematic_q2_smooth_short_distance_proxy_no_lhcb_nonresonant_window_"
    "or_resonance_acceptance"
)
_NEEDS_HUMAN_PHYSICS = (
    "NEEDS-HUMAN-PHYSICS: C007 uses rigorous Phase-3a light-Z "
    "rs_semileptonic_wilsons.c_to_u_ll vector/axial C9/C10/C9'/C10' NP; "
    "long-distance resonance amplitudes, heavy-neutral/lepton-tower completion, "
    "scalar/tensor terms and exact dimuon-window acceptance remain deferred."
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
    companion_95cl_limit: float | None = None
    integrated_luminosity_fb_inv: float | None = None
    dataset_year: int | None = None
    collision_energy_tev: float | None = None


@dataclass(frozen=True)
class NumericContextAnchor:
    """Typed numeric context entry loaded through ``load_anchor``."""

    block_key: str
    source: str | None
    year: int | None
    observable: str | None
    value: float
    units: str | None
    qualifier: str | None
    source_url: str | None
    snapshot_path: str | None


@dataclass(frozen=True)
class RareCharmTheoryContext:
    """Typed rare-charm theory provenance from the C007 sidecar."""

    block_key: str
    source: str | None
    year: int | None
    observable: str | None
    wilson_coefficients: tuple[str, ...]
    long_distance_note: str | None
    source_url: str | None
    snapshot_path: str | None


@dataclass(frozen=True)
class RSBaselineContext:
    """Typed paper-era RS reference block carried for proxy provenance."""

    block_key: str
    source: str | None
    year: int | None
    generic_rs_kk_gluon_scale_tev: float
    composite_pseudo_goldstone_kk_gluon_scale_tev: float
    source_url: str | None
    snapshot_path: str | None


@dataclass(frozen=True)
class DToPiFormFactorContext:
    """Typed YAML-loaded exclusive ``D -> pi`` form-factor inputs."""

    fplus0: Anchor
    pole_mass: Anchor
    fplus_shape_a: Anchor
    fzero_shape_b: Anchor

    @property
    def normalization_uncertainty_fraction(self) -> float:
        return float(2.0 * self.fplus0.uncertainty / self.fplus0.value)

    @property
    def source_summary(self) -> str:
        return (
            f"{self.fplus0.source}: f_+(0)={self.fplus0.value}"
            f"({self.fplus0.uncertainty}), m_pole={self.pole_mass.value}"
            f"({self.pole_mass.uncertainty}) GeV, a={self.fplus_shape_a.value}"
            f"({self.fplus_shape_a.uncertainty}), b={self.fzero_shape_b.value}"
            f"({self.fzero_shape_b.uncertainty}); snapshot="
            f"{self.fplus0.snapshot_path}"
        )

    def to_inputs(self) -> RareCharmDToPiFormFactorInputs:
        return RareCharmDToPiFormFactorInputs(
            fplus_0=float(self.fplus0.value),
            fplus_shape_a=float(self.fplus_shape_a.value),
            fzero_shape_b=float(self.fzero_shape_b.value),
            pole_mass_gev=float(self.pole_mass.value),
            source=self.source_summary,
        )


@dataclass(frozen=True)
class C007Anchor:
    """Typed C007 anchor: limits, SD SM scale, theory context, and budget."""

    current_limit: BranchingLimitAnchor
    lhcb_current_limit: BranchingLimitAnchor
    previous_nonresonant_limit: BranchingLimitAnchor
    search_scope: NumericContextAnchor
    short_distance_sm_scale: NumericContextAnchor
    theory_context: RareCharmTheoryContext
    rs_baseline: RSBaselineContext
    form_factor: DToPiFormFactorContext

    @property
    def value(self) -> float:
        """Current PDG/LHCb upper limit used as the experimental value."""
        return self.current_limit.value

    @property
    def sm_value(self) -> float:
        """LHCb-quoted order-of-magnitude SM short-distance scale."""
        return self.short_distance_sm_scale.value

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


def _optional_float(value: Any, *, process_id: str, field_name: str) -> float | None:
    if value is None:
        return None
    return _required_float(value, process_id=process_id, field_name=field_name)


def _required_float(value: Any, *, process_id: str, field_name: str) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError) as exc:
        raise AnchorError(
            f"{process_id}: C007 anchor field {field_name!r}={value!r} is not numeric"
        ) from exc
    if not math.isfinite(out):
        raise AnchorError(f"{process_id}: C007 anchor field {field_name!r} is not finite")
    return out


def _positive_float(value: Any, *, process_id: str, field_name: str) -> float:
    out = _required_float(value, process_id=process_id, field_name=field_name)
    if out <= 0.0:
        raise AnchorError(f"{process_id}: C007 anchor field {field_name!r} must be positive")
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


def _pdg_subblock(process_id: str, block_key: str) -> Mapping[str, Any]:
    data = load_full_yaml(process_id, family=_FAMILY)
    pdg_block = data.get("pdg_or_equivalent")
    if not isinstance(pdg_block, Mapping):
        raise AnchorError(f"{process_id}: expected mapping-shaped pdg_or_equivalent")
    sub = pdg_block.get(block_key)
    if not isinstance(sub, Mapping):
        present = sorted(str(key) for key in pdg_block)
        raise AnchorError(
            f"{process_id}: pdg_or_equivalent.{block_key} missing or not a mapping "
            f"(present keys: {present})"
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
    if scaffold_anchor.units != _EXPECTED_BRANCHING_UNITS:
        raise AnchorError(
            f"{process_id}: {scaffold_anchor.block_key} must use units "
            f"{_EXPECTED_BRANCHING_UNITS!r}, got {scaffold_anchor.units!r}"
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
        companion_95cl_limit=_optional_float(
            sub.get("companion_95cl_limit"),
            process_id=process_id,
            field_name=f"{scaffold_anchor.block_key}.companion_95cl_limit",
        ),
        integrated_luminosity_fb_inv=_optional_float(
            sub.get("integrated_luminosity_fb_inv"),
            process_id=process_id,
            field_name=f"{scaffold_anchor.block_key}.integrated_luminosity_fb_inv",
        ),
        dataset_year=_optional_int(sub.get("dataset_year")),
        collision_energy_tev=_optional_float(
            sub.get("collision_energy_TeV"),
            process_id=process_id,
            field_name=f"{scaffold_anchor.block_key}.collision_energy_TeV",
        ),
    )


def _load_numeric_context(
    candidates: tuple[str, ...],
    *,
    process_id: str,
    expected_units: str,
) -> NumericContextAnchor:
    anchor = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=candidates,
    )
    sub = _pdg_subblock_for_anchor(anchor, process_id=process_id)
    if anchor.units != expected_units:
        raise AnchorError(
            f"{process_id}: {anchor.block_key} must use units "
            f"{expected_units!r}, got {anchor.units!r}"
        )
    value = _positive_float(
        anchor.value,
        process_id=process_id,
        field_name=f"{anchor.block_key}.value",
    )
    return NumericContextAnchor(
        block_key=anchor.block_key,
        source=_optional_str(anchor.source),
        year=anchor.year,
        observable=_optional_str(anchor.observable),
        value=float(value),
        units=anchor.units,
        qualifier=_optional_str(sub.get("qualifier")),
        source_url=_optional_str(anchor.source_url),
        snapshot_path=_optional_str(anchor.snapshot_path),
    )


def _load_theory_context(process_id: str) -> RareCharmTheoryContext:
    sub = _pdg_subblock(process_id, _THEORY_CONTEXT_KEY)
    wilsons = sub.get("wilson_coefficients", ())
    if not isinstance(wilsons, (list, tuple)):
        raise AnchorError(f"{process_id}: theory_context.wilson_coefficients invalid")
    return RareCharmTheoryContext(
        block_key=_THEORY_CONTEXT_KEY,
        source=_optional_str(sub.get("source")),
        year=_optional_int(sub.get("year")),
        observable=_optional_str(sub.get("observable")),
        wilson_coefficients=tuple(str(item) for item in wilsons),
        long_distance_note=_optional_str(sub.get("long_distance_note")),
        source_url=_optional_str(sub.get("source_url")),
        snapshot_path=_optional_str(sub.get("snapshot_path")),
    )


def _load_rs_baseline(process_id: str) -> RSBaselineContext:
    sub = _pdg_subblock(process_id, _RS_BASELINE_KEY)
    return RSBaselineContext(
        block_key=_RS_BASELINE_KEY,
        source=_optional_str(sub.get("source")),
        year=_optional_int(sub.get("year")),
        generic_rs_kk_gluon_scale_tev=_positive_float(
            sub.get("generic_rs_kk_gluon_scale_TeV"),
            process_id=process_id,
            field_name="rs_baseline.generic_rs_kk_gluon_scale_TeV",
        ),
        composite_pseudo_goldstone_kk_gluon_scale_tev=_positive_float(
            sub.get("composite_pseudo_goldstone_kk_gluon_scale_TeV"),
            process_id=process_id,
            field_name="rs_baseline.composite_pseudo_goldstone_kk_gluon_scale_TeV",
        ),
        source_url=_optional_str(sub.get("source_url")),
        snapshot_path=_optional_str(sub.get("snapshot_path")),
    )


def _load_form_factor_value(
    candidates: tuple[str, ...],
    *,
    process_id: str,
    expected_units: str,
) -> Anchor:
    anchor = load_anchor(process_id, family=_FAMILY, candidates=candidates)
    if anchor.units != expected_units:
        raise AnchorError(
            f"{process_id}: {anchor.block_key} must use units "
            f"{expected_units!r}, got {anchor.units!r}"
        )
    _positive_float(
        anchor.value,
        process_id=process_id,
        field_name=f"{anchor.block_key}.value",
    )
    if anchor.uncertainty is None:
        raise AnchorError(f"{process_id}: {anchor.block_key}.uncertainty is required")
    _positive_float(
        anchor.uncertainty,
        process_id=process_id,
        field_name=f"{anchor.block_key}.uncertainty",
    )
    return anchor


def _load_form_factor_context(process_id: str) -> DToPiFormFactorContext:
    return DToPiFormFactorContext(
        fplus0=_load_form_factor_value(
            _FF_FPLUS0_CANDIDATES,
            process_id=process_id,
            expected_units="dimensionless",
        ),
        pole_mass=_load_form_factor_value(
            _FF_POLE_MASS_CANDIDATES,
            process_id=process_id,
            expected_units="GeV",
        ),
        fplus_shape_a=_load_form_factor_value(
            _FF_FPLUS_SHAPE_A_CANDIDATES,
            process_id=process_id,
            expected_units="dimensionless",
        ),
        fzero_shape_b=_load_form_factor_value(
            _FF_FZERO_SHAPE_B_CANDIDATES,
            process_id=process_id,
            expected_units="dimensionless",
        ),
    )


def _load_c007_anchor(process_id: str) -> C007Anchor:
    anchor = C007Anchor(
        current_limit=_load_limit_anchor(
            _CURRENT_LIMIT_CANDIDATES,
            process_id=process_id,
        ),
        lhcb_current_limit=_load_limit_anchor(
            _LHCB_CURRENT_CANDIDATES,
            process_id=process_id,
        ),
        previous_nonresonant_limit=_load_limit_anchor(
            _PREVIOUS_LIMIT_CANDIDATES,
            process_id=process_id,
        ),
        search_scope=_load_numeric_context(
            _SEARCH_SCOPE_CANDIDATES,
            process_id=process_id,
            expected_units="decay modes",
        ),
        short_distance_sm_scale=_load_numeric_context(
            _SM_SD_SCALE_CANDIDATES,
            process_id=process_id,
            expected_units="branching fraction, order of magnitude",
        ),
        theory_context=_load_theory_context(process_id),
        rs_baseline=_load_rs_baseline(process_id),
        form_factor=_load_form_factor_context(process_id),
    )
    if anchor.budget <= 0.0 or anchor.sm_value <= 0.0:
        raise AnchorError(f"{process_id}: C007 budget and SM SD scale must be positive")
    return anchor


def _complex_wilsons(result: Any) -> dict[str, complex]:
    if result.wilsons is None:
        return {}
    return {key: complex(value) for key, value in result.wilsons.wilsons.items()}


@register
class Constraint:
    """Catalogued short-distance ``D+ -> pi+ mu+ mu-`` constraint (C007)."""

    process_id = "C007"
    severity = Severity.HARD
    observable = "BR(D+ -> pi+ mu+ mu-) smooth full-q2 proxy"

    def __init__(self) -> None:
        self.anchor = _load_c007_anchor(self.process_id)
        base_inputs = rare_charm_dtopi_mumu_default_inputs()
        self.sd_inputs = replace(
            base_inputs,
            form_factor=self.anchor.form_factor.to_inputs(),
        )
        self.smooth_sm_result = dplus_piplus_mumu_sm(self.sd_inputs)

    def _base_diagnostics(self) -> dict[str, Any]:
        theory = self.anchor.theory_context
        rs = self.anchor.rs_baseline
        ff = self.anchor.form_factor
        return {
            "experimental_block": self.anchor.current_limit.block_key,
            "experimental_confidence_level": float(
                self.anchor.current_limit.confidence_level
            ),
            "lhcb_current_95cl_limit": float(
                self.anchor.lhcb_current_limit.companion_95cl_limit or 0.0
            ),
            "lhcb_current_integrated_luminosity_fb_inv": float(
                self.anchor.lhcb_current_limit.integrated_luminosity_fb_inv or 0.0
            ),
            "lhcb_previous_nonresonant_90cl_limit": float(
                self.anchor.previous_nonresonant_limit.value
            ),
            "lhcb_previous_nonresonant_95cl_limit": float(
                self.anchor.previous_nonresonant_limit.companion_95cl_limit or 0.0
            ),
            "lhcb_2021_search_scope_decay_modes": float(self.anchor.search_scope.value),
            "catalog_sm_short_distance_scale": float(self.anchor.sm_value),
            "catalog_sm_short_distance_scale_qualifier": (
                self.anchor.short_distance_sm_scale.qualifier
            ),
            "smooth_formula_sm_branching_fraction": float(
                self.smooth_sm_result.branching_fraction
            ),
            "smooth_formula_sm_policy": (
                "The shared rare_charm_dilepton core has no nonzero SM C9/C10 "
                "input for C007; the LHCb YAML O(1e-12) short-distance SM "
                "scale is carried as the SM prediction and added as a tiny "
                "incoherent floor to the NP proxy."
            ),
            "budget_source": _BUDGET_SOURCE,
            "budget_is_upper_limit": True,
            "full_q2_proxy_constrained": True,
            "q2_treatment": _FULL_Q2_PROXY_TREATMENT,
            "lhcb_nonresonant_dimuon_window_applied": False,
            "long_distance_resonance_dominated": True,
            "long_distance_not_subtracted": True,
            "resonance_window_recast_available": False,
            "budget_form_factor_uncertainty_policy": (
                "f_+(0) normalization uncertainty is diagnostic-only; the HARD "
                "budget remains the YAML 90% CL upper limit."
            ),
            "form_factor_fplus0_anchor_block": ff.fplus0.block_key,
            "form_factor_fplus0_uncertainty": float(ff.fplus0.uncertainty),
            "form_factor_pole_mass_anchor_block": ff.pole_mass.block_key,
            "form_factor_pole_mass_uncertainty_gev": float(ff.pole_mass.uncertainty),
            "form_factor_fplus_shape_a_anchor_block": ff.fplus_shape_a.block_key,
            "form_factor_fplus_shape_a_uncertainty": float(
                ff.fplus_shape_a.uncertainty
            ),
            "form_factor_fzero_shape_b_anchor_block": ff.fzero_shape_b.block_key,
            "form_factor_fzero_shape_b_uncertainty": float(
                ff.fzero_shape_b.uncertainty
            ),
            "form_factor_normalization_uncertainty_fraction": (
                ff.normalization_uncertainty_fraction
            ),
            "form_factor_normalization_uncertainty_percent": (
                100.0 * ff.normalization_uncertainty_fraction
            ),
            "theory_context_block": theory.block_key,
            "theory_context_source": theory.source,
            "theory_context_year": theory.year,
            "theory_context_observable": theory.observable,
            "theory_context_wilson_coefficients": theory.wilson_coefficients,
            "theory_context_long_distance_note": theory.long_distance_note,
            "rs_baseline_block": rs.block_key,
            "rs_baseline_source": rs.source,
            "rs_baseline_year": rs.year,
            "rs_baseline_generic_kk_gluon_scale_tev": float(
                rs.generic_rs_kk_gluon_scale_tev
            ),
            "rs_baseline_composite_pseudo_goldstone_kk_gluon_scale_tev": float(
                rs.composite_pseudo_goldstone_kk_gluon_scale_tev
            ),
            "parametrization_citation": RARE_CHARM_DTOPI_MUMU_PARAMETRIZATION_CITATION,
            "resonance_limitation": RARE_CHARM_DTOPI_MUMU_RESONANCE_LIMITATION_V1,
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
                sm_prediction=float(self.anchor.sm_value),
                experimental=float(self.anchor.value),
                budget=float(self.anchor.budget),
                notes=(
                    f"extra {_REQUIRED_EXTRA!r} absent; D+ -> pi+ mu+ mu- "
                    "smooth full-q2 proxy constraint was not evaluated."
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
            result = dplus_piplus_mumu_from_rs_semileptonic_wilsons(
                rs_wilsons,
                lepton="mu",
                matching_scale_gev=None if kk_ew_mass is None else float(kk_ew_mass),
                inputs=self.sd_inputs,
            )
        except Exception as exc:  # noqa: BLE001 - constraints degrade cleanly
            return ConstraintResult(
                process_id=self.process_id,
                severity=self.severity,
                passes=True,
                sm_prediction=float(self.anchor.sm_value),
                experimental=float(self.anchor.value),
                budget=float(self.anchor.budget),
                notes=(
                    "NOT EVALUATED - invalid rs_semileptonic_wilsons for "
                    "D+ -> pi+ mu+ mu-"
                ),
                diagnostics={
                    "evaluated": False,
                    "invalid_extra": _REQUIRED_EXTRA,
                    "exception_type": type(exc).__name__,
                    "exception_message": str(exc),
                    **self._base_diagnostics(),
                },
            )
        smooth_short_distance = float(result.branching_fraction)
        predicted = float(self.anchor.sm_value + smooth_short_distance)
        budget = float(self.anchor.budget)
        ratio = predicted / budget if budget > 0.0 else float("inf")

        diagnostics = dict(result.diagnostics)
        diagnostics.update(
            {
                **self._base_diagnostics(),
                "evaluated": True,
                "smooth_full_q2_np_proxy_branching_fraction": (
                    smooth_short_distance
                ),
                "catalog_sm_short_distance_scale_added_incoherently": True,
                "predicted_full_q2_proxy_total": float(predicted),
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
            sm_prediction=float(self.anchor.sm_value),
            experimental=float(self.anchor.value),
            ratio=float(ratio),
            budget=budget,
            notes=(
                "BR_SD(D+ -> pi+ mu+ mu-) consumes Phase-3a "
                "rs_semileptonic_wilsons.c_to_u_ll C9/C10/C9'/C10' additively "
                "and adds the existing exclusive D->pi form-factor rate. The "
                "HARD ratio is the smooth full-q2 proxy prediction over the "
                "C007 YAML PDG/LHCb 90% CL upper limit. Long-distance "
                "rho/omega/phi resonance physics and exact LHCb dimuon-window "
                "acceptance are documented but not applied."
            ),
            diagnostics=diagnostics,
        )
