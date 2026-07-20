"""B033 - mixing-induced CP asymmetry in ``B0 -> phi K_S``.

Physics
-------
``S_phiK_S`` is a penguin-dominated, time-dependent CP observable in the
``b -> s sbar s`` channel.  It probes the same clean ``B_d`` mixing phase
that appears in charmonium modes, but its Standard-Model interpretation is
polluted by mode-dependent hadronic penguin amplitudes and strong phases.  A
first-principles SM number would require a QCDF/SCET or phenomenological
nonleptonic amplitude treatment.  A rigorous RS prediction would require
``Delta B = 1`` penguin matching, running, and ``B -> phi K`` matrix elements.

This file is therefore an honest non-vetoing stub.  It loads the HFLAV Summer
2025 ``S_phiK_S`` measurement and the local HFLAV ``sin(2 beta)`` comparison
reference from ``B033.yaml``.  It does not import or reuse B002's clean
``B_d`` mixing-phase machinery as a full ``S_phiK_S`` prediction.

Severity
--------
INFO.  Both the SM hadronic-penguin side and the RS ``Delta B = 1`` penguin
side are flagged ``NEEDS-HUMAN-PHYSICS``.  The returned result is
non-vetoing and must not be used as a hard exclusion.

Catalog sidecar
---------------
``flavor_catalog/processes/beauty/B033.yaml`` stores the relevant HFLAV values
inside list-shaped ``observables`` blocks.  Selected list entries are
virtualized as one-entry mapping blocks and routed through the scaffold
``load_anchor`` path, then checked with typed, loud validation.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Mapping

from flavor_catalog_constraints import anchors as anchor_scaffold
from flavor_catalog_constraints.anchors import Anchor, AnchorError, load_anchor, load_full_yaml
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.charmless_b_phiks import (
    B_TO_PHI_KS_NONLEPTONIC_STUB_MODEL_V1,
    B_TO_PHI_KS_RS_NEEDS_HUMAN_PHYSICS,
    B_TO_PHI_KS_SM_NEEDS_HUMAN_PHYSICS,
    SphiKsReferenceComparison,
    compare_sphiks_to_sin2beta_reference,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "beauty"
_MEASUREMENT_BLOCK = "canonical_hflav_phiK0"
_REFERENCE_BLOCK = "hflav_comparison_sin2beta"

_S_PHI_KS_OBSERVABLE = "sin(2 beta_eff)(phi K0)"
_C_PHI_K0_OBSERVABLE = "C_CP(phi K0)"
_SIN2BETA_ALL_CHARMONIUM_OBSERVABLE = "sin(2 beta), all charmonium"
_SIN2BETA_JPSI_KS_OBSERVABLE = "sin(2 beta), J/psi K_S mode"
_DELTA_S_OBSERVABLE = (
    "difference sin(2 beta_eff)(phi K0) - sin(2 beta all charmonium)"
)

_EXPECTED_UNITS = "dimensionless"
_BUDGET_SOURCE = (
    "flavor_catalog/processes/beauty/B033.yaml "
    "canonical_hflav_phiK0 + hflav_comparison_sin2beta; "
    "budget is the non-vetoing quadrature scale for the measured Delta S "
    "reference comparison"
)
_PARAMETRIZATION_CITATION = (
    "HFLAV CPV & Unitarity Triangle Parameters Summer 2025; "
    "S_phiK_S is recorded as sin(2 beta_eff)(phi K0), with sin(2 beta) "
    "kept only as the clean charmonium reference"
)
_NEEDS_HUMAN_PHYSICS = (
    B_TO_PHI_KS_SM_NEEDS_HUMAN_PHYSICS,
    B_TO_PHI_KS_RS_NEEDS_HUMAN_PHYSICS,
)


@dataclass(frozen=True)
class HflavObservableAnchor:
    """Typed HFLAV observable entry from a B033 list block."""

    block_key: str
    parent_block_key: str
    observable: str
    latex: str | None
    source: str | None
    year: int | None
    value: float
    uncertainty: float | None
    units: str | None
    chi2: float | None
    dof: int | None
    confidence_level: float | None
    pull_sigma: float | None
    stat_only_uncertainty: float | None
    uncertainty_note: str | None
    convention: str | None
    source_url: str | None
    access_date: str | None
    snapshot_path: str | None
    sha256: str | None


@dataclass(frozen=True)
class B033Anchor:
    """Typed B033 HFLAV measurement, reference, and Delta S bookkeeping."""

    s_phi_ks: HflavObservableAnchor
    c_phi_k0: HflavObservableAnchor
    sin2beta_all_charmonium: HflavObservableAnchor
    sin2beta_jpsi_ks: HflavObservableAnchor
    delta_s_yaml: HflavObservableAnchor
    comparison: SphiKsReferenceComparison

    @property
    def value(self) -> float:
        """HFLAV ``S_phiK_S`` central value."""
        return self.s_phi_ks.value

    @property
    def uncertainty(self) -> float:
        """HFLAV ``S_phiK_S`` one-sigma uncertainty."""
        assert self.s_phi_ks.uncertainty is not None
        return self.s_phi_ks.uncertainty

    @property
    def sm_value(self) -> float:
        """HFLAV ``sin(2 beta)`` reference, not a full SM ``S_phiK_S`` prediction."""
        return self.sin2beta_all_charmonium.value

    @property
    def sm_uncertainty(self) -> float:
        """HFLAV ``sin(2 beta)`` reference uncertainty."""
        assert self.sin2beta_all_charmonium.uncertainty is not None
        return self.sin2beta_all_charmonium.uncertainty

    @property
    def budget(self) -> float:
        """Non-vetoing ``Delta S`` reference-comparison uncertainty scale."""
        return self.comparison.delta_s_uncertainty


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
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise AnchorError(
            f"{process_id}: B033 anchor field {field_name!r}={value!r} "
            "is not numeric"
        ) from exc
    if not math.isfinite(number):
        raise AnchorError(
            f"{process_id}: B033 anchor field {field_name!r}={value!r} "
            "is not finite"
        )
    return number


def _positive_float(value: Any, *, process_id: str, field_name: str) -> float:
    number = _required_float(value, process_id=process_id, field_name=field_name)
    if number <= 0.0:
        raise AnchorError(f"{process_id}: B033 anchor field {field_name!r} <= 0")
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


def _mapping_block(
    pdg_block: Mapping[str, Any],
    block_key: str,
    *,
    process_id: str,
) -> Mapping[str, Any]:
    block = pdg_block.get(block_key)
    if not isinstance(block, Mapping):
        raise AnchorError(
            f"{process_id}: pdg_or_equivalent.{block_key} is not a mapping"
        )
    return block


def _observable_entries(
    block: Mapping[str, Any],
    *,
    process_id: str,
    block_key: str,
) -> tuple[Mapping[str, Any], ...]:
    raw = block.get("observables")
    if not isinstance(raw, list) or not raw:
        raise AnchorError(f"{process_id}: {block_key}.observables is not a list")
    entries: list[Mapping[str, Any]] = []
    for index, entry in enumerate(raw):
        if not isinstance(entry, Mapping):
            raise AnchorError(
                f"{process_id}: {block_key}.observables[{index}] is not a mapping "
                f"(got {type(entry).__name__})"
            )
        entries.append(entry)
    return tuple(entries)


def _entry_name(
    entry: Mapping[str, Any],
    *,
    process_id: str,
    block_key: str,
    index: int,
) -> str:
    name = entry.get("name")
    if not isinstance(name, str) or not name:
        raise AnchorError(
            f"{process_id}: {block_key}.observables[{index}] has no non-empty 'name'"
        )
    return name


def _merged_observable_entry(
    pdg_block: Mapping[str, Any],
    block_key: str,
    observable: str,
    *,
    process_id: str,
) -> tuple[int, Mapping[str, Any]]:
    block = _mapping_block(pdg_block, block_key, process_id=process_id)
    entries = _observable_entries(block, process_id=process_id, block_key=block_key)
    for index, entry in enumerate(entries):
        name = _entry_name(
            entry,
            process_id=process_id,
            block_key=block_key,
            index=index,
        )
        if name == observable:
            merged = {key: value for key, value in block.items() if key != "observables"}
            merged.update(entry)
            merged["observable"] = name
            return index, merged
    raise AnchorError(
        f"{process_id}: {block_key}.observables has no observable {observable!r}"
    )


def _scaffold_anchor_for_observable(
    pdg_block: Mapping[str, Any],
    block_key: str,
    observable: str,
    *,
    process_id: str,
) -> tuple[Anchor, Mapping[str, Any]]:
    index, merged = _merged_observable_entry(
        pdg_block,
        block_key,
        observable,
        process_id=process_id,
    )
    virtual_key = f"{block_key}.observables[{index}]"
    virtual_block = {virtual_key: dict(merged)}
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
            candidates=(virtual_key,),
        )
    finally:
        anchor_scaffold.load_pdg_block = original_load_pdg_block

    if scaffold_anchor.block_key != virtual_key:
        raise AnchorError(
            f"{process_id}: load_anchor selected {scaffold_anchor.block_key!r}, "
            f"expected {virtual_key!r} for B033 observable {observable!r}"
        )
    return scaffold_anchor, merged


def _load_hflav_observable(
    pdg_block: Mapping[str, Any],
    block_key: str,
    observable: str,
    *,
    process_id: str,
    require_uncertainty: bool = True,
) -> HflavObservableAnchor:
    anchor, merged = _scaffold_anchor_for_observable(
        pdg_block,
        block_key,
        observable,
        process_id=process_id,
    )
    if anchor.units != _EXPECTED_UNITS:
        raise AnchorError(
            f"{process_id}: {anchor.block_key} expected units "
            f"{_EXPECTED_UNITS!r}, got {anchor.units!r}"
        )
    uncertainty = (
        None
        if anchor.uncertainty is None
        else _positive_float(
            anchor.uncertainty,
            process_id=process_id,
            field_name=f"{anchor.block_key}.uncertainty",
        )
    )
    if require_uncertainty and uncertainty is None:
        raise AnchorError(f"{process_id}: {anchor.block_key} lacks uncertainty")
    return HflavObservableAnchor(
        block_key=anchor.block_key,
        parent_block_key=block_key,
        observable=str(anchor.observable),
        latex=_optional_str(merged.get("latex")),
        source=_optional_str(anchor.source),
        year=anchor.year,
        value=_required_float(
            anchor.value,
            process_id=process_id,
            field_name=f"{anchor.block_key}.value",
        ),
        uncertainty=uncertainty,
        units=anchor.units,
        chi2=_optional_float(
            merged.get("chi2"),
            process_id=process_id,
            field_name=f"{anchor.block_key}.chi2",
        ),
        dof=_optional_int(merged.get("dof")),
        confidence_level=_optional_float(
            merged.get("confidence_level"),
            process_id=process_id,
            field_name=f"{anchor.block_key}.confidence_level",
        ),
        pull_sigma=_optional_float(
            merged.get("pull_sigma"),
            process_id=process_id,
            field_name=f"{anchor.block_key}.pull_sigma",
        ),
        stat_only_uncertainty=_optional_float(
            merged.get("stat_only_uncertainty"),
            process_id=process_id,
            field_name=f"{anchor.block_key}.stat_only_uncertainty",
        ),
        uncertainty_note=_optional_str(merged.get("uncertainty_note")),
        convention=_optional_str(merged.get("convention")),
        source_url=_optional_str(anchor.source_url),
        access_date=_optional_str(merged.get("access_date")),
        snapshot_path=_optional_str(anchor.snapshot_path),
        sha256=_optional_str(merged.get("sha256_of_text_snapshot")),
    )


def _load_b033_anchor(process_id: str) -> B033Anchor:
    pdg_block = _pdg_block(process_id)
    s_phi_ks = _load_hflav_observable(
        pdg_block,
        _MEASUREMENT_BLOCK,
        _S_PHI_KS_OBSERVABLE,
        process_id=process_id,
    )
    c_phi_k0 = _load_hflav_observable(
        pdg_block,
        _MEASUREMENT_BLOCK,
        _C_PHI_K0_OBSERVABLE,
        process_id=process_id,
    )
    sin2beta_all = _load_hflav_observable(
        pdg_block,
        _REFERENCE_BLOCK,
        _SIN2BETA_ALL_CHARMONIUM_OBSERVABLE,
        process_id=process_id,
    )
    sin2beta_jpsi = _load_hflav_observable(
        pdg_block,
        _REFERENCE_BLOCK,
        _SIN2BETA_JPSI_KS_OBSERVABLE,
        process_id=process_id,
    )
    delta_s_yaml = _load_hflav_observable(
        pdg_block,
        _REFERENCE_BLOCK,
        _DELTA_S_OBSERVABLE,
        process_id=process_id,
        require_uncertainty=False,
    )
    assert s_phi_ks.uncertainty is not None
    assert sin2beta_all.uncertainty is not None
    comparison = compare_sphiks_to_sin2beta_reference(
        measured_s_phi_ks=s_phi_ks.value,
        measured_uncertainty=s_phi_ks.uncertainty,
        sin2beta_reference=sin2beta_all.value,
        sin2beta_uncertainty=sin2beta_all.uncertainty,
    )
    if not math.isclose(
        delta_s_yaml.value,
        comparison.delta_s,
        rel_tol=0.0,
        abs_tol=1.0e-12,
    ):
        raise AnchorError(
            f"{process_id}: YAML Delta S={delta_s_yaml.value} does not match "
            f"S_phiK_S - sin2beta={comparison.delta_s}"
        )
    if comparison.delta_s_uncertainty <= 0.0:
        raise AnchorError(f"{process_id}: Delta S reference scale must be positive")
    return B033Anchor(
        s_phi_ks=s_phi_ks,
        c_phi_k0=c_phi_k0,
        sin2beta_all_charmonium=sin2beta_all,
        sin2beta_jpsi_ks=sin2beta_jpsi,
        delta_s_yaml=delta_s_yaml,
        comparison=comparison,
    )


@register
class Constraint:
    """Catalogued non-vetoing ``S_phiK_S`` penguin-CP stub (process_id B033)."""

    process_id = "B033"
    severity = Severity.INFO
    observable = "S_phiK_S"

    def __init__(self) -> None:
        self.anchor = _load_b033_anchor(self.process_id)
        self.sm_value = self.anchor.sm_value

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        del point
        comparison = self.anchor.comparison
        return ConstraintResult(
            process_id=self.process_id,
            severity=self.severity,
            passes=bool(comparison.passes),
            predicted=None,
            sm_prediction=float(self.anchor.sm_value),
            experimental=float(self.anchor.value),
            ratio=float(comparison.ratio_to_reference_uncertainty),
            budget=float(comparison.delta_s_uncertainty),
            notes=(
                "INFO-only B033 stub: loads HFLAV S_phiK_S and reports "
                "Delta S relative to the HFLAV sin(2 beta) charmonium "
                "reference. The sm_prediction field is that reference only, "
                "not a QCDF/SCET SM S_phiK_S prediction. No RS Delta B=1 "
                "penguin amplitude is computed; both missing pieces are "
                "flagged NEEDS-HUMAN-PHYSICS."
            ),
            diagnostics={
                "non_vetoing": True,
                "stub_model": B_TO_PHI_KS_NONLEPTONIC_STUB_MODEL_V1,
                "budget_source": _BUDGET_SOURCE,
                "parametrization_citation": _PARAMETRIZATION_CITATION,
                "needs_human_physics": _NEEDS_HUMAN_PHYSICS,
                "needs_human_physics_sm": B_TO_PHI_KS_SM_NEEDS_HUMAN_PHYSICS,
                "needs_human_physics_np": B_TO_PHI_KS_RS_NEEDS_HUMAN_PHYSICS,
                "no_hadronic_penguin_calculation": True,
                "no_qcdf_scet_prediction": True,
                "no_delta_b1_penguin_matching": True,
                "rs_delta_b1_penguin_amplitude_available": False,
                "full_sm_sphi_ks_prediction_available": False,
                "sin2beta_reference_is_not_full_sm_prediction": True,
                "sm_prediction_field_semantics": (
                    "HFLAV sin(2 beta), all charmonium reference only; not a "
                    "full SM S_phiK_S prediction."
                ),
                "b002_clean_mixing_phase_reused": False,
                "bd_mixing_phase_used_as_full_sphiks_prediction": False,
                "parameter_point_inputs_used": (),
                "hflav_convention": self.anchor.s_phi_ks.convention,
                "experimental_block": self.anchor.s_phi_ks.parent_block_key,
                "reference_block": self.anchor.sin2beta_all_charmonium.parent_block_key,
                "s_phi_ks_value": float(comparison.measured_s_phi_ks),
                "s_phi_ks_uncertainty": float(comparison.measured_uncertainty),
                "sin2beta_reference_value": float(comparison.sin2beta_reference),
                "sin2beta_reference_uncertainty": float(
                    comparison.sin2beta_uncertainty
                ),
                "delta_s_value": float(comparison.delta_s),
                "delta_s_yaml_central": float(self.anchor.delta_s_yaml.value),
                "delta_s_uncertainty": float(comparison.delta_s_uncertainty),
                "delta_s_reference_ratio": float(
                    comparison.ratio_to_reference_uncertainty
                ),
                "c_phi_k0_value": float(self.anchor.c_phi_k0.value),
                "c_phi_k0_uncertainty": float(self.anchor.c_phi_k0.uncertainty or 0.0),
                "sin2beta_jpsi_ks_value": float(self.anchor.sin2beta_jpsi_ks.value),
                "sin2beta_jpsi_ks_uncertainty": float(
                    self.anchor.sin2beta_jpsi_ks.uncertainty or 0.0
                ),
                "s_phi_ks_chi2": self.anchor.s_phi_ks.chi2,
                "s_phi_ks_dof": self.anchor.s_phi_ks.dof,
                "s_phi_ks_confidence_level": self.anchor.s_phi_ks.confidence_level,
                "s_phi_ks_pull_sigma": self.anchor.s_phi_ks.pull_sigma,
                "source_url": self.anchor.s_phi_ks.source_url,
                "snapshot_path": self.anchor.s_phi_ks.snapshot_path,
            },
        )
