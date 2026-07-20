"""B034 - penguin CP phase in ``B_s -> phi phi``.

Physics
-------
``B_s -> phi phi`` is a pure charmless ``b -> s sbar s`` penguin mode.  The
measured phase ``phi_s^{s sbar s}`` is sensitive to the same small ``B_s``
mixing phase convention used in cleaner modes, but its SM interpretation is
polluted by mode- and helicity-dependent hadronic penguin amplitudes and
strong phases.  A first-principles SM number would require a QCDF/SCET or
phenomenological nonleptonic amplitude framework.  A rigorous RS prediction
would additionally require ``Delta B = 1`` penguin matching, running, and
``B_s -> phi phi`` matrix elements.

This file is therefore an honest non-vetoing stub, like B033.  It loads the
measured LHCb-combination ``phi_s^{s sbar s}`` entry and companion
``|lambda|`` value from ``B034.yaml``.  It also records the PDG ``beta_s``
convention entry as a sign-convention cross-check, but does not use the clean
``B_s`` mixing machinery as a full nonleptonic prediction.

Severity
--------
INFO.  Both the SM hadronic-penguin side and the RS ``Delta B = 1`` penguin
side are flagged ``NEEDS-HUMAN-PHYSICS``.  The returned result is
non-vetoing and must not be used as a hard exclusion.

Catalog sidecar
---------------
``flavor_catalog/processes/beauty/B034.yaml`` stores the relevant PDG/LHCb
values inside list-shaped ``observables`` blocks.  Selected list entries are
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
from flavor_catalog_constraints.physics_adapters.charmless_b_phiphi import (
    BS_TO_PHI_PHI_NONLEPTONIC_STUB_MODEL_V1,
    BS_TO_PHI_PHI_RS_NEEDS_HUMAN_PHYSICS,
    BS_TO_PHI_PHI_SM_NEEDS_HUMAN_PHYSICS,
    BsPhiPhiPhaseRoomComparison,
    compare_bs_phiphi_phase_to_room,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "beauty"
_PDG_MEASUREMENT_BLOCK = "pdg_beta_s_b_to_ssbars"
_LHCB_2023_BLOCK = "lhcb2023_precision_and_combination"

_BETA_S_OBSERVABLE = "beta_s(b -> s sbar s)"
_PHI_S_OBSERVABLE = "phi_s^{s sbar s}, LHCb combined note"
_LAMBDA_OBSERVABLE = "|lambda|, LHCb combined note"
_LHCB_RUN12_PHI_S_OBSERVABLE = (
    "phi_s^{s sbar s}, LHCb Run 1 + Run 2 combination"
)
_LHCB_RUN12_LAMBDA_OBSERVABLE = "|lambda|, LHCb Run 1 + Run 2 combination"

_RAD_UNITS = "rad"
_BETA_S_UNITS = "10^-2 rad"
_DIMENSIONLESS_UNITS = "dimensionless"
_BETA_TO_PHI_SCALE = -2.0e-2
_BETA_TO_PHI_TOLERANCE_RAD = 1.0e-12

_BUDGET_SOURCE = (
    "flavor_catalog/processes/beauty/B034.yaml pdg_beta_s_b_to_ssbars; "
    "budget is the one-sigma experimental uncertainty on phi_s^{s sbar s}, "
    "reported as a non-vetoing INFO bookkeeping scale because no reliable SM "
    "penguin subtraction exists"
)
_PARAMETRIZATION_CITATION = (
    "PDG pdgLive B_s^0 listing S086 and LHCb Phys. Rev. Lett. 131 (2023) "
    "171802; phi_s^{s sbar s} is recorded as the LHCb Run 1 + Run 2 "
    "combination, with beta_s retained only as the PDG convention cross-check"
)
_NEEDS_HUMAN_PHYSICS = (
    BS_TO_PHI_PHI_SM_NEEDS_HUMAN_PHYSICS,
    BS_TO_PHI_PHI_RS_NEEDS_HUMAN_PHYSICS,
)


@dataclass(frozen=True)
class B034ObservableAnchor:
    """Typed PDG/LHCb observable entry from a B034 list block."""

    block_key: str
    parent_block_key: str
    observable: str
    latex: str | None
    source: str | None
    year: int | None
    value: float
    uncertainty: float | None
    stat_uncertainty: float | None
    syst_uncertainty: float | None
    units: str | None
    convention: str | None
    arxiv_id: str | None
    doi: str | None
    source_url: str | None
    access_date: str | None
    snapshot_path: str | None
    sha256: str | None


@dataclass(frozen=True)
class B034Anchor:
    """Typed B034 measurement bundle and non-vetoing room convention."""

    phi_s_sss: B034ObservableAnchor
    beta_s: B034ObservableAnchor
    lambda_abs: B034ObservableAnchor
    lhcb_run12_phi_s_sss: B034ObservableAnchor
    lhcb_run12_lambda_abs: B034ObservableAnchor
    comparison: BsPhiPhiPhaseRoomComparison
    phi_s_from_beta_convention: float
    phi_s_uncertainty_from_beta_convention: float

    @property
    def value(self) -> float:
        """Measured ``phi_s^{s sbar s}`` central value in radians."""
        return self.phi_s_sss.value

    @property
    def uncertainty(self) -> float:
        """Measured ``phi_s^{s sbar s}`` one-sigma uncertainty in radians."""
        assert self.phi_s_sss.uncertainty is not None
        return self.phi_s_sss.uncertainty

    @property
    def budget(self) -> float:
        """Non-vetoing one-sigma phase room used only for scalar bookkeeping."""
        return self.comparison.documented_np_room_abs


def _optional_str(value: Any) -> str | None:
    return None if value is None else str(value)


def _optional_float(value: Any, *, process_id: str, field_name: str) -> float | None:
    if value is None:
        return None
    return _required_float(value, process_id=process_id, field_name=field_name)


def _required_float(value: Any, *, process_id: str, field_name: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise AnchorError(
            f"{process_id}: B034 anchor field {field_name!r}={value!r} "
            "is not numeric"
        ) from exc
    if not math.isfinite(number):
        raise AnchorError(
            f"{process_id}: B034 anchor field {field_name!r}={value!r} "
            "is not finite"
        )
    return number


def _positive_float(value: Any, *, process_id: str, field_name: str) -> float:
    number = _required_float(value, process_id=process_id, field_name=field_name)
    if number <= 0.0:
        raise AnchorError(f"{process_id}: B034 anchor field {field_name!r} <= 0")
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
            f"expected {virtual_key!r} for B034 observable {observable!r}"
        )
    return scaffold_anchor, merged


def _load_b034_observable(
    pdg_block: Mapping[str, Any],
    block_key: str,
    observable: str,
    *,
    process_id: str,
    expected_units: str,
    require_uncertainty: bool = True,
) -> B034ObservableAnchor:
    anchor, merged = _scaffold_anchor_for_observable(
        pdg_block,
        block_key,
        observable,
        process_id=process_id,
    )
    if anchor.units != expected_units:
        raise AnchorError(
            f"{process_id}: {anchor.block_key} expected units "
            f"{expected_units!r}, got {anchor.units!r}"
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
    return B034ObservableAnchor(
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
        stat_uncertainty=_optional_float(
            merged.get("stat_uncertainty"),
            process_id=process_id,
            field_name=f"{anchor.block_key}.stat_uncertainty",
        ),
        syst_uncertainty=_optional_float(
            merged.get("syst_uncertainty"),
            process_id=process_id,
            field_name=f"{anchor.block_key}.syst_uncertainty",
        ),
        units=anchor.units,
        convention=_optional_str(merged.get("convention")),
        arxiv_id=_optional_str(merged.get("arxiv_id")),
        doi=_optional_str(merged.get("doi")),
        source_url=_optional_str(anchor.source_url),
        access_date=_optional_str(merged.get("access_date")),
        snapshot_path=_optional_str(anchor.snapshot_path),
        sha256=_optional_str(merged.get("sha256_of_text_snapshot")),
    )


def _load_b034_anchor(process_id: str) -> B034Anchor:
    pdg_block = _pdg_block(process_id)
    phi_s = _load_b034_observable(
        pdg_block,
        _PDG_MEASUREMENT_BLOCK,
        _PHI_S_OBSERVABLE,
        process_id=process_id,
        expected_units=_RAD_UNITS,
    )
    beta_s = _load_b034_observable(
        pdg_block,
        _PDG_MEASUREMENT_BLOCK,
        _BETA_S_OBSERVABLE,
        process_id=process_id,
        expected_units=_BETA_S_UNITS,
    )
    lambda_abs = _load_b034_observable(
        pdg_block,
        _PDG_MEASUREMENT_BLOCK,
        _LAMBDA_OBSERVABLE,
        process_id=process_id,
        expected_units=_DIMENSIONLESS_UNITS,
    )
    lhcb_phi_s = _load_b034_observable(
        pdg_block,
        _LHCB_2023_BLOCK,
        _LHCB_RUN12_PHI_S_OBSERVABLE,
        process_id=process_id,
        expected_units=_RAD_UNITS,
    )
    lhcb_lambda_abs = _load_b034_observable(
        pdg_block,
        _LHCB_2023_BLOCK,
        _LHCB_RUN12_LAMBDA_OBSERVABLE,
        process_id=process_id,
        expected_units=_DIMENSIONLESS_UNITS,
    )
    assert phi_s.uncertainty is not None
    assert beta_s.uncertainty is not None
    phi_s_from_beta = float(_BETA_TO_PHI_SCALE * beta_s.value)
    phi_s_uncertainty_from_beta = float(abs(_BETA_TO_PHI_SCALE) * beta_s.uncertainty)
    if not math.isclose(
        phi_s.value,
        phi_s_from_beta,
        rel_tol=0.0,
        abs_tol=_BETA_TO_PHI_TOLERANCE_RAD,
    ):
        raise AnchorError(
            f"{process_id}: phi_s={phi_s.value} does not match the PDG "
            f"-2 beta_s convention value {phi_s_from_beta}"
        )
    if not math.isclose(
        lhcb_phi_s.value,
        phi_s.value,
        rel_tol=0.0,
        abs_tol=_BETA_TO_PHI_TOLERANCE_RAD,
    ):
        raise AnchorError(
            f"{process_id}: LHCb Run 1 + Run 2 phi_s={lhcb_phi_s.value} does "
            f"not match canonical phi_s={phi_s.value}"
        )
    comparison = compare_bs_phiphi_phase_to_room(
        measured_phi_s=phi_s.value,
        experimental_uncertainty=phi_s.uncertainty,
        documented_np_room_abs=phi_s.uncertainty,
    )
    return B034Anchor(
        phi_s_sss=phi_s,
        beta_s=beta_s,
        lambda_abs=lambda_abs,
        lhcb_run12_phi_s_sss=lhcb_phi_s,
        lhcb_run12_lambda_abs=lhcb_lambda_abs,
        comparison=comparison,
        phi_s_from_beta_convention=phi_s_from_beta,
        phi_s_uncertainty_from_beta_convention=phi_s_uncertainty_from_beta,
    )


@register
class Constraint:
    """Catalogued non-vetoing ``phi_s^{s sbar s}`` penguin-CP stub."""

    process_id = "B034"
    severity = Severity.INFO
    observable = "phi_s^sss(Bs -> phi phi)"

    def __init__(self) -> None:
        self.anchor = _load_b034_anchor(self.process_id)

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        del point
        comparison = self.anchor.comparison
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
                "INFO-only B034 stub: loads the measured LHCb-combination "
                "phi_s^{s sbar s} in B_s -> phi phi and records the "
                "one-sigma experimental phase uncertainty as non-vetoing "
                "INFO room. No QCDF/SCET SM amplitude and no RS Delta B=1 "
                "penguin amplitude are computed; both missing pieces are "
                "flagged NEEDS-HUMAN-PHYSICS."
            ),
            diagnostics={
                "non_vetoing": True,
                "severity_policy": "INFO/non-vetoing even when advisory passes=False",
                "stub_model": BS_TO_PHI_PHI_NONLEPTONIC_STUB_MODEL_V1,
                "budget_source": _BUDGET_SOURCE,
                "budget_is_measurement_uncertainty": True,
                "parametrization_citation": _PARAMETRIZATION_CITATION,
                "needs_human_physics": _NEEDS_HUMAN_PHYSICS,
                "needs_human_physics_sm": BS_TO_PHI_PHI_SM_NEEDS_HUMAN_PHYSICS,
                "needs_human_physics_np": BS_TO_PHI_PHI_RS_NEEDS_HUMAN_PHYSICS,
                "no_hadronic_penguin_calculation": True,
                "no_qcdf_scet_prediction": True,
                "no_delta_b1_penguin_matching": True,
                "rs_delta_b1_penguin_amplitude_available": False,
                "full_sm_phi_s_sss_prediction_available": False,
                "bs_mixing_phase_reused": False,
                "bs_mixing_phase_used_as_full_phi_s_sss_prediction": False,
                "parameter_point_inputs_used": (),
                "hflav_lhcb_convention": self.anchor.phi_s_sss.convention,
                "experimental_block": self.anchor.phi_s_sss.parent_block_key,
                "lhcb_combination_block": self.anchor.lhcb_run12_phi_s_sss.parent_block_key,
                "phi_s_sss_value": float(comparison.measured_phi_s),
                "phi_s_sss_uncertainty": float(comparison.experimental_uncertainty),
                "phi_s_sss_abs": float(comparison.measured_abs),
                "phi_s_room_abs": float(comparison.documented_np_room_abs),
                "beta_s_value_10minus2_rad": float(self.anchor.beta_s.value),
                "beta_s_uncertainty_10minus2_rad": float(
                    self.anchor.beta_s.uncertainty or 0.0
                ),
                "phi_s_from_beta_convention_rad": float(
                    self.anchor.phi_s_from_beta_convention
                ),
                "phi_s_uncertainty_from_beta_convention_rad": float(
                    self.anchor.phi_s_uncertainty_from_beta_convention
                ),
                "lambda_abs_value": float(self.anchor.lambda_abs.value),
                "lambda_abs_uncertainty": float(
                    self.anchor.lambda_abs.uncertainty or 0.0
                ),
                "lhcb_run12_phi_s_value": float(
                    self.anchor.lhcb_run12_phi_s_sss.value
                ),
                "lhcb_run12_phi_s_uncertainty": float(
                    self.anchor.lhcb_run12_phi_s_sss.uncertainty or 0.0
                ),
                "lhcb_run12_lambda_abs_value": float(
                    self.anchor.lhcb_run12_lambda_abs.value
                ),
                "lhcb_run12_lambda_abs_uncertainty": float(
                    self.anchor.lhcb_run12_lambda_abs.uncertainty or 0.0
                ),
                "source_url": self.anchor.phi_s_sss.source_url,
                "lhcb_source_url": self.anchor.lhcb_run12_phi_s_sss.source_url,
                "snapshot_path": self.anchor.phi_s_sss.snapshot_path,
            },
        )
