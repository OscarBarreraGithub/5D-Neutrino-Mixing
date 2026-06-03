"""T010 - ``Z -> b bbar`` pole pseudo-observables ``R_b`` and ``A_b``.

Physics
-------
This constraint uses the shared Z-pole pseudo-observable machinery introduced
in ``quarkConstraints.zpole`` and reached only through the
``flavor_catalog_constraints.physics_adapters.zpole`` boundary.  Effective
chiral couplings are normalized as

    L_Z = g_Z Z_mu bbar gamma^mu (g_L P_L + g_R P_R) b,

with

    R_b = Gamma_b / Gamma_had,
    A_b = (|g_L|^2 - |g_R|^2) / (|g_L|^2 + |g_R|^2).

The SM-limit width radiator for the bottom channel is calibrated to the
LEP/SLC SM-fit ``R_b`` value in the YAML-referenced local source snapshot.

RS matching status
------------------
PARTIAL / NEEDS-HUMAN-PHYSICS.  The minimal-RS gauge neutral-current piece is
read from the Phase-3a ``rs_ew_couplings`` extra.  The classic ``Zbb``
fermion-KK/custodial/BKT completion is not included until Phase 6, so the
diagnostics keep an explicit partial-status flag.

Severity
--------
HARD.  ``R_b^0`` and ``A_b`` are observed Z-pole precision measurements; a
point-specific bottom-coupling shift must fit inside their LEP/SLC
uncertainty-aware budgets.  The legacy ``A_FB^{0,b}`` 2.8 sigma SM tension is
loaded and reported as context but is not used as the T010 veto scalar.

Catalog sidecar
---------------
``flavor_catalog/processes/top_higgs_ew/T010.yaml`` is the source of truth for
the measured anchors and provenance.  Its ``pdg_or_equivalent`` block is a
list, so selected entries are routed through the scaffold ``load_anchor`` path
and fail loudly on missing observables.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
import re
from typing import Any, Mapping, Sequence

from flavor_catalog_constraints import anchors as anchor_scaffold
from flavor_catalog_constraints.anchors import Anchor, AnchorError, load_anchor, load_full_yaml
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.zpole import (
    ZPoleQuarkObservables,
    zpole_default_sm_inputs,
    zpole_evaluate_quark,
    zpole_inputs_with_bottom_radiator,
    zpole_shifted_couplings,
    zpole_sm_couplings,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "top_higgs_ew"
_REQUIRED_EXTRA = "rs_ew_couplings"
_REPO_ROOT = Path(__file__).resolve().parents[3]

_RB_OBSERVABLE = "R_b^0"
_AFB_OBSERVABLE = "A_FB^{0,b}"
_AB_OBSERVABLE = "A_b"
_AFB_PULL_OBSERVABLE = "LEP/SLC final-combination pull for A_FB^{0,b}"
_FCC_PROJECTION_OBSERVABLE = "FCC-ee projected relative uncertainty for R_b and A_FB^b"
_VETO_OBSERVABLES = (_RB_OBSERVABLE, _AB_OBSERVABLE)

_SM_SNAPSHOT_LABELS = {
    _RB_OBSERVABLE: "R_b^0",
    _AFB_OBSERVABLE: "A_FB^(0,b)",
    _AB_OBSERVABLE: "A_b",
}
_NUMBER_RE = r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?"
_NEEDS_HUMAN_PHYSICS = (
    "PARTIAL/NEEDS-HUMAN-PHYSICS: Phase-3b uses rigorous minimal-RS gauge "
    "neutral-current z_delta_g_L/R_d[2,2] from rs_ew_couplings; the classic "
    "Zbb fermion-KK, custodial-representation, and brane-kinetic-term "
    "completion is deferred to Phase 6."
)
_UNEVALUATED_REASON = "rs_ew_couplings not provided"
_UNEVALUATED_NOTES = (
    "T010 not evaluated: rs_ew_couplings not provided on ParameterPoint. "
    "Result is non-vetoing and diagnostic-only; build points with "
    "build_from_rs_ew_inputs for the rigorous Phase-3b path."
)


@dataclass(frozen=True)
class SMFitValue:
    """SM-fit value parsed from the YAML-referenced LEP/SLC snapshot."""

    observable: str
    value: float
    uncertainty: float
    pull: float
    source_path: str
    source_url: str | None


@dataclass(frozen=True)
class ObservableBudget:
    """Uncertainty-aware budget for one Z-pole pseudo-observable."""

    observable: str
    experimental_sigma: float
    sm_fit_sigma: float
    combined_sigma: float
    source: str


@dataclass(frozen=True)
class T010Anchor:
    """Typed T010 anchor bundle."""

    r_b: Anchor
    a_fb: Anchor
    a_b: Anchor
    a_fb_pull: Anchor
    fcc_projection: Anchor
    sm_fit_values: Mapping[str, SMFitValue]
    budgets: Mapping[str, ObservableBudget]

    @property
    def value(self) -> float:
        return self.r_b.value

    @property
    def budget(self) -> float:
        return self.budgets[_RB_OBSERVABLE].combined_sigma


@dataclass(frozen=True)
class ObservablePull:
    """One scalar pull entering the two-observable T010 max-pull result."""

    observable: str
    predicted: float
    sm_prediction: float
    experimental: float
    budget: float
    pull: float
    sm_pull: float

    @property
    def abs_pull(self) -> float:
        return float(abs(self.pull))


def _required_float(value: Any, *, process_id: str, field_name: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise AnchorError(
            f"{process_id}: T010 field {field_name!r}={value!r} is not numeric"
        ) from exc
    if not math.isfinite(number):
        raise AnchorError(f"{process_id}: T010 field {field_name!r} is not finite")
    return number


def _positive_float(value: Any, *, process_id: str, field_name: str) -> float:
    number = _required_float(value, process_id=process_id, field_name=field_name)
    if number <= 0.0:
        raise AnchorError(f"{process_id}: T010 field {field_name!r} must be positive")
    return number


def _pdg_entries(process_id: str) -> Sequence[Mapping[str, Any]]:
    data = load_full_yaml(process_id, family=_FAMILY)
    entries = data.get("pdg_or_equivalent")
    if not isinstance(entries, list) or not entries:
        raise AnchorError(f"{process_id}: expected non-empty pdg_or_equivalent list")
    for index, entry in enumerate(entries):
        if not isinstance(entry, Mapping):
            raise AnchorError(
                f"{process_id}: pdg_or_equivalent[{index}] is not a mapping"
            )
    return entries


def _find_entry(
    process_id: str,
    observable: str,
) -> tuple[int, Mapping[str, Any]]:
    for index, entry in enumerate(_pdg_entries(process_id)):
        if entry.get("observable") == observable:
            return index, entry
    present = [str(entry.get("observable")) for entry in _pdg_entries(process_id)]
    raise AnchorError(
        f"{process_id}: observable {observable!r} not found in "
        f"pdg_or_equivalent list (present: {present})"
    )


def _load_scaffold_list_anchor(
    observable: str,
    *,
    process_id: str,
) -> Anchor:
    index, entry = _find_entry(process_id, observable)
    block_key = f"pdg_or_equivalent[{index}]"
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
        )
    finally:
        anchor_scaffold.load_pdg_block = original_load_pdg_block

    if scaffold_anchor.block_key != block_key:
        raise AnchorError(
            f"{process_id}: load_anchor selected {scaffold_anchor.block_key!r}, "
            f"expected {block_key!r} for T010 observable {observable!r}"
        )
    return scaffold_anchor


def _sm_snapshot_path(anchor: Anchor, *, process_id: str) -> Path:
    if anchor.snapshot_path is None:
        raise AnchorError(f"{process_id}: LEP/SLC SM-fit anchor lacks snapshot_path")
    path = _REPO_ROOT / anchor.snapshot_path
    if not path.is_file():
        raise AnchorError(f"{process_id}: LEP/SLC SM-fit snapshot not found: {path}")
    return path


def _parse_sm_fit_line(
    text: str,
    *,
    process_id: str,
    observable: str,
    source_path: Path,
    source_url: str | None,
) -> SMFitValue:
    label = _SM_SNAPSHOT_LABELS[observable]
    pattern = re.compile(
        rf"{re.escape(label)}\s*=\s*"
        rf"(?P<exp>{_NUMBER_RE})\s*\+/-\s*(?P<exp_unc>{_NUMBER_RE})"
        rf";\s*SM fit value shown there:\s*"
        rf"(?P<sm>{_NUMBER_RE})\s*\+/-\s*(?P<sm_unc>{_NUMBER_RE})"
        rf";\s*pull\s*(?P<pull>{_NUMBER_RE})",
    )
    match = pattern.search(text)
    if match is None:
        raise AnchorError(
            f"{process_id}: could not parse SM-fit value for {observable!r} "
            f"from {source_path}"
        )
    return SMFitValue(
        observable=observable,
        value=_required_float(
            match.group("sm"),
            process_id=process_id,
            field_name=f"{observable}.sm_fit_value",
        ),
        uncertainty=_positive_float(
            match.group("sm_unc"),
            process_id=process_id,
            field_name=f"{observable}.sm_fit_uncertainty",
        ),
        pull=_required_float(
            match.group("pull"),
            process_id=process_id,
            field_name=f"{observable}.sm_fit_pull",
        ),
        source_path=str(source_path.relative_to(_REPO_ROOT)),
        source_url=source_url,
    )


def _load_sm_fit_values(
    *,
    process_id: str,
    lep_slc_anchor: Anchor,
) -> Mapping[str, SMFitValue]:
    path = _sm_snapshot_path(lep_slc_anchor, process_id=process_id)
    text = path.read_text()
    return {
        observable: _parse_sm_fit_line(
            text,
            process_id=process_id,
            observable=observable,
            source_path=path,
            source_url=lep_slc_anchor.source_url,
        )
        for observable in (_RB_OBSERVABLE, _AFB_OBSERVABLE, _AB_OBSERVABLE)
    }


def _build_budget(
    *,
    process_id: str,
    observable: str,
    experimental: Anchor,
    sm_fit: SMFitValue,
) -> ObservableBudget:
    if experimental.uncertainty is None:
        raise AnchorError(f"{process_id}: {observable} must have an uncertainty")
    exp_sigma = _positive_float(
        experimental.uncertainty,
        process_id=process_id,
        field_name=f"{observable}.experimental_uncertainty",
    )
    sm_sigma = _positive_float(
        sm_fit.uncertainty,
        process_id=process_id,
        field_name=f"{observable}.sm_fit_uncertainty",
    )
    combined = math.sqrt(exp_sigma * exp_sigma + sm_sigma * sm_sigma)
    return ObservableBudget(
        observable=observable,
        experimental_sigma=float(exp_sigma),
        sm_fit_sigma=float(sm_sigma),
        combined_sigma=float(combined),
        source=(
            "flavor_catalog/processes/top_higgs_ew/T010.yaml measured anchor "
            f"+ {sm_fit.source_path} SM-fit uncertainty"
        ),
    )


def _load_t010_anchor(process_id: str) -> T010Anchor:
    r_b = _load_scaffold_list_anchor(_RB_OBSERVABLE, process_id=process_id)
    a_fb = _load_scaffold_list_anchor(_AFB_OBSERVABLE, process_id=process_id)
    a_b = _load_scaffold_list_anchor(_AB_OBSERVABLE, process_id=process_id)
    a_fb_pull = _load_scaffold_list_anchor(_AFB_PULL_OBSERVABLE, process_id=process_id)
    fcc_projection = _load_scaffold_list_anchor(
        _FCC_PROJECTION_OBSERVABLE,
        process_id=process_id,
    )
    sm_fit = _load_sm_fit_values(process_id=process_id, lep_slc_anchor=a_fb_pull)
    budgets = {
        _RB_OBSERVABLE: _build_budget(
            process_id=process_id,
            observable=_RB_OBSERVABLE,
            experimental=r_b,
            sm_fit=sm_fit[_RB_OBSERVABLE],
        ),
        _AB_OBSERVABLE: _build_budget(
            process_id=process_id,
            observable=_AB_OBSERVABLE,
            experimental=a_b,
            sm_fit=sm_fit[_AB_OBSERVABLE],
        ),
    }
    return T010Anchor(
        r_b=r_b,
        a_fb=a_fb,
        a_b=a_b,
        a_fb_pull=a_fb_pull,
        fcc_projection=fcc_projection,
        sm_fit_values=sm_fit,
        budgets=budgets,
    )


def _observable_pull(
    *,
    observable: str,
    predicted: float,
    sm_prediction: float,
    experimental: Anchor,
    budget: ObservableBudget,
) -> ObservablePull:
    pull = float((predicted - experimental.value) / budget.combined_sigma)
    sm_pull = float((sm_prediction - experimental.value) / budget.combined_sigma)
    return ObservablePull(
        observable=observable,
        predicted=float(predicted),
        sm_prediction=float(sm_prediction),
        experimental=float(experimental.value),
        budget=float(budget.combined_sigma),
        pull=pull,
        sm_pull=sm_pull,
    )


def _pulls_for_observables(
    *,
    anchor: T010Anchor,
    prediction: ZPoleQuarkObservables,
    sm_prediction: ZPoleQuarkObservables,
) -> Mapping[str, ObservablePull]:
    return {
        _RB_OBSERVABLE: _observable_pull(
            observable=_RB_OBSERVABLE,
            predicted=prediction.r_q,
            sm_prediction=sm_prediction.r_q,
            experimental=anchor.r_b,
            budget=anchor.budgets[_RB_OBSERVABLE],
        ),
        _AB_OBSERVABLE: _observable_pull(
            observable=_AB_OBSERVABLE,
            predicted=prediction.a_q,
            sm_prediction=sm_prediction.a_q,
            experimental=anchor.a_b,
            budget=anchor.budgets[_AB_OBSERVABLE],
        ),
    }


def _pull_diagnostics(pull: ObservablePull, budget: ObservableBudget) -> dict[str, float | str]:
    return {
        "predicted": float(pull.predicted),
        "sm_prediction": float(pull.sm_prediction),
        "experimental": float(pull.experimental),
        "budget": float(pull.budget),
        "pull": float(pull.pull),
        "abs_pull": float(pull.abs_pull),
        "sm_pull": float(pull.sm_pull),
        "experimental_sigma": float(budget.experimental_sigma),
        "sm_fit_sigma": float(budget.sm_fit_sigma),
        "budget_source": budget.source,
    }


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
    """Catalogued ``Z -> b bbar`` precision-electroweak constraint."""

    process_id = "T010"
    severity = Severity.HARD
    observable = "max Zbb pull from R_b^0 and A_b"

    def __init__(self) -> None:
        self.anchor = _load_t010_anchor(self.process_id)
        base_inputs = zpole_default_sm_inputs()
        self.sm_inputs = zpole_inputs_with_bottom_radiator(
            self.anchor.sm_fit_values[_RB_OBSERVABLE].value,
            base_inputs,
        )
        self.sm_observables = zpole_evaluate_quark("b", inputs=self.sm_inputs)

    def _unevaluated_result(self, diagnostics: Mapping[str, Any]) -> ConstraintResult:
        return ConstraintResult(
            process_id=self.process_id,
            severity=self.severity,
            passes=True,
            predicted=None,
            sm_prediction=float(self.sm_observables.r_q),
            experimental=float(self.anchor.value),
            ratio=None,
            budget=float(self.anchor.budget),
            notes=_UNEVALUATED_NOTES,
            diagnostics={
                "evaluated": False,
                "unevaluated_reason": _UNEVALUATED_REASON,
                "passes_semantics": (
                    "non-vetoing only; no Zbb Phase-3b NP prediction was evaluated"
                ),
                "required_parameter_point_extras": [_REQUIRED_EXTRA],
                "needs_human_physics": _NEEDS_HUMAN_PHYSICS,
                **dict(diagnostics),
            },
        )

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        rs_ew_couplings = point.get_extra(_REQUIRED_EXTRA)
        if rs_ew_couplings is None:
            return self._unevaluated_result(
                {
                    "missing_extra": _REQUIRED_EXTRA,
                    "legacy_quark_mass_basis_couplings_present": (
                        point.get_extra("quark_mass_basis_couplings") is not None
                    ),
                }
            )

        try:
            delta_g_left_b = _coupling_entry(
                rs_ew_couplings, "z_delta_g_L_d", 2, 2
            )
            delta_g_right_b = _coupling_entry(
                rs_ew_couplings, "z_delta_g_R_d", 2, 2
            )
            shifted_bottom = zpole_shifted_couplings(
                zpole_sm_couplings("b", self.sm_inputs),
                delta_g_left=delta_g_left_b,
                delta_g_right=delta_g_right_b,
            )
            prediction = zpole_evaluate_quark(
                "b", {"b": shifted_bottom}, inputs=self.sm_inputs
            )
        except (AttributeError, TypeError, ValueError) as exc:
            return self._unevaluated_result(
                {
                    "invalid_extra": _REQUIRED_EXTRA,
                    "exception_type": type(exc).__name__,
                    "exception": str(exc),
                }
            )

        pulls = _pulls_for_observables(
            anchor=self.anchor,
            prediction=prediction,
            sm_prediction=self.sm_observables,
        )
        selected = max(pulls.values(), key=lambda item: item.abs_pull)
        ratio = float(selected.abs_pull)
        passes = bool(ratio <= 1.0)

        diagnostics: dict[str, Any] = {
            "evaluated": True,
            "selected_observable": selected.observable,
            "observables": {
                key: _pull_diagnostics(pull, self.anchor.budgets[key])
                for key, pull in pulls.items()
            },
            "sm_validation": {
                "r_b_formula": float(self.sm_observables.r_q),
                "r_b_lep_slc_sm_fit": float(
                    self.anchor.sm_fit_values[_RB_OBSERVABLE].value
                ),
                "a_b_formula": float(self.sm_observables.a_q),
                "a_b_lep_slc_sm_fit": float(
                    self.anchor.sm_fit_values[_AB_OBSERVABLE].value
                ),
                "a_fb_formula": float(self.sm_observables.a_fb),
                "a_fb_lep_slc_sm_fit": float(
                    self.anchor.sm_fit_values[_AFB_OBSERVABLE].value
                ),
                "sin2_theta_eff": float(self.sm_inputs.sin2_theta_eff),
                "bottom_radiator": float(self.sm_inputs.radiator_for("b")),
            },
            "a_fb_legacy_context": {
                "experimental": float(self.anchor.a_fb.value),
                "experimental_uncertainty": float(self.anchor.a_fb.uncertainty),
                "lep_slc_sm_fit": float(
                    self.anchor.sm_fit_values[_AFB_OBSERVABLE].value
                ),
                "lep_slc_pull_sigma": float(self.anchor.a_fb_pull.value),
                "source_url": self.anchor.a_fb_pull.source_url,
            },
            "fcc_ee_projection_relative_uncertainty_percent": float(
                self.anchor.fcc_projection.value
            ),
            "needs_human_physics": _NEEDS_HUMAN_PHYSICS,
            "rs_matching_assumption": getattr(
                rs_ew_couplings, "matching_assumption", None
            ),
            "rs_ew_model_label": getattr(rs_ew_couplings, "model_label", None),
            "rs_ew_kk_mass_gev": float(getattr(rs_ew_couplings, "kk_ew_mass_gev")),
            "delta_g_left_b": complex(delta_g_left_b),
            "delta_g_right_b": complex(delta_g_right_b),
            "shifted_zbb_couplings": {
                "g_left": complex(shifted_bottom.g_left),
                "g_right": complex(shifted_bottom.g_right),
            },
            "required_parameter_point_extras": [_REQUIRED_EXTRA],
        }

        return ConstraintResult(
            process_id=self.process_id,
            severity=self.severity,
            passes=passes,
            predicted=float(selected.predicted),
            sm_prediction=float(selected.sm_prediction),
            experimental=float(selected.experimental),
            ratio=ratio,
            budget=float(selected.budget),
            notes=(
                "Max one-sigma pull over R_b^0 and A_b.  SM-limit "
                "pseudo-observables use effective Zbb couplings with the "
                "bottom width radiator calibrated to the YAML-referenced "
                "LEP/SLC SM-fit R_b value. Point-specific RS shifts use "
                "Phase-3a minimal-RS gauge neutral-current z_delta_g_L/R_d[2,2]; "
                "Zbb fermion-KK/custodial/BKT completion remains "
                "PARTIAL/NEEDS-HUMAN-PHYSICS."
            ),
            diagnostics=diagnostics,
        )
