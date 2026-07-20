"""T012 - ``Z -> c cbar`` pole pseudo-observables ``R_c`` and ``A_c``.

Physics
-------
This constraint reuses the shared Z-pole pseudo-observable machinery introduced
for T010 in ``quarkConstraints.zpole`` and reached only through the
``flavor_catalog_constraints.physics_adapters.zpole_charm`` boundary.
Effective chiral couplings are normalized as

    L_Z = g_Z Z_mu cbar gamma^mu (g_L P_L + g_R P_R) c,

with

    R_c = Gamma_c / Gamma_had,
    A_c = (|g_L|^2 - |g_R|^2) / (|g_L|^2 + |g_R|^2).

The SM-limit charm width radiator is calibrated to the YAML ``R_c^0`` anchor;
the independent chiral-coupling formula then gives the expected ``A_c`` value.

RS matching status
------------------
The minimal-RS gauge neutral-current piece is read from the Phase-3a
``rs_ew_couplings`` extra.  The point-specific ``Zcc`` prediction uses
``z_delta_g_L/R_u[1,1]`` directly, with no quark-overlap proxy.

Severity
--------
HARD.  ``R_c^0`` and ``A_c`` are observed LEP/SLC/PDG Z-pole measurements; a
point-specific charm-coupling shift must fit inside the one-sigma budgets loaded
from ``T012.yaml``.  ``A_FB^{0,c}`` is loaded and reported as context but is not
the T012 veto scalar requested here.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Mapping, Sequence

from flavor_catalog_constraints import anchors as anchor_scaffold
from flavor_catalog_constraints.anchors import Anchor, AnchorError, load_anchor, load_full_yaml
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.zpole_charm import (
    ZPoleQuarkObservables,
    zpole_default_sm_inputs,
    zpole_evaluate_quark,
    zpole_inputs_with_charm_radiator,
    zpole_shifted_couplings,
    zpole_sm_couplings,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "top_higgs_ew"
_REQUIRED_EXTRA = "rs_ew_couplings"

_RC_OBSERVABLE = "R_c^0"
_AFB_OBSERVABLE = "A_FB^{0,c}"
_AC_OBSERVABLE = "A_c"
_VETO_OBSERVABLES = (_RC_OBSERVABLE, _AC_OBSERVABLE)

_UNEVALUATED_REASON = "rs_ew_couplings not provided"
_UNEVALUATED_NOTES = (
    "T012 not evaluated: rs_ew_couplings not provided on ParameterPoint. "
    "Result is non-vetoing and diagnostic-only; build points with "
    "build_from_rs_ew_inputs for the rigorous Phase-3b path."
)


@dataclass(frozen=True)
class ObservableBudget:
    """One-sigma budget for one charm Z-pole pseudo-observable."""

    observable: str
    experimental_sigma: float
    sm_reference_sigma: float
    combined_sigma: float
    source: str


@dataclass(frozen=True)
class T012Anchor:
    """Typed T012 anchor bundle."""

    r_c: Anchor
    a_fb: Anchor
    a_c: Anchor
    budgets: Mapping[str, ObservableBudget]

    @property
    def value(self) -> float:
        return self.r_c.value

    @property
    def budget(self) -> float:
        return self.budgets[_RC_OBSERVABLE].combined_sigma


@dataclass(frozen=True)
class ObservablePull:
    """One scalar pull entering the T012 max-pull result."""

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
            f"{process_id}: T012 field {field_name!r}={value!r} is not numeric"
        ) from exc
    if not math.isfinite(number):
        raise AnchorError(f"{process_id}: T012 field {field_name!r} is not finite")
    return number


def _positive_float(value: Any, *, process_id: str, field_name: str) -> float:
    number = _required_float(value, process_id=process_id, field_name=field_name)
    if number <= 0.0:
        raise AnchorError(f"{process_id}: T012 field {field_name!r} must be positive")
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
            f"expected {block_key!r} for T012 observable {observable!r}"
        )
    return scaffold_anchor


def _build_budget(
    *,
    process_id: str,
    observable: str,
    experimental: Anchor,
) -> ObservableBudget:
    if experimental.uncertainty is None:
        raise AnchorError(f"{process_id}: {observable} must have an uncertainty")
    exp_sigma = _positive_float(
        experimental.uncertainty,
        process_id=process_id,
        field_name=f"{observable}.experimental_uncertainty",
    )
    return ObservableBudget(
        observable=observable,
        experimental_sigma=float(exp_sigma),
        sm_reference_sigma=0.0,
        combined_sigma=float(exp_sigma),
        source=(
            "flavor_catalog/processes/top_higgs_ew/T012.yaml measured "
            f"{observable} uncertainty; no separate SM-fit uncertainty block "
            "is present in the sidecar"
        ),
    )


def _load_t012_anchor(process_id: str) -> T012Anchor:
    r_c = _load_scaffold_list_anchor(_RC_OBSERVABLE, process_id=process_id)
    a_fb = _load_scaffold_list_anchor(_AFB_OBSERVABLE, process_id=process_id)
    a_c = _load_scaffold_list_anchor(_AC_OBSERVABLE, process_id=process_id)
    budgets = {
        _RC_OBSERVABLE: _build_budget(
            process_id=process_id,
            observable=_RC_OBSERVABLE,
            experimental=r_c,
        ),
        _AC_OBSERVABLE: _build_budget(
            process_id=process_id,
            observable=_AC_OBSERVABLE,
            experimental=a_c,
        ),
    }
    return T012Anchor(r_c=r_c, a_fb=a_fb, a_c=a_c, budgets=budgets)


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
    anchor: T012Anchor,
    prediction: ZPoleQuarkObservables,
    sm_prediction: ZPoleQuarkObservables,
) -> Mapping[str, ObservablePull]:
    return {
        _RC_OBSERVABLE: _observable_pull(
            observable=_RC_OBSERVABLE,
            predicted=prediction.r_q,
            sm_prediction=sm_prediction.r_q,
            experimental=anchor.r_c,
            budget=anchor.budgets[_RC_OBSERVABLE],
        ),
        _AC_OBSERVABLE: _observable_pull(
            observable=_AC_OBSERVABLE,
            predicted=prediction.a_q,
            sm_prediction=sm_prediction.a_q,
            experimental=anchor.a_c,
            budget=anchor.budgets[_AC_OBSERVABLE],
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
        "sm_reference_sigma": float(budget.sm_reference_sigma),
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
    """Catalogued ``Z -> c cbar`` precision-electroweak constraint."""

    process_id = "T012"
    severity = Severity.HARD
    observable = "max Zcc pull from R_c^0 and A_c"

    def __init__(self) -> None:
        self.anchor = _load_t012_anchor(self.process_id)
        base_inputs = zpole_default_sm_inputs()
        self.sm_inputs = zpole_inputs_with_charm_radiator(
            self.anchor.r_c.value,
            base_inputs,
        )
        self.sm_observables = zpole_evaluate_quark("c", inputs=self.sm_inputs)

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
                    "non-vetoing only; no Zcc Phase-3b NP prediction was evaluated"
                ),
                "required_parameter_point_extras": [_REQUIRED_EXTRA],
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
            delta_g_left_c = _coupling_entry(
                rs_ew_couplings, "z_delta_g_L_u", 1, 1
            )
            delta_g_right_c = _coupling_entry(
                rs_ew_couplings, "z_delta_g_R_u", 1, 1
            )
            shifted_charm = zpole_shifted_couplings(
                zpole_sm_couplings("c", self.sm_inputs),
                delta_g_left=delta_g_left_c,
                delta_g_right=delta_g_right_c,
            )
            prediction = zpole_evaluate_quark(
                "c", {"c": shifted_charm}, inputs=self.sm_inputs
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
                "r_c_formula": float(self.sm_observables.r_q),
                "r_c_yaml_anchor": float(self.anchor.r_c.value),
                "a_c_formula": float(self.sm_observables.a_q),
                "a_c_yaml_anchor": float(self.anchor.a_c.value),
                "a_fb_formula": float(self.sm_observables.a_fb),
                "a_fb_yaml_anchor": float(self.anchor.a_fb.value),
                "sin2_theta_eff": float(self.sm_inputs.sin2_theta_eff),
                "charm_radiator": float(self.sm_inputs.radiator_for("c")),
            },
            "a_fb_context": {
                "experimental": float(self.anchor.a_fb.value),
                "experimental_uncertainty": float(self.anchor.a_fb.uncertainty),
                "sm_formula": float(self.sm_observables.a_fb),
                "source_url": self.anchor.a_fb.source_url,
            },
            "rs_matching_assumption": getattr(
                rs_ew_couplings, "matching_assumption", None
            ),
            "rs_ew_model_label": getattr(rs_ew_couplings, "model_label", None),
            "rs_ew_kk_mass_gev": float(rs_ew_couplings.kk_ew_mass_gev),
            "delta_g_left_c": complex(delta_g_left_c),
            "delta_g_right_c": complex(delta_g_right_c),
            "shifted_zcc_couplings": {
                "g_left": complex(shifted_charm.g_left),
                "g_right": complex(shifted_charm.g_right),
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
                "Max one-sigma pull over R_c^0 and A_c. SM-limit "
                "pseudo-observables use effective Zcc couplings with the "
                "charm width radiator calibrated to the YAML R_c^0 anchor. "
                "Point-specific RS shifts use Phase-3a minimal-RS gauge "
                "neutral-current z_delta_g_L/R_u[1,1]."
            ),
            diagnostics=diagnostics,
        )
