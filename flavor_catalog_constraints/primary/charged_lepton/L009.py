"""L009 - charged-LFV three-body tau decay ``tau -> 3mu``.

Physics
-------
The Standard Model rate is negligible for catalog purposes, so L009 is a
pure-new-physics branching-fraction upper bound.  The prediction reuses the
audited L002 low-energy machinery in ``quarkConstraints.lfv_three_body`` via
the tau-specific adapter:

* off-shell ``tau -> mu gamma`` dipole conversion when an explicit dipole
  parent branching fraction is supplied,
* tree-level light-Z vector-contact amplitudes from Phase-4a lepton Z matrices,
* deferred box vector-contact inputs,
* the same documented constructive dipole-contact interference envelope used
  by L002 when the chiral phase split is not supplied.

NEEDS-HUMAN-PHYSICS
-------------------
The tree-level light-Z contact is rigorous when ``rs_ew_couplings`` is present
and is zero for the Phase-4a diagonal charged-lepton fit.  Tau dipole matching,
dipole-contact phase, heavy neutral exchange, and box matching remain deferred.

Catalog sidecar
---------------
``flavor_catalog/processes/charged_lepton/L009.yaml`` is the source of truth
for the Belle II/PDG branching-fraction limit.  Numeric values below are loaded
through the scaffold anchor loader, not hardcoded in this constraint.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Mapping

import numpy as np

from flavor_catalog_constraints.anchors import Anchor, AnchorError, load_anchor
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.lfv_three_body_tau import (
    LFV_THREE_BODY_DEFERRED_PIECES_V1,
    LFV_THREE_BODY_DIPOLE_CONTACT_INTERFERENCE_CONVENTION,
    LFV_THREE_BODY_OPERATOR_CONVENTION,
    LFV_THREE_BODY_TREE_CONTACT_RIGOROUS_V1,
    TAU_TO_3MU_PROXY_V1,
    tau_to_3mu_default_sm_inputs,
    tau_to_3mu_from_lepton_input,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "charged_lepton"
_REQUIRED_EXTRA = "lepton_mass_basis_couplings"
_RS_EW_EXTRA = "rs_ew_couplings"
_OPTIONAL_EW_MASS_EXTRA = "kk_ew_mass_gev"
_LIMIT_ANCHOR_CANDIDATES = ("primary_current_limit",)
_PRIMARY_EXPERIMENT_CANDIDATES = ("primary_experiment",)
_POST_PDG_UPDATE_CANDIDATES = ("post_pdg_update",)
_UNEVALUATED_REASON = (
    "no tau -> 3mu prediction available "
    "(lepton-sector RS contact/dipole inputs not on ParameterPoint)"
)
_UNEVALUATED_NOTES = f"NOT EVALUATED - {_UNEVALUATED_REASON}"


@dataclass(frozen=True)
class L009Anchor:
    """Typed L009 inputs: PDG/Belle II bound plus post-PDG context."""

    experimental: Anchor
    primary_experiment: Anchor
    post_pdg_update_90cl: Anchor

    @property
    def value(self) -> float:
        return self.experimental.value

    @property
    def budget(self) -> float:
        return self.experimental.value

    @property
    def source_url(self) -> str | None:
        return self.experimental.source_url

    @property
    def belle_ii_source_url(self) -> str | None:
        return self.primary_experiment.source_url


def _assert_positive_limit(anchor: Anchor, *, process_id: str, label: str) -> None:
    if not math.isfinite(anchor.value) or anchor.value <= 0.0:
        raise AnchorError(f"{process_id}: {label} branching-ratio limit must be positive")


def _load_l009_anchor(process_id: str) -> L009Anchor:
    experimental = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=_LIMIT_ANCHOR_CANDIDATES,
    )
    primary_experiment = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=_PRIMARY_EXPERIMENT_CANDIDATES,
    )
    post_pdg_update = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=_POST_PDG_UPDATE_CANDIDATES,
        value_key="value_90cl",
    )
    _assert_positive_limit(
        experimental,
        process_id=process_id,
        label="primary_current_limit",
    )
    _assert_positive_limit(
        primary_experiment,
        process_id=process_id,
        label="primary_experiment",
    )
    _assert_positive_limit(
        post_pdg_update,
        process_id=process_id,
        label="post_pdg_update.value_90cl",
    )
    if not math.isclose(experimental.value, primary_experiment.value, rel_tol=0.0, abs_tol=0.0):
        raise AnchorError(
            f"{process_id}: primary_current_limit and primary_experiment differ"
        )
    return L009Anchor(
        experimental=experimental,
        primary_experiment=primary_experiment,
        post_pdg_update_90cl=post_pdg_update,
    )


@register
class Constraint:
    """Catalogued ``tau -> 3mu`` charged-LFV three-body constraint."""

    process_id = "L009"
    severity = Severity.HARD
    observable = "BR(tau -> 3mu)"

    def __init__(self) -> None:
        self.anchor = _load_l009_anchor(self.process_id)
        self.sm_inputs = tau_to_3mu_default_sm_inputs()

    def _unevaluated_result(
        self,
        *,
        diagnostics: dict[str, object],
    ) -> ConstraintResult:
        merged = {
            "evaluated": False,
            "unevaluated_reason": _UNEVALUATED_REASON,
            "passes_semantics": (
                "non-vetoing only; no BR(tau -> 3mu) prediction was evaluated"
            ),
            "sm_prediction_policy": (
                "Charged-LFV SM rate is negligible; L009 is applied as a "
                "pure-NP branching-fraction upper bound."
            ),
            "budget_source": self.anchor.source_url,
            "belle_ii_source": self.anchor.belle_ii_source_url,
            "experimental_block": self.anchor.experimental.block_key,
            "primary_experiment_block": self.anchor.primary_experiment.block_key,
            "post_pdg_update_90cl": float(self.anchor.post_pdg_update_90cl.value),
            "post_pdg_update_source": self.anchor.post_pdg_update_90cl.source_url,
            **dict(diagnostics),
        }
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
            diagnostics=merged,
        )

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        lepton_input = point.get_extra(_REQUIRED_EXTRA)
        rs_ew_couplings = point.get_extra(_RS_EW_EXTRA)
        if lepton_input is None:
            if rs_ew_couplings is None:
                return self._unevaluated_result(
                    diagnostics={"missing_extra": _REQUIRED_EXTRA},
                )
            lepton_input = _tree_only_lepton_input()

        kk_ew_mass = point.get_extra(_OPTIONAL_EW_MASS_EXTRA)
        try:
            result = tau_to_3mu_from_lepton_input(
                _adapter_input(lepton_input, rs_ew_couplings),
                br_limit=self.anchor.budget,
                m_kk_gev=None if kk_ew_mass is None else float(kk_ew_mass),
                inputs=self.sm_inputs,
            )
        except (AttributeError, KeyError, TypeError, ValueError) as exc:
            diagnostics: dict[str, object] = {
                "exception_type": type(exc).__name__,
                "exception": str(exc),
            }
            if _is_unevaluated_legacy_proxy_exception(exc, lepton_input):
                diagnostics["legacy_overlap_tree_proxy_ignored"] = True
            else:
                diagnostics["invalid_extra"] = _REQUIRED_EXTRA
            return self._unevaluated_result(diagnostics=diagnostics)

        diagnostics = dict(result.diagnostics)
        diagnostics.update(
            {
                "evaluated": True,
                "needs_human_physics": LFV_THREE_BODY_DEFERRED_PIECES_V1,
                "tree_contact_matching": diagnostics.get(
                    "matching_assumption",
                    LFV_THREE_BODY_TREE_CONTACT_RIGOROUS_V1,
                ),
                "shared_contact_policy": TAU_TO_3MU_PROXY_V1,
                "operator_convention": LFV_THREE_BODY_OPERATOR_CONVENTION,
                "dipole_contact_interference_convention": (
                    LFV_THREE_BODY_DIPOLE_CONTACT_INTERFERENCE_CONVENTION
                ),
                "sm_prediction_policy": (
                    "Charged-LFV SM rate is negligible; L009 is applied as a "
                    "pure-NP branching-fraction upper bound."
                ),
                "budget_source": self.anchor.source_url,
                "belle_ii_source": self.anchor.belle_ii_source_url,
                "experimental_block": self.anchor.experimental.block_key,
                "primary_experiment_block": self.anchor.primary_experiment.block_key,
                "post_pdg_update_90cl": float(self.anchor.post_pdg_update_90cl.value),
                "post_pdg_update_source": self.anchor.post_pdg_update_90cl.source_url,
                "initial_flavor": result.initial_flavor,
                "final_flavor": result.final_flavor,
                "kk_ew_mass_extra_used": kk_ew_mass is not None,
                "rs_ew_couplings_extra_present": rs_ew_couplings is not None,
                "reused_physics_module": "quarkConstraints.lfv_three_body",
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
                "Pure-NP BR(tau -> 3mu) from the shared L002 "
                "lfv_three_body dipole/contact formula pinned to "
                "initial=tau and final=mu. The HARD budget is the Belle II/PDG "
                "branching-fraction limit from L009.yaml; tree light-Z contact "
                "is rigorous when present, while box, heavy-neutral, and tau "
                "dipole matching remain NEEDS-HUMAN-PHYSICS."
            ),
            diagnostics=diagnostics,
        )


def _adapter_input(lepton_input: Any, rs_ew_couplings: Any | None) -> Any:
    if rs_ew_couplings is None:
        return lepton_input
    if isinstance(lepton_input, Mapping):
        return {
            **dict(lepton_input),
            "lepton_mass_basis_couplings": lepton_input,
            "rs_ew_couplings": rs_ew_couplings,
        }
    return {
        "lepton_mass_basis_couplings": lepton_input,
        "dipole": lepton_input,
        "rs_ew_couplings": rs_ew_couplings,
    }


def _is_unevaluated_legacy_proxy_exception(exc: Exception, lepton_input: Any) -> bool:
    """Return true for stale overlap/flavor probes that carry no L009 prediction."""

    message = str(exc)
    if isinstance(exc, ValueError):
        return (
            "pinned to initial_flavor='tau'" in message
            or "mu->e/e-mu overlap aliases" in message
        )
    if not isinstance(exc, TypeError):
        return False
    if "must provide tau->3mu contact proxy inputs" not in message:
        return False
    return _has_finite_legacy_overlap_like_input(lepton_input)


_LEGACY_SCALAR_OVERLAP_KEYS = (
    "left_lfv_overlap",
    "right_lfv_overlap",
    "left_emu_overlap",
    "right_emu_overlap",
    "left_emu",
    "right_emu",
)
_LEGACY_MATRIX_OVERLAP_KEYS = (
    "left_charged_lepton_overlap",
    "right_charged_lepton_overlap",
    "left_lepton_overlap",
    "right_lepton_overlap",
    "left_overlap",
    "right_overlap",
)


def _has_finite_legacy_overlap_like_input(lepton_input: Any) -> bool:
    return any(
        _finite_complex_like(value)
        for value in _present_values(lepton_input, _LEGACY_SCALAR_OVERLAP_KEYS)
    ) or any(
        _finite_matrix_like(value)
        for value in _present_values(lepton_input, _LEGACY_MATRIX_OVERLAP_KEYS)
    )


def _present_values(value: Any, keys: tuple[str, ...]):
    if isinstance(value, Mapping):
        for key in keys:
            if key in value:
                yield value[key]
        return
    for key in keys:
        if hasattr(value, key):
            yield getattr(value, key)


def _finite_complex_like(value: Any) -> bool:
    try:
        number = complex(value)
    except (TypeError, ValueError):
        return False
    return math.isfinite(number.real) and math.isfinite(number.imag)


def _finite_matrix_like(value: Any) -> bool:
    try:
        matrix = np.asarray(value, dtype=np.complex128)
    except (TypeError, ValueError):
        return False
    return (
        matrix.shape == (3, 3)
        and bool(np.all(np.isfinite(matrix.real)))
        and bool(np.all(np.isfinite(matrix.imag)))
    )


def _tree_only_lepton_input() -> dict[str, str]:
    return {
        "source": "tree-only rs_ew_couplings input; lepton_mass_basis_couplings absent"
    }
