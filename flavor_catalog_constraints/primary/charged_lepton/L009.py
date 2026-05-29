"""L009 - charged-LFV three-body tau decay ``tau -> 3mu``.

Physics
-------
The Standard Model rate is negligible for catalog purposes, so L009 is a
pure-new-physics branching-fraction upper bound.  The prediction reuses the
audited L002 low-energy machinery in ``quarkConstraints.lfv_three_body`` via
the tau-specific adapter:

* off-shell ``tau -> mu gamma`` dipole conversion when an explicit dipole
  parent branching fraction is supplied,
* Z-penguin and box vector-contact proxy amplitudes,
* the same documented constructive dipole-contact interference envelope used
  by L002 when the chiral phase split is not supplied.

NEEDS-HUMAN-PHYSICS
-------------------
The current ``ParameterPoint`` does not carry the charged-lepton RS
neutral-current, EW KK/Z/Z', loop-level tau dipole, or box-matching inputs
needed for rigorous ``tau -> 3mu``.  The contact and dipole inputs are therefore
explicit proxies and are flagged in diagnostics.

Catalog sidecar
---------------
``flavor_catalog/processes/charged_lepton/L009.yaml`` is the source of truth
for the Belle II/PDG branching-fraction limit.  Numeric values below are loaded
through the scaffold anchor loader, not hardcoded in this constraint.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

from flavor_catalog_constraints.anchors import Anchor, AnchorError, load_anchor
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.lfv_three_body_tau import (
    LFV_THREE_BODY_DIPOLE_CONTACT_INTERFERENCE_CONVENTION,
    LFV_THREE_BODY_OPERATOR_CONVENTION,
    LFV_THREE_BODY_PROXY_V1,
    TAU_TO_3MU_PROXY_V1,
    tau_to_3mu_default_sm_inputs,
    tau_to_3mu_from_lepton_input,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "charged_lepton"
_REQUIRED_EXTRA = "lepton_mass_basis_couplings"
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
        if lepton_input is None:
            return self._unevaluated_result(
                diagnostics={"missing_extra": _REQUIRED_EXTRA},
            )

        kk_ew_mass = point.get_extra(_OPTIONAL_EW_MASS_EXTRA)
        try:
            result = tau_to_3mu_from_lepton_input(
                lepton_input,
                br_limit=self.anchor.budget,
                m_kk_gev=None if kk_ew_mass is None else float(kk_ew_mass),
                inputs=self.sm_inputs,
            )
        except (AttributeError, KeyError, TypeError, ValueError) as exc:
            return self._unevaluated_result(
                diagnostics={
                    "invalid_extra": _REQUIRED_EXTRA,
                    "exception_type": type(exc).__name__,
                    "exception": str(exc),
                },
            )

        diagnostics = dict(result.diagnostics)
        diagnostics.update(
            {
                "evaluated": True,
                "needs_human_physics": TAU_TO_3MU_PROXY_V1,
                "shared_contact_proxy": LFV_THREE_BODY_PROXY_V1,
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
                "branching-fraction limit from L009.yaml; RS contact, box, "
                "and tau dipole matching are flagged NEEDS-HUMAN-PHYSICS."
            ),
            diagnostics=diagnostics,
        )
