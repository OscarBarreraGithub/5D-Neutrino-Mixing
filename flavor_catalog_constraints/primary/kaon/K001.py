"""K001 — Indirect CP violation in K0-K0bar mixing (epsilon_K).

THROWAWAY EXAMPLE CONSTRAINT
============================
This file exists ONLY to prove the scaffold works end to end (anchor
loading + adapter call + result shape + registry discovery + the
example test). The orchestrator may delete it before the real,
physics-reviewed constraints land. It wraps existing, audited physics
(``quarkConstraints.deltaf2.evaluate_epsilon_k``) via the adapter, so it
is a faithful demonstration rather than a fake.

Severity: HARD. epsilon_K is an observed bound (~0.5% experimental
error); NP must fit inside |exp - SM|. A failure vetoes the point.

Sidecar: flavor_catalog/processes/kaon/K001.yaml
(uses the ``canonical_experimental_value`` block — the schema-flex
loader resolves it from the candidate list below).
"""

from __future__ import annotations

from flavor_catalog_constraints.anchors import load_anchor, load_pdg_block
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.deltaf2 import epsilon_k_from_couplings
from flavor_catalog_constraints.registry import register

_FAMILY = "kaon"
_REQUIRED_EXTRA = "quark_mass_basis_couplings"


@register
class Constraint:
    """Catalogued epsilon_K Delta F=2 constraint (process_id K001)."""

    process_id = "K001"
    severity = Severity.HARD
    observable = "epsilon_K"

    def __init__(self) -> None:
        # Experimental anchor (typed, loud on miss).
        self.anchor = load_anchor(
            self.process_id,
            family=_FAMILY,
            candidates=("canonical_experimental_value",),
        )
        # SM reference is in a sibling block; pull it explicitly and
        # type it so the NP budget |exp - SM| is reproducible.
        block = load_pdg_block(self.process_id, family=_FAMILY)
        self.sm_value = float(block["standard_model_reference"]["value"])

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        couplings = point.get_extra(_REQUIRED_EXTRA)
        if couplings is None:
            return ConstraintResult(
                process_id=self.process_id,
                severity=self.severity,
                passes=True,
                sm_prediction=self.sm_value,
                experimental=self.anchor.value,
                notes=f"extra {_REQUIRED_EXTRA!r} absent; epsilon_K not evaluated.",
                diagnostics={"missing_extra": _REQUIRED_EXTRA},
            )

        result = epsilon_k_from_couplings(couplings)
        return ConstraintResult(
            process_id=self.process_id,
            severity=self.severity,
            passes=bool(result.passes),
            predicted=float(result.epsilon_k_np),
            sm_prediction=self.sm_value,
            experimental=self.anchor.value,
            ratio=float(result.ratio_to_budget),
            budget=float(result.epsilon_k_np_budget),
            notes="|epsilon_K^NP| via kappa_epsilon convention; budget = |exp - SM|.",
            diagnostics={"im_m12_np": float(result.im_m12_np)},
        )
