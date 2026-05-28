"""K001 — Indirect CP violation in K0-K0bar mixing (``epsilon_K``).

Severity rationale
------------------
HARD. The experimental value is PDG-bounded to ~0.5% relative
uncertainty; the SM prediction (Brod, Gorbahn, Stamou 2020) has a
~8% theory uncertainty dominated by hadronic and parametric inputs.
The NP budget is taken as ``|epsilon_K^exp - epsilon_K^SM|`` and a
point fails if the predicted ``|epsilon_K^NP|`` exceeds it. This is
a veto-grade constraint in any RS flavor scan.

Catalog sidecar
---------------
``flavor_catalog/processes/kaon/K001.yaml`` (Pattern A of the plan —
the *only* yaml that uses ``canonical_experimental_value:``).

Physics core
------------
``quarkConstraints/deltaf2.py``:

- :func:`compute_delta_f2_wilsons` matches tree-level KK-gluon Wilson
  coefficients onto the four-quark Delta F=2 operator basis.
- :func:`evaluate_epsilon_k` plugs the Wilsons into the
  ``kappa_epsilon`` convention and returns
  :class:`EpsilonKResult`.

This module reaches the physics core via
:mod:`flavor_catalog_constraints.physics_adapters.deltaf2` so that
any future signature change in the physics module is absorbed at the
adapter boundary.
"""

from __future__ import annotations

from dataclasses import dataclass

from flavor_catalog_constraints.anchors import load_raw
from flavor_catalog_constraints.base import (
    ConstraintResult,
    ParameterPoint,
    Severity,
)
from flavor_catalog_constraints.physics_adapters.deltaf2 import (
    evaluate_epsilon_k_from_couplings,
)
from flavor_catalog_constraints.registry import register_constraint


@dataclass(frozen=True)
class _K001Anchor:
    """Typed view over the subset of ``K001.yaml`` this constraint reads."""

    eps_K_exp_value: float
    eps_K_exp_uncertainty: float
    eps_K_sm_value: float
    snapshot_path_exp: str
    sha256_exp: str
    snapshot_path_sm: str
    sha256_sm: str


def _build_anchor(raw) -> _K001Anchor:
    """Pattern A — K001 uses ``canonical_experimental_value:`` (only yaml that does)."""
    exp = raw["canonical_experimental_value"]
    sm = raw["standard_model_reference"]
    return _K001Anchor(
        eps_K_exp_value=float(exp["value"]),
        eps_K_exp_uncertainty=float(exp["uncertainty"]),
        eps_K_sm_value=float(sm["value"]),
        snapshot_path_exp=str(exp["snapshot_path"]),
        sha256_exp=str(exp["sha256_of_local_snapshot"]),
        snapshot_path_sm=str(sm["snapshot_path"]),
        sha256_sm=str(sm["sha256_of_local_snapshot"]),
    )


@register_constraint
class Constraint:
    """Catalogued epsilon_K constraint (process_id ``K001``)."""

    process_id = "K001"
    severity = Severity.HARD
    observable = "epsilon_K"

    def __init__(self) -> None:
        raw = load_raw(self.process_id, family="kaon", tier="primary")
        self.anchor = _build_anchor(raw)
        self.references = (
            self.anchor.snapshot_path_exp,
            self.anchor.snapshot_path_sm,
        )

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        couplings = point.extras.get("quark_mass_basis_couplings")
        if couplings is None:
            return ConstraintResult(
                process_id=self.process_id,
                passes=False,
                predicted=None,
                sm_prediction=self.anchor.eps_K_sm_value,
                experimental=self.anchor.eps_K_exp_value,
                ratio=None,
                budget=None,
                severity=self.severity,
                notes=(
                    "ParameterPoint.extras['quark_mass_basis_couplings'] is "
                    "missing. point_builder.build_parameter_point() must "
                    "populate it before K001 can be evaluated."
                ),
                diagnostics={"missing_extra": "quark_mass_basis_couplings"},
            )

        # Adapter returns the existing EpsilonKResult dataclass unchanged.
        result = evaluate_epsilon_k_from_couplings(couplings)
        return ConstraintResult(
            process_id=self.process_id,
            passes=bool(result.passes),
            predicted=float(result.epsilon_k_np),
            sm_prediction=self.anchor.eps_K_sm_value,
            experimental=self.anchor.eps_K_exp_value,
            ratio=float(result.ratio_to_budget),
            budget=float(result.epsilon_k_np_budget),
            severity=self.severity,
            notes=(
                "|epsilon_K^NP| from Im(M12^NP) via the kappa_epsilon "
                "convention; budget = |epsilon_K^exp - epsilon_K^SM|."
            ),
            diagnostics={"im_m12_np": float(result.im_m12_np)},
        )
