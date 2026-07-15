"""TEMPLATE — copy this file to create a new constraint.

HOW TO USE
==========
1. Copy this file to:
       flavor_catalog_constraints/primary/<family>/<ID>.py
   (or secondary/<family>/<ID>.py for SECONDARY-tier).
   <family> must be one of registry.VALID_FAMILIES.
   <ID> must match ^[A-Z]+[0-9]+$ and have a sidecar at
       flavor_catalog/processes/<family>/<ID>.yaml
   (or .../secondary/<family>/<ID>.yaml).

2. Set ``process_id``, ``severity``, ``observable`` below.
   Do NOT set ``level`` or ``family`` — the registry derives them from
   this file's path.

3. List the candidate anchor keys your sidecar's ``pdg_or_equivalent``
   block actually uses (see the yaml). The loader tries them in order
   and raises loudly if none match — no silent default.

4. Reach physics ONLY through flavor_catalog_constraints.physics_adapters.*
   Never import quarkConstraints / flavorConstraints / qcd directly.

5. Add a per-constraint test by copying
   tests/constraints/test_constraint_template.py.

THE COMMON CASE (load an anchor, call an adapter, return a result) is
designed to be ~30 lines. Everything below the dashed line is fill-in.
"""

from __future__ import annotations

from flavor_catalog_constraints.anchors import load_anchor
from flavor_catalog_constraints.base import (
    ConstraintResult,
    ParameterPoint,
    Severity,
)
from flavor_catalog_constraints.registry import register

# Import the adapter wrapper(s) you need, e.g.:
# from flavor_catalog_constraints.physics_adapters.deltaf2 import epsilon_k_from_couplings

# --------------------------------------------------------------------- #
# FILL IN BELOW
# --------------------------------------------------------------------- #

# The family this constraint belongs to — must equal the directory name.
# (Used only for the anchor lookup; the registry independently re-derives
#  level+family from the module path and would reject a mismatch.)
_FAMILY = "kaon"

# Ordered list of keys your sidecar's pdg_or_equivalent block may use.
# The loader returns the first match and raises AnchorError if none exist.
_ANCHOR_CANDIDATES = (
    "canonical_experimental_value",
    "pdg_fit_assuming_cpt",
    "canonical_limit",
    "primary_current_limit",
)

# Optional defense-in-depth checks for the exact YAML value this constraint
# intends to consume. Set to None only when the sidecar genuinely lacks the
# field.
_EXPECTED_ANCHOR_VALUE_ID = "PDG2026:XX000:canonical_value"
_EXPECTED_ANCHOR_BLOCK_KEY = "canonical_experimental_value"
_EXPECTED_ANCHOR_UNITS = "dimensionless"
_EXPECTED_ANCHOR_CONFIDENCE_LEVEL = "95% CL"

# The extra (if any) this constraint reads off the ParameterPoint. Must be
# declared in point_builder.KNOWN_EXTRA_KEYS.
_REQUIRED_EXTRA = "quark_mass_basis_couplings"


@register
class Constraint:
    """One-line description of the observable and its RS sensitivity."""

    process_id = "XX000"           # <-- set to the catalog ID
    severity = Severity.HARD       # HARD / SOFT / INFO
    observable = "describe me"     # human-readable label

    def __init__(self) -> None:
        # Load + type the anchor at construction so a bad sidecar fails
        # at import time (and is isolated by the registry) rather than
        # mid-sweep.
        self.anchor = load_anchor(
            self.process_id,
            family=_FAMILY,
            candidates=_ANCHOR_CANDIDATES,
            expected_value_id=_EXPECTED_ANCHOR_VALUE_ID,
            expected_block_key=_EXPECTED_ANCHOR_BLOCK_KEY,
            expected_units=_EXPECTED_ANCHOR_UNITS,
            expected_confidence_level=_EXPECTED_ANCHOR_CONFIDENCE_LEVEL,
        )

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        inputs = point.get_extra(_REQUIRED_EXTRA)
        if inputs is None:
            # Legitimately absent extra: no input -> cannot predict. Report
            # INFO-style (no veto via ratio) with a clear note.
            return ConstraintResult(
                process_id=self.process_id,
                severity=self.severity,
                passes=True,
                experimental=self.anchor.value,
                notes=f"extra {_REQUIRED_EXTRA!r} absent; constraint not evaluated.",
                diagnostics={"missing_extra": _REQUIRED_EXTRA},
            )

        # ---- call the physics adapter and build the result --------------
        # M-18: if the extra is present but malformed/NaN, either let the
        # exception propagate to the registry or return diagnostics with
        # {"invalid_extra": _REQUIRED_EXTRA}; HARD invalid extras fail closed
        # under the base.ConstraintResult contract.
        # result = epsilon_k_from_couplings(inputs)
        # predicted = float(result.epsilon_k_np)
        # budget = float(result.epsilon_k_np_budget)
        # ratio = float(result.ratio_to_budget)
        # return ConstraintResult(
        #     process_id=self.process_id,
        #     severity=self.severity,
        #     passes=bool(result.passes),
        #     predicted=predicted,
        #     experimental=self.anchor.value,
        #     ratio=ratio,
        #     budget=budget,
        #     notes="…formula + budget definition…",
        #     diagnostics={"im_m12_np": float(result.im_m12_np)},
        # )
        raise NotImplementedError("fill in the physics evaluation")
