"""C001 - neutral charm mixing.

Physics
-------
``D0-D0bar`` mixing is long-distance dominated in the Standard Model, so a
short-distance SM subtraction is not defensible for this catalog veto.  The
new-physics contribution is evaluated as ``|M12^NP|`` with the audited
Delta F = 2 running path in ``quarkConstraints.deltaf2``:

    |M12^NP| <= Delta m_D^exp / 2.

The Wilson coefficients are QCD-evolved to ``mu_had = 3 GeV`` before applying
the D0 hadronic matrix elements.  This module reaches the core only through
``flavor_catalog_constraints.physics_adapters.deltaf2``.

Severity
--------
HARD.  Neutral-D mixing is observed, and a new short-distance contribution
larger than the measured splitting would overfill the conservative NP room.
The budget is intentionally the full measured splitting, not a central
SM-minus-experiment residual, because the charm SM prediction is theoretically
dirty.  Reviewers should scrutinize this physics policy if a future charm-SM
budget audit supersedes the current core convention.

Catalog sidecar
---------------
``flavor_catalog/processes/charm/C001.yaml`` is the source of truth for
``x_D``, ``y_D``, and ``Delta m_D`` provenance.  The HARD veto budget uses the
``Delta_m_D.value_GeV`` anchor loaded through the scaffold anchor loader:
``6.562e-15 GeV / 2 = 3.281e-15 GeV``.  This is an explicit catalog override of
the older hardcoded D0 value still present in the Delta F = 2 core
(``6.25e-15 GeV / 2 = 3.125e-15 GeV``).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

from flavor_catalog_constraints.anchors import Anchor, AnchorError, load_anchor
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.deltaf2 import (
    d0_mixing_from_wilsons_with_running,
    d0_mixing_wilsons_from_couplings,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "charm"
_REQUIRED_EXTRA = "quark_mass_basis_couplings"

_X_D_ANCHOR_CANDIDATES = ("x_D",)
_Y_D_ANCHOR_CANDIDATES = ("y_D",)
_DELTA_M_D_ANCHOR_CANDIDATES = ("Delta_m_D",)
_BUDGET_DOC_CITATION = (
    "quarkConstraints/deltaf2.py:959-961; "
    "docs/audits/bag_param_inventory.md:53; "
    "flavor_catalog/processes/charm/C001.yaml:145-157"
)
_MU_HAD_GEV = 3.0


@dataclass(frozen=True)
class D0MixingAnchor:
    """Typed C001 anchor: x_D, y_D, and Delta m_D in native/GeV units."""

    x_d: Anchor
    y_d: Anchor
    delta_m_d_native: Anchor
    delta_m_d_gev: Anchor

    @property
    def value(self) -> float:
        """Experimental ``Delta m_D`` central value in GeV."""
        return self.delta_m_d_gev.value

    @property
    def uncertainty(self) -> float | None:
        """Experimental ``Delta m_D`` uncertainty in GeV."""
        return self.delta_m_d_gev.uncertainty

    @property
    def source_url(self) -> str | None:
        """Experimental ``Delta m_D`` source URL."""
        return self.delta_m_d_gev.source_url

    @property
    def budget(self) -> float:
        """Conservative NP room ``Delta m_D^exp / 2`` in GeV."""
        return 0.5 * self.delta_m_d_gev.value


def _load_d0_mixing_anchor(process_id: str) -> D0MixingAnchor:
    anchor = D0MixingAnchor(
        x_d=load_anchor(
            process_id,
            family=_FAMILY,
            candidates=_X_D_ANCHOR_CANDIDATES,
        ),
        y_d=load_anchor(
            process_id,
            family=_FAMILY,
            candidates=_Y_D_ANCHOR_CANDIDATES,
        ),
        delta_m_d_native=load_anchor(
            process_id,
            family=_FAMILY,
            candidates=_DELTA_M_D_ANCHOR_CANDIDATES,
        ),
        delta_m_d_gev=load_anchor(
            process_id,
            family=_FAMILY,
            candidates=_DELTA_M_D_ANCHOR_CANDIDATES,
            value_key="value_GeV",
            uncertainty_key="uncertainty_GeV",
        ),
    )
    if anchor.value <= 0.0 or anchor.budget <= 0.0:
        raise AnchorError(
            f"{process_id}: Delta m_D anchor and NP budget must be positive "
            f"(delta_m_d={anchor.value}, budget={anchor.budget})"
        )
    return anchor


def _complex_mapping(values: Mapping[str, complex]) -> dict[str, complex]:
    return {str(key): complex(value) for key, value in values.items()}


@register
class Constraint:
    """Catalogued D0-D0bar Delta F=2 mixing constraint (process_id C001)."""

    process_id = "C001"
    severity = Severity.HARD
    observable = "D0-D0bar mixing"

    def __init__(self) -> None:
        self.anchor = _load_d0_mixing_anchor(self.process_id)

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        couplings = point.get_extra(_REQUIRED_EXTRA)
        if couplings is None:
            return ConstraintResult(
                process_id=self.process_id,
                severity=self.severity,
                passes=True,
                experimental=float(self.anchor.value),
                budget=float(self.anchor.budget),
                notes=(
                    f"extra {_REQUIRED_EXTRA!r} absent; D0 mixing constraint "
                    "was not evaluated."
                ),
                diagnostics={"missing_extra": _REQUIRED_EXTRA},
            )

        wilsons = d0_mixing_wilsons_from_couplings(couplings)
        # Budget policy: use the current catalog YAML PDG/HFLAV Delta_m_D
        # anchor, not the older hardcoded D0 default in quarkConstraints.
        result = d0_mixing_from_wilsons_with_running(
            wilsons,
            mu_had=_MU_HAD_GEV,
            m12_np_budget=self.anchor.budget,
        )

        predicted = float(result.abs_m12_np)
        ratio = float(result.ratio_to_budget)
        budget = float(result.budget)

        return ConstraintResult(
            process_id=self.process_id,
            severity=self.severity,
            passes=bool(result.passes),
            predicted=predicted,
            experimental=float(self.anchor.value),
            ratio=ratio,
            budget=budget,
            notes=(
                "|M12^NP| for D0 mixing; Wilsons are QCD-evolved to 3 GeV; "
                "HARD budget uses the conservative long-distance-dominated "
                "charm convention |M12^NP| <= Delta m_D^exp/2 from "
                f"{_BUDGET_DOC_CITATION}."
            ),
            diagnostics={
                "abs_m12_np_gev": predicted,
                "qcd_running_applied": True,
                "hadronic_scale_gev": _MU_HAD_GEV,
                "matching_scale_gev": float(wilsons.matching_scale),
                "m_kk_gev": float(wilsons.M_KK),
                "left_uc_coupling": complex(wilsons.left_coupling),
                "right_uc_coupling": complex(wilsons.right_coupling),
                "wilson_coefficients": _complex_mapping(wilsons.wilsons),
                "budget_doc_citation": _BUDGET_DOC_CITATION,
                "budget_policy": (
                    "Full measured D0 mass splitting is treated as NP room; "
                    "no short-distance SM subtraction is applied."
                ),
                "physics_uncertainty": (
                    "Charm mixing is long-distance dominated and theoretically "
                    "dirty; this is a conservative observed-amplitude veto."
                ),
                "delta_m_d_budget_gev": budget,
                "delta_m_d_experimental_gev": float(self.anchor.value),
                "delta_m_d_uncertainty_gev": (
                    None
                    if self.anchor.uncertainty is None
                    else float(self.anchor.uncertainty)
                ),
                "delta_m_d_native_value": float(self.anchor.delta_m_d_native.value),
                "delta_m_d_native_uncertainty": (
                    None
                    if self.anchor.delta_m_d_native.uncertainty is None
                    else float(self.anchor.delta_m_d_native.uncertainty)
                ),
                "delta_m_d_native_units": self.anchor.delta_m_d_native.units,
                "x_d_percent": float(self.anchor.x_d.value),
                "x_d_uncertainty_percent": (
                    None
                    if self.anchor.x_d.uncertainty is None
                    else float(self.anchor.x_d.uncertainty)
                ),
                "y_d_percent": float(self.anchor.y_d.value),
                "y_d_uncertainty_percent": (
                    None
                    if self.anchor.y_d.uncertainty is None
                    else float(self.anchor.y_d.uncertainty)
                ),
                "experimental_block": self.anchor.delta_m_d_gev.block_key,
                "x_d_block": self.anchor.x_d.block_key,
                "y_d_block": self.anchor.y_d.block_key,
            },
        )
