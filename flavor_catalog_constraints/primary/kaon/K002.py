"""K002 — K_L-K_S mass difference (``Delta m_K``).

Severity rationale
------------------
HARD. The short-distance SM prediction is not used as a subtraction
because ``Delta m_K`` is long-distance dominated. Following the
existing Delta F = 2 implementation, the full experimental mass
difference is treated as the allowed NP room, i.e.
``|M_12^NP| <= Delta m_K^exp / 2``.

Catalog sidecar
---------------
``flavor_catalog/processes/kaon/K002.yaml`` (schema-flex anchor:
``pdg_or_equivalent.pdg_fit_assuming_cpt``).

Physics core
------------
``quarkConstraints/deltaf2.py``:

- :func:`compute_delta_f2_wilsons` matches tree-level KK-gluon Wilson
  coefficients onto the four-quark Delta F=2 operator basis.
- :func:`evaluate_delta_mk` plugs the kaon Wilsons into the
  ``|M_12^NP| / (Delta m_K / 2)`` convention and returns
  :class:`DeltaMKResult`.

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
    evaluate_delta_mk_from_couplings,
)
from flavor_catalog_constraints.registry import register_constraint

_HBAR_GEV_S = 6.582119569e-25


@dataclass(frozen=True)
class _Anchor:
    """Typed view over ``K002.yaml``'s ``pdg_fit_assuming_cpt`` block."""

    year: int
    source: str
    observable: str
    value: float
    uncertainty: float
    units: str
    display_value: str
    assumptions: tuple[str, ...]
    source_url: str
    access_date: str
    snapshot_path: str
    sha256: str
    sha256_of_text_snapshot: str

    @property
    def delta_m_k_exp_gev(self) -> float:
        """Convert ``10^10 hbar s^-1`` catalog units to GeV."""
        if self.units != "10^10 hbar s^-1":
            raise ValueError(f"Unsupported Delta m_K units: {self.units!r}")
        return self.value * 1.0e10 * _HBAR_GEV_S

    @property
    def m12_np_budget_gev(self) -> float:
        """The implemented bound is ``|M_12^NP| <= Delta m_K^exp / 2``."""
        return self.delta_m_k_exp_gev / 2.0


def _build_anchor(raw) -> _Anchor:
    """K002 uses ``pdg_fit_assuming_cpt``, not K001's canonical schema."""
    fit = raw["pdg_fit_assuming_cpt"]
    return _Anchor(
        year=int(fit["year"]),
        source=str(fit["source"]),
        observable=str(fit["observable"]),
        value=float(fit["value"]),
        uncertainty=float(fit["uncertainty"]),
        units=str(fit["units"]),
        display_value=str(fit["display_value"]),
        assumptions=tuple(str(item) for item in fit["assumptions"]),
        source_url=str(fit["source_url"]),
        access_date=str(fit["access_date"]),
        snapshot_path=str(fit["snapshot_path"]),
        sha256=str(fit["sha256"]),
        sha256_of_text_snapshot=str(fit["sha256_of_text_snapshot"]),
    )


@register_constraint
class Constraint:
    """Catalogued Delta m_K constraint (process_id ``K002``)."""

    process_id = "K002"
    severity = Severity.HARD
    observable = "Delta_m_K"

    def __init__(self) -> None:
        raw = load_raw(self.process_id, family="kaon", tier="primary")
        self.anchor = _build_anchor(raw)
        self.references = (self.anchor.snapshot_path,)

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        couplings = point.extras.get("quark_mass_basis_couplings")
        if couplings is None:
            return ConstraintResult(
                process_id=self.process_id,
                passes=False,
                predicted=None,
                sm_prediction=None,
                experimental=self.anchor.m12_np_budget_gev,
                ratio=None,
                budget=self.anchor.m12_np_budget_gev,
                severity=self.severity,
                notes=(
                    "ParameterPoint.extras['quark_mass_basis_couplings'] is "
                    "missing. point_builder.build_parameter_point() must "
                    "populate it before K002 can be evaluated."
                ),
                diagnostics={"missing_extra": "quark_mass_basis_couplings"},
            )

        # Adapter returns the existing DeltaMKResult dataclass unchanged.
        result = evaluate_delta_mk_from_couplings(couplings)
        return ConstraintResult(
            process_id=self.process_id,
            passes=bool(result.passes),
            predicted=float(result.abs_m12_np),
            sm_prediction=None,
            experimental=self.anchor.m12_np_budget_gev,
            ratio=float(result.ratio_to_exp),
            budget=self.anchor.m12_np_budget_gev,
            severity=self.severity,
            notes=(
                "|M_12^NP| from kaon Delta F=2 Wilsons; budget = "
                "Delta m_K^exp / 2 because the SM long-distance "
                "prediction is not subtracted."
            ),
            diagnostics={
                "delta_m_k_exp_gev": self.anchor.delta_m_k_exp_gev,
                "ratio_to_exp_from_core": float(result.ratio_to_exp),
            },
        )
