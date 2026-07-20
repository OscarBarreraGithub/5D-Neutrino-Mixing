"""E001 - electron electric dipole moment ``|d_e|``.

Physics
-------
The Standard Model electron EDM is negligible at the precision of the current
bound, so E001 is a pure-new-physics HARD constraint:

    |d_e^NP| <= d_e^limit.

The observable-side conversion is the low-energy CP-odd dipole convention
implemented in ``quarkConstraints.edm`` and reached only through
``flavor_catalog_constraints.physics_adapters.edm``:

    d_e [e cm] = c_e^CP-odd [GeV^-1] * hbar*c [GeV cm].

NEEDS-HUMAN-PHYSICS
-------------------
The RS contribution is a one-loop CP-violating dipole that needs complex
lepton, KK-fermion, Higgs, and electroweak couplings.  Those inputs are not
available on ``ParameterPoint``.  If the caller supplies an explicit electron
EDM proxy coefficient through ``lepton_mass_basis_couplings``, the result is
evaluated but flagged ``NEEDS-HUMAN-PHYSICS`` in diagnostics.

Catalog sidecar
---------------
``flavor_catalog/processes/edm_neutrino/E001.yaml`` is the source of truth for
the PDG/Roussy electron-EDM upper limit.  Numeric limit values are loaded from
that sidecar, not hardcoded here.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Mapping

from flavor_catalog_constraints.anchors import Anchor, AnchorError, load_anchor, load_pdg_block
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.edm import (
    EDM_COEFFICIENT_CONVENTION_V1,
    EDM_RS_MATCHING_GAP_V1,
    electron_edm_from_lepton_input,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "edm_neutrino"
_REQUIRED_EXTRA = "lepton_mass_basis_couplings"
_OPTIONAL_EW_MASS_EXTRA = "kk_ew_mass_gev"
_CURRENT_LIMIT_BLOCK = "canonical_limit"
_LIMIT_ANCHOR_CANDIDATES = (_CURRENT_LIMIT_BLOCK,)
_EXPECTED_UNITS = "e cm"
_SM_ELECTRON_EDM_E_CM = 0.0
_UNEVALUATED_REASON = (
    "no electron EDM dipole prediction available "
    "(complex RS lepton dipole matching inputs not on ParameterPoint)"
)
_UNEVALUATED_NOTES = f"NOT EVALUATED - {_UNEVALUATED_REASON}"


@dataclass(frozen=True)
class ElectronEDMAnchor:
    """Typed E001 anchor: current upper limit plus limit metadata."""

    experimental: Anchor
    limit_operator: str | None
    confidence_level: str | None
    used_measurement: str | None
    table_value: float | None
    table_units: str | None
    value_summary: str | None

    @property
    def value(self) -> float:
        """Upper limit on ``|d_e|`` in e cm."""
        return self.experimental.value

    @property
    def budget(self) -> float:
        """HARD veto budget for the pure-NP electron EDM."""
        return self.experimental.value

    @property
    def source_url(self) -> str | None:
        """Primary source URL for the limit."""
        return self.experimental.source_url


def _optional_str(value: Any) -> str | None:
    return None if value is None else str(value)


def _optional_float(value: Any, *, process_id: str, field_name: str) -> float | None:
    if value is None:
        return None
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise AnchorError(
            f"{process_id}: YAML field {field_name!r}={value!r} is not numeric"
        ) from exc
    if not math.isfinite(number):
        raise AnchorError(f"{process_id}: YAML field {field_name!r}={value!r} is not finite")
    return number


def _load_e001_anchor(process_id: str) -> ElectronEDMAnchor:
    experimental = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=_LIMIT_ANCHOR_CANDIDATES,
    )
    if experimental.block_key != _CURRENT_LIMIT_BLOCK:
        raise AnchorError(
            f"{process_id}: load_anchor selected {experimental.block_key!r}, "
            f"expected {_CURRENT_LIMIT_BLOCK!r} for the current electron EDM limit"
        )
    if experimental.value <= 0.0:
        raise AnchorError(f"{process_id}: electron EDM limit must be positive")
    if experimental.units != _EXPECTED_UNITS:
        raise AnchorError(
            f"{process_id}: expected units {_EXPECTED_UNITS!r}, got {experimental.units!r}"
        )

    pdg = load_pdg_block(process_id, family=_FAMILY)
    sub = pdg.get(experimental.block_key)
    if not isinstance(sub, Mapping):
        raise AnchorError(
            f"{process_id}: selected anchor block {experimental.block_key!r} "
            "is not available as a mapping"
        )
    return ElectronEDMAnchor(
        experimental=experimental,
        limit_operator=_optional_str(sub.get("limit_operator")),
        confidence_level=_optional_str(sub.get("confidence_level")),
        used_measurement=_optional_str(sub.get("used_measurement")),
        table_value=_optional_float(
            sub.get("table_value"),
            process_id=process_id,
            field_name=f"{experimental.block_key}.table_value",
        ),
        table_units=_optional_str(sub.get("table_units")),
        value_summary=_optional_str(sub.get("value_summary")),
    )


@register
class Constraint:
    """Catalogued electron EDM bound."""

    process_id = "E001"
    severity = Severity.HARD
    observable = "|d_e|"

    def __init__(self) -> None:
        self.anchor = _load_e001_anchor(self.process_id)
        self.sm_value = _SM_ELECTRON_EDM_E_CM

    def _unevaluated_result(
        self,
        *,
        diagnostics: Mapping[str, Any],
    ) -> ConstraintResult:
        merged_diagnostics = {
            "evaluated": False,
            "unevaluated_reason": _UNEVALUATED_REASON,
            "passes_semantics": (
                "non-vetoing only; no electron EDM prediction was evaluated"
            ),
            "budget_source": self.anchor.source_url,
            "experimental_block": self.anchor.experimental.block_key,
            "sm_edm_policy": (
                "SM electron EDM is negligible at catalog precision; "
                "reported sm_prediction is 0.0 for a pure-NP bound."
            ),
            "needs_human_physics": EDM_RS_MATCHING_GAP_V1,
            **dict(diagnostics),
        }
        return ConstraintResult(
            process_id=self.process_id,
            severity=self.severity,
            passes=True,
            predicted=None,
            sm_prediction=float(self.sm_value),
            experimental=float(self.anchor.value),
            ratio=None,
            budget=float(self.anchor.budget),
            notes=_UNEVALUATED_NOTES,
            diagnostics=merged_diagnostics,
        )

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        lepton_input = point.get_extra(_REQUIRED_EXTRA)
        if lepton_input is None:
            return self._unevaluated_result(
                diagnostics={
                    "missing_extra": _REQUIRED_EXTRA,
                },
            )

        kk_ew_mass = point.get_extra(_OPTIONAL_EW_MASS_EXTRA)
        try:
            result = electron_edm_from_lepton_input(
                lepton_input,
                experimental_limit_e_cm=self.anchor.budget,
                sm_edm_e_cm=self.sm_value,
                m_kk_gev=None if kk_ew_mass is None else float(kk_ew_mass),
            )
        except (KeyError, TypeError, ValueError) as exc:
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
                "edm_e_cm_signed": float(result.edm_e_cm),
                "abs_edm_e_cm": float(result.abs_edm_e_cm),
                "cp_odd_dipole_coefficient_gev_inv": float(
                    result.cp_odd_dipole_coefficient_gev_inv
                ),
                "hbarc_gev_cm": float(result.hbarc_gev_cm),
                "coefficient_source": result.coefficient_source,
                "coefficient_convention": EDM_COEFFICIENT_CONVENTION_V1,
                "rs_matching_gap": EDM_RS_MATCHING_GAP_V1,
                "sm_prediction_policy": (
                    "SM d_e is negligible compared with the current bound; "
                    "the limit is applied to the pure NP EDM."
                ),
                "budget_source": self.anchor.source_url,
                "experimental_block": self.anchor.experimental.block_key,
                "limit_operator": self.anchor.limit_operator,
                "confidence_level": self.anchor.confidence_level,
                "used_measurement": self.anchor.used_measurement,
                "value_summary": self.anchor.value_summary,
                "kk_ew_mass_extra_used": kk_ew_mass is not None,
            }
        )

        return ConstraintResult(
            process_id=self.process_id,
            severity=self.severity,
            passes=bool(result.passes),
            predicted=float(result.abs_edm_e_cm),
            sm_prediction=float(result.sm_edm_e_cm),
            experimental=float(self.anchor.value),
            ratio=float(result.ratio_to_limit),
            budget=float(result.experimental_limit_e_cm),
            notes=(
                "Pure-NP |d_e| bound using d_e[e cm] = c_CPodd[GeV^-1] "
                "* hbar*c. HARD budget is the E001.yaml canonical limit. "
                "The RS contribution is an explicit proxy coefficient and is "
                "flagged NEEDS-HUMAN-PHYSICS."
            ),
            diagnostics=diagnostics,
        )
