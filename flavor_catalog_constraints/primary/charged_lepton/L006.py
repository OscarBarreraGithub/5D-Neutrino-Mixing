"""L006 - muonium-antimuonium conversion.

Physics
-------
Muonium-antimuonium conversion is a Delta L=2 charged-lepton process,
``mu+ e- <-> mu- e+``.  The SM rate is negligible for catalog purposes, so
L006 is a pure-new-physics bound on the conversion probability,

    P_proxy = P_limit * |(G_C/G_F)_proxy / (G_C/G_F)_limit|^2,

when a caller supplies an explicit PDG-convention ``G_C/G_F`` proxy.  A direct
probability proxy is also accepted.  The low-energy probability mapper lives
in ``flavor_catalog_constraints.physics_adapters.muonium_conversion``; this
constraint imports no physics module directly.

NEEDS-HUMAN-PHYSICS
-------------------
The current ``ParameterPoint`` does not carry RS Delta L=2 four-lepton Wilson
coefficients, EW KK/Z/Z' lepton-current matching, spin-state normalization, or
magnetic-field response.  The v1 prediction is therefore a documented proxy
and is flagged in diagnostics.

Catalog sidecar
---------------
``flavor_catalog/processes/charged_lepton/L006.yaml`` is the source of truth
for the MACS/PSI probability limit and the PDG ``G_C/G_F`` limit.  Numeric
anchors are loaded from that sidecar, not hardcoded here.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Mapping

from flavor_catalog_constraints.anchors import (
    Anchor,
    AnchorError,
    load_anchor,
    load_pdg_block,
)
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.muonium_conversion import (
    MUONIUM_CONVERSION_OPERATOR_CONVENTION,
    MUONIUM_CONVERSION_PROXY_V1,
    muonium_conversion_from_lepton_input,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "charged_lepton"
_REQUIRED_EXTRA = "lepton_mass_basis_couplings"
_PROBABILITY_LIMIT_CANDIDATES = ("primary_current_probability_limit",)
_COUPLING_LIMIT_CANDIDATES = ("pdg_effective_coupling_limit",)
_ORIGINAL_EXPERIMENT_CANDIDATES = ("original_experiment",)
_UNEVALUATED_REASON = (
    "no muonium-antimuonium Delta L=2 prediction available "
    "(four-lepton RS matching inputs not on ParameterPoint)"
)
_UNEVALUATED_NOTES = f"NOT EVALUATED -- {_UNEVALUATED_REASON}"
_BUDGET_POLICY = (
    "Observed MACS/PSI 90% C.L. probability upper limit from L006.yaml; HARD "
    "veto uses the pure-NP M -> Mbar conversion probability proxy."
)


@dataclass(frozen=True)
class LimitMetadata:
    """Small provenance view over one L006 YAML limit block."""

    block_key: str
    relation: str | None
    confidence_level: str | None
    conditions: str | None


@dataclass(frozen=True)
class L006Anchor:
    """Typed L006 inputs: probability limit, coupling limit, MACS source."""

    probability_limit: Anchor
    effective_coupling_limit: Anchor
    original_experiment: Anchor
    probability_metadata: LimitMetadata
    coupling_metadata: LimitMetadata
    original_metadata: LimitMetadata

    @property
    def value(self) -> float:
        """Current upper limit on the conversion probability."""

        return self.probability_limit.value

    @property
    def budget(self) -> float:
        """HARD veto budget for the pure-NP conversion probability."""

        return self.probability_limit.value

    @property
    def coupling_limit_over_gf(self) -> float:
        """PDG effective-coupling limit in the stated convention."""

        return self.effective_coupling_limit.value

    @property
    def source_url(self) -> str | None:
        """Primary PDG probability-limit source URL."""

        return self.probability_limit.source_url


def _optional_str(value: Any) -> str | None:
    return None if value is None else str(value)


def _required_mapping(value: Any, *, process_id: str, path: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise AnchorError(f"{process_id}: YAML block {path!r} is not a mapping")
    return value


def _metadata_for_block(
    root: Mapping[str, Any],
    block_key: str,
    *,
    process_id: str,
) -> LimitMetadata:
    sub = _required_mapping(
        root.get(block_key),
        process_id=process_id,
        path=block_key,
    )
    return LimitMetadata(
        block_key=block_key,
        relation=_optional_str(sub.get("relation")),
        confidence_level=_optional_str(sub.get("confidence_level")),
        conditions=_optional_str(sub.get("conditions")),
    )


def _assert_positive(anchor: Anchor, *, process_id: str, label: str) -> None:
    if not math.isfinite(anchor.value) or anchor.value <= 0.0:
        raise AnchorError(f"{process_id}: {label} must be positive and finite")


def _load_l006_anchor(process_id: str) -> L006Anchor:
    probability = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=_PROBABILITY_LIMIT_CANDIDATES,
    )
    coupling = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=_COUPLING_LIMIT_CANDIDATES,
    )
    original = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=_ORIGINAL_EXPERIMENT_CANDIDATES,
    )
    _assert_positive(
        probability,
        process_id=process_id,
        label="primary_current_probability_limit",
    )
    _assert_positive(
        coupling,
        process_id=process_id,
        label="pdg_effective_coupling_limit",
    )
    _assert_positive(
        original,
        process_id=process_id,
        label="original_experiment",
    )
    if not math.isclose(probability.value, original.value, rel_tol=1.0e-12):
        raise AnchorError(
            f"{process_id}: PDG probability limit and original MACS value differ "
            f"({probability.value} vs {original.value})"
        )

    pdg = load_pdg_block(process_id, family=_FAMILY)
    return L006Anchor(
        probability_limit=probability,
        effective_coupling_limit=coupling,
        original_experiment=original,
        probability_metadata=_metadata_for_block(
            pdg,
            probability.block_key,
            process_id=process_id,
        ),
        coupling_metadata=_metadata_for_block(
            pdg,
            coupling.block_key,
            process_id=process_id,
        ),
        original_metadata=_metadata_for_block(
            pdg,
            original.block_key,
            process_id=process_id,
        ),
    )


@register
class Constraint:
    """Catalogued muonium-antimuonium conversion probability constraint."""

    process_id = "L006"
    severity = Severity.HARD
    observable = "P(M -> Mbar)"

    def __init__(self) -> None:
        self.anchor = _load_l006_anchor(self.process_id)

    def _unevaluated_result(
        self,
        *,
        diagnostics: Mapping[str, Any],
    ) -> ConstraintResult:
        merged = {
            "evaluated": False,
            "unevaluated_reason": _UNEVALUATED_REASON,
            "passes_semantics": (
                "non-vetoing only; no M -> Mbar probability prediction was "
                "evaluated"
            ),
            "sm_prediction_policy": (
                "SM muonium-antimuonium conversion is negligible for catalog "
                "purposes; L006 is a pure-NP upper-limit constraint."
            ),
            "budget_source": self.anchor.source_url,
            "budget_policy": _BUDGET_POLICY,
            "budget_limit_status": "observed_experimental_bound",
            "budget_verdict_role": "HARD observed upper-limit veto when evaluated",
            "probability_limit_block": self.anchor.probability_limit.block_key,
            "effective_coupling_limit_block": (
                self.anchor.effective_coupling_limit.block_key
            ),
            "effective_coupling_limit_over_gf": float(
                self.anchor.coupling_limit_over_gf
            ),
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

        try:
            result = muonium_conversion_from_lepton_input(
                lepton_input,
                probability_limit=self.anchor.budget,
                coupling_limit_over_gf=self.anchor.coupling_limit_over_gf,
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
                "needs_human_physics": MUONIUM_CONVERSION_PROXY_V1,
                "operator_convention": MUONIUM_CONVERSION_OPERATOR_CONVENTION,
                "sm_prediction_policy": (
                    "SM muonium-antimuonium conversion is negligible for "
                    "catalog purposes; budget is applied to the pure NP "
                    "conversion probability."
                ),
                "budget_source": self.anchor.source_url,
                "budget_policy": _BUDGET_POLICY,
                "budget_limit_status": "observed_experimental_bound",
                "budget_verdict_role": "HARD observed upper-limit veto",
                "probability_limit_block": self.anchor.probability_limit.block_key,
                "probability_limit_relation": (
                    self.anchor.probability_metadata.relation
                ),
                "probability_limit_confidence_level": (
                    self.anchor.probability_metadata.confidence_level
                ),
                "probability_limit_conditions": (
                    self.anchor.probability_metadata.conditions
                ),
                "effective_coupling_limit_block": (
                    self.anchor.effective_coupling_limit.block_key
                ),
                "effective_coupling_limit_relation": (
                    self.anchor.coupling_metadata.relation
                ),
                "effective_coupling_limit_confidence_level": (
                    self.anchor.coupling_metadata.confidence_level
                ),
                "effective_coupling_limit_conditions": (
                    self.anchor.coupling_metadata.conditions
                ),
                "original_experiment_block": self.anchor.original_experiment.block_key,
                "original_experiment_source": self.anchor.original_experiment.source,
                "original_experiment_source_url": (
                    self.anchor.original_experiment.source_url
                ),
            }
        )

        return ConstraintResult(
            process_id=self.process_id,
            severity=self.severity,
            passes=bool(result.passes),
            predicted=float(result.conversion_probability),
            sm_prediction=0.0,
            experimental=float(self.anchor.value),
            ratio=float(result.ratio_to_limit),
            budget=float(result.probability_limit),
            notes=(
                "Pure-NP muonium-antimuonium conversion probability proxy. "
                "A caller-supplied low-energy G_C/G_F value is calibrated "
                "with the L006 PDG/MACS probability and PDG effective-coupling "
                "anchors; full RS Delta L=2 four-lepton matching is flagged "
                "NEEDS-HUMAN-PHYSICS. The observed MACS probability limit is "
                "a HARD veto."
            ),
            diagnostics=diagnostics,
        )
