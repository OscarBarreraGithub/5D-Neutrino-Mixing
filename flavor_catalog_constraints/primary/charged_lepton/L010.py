"""L010 - charged-LFV three-body tau decay ``tau -> 3e``.

Physics
-------
The Standard Model rate is negligible for catalog purposes, so L010 is a
pure-new-physics branching-fraction upper bound.  The prediction reuses the
audited L002 low-energy machinery in ``quarkConstraints.lfv_three_body`` via
the tau-electron adapter:

* off-shell ``tau -> e gamma`` dipole conversion when an explicit dipole
  parent branching fraction is supplied,
* Z-penguin and box vector-contact proxy amplitudes,
* the same documented constructive dipole-contact interference envelope used
  by L002 when the chiral phase split is not supplied.

NEEDS-HUMAN-PHYSICS
-------------------
The current ``ParameterPoint`` does not carry the charged-lepton RS
neutral-current, EW KK/Z/Z', loop-level tau dipole, or box-matching inputs
needed for rigorous ``tau -> 3e``.  The contact and dipole inputs are therefore
explicit proxies and are flagged in diagnostics.

Catalog sidecar
---------------
``flavor_catalog/processes/charged_lepton/L010.yaml`` is the source of truth
for the Belle II branching-fraction limit.  Numeric values below are loaded
through the scaffold anchor loader, not hardcoded in this constraint.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Mapping

from flavor_catalog_constraints.anchors import (
    Anchor,
    AnchorError,
    load_anchor,
    load_pdg_block,
)
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.lfv_three_body_taue import (
    LFV_THREE_BODY_DIPOLE_CONTACT_INTERFERENCE_CONVENTION,
    LFV_THREE_BODY_OPERATOR_CONVENTION,
    LFV_THREE_BODY_PROXY_V1,
    TAU_TO_3E_PROXY_V1,
    tau_to_3e_default_sm_inputs,
    tau_to_3e_from_lepton_input,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "charged_lepton"
_REQUIRED_EXTRA = "lepton_mass_basis_couplings"
_OPTIONAL_EW_MASS_EXTRA = "kk_ew_mass_gev"
_LIMIT_ANCHOR_CANDIDATES = ("primary_current_limit",)
_PRIMARY_EXPERIMENT_CANDIDATES = ("primary_experiment",)
_SUPPORTING_EXPERIMENT_CANDIDATES = ("supporting_experiment",)
_UNEVALUATED_REASON = (
    "no tau -> 3e prediction available "
    "(lepton-sector RS contact/dipole inputs not on ParameterPoint)"
)
_UNEVALUATED_NOTES = f"NOT EVALUATED - {_UNEVALUATED_REASON}"


@dataclass(frozen=True)
class OptionalYAMLBlock:
    """Minimal typed provenance for a non-numeric YAML context block."""

    block_key: str
    source: str | None
    source_url: str | None
    statement: str | None
    snapshot_path: str | None


@dataclass(frozen=True)
class L010Anchor:
    """Typed L010 inputs: Belle II bound plus supporting context."""

    experimental: Anchor
    primary_experiment: Anchor
    supporting_experiment: Anchor
    prospects: OptionalYAMLBlock | None
    lhcb_context: OptionalYAMLBlock | None

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
    def belle_source_url(self) -> str | None:
        return self.primary_experiment.source_url


def _optional_str(value: Any) -> str | None:
    return None if value is None else str(value)


def _assert_positive_limit(anchor: Anchor, *, process_id: str, label: str) -> None:
    if not math.isfinite(anchor.value) or anchor.value <= 0.0:
        raise AnchorError(f"{process_id}: {label} branching-ratio limit must be positive")


def _optional_context_block(
    pdg: Mapping[str, Any],
    key: str,
    *,
    process_id: str,
) -> OptionalYAMLBlock | None:
    raw = pdg.get(key)
    if raw is None:
        return None
    if not isinstance(raw, Mapping):
        raise AnchorError(f"{process_id}: pdg_or_equivalent.{key} is not a mapping")
    return OptionalYAMLBlock(
        block_key=key,
        source=_optional_str(raw.get("source")),
        source_url=_optional_str(raw.get("source_url")),
        statement=_optional_str(raw.get("statement")),
        snapshot_path=_optional_str(raw.get("snapshot_path")),
    )


def _load_l010_anchor(process_id: str) -> L010Anchor:
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
    supporting_experiment = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=_SUPPORTING_EXPERIMENT_CANDIDATES,
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
        supporting_experiment,
        process_id=process_id,
        label="supporting_experiment",
    )
    if not math.isclose(experimental.value, primary_experiment.value, rel_tol=0.0, abs_tol=0.0):
        raise AnchorError(
            f"{process_id}: primary_current_limit and primary_experiment differ"
        )
    pdg = load_pdg_block(process_id, family=_FAMILY)
    return L010Anchor(
        experimental=experimental,
        primary_experiment=primary_experiment,
        supporting_experiment=supporting_experiment,
        prospects=_optional_context_block(pdg, "prospects", process_id=process_id),
        lhcb_context=_optional_context_block(pdg, "lhcb_context", process_id=process_id),
    )


@register
class Constraint:
    """Catalogued ``tau -> 3e`` charged-LFV three-body constraint."""

    process_id = "L010"
    severity = Severity.HARD
    observable = "BR(tau -> 3e)"

    def __init__(self) -> None:
        self.anchor = _load_l010_anchor(self.process_id)
        self.sm_inputs = tau_to_3e_default_sm_inputs()

    def _unevaluated_result(
        self,
        *,
        diagnostics: Mapping[str, object],
    ) -> ConstraintResult:
        merged = {
            "evaluated": False,
            "unevaluated_reason": _UNEVALUATED_REASON,
            "passes_semantics": (
                "non-vetoing only; no BR(tau -> 3e) prediction was evaluated"
            ),
            "sm_prediction_policy": (
                "Charged-LFV SM rate is negligible; L010 is applied as a "
                "pure-NP branching-fraction upper bound."
            ),
            "budget_source": self.anchor.source_url,
            "belle_source": self.anchor.belle_source_url,
            "experimental_block": self.anchor.experimental.block_key,
            "primary_experiment_block": self.anchor.primary_experiment.block_key,
            "supporting_experiment_block": self.anchor.supporting_experiment.block_key,
            "supporting_experiment_limit": float(self.anchor.supporting_experiment.value),
            "supporting_experiment_source": self.anchor.supporting_experiment.source_url,
            "prospects_source": (
                None if self.anchor.prospects is None else self.anchor.prospects.source_url
            ),
            "lhcb_context_source": (
                None
                if self.anchor.lhcb_context is None
                else self.anchor.lhcb_context.source_url
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

        kk_ew_mass = point.get_extra(_OPTIONAL_EW_MASS_EXTRA)
        try:
            result = tau_to_3e_from_lepton_input(
                lepton_input,
                br_limit=self.anchor.budget,
                m_kk_gev=None if kk_ew_mass is None else float(kk_ew_mass),
                inputs=self.sm_inputs,
            )
        except (AttributeError, KeyError, TypeError, ValueError) as exc:
            diagnostics: dict[str, object] = {
                "exception_type": type(exc).__name__,
                "exception": str(exc),
            }
            if _is_unevaluated_rejection(exc, lepton_input):
                diagnostics["caller_flavor_or_alias_rejected"] = True
            else:
                diagnostics["invalid_extra"] = _REQUIRED_EXTRA
            return self._unevaluated_result(diagnostics=diagnostics)

        diagnostics = dict(result.diagnostics)
        diagnostics.update(
            {
                "evaluated": True,
                "needs_human_physics": TAU_TO_3E_PROXY_V1,
                "shared_contact_proxy": LFV_THREE_BODY_PROXY_V1,
                "operator_convention": LFV_THREE_BODY_OPERATOR_CONVENTION,
                "dipole_contact_interference_convention": (
                    LFV_THREE_BODY_DIPOLE_CONTACT_INTERFERENCE_CONVENTION
                ),
                "sm_prediction_policy": (
                    "Charged-LFV SM rate is negligible; L010 is applied as a "
                    "pure-NP branching-fraction upper bound."
                ),
                "budget_source": self.anchor.source_url,
                "belle_source": self.anchor.belle_source_url,
                "experimental_block": self.anchor.experimental.block_key,
                "primary_experiment_block": self.anchor.primary_experiment.block_key,
                "supporting_experiment_block": self.anchor.supporting_experiment.block_key,
                "supporting_experiment_limit": float(
                    self.anchor.supporting_experiment.value
                ),
                "supporting_experiment_source": (
                    self.anchor.supporting_experiment.source_url
                ),
                "prospects_source": (
                    None
                    if self.anchor.prospects is None
                    else self.anchor.prospects.source_url
                ),
                "lhcb_context_source": (
                    None
                    if self.anchor.lhcb_context is None
                    else self.anchor.lhcb_context.source_url
                ),
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
                "Pure-NP BR(tau -> 3e) from the shared L002 "
                "lfv_three_body dipole/contact formula pinned to "
                "initial=tau and final=e. The HARD budget is the Belle II "
                "observed branching-fraction limit from L010.yaml; RS contact, box, "
                "and tau dipole matching are flagged NEEDS-HUMAN-PHYSICS."
            ),
            diagnostics=diagnostics,
        )


_MISMATCHED_SPURION_KEYS = (
    "left_emu_overlap",
    "right_emu_overlap",
    "left_emu",
    "right_emu",
    "left_mue_overlap",
    "right_mue_overlap",
    "left_mu_e_overlap",
    "right_mu_e_overlap",
    "left_e_mu_overlap",
    "right_e_mu_overlap",
    "left_taumu_overlap",
    "right_taumu_overlap",
    "left_mutau_overlap",
    "right_mutau_overlap",
    "left_tau_mu_overlap",
    "right_tau_mu_overlap",
    "left_mu_tau_overlap",
    "right_mu_tau_overlap",
)


def _is_unevaluated_rejection(exc: Exception, lepton_input: Any) -> bool:
    message = str(exc)
    if not isinstance(exc, ValueError):
        return False
    if (
        "tau_to_3e_from_lepton_input is pinned to initial_flavor='tau'" in message
        or "tau_to_3e_from_lepton_input is pinned to final_flavor='e'" in message
    ):
        return True
    if "tau->3e does not accept mu-e or tau-mu overlap aliases" not in message:
        return False
    return _has_finite_mismatched_spurion_alias(lepton_input)


def _has_finite_mismatched_spurion_alias(lepton_input: Any) -> bool:
    return any(
        _finite_complex_like(value)
        for value in _present_values(lepton_input, _MISMATCHED_SPURION_KEYS)
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
