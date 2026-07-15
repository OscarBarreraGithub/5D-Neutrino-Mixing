"""L002 - charged-LFV three-body muon decay ``mu -> 3e``.

Physics
-------
The Standard Model rate is negligible for catalog purposes, so L002 is a
pure-new-physics branching-fraction bound.  The v1 prediction combines

* the off-shell dipole contribution reused from the L001 lepton-dipole adapter,
* tree-level light-Z vector-contact amplitudes from Phase-4a lepton Z matrices,
* caller-supplied box vector-contact amplitudes in the same convention.

The lower-level formula and proxy matching live in
``quarkConstraints.lfv_three_body`` and are reached only through the
``flavor_catalog_constraints.physics_adapters.lfv_three_body`` boundary.

NEEDS-HUMAN-PHYSICS
-------------------
The tree-level light-Z contact is rigorous when ``rs_ew_couplings`` is present
and is zero for the Phase-4a diagonal charged-lepton fit.  Dipole matching,
dipole-contact phase, heavy neutral exchange, and box matching remain deferred.

Catalog sidecar
---------------
``flavor_catalog/processes/charged_lepton/L002.yaml`` is the source of truth
for the SINDRUM/PDG branching-fraction limit.  The L001 dipole prefactor and
normalization limit are loaded from ``L001.yaml`` so no experimental or
normalization numbers are hardcoded here.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Mapping

import numpy as np

from flavor_catalog_constraints.anchors import (
    Anchor,
    AnchorError,
    load_anchor,
    load_pdg_block,
)
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.lfv_three_body import (
    LFV_THREE_BODY_DEFERRED_PIECES_V1,
    LFV_THREE_BODY_OPERATOR_CONVENTION,
    LFV_THREE_BODY_TREE_CONTACT_RIGOROUS_V1,
    lfv_three_body_default_sm_inputs,
    mu_to_3e_from_lepton_input,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "charged_lepton"
_DIPOLE_PROCESS_ID = "L001"
_REQUIRED_EXTRA = "lepton_mass_basis_couplings"
_RS_EW_EXTRA = "rs_ew_couplings"
_OPTIONAL_EW_MASS_EXTRA = "kk_ew_mass_gev"
_LIMIT_ANCHOR_CANDIDATES = ("primary_current_limit",)
_ORIGINAL_EXPERIMENT_CANDIDATES = ("original_experiment",)
_REFERENCE_SCALE_GEV = 3000.0
_UNEVALUATED_REASON = (
    "no mu -> 3e prediction available "
    "(lepton-sector RS contact/dipole inputs not on ParameterPoint)"
)
_UNEVALUATED_NOTES = f"NOT EVALUATED -- {_UNEVALUATED_REASON}"


@dataclass(frozen=True)
class ScalarYAMLInput:
    """Typed scalar loaded from a nested YAML provenance block."""

    block_path: str
    value: float
    source_url: str | None
    supporting_source_url: str | None
    snapshot_path: str | None
    sha256: str | None


@dataclass(frozen=True)
class L002Anchor:
    """Typed L002 inputs: SINDRUM limit plus L001 dipole normalization."""

    experimental: Anchor
    original_experiment: Anchor
    dipole_br_limit: Anchor
    dipole_prefactor_br: ScalarYAMLInput

    @property
    def value(self) -> float:
        return self.experimental.value

    @property
    def budget(self) -> float:
        return self.experimental.value

    @property
    def source_url(self) -> str | None:
        return self.experimental.source_url


def _optional_str(value: Any) -> str | None:
    return None if value is None else str(value)


def _required_mapping(value: Any, *, process_id: str, path: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise AnchorError(f"{process_id}: YAML block {path!r} is not a mapping")
    return value


def _required_float(value: Any, *, process_id: str, path: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise AnchorError(f"{process_id}: YAML field {path!r}={value!r} is not numeric") from exc
    if not math.isfinite(number):
        raise AnchorError(f"{process_id}: YAML field {path!r}={value!r} is not finite")
    return number


def _load_nested_scalar(
    root: Mapping[str, Any],
    path: tuple[str, ...],
    *,
    process_id: str,
) -> ScalarYAMLInput:
    current: Any = root
    walked: list[str] = []
    for part in path:
        walked.append(part)
        current = _required_mapping(
            current,
            process_id=process_id,
            path=".".join(walked[:-1]) or "pdg_or_equivalent",
        ).get(part)
        if current is None:
            raise AnchorError(f"{process_id}: missing YAML block {'.'.join(walked)!r}")
    sub = _required_mapping(current, process_id=process_id, path=".".join(path))
    if "value" not in sub:
        raise AnchorError(f"{process_id}: YAML block {'.'.join(path)!r} has no 'value'")
    return ScalarYAMLInput(
        block_path=".".join(path),
        value=_required_float(
            sub["value"],
            process_id=process_id,
            path=f"{'.'.join(path)}.value",
        ),
        source_url=_optional_str(sub.get("source_url")),
        supporting_source_url=_optional_str(sub.get("supporting_source_url")),
        snapshot_path=_optional_str(sub.get("snapshot_path")),
        sha256=_optional_str(sub.get("sha256")),
    )


def _load_l002_anchor(process_id: str) -> L002Anchor:
    experimental = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=_LIMIT_ANCHOR_CANDIDATES,
    )
    original = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=_ORIGINAL_EXPERIMENT_CANDIDATES,
    )
    if experimental.value <= 0.0:
        raise AnchorError(f"{process_id}: mu -> 3e branching-ratio limit must be positive")
    if not math.isclose(experimental.value, original.value, rel_tol=0.0, abs_tol=0.0):
        raise AnchorError(
            f"{process_id}: primary_current_limit and original_experiment differ"
        )

    dipole_limit = load_anchor(
        _DIPOLE_PROCESS_ID,
        family=_FAMILY,
        candidates=_LIMIT_ANCHOR_CANDIDATES,
    )
    l001_pdg = load_pdg_block(_DIPOLE_PROCESS_ID, family=_FAMILY)
    dipole_prefactor = _load_nested_scalar(
        l001_pdg,
        ("repo_default", "prefac_br"),
        process_id=_DIPOLE_PROCESS_ID,
    )
    if dipole_limit.value <= 0.0 or dipole_prefactor.value <= 0.0:
        raise AnchorError("L002: L001 dipole normalization inputs must be positive")
    return L002Anchor(
        experimental=experimental,
        original_experiment=original,
        dipole_br_limit=dipole_limit,
        dipole_prefactor_br=dipole_prefactor,
    )


@register
class Constraint:
    """Catalogued ``mu -> 3e`` charged-LFV three-body constraint."""

    process_id = "L002"
    severity = Severity.HARD
    observable = "BR(mu -> 3e)"

    def __init__(self) -> None:
        self.anchor = _load_l002_anchor(self.process_id)
        self.sm_inputs = lfv_three_body_default_sm_inputs()

    def _unevaluated_result(
        self,
        *,
        diagnostics: Mapping[str, Any],
    ) -> ConstraintResult:
        merged = {
            "evaluated": False,
            "unevaluated_reason": _UNEVALUATED_REASON,
            "passes_semantics": (
                "non-vetoing only; no BR(mu -> 3e) prediction was evaluated"
            ),
            "budget_source": self.anchor.source_url,
            "experimental_block": self.anchor.experimental.block_key,
            "original_experiment_block": self.anchor.original_experiment.block_key,
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
            result = mu_to_3e_from_lepton_input(
                _adapter_input(lepton_input, rs_ew_couplings),
                br_limit=self.anchor.budget,
                dipole_br_limit=self.anchor.dipole_br_limit.value,
                dipole_prefactor_br=self.anchor.dipole_prefactor_br.value,
                reference_scale_gev=_REFERENCE_SCALE_GEV,
                m_kk_gev=None if kk_ew_mass is None else float(kk_ew_mass),
                inputs=self.sm_inputs,
            )
        except (KeyError, TypeError, ValueError) as exc:
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
                "operator_convention": LFV_THREE_BODY_OPERATOR_CONVENTION,
                "sm_prediction_policy": (
                    "Charged-LFV SM rate is negligible; L002 is applied as a "
                    "pure-NP branching-fraction upper bound."
                ),
                "budget_source": self.anchor.source_url,
                "experimental_block": self.anchor.experimental.block_key,
                "original_experiment_block": self.anchor.original_experiment.block_key,
                "dipole_reused_process_id": _DIPOLE_PROCESS_ID,
                "dipole_limit_source": self.anchor.dipole_br_limit.source_url,
                "dipole_prefactor_block": self.anchor.dipole_prefactor_br.block_path,
                "reference_scale_gev": float(_REFERENCE_SCALE_GEV),
                "kk_ew_mass_extra_used": kk_ew_mass is not None,
                "rs_ew_couplings_extra_present": rs_ew_couplings is not None,
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
                "Pure-NP BR(mu -> 3e) from L001 dipole reuse plus documented "
                "rigorous tree light-Z contact and deferred box inputs, including the "
                "documented dipole-contact interference envelope; HARD budget "
                "is the SINDRUM/PDG branching-fraction limit from L002.yaml. "
                "Dipole, box, heavy-neutral, and unresolved interference pieces "
                "remain NEEDS-HUMAN-PHYSICS."
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


def _is_unevaluated_legacy_proxy_exception(exc: Exception, lepton_input: Any) -> bool:
    message = str(exc)
    if isinstance(exc, ValueError):
        return (
            "pinned to initial_flavor='mu'" in message
            or "pinned to final_flavor='e'" in message
        )
    if not isinstance(exc, TypeError):
        return False
    if (
        "must provide a dipole source and/or lfv three-body tree contact inputs"
        not in message
    ):
        return False
    return _has_finite_legacy_overlap_like_input(lepton_input)


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
