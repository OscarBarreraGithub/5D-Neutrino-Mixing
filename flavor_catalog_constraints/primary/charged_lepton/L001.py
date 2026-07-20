"""L001 - charged-LFV radiative muon decay ``mu -> e gamma``.

Physics
-------
The Standard Model rate is negligible for catalog purposes, so L001 is a pure
new-physics dipole bound:

    BR_NP = prefactor_br * |(Y_N_bar Y_N_bar^dagger)_{e mu}|^2
            * (3 TeV / M_KK)^4.

The lower-level dipole check is the existing
``flavorConstraints.muToEGamma.check_mu_to_e_gamma_raw`` machinery, reached only
through ``flavor_catalog_constraints.physics_adapters.lepton``.  The catalog
hard veto is the branching fraction; the paper coefficient ``C=0.02`` is kept
only as a dipole-RHS diagnostic.

Catalog sidecar
---------------
``flavor_catalog/processes/charged_lepton/L001.yaml`` is the source of truth
for the MEG II branching-fraction limit and the repo-local NDA prefactor /
coefficient provenance.  Numeric anchors are loaded and consistency-checked
from that sidecar, not hardcoded here.
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
from flavor_catalog_constraints.physics_adapters.lepton import (
    LMFVLeptonParameters,
    mu_to_e_gamma_coefficient_from_limit,
    mu_to_e_gamma_from_lepton_input,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "charged_lepton"
_REQUIRED_EXTRA = "lepton_lmfv_parameters"
_OPTIONAL_EW_MASS_EXTRA = "kk_ew_mass_gev"
_LIMIT_ANCHOR_CANDIDATES = ("primary_current_limit",)
_REFERENCE_SCALE_GEV = 3000.0
_REFERENCE_SCALE_SOURCE = (
    "Perez-Randall NDA convention in flavor_catalog/processes/charged_lepton/"
    "L001.tex and L001.yaml paper_era_reference"
)
_UNEVALUATED_REASON = (
    "no lepton dipole prediction available "
    "(LMFV lepton carrier not on ParameterPoint)"
)
_UNEVALUATED_NOTES = f"NOT EVALUATED — {_UNEVALUATED_REASON}"


@dataclass(frozen=True)
class ScalarYAMLInput:
    """Typed scalar loaded from a nested L001 YAML provenance block."""

    block_path: str
    value: float
    source_url: str | None
    supporting_source_url: str | None
    snapshot_path: str | None
    sha256: str | None


@dataclass(frozen=True)
class L001Anchor:
    """Typed L001 inputs: MEG II limit plus NDA prefactor/coefficient."""

    experimental: Anchor
    repo_default_br_limit: ScalarYAMLInput
    prefactor_br: ScalarYAMLInput
    lfv_c: ScalarYAMLInput
    c_paper: ScalarYAMLInput

    @property
    def value(self) -> float:
        """MEG II upper limit on the branching fraction."""
        return self.experimental.value

    @property
    def budget(self) -> float:
        """HARD veto budget for the pure-NP branching fraction."""
        return self.experimental.value

    @property
    def source_url(self) -> str | None:
        """Primary experimental source URL."""
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


def _assert_close(
    *,
    process_id: str,
    label: str,
    actual: float,
    expected: float,
    rel_tol: float = 1.0e-12,
) -> None:
    if not math.isclose(actual, expected, rel_tol=rel_tol, abs_tol=0.0):
        raise AnchorError(
            f"{process_id}: {label} mismatch, got {actual:.17g}, "
            f"expected {expected:.17g}"
        )


def _load_l001_anchor(process_id: str) -> L001Anchor:
    experimental = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=_LIMIT_ANCHOR_CANDIDATES,
    )
    if experimental.value <= 0.0:
        raise AnchorError(f"{process_id}: MEG II branching-ratio limit must be positive")

    pdg = load_pdg_block(process_id, family=_FAMILY)
    br_limit = _load_nested_scalar(
        pdg,
        ("repo_default", "br_limit"),
        process_id=process_id,
    )
    prefactor = _load_nested_scalar(
        pdg,
        ("repo_default", "prefac_br"),
        process_id=process_id,
    )
    lfv_c = _load_nested_scalar(
        pdg,
        ("repo_default", "lfv_C"),
        process_id=process_id,
    )
    c_paper = _load_nested_scalar(
        pdg,
        ("repo_default", "c_paper"),
        process_id=process_id,
    )

    if prefactor.value <= 0.0 or lfv_c.value <= 0.0 or c_paper.value <= 0.0:
        raise AnchorError(f"{process_id}: L001 prefactor and LFV coefficients must be positive")
    _assert_close(
        process_id=process_id,
        label="repo_default.br_limit vs primary_current_limit",
        actual=br_limit.value,
        expected=experimental.value,
    )
    expected_c = mu_to_e_gamma_coefficient_from_limit(
        experimental.value,
        prefactor.value,
    )
    _assert_close(
        process_id=process_id,
        label="repo_default.lfv_C vs sqrt(BR_limit/prefac_br)",
        actual=lfv_c.value,
        expected=expected_c,
    )

    return L001Anchor(
        experimental=experimental,
        repo_default_br_limit=br_limit,
        prefactor_br=prefactor,
        lfv_c=lfv_c,
        c_paper=c_paper,
    )


@register
class Constraint:
    """Catalogued ``mu -> e gamma`` charged-LFV dipole constraint."""

    process_id = "L001"
    severity = Severity.HARD
    observable = "BR(mu -> e gamma)"

    def __init__(self) -> None:
        self.anchor = _load_l001_anchor(self.process_id)

    def _unevaluated_result(
        self,
        *,
        diagnostics: Mapping[str, Any],
    ) -> ConstraintResult:
        merged_diagnostics = {
            "evaluated": False,
            "unevaluated_reason": _UNEVALUATED_REASON,
            "passes_semantics": (
                "non-vetoing only; no BR(mu -> e gamma) prediction was evaluated"
            ),
            "budget_source": self.anchor.source_url,
            "experimental_block": self.anchor.experimental.block_key,
            **dict(diagnostics),
        }
        return ConstraintResult(
            process_id=self.process_id,
            severity=self.severity,
            passes=True,
            predicted=None,
            sm_prediction=None,
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
        if not isinstance(lepton_input, LMFVLeptonParameters):
            return self._unevaluated_result(
                diagnostics={
                    "invalid_extra": _REQUIRED_EXTRA,
                    "exception_type": "TypeError",
                    "exception": (
                        f"{_REQUIRED_EXTRA} must be an LMFVLeptonParameters "
                        "carrier"
                    ),
                },
            )

        kk_ew_mass = point.get_extra(_OPTIONAL_EW_MASS_EXTRA)
        kk_ew_mass_value = None if kk_ew_mass is None else float(kk_ew_mass)
        kk_ew_mass_consistent = (
            None
            if kk_ew_mass_value is None
            else math.isclose(
                kk_ew_mass_value,
                float(lepton_input.M_KK_gev),
                rel_tol=1.0e-12,
                abs_tol=1.0e-9,
            )
        )
        try:
            result = mu_to_e_gamma_from_lepton_input(
                lepton_input,
                br_limit=self.anchor.budget,
                prefactor_br=self.anchor.prefactor_br.value,
                c_lfv=self.anchor.c_paper.value,
                reference_scale_gev=_REFERENCE_SCALE_GEV,
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
                "lmfv_model": "Perez-Randall LMFV NDA",
                "dipole_lhs": float(result.dipole_lhs),
                "dipole_rhs": float(result.dipole_rhs),
                "dipole_ratio_to_bound": float(result.dipole_ratio_to_bound),
                "lfv_coefficient": float(result.c_lfv),
                "c_lfv_role": "dipole_rhs_diagnostic_only",
                "prefactor_br": float(result.prefactor_br),
                "br_limit": float(result.br_limit),
                "br_ratio_to_limit": float(result.ratio_to_limit),
                "m_kk_gev": float(result.m_kk_gev),
                "reference_scale_gev": float(result.reference_scale_gev),
                "reference_scale_source": _REFERENCE_SCALE_SOURCE,
                "off_diagonal_12": complex(result.off_diagonal_12),
                "product_matrix": result.product_matrix,
                "input_kind": result.input_kind,
                "extra_used": _REQUIRED_EXTRA,
                "used_proxy": bool(result.used_proxy),
                "sm_prediction_policy": (
                    "Charged-LFV SM rate is negligible; budget is applied to "
                    "the pure NP dipole branching fraction."
                ),
                "budget_source": self.anchor.source_url,
                "experimental_block": self.anchor.experimental.block_key,
                "repo_default_br_limit_block": self.anchor.repo_default_br_limit.block_path,
                "prefactor_block": self.anchor.prefactor_br.block_path,
                "lfv_c_block": self.anchor.lfv_c.block_path,
                "paper_c_lfv": float(self.anchor.c_paper.value),
                "kk_ew_mass_extra_present": kk_ew_mass is not None,
                "kk_ew_mass_extra_used": False,
                "kk_ew_mass_extra_value": kk_ew_mass_value,
                "kk_ew_mass_extra_consistent": kk_ew_mass_consistent,
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
                "Pure-NP BR(mu -> e gamma) LMFV NDA prediction using the "
                "existing flavorConstraints.muToEGamma machinery through the "
                "lepton adapter; HARD budget is the MEG II branching-fraction "
                "limit from L001.yaml."
            ),
            diagnostics=diagnostics,
        )
