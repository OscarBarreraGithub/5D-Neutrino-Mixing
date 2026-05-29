"""L007 - charged-LFV radiative tau decay ``tau -> mu gamma``.

Physics
-------
The Standard Model charged-LFV rate is negligible for catalog purposes, so
L007 is a pure-new-physics dipole bound:

    BR_NP = prefactor_br * |(Y_N_bar Y_N_bar^dagger)_{mu tau}|^2
            * (3 TeV / M_KK)^4.

The calculation reuses the L001 lepton dipole machinery in
``flavorConstraints.muToEGamma`` through the tau-mu adapter
``flavor_catalog_constraints.physics_adapters.lepton_tau_mu``.  The adapter
maps the L001 core's hardcoded ``(e,mu)`` element to the ``(mu,tau)`` element by
permuting charged-lepton rows before calling the existing core.

NEEDS-HUMAN-PHYSICS
-------------------
The current ``ParameterPoint`` can be quark-only and does not carry a full
lepton-sector RS coupling object or tau-specific radiative matching inputs.  If
a caller supplies the explicit ``lepton_mass_basis_couplings`` proxy accepted
by the adapter, the result is flagged ``NEEDS-HUMAN-PHYSICS`` because that
proxy is not a loop-level RS ``tau -> mu gamma`` calculation.

Catalog sidecar
---------------
``flavor_catalog/processes/charged_lepton/L007.yaml`` is the source of truth
for the PDG/Belle branching-fraction limit.  The shared Perez-Randall NDA
prefactor is loaded from the L001 sidecar, as in the existing mu->e conversion
constraints that reuse the L001 dipole normalization.
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
from flavor_catalog_constraints.physics_adapters.lepton_tau_mu import (
    TAU_TO_MU_GAMMA_PROXY_V1,
    tau_to_mu_gamma_from_lepton_input,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "charged_lepton"
_DIPOLE_PROCESS_ID = "L001"
_REQUIRED_EXTRA = "lepton_mass_basis_couplings"
_OPTIONAL_EW_MASS_EXTRA = "kk_ew_mass_gev"
_LIMIT_ANCHOR_CANDIDATES = ("primary_current_limit",)
_PRIMARY_EXPERIMENT_CANDIDATES = ("primary_experiment",)
_PRIOR_EXPERIMENT_CANDIDATES = ("prior_experimental_limit",)
_BELLE_II_PROJECTION_CANDIDATES = ("belle_ii_projection",)
_REFERENCE_SCALE_GEV = 3000.0
_REFERENCE_SCALE_SOURCE = (
    "Perez-Randall NDA convention reused from the L001 lepton dipole "
    "machinery; see flavor_catalog/processes/charged_lepton/L001.yaml "
    "paper_era_reference"
)
_UNEVALUATED_REASON = (
    "no tau->mu gamma dipole prediction available "
    "(lepton-sector RS couplings not on ParameterPoint)"
)
_UNEVALUATED_NOTES = f"NOT EVALUATED - {_UNEVALUATED_REASON}"


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
class L007Anchor:
    """Typed L007 inputs: current tau->mu gamma limit plus dipole prefactor."""

    experimental: Anchor
    primary_experiment: Anchor
    prior_experimental_limit: Anchor
    belle_ii_projection: Anchor
    dipole_prefactor_br: ScalarYAMLInput

    @property
    def value(self) -> float:
        """Current observed upper limit on the branching fraction."""
        return self.experimental.value

    @property
    def budget(self) -> float:
        """HARD veto budget for the pure-NP branching fraction."""
        return self.experimental.value

    @property
    def source_url(self) -> str | None:
        """Primary PDG/equivalent source URL."""
        return self.experimental.source_url

    @property
    def belle_source_url(self) -> str | None:
        """Belle full-data source URL for the current limit."""
        return self.primary_experiment.source_url


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


def _assert_positive_limit(anchor: Anchor, *, process_id: str, label: str) -> None:
    if not math.isfinite(anchor.value) or anchor.value <= 0.0:
        raise AnchorError(f"{process_id}: {label} branching-ratio limit must be positive")


def _load_l007_anchor(process_id: str) -> L007Anchor:
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
    prior_experimental_limit = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=_PRIOR_EXPERIMENT_CANDIDATES,
    )
    belle_ii_projection = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=_BELLE_II_PROJECTION_CANDIDATES,
        value_key="projected_limit",
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
        prior_experimental_limit,
        process_id=process_id,
        label="prior_experimental_limit",
    )
    _assert_positive_limit(
        belle_ii_projection,
        process_id=process_id,
        label="belle_ii_projection.projected_limit",
    )
    if not math.isclose(experimental.value, primary_experiment.value, rel_tol=0.0, abs_tol=0.0):
        raise AnchorError(
            f"{process_id}: primary_current_limit and primary_experiment differ"
        )

    l001_pdg = load_pdg_block(_DIPOLE_PROCESS_ID, family=_FAMILY)
    dipole_prefactor = _load_nested_scalar(
        l001_pdg,
        ("repo_default", "prefac_br"),
        process_id=_DIPOLE_PROCESS_ID,
    )
    if dipole_prefactor.value <= 0.0:
        raise AnchorError("L007: L001 dipole prefactor must be positive")

    return L007Anchor(
        experimental=experimental,
        primary_experiment=primary_experiment,
        prior_experimental_limit=prior_experimental_limit,
        belle_ii_projection=belle_ii_projection,
        dipole_prefactor_br=dipole_prefactor,
    )


@register
class Constraint:
    """Catalogued ``tau -> mu gamma`` charged-LFV dipole constraint."""

    process_id = "L007"
    severity = Severity.HARD
    observable = "BR(tau -> mu gamma)"

    def __init__(self) -> None:
        self.anchor = _load_l007_anchor(self.process_id)

    def _unevaluated_result(
        self,
        *,
        diagnostics: Mapping[str, Any],
    ) -> ConstraintResult:
        merged_diagnostics = {
            "evaluated": False,
            "unevaluated_reason": _UNEVALUATED_REASON,
            "passes_semantics": (
                "non-vetoing only; no BR(tau -> mu gamma) prediction was evaluated"
            ),
            "sm_prediction_policy": (
                "Charged-LFV SM rate is negligible; L007 is applied as a "
                "pure-NP branching-fraction upper bound when proxy inputs exist."
            ),
            "budget_source": self.anchor.source_url,
            "belle_source": self.anchor.belle_source_url,
            "experimental_block": self.anchor.experimental.block_key,
            "primary_experiment_block": self.anchor.primary_experiment.block_key,
            "prior_experimental_limit": float(self.anchor.prior_experimental_limit.value),
            "prior_experimental_source": self.anchor.prior_experimental_limit.source_url,
            "belle_ii_projection": float(self.anchor.belle_ii_projection.value),
            "belle_ii_projection_source": self.anchor.belle_ii_projection.source_url,
            "belle_ii_projection_not_used_for_hard_veto": True,
            "prefactor_block": self.anchor.dipole_prefactor_br.block_path,
            "prefactor_source": self.anchor.dipole_prefactor_br.source_url,
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
                diagnostics={"missing_extra": _REQUIRED_EXTRA},
            )

        kk_ew_mass = point.get_extra(_OPTIONAL_EW_MASS_EXTRA)
        try:
            result = tau_to_mu_gamma_from_lepton_input(
                lepton_input,
                br_limit=self.anchor.budget,
                prefactor_br=self.anchor.dipole_prefactor_br.value,
                reference_scale_gev=_REFERENCE_SCALE_GEV,
                m_kk_gev=None if kk_ew_mass is None else float(kk_ew_mass),
            )
        except (AttributeError, KeyError, TypeError, ValueError) as exc:
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
                "needs_human_physics": TAU_TO_MU_GAMMA_PROXY_V1,
                "dipole_lhs": float(result.dipole_lhs),
                "dipole_rhs": float(result.dipole_rhs),
                "dipole_ratio_to_bound": float(result.dipole_ratio_to_bound),
                "lfv_coefficient": float(result.c_lfv),
                "prefactor_br": float(result.prefactor_br),
                "m_kk_gev": float(result.m_kk_gev),
                "reference_scale_gev": float(result.reference_scale_gev),
                "reference_scale_source": _REFERENCE_SCALE_SOURCE,
                "off_diagonal_23": complex(result.off_diagonal_23),
                "off_diagonal_32": complex(result.off_diagonal_32),
                "product_matrix": result.product_matrix,
                "row_permutation": result.row_permutation,
                "input_kind": result.input_kind,
                "used_proxy": bool(result.used_proxy),
                "sm_prediction_policy": (
                    "Charged-LFV SM rate is negligible; budget is applied to "
                    "the pure NP tau->mu dipole branching fraction."
                ),
                "budget_source": self.anchor.source_url,
                "belle_source": self.anchor.belle_source_url,
                "experimental_block": self.anchor.experimental.block_key,
                "primary_experiment_block": self.anchor.primary_experiment.block_key,
                "prior_experimental_limit": float(
                    self.anchor.prior_experimental_limit.value
                ),
                "prior_experimental_source": (
                    self.anchor.prior_experimental_limit.source_url
                ),
                "belle_ii_projection": float(self.anchor.belle_ii_projection.value),
                "belle_ii_projection_source": self.anchor.belle_ii_projection.source_url,
                "belle_ii_projection_not_used_for_hard_veto": True,
                "prefactor_block": self.anchor.dipole_prefactor_br.block_path,
                "prefactor_source": self.anchor.dipole_prefactor_br.source_url,
                "kk_ew_mass_extra_used": kk_ew_mass is not None,
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
                "Pure-NP BR(tau -> mu gamma) dipole check using the existing "
                "flavorConstraints.muToEGamma machinery through the tau-mu "
                "lepton adapter. The HARD budget is the PDG/Belle observed "
                "branching-fraction limit from L007.yaml; the tau-mu dipole "
                "proxy is flagged NEEDS-HUMAN-PHYSICS."
            ),
            diagnostics=diagnostics,
        )
