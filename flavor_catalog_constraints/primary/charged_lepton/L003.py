"""L003 - coherent mu->e conversion in aluminum.

Physics
-------
The Standard Model rate is negligible for catalog purposes, so L003 is a
pure-new-physics conversion-rate bound,

    CR(mu Al -> e Al) = omega_conv / Gamma_capture,

with the Kitano-Koike-Okada overlap-integral formula.  The prediction combines

* the L001 ``mu -> e gamma`` dipole adapter,
* scalar/vector proton and neutron coefficient proxies in the KKO convention,
* aluminum nuclear overlaps and capture rate from the reusable
  ``quarkConstraints.mu_e_conversion`` core.

The lower-level formula is reached only through
``flavor_catalog_constraints.physics_adapters.mu_e_conversion``.

NEEDS-HUMAN-PHYSICS
-------------------
The current ``ParameterPoint`` does not carry full RS charged-lepton
neutral-current, Higgs/scalar, EW KK/Z/Z', or lepton-quark matching data.  The
scalar/vector terms therefore use explicit low-energy coefficient proxies and
are flagged in diagnostics.

Catalog sidecar
---------------
``flavor_catalog/processes/charged_lepton/L003.yaml`` is the source of truth
for the SINDRUM II current benchmark and the Mu2e aluminum projected limit.
There is no published direct aluminum limit in the sidecar.  This constraint
therefore reports the Mu2e full-exposure expected 90% C.L. aluminum projection
as a SOFT sensitivity-reach comparison, while reporting the observed SINDRUM II
Au benchmark separately in diagnostics.
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
from flavor_catalog_constraints.physics_adapters.mu_e_conversion import (
    MU_E_CONVERSION_OPERATOR_CONVENTION,
    MU_E_CONVERSION_PROXY_V1,
    mu_e_conversion_aluminum_nuclear_inputs,
    mu_e_conversion_from_lepton_input,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "charged_lepton"
_DIPOLE_PROCESS_ID = "L001"
_REQUIRED_EXTRA = "lepton_mass_basis_couplings"
_OPTIONAL_EW_MASS_EXTRA = "kk_ew_mass_gev"
_CURRENT_WORLD_LIMIT_CANDIDATES = ("current_world_limit",)
_DIPOLE_LIMIT_CANDIDATES = ("primary_current_limit",)
_REFERENCE_SCALE_GEV = 3000.0
_BUDGET_PROJECTION_SOURCE = "Mu2e experiment at Fermilab"
_BUDGET_PROJECTION_VALUE_KEY = "expected_upper_limit_90cl"
_UNEVALUATED_REASON = (
    "no mu->e conversion prediction available "
    "(lepton-sector scalar/vector/dipole inputs not on ParameterPoint)"
)
_UNEVALUATED_NOTES = f"NOT EVALUATED -- {_UNEVALUATED_REASON}"
_BUDGET_POLICY = (
    "Mu2e full-exposure expected 90% C.L. aluminum projection from L003.yaml "
    "is used only as a SOFT sensitivity-reach comparison; SINDRUM II Au is "
    "loaded as the current observed non-aluminum benchmark and is not treated "
    "as a direct aluminum limit or hard veto."
)


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
class AluminumProjectionLimit:
    """Typed view over one L003 aluminum projection item."""

    list_path: str
    source: str
    year: int | None
    target_material: str
    observable: str
    value_key: str
    value: float
    single_event_sensitivity: float | None
    relation: str | None
    measurement_status: str | None
    source_url: str | None
    snapshot_path: str | None
    sha256: str | None


@dataclass(frozen=True)
class L003Anchor:
    """Typed L003 inputs: observed benchmark, Mu2e budget, and L001 dipole."""

    current_world_limit: Anchor
    aluminum_budget: AluminumProjectionLimit
    dipole_br_limit: Anchor
    dipole_prefactor_br: ScalarYAMLInput
    direct_aluminum_limit_available: bool
    aluminum_status_note: str | None

    @property
    def value(self) -> float:
        return self.aluminum_budget.value

    @property
    def budget(self) -> float:
        return self.aluminum_budget.value

    @property
    def source_url(self) -> str | None:
        return self.aluminum_budget.source_url

    @property
    def current_benchmark_value(self) -> float:
        return self.current_world_limit.value


def _optional_str(value: Any) -> str | None:
    return None if value is None else str(value)


def _optional_int(value: Any) -> int | None:
    if value is None:
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


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


def _optional_float(value: Any, *, process_id: str, path: str) -> float | None:
    if value is None:
        return None
    return _required_float(value, process_id=process_id, path=path)


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


def _load_aluminum_budget_projection(process_id: str) -> AluminumProjectionLimit:
    pdg = load_pdg_block(process_id, family=_FAMILY)
    projections = pdg.get("aluminum_projections")
    if not isinstance(projections, list) or not projections:
        raise AnchorError(f"{process_id}: missing nonempty aluminum_projections list")

    for index, item in enumerate(projections):
        sub = _required_mapping(
            item,
            process_id=process_id,
            path=f"aluminum_projections[{index}]",
        )
        if sub.get("source") != _BUDGET_PROJECTION_SOURCE:
            continue
        if _BUDGET_PROJECTION_VALUE_KEY not in sub:
            raise AnchorError(
                f"{process_id}: {_BUDGET_PROJECTION_SOURCE} projection has no "
                f"{_BUDGET_PROJECTION_VALUE_KEY!r}"
            )
        value = _required_float(
            sub[_BUDGET_PROJECTION_VALUE_KEY],
            process_id=process_id,
            path=f"aluminum_projections[{index}].{_BUDGET_PROJECTION_VALUE_KEY}",
        )
        if value <= 0.0:
            raise AnchorError(f"{process_id}: Mu2e projected limit must be positive")
        target = _optional_str(sub.get("target_material"))
        if target != "Al":
            raise AnchorError(
                f"{process_id}: Mu2e projection target must be 'Al', got {target!r}"
            )
        return AluminumProjectionLimit(
            list_path=f"aluminum_projections[{index}]",
            source=str(sub["source"]),
            year=_optional_int(sub.get("year")),
            target_material=target,
            observable=str(sub.get("observable", "R_mue")),
            value_key=_BUDGET_PROJECTION_VALUE_KEY,
            value=float(value),
            single_event_sensitivity=_optional_float(
                sub.get("single_event_sensitivity"),
                process_id=process_id,
                path=f"aluminum_projections[{index}].single_event_sensitivity",
            ),
            relation=_optional_str(sub.get("relation")),
            measurement_status=_optional_str(sub.get("measurement_status")),
            source_url=_optional_str(sub.get("source_url")),
            snapshot_path=_optional_str(sub.get("snapshot_path")),
            sha256=_optional_str(sub.get("sha256")),
        )

    raise AnchorError(
        f"{process_id}: no aluminum projection from {_BUDGET_PROJECTION_SOURCE!r}"
    )


def _load_aluminum_status(process_id: str) -> tuple[bool, str | None]:
    pdg = load_pdg_block(process_id, family=_FAMILY)
    status = _required_mapping(
        pdg.get("aluminum_status"),
        process_id=process_id,
        path="aluminum_status",
    )
    return (
        status.get("published_direct_aluminum_limit") is not None,
        _optional_str(status.get("note")),
    )


def _load_l003_anchor(process_id: str) -> L003Anchor:
    current_world_limit = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=_CURRENT_WORLD_LIMIT_CANDIDATES,
    )
    if current_world_limit.value <= 0.0:
        raise AnchorError(f"{process_id}: current world benchmark must be positive")

    dipole_limit = load_anchor(
        _DIPOLE_PROCESS_ID,
        family=_FAMILY,
        candidates=_DIPOLE_LIMIT_CANDIDATES,
    )
    l001_pdg = load_pdg_block(_DIPOLE_PROCESS_ID, family=_FAMILY)
    dipole_prefactor = _load_nested_scalar(
        l001_pdg,
        ("repo_default", "prefac_br"),
        process_id=_DIPOLE_PROCESS_ID,
    )
    if dipole_limit.value <= 0.0 or dipole_prefactor.value <= 0.0:
        raise AnchorError("L003: L001 dipole normalization inputs must be positive")

    direct_available, status_note = _load_aluminum_status(process_id)
    return L003Anchor(
        current_world_limit=current_world_limit,
        aluminum_budget=_load_aluminum_budget_projection(process_id),
        dipole_br_limit=dipole_limit,
        dipole_prefactor_br=dipole_prefactor,
        direct_aluminum_limit_available=direct_available,
        aluminum_status_note=status_note,
    )


@register
class Constraint:
    """Catalogued coherent mu->e conversion constraint for aluminum."""

    process_id = "L003"
    severity = Severity.SOFT
    observable = "CR(mu Al -> e Al)"

    def __init__(self) -> None:
        self.anchor = _load_l003_anchor(self.process_id)
        self.nuclear_inputs = mu_e_conversion_aluminum_nuclear_inputs()

    def _unevaluated_result(
        self,
        *,
        diagnostics: Mapping[str, Any],
    ) -> ConstraintResult:
        merged = {
            "evaluated": False,
            "unevaluated_reason": _UNEVALUATED_REASON,
            "passes_semantics": (
                "non-vetoing only; no CR(mu Al -> e Al) prediction was evaluated"
            ),
            "budget_source": self.anchor.source_url,
            "budget_policy": _BUDGET_POLICY,
            "budget_limit_status": "projected_expected_sensitivity",
            "budget_verdict_role": "SOFT sensitivity reach; no hard veto",
            "current_world_limit": float(self.anchor.current_benchmark_value),
            "current_world_limit_target": (
                self.anchor.current_world_limit.source
                or "SINDRUM II Au benchmark"
            ),
            "current_world_limit_status": "observed_experimental_bound",
            "current_world_limit_verdict_role": (
                "diagnostic observed benchmark only; different Au target"
            ),
            "direct_aluminum_limit_available": (
                self.anchor.direct_aluminum_limit_available
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
            result = mu_e_conversion_from_lepton_input(
                lepton_input,
                conversion_rate_limit=self.anchor.budget,
                dipole_br_limit=self.anchor.dipole_br_limit.value,
                dipole_prefactor_br=self.anchor.dipole_prefactor_br.value,
                reference_scale_gev=_REFERENCE_SCALE_GEV,
                m_kk_gev=None if kk_ew_mass is None else float(kk_ew_mass),
                nuclear_inputs=self.nuclear_inputs,
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
                "needs_human_physics": MU_E_CONVERSION_PROXY_V1,
                "operator_convention": MU_E_CONVERSION_OPERATOR_CONVENTION,
                "sm_prediction_policy": (
                    "Charged-LFV mu->e conversion has negligible SM rate; L003 "
                    "is applied as a pure-NP conversion-rate upper bound."
                ),
                "budget_source": self.anchor.source_url,
                "budget_policy": _BUDGET_POLICY,
                "budget_projection_block": self.anchor.aluminum_budget.list_path,
                "budget_projection_source": self.anchor.aluminum_budget.source,
                "budget_projection_value_key": self.anchor.aluminum_budget.value_key,
                "budget_limit_status": (
                    self.anchor.aluminum_budget.measurement_status
                    or "projected_expected_sensitivity"
                ),
                "budget_verdict_role": "SOFT sensitivity reach; no hard veto",
                "budget_projection_single_event_sensitivity": (
                    self.anchor.aluminum_budget.single_event_sensitivity
                ),
                "current_world_limit": float(self.anchor.current_benchmark_value),
                "current_world_limit_block": (
                    self.anchor.current_world_limit.block_key
                ),
                "current_world_limit_target_material": "Au",
                "current_world_limit_source": self.anchor.current_world_limit.source,
                "current_world_limit_status": "observed_experimental_bound",
                "current_world_limit_verdict_role": (
                    "diagnostic observed benchmark only; different Au target"
                ),
                "direct_aluminum_limit_available": (
                    self.anchor.direct_aluminum_limit_available
                ),
                "aluminum_status_note": self.anchor.aluminum_status_note,
                "dipole_reused_process_id": _DIPOLE_PROCESS_ID,
                "dipole_limit_source": self.anchor.dipole_br_limit.source_url,
                "dipole_prefactor_block": self.anchor.dipole_prefactor_br.block_path,
                "reference_scale_gev": float(_REFERENCE_SCALE_GEV),
                "kk_ew_mass_extra_used": kk_ew_mass is not None,
            }
        )

        return ConstraintResult(
            process_id=self.process_id,
            severity=self.severity,
            passes=bool(result.passes),
            predicted=float(result.conversion_rate),
            sm_prediction=0.0,
            experimental=float(self.anchor.value),
            ratio=float(result.ratio_to_limit),
            budget=float(result.conversion_rate_limit),
            notes=(
                "Pure-NP CR(mu Al -> e Al) from L001 dipole reuse plus "
                "documented scalar/vector nucleon coefficient proxies and "
                "aluminum nuclear overlap integrals; Mu2e aluminum projected "
                "expected upper limit from L003.yaml is a SOFT sensitivity "
                "reach, not an observed veto. Ambiguous dipole-contact phase "
                "uses the lower-envelope verdict and reports the full interval; "
                "scalar/vector matching is flagged NEEDS-HUMAN-PHYSICS."
            ),
            diagnostics=diagnostics,
        )
