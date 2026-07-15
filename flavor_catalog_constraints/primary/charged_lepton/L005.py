"""L005 - coherent mu->e conversion in titanium.

Physics
-------
The Standard Model charged-LFV conversion rate is negligible for catalog
purposes, so L005 is a pure-new-physics bound,

    CR(mu Ti -> e Ti) = omega_conv / Gamma_capture.

The prediction reuses the Kitano-Koike-Okada overlap-integral machinery built
for L003/L004, reached only through
``flavor_catalog_constraints.physics_adapters.mu_e_conversion``.  This file
selects the titanium target nuclear inputs and the observed SINDRUM II Ti
limit from the L005 catalog sidecar.

NEEDS-HUMAN-PHYSICS
-------------------
The tree-level light-Z vector contact is rigorous when ``rs_ew_couplings`` is
present and is zero for the Phase-4a diagonal charged-lepton fit.  Scalar,
dipole, heavy-neutral, and interference pieces remain deferred.  Unknown
dipole-contact phases use the conservative destructive-interference envelope
for the HARD verdict and report the full interval.

Catalog sidecar
---------------
``flavor_catalog/processes/charged_lepton/L005.yaml`` is the source of truth
for the observed SINDRUM II titanium conversion bound.  L001.yaml supplies the
dipole-parent normalization reused by the mu->e conversion adapter.
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
    MU_E_CONVERSION_DEFERRED_PIECES_V1,
    MU_E_CONVERSION_OPERATOR_CONVENTION,
    MU_E_CONVERSION_VECTOR_TREE_RIGOROUS_V1,
    mu_e_conversion_from_lepton_input,
    mu_e_conversion_titanium_nuclear_inputs,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "charged_lepton"
_DIPOLE_PROCESS_ID = "L001"
_REQUIRED_EXTRA = "lepton_mass_basis_couplings"
_RS_EW_EXTRA = "rs_ew_couplings"
_OPTIONAL_EW_MASS_EXTRA = "kk_ew_mass_gev"
_LIMIT_ANCHOR_CANDIDATES = ("primary_current_limit",)
_DIPOLE_LIMIT_CANDIDATES = ("primary_current_limit",)
_REFERENCE_SCALE_GEV = 3000.0
_UNEVALUATED_REASON = (
    "no mu->e conversion prediction available "
    "(lepton-sector scalar/vector/dipole inputs not on ParameterPoint)"
)
_UNEVALUATED_NOTES = f"NOT EVALUATED -- {_UNEVALUATED_REASON}"
_BUDGET_POLICY = (
    "Observed SINDRUM II Ti 90% C.L. upper limit from L005.yaml; HARD veto "
    "uses the pure-NP CR(mu Ti -> e Ti) prediction normalized to the titanium "
    "capture rate."
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
class PublicationYAMLInput:
    """Typed provenance view over the SINDRUM II publication block."""

    block_path: str
    source: str | None
    year: int | None
    title: str | None
    journal: str | None
    doi: str | None
    source_url: str | None
    snapshot_path: str | None
    sha256: str | None


@dataclass(frozen=True)
class L005Anchor:
    """Typed L005 inputs: SINDRUM II Ti limit plus L001 dipole normalization."""

    experimental: Anchor
    experimental_publication: PublicationYAMLInput
    target_material: str
    relation: str | None
    confidence_level: str | None
    experiment: str | None
    technique: str | None
    pdg_note: str | None
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


def _load_publication(
    root: Mapping[str, Any],
    *,
    process_id: str,
) -> PublicationYAMLInput:
    block_key = "experimental_publication"
    sub = _required_mapping(
        root.get(block_key),
        process_id=process_id,
        path=block_key,
    )
    return PublicationYAMLInput(
        block_path=block_key,
        source=_optional_str(sub.get("source")),
        year=_optional_int(sub.get("year")),
        title=_optional_str(sub.get("title")),
        journal=_optional_str(sub.get("journal")),
        doi=_optional_str(sub.get("doi")),
        source_url=_optional_str(sub.get("source_url")),
        snapshot_path=_optional_str(sub.get("snapshot_path")),
        sha256=_optional_str(sub.get("sha256")),
    )


def _load_l005_anchor(process_id: str) -> L005Anchor:
    experimental = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=_LIMIT_ANCHOR_CANDIDATES,
    )
    if experimental.value <= 0.0:
        raise AnchorError(f"{process_id}: Ti conversion limit must be positive")

    pdg = load_pdg_block(process_id, family=_FAMILY)
    experimental_sub = _required_mapping(
        pdg.get(experimental.block_key),
        process_id=process_id,
        path=experimental.block_key,
    )
    target = _optional_str(experimental_sub.get("target_material"))
    if target != "Ti":
        raise AnchorError(f"{process_id}: primary limit target must be 'Ti', got {target!r}")

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
        raise AnchorError("L005: L001 dipole normalization inputs must be positive")

    return L005Anchor(
        experimental=experimental,
        experimental_publication=_load_publication(pdg, process_id=process_id),
        target_material=target,
        relation=_optional_str(experimental_sub.get("relation")),
        confidence_level=_optional_str(experimental_sub.get("confidence_level")),
        experiment=_optional_str(experimental_sub.get("experiment")),
        technique=_optional_str(experimental_sub.get("technique")),
        pdg_note=_optional_str(experimental_sub.get("pdg_note")),
        dipole_br_limit=dipole_limit,
        dipole_prefactor_br=dipole_prefactor,
    )


@register
class Constraint:
    """Catalogued coherent mu->e conversion constraint for titanium."""

    process_id = "L005"
    severity = Severity.HARD
    observable = "CR(mu Ti -> e Ti)"

    def __init__(self) -> None:
        self.anchor = _load_l005_anchor(self.process_id)
        self.nuclear_inputs = mu_e_conversion_titanium_nuclear_inputs()

    def _unevaluated_result(
        self,
        *,
        diagnostics: Mapping[str, Any],
    ) -> ConstraintResult:
        merged = {
            "evaluated": False,
            "unevaluated_reason": _UNEVALUATED_REASON,
            "passes_semantics": (
                "non-vetoing only; no CR(mu Ti -> e Ti) prediction was evaluated"
            ),
            "budget_source": self.anchor.source_url,
            "budget_policy": _BUDGET_POLICY,
            "budget_limit_status": "observed_experimental_bound",
            "budget_verdict_role": "HARD observed upper-limit veto when evaluated",
            "target_material": self.anchor.target_material,
            "experimental_block": self.anchor.experimental.block_key,
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
            result = mu_e_conversion_from_lepton_input(
                _adapter_input(lepton_input, rs_ew_couplings),
                conversion_rate_limit=self.anchor.budget,
                dipole_br_limit=self.anchor.dipole_br_limit.value,
                dipole_prefactor_br=self.anchor.dipole_prefactor_br.value,
                reference_scale_gev=_REFERENCE_SCALE_GEV,
                m_kk_gev=None if kk_ew_mass is None else float(kk_ew_mass),
                nuclear_inputs=self.nuclear_inputs,
            )
        except (KeyError, TypeError, ValueError) as exc:
            diagnostics: dict[str, object] = {
                "exception_type": type(exc).__name__,
                "exception": str(exc),
            }
            if _is_unevaluated_legacy_vector_exception(exc, lepton_input):
                diagnostics["legacy_vector_proxy_ignored"] = True
            else:
                diagnostics["invalid_extra"] = _REQUIRED_EXTRA
            return self._unevaluated_result(diagnostics=diagnostics)

        diagnostics = dict(result.diagnostics)
        diagnostics.update(
            {
                "evaluated": True,
                "needs_human_physics": MU_E_CONVERSION_DEFERRED_PIECES_V1,
                "vector_tree_matching": diagnostics.get(
                    "vector_tree_matching",
                    MU_E_CONVERSION_VECTOR_TREE_RIGOROUS_V1,
                ),
                "operator_convention": MU_E_CONVERSION_OPERATOR_CONVENTION,
                "sm_prediction_policy": (
                    "Charged-LFV mu->e conversion has negligible SM rate; L005 "
                    "is applied as a pure-NP conversion-rate upper bound."
                ),
                "budget_source": self.anchor.source_url,
                "budget_policy": _BUDGET_POLICY,
                "budget_limit_status": "observed_experimental_bound",
                "budget_verdict_role": "HARD observed upper-limit veto",
                "experimental_block": self.anchor.experimental.block_key,
                "experimental_relation": self.anchor.relation,
                "experimental_confidence_level": self.anchor.confidence_level,
                "experimental_source": self.anchor.experimental.source,
                "experimental_experiment": self.anchor.experiment,
                "experimental_technique": self.anchor.technique,
                "experimental_pdg_note": self.anchor.pdg_note,
                "experimental_publication_block": (
                    self.anchor.experimental_publication.block_path
                ),
                "experimental_publication_doi": (
                    self.anchor.experimental_publication.doi
                ),
                "target_material": self.anchor.target_material,
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
            predicted=float(result.conversion_rate),
            sm_prediction=0.0,
            experimental=float(self.anchor.value),
            ratio=float(result.ratio_to_limit),
            budget=float(result.conversion_rate_limit),
            notes=(
                "Pure-NP CR(mu Ti -> e Ti) from L001 dipole reuse plus "
                "rigorous tree light-Z vector coefficients, deferred scalar inputs, and "
                "titanium nuclear overlap integrals; observed SINDRUM II Ti "
                "limit from L005.yaml is a HARD veto. Ambiguous dipole-contact "
                "phase uses the destructive-interference lower envelope for "
                "the verdict and reports the full interval; scalar, dipole, "
                "and interference matching remain NEEDS-HUMAN-PHYSICS."
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


_LEGACY_VECTOR_COEFFICIENT_KEYS = (
    "g_lv_p",
    "g_LV_p",
    "g_lv_proton",
    "g_LV_proton",
    "g_lv_n",
    "g_LV_n",
    "g_lv_neutron",
    "g_LV_neutron",
    "g_rv_p",
    "g_RV_p",
    "g_rv_proton",
    "g_RV_proton",
    "g_rv_n",
    "g_RV_n",
    "g_rv_neutron",
    "g_RV_neutron",
)


def _is_unevaluated_legacy_vector_exception(exc: Exception, lepton_input: Any) -> bool:
    if not isinstance(exc, TypeError):
        return False
    if (
        "must provide a dipole source and/or mu-e conversion scalar/dipole "
        "coefficient or rs_ew_couplings vector inputs"
        not in str(exc)
    ):
        return False
    return _has_finite_legacy_vector_like_input(lepton_input)


def _has_finite_legacy_vector_like_input(lepton_input: Any) -> bool:
    return any(
        _finite_complex_like(value)
        for value in _present_values(lepton_input, _LEGACY_VECTOR_COEFFICIENT_KEYS)
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


def _tree_only_lepton_input() -> dict[str, str]:
    return {
        "source": "tree-only rs_ew_couplings input; lepton_mass_basis_couplings absent"
    }
