"""C002 - CP violation in neutral charm mixing.

Physics
-------
``|q/p|`` and ``phi_D = Arg(q/p)`` constrain indirect CP violation in
``D0-D0bar`` mixing.  The repo has an audited short-distance Delta F = 2
D0-amplitude path, but no full HFLAV ``(|q/p|, phi_D)`` likelihood and no
grounded Standard-Model long-distance charm phase.  This constraint therefore
uses a conservative CP-odd amplitude proxy:

    |Im M12^NP| / (Delta m_D^exp / 2).

The Wilson coefficients are QCD-evolved to ``mu_had = 2 GeV`` before applying
the D0 matrix elements.  The HARD budget is set by the HFLAV all-CPV
``phiM = Arg(M12)`` 95% CL interval in the C002 HFLAV snapshot,
``[-1.48, 1.35]`` degrees, mapped to ``max|sin(phiM)|``.

Severity
--------
HARD.  The no-indirect-CPV point ``(|q/p|, phi) = (1, 0)`` is experimentally
allowed, so a new short-distance CP-odd D0 mixing component must fit inside
the current HFLAV near-zero CPV room.  The mapping is intentionally marked
``NEEDS-HUMAN-PHYSICS`` because a full treatment needs the HFLAV contour,
Gamma12, and the long-distance SM phase.

Catalog sidecar
---------------
``flavor_catalog/processes/charm/C002.yaml`` is the source of truth for the
HFLAV ``|q/p|``, ``phi_D``, no-indirect-CPV anchors, and local HFLAV snapshot
path used for ``phiM``.  The D0 amplitude normalization uses the same
``Delta_m_D.value_GeV / 2`` YAML anchor as C001.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
import re
from typing import Any, Mapping, Sequence

from flavor_catalog_constraints.anchors import (
    Anchor,
    AnchorError,
    load_full_yaml,
    load_anchor,
)
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.deltaf2 import (
    d0_mixing_from_wilsons_with_running,
    d0_mixing_m12_np_from_wilsons_with_running,
    d0_mixing_wilsons_from_couplings,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "charm"
_REQUIRED_EXTRA = "quark_mass_basis_couplings"

_Q_OVER_P_ANCHOR_CANDIDATES = ("q_over_p",)
_PHI_D_ANCHOR_CANDIDATES = ("phi_D",)
_NO_INDIRECT_CPV_CANDIDATES = ("no_indirect_cpv_test",)
_DELTA_M_D_PROCESS_ID = "C001"
_DELTA_M_D_ANCHOR_CANDIDATES = ("Delta_m_D",)
_CP_CONSERVING_Q_OVER_P = 1.0
_CP_CONSERVING_PHI_D_DEG = 0.0
_CP_CONSERVING_PHI_M_DEG = 0.0
_MU_HAD_GEV = 2.0
_REPO_ROOT = Path(__file__).resolve().parents[3]
_FLOAT_RE = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
_PLUS_MINUS_RE = re.compile(
    rf"(?P<value>{_FLOAT_RE})\s+\+(?P<plus>{_FLOAT_RE})\/-(?P<minus>{_FLOAT_RE})"
)
_CI_RE = re.compile(rf"\[(?P<low>{_FLOAT_RE}),\s*(?P<high>{_FLOAT_RE})\]")
_BUDGET_DOC_CITATION = (
    "flavor_catalog/processes/charm/C002.yaml:q_over_p,phi_D snapshot_path -> "
    "flavor_catalog/references/C002/hflav_ckm25_dmixing_cpv_fit.txt:phiM; "
    "flavor_catalog/processes/charm/C001.yaml:Delta_m_D; "
    "docs/quark_scan_assumptions_compact.tex:464-470"
)


@dataclass(frozen=True)
class AsymmetricCPVAnchor:
    """Typed HFLAV value with asymmetric one-sigma and 95% CL interval."""

    process_id: str
    block_key: str
    value: float
    uncertainty_plus: float
    uncertainty_minus: float
    confidence_interval_95: tuple[float, float]
    observable: str | None
    units: str | None
    source: str | None
    source_url: str | None
    year: int | None
    snapshot_path: str | None


@dataclass(frozen=True)
class NoIndirectCPVAnchor:
    """Typed view over the HFLAV no-indirect-CPV test block."""

    block_key: str
    delta_chi2: float
    significance_sigma: float
    interpretation: str | None
    source: str | None
    source_url: str | None
    year: int | None
    snapshot_path: str | None


@dataclass(frozen=True)
class CPVBudget:
    """Uncertainty-aware C002 CPV budget from the HFLAV ``phiM`` 95% room."""

    doc_citation: str
    q_over_p_central_deviation: float
    phi_d_central_deviation_sine: float
    phi_m_central_deviation_sine: float
    central_observed_deviation_proxy: float
    q_over_p_95_deviation_budget: float
    phi_d_95_deviation_degrees: float
    phi_d_95_deviation_sine_budget: float
    phi_m_95_deviation_degrees: float
    phi_m_95_deviation_sine_budget: float
    hard_veto_budget: float
    active_budget_source: str


@dataclass(frozen=True)
class D0CPVAnchor:
    """Typed C002 anchor: HFLAV CPV values plus D0 amplitude normalization."""

    q_over_p: AsymmetricCPVAnchor
    phi_d: AsymmetricCPVAnchor
    phi_m: AsymmetricCPVAnchor
    no_indirect_cpv: NoIndirectCPVAnchor
    delta_m_d_gev: Anchor
    budget_band: CPVBudget

    @property
    def value(self) -> float:
        """Central observed CPV deviation proxy; zero is CP conservation."""
        return self.budget_band.central_observed_deviation_proxy

    @property
    def sm_value(self) -> float:
        """CP-conserving SM proxy used when no LD charm phase is grounded."""
        return 0.0

    @property
    def budget(self) -> float:
        """Dimensionless HARD budget for the CP-odd D0 amplitude fraction."""
        return self.budget_band.hard_veto_budget

    @property
    def m12_budget_gev(self) -> float:
        """Conservative D0 amplitude room ``Delta m_D^exp / 2`` in GeV."""
        return 0.5 * self.delta_m_d_gev.value

    @property
    def cp_odd_m12_budget_gev(self) -> float:
        """Allowed CP-odd ``|Im M12^NP|`` room in GeV."""
        return self.m12_budget_gev * self.budget


def _optional_str(value: Any) -> str | None:
    return None if value is None else str(value)


def _optional_int(value: Any) -> int | None:
    if value is None:
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _required_float(value: Any, *, process_id: str, field_name: str) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError) as exc:
        raise AnchorError(
            f"{process_id}: C002 anchor field {field_name!r}={value!r} is not numeric"
        ) from exc
    if not math.isfinite(out):
        raise AnchorError(f"{process_id}: C002 anchor field {field_name!r} is not finite")
    return out


def _required_ci95(
    value: Any,
    *,
    process_id: str,
    field_name: str,
) -> tuple[float, float]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)) or len(value) != 2:
        raise AnchorError(
            f"{process_id}: {field_name!r} must be a two-entry 95% confidence interval"
        )
    low = _required_float(value[0], process_id=process_id, field_name=f"{field_name}[0]")
    high = _required_float(value[1], process_id=process_id, field_name=f"{field_name}[1]")
    if low >= high:
        raise AnchorError(f"{process_id}: {field_name!r} must satisfy low < high")
    return (low, high)


def _load_pdg_or_equivalent_block(process_id: str) -> Mapping[str, Any]:
    data = load_full_yaml(process_id, family=_FAMILY)
    block = data.get("pdg_or_equivalent")
    if not isinstance(block, Mapping):
        raise AnchorError(
            f"{process_id}: missing or invalid 'pdg_or_equivalent' mapping in sidecar"
        )
    return block


def _subblock_from_anchor(anchor: Anchor) -> Mapping[str, Any]:
    block = _load_pdg_or_equivalent_block(anchor.process_id)
    sub = block.get(anchor.block_key)
    if not isinstance(sub, Mapping):
        raise AnchorError(
            f"{anchor.process_id}: anchor block {anchor.block_key!r} is not a mapping"
        )
    return sub


def _load_asymmetric_cpv_anchor(
    process_id: str,
    candidates: Sequence[str],
) -> AsymmetricCPVAnchor:
    value_anchor = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=candidates,
        value_key="value",
        uncertainty_key="uncertainty_plus",
    )
    sub = _subblock_from_anchor(value_anchor)
    return AsymmetricCPVAnchor(
        process_id=process_id,
        block_key=value_anchor.block_key,
        value=value_anchor.value,
        uncertainty_plus=_required_float(
            sub.get("uncertainty_plus"),
            process_id=process_id,
            field_name="uncertainty_plus",
        ),
        uncertainty_minus=_required_float(
            sub.get("uncertainty_minus"),
            process_id=process_id,
            field_name="uncertainty_minus",
        ),
        confidence_interval_95=_required_ci95(
            sub.get("confidence_interval_95_percent"),
            process_id=process_id,
            field_name="confidence_interval_95_percent",
        ),
        observable=_optional_str(sub.get("observable")),
        units=_optional_str(sub.get("units")),
        source=_optional_str(sub.get("source")),
        source_url=_optional_str(sub.get("source_url")),
        year=_optional_int(sub.get("year")),
        snapshot_path=_optional_str(sub.get("snapshot_path")),
    )


def _load_no_indirect_cpv_anchor(process_id: str) -> NoIndirectCPVAnchor:
    pdg_block = _load_pdg_or_equivalent_block(process_id)
    block_key = next(
        (key for key in _NO_INDIRECT_CPV_CANDIDATES if key in pdg_block),
        None,
    )
    if block_key is None:
        raise AnchorError(
            f"{process_id}: none of the expected anchor keys "
            f"{list(_NO_INDIRECT_CPV_CANDIDATES)} found in pdg_or_equivalent"
        )
    sub = pdg_block[block_key]
    if not isinstance(sub, Mapping):
        raise AnchorError(f"{process_id}: anchor block {block_key!r} is not a mapping")
    return NoIndirectCPVAnchor(
        block_key=block_key,
        delta_chi2=_required_float(
            sub.get("delta_chi2"),
            process_id=process_id,
            field_name="delta_chi2",
        ),
        significance_sigma=_required_float(
            sub.get("significance_sigma"),
            process_id=process_id,
            field_name="significance_sigma",
        ),
        interpretation=_optional_str(sub.get("interpretation")),
        source=_optional_str(sub.get("source")),
        source_url=_optional_str(sub.get("source_url")),
        year=_optional_int(sub.get("year")),
        snapshot_path=_optional_str(sub.get("snapshot_path")),
    )


def _snapshot_path_from_anchor(process_id: str, anchor: AsymmetricCPVAnchor) -> Path:
    if anchor.snapshot_path is None:
        raise AnchorError(f"{process_id}: C002 HFLAV anchor has no snapshot_path")
    path = Path(anchor.snapshot_path)
    if not path.is_absolute():
        path = _REPO_ROOT / path
    if not path.is_file():
        raise AnchorError(f"{process_id}: C002 HFLAV snapshot not found: {path}")
    return path


def _load_phi_m_anchor_from_snapshot(
    process_id: str,
    reference_anchor: AsymmetricCPVAnchor,
) -> AsymmetricCPVAnchor:
    snapshot_path = _snapshot_path_from_anchor(process_id, reference_anchor)
    try:
        text = snapshot_path.read_text(encoding="utf-8")
    except OSError as exc:
        raise AnchorError(f"{process_id}: could not read C002 HFLAV snapshot") from exc

    for line in text.splitlines():
        stripped = line.strip()
        if not stripped.startswith("phiM") or "(degrees)" not in stripped:
            continue
        value_matches = list(_PLUS_MINUS_RE.finditer(stripped))
        ci_match = _CI_RE.search(stripped)
        if not value_matches or ci_match is None:
            raise AnchorError(
                f"{process_id}: malformed phiM row in C002 HFLAV snapshot"
            )
        all_cpv_match = value_matches[-1]
        return AsymmetricCPVAnchor(
            process_id=process_id,
            block_key="phiM_snapshot",
            value=_required_float(
                all_cpv_match.group("value"),
                process_id=process_id,
                field_name="phiM.value",
            ),
            uncertainty_plus=_required_float(
                all_cpv_match.group("plus"),
                process_id=process_id,
                field_name="phiM.uncertainty_plus",
            ),
            uncertainty_minus=_required_float(
                all_cpv_match.group("minus"),
                process_id=process_id,
                field_name="phiM.uncertainty_minus",
            ),
            confidence_interval_95=_required_ci95(
                (ci_match.group("low"), ci_match.group("high")),
                process_id=process_id,
                field_name="phiM.confidence_interval_95_percent",
            ),
            observable="phiM = Arg(M12) in D0-D0bar mixing",
            units="degrees",
            source=reference_anchor.source,
            source_url=reference_anchor.source_url,
            year=reference_anchor.year,
            snapshot_path=reference_anchor.snapshot_path,
        )

    raise AnchorError(f"{process_id}: phiM row missing from C002 HFLAV snapshot")


def _build_budget_band(
    *,
    process_id: str,
    q_over_p: AsymmetricCPVAnchor,
    phi_d: AsymmetricCPVAnchor,
    phi_m: AsymmetricCPVAnchor,
) -> CPVBudget:
    q_low, q_high = q_over_p.confidence_interval_95
    phi_low, phi_high = phi_d.confidence_interval_95
    phi_m_low, phi_m_high = phi_m.confidence_interval_95
    q_budget = max(
        abs(q_low - _CP_CONSERVING_Q_OVER_P),
        abs(q_high - _CP_CONSERVING_Q_OVER_P),
    )
    phi_budget_degrees = max(
        abs(phi_low - _CP_CONSERVING_PHI_D_DEG),
        abs(phi_high - _CP_CONSERVING_PHI_D_DEG),
    )
    phi_budget_sine = abs(math.sin(math.radians(phi_budget_degrees)))
    phi_m_budget_degrees = max(
        abs(phi_m_low - _CP_CONSERVING_PHI_M_DEG),
        abs(phi_m_high - _CP_CONSERVING_PHI_M_DEG),
    )
    phi_m_budget_sine = abs(math.sin(math.radians(phi_m_budget_degrees)))
    hard_budget = phi_m_budget_sine
    if hard_budget <= 0.0:
        raise AnchorError(f"{process_id}: C002 CPV budget must be positive")

    q_central_deviation = abs(q_over_p.value - _CP_CONSERVING_Q_OVER_P)
    phi_central_sine = abs(
        math.sin(math.radians(phi_d.value - _CP_CONSERVING_PHI_D_DEG))
    )
    phi_m_central_sine = abs(
        math.sin(math.radians(phi_m.value - _CP_CONSERVING_PHI_M_DEG))
    )
    return CPVBudget(
        doc_citation=_BUDGET_DOC_CITATION,
        q_over_p_central_deviation=q_central_deviation,
        phi_d_central_deviation_sine=phi_central_sine,
        phi_m_central_deviation_sine=phi_m_central_sine,
        central_observed_deviation_proxy=phi_m_central_sine,
        q_over_p_95_deviation_budget=q_budget,
        phi_d_95_deviation_degrees=phi_budget_degrees,
        phi_d_95_deviation_sine_budget=phi_budget_sine,
        phi_m_95_deviation_degrees=phi_m_budget_degrees,
        phi_m_95_deviation_sine_budget=phi_m_budget_sine,
        hard_veto_budget=hard_budget,
        active_budget_source="phi_m_95_sine_deviation",
    )


def _load_d0_cpv_anchor(process_id: str) -> D0CPVAnchor:
    q_over_p = _load_asymmetric_cpv_anchor(process_id, _Q_OVER_P_ANCHOR_CANDIDATES)
    phi_d = _load_asymmetric_cpv_anchor(process_id, _PHI_D_ANCHOR_CANDIDATES)
    phi_m = _load_phi_m_anchor_from_snapshot(process_id, phi_d)
    delta_m_d_gev = load_anchor(
        _DELTA_M_D_PROCESS_ID,
        family=_FAMILY,
        candidates=_DELTA_M_D_ANCHOR_CANDIDATES,
        value_key="value_GeV",
        uncertainty_key="uncertainty_GeV",
    )
    anchor = D0CPVAnchor(
        q_over_p=q_over_p,
        phi_d=phi_d,
        phi_m=phi_m,
        no_indirect_cpv=_load_no_indirect_cpv_anchor(process_id),
        delta_m_d_gev=delta_m_d_gev,
        budget_band=_build_budget_band(
            process_id=process_id,
            q_over_p=q_over_p,
            phi_d=phi_d,
            phi_m=phi_m,
        ),
    )
    if anchor.m12_budget_gev <= 0.0:
        raise AnchorError(
            f"{process_id}: Delta m_D normalization budget must be positive, "
            f"got {anchor.m12_budget_gev}"
        )
    return anchor


def _complex_mapping(values: Mapping[str, complex]) -> dict[str, complex]:
    return {str(key): complex(value) for key, value in values.items()}


@register
class Constraint:
    """Catalogued D0 indirect-CPV Delta F=2 constraint (process_id C002)."""

    process_id = "C002"
    severity = Severity.HARD
    observable = "D0-D0bar indirect CP violation"

    def __init__(self) -> None:
        self.anchor = _load_d0_cpv_anchor(self.process_id)
        self.sm_value = self.anchor.sm_value

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        couplings = point.get_extra(_REQUIRED_EXTRA)
        if couplings is None:
            return ConstraintResult(
                process_id=self.process_id,
                severity=self.severity,
                passes=True,
                sm_prediction=float(self.anchor.sm_value),
                experimental=float(self.anchor.value),
                budget=float(self.anchor.budget),
                notes=(
                    f"extra {_REQUIRED_EXTRA!r} absent; D0 indirect-CPV "
                    "constraint was not evaluated."
                ),
                diagnostics={"missing_extra": _REQUIRED_EXTRA},
            )

        wilsons = d0_mixing_wilsons_from_couplings(couplings)
        mixing = d0_mixing_from_wilsons_with_running(
            wilsons,
            mu_had=_MU_HAD_GEV,
            m12_np_budget=self.anchor.m12_budget_gev,
        )
        m12_np = d0_mixing_m12_np_from_wilsons_with_running(
            wilsons,
            mu_had=_MU_HAD_GEV,
        )

        cp_odd_fraction = abs(float(m12_np.imag)) / self.anchor.m12_budget_gev
        ratio = cp_odd_fraction / self.anchor.budget
        m12_raw_argument_degrees = math.degrees(math.atan2(m12_np.imag, m12_np.real))
        abs_m12_np = abs(m12_np)
        amplitude_ratio = abs_m12_np / self.anchor.m12_budget_gev
        phase_sine = 0.0 if abs_m12_np == 0.0 else abs(m12_np.imag) / abs_m12_np

        return ConstraintResult(
            process_id=self.process_id,
            severity=self.severity,
            passes=ratio <= 1.0,
            predicted=float(cp_odd_fraction),
            sm_prediction=float(self.anchor.sm_value),
            experimental=float(self.anchor.value),
            ratio=float(ratio),
            budget=float(self.anchor.budget),
            notes=(
                "|Im(M12^NP)|/(Delta m_D^exp/2) uses QCD-evolved D0 Wilsons "
                "at 2 GeV.  HARD budget uses the C002 HFLAV all-CPV 95% "
                "phiM=Arg(M12) room from the local snapshot; full charm LD "
                "phase and HFLAV contour mapping are marked NEEDS-HUMAN-PHYSICS."
            ),
            diagnostics={
                "m12_np_gev": complex(m12_np),
                "re_m12_np_gev": float(m12_np.real),
                "im_m12_np_gev": float(m12_np.imag),
                "abs_m12_np_gev": float(abs_m12_np),
                "abs_m12_np_from_magnitude_evaluator_gev": float(mixing.abs_m12_np),
                "m12_abs_consistency_delta_gev": float(abs(abs_m12_np - mixing.abs_m12_np)),
                "m12_np_phase_degrees": float(m12_raw_argument_degrees),
                "m12_np_raw_argument_degrees": float(m12_raw_argument_degrees),
                "m12_np_raw_argument_note": (
                    "Raw Arg(M12^NP), not a q/p or phi_D prediction; pass/fail "
                    "uses only |Im(M12^NP)|."
                ),
                "m12_np_phase_sine_abs": float(phase_sine),
                "amplitude_ratio_to_delta_m_d_half": float(amplitude_ratio),
                "cp_odd_fraction": float(cp_odd_fraction),
                "q_over_p_deviation_proxy": float(cp_odd_fraction),
                "phi_d_sine_deviation_proxy": float(cp_odd_fraction),
                "qcd_running_applied": True,
                "hadronic_scale_gev": _MU_HAD_GEV,
                "matching_scale_gev": float(wilsons.matching_scale),
                "m_kk_gev": float(wilsons.M_KK),
                "left_uc_coupling": complex(wilsons.left_coupling),
                "right_uc_coupling": complex(wilsons.right_coupling),
                "wilson_coefficients": _complex_mapping(wilsons.wilsons),
                "m12_normalization_budget_gev": float(self.anchor.m12_budget_gev),
                "cp_odd_m12_budget_gev": float(self.anchor.cp_odd_m12_budget_gev),
                "delta_m_d_experimental_gev": float(self.anchor.delta_m_d_gev.value),
                "delta_m_d_uncertainty_gev": (
                    None
                    if self.anchor.delta_m_d_gev.uncertainty is None
                    else float(self.anchor.delta_m_d_gev.uncertainty)
                ),
                "q_over_p_value": float(self.anchor.q_over_p.value),
                "q_over_p_uncertainty_plus": float(self.anchor.q_over_p.uncertainty_plus),
                "q_over_p_uncertainty_minus": float(
                    self.anchor.q_over_p.uncertainty_minus
                ),
                "q_over_p_95_low": float(self.anchor.q_over_p.confidence_interval_95[0]),
                "q_over_p_95_high": float(self.anchor.q_over_p.confidence_interval_95[1]),
                "q_over_p_95_deviation_budget": float(
                    self.anchor.budget_band.q_over_p_95_deviation_budget
                ),
                "phi_d_value_degrees": float(self.anchor.phi_d.value),
                "phi_d_uncertainty_plus_degrees": float(
                    self.anchor.phi_d.uncertainty_plus
                ),
                "phi_d_uncertainty_minus_degrees": float(
                    self.anchor.phi_d.uncertainty_minus
                ),
                "phi_d_95_low_degrees": float(
                    self.anchor.phi_d.confidence_interval_95[0]
                ),
                "phi_d_95_high_degrees": float(
                    self.anchor.phi_d.confidence_interval_95[1]
                ),
                "phi_d_95_deviation_degrees": float(
                    self.anchor.budget_band.phi_d_95_deviation_degrees
                ),
                "phi_d_95_sine_budget": float(
                    self.anchor.budget_band.phi_d_95_deviation_sine_budget
                ),
                "phi_m_value_degrees": float(self.anchor.phi_m.value),
                "phi_m_uncertainty_plus_degrees": float(
                    self.anchor.phi_m.uncertainty_plus
                ),
                "phi_m_uncertainty_minus_degrees": float(
                    self.anchor.phi_m.uncertainty_minus
                ),
                "phi_m_95_low_degrees": float(
                    self.anchor.phi_m.confidence_interval_95[0]
                ),
                "phi_m_95_high_degrees": float(
                    self.anchor.phi_m.confidence_interval_95[1]
                ),
                "phi_m_95_deviation_degrees": float(
                    self.anchor.budget_band.phi_m_95_deviation_degrees
                ),
                "phi_m_95_sine_budget": float(
                    self.anchor.budget_band.phi_m_95_deviation_sine_budget
                ),
                "hard_veto_budget_source": self.anchor.budget_band.active_budget_source,
                "budget_doc_citation": self.anchor.budget_band.doc_citation,
                "central_observed_deviation_proxy": float(self.anchor.value),
                "q_over_p_central_deviation": float(
                    self.anchor.budget_band.q_over_p_central_deviation
                ),
                "phi_d_central_deviation_sine": float(
                    self.anchor.budget_band.phi_d_central_deviation_sine
                ),
                "phi_m_central_deviation_sine": float(
                    self.anchor.budget_band.phi_m_central_deviation_sine
                ),
                "cp_conserving_q_over_p": _CP_CONSERVING_Q_OVER_P,
                "cp_conserving_phi_d_degrees": _CP_CONSERVING_PHI_D_DEG,
                "cp_conserving_phi_m_degrees": _CP_CONSERVING_PHI_M_DEG,
                "no_indirect_cpv_delta_chi2": float(
                    self.anchor.no_indirect_cpv.delta_chi2
                ),
                "no_indirect_cpv_significance_sigma": float(
                    self.anchor.no_indirect_cpv.significance_sigma
                ),
                "sm_long_distance_phase_grounded": False,
                "full_hflav_likelihood_used": False,
                "needs_human_physics": (
                    "NEEDS-HUMAN-PHYSICS: full |q/p|-phi_D contour/covariance, "
                    "Gamma12, and the SM long-distance charm phase are not "
                    "available in the current adapter; v1 constrains the "
                    "QCD-evolved CP-odd M12^NP fraction with one-dimensional "
                    "HFLAV phiM 95% room."
                ),
                "experimental_blocks": {
                    "q_over_p": self.anchor.q_over_p.block_key,
                    "phi_d": self.anchor.phi_d.block_key,
                    "phi_m": self.anchor.phi_m.block_key,
                    "no_indirect_cpv": self.anchor.no_indirect_cpv.block_key,
                    "delta_m_d": self.anchor.delta_m_d_gev.block_key,
                },
            },
        )
