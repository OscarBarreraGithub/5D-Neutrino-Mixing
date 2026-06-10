"""EW001 - Peskin-Takeuchi oblique electroweak parameters.

Physics
-------
EW001 compares the RS oblique-parameter prediction to the PDG 2025 global
``(S,T)`` fit with ``U`` fixed to zero.  The rigorous part is the correlated
Gaussian ellipse

    chi2 = (x - x0)^T Cov(S,T)^-1 (x - x0),

where ``x = (S,T)`` and the central values, uncertainties, and correlation are
loaded from ``EW001.yaml``.  The SM reference is ``S = T = U = 0`` by
construction.

RS matching status
------------------
NEEDS-HUMAN-PHYSICS.  The current ``ParameterPoint`` does not carry the full
electroweak KK gauge sector.  The NP side therefore uses the documented
minimal-RS proxy from ``quarkConstraints.oblique_stu``:
``Delta S = c_S v^2/M_KK^2`` with ``c_S`` loaded from the PDG warped-context
anchor, and ``Delta T = pi L/(2 c_W^2) v^2/M_KK^2`` with ``U = 0``.  This is
the classic RS T-parameter problem proxy, not a full custodial/EW fit.

Severity
--------
HARD.  Oblique parameters are observed electroweak-precision bounds on new
physics, and the RS proxy point must lie inside the chosen 95% ``(S,T)`` fit
contour.  The scalar ``predicted`` field is the ST chi2; vector components are
reported in diagnostics.

Catalog sidecar
---------------
``flavor_catalog/processes/top_higgs_ew/EW001.yaml`` is the source of truth
for the PDG central values, uncertainties, correlation, and warped-S
coefficient.  Its ``pdg_or_equivalent`` block is list-shaped, so this module
routes selected list entries through the scaffold ``load_anchor`` helper and
fails loudly on missing observables.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Mapping, Sequence

from flavor_catalog_constraints import anchors as anchor_scaffold
from flavor_catalog_constraints.anchors import Anchor, AnchorError, load_anchor, load_full_yaml
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.oblique_stu import (
    CHI2_2DOF_95,
    DEFAULT_HIGGS_VEV_GEV,
    DEFAULT_RS_VOLUME_LOG,
    DEFAULT_SIN2_THETA_W,
    OBLIQUE_STU_LIKELIHOOD_V1,
    OBLIQUE_STU_RS_PROXY_V1,
    ObliqueSTFit,
    compare_oblique_st_to_fit,
    evaluate_rs_oblique_proxy,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "top_higgs_ew"
_EXPECTED_UNITS = "dimensionless"
_MASS_EXTRAS = ("kk_ew_mass_gev", "kk_gluon_mass_gev")
_COUPLINGS_EXTRA = "quark_mass_basis_couplings"
_RS_EW_COUPLINGS_EXTRA = "rs_ew_couplings"
_S_U_FIXED = "PDG 2025 global oblique fit, S with U fixed to zero"
_T_U_FIXED = "PDG 2025 global oblique fit, T with U fixed to zero"
_U_FIXED = "PDG 2025 global oblique fit, U fixed value"
_RHO_ST_U_FIXED = "PDG 2025 correlation rho(S,T), U fixed to zero"
_WARPED_S_COEFF = "PDG 2025 warped-extra-dimension S coefficient"
_WARPED_S_UPPER = "PDG 2025 one-sided S upper value used for warped M_KK context"
_WARPED_MKK_CONTEXT = "PDG 2025 warped M_KK lower-bound context"
_FIT_CONTOUR_SOURCE = (
    "PDG 2025 EW001.yaml U-fixed S/T central values, uncertainties, "
    "rho(S,T), and standard chi2_2 95% contour"
)


@dataclass(frozen=True)
class EW001ValueAnchor:
    """Typed view over one EW001 list entry selected via ``load_anchor``."""

    anchor: Anchor
    observable: str
    entry_index: int
    cl: str | None
    access_date: str | None
    sha256: str | None

    @property
    def value(self) -> float:
        return self.anchor.value

    @property
    def uncertainty(self) -> float | None:
        return self.anchor.uncertainty

    @property
    def units(self) -> str | None:
        return self.anchor.units

    @property
    def source(self) -> str | None:
        return self.anchor.source

    @property
    def source_url(self) -> str | None:
        return self.anchor.source_url

    @property
    def snapshot_path(self) -> str | None:
        return self.anchor.snapshot_path

    @property
    def block_key(self) -> str:
        return self.anchor.block_key


@dataclass(frozen=True)
class EW001Anchor:
    """YAML-loaded ST fit and warped-context anchors for EW001."""

    s_u_fixed: EW001ValueAnchor
    t_u_fixed: EW001ValueAnchor
    u_fixed: EW001ValueAnchor
    rho_st_u_fixed: EW001ValueAnchor
    warped_s_coefficient: EW001ValueAnchor
    warped_s_upper_95: EW001ValueAnchor
    warped_mkk_context: EW001ValueAnchor
    fit: ObliqueSTFit

    @property
    def value(self) -> float:
        return self.fit.chi2_budget

    @property
    def budget(self) -> float:
        return self.fit.chi2_budget


def _required_float(value: Any, *, process_id: str, field_name: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise AnchorError(
            f"{process_id}: EW001 field {field_name!r}={value!r} is not numeric"
        ) from exc
    if not math.isfinite(number):
        raise AnchorError(
            f"{process_id}: EW001 field {field_name!r}={value!r} is not finite"
        )
    return number


def _positive_float(value: Any, *, process_id: str, field_name: str) -> float:
    number = _required_float(value, process_id=process_id, field_name=field_name)
    if number <= 0.0:
        raise AnchorError(f"{process_id}: EW001 field {field_name!r} must be positive")
    return number


def _optional_str(value: Any) -> str | None:
    return None if value is None else str(value)


def _pdg_entries(process_id: str) -> Sequence[Mapping[str, Any]]:
    data = load_full_yaml(process_id, family=_FAMILY)
    entries = data.get("pdg_or_equivalent")
    if not isinstance(entries, list) or not entries:
        raise AnchorError(f"{process_id}: expected non-empty pdg_or_equivalent list")
    for index, entry in enumerate(entries):
        if not isinstance(entry, Mapping):
            raise AnchorError(
                f"{process_id}: pdg_or_equivalent[{index}] is not a mapping"
            )
    return entries


def _find_entry(
    process_id: str,
    observable: str,
) -> tuple[int, Mapping[str, Any]]:
    for index, entry in enumerate(_pdg_entries(process_id)):
        if entry.get("observable") == observable:
            return index, entry
    present = [str(entry.get("observable")) for entry in _pdg_entries(process_id)]
    raise AnchorError(
        f"{process_id}: observable {observable!r} not found in "
        f"pdg_or_equivalent list (present: {present})"
    )


def _load_scaffold_list_anchor(
    observable: str,
    *,
    process_id: str,
    uncertainty_key: str = "uncertainty",
) -> tuple[Anchor, Mapping[str, Any], int]:
    index, entry = _find_entry(process_id, observable)
    block_key = f"pdg_or_equivalent[{index}]"
    virtual_block = {block_key: dict(entry)}
    original_load_pdg_block = anchor_scaffold.load_pdg_block

    def _load_virtual_pdg_block(
        request_process_id: str,
        **kwargs: Any,
    ) -> Mapping[str, Any]:
        if request_process_id == process_id and kwargs.get("family") == _FAMILY:
            return virtual_block
        return original_load_pdg_block(request_process_id, **kwargs)

    anchor_scaffold.load_pdg_block = _load_virtual_pdg_block
    try:
        scaffold_anchor = load_anchor(
            process_id,
            family=_FAMILY,
            candidates=(block_key,),
            uncertainty_key=uncertainty_key,
        )
    finally:
        anchor_scaffold.load_pdg_block = original_load_pdg_block

    if scaffold_anchor.block_key != block_key:
        raise AnchorError(
            f"{process_id}: load_anchor selected {scaffold_anchor.block_key!r}, "
            f"expected {block_key!r} for EW001 observable {observable!r}"
        )
    return scaffold_anchor, entry, index


def _load_value_anchor(
    observable: str,
    *,
    process_id: str,
    expected_units: str = _EXPECTED_UNITS,
    uncertainty_key: str = "uncertainty",
) -> EW001ValueAnchor:
    scaffold_anchor, entry, index = _load_scaffold_list_anchor(
        observable,
        process_id=process_id,
        uncertainty_key=uncertainty_key,
    )
    if scaffold_anchor.units != expected_units:
        raise AnchorError(
            f"{process_id}: expected units {expected_units!r} for {observable}, "
            f"got {scaffold_anchor.units!r}"
        )
    return EW001ValueAnchor(
        anchor=scaffold_anchor,
        observable=observable,
        entry_index=index,
        cl=_optional_str(entry.get("cl")),
        access_date=_optional_str(entry.get("access_date")),
        sha256=_optional_str(entry.get("sha256")),
    )


def _build_st_fit(
    *,
    process_id: str,
    s_anchor: EW001ValueAnchor,
    t_anchor: EW001ValueAnchor,
    u_anchor: EW001ValueAnchor,
    rho_anchor: EW001ValueAnchor,
) -> ObliqueSTFit:
    if s_anchor.uncertainty is None or t_anchor.uncertainty is None:
        raise AnchorError(f"{process_id}: S and T anchors must carry uncertainties")
    sigma_s = _positive_float(
        s_anchor.uncertainty,
        process_id=process_id,
        field_name=f"{s_anchor.observable}.uncertainty",
    )
    sigma_t = _positive_float(
        t_anchor.uncertainty,
        process_id=process_id,
        field_name=f"{t_anchor.observable}.uncertainty",
    )
    return ObliqueSTFit(
        s_central=float(s_anchor.value),
        t_central=float(t_anchor.value),
        u_fixed=float(u_anchor.value),
        sigma_s=float(sigma_s),
        sigma_t=float(sigma_t),
        rho_st=float(rho_anchor.value),
        chi2_budget=CHI2_2DOF_95,
        confidence_level=0.95,
        source=s_anchor.source,
        source_url=s_anchor.source_url,
        covariance_source=_FIT_CONTOUR_SOURCE,
    )


def _load_ew001_anchor(process_id: str) -> EW001Anchor:
    s_anchor = _load_value_anchor(_S_U_FIXED, process_id=process_id)
    t_anchor = _load_value_anchor(_T_U_FIXED, process_id=process_id)
    u_anchor = _load_value_anchor(_U_FIXED, process_id=process_id)
    rho_anchor = _load_value_anchor(
        _RHO_ST_U_FIXED,
        process_id=process_id,
        uncertainty_key="__ew001_no_scalar_uncertainty__",
    )
    warped_s_coefficient = _load_value_anchor(
        _WARPED_S_COEFF,
        process_id=process_id,
        expected_units="coefficient in S approximately coefficient * v^2 / M_KK^2",
        uncertainty_key="__ew001_no_scalar_uncertainty__",
    )
    warped_s_upper_95 = _load_value_anchor(
        _WARPED_S_UPPER,
        process_id=process_id,
        uncertainty_key="__ew001_no_scalar_uncertainty__",
    )
    warped_mkk_context = _load_value_anchor(
        _WARPED_MKK_CONTEXT,
        process_id=process_id,
        expected_units="TeV",
        uncertainty_key="__ew001_no_scalar_uncertainty__",
    )
    return EW001Anchor(
        s_u_fixed=s_anchor,
        t_u_fixed=t_anchor,
        u_fixed=u_anchor,
        rho_st_u_fixed=rho_anchor,
        warped_s_coefficient=warped_s_coefficient,
        warped_s_upper_95=warped_s_upper_95,
        warped_mkk_context=warped_mkk_context,
        fit=_build_st_fit(
            process_id=process_id,
            s_anchor=s_anchor,
            t_anchor=t_anchor,
            u_anchor=u_anchor,
            rho_anchor=rho_anchor,
        ),
    )


def _resolve_m_kk_gev(point: ParameterPoint) -> tuple[float | None, str | None]:
    for key in _MASS_EXTRAS:
        value = point.get_extra(key)
        if value is not None:
            return float(value), key
    couplings = point.get_extra(_COUPLINGS_EXTRA)
    if couplings is not None and getattr(couplings, "M_KK", None) is not None:
        return float(getattr(couplings, "M_KK")), f"{_COUPLINGS_EXTRA}.M_KK"
    return None, None


def _resolve_ew_model(point: ParameterPoint) -> str:
    rs_ew_couplings = point.get_extra(_RS_EW_COUPLINGS_EXTRA)
    metadata = getattr(rs_ew_couplings, "metadata", None)
    if isinstance(metadata, Mapping) and metadata.get("ew_model") is not None:
        return str(metadata["ew_model"])

    raw = point.raw
    if isinstance(raw, Mapping) and raw.get("ew_model") is not None:
        return str(raw["ew_model"])
    raw_model = getattr(raw, "ew_model", None)
    if raw_model is not None:
        return str(raw_model)

    point_model = getattr(point, "ew_model", None)
    if point_model is not None:
        return str(point_model)
    return "minimal_rs"


def _resolve_rs_ew_metadata(point: ParameterPoint) -> Mapping[str, Any]:
    rs_ew_couplings = point.get_extra(_RS_EW_COUPLINGS_EXTRA)
    metadata = getattr(rs_ew_couplings, "metadata", None)
    if isinstance(metadata, Mapping):
        return metadata
    return {}


def _sm_reference_chi2(fit: ObliqueSTFit) -> tuple[float, float, float, bool]:
    return compare_oblique_st_to_fit(s=0.0, t=0.0, fit=fit)


@register
class Constraint:
    """Catalogued oblique ``S,T,U`` electroweak-precision constraint."""

    process_id = "EW001"
    severity = Severity.HARD
    observable = "S,T,U oblique electroweak parameters"

    def __init__(self) -> None:
        self.anchor = _load_ew001_anchor(self.process_id)
        self.sm_reference = _sm_reference_chi2(self.anchor.fit)

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        sm_chi2, _, _, _ = self.sm_reference
        ew_model = _resolve_ew_model(point)
        rs_ew_metadata = _resolve_rs_ew_metadata(point)
        m_kk_gev, mass_source = _resolve_m_kk_gev(point)
        if m_kk_gev is None:
            return ConstraintResult(
                process_id=self.process_id,
                severity=self.severity,
                passes=True,
                sm_prediction=float(sm_chi2),
                experimental=0.0,
                budget=float(self.anchor.budget),
                notes=(
                    "No KK electroweak mass was available; EW001 ST ellipse "
                    "constraint was not evaluated."
                ),
                diagnostics={
                    "missing_extras": _MASS_EXTRAS + (_COUPLINGS_EXTRA,),
                    "ew_model": ew_model,
                    "needs_human_physics": OBLIQUE_STU_RS_PROXY_V1,
                    "sm_reference_s": 0.0,
                    "sm_reference_t": 0.0,
                    "sm_reference_u": 0.0,
                    "sm_reference_chi2": float(sm_chi2),
                    "budget_source": _FIT_CONTOUR_SOURCE,
                },
            )

        loop_metadata: Mapping[str, Any] | None = None
        delta_t_loop: Any = 0.0
        if (
            ew_model == "custodial_rs_plr"
            and bool(rs_ew_metadata.get("top_partner_loop_numerics_included", False))
        ):
            loop_metadata = dict(rs_ew_metadata)
            delta_t_loop = rs_ew_metadata.get("top_partner_delta_t_loop_applied")

        try:
            comparison = evaluate_rs_oblique_proxy(
                m_kk_gev=float(m_kk_gev),
                fit=self.anchor.fit,
                s_coefficient=float(self.anchor.warped_s_coefficient.value),
                ew_model=ew_model,
                delta_t_loop=delta_t_loop,
                loop_metadata=loop_metadata,
            )
        except ValueError as exc:
            return ConstraintResult(
                process_id=self.process_id,
                severity=self.severity,
                passes=False,
                sm_prediction=float(sm_chi2),
                experimental=0.0,
                budget=float(self.anchor.budget),
                notes=f"Invalid EW001 mass/proxy input: {exc}",
                diagnostics={
                    "invalid_m_kk_gev": m_kk_gev,
                    "mass_source": mass_source,
                    "ew_model": ew_model,
                    "needs_human_physics": OBLIQUE_STU_RS_PROXY_V1,
                },
            )

        prediction = comparison.prediction
        diagnostics = dict(comparison.diagnostics)
        diagnostics.update(
            {
                "needs_human_physics": OBLIQUE_STU_RS_PROXY_V1,
                "likelihood_model": OBLIQUE_STU_LIKELIHOOD_V1,
                "ew_model": ew_model,
                "m_kk_gev": float(prediction.m_kk_gev),
                "mass_source": mass_source,
                "s_prediction": float(prediction.s),
                "t_prediction": float(prediction.t),
                "u_prediction": float(prediction.u),
                "sm_reference_s": 0.0,
                "sm_reference_t": 0.0,
                "sm_reference_u": 0.0,
                "sm_reference_chi2": float(sm_chi2),
                "delta_s_from_fit_center": float(comparison.delta_s_from_fit_center),
                "delta_t_from_fit_center": float(comparison.delta_t_from_fit_center),
                "s_coefficient": float(prediction.s_coefficient),
                "t_coefficient": float(prediction.t_coefficient),
                "higgs_vev_gev": float(DEFAULT_HIGGS_VEV_GEV),
                "sin2_theta_w": float(DEFAULT_SIN2_THETA_W),
                "rs_volume_log": float(DEFAULT_RS_VOLUME_LOG),
                "fit_contour_chi2": float(comparison.chi2_budget),
                "fit_contour_confidence_level": float(self.anchor.fit.confidence_level),
                "budget_source": _FIT_CONTOUR_SOURCE,
                "experimental_s_central": float(self.anchor.fit.s_central),
                "experimental_t_central": float(self.anchor.fit.t_central),
                "experimental_u_fixed": float(self.anchor.fit.u_fixed),
                "pdg_warped_s_upper_95": float(self.anchor.warped_s_upper_95.value),
                "pdg_warped_mkk_context_tev": float(
                    self.anchor.warped_mkk_context.value
                ),
                "pdg_warped_s_coefficient_source": (
                    self.anchor.warped_s_coefficient.source
                ),
                "pdg_s_block": self.anchor.s_u_fixed.block_key,
                "pdg_t_block": self.anchor.t_u_fixed.block_key,
                "pdg_rho_block": self.anchor.rho_st_u_fixed.block_key,
            }
        )

        return ConstraintResult(
            process_id=self.process_id,
            severity=self.severity,
            passes=bool(comparison.passes),
            predicted=float(comparison.chi2),
            sm_prediction=float(sm_chi2),
            experimental=0.0,
            ratio=float(comparison.ratio_to_budget),
            budget=float(comparison.chi2_budget),
            notes=(
                "EW001 uses the PDG U-fixed correlated (S,T) ellipse. "
                "The RS point is a documented minimal-RS v^2/M_KK^2 proxy "
                "for Delta S and the volume-enhanced Delta T term; full "
                "EW KK matching is flagged NEEDS-HUMAN-PHYSICS."
            ),
            diagnostics=diagnostics,
        )
