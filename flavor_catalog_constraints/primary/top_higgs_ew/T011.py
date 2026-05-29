"""T011 - ``Z -> b bbar`` pole asymmetries ``A_b`` and ``A_FB^{0,b}``.

Physics
-------
This constraint reuses the shared Z-pole pseudo-observable machinery from
``quarkConstraints.zpole`` through the
``flavor_catalog_constraints.physics_adapters.zpole`` boundary.  Effective
chiral couplings are normalized as

    L_Z = g_Z Z_mu bbar gamma^mu (g_L P_L + g_R P_R) b,

with

    A_b = (|g_L|^2 - |g_R|^2) / (|g_L|^2 + |g_R|^2),
    A_FB^{0,b} = 3/4 A_e A_b.

The SM-limit values are evaluated directly from the effective couplings; the
LEP/SLC SM-fit values and uncertainty budget are loaded from the catalog
sidecar/snapshot provenance.

RS matching status
------------------
NEEDS-HUMAN-PHYSICS.  A rigorous RS prediction for ``Z b_L b_L`` and
``Z b_R b_R`` shifts requires the electroweak KK/Z/Z' spectrum, custodial
representations, fermion embeddings, brane kinetic terms, and Z-mixing data.
Those inputs are not present on ``ParameterPoint``.  The NP part is therefore
the documented T010 proxy: available bottom-vs-light overlap
non-universality is mapped onto ``delta g_b ~ (m_Z/M_KK)^2 Delta overlap``.

Severity
--------
HARD.  The asymmetry measurements constrain point-dependent shifts in the
bottom chiral Z couplings.  Because ``A_FB^{0,b}`` is also a long-standing
SM tension, the veto scalar is the absolute NP shift relative to the SM-limit
pseudo-observable divided by the YAML/snapshot LEP-SLC-vs-SM loose-edge
budget, ``|exp - SM_fit| + sqrt(sigma_exp^2 + sigma_SM^2)``.  Thus the SM
point is not vetoed merely because of the legacy anomaly, while a large RS
shift is excluded.

Catalog sidecar
---------------
``flavor_catalog/processes/top_higgs_ew/T011.yaml`` is the source of truth for
T011.  Historical ``canonical_home``/T010 aliases are rejected: the asymmetry
anchors must live under T011 and load through the scaffold ``load_anchor`` path.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
import re
from typing import Any, Mapping, Sequence

from flavor_catalog_constraints import anchors as anchor_scaffold
from flavor_catalog_constraints.anchors import Anchor, AnchorError, load_anchor, load_full_yaml
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.zpole import (
    ZPOLE_RS_ZBB_PROXY_V1,
    ZPoleQuarkObservables,
    zpole_default_sm_inputs,
    zpole_evaluate_quark,
    zpole_evaluate_zbb_with_proxy,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "top_higgs_ew"
_REQUIRED_EXTRA = "quark_mass_basis_couplings"
_OPTIONAL_EW_MASS_EXTRA = "kk_ew_mass_gev"
_REPO_ROOT = Path(__file__).resolve().parents[3]

_AFB_OBSERVABLE = "A_FB^{0,b}"
_AB_OBSERVABLE = "A_b"
_AFB_PULL_OBSERVABLE = "LEP/SLC final-combination pull for A_FB^{0,b}"
_VETO_OBSERVABLES = (_AFB_OBSERVABLE, _AB_OBSERVABLE)
_SM_SNAPSHOT_LABELS = {
    _AFB_OBSERVABLE: "A_FB^(0,b)",
    _AB_OBSERVABLE: "A_b",
}
_NUMBER_RE = r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?"
_NEEDS_HUMAN_PHYSICS = (
    "NEEDS-HUMAN-PHYSICS: rigorous RS Zbb asymmetry matching requires EW "
    "KK/Z/Z' coupling shifts, custodial representations, fermion embeddings, "
    "brane kinetic terms, and Z-mixing data not present on ParameterPoint; "
    "this v1 uses the documented bottom-overlap non-universality proxy."
)


@dataclass(frozen=True)
class SMFitValue:
    """LEP/SLC SM-fit reference for one bottom asymmetry."""

    observable: str
    value: float
    uncertainty: float
    pull: float | None
    source_path: str | None
    source_url: str | None


@dataclass(frozen=True)
class AsymmetryBudget:
    """Loose-edge NP-shift budget for one bottom asymmetry."""

    observable: str
    central_residual: float
    experimental_sigma: float
    sm_fit_sigma: float
    combined_sigma: float
    hard_veto_budget: float
    source: str


@dataclass(frozen=True)
class T011Anchor:
    """Typed T011 anchor bundle."""

    a_fb: Anchor
    a_b: Anchor
    a_fb_pull: Anchor | None
    sm_fit_values: Mapping[str, SMFitValue]
    budgets: Mapping[str, AsymmetryBudget]
    source_process_id: str
    canonical_home_fallback: bool

    @property
    def value(self) -> float:
        return self.a_fb.value

    @property
    def budget(self) -> float:
        return self.budgets[_AFB_OBSERVABLE].hard_veto_budget


@dataclass(frozen=True)
class AsymmetryShift:
    """One scalar NP-shift comparison entering the T011 max-ratio result."""

    observable: str
    predicted: float
    sm_prediction: float
    experimental: float
    budget: float
    np_shift: float
    ratio: float
    total_exp_pull: float
    total_exp_pull_ratio: float


def _required_float(value: Any, *, process_id: str, field_name: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise AnchorError(
            f"{process_id}: T011 field {field_name!r}={value!r} is not numeric"
        ) from exc
    if not math.isfinite(number):
        raise AnchorError(f"{process_id}: T011 field {field_name!r} is not finite")
    return number


def _positive_float(value: Any, *, process_id: str, field_name: str) -> float:
    number = _required_float(value, process_id=process_id, field_name=field_name)
    if number <= 0.0:
        raise AnchorError(f"{process_id}: T011 field {field_name!r} must be positive")
    return number


def _standalone_sidecar(process_id: str) -> Mapping[str, Any]:
    data = load_full_yaml(process_id, family=_FAMILY)
    legacy_markers = [
        key for key in ("canonical_home", "merged_into") if data.get(key) is not None
    ]
    if legacy_markers:
        raise AnchorError(
            f"{process_id}: standalone T011 anchors are required; legacy "
            f"T010 fallback marker(s) present: {legacy_markers}"
        )
    if "pdg_or_equivalent" in data:
        return data
    raise AnchorError(
        f"{process_id}: missing standalone pdg_or_equivalent anchors for "
        f"{_AFB_OBSERVABLE} and {_AB_OBSERVABLE}; T010 fallback is not allowed"
    )


def _pdg_container(process_id: str) -> Mapping[str, Any] | Sequence[Mapping[str, Any]]:
    data = _standalone_sidecar(process_id)
    entries = data.get("pdg_or_equivalent")
    if isinstance(entries, Mapping):
        return entries
    if isinstance(entries, list) and entries:
        for index, entry in enumerate(entries):
            if not isinstance(entry, Mapping):
                raise AnchorError(
                    f"{process_id}: pdg_or_equivalent[{index}] is not a mapping"
                )
        return entries
    raise AnchorError(f"{process_id}: expected non-empty pdg_or_equivalent anchors")


def _find_entry(
    process_id: str,
    observable: str,
) -> tuple[str, Mapping[str, Any]]:
    entries = _pdg_container(process_id)
    if isinstance(entries, Mapping):
        direct = entries.get(observable)
        if isinstance(direct, Mapping):
            return observable, direct
        for key, entry in entries.items():
            if isinstance(entry, Mapping) and entry.get("observable") == observable:
                return str(key), entry
        present = [
            str(key)
            for key, entry in entries.items()
            if isinstance(entry, Mapping) and entry.get("observable") is not None
        ]
    else:
        for index, entry in enumerate(entries):
            if entry.get("observable") == observable:
                return f"pdg_or_equivalent[{index}]", entry
        present = [str(entry.get("observable")) for entry in entries]
    raise AnchorError(
        f"{process_id}: observable {observable!r} not found in pdg_or_equivalent "
        f"(present: {present})"
    )


def _load_scaffold_entry_anchor(
    observable: str,
    *,
    process_id: str,
) -> Anchor:
    block_key, entry = _find_entry(process_id, observable)
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
        )
    finally:
        anchor_scaffold.load_pdg_block = original_load_pdg_block

    if scaffold_anchor.process_id != process_id:
        raise AnchorError(
            f"{process_id}: load_anchor selected process {scaffold_anchor.process_id!r}, "
            f"expected standalone T011 process {process_id!r}"
        )
    if scaffold_anchor.block_key != block_key:
        raise AnchorError(
            f"{process_id}: load_anchor selected {scaffold_anchor.block_key!r}, "
            f"expected {block_key!r} for T011 observable {observable!r}"
        )
    return scaffold_anchor


def _optional_entry_anchor(observable: str, *, process_id: str) -> Anchor | None:
    try:
        return _load_scaffold_entry_anchor(observable, process_id=process_id)
    except AnchorError:
        return None


def _entry_for_anchor(anchor: Anchor, *, observable: str) -> Mapping[str, Any]:
    _, entry = _find_entry(anchor.process_id, observable)
    return entry


def _numeric_from_mapping(value: Any, *, key: str) -> Any:
    if isinstance(value, Mapping):
        return value.get(key)
    return value


def _sm_fit_from_entry(
    entry: Mapping[str, Any],
    *,
    process_id: str,
    observable: str,
    source_url: str | None,
) -> SMFitValue | None:
    value_keys = (
        "sm_fit_value",
        "sm_prediction",
        "standard_model_prediction",
        "standard_model_value",
        "sm_value",
    )
    uncertainty_keys = (
        "sm_fit_uncertainty",
        "sm_prediction_uncertainty",
        "standard_model_uncertainty",
        "sm_uncertainty",
    )
    value = next(
        (
            _numeric_from_mapping(entry[key], key="value")
            for key in value_keys
            if key in entry
        ),
        None,
    )
    uncertainty = next(
        (
            _numeric_from_mapping(entry[key], key="uncertainty")
            for key in uncertainty_keys
            if key in entry
        ),
        None,
    )
    if value is None:
        return None
    if uncertainty is None:
        raise AnchorError(
            f"{process_id}: {observable} has SM value in YAML but no SM uncertainty"
        )
    pull_value = entry.get("pull", entry.get("sm_pull", entry.get("pull_sigma")))
    return SMFitValue(
        observable=observable,
        value=_required_float(
            value,
            process_id=process_id,
            field_name=f"{observable}.sm_fit_value",
        ),
        uncertainty=_positive_float(
            uncertainty,
            process_id=process_id,
            field_name=f"{observable}.sm_fit_uncertainty",
        ),
        pull=(
            None
            if pull_value is None
            else _required_float(
                pull_value,
                process_id=process_id,
                field_name=f"{observable}.pull",
            )
        ),
        source_path=None,
        source_url=source_url,
    )


def _snapshot_path(path_value: str | None, *, process_id: str) -> Path | None:
    if path_value is None:
        return None
    path = _REPO_ROOT / path_value
    if not path.is_file():
        raise AnchorError(f"{process_id}: T011 snapshot not found: {path}")
    return path


def _snapshot_candidates(
    *,
    process_id: str,
    anchors: Sequence[Anchor | None],
) -> list[tuple[Path, str | None]]:
    seen: set[Path] = set()
    candidates: list[tuple[Path, str | None]] = []
    for anchor in anchors:
        if anchor is None:
            continue
        path = _snapshot_path(anchor.snapshot_path, process_id=process_id)
        if path is not None and path not in seen:
            seen.add(path)
            candidates.append((path, anchor.source_url))
    return candidates


def _parse_sm_fit_line(
    text: str,
    *,
    process_id: str,
    observable: str,
    source_path: Path,
    source_url: str | None,
) -> SMFitValue:
    label = _SM_SNAPSHOT_LABELS[observable]
    pattern = re.compile(
        rf"{re.escape(label)}\s*=\s*"
        rf"(?P<exp>{_NUMBER_RE})\s*\+/-\s*(?P<exp_unc>{_NUMBER_RE})"
        rf";\s*SM fit value shown there:\s*"
        rf"(?P<sm>{_NUMBER_RE})\s*\+/-\s*(?P<sm_unc>{_NUMBER_RE})"
        rf";\s*pull\s*(?P<pull>{_NUMBER_RE})",
    )
    match = pattern.search(text)
    if match is None:
        raise AnchorError(
            f"{process_id}: could not parse SM-fit value for {observable!r} "
            f"from {source_path}"
        )
    return SMFitValue(
        observable=observable,
        value=_required_float(
            match.group("sm"),
            process_id=process_id,
            field_name=f"{observable}.sm_fit_value",
        ),
        uncertainty=_positive_float(
            match.group("sm_unc"),
            process_id=process_id,
            field_name=f"{observable}.sm_fit_uncertainty",
        ),
        pull=_required_float(
            match.group("pull"),
            process_id=process_id,
            field_name=f"{observable}.sm_fit_pull",
        ),
        source_path=str(source_path.relative_to(_REPO_ROOT)),
        source_url=source_url,
    )


def _sm_fit_value(
    *,
    process_id: str,
    observable: str,
    anchor: Anchor,
    snapshot_candidates: Sequence[tuple[Path, str | None]],
) -> SMFitValue:
    entry = _entry_for_anchor(anchor, observable=observable)
    from_entry = _sm_fit_from_entry(
        entry,
        process_id=process_id,
        observable=observable,
        source_url=anchor.source_url,
    )
    if from_entry is not None:
        return from_entry

    failures: list[str] = []
    for path, source_url in snapshot_candidates:
        try:
            return _parse_sm_fit_line(
                path.read_text(),
                process_id=process_id,
                observable=observable,
                source_path=path,
                source_url=source_url,
            )
        except AnchorError as exc:
            failures.append(str(exc))
    raise AnchorError(
        f"{process_id}: no YAML SM-fit fields or parseable snapshot found for "
        f"{observable!r}; failures: {failures}"
    )


def _build_budget(
    *,
    process_id: str,
    observable: str,
    experimental: Anchor,
    sm_fit: SMFitValue,
) -> AsymmetryBudget:
    if experimental.uncertainty is None:
        raise AnchorError(f"{process_id}: {observable} must have an uncertainty")
    exp_sigma = _positive_float(
        experimental.uncertainty,
        process_id=process_id,
        field_name=f"{observable}.experimental_uncertainty",
    )
    sm_sigma = _positive_float(
        sm_fit.uncertainty,
        process_id=process_id,
        field_name=f"{observable}.sm_fit_uncertainty",
    )
    central = abs(float(experimental.value) - float(sm_fit.value))
    combined = math.sqrt(exp_sigma * exp_sigma + sm_sigma * sm_sigma)
    budget = central + combined
    source_note = (
        f"flavor_catalog/processes/{_FAMILY}/{process_id}.yaml "
        f"{observable} experimental anchor + LEP/SLC SM-fit reference"
    )
    return AsymmetryBudget(
        observable=observable,
        central_residual=float(central),
        experimental_sigma=float(exp_sigma),
        sm_fit_sigma=float(sm_sigma),
        combined_sigma=float(combined),
        hard_veto_budget=float(budget),
        source=source_note,
    )


def _load_t011_anchor(process_id: str) -> T011Anchor:
    a_fb = _load_scaffold_entry_anchor(_AFB_OBSERVABLE, process_id=process_id)
    a_b = _load_scaffold_entry_anchor(_AB_OBSERVABLE, process_id=process_id)
    a_fb_pull = _optional_entry_anchor(_AFB_PULL_OBSERVABLE, process_id=process_id)
    snapshots = _snapshot_candidates(
        process_id=process_id,
        anchors=(a_fb_pull, a_fb, a_b),
    )
    sm_fit = {
        _AFB_OBSERVABLE: _sm_fit_value(
            process_id=process_id,
            observable=_AFB_OBSERVABLE,
            anchor=a_fb,
            snapshot_candidates=snapshots,
        ),
        _AB_OBSERVABLE: _sm_fit_value(
            process_id=process_id,
            observable=_AB_OBSERVABLE,
            anchor=a_b,
            snapshot_candidates=snapshots,
        ),
    }
    budgets = {
        observable: _build_budget(
            process_id=process_id,
            observable=observable,
            experimental=a_fb if observable == _AFB_OBSERVABLE else a_b,
            sm_fit=sm_fit[observable],
        )
        for observable in _VETO_OBSERVABLES
    }
    return T011Anchor(
        a_fb=a_fb,
        a_b=a_b,
        a_fb_pull=a_fb_pull,
        sm_fit_values=sm_fit,
        budgets=budgets,
        source_process_id=process_id,
        canonical_home_fallback=False,
    )


def _asymmetry_shift(
    *,
    observable: str,
    predicted: float,
    sm_prediction: float,
    experimental: Anchor,
    budget: AsymmetryBudget,
) -> AsymmetryShift:
    np_shift = float(predicted - sm_prediction)
    ratio = abs(np_shift) / budget.hard_veto_budget
    total_exp_pull = float(predicted - experimental.value)
    total_exp_pull_ratio = abs(total_exp_pull) / budget.combined_sigma
    return AsymmetryShift(
        observable=observable,
        predicted=float(predicted),
        sm_prediction=float(sm_prediction),
        experimental=float(experimental.value),
        budget=float(budget.hard_veto_budget),
        np_shift=np_shift,
        ratio=float(ratio),
        total_exp_pull=total_exp_pull,
        total_exp_pull_ratio=float(total_exp_pull_ratio),
    )


def _shifts_for_observables(
    *,
    anchor: T011Anchor,
    prediction: ZPoleQuarkObservables,
    sm_prediction: ZPoleQuarkObservables,
) -> Mapping[str, AsymmetryShift]:
    return {
        _AFB_OBSERVABLE: _asymmetry_shift(
            observable=_AFB_OBSERVABLE,
            predicted=prediction.a_fb,
            sm_prediction=sm_prediction.a_fb,
            experimental=anchor.a_fb,
            budget=anchor.budgets[_AFB_OBSERVABLE],
        ),
        _AB_OBSERVABLE: _asymmetry_shift(
            observable=_AB_OBSERVABLE,
            predicted=prediction.a_q,
            sm_prediction=sm_prediction.a_q,
            experimental=anchor.a_b,
            budget=anchor.budgets[_AB_OBSERVABLE],
        ),
    }


def _shift_diagnostics(
    shift: AsymmetryShift,
    budget: AsymmetryBudget,
) -> dict[str, float | str]:
    return {
        "predicted": float(shift.predicted),
        "sm_prediction": float(shift.sm_prediction),
        "experimental": float(shift.experimental),
        "np_shift": float(shift.np_shift),
        "ratio": float(shift.ratio),
        "budget": float(shift.budget),
        "total_minus_experiment": float(shift.total_exp_pull),
        "total_minus_experiment_ratio_to_combined_sigma": float(
            shift.total_exp_pull_ratio
        ),
        "central_exp_minus_sm_fit": float(budget.central_residual),
        "experimental_sigma": float(budget.experimental_sigma),
        "sm_fit_sigma": float(budget.sm_fit_sigma),
        "combined_sigma": float(budget.combined_sigma),
        "budget_source": budget.source,
    }


@register
class Constraint:
    """Catalogued ``Z -> b bbar`` pole-asymmetry constraint."""

    process_id = "T011"
    severity = Severity.HARD
    observable = "max Zbb asymmetry NP shift from A_FB^0,b and A_b"

    def __init__(self) -> None:
        self.anchor = _load_t011_anchor(self.process_id)
        self.sm_inputs = zpole_default_sm_inputs()
        self.sm_observables = zpole_evaluate_quark("b", inputs=self.sm_inputs)

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        couplings = point.get_extra(_REQUIRED_EXTRA)
        kk_ew_mass = point.get_extra(_OPTIONAL_EW_MASS_EXTRA)
        proxy = None
        missing_extra = couplings is None
        if missing_extra:
            prediction = self.sm_observables
        else:
            prediction, proxy = zpole_evaluate_zbb_with_proxy(
                couplings,
                m_kk_gev=None if kk_ew_mass is None else float(kk_ew_mass),
                inputs=self.sm_inputs,
            )

        shifts = _shifts_for_observables(
            anchor=self.anchor,
            prediction=prediction,
            sm_prediction=self.sm_observables,
        )
        selected = max(shifts.values(), key=lambda item: item.ratio)
        passes = bool(selected.ratio <= 1.0)

        diagnostics: dict[str, Any] = {
            "selected_observable": selected.observable,
            "observables": {
                key: _shift_diagnostics(shift, self.anchor.budgets[key])
                for key, shift in shifts.items()
            },
            "sm_validation": {
                "a_b_formula": float(self.sm_observables.a_q),
                "a_b_lep_slc_sm_fit": float(
                    self.anchor.sm_fit_values[_AB_OBSERVABLE].value
                ),
                "a_fb_formula": float(self.sm_observables.a_fb),
                "a_fb_lep_slc_sm_fit": float(
                    self.anchor.sm_fit_values[_AFB_OBSERVABLE].value
                ),
                "a_e_formula": float(self.sm_observables.a_e),
                "sin2_theta_eff": float(self.sm_inputs.sin2_theta_eff),
            },
            "legacy_a_fb_context": {
                "experimental": float(self.anchor.a_fb.value),
                "experimental_uncertainty": float(self.anchor.a_fb.uncertainty),
                "lep_slc_sm_fit": float(
                    self.anchor.sm_fit_values[_AFB_OBSERVABLE].value
                ),
                "lep_slc_pull_sigma": self.anchor.sm_fit_values[
                    _AFB_OBSERVABLE
                ].pull,
            },
            "source_process_id": self.anchor.source_process_id,
            "canonical_home_fallback": self.anchor.canonical_home_fallback,
            "needs_human_physics": _NEEDS_HUMAN_PHYSICS,
            "rs_matching_assumption": ZPOLE_RS_ZBB_PROXY_V1,
            "required_parameter_point_extras": [_REQUIRED_EXTRA],
            "kk_ew_mass_extra_used": kk_ew_mass is not None,
        }
        if missing_extra:
            diagnostics["missing_extra"] = _REQUIRED_EXTRA
        if proxy is not None:
            diagnostics["zbb_proxy"] = dict(proxy.diagnostics)
            diagnostics["delta_g_left_b"] = float(proxy.delta_g_left_b)
            diagnostics["delta_g_right_b"] = float(proxy.delta_g_right_b)

        return ConstraintResult(
            process_id=self.process_id,
            severity=self.severity,
            passes=passes,
            predicted=float(selected.predicted),
            sm_prediction=float(selected.sm_prediction),
            experimental=float(selected.experimental),
            ratio=float(selected.ratio),
            budget=float(selected.budget),
            notes=(
                "Max NP-shift ratio over A_FB^{0,b} and A_b. SM-limit "
                "pseudo-observables use effective Zbb and Zee couplings; "
                "point-specific RS shifts use the documented Zbb "
                "coupling-shift proxy and are flagged NEEDS-HUMAN-PHYSICS. "
                "The budget is the LEP/SLC exp-vs-SM loose edge from YAML "
                "provenance."
            ),
            diagnostics=diagnostics,
        )
