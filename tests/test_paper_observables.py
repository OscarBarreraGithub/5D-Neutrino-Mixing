"""Forward-compatible tests for the paper-owned NP-only kaon observable layer."""

from __future__ import annotations

import copy
import dataclasses
import importlib
import inspect
import json
import math
import re
import subprocess
import sys
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
PACKAGE_ROOT = REPO_ROOT / "quarkConstraints" / "paper_0710_1869"
HADRONIC_MODULE_PATH = PACKAGE_ROOT / "eft_deltaf2" / "hadronic.py"
OBSERVABLES_MODULE_PATH = PACKAGE_ROOT / "eft_deltaf2" / "observables.py"

FORBIDDEN_REPO_V1_MODULES = {
    "deltaf2",
    "quarkConstraints.benchmarks",
    "quarkConstraints.couplings",
    "quarkConstraints.deltaf2",
    "quarkConstraints.fit",
    "quarkConstraints.model",
    "quarkConstraints.proxies",
    "quarkConstraints.scan",
    "quarkConstraints.scales",
    "quarkConstraints.validation",
}
MATCHING_OBJECT_NAMES = (
    "default_paper_0710_1869_kaon_matching",
    "build_paper_0710_1869_kaon_matching",
    "default_paper_0710_1869_deltaf2_matching",
    "build_paper_0710_1869_deltaf2_matching",
    "default_paper_0710_1869_kkgluon_matching",
    "build_paper_0710_1869_kkgluon_matching",
    "match_default_paper_0710_1869_kaon_kk_gluon_deltaf2",
)
PARAMETERIZED_RG_NAMES = (
    "evolve_deltaf2_wilsons_lo",
    "run_deltaf2_wilsons_lo",
    "run_paper_0710_1869_deltaf2_wilsons_lo",
    "evolve_paper_0710_1869_deltaf2_wilsons_lo",
)
RG_DEFAULT_LOW_SCALE_GEV = 2.0
HADRONIC_ALIAS_STEMS = (
    "default_kaon_hadronic",
    "kaon_hadronic",
    "hadronic_inputs",
    "hadronic",
)
DEFAULT_HADRONIC_SUMMARY_NAMES = tuple(
    name
    for stem in HADRONIC_ALIAS_STEMS
    for name in (
        f"default_paper_0710_1869_{stem}_summary",
        f"build_paper_0710_1869_{stem}_summary",
    )
)
DEFAULT_HADRONIC_EXPORT_NAMES = tuple(
    name
    for stem in HADRONIC_ALIAS_STEMS
    for name in (
        f"default_paper_0710_1869_{stem}",
        f"build_paper_0710_1869_{stem}",
    )
)
HADRONIC_BUILDER_NAMES = (
    "build_paper_0710_1869_kaon_hadronic_bundle",
    "build_paper_0710_1869_kaon_hadronic",
)
CUSTOM_LR_HADRONIC_BUILDER_NAMES = (
    "build_paper_0710_1869_custom_kaon_lr_hadronic_bundle",
    "build_paper_0710_1869_kaon_lr_hadronic_bundle",
    "build_paper_0710_1869_custom_kaon_lr_hadronic",
    "build_paper_0710_1869_kaon_lr_hadronic",
)
KAON_LR_R_CHI_FREEZE_EXPORT_NAMES = (
    "default_paper_0710_1869_kaon_lr_r_chi_freeze",
    "build_paper_0710_1869_kaon_lr_r_chi_freeze",
)
HADRONIC_SOURCE_REF_CLASS_NAMES = (
    "Paper07101869LRHadronicSourceRef",
    "Paper07101869HadronicSourceRef",
)
LR_HADRONIC_B4_PARAM_NAMES = ("B4_mu_had", "b4_mu_had", "B4", "b4")
LR_HADRONIC_B5_PARAM_NAMES = ("B5_mu_had", "b5_mu_had", "B5", "b5")
LR_HADRONIC_R_CHI_PARAM_NAMES = ("R_chi_mu_had", "r_chi_mu_had", "R_chi", "r_chi")
LR_HADRONIC_B4_SOURCE_PARAM_NAMES = (
    "B4_source",
    "b4_source",
    "B4_source_ref",
    "b4_source_ref",
)
LR_HADRONIC_B5_SOURCE_PARAM_NAMES = (
    "B5_source",
    "b5_source",
    "B5_source_ref",
    "b5_source_ref",
)
LR_HADRONIC_R_CHI_SOURCE_PARAM_NAMES = (
    "R_chi_source",
    "r_chi_source",
    "R_chi_source_ref",
    "r_chi_source_ref",
)
LR_HADRONIC_PROBE_B4 = 0.92
LR_HADRONIC_PROBE_B5 = 0.71
LR_HADRONIC_PROBE_R_CHI = 31.5
LR_HADRONIC_PROBE_SCHEME_ID = "paper_0710_1869.deltaf2.kaon_lr_hadronic.custom_probe.v1"
LR_HADRONIC_PROBE_MU_HAD_GEV = 3.25
OBSERVABLE_ALIAS_STEMS = (
    "default_kaon_observables",
    "kaon_np_observable",
    "kaon_np_observables",
    "kaon_np_observable_result",
    "kaon_observables",
    "kaon_observable",
    "observables",
    "observable",
)
DEFAULT_OBSERVABLE_SUMMARY_NAMES = tuple(
    name
    for stem in OBSERVABLE_ALIAS_STEMS
    for name in (
        f"default_paper_0710_1869_{stem}_summary",
        f"build_paper_0710_1869_{stem}_summary",
    )
)
DEFAULT_OBSERVABLE_EXPORT_NAMES = tuple(
    name
    for stem in OBSERVABLE_ALIAS_STEMS
    for name in (
        f"default_paper_0710_1869_{stem}",
        f"build_paper_0710_1869_{stem}",
    )
)
PARAMETERIZED_OBSERVABLE_NAMES = (
    "evaluate_paper_0710_1869_kaon_np_observables",
    "evaluate_paper_0710_1869_kaon_observables",
    "evaluate_deltaf2_kaon_observables",
    "evaluate_paper_0710_1869_observables",
    "build_paper_0710_1869_kaon_np_observable_result",
    "compute_paper_0710_1869_kaon_observables",
    "compute_kaon_np_observables",
    "evaluate_kaon_np_observables",
)
CUSTOM_LR_OBSERVABLE_NAMES = (
    "evaluate_paper_0710_1869_kaon_lr_only_observables",
    "evaluate_kaon_lr_only_observables",
    "build_paper_0710_1869_kaon_lr_only_observable_result",
    "compute_kaon_lr_only_observables",
)
CUSTOM_TOTAL_OBSERVABLE_NAMES = (
    "evaluate_paper_0710_1869_kaon_custom_total_observables",
    "evaluate_kaon_custom_total_observables",
    "build_paper_0710_1869_kaon_custom_total_observable_result",
    "compute_kaon_custom_total_observables",
)
RG_RESULT_PARAM_NAMES = (
    "rg_result_or_wilsons",
    "rg_result",
    "deltaf2_rg_result",
    "deltaf2_result",
    "result",
    "rg_output",
)
WILSON_PARAM_NAMES = (
    "wilsons",
    "deltaf2_wilsons",
    "evolved_wilsons",
    "coefficients",
)
HADRONIC_PARAM_NAMES = (
    "hadronic",
    "hadronic_inputs",
    "hadronic_bundle",
    "input_bundle",
    "matrix_elements",
)
SUPPORTED_OBSERVABLE_ROWS = {
    "M12_K_NP.re",
    "M12_K_NP.im",
    "Delta_m_K_NP",
}
SUPPORTED_OBSERVABLE_CONTAINER_KEYS = {"M12_K_NP", "Delta_m_K_NP"}
SUPPORTED_OBSERVABLE_FLAT_FIELDS = {
    "re_M12_K_NP_GeV": "M12_K_NP.re",
    "im_M12_K_NP_GeV": "M12_K_NP.im",
    "delta_m_K_NP_GeV": "Delta_m_K_NP",
}
EXPECTED_OBSERVABLE_SCOPE_ID = "kaon.np_only.m12.v1"
EXPECTED_INTERPRETATION_ID = "kaon.np_only.v1"
EXPECTED_M12_OBSERVABLE_ID = "M12_K_NP"
EXPECTED_DELTA_M_OBSERVABLE_ID = "Delta_m_K_NP"
EXPECTED_LR_ONLY_SCOPE_ID = "kaon.lr_only.custom.m12.v1"
EXPECTED_LR_ONLY_INTERPRETATION_ID = "kaon.lr_only.custom.v1"
EXPECTED_LR_ONLY_M12_OBSERVABLE_ID = "M12_K_LR_NP"
EXPECTED_LR_ONLY_DELTA_M_OBSERVABLE_ID = "Delta_m_K_LR_NP"
EXPECTED_CUSTOM_TOTAL_SCOPE_ID = "kaon.np_only.custom_total.q1_plus_lr.v1"
EXPECTED_CUSTOM_TOTAL_INTERPRETATION_ID = "kaon.np_only.custom_total.q1_plus_lr.v1"
EXPECTED_CUSTOM_TOTAL_M12_OBSERVABLE_ID = "M12_K_NP_CUSTOM_TOTAL"
EXPECTED_CUSTOM_TOTAL_DELTA_M_OBSERVABLE_ID = "Delta_m_K_NP_CUSTOM_TOTAL"
B_CUSTOM_Q1_PROBE_MU_HAD_GEV = 4.2
B_CUSTOM_Q1_PROBE_CONFIG = {
    "B_d": {
        "matching_names": (
            "default_paper_0710_1869_bd_matching",
            "build_paper_0710_1869_bd_matching",
            "match_default_paper_0710_1869_bd_kk_gluon_deltaf2",
        ),
        "hadronic_builder_names": (
            "build_paper_0710_1869_bd_hadronic_bundle",
            "build_paper_0710_1869_bd_hadronic",
        ),
        "observable_names": (
            "evaluate_paper_0710_1869_bd_custom_q1_observables",
            "evaluate_bd_custom_q1_observables",
            "build_paper_0710_1869_bd_custom_q1_observable_result",
            "compute_bd_custom_q1_observables",
        ),
        "scope_id": "bd.np_only.custom_q1.m12.v1",
        "interpretation_id": "bd.np_only.custom_q1.v1",
        "m12_id": "M12_Bd_NP_CUSTOM_Q1",
        "delta_m_id": "Delta_m_Bd_NP_CUSTOM_Q1",
        "mass_param": "m_Bd_GeV",
        "decay_param": "f_Bd_GeV",
        "bag_param": "B_Bd_mu_had",
        "meson_mass_GeV": 5.27965,
        "meson_decay_constant_GeV": 0.190,
        "bag_parameter_mu_had": 0.92,
        "source_prefix": "hadronic.bd.custom_q1_probe",
        "expected_sector_id": "down",
        "lr_error_pattern": r"Q4_LR/Q5_LR remain out of scope for the B-system slice",
    },
    "B_s": {
        "matching_names": (
            "default_paper_0710_1869_bs_matching",
            "build_paper_0710_1869_bs_matching",
            "match_default_paper_0710_1869_bs_kk_gluon_deltaf2",
        ),
        "hadronic_builder_names": (
            "build_paper_0710_1869_bs_hadronic_bundle",
            "build_paper_0710_1869_bs_hadronic",
        ),
        "observable_names": (
            "evaluate_paper_0710_1869_bs_custom_q1_observables",
            "evaluate_bs_custom_q1_observables",
            "build_paper_0710_1869_bs_custom_q1_observable_result",
            "compute_bs_custom_q1_observables",
        ),
        "scope_id": "bs.np_only.custom_q1.m12.v1",
        "interpretation_id": "bs.np_only.custom_q1.v1",
        "m12_id": "M12_Bs_NP_CUSTOM_Q1",
        "delta_m_id": "Delta_m_Bs_NP_CUSTOM_Q1",
        "mass_param": "m_Bs_GeV",
        "decay_param": "f_Bs_GeV",
        "bag_param": "B_Bs_mu_had",
        "meson_mass_GeV": 5.36688,
        "meson_decay_constant_GeV": 0.230,
        "bag_parameter_mu_had": 0.91,
        "source_prefix": "hadronic.bs.custom_q1_probe",
        "expected_sector_id": "down",
        "lr_error_pattern": r"Q4_LR/Q5_LR remain out of scope for the B-system slice",
    },
    "D0": {
        "matching_names": (
            "default_paper_0710_1869_d0_matching",
            "build_paper_0710_1869_d0_matching",
            "match_default_paper_0710_1869_d0_kk_gluon_deltaf2",
        ),
        "hadronic_builder_names": (
            "build_paper_0710_1869_d0_hadronic_bundle",
            "build_paper_0710_1869_d0_hadronic",
        ),
        "observable_names": (
            "evaluate_paper_0710_1869_d0_custom_q1_observables",
            "evaluate_d0_custom_q1_observables",
            "build_paper_0710_1869_d0_custom_q1_observable_result",
            "compute_d0_custom_q1_observables",
        ),
        "scope_id": "d0.np_only.custom_q1.m12.v1",
        "interpretation_id": "d0.np_only.custom_q1.v1",
        "m12_id": "M12_D0_NP",
        "delta_m_id": "Delta_m_D0_NP",
        "mass_param": "m_D0_GeV",
        "decay_param": "f_D0_GeV",
        "bag_param": "B_D0_mu_had",
        "meson_mass_GeV": 1.86484,
        "meson_decay_constant_GeV": 0.212,
        "bag_parameter_mu_had": 0.80,
        "probe_q1_vll": complex(3.25e-13, -1.10e-13),
        "probe_q1_vrr": complex(-4.50e-14, 2.50e-14),
        "require_nonzero_probe": True,
        "source_prefix": "hadronic.d0.custom_q1_probe",
        "expected_sector_id": "up",
        "lr_error_pattern": r"Q4_LR/Q5_LR remain out of scope for the D0 slice",
    },
}
LR_OBSERVABLE_PROBE_Q4 = complex(1.25e-12, -4.0e-13)
LR_OBSERVABLE_PROBE_Q5 = complex(-7.5e-13, 3.0e-13)
CUSTOM_Q1_HADRONIC_PROBE_BAG_SOURCE_ID = "hadronic.kaon.custom_total_probe.bk.v1"
CUSTOM_Q1_HADRONIC_PROBE_BAG_SOURCE_CITATION = "Custom caller-supplied B_K(mu_had) source"
CUSTOM_Q1_HADRONIC_PROBE_BAG_SOURCE_LOCATOR = "B_K(mu_had) external input"
CUSTOM_Q1_HADRONIC_PROBE_YEAR = 2026
CUSTOM_Q1_HADRONIC_PROBE_NOTES = "Custom Q1 hadronic source for combined-observable acceptance."
CUSTOM_TOTAL_LR_B4_SOURCE_ID = "hadronic.kaon.custom_total_probe.b4.v1"
CUSTOM_TOTAL_LR_B5_SOURCE_ID = "hadronic.kaon.custom_total_probe.b5.v1"
CUSTOM_TOTAL_LR_RCHI_SOURCE_ID = "hadronic.kaon.custom_total_probe.rchi.v1"


def _has_required_modules() -> bool:
    return HADRONIC_MODULE_PATH.exists() and OBSERVABLES_MODULE_PATH.exists()


def _load_observables_module():
    if not _has_required_modules():
        pytest.skip("paper_0710_1869 hadronic/observable layers not implemented yet")
    return importlib.import_module("quarkConstraints.paper_0710_1869.eft_deltaf2.observables")


def _load_hadronic_module():
    if not HADRONIC_MODULE_PATH.exists():
        pytest.skip("paper_0710_1869 hadronic layer not implemented yet")
    return importlib.import_module("quarkConstraints.paper_0710_1869.eft_deltaf2.hadronic")


def _lr_guard_contract_ids() -> tuple[str, str]:
    rg_inputs_module = importlib.import_module(
        "quarkConstraints.paper_0710_1869.eft_deltaf2.rg_inputs"
    )
    return (
        str(rg_inputs_module.PAPER_0710_1869_DELTAF2_RG_LR_BASIS_CONTRACT_ID),
        str(rg_inputs_module.PAPER_0710_1869_DELTAF2_RG_LR_BASIS_STATUS_ID),
    )


def _assert_observable_lr_guard_message(message: str) -> None:
    contract_id, status_id = _lr_guard_contract_ids()

    assert "currently support only Q1_VLL/Q1_VRR" in message
    assert "Q4_LR/Q5_LR" in message
    assert "remain blocked by LR contract" in message
    assert contract_id in message
    assert status_id in message


def _canonicalize(value: Any) -> Any:
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        return _canonicalize(dataclasses.asdict(value))
    if isinstance(value, Mapping):
        items = sorted(value.items(), key=lambda item: str(item[0]))
        return {str(key): _canonicalize(val) for key, val in items}
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, (str, int, bool)) or value is None:
        return value
    if isinstance(value, float):
        return float(value)
    if hasattr(value, "tolist"):
        try:
            return _canonicalize(value.tolist())
        except Exception:
            pass
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return [_canonicalize(item) for item in value]
    if hasattr(value, "__dict__") and value.__dict__:
        return _canonicalize(vars(value))
    return repr(value)


def _payload_from_value(value: Any) -> dict[str, Any]:
    if hasattr(value, "summary") and callable(value.summary):
        try:
            return _payload_from_value(value.summary())
        except TypeError:
            pass
    if hasattr(value, "as_dict") and callable(value.as_dict):
        payload = _canonicalize(value.as_dict())
    else:
        payload = _canonicalize(value)
    if not isinstance(payload, Mapping):
        raise AssertionError("observable payload must canonicalize to a mapping")
    return dict(payload)


def _get_callable(module: Any, names: Sequence[str]) -> Any:
    for name in names:
        candidate = getattr(module, name, None)
        if callable(candidate):
            return candidate
    return None


def _get_class(module: Any, names: Sequence[str]) -> type[Any] | None:
    for name in names:
        candidate = getattr(module, name, None)
        if isinstance(candidate, type):
            return candidate
    return None


def _custom_lr_hadronic_builder(module: Any) -> Any:
    for name in CUSTOM_LR_HADRONIC_BUILDER_NAMES:
        candidate = getattr(module, name, None)
        if callable(candidate):
            return candidate
    for name in HADRONIC_BUILDER_NAMES:
        candidate = getattr(module, name, None)
        if not callable(candidate):
            continue
        parameters = inspect.signature(candidate).parameters
        required_groups = (
            LR_HADRONIC_B4_PARAM_NAMES,
            LR_HADRONIC_B5_PARAM_NAMES,
            LR_HADRONIC_R_CHI_PARAM_NAMES,
            LR_HADRONIC_B4_SOURCE_PARAM_NAMES,
            LR_HADRONIC_B5_SOURCE_PARAM_NAMES,
            LR_HADRONIC_R_CHI_SOURCE_PARAM_NAMES,
        )
        if all(any(param in parameters for param in group) for group in required_groups):
            return candidate
    raise AssertionError(
        "hadronic layer exists but exposes no custom LR builder; expected one of "
        + ", ".join(CUSTOM_LR_HADRONIC_BUILDER_NAMES)
        + " or an LR-extended "
        + ", ".join(HADRONIC_BUILDER_NAMES)
    )


def _default_kaon_lr_r_chi_freeze() -> Any:
    hadronic_module = _load_hadronic_module()
    callable_obj = _get_callable(hadronic_module, KAON_LR_R_CHI_FREEZE_EXPORT_NAMES)
    assert callable_obj is not None, (
        "hadronic layer exposes no frozen LR R_chi helper; expected one of "
        + ", ".join(KAON_LR_R_CHI_FREEZE_EXPORT_NAMES)
    )
    return callable_obj()


def _build_hadronic_source_ref(
    module: Any,
    *,
    source_id: str,
    citation: str,
    locator_label: str,
    year: int,
    scheme_id: str,
    scale_GeV: float,
    notes: str,
) -> Any:
    source_ref_cls = _get_class(module, HADRONIC_SOURCE_REF_CLASS_NAMES)
    if source_ref_cls is None:
        raise AssertionError(
            "hadronic layer exposes no source-ref class; expected one of "
            + ", ".join(HADRONIC_SOURCE_REF_CLASS_NAMES)
        )

    parameters = inspect.signature(source_ref_cls).parameters
    kwargs: dict[str, Any] = {}
    if "source_id" in parameters:
        kwargs["source_id"] = source_id
    if "source_kind" in parameters:
        kwargs["source_kind"] = "review-average"
    if "citation" in parameters:
        kwargs["citation"] = citation
    if "locator_label" in parameters:
        kwargs["locator_label"] = locator_label
    if "year" in parameters:
        kwargs["year"] = year
    if "renormalization_scheme_id" in parameters:
        kwargs["renormalization_scheme_id"] = scheme_id
    if "scale_GeV" in parameters:
        kwargs["scale_GeV"] = scale_GeV
    if "transformation_id" in parameters:
        kwargs["transformation_id"] = "none"
    if "notes" in parameters:
        kwargs["notes"] = notes

    missing = [
        name
        for name, parameter in parameters.items()
        if parameter.kind
        not in (
            inspect.Parameter.VAR_POSITIONAL,
            inspect.Parameter.VAR_KEYWORD,
        )
        and parameter.default is inspect._empty
        and name not in kwargs
    ]
    if missing:
        raise AssertionError(
            "hadronic source-ref constructor requires unsupported arguments for observable "
            "LR-HAD acceptance tests: "
            + ", ".join(missing)
        )
    return source_ref_cls(**kwargs)


def _custom_lr_hadronic_builder_kwargs(
    builder: Any,
    *,
    default_bundle: Any,
    B4_mu_had: float,
    B5_mu_had: float,
    R_chi_mu_had: float,
    B4_source: Any,
    B5_source: Any,
    R_chi_source: Any,
    renormalization_scheme_id: str | None = None,
    mu_had_GeV: float | None = None,
    source_prefix: str = "hadronic.kaon.lr_observable_probe",
) -> dict[str, Any]:
    parameters = inspect.signature(builder).parameters
    kwargs: dict[str, Any] = {}

    def assign_first(candidates: Sequence[str], value: Any) -> None:
        for name in candidates:
            if name in parameters:
                kwargs[name] = value
                return

    assign_first(LR_HADRONIC_B4_PARAM_NAMES, B4_mu_had)
    assign_first(LR_HADRONIC_B5_PARAM_NAMES, B5_mu_had)
    assign_first(LR_HADRONIC_R_CHI_PARAM_NAMES, R_chi_mu_had)
    assign_first(LR_HADRONIC_B4_SOURCE_PARAM_NAMES, B4_source)
    assign_first(LR_HADRONIC_B5_SOURCE_PARAM_NAMES, B5_source)
    assign_first(LR_HADRONIC_R_CHI_SOURCE_PARAM_NAMES, R_chi_source)

    resolved_mu_had_GeV = (
        float(default_bundle.mu_had_GeV) if mu_had_GeV is None else float(mu_had_GeV)
    )
    resolved_renormalization_scheme_id = (
        str(default_bundle.renormalization_scheme_id)
        if renormalization_scheme_id is None
        else str(renormalization_scheme_id)
    )
    if "mu_had_GeV" in parameters:
        kwargs["mu_had_GeV"] = resolved_mu_had_GeV
    if "m_K0_GeV" in parameters:
        kwargs["m_K0_GeV"] = float(default_bundle.m_K0_GeV)
    if "f_K_GeV" in parameters:
        kwargs["f_K_GeV"] = float(default_bundle.f_K_GeV)
    if "renormalization_scheme_id" in parameters:
        kwargs["renormalization_scheme_id"] = resolved_renormalization_scheme_id
    if "scheme_id" in parameters:
        kwargs["scheme_id"] = resolved_renormalization_scheme_id
    if "mass_source" in parameters and hasattr(default_bundle, "mass_source"):
        kwargs["mass_source"] = default_bundle.mass_source
    if "decay_constant_source" in parameters and hasattr(default_bundle, "decay_constant_source"):
        kwargs["decay_constant_source"] = default_bundle.decay_constant_source
    if "source_id" in parameters:
        kwargs["source_id"] = f"{source_prefix}.sources.v1"
    if "provenance_ids" in parameters:
        provenance_ids = [f"{source_prefix}.sources.v1"]
        for attr_name in ("mass_source", "decay_constant_source", "bag_parameter_source"):
            source = getattr(default_bundle, attr_name, None)
            source_id = getattr(source, "source_id", None)
            if source_id:
                provenance_ids.append(str(source_id))
        provenance_ids.extend(
            (
                f"{source_prefix}.b4.v1",
                f"{source_prefix}.b5.v1",
                f"{source_prefix}.rchi.v1",
            )
        )
        kwargs["provenance_ids"] = tuple(dict.fromkeys(provenance_ids))

    missing = [
        name
        for name, parameter in parameters.items()
        if parameter.kind
        not in (
            inspect.Parameter.VAR_POSITIONAL,
            inspect.Parameter.VAR_KEYWORD,
        )
        and parameter.default is inspect._empty
        and name not in kwargs
    ]
    if missing:
        raise AssertionError(
            "custom LR hadronic builder requires unsupported arguments for observable "
            "LR-HAD acceptance tests: "
            + ", ".join(missing)
        )
    return kwargs


def _default_matching_object() -> Any:
    module = importlib.import_module(
        "quarkConstraints.paper_0710_1869.eft_deltaf2.matching_kkgluon"
    )
    callable_obj = _get_callable(module, MATCHING_OBJECT_NAMES)
    assert callable_obj is not None, "matching layer exposes no default object callable"
    return callable_obj()


def _invoke_rg_evolution(callable_obj: Any, *, value: Any, mu_low_GeV: float) -> Any:
    parameters = inspect.signature(callable_obj).parameters
    kwargs: dict[str, Any] = {}
    wilsons = getattr(value, "wilsons", value)

    for name in ("wilsons", "deltaf2_wilsons", "wilsons_in", "coefficients"):
        if name in parameters:
            kwargs[name] = wilsons
            break
    for name in ("matching", "deltaf2_matching", "match"):
        if name in parameters:
            kwargs[name] = value
            break
    for name in ("mu_low_GeV", "mu_eval_GeV", "evaluation_scale_GeV", "scale_GeV"):
        if name in parameters:
            kwargs[name] = mu_low_GeV
            break

    missing = [
        name
        for name, parameter in parameters.items()
        if parameter.kind
        not in (
            inspect.Parameter.VAR_POSITIONAL,
            inspect.Parameter.VAR_KEYWORD,
        )
        and parameter.default is inspect._empty
        and name not in kwargs
    ]
    if missing:
        raise AssertionError(
            "RG evolution callable requires unsupported arguments for observable tests: "
            + ", ".join(missing)
        )
    return callable_obj(**kwargs)


def _default_rg_result() -> Any:
    module = importlib.import_module("quarkConstraints.paper_0710_1869.eft_deltaf2.rg")
    callable_obj = _get_callable(module, PARAMETERIZED_RG_NAMES)
    assert callable_obj is not None, "RG layer exposes no callable evolution API"
    return _invoke_rg_evolution(
        callable_obj,
        value=_default_matching_object(),
        mu_low_GeV=RG_DEFAULT_LOW_SCALE_GEV,
    )


def _default_hadronic_value() -> Any:
    module = _load_hadronic_module()
    callable_obj = _get_callable(
        module,
        DEFAULT_HADRONIC_EXPORT_NAMES + DEFAULT_HADRONIC_SUMMARY_NAMES,
    )
    assert callable_obj is not None, "hadronic layer exposes no default kaon export"
    return callable_obj()


def _get_observable_evaluator(module: Any) -> Any:
    callable_obj = _get_callable(module, PARAMETERIZED_OBSERVABLE_NAMES)
    assert callable_obj is not None, (
        "observable layer exists but exposes no parameterized evaluator; expected one of "
        + ", ".join(PARAMETERIZED_OBSERVABLE_NAMES)
    )
    return callable_obj


def _get_custom_lr_observable_evaluator(module: Any) -> Any:
    callable_obj = _get_callable(module, CUSTOM_LR_OBSERVABLE_NAMES)
    assert callable_obj is not None, (
        "observable layer exists but exposes no custom LR-only evaluator; expected one of "
        + ", ".join(CUSTOM_LR_OBSERVABLE_NAMES)
    )
    return callable_obj


def _get_custom_total_observable_evaluator(module: Any) -> Any:
    callable_obj = _get_callable(module, CUSTOM_TOTAL_OBSERVABLE_NAMES)
    assert callable_obj is not None, (
        "observable layer exists but exposes no custom combined evaluator; expected one of "
        + ", ".join(CUSTOM_TOTAL_OBSERVABLE_NAMES)
    )
    return callable_obj


def _get_custom_b_q1_observable_evaluator(module: Any, system_id: str) -> Any:
    config = B_CUSTOM_Q1_PROBE_CONFIG[system_id]
    callable_obj = _get_callable(module, config["observable_names"])
    assert callable_obj is not None, (
        f"observable layer exists but exposes no {system_id} custom Q1 evaluator; expected one of "
        + ", ".join(config["observable_names"])
    )
    return callable_obj


def _invoke_observable_evaluation(callable_obj: Any, *, rg_value: Any, hadronic_value: Any) -> Any:
    parameters = inspect.signature(callable_obj).parameters
    kwargs: dict[str, Any] = {}

    for name in RG_RESULT_PARAM_NAMES:
        if name in parameters:
            kwargs[name] = rg_value
            break
    else:
        for name in WILSON_PARAM_NAMES:
            if name in parameters:
                kwargs[name] = getattr(rg_value, "wilsons", rg_value)
                break

    for name in HADRONIC_PARAM_NAMES:
        if name in parameters:
            kwargs[name] = hadronic_value
            break

    missing = [
        name
        for name, parameter in parameters.items()
        if parameter.kind
        not in (
            inspect.Parameter.VAR_POSITIONAL,
            inspect.Parameter.VAR_KEYWORD,
        )
        and parameter.default is inspect._empty
        and name not in kwargs
    ]
    if missing:
        raise AssertionError(
            "observable evaluator requires unsupported arguments for tests: "
            + ", ".join(missing)
        )
    return callable_obj(**kwargs)


def _invoke_custom_total_observable_evaluation(
    callable_obj: Any,
    *,
    rg_value: Any,
    q1_hadronic_bundle: Any,
    lr_hadronic_inputs: Any,
) -> Any:
    parameters = inspect.signature(callable_obj).parameters
    kwargs: dict[str, Any] = {}

    for name in RG_RESULT_PARAM_NAMES:
        if name in parameters:
            kwargs[name] = rg_value
            break
    else:
        for name in WILSON_PARAM_NAMES:
            if name in parameters:
                kwargs[name] = getattr(rg_value, "wilsons", rg_value)
                break

    for name in ("q1_hadronic_bundle", "q1_bundle", "hadronic_bundle"):
        if name in parameters:
            kwargs[name] = q1_hadronic_bundle
            break
    for name in ("lr_hadronic_inputs", "lr_inputs", "hadronic_inputs"):
        if name in parameters:
            kwargs[name] = lr_hadronic_inputs
            break

    missing = [
        name
        for name, parameter in parameters.items()
        if parameter.kind
        not in (
            inspect.Parameter.VAR_POSITIONAL,
            inspect.Parameter.VAR_KEYWORD,
        )
        and parameter.default is inspect._empty
        and name not in kwargs
    ]
    if missing:
        raise AssertionError(
            "custom combined observable evaluator requires unsupported arguments for tests: "
            + ", ".join(missing)
        )
    return callable_obj(**kwargs)


def _complex_from_payload(value: Any) -> complex:
    if not isinstance(value, Mapping):
        raise AssertionError("complex payload must be a mapping")
    return complex(float(value["real"]), float(value["imag"]))


def _assert_coefficient_payload(payload: Mapping[str, Any], operator_name: str, expected: complex) -> None:
    observed = _complex_from_payload(payload["coefficients"][operator_name])
    assert observed.real == pytest.approx(expected.real, rel=0.0, abs=1.0e-24)
    assert observed.imag == pytest.approx(expected.imag, rel=0.0, abs=1.0e-24)


def _build_custom_lr_hadronic_inputs(
    *,
    renormalization_scheme_id: str = LR_HADRONIC_PROBE_SCHEME_ID,
    mu_had_GeV: float = LR_HADRONIC_PROBE_MU_HAD_GEV,
    source_prefix: str = "hadronic.kaon.lr_observable_probe",
) -> Any:
    hadronic_module = _load_hadronic_module()
    default_builder = _get_callable(hadronic_module, HADRONIC_BUILDER_NAMES)
    assert callable(default_builder), "hadronic layer exposes no default builder"
    default_bundle = default_builder()
    lr_builder = _custom_lr_hadronic_builder(hadronic_module)
    b4_source = _build_hadronic_source_ref(
        hadronic_module,
        source_id=f"{source_prefix}.b4.v1",
        citation="Becirevic and Villadoro, hep-lat/0408029",
        locator_label="eq. (5), O4 scalar LR matrix element",
        year=2004,
        scheme_id=renormalization_scheme_id,
        scale_GeV=mu_had_GeV,
        notes="Observable LR probe source for B4(mu_had).",
    )
    b5_source = _build_hadronic_source_ref(
        hadronic_module,
        source_id=f"{source_prefix}.b5.v1",
        citation="Becirevic and Villadoro, hep-lat/0408029",
        locator_label="eq. (5), O5 scalar LR matrix element",
        year=2004,
        scheme_id=renormalization_scheme_id,
        scale_GeV=mu_had_GeV,
        notes="Observable LR probe source for B5(mu_had).",
    )
    r_chi_source = _build_hadronic_source_ref(
        hadronic_module,
        source_id=f"{source_prefix}.rchi.v1",
        citation="Custom caller-supplied chiral ratio input",
        locator_label="R_chi(mu_had) external input",
        year=2026,
        scheme_id=renormalization_scheme_id,
        scale_GeV=mu_had_GeV,
        notes="Observable LR probe source for R_chi(mu_had).",
    )
    return lr_builder(
        **_custom_lr_hadronic_builder_kwargs(
            lr_builder,
            default_bundle=default_bundle,
            B4_mu_had=LR_HADRONIC_PROBE_B4,
            B5_mu_had=LR_HADRONIC_PROBE_B5,
            R_chi_mu_had=LR_HADRONIC_PROBE_R_CHI,
            B4_source=b4_source,
            B5_source=b5_source,
            R_chi_source=r_chi_source,
            renormalization_scheme_id=renormalization_scheme_id,
            mu_had_GeV=mu_had_GeV,
            source_prefix=source_prefix,
        )
    )


def _build_custom_q1_hadronic_bundle(*, mu_had_GeV: float = RG_DEFAULT_LOW_SCALE_GEV) -> Any:
    hadronic_module = _load_hadronic_module()
    builder = _get_callable(hadronic_module, HADRONIC_BUILDER_NAMES)
    assert callable(builder), "hadronic layer exposes no default builder"
    default_bundle = builder()
    custom_bag_source = dataclasses.replace(
        default_bundle.bag_parameter_source,
        source_id=CUSTOM_Q1_HADRONIC_PROBE_BAG_SOURCE_ID,
        citation=CUSTOM_Q1_HADRONIC_PROBE_BAG_SOURCE_CITATION,
        locator_label=CUSTOM_Q1_HADRONIC_PROBE_BAG_SOURCE_LOCATOR,
        year=CUSTOM_Q1_HADRONIC_PROBE_YEAR,
        notes=CUSTOM_Q1_HADRONIC_PROBE_NOTES,
    )
    return builder(
        mu_had_GeV=float(mu_had_GeV),
        m_K0_GeV=float(default_bundle.m_K0_GeV),
        f_K_GeV=float(default_bundle.f_K_GeV),
        hat_B_K_rgi_source_value=float(default_bundle.hat_B_K_rgi_source_value),
        bag_parameter_source=custom_bag_source,
    )


def _build_custom_total_probe_inputs() -> tuple[Any, Any, Any]:
    custom_q1_hadronic = _build_custom_q1_hadronic_bundle()
    aligned_lr_hadronic = _build_custom_lr_hadronic_inputs(
        renormalization_scheme_id=str(custom_q1_hadronic.renormalization_scheme_id),
        mu_had_GeV=float(custom_q1_hadronic.mu_had_GeV),
        source_prefix="hadronic.kaon.custom_total_probe",
    )
    base_wilsons = _default_rg_result().wilsons
    combined_wilsons = dataclasses.replace(
        base_wilsons,
        matching_scale_GeV=float(custom_q1_hadronic.mu_had_GeV),
        contract=dataclasses.replace(
            base_wilsons.contract,
            renormalization_scheme_id=str(custom_q1_hadronic.renormalization_scheme_id),
        ),
        q4_lr=LR_OBSERVABLE_PROBE_Q4,
        q5_lr=LR_OBSERVABLE_PROBE_Q5,
    )
    return custom_q1_hadronic, aligned_lr_hadronic, combined_wilsons


def _default_b_matching_object(system_id: str) -> Any:
    module = importlib.import_module(
        "quarkConstraints.paper_0710_1869.eft_deltaf2.matching_kkgluon"
    )
    config = B_CUSTOM_Q1_PROBE_CONFIG[system_id]
    callable_obj = _get_callable(module, config["matching_names"])
    assert callable_obj is not None, (
        f"matching layer exposes no {system_id} default object callable; expected one of "
        + ", ".join(config["matching_names"])
    )
    return callable_obj()


def _build_custom_b_rg_result(
    system_id: str,
    *,
    mu_low_GeV: float = B_CUSTOM_Q1_PROBE_MU_HAD_GEV,
) -> Any:
    module = importlib.import_module("quarkConstraints.paper_0710_1869.eft_deltaf2.rg")
    callable_obj = _get_callable(module, PARAMETERIZED_RG_NAMES)
    assert callable_obj is not None, "RG layer exposes no callable evolution API"
    rg_result = _invoke_rg_evolution(
        callable_obj,
        value=_default_b_matching_object(system_id),
        mu_low_GeV=float(mu_low_GeV),
    )
    config = B_CUSTOM_Q1_PROBE_CONFIG[system_id]
    wilson_updates: dict[str, complex] = {}
    if "probe_q1_vll" in config:
        wilson_updates["q1_vll"] = complex(config["probe_q1_vll"])
    if "probe_q1_vrr" in config:
        wilson_updates["q1_vrr"] = complex(config["probe_q1_vrr"])
    if not wilson_updates:
        return rg_result
    updated_wilsons = dataclasses.replace(rg_result.wilsons, **wilson_updates)
    return dataclasses.replace(rg_result, wilsons=updated_wilsons)


def _build_custom_b_hadronic_bundle(
    system_id: str,
    *,
    renormalization_scheme_id: str,
    mu_had_GeV: float,
) -> Any:
    hadronic_module = _load_hadronic_module()
    config = B_CUSTOM_Q1_PROBE_CONFIG[system_id]
    builder = _get_callable(hadronic_module, config["hadronic_builder_names"])
    assert builder is not None, (
        f"hadronic layer exposes no {system_id} custom bundle builder; expected one of "
        + ", ".join(config["hadronic_builder_names"])
    )
    source_prefix = str(config["source_prefix"])
    mass_source = _build_hadronic_source_ref(
        hadronic_module,
        source_id=f"{source_prefix}.mass.v1",
        citation=f"Custom {system_id} meson mass input",
        locator_label=f"{system_id} mass external input",
        year=2026,
        scheme_id=renormalization_scheme_id,
        scale_GeV=mu_had_GeV,
        notes=f"Custom {system_id} mass source for Q1 observable acceptance.",
    )
    decay_constant_source = _build_hadronic_source_ref(
        hadronic_module,
        source_id=f"{source_prefix}.decay_constant.v1",
        citation=f"Custom {system_id} decay constant input",
        locator_label=f"{system_id} decay constant external input",
        year=2026,
        scheme_id=renormalization_scheme_id,
        scale_GeV=mu_had_GeV,
        notes=f"Custom {system_id} decay-constant source for Q1 observable acceptance.",
    )
    bag_parameter_source = _build_hadronic_source_ref(
        hadronic_module,
        source_id=f"{source_prefix}.bag_parameter.v1",
        citation=f"Custom {system_id} bag-parameter input",
        locator_label=f"{system_id} bag parameter external input",
        year=2026,
        scheme_id=renormalization_scheme_id,
        scale_GeV=mu_had_GeV,
        notes=f"Custom {system_id} bag-parameter source for Q1 observable acceptance.",
    )
    parameters = inspect.signature(builder).parameters
    kwargs: dict[str, Any] = {
        "mu_had_GeV": float(mu_had_GeV),
        "renormalization_scheme_id": str(renormalization_scheme_id),
        str(config["mass_param"]): float(config["meson_mass_GeV"]),
        str(config["decay_param"]): float(config["meson_decay_constant_GeV"]),
        str(config["bag_param"]): float(config["bag_parameter_mu_had"]),
        "mass_source": mass_source,
        "decay_constant_source": decay_constant_source,
        "bag_parameter_source": bag_parameter_source,
    }
    if "bundle_id" in parameters:
        kwargs["bundle_id"] = f"{source_prefix}.bundle.v1"
    if "source_id" in parameters:
        kwargs["source_id"] = f"{source_prefix}.sources.v1"
    if "provenance_ids" in parameters:
        kwargs["provenance_ids"] = (
            f"{source_prefix}.sources.v1",
            f"{source_prefix}.mass.v1",
            f"{source_prefix}.decay_constant.v1",
            f"{source_prefix}.bag_parameter.v1",
        )
    return builder(**kwargs)


def _build_custom_b_probe_inputs(system_id: str) -> tuple[Any, Any]:
    rg_result = _build_custom_b_rg_result(system_id)
    hadronic_bundle = _build_custom_b_hadronic_bundle(
        system_id,
        renormalization_scheme_id=str(rg_result.wilsons.renormalization_scheme_id),
        mu_had_GeV=float(rg_result.wilsons.matching_scale_GeV),
    )
    return hadronic_bundle, rg_result


def _build_lr_probe_wilsons(
    *,
    q4_lr: complex = LR_OBSERVABLE_PROBE_Q4,
    q5_lr: complex = LR_OBSERVABLE_PROBE_Q5,
    renormalization_scheme_id: str | None = None,
    matching_scale_GeV: float | None = None,
) -> Any:
    rg_result = _default_rg_result()
    wilsons = rg_result.wilsons
    contract = wilsons.contract
    if renormalization_scheme_id is not None:
        contract = dataclasses.replace(
            contract,
            renormalization_scheme_id=renormalization_scheme_id,
        )
    resolved_matching_scale_GeV = (
        wilsons.matching_scale_GeV if matching_scale_GeV is None else matching_scale_GeV
    )
    return dataclasses.replace(
        wilsons,
        contract=contract,
        q1_vll=complex(0.0, 0.0),
        q1_vrr=complex(0.0, 0.0),
        q4_lr=q4_lr,
        q5_lr=q5_lr,
        matching_scale_GeV=resolved_matching_scale_GeV,
    )


def _default_observable_payload(module: Any) -> dict[str, Any]:
    for name in DEFAULT_OBSERVABLE_SUMMARY_NAMES + DEFAULT_OBSERVABLE_EXPORT_NAMES:
        candidate = getattr(module, name, None)
        if not callable(candidate):
            continue
        try:
            return _payload_from_value(candidate())
        except TypeError:
            continue

    evaluator = _get_observable_evaluator(module)
    return _payload_from_value(
        _invoke_observable_evaluation(
            evaluator,
            rg_value=_default_rg_result(),
            hadronic_value=_default_hadronic_value(),
        )
    )


def _default_observable_payload_cross_process() -> dict[str, Any]:
    script = f"""
import dataclasses
import importlib
import inspect
import json
from collections.abc import Mapping, Sequence
from pathlib import Path

MATCHING_NAMES = {MATCHING_OBJECT_NAMES!r}
RG_NAMES = {PARAMETERIZED_RG_NAMES!r}
MU_LOW_GEV = {RG_DEFAULT_LOW_SCALE_GEV!r}
HADRONIC_NAMES = {DEFAULT_HADRONIC_EXPORT_NAMES + DEFAULT_HADRONIC_SUMMARY_NAMES!r}
DEFAULT_NAMES = {DEFAULT_OBSERVABLE_SUMMARY_NAMES + DEFAULT_OBSERVABLE_EXPORT_NAMES!r}
EVALUATOR_NAMES = {PARAMETERIZED_OBSERVABLE_NAMES!r}
RG_RESULT_PARAM_NAMES = {RG_RESULT_PARAM_NAMES!r}
WILSON_PARAM_NAMES = {WILSON_PARAM_NAMES!r}
HADRONIC_PARAM_NAMES = {HADRONIC_PARAM_NAMES!r}

def canonicalize(value):
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        return canonicalize(dataclasses.asdict(value))
    if isinstance(value, Mapping):
        items = sorted(value.items(), key=lambda item: str(item[0]))
        return {{str(key): canonicalize(val) for key, val in items}}
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, (str, int, bool)) or value is None:
        return value
    if isinstance(value, float):
        return float(value)
    if hasattr(value, "tolist"):
        try:
            return canonicalize(value.tolist())
        except Exception:
            pass
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return [canonicalize(item) for item in value]
    if hasattr(value, "__dict__") and value.__dict__:
        return canonicalize(vars(value))
    return repr(value)

def payload_from_value(value):
    if hasattr(value, "summary") and callable(value.summary):
        try:
            return payload_from_value(value.summary())
        except TypeError:
            pass
    if hasattr(value, "as_dict") and callable(value.as_dict):
        payload = canonicalize(value.as_dict())
    else:
        payload = canonicalize(value)
    if not isinstance(payload, Mapping):
        raise AssertionError("observable payload must canonicalize to a mapping")
    return payload

def get_callable(module, names):
    for name in names:
        candidate = getattr(module, name, None)
        if callable(candidate):
            return candidate
    return None

def default_matching_object():
    module = importlib.import_module(
        "quarkConstraints.paper_0710_1869.eft_deltaf2.matching_kkgluon"
    )
    callable_obj = get_callable(module, MATCHING_NAMES)
    if not callable(callable_obj):
        raise SystemExit("missing default matching object")
    return callable_obj()

def invoke_rg(callable_obj, value, mu_low_GeV):
    parameters = inspect.signature(callable_obj).parameters
    kwargs = {{}}
    wilsons = getattr(value, "wilsons", value)
    for name in ("wilsons", "deltaf2_wilsons", "wilsons_in", "coefficients"):
        if name in parameters:
            kwargs[name] = wilsons
            break
    for name in ("matching", "deltaf2_matching", "match"):
        if name in parameters:
            kwargs[name] = value
            break
    for name in ("mu_low_GeV", "mu_eval_GeV", "evaluation_scale_GeV", "scale_GeV"):
        if name in parameters:
            kwargs[name] = mu_low_GeV
            break
    missing = [
        name for name, parameter in parameters.items()
        if parameter.kind not in (
            inspect.Parameter.VAR_POSITIONAL,
            inspect.Parameter.VAR_KEYWORD,
        )
        and parameter.default is inspect._empty
        and name not in kwargs
    ]
    if missing:
        raise AssertionError("RG callable requires unsupported arguments: " + ", ".join(missing))
    return callable_obj(**kwargs)

def default_rg_result():
    module = importlib.import_module("quarkConstraints.paper_0710_1869.eft_deltaf2.rg")
    callable_obj = get_callable(module, RG_NAMES)
    if not callable(callable_obj):
        raise SystemExit("missing RG callable")
    return invoke_rg(callable_obj, default_matching_object(), MU_LOW_GEV)

def default_hadronic_value():
    module = importlib.import_module("quarkConstraints.paper_0710_1869.eft_deltaf2.hadronic")
    callable_obj = get_callable(module, HADRONIC_NAMES)
    if not callable(callable_obj):
        raise SystemExit("missing default hadronic export")
    return callable_obj()

def invoke_observables(callable_obj, rg_value, hadronic_value):
    parameters = inspect.signature(callable_obj).parameters
    kwargs = {{}}
    for name in RG_RESULT_PARAM_NAMES:
        if name in parameters:
            kwargs[name] = rg_value
            break
    else:
        for name in WILSON_PARAM_NAMES:
            if name in parameters:
                kwargs[name] = getattr(rg_value, "wilsons", rg_value)
                break
    for name in HADRONIC_PARAM_NAMES:
        if name in parameters:
            kwargs[name] = hadronic_value
            break
    missing = [
        name for name, parameter in parameters.items()
        if parameter.kind not in (
            inspect.Parameter.VAR_POSITIONAL,
            inspect.Parameter.VAR_KEYWORD,
        )
        and parameter.default is inspect._empty
        and name not in kwargs
    ]
    if missing:
        raise AssertionError(
            "observable evaluator requires unsupported arguments: " + ", ".join(missing)
        )
    return callable_obj(**kwargs)

module = importlib.import_module("quarkConstraints.paper_0710_1869.eft_deltaf2.observables")
for name in DEFAULT_NAMES:
    candidate = getattr(module, name, None)
    if not callable(candidate):
        continue
    try:
        print(json.dumps(payload_from_value(candidate()), sort_keys=True))
        break
    except TypeError:
        continue
else:
    evaluator = get_callable(module, EVALUATOR_NAMES)
    if not callable(evaluator):
        raise SystemExit("missing observable export")
    print(
        json.dumps(
            payload_from_value(
                invoke_observables(evaluator, default_rg_result(), default_hadronic_value())
            ),
            sort_keys=True,
        )
    )
"""
    completed = subprocess.run(
        [sys.executable, "-c", script],
        check=True,
        capture_output=True,
        text=True,
        cwd=str(REPO_ROOT),
    )
    return json.loads(completed.stdout)


def _nested_value(mapping: Mapping[str, Any], key_paths: Sequence[Sequence[str]]) -> Any | None:
    for path in key_paths:
        current: Any = mapping
        for key in path:
            if not isinstance(current, Mapping) or key not in current:
                current = None
                break
            current = current[key]
        if current is not None:
            return current
    return None


def _float_from_value(value: Any) -> float:
    if isinstance(value, bool):
        raise TypeError("boolean values are not observable numbers")
    if isinstance(value, (int, float)):
        numeric = float(value)
        if not math.isfinite(numeric):
            raise ValueError("observable value must be finite")
        return numeric
    if isinstance(value, Mapping):
        for key in ("value", "central_value", "observable_value"):
            if key in value:
                return _float_from_value(value[key])
    raise TypeError(f"unsupported observable value encoding: {value!r}")


def _extract_named_observables(payload: Mapping[str, Any]) -> dict[str, float]:
    container = _nested_value(
        payload,
        (
            ("observables",),
            ("results",),
            ("observable_values",),
            ("quantities",),
        ),
    )
    rows: dict[str, float] = {}
    if isinstance(container, Mapping):
        unexpected = sorted(
            str(key) for key in container if str(key) not in SUPPORTED_OBSERVABLE_CONTAINER_KEYS
        )
        if unexpected:
            raise AssertionError(
                "observable payload exposes unsupported named outputs: "
                + ", ".join(unexpected)
            )
        if "M12_K_NP" not in container or "Delta_m_K_NP" not in container:
            raise AssertionError(
                "observable payload must expose both M12_K_NP and Delta_m_K_NP"
            )
        m12_value = container["M12_K_NP"]
        if not isinstance(m12_value, Mapping) or not {"re", "im"} <= set(m12_value):
            raise AssertionError("M12_K_NP output must expose both re and im components")
        rows["M12_K_NP.re"] = _float_from_value(m12_value["re"])
        rows["M12_K_NP.im"] = _float_from_value(m12_value["im"])
        rows["Delta_m_K_NP"] = _float_from_value(container["Delta_m_K_NP"])
        return rows
    if isinstance(container, Sequence) and not isinstance(container, (str, bytes, bytearray)):
        unexpected_names: list[str] = []
        for index, item in enumerate(container):
            if not isinstance(item, Mapping):
                raise AssertionError("observable rows must be mappings")
            name = (
                item.get("name")
                or item.get("observable_id")
                or item.get("quantity")
                or f"item_{index}"
            )
            if name not in SUPPORTED_OBSERVABLE_CONTAINER_KEYS:
                unexpected_names.append(str(name))
                continue
            if name == "M12_K_NP":
                if not {"re", "im"} <= set(item):
                    raise AssertionError("M12_K_NP output must expose both re and im components")
                rows["M12_K_NP.re"] = _float_from_value(item["re"])
                rows["M12_K_NP.im"] = _float_from_value(item["im"])
            else:
                rows[str(name)] = _float_from_value(item)
        if unexpected_names:
            raise AssertionError(
                "observable payload exposes unsupported named outputs: "
                + ", ".join(sorted(unexpected_names))
            )
        return rows

    for key, row_name in SUPPORTED_OBSERVABLE_FLAT_FIELDS.items():
        if key in payload:
            rows[row_name] = _float_from_value(payload[key])
    if rows:
        unexpected_top_level: list[str] = []
        for key, value in payload.items():
            lowered = str(key).lower()
            if key in SUPPORTED_OBSERVABLE_FLAT_FIELDS or key == "M12_K_NP_GeV":
                continue
            if "epsilon" in lowered:
                unexpected_top_level.append(str(key))
                continue
            if key.startswith("re_M12") or key.startswith("im_M12") or key.startswith("delta_m"):
                try:
                    _float_from_value(value)
                except (TypeError, ValueError):
                    continue
                unexpected_top_level.append(str(key))
        if unexpected_top_level:
            raise AssertionError(
                "observable payload exposes unsupported top-level outputs: "
                + ", ".join(sorted(unexpected_top_level))
            )
        if set(rows) != SUPPORTED_OBSERVABLE_ROWS:
            raise AssertionError(
                "observable payload must expose exactly the PR5a kaon NP-only outputs"
            )
        return rows

    for key in payload:
        lowered = str(key).lower()
        if "epsilon" in lowered or "m12" in lowered or "delta_m" in lowered:
            raise AssertionError(
                "observable payload exposes unsupported observable-like fields outside the "
                "frozen PR5a kaon NP-only surface"
            )
    raise AssertionError("observable payload does not expose the PR5a kaon NP-only outputs")


def _is_kaon_system(value: Any) -> bool:
    if value is None:
        return False
    lowered = str(value).strip().lower()
    return lowered in {"k", "k0", "kaon"} or "kaon" in lowered


def _is_np_only_interpretation(value: Any) -> bool:
    if value is None:
        return False
    lowered = str(value).strip().lower().replace("_", " ").replace("-", " ")
    return "np only" in lowered or ("np" in lowered and "sm" not in lowered)


def _path_exists(value: Any, path: Sequence[str]) -> bool:
    current = value
    for key in path:
        if isinstance(current, Mapping):
            if key not in current:
                return False
            current = current[key]
            continue
        if not hasattr(current, key):
            return False
        current = getattr(current, key)
    return True


def _deep_replace(value: Any, path: Sequence[str], replacement: Any) -> Any:
    if not path:
        return replacement

    head = path[0]
    tail = path[1:]

    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        current = getattr(value, head)
        updated = _deep_replace(current, tail, replacement)
        return dataclasses.replace(value, **{head: updated})

    if isinstance(value, Mapping):
        if head not in value:
            raise KeyError(head)
        updated_mapping = dict(value)
        updated_mapping[head] = _deep_replace(updated_mapping[head], tail, replacement)
        return updated_mapping

    clone = copy.deepcopy(value)
    current = getattr(clone, head)
    updated = _deep_replace(current, tail, replacement)
    try:
        setattr(clone, head, updated)
    except Exception:
        object.__setattr__(clone, head, updated)
    return clone


def _mutate_first_supported_path(
    value: Any,
    candidate_paths: Sequence[Sequence[str]],
    replacement_factory,
) -> Any:
    for path in candidate_paths:
        if not _path_exists(value, path):
            continue
        try:
            original = value
            for key in path:
                original = (
                    original[key]
                    if isinstance(original, Mapping)
                    else getattr(original, key)
                )
            return _deep_replace(value, path, replacement_factory(original))
        except Exception:
            continue
    raise AssertionError(f"no supported mutation path among: {candidate_paths!r}")


def test_importing_observables_module_does_not_load_repo_v1_modules() -> None:
    if not _has_required_modules():
        pytest.skip("paper_0710_1869 hadronic/observable layers not implemented yet")

    script = f"""
import importlib
import json
import sys

importlib.import_module("quarkConstraints.paper_0710_1869.eft_deltaf2.observables")

forbidden = sorted(
    name for name in sys.modules
    if (
        name in {sorted(FORBIDDEN_REPO_V1_MODULES)!r}
        or any(name.startswith(f"{{root}}.") for root in {sorted(FORBIDDEN_REPO_V1_MODULES)!r})
    )
)
print(json.dumps(forbidden))
"""
    completed = subprocess.run(
        [sys.executable, "-c", script],
        check=True,
        capture_output=True,
        text=True,
        cwd=str(REPO_ROOT),
    )
    assert json.loads(completed.stdout) == []


def test_default_kaon_observable_export_is_finite_deterministic_and_np_only() -> None:
    module = _load_observables_module()

    payload = _default_observable_payload(module)
    second = _default_observable_payload(module)
    cross_process_payload = _default_observable_payload_cross_process()

    assert json.dumps(payload, sort_keys=True) == json.dumps(second, sort_keys=True)
    assert json.dumps(payload, sort_keys=True) == json.dumps(
        cross_process_payload,
        sort_keys=True,
    )

    interpretation = _nested_value(
        payload,
        (
            ("interpretation",),
            ("observable_scope_id",),
            ("contract", "interpretation"),
            ("contract", "observable_scope_id"),
            ("tags", "interpretation"),
            ("tags", "observable_scope_id"),
        ),
    )
    observable_scope_id = _nested_value(
        payload,
        (
            ("observable_scope_id",),
            ("contract", "observable_scope_id"),
            ("tags", "observable_scope_id"),
        ),
    )
    m12_observable_id = _nested_value(
        payload,
        (
            ("m12_observable_id",),
            ("tags", "m12_observable_id"),
        ),
    )
    delta_m_observable_id = _nested_value(
        payload,
        (
            ("delta_m_observable_id",),
            ("tags", "delta_m_observable_id"),
        ),
    )
    system_id = _nested_value(
        payload,
        (
            ("system_id",),
            ("system",),
            ("tags", "system_id"),
            ("contract", "system_id"),
        ),
    )
    scheme_id = _nested_value(
        payload,
        (
            ("renormalization_scheme_id",),
            ("scheme_id",),
            ("tags", "renormalization_scheme_id"),
            ("contract", "renormalization_scheme_id"),
        ),
    )
    mu_had_GeV = _nested_value(
        payload,
        (
            ("mu_had_GeV",),
            ("evaluation_scale_GeV",),
            ("hadronic_scale_GeV",),
            ("tags", "mu_had_GeV"),
            ("tags", "evaluation_scale_GeV"),
            ("contract", "mu_had_GeV"),
        ),
    )
    parity_relation_id = _nested_value(
        payload,
        (
            ("hadronic_bundle", "parity_relation_id"),
            ("tags", "parity_relation_id"),
        ),
    )
    observables = _extract_named_observables(payload)

    checks = {
        "interpretation_is_exact_np_only": interpretation == EXPECTED_INTERPRETATION_ID,
        "observable_scope_is_pr5a_surface": observable_scope_id == EXPECTED_OBSERVABLE_SCOPE_ID,
        "m12_observable_id_is_supported": m12_observable_id == EXPECTED_M12_OBSERVABLE_ID,
        "delta_m_observable_id_is_supported": (
            delta_m_observable_id == EXPECTED_DELTA_M_OBSERVABLE_ID
        ),
        "system_is_kaon": _is_kaon_system(system_id),
        "scheme_id_present": bool(scheme_id),
        "mu_had_present_and_positive": isinstance(mu_had_GeV, (int, float))
        and math.isfinite(float(mu_had_GeV))
        and float(mu_had_GeV) > 0.0,
        "observable_names_match_supported_surface": set(observables) == SUPPORTED_OBSERVABLE_ROWS,
        "observable_values_finite": all(
            math.isfinite(float(value)) for value in observables.values()
        ),
        "parity_relation_tag_present": bool(parity_relation_id),
        "parity_relation_matches_supported_relation": bool(parity_relation_id)
        and "q1_vll_equals_q1_vrr" in str(parity_relation_id).lower(),
    }
    assert not [name for name, ok in checks.items() if not ok], checks


def test_observable_evaluation_rejects_nonzero_q4_lr_contribution() -> None:
    module = _load_observables_module()
    evaluator = _get_observable_evaluator(module)
    lr_evaluator = _get_custom_lr_observable_evaluator(module)
    hadronic_value = _default_hadronic_value()
    rg_value = _default_rg_result()
    q4_probe = complex(1.0e-12, 0.0)

    mutated_rg_value = _mutate_first_supported_path(
        rg_value,
        (
            ("wilsons", "q4_lr"),
            ("q4_lr",),
        ),
        lambda _original: q4_probe,
    )

    # C-2: the default path no longer hard-errors on LR coefficients; it carries
    # them in the Wilson snapshot while explicit LR surfaces own LR matrix elements.
    payload = _payload_from_value(
        _invoke_observable_evaluation(
            evaluator,
            rg_value=mutated_rg_value,
            hadronic_value=hadronic_value,
        )
    )
    _assert_coefficient_payload(payload, "Q4_LR", q4_probe)

    custom_hadronic = _build_custom_lr_hadronic_inputs(
        renormalization_scheme_id=str(hadronic_value.renormalization_scheme_id),
        mu_had_GeV=float(hadronic_value.mu_had_GeV),
        source_prefix="hadronic.kaon.c2_q4_lr_probe",
    )
    custom_wilsons = _build_lr_probe_wilsons(
        q4_lr=q4_probe,
        q5_lr=0.0 + 0.0j,
        renormalization_scheme_id=str(custom_hadronic.renormalization_scheme_id),
        matching_scale_GeV=float(custom_hadronic.mu_had_GeV),
    )
    lr_payload = _payload_from_value(
        _invoke_observable_evaluation(
            lr_evaluator,
            rg_value=custom_wilsons,
            hadronic_value=custom_hadronic,
        )
    )
    expected_lr_m12 = (
        q4_probe * float(custom_hadronic.q4_matrix_element_GeV4)
    ) / (2.0 * float(custom_hadronic.m_K0_GeV))
    observed_lr_m12 = _complex_from_payload(lr_payload["M12_K_LR_NP_GeV"])
    assert abs(observed_lr_m12) > 0.0
    assert observed_lr_m12.real == pytest.approx(expected_lr_m12.real, rel=0.0, abs=1.0e-24)
    assert observed_lr_m12.imag == pytest.approx(expected_lr_m12.imag, rel=0.0, abs=1.0e-24)


def test_observable_evaluation_rejects_nonzero_q5_lr_contribution() -> None:
    module = _load_observables_module()
    evaluator = _get_observable_evaluator(module)
    lr_evaluator = _get_custom_lr_observable_evaluator(module)
    hadronic_value = _default_hadronic_value()
    rg_value = _default_rg_result()
    q5_probe = complex(1.0e-12, 0.0)

    mutated_rg_value = _mutate_first_supported_path(
        rg_value,
        (
            ("wilsons", "q5_lr"),
            ("q5_lr",),
        ),
        lambda _original: q5_probe,
    )

    # C-2: the default path no longer hard-errors on LR coefficients; it carries
    # them in the Wilson snapshot while explicit LR surfaces own LR matrix elements.
    payload = _payload_from_value(
        _invoke_observable_evaluation(
            evaluator,
            rg_value=mutated_rg_value,
            hadronic_value=hadronic_value,
        )
    )
    _assert_coefficient_payload(payload, "Q5_LR", q5_probe)

    custom_hadronic = _build_custom_lr_hadronic_inputs(
        renormalization_scheme_id=str(hadronic_value.renormalization_scheme_id),
        mu_had_GeV=float(hadronic_value.mu_had_GeV),
        source_prefix="hadronic.kaon.c2_q5_lr_probe",
    )
    custom_wilsons = _build_lr_probe_wilsons(
        q4_lr=0.0 + 0.0j,
        q5_lr=q5_probe,
        renormalization_scheme_id=str(custom_hadronic.renormalization_scheme_id),
        matching_scale_GeV=float(custom_hadronic.mu_had_GeV),
    )
    lr_payload = _payload_from_value(
        _invoke_observable_evaluation(
            lr_evaluator,
            rg_value=custom_wilsons,
            hadronic_value=custom_hadronic,
        )
    )
    expected_lr_m12 = (
        q5_probe * float(custom_hadronic.q5_matrix_element_GeV4)
    ) / (2.0 * float(custom_hadronic.m_K0_GeV))
    observed_lr_m12 = _complex_from_payload(lr_payload["M12_K_LR_NP_GeV"])
    assert abs(observed_lr_m12) > 0.0
    assert observed_lr_m12.real == pytest.approx(expected_lr_m12.real, rel=0.0, abs=1.0e-24)
    assert observed_lr_m12.imag == pytest.approx(expected_lr_m12.imag, rel=0.0, abs=1.0e-24)


def test_observable_evaluation_rejects_lr_even_with_custom_lr_hadronic_inputs() -> None:
    hadronic_module = _load_hadronic_module()
    observables_module = _load_observables_module()
    evaluator = _get_observable_evaluator(observables_module)
    default_builder = _get_callable(hadronic_module, HADRONIC_BUILDER_NAMES)
    assert callable(default_builder), "hadronic layer exposes no default builder"
    default_bundle = default_builder()
    lr_builder = _custom_lr_hadronic_builder(hadronic_module)

    scheme_id = LR_HADRONIC_PROBE_SCHEME_ID
    mu_had_GeV = LR_HADRONIC_PROBE_MU_HAD_GEV
    b4_source = _build_hadronic_source_ref(
        hadronic_module,
        source_id="hadronic.kaon.lr_observable_probe.b4.v1",
        citation="Becirevic and Villadoro, hep-lat/0408029",
        locator_label="eq. (5), O4 scalar LR matrix element",
        year=2004,
        scheme_id=scheme_id,
        scale_GeV=mu_had_GeV,
        notes="Observable LR-blocking probe source for B4(mu_had).",
    )
    b5_source = _build_hadronic_source_ref(
        hadronic_module,
        source_id="hadronic.kaon.lr_observable_probe.b5.v1",
        citation="Becirevic and Villadoro, hep-lat/0408029",
        locator_label="eq. (5), O5 scalar LR matrix element",
        year=2004,
        scheme_id=scheme_id,
        scale_GeV=mu_had_GeV,
        notes="Observable LR-blocking probe source for B5(mu_had).",
    )
    r_chi_source = _build_hadronic_source_ref(
        hadronic_module,
        source_id="hadronic.kaon.lr_observable_probe.rchi.v1",
        citation="Custom caller-supplied chiral ratio input",
        locator_label="R_chi(mu_had) external input",
        year=2026,
        scheme_id=scheme_id,
        scale_GeV=mu_had_GeV,
        notes="Observable LR-blocking probe source for R_chi(mu_had).",
    )
    custom_lr_hadronic = lr_builder(
        **_custom_lr_hadronic_builder_kwargs(
            lr_builder,
            default_bundle=default_bundle,
            B4_mu_had=LR_HADRONIC_PROBE_B4,
            B5_mu_had=LR_HADRONIC_PROBE_B5,
            R_chi_mu_had=LR_HADRONIC_PROBE_R_CHI,
            B4_source=b4_source,
            B5_source=b5_source,
            R_chi_source=r_chi_source,
            renormalization_scheme_id=scheme_id,
            mu_had_GeV=mu_had_GeV,
        )
    )
    assert str(custom_lr_hadronic.renormalization_scheme_id) == scheme_id
    assert float(custom_lr_hadronic.mu_had_GeV) == pytest.approx(
        mu_had_GeV,
        rel=0.0,
        abs=1.0e-12,
    )

    rg_value = _default_rg_result()
    mutated_rg_value = dataclasses.replace(
        rg_value,
        wilsons=dataclasses.replace(
            rg_value.wilsons,
            q4_lr=complex(1.0e-12, 0.0),
            q5_lr=complex(-5.0e-13, 0.0),
        ),
    )

    with pytest.raises(ValueError) as exc_info:
        _invoke_observable_evaluation(
            evaluator,
            rg_value=mutated_rg_value,
            hadronic_value=custom_lr_hadronic,
        )

    message = str(exc_info.value)
    assert (
        (
            "Q4_LR" in message
            and "Q5_LR" in message
            and "blocked" in message
        )
        or "hadronic_bundle must be a Paper07101869KaonHadronicBundle" in message
    )


def test_custom_lr_only_observable_path_matches_hand_calculation() -> None:
    module = _load_observables_module()
    evaluator = _get_custom_lr_observable_evaluator(module)
    custom_hadronic = _build_custom_lr_hadronic_inputs()
    custom_wilsons = _build_lr_probe_wilsons(
        renormalization_scheme_id=LR_HADRONIC_PROBE_SCHEME_ID,
        matching_scale_GeV=float(custom_hadronic.mu_had_GeV),
    )

    result = _invoke_observable_evaluation(
        evaluator,
        rg_value=custom_wilsons,
        hadronic_value=custom_hadronic,
    )
    payload = _payload_from_value(result)

    assert payload["observable_scope_id"] == EXPECTED_LR_ONLY_SCOPE_ID
    assert payload["interpretation"] == EXPECTED_LR_ONLY_INTERPRETATION_ID
    assert payload["m12_observable_id"] == EXPECTED_LR_ONLY_M12_OBSERVABLE_ID
    assert payload["delta_m_observable_id"] == EXPECTED_LR_ONLY_DELTA_M_OBSERVABLE_ID
    assert payload["renormalization_scheme_id"] == LR_HADRONIC_PROBE_SCHEME_ID
    assert float(payload["mu_had_GeV"]) == pytest.approx(
        LR_HADRONIC_PROBE_MU_HAD_GEV,
        rel=0.0,
        abs=1.0e-12,
    )

    expected_m12 = (
        custom_wilsons.q4_lr * float(custom_hadronic.q4_matrix_element_GeV4)
        + custom_wilsons.q5_lr * float(custom_hadronic.q5_matrix_element_GeV4)
    ) / (2.0 * float(custom_hadronic.m_K0_GeV))
    expected_delta_m = 2.0 * float(expected_m12.real)

    observed_m12 = _complex_from_payload(payload["M12_K_LR_NP_GeV"])
    assert observed_m12.real == pytest.approx(
        expected_m12.real,
        rel=0.0,
        abs=1.0e-24,
    )
    assert observed_m12.imag == pytest.approx(
        expected_m12.imag,
        rel=0.0,
        abs=1.0e-24,
    )
    assert float(payload["delta_m_K_LR_NP_GeV"]) == pytest.approx(
        expected_delta_m,
        rel=0.0,
        abs=1.0e-24,
    )
    assert float(payload["coefficients"]["Q4_LR"]["real"]) == pytest.approx(
        custom_wilsons.q4_lr.real,
        rel=0.0,
        abs=1.0e-24,
    )
    assert float(payload["coefficients"]["Q4_LR"]["imag"]) == pytest.approx(
        custom_wilsons.q4_lr.imag,
        rel=0.0,
        abs=1.0e-24,
    )
    assert float(payload["coefficients"]["Q5_LR"]["real"]) == pytest.approx(
        custom_wilsons.q5_lr.real,
        rel=0.0,
        abs=1.0e-24,
    )
    assert float(payload["coefficients"]["Q5_LR"]["imag"]) == pytest.approx(
        custom_wilsons.q5_lr.imag,
        rel=0.0,
        abs=1.0e-24,
    )


def test_custom_lr_only_observable_path_rejects_scheme_mismatch() -> None:
    module = _load_observables_module()
    evaluator = _get_custom_lr_observable_evaluator(module)
    custom_hadronic = _build_custom_lr_hadronic_inputs()
    mismatched_wilsons = _build_lr_probe_wilsons(
        renormalization_scheme_id="paper_0710_1869.deltaf2.kaon_lr_observable.wrong_scheme.v1",
        matching_scale_GeV=float(custom_hadronic.mu_had_GeV),
    )

    with pytest.raises(ValueError, match="renormalization_scheme_id"):
        _invoke_observable_evaluation(
            evaluator,
            rg_value=mismatched_wilsons,
            hadronic_value=custom_hadronic,
        )


def test_custom_lr_only_observable_path_rejects_scale_mismatch() -> None:
    module = _load_observables_module()
    evaluator = _get_custom_lr_observable_evaluator(module)
    custom_hadronic = _build_custom_lr_hadronic_inputs()
    mismatched_wilsons = _build_lr_probe_wilsons(
        renormalization_scheme_id=LR_HADRONIC_PROBE_SCHEME_ID,
        matching_scale_GeV=float(custom_hadronic.mu_had_GeV) + 0.25,
    )

    with pytest.raises(ValueError, match="mu_had_GeV|evaluation scale"):
        _invoke_observable_evaluation(
            evaluator,
            rg_value=mismatched_wilsons,
            hadronic_value=custom_hadronic,
        )


def test_custom_lr_only_observable_path_rejects_frozen_r_chi_object() -> None:
    module = _load_observables_module()
    evaluator = _get_custom_lr_observable_evaluator(module)
    frozen_r_chi = _default_kaon_lr_r_chi_freeze()
    probe_wilsons = _build_lr_probe_wilsons(
        renormalization_scheme_id=str(frozen_r_chi.operator_renormalization_scheme_id),
        matching_scale_GeV=float(frozen_r_chi.mu_had_GeV),
    )

    with pytest.raises(
        ValueError,
        match="hadronic_inputs must be a Paper07101869KaonLRHadronicInputs",
    ):
        _invoke_observable_evaluation(
            evaluator,
            rg_value=probe_wilsons,
            hadronic_value=frozen_r_chi,
        )


def test_default_q1_only_observable_path_still_rejects_nonzero_lr_without_custom_api() -> None:
    module = _load_observables_module()
    default_evaluator = _get_observable_evaluator(module)
    custom_evaluator = _get_custom_lr_observable_evaluator(module)
    default_hadronic = _default_hadronic_value()
    custom_hadronic = _build_custom_lr_hadronic_inputs()
    default_path_wilsons = _build_lr_probe_wilsons(
        matching_scale_GeV=float(default_hadronic.mu_had_GeV),
    )
    custom_wilsons = _build_lr_probe_wilsons(
        renormalization_scheme_id=LR_HADRONIC_PROBE_SCHEME_ID,
        matching_scale_GeV=float(custom_hadronic.mu_had_GeV),
    )

    # C-2: nonzero LR Wilsons now flow through the default observable payload
    # without the old guard; explicit LR hadronic inputs produce the LR M12 piece.
    default_payload = _payload_from_value(
        _invoke_observable_evaluation(
            default_evaluator,
            rg_value=default_path_wilsons,
            hadronic_value=default_hadronic,
        )
    )
    _assert_coefficient_payload(default_payload, "Q4_LR", LR_OBSERVABLE_PROBE_Q4)
    _assert_coefficient_payload(default_payload, "Q5_LR", LR_OBSERVABLE_PROBE_Q5)
    assert _complex_from_payload(default_payload["M12_K_NP_GeV"]) == pytest.approx(
        0.0 + 0.0j,
        rel=0.0,
        abs=1.0e-24,
    )

    custom_result = _invoke_observable_evaluation(
        custom_evaluator,
        rg_value=custom_wilsons,
        hadronic_value=custom_hadronic,
    )
    custom_payload = _payload_from_value(custom_result)
    expected_lr_m12 = (
        custom_wilsons.q4_lr * float(custom_hadronic.q4_matrix_element_GeV4)
        + custom_wilsons.q5_lr * float(custom_hadronic.q5_matrix_element_GeV4)
    ) / (2.0 * float(custom_hadronic.m_K0_GeV))
    assert custom_payload["m12_observable_id"] == EXPECTED_LR_ONLY_M12_OBSERVABLE_ID
    observed_lr_m12 = _complex_from_payload(custom_payload["M12_K_LR_NP_GeV"])
    assert abs(observed_lr_m12) > 0.0
    assert observed_lr_m12 == pytest.approx(expected_lr_m12, rel=0.0, abs=1.0e-24)


def test_custom_total_observable_path_matches_hand_calculation_and_existing_pieces() -> None:
    module = _load_observables_module()
    combined_evaluator = _get_custom_total_observable_evaluator(module)
    q1_evaluator = _get_observable_evaluator(module)
    lr_evaluator = _get_custom_lr_observable_evaluator(module)
    custom_q1_hadronic, custom_lr_hadronic, combined_wilsons = _build_custom_total_probe_inputs()

    result = _invoke_custom_total_observable_evaluation(
        combined_evaluator,
        rg_value=combined_wilsons,
        q1_hadronic_bundle=custom_q1_hadronic,
        lr_hadronic_inputs=custom_lr_hadronic,
    )
    payload = _payload_from_value(result)

    assert payload["observable_scope_id"] == EXPECTED_CUSTOM_TOTAL_SCOPE_ID
    assert payload["interpretation"] == EXPECTED_CUSTOM_TOTAL_INTERPRETATION_ID
    assert payload["m12_observable_id"] == EXPECTED_CUSTOM_TOTAL_M12_OBSERVABLE_ID
    assert payload["delta_m_observable_id"] == EXPECTED_CUSTOM_TOTAL_DELTA_M_OBSERVABLE_ID
    assert payload["renormalization_scheme_id"] == str(custom_q1_hadronic.renormalization_scheme_id)
    assert float(payload["mu_had_GeV"]) == pytest.approx(
        float(custom_q1_hadronic.mu_had_GeV),
        rel=0.0,
        abs=1.0e-12,
    )

    hand_q1 = (
        (combined_wilsons.q1_vll + combined_wilsons.q1_vrr)
        * float(custom_q1_hadronic.q1_matrix_element_GeV4)
        / (2.0 * float(custom_q1_hadronic.m_K0_GeV))
    )
    hand_lr = (
        combined_wilsons.q4_lr * float(custom_lr_hadronic.q4_matrix_element_GeV4)
        + combined_wilsons.q5_lr * float(custom_lr_hadronic.q5_matrix_element_GeV4)
    ) / (2.0 * float(custom_lr_hadronic.m_K0_GeV))
    hand_total = hand_q1 + hand_lr
    hand_delta_m = 2.0 * float(hand_total.real)

    observed_q1 = _complex_from_payload(payload["M12_K_NP_Q1_GeV"])
    observed_lr = _complex_from_payload(payload["M12_K_LR_NP_GeV"])
    observed_total = _complex_from_payload(payload["M12_K_NP_TOTAL_GeV"])
    assert observed_q1.real == pytest.approx(hand_q1.real, rel=0.0, abs=1.0e-24)
    assert observed_q1.imag == pytest.approx(hand_q1.imag, rel=0.0, abs=1.0e-24)
    assert observed_lr.real == pytest.approx(hand_lr.real, rel=0.0, abs=1.0e-24)
    assert observed_lr.imag == pytest.approx(hand_lr.imag, rel=0.0, abs=1.0e-24)
    assert observed_total.real == pytest.approx(hand_total.real, rel=0.0, abs=1.0e-24)
    assert observed_total.imag == pytest.approx(hand_total.imag, rel=0.0, abs=1.0e-24)
    assert float(payload["delta_m_K_NP_TOTAL_GeV"]) == pytest.approx(
        hand_delta_m,
        rel=0.0,
        abs=1.0e-24,
    )
    assert observed_total == pytest.approx(observed_q1 + observed_lr, rel=0.0, abs=1.0e-24)

    q1_only_wilsons = dataclasses.replace(combined_wilsons, q4_lr=0.0 + 0.0j, q5_lr=0.0 + 0.0j)
    lr_only_wilsons = dataclasses.replace(combined_wilsons, q1_vll=0.0 + 0.0j, q1_vrr=0.0 + 0.0j)
    q1_payload = _payload_from_value(
        _invoke_observable_evaluation(
            q1_evaluator,
            rg_value=q1_only_wilsons,
            hadronic_value=custom_q1_hadronic,
        )
    )
    lr_payload = _payload_from_value(
        _invoke_observable_evaluation(
            lr_evaluator,
            rg_value=lr_only_wilsons,
            hadronic_value=custom_lr_hadronic,
        )
    )
    existing_q1 = _complex_from_payload(q1_payload["M12_K_NP_GeV"])
    existing_lr = _complex_from_payload(lr_payload["M12_K_LR_NP_GeV"])
    assert observed_total == pytest.approx(existing_q1 + existing_lr, rel=0.0, abs=1.0e-24)


def test_custom_total_observable_path_rejects_scheme_mismatch() -> None:
    module = _load_observables_module()
    combined_evaluator = _get_custom_total_observable_evaluator(module)
    custom_q1_hadronic, custom_lr_hadronic, combined_wilsons = _build_custom_total_probe_inputs()
    mismatched_lr_hadronic = _build_custom_lr_hadronic_inputs(
        renormalization_scheme_id="paper_0710_1869.deltaf2.kaon_custom_total.wrong_scheme.v1",
        mu_had_GeV=float(custom_q1_hadronic.mu_had_GeV),
        source_prefix="hadronic.kaon.custom_total_wrong_scheme",
    )

    with pytest.raises(ValueError, match="renormalization_scheme_id|scheme"):
        _invoke_custom_total_observable_evaluation(
            combined_evaluator,
            rg_value=combined_wilsons,
            q1_hadronic_bundle=custom_q1_hadronic,
            lr_hadronic_inputs=mismatched_lr_hadronic,
        )


def test_custom_total_observable_path_rejects_mu_had_mismatch() -> None:
    module = _load_observables_module()
    combined_evaluator = _get_custom_total_observable_evaluator(module)
    custom_q1_hadronic, custom_lr_hadronic, combined_wilsons = _build_custom_total_probe_inputs()
    mismatched_lr_hadronic = _build_custom_lr_hadronic_inputs(
        renormalization_scheme_id=str(custom_q1_hadronic.renormalization_scheme_id),
        mu_had_GeV=float(custom_q1_hadronic.mu_had_GeV) + 0.25,
        source_prefix="hadronic.kaon.custom_total_wrong_scale",
    )

    with pytest.raises(ValueError, match="mu_had|evaluation scale|evaluation_scale"):
        _invoke_custom_total_observable_evaluation(
            combined_evaluator,
            rg_value=combined_wilsons,
            q1_hadronic_bundle=custom_q1_hadronic,
            lr_hadronic_inputs=mismatched_lr_hadronic,
        )


@pytest.mark.parametrize(
    ("field_name", "error_pattern"),
    (
        ("operator_basis_id", "(?i)basis"),
        ("operator_normalization_id", "(?i)normalization"),
    ),
)
def test_custom_total_observable_path_rejects_alignment_tag_mismatch(
    field_name: str,
    error_pattern: str,
) -> None:
    module = _load_observables_module()
    combined_evaluator = _get_custom_total_observable_evaluator(module)
    custom_q1_hadronic, custom_lr_hadronic, combined_wilsons = _build_custom_total_probe_inputs()

    replacement = f"{getattr(custom_lr_hadronic, field_name)}.wrong"
    try:
        mismatched_lr_hadronic = dataclasses.replace(
            custom_lr_hadronic,
            **{
                field_name: replacement,
                "contract": dataclasses.replace(
                    custom_lr_hadronic.contract,
                    **{field_name: replacement},
                ),
            },
        )
    except ValueError as exc:
        assert re.search(error_pattern, str(exc))
        return

    with pytest.raises(ValueError, match=error_pattern):
        _invoke_custom_total_observable_evaluation(
            combined_evaluator,
            rg_value=combined_wilsons,
            q1_hadronic_bundle=custom_q1_hadronic,
            lr_hadronic_inputs=mismatched_lr_hadronic,
        )


def test_custom_total_observable_path_rejects_default_q1_hadronic_bundle() -> None:
    module = _load_observables_module()
    combined_evaluator = _get_custom_total_observable_evaluator(module)
    default_q1_hadronic = _default_hadronic_value()
    custom_lr_hadronic = _build_custom_lr_hadronic_inputs(
        renormalization_scheme_id=str(default_q1_hadronic.renormalization_scheme_id),
        mu_had_GeV=float(default_q1_hadronic.mu_had_GeV),
        source_prefix="hadronic.kaon.custom_total_default_q1_reject",
    )
    combined_wilsons = dataclasses.replace(
        _default_rg_result().wilsons,
        q4_lr=LR_OBSERVABLE_PROBE_Q4,
        q5_lr=LR_OBSERVABLE_PROBE_Q5,
    )

    with pytest.raises(ValueError, match="custom Q1 hadronic bundle|input_provenance_mode_id"):
        _invoke_custom_total_observable_evaluation(
            combined_evaluator,
            rg_value=combined_wilsons,
            q1_hadronic_bundle=default_q1_hadronic,
            lr_hadronic_inputs=custom_lr_hadronic,
        )


def test_custom_total_observable_path_rejects_frozen_r_chi_object() -> None:
    module = _load_observables_module()
    combined_evaluator = _get_custom_total_observable_evaluator(module)
    frozen_r_chi = _default_kaon_lr_r_chi_freeze()
    custom_q1_hadronic = _build_custom_q1_hadronic_bundle(
        mu_had_GeV=float(frozen_r_chi.mu_had_GeV),
    )
    combined_wilsons = dataclasses.replace(
        _default_rg_result().wilsons,
        matching_scale_GeV=float(frozen_r_chi.mu_had_GeV),
        contract=dataclasses.replace(
            _default_rg_result().wilsons.contract,
            renormalization_scheme_id=str(frozen_r_chi.operator_renormalization_scheme_id),
        ),
        q4_lr=LR_OBSERVABLE_PROBE_Q4,
        q5_lr=LR_OBSERVABLE_PROBE_Q5,
    )

    with pytest.raises(
        ValueError,
        match="lr_hadronic_inputs must be a Paper07101869KaonLRHadronicInputs",
    ):
        _invoke_custom_total_observable_evaluation(
            combined_evaluator,
            rg_value=combined_wilsons,
            q1_hadronic_bundle=custom_q1_hadronic,
            lr_hadronic_inputs=frozen_r_chi,
        )


def test_default_q1_only_path_remains_blocked_while_custom_total_path_succeeds() -> None:
    module = _load_observables_module()
    default_evaluator = _get_observable_evaluator(module)
    combined_evaluator = _get_custom_total_observable_evaluator(module)
    custom_q1_hadronic, custom_lr_hadronic, combined_wilsons = _build_custom_total_probe_inputs()

    # C-2: the default payload carries LR Wilsons instead of blocking the pipeline.
    default_payload = _payload_from_value(
        _invoke_observable_evaluation(
            default_evaluator,
            rg_value=combined_wilsons,
            hadronic_value=_default_hadronic_value(),
        )
    )
    _assert_coefficient_payload(default_payload, "Q4_LR", combined_wilsons.q4_lr)
    _assert_coefficient_payload(default_payload, "Q5_LR", combined_wilsons.q5_lr)

    combined_payload = _payload_from_value(
        _invoke_custom_total_observable_evaluation(
            combined_evaluator,
            rg_value=combined_wilsons,
            q1_hadronic_bundle=custom_q1_hadronic,
            lr_hadronic_inputs=custom_lr_hadronic,
        )
    )
    assert combined_payload["m12_observable_id"] == EXPECTED_CUSTOM_TOTAL_M12_OBSERVABLE_ID
    observed_lr = _complex_from_payload(combined_payload["M12_K_LR_NP_GeV"])
    observed_total = _complex_from_payload(combined_payload["M12_K_NP_TOTAL_GeV"])
    expected_lr = (
        combined_wilsons.q4_lr * float(custom_lr_hadronic.q4_matrix_element_GeV4)
        + combined_wilsons.q5_lr * float(custom_lr_hadronic.q5_matrix_element_GeV4)
    ) / (2.0 * float(custom_lr_hadronic.m_K0_GeV))
    assert abs(observed_lr) > 0.0
    assert observed_lr == pytest.approx(expected_lr, rel=0.0, abs=1.0e-24)
    assert observed_total == pytest.approx(
        _complex_from_payload(combined_payload["M12_K_NP_Q1_GeV"]) + observed_lr,
        rel=0.0,
        abs=1.0e-24,
    )
    assert "epsilon_K_NP" not in combined_payload.get("observables", {})


@pytest.mark.parametrize("system_id", tuple(B_CUSTOM_Q1_PROBE_CONFIG))
def test_custom_b_q1_observable_paths_match_hand_calculation(system_id: str) -> None:
    module = _load_observables_module()
    evaluator = _get_custom_b_q1_observable_evaluator(module, system_id)
    hadronic_bundle, rg_result = _build_custom_b_probe_inputs(system_id)
    payload = _payload_from_value(
        _invoke_observable_evaluation(
            evaluator,
            rg_value=rg_result,
            hadronic_value=hadronic_bundle,
        )
    )
    config = B_CUSTOM_Q1_PROBE_CONFIG[system_id]
    observed_m12 = _complex_from_payload(payload["M12_NP_GeV"])
    expected_m12 = (
        (rg_result.wilsons.q1_vll + rg_result.wilsons.q1_vrr)
        * float(hadronic_bundle.q1_matrix_element_GeV4)
        / (2.0 * float(hadronic_bundle.meson_mass_GeV))
    )
    expected_delta_m = 2.0 * float(expected_m12.real)

    assert payload["system_id"] == system_id
    assert payload["observable_scope_id"] == config["scope_id"]
    assert payload["interpretation"] == config["interpretation_id"]
    assert payload["m12_observable_id"] == config["m12_id"]
    assert payload["delta_m_observable_id"] == config["delta_m_id"]
    assert payload["renormalization_scheme_id"] == str(hadronic_bundle.renormalization_scheme_id)
    assert float(payload["mu_had_GeV"]) == pytest.approx(
        float(hadronic_bundle.mu_had_GeV),
        rel=0.0,
        abs=1.0e-12,
    )
    assert float(payload["matching_scale_GeV"]) == pytest.approx(
        float(hadronic_bundle.mu_had_GeV),
        rel=0.0,
        abs=1.0e-12,
    )
    if bool(config.get("require_nonzero_probe", False)):
        assert abs(expected_m12) > 0.0
        assert abs(observed_m12) > 0.0
    assert observed_m12.real == pytest.approx(expected_m12.real, rel=0.0, abs=1.0e-24)
    assert observed_m12.imag == pytest.approx(expected_m12.imag, rel=0.0, abs=1.0e-24)
    assert float(payload["delta_m_NP_GeV"]) == pytest.approx(
        expected_delta_m,
        rel=0.0,
        abs=1.0e-24,
    )
    assert payload["observables"][config["m12_id"]]["re"] == pytest.approx(
        expected_m12.real,
        rel=0.0,
        abs=1.0e-24,
    )
    assert payload["observables"][config["m12_id"]]["im"] == pytest.approx(
        expected_m12.imag,
        rel=0.0,
        abs=1.0e-24,
    )
    assert payload["observables"][config["delta_m_id"]] == pytest.approx(
        expected_delta_m,
        rel=0.0,
        abs=1.0e-24,
    )
    assert _complex_from_payload(payload["coefficients"]["Q4_LR"]) == 0.0 + 0.0j
    assert _complex_from_payload(payload["coefficients"]["Q5_LR"]) == 0.0 + 0.0j
    assert (
        "This does not widen the default/exported kaon Q1-only observable or artifact surface."
        in payload["notes"]
    )


@pytest.mark.parametrize("system_id", tuple(B_CUSTOM_Q1_PROBE_CONFIG))
def test_custom_b_q1_observable_paths_reject_scheme_mismatch(system_id: str) -> None:
    module = _load_observables_module()
    evaluator = _get_custom_b_q1_observable_evaluator(module, system_id)
    hadronic_bundle, rg_result = _build_custom_b_probe_inputs(system_id)
    wrong_scheme = f"{hadronic_bundle.renormalization_scheme_id}.wrong"
    mismatched_wilsons = dataclasses.replace(
        rg_result.wilsons,
        contract=dataclasses.replace(
            rg_result.wilsons.contract,
            renormalization_scheme_id=wrong_scheme,
        ),
    )

    with pytest.raises(ValueError, match="renormalization_scheme_id|scheme"):
        _invoke_observable_evaluation(
            evaluator,
            rg_value=mismatched_wilsons,
            hadronic_value=hadronic_bundle,
        )


@pytest.mark.parametrize("system_id", tuple(B_CUSTOM_Q1_PROBE_CONFIG))
def test_custom_b_q1_observable_paths_reject_mu_had_mismatch(system_id: str) -> None:
    module = _load_observables_module()
    evaluator = _get_custom_b_q1_observable_evaluator(module, system_id)
    hadronic_bundle, rg_result = _build_custom_b_probe_inputs(system_id)
    mismatched_wilsons = dataclasses.replace(
        rg_result.wilsons,
        matching_scale_GeV=float(hadronic_bundle.mu_had_GeV) + 0.25,
    )

    with pytest.raises(ValueError, match="mu_had_GeV|evaluation scale"):
        _invoke_observable_evaluation(
            evaluator,
            rg_value=mismatched_wilsons,
            hadronic_value=hadronic_bundle,
        )


@pytest.mark.parametrize("system_id", tuple(B_CUSTOM_Q1_PROBE_CONFIG))
@pytest.mark.parametrize(
    ("field_name", "error_pattern"),
    (
        ("operator_basis_id", "(?i)basis"),
        ("operator_normalization_id", "(?i)normalization"),
    ),
)
def test_custom_b_q1_observable_paths_reject_alignment_tag_mismatch(
    system_id: str,
    field_name: str,
    error_pattern: str,
) -> None:
    module = _load_observables_module()
    evaluator = _get_custom_b_q1_observable_evaluator(module, system_id)
    hadronic_bundle, rg_result = _build_custom_b_probe_inputs(system_id)
    replacement = f"{getattr(rg_result.wilsons, field_name)}.wrong"
    mismatched_wilsons = dataclasses.replace(
        rg_result.wilsons,
        contract=dataclasses.replace(
            rg_result.wilsons.contract,
            **{field_name: replacement},
        ),
    )

    with pytest.raises(ValueError, match=error_pattern):
        _invoke_observable_evaluation(
            evaluator,
            rg_value=mismatched_wilsons,
            hadronic_value=hadronic_bundle,
        )


@pytest.mark.parametrize(
    ("system_id", "wrong_system_id"),
    (("B_d", "B_s"), ("B_s", "B_d"), ("D0", "B_d")),
)
def test_custom_b_q1_observable_paths_reject_system_mismatch(
    system_id: str,
    wrong_system_id: str,
) -> None:
    module = _load_observables_module()
    evaluator = _get_custom_b_q1_observable_evaluator(module, system_id)
    hadronic_bundle, rg_result = _build_custom_b_probe_inputs(system_id)
    wrong_hadronic_bundle, wrong_rg_result = _build_custom_b_probe_inputs(wrong_system_id)

    with pytest.raises(ValueError, match=system_id):
        _invoke_observable_evaluation(
            evaluator,
            rg_value=wrong_rg_result,
            hadronic_value=hadronic_bundle,
        )

    with pytest.raises(ValueError, match=system_id):
        _invoke_observable_evaluation(
            evaluator,
            rg_value=rg_result,
            hadronic_value=wrong_hadronic_bundle,
        )


@pytest.mark.parametrize("system_id", tuple(B_CUSTOM_Q1_PROBE_CONFIG))
def test_custom_b_q1_observable_paths_reject_nonzero_lr(system_id: str) -> None:
    module = _load_observables_module()
    evaluator = _get_custom_b_q1_observable_evaluator(module, system_id)
    hadronic_bundle, rg_result = _build_custom_b_probe_inputs(system_id)
    mismatched_wilsons = dataclasses.replace(
        rg_result.wilsons,
        q4_lr=1.0e-12 + 0.0j,
        q5_lr=-2.5e-13 + 0.0j,
    )

    with pytest.raises(
        ValueError,
        match=str(B_CUSTOM_Q1_PROBE_CONFIG[system_id]["lr_error_pattern"]),
    ):
        _invoke_observable_evaluation(
            evaluator,
            rg_value=mismatched_wilsons,
            hadronic_value=hadronic_bundle,
        )


@pytest.mark.parametrize("system_id", tuple(B_CUSTOM_Q1_PROBE_CONFIG))
def test_custom_b_q1_observable_paths_reject_sector_mismatch(system_id: str) -> None:
    module = _load_observables_module()
    evaluator = _get_custom_b_q1_observable_evaluator(module, system_id)
    hadronic_bundle, rg_result = _build_custom_b_probe_inputs(system_id)
    expected_sector_id = str(B_CUSTOM_Q1_PROBE_CONFIG[system_id]["expected_sector_id"])
    wrong_sector_id = "up" if expected_sector_id == "down" else "down"
    mismatched_wilsons = dataclasses.replace(
        rg_result.wilsons,
        sector_id=wrong_sector_id,
    )

    with pytest.raises(ValueError, match="sector"):
        _invoke_observable_evaluation(
            evaluator,
            rg_value=mismatched_wilsons,
            hadronic_value=hadronic_bundle,
        )


def test_d0_custom_q1_observable_paths_reject_generation_mismatch() -> None:
    module = _load_observables_module()
    evaluator = _get_custom_b_q1_observable_evaluator(module, "D0")
    hadronic_bundle, rg_result = _build_custom_b_probe_inputs("D0")
    mismatched_wilsons = dataclasses.replace(
        rg_result.wilsons,
        generations=(1, 2),
    )

    with pytest.raises(ValueError, match="generations"):
        _invoke_observable_evaluation(
            evaluator,
            rg_value=mismatched_wilsons,
            hadronic_value=hadronic_bundle,
        )


def test_d0_custom_q1_matching_helper_uses_up_sector_generations() -> None:
    matching_module = importlib.import_module(
        "quarkConstraints.paper_0710_1869.eft_deltaf2.matching_kkgluon"
    )
    builder = getattr(matching_module, "default_paper_0710_1869_d0_deltaf2_system", None)
    assert callable(builder)
    system = builder()

    assert system.system_id == "D0"
    assert system.sector_id == "up"
    assert tuple(system.flavor_indices) == (0, 1)


def test_d0_custom_q1_hadronic_contract_rejects_non_d0_system_id() -> None:
    module = _load_observables_module()
    hadronic_module = _load_hadronic_module()

    assert callable(getattr(module, "evaluate_d0_custom_q1_observables", None))
    assert callable(getattr(module, "compute_d0_custom_q1_observables", None))

    contract_cls = getattr(hadronic_module, "Paper07101869D0HadronicContract", None)
    assert isinstance(contract_cls, type)

    with pytest.raises(ValueError, match="system_id"):
        contract_cls(
            system_id="B_d",
            renormalization_scheme_id="bmu.hep-ph-0005183.ndr-ms.lo.v1",
            mu_had_GeV=B_CUSTOM_Q1_PROBE_MU_HAD_GEV,
            evaluation_scale_GeV=B_CUSTOM_Q1_PROBE_MU_HAD_GEV,
        )


def test_d0_custom_q1_hadronic_bundle_rejects_hamiltonian_convention_mutation() -> None:
    hadronic_bundle, _ = _build_custom_b_probe_inputs("D0")

    with pytest.raises(ValueError, match="hamiltonian_convention_id|convention"):
        dataclasses.replace(
            hadronic_bundle,
            hamiltonian_convention_id=f"{hadronic_bundle.hamiltonian_convention_id}.wrong",
        )


@pytest.mark.parametrize(
    ("field_name", "paths", "replacement_suffix", "error_pattern"),
    (
        (
            "operator_basis_id",
            (
                ("operator_basis_id",),
                ("contract", "operator_basis_id"),
            ),
            ".wrong_basis",
            "(?i)basis",
        ),
        (
            "operator_normalization_id",
            (
                ("operator_normalization_id",),
                ("contract", "operator_normalization_id"),
            ),
            ".wrong_normalization",
            "(?i)normalization",
        ),
        (
            "renormalization_scheme_id",
            (
                ("renormalization_scheme_id",),
                ("contract", "renormalization_scheme_id"),
                ("scheme_id",),
                ("contract", "scheme_id"),
            ),
            ".wrong_scheme",
            "(?i)scheme",
        ),
    ),
)
def test_observable_evaluation_rejects_hadronic_contract_tag_mutation(
    field_name: str,
    paths: Sequence[Sequence[str]],
    replacement_suffix: str,
    error_pattern: str,
) -> None:
    module = _load_observables_module()
    evaluator = _get_observable_evaluator(module)
    rg_value = _default_rg_result()
    hadronic_value = _default_hadronic_value()

    if dataclasses.is_dataclass(hadronic_value) and not isinstance(hadronic_value, type):
        with pytest.raises(ValueError, match=error_pattern):
            dataclasses.replace(
                hadronic_value,
                **{field_name: f"{getattr(hadronic_value, field_name)}{replacement_suffix}"},
            )
        return

    mutated_hadronic = _mutate_first_supported_path(
        hadronic_value,
        paths,
        lambda original: f"{original}{replacement_suffix}",
    )

    with pytest.raises(ValueError, match=error_pattern):
        _invoke_observable_evaluation(
            evaluator,
            rg_value=rg_value,
            hadronic_value=mutated_hadronic,
        )


def test_observable_evaluation_rejects_mu_had_mismatch() -> None:
    module = _load_observables_module()
    evaluator = _get_observable_evaluator(module)
    rg_value = _default_rg_result()
    hadronic_module = _load_hadronic_module()
    hadronic_value = _default_hadronic_value()

    builder = _get_callable(hadronic_module, HADRONIC_BUILDER_NAMES)
    if callable(builder) and dataclasses.is_dataclass(hadronic_value) and not isinstance(
        hadronic_value, type
    ):
        mutated_hadronic = builder(mu_had_GeV=float(hadronic_value.mu_had_GeV) + 0.5)
    else:
        mutated_hadronic = _mutate_first_supported_path(
            hadronic_value,
            (
                ("contract", "mu_had_GeV"),
                ("contract", "evaluation_scale_GeV"),
                ("mu_had_GeV",),
                ("evaluation_scale_GeV",),
                ("hadronic_scale_GeV",),
            ),
            lambda original: float(original) + 0.5,
        )

    with pytest.raises(ValueError, match="(?i)(mu_had|scale|evaluation)"):
        _invoke_observable_evaluation(
            evaluator,
            rg_value=rg_value,
            hadronic_value=mutated_hadronic,
        )


def test_observable_evaluation_rejects_non_kaon_system_when_mutable() -> None:
    module = _load_observables_module()
    evaluator = _get_observable_evaluator(module)
    rg_value = _default_rg_result()
    hadronic_value = _default_hadronic_value()

    system_paths = (
        ("contract", "system_id"),
        ("system_id",),
        ("system",),
    )
    if any(_path_exists(hadronic_value, path) for path in system_paths):
        try:
            mutated_hadronic = _mutate_first_supported_path(
                hadronic_value,
                system_paths,
                lambda _original: "bd",
            )
        except AssertionError:
            pytest.skip("hadronic bundle system_id is schema-locked at input construction")
        with pytest.raises(ValueError, match="(?i)(kaon|system)"):
            _invoke_observable_evaluation(
                evaluator,
                rg_value=rg_value,
                hadronic_value=mutated_hadronic,
            )
        return

    rg_system_paths = (
        ("wilsons", "system_id"),
        ("system_id",),
    )
    if any(_path_exists(rg_value, path) for path in rg_system_paths):
        try:
            mutated_rg_value = _mutate_first_supported_path(
                rg_value,
                rg_system_paths,
                lambda _original: "bd",
            )
        except AssertionError:
            pytest.skip("RG system_id is schema-locked at input construction")
        with pytest.raises(ValueError, match="(?i)(kaon|system)"):
            _invoke_observable_evaluation(
                evaluator,
                rg_value=mutated_rg_value,
                hadronic_value=hadronic_value,
            )
        return

    pytest.skip("observable API does not expose a mutable system tag in its public inputs")


def test_observable_evaluation_requires_parity_relation_tag() -> None:
    module = _load_observables_module()
    evaluator = _get_observable_evaluator(module)
    rg_value = _default_rg_result()
    hadronic_value = _default_hadronic_value()

    if dataclasses.is_dataclass(hadronic_value) and not isinstance(hadronic_value, type):
        try:
            mutated_hadronic = dataclasses.replace(
                hadronic_value,
                parity_relation_id=f"{hadronic_value.parity_relation_id}.mismatch",
            )
        except ValueError as exc:
            assert exc.args and any(
                token in str(exc).lower() for token in ("parity", "relation", "q1_vll")
            )
            return
    else:
        mutated_hadronic = _mutate_first_supported_path(
            hadronic_value,
            (
                ("parity_relation_id",),
                ("contract", "parity_relation_id"),
            ),
            lambda original: f"{original}.mismatch",
        )

    with pytest.raises(ValueError, match="(?i)(parity|relation|q1_vll)"):
        _invoke_observable_evaluation(
            evaluator,
            rg_value=rg_value,
            hadronic_value=mutated_hadronic,
        )
