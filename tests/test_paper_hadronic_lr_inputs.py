from __future__ import annotations

import dataclasses
import importlib
import inspect
import math
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
PACKAGE_ROOT = REPO_ROOT / "quarkConstraints" / "paper_0710_1869"
HADRONIC_MODULE_PATH = PACKAGE_ROOT / "eft_deltaf2" / "hadronic.py"

DEFAULT_HADRONIC_EXPORT_NAMES = (
    "default_paper_0710_1869_default_kaon_hadronic",
    "default_paper_0710_1869_kaon_hadronic",
    "default_paper_0710_1869_hadronic_inputs",
    "default_paper_0710_1869_hadronic",
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
DEFAULT_LR_HADRONIC_BUILDER_NAMES = (
    "build_paper_0710_1869_default_kaon_lr_hadronic_inputs",
)
DEFAULT_LR_HADRONIC_EXPORT_NAMES = (
    "default_paper_0710_1869_kaon_lr_hadronic_inputs",
    "default_paper_0710_1869_kaon_lr_hadronic_bundle",
    "default_paper_0710_1869_kaon_lr_hadronic",
    "default_paper_0710_1869_kaon_lr_hadronic_summary",
)
DEFAULT_LR_HADRONIC_VALUE_EXPORT_NAMES = DEFAULT_LR_HADRONIC_EXPORT_NAMES[:-1]
KAON_LR_R_CHI_FREEZE_EXPORT_NAMES = (
    "default_paper_0710_1869_kaon_lr_r_chi_freeze",
    "build_paper_0710_1869_kaon_lr_r_chi_freeze",
)
KAON_LR_R_CHI_SUMMARY_EXPORT_NAMES = (
    "default_paper_0710_1869_kaon_lr_r_chi_summary",
    "build_paper_0710_1869_kaon_lr_r_chi_summary",
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
LR_HADRONIC_PROBE_BUNDLE_ID = "hadronic.kaon.lr_custom_probe.v1"
LR_HADRONIC_PROBE_SOURCE_ID = "hadronic.kaon.lr_custom_probe.sources.v1"
LR_HADRONIC_PROBE_B4_SOURCE_ID = "hadronic.kaon.lr_custom_probe.b4.v1"
LR_HADRONIC_PROBE_B5_SOURCE_ID = "hadronic.kaon.lr_custom_probe.b5.v1"
LR_HADRONIC_PROBE_R_CHI_SOURCE_ID = "hadronic.kaon.lr_custom_probe.rchi.v1"
EXPECTED_KAON_LR_R_CHI_MASS_SCHEME_ID = "pdg.2024.msbar.running_masses.at_2gev.v1"
EXPECTED_KAON_LR_R_CHI_ACTIVE_FLAVOR_POLICY_ID = (
    "pdg.2024.quark_masses.n_l_4.at_2gev.explicit.v1"
)
EXPECTED_KAON_LR_R_CHI_NO_HIDDEN_CONVERSION_POLICY_ID = (
    "none.freeze_pdg2024_msbar_nl4_inputs_at_2gev.v1"
)
EXPECTED_KAON_LR_R_CHI_M_S_2GEV_GEV = 0.09274
EXPECTED_KAON_LR_R_CHI_M_D_2GEV_GEV = 0.00469
EXPECTED_KAON_LR_R_CHI_M_K0_GEV = 0.497611
EXPECTED_KAON_LR_R_CHI_EXACT_VALUE = 26.085222120747908
EXPECTED_KAON_LR_R_CHI_FREEZE_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.eft_deltaf2.kaon_lr_r_chi_freeze.v1"
)
EXPECTED_KAON_LR_R_CHI_SUMMARY_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.eft_deltaf2.kaon_lr_r_chi_summary.v1"
)
EXPECTED_DEFAULT_LR_HADRONIC_BUNDLE_ID = "hadronic.kaon.lr.default.etm2013_ms_2gev.v1"
EXPECTED_DEFAULT_LR_HADRONIC_SOURCE_ID = (
    "hadronic.kaon.lr.default.etm2013_ms_2gev.aggregate.v1"
)
EXPECTED_DEFAULT_LR_HADRONIC_INPUT_POLICY_ID = (
    "default_source.etm2013.table1.ms_2gev.no_hidden_conversion.v1"
)
EXPECTED_CUSTOM_LR_HADRONIC_INPUT_POLICY_ID = (
    "custom_input_only.default_lr_bundle_frozen_separately.no_auto_consumption.v1"
)
EXPECTED_DEFAULT_LR_HADRONIC_B4_SOURCE_ID = (
    "hadronic.kaon.lr.b4.etm2013.table1.ms_2gev.v1"
)
EXPECTED_DEFAULT_LR_HADRONIC_B5_SOURCE_ID = (
    "hadronic.kaon.lr.b5.etm2013.table1.ms_2gev.v1"
)
EXPECTED_DEFAULT_LR_HADRONIC_R_CHI_SOURCE_ID = "hadronic.kaon.lr.r_chi.pdg2024_msbar_nl4.v1"
EXPECTED_DEFAULT_LR_HADRONIC_MASS_SOURCE_ID = "pdg.2024.k0.mass.v1"
EXPECTED_DEFAULT_LR_HADRONIC_DECAY_CONSTANT_SOURCE_ID = "pdg.2024.fkplus.eq72.14.v1"
EXPECTED_DEFAULT_LR_HADRONIC_ETM_CITATION = (
    "ETM Collaboration, JHEP 03 (2013) 089, arXiv:1207.1287"
)
EXPECTED_DEFAULT_LR_HADRONIC_B4_VALUE = 0.78
EXPECTED_DEFAULT_LR_HADRONIC_B5_VALUE = 0.57
EXPECTED_DEFAULT_LR_HADRONIC_OPERATOR_BASIS_ID = "kk_gluon_tree_np_only.v1"
EXPECTED_DEFAULT_LR_HADRONIC_OPERATOR_NORMALIZATION_ID = (
    "paper_0710_1869.deltaf2.kk_gluon_tree_color_normalization.v1"
)
EXPECTED_DEFAULT_LR_HADRONIC_HAMILTONIAN_CONVENTION_ID = (
    "heff.sum_ci_qi.no_hc_factor.v1"
)
EXPECTED_DEFAULT_LR_HADRONIC_Q4_MATRIX_ELEMENT_FORMULA_ID = (
    "kaon.q4_lr.o4_scalar_lr.1over2_rchi_fk2_mk2_b4_mu.plpr_projectors.v2"
)
EXPECTED_DEFAULT_LR_HADRONIC_Q5_MATRIX_ELEMENT_FORMULA_ID = (
    "kaon.q5_lr.o5_scalar_lr.1over6_rchi_fk2_mk2_b5_mu.plpr_projectors.v2"
)
EXPECTED_DEFAULT_LR_HADRONIC_PROVENANCE_IDS = [
    EXPECTED_DEFAULT_LR_HADRONIC_SOURCE_ID,
    EXPECTED_DEFAULT_LR_HADRONIC_B4_SOURCE_ID,
    EXPECTED_DEFAULT_LR_HADRONIC_B5_SOURCE_ID,
    EXPECTED_DEFAULT_LR_HADRONIC_R_CHI_SOURCE_ID,
    EXPECTED_DEFAULT_LR_HADRONIC_MASS_SOURCE_ID,
    EXPECTED_DEFAULT_LR_HADRONIC_DECAY_CONSTANT_SOURCE_ID,
]
EXPECTED_DEFAULT_LR_HADRONIC_SUPPORTED_OPERATORS = ["Q4_LR", "Q5_LR"]
EXPECTED_DEFAULT_LR_HADRONIC_UNSUPPORTED_OPERATORS = ["Q1_VLL", "Q1_VRR"]


def _load_hadronic_module():
    if not HADRONIC_MODULE_PATH.exists():
        pytest.skip("paper_0710_1869 hadronic layer not implemented yet")
    return importlib.import_module("quarkConstraints.paper_0710_1869.eft_deltaf2.hadronic")


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
        raise AssertionError("hadronic payload must canonicalize to a mapping")
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


def _default_kaon_hadronic(module: Any) -> Any:
    callable_obj = _get_callable(module, DEFAULT_HADRONIC_EXPORT_NAMES)
    if not callable(callable_obj):
        raise AssertionError(
            "hadronic layer exists but exposes no default kaon hadronic export; expected one "
            "of " + ", ".join(DEFAULT_HADRONIC_EXPORT_NAMES)
        )
    return callable_obj()


def _default_kaon_hadronic_builder(module: Any) -> Any:
    builder = _get_callable(module, HADRONIC_BUILDER_NAMES)
    if not callable(builder):
        raise AssertionError(
            "hadronic layer exists but exposes no default builder callable; expected one of "
            + ", ".join(HADRONIC_BUILDER_NAMES)
        )
    return builder()


def _default_hadronic_bundle(module: Any) -> Any:
    return _default_kaon_hadronic_builder(module)


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


def _default_kaon_lr_r_chi_freeze(module: Any) -> Any:
    callable_obj = _get_callable(module, KAON_LR_R_CHI_FREEZE_EXPORT_NAMES)
    if not callable(callable_obj):
        raise AssertionError(
            "hadronic layer exposes no frozen LR R_chi helper; expected one of "
            + ", ".join(KAON_LR_R_CHI_FREEZE_EXPORT_NAMES)
        )
    return callable_obj()


def _default_kaon_lr_hadronic(module: Any) -> Any:
    callable_obj = _get_callable(module, DEFAULT_LR_HADRONIC_VALUE_EXPORT_NAMES)
    if not callable(callable_obj):
        raise AssertionError(
            "hadronic layer exposes no default LR bundle helper; expected one of "
            + ", ".join(DEFAULT_LR_HADRONIC_VALUE_EXPORT_NAMES)
        )
    return callable_obj()


def _default_kaon_lr_hadronic_builder(module: Any) -> Any:
    callable_obj = _get_callable(module, DEFAULT_LR_HADRONIC_BUILDER_NAMES)
    if not callable(callable_obj):
        raise AssertionError(
            "hadronic layer exposes no default LR builder helper; expected one of "
            + ", ".join(DEFAULT_LR_HADRONIC_BUILDER_NAMES)
        )
    return callable_obj()


def _default_kaon_lr_r_chi_summary(module: Any) -> dict[str, Any]:
    callable_obj = _get_callable(module, KAON_LR_R_CHI_SUMMARY_EXPORT_NAMES)
    if not callable(callable_obj):
        raise AssertionError(
            "hadronic layer exposes no frozen LR R_chi summary helper; expected one of "
            + ", ".join(KAON_LR_R_CHI_SUMMARY_EXPORT_NAMES)
        )
    return _payload_from_value(callable_obj())


def _build_source_ref(
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
            "hadronic source-ref constructor requires unsupported arguments for LR-HAD-1 "
            "acceptance tests: "
            + ", ".join(missing)
        )
    return source_ref_cls(**kwargs)


def _assign_first(
    parameters: Mapping[str, inspect.Parameter],
    kwargs: dict[str, Any],
    assigned: dict[str, str],
    canonical_name: str,
    candidate_names: Sequence[str],
    value: Any,
) -> None:
    for name in candidate_names:
        if name in parameters:
            kwargs[name] = value
            assigned[canonical_name] = name
            return


def _custom_lr_builder_kwargs(
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
) -> tuple[dict[str, Any], dict[str, str]]:
    parameters = inspect.signature(builder).parameters
    kwargs: dict[str, Any] = {}
    assigned: dict[str, str] = {}

    _assign_first(parameters, kwargs, assigned, "B4_mu_had", LR_HADRONIC_B4_PARAM_NAMES, B4_mu_had)
    _assign_first(parameters, kwargs, assigned, "B5_mu_had", LR_HADRONIC_B5_PARAM_NAMES, B5_mu_had)
    _assign_first(
        parameters,
        kwargs,
        assigned,
        "R_chi_mu_had",
        LR_HADRONIC_R_CHI_PARAM_NAMES,
        R_chi_mu_had,
    )
    _assign_first(
        parameters,
        kwargs,
        assigned,
        "B4_source",
        LR_HADRONIC_B4_SOURCE_PARAM_NAMES,
        B4_source,
    )
    _assign_first(
        parameters,
        kwargs,
        assigned,
        "B5_source",
        LR_HADRONIC_B5_SOURCE_PARAM_NAMES,
        B5_source,
    )
    _assign_first(
        parameters,
        kwargs,
        assigned,
        "R_chi_source",
        LR_HADRONIC_R_CHI_SOURCE_PARAM_NAMES,
        R_chi_source,
    )

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
    if "bundle_id" in parameters:
        kwargs["bundle_id"] = LR_HADRONIC_PROBE_BUNDLE_ID
    if "source_id" in parameters:
        kwargs["source_id"] = LR_HADRONIC_PROBE_SOURCE_ID
    if "provenance_ids" in parameters:
        provenance_ids = [LR_HADRONIC_PROBE_SOURCE_ID]
        for attr_name in ("mass_source", "decay_constant_source", "bag_parameter_source"):
            source = getattr(default_bundle, attr_name, None)
            source_id = getattr(source, "source_id", None)
            if source_id:
                provenance_ids.append(str(source_id))
        provenance_ids.extend(
            (
                LR_HADRONIC_PROBE_B4_SOURCE_ID,
                LR_HADRONIC_PROBE_B5_SOURCE_ID,
                LR_HADRONIC_PROBE_R_CHI_SOURCE_ID,
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
            "custom LR hadronic builder requires unsupported arguments for LR-HAD-1 "
            "acceptance tests: "
            + ", ".join(missing)
        )
    return kwargs, assigned


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


def _extract_lr_matrix_element(payload: Mapping[str, Any], operator_name: str) -> float | None:
    short_name = operator_name.lower().split("_", 1)[0]
    candidate = _nested_value(
        payload,
        (
            (f"{short_name}_matrix_element_GeV4",),
            (f"{short_name}_matrix_element_gev4",),
            (f"{operator_name.lower()}_matrix_element_GeV4",),
            (f"{operator_name.lower()}_matrix_element_gev4",),
            ("matrix_elements_GeV4", operator_name),
            ("matrix_elements", operator_name, "GeV4"),
            ("matrix_elements", operator_name, "value_GeV4"),
            ("matrix_elements", operator_name),
            ("lr_matrix_elements_GeV4", operator_name),
            ("lr_matrix_elements", operator_name, "GeV4"),
            ("lr_matrix_elements", operator_name),
        ),
    )
    if candidate is None:
        return None
    try:
        return float(candidate)
    except (TypeError, ValueError):
        return None


def _source_metadata_text(payload: Mapping[str, Any], source_key: str) -> str:
    source_payload = payload.get(source_key)
    if not isinstance(source_payload, Mapping):
        raise AssertionError(f"{source_key} must be present in the LR hadronic payload")
    return " ".join(
        str(source_payload.get(field_name, ""))
        for field_name in ("citation", "locator_label", "notes")
    )


def _normalized_metadata_text(value: object) -> str:
    lowered = str(value).strip().lower().replace("-", " ").replace("_", " ")
    return " ".join(lowered.split())


def _is_current_custom_lr_input_policy(policy_id: object) -> bool:
    lowered = str(policy_id).strip().lower()
    return bool(lowered) and lowered == EXPECTED_CUSTOM_LR_HADRONIC_INPUT_POLICY_ID.lower()


def _has_current_custom_lr_note_core(notes: object) -> bool:
    lowered = _normalized_metadata_text(notes)
    return bool(lowered) and "defaults not frozen" not in lowered


def _is_current_custom_lr_contract_notes(notes: object) -> bool:
    lowered = _normalized_metadata_text(notes)
    return (
        _has_current_custom_lr_note_core(notes)
        and "frozen default lr hadronic bundle remains separate" in lowered
        and "not auto consumed on this custom surface" in lowered
    )


def _is_current_custom_lr_bundle_notes(notes: object) -> bool:
    lowered = _normalized_metadata_text(notes)
    return (
        _has_current_custom_lr_note_core(notes)
        and "frozen default lr hadronic bundle remains separate" in lowered
        and (
            "lr only and custom combined observable surfaces still "
            "require explicit custom lr inputs"
            in lowered
        )
    )


def test_default_kaon_hadronic_export_matches_builder_default_payload() -> None:
    module = _load_hadronic_module()

    assert _payload_from_value(_default_kaon_hadronic_builder(module)) == _payload_from_value(
        _default_kaon_hadronic(module)
    )


def test_custom_lr_hadronic_builder_is_custom_only_and_matches_bv2004_eq5() -> None:
    module = _load_hadronic_module()
    default_bundle = _default_hadronic_bundle(module)
    builder = _custom_lr_hadronic_builder(module)

    scheme_id = LR_HADRONIC_PROBE_SCHEME_ID
    mu_had_GeV = LR_HADRONIC_PROBE_MU_HAD_GEV
    m_K0_GeV = float(default_bundle.m_K0_GeV)
    f_K_GeV = float(default_bundle.f_K_GeV)

    b4_source = _build_source_ref(
        module,
        source_id=LR_HADRONIC_PROBE_B4_SOURCE_ID,
        citation="Becirevic and Villadoro, hep-lat/0408029",
        locator_label="eq. (5), O4 scalar LR matrix element",
        year=2004,
        scheme_id=scheme_id,
        scale_GeV=mu_had_GeV,
        notes="Custom LR hadronic acceptance probe source for B4(mu_had).",
    )
    b5_source = _build_source_ref(
        module,
        source_id=LR_HADRONIC_PROBE_B5_SOURCE_ID,
        citation="Becirevic and Villadoro, hep-lat/0408029",
        locator_label="eq. (5), O5 scalar LR matrix element",
        year=2004,
        scheme_id=scheme_id,
        scale_GeV=mu_had_GeV,
        notes="Custom LR hadronic acceptance probe source for B5(mu_had).",
    )
    r_chi_source = _build_source_ref(
        module,
        source_id=LR_HADRONIC_PROBE_R_CHI_SOURCE_ID,
        citation="Custom caller-supplied chiral ratio input",
        locator_label="R_chi(mu_had) external input",
        year=2026,
        scheme_id=scheme_id,
        scale_GeV=mu_had_GeV,
        notes="Custom LR hadronic acceptance probe source for R_chi(mu_had).",
    )
    kwargs, _ = _custom_lr_builder_kwargs(
        builder,
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

    built = builder(**kwargs)
    second = builder(**kwargs)
    payload = _payload_from_value(built)
    second_payload = _payload_from_value(second)
    assert payload == second_payload
    assert str(payload["renormalization_scheme_id"]) == scheme_id
    assert float(payload["mu_had_GeV"]) == pytest.approx(
        mu_had_GeV,
        rel=0.0,
        abs=1.0e-12,
    )
    assert str(payload["contract"]["renormalization_scheme_id"]) == scheme_id
    assert float(payload["contract"]["mu_had_GeV"]) == pytest.approx(
        mu_had_GeV,
        rel=0.0,
        abs=1.0e-12,
    )

    input_provenance_mode_id = _nested_value(
        payload,
        (
            ("input_provenance_mode_id",),
            ("tags", "input_provenance_mode_id"),
        ),
    )
    assert input_provenance_mode_id is not None
    assert "custom" in str(input_provenance_mode_id).lower()
    input_policy_id = _nested_value(payload, (("input_policy_id",),))
    contract_input_policy_id = _nested_value(payload, (("contract", "input_policy_id"),))
    bundle_notes = _nested_value(payload, (("notes",),))
    contract_notes = _nested_value(payload, (("contract", "notes"),))
    assert _is_current_custom_lr_input_policy(input_policy_id)
    assert _is_current_custom_lr_input_policy(contract_input_policy_id)
    assert input_policy_id == contract_input_policy_id
    assert _is_current_custom_lr_bundle_notes(bundle_notes)
    assert _is_current_custom_lr_contract_notes(contract_notes)

    contract_id = _nested_value(
        payload,
        (
            ("lr_hadronic_contract_id",),
            ("contract", "lr_hadronic_contract_id"),
            ("contract", "contract_id"),
            ("contract", "schema_id"),
            ("contract", "input_policy_id"),
            ("contract_id",),
        ),
    )
    assert contract_id
    q4_formula_id = _nested_value(
        payload,
        (
            ("q4_matrix_element_formula_id",),
            ("formula_ids", "Q4_LR"),
            ("contract", "q4_matrix_element_formula_id"),
            ("matrix_element_formula_id",),
        ),
    )
    q5_formula_id = _nested_value(
        payload,
        (
            ("q5_matrix_element_formula_id",),
            ("formula_ids", "Q5_LR"),
            ("contract", "q5_matrix_element_formula_id"),
            ("matrix_element_formula_id",),
        ),
    )
    assert q4_formula_id
    assert q5_formula_id

    q4_matrix_element_GeV4 = _extract_lr_matrix_element(payload, "Q4_LR")
    q5_matrix_element_GeV4 = _extract_lr_matrix_element(payload, "Q5_LR")
    assert q4_matrix_element_GeV4 is not None
    assert q5_matrix_element_GeV4 is not None

    # M-27: PL/PR projectors add the audited /4 relative to the old Eq. (5) pin.
    expected_q4 = (
        0.5
        * LR_HADRONIC_PROBE_R_CHI
        * (m_K0_GeV**2)
        * (f_K_GeV**2)
        * LR_HADRONIC_PROBE_B4
    )
    expected_q5 = (
        (1.0 / 6.0)
        * LR_HADRONIC_PROBE_R_CHI
        * (m_K0_GeV**2)
        * (f_K_GeV**2)
        * LR_HADRONIC_PROBE_B5
    )
    assert q4_matrix_element_GeV4 == pytest.approx(expected_q4, rel=0.0, abs=1.0e-15)
    assert q5_matrix_element_GeV4 == pytest.approx(expected_q5, rel=0.0, abs=1.0e-15)

    assert str(b4_source.renormalization_scheme_id) == scheme_id
    assert str(b5_source.renormalization_scheme_id) == scheme_id
    assert str(r_chi_source.renormalization_scheme_id) == scheme_id
    assert math.isclose(
        float(b4_source.scale_GeV),
        mu_had_GeV,
        rel_tol=0.0,
        abs_tol=1.0e-12,
    )
    assert math.isclose(
        float(b5_source.scale_GeV),
        mu_had_GeV,
        rel_tol=0.0,
        abs_tol=1.0e-12,
    )
    assert math.isclose(
        float(r_chi_source.scale_GeV),
        mu_had_GeV,
        rel_tol=0.0,
        abs_tol=1.0e-12,
    )


@pytest.mark.parametrize(
    ("missing_key", "error_pattern"),
    (
        ("B4_source", "(?i)(B4|source)"),
        ("B5_source", "(?i)(B5|source)"),
        ("R_chi_source", "(?i)(R_chi|source)"),
        ("R_chi_mu_had", "(?i)(R_chi)"),
    ),
)
def test_custom_lr_hadronic_builder_rejects_missing_explicit_inputs(
    missing_key: str,
    error_pattern: str,
) -> None:
    module = _load_hadronic_module()
    default_bundle = _default_hadronic_bundle(module)
    builder = _custom_lr_hadronic_builder(module)

    scheme_id = LR_HADRONIC_PROBE_SCHEME_ID
    mu_had_GeV = LR_HADRONIC_PROBE_MU_HAD_GEV
    b4_source = _build_source_ref(
        module,
        source_id=LR_HADRONIC_PROBE_B4_SOURCE_ID,
        citation="Becirevic and Villadoro, hep-lat/0408029",
        locator_label="eq. (5), O4 scalar LR matrix element",
        year=2004,
        scheme_id=scheme_id,
        scale_GeV=mu_had_GeV,
        notes="Custom LR hadronic acceptance probe source for B4(mu_had).",
    )
    b5_source = _build_source_ref(
        module,
        source_id=LR_HADRONIC_PROBE_B5_SOURCE_ID,
        citation="Becirevic and Villadoro, hep-lat/0408029",
        locator_label="eq. (5), O5 scalar LR matrix element",
        year=2004,
        scheme_id=scheme_id,
        scale_GeV=mu_had_GeV,
        notes="Custom LR hadronic acceptance probe source for B5(mu_had).",
    )
    r_chi_source = _build_source_ref(
        module,
        source_id=LR_HADRONIC_PROBE_R_CHI_SOURCE_ID,
        citation="Custom caller-supplied chiral ratio input",
        locator_label="R_chi(mu_had) external input",
        year=2026,
        scheme_id=scheme_id,
        scale_GeV=mu_had_GeV,
        notes="Custom LR hadronic acceptance probe source for R_chi(mu_had).",
    )
    kwargs, assigned = _custom_lr_builder_kwargs(
        builder,
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
    assert missing_key in assigned
    kwargs.pop(assigned[missing_key])

    with pytest.raises((TypeError, ValueError, AssertionError), match=error_pattern):
        builder(**kwargs)


@pytest.mark.parametrize(
    ("source_key", "bad_scheme", "bad_scale_GeV", "error_pattern"),
    (
        (
            "B4_source",
            "wrong.scheme.v1",
            LR_HADRONIC_PROBE_MU_HAD_GEV,
            "(?i)(scheme|B4|source)",
        ),
        (
            "B5_source",
            LR_HADRONIC_PROBE_SCHEME_ID,
            LR_HADRONIC_PROBE_MU_HAD_GEV + 0.25,
            "(?i)(scale|mu_had|B5|source)",
        ),
        (
            "R_chi_source",
            "wrong.scheme.v1",
            LR_HADRONIC_PROBE_MU_HAD_GEV,
            "(?i)(scheme|R_chi|source)",
        ),
        (
            "R_chi_source",
            LR_HADRONIC_PROBE_SCHEME_ID,
            LR_HADRONIC_PROBE_MU_HAD_GEV + 0.25,
            "(?i)(scale|mu_had|R_chi|source)",
        ),
    ),
)
def test_custom_lr_hadronic_builder_rejects_source_scheme_or_scale_mismatch(
    source_key: str,
    bad_scheme: str,
    bad_scale_GeV: float,
    error_pattern: str,
) -> None:
    module = _load_hadronic_module()
    default_bundle = _default_hadronic_bundle(module)
    builder = _custom_lr_hadronic_builder(module)

    scheme_id = LR_HADRONIC_PROBE_SCHEME_ID
    mu_had_GeV = LR_HADRONIC_PROBE_MU_HAD_GEV
    b4_source = _build_source_ref(
        module,
        source_id=LR_HADRONIC_PROBE_B4_SOURCE_ID,
        citation="Becirevic and Villadoro, hep-lat/0408029",
        locator_label="eq. (5), O4 scalar LR matrix element",
        year=2004,
        scheme_id=scheme_id,
        scale_GeV=mu_had_GeV,
        notes="Custom LR hadronic acceptance probe source for B4(mu_had).",
    )
    b5_source = _build_source_ref(
        module,
        source_id=LR_HADRONIC_PROBE_B5_SOURCE_ID,
        citation="Becirevic and Villadoro, hep-lat/0408029",
        locator_label="eq. (5), O5 scalar LR matrix element",
        year=2004,
        scheme_id=scheme_id,
        scale_GeV=mu_had_GeV,
        notes="Custom LR hadronic acceptance probe source for B5(mu_had).",
    )
    r_chi_source = _build_source_ref(
        module,
        source_id=LR_HADRONIC_PROBE_R_CHI_SOURCE_ID,
        citation="Custom caller-supplied chiral ratio input",
        locator_label="R_chi(mu_had) external input",
        year=2026,
        scheme_id=scheme_id,
        scale_GeV=mu_had_GeV,
        notes="Custom LR hadronic acceptance probe source for R_chi(mu_had).",
    )
    if source_key == "B4_source":
        b4_source = _build_source_ref(
            module,
            source_id=LR_HADRONIC_PROBE_B4_SOURCE_ID,
            citation="Becirevic and Villadoro, hep-lat/0408029",
            locator_label="eq. (5), O4 scalar LR matrix element",
            year=2004,
            scheme_id=bad_scheme,
            scale_GeV=bad_scale_GeV,
            notes="Mismatched custom LR hadronic acceptance probe source for B4(mu_had).",
        )
    elif source_key == "B5_source":
        b5_source = _build_source_ref(
            module,
            source_id=LR_HADRONIC_PROBE_B5_SOURCE_ID,
            citation="Becirevic and Villadoro, hep-lat/0408029",
            locator_label="eq. (5), O5 scalar LR matrix element",
            year=2004,
            scheme_id=bad_scheme,
            scale_GeV=bad_scale_GeV,
            notes="Mismatched custom LR hadronic acceptance probe source for B5(mu_had).",
        )
    elif source_key == "R_chi_source":
        r_chi_source = _build_source_ref(
            module,
            source_id=LR_HADRONIC_PROBE_R_CHI_SOURCE_ID,
            citation="Custom caller-supplied chiral ratio input",
            locator_label="R_chi(mu_had) external input",
            year=2026,
            scheme_id=bad_scheme,
            scale_GeV=bad_scale_GeV,
            notes="Mismatched custom LR hadronic acceptance probe source for R_chi(mu_had).",
        )
    else:  # pragma: no cover - parameterized values are fixed above
        raise AssertionError(f"unexpected source_key: {source_key}")
    kwargs, _ = _custom_lr_builder_kwargs(
        builder,
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

    with pytest.raises((TypeError, ValueError, AssertionError), match=error_pattern):
        builder(**kwargs)


def test_default_kaon_lr_r_chi_freeze_is_deterministic_and_matches_bv2004_definition() -> None:
    module = _load_hadronic_module()

    freeze_value = _default_kaon_lr_r_chi_freeze(module)
    second_value = _default_kaon_lr_r_chi_freeze(module)
    payload = _payload_from_value(freeze_value.as_dict())
    second_payload = _payload_from_value(second_value.as_dict())
    summary_payload = _default_kaon_lr_r_chi_summary(module)

    assert payload == second_payload
    assert payload["schema_id"] == EXPECTED_KAON_LR_R_CHI_FREEZE_SCHEMA_ID
    assert payload["system_id"] == "kaon"
    assert float(payload["mu_had_GeV"]) == pytest.approx(2.0, rel=0.0, abs=1.0e-12)
    assert payload["freeze_id"]
    assert payload["source_id"]
    assert "derived" in str(payload["input_provenance_mode_id"]).lower()
    assert payload["operator_renormalization_scheme_id"] == payload["operator_scheme_id"]
    assert payload["mass_renormalization_scheme_id"] == payload["mass_scheme_id"]
    assert payload["mass_renormalization_scheme_id"] == EXPECTED_KAON_LR_R_CHI_MASS_SCHEME_ID
    assert payload["operator_scheme_id"] != payload["mass_scheme_id"]
    assert (
        payload["mass_active_flavor_policy_id"]
        == EXPECTED_KAON_LR_R_CHI_ACTIVE_FLAVOR_POLICY_ID
    )
    assert payload["derivation_formula_id"]
    assert payload["derivation_formula_source_id"]
    assert (
        payload["no_hidden_conversion_policy_id"]
        == EXPECTED_KAON_LR_R_CHI_NO_HIDDEN_CONVERSION_POLICY_ID
    )
    assert float(payload["m_K0_GeV"]) == pytest.approx(
        EXPECTED_KAON_LR_R_CHI_M_K0_GEV,
        rel=0.0,
        abs=1.0e-15,
    )
    assert float(payload["m_s_mu_had_GeV"]) == pytest.approx(
        EXPECTED_KAON_LR_R_CHI_M_S_2GEV_GEV,
        rel=0.0,
        abs=1.0e-15,
    )
    assert float(payload["m_d_mu_had_GeV"]) == pytest.approx(
        EXPECTED_KAON_LR_R_CHI_M_D_2GEV_GEV,
        rel=0.0,
        abs=1.0e-15,
    )
    assert payload["provenance_ids"] == [
        payload["source_id"],
        payload["kaon_mass_source"]["source_id"],
        payload["strange_mass_source"]["source_id"],
        payload["down_mass_source"]["source_id"],
    ]
    assert (
        payload["strange_mass_source"]["renormalization_scheme_id"]
        == payload["mass_renormalization_scheme_id"]
    )
    assert (
        payload["down_mass_source"]["renormalization_scheme_id"]
        == payload["mass_renormalization_scheme_id"]
    )
    assert float(payload["strange_mass_source"]["scale_GeV"]) == pytest.approx(
        float(payload["mu_had_GeV"]),
        rel=0.0,
        abs=1.0e-12,
    )
    assert float(payload["down_mass_source"]["scale_GeV"]) == pytest.approx(
        float(payload["mu_had_GeV"]),
        rel=0.0,
        abs=1.0e-12,
    )

    expected_r_chi = (float(payload["m_K0_GeV"]) / (
        float(payload["m_s_mu_had_GeV"]) + float(payload["m_d_mu_had_GeV"])
    )) ** 2
    assert expected_r_chi == pytest.approx(
        EXPECTED_KAON_LR_R_CHI_EXACT_VALUE,
        rel=0.0,
        abs=1.0e-15,
    )
    assert float(payload["R_chi_mu_had"]) == pytest.approx(
        EXPECTED_KAON_LR_R_CHI_EXACT_VALUE,
        rel=0.0,
        abs=1.0e-15,
    )
    assert float(summary_payload["R_chi_mu_had"]) == pytest.approx(
        EXPECTED_KAON_LR_R_CHI_EXACT_VALUE,
        rel=0.0,
        abs=1.0e-15,
    )
    assert summary_payload["schema_id"] == EXPECTED_KAON_LR_R_CHI_SUMMARY_SCHEMA_ID
    assert summary_payload["schema_id"] != payload["schema_id"]


@pytest.mark.parametrize(
    ("mutator", "error_pattern"),
    (
        (
            lambda freeze: {"provenance_ids": tuple(freeze.provenance_ids[:-1])},
            "(?i)provenance_ids",
        ),
        (
            lambda _freeze: {"derivation_formula_id": ""},
            "(?i)derivation_formula_id",
        ),
        (
            lambda freeze: {"mu_had_GeV": float(freeze.mu_had_GeV) + 0.25},
            "(?i)mu_had_GeV",
        ),
        (
            lambda freeze: {"source_id": f"{freeze.source_id}.drift"},
            "(?i)source_id",
        ),
        (
            lambda _freeze: {
                "no_hidden_conversion_policy_id": "freeze.hidden_conversion.policy_drift.v1"
            },
            "(?i)no_hidden_conversion_policy_id",
        ),
        (
            lambda _freeze: {"input_policy_id": "freeze.hidden_conversion.input_policy_drift.v1"},
            "(?i)input_policy_id",
        ),
        (
            lambda _freeze: {"mass_active_flavor_policy_id": "wrong.active_flavor.v1"},
            "(?i)mass_active_flavor_policy_id",
        ),
        (
            lambda freeze: {
                "mass_renormalization_scheme_id": str(
                    freeze.operator_renormalization_scheme_id
                )
            },
            "(?i)(mass_renormalization_scheme_id|distinct)",
        ),
    ),
)
def test_default_kaon_lr_r_chi_freeze_rejects_metadata_drift(
    mutator: Any,
    error_pattern: str,
) -> None:
    module = _load_hadronic_module()
    freeze_value = _default_kaon_lr_r_chi_freeze(module)

    with pytest.raises(ValueError, match=error_pattern):
        dataclasses.replace(freeze_value, **mutator(freeze_value))


def test_default_kaon_lr_hadronic_bundle_is_deterministic_and_carries_frozen_etm_metadata() -> None:
    module = _load_hadronic_module()
    payload = _payload_from_value(_default_kaon_lr_hadronic(module))
    second_payload = _payload_from_value(_default_kaon_lr_hadronic(module))
    contract = payload["contract"]

    assert payload == second_payload
    assert payload["system_id"] == "kaon"
    assert float(payload["mu_had_GeV"]) == pytest.approx(2.0, rel=0.0, abs=1.0e-12)
    assert payload["bundle_id"] == EXPECTED_DEFAULT_LR_HADRONIC_BUNDLE_ID
    assert payload["source_id"] == EXPECTED_DEFAULT_LR_HADRONIC_SOURCE_ID
    assert payload["input_policy_id"] == EXPECTED_DEFAULT_LR_HADRONIC_INPUT_POLICY_ID
    assert (
        payload["operator_basis_id"]
        == EXPECTED_DEFAULT_LR_HADRONIC_OPERATOR_BASIS_ID
    )
    assert (
        payload["operator_normalization_id"]
        == EXPECTED_DEFAULT_LR_HADRONIC_OPERATOR_NORMALIZATION_ID
    )
    assert (
        payload["hamiltonian_convention_id"]
        == EXPECTED_DEFAULT_LR_HADRONIC_HAMILTONIAN_CONVENTION_ID
    )
    assert (
        payload["q4_matrix_element_formula_id"]
        == EXPECTED_DEFAULT_LR_HADRONIC_Q4_MATRIX_ELEMENT_FORMULA_ID
    )
    assert (
        payload["q5_matrix_element_formula_id"]
        == EXPECTED_DEFAULT_LR_HADRONIC_Q5_MATRIX_ELEMENT_FORMULA_ID
    )
    assert contract["input_policy_id"] == EXPECTED_DEFAULT_LR_HADRONIC_INPUT_POLICY_ID
    assert (
        contract["operator_basis_id"]
        == EXPECTED_DEFAULT_LR_HADRONIC_OPERATOR_BASIS_ID
    )
    assert (
        contract["operator_normalization_id"]
        == EXPECTED_DEFAULT_LR_HADRONIC_OPERATOR_NORMALIZATION_ID
    )
    assert (
        contract["hamiltonian_convention_id"]
        == EXPECTED_DEFAULT_LR_HADRONIC_HAMILTONIAN_CONVENTION_ID
    )
    assert payload["provenance_ids"] == EXPECTED_DEFAULT_LR_HADRONIC_PROVENANCE_IDS
    assert payload["supported_operator_names"] == EXPECTED_DEFAULT_LR_HADRONIC_SUPPORTED_OPERATORS
    assert (
        payload["unsupported_operator_names"]
        == EXPECTED_DEFAULT_LR_HADRONIC_UNSUPPORTED_OPERATORS
    )
    assert float(payload["B4_mu_had"]) == pytest.approx(
        EXPECTED_DEFAULT_LR_HADRONIC_B4_VALUE,
        rel=0.0,
        abs=1.0e-15,
    )
    assert float(payload["B5_mu_had"]) == pytest.approx(
        EXPECTED_DEFAULT_LR_HADRONIC_B5_VALUE,
        rel=0.0,
        abs=1.0e-15,
    )
    assert float(payload["R_chi_mu_had"]) == pytest.approx(
        EXPECTED_KAON_LR_R_CHI_EXACT_VALUE,
        rel=0.0,
        abs=1.0e-15,
    )
    assert payload["b4_source"]["source_id"] == EXPECTED_DEFAULT_LR_HADRONIC_B4_SOURCE_ID
    assert payload["b5_source"]["source_id"] == EXPECTED_DEFAULT_LR_HADRONIC_B5_SOURCE_ID
    assert payload["b4_source"]["transformation_id"] == "none"
    assert payload["b5_source"]["transformation_id"] == "none"
    assert payload["r_chi_source"]["source_id"] == EXPECTED_DEFAULT_LR_HADRONIC_R_CHI_SOURCE_ID
    assert payload["mass_source"]["source_id"] == EXPECTED_DEFAULT_LR_HADRONIC_MASS_SOURCE_ID
    assert (
        payload["decay_constant_source"]["source_id"]
        == EXPECTED_DEFAULT_LR_HADRONIC_DECAY_CONSTANT_SOURCE_ID
    )
    for source_key in ("b4_source", "b5_source"):
        assert (
            payload[source_key]["renormalization_scheme_id"]
            == payload["renormalization_scheme_id"]
        )
        assert float(payload[source_key]["scale_GeV"]) == pytest.approx(
            float(payload["mu_had_GeV"]),
            rel=0.0,
            abs=1.0e-12,
        )
    assert (
        payload["r_chi_source"]["renormalization_scheme_id"]
        == EXPECTED_KAON_LR_R_CHI_MASS_SCHEME_ID
    )
    assert (
        payload["r_chi_source"]["renormalization_scheme_id"]
        != payload["renormalization_scheme_id"]
    )
    assert float(payload["r_chi_source"]["scale_GeV"]) == pytest.approx(
        float(payload["mu_had_GeV"]),
        rel=0.0,
        abs=1.0e-12,
    )

    b4_source_text = _source_metadata_text(payload, "b4_source")
    b5_source_text = _source_metadata_text(payload, "b5_source")
    for source_text, operator_label, bag_value in (
        (b4_source_text, "B4", EXPECTED_DEFAULT_LR_HADRONIC_B4_VALUE),
        (b5_source_text, "B5", EXPECTED_DEFAULT_LR_HADRONIC_B5_VALUE),
    ):
        assert EXPECTED_DEFAULT_LR_HADRONIC_ETM_CITATION in source_text
        assert "Table 1" in source_text
        assert "Buras" in source_text
        assert "2 GeV" in source_text
        assert operator_label in source_text
        assert f"{bag_value:.2f}" in source_text

    # M-27: PL/PR projectors add the audited /4 relative to the old Eq. (5) pin.
    expected_q4 = (
        0.5
        * float(payload["R_chi_mu_had"])
        * (float(payload["m_K0_GeV"]) ** 2)
        * (float(payload["f_K_GeV"]) ** 2)
        * float(payload["B4_mu_had"])
    )
    expected_q5 = (
        (1.0 / 6.0)
        * float(payload["R_chi_mu_had"])
        * (float(payload["m_K0_GeV"]) ** 2)
        * (float(payload["f_K_GeV"]) ** 2)
        * float(payload["B5_mu_had"])
    )
    assert float(payload["q4_matrix_element_GeV4"]) == pytest.approx(
        expected_q4,
        rel=0.0,
        abs=1.0e-15,
    )
    assert float(payload["q5_matrix_element_GeV4"]) == pytest.approx(
        expected_q5,
        rel=0.0,
        abs=1.0e-15,
    )


def test_default_kaon_lr_hadronic_export_matches_builder_default_payload() -> None:
    module = _load_hadronic_module()
    export_payload = _payload_from_value(_default_kaon_lr_hadronic(module))
    builder_payload = _payload_from_value(_default_kaon_lr_hadronic_builder(module))

    assert export_payload == builder_payload


@pytest.mark.parametrize(
    ("mutator", "error_pattern"),
    (
        (
            lambda bundle: {"B4_mu_had": float(bundle.B4_mu_had) + 0.01},
            "(?i)(B4\\(2 GeV\\)|B4_mu_had)",
        ),
        (
            lambda bundle: {"B5_mu_had": float(bundle.B5_mu_had) + 0.01},
            "(?i)(B5\\(2 GeV\\)|B5_mu_had)",
        ),
        (
            lambda bundle: {"source_id": f"{bundle.source_id}.drift"},
            "(?i)source_id",
        ),
        (
            lambda bundle: {
                "b4_source": dataclasses.replace(
                    bundle.b4_source,
                    source_id=f"{bundle.b4_source.source_id}.drift",
                )
            },
            "(?i)(B4 source|b4_source|source_id)",
        ),
        (
            lambda bundle: {
                "b5_source": dataclasses.replace(
                    bundle.b5_source,
                    source_id=f"{bundle.b5_source.source_id}.drift",
                )
            },
            "(?i)(B5 source|b5_source|source_id)",
        ),
        (
            lambda bundle: {
                "r_chi_source": dataclasses.replace(
                    bundle.r_chi_source,
                    source_id=f"{bundle.r_chi_source.source_id}.drift",
                )
            },
            "(?i)(R_chi source|r_chi_source|source_id)",
        ),
        (
            lambda bundle: {
                "renormalization_scheme_id": f"{bundle.renormalization_scheme_id}.drift"
            },
            "(?i)(renormalization_scheme_id|scheme)",
        ),
        (
            lambda bundle: {"operator_basis_id": f"{bundle.operator_basis_id}.drift"},
            "(?i)operator_basis_id",
        ),
        (
            lambda bundle: {
                "operator_normalization_id": f"{bundle.operator_normalization_id}.drift"
            },
            "(?i)operator_normalization_id",
        ),
        (
            lambda bundle: {
                "hamiltonian_convention_id": (
                    f"{bundle.hamiltonian_convention_id}.drift"
                )
            },
            "(?i)hamiltonian_convention_id",
        ),
        (
            lambda bundle: {
                "q4_matrix_element_formula_id": (
                    f"{bundle.q4_matrix_element_formula_id}.drift"
                )
            },
            "(?i)q4_matrix_element_formula_id",
        ),
        (
            lambda bundle: {
                "q5_matrix_element_formula_id": (
                    f"{bundle.q5_matrix_element_formula_id}.drift"
                )
            },
            "(?i)q5_matrix_element_formula_id",
        ),
        (
            lambda bundle: {
                "r_chi_source": dataclasses.replace(
                    bundle.r_chi_source,
                    renormalization_scheme_id=str(bundle.renormalization_scheme_id),
                )
            },
            "(?i)(r_chi_source\\.renormalization_scheme_id|frozen default mass scheme)",
        ),
        (
            lambda bundle: {"mu_had_GeV": float(bundle.mu_had_GeV) + 0.25},
            "(?i)(mu_had|scale)",
        ),
        (
            lambda _bundle: {
                "input_policy_id": "default_source.etm2013.table1.ms_2gev.hidden_conversion.v1"
            },
            "(?i)input_policy_id",
        ),
        (
            lambda bundle: {
                "contract": dataclasses.replace(
                    bundle.contract,
                    input_policy_id=f"{bundle.contract.input_policy_id}.drift",
                )
            },
            "(?i)contract\\.input_policy_id",
        ),
        (
            lambda bundle: {
                "contract": dataclasses.replace(
                    bundle.contract,
                    operator_basis_id=f"{bundle.contract.operator_basis_id}.drift",
                )
            },
            "(?i)contract\\.operator_basis_id",
        ),
        (
            lambda bundle: {
                "contract": dataclasses.replace(
                    bundle.contract,
                    operator_normalization_id=(
                        f"{bundle.contract.operator_normalization_id}.drift"
                    ),
                )
            },
            "(?i)contract\\.operator_normalization_id",
        ),
        (
            lambda bundle: {
                "contract": dataclasses.replace(
                    bundle.contract,
                    hamiltonian_convention_id=(
                        f"{bundle.contract.hamiltonian_convention_id}.drift"
                    ),
                )
            },
            "(?i)contract\\.hamiltonian_convention_id",
        ),
        (
            lambda bundle: {
                "contract": dataclasses.replace(
                    bundle.contract,
                    q4_matrix_element_formula_id=(
                        f"{bundle.contract.q4_matrix_element_formula_id}.drift"
                    ),
                )
            },
            "(?i)contract\\.q4_matrix_element_formula_id",
        ),
        (
            lambda bundle: {
                "contract": dataclasses.replace(
                    bundle.contract,
                    q5_matrix_element_formula_id=(
                        f"{bundle.contract.q5_matrix_element_formula_id}.drift"
                    ),
                )
            },
            "(?i)contract\\.q5_matrix_element_formula_id",
        ),
        (
            lambda bundle: {
                "contract": dataclasses.replace(
                    bundle.contract,
                    formula_source_id=f"{bundle.contract.formula_source_id}.drift",
                )
            },
            "(?i)contract\\.formula_source_id",
        ),
        (
            lambda bundle: {
                "r_chi_source": dataclasses.replace(
                    bundle.r_chi_source,
                    transformation_id="ri_mom_to_ms.hidden_conversion_probe.v1",
                )
            },
            "(?i)(R_chi source|r_chi_source|frozen default)",
        ),
    ),
)
def test_default_kaon_lr_hadronic_bundle_rejects_frozen_default_drift(
    mutator: Any,
    error_pattern: str,
) -> None:
    module = _load_hadronic_module()
    default_bundle = _default_kaon_lr_hadronic(module)

    with pytest.raises(ValueError, match=error_pattern):
        dataclasses.replace(default_bundle, **mutator(default_bundle))
