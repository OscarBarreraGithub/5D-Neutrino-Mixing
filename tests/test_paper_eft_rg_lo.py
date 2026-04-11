"""Forward-compatible tests for the paper-owned LO Delta F=2 RG layer."""

from __future__ import annotations

import dataclasses
import importlib
import inspect
import json
import subprocess
import sys
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import pytest

from qcd.constants import M_BOTTOM, M_CHARM, M_TOP_MS
from qcd.running import alpha_s

REPO_ROOT = Path(__file__).resolve().parents[1]
PACKAGE_ROOT = REPO_ROOT / "quarkConstraints" / "paper_0710_1869"
RG_MODULE_PATH = PACKAGE_ROOT / "eft_deltaf2" / "rg.py"

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
IMPLEMENTATION_DEFAULT_RG_SUMMARY_NAMES = (
    "default_paper_0710_1869_default_kaon_rg_summary",
    "build_paper_0710_1869_default_kaon_rg_summary",
)
RG_ALIAS_STEMS = ("kaon_rg", "deltaf2_rg", "lo_rg", "rg_lo")
DEFAULT_RG_SUMMARY_NAMES = tuple(
    name
    for stem in RG_ALIAS_STEMS
    for name in (
        f"default_paper_0710_1869_{stem}_summary",
        f"build_paper_0710_1869_{stem}_summary",
    )
) + IMPLEMENTATION_DEFAULT_RG_SUMMARY_NAMES
DEFAULT_RG_EXPORT_NAMES = tuple(
    name
    for stem in RG_ALIAS_STEMS
    for name in (
        f"default_paper_0710_1869_{stem}",
        f"build_paper_0710_1869_{stem}",
    )
)
PARAMETERIZED_RG_NAMES = (
    "evolve_deltaf2_wilsons_lo",
    "run_deltaf2_wilsons_lo",
    "run_paper_0710_1869_deltaf2_wilsons_lo",
    "evolve_paper_0710_1869_deltaf2_wilsons_lo",
)
RG_DEFAULT_LOW_SCALE_GEV = 2.0
EXPECTED_RG_LR_BASIS_CONTRACT_ID = (
    "paper_q4q5_to_bmu_lr_basis.map_frozen.audit_ready.v2"
)
EXPECTED_RG_LR_BASIS_DIRECTION_ID = (
    "paper_q4q5_operator_order.to_bmu_q1lr_q2lr_operator_order.map_frozen.v2"
)
EXPECTED_RG_LR_BMU_BASIS_ID = "bmu.lr_basis.q1lr_q2lr.ndr_ms.lo.v1"
EXPECTED_PAPER_LR_OPERATOR_DEFINITION_IDS = (
    "paper_0710_1869.deltaf2.q4_lr.susy_o4.scalar_lr.color_singlet.v1",
    "paper_0710_1869.deltaf2.q5_lr.susy_o5.scalar_lr.color_mixed.v1",
)
EXPECTED_FULL_RG_SUPPORTED_OPERATORS = (
    "Q1_VLL",
    "Q1_VRR",
    "Q4_LR",
    "Q5_LR",
)


def _has_rg_module() -> bool:
    return RG_MODULE_PATH.exists()


def _load_rg_module():
    if not _has_rg_module():
        pytest.skip("paper_0710_1869 LO RG layer not implemented yet")
    return importlib.import_module("quarkConstraints.paper_0710_1869.eft_deltaf2.rg")


def _rg_inputs_ids() -> tuple[str, str, str]:
    module = importlib.import_module("quarkConstraints.paper_0710_1869.eft_deltaf2.rg_inputs")
    return (
        str(module.PAPER_0710_1869_DELTAF2_RG_LR_BASIS_CONTRACT_ID),
        str(module.PAPER_0710_1869_DELTAF2_RG_LR_BASIS_STATUS_ID),
        str(module.PAPER_0710_1869_DELTAF2_RG_SUPPORTED_OPERATOR_SUBSET_ID),
    )


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
        raise AssertionError("RG payload must canonicalize to a mapping")
    return dict(payload)


def _get_callable(module: Any, names: Sequence[str]) -> Any:
    for name in names:
        candidate = getattr(module, name, None)
        if callable(candidate):
            return candidate
    return None


def _first_nested_value(
    mapping: Mapping[str, Any],
    key_paths: Sequence[Sequence[str]],
) -> Any | None:
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


def _require_nested_value(
    mapping: Mapping[str, Any],
    key_paths: Sequence[Sequence[str]],
    description: str,
) -> Any:
    value = _first_nested_value(mapping, key_paths)
    assert value is not None, description
    return value


def _complex_from_value(value: Any) -> complex:
    if isinstance(value, bool):
        raise TypeError("boolean values do not represent Wilson coefficients")
    if isinstance(value, (int, float)):
        return complex(float(value), 0.0)
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        if len(value) == 2 and all(isinstance(item, (int, float)) for item in value):
            return complex(float(value[0]), float(value[1]))
    if isinstance(value, Mapping):
        if {"re", "im"} <= set(value):
            return complex(float(value["re"]), float(value["im"]))
        if {"real", "imag"} <= set(value):
            return complex(float(value["real"]), float(value["imag"]))
        for key in ("value", "coefficient", "wilson"):
            if key in value:
                return _complex_from_value(value[key])
    raise TypeError(f"unsupported coefficient encoding: {value!r}")


def _coefficient_complex_map_from_payload(payload: Mapping[str, Any]) -> dict[str, complex]:
    data = _require_nested_value(
        payload,
        (
            ("coefficients",),
            ("evolved_coefficients",),
            ("wilsons", "coefficients"),
            ("output_wilsons", "coefficients"),
        ),
        "RG payload does not expose coefficients",
    )
    if not isinstance(data, Mapping):
        raise AssertionError("RG coefficient payload must be a mapping")
    return {str(key): _complex_from_value(value) for key, value in data.items()}


def _coefficient_complex_map(value: Any) -> dict[str, complex]:
    coefficients = getattr(value, "coefficients", None)
    if isinstance(coefficients, Mapping):
        return {str(key): complex(item) for key, item in coefficients.items()}
    return _coefficient_complex_map_from_payload(_payload_from_value(value))


def _rg_contract_wilson_scheme_id(payload: Mapping[str, Any]) -> str:
    return str(
        _require_nested_value(
            payload,
            (("contract", "renormalization_scheme_id"),),
            "missing RG contract Wilson renormalization scheme tag",
        )
    )


def _rg_contract_rg_scheme_id(payload: Mapping[str, Any]) -> str:
    return str(
        _require_nested_value(
            payload,
            (("contract", "rg_scheme_id"),),
            "missing RG contract RG-scheme tag",
        )
    )


def _output_wilson_scheme_id(payload: Mapping[str, Any]) -> str:
    return str(
        _require_nested_value(
            payload,
            (
                ("scheme_id",),
                ("renormalization_scheme_id",),
                ("matching_summary", "scheme_id"),
                ("matching_summary", "renormalization_scheme_id"),
                ("matching_summary", "tags", "renormalization_scheme_id"),
                ("tags", "renormalization_scheme_id"),
                ("wilsons", "tags", "renormalization_scheme_id"),
            ),
            "missing RG output Wilson renormalization scheme tag",
        )
    )


def _output_rg_scheme_id(payload: Mapping[str, Any]) -> str:
    return str(
        _require_nested_value(
            payload,
            (
                ("rg_scheme_id",),
                ("evolution", "rg_scheme_id"),
            ),
            "missing RG output RG-scheme tag",
        )
    )


def _rg_contract_matching_input_scheme_id(payload: Mapping[str, Any]) -> str:
    return str(
        _require_nested_value(
            payload,
            (("contract", "matching_input_scheme_id"),),
            "missing RG contract matching-input scheme tag",
        )
    )


def _rg_contract_matching_input_matching_id(payload: Mapping[str, Any]) -> str:
    return str(
        _require_nested_value(
            payload,
            (("contract", "matching_input_matching_id"),),
            "missing RG contract matching-input matching_id tag",
        )
    )


def _rg_contract_matching_input_wilson_schema_id(payload: Mapping[str, Any]) -> str:
    return str(
        _require_nested_value(
            payload,
            (("contract", "matching_input_wilson_schema_id"),),
            "missing RG contract matching-input Wilson schema tag",
        )
    )


def _rg_contract_lr_basis_contract_id(payload: Mapping[str, Any]) -> str:
    return str(
        _require_nested_value(
            payload,
            (("contract", "lr_basis_contract_id"),),
            "missing RG contract LR-basis contract tag",
        )
    )


def _rg_contract_lr_basis_status_id(payload: Mapping[str, Any]) -> str:
    return str(
        _require_nested_value(
            payload,
            (("contract", "lr_basis_status_id"),),
            "missing RG contract LR-basis status tag",
        )
    )


def _output_lr_basis_status_id(payload: Mapping[str, Any]) -> str:
    return str(
        _require_nested_value(
            payload,
            (
                ("lr_basis_status_id",),
                ("evolution", "lr_basis_status_id"),
            ),
            "missing RG output LR-basis status tag",
        )
    )


def _rg_contract_supported_operator_subset_id(payload: Mapping[str, Any]) -> str:
    return str(
        _require_nested_value(
            payload,
            (("contract", "supported_operator_subset_id"),),
            "missing RG contract supported-operator-subset tag",
        )
    )


def _rg_supported_operator_names(payload: Mapping[str, Any]) -> tuple[str, ...]:
    values = _require_nested_value(
        payload,
        (
            ("contract", "supported_operator_names"),
            ("evolution", "supported_operator_names"),
        ),
        "missing RG supported-operator names",
    )
    assert isinstance(values, Sequence) and not isinstance(
        values,
        (str, bytes, bytearray),
    ), "RG supported_operator_names must be a sequence"
    return tuple(str(value) for value in values)


def _rg_unsupported_operator_names(payload: Mapping[str, Any]) -> tuple[str, ...]:
    values = _require_nested_value(
        payload,
        (
            ("contract", "unsupported_operator_names"),
            ("evolution", "unsupported_operator_names"),
        ),
        "missing RG unsupported-operator names",
    )
    assert isinstance(values, Sequence) and not isinstance(
        values,
        (str, bytes, bytearray),
    ), "RG unsupported_operator_names must be a sequence"
    return tuple(str(value) for value in values)


def _lr_basis_contract_map_direction_id(payload: Mapping[str, Any]) -> str:
    return str(
        _require_nested_value(
            payload,
            (("map_direction_id",),),
            "missing LR-basis map-direction tag",
        )
    )


def _lr_basis_contract_bmu_basis_id(payload: Mapping[str, Any]) -> str:
    return str(
        _require_nested_value(
            payload,
            (("bmu_lr_basis_id",),),
            "missing LR-basis BMU-basis tag",
        )
    )


def _lr_basis_contract_definition_ids(payload: Mapping[str, Any]) -> tuple[str, ...]:
    values = _require_nested_value(
        payload,
        (("paper_operator_definition_ids",),),
        "missing LR-basis definition-id tags",
    )
    assert isinstance(values, Sequence) and not isinstance(
        values,
        (str, bytes, bytearray),
    ), "LR-basis paper_operator_definition_ids must be a sequence"
    return tuple(str(value) for value in values)


def _lr_basis_contract_definitions_frozen(payload: Mapping[str, Any]) -> bool:
    value = _require_nested_value(
        payload,
        (("lr_definitions_frozen",),),
        "missing LR-basis definitions-frozen flag",
    )
    assert isinstance(value, bool), "LR-basis lr_definitions_frozen must be boolean"
    return value


def _lr_basis_contract_mapping_matrix_frozen(payload: Mapping[str, Any]) -> bool:
    value = _require_nested_value(
        payload,
        (("mapping_matrix_frozen",),),
        "missing LR-basis mapping-matrix-frozen flag",
    )
    assert isinstance(value, bool), "LR-basis mapping_matrix_frozen must be boolean"
    return value


def _lr_basis_contract_running_activated(payload: Mapping[str, Any]) -> bool:
    value = _require_nested_value(
        payload,
        (("lr_running_activated",),),
        "missing LR-basis running-activated flag",
    )
    assert isinstance(value, bool), "LR-basis lr_running_activated must be boolean"
    return value


def _matching_summary_scheme_id(payload: Mapping[str, Any]) -> str:
    return str(
        _require_nested_value(
            payload,
            (
                ("matching_summary", "scheme_id"),
                ("matching_summary", "renormalization_scheme_id"),
                ("matching_summary", "tags", "renormalization_scheme_id"),
                ("input_matching_summary", "scheme_id"),
                ("input_matching_summary", "renormalization_scheme_id"),
                ("input_matching_summary", "tags", "renormalization_scheme_id"),
                ("source_renormalization_scheme_id",),
                ("tags", "source_renormalization_scheme_id"),
                ("wilsons", "source_renormalization_scheme_id"),
                ("wilsons", "tags", "source_renormalization_scheme_id"),
            ),
            "missing matching-summary renormalization scheme tag",
        )
    )


def _matching_summary_matching_id(payload: Mapping[str, Any]) -> str:
    return str(
        _require_nested_value(
            payload,
            (
                ("matching_summary", "matching_id"),
                ("input_matching_summary", "matching_id"),
                ("tags", "matching_input_matching_id"),
                ("wilsons", "tags", "matching_input_matching_id"),
            ),
            "missing matching-summary matching_id tag",
        )
    )


def _matching_summary_wilson_schema_id(payload: Mapping[str, Any]) -> str:
    return str(
        _require_nested_value(
            payload,
            (
                ("matching_summary", "wilson_schema_id"),
                ("input_matching_summary", "wilson_schema_id"),
                ("source_wilson_schema_id",),
                ("tags", "source_wilson_schema_id"),
                ("wilsons", "source_wilson_schema_id"),
                ("wilsons", "tags", "source_wilson_schema_id"),
            ),
            "missing matching-summary Wilson schema tag",
        )
    )


def _assert_complex_maps_close(
    observed: Mapping[str, complex],
    expected: Mapping[str, complex],
    *,
    rel: float = 1.0e-10,
    abs_tol: float = 1.0e-12,
) -> None:
    assert set(observed) == set(expected)
    for name in sorted(expected):
        assert observed[name].real == pytest.approx(expected[name].real, rel=rel, abs=abs_tol)
        assert observed[name].imag == pytest.approx(expected[name].imag, rel=rel, abs=abs_tol)


def _default_matching_object() -> Any:
    matching_module = importlib.import_module(
        "quarkConstraints.paper_0710_1869.eft_deltaf2.matching_kkgluon"
    )
    callable_obj = _get_callable(matching_module, MATCHING_OBJECT_NAMES)
    assert callable_obj is not None, "matching layer exposes no default object callable"
    return callable_obj()


def _matching_wilsons(value: Any) -> Any:
    return getattr(value, "wilsons", value)


def _matching_scale_GeV(value: Any) -> float:
    if hasattr(value, "matching_scale_GeV"):
        return float(value.matching_scale_GeV)
    payload = _payload_from_value(value)
    scale = _require_nested_value(
        payload,
        (
            ("matching_scale_GeV",),
            ("mu_match_GeV",),
            ("wilsons", "matching_scale_GeV"),
        ),
        "missing matching-scale tag",
    )
    return float(scale)


def _invoke_rg_evolution(callable_obj: Any, *, value: Any, mu_low_GeV: float) -> Any:
    parameters = inspect.signature(callable_obj).parameters
    kwargs: dict[str, Any] = {}
    wilsons = _matching_wilsons(value)

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
            "RG evolution callable requires unsupported arguments for acceptance tests: "
            + ", ".join(missing)
        )
    return callable_obj(**kwargs)


def _default_rg_payload(module: Any) -> dict[str, Any]:
    summary_fn = _get_callable(module, DEFAULT_RG_SUMMARY_NAMES)
    if callable(summary_fn):
        return _payload_from_value(summary_fn())

    export_fn = _get_callable(module, DEFAULT_RG_EXPORT_NAMES)
    if callable(export_fn):
        return _payload_from_value(export_fn())

    evolve_fn = _get_callable(module, PARAMETERIZED_RG_NAMES)
    if callable(evolve_fn):
        return _payload_from_value(
            _invoke_rg_evolution(
                evolve_fn,
                value=_default_matching_object(),
                mu_low_GeV=RG_DEFAULT_LOW_SCALE_GEV,
            )
        )

    raise AssertionError(
        "RG layer exists but exposes no default export; expected one of "
        + ", ".join(DEFAULT_RG_SUMMARY_NAMES + DEFAULT_RG_EXPORT_NAMES + PARAMETERIZED_RG_NAMES)
    )


def _default_lr_basis_contract_payload(module: Any) -> dict[str, Any]:
    contract_fn = getattr(module, "default_paper_0710_1869_deltaf2_lr_basis_contract", None)
    assert callable(contract_fn), "RG module exposes no default LR-basis contract export"
    return _payload_from_value(contract_fn())


def _default_rg_result_payload(module: Any) -> dict[str, Any]:
    evolve_fn = _get_callable(module, PARAMETERIZED_RG_NAMES)
    assert callable(evolve_fn), "RG module exposes no evolution callable"
    return _payload_from_value(
        _invoke_rg_evolution(
            evolve_fn,
            value=_default_matching_object(),
            mu_low_GeV=RG_DEFAULT_LOW_SCALE_GEV,
        )
    )


def _default_rg_payload_cross_process() -> dict[str, Any]:
    script = f"""
import dataclasses
import importlib
import inspect
import json
from collections.abc import Mapping, Sequence
from pathlib import Path

SUMMARY_NAMES = {DEFAULT_RG_SUMMARY_NAMES!r}
EXPORT_NAMES = {DEFAULT_RG_EXPORT_NAMES!r}
EVOLUTION_NAMES = {PARAMETERIZED_RG_NAMES!r}
MATCHING_NAMES = {MATCHING_OBJECT_NAMES!r}
MU_LOW_GEV = {RG_DEFAULT_LOW_SCALE_GEV!r}

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
    if not isinstance(payload, dict):
        raise AssertionError("RG payload must canonicalize to a mapping")
    return payload

def get_callable(module, names):
    for name in names:
        candidate = getattr(module, name, None)
        if callable(candidate):
            return candidate
    return None

def matching_wilsons(value):
    return getattr(value, "wilsons", value)

def invoke_evolution(callable_obj, value, mu_low_GeV):
    parameters = inspect.signature(callable_obj).parameters
    kwargs = {{}}
    wilsons = matching_wilsons(value)
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
        raise AssertionError(
            "RG evolution callable requires unsupported arguments: " + ", ".join(missing)
        )
    return callable_obj(**kwargs)

rg_module = importlib.import_module("quarkConstraints.paper_0710_1869.eft_deltaf2.rg")
summary_fn = get_callable(rg_module, SUMMARY_NAMES)
if callable(summary_fn):
    print(json.dumps(payload_from_value(summary_fn()), sort_keys=True))
else:
    export_fn = get_callable(rg_module, EXPORT_NAMES)
    if callable(export_fn):
        print(json.dumps(payload_from_value(export_fn()), sort_keys=True))
    else:
        evolve_fn = get_callable(rg_module, EVOLUTION_NAMES)
        if not callable(evolve_fn):
            raise SystemExit("missing default RG export")
        matching_module = importlib.import_module(
            "quarkConstraints.paper_0710_1869.eft_deltaf2.matching_kkgluon"
        )
        matching_fn = get_callable(matching_module, MATCHING_NAMES)
        if not callable(matching_fn):
            raise SystemExit("missing default matching object export")
        print(
            json.dumps(
                payload_from_value(invoke_evolution(evolve_fn, matching_fn(), MU_LOW_GEV)),
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


def _toy_wilsons(
    *,
    matching_scale_GeV: float,
    q1_vll: complex = 0.0,
    q1_vrr: complex = 0.0,
    q4_lr: complex = 0.0,
    q5_lr: complex = 0.0,
):
    operators_module = importlib.import_module(
        "quarkConstraints.paper_0710_1869.eft_deltaf2.operators"
    )
    return operators_module.Paper07101869DeltaF2WilsonCoefficients(
        contract=operators_module.default_paper_0710_1869_deltaf2_wilson_contract(),
        benchmark_id="toy_rg_contract",
        scale_label=f"toy_mu_{matching_scale_GeV:g}",
        system_id="kaon",
        sector_id="down",
        generations=(0, 1),
        matching_scale_GeV=float(matching_scale_GeV),
        propagator_mass_GeV=3000.0,
        left_coupling=0.0,
        right_coupling=0.0,
        q1_vll=q1_vll,
        q1_vrr=q1_vrr,
        q4_lr=q4_lr,
        q5_lr=q5_lr,
    )


def _vll_lo_factor(mu_high_GeV: float, mu_low_GeV: float, *, n_f: int) -> float:
    exponent = 6.0 / (33.0 - (2.0 * float(n_f)))
    alpha_high = alpha_s(mu_high_GeV, n_loops=1, matching_loops=0)
    alpha_low = alpha_s(mu_low_GeV, n_loops=1, matching_loops=0)
    return float((alpha_low / alpha_high) ** exponent)


def test_importing_rg_module_does_not_load_repo_v1_modules() -> None:
    if not _has_rg_module():
        pytest.skip("paper_0710_1869 LO RG layer not implemented yet")

    script = f"""
import importlib
import json
import sys

importlib.import_module("quarkConstraints.paper_0710_1869.eft_deltaf2.rg")

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


def test_default_rg_export_is_deterministic_and_tagged() -> None:
    module = _load_rg_module()

    first = _default_rg_payload(module)
    second = _default_rg_payload(module)
    third = _default_rg_payload_cross_process()
    assert json.dumps(first, sort_keys=True) == json.dumps(second, sort_keys=True)
    assert json.dumps(first, sort_keys=True) == json.dumps(third, sort_keys=True)

    basis_id = _require_nested_value(
        first,
        (
            ("basis_id",),
            ("operator_basis_id",),
            ("contract", "operator_basis_id"),
            ("tags", "operator_basis_id"),
            ("wilsons", "tags", "operator_basis_id"),
        ),
        "missing operator-basis tag",
    )
    scheme_id = _output_wilson_scheme_id(first)
    contract_wilson_scheme_id = _rg_contract_wilson_scheme_id(first)
    rg_scheme_id = _output_rg_scheme_id(first)
    contract_rg_scheme_id = _rg_contract_rg_scheme_id(first)
    rg_order = _require_nested_value(
        first,
        (
            ("rg_order_id",),
            ("rg_order",),
            ("contract", "rg_order_id"),
            ("evolution", "rg_order_id"),
        ),
        "missing RG-order tag",
    )
    alpha_policy_id = _require_nested_value(
        first,
        (
            ("alpha_s_policy_id",),
            ("contract", "alpha_s_policy_id"),
            ("evolution", "alpha_s_policy_id"),
        ),
        "missing alpha_s policy tag",
    )
    threshold_policy_id = _require_nested_value(
        first,
        (
            ("threshold_policy_id",),
            ("contract", "threshold_policy_id"),
            ("evolution", "threshold_policy_id"),
        ),
        "missing threshold policy tag",
    )
    low_scale = float(
        _require_nested_value(
            first,
            (
                ("mu_low_GeV",),
                ("evaluation_scale_GeV",),
                ("running_scale_GeV",),
                ("scale_GeV",),
                ("wilsons", "matching_scale_GeV"),
                ("matching_scale_GeV",),
            ),
            "missing RG low-scale tag",
        )
    )
    segments = _require_nested_value(
        first,
        (
            ("segments",),
            ("evolution", "segments"),
        ),
        "missing RG threshold segments",
    )
    coefficient_map = _coefficient_complex_map_from_payload(first)

    assert str(basis_id)
    assert str(scheme_id)
    assert str(contract_wilson_scheme_id)
    assert str(rg_scheme_id)
    assert str(contract_rg_scheme_id)
    assert str(alpha_policy_id)
    assert str(threshold_policy_id)
    assert str(rg_order).lower() == "lo"
    assert low_scale > 0.0
    assert isinstance(segments, Sequence)
    assert not isinstance(segments, (str, bytes, bytearray))
    assert coefficient_map


def test_default_rg_export_scheme_tags_are_contract_coherent() -> None:
    payload = _default_rg_payload(_load_rg_module())

    assert _output_wilson_scheme_id(payload) == _rg_contract_wilson_scheme_id(payload)
    assert _output_rg_scheme_id(payload) == _rg_contract_rg_scheme_id(payload)


def test_default_rg_export_preserves_frozen_pr3_input_source_tags() -> None:
    payload = _default_rg_payload(_load_rg_module())

    assert _matching_summary_scheme_id(payload) == _rg_contract_matching_input_scheme_id(payload)
    assert _matching_summary_matching_id(payload) == _rg_contract_matching_input_matching_id(
        payload
    )
    assert _matching_summary_wilson_schema_id(payload) == (
        _rg_contract_matching_input_wilson_schema_id(payload)
    )


def test_default_rg_export_carries_frozen_lr_contract_ids() -> None:
    module = _load_rg_module()
    payload = _default_rg_result_payload(module)
    lr_contract_payload = _default_lr_basis_contract_payload(module)
    expected_contract_id, expected_status_id, expected_subset_id = _rg_inputs_ids()

    assert _rg_contract_lr_basis_contract_id(payload) == EXPECTED_RG_LR_BASIS_CONTRACT_ID
    assert _rg_contract_lr_basis_contract_id(payload) == expected_contract_id
    assert _rg_contract_lr_basis_contract_id(payload) == str(lr_contract_payload["contract_id"])
    assert _rg_contract_lr_basis_status_id(payload) == expected_status_id
    assert _rg_contract_lr_basis_status_id(payload) == str(lr_contract_payload["status_id"])
    assert _output_lr_basis_status_id(payload) == expected_status_id
    assert _rg_contract_supported_operator_subset_id(payload) == expected_subset_id
    assert _rg_supported_operator_names(payload) == EXPECTED_FULL_RG_SUPPORTED_OPERATORS
    assert _rg_unsupported_operator_names(payload) == tuple()
    assert _first_nested_value(
        payload,
        (
            ("lr_basis_map_supported",),
            ("evolution", "lr_basis_map_supported"),
        ),
    ) is True
    assert _lr_basis_contract_map_direction_id(lr_contract_payload) == (
        EXPECTED_RG_LR_BASIS_DIRECTION_ID
    )
    assert _lr_basis_contract_bmu_basis_id(lr_contract_payload) == (
        EXPECTED_RG_LR_BMU_BASIS_ID
    )
    assert _lr_basis_contract_definition_ids(lr_contract_payload) == (
        EXPECTED_PAPER_LR_OPERATOR_DEFINITION_IDS
    )
    assert _lr_basis_contract_definitions_frozen(lr_contract_payload) is True
    assert _lr_basis_contract_mapping_matrix_frozen(lr_contract_payload) is True
    assert _lr_basis_contract_running_activated(lr_contract_payload) is True


def test_lo_rg_identity_holds_at_same_scale() -> None:
    module = _load_rg_module()
    evolve_fn = _get_callable(module, PARAMETERIZED_RG_NAMES)
    assert callable(evolve_fn), "RG module exposes no evolution callable"

    matching = _default_matching_object()
    wilsons = _matching_wilsons(matching)
    mu_same = _matching_scale_GeV(wilsons)
    evolved = _invoke_rg_evolution(evolve_fn, value=matching, mu_low_GeV=mu_same)

    _assert_complex_maps_close(
        _coefficient_complex_map(evolved),
        _coefficient_complex_map(wilsons),
        rel=0.0,
        abs_tol=1.0e-14,
    )


def test_lo_rg_composition_matches_two_step_running() -> None:
    module = _load_rg_module()
    evolve_fn = _get_callable(module, PARAMETERIZED_RG_NAMES)
    assert callable(evolve_fn), "RG module exposes no evolution callable"

    matching = _default_matching_object()
    direct = _invoke_rg_evolution(evolve_fn, value=matching, mu_low_GeV=2.0)
    mid = _invoke_rg_evolution(evolve_fn, value=matching, mu_low_GeV=25.0)
    composed = _invoke_rg_evolution(evolve_fn, value=mid, mu_low_GeV=2.0)

    _assert_complex_maps_close(
        _coefficient_complex_map(composed),
        _coefficient_complex_map(direct),
        rel=1.0e-9,
        abs_tol=1.0e-12,
    )


def test_lo_rg_rejects_mutated_pr3_matching_id() -> None:
    module = _load_rg_module()
    evolve_fn = _get_callable(module, PARAMETERIZED_RG_NAMES)
    assert callable(evolve_fn), "RG module exposes no evolution callable"

    operators_module = importlib.import_module(
        "quarkConstraints.paper_0710_1869.eft_deltaf2.operators"
    )
    matching = _default_matching_object()
    wilsons = _matching_wilsons(matching)
    bad_basis = dataclasses.replace(
        wilsons.contract.operator_basis,
        matching_id="eft_matching.lo.mutated_for_negative_test.v1",
    )
    bad_contract = dataclasses.replace(wilsons.contract, operator_basis=bad_basis)
    bad_wilsons = dataclasses.replace(wilsons, contract=bad_contract)

    assert bad_wilsons.contract.operator_basis_id == wilsons.contract.operator_basis_id
    assert (
        bad_wilsons.contract.renormalization_scheme_id
        == wilsons.contract.renormalization_scheme_id
    )
    assert (
        bad_wilsons.contract.operator_normalization_id
        == wilsons.contract.operator_normalization_id
    )
    assert bad_wilsons.contract.operator_order == wilsons.contract.operator_order
    assert bad_wilsons.contract.matching_id != wilsons.contract.matching_id
    assert bad_wilsons.schema_id == operators_module.PAPER_0710_1869_DELTAF2_WILSON_SCHEMA_ID

    with pytest.raises(ValueError, match="matching_id is incompatible"):
        _invoke_rg_evolution(evolve_fn, value=bad_wilsons, mu_low_GeV=2.0)


def test_lo_rg_threshold_segmentation_matches_piecewise_running() -> None:
    module = _load_rg_module()
    evolve_fn = _get_callable(module, PARAMETERIZED_RG_NAMES)
    assert callable(evolve_fn), "RG module exposes no evolution callable"

    wilsons = _toy_wilsons(matching_scale_GeV=3000.0, q1_vll=1.0)
    direct = _invoke_rg_evolution(evolve_fn, value=wilsons, mu_low_GeV=1.0)
    step_top = _invoke_rg_evolution(evolve_fn, value=wilsons, mu_low_GeV=M_TOP_MS)
    step_bottom = _invoke_rg_evolution(evolve_fn, value=step_top, mu_low_GeV=M_BOTTOM)
    step_charm = _invoke_rg_evolution(evolve_fn, value=step_bottom, mu_low_GeV=M_CHARM)
    piecewise = _invoke_rg_evolution(evolve_fn, value=step_charm, mu_low_GeV=1.0)

    _assert_complex_maps_close(
        _coefficient_complex_map(piecewise),
        _coefficient_complex_map(direct),
        rel=2.0e-8,
        abs_tol=1.0e-12,
    )


def test_vll_vrr_fixed_nf_lo_regression_matches_bmu_exponent() -> None:
    module = _load_rg_module()
    evolve_fn = _get_callable(module, PARAMETERIZED_RG_NAMES)
    assert callable(evolve_fn), "RG module exposes no evolution callable"

    mu_high = 100.0
    mu_low = 10.0
    factor = _vll_lo_factor(mu_high, mu_low, n_f=5)
    wilsons = _toy_wilsons(
        matching_scale_GeV=mu_high,
        q1_vll=1.0,
        q1_vrr=1.0,
    )
    evolved = _invoke_rg_evolution(evolve_fn, value=wilsons, mu_low_GeV=mu_low)
    coefficients = _coefficient_complex_map(evolved)

    assert coefficients["Q1_VLL"].real == pytest.approx(factor, rel=1.0e-9, abs=1.0e-12)
    assert coefficients["Q1_VLL"].imag == pytest.approx(0.0, abs=1.0e-12)
    assert coefficients["Q1_VRR"].real == pytest.approx(factor, rel=1.0e-9, abs=1.0e-12)
    assert coefficients["Q1_VRR"].imag == pytest.approx(0.0, abs=1.0e-12)
    assert coefficients["Q4_LR"].real == pytest.approx(0.0, abs=1.0e-12)
    assert coefficients["Q4_LR"].imag == pytest.approx(0.0, abs=1.0e-12)
    assert coefficients["Q5_LR"].real == pytest.approx(0.0, abs=1.0e-12)
    assert coefficients["Q5_LR"].imag == pytest.approx(0.0, abs=1.0e-12)


def test_lo_rg_accepts_nonzero_lr_and_evolves_it() -> None:
    module = _load_rg_module()
    evolve_fn = _get_callable(module, PARAMETERIZED_RG_NAMES)
    assert callable(evolve_fn), "RG module exposes no evolution callable"

    wilsons = _toy_wilsons(
        matching_scale_GeV=100.0,
        q4_lr=1.0 + 0.25j,
        q5_lr=-0.5 + 0.125j,
    )
    evolved = _invoke_rg_evolution(evolve_fn, value=wilsons, mu_low_GeV=10.0)
    coefficients = _coefficient_complex_map(evolved)

    assert set(coefficients) == {"Q1_VLL", "Q1_VRR", "Q4_LR", "Q5_LR"}
    assert abs(coefficients["Q1_VLL"]) == pytest.approx(0.0, abs=1.0e-18)
    assert abs(coefficients["Q1_VRR"]) == pytest.approx(0.0, abs=1.0e-18)
    assert abs(coefficients["Q4_LR"]) > 0.0
    assert abs(coefficients["Q5_LR"]) > 0.0
