"""Forward-compatible tests for the paper-owned Delta F=2 matching layer."""

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

import numpy as np
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
PACKAGE_ROOT = REPO_ROOT / "quarkConstraints" / "paper_0710_1869"
OPERATORS_MODULE_PATH = PACKAGE_ROOT / "eft_deltaf2" / "operators.py"
MATCHING_MODULE_PATH = PACKAGE_ROOT / "eft_deltaf2" / "matching_kkgluon.py"

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
IMPLEMENTATION_DEFAULT_MATCHING_EXPORT_NAMES = (
    "match_default_paper_0710_1869_kaon_kk_gluon_deltaf2",
)
MATCHING_ALIAS_STEMS = ("kaon_matching", "deltaf2_matching", "kkgluon_matching")
DEFAULT_MATCHING_SUMMARY_NAMES = tuple(
    name
    for stem in MATCHING_ALIAS_STEMS
    for name in (
        f"default_paper_0710_1869_{stem}_summary",
        f"build_paper_0710_1869_{stem}_summary",
    )
)
DEFAULT_MATCHING_EXPORT_NAMES = IMPLEMENTATION_DEFAULT_MATCHING_EXPORT_NAMES + tuple(
    name
    for stem in MATCHING_ALIAS_STEMS
    for name in (
        f"default_paper_0710_1869_{stem}",
        f"build_paper_0710_1869_{stem}",
    )
)
PARAMETERIZED_MATCHING_NAMES = (
    "match_paper_0710_1869_kk_gluon_deltaf2",
    "match_kkgluon_to_deltaf2",
    "match_paper_0710_1869_kkgluon_to_deltaf2",
)
EXPECTED_Q4_LR_DEFINITION_ID = (
    "paper_0710_1869.deltaf2.q4_lr.susy_o4.scalar_lr.color_singlet.v1"
)
EXPECTED_Q5_LR_DEFINITION_ID = (
    "paper_0710_1869.deltaf2.q5_lr.susy_o5.scalar_lr.color_mixed.v1"
)
EXPECTED_PROJECTOR_NORMALIZATION_ID = (
    "paper_0710_1869.deltaf2.projectors.pl_pr_equals_1mp_gamma5_over_2.v1"
)
EXPECTED_PROJECTOR_NORMALIZATION_NOTE = (
    "All chiral bilinears use P_L = (1-gamma5)/2 and P_R = (1+gamma5)/2. "
    "The legacy susy_o4/susy_o5 identifiers are retained for compatibility only. "
    "They label the O4/O5 scalar LR basis used in the DeltaF=2 literature and do "
    "not imply a bare (1-gamma5) or (1+gamma5) normalization without the "
    "factor-of-two projector convention."
)


def _has_matching_module() -> bool:
    return MATCHING_MODULE_PATH.exists()


def _load_matching_module():
    if not _has_matching_module():
        pytest.skip("paper_0710_1869 EFT matching layer not implemented yet")
    assert OPERATORS_MODULE_PATH.exists(), (
        "matching_kkgluon.py exists but eft_deltaf2/operators.py is missing"
    )
    return importlib.import_module(
        "quarkConstraints.paper_0710_1869.eft_deltaf2.matching_kkgluon"
    )


def _load_operators_module():
    if not OPERATORS_MODULE_PATH.exists():
        pytest.skip("paper_0710_1869 EFT operator basis not implemented yet")
    return importlib.import_module("quarkConstraints.paper_0710_1869.eft_deltaf2.operators")


def _physical_bulk_state() -> Any:
    benchmarks_module = importlib.import_module("quarkConstraints.paper_0710_1869.benchmarks")
    model_module = importlib.import_module("quarkConstraints.paper_0710_1869.model")

    seed = benchmarks_module.Paper07101869BenchmarkSpurionSeed(
        up_singular_values=(0.45, 1.2, 3.4),
        down_singular_values=(0.15, 0.55, 1.7),
        overall_scale=0.8,
        up_left=model_module.Paper07101869RotationParameters.from_degrees(
            theta12_deg=19.0,
            theta13_deg=4.0,
            theta23_deg=11.0,
            delta=0.3,
        ),
        up_right=model_module.Paper07101869RotationParameters.from_degrees(
            theta12_deg=13.0,
            theta13_deg=7.0,
            theta23_deg=5.0,
            delta=0.1,
        ),
        down_left=model_module.Paper07101869RotationParameters.from_degrees(
            theta12_deg=7.0,
            theta13_deg=2.0,
            theta23_deg=9.0,
            delta=0.5,
        ),
        down_right=model_module.Paper07101869RotationParameters.from_degrees(
            theta12_deg=17.0,
            theta13_deg=3.0,
            theta23_deg=6.0,
            delta=0.2,
        ),
        notes="matching-test-seed",
    )
    benchmark = benchmarks_module.default_paper_0710_1869_pr1_benchmark()
    physical_point = benchmarks_module.build_paper_0710_1869_seeded_physical_point(
        benchmark,
        seed,
        metadata={"test_case": "matching_acceptance"},
    )
    return model_module.derive_paper_0710_1869_physical_bulk_state(physical_point)


def _physical_kk_gluon_adapter_callable(module: Any) -> Any:
    for name in (
        "build_paper_0710_1869_physical_kk_gluon_couplings",
        "build_paper_0710_1869_kk_gluon_couplings_from_physical_bulk_state",
        "build_paper_0710_1869_kk_gluon_couplings_from_physical_state",
    ):
        candidate = getattr(module, name, None)
        if callable(candidate):
            return candidate

    candidate = getattr(module, "build_paper_0710_1869_kk_gluon_couplings", None)
    if callable(candidate):
        try:
            parameters = inspect.signature(candidate).parameters
        except (TypeError, ValueError):
            parameters = {}
        if any(name in parameters for name in ("physical_bulk_state", "bulk_state")):
            return candidate
    return None


def _build_physical_kk_gluon_couplings() -> Any:
    kkgluon_module = importlib.import_module("quarkConstraints.paper_0710_1869.kkgluon")
    adapter = _physical_kk_gluon_adapter_callable(kkgluon_module)
    if not callable(adapter):
        pytest.skip("physical bulk-state KK-gluon adapter not exposed yet")

    physical_bulk_state = _physical_bulk_state()
    try:
        parameters = inspect.signature(adapter).parameters
    except (TypeError, ValueError):
        parameters = {}
    for name in ("physical_bulk_state", "bulk_state", "physical_state"):
        if name in parameters:
            return adapter(**{name: physical_bulk_state})
    pytest.skip("physical bulk-state KK-gluon adapter exposes no supported keyword parameter")


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
        return _canonicalize(value.as_dict())
    payload = _canonicalize(value)
    if not isinstance(payload, dict):
        raise AssertionError("matching payload must canonicalize to a mapping")
    return payload


def _get_callable(module: Any, names: Sequence[str]) -> Any:
    for name in names:
        candidate = getattr(module, name, None)
        if callable(candidate):
            return candidate
    return None


def _first_mapping_value(mapping: Mapping[str, Any], keys: Sequence[str]) -> Any | None:
    for key in keys:
        if key in mapping:
            return mapping[key]
    return None


def _require_mapping_value(
    mapping: Mapping[str, Any],
    keys: Sequence[str],
    description: str,
) -> Any:
    found = _first_mapping_value(mapping, keys)
    assert found is not None, f"matching payload missing {description}: {', '.join(keys)}"
    return found


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


def _coefficient_abs_map(payload: Mapping[str, Any]) -> dict[str, float]:
    data = _first_mapping_value(payload, ("coefficients", "wilson_coefficients", "C_i", "C"))
    assert data is not None, "matching payload does not expose C_i(mu_match)"
    if isinstance(data, Mapping):
        return {str(key): abs(_complex_from_value(value)) for key, value in data.items()}
    if isinstance(data, Sequence) and not isinstance(data, (str, bytes, bytearray)):
        rows: dict[str, float] = {}
        for item in data:
            if not isinstance(item, Mapping):
                raise AssertionError("coefficient row payload must contain mappings")
            label = _require_mapping_value(
                item,
                ("operator_id", "name", "label"),
                "operator identifier",
            )
            coefficient = _require_mapping_value(
                item,
                ("value", "coefficient", "wilson"),
                "coefficient value",
            )
            rows[str(label)] = abs(_complex_from_value(coefficient))
        return rows
    raise AssertionError("unsupported C_i(mu_match) payload shape")


def _default_matching_payload(module: Any) -> dict[str, Any]:
    summary_fn = _get_callable(module, DEFAULT_MATCHING_SUMMARY_NAMES)
    if callable(summary_fn):
        return _payload_from_value(summary_fn())
    export_fn = _get_callable(module, DEFAULT_MATCHING_EXPORT_NAMES)
    if callable(export_fn):
        return _payload_from_value(export_fn())
    raise AssertionError(
        "matching layer exists but exposes no default benchmark export; expected one of "
        f"{', '.join(DEFAULT_MATCHING_SUMMARY_NAMES + DEFAULT_MATCHING_EXPORT_NAMES)}"
    )


def _default_matching_payload_cross_process() -> dict[str, Any]:
    script = f"""
import dataclasses
import importlib
import json
from collections.abc import Mapping, Sequence
from pathlib import Path

SUMMARY_NAMES = {DEFAULT_MATCHING_SUMMARY_NAMES!r}
EXPORT_NAMES = {DEFAULT_MATCHING_EXPORT_NAMES!r}

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
        return canonicalize(value.as_dict())
    payload = canonicalize(value)
    if not isinstance(payload, dict):
        raise AssertionError("matching payload must canonicalize to a mapping")
    return payload

module = importlib.import_module("quarkConstraints.paper_0710_1869.eft_deltaf2.matching_kkgluon")
for name in SUMMARY_NAMES + EXPORT_NAMES:
    candidate = getattr(module, name, None)
    if callable(candidate):
        print(json.dumps(payload_from_value(candidate()), sort_keys=True))
        break
else:
    raise SystemExit("missing default matching export")
"""
    completed = subprocess.run(
        [sys.executable, "-c", script],
        check=True,
        capture_output=True,
        text=True,
        cwd=str(REPO_ROOT),
    )
    return json.loads(completed.stdout)


def _invoke_matching(callable_obj: Any, *, kk_couplings: Any) -> Any:
    parameters = inspect.signature(callable_obj).parameters
    kwargs: dict[str, Any] = {}

    for name in ("kk_couplings", "kk_gluon_couplings", "couplings"):
        if name in parameters:
            kwargs[name] = kk_couplings
            break

    for name in ("system", "system_id", "meson", "meson_id"):
        if name in parameters and parameters[name].default is inspect._empty:
            kwargs[name] = "kaon"

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
            "matching callable requires unsupported arguments for acceptance tests: "
            + ", ".join(missing)
        )
    return callable_obj(**kwargs)


def test_importing_matching_module_does_not_load_repo_v1_modules() -> None:
    if not _has_matching_module():
        pytest.skip("paper_0710_1869 EFT matching layer not implemented yet")

    script = f"""
import importlib
import json
import sys

importlib.import_module("quarkConstraints.paper_0710_1869.eft_deltaf2.matching_kkgluon")

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


def test_default_matching_export_is_deterministic_and_tagged() -> None:
    module = _load_matching_module()

    first = _default_matching_payload(module)
    second = _default_matching_payload(module)
    third = _default_matching_payload_cross_process()
    assert json.dumps(first, sort_keys=True) == json.dumps(second, sort_keys=True)
    assert json.dumps(first, sort_keys=True) == json.dumps(third, sort_keys=True)

    basis_id = _require_mapping_value(first, ("basis_id", "operator_basis_id"), "basis tag")
    scheme_id = _require_mapping_value(
        first,
        ("scheme_id", "renormalization_scheme_id"),
        "renormalization scheme tag",
    )
    mu_match_GeV = float(
        _require_mapping_value(first, ("mu_match_GeV", "matching_scale_GeV"), "mu_match tag")
    )
    propagator_mass_GeV = float(
        _require_mapping_value(first, ("propagator_mass_GeV",), "propagator-mass tag")
    )
    coefficient_map = _coefficient_abs_map(first)

    assert str(basis_id)
    assert str(scheme_id)
    assert mu_match_GeV > 0.0
    assert propagator_mass_GeV > 0.0
    assert coefficient_map


def test_matching_scale_tags_track_mu_match_not_propagator_mass() -> None:
    module = _load_matching_module()
    callable_obj = _get_callable(module, PARAMETERIZED_MATCHING_NAMES)
    if not callable(callable_obj):
        pytest.skip("parameterized matching callable is not exposed yet")

    kkgluon_module = importlib.import_module("quarkConstraints.paper_0710_1869.kkgluon")
    scales_module = importlib.import_module("quarkConstraints.paper_0710_1869.scales")

    build_couplings = kkgluon_module.build_paper_0710_1869_kk_gluon_couplings
    Paper07101869ScalePoint = scales_module.Paper07101869ScalePoint

    base_scale = Paper07101869ScalePoint(
        label="kaon_base",
        Lambda_IR_GeV=3000.0,
        m_g1_GeV=4200.0,
        mu_match_GeV=3000.0,
        mu_gs_GeV=4200.0,
    )
    heavier_scale = Paper07101869ScalePoint(
        label="kaon_heavier",
        Lambda_IR_GeV=3000.0,
        m_g1_GeV=4200.0,
        mu_match_GeV=3000.0,
        mu_gs_GeV=4200.0,
        m_KK_eff_GeV=8400.0,
    )

    base_payload = _payload_from_value(
        _invoke_matching(callable_obj, kk_couplings=build_couplings(scale_point=base_scale))
    )
    heavier_payload = _payload_from_value(
        _invoke_matching(callable_obj, kk_couplings=build_couplings(scale_point=heavier_scale))
    )

    base_mu_match = float(
        _require_mapping_value(
            base_payload,
            ("mu_match_GeV", "matching_scale_GeV"),
            "mu_match tag",
        )
    )
    heavier_mu_match = float(
        _require_mapping_value(
            heavier_payload,
            ("mu_match_GeV", "matching_scale_GeV"),
            "mu_match tag",
        )
    )
    base_propagator_mass = float(
        _require_mapping_value(
            base_payload,
            ("propagator_mass_GeV",),
            "propagator-mass tag",
        )
    )
    heavier_propagator_mass = float(
        _require_mapping_value(
            heavier_payload,
            ("propagator_mass_GeV",),
            "propagator-mass tag",
        )
    )

    assert base_mu_match == pytest.approx(base_scale.mu_match_GeV)
    assert heavier_mu_match == pytest.approx(heavier_scale.mu_match_GeV)
    assert base_mu_match == pytest.approx(heavier_mu_match)
    assert base_propagator_mass == pytest.approx(base_scale.propagator_mass_GeV)
    assert heavier_propagator_mass == pytest.approx(heavier_scale.propagator_mass_GeV)
    assert base_propagator_mass != pytest.approx(heavier_propagator_mass)


def test_matching_accepts_physical_bulk_state_couplings_without_new_logic() -> None:
    module = _load_matching_module()
    callable_obj = _get_callable(module, PARAMETERIZED_MATCHING_NAMES)
    if not callable(callable_obj):
        pytest.skip("parameterized matching callable is not exposed yet")

    kkgluon_module = importlib.import_module("quarkConstraints.paper_0710_1869.kkgluon")
    if _physical_kk_gluon_adapter_callable(kkgluon_module) is None:
        pytest.skip("physical bulk-state KK-gluon adapter not exposed yet")

    physical_couplings = _build_physical_kk_gluon_couplings()
    assert isinstance(
        physical_couplings,
        kkgluon_module.Paper07101869PhysicalKKGluonCouplings,
    )
    assert (
        physical_couplings.physical_profile_status_id
        == kkgluon_module.PAPER_0710_1869_PHYSICAL_KK_GLUON_STATUS_ID
    )
    payload = _payload_from_value(_invoke_matching(callable_obj, kk_couplings=physical_couplings))

    assert float(
        _require_mapping_value(payload, ("mu_match_GeV", "matching_scale_GeV"), "mu_match tag")
    ) == pytest.approx(float(physical_couplings.normalization.mu_match_GeV))
    assert float(
        _require_mapping_value(payload, ("propagator_mass_GeV",), "propagator-mass tag")
    ) == pytest.approx(float(physical_couplings.normalization.propagator_mass_GeV))
    assert _coefficient_abs_map(payload)
    assert _require_mapping_value(payload, ("benchmark_id",), "benchmark identifier") == physical_couplings.benchmark_id
    assert _require_mapping_value(payload, ("scale_label",), "scale label") == physical_couplings.normalization.scale_label
    assert _require_mapping_value(payload, ("left_basis_id",), "left basis identifier") == physical_couplings.left_down_aligned.basis_id
    assert _require_mapping_value(payload, ("right_basis_id",), "right basis identifier") == physical_couplings.right_down.basis_id


def test_lr_operator_descriptors_are_explicit_and_non_placeholder() -> None:
    _load_matching_module()
    operators_module = _load_operators_module()

    basis = operators_module.default_paper_0710_1869_deltaf2_operator_basis()
    operator_map = {operator.name: operator for operator in basis.operators}
    q4_lr = operator_map["Q4_LR"]
    q5_lr = operator_map["Q5_LR"]

    assert q4_lr.definition_id == EXPECTED_Q4_LR_DEFINITION_ID
    assert q5_lr.definition_id == EXPECTED_Q5_LR_DEFINITION_ID
    assert q4_lr.lorentz_structure == "(S-P)x(S+P)"
    assert q5_lr.lorentz_structure == "(S-P)x(S+P)"
    assert q4_lr.color_structure == "color_singlet_density_density"
    assert q5_lr.color_structure == "color_mixed_density_density"
    assert "P_L" in q4_lr.operator_formula and "P_R" in q4_lr.operator_formula
    assert "P_L" in q5_lr.operator_formula and "P_R" in q5_lr.operator_formula

    for descriptor in (
        q4_lr.definition_id,
        q4_lr.lorentz_structure,
        q4_lr.color_structure,
        q4_lr.operator_formula,
        q5_lr.definition_id,
        q5_lr.lorentz_structure,
        q5_lr.color_structure,
        q5_lr.operator_formula,
    ):
        assert "placeholder" not in descriptor.lower()

    assert q5_lr.lorentz_structure != "tensor_mixed_fierz_partner"


def test_zero_fcnc_toy_couplings_give_zero_wilsons() -> None:
    module = _load_matching_module()
    callable_obj = _get_callable(module, PARAMETERIZED_MATCHING_NAMES)
    if not callable(callable_obj):
        pytest.skip("parameterized matching callable is not exposed yet")

    kkgluon_module = importlib.import_module("quarkConstraints.paper_0710_1869.kkgluon")
    couplings = kkgluon_module.default_paper_0710_1869_kk_gluon_couplings()
    original_raw = np.asarray(couplings.left_down_aligned.raw_dimensionless, dtype=np.complex128)
    zero_left_down_raw = np.diag(np.diag(original_raw)).astype(np.complex128, copy=False)
    zero_left_down_universal_component = float(np.trace(zero_left_down_raw).real / 3.0)
    zero_left_down_subtracted = zero_left_down_raw - (
        zero_left_down_universal_component * np.eye(3, dtype=np.complex128)
    )
    g_s_mu_gs = float(couplings.left_down_aligned.g_s_mu_gs)

    assert np.allclose(
        zero_left_down_raw - np.diag(np.diag(zero_left_down_raw)),
        0.0,
        atol=0.0,
    )

    zero_left_down = dataclasses.replace(
        couplings.left_down_aligned,
        label="left_down_aligned_zero_fcnc",
        basis_id=couplings.left_down_aligned.basis_id,
        raw_dimensionless=zero_left_down_raw,
        universal_component=zero_left_down_universal_component,
        universal_subtracted_dimensionless=zero_left_down_subtracted,
        raw_gs_normalized=g_s_mu_gs * zero_left_down_raw,
        universal_subtracted_gs_normalized=g_s_mu_gs * zero_left_down_subtracted,
    )
    zero_fcnc_couplings = dataclasses.replace(couplings, left_down_aligned=zero_left_down)
    payload = _payload_from_value(_invoke_matching(callable_obj, kk_couplings=zero_fcnc_couplings))

    coefficient_abs = _coefficient_abs_map(payload)
    assert coefficient_abs
    assert max(coefficient_abs.values()) <= 1.0e-14


def test_toy_formula_regression_scales_with_inverse_propagator_mass_squared() -> None:
    module = _load_matching_module()
    callable_obj = _get_callable(module, PARAMETERIZED_MATCHING_NAMES)
    if not callable(callable_obj):
        pytest.skip("parameterized matching callable is not exposed yet")

    kkgluon_module = importlib.import_module("quarkConstraints.paper_0710_1869.kkgluon")
    scales_module = importlib.import_module("quarkConstraints.paper_0710_1869.scales")

    build_couplings = kkgluon_module.build_paper_0710_1869_kk_gluon_couplings
    Paper07101869ScalePoint = scales_module.Paper07101869ScalePoint

    base_scale = Paper07101869ScalePoint(
        label="kaon_formula_base",
        Lambda_IR_GeV=3000.0,
        m_g1_GeV=4200.0,
        mu_match_GeV=3000.0,
        mu_gs_GeV=4200.0,
    )
    heavier_scale = Paper07101869ScalePoint(
        label="kaon_formula_heavier",
        Lambda_IR_GeV=3000.0,
        m_g1_GeV=4200.0,
        mu_match_GeV=3000.0,
        mu_gs_GeV=4200.0,
        m_KK_eff_GeV=8400.0,
    )

    base_payload = _payload_from_value(
        _invoke_matching(callable_obj, kk_couplings=build_couplings(scale_point=base_scale))
    )
    heavier_payload = _payload_from_value(
        _invoke_matching(callable_obj, kk_couplings=build_couplings(scale_point=heavier_scale))
    )
    base_coefficients = _coefficient_abs_map(base_payload)
    heavier_coefficients = _coefficient_abs_map(heavier_payload)

    expected_ratio = (
        base_scale.propagator_mass_GeV / heavier_scale.propagator_mass_GeV
    ) ** 2
    informative = {
        name: value
        for name, value in base_coefficients.items()
        if value > 1.0e-18 and name in heavier_coefficients
    }
    assert informative, "matching export has no nonzero coefficients for the default kaon point"

    for name, base_value in informative.items():
        heavier_value = heavier_coefficients[name]
        assert heavier_value == pytest.approx(base_value * expected_ratio, rel=1.0e-10)


def test_lr_operator_descriptors_are_frozen_and_not_placeholders() -> None:
    operators_module = _load_operators_module()
    basis = operators_module.default_paper_0710_1869_deltaf2_operator_basis()
    operators_by_name = {operator.name: operator for operator in basis.operators}

    assert basis.projector_normalization_id == EXPECTED_PROJECTOR_NORMALIZATION_ID
    assert basis.projector_normalization_note == EXPECTED_PROJECTOR_NORMALIZATION_NOTE

    q4_lr = operators_by_name[operators_module.PAPER_0710_1869_DELTAF2_Q4_LR]
    q5_lr = operators_by_name[operators_module.PAPER_0710_1869_DELTAF2_Q5_LR]

    assert (
        q4_lr.definition_id
        == operators_module.PAPER_0710_1869_DELTAF2_OPERATOR_DEFINITION_Q4_LR_ID
    )
    assert (
        q5_lr.definition_id
        == operators_module.PAPER_0710_1869_DELTAF2_OPERATOR_DEFINITION_Q5_LR_ID
    )
    assert q4_lr.chirality == "LR"
    assert q5_lr.chirality == "LR"
    assert q4_lr.lorentz_structure == "(S-P)x(S+P)"
    assert q5_lr.lorentz_structure == "(S-P)x(S+P)"
    assert q4_lr.color_structure == "color_singlet_density_density"
    assert q5_lr.color_structure == "color_mixed_density_density"
    assert "P_L" in q4_lr.operator_formula
    assert "P_R" in q4_lr.operator_formula
    assert "P_L" in q5_lr.operator_formula
    assert "P_R" in q5_lr.operator_formula

    for descriptor in (
        q4_lr.definition_id,
        q5_lr.definition_id,
        q4_lr.lorentz_structure,
        q5_lr.lorentz_structure,
        q4_lr.color_structure,
        q5_lr.color_structure,
        q4_lr.operator_formula,
        q5_lr.operator_formula,
    ):
        lowered = descriptor.lower()
        assert "placeholder" not in lowered
        assert "todo" not in lowered
        assert "tbd" not in lowered


def test_human_facing_lr_notes_use_neutral_basis_language() -> None:
    operators_module = _load_operators_module()
    matching_module = _load_matching_module()
    matching_summary_fn = _get_callable(matching_module, DEFAULT_MATCHING_SUMMARY_NAMES)
    assert callable(matching_summary_fn), "matching layer exposes no default summary callable"

    basis = operators_module.default_paper_0710_1869_deltaf2_operator_basis()
    matching_summary = _payload_from_value(matching_summary_fn())
    support_note = str(matching_summary["lr_observable_support_note"])

    assert "legacy susy_o4/susy_o5 identifiers" in basis.projector_normalization_note
    assert "stable machine ids" in basis.metadata_compatibility_note
    assert "O4/O5 scalar LR operators" in basis.notes
    assert "paper O4/O5 scalar LR basis used in the DeltaF=2 literature" in support_note

    for text in (
        basis.projector_normalization_note,
        basis.metadata_compatibility_note,
        basis.notes,
        support_note,
    ):
        lowered = text.lower()
        assert "supersym" not in lowered
        assert "susy " not in lowered
