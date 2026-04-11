"""Forward-compatible tests for the paper-owned hadronic input layer."""

from __future__ import annotations

import dataclasses
import importlib
import json
import math
import subprocess
import sys
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
PACKAGE_ROOT = REPO_ROOT / "quarkConstraints" / "paper_0710_1869"
HADRONIC_MODULE_PATH = PACKAGE_ROOT / "eft_deltaf2" / "hadronic.py"

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
EXPECTED_SUPPORTED_OPERATORS = {"Q1_VLL", "Q1_VRR"}
EXPECTED_UNSUPPORTED_OPERATORS = {"Q4_LR", "Q5_LR"}
FORBIDDEN_DEFAULT_LR_FIELD_NAMES = (
    "q4_matrix_element_GeV4",
    "q5_matrix_element_GeV4",
    "B4_mu_had",
    "B5_mu_had",
    "R_chi_mu_had",
)


def _has_hadronic_module() -> bool:
    return HADRONIC_MODULE_PATH.exists()


def _load_hadronic_module():
    if not _has_hadronic_module():
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


def _default_hadronic_payload(module: Any) -> dict[str, Any]:
    callable_obj = _get_callable(
        module,
        DEFAULT_HADRONIC_SUMMARY_NAMES + DEFAULT_HADRONIC_EXPORT_NAMES,
    )
    if not callable(callable_obj):
        raise AssertionError(
            "hadronic layer exists but exposes no default kaon export; expected one of "
            + ", ".join(DEFAULT_HADRONIC_SUMMARY_NAMES + DEFAULT_HADRONIC_EXPORT_NAMES)
        )
    return _payload_from_value(callable_obj())


def _hadronic_builder(module: Any) -> Any:
    callable_obj = _get_callable(module, HADRONIC_BUILDER_NAMES)
    if not callable(callable_obj):
        raise AssertionError(
            "hadronic layer exists but exposes no builder callable; expected one of "
            + ", ".join(HADRONIC_BUILDER_NAMES)
        )
    return callable_obj


def _default_hadronic_payload_cross_process() -> dict[str, Any]:
    script = f"""
import dataclasses
import importlib
import json
from collections.abc import Mapping, Sequence
from pathlib import Path

NAMES = {DEFAULT_HADRONIC_SUMMARY_NAMES + DEFAULT_HADRONIC_EXPORT_NAMES!r}

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
        raise AssertionError("hadronic payload must canonicalize to a mapping")
    return payload

module = importlib.import_module("quarkConstraints.paper_0710_1869.eft_deltaf2.hadronic")
for name in NAMES:
    candidate = getattr(module, name, None)
    if callable(candidate):
        print(json.dumps(payload_from_value(candidate()), sort_keys=True))
        break
else:
    raise SystemExit("missing default hadronic export")
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


def _mapping_contains_key_deep(value: Any, target_key: str) -> bool:
    if isinstance(value, Mapping):
        if target_key in value:
            return True
        return any(_mapping_contains_key_deep(item, target_key) for item in value.values())
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return any(_mapping_contains_key_deep(item, target_key) for item in value)
    return False


def _positive_finite(value: Any) -> bool:
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return False
    return math.isfinite(numeric) and numeric > 0.0


def _is_kaon_system(value: Any) -> bool:
    if value is None:
        return False
    lowered = str(value).strip().lower()
    return lowered in {"k", "k0", "kaon"} or "kaon" in lowered


def test_importing_hadronic_module_does_not_load_repo_v1_modules() -> None:
    if not _has_hadronic_module():
        pytest.skip("paper_0710_1869 hadronic layer not implemented yet")

    script = f"""
import importlib
import json
import sys

importlib.import_module("quarkConstraints.paper_0710_1869.eft_deltaf2.hadronic")

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


def test_default_hadronic_export_is_deterministic_and_contract_tagged() -> None:
    module = _load_hadronic_module()

    payload = _default_hadronic_payload(module)
    second = _default_hadronic_payload(module)
    cross_process_payload = _default_hadronic_payload_cross_process()

    assert json.dumps(payload, sort_keys=True) == json.dumps(second, sort_keys=True)
    assert json.dumps(payload, sort_keys=True) == json.dumps(
        cross_process_payload,
        sort_keys=True,
    )

    schema_tag = _nested_value(
        payload,
        (
            ("schema_id",),
            ("metadata", "schema_name"),
            ("contract", "schema_id"),
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
    normalization_id = _nested_value(
        payload,
        (
            ("operator_normalization_id",),
            ("normalization_id",),
            ("tags", "operator_normalization_id"),
            ("contract", "operator_normalization_id"),
        ),
    )
    basis_id = _nested_value(
        payload,
        (
            ("operator_basis_id",),
            ("basis_id",),
            ("tags", "operator_basis_id"),
            ("contract", "operator_basis_id"),
        ),
    )
    contract_system_id = _nested_value(
        payload,
        (
            ("contract", "system_id"),
        ),
    )
    contract_scheme_id = _nested_value(
        payload,
        (
            ("contract", "renormalization_scheme_id"),
            ("contract", "scheme_id"),
        ),
    )
    contract_normalization_id = _nested_value(
        payload,
        (
            ("contract", "operator_normalization_id"),
        ),
    )
    contract_basis_id = _nested_value(
        payload,
        (
            ("contract", "operator_basis_id"),
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
    contract_mu_had_GeV = _nested_value(
        payload,
        (
            ("contract", "mu_had_GeV"),
            ("contract", "evaluation_scale_GeV"),
        ),
    )
    mass_source_id = _nested_value(
        payload,
        (
            ("mass_source", "source_id"),
        ),
    )
    decay_constant_source_id = _nested_value(
        payload,
        (
            ("decay_constant_source", "source_id"),
        ),
    )
    bag_parameter_source_id = _nested_value(
        payload,
        (
            ("bag_parameter_source", "source_id"),
        ),
    )
    bundle_source_id = _nested_value(
        payload,
        (
            ("source_id",),
            ("tags", "source_id"),
        ),
    )
    provenance_ids = _nested_value(
        payload,
        (
            ("provenance_ids",),
            ("tags", "provenance_ids"),
        ),
    )
    parity_relation_id = _nested_value(
        payload,
        (
            ("parity_relation_id",),
            ("tags", "parity_relation_id"),
            ("contract", "parity_relation_id"),
        ),
    )
    contract_parity_relation_id = _nested_value(
        payload,
        (
            ("contract", "parity_relation_id"),
        ),
    )
    operator_names = _nested_value(
        payload,
        (
            ("supported_operator_names",),
            ("contract", "supported_operator_names"),
        ),
    )
    unsupported_operator_names = _nested_value(
        payload,
        (
            ("unsupported_operator_names",),
            ("contract", "unsupported_operator_names"),
        ),
    )

    checks = {
        "schema_tag_present": bool(schema_tag),
        "schema_tag_mentions_hadronic": schema_tag is not None
        and "hadronic" in str(schema_tag).lower(),
        "system_is_kaon": _is_kaon_system(system_id),
        "scheme_id_present": bool(scheme_id),
        "operator_basis_id_present": bool(basis_id),
        "operator_normalization_id_present": bool(normalization_id),
        "contract_system_matches_top_level": contract_system_id == system_id,
        "contract_basis_matches_top_level": contract_basis_id == basis_id,
        "contract_normalization_matches_top_level": (
            contract_normalization_id == normalization_id
        ),
        "contract_scheme_matches_top_level": contract_scheme_id == scheme_id,
        "contract_mu_had_matches_top_level": (
            _positive_finite(contract_mu_had_GeV)
            and float(contract_mu_had_GeV) == float(mu_had_GeV)
        ),
        "mu_had_present_and_positive": _positive_finite(mu_had_GeV),
        "source_ids_present": all(
            bool(source_id)
            for source_id in (
                mass_source_id,
                decay_constant_source_id,
                bag_parameter_source_id,
            )
        ),
        "bundle_source_present": bool(bundle_source_id),
        "bundle_source_in_provenance_ids": bool(bundle_source_id)
        and isinstance(provenance_ids, Sequence)
        and not isinstance(provenance_ids, (str, bytes, bytearray))
        and bundle_source_id in {str(item) for item in provenance_ids},
        "parity_relation_tag_present": bool(parity_relation_id),
        "parity_relation_tag_mentions_supported_relation": bool(parity_relation_id)
        and "q1_vll_equals_q1_vrr" in str(parity_relation_id).lower(),
        "contract_parity_relation_matches_top_level": (
            contract_parity_relation_id == parity_relation_id
        ),
    }
    assert not [name for name, ok in checks.items() if not ok], checks

    if isinstance(operator_names, Sequence) and not isinstance(
        operator_names,
        (str, bytes, bytearray),
    ):
        assert {str(item) for item in operator_names} == EXPECTED_SUPPORTED_OPERATORS
    if isinstance(unsupported_operator_names, Sequence) and not isinstance(
        unsupported_operator_names,
        (str, bytes, bytearray),
    ):
        assert {str(item) for item in unsupported_operator_names} == EXPECTED_UNSUPPORTED_OPERATORS
    for field_name in FORBIDDEN_DEFAULT_LR_FIELD_NAMES:
        assert not _mapping_contains_key_deep(payload, field_name), field_name


def test_hadronic_builder_preserves_override_sources_and_bundle_provenance_contract() -> None:
    module = _load_hadronic_module()
    builder = _hadronic_builder(module)

    default_bundle = builder()
    override_mass_source = dataclasses.replace(
        default_bundle.mass_source,
        source_id="override.mass.source.v1",
    )
    override_decay_constant_source = dataclasses.replace(
        default_bundle.decay_constant_source,
        source_id="override.decay_constant.source.v1",
    )
    override_bag_parameter_source = dataclasses.replace(
        default_bundle.bag_parameter_source,
        source_id="override.bag_parameter.source.v1",
    )

    overridden_bundle = builder(
        mass_source=override_mass_source,
        decay_constant_source=override_decay_constant_source,
        bag_parameter_source=override_bag_parameter_source,
    )
    second_bundle = builder(
        mass_source=override_mass_source,
        decay_constant_source=override_decay_constant_source,
        bag_parameter_source=override_bag_parameter_source,
    )

    overridden_payload = _payload_from_value(overridden_bundle)
    second_payload = _payload_from_value(second_bundle)
    assert json.dumps(overridden_payload, sort_keys=True) == json.dumps(
        second_payload,
        sort_keys=True,
    )
    assert overridden_bundle.mass_source == override_mass_source
    assert overridden_bundle.decay_constant_source == override_decay_constant_source
    assert overridden_bundle.bag_parameter_source == override_bag_parameter_source
    assert overridden_payload["mass_source"]["source_id"] == override_mass_source.source_id
    assert (
        overridden_payload["decay_constant_source"]["source_id"]
        == override_decay_constant_source.source_id
    )
    assert (
        overridden_payload["bag_parameter_source"]["source_id"]
        == override_bag_parameter_source.source_id
    )
    provenance_ids = overridden_payload.get("provenance_ids")
    assert isinstance(provenance_ids, list)
    assert overridden_payload["source_id"] in provenance_ids


@pytest.mark.parametrize(
    ("field_name", "replacement", "error_pattern"),
    (
        ("operator_basis_id", "kk_gluon_tree_np_only.v1.inconsistent", "(?i)(contract|basis)"),
        (
            "operator_normalization_id",
            "paper_0710_1869.deltaf2.kk_gluon_tree_color_normalization.v1.inconsistent",
            "(?i)(contract|normalization)",
        ),
        (
            "renormalization_scheme_id",
            "bmu.hep-ph-0005183.ndr-ms.lo.v1.inconsistent",
            "(?i)(contract|scheme)",
        ),
        ("mu_had_GeV", 2.5, "(?i)(contract|mu_had|evaluation)"),
        (
            "parity_relation_id",
            "kaon.q1_vll_equals_q1_vrr.by_parity.v1.inconsistent",
            "(?i)(contract|parity|relation)",
        ),
    ),
)
def test_hadronic_bundle_rejects_nested_contract_mismatch(
    field_name: str,
    replacement: str | float,
    error_pattern: str,
) -> None:
    module = _load_hadronic_module()
    builder = _hadronic_builder(module)
    bundle = builder()

    mismatched_contract = dataclasses.replace(bundle.contract, **{field_name: replacement})

    with pytest.raises(ValueError, match=error_pattern):
        dataclasses.replace(bundle, contract=mismatched_contract)


@pytest.mark.parametrize(
    ("field_name", "replacement", "error_pattern"),
    (
        (
            "renormalization_scheme_id",
            "arbitrary.scheme.v1",
            "(?i)(bag_parameter_source|scheme|derived)",
        ),
        (
            "transformation_id",
            "arbitrary.transformation.v1",
            "(?i)(bag_parameter_source|transformation|derived)",
        ),
    ),
)
def test_hadronic_builder_rejects_arbitrary_custom_bag_source_metadata_for_derived_bk(
    field_name: str,
    replacement: str,
    error_pattern: str,
) -> None:
    module = _load_hadronic_module()
    builder = _hadronic_builder(module)
    default_bundle = builder()
    custom_bag_source = dataclasses.replace(
        default_bundle.bag_parameter_source,
        source_id="custom.bag.source.v1",
        **{field_name: replacement},
    )

    with pytest.raises(ValueError, match=error_pattern):
        builder(
            bag_parameter_source=custom_bag_source,
        )


def test_hadronic_builder_rejects_custom_source_id_missing_from_provenance_ids() -> None:
    module = _load_hadronic_module()
    builder = _hadronic_builder(module)
    default_bundle = builder()
    custom_bag_source = dataclasses.replace(
        default_bundle.bag_parameter_source,
        source_id="custom.bag.source.v1",
    )

    with pytest.raises(ValueError, match="(?i)(source_id|provenance)"):
        builder(
            bag_parameter_source=custom_bag_source,
            source_id="custom.hadronic.source.v1",
            provenance_ids=("other.record.v1",),
        )
