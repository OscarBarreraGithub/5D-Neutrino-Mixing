"""Contract tests for the modern quark input registry fence."""

from __future__ import annotations

from dataclasses import replace
import json
import subprocess
import sys
from pathlib import Path

import pytest

from quarkConstraints.modern.conventions import (
    MODERN_BUNDLE_ID_TEMPLATE,
    MODERN_CONVENTIONS_SCHEMA_ID,
    MODERN_DEFAULT_BUNDLE_FAMILY_ID,
    MODERN_DEFAULT_BUNDLE_SLUG,
    MODERN_DEFAULT_INPUT_BUNDLE_ID,
    MODERN_DEFAULT_INPUT_PROVENANCE_ID,
    MODERN_INPUT_NAMESPACE,
    MODERN_INPUT_REGISTRY_SCHEMA_ID,
    MODERN_LANE_ID,
    MODERN_PROVENANCE_ID_TEMPLATE,
    MODERN_PROVENANCE_NAMESPACE,
    MODERN_STRICT_PAPER_BUNDLE_FAMILY_ID,
    MODERN_STRICT_PAPER_BUNDLE_SLUG,
    MODERN_STRICT_PAPER_INPUT_BUNDLE_ID,
    MODERN_STRICT_PAPER_INPUT_PROVENANCE_ID,
    MODERN_STRICT_PAPER_RESOLUTION_POLICY_ID,
    MODERN_STRICT_PAPER_SOURCE_SNAPSHOT_ID,
    MODERN_STRICT_PAPER_SOURCE_SNAPSHOT_NAMES,
    MODERN_STRICT_PAPER_SOURCE_SCHEMA_IDS,
    ModernLaneConventions,
    build_modern_bundle_id,
    build_modern_provenance_id,
)
from quarkConstraints.modern.inputs import (
    MODERN_DEFAULT_ALPHA_S_POLICY_ID,
    MODERN_DEFAULT_INPUTS_SCHEMA_ID,
    MODERN_DEFAULT_NEUTRAL_MESON_SYSTEM_IDS,
    MODERN_DEFAULT_PROVENANCE_RECORD_IDS,
    MODERN_DEFAULT_REFERENCE_ALPHA_S_3TEV,
    MODERN_DEFAULT_RESOLUTION_POLICY_ID,
    MODERN_DEFAULT_TARGET_SCALE_GEV,
    MODERN_STRICT_PAPER_INPUTS_SCHEMA_ID,
    ModernDefaultInputs,
    ModernInputRegistry,
    ModernStrictPaperInputs,
)
from quarkConstraints.paper_0710_1869.conventions import PAPER_0710_1869_MODE_ID
from quarkConstraints.paper_0710_1869.validation import module_has_forbidden_import

REPO_ROOT = Path(__file__).resolve().parents[1]
PACKAGE_ROOT = REPO_ROOT / "quarkConstraints" / "modern"

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

FORBIDDEN_PAPER_CLAIM_MODULES = {
    "quarkConstraints.paper_0710_1869.artifacts",
    "quarkConstraints.paper_0710_1869.benchmarks",
    "quarkConstraints.paper_0710_1869.couplings",
    "quarkConstraints.paper_0710_1869.deltaf2",
    "quarkConstraints.paper_0710_1869.eft_deltaf2",
    "quarkConstraints.paper_0710_1869.fit",
    "quarkConstraints.paper_0710_1869.model",
    "quarkConstraints.paper_0710_1869.proxies",
    "quarkConstraints.paper_0710_1869.scan",
    "quarkConstraints.paper_0710_1869.scales",
    "quarkConstraints.paper_0710_1869.verifier",
}


def test_modern_registry_freezes_lane_and_schema_ids():
    conventions = ModernLaneConventions()
    registry = ModernInputRegistry(conventions=conventions)

    strict = registry.strict_paper_inputs
    modern_default = registry.modern_default_inputs

    assert conventions.schema_id == MODERN_CONVENTIONS_SCHEMA_ID
    assert conventions.lane_id == MODERN_LANE_ID
    assert conventions.input_registry_schema_id == MODERN_INPUT_REGISTRY_SCHEMA_ID
    assert conventions.input_namespace == MODERN_INPUT_NAMESPACE
    assert conventions.provenance_namespace == MODERN_PROVENANCE_NAMESPACE
    assert conventions.bundle_id_template == MODERN_BUNDLE_ID_TEMPLATE
    assert conventions.provenance_id_template == MODERN_PROVENANCE_ID_TEMPLATE
    assert conventions.strict_paper_bundle_family_id == MODERN_STRICT_PAPER_BUNDLE_FAMILY_ID
    assert conventions.modern_default_bundle_family_id == MODERN_DEFAULT_BUNDLE_FAMILY_ID

    assert registry.schema_id == MODERN_INPUT_REGISTRY_SCHEMA_ID
    assert registry.lane_id == MODERN_LANE_ID

    assert strict.schema_id == MODERN_STRICT_PAPER_INPUTS_SCHEMA_ID
    assert strict.bundle_id == MODERN_STRICT_PAPER_INPUT_BUNDLE_ID
    assert strict.provenance_id == MODERN_STRICT_PAPER_INPUT_PROVENANCE_ID
    assert strict.source_snapshot_id == MODERN_STRICT_PAPER_SOURCE_SNAPSHOT_ID
    assert strict.source_snapshot_names == MODERN_STRICT_PAPER_SOURCE_SNAPSHOT_NAMES
    assert strict.source_schema_ids == MODERN_STRICT_PAPER_SOURCE_SCHEMA_IDS
    assert strict.source_resolution_policy_id == MODERN_STRICT_PAPER_RESOLUTION_POLICY_ID
    assert modern_default.schema_id == MODERN_DEFAULT_INPUTS_SCHEMA_ID
    assert modern_default.bundle_id == MODERN_DEFAULT_INPUT_BUNDLE_ID
    assert modern_default.provenance_id == MODERN_DEFAULT_INPUT_PROVENANCE_ID
    assert modern_default.source_resolution_policy_id == MODERN_DEFAULT_RESOLUTION_POLICY_ID
    assert tuple(item.system_id for item in modern_default.neutral_meson_inputs) == (
        MODERN_DEFAULT_NEUTRAL_MESON_SYSTEM_IDS
    )
    assert modern_default.operator_weight_policy.reference_scale_GeV == pytest.approx(
        MODERN_DEFAULT_TARGET_SCALE_GEV
    )
    assert modern_default.operator_weight_policy.lr1_weight == pytest.approx(7.0)
    assert modern_default.ckm_target.target_id
    assert modern_default.ckm_target.theta12 == pytest.approx(0.2274)
    assert modern_default.ckm_target.theta13 == pytest.approx(0.00368)
    assert modern_default.ckm_target.theta23 == pytest.approx(0.0415)
    assert modern_default.ckm_target.delta == pytest.approx(1.196)
    assert modern_default.quark_mass_target.scale_GeV == pytest.approx(
        MODERN_DEFAULT_TARGET_SCALE_GEV
    )
    assert modern_default.quark_mass_target.up_masses_GeV == pytest.approx((0.0013, 0.62, 172.0))
    assert modern_default.quark_mass_target.down_masses_GeV == pytest.approx((0.0028, 0.057, 2.86))
    assert modern_default.qcd_metadata.matching_scale_GeV == pytest.approx(
        MODERN_DEFAULT_TARGET_SCALE_GEV
    )
    assert modern_default.qcd_metadata.xi_KK == pytest.approx(1.0)
    assert modern_default.qcd_metadata.alpha_s_reference_value == pytest.approx(
        MODERN_DEFAULT_REFERENCE_ALPHA_S_3TEV
    )
    assert modern_default.qcd_metadata.alpha_s_policy_id == MODERN_DEFAULT_ALPHA_S_POLICY_ID
    assert tuple(record.record_id for record in modern_default.provenance_records) == (
        MODERN_DEFAULT_PROVENANCE_RECORD_IDS
    )
    assert modern_default.paper_inputs == ()

    assert build_modern_bundle_id(MODERN_STRICT_PAPER_BUNDLE_FAMILY_ID, MODERN_STRICT_PAPER_BUNDLE_SLUG) == (
        MODERN_STRICT_PAPER_INPUT_BUNDLE_ID
    )
    assert build_modern_provenance_id(
        MODERN_STRICT_PAPER_BUNDLE_FAMILY_ID,
        MODERN_STRICT_PAPER_BUNDLE_SLUG,
    ) == MODERN_STRICT_PAPER_INPUT_PROVENANCE_ID
    assert build_modern_bundle_id(MODERN_DEFAULT_BUNDLE_FAMILY_ID, MODERN_DEFAULT_BUNDLE_SLUG) == (
        MODERN_DEFAULT_INPUT_BUNDLE_ID
    )
    assert build_modern_provenance_id(
        MODERN_DEFAULT_BUNDLE_FAMILY_ID,
        MODERN_DEFAULT_BUNDLE_SLUG,
    ) == MODERN_DEFAULT_INPUT_PROVENANCE_ID


def test_modern_registry_separates_strict_paper_and_modern_default_bundle_families():
    registry = ModernInputRegistry()
    strict = registry.strict_paper_inputs
    modern_default = registry.modern_default_inputs

    assert strict.family_id == MODERN_STRICT_PAPER_BUNDLE_FAMILY_ID
    assert modern_default.family_id == MODERN_DEFAULT_BUNDLE_FAMILY_ID
    assert strict.bundle_slug == MODERN_STRICT_PAPER_BUNDLE_SLUG
    assert modern_default.bundle_slug == MODERN_DEFAULT_BUNDLE_SLUG
    assert strict.bundle_id != modern_default.bundle_id
    assert strict.provenance_id != modern_default.provenance_id
    assert strict.source_lane_id == PAPER_0710_1869_MODE_ID
    assert modern_default.source_lane_id == MODERN_LANE_ID
    resolved = strict.resolve_paper_inputs()
    assert len(resolved) == 2
    assert strict.validate_resolved_paper_inputs(resolved) == resolved
    assert modern_default.neutral_meson_inputs[0].system_id == "epsilon_K"
    assert modern_default.neutral_meson_inputs[-1].system_id == "D0"
    assert modern_default.operator_weight_policy.policy_id == (
        modern_default.neutral_meson_inputs[0].weight_policy_id
    )
    assert modern_default.paper_inputs == ()


@pytest.mark.parametrize(
    ("factory", "kwargs", "expected_message"),
    [
        (ModernLaneConventions, {"lane_id": "paper_0710_1869"}, "lane_id"),
        (ModernLaneConventions, {"input_namespace": "widened"}, "input_namespace"),
        (ModernLaneConventions, {"provenance_namespace": "widened"}, "provenance_namespace"),
        (ModernStrictPaperInputs, {"family_id": MODERN_DEFAULT_BUNDLE_FAMILY_ID}, "family_id"),
        (ModernStrictPaperInputs, {"bundle_id": MODERN_DEFAULT_INPUT_BUNDLE_ID}, "bundle_id"),
        (
            ModernStrictPaperInputs,
            {"provenance_id": MODERN_DEFAULT_INPUT_PROVENANCE_ID},
            "provenance_id",
        ),
        (
            ModernStrictPaperInputs,
            {"source_snapshot_id": "widened"},
            "source_snapshot_id",
        ),
        (
            ModernStrictPaperInputs,
            {"source_snapshot_names": ("other", "loader")},
            "source_snapshot_names",
        ),
        (
            ModernStrictPaperInputs,
            {"source_schema_ids": ("other.schema", "loader.schema")},
            "source_schema_ids",
        ),
        (
            ModernStrictPaperInputs,
            {"source_resolution_policy_id": "widened"},
            "source_resolution_policy_id",
        ),
        (ModernDefaultInputs, {"family_id": MODERN_STRICT_PAPER_BUNDLE_FAMILY_ID}, "family_id"),
        (ModernDefaultInputs, {"bundle_id": MODERN_STRICT_PAPER_INPUT_BUNDLE_ID}, "bundle_id"),
        (
            ModernDefaultInputs,
            {"provenance_id": MODERN_STRICT_PAPER_INPUT_PROVENANCE_ID},
            "provenance_id",
        ),
        (
            ModernDefaultInputs,
            {"source_resolution_policy_id": "widened"},
            "source_resolution_policy_id",
        ),
        (ModernDefaultInputs, {"neutral_meson_inputs": []}, "neutral_meson_inputs"),
        (ModernDefaultInputs, {"provenance_records": []}, "provenance_records"),
        (ModernDefaultInputs, {"paper_inputs": []}, "paper_inputs"),
        (ModernInputRegistry, {"schema_id": "widened"}, "schema_id"),
    ],
)
def test_modern_registry_rejects_bundle_aliasing_or_widened_ids(
    factory, kwargs, expected_message
):
    with pytest.raises(ValueError, match=expected_message):
        factory(**kwargs)


def test_modern_registry_modules_do_not_import_repo_v1_or_forbidden_paper_claim_modules():
    assert not module_has_forbidden_import(
        PACKAGE_ROOT / "__init__.py",
        FORBIDDEN_REPO_V1_MODULES | FORBIDDEN_PAPER_CLAIM_MODULES,
    )
    assert not module_has_forbidden_import(
        PACKAGE_ROOT / "conventions.py",
        FORBIDDEN_REPO_V1_MODULES | FORBIDDEN_PAPER_CLAIM_MODULES,
    )
    assert not module_has_forbidden_import(
        PACKAGE_ROOT / "inputs.py",
        FORBIDDEN_REPO_V1_MODULES | FORBIDDEN_PAPER_CLAIM_MODULES,
    )


def test_importing_modern_package_does_not_load_repo_v1_or_paper_modules():
    script = """
import importlib
import json
import sys

importlib.import_module("quarkConstraints.modern")

forbidden = sorted(
    name
    for name in sys.modules
    if (
        name in {
            "deltaf2",
            "quarkConstraints.benchmarks",
            "quarkConstraints.couplings",
            "quarkConstraints.deltaf2",
            "quarkConstraints.fit",
            "quarkConstraints.model",
            "quarkConstraints.paper_0710_1869",
            "quarkConstraints.paper_0710_1869.artifacts",
            "quarkConstraints.paper_0710_1869.benchmarks",
            "quarkConstraints.paper_0710_1869.couplings",
            "quarkConstraints.paper_0710_1869.deltaf2",
            "quarkConstraints.paper_0710_1869.fit",
            "quarkConstraints.paper_0710_1869.model",
            "quarkConstraints.paper_0710_1869.proxies",
            "quarkConstraints.paper_0710_1869.scan",
            "quarkConstraints.paper_0710_1869.scales",
            "quarkConstraints.paper_0710_1869.verifier",
            "quarkConstraints.proxies",
            "quarkConstraints.scan",
            "quarkConstraints.scales",
            "quarkConstraints.validation",
        }
        or name.startswith("quarkConstraints.benchmarks.")
        or name.startswith("quarkConstraints.couplings.")
        or name.startswith("quarkConstraints.deltaf2.")
        or name.startswith("quarkConstraints.fit.")
        or name.startswith("quarkConstraints.model.")
        or name.startswith("quarkConstraints.paper_0710_1869.artifacts.")
        or name.startswith("quarkConstraints.paper_0710_1869.benchmarks.")
        or name.startswith("quarkConstraints.paper_0710_1869.couplings.")
        or name.startswith("quarkConstraints.paper_0710_1869.deltaf2.")
        or name.startswith("quarkConstraints.paper_0710_1869.fit.")
        or name.startswith("quarkConstraints.paper_0710_1869.model.")
        or name.startswith("quarkConstraints.paper_0710_1869.proxies.")
        or name.startswith("quarkConstraints.paper_0710_1869.scan.")
        or name.startswith("quarkConstraints.paper_0710_1869.scales.")
        or name.startswith("quarkConstraints.paper_0710_1869.verifier.")
        or name.startswith("quarkConstraints.proxies.")
        or name.startswith("quarkConstraints.scan.")
        or name.startswith("quarkConstraints.scales.")
        or name.startswith("quarkConstraints.validation.")
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


def test_importing_modern_inputs_stays_within_thin_paper_adapter_boundary():
    script = """
import importlib
import json
import sys

importlib.import_module("quarkConstraints.modern.inputs")

loaded = sorted(name for name in sys.modules if name.startswith("quarkConstraints.paper_0710_1869"))
print(json.dumps(loaded))
"""
    completed = subprocess.run(
        [sys.executable, "-c", script],
        check=True,
        capture_output=True,
        text=True,
        cwd=str(REPO_ROOT),
    )
    loaded = json.loads(completed.stdout)

    assert loaded == []


def test_modern_registry_construction_does_not_load_paper_modules():
    script = """
import json
import sys

from quarkConstraints.modern.inputs import ModernInputRegistry

ModernInputRegistry()

loaded = sorted(
    name for name in sys.modules if name.startswith("quarkConstraints.paper_0710_1869")
)
print(json.dumps(loaded))
"""
    completed = subprocess.run(
        [sys.executable, "-c", script],
        check=True,
        capture_output=True,
        text=True,
        cwd=str(REPO_ROOT),
    )
    loaded = json.loads(completed.stdout)

    assert loaded == []


def test_modern_strict_paper_inputs_resolve_path_does_not_load_paper_modules():
    script = """
import json
import sys

from quarkConstraints.modern.inputs import ModernStrictPaperInputs

strict = ModernStrictPaperInputs()
resolved = strict.resolve_paper_inputs()
strict.validate_resolved_paper_inputs(resolved)

loaded = sorted(
    name for name in sys.modules if name.startswith("quarkConstraints.paper_0710_1869")
)
print(json.dumps(loaded))
"""
    completed = subprocess.run(
        [sys.executable, "-c", script],
        check=True,
        capture_output=True,
        text=True,
        cwd=str(REPO_ROOT),
    )
    loaded = json.loads(completed.stdout)

    assert loaded == []


def test_modern_strict_paper_inputs_rejects_mutated_canonical_payload_aliasing():
    strict = ModernStrictPaperInputs()
    canonical = strict.resolve_paper_inputs()
    table, eq3 = canonical
    object.__setattr__(table, "schema_id", "widened.schema")

    with pytest.raises(ValueError, match="Table I schema_id"):
        strict.validate_resolved_paper_inputs(canonical)

    object.__setattr__(table, "schema_id", MODERN_STRICT_PAPER_SOURCE_SCHEMA_IDS[0])
    object.__setattr__(
        eq3,
        "payload_signature",
        (("theta12_deg", 999.0),),
    )

    with pytest.raises(ValueError, match="Eq\\. \\(3\\) snapshot"):
        strict.validate_resolved_paper_inputs(canonical)


def test_modern_default_inputs_export_explicit_versioned_numeric_payloads_deterministically():
    default_inputs = ModernDefaultInputs()

    payload_a = default_inputs.as_dict()
    payload_b = ModernDefaultInputs().as_dict()

    assert payload_a == payload_b
    assert json.dumps(payload_a, sort_keys=True) == json.dumps(payload_b, sort_keys=True)
    assert payload_a["schema_id"] == MODERN_DEFAULT_INPUTS_SCHEMA_ID
    assert payload_a["source_resolution_policy_id"] == MODERN_DEFAULT_RESOLUTION_POLICY_ID
    assert [item["system_id"] for item in payload_a["neutral_meson_inputs"]] == list(
        MODERN_DEFAULT_NEUTRAL_MESON_SYSTEM_IDS
    )
    assert payload_a["neutral_meson_inputs"][0]["schema_id"].endswith(".v1")
    assert payload_a["operator_weight_policy"]["schema_id"].endswith(".v1")
    assert payload_a["operator_weight_policy"]["reference_scale_GeV"] == pytest.approx(
        MODERN_DEFAULT_TARGET_SCALE_GEV
    )
    assert payload_a["ckm_target"]["schema_id"].endswith(".v1")
    assert payload_a["ckm_target"]["theta12"] == pytest.approx(0.2274)
    assert payload_a["quark_mass_target"]["schema_id"].endswith(".v1")
    assert payload_a["quark_mass_target"]["up_masses_GeV"] == pytest.approx(
        [0.0013, 0.62, 172.0]
    )
    assert payload_a["qcd_metadata"]["schema_id"].endswith(".v1")
    assert payload_a["qcd_metadata"]["alpha_s_policy_id"] == MODERN_DEFAULT_ALPHA_S_POLICY_ID
    assert payload_a["qcd_metadata"]["alpha_s_reference_value"] == pytest.approx(
        MODERN_DEFAULT_REFERENCE_ALPHA_S_3TEV
    )
    assert [record["record_id"] for record in payload_a["provenance_records"]] == list(
        MODERN_DEFAULT_PROVENANCE_RECORD_IDS
    )
    assert payload_a["paper_inputs"] == []


def test_modern_default_inputs_reject_mutated_frozen_numeric_payloads():
    default_inputs = ModernDefaultInputs()
    mutated_system = replace(
        default_inputs.neutral_meson_inputs[0],
        bound=default_inputs.neutral_meson_inputs[0].bound * 2.0,
    )
    with pytest.raises(ValueError, match="neutral_meson_inputs must match the frozen default payload"):
        ModernDefaultInputs(
            neutral_meson_inputs=(mutated_system,) + default_inputs.neutral_meson_inputs[1:]
        )

    mutated_qcd = replace(
        default_inputs.qcd_metadata,
        alpha_s_reference_value=default_inputs.qcd_metadata.alpha_s_reference_value * 1.01,
    )
    with pytest.raises(ValueError, match="qcd_metadata must match the frozen default payload"):
        ModernDefaultInputs(qcd_metadata=mutated_qcd)


def test_modern_default_inputs_reject_non_tuple_paper_inputs():
    with pytest.raises(ValueError, match="paper_inputs must be a tuple"):
        ModernDefaultInputs(paper_inputs=[])  # type: ignore[arg-type]

    default_inputs = ModernDefaultInputs()
    assert isinstance(default_inputs.paper_inputs, tuple)
    assert default_inputs.paper_inputs == ()
    assert tuple(item.system_id for item in default_inputs.neutral_meson_inputs) == (
        MODERN_DEFAULT_NEUTRAL_MESON_SYSTEM_IDS
    )
