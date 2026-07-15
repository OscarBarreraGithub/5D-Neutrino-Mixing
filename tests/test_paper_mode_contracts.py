"""Contract tests for the dedicated paper-facing 0710.1869 package."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

from quarkConstraints.paper_0710_1869.conventions import (
    PAPER_0710_1869_CONVENTIONS_SCHEMA_ID,
    PAPER_0710_1869_MODE_ID,
    PAPER_0710_1869_KK_GLUON_NORMALIZATION_ID,
    PAPER_0710_1869_MATCHING_ID,
    PAPER_0710_1869_OPERATOR_BASIS_ID,
    PAPER_0710_1869_PROFILE_DERIVATION_POLICY_ID,
    PAPER_0710_1869_PROVENANCE_POLICY_ID,
    PAPER_0710_1869_SEED_TO_PROFILE_MAPPING_POLICY_ID,
    PAPER_0710_1869_UNIVERSAL_TERM_COEFFICIENT_POLICY_ID,
    PAPER_0710_1869_VERIFIER_POLICY_ID,
    Paper07101869Conventions,
)
from quarkConstraints.paper_0710_1869.inputs import (
    PAPER_0710_1869_AFFINE_BULK_MASS_LEADING_TERM_COEFFICIENT,
    PAPER_0710_1869_AFFINE_BULK_MASS_SECTOR_POLICY_SCHEMA_ID,
    PAPER_0710_1869_AFFINE_BULK_MASS_UNIVERSAL_OFFSET,
    PAPER_0710_1869_BULK_MASS_PARAMETERIZATION_ID,
    PAPER_0710_1869_EIGENVALUE_ORDERING_ID,
    PAPER_0710_1869_PHYSICAL_SEED_TO_PROFILE_CONTRACT_SCHEMA_ID,
    PAPER_0710_1869_SEED_TO_PROFILE_MAPPING_SCHEMA_ID,
    PAPER_0710_1869_UNIVERSAL_TERM_POLICY_SCHEMA_ID,
    Paper07101869AffineBulkMassSectorPolicy,
    Paper07101869PhysicalSeedToProfileContract,
    Paper07101869SeedToProfileMappingPolicy,
    Paper07101869UniversalTermPolicy,
)
from quarkConstraints.paper_0710_1869.scales import (
    PAPER_0710_1869_SCALES_SCHEMA_ID,
    Paper07101869ScalePoint,
)
from quarkConstraints.paper_0710_1869.scan import (
    PAPER_0710_1869_SCAN_SCHEMA_ID,
    Paper07101869ScanRequest,
    build_structural_scan_rows,
)
from quarkConstraints.paper_0710_1869.validation import module_has_forbidden_import

REPO_ROOT = Path(__file__).resolve().parents[1]
PACKAGE_ROOT = REPO_ROOT / "quarkConstraints" / "paper_0710_1869"
CANONICAL_MODULES = (
    "__init__.py",
    "conventions.py",
    "couplings.py",
    "kkgluon.py",
    "scales.py",
    "scan.py",
    "validation.py",
)
OPTIONAL_CANONICAL_MODULES = tuple(
    str(path.relative_to(PACKAGE_ROOT))
    for path in (
        PACKAGE_ROOT / "eft_deltaf2" / "operators.py",
        PACKAGE_ROOT / "eft_deltaf2" / "matching_kkgluon.py",
        PACKAGE_ROOT / "eft_deltaf2" / "rg.py",
        PACKAGE_ROOT / "eft_deltaf2" / "hadronic.py",
        PACKAGE_ROOT / "eft_deltaf2" / "observables.py",
    )
    if path.exists()
)
OPTIONAL_RUNTIME_IMPORTS = tuple(
    (module_path, f"quarkConstraints.paper_0710_1869.{module_path[:-3].replace('/', '.')}")
    for module_path in (
        "eft_deltaf2/rg.py",
        "eft_deltaf2/hadronic.py",
        "eft_deltaf2/observables.py",
    )
    if (PACKAGE_ROOT / module_path).exists()
)
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


def test_conventions_validate_default_canonical_contract():
    conventions = Paper07101869Conventions()

    assert conventions.schema_id == PAPER_0710_1869_CONVENTIONS_SCHEMA_ID
    assert conventions.mode_id == PAPER_0710_1869_MODE_ID
    assert conventions.seed_to_profile_mapping_policy_id == (
        PAPER_0710_1869_SEED_TO_PROFILE_MAPPING_POLICY_ID
    )
    assert conventions.universal_term_coefficient_policy_id == (
        PAPER_0710_1869_UNIVERSAL_TERM_COEFFICIENT_POLICY_ID
    )
    assert conventions.profile_derivation_policy_id == (
        PAPER_0710_1869_PROFILE_DERIVATION_POLICY_ID
    )
    assert conventions.operator_basis_id == PAPER_0710_1869_OPERATOR_BASIS_ID
    assert conventions.matching_id == PAPER_0710_1869_MATCHING_ID
    assert conventions.kk_gluon_normalization_id == PAPER_0710_1869_KK_GLUON_NORMALIZATION_ID
    assert conventions.provenance_policy_id == PAPER_0710_1869_PROVENANCE_POLICY_ID
    assert conventions.verifier_policy_id == PAPER_0710_1869_VERIFIER_POLICY_ID
    assert conventions.rg_order == "lo"
    assert conventions.observable_scope == "np_only"


@pytest.mark.parametrize(
    ("overrides", "expected_message"),
    [
        ({"schema_id": "wrong.schema"}, "schema_id"),
        ({"mode_id": "proxy"}, "mode_id"),
        ({"rg_order": "nnlo"}, "rg_order"),
        ({"observable_scope": "full"}, "observable_scope"),
    ],
)
def test_conventions_reject_invalid_identifiers(overrides, expected_message):
    with pytest.raises(ValueError, match=expected_message):
        Paper07101869Conventions(**overrides)


def test_qs1_physical_seed_to_profile_policies_default_to_frozen_contract():
    sector_policy = Paper07101869AffineBulkMassSectorPolicy()
    universal_policy = Paper07101869UniversalTermPolicy()
    mapping_policy = Paper07101869SeedToProfileMappingPolicy()
    contract = Paper07101869PhysicalSeedToProfileContract()
    expected_affine_coefficients = [
        (
            PAPER_0710_1869_AFFINE_BULK_MASS_LEADING_TERM_COEFFICIENT,
            PAPER_0710_1869_AFFINE_BULK_MASS_UNIVERSAL_OFFSET,
        )
    ] * 3

    assert sector_policy.schema_id == PAPER_0710_1869_AFFINE_BULK_MASS_SECTOR_POLICY_SCHEMA_ID
    assert sector_policy.leading_term_coefficient == (
        PAPER_0710_1869_AFFINE_BULK_MASS_LEADING_TERM_COEFFICIENT
    )
    assert sector_policy.universal_offset == PAPER_0710_1869_AFFINE_BULK_MASS_UNIVERSAL_OFFSET

    assert universal_policy.schema_id == PAPER_0710_1869_UNIVERSAL_TERM_POLICY_SCHEMA_ID
    assert universal_policy.policy_id == PAPER_0710_1869_UNIVERSAL_TERM_COEFFICIENT_POLICY_ID
    assert [policy.sector_id for policy in universal_policy.sector_policies] == ["Q", "u", "d"]
    # C-6: the audited seed->profile policy is negative-slope, c = -lambda + 0.6.
    assert [
        (policy.leading_term_coefficient, policy.universal_offset)
        for policy in universal_policy.sector_policies
    ] == expected_affine_coefficients

    assert mapping_policy.schema_id == PAPER_0710_1869_SEED_TO_PROFILE_MAPPING_SCHEMA_ID
    assert mapping_policy.policy_id == PAPER_0710_1869_SEED_TO_PROFILE_MAPPING_POLICY_ID
    assert mapping_policy.bulk_mass_parameterization_id == PAPER_0710_1869_BULK_MASS_PARAMETERIZATION_ID
    assert mapping_policy.eigenvalue_ordering_id == PAPER_0710_1869_EIGENVALUE_ORDERING_ID
    assert mapping_policy.profile_derivation_policy_id == (
        PAPER_0710_1869_PROFILE_DERIVATION_POLICY_ID
    )
    assert mapping_policy.uses_hidden_bulk_mass_map_surrogate is False

    assert contract.schema_id == PAPER_0710_1869_PHYSICAL_SEED_TO_PROFILE_CONTRACT_SCHEMA_ID
    assert contract.mapping_policy.policy_id == PAPER_0710_1869_SEED_TO_PROFILE_MAPPING_POLICY_ID
    assert contract.universal_term_policy.policy_id == (
        PAPER_0710_1869_UNIVERSAL_TERM_COEFFICIENT_POLICY_ID
    )
    assert contract.mapping_policy.bulk_mass_parameterization_id == (
        PAPER_0710_1869_BULK_MASS_PARAMETERIZATION_ID
    )
    assert contract.mapping_policy.eigenvalue_ordering_id == PAPER_0710_1869_EIGENVALUE_ORDERING_ID
    assert [
        (policy.leading_term_coefficient, policy.universal_offset)
        for policy in contract.universal_term_policy.sector_policies
    ] == expected_affine_coefficients


@pytest.mark.parametrize(
    ("overrides", "expected_message"),
    [
        ({"leading_term_coefficient": 1.5}, "leading_term_coefficient"),
    ],
)
def test_affine_sector_policy_rejects_widened_numerics(overrides, expected_message):
    with pytest.raises(ValueError, match=expected_message):
        Paper07101869AffineBulkMassSectorPolicy(**overrides)


def test_affine_sector_policy_accepts_explicit_finite_offset_override():
    # C-6 freezes the default offset at 0.6, while non-default sector policies only
    # require a finite explicit offset and a negative leading coefficient.
    policy = Paper07101869AffineBulkMassSectorPolicy(universal_offset=0.25)

    assert policy.leading_term_coefficient == PAPER_0710_1869_AFFINE_BULK_MASS_LEADING_TERM_COEFFICIENT
    assert policy.universal_offset == 0.25


@pytest.mark.parametrize(
    ("overrides", "expected_message"),
    [
        ({"seed_to_profile_mapping_policy_id": "widened"}, "seed_to_profile_mapping_policy_id"),
        ({"universal_term_coefficient_policy_id": "widened"}, "universal_term_coefficient_policy_id"),
        ({"profile_derivation_policy_id": "widened"}, "profile_derivation_policy_id"),
    ],
)
def test_conventions_reject_widened_qs1_policy_ids(overrides, expected_message):
    with pytest.raises(ValueError, match=expected_message):
        Paper07101869Conventions(**overrides)


@pytest.mark.parametrize(
    ("overrides", "expected_message"),
    [
        ({"operator_basis_id": "widened"}, "operator_basis_id"),
        ({"matching_id": "widened"}, "matching_id"),
        ({"kk_gluon_normalization_id": "widened"}, "kk_gluon_normalization_id"),
        ({"provenance_policy_id": "widened"}, "provenance_policy_id"),
        ({"verifier_policy_id": "widened"}, "verifier_policy_id"),
    ],
)
def test_conventions_reject_widened_canonical_bundle_ids(overrides, expected_message):
    with pytest.raises(ValueError, match=expected_message):
        Paper07101869Conventions(**overrides)


@pytest.mark.parametrize(
    ("overrides", "expected_message"),
    [
        ({"policy_id": "widened"}, "policy_id"),
        ({"bulk_mass_parameterization_id": "widened"}, "bulk_mass_parameterization_id"),
        ({"eigenvalue_ordering_id": "widened"}, "eigenvalue_ordering_id"),
        ({"profile_derivation_policy_id": "widened"}, "profile_derivation_policy_id"),
    ],
)
def test_seed_to_profile_mapping_policy_rejects_widened_overrides(overrides, expected_message):
    with pytest.raises(ValueError, match=expected_message):
        Paper07101869SeedToProfileMappingPolicy(**overrides)


def test_universal_term_policy_rejects_widened_sector_coefficients():
    widened = Paper07101869AffineBulkMassSectorPolicy()
    object.__setattr__(widened, "leading_term_coefficient", 2.0)

    with pytest.raises(ValueError, match="sector_policies must use negative leading_term_coefficient"):
        Paper07101869UniversalTermPolicy(
            sector_policies=(
                widened,
                Paper07101869AffineBulkMassSectorPolicy(sector_id="u"),
                Paper07101869AffineBulkMassSectorPolicy(sector_id="d"),
            )
        )


def test_physical_seed_to_profile_contract_rejects_widened_nested_policy_ids():
    mapping_policy = Paper07101869SeedToProfileMappingPolicy()
    object.__setattr__(mapping_policy, "policy_id", "widened")

    with pytest.raises(ValueError, match="policy_id must be exactly"):
        Paper07101869PhysicalSeedToProfileContract(mapping_policy=mapping_policy)


def test_physical_seed_to_profile_contract_revalidates_nested_policy_contents():
    universal_policy = Paper07101869UniversalTermPolicy()
    object.__setattr__(
        universal_policy.sector_policies[0],
        "leading_term_coefficient",
        2.0,
    )

    with pytest.raises(ValueError, match="leading_term_coefficient"):
        Paper07101869PhysicalSeedToProfileContract(universal_term_policy=universal_policy)


def test_physical_seed_to_profile_contract_revalidates_mapping_policy_contents():
    mapping_policy = Paper07101869SeedToProfileMappingPolicy()
    object.__setattr__(mapping_policy, "eigenvalue_ordering_id", "widened")

    with pytest.raises(ValueError, match="eigenvalue_ordering_id"):
        Paper07101869PhysicalSeedToProfileContract(mapping_policy=mapping_policy)


def test_scale_point_keeps_matching_and_coupling_scales_explicit():
    point = Paper07101869ScalePoint(
        label="kaon_central",
        Lambda_IR_GeV=2500.0,
        m_g1_GeV=4100.0,
        xi_g=4100.0 / 2500.0,
        mu_match_GeV=3500.0,
        mu_gs_GeV=4100.0,
        m_KK_eff_GeV=3900.0,
    )

    assert point.schema_id == PAPER_0710_1869_SCALES_SCHEMA_ID
    assert point.xi_g == pytest.approx(4100.0 / 2500.0)
    assert point.mu_match_GeV == 3500.0
    assert point.mu_gs_GeV == 4100.0
    assert point.propagator_mass_GeV == 3900.0
    assert point.has_explicit_effective_kk_scale


@pytest.mark.parametrize(
    ("overrides", "expected_message"),
    [
        ({"label": ""}, "label"),
        ({"Lambda_IR_GeV": 0.0}, "Lambda_IR_GeV"),
        ({"m_g1_GeV": -1.0}, "m_g1_GeV"),
        ({"xi_g": 0.0}, "xi_g"),
        ({"Lambda_IR_GeV": 2500.0, "m_g1_GeV": 4100.0, "xi_g": 1.0}, "xi_g"),
        ({"mu_match_GeV": float("nan")}, "mu_match_GeV"),
        ({"m_KK_eff_GeV": -10.0}, "m_KK_eff_GeV"),
    ],
)
def test_scale_point_rejects_invalid_values(overrides, expected_message):
    with pytest.raises(ValueError, match=expected_message):
        Paper07101869ScalePoint(**overrides)


def test_scan_request_emits_structural_rows_without_overloaded_mass_contract():
    request = Paper07101869ScanRequest(
        scale_points=(
            Paper07101869ScalePoint(
                label="kaon_lo",
                Lambda_IR_GeV=3000.0,
                m_g1_GeV=4300.0,
                mu_match_GeV=3000.0,
                mu_gs_GeV=4300.0,
            ),
            Paper07101869ScalePoint(
                label="kaon_lo_eff",
                Lambda_IR_GeV=3000.0,
                m_g1_GeV=4300.0,
                mu_match_GeV=3000.0,
                mu_gs_GeV=4300.0,
                m_KK_eff_GeV=4050.0,
            ),
        )
    )

    rows = build_structural_scan_rows(request)

    assert request.schema_id == PAPER_0710_1869_SCAN_SCHEMA_ID
    assert [row.point_id for row in rows] == ["kaon_lo", "kaon_lo_eff"]
    assert rows[0].propagator_mass_GeV == 4300.0
    assert rows[1].propagator_mass_GeV == 4050.0
    assert all(row.status == "structural_only" for row in rows)


def test_scan_request_rejects_duplicate_labels():
    first = Paper07101869ScalePoint(label="duplicate")
    second = Paper07101869ScalePoint(label="duplicate")

    with pytest.raises(ValueError, match="labels must be unique"):
        Paper07101869ScanRequest(scale_points=(first, second))


@pytest.mark.parametrize("module_name", CANONICAL_MODULES + OPTIONAL_CANONICAL_MODULES)
def test_canonical_paper_modules_do_not_import_proxy_deltaf2(module_name):
    module_path = PACKAGE_ROOT / module_name

    assert module_path.exists()
    assert not module_has_forbidden_import(
        module_path,
        FORBIDDEN_REPO_V1_MODULES,
    )


def test_importing_paper_package_does_not_load_repo_v1_proxy_modules():
    script = """
import importlib
import json
import sys

importlib.import_module("quarkConstraints.paper_0710_1869")

forbidden_roots = {
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

forbidden = sorted(
        name
        for name in sys.modules
        if (
            name in forbidden_roots
            or any(name.startswith(f"{root}.") for root in forbidden_roots)
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


@pytest.mark.parametrize(("module_path", "module_import"), OPTIONAL_RUNTIME_IMPORTS)
def test_importing_optional_eft_module_does_not_load_repo_v1_proxy_modules(
    module_path: str,
    module_import: str,
):
    if not (PACKAGE_ROOT / module_path).exists():
        pytest.skip(f"{module_path} is not implemented in this checkout")

    script = f"""
import importlib
import json
import sys

importlib.import_module({module_import!r})

forbidden_roots = {
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

forbidden = sorted(
    name
    for name in sys.modules
    if (
        name in forbidden_roots
        or any(name.startswith(f"{{root}}.") for root in forbidden_roots)
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
