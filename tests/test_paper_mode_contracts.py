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
    Paper07101869Conventions,
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
    assert conventions.rg_order == "lo"
    assert conventions.observable_scope == "np_only"
    assert conventions.verifier_policy_id == "independent_verifier.required.v1"


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
