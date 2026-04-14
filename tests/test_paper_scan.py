"""Acceptance tests for the paper-facing 0710.1869 QS1 scan surface."""

from __future__ import annotations

import inspect
from dataclasses import asdict

import numpy as np
import pytest

from quarkConstraints.paper_0710_1869.conventions import (
    PAPER_0710_1869_MODE_ID,
    PAPER_0710_1869_CONVENTIONS_SCHEMA_ID,
    Paper07101869Conventions,
)
from quarkConstraints.paper_0710_1869.scales import Paper07101869ScalePoint
from quarkConstraints.paper_0710_1869.scan import (
    PAPER_0710_1869_SCAN_SCHEMA_ID,
    PAPER_0710_1869_STRICT_PAPER_BENCHMARK_ID,
    PAPER_0710_1869_STRICT_PAPER_CLAIM_LEVEL_ID,
    PAPER_0710_1869_STRICT_PAPER_INPUT_BUNDLE_ID,
    PAPER_0710_1869_STRICT_PAPER_INPUT_PROVENANCE_ID,
    PAPER_0710_1869_STRICT_PAPER_SCAN_STATUS,
    Paper07101869ScanRequest,
    Paper07101869StrictPaperScanRow,
    build_structural_scan_rows,
)


def test_qs1_structural_scan_rows_are_deterministic_and_noncanonical():
    request = Paper07101869ScanRequest(
        conventions=Paper07101869Conventions(),
        scale_points=(
            Paper07101869ScalePoint(
                label="qs1_lo",
                Lambda_IR_GeV=3000.0,
                m_g1_GeV=4300.0,
                mu_match_GeV=3000.0,
                mu_gs_GeV=4300.0,
            ),
            Paper07101869ScalePoint(
                label="qs1_lo_eff",
                Lambda_IR_GeV=3000.0,
                m_g1_GeV=4300.0,
                mu_match_GeV=3000.0,
                mu_gs_GeV=4300.0,
                m_KK_eff_GeV=4050.0,
            ),
        ),
    )

    rows = build_structural_scan_rows(request)
    rows_again = build_structural_scan_rows(request)

    expected_keys = {
        "point_id",
        "mode_id",
        "conventions_schema_id",
        "scales_schema_id",
        "Lambda_IR_GeV",
        "m_g1_GeV",
        "mu_match_GeV",
        "mu_gs_GeV",
        "m_KK_eff_GeV",
        "propagator_mass_GeV",
        "status",
        "note",
    }

    assert request.schema_id == PAPER_0710_1869_SCAN_SCHEMA_ID
    assert request.mode_id == PAPER_0710_1869_MODE_ID
    assert request.conventions.schema_id == PAPER_0710_1869_CONVENTIONS_SCHEMA_ID
    assert rows == rows_again
    assert [row.point_id for row in rows] == ["qs1_lo", "qs1_lo_eff"]
    assert all(row.mode_id == PAPER_0710_1869_MODE_ID for row in rows)
    assert all(row.status == "structural_only" for row in rows)
    assert all(
        "observable evaluation is intentionally not implemented" in row.note.lower()
        for row in rows
    )
    assert all(set(asdict(row)) == expected_keys for row in rows)
    assert rows[0].m_KK_eff_GeV is None
    assert rows[0].propagator_mass_GeV == pytest.approx(4300.0)
    assert rows[1].m_KK_eff_GeV == pytest.approx(4050.0)
    assert rows[1].propagator_mass_GeV == pytest.approx(4050.0)
    np.testing.assert_allclose(rows[0].Lambda_IR_GeV, rows_again[0].Lambda_IR_GeV)
    np.testing.assert_allclose(rows[1].mu_gs_GeV, rows_again[1].mu_gs_GeV)


def _strict_paper_scan_builder(module):
    for name in (
        "build_strict_paper_scan_rows",
        "build_paper_0710_1869_strict_paper_scan_rows",
        "build_paper_0710_1869_strict_scan_rows",
    ):
        candidate = getattr(module, name, None)
        if callable(candidate):
            return candidate
    pytest.fail("strict-paper scan path not exposed yet")


def _invoke_strict_paper_scan_builder(builder, request):
    parameters = inspect.signature(builder).parameters
    kwargs = {}
    for name, parameter in parameters.items():
        if parameter.kind in (inspect.Parameter.VAR_POSITIONAL, inspect.Parameter.VAR_KEYWORD):
            continue
        if parameter.default is inspect._empty:
            if name in ("request", "scan_request", "strict_paper_request"):
                kwargs[name] = request
            else:
                pytest.fail(
                    "strict-paper scan builder exposes an unsupported required parameter: "
                    f"{name}"
                )
    return builder(**kwargs)


def _row_payload(row):
    if hasattr(row, "as_dict") and callable(row.as_dict):
        return row.as_dict()
    return asdict(row)


def test_qs3_strict_paper_scan_rows_require_explicit_strict_claim_and_provenance() -> None:
    module = __import__("quarkConstraints.paper_0710_1869.scan", fromlist=["*"])
    builder = _strict_paper_scan_builder(module)
    request = Paper07101869ScanRequest(
        conventions=Paper07101869Conventions(),
        scale_points=(
            Paper07101869ScalePoint(
                label="qs3_strict_lo",
                Lambda_IR_GeV=3000.0,
                m_g1_GeV=4300.0,
                mu_match_GeV=3000.0,
                mu_gs_GeV=4300.0,
            ),
        ),
    )

    rows = _invoke_strict_paper_scan_builder(builder, request)
    rows_again = _invoke_strict_paper_scan_builder(builder, request)

    assert rows == rows_again
    first = _row_payload(rows[0])

    assert first["mode_id"] == PAPER_0710_1869_MODE_ID
    assert first["paper_id"] == "arXiv:0710.1869"
    assert first["claim_level_id"] != "structural_only"
    assert first["status"] == "strict_paper"
    assert first["provenance_policy_id"] == request.conventions.provenance_policy_id
    assert first["verifier_policy_id"] == request.conventions.verifier_policy_id
    assert first["benchmark_id"]
    assert first["benchmark_id"] != request.scale_points[0].label

    input_bundle_id = first.get("strict_paper_input_bundle_id") or first.get("input_bundle_id")
    input_provenance_id = first.get("strict_paper_input_provenance_id") or first.get(
        "input_provenance_id"
    )
    assert input_bundle_id
    assert input_provenance_id
    assert "strict_paper" in str(input_bundle_id)
    assert "strict_paper" in str(input_provenance_id)


def test_qs3_strict_paper_row_constructor_rejects_arbitrary_provenance_ids() -> None:
    request = Paper07101869ScanRequest(
        conventions=Paper07101869Conventions(),
        scale_points=(
            Paper07101869ScalePoint(
                label="qs3_strict_row",
                Lambda_IR_GeV=3000.0,
                m_g1_GeV=4300.0,
                mu_match_GeV=3000.0,
                mu_gs_GeV=4300.0,
            ),
        ),
    )
    point = request.scale_points[0]

    with pytest.raises(ValueError, match="benchmark_id"):
        Paper07101869StrictPaperScanRow(
            point_id=point.label,
            mode_id=request.mode_id,
            paper_id=request.conventions.paper_id,
            claim_level_id=PAPER_0710_1869_STRICT_PAPER_CLAIM_LEVEL_ID,
            conventions_schema_id=request.conventions.schema_id,
            scales_schema_id=point.schema_id,
            provenance_policy_id=request.conventions.provenance_policy_id,
            verifier_policy_id=request.conventions.verifier_policy_id,
            seed_to_profile_mapping_policy_id=request.conventions.seed_to_profile_mapping_policy_id,
            universal_term_coefficient_policy_id=(
                request.conventions.universal_term_coefficient_policy_id
            ),
            profile_derivation_policy_id=request.conventions.profile_derivation_policy_id,
            benchmark_id="arbitrary.benchmark.v999",
            input_bundle_id=PAPER_0710_1869_STRICT_PAPER_INPUT_BUNDLE_ID,
            input_provenance_id=PAPER_0710_1869_STRICT_PAPER_INPUT_PROVENANCE_ID,
            Lambda_IR_GeV=point.Lambda_IR_GeV,
            m_g1_GeV=point.m_g1_GeV,
            mu_match_GeV=point.mu_match_GeV,
            mu_gs_GeV=point.mu_gs_GeV,
            m_KK_eff_GeV=point.m_KK_eff_GeV,
            propagator_mass_GeV=point.propagator_mass_GeV,
        )

    with pytest.raises(ValueError, match="input_bundle_id"):
        Paper07101869StrictPaperScanRow(
            point_id=point.label,
            mode_id=request.mode_id,
            paper_id=request.conventions.paper_id,
            claim_level_id=PAPER_0710_1869_STRICT_PAPER_CLAIM_LEVEL_ID,
            conventions_schema_id=request.conventions.schema_id,
            scales_schema_id=point.schema_id,
            provenance_policy_id=request.conventions.provenance_policy_id,
            verifier_policy_id=request.conventions.verifier_policy_id,
            seed_to_profile_mapping_policy_id=request.conventions.seed_to_profile_mapping_policy_id,
            universal_term_coefficient_policy_id=(
                request.conventions.universal_term_coefficient_policy_id
            ),
            profile_derivation_policy_id=request.conventions.profile_derivation_policy_id,
            benchmark_id=PAPER_0710_1869_STRICT_PAPER_BENCHMARK_ID,
            input_bundle_id="arbitrary.input_bundle.v999",
            input_provenance_id=PAPER_0710_1869_STRICT_PAPER_INPUT_PROVENANCE_ID,
            Lambda_IR_GeV=point.Lambda_IR_GeV,
            m_g1_GeV=point.m_g1_GeV,
            mu_match_GeV=point.mu_match_GeV,
            mu_gs_GeV=point.mu_gs_GeV,
            m_KK_eff_GeV=point.m_KK_eff_GeV,
            propagator_mass_GeV=point.propagator_mass_GeV,
        )

    with pytest.raises(ValueError, match="input_provenance_id"):
        Paper07101869StrictPaperScanRow(
            point_id=point.label,
            mode_id=request.mode_id,
            paper_id=request.conventions.paper_id,
            claim_level_id=PAPER_0710_1869_STRICT_PAPER_CLAIM_LEVEL_ID,
            conventions_schema_id=request.conventions.schema_id,
            scales_schema_id=point.schema_id,
            provenance_policy_id=request.conventions.provenance_policy_id,
            verifier_policy_id=request.conventions.verifier_policy_id,
            seed_to_profile_mapping_policy_id=request.conventions.seed_to_profile_mapping_policy_id,
            universal_term_coefficient_policy_id=(
                request.conventions.universal_term_coefficient_policy_id
            ),
            profile_derivation_policy_id=request.conventions.profile_derivation_policy_id,
            benchmark_id=PAPER_0710_1869_STRICT_PAPER_BENCHMARK_ID,
            input_bundle_id=PAPER_0710_1869_STRICT_PAPER_INPUT_BUNDLE_ID,
            input_provenance_id="arbitrary.input_provenance.v999",
            Lambda_IR_GeV=point.Lambda_IR_GeV,
            m_g1_GeV=point.m_g1_GeV,
            mu_match_GeV=point.mu_match_GeV,
            mu_gs_GeV=point.mu_gs_GeV,
            m_KK_eff_GeV=point.m_KK_eff_GeV,
            propagator_mass_GeV=point.propagator_mass_GeV,
        )
