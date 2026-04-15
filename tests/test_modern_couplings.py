"""Focused tests for the modern QS2 coupling bridge."""

from __future__ import annotations

import numpy as np

from quarkConstraints.benchmarks import evaluate_default_benchmark
from quarkConstraints.couplings import compute_quark_kk_gluon_couplings
from quarkConstraints.modern.couplings import (
    MODERN_POINT_COUPLINGS_SCHEMA_ID,
    ModernPointCouplings,
    build_modern_point_couplings,
)
from quarkConstraints.modern.inputs import default_modern_default_inputs


def test_modern_point_couplings_matches_repo_formula_for_benchmark_point() -> None:
    fit_result = evaluate_default_benchmark().result
    inputs = default_modern_default_inputs()
    modern = build_modern_point_couplings(
        fit_result,
        point_id="benchmark-point",
        point_label="benchmark-label",
    )
    # The repo formula must be called with the same g_s_star that the modern
    # bundle uses (default 3.0) so the coupling matrices agree.
    repo = compute_quark_kk_gluon_couplings(
        fit_result, g_s_star=inputs.qcd_metadata.g_s_star,
    )

    assert modern.schema_id == MODERN_POINT_COUPLINGS_SCHEMA_ID
    assert modern.point_id == "benchmark-point"
    assert modern.point_label == "benchmark-label"
    assert modern.input_bundle_id == inputs.bundle_id
    assert modern.input_provenance_id == inputs.provenance_id
    assert modern.qcd_metadata_id == inputs.qcd_metadata.metadata_id
    assert modern.alpha_s_policy_id == inputs.qcd_metadata.alpha_s_policy_id
    assert modern.M_KK == repo.M_KK
    assert modern.alpha_s == repo.alpha_s
    assert modern.g_s == repo.g_s
    assert np.allclose(np.asarray(modern.left_up), repo.left_up)
    assert np.allclose(np.asarray(modern.left_down), repo.left_down)
    assert np.allclose(np.asarray(modern.right_up), repo.right_up)
    assert np.allclose(np.asarray(modern.right_down), repo.right_down)


def test_modern_point_couplings_round_trip_is_deterministic() -> None:
    fit_result = evaluate_default_benchmark().result
    modern = build_modern_point_couplings(fit_result)

    payload = modern.as_dict()
    assert ModernPointCouplings.from_dict(payload) == modern
    assert ModernPointCouplings.from_dict(payload).as_dict() == payload
