"""Focused tests for the modern QS2 matching bridge."""

from __future__ import annotations

import pytest

from quarkConstraints.benchmarks import evaluate_default_benchmark
from quarkConstraints.couplings import compute_quark_kk_gluon_couplings
from quarkConstraints.deltaf2 import DeltaF2Input, compute_delta_f2_wilsons
from quarkConstraints.modern.couplings import build_modern_point_couplings
from quarkConstraints.modern.inputs import ModernDefaultInputs
from quarkConstraints.modern.matching import (
    MODERN_POINT_MATCHING_SCHEMA_ID,
    ModernPointMatching,
    build_modern_point_matching,
)


def _backend_inputs_from_modern_inputs(inputs: ModernDefaultInputs) -> tuple[DeltaF2Input, ...]:
    weights = inputs.operator_weight_policy
    return tuple(
        DeltaF2Input(
            key=system.backend_key,
            display_name=system.display_name,
            column_name=system.column_name,
            reject_reason=system.reject_reason,
            sector=system.sector_id,
            generations=system.generations,
            bound=system.bound,
            ll_weight=weights.ll_weight,
            rr_weight=weights.rr_weight,
            lr1_weight=weights.lr1_weight,
            lr2_weight=weights.lr2_weight,
            reference_scale=weights.reference_scale_GeV,
            note=system.note,
        )
        for system in inputs.neutral_meson_inputs
    )


def test_modern_point_matching_matches_repo_formulas_for_benchmark_point() -> None:
    fit_result = evaluate_default_benchmark().result
    inputs = ModernDefaultInputs()
    modern_couplings = build_modern_point_couplings(fit_result)
    modern_matching = build_modern_point_matching(modern_couplings, inputs=inputs)
    repo_couplings = compute_quark_kk_gluon_couplings(fit_result)
    repo_matching = compute_delta_f2_wilsons(
        repo_couplings,
        inputs=_backend_inputs_from_modern_inputs(inputs),
    )

    assert modern_matching.schema_id == MODERN_POINT_MATCHING_SCHEMA_ID
    assert modern_matching.point_id == modern_couplings.point_id
    assert modern_matching.input_bundle_id == inputs.bundle_id
    assert modern_matching.input_provenance_id == inputs.provenance_id
    assert modern_matching.operator_basis_id == inputs.operator_weight_policy.operator_convention_id

    for modern_system, repo_system in zip(modern_matching.systems, repo_matching):
        assert modern_system.backend_key == repo_system.input.key
        assert modern_system.generations == repo_system.input.generations
        assert modern_system.left_coupling == pytest.approx(repo_system.left_coupling)
        assert modern_system.right_coupling == pytest.approx(repo_system.right_coupling)
        assert modern_system.c1_vll == pytest.approx(repo_system.c1_vll)
        assert modern_system.c1_vrr == pytest.approx(repo_system.c1_vrr)
        assert modern_system.c4_lr == pytest.approx(repo_system.c4_lr)
        assert modern_system.c5_lr == pytest.approx(repo_system.c5_lr)

    kaon = modern_matching.system("K")
    assert kaon.observable_id == "epsilon_K"
    assert kaon.backend_system_id == "K"
    assert kaon.backend_key == "epsilon_k"
    assert "epsilon_K" in kaon.note

    d0 = modern_matching.system("D0")
    assert d0.observable_id == "D0"
    assert d0.backend_system_id == "D"
    assert d0.backend_key == "d"
    assert "conservative" in d0.note


def test_modern_point_matching_round_trip_is_deterministic() -> None:
    fit_result = evaluate_default_benchmark().result
    matching = build_modern_point_matching(build_modern_point_couplings(fit_result))
    payload = matching.as_dict()

    assert ModernPointMatching.from_dict(payload) == matching
    assert ModernPointMatching.from_dict(payload).as_dict() == payload
