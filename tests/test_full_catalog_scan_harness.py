from __future__ import annotations

import importlib.util
import ast
import json
import sys
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from flavor_catalog_constraints.base import ConstraintResult, Severity


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = REPO_ROOT / "scripts" / "run_full_catalog_scan.py"


def _load_harness():
    spec = importlib.util.spec_from_file_location("run_full_catalog_scan", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_hard_veto_policy_separates_rigorous_proxy_partial_stub_and_advisory():
    harness = _load_harness()
    payload = harness._classify_results(
        {
            "R001": ConstraintResult(
                process_id="R001",
                severity=Severity.HARD,
                passes=False,
                ratio=2.0,
                diagnostics={"evaluated": True, "matching_status": "tree_level_exact"},
            ),
            "P001": ConstraintResult(
                process_id="P001",
                severity=Severity.HARD,
                passes=False,
                ratio=3.0,
                diagnostics={
                    "evaluated": True,
                    "used_proxy": True,
                    "needs_human_physics": "documented mass-limit proxy",
                },
            ),
            "H001": ConstraintResult(
                process_id="H001",
                severity=Severity.HARD,
                passes=False,
                ratio=4.0,
                diagnostics={
                    "evaluated": True,
                    "needs_human_physics": "NEEDS-HUMAN-PHYSICS: partial EFT matching",
                },
            ),
            "S001": ConstraintResult(
                process_id="S001",
                severity=Severity.HARD,
                passes=False,
                diagnostics={"evaluated": False, "missing_extra": "rs_loop_input"},
            ),
            "SOFT1": ConstraintResult(
                process_id="SOFT1",
                severity=Severity.SOFT,
                passes=False,
                ratio=9.0,
            ),
            "INFO1": ConstraintResult(
                process_id="INFO1",
                severity=Severity.INFO,
                passes=False,
                ratio=9.0,
            ),
        }
    )

    assert payload["excluded_by_rigorous"] == ["R001"]
    assert payload["excluded_by_proxy"] == ["P001"]
    assert payload["hard_not_evaluated"] == ["H001", "S001"]
    assert payload["survives_all_HARD_strict"] is False
    assert payload["survives_all_HARD_inclusive"] is False
    assert "SOFT:SOFT1" in payload["advisory_flags"]
    assert "INFO:INFO1" in payload["advisory_flags"]
    assert payload["constraints"]["H001"]["tag"] == "partial"
    assert payload["constraints"]["S001"]["tag"] == "stub"


def test_evaluated_proxy_with_unevaluated_subobservable_note_stays_proxy():
    harness = _load_harness()
    result = ConstraintResult(
        process_id="B013",
        severity=Severity.HARD,
        passes=False,
        ratio=2.0,
        notes="A_Delta and S_phi gamma are not evaluated by this branching proxy.",
        diagnostics={
            "evaluated": True,
            "mass_proxy": "exclusive C7 normalization proxy",
            "needs_human_physics": "documented C7 proxy",
        },
    )

    tag, _, _, proxy_flags = harness.tag_result(result)

    assert tag == "proxy"
    assert proxy_flags["mass_proxy"] == "exclusive C7 normalization proxy"


def test_plural_proxies_in_needs_human_is_tagged_proxy_not_partial():
    """M5 regression: the legacy ``"proxy" in needs_text`` missed the plural
    "proxies" (the T001/T002 t->qZ overlap-proxy wording), mis-tagging an
    evaluated, passing proxy as ``partial`` on 100% of rows.  The plural-robust
    match must now tag it ``proxy``."""
    harness = _load_harness()
    result = ConstraintResult(
        process_id="T001",
        severity=Severity.HARD,
        passes=True,
        ratio=1.0e-8,
        diagnostics={
            "evaluated": True,
            "needs_human_physics": (
                "NEEDS-HUMAN-PHYSICS: up-sector KK-gluon couplings are used as "
                "neutral-current flavor-overlap proxies; ..."
            ),
        },
    )

    tag, _, needs_human, _ = harness.tag_result(result)
    assert tag == "proxy"
    assert needs_human is not None


def test_structured_tag_class_hint_overrides_prose_matching():
    """M5: a constraint declaring ``diagnostics["tag_class"]`` is tagged by the
    structured hint, independent of any prose phrasing (the durable fix)."""
    harness = _load_harness()
    result = ConstraintResult(
        process_id="T002",
        severity=Severity.HARD,
        passes=True,
        ratio=1.0e-9,
        diagnostics={
            "evaluated": True,
            "tag_class": "proxy",
            # No "proxy"/"proxies" anywhere in the prose -- structured hint wins.
            "needs_human_physics": "documented overlap stand-in for the full match",
        },
    )

    tag, _, _, _ = harness.tag_result(result)
    assert tag == "proxy"

    # An unevaluated point still reports "stub" even with a tag_class hint.
    unevaluated = ConstraintResult(
        process_id="T002",
        severity=Severity.HARD,
        passes=True,
        diagnostics={"missing_extra": "quark_mass_basis_couplings", "tag_class": "proxy"},
    )
    tag_uneval, _, _, _ = harness.tag_result(unevaluated)
    assert tag_uneval == "stub"


def test_resolved_proxy_wording_in_matching_status_remains_rigorous():
    harness = _load_harness()
    result = ConstraintResult(
        process_id="B022",
        severity=Severity.HARD,
        passes=False,
        ratio=1.1,
        diagnostics={
            "evaluated": True,
            "rs_semileptonic_nunu_matching_status": (
                "Phase-4d rigorous active-neutrino Wilson; old one-Z-like "
                "proxy resolved and proxy is not used."
            ),
        },
    )

    tag, matching_status, needs_human, proxy_flags = harness.tag_result(result)

    assert tag == "rigorous"
    assert "proxy resolved" in matching_status
    assert needs_human is None
    assert proxy_flags == {}


def test_deferred_custodial_refinement_without_needs_human_stays_rigorous():
    harness = _load_harness()
    result = ConstraintResult(
        process_id="T010",
        severity=Severity.HARD,
        passes=False,
        ratio=2.0,
        diagnostics={
            "evaluated": True,
            "minimal_rs_tree_complete": True,
            "minimal_rs_tree_veto_ready": True,
            "custodial_variant_deferred": True,
            "custodial_variant_deferred_note": "Deferred refinement: custodial variant.",
        },
    )

    tag, matching_status, needs_human, proxy_flags = harness.tag_result(result)

    assert tag == "rigorous"
    assert matching_status is None
    assert needs_human is None
    assert proxy_flags == {}


def test_tile_seed_stride_and_config_hash_are_deterministic():
    harness = _load_harness()
    cfg = harness.ScanConfig(
        mkk_values_gev=(1000.0, 3000.0),
        n_draws_per_tile=5,
        xi_kk=2.0,
        base_seed=11,
        tile_seed_stride=17,
    )

    tiles = harness._build_tiles(cfg)
    same_hash = harness._config_hash(cfg)
    changed_hash = harness._config_hash(
        harness.ScanConfig(mkk_values_gev=(1000.0, 3000.0), n_draws_per_tile=6)
    )

    assert [tile.seed for tile in tiles] == [11, 28]
    assert [tile.lambda_ir_gev for tile in tiles] == [500.0, 1500.0]
    assert same_hash == harness._config_hash(cfg)
    assert same_hash != changed_hash


def test_default_tile_seed_stride_covers_large_smoke_tiles():
    harness = _load_harness()
    cfg = harness.ScanConfig(
        mkk_values_gev=(1000.0, 3000.0, 6000.0, 10000.0),
        n_draws_per_tile=10_000,
    )
    tiles = harness._build_tiles(cfg)

    assert harness.DEFAULT_TILE_SEED_STRIDE > cfg.n_draws_per_tile
    ranges = [
        set(range(tile.seed, tile.seed + tile.n_draws))
        for tile in tiles
    ]
    assert sum(len(seed_range) for seed_range in ranges) == len(set().union(*ranges))


def test_completed_summary_resume_requires_complete_matching_hash_and_draws(tmp_path):
    harness = _load_harness()
    path = tmp_path / "tile-00000.summary.json"
    path.write_text(
        json.dumps({"complete": True, "config_hash": "abc", "n_requested": 10, "n_rows": 10}),
        encoding="utf-8",
    )

    assert harness._completed_summary(path, expected_hash="abc", expected_draws=10) is not None
    assert harness._completed_summary(path, expected_hash="def", expected_draws=10) is None
    assert harness._completed_summary(path, expected_hash="abc", expected_draws=9) is None

    path.write_text(
        json.dumps({"complete": True, "config_hash": "abc", "n_requested": 10, "n_rows": 9}),
        encoding="utf-8",
    )
    assert harness._completed_summary(path, expected_hash="abc", expected_draws=10) is None

    path.write_text(
        json.dumps({"complete": False, "config_hash": "abc", "n_requested": 10, "n_rows": 10}),
        encoding="utf-8",
    )
    assert harness._completed_summary(path, expected_hash="abc", expected_draws=10) is None


def test_quark_only_config_payload_preserves_full_mode_hash_surface():
    harness = _load_harness()
    full_cfg = harness.ScanConfig(
        mkk_values_gev=(1000.0, 3000.0),
        n_draws_per_tile=5,
        xi_kk=2.0,
        base_seed=11,
        tile_seed_stride=17,
    )
    quark_cfg = harness.ScanConfig(
        mkk_values_gev=(1000.0, 3000.0),
        n_draws_per_tile=5,
        quark_only=True,
        xi_kk=2.0,
        base_seed=11,
        tile_seed_stride=17,
    )

    assert "quark_only" not in harness._config_payload(full_cfg)
    assert harness._config_payload(quark_cfg)["quark_only"] is True
    assert harness._config_hash(full_cfg) == "45e21a07585f7489"
    assert harness._config_hash(quark_cfg) == "d96cb734f724aedb"
    assert harness._config_hash(full_cfg) != harness._config_hash(quark_cfg)


def test_ew_model_config_payload_round_trip_and_pinned_minimal_hashes():
    harness = _load_harness()
    wq_cfg = harness.ScanConfig(
        mkk_values_gev=(
            1000.0,
            2000.0,
            3000.0,
            5000.0,
            7000.0,
            10000.0,
            15000.0,
            20000.0,
            30000.0,
            50000.0,
        ),
        n_draws_per_tile=2000,
        quark_only=True,
        base_seed=202606040000,
        tile_seed_stride=1000003,
        quark_fit_r=0.05,
    )
    custodial_cfg = harness.ScanConfig(
        mkk_values_gev=(1000.0, 3000.0),
        n_draws_per_tile=5,
        quark_only=True,
        ew_model=harness.CUSTODIAL_RS_PLR_EW_MODEL,
    )

    assert "ew_model" not in harness._config_payload(wq_cfg)
    assert harness._config_hash(wq_cfg) == "c6939cc65d71f86a"
    assert harness._config_from_payload(harness._config_payload(wq_cfg)).ew_model == (
        harness.MINIMAL_RS_EW_MODEL
    )
    custodial_payload = harness._config_payload(custodial_cfg)
    assert custodial_payload["ew_model"] == harness.CUSTODIAL_RS_PLR_EW_MODEL
    assert harness._config_from_payload(custodial_payload).ew_model == (
        harness.CUSTODIAL_RS_PLR_EW_MODEL
    )


def test_ew_model_argparse_choices():
    harness = _load_harness()
    parser = harness._build_argparser()

    args = parser.parse_args(
        [
            "--output-dir",
            "out",
            "--ew-model",
            harness.CUSTODIAL_RS_PLR_EW_MODEL,
        ]
    )

    assert args.ew_model == harness.CUSTODIAL_RS_PLR_EW_MODEL
    with pytest.raises(SystemExit):
        parser.parse_args(["--output-dir", "out", "--ew-model", "nonsense"])


def test_quark_only_allowlist_matches_candidate_verification_and_drops_deferred_leptons():
    harness = _load_harness()
    expected_collider_in = (
        "CR001",
        "CR002",
        "CR003",
        "CR004",
        "CR007",
        "CR008",
        "CR010",
        "CR012",
        "CR013",
    )
    expected_collider_deferred = ("CR005", "CR006", "CR009", "CR011", "CR014")

    assert set(harness.QUARK_ONLY_CANDIDATE_IDS) == (
        set(harness.QUARK_ONLY_ALLOWLIST_IDS)
        | set(harness.QUARK_ONLY_DEFERRED_LEPTON_FOLLOWUP)
    )
    assert len(harness.QUARK_ONLY_ALLOWLIST_IDS) == 46
    assert set(expected_collider_in) <= harness.QUARK_ONLY_ALLOWLIST_SET
    assert set(expected_collider_deferred).isdisjoint(harness.QUARK_ONLY_ALLOWLIST_SET)
    assert harness.QUARK_ONLY_DEFERRED_LEPTON_FOLLOWUP == (
        "EW002",
        *expected_collider_deferred,
    )
    assert "EW002" not in harness.QUARK_ONLY_ALLOWLIST_SET
    assert "lepton_lmfv_parameters" in harness.QUARK_ONLY_FORBIDDEN_EXTRAS

    actual = _candidate_get_extra_usage(harness)
    for pid in harness.QUARK_ONLY_ALLOWLIST_IDS:
        assert tuple(actual[pid]) == harness.QUARK_ONLY_ALLOWLIST_EXTRAS[pid]
        forbidden = set(actual[pid]) & harness.QUARK_ONLY_FORBIDDEN_EXTRAS
        if pid not in harness.QUARK_ONLY_OPTIONAL_LEPTON_DIAGNOSTIC_IDS:
            assert forbidden == set()
    for pid in harness.QUARK_ONLY_DEFERRED_LEPTON_FOLLOWUP:
        assert actual[pid] == harness.QUARK_ONLY_DEFERRED_EXTRAS[pid]
    for pid in expected_collider_in:
        assert not set(actual[pid]) & harness.QUARK_ONLY_FORBIDDEN_EXTRAS


def test_quark_only_collider_allowlist_is_active_proxy_and_monotonic():
    from flavor_catalog_constraints import point_builder

    harness = _load_harness()
    collider_ids = (
        "CR001",
        "CR002",
        "CR003",
        "CR004",
        "CR007",
        "CR008",
        "CR010",
        "CR012",
        "CR013",
    )

    def point(m_kk_gev: float):
        return point_builder.make_point(
            kk_ew_mass_gev=m_kk_gev,
            kk_gluon_mass_gev=m_kk_gev,
            quark_mass_basis_couplings=SimpleNamespace(M_KK=m_kk_gev, xi_KK=1.0),
        )

    low_results = harness._evaluate_constraint_ids(point(1000.0), collider_ids)
    high_results = harness._evaluate_constraint_ids(point(50000.0), collider_ids)
    low_payload = harness._classify_results(low_results)
    high_payload = harness._classify_results(high_results)

    assert low_payload["excluded_by_proxy"]
    assert low_payload["survives_all_HARD_inclusive"] is False
    assert high_payload["excluded_by_proxy"] == []
    assert high_payload["survives_all_HARD_inclusive"] is True

    for pid in collider_ids:
        result = low_results[pid]
        high_result = high_results[pid]
        diagnostics = dict(result.diagnostics)
        assert result.severity is Severity.HARD
        assert low_payload["constraints"][pid]["evaluated"] is True
        assert low_payload["constraints"][pid]["active"] is True
        assert low_payload["constraints"][pid]["tag"] == "proxy"
        assert "missing_extra" not in diagnostics
        assert "missing_extras" not in diagnostics
        assert result.ratio is not None
        assert high_result.ratio is not None
        assert result.ratio > high_result.ratio
        assert high_result.passes is True
        if pid in {"CR001", "CR007", "CR012", "CR013"}:
            assert diagnostics["mass_source"] in {
                "kk_gluon_mass_gev",
                "kk_ew_mass_gev",
            }
        else:
            assert diagnostics["mass_source"] == "quark_mass_basis_couplings.M_KK"


def test_quark_only_evaluate_draw_skips_leptons_filters_allowlist_and_is_deterministic(monkeypatch):
    harness = _load_harness()
    cfg = harness.ScanConfig(
        mkk_values_gev=(3000.0,),
        n_draws_per_tile=1,
        quark_only=True,
        quark_fit_max_nfev=1,
    )
    tile = harness.TileSpec(
        tile_id=0,
        mkk_gev=3000.0,
        lambda_ir_gev=3000.0 / cfg.xi_kk,
        n_draws=1,
        seed=123,
    )
    spectrum = object()

    @dataclass
    class _Cache:
        max_a_rel_err: float = 0.0

    @dataclass
    class _BulkState:
        c_Q: np.ndarray
        c_u: np.ndarray
        c_d: np.ndarray

    @dataclass
    class _Point:
        Y_u: np.ndarray
        Y_d: np.ndarray

    @dataclass
    class _FitResult:
        bulk_state: _BulkState
        point: _Point
        score: float
        residual_norm: float
        mass_residuals_up: np.ndarray
        mass_residuals_down: np.ndarray
        ckm_residuals: np.ndarray

    @dataclass
    class _Solution:
        result: _FitResult
        success: bool = True
        message: str = "ok"
        nfev: int = 1
        initial_score: float = 1.0

    def fake_fit(*args, **kwargs):
        assert kwargs["fit_orientation"] is True
        assert kwargs["r"] == cfg.quark_fit_r
        fit_result = _FitResult(
            bulk_state=_BulkState(
                c_Q=np.array([0.61, 0.62, 0.63]),
                c_u=np.array([0.64, 0.65, 0.66]),
                c_d=np.array([0.67, 0.68, 0.69]),
            ),
            point=_Point(
                Y_u=np.ones((3, 3), dtype=np.complex128),
                Y_d=np.ones((3, 3), dtype=np.complex128),
            ),
            score=0.0,
            residual_norm=0.0,
            mass_residuals_up=np.zeros(3),
            mass_residuals_down=np.zeros(3),
            ckm_residuals=np.zeros(4),
        )
        return _Solution(result=fit_result)

    def forbidden(*_args, **_kwargs):
        raise AssertionError("lepton stage must not run in quark-only mode")

    def fake_build(fit_result, **kwargs):
        assert kwargs["spectrum"] is spectrum
        assert kwargs["rs_ew_cache"] is cache
        assert kwargs["lepton_yukawa_result"] is None
        assert kwargs["include_charged_current"] is False
        assert kwargs["include_fermion_kk_mixing"] is True
        assert kwargs["include_higgs_yukawas"] is False
        point = harness.point_builder.make_point(
            raw=kwargs["raw"],
            kk_ew_mass_gev=3000.0,
            rs_ew_spectrum=spectrum,
            rs_ew_couplings=object(),
        )
        assert "lepton_lmfv_parameters" not in point.extras
        return point

    @dataclass
    class _Couplings:
        M_KK: float

    def fake_couplings(_fit_result, **kwargs):
        return _Couplings(M_KK=float(kwargs["M_KK"]))

    captured_ids: list[tuple[str, ...]] = []

    def fake_evaluate_ids(_point, process_ids):
        ids = tuple(process_ids)
        captured_ids.append(ids)
        return {
            pid: ConstraintResult(
                process_id=pid,
                severity=Severity.HARD,
                passes=pid != "K001",
                ratio=2.0 if pid == "K001" else 0.1,
                diagnostics={"evaluated": True},
            )
            for pid in ids
        }

    cache = _Cache()
    monkeypatch.setattr(harness, "fit_quark_sector", fake_fit)
    monkeypatch.setattr(harness, "_draw_lepton_inputs", forbidden)
    monkeypatch.setattr(harness, "compute_all_yukawas", forbidden)
    monkeypatch.setattr(harness, "_require_perturbative_leptons", forbidden)
    monkeypatch.setattr(harness.point_builder, "build_from_rs_ew_inputs", fake_build)
    monkeypatch.setattr(harness, "compute_quark_kk_gluon_couplings", fake_couplings)
    monkeypatch.setattr(harness, "_evaluate_constraint_ids", fake_evaluate_ids)

    kwargs = dict(
        tile=tile,
        draw_idx=0,
        draw_seed=123,
        cfg=cfg,
        spectrum=spectrum,
        overlap_cache=cache,
        provenance={"mode": "quark_only", "quark_fit_r": cfg.quark_fit_r},
        registry_count=harness.EXPECTED_REGISTRY_COUNT,
        config_hash="abc",
    )
    row1 = harness._evaluate_draw(np.random.default_rng(123), **kwargs)
    row2 = harness._evaluate_draw(np.random.default_rng(123), **kwargs)

    assert row1["skipped"] is False
    assert row1["mode"] == "quark_only"
    assert row1["quark_fit_r"] == cfg.quark_fit_r
    assert row1["params"]["quark_fit_r"] == cfg.quark_fit_r
    assert row1["provenance"]["quark_fit_r"] == cfg.quark_fit_r
    assert row1["provenance"]["config_hash"] == "abc"
    assert "lepton_inputs" not in row1["params"]
    assert row1["lepton_sector"] == "dropped (not rigorous)"
    assert "lepton_lmfv_parameters" not in row1["params"]
    assert "fitted_up_yukawa_singular_values" in row1["fit_diagnostics"]
    assert "fitted_down_yukawa_singular_values" in row1["fit_diagnostics"]
    np.testing.assert_allclose(
        row1["fit_diagnostics"]["fitted_up_yukawa_singular_values"],
        [0.0, 0.0, 3.0],
        atol=1.0e-12,
    )
    assert set(row1["constraints"]) == harness.QUARK_ONLY_ALLOWLIST_SET
    assert "EW002" not in row1["constraints"]
    assert captured_ids[-1] == harness.QUARK_ONLY_ALLOWLIST_IDS
    assert row1["params"] == row2["params"]
    assert row1["constraints"] == row2["constraints"]
    serialized1 = json.dumps(
        harness._json_sanitize(row1),
        sort_keys=True,
        separators=(",", ":"),
    ).encode()
    serialized2 = json.dumps(
        harness._json_sanitize(row2),
        sort_keys=True,
        separators=(",", ":"),
    ).encode()
    expected_serialized = '{"advisory_flags":["deferred_lepton_followup"],"allowlist":["K001","K002","B001","B002","B003","B004","C001","C002","B011","B012","B013","B014","B032","B033","B034","C003","K003","K013","T001","T002","T003","T004","T005","T006","T007","T008","T010","T011","T012","T014","EW001","EW003","CR001","CR002","CR003","CR004","CR007","CR008","CR010","CR012","CR013","E004","E006","E007","E008","E009"],"cache_metrics":{"rs_ew_cache_injected":true,"spectrum_injected":true,"spline_max_a_rel_err":0.0},"constraints":{"B001":{"active":true,"evaluated":true,"match_status":null,"needs_human_physics":null,"passes":true,"proxy_flags":{},"ratio":0.1,"severity":"HARD","tag":"rigorous"},"B002":{"active":true,"evaluated":true,"match_status":null,"needs_human_physics":null,"passes":true,"proxy_flags":{},"ratio":0.1,"severity":"HARD","tag":"rigorous"},"B003":{"active":true,"evaluated":true,"match_status":null,"needs_human_physics":null,"passes":true,"proxy_flags":{},"ratio":0.1,"severity":"HARD","tag":"rigorous"},"B004":{"active":true,"evaluated":true,"match_status":null,"needs_human_physics":null,"passes":true,"proxy_flags":{},"ratio":0.1,"severity":"HARD","tag":"rigorous"},"B011":{"active":true,"evaluated":true,"match_status":null,"needs_human_physics":null,"passes":true,"proxy_flags":{},"ratio":0.1,"severity":"HARD","tag":"rigorous"},"B012":{"active":true,"evaluated":true,"match_status":null,"needs_human_physics":null,"passes":true,"proxy_flags":{},"ratio":0.1,"severity":"HARD","tag":"rigorous"},"B013":{"active":true,"evaluated":true,"match_status":null,"needs_human_physics":null,"passes":true,"proxy_flags":{},"ratio":0.1,"severity":"HARD","tag":"rigorous"},"B014":{"active":true,"evaluated":true,"match_status":null,"needs_human_physics":null,"passes":true,"proxy_flags":{},"ratio":0.1,"severity":"HARD","tag":"rigorous"},"B032":{"active":true,"evaluated":true,"match_status":null,"needs_human_physics":null,"passes":true,"proxy_flags":{},"ratio":0.1,"severity":"HARD","tag":"rigorous"},"B033":{"active":true,"evaluated":true,"match_status":null,"needs_human_physics":null,"passes":true,"proxy_flags":{},"ratio":0.1,"severity":"HARD","tag":"rigorous"},"B034":{"active":true,"evaluated":true,"match_status":null,"needs_human_physics":null,"passes":true,"proxy_flags":{},"ratio":0.1,"severity":"HARD","tag":"rigorous"},"C001":{"active":true,"evaluated":true,"match_status":null,"needs_human_physics":null,"passes":true,"proxy_flags":{},"ratio":0.1,"severity":"HARD","tag":"rigorous"},"C002":{"active":true,"evaluated":true,"match_status":null,"needs_human_physics":null,"passes":true,"proxy_flags":{},"ratio":0.1,"severity":"HARD","tag":"rigorous"},"C003":{"active":true,"evaluated":true,"match_status":null,"needs_human_physics":null,"passes":true,"proxy_flags":{},"ratio":0.1,"severity":"HARD","tag":"rigorous"},"CR001":{"active":true,"evaluated":true,"match_status":null,"needs_human_physics":null,"passes":true,"proxy_flags":{},"ratio":0.1,"severity":"HARD","tag":"rigorous"},"CR002":{"active":true,"evaluated":true,"match_status":null,"needs_human_physics":null,"passes":true,"proxy_flags":{},"ratio":0.1,"severity":"HARD","tag":"rigorous"},"CR003":{"active":true,"evaluated":true,"match_status":null,"needs_human_physics":null,"passes":true,"proxy_flags":{},"ratio":0.1,"severity":"HARD","tag":"rigorous"},"CR004":{"active":true,"evaluated":true,"match_status":null,"needs_human_physics":null,"passes":true,"proxy_flags":{},"ratio":0.1,"severity":"HARD","tag":"rigorous"},"CR007":{"active":true,"evaluated":true,"match_status":null,"needs_human_physics":null,"passes":true,"proxy_flags":{},"ratio":0.1,"severity":"HARD","tag":"rigorous"},"CR008":{"active":true,"evaluated":true,"match_status":null,"needs_human_physics":null,"passes":true,"proxy_flags":{},"ratio":0.1,"severity":"HARD","tag":"rigorous"},"CR010":{"active":true,"evaluated":true,"match_status":null,"needs_human_physics":null,"passes":true,"proxy_flags":{},"ratio":0.1,"severity":"HARD","tag":"rigorous"},"CR012":{"active":true,"evaluated":true,"match_status":null,"needs_human_physics":null,"passes":true,"proxy_flags":{},"ratio":0.1,"severity":"HARD","tag":"rigorous"},"CR013":{"active":true,"evaluated":true,"match_status":null,"needs_human_physics":null,"passes":true,"proxy_flags":{},"ratio":0.1,"severity":"HARD","tag":"rigorous"},"E004":{"active":true,"evaluated":true,"match_status":null,"needs_human_physics":null,"passes":true,"proxy_flags":{},"ratio":0.1,"severity":"HARD","tag":"rigorous"},"E006":{"active":true,"evaluated":true,"match_status":null,"needs_human_physics":null,"passes":true,"proxy_flags":{},"ratio":0.1,"severity":"HARD","tag":"rigorous"},"E007":{"active":true,"evaluated":true,"match_status":null,"needs_human_physics":null,"passes":true,"proxy_flags":{},"ratio":0.1,"severity":"HARD","tag":"rigorous"},"E008":{"active":true,"evaluated":true,"match_status":null,"needs_human_physics":null,"passes":true,"proxy_flags":{},"ratio":0.1,"severity":"HARD","tag":"rigorous"},"E009":{"active":true,"evaluated":true,"match_status":null,"needs_human_physics":null,"passes":true,"proxy_flags":{},"ratio":0.1,"severity":"HARD","tag":"rigorous"},"EW001":{"active":true,"evaluated":true,"match_status":null,"needs_human_physics":null,"passes":true,"proxy_flags":{},"ratio":0.1,"severity":"HARD","tag":"rigorous"},"EW003":{"active":true,"evaluated":true,"match_status":null,"needs_human_physics":null,"passes":true,"proxy_flags":{},"ratio":0.1,"severity":"HARD","tag":"rigorous"},"K001":{"active":true,"evaluated":true,"match_status":null,"needs_human_physics":null,"passes":false,"proxy_flags":{},"ratio":2.0,"severity":"HARD","tag":"rigorous"},"K002":{"active":true,"evaluated":true,"match_status":null,"needs_human_physics":null,"passes":true,"proxy_flags":{},"ratio":0.1,"severity":"HARD","tag":"rigorous"},"K003":{"active":true,"evaluated":true,"match_status":null,"needs_human_physics":null,"passes":true,"proxy_flags":{},"ratio":0.1,"severity":"HARD","tag":"rigorous"},"K013":{"active":true,"evaluated":true,"match_status":null,"needs_human_physics":null,"passes":true,"proxy_flags":{},"ratio":0.1,"severity":"HARD","tag":"rigorous"},"T001":{"active":true,"evaluated":true,"match_status":null,"needs_human_physics":null,"passes":true,"proxy_flags":{},"ratio":0.1,"severity":"HARD","tag":"rigorous"},"T002":{"active":true,"evaluated":true,"match_status":null,"needs_human_physics":null,"passes":true,"proxy_flags":{},"ratio":0.1,"severity":"HARD","tag":"rigorous"},"T003":{"active":true,"evaluated":true,"match_status":null,"needs_human_physics":null,"passes":true,"proxy_flags":{},"ratio":0.1,"severity":"HARD","tag":"rigorous"},"T004":{"active":true,"evaluated":true,"match_status":null,"needs_human_physics":null,"passes":true,"proxy_flags":{},"ratio":0.1,"severity":"HARD","tag":"rigorous"},"T005":{"active":true,"evaluated":true,"match_status":null,"needs_human_physics":null,"passes":true,"proxy_flags":{},"ratio":0.1,"severity":"HARD","tag":"rigorous"},"T006":{"active":true,"evaluated":true,"match_status":null,"needs_human_physics":null,"passes":true,"proxy_flags":{},"ratio":0.1,"severity":"HARD","tag":"rigorous"},"T007":{"active":true,"evaluated":true,"match_status":null,"needs_human_physics":null,"passes":true,"proxy_flags":{},"ratio":0.1,"severity":"HARD","tag":"rigorous"},"T008":{"active":true,"evaluated":true,"match_status":null,"needs_human_physics":null,"passes":true,"proxy_flags":{},"ratio":0.1,"severity":"HARD","tag":"rigorous"},"T010":{"active":true,"evaluated":true,"match_status":null,"needs_human_physics":null,"passes":true,"proxy_flags":{},"ratio":0.1,"severity":"HARD","tag":"rigorous"},"T011":{"active":true,"evaluated":true,"match_status":null,"needs_human_physics":null,"passes":true,"proxy_flags":{},"ratio":0.1,"severity":"HARD","tag":"rigorous"},"T012":{"active":true,"evaluated":true,"match_status":null,"needs_human_physics":null,"passes":true,"proxy_flags":{},"ratio":0.1,"severity":"HARD","tag":"rigorous"},"T014":{"active":true,"evaluated":true,"match_status":null,"needs_human_physics":null,"passes":true,"proxy_flags":{},"ratio":0.1,"severity":"HARD","tag":"rigorous"}},"coverage_complete":true,"deferred_lepton_followup":["EW002","CR005","CR006","CR009","CR011","CR014"],"draw_id":0,"excluded_by_proxy":[],"excluded_by_rigorous":["K001"],"fit_diagnostics":{"bulk_c_Q":[0.61,0.62,0.63],"bulk_c_d":[0.67,0.68,0.69],"bulk_c_u":[0.64,0.65,0.66],"fitted_down_yukawa_singular_values":[0.0,0.0,3.0],"fitted_up_yukawa_singular_values":[0.0,0.0,3.0],"initial_score":1.0,"max_abs_quark_yukawa":1.0,"max_ckm_residual":0.0,"max_down_log_residual":0.0,"max_up_log_residual":0.0,"message":"ok","nfev":1,"residual_norm":0.0,"score":0.0,"success":true},"hard_not_evaluated":[],"lepton_sector":"dropped (not rigorous)","mode":"quark_only","params":{"Lambda_IR":1225.1398701351736,"M_KK":3000.0,"k":1.2209e+19,"quark_fit_r":0.25,"quark_yukawa_seed":{"Y_d_im":[[0.6847299845497847,1.25580147532063,0.37660201720639197],[1.251367717642807,1.0940707532813247,-0.8455713802500362],[1.0983822922147288,0.6922558091137612,-0.6664041291032166]],"Y_d_re":[[0.8973753858602489,0.0544951105581426,-0.8053331255487975],[-1.0022880202777662,-0.00663309450661842,0.24817392184595977],[-0.9469860377145608,-1.4553152497193032,-0.0866003133286175]],"Y_u_im":[[1.1696780793335577,0.038911365688595545,-0.7651061967936106],[0.9727247882922341,-0.8587111098747136,0.7244011567041291],[0.38982061376904253,1.282221775575501,-0.8042754341807435]],"Y_u_re":[[0.5470555897444305,-1.3385369435933319,-0.8389203816821659],[-0.9468845679039909,-0.972282296744909,0.9362835199673212],[1.2700349940811693,-0.6702768066086813,0.9592636847790064]]},"xi_KK":2.4487},"provenance":{"config_hash":"abc","mode":"quark_only","quark_fit_r":0.25,"registry_count":103},"quark_fit_r":0.25,"seed":123,"skipped":false,"survives_all_HARD_inclusive":false,"survives_all_HARD_strict":false,"tile_id":0}'.encode()
    assert serialized1 == serialized2
    assert serialized1 == expected_serialized


def test_evaluate_draw_threads_custodial_ew_model_to_quark_only_and_full_builders(monkeypatch):
    harness = _load_harness()
    spectrum = object()

    @dataclass
    class _Cache:
        max_a_rel_err: float = 0.0

    @dataclass
    class _BulkState:
        c_Q: np.ndarray
        c_u: np.ndarray
        c_d: np.ndarray

    @dataclass
    class _Point:
        Y_u: np.ndarray
        Y_d: np.ndarray

    @dataclass
    class _FitResult:
        bulk_state: _BulkState
        point: _Point
        score: float
        residual_norm: float
        mass_residuals_up: np.ndarray
        mass_residuals_down: np.ndarray
        ckm_residuals: np.ndarray

    @dataclass
    class _Solution:
        result: _FitResult
        success: bool = True
        message: str = "ok"
        nfev: int = 1
        initial_score: float = 1.0

    @dataclass
    class _Leptons:
        Y_E_bar: np.ndarray
        Y_N_bar: np.ndarray

        def is_perturbative(self, _max_ybar):
            return True

    @dataclass
    class _Couplings:
        M_KK: float

    def fake_fit(*args, **kwargs):
        fit_result = _FitResult(
            bulk_state=_BulkState(
                c_Q=np.array([0.61, 0.62, 0.63]),
                c_u=np.array([0.64, 0.65, 0.66]),
                c_d=np.array([0.67, 0.68, 0.69]),
            ),
            point=_Point(
                Y_u=np.ones((3, 3), dtype=np.complex128),
                Y_d=np.ones((3, 3), dtype=np.complex128),
            ),
            score=0.0,
            residual_norm=0.0,
            mass_residuals_up=np.zeros(3),
            mass_residuals_down=np.zeros(3),
            ckm_residuals=np.zeros(4),
        )
        return _Solution(result=fit_result)

    captured_builds: list[dict[str, object]] = []

    def fake_build(_fit_result, **kwargs):
        captured_builds.append(dict(kwargs))
        return harness.point_builder.make_point(
            raw=kwargs["raw"],
            kk_ew_mass_gev=3000.0,
            rs_ew_spectrum=spectrum,
            rs_ew_couplings=object(),
        )

    def fake_couplings(_fit_result, **kwargs):
        return _Couplings(M_KK=float(kwargs["M_KK"]))

    def fake_results(process_ids):
        return {
            pid: ConstraintResult(
                process_id=pid,
                severity=Severity.HARD,
                passes=True,
                ratio=0.1,
                diagnostics={"evaluated": True},
            )
            for pid in process_ids
        }

    monkeypatch.setattr(harness, "fit_quark_sector", fake_fit)
    monkeypatch.setattr(harness, "compute_all_yukawas", lambda *a, **k: _Leptons(np.ones(3), np.ones(3)))
    monkeypatch.setattr(harness.point_builder, "build_from_rs_ew_inputs", fake_build)
    monkeypatch.setattr(harness, "compute_quark_kk_gluon_couplings", fake_couplings)
    monkeypatch.setattr(
        harness,
        "_draw_lepton_inputs",
        lambda _rng, _cfg: {
            "c_L": 0.61,
            "c_E": [0.62, 0.63, 0.64],
            "c_N": 0.65,
            "M_N": 1.0e12,
            "lightest_nu_mass": 0.001,
            "ordering": "normal",
            "majorana_alpha": 0.0,
            "majorana_beta": 0.0,
        },
    )
    monkeypatch.setattr(
        harness,
        "_evaluate_constraint_ids",
        lambda _point, ids: fake_results(ids),
    )
    monkeypatch.setattr(
        harness.registry,
        "evaluate_all",
        lambda _point: fake_results(("EW001",)),
    )

    for idx, quark_only in enumerate((True, False)):
        cfg = harness.ScanConfig(
            mkk_values_gev=(3000.0,),
            n_draws_per_tile=1,
            quark_only=quark_only,
            ew_model=harness.CUSTODIAL_RS_PLR_EW_MODEL,
            quark_fit_max_nfev=1,
        )
        tile = harness.TileSpec(
            tile_id=0,
            mkk_gev=3000.0,
            lambda_ir_gev=3000.0 / cfg.xi_kk,
            n_draws=1,
            seed=900 + idx,
        )
        row = harness._evaluate_draw(
            np.random.default_rng(900 + idx),
            tile=tile,
            draw_idx=0,
            draw_seed=900 + idx,
            cfg=cfg,
            spectrum=spectrum,
            overlap_cache=_Cache(),
            provenance={},
            registry_count=harness.EXPECTED_REGISTRY_COUNT,
            config_hash="cust",
        )

        assert row["skipped"] is False
        assert row["ew_model"] == harness.CUSTODIAL_RS_PLR_EW_MODEL
        assert row["provenance"]["ew_model"] == harness.CUSTODIAL_RS_PLR_EW_MODEL

    assert [kwargs["ew_model"] for kwargs in captured_builds] == [
        harness.CUSTODIAL_RS_PLR_EW_MODEL,
        harness.CUSTODIAL_RS_PLR_EW_MODEL,
    ]
    assert all(kwargs["spectrum"] is spectrum for kwargs in captured_builds)


def test_nonperturbative_leptons_skip_before_l001_evaluation(monkeypatch):
    harness = _load_harness()
    cfg = harness.ScanConfig(
        mkk_values_gev=(3000.0,),
        n_draws_per_tile=1,
        quark_fit_max_nfev=1,
        perturbative_ybar_max=4.0,
    )
    tile = harness.TileSpec(
        tile_id=0,
        mkk_gev=3000.0,
        lambda_ir_gev=3000.0 / cfg.xi_kk,
        n_draws=1,
        seed=456,
    )
    spectrum = object()

    @dataclass
    class _Cache:
        max_a_rel_err: float = 0.0

    @dataclass
    class _BulkState:
        c_Q: np.ndarray
        c_u: np.ndarray
        c_d: np.ndarray

    @dataclass
    class _Point:
        Y_u: np.ndarray
        Y_d: np.ndarray

    @dataclass
    class _FitResult:
        bulk_state: _BulkState
        point: _Point
        score: float
        residual_norm: float
        mass_residuals_up: np.ndarray
        mass_residuals_down: np.ndarray
        ckm_residuals: np.ndarray

    @dataclass
    class _Solution:
        result: _FitResult
        success: bool = True
        message: str = "ok"
        nfev: int = 1
        initial_score: float = 1.0

    @dataclass
    class _NonperturbativeLeptons:
        Y_E_bar: np.ndarray
        Y_N_bar: np.ndarray

        def is_perturbative(self, max_ybar):
            assert max_ybar == cfg.perturbative_ybar_max
            return False

    def fake_fit(*args, **kwargs):
        fit_result = _FitResult(
            bulk_state=_BulkState(
                c_Q=np.array([0.61, 0.62, 0.63]),
                c_u=np.array([0.64, 0.65, 0.66]),
                c_d=np.array([0.67, 0.68, 0.69]),
            ),
            point=_Point(
                Y_u=np.ones((3, 3), dtype=np.complex128),
                Y_d=np.ones((3, 3), dtype=np.complex128),
            ),
            score=0.0,
            residual_norm=0.0,
            mass_residuals_up=np.zeros(3),
            mass_residuals_down=np.zeros(3),
            ckm_residuals=np.zeros(4),
        )
        return _Solution(result=fit_result)

    def fake_compute_all_yukawas(*args, **kwargs):
        return _NonperturbativeLeptons(
            Y_E_bar=np.array([0.1, 0.2, 0.3]),
            Y_N_bar=np.array([0.2, 4.5, 0.7]),
        )

    def forbidden(*_args, **_kwargs):
        raise AssertionError("point building and L001 evaluation must not run")

    monkeypatch.setattr(harness, "fit_quark_sector", fake_fit)
    monkeypatch.setattr(harness, "compute_all_yukawas", fake_compute_all_yukawas)
    monkeypatch.setattr(harness.point_builder, "build_from_rs_ew_inputs", forbidden)
    monkeypatch.setattr(harness.registry, "evaluate_all", forbidden)

    row = harness._evaluate_draw(
        np.random.default_rng(456),
        tile=tile,
        draw_idx=0,
        draw_seed=456,
        cfg=cfg,
        spectrum=spectrum,
        overlap_cache=_Cache(),
        provenance={"mode": "full"},
        registry_count=harness.EXPECTED_REGISTRY_COUNT,
        config_hash="def",
    )

    assert row["skipped"] is True
    assert row["skip_reason"] == "nonperturbative_lepton_yukawa"
    assert "nonperturbative_lepton_yukawa" in row["skip_error"]
    assert row["constraints"] == {}
    assert "lepton_inputs" in row["params"]


def test_universal_c_sanity_uses_cfg_ew_model_for_spectrum_and_builder(monkeypatch):
    harness = _load_harness()
    captured = {"model_label": None, "ew_model": None}

    @dataclass
    class _Spectrum:
        pass

    @dataclass
    class _Cache:
        pass

    @dataclass
    class _BulkState:
        c_Q: np.ndarray
        c_u: np.ndarray
        c_d: np.ndarray
        F_Q: np.ndarray
        F_u: np.ndarray
        F_d: np.ndarray

    @dataclass
    class _FitResult:
        bulk_state: _BulkState
        U_L_u: np.ndarray
        U_L_d: np.ndarray
        U_R_u: np.ndarray
        U_R_d: np.ndarray
        ckm: np.ndarray

    @dataclass
    class _Solution:
        result: _FitResult

    fit = _FitResult(
        bulk_state=_BulkState(
            c_Q=np.ones(3),
            c_u=np.ones(3),
            c_d=np.ones(3),
            F_Q=np.ones(3),
            F_u=np.ones(3),
            F_d=np.ones(3),
        ),
        U_L_u=np.eye(3),
        U_L_d=np.eye(3),
        U_R_u=np.eye(3),
        U_R_d=np.eye(3),
        ckm=np.eye(3),
    )
    spectrum = _Spectrum()

    def fake_spectrum_build(**kwargs):
        captured["model_label"] = kwargs["model_label"]
        return spectrum

    def fake_builder(_fit, **kwargs):
        captured["ew_model"] = kwargs["ew_model"]
        assert kwargs["spectrum"] is spectrum
        return harness.point_builder.make_point(raw=kwargs["raw"], rs_ew_spectrum=spectrum)

    monkeypatch.setattr(harness.registry, "discover", lambda: None)
    monkeypatch.setattr(harness.registry, "import_failures", lambda: {})
    monkeypatch.setattr(harness.registry, "all_constraints", lambda: {})
    monkeypatch.setattr(harness.registry, "evaluate_all", lambda _point: {})
    monkeypatch.setattr(harness.RSEWSpectrum, "build", fake_spectrum_build)
    monkeypatch.setattr(harness.RSEWOverlapSplineCache, "build", lambda *a, **k: _Cache())
    monkeypatch.setattr(harness, "fit_quark_sector", lambda *a, **k: _Solution(result=fit))
    monkeypatch.setattr(harness, "f_IR", lambda *a, **k: 1.0)
    monkeypatch.setattr(
        harness,
        "compute_all_yukawas",
        lambda *a, **k: SimpleNamespace(Y_N_bar=np.ones(3)),
    )
    monkeypatch.setattr(harness, "_force_degenerate_neutrino_yukawas", lambda *a, **k: None)
    monkeypatch.setattr(harness.point_builder, "build_from_rs_ew_inputs", fake_builder)
    monkeypatch.setattr(harness, "compute_quark_kk_gluon_couplings", lambda *a, **k: object())
    cfg = harness.ScanConfig(
        mkk_values_gev=(3000.0,),
        n_draws_per_tile=1,
        ew_model=harness.CUSTODIAL_RS_PLR_EW_MODEL,
        expected_registry_count=0,
    )

    result = harness.run_universal_c_sanity(cfg, provenance={}, config_hash="cust")

    assert captured == {
        "model_label": harness.CUSTODIAL_RS_PLR_EW_MODEL,
        "ew_model": harness.CUSTODIAL_RS_PLR_EW_MODEL,
    }
    assert result["passes_no_spurious_hard_exclusions"] is True
    assert result["ew_model"] == harness.CUSTODIAL_RS_PLR_EW_MODEL


def test_quark_only_constraint_tallies_count_evaluated_active_failed_and_vetoed():
    harness = _load_harness()
    tallies = harness._empty_constraint_tallies(("K001", "EW003"))
    counters = {}
    from collections import Counter

    row = {
        "skipped": False,
        "survives_all_HARD_strict": False,
        "survives_all_HARD_inclusive": False,
        "excluded_by_rigorous": ["K001"],
        "excluded_by_proxy": [],
        "hard_not_evaluated": [],
        "constraints": {
            "K001": {
                "passes": False,
                "severity": "HARD",
                "ratio": 2.0,
                "active": True,
                "evaluated": True,
                "tag": "rigorous",
            },
            "EW003": {
                "passes": False,
                "severity": "SOFT",
                "ratio": 1.2,
                "active": True,
                "evaluated": True,
                "tag": "partial",
            },
        },
    }
    harness._accumulate_row(
        row,
        counters=Counter(counters),
        hard_veto_rigorous=Counter(),
        hard_veto_proxy=Counter(),
        hard_not_evaluated=Counter(),
        tag_counts=Counter(),
        exception_ids=Counter(),
        constraint_tallies=tallies,
    )
    final = harness._finalize_constraint_tallies(tallies)

    assert final["K001"]["evaluated"] == 1
    assert final["K001"]["active"] == 1
    assert final["K001"]["failed"] == 1
    assert final["K001"]["vetoed"] == 1
    assert final["EW003"]["failed"] == 1
    assert final["EW003"]["vetoed"] == 0


def test_deltaf2_bucket1_ratios_relax_with_high_mkk():
    from flavor_catalog_constraints import point_builder, registry
    from quarkConstraints.benchmarks import default_quark_targets, default_spurion_seed
    from quarkConstraints.couplings import compute_quark_kk_gluon_couplings
    from quarkConstraints.fit import fit_quark_sector

    fit = fit_quark_sector(
        default_quark_targets(),
        seed=default_spurion_seed(),
        overall_scale=3.0,
        max_nfev=120,
    ).result
    low = compute_quark_kk_gluon_couplings(fit, M_KK=1500.0, xi_KK=1.0, g_s_star=None)
    high = compute_quark_kk_gluon_couplings(fit, M_KK=15000.0, xi_KK=1.0, g_s_star=None)
    low_point = point_builder.make_point(
        quark_mass_basis_couplings=low,
        kk_gluon_mass_gev=1500.0,
    )
    high_point = point_builder.make_point(
        quark_mass_basis_couplings=high,
        kk_gluon_mass_gev=15000.0,
    )

    k001_low = registry.get("K001").evaluate(low_point)
    k001_high = registry.get("K001").evaluate(high_point)
    k002_low = registry.get("K002").evaluate(low_point)
    k002_high = registry.get("K002").evaluate(high_point)

    assert k001_low.ratio > k001_high.ratio
    assert k002_low.ratio > k002_high.ratio
    assert k001_low.passes is False
    assert k001_high.passes is True


def _candidate_get_extra_usage(harness):
    constants_by_id = {}
    usage = {}
    for pid in harness.QUARK_ONLY_CANDIDATE_IDS:
        paths = list((REPO_ROOT / "flavor_catalog_constraints").glob(f"**/{pid}.py"))
        assert len(paths) == 1
        tree = ast.parse(paths[0].read_text(encoding="utf-8"))
        constants = {}
        for node in tree.body:
            if isinstance(node, ast.Assign):
                for target in node.targets:
                    if isinstance(target, ast.Name):
                        constants[target.id] = _literal_or_none(node.value)
        constants_by_id[pid] = constants
        extras: list[str] = []

        class _Visitor(ast.NodeVisitor):
            def __init__(self):
                self.env = {}

            def visit_For(self, node):
                if isinstance(node.target, ast.Name) and isinstance(node.iter, ast.Name):
                    values = constants.get(node.iter.id)
                    if isinstance(values, tuple):
                        old = self.env.get(node.target.id, None)
                        had_old = node.target.id in self.env
                        for value in values:
                            self.env[node.target.id] = (str(value),)
                            for child in node.body:
                                self.visit(child)
                        if had_old:
                            self.env[node.target.id] = old
                        else:
                            self.env.pop(node.target.id, None)
                        for child in node.orelse:
                            self.visit(child)
                        return
                self.generic_visit(node)

            def visit_Call(self, node):
                if isinstance(node.func, ast.Attribute) and node.func.attr == "get_extra":
                    if node.args:
                        extras.extend(_resolve_extra_arg(node.args[0], constants, self.env))
                self.generic_visit(node)

        _Visitor().visit(tree)
        usage[pid] = tuple(dict.fromkeys(extras))
    return usage


def _literal_or_none(node):
    try:
        return ast.literal_eval(node)
    except (ValueError, SyntaxError):
        return None


def _resolve_extra_arg(node, constants, env):
    if isinstance(node, ast.Constant):
        return [str(node.value)]
    if isinstance(node, ast.Name):
        if node.id in env:
            return list(env[node.id])
        value = constants.get(node.id)
        if isinstance(value, str):
            return [value]
        if isinstance(value, tuple):
            return [str(item) for item in value]
        return [node.id]
    return [ast.unparse(node)]
