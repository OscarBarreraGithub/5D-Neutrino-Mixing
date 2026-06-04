from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

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
