from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = REPO_ROOT / "scripts" / "analyze_wq_quarkonly.py"
PLAN_PATH = REPO_ROOT / "scripts" / "wq_quarkonly_1m_plan.py"
SMOKE_DIR = REPO_ROOT / ".orchestration" / "runs" / "WQ-QUARKONLY" / "smoke-1k"


def _load_analysis():
    spec = importlib.util.spec_from_file_location("analyze_wq_quarkonly", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _load_plan():
    spec = importlib.util.spec_from_file_location("wq_quarkonly_1m_plan", PLAN_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_wq_quarkonly_analysis_emits_smoke_plots_and_report(tmp_path):
    if not SMOKE_DIR.is_dir():
        pytest.skip(f"missing smoke output fixture: {SMOKE_DIR}")
    module = _load_analysis()

    analysis = module.analyze(SMOKE_DIR, max_samples_per_group=1000)
    plot_paths = module.write_outputs(analysis, tmp_path)
    report = tmp_path / "analysis_report.md"

    assert report.is_file()
    assert report.stat().st_size > 0
    assert "grouped by exact (quark_fit_r, M_KK)" in report.read_text(encoding="utf-8")
    assert {path.name for path in plot_paths} == set(module.EXPECTED_PLOTS)
    for name in module.EXPECTED_PLOTS:
        path = tmp_path / "plots" / name
        assert path.is_file()
        assert path.stat().st_size > 0


def test_wq_quarkonly_1m_seed_intervals_are_disjoint():
    plan = _load_plan()

    intervals = []
    for task_id in range(plan.total_tasks()):
        task = plan.task_plan(task_id)
        for tile_id in range(task.n_mkk_tiles):
            start = task.base_seed + task.tile_seed_stride * tile_id
            stop = start + task.draws_per_mkk
            intervals.append((start, stop, task_id, tile_id))

    intervals.sort()
    for left, right in zip(intervals, intervals[1:]):
        assert left[1] <= right[0], (left, right)

    proof = plan.summary_payload()["seed_disjointness"]
    assert proof["passes"] is True
    assert plan.SHARD_SEED_BLOCK > (
        plan.TILE_SEED_STRIDE * len(plan.M_KK_TEV)
        + plan.DRAWS_PER_MKK_PER_SHARD
    )
