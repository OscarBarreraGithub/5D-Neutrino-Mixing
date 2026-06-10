from __future__ import annotations

import csv
import importlib.util
import json
import sys
from pathlib import Path

import pyarrow.parquet as pq


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = REPO_ROOT / "scripts" / "build_wq_quarkonly_comparison.py"


def _load_comparison():
    spec = importlib.util.spec_from_file_location("build_wq_quarkonly_comparison", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_wq_quarkonly_comparison_pairs_seed_normalizes_survival_and_veto_enums(tmp_path):
    module = _load_comparison()
    minimal = tmp_path / "wq_quarkonly_1M_20128400"
    custodial = tmp_path / "wq_quarkonly_1M_custodial_999"
    output = custodial / "comparison"
    _write_run(
        minimal,
        [
            _row(
                seed=101,
                strict=False,
                inclusive=False,
                excluded_by_rigorous=["K001"],
                hard_not_evaluated=["T001"],
                constraints={
                    "K001": _constraint(False, "rigorous", ratio=2.0),
                    "T001": _constraint(False, "stub", evaluated=False, active=False),
                },
            ),
            _row(seed=102, skipped=True, skip_reason="quark_fit_failed"),
        ],
        config_hash="minhash",
        git_sha="minsha",
    )
    _write_run(
        custodial,
        [
            _row(
                seed=101,
                strict=False,
                inclusive=False,
                excluded_by_proxy=["B011"],
                constraints={"B011": _constraint(False, "proxy", ratio=3.0)},
            ),
            _row(seed=102, strict=True, inclusive=True),
        ],
        config_hash="custhash",
        git_sha="custsha",
    )

    manifest = module.build_comparison(minimal, custodial, output)

    assert manifest["validation"]["paired_rows"] == 2
    assert manifest["pairing_key"] == ["r", "mkk_tev", "draw_seed"]
    assert manifest["physics_tags"]["custodial_top_partner_loop_status"] == "deferred"
    paired = pq.read_table(output / "paired_draws.parquet").to_pylist()
    assert [(row["r"], row["mkk_tev"], row["draw_seed"]) for row in paired] == [
        (0.05, 2.0, 101),
        (0.05, 2.0, 102),
    ]
    assert paired[0]["minimal_excluded_by_rigorous"] == ["K001"]
    assert paired[1]["minimal_skipped"] is True
    assert paired[1]["minimal_skip_reason"] == "quark_fit_failed"

    vetoes = pq.read_table(output / "paired_vetoes.parquet").to_pylist()
    veto_keys = {
        (row["ew_model"], row["constraint_id"], row["veto_class"])
        for row in vetoes
    }
    assert (
        module.MINIMAL_RS_EW_MODEL,
        "K001",
        "rigorous",
    ) in veto_keys
    assert (
        module.MINIMAL_RS_EW_MODEL,
        "T001",
        "not_evaluated",
    ) in veto_keys
    assert (
        module.CUSTODIAL_RS_PLR_EW_MODEL,
        "B011",
        "proxy",
    ) in veto_keys

    with (output / "survival_by_r_mkk.csv").open(encoding="utf-8") as fh:
        survival = list(csv.DictReader(fh))
    assert len(survival) == 1
    assert survival[0]["minimal_rows"] == "2"
    assert survival[0]["minimal_evaluated"] == "1"
    assert survival[0]["minimal_skipped"] == "1"
    assert survival[0]["minimal_strict_survival_frac"] == "0.0"
    assert survival[0]["custodial_strict_survival_frac"] == "0.5"

    readme = (output / "README.md").read_text(encoding="utf-8")
    assert "Raw scan rows use `seed`, not `draw_seed`" in readme
    assert (output / "manifest.json").is_file()
    assert (output / "schema.json").is_file()
    assert (output / "run_index.json").is_file()


def _write_run(root: Path, rows, *, config_hash: str, git_sha: str) -> None:
    scan_plan = {
        "r_grid": [0.05],
        "mkk_tev": [2.0],
        "base_seed": 100,
        "tile_seed_stride": 1000003,
        "shard_seed_block": 20000000,
        "seed_disjointness": {
            "draw_seed_formula": "tile_seed + draw_idx",
            "tile_seed_formula": "base_seed + tile_seed_stride * tile_id",
            "per_shard_base_formula": "BASE_SEED + SHARD_SEED_BLOCK * task_id",
        },
    }
    shard = root / "r0p05" / "shard-00"
    shard.mkdir(parents=True)
    (root / "scan_plan.json").write_text(json.dumps(scan_plan, indent=2) + "\n", encoding="utf-8")
    path = shard / "tile-00000.jsonl"
    with path.open("w", encoding="utf-8") as fh:
        for draw_id, row in enumerate(rows):
            payload = dict(row)
            payload["draw_id"] = draw_id
            payload["provenance"] = {
                "config_hash": config_hash,
                "git_sha": git_sha,
            }
            fh.write(json.dumps(payload, sort_keys=True) + "\n")
    summary = {
        "complete": True,
        "config_hash": config_hash,
        "provenance": {"git_sha": git_sha},
    }
    (shard / "tile-00000.summary.json").write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )


def _row(
    *,
    seed: int,
    strict: bool = True,
    inclusive: bool = True,
    skipped: bool = False,
    skip_reason: str | None = None,
    excluded_by_rigorous=None,
    excluded_by_proxy=None,
    hard_not_evaluated=None,
    constraints=None,
):
    return {
        "tile_id": 0,
        "seed": seed,
        "skipped": skipped,
        "skip_reason": skip_reason,
        "quark_fit_r": 0.05,
        "params": {
            "M_KK": 2000.0,
            "quark_fit_r": 0.05,
        },
        "survives_all_HARD_strict": strict,
        "survives_all_HARD_inclusive": inclusive,
        "excluded_by_rigorous": list(excluded_by_rigorous or []),
        "excluded_by_proxy": list(excluded_by_proxy or []),
        "hard_not_evaluated": list(hard_not_evaluated or []),
        "constraints": dict(constraints or {}),
    }


def _constraint(
    passes: bool,
    tag: str,
    *,
    ratio: float | None = None,
    severity: str = "HARD",
    evaluated: bool = True,
    active: bool = True,
):
    return {
        "passes": passes,
        "tag": tag,
        "ratio": ratio,
        "severity": severity,
        "evaluated": evaluated,
        "active": active,
    }
