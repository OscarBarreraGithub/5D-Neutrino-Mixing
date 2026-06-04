#!/usr/bin/env python3
"""Plan helper for the WQ quark-only 1M r x M_KK scan."""

from __future__ import annotations

import argparse
import json
import shlex
from dataclasses import asdict, dataclass
from typing import Sequence


R_GRID = (0.05, 0.10, 0.25, 0.50, 1.00)
M_KK_TEV = (1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0, 30.0, 50.0)
DRAWS_PER_MKK_PER_R = 20_000
DRAW_SHARDS_PER_R = 10
DRAWS_PER_MKK_PER_SHARD = DRAWS_PER_MKK_PER_R // DRAW_SHARDS_PER_R
BASE_SEED = 202_606_040_000
TILE_SEED_STRIDE = 1_000_003
SHARD_SEED_BLOCK = 20_000_000


@dataclass(frozen=True)
class TaskPlan:
    task_id: int
    r_index: int
    r_value: float
    r_label: str
    draw_shard: int
    draw_start: int
    draw_stop: int
    draws_per_mkk: int
    base_seed: int
    tile_seed_stride: int
    shard_seed_block: int
    mkk_tev_csv: str
    n_mkk_tiles: int


def total_tasks() -> int:
    return len(R_GRID) * DRAW_SHARDS_PER_R


def r_label(value: float) -> str:
    text = f"{value:.2f}".rstrip("0").rstrip(".")
    return "r" + text.replace(".", "p")


def task_plan(task_id: int) -> TaskPlan:
    if task_id < 0 or task_id >= total_tasks():
        raise ValueError(f"task_id must be in [0,{total_tasks() - 1}]")
    r_index = task_id // DRAW_SHARDS_PER_R
    draw_shard = task_id % DRAW_SHARDS_PER_R
    draw_start = draw_shard * DRAWS_PER_MKK_PER_SHARD
    draw_stop = draw_start + DRAWS_PER_MKK_PER_SHARD
    r_value = float(R_GRID[r_index])
    return TaskPlan(
        task_id=task_id,
        r_index=r_index,
        r_value=r_value,
        r_label=r_label(r_value),
        draw_shard=draw_shard,
        draw_start=draw_start,
        draw_stop=draw_stop,
        draws_per_mkk=DRAWS_PER_MKK_PER_SHARD,
        base_seed=BASE_SEED + SHARD_SEED_BLOCK * task_id,
        tile_seed_stride=TILE_SEED_STRIDE,
        shard_seed_block=SHARD_SEED_BLOCK,
        mkk_tev_csv=",".join(_format_float(x) for x in M_KK_TEV),
        n_mkk_tiles=len(M_KK_TEV),
    )


def summary_payload() -> dict[str, object]:
    n_mkk = len(M_KK_TEV)
    n_r = len(R_GRID)
    total_draws = n_mkk * n_r * DRAWS_PER_MKK_PER_R
    max_tile_span = TILE_SEED_STRIDE * n_mkk + DRAWS_PER_MKK_PER_SHARD
    return {
        "r_grid": list(R_GRID),
        "mkk_tev": list(M_KK_TEV),
        "draws_per_mkk_per_r": DRAWS_PER_MKK_PER_R,
        "draw_shards_per_r": DRAW_SHARDS_PER_R,
        "draws_per_mkk_per_shard": DRAWS_PER_MKK_PER_SHARD,
        "total_tasks": total_tasks(),
        "total_draws": total_draws,
        "base_seed": BASE_SEED,
        "tile_seed_stride": TILE_SEED_STRIDE,
        "shard_seed_block": SHARD_SEED_BLOCK,
        "seed_disjointness": {
            "tile_seed_formula": "base_seed + tile_seed_stride * tile_id",
            "draw_seed_formula": "tile_seed + draw_idx",
            "per_shard_base_formula": "BASE_SEED + SHARD_SEED_BLOCK * task_id",
            "required_block_gt": "tile_seed_stride * n_mkk_tiles + draws_per_mkk_per_shard",
            "required_block_value": max_tile_span,
            "chosen_block_value": SHARD_SEED_BLOCK,
            "passes": SHARD_SEED_BLOCK > max_tile_span,
        },
        "physics_note": (
            "r weights the up-type MFV spurion in C_Q = r Yu Yu^dagger + "
            "Yd Yd^dagger; the grid spans down-dominated, default-near, and "
            "up-dominated LH-doublet localization."
        ),
    }


def seed_proof_text() -> str:
    payload = summary_payload()
    proof = payload["seed_disjointness"]
    return "\n".join(
        [
            "# WQ Quark-Only 1M Seed Disjointness",
            "",
            f"- tasks: {payload['total_tasks']}",
            f"- M_KK tiles/task: {len(M_KK_TEV)}",
            f"- draws/M_KK/task: {DRAWS_PER_MKK_PER_SHARD}",
            f"- tile seed stride: {TILE_SEED_STRIDE}",
            f"- shard seed block: {SHARD_SEED_BLOCK}",
            f"- per-shard base: {BASE_SEED} + {SHARD_SEED_BLOCK} * task_id",
            (
                "- tile seed: base + "
                f"{TILE_SEED_STRIDE} * tile_id, for tile_id=0..{len(M_KK_TEV) - 1}"
            ),
            "- draw seed: tile_seed + draw_idx",
            (
                "- required block > stride*n_tiles + n_draws = "
                f"{proof['required_block_value']}; chosen block = {SHARD_SEED_BLOCK}"
            ),
            "- conclusion: draw-seed intervals are disjoint within and across all shards.",
            "",
        ]
    )


def emit_task_env(plan: TaskPlan) -> str:
    values = {
        "WQ_TASK_ID": plan.task_id,
        "WQ_R_INDEX": plan.r_index,
        "WQ_R_VALUE": _format_float(plan.r_value),
        "WQ_R_LABEL": plan.r_label,
        "WQ_DRAW_SHARD": plan.draw_shard,
        "WQ_DRAW_SHARD_PADDED": f"{plan.draw_shard:02d}",
        "WQ_DRAW_START": plan.draw_start,
        "WQ_DRAW_STOP": plan.draw_stop,
        "WQ_DRAWS_PER_MKK_PER_SHARD": plan.draws_per_mkk,
        "WQ_BASE_SEED": plan.base_seed,
        "WQ_TILE_SEED_STRIDE": plan.tile_seed_stride,
        "WQ_SHARD_SEED_BLOCK": plan.shard_seed_block,
        "WQ_M_KK_TEV": plan.mkk_tev_csv,
        "WQ_N_MKK_TILES": plan.n_mkk_tiles,
        "WQ_TOTAL_TASKS": total_tasks(),
    }
    return "\n".join(
        f"export {key}={shlex.quote(str(value))}"
        for key, value in values.items()
    )


def _format_float(value: float) -> str:
    return f"{float(value):g}"


def _build_argparser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("summary")
    sub.add_parser("seed-proof")
    task_env = sub.add_parser("task-env")
    task_env.add_argument("--task-id", type=int, required=True)
    task_json = sub.add_parser("task-json")
    task_json.add_argument("--task-id", type=int, required=True)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _build_argparser().parse_args(argv)
    if args.command == "summary":
        print(json.dumps(summary_payload(), indent=2, sort_keys=True))
        return 0
    if args.command == "seed-proof":
        print(seed_proof_text(), end="")
        return 0
    if args.command == "task-env":
        print(emit_task_env(task_plan(int(args.task_id))))
        return 0
    if args.command == "task-json":
        print(json.dumps(asdict(task_plan(int(args.task_id))), indent=2, sort_keys=True))
        return 0
    raise AssertionError(args.command)


if __name__ == "__main__":
    raise SystemExit(main())
