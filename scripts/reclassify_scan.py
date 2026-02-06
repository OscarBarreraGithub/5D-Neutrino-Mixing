#!/usr/bin/env python3
"""Reclassify existing scan CSV outputs without rerunning the physics scan."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from scanParams import AnarchyConfig, ReclassifyConfig, classify_row


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Recompute pass/fail categorization from an existing scan CSV. "
            "This lets you change selection policy (including anarchic criteria) "
            "without rerunning compute_all_yukawas over the full grid."
        )
    )
    parser.add_argument("input_csv", help="Path to original scan CSV")
    parser.add_argument("output_csv", help="Path to write reclassified CSV")

    parser.add_argument("--max-y-bar", type=float, default=4.0, help="Perturbativity bound on max |Y_bar|")
    parser.add_argument(
        "--natural-min",
        type=float,
        default=0.1,
        help="Lower bound for naturalness window on |Y_bar|",
    )
    parser.add_argument(
        "--natural-max",
        type=float,
        default=4.0,
        help="Upper bound for naturalness window on |Y_bar|",
    )
    parser.add_argument(
        "--skip-lfv",
        action="store_true",
        help="Do not include LFV pass/fail in reclassification",
    )

    parser.add_argument(
        "--compute-anarchy",
        action="store_true",
        help="Recompute anarchy score from stored Y_N_bar_i + PMNS metadata",
    )
    parser.add_argument(
        "--anarchy-min-score",
        type=float,
        default=None,
        help="Optional minimum anarchy score cut (requires --compute-anarchy)",
    )
    parser.add_argument("--anarchy-magnitude-min", type=float, default=1.0 / 3.0)
    parser.add_argument("--anarchy-magnitude-max", type=float, default=3.0)
    parser.add_argument("--anarchy-phase-min", type=float, default=0.0)
    parser.add_argument("--anarchy-phase-max", type=float, default=2.0 * 3.141592653589793)
    parser.add_argument("--anarchy-yn-overall-min", type=float, default=0.01)
    parser.add_argument("--anarchy-yn-overall-max", type=float, default=0.2)
    parser.add_argument("--anarchy-w-band", type=float, default=1.0)
    parser.add_argument("--anarchy-w-cond", type=float, default=1.0)
    parser.add_argument("--anarchy-w-fit", type=float, default=0.0)
    return parser


def main() -> int:
    args = _build_parser().parse_args()

    if args.anarchy_min_score is not None and not args.compute_anarchy:
        raise ValueError("--anarchy-min-score requires --compute-anarchy")

    anarchy_config = None
    if args.compute_anarchy:
        anarchy_config = AnarchyConfig(
            magnitude_min=args.anarchy_magnitude_min,
            magnitude_max=args.anarchy_magnitude_max,
            phase_min=args.anarchy_phase_min,
            phase_max=args.anarchy_phase_max,
            yN_overall_min=args.anarchy_yn_overall_min,
            yN_overall_max=args.anarchy_yn_overall_max,
            w_band=args.anarchy_w_band,
            w_cond=args.anarchy_w_cond,
            w_fit=args.anarchy_w_fit,
        )

    reclass_config = ReclassifyConfig(
        max_Y_bar=args.max_y_bar,
        naturalness_range=(args.natural_min, args.natural_max),
        require_lfv=not args.skip_lfv,
        anarchy=anarchy_config,
        anarchy_min_score=args.anarchy_min_score,
    )

    input_path = Path(args.input_csv)
    output_path = Path(args.output_csv)
    if not input_path.exists():
        raise FileNotFoundError(f"Input CSV not found: {input_path}")

    with input_path.open("r", encoding="utf-8", newline="") as in_handle:
        reader = csv.DictReader(in_handle)
        if reader.fieldnames is None:
            raise ValueError(f"Input CSV has no header: {input_path}")
        in_fields = list(reader.fieldnames)
        rows = list(reader)

    reclass_columns = [
        "reclass_max_Y_bar_observed",
        "reclass_min_Y_bar_observed",
        "reclass_perturbative",
        "reclass_natural",
        "reclass_lfv_passes",
        "reclass_anarchy_score",
        "reclass_anarchy_band_penalty",
        "reclass_anarchy_condition_penalty",
        "reclass_anarchy_yN_overall",
        "reclass_passes_all",
        "reclass_reject_reason",
    ]
    out_fields = in_fields + [c for c in reclass_columns if c not in in_fields]

    n_pass = 0
    out_rows = []
    for row in rows:
        rec = classify_row(row, reclass_config)
        merged = dict(row)
        merged.update(rec)
        out_rows.append(merged)
        if rec["reclass_passes_all"]:
            n_pass += 1

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as out_handle:
        writer = csv.DictWriter(out_handle, fieldnames=out_fields)
        writer.writeheader()
        writer.writerows(out_rows)

    print(f"Wrote {len(out_rows)} rows to {output_path}")
    print(f"Reclassified acceptance: {n_pass}/{len(out_rows)} ({100.0 * n_pass / max(len(out_rows), 1):.1f}%)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
