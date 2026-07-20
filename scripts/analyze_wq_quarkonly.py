#!/usr/bin/env python3
"""Analyze WQ quark-only scan JSONL shards without mixing r or M_KK groups."""

from __future__ import annotations

import argparse
import json
import math
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

EXPECTED_PLOTS = (
    "survival_vs_mkk.png",
    "constraint_veto_fraction_rigorous.png",
    "constraint_veto_fraction_proxy.png",
    "constraining_power_ranking_rigorous.png",
    "constraining_power_ranking_proxy.png",
    "yukawa_singular_values.png",
    "bulk_c_distributions.png",
    "max_abs_quark_yukawa.png",
    "per_r_survival_by_mkk.png",
    "per_r_cq_localization_by_mkk.png",
    "per_r_yukawa_singular_values_by_mkk.png",
)


@dataclass
class GroupStats:
    rows: int = 0
    skipped: int = 0
    evaluated: int = 0
    survives_strict: int = 0
    survives_inclusive: int = 0


@dataclass
class ConstraintStats:
    points: int = 0
    evaluated: int = 0
    active: int = 0
    failed: int = 0
    vetoed: int = 0
    severity_counts: Counter[str] = field(default_factory=Counter)


@dataclass(frozen=True)
class SampleRecord:
    up_singular: tuple[float, ...]
    down_singular: tuple[float, ...]
    c_q: tuple[float, ...]
    c_u: tuple[float, ...]
    c_d: tuple[float, ...]
    max_abs_y: float | None
    singular_source: str


@dataclass
class Reservoir:
    limit: int
    seen: int = 0
    records: list[SampleRecord] = field(default_factory=list)

    def add(self, record: SampleRecord, rng: np.random.Generator) -> None:
        self.seen += 1
        if len(self.records) < self.limit:
            self.records.append(record)
            return
        index = int(rng.integers(0, self.seen))
        if index < self.limit:
            self.records[index] = record


@dataclass
class Analysis:
    input_paths: list[Path]
    group_stats: dict[tuple[float, float], GroupStats] = field(default_factory=dict)
    constraint_stats: dict[tuple[float, float, str, str], ConstraintStats] = field(
        default_factory=dict
    )
    samples: dict[tuple[float, float], Reservoir] = field(default_factory=dict)
    singular_source_counts: Counter[str] = field(default_factory=Counter)
    malformed_rows: int = 0
    default_r: float | None = None


def analyze(input_path: Path, *, max_samples_per_group: int = 5000) -> Analysis:
    jsonl_paths = _discover_jsonl(input_path)
    if not jsonl_paths:
        raise FileNotFoundError(f"no JSONL files found under {input_path}")
    default_r = _discover_default_r(input_path)
    analysis = Analysis(input_paths=jsonl_paths, default_r=default_r)
    rng = np.random.default_rng(20260604)

    for path in jsonl_paths:
        with path.open("r", encoding="utf-8") as fh:
            for line_number, line in enumerate(fh, start=1):
                if not line.strip():
                    continue
                try:
                    row = json.loads(line)
                    key = _row_group_key(row, default_r=default_r)
                except Exception as exc:  # noqa: BLE001 - report and continue analysis
                    analysis.malformed_rows += 1
                    print(f"[analyze_wq] skipped malformed {path}:{line_number}: {exc}", file=sys.stderr)
                    continue
                stats = analysis.group_stats.setdefault(key, GroupStats())
                stats.rows += 1
                if bool(row.get("skipped")):
                    stats.skipped += 1
                    continue
                stats.evaluated += 1
                if bool(row.get("survives_all_HARD_strict")):
                    stats.survives_strict += 1
                if bool(row.get("survives_all_HARD_inclusive")):
                    stats.survives_inclusive += 1
                _accumulate_constraints(analysis, key, row)
                record = _sample_from_row(row)
                analysis.singular_source_counts[record.singular_source] += 1
                reservoir = analysis.samples.setdefault(
                    key,
                    Reservoir(limit=max(0, int(max_samples_per_group))),
                )
                if reservoir.limit > 0:
                    reservoir.add(record, rng)
    return analysis


def write_outputs(analysis: Analysis, output_root: Path) -> list[Path]:
    plots_dir = output_root / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    output_root.mkdir(parents=True, exist_ok=True)
    plot_paths = [
        _plot_survival_vs_mkk(analysis, plots_dir / "survival_vs_mkk.png"),
        _plot_constraint_veto_by_mkk(
            analysis,
            tag="rigorous",
            path=plots_dir / "constraint_veto_fraction_rigorous.png",
        ),
        _plot_constraint_veto_by_mkk(
            analysis,
            tag="proxy",
            path=plots_dir / "constraint_veto_fraction_proxy.png",
        ),
        _plot_constraining_power_heatmap(
            analysis,
            tag="rigorous",
            path=plots_dir / "constraining_power_ranking_rigorous.png",
        ),
        _plot_constraining_power_heatmap(
            analysis,
            tag="proxy",
            path=plots_dir / "constraining_power_ranking_proxy.png",
        ),
        _plot_yukawa_singular_values(analysis, plots_dir / "yukawa_singular_values.png"),
        _plot_bulk_c_distributions(analysis, plots_dir / "bulk_c_distributions.png"),
        _plot_max_abs_yukawa(analysis, plots_dir / "max_abs_quark_yukawa.png"),
        _plot_per_r_survival(analysis, plots_dir / "per_r_survival_by_mkk.png"),
        _plot_per_r_cq(analysis, plots_dir / "per_r_cq_localization_by_mkk.png"),
        _plot_per_r_yukawa(
            analysis,
            plots_dir / "per_r_yukawa_singular_values_by_mkk.png",
        ),
    ]
    _write_report(analysis, output_root / "analysis_report.md", plot_paths)
    return plot_paths


def _discover_jsonl(input_path: Path) -> list[Path]:
    if input_path.is_file():
        return [input_path] if input_path.suffix == ".jsonl" else []
    paths = [
        path
        for path in input_path.rglob("*.jsonl")
        if path.is_file() and ".tmp." not in path.name
    ]
    return sorted(paths)


def _discover_default_r(input_path: Path) -> float | None:
    roots = [input_path] if input_path.is_dir() else [input_path.parent]
    candidates: set[float] = set()
    for root in roots:
        for summary_path in root.rglob("run_summary.json"):
            try:
                summary = json.loads(summary_path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                continue
            value = _finite_or_none(
                dict(summary.get("config", {})).get("quark_fit_r")
            )
            if value is not None:
                candidates.add(float(value))
    if len(candidates) == 1:
        return next(iter(candidates))
    return None


def _row_group_key(row: Mapping[str, Any], *, default_r: float | None) -> tuple[float, float]:
    params = dict(row.get("params") or {})
    provenance = dict(row.get("provenance") or {})
    r_value = _finite_or_none(
        row.get("quark_fit_r", params.get("quark_fit_r", provenance.get("quark_fit_r")))
    )
    if r_value is None:
        r_value = default_r
    if r_value is None:
        raise ValueError("row has no quark_fit_r and no unique run_summary fallback")
    mkk = _finite_or_none(params.get("M_KK"))
    if mkk is None:
        raise ValueError("row has no params.M_KK")
    return (float(r_value), float(mkk))


def _accumulate_constraints(
    analysis: Analysis,
    key: tuple[float, float],
    row: Mapping[str, Any],
) -> None:
    for pid, item_raw in dict(row.get("constraints") or {}).items():
        item = dict(item_raw)
        tag = str(item.get("tag", "unknown"))
        stats_key = (key[0], key[1], str(pid), tag)
        stats = analysis.constraint_stats.setdefault(stats_key, ConstraintStats())
        stats.points += 1
        if bool(item.get("evaluated")):
            stats.evaluated += 1
        if bool(item.get("active")):
            stats.active += 1
        severity = str(item.get("severity", "unknown"))
        stats.severity_counts[severity] += 1
        failed = not bool(item.get("passes", True))
        if failed:
            stats.failed += 1
        if (
            failed
            and severity == "HARD"
            and bool(item.get("evaluated"))
            and tag in {"rigorous", "proxy"}
        ):
            stats.vetoed += 1


def _sample_from_row(row: Mapping[str, Any]) -> SampleRecord:
    fit = dict(row.get("fit_diagnostics") or {})
    up = _float_tuple(fit.get("fitted_up_yukawa_singular_values"))
    down = _float_tuple(fit.get("fitted_down_yukawa_singular_values"))
    source = "fitted"
    if not up or not down:
        up, down = _seed_singular_values(row)
        source = "seed_fallback" if up and down else "missing"
    return SampleRecord(
        up_singular=up,
        down_singular=down,
        c_q=_float_tuple(fit.get("bulk_c_Q")),
        c_u=_float_tuple(fit.get("bulk_c_u")),
        c_d=_float_tuple(fit.get("bulk_c_d")),
        max_abs_y=_finite_or_none(fit.get("max_abs_quark_yukawa")),
        singular_source=source,
    )


def _seed_singular_values(row: Mapping[str, Any]) -> tuple[tuple[float, ...], tuple[float, ...]]:
    seed = dict(dict(row.get("params") or {}).get("quark_yukawa_seed") or {})
    try:
        y_u = np.asarray(seed["Y_u_re"], dtype=float) + 1j * np.asarray(
            seed["Y_u_im"], dtype=float
        )
        y_d = np.asarray(seed["Y_d_re"], dtype=float) + 1j * np.asarray(
            seed["Y_d_im"], dtype=float
        )
    except (KeyError, TypeError, ValueError):
        return (), ()
    return (_singular_values(y_u), _singular_values(y_d))


def _singular_values(matrix: np.ndarray) -> tuple[float, ...]:
    values = np.linalg.svd(np.asarray(matrix, dtype=np.complex128), compute_uv=False)
    return tuple(float(x) for x in np.sort(np.asarray(values, dtype=float)))


def _plot_survival_vs_mkk(analysis: Analysis, path: Path) -> Path:
    fig, ax = plt.subplots(figsize=(8.8, 5.4))
    r_values = _r_values(analysis)
    if not r_values:
        _blank_axis(ax, "No grouped rows found")
    for r_value in r_values:
        mkk, strict, inclusive = [], [], []
        for mkk_value in _mkk_values(analysis):
            stats = analysis.group_stats.get((r_value, mkk_value))
            if not stats or stats.evaluated == 0:
                continue
            mkk.append(mkk_value / 1000.0)
            strict.append(stats.survives_strict / stats.evaluated)
            inclusive.append(stats.survives_inclusive / stats.evaluated)
        if mkk:
            ax.plot(mkk, strict, marker="o", label=f"r={r_value:g} strict")
            ax.plot(mkk, inclusive, marker="s", linestyle="--", label=f"r={r_value:g} incl.")
    ax.set_xscale("log")
    ax.set_xlabel("M_KK [TeV]")
    ax.set_ylabel("Survival fraction")
    ax.set_title("Survival vs M_KK, grouped separately by r")
    ax.set_ylim(-0.03, 1.03)
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8, ncols=2)
    _save(fig, path)
    return path


def _plot_constraint_veto_by_mkk(analysis: Analysis, *, tag: str, path: Path) -> Path:
    r_values = _r_values(analysis)
    fig, axes = _facet_axes(len(r_values), figsize=(12.0, max(3.2, 2.8 * math.ceil(max(len(r_values), 1) / 2))))
    top_pids = _top_constraint_ids(analysis, tag=tag, limit=10)
    if not r_values or not top_pids:
        _blank_axis(axes[0], f"No {tag} HARD vetoes found")
    for ax, r_value in zip(axes, r_values):
        for pid in top_pids:
            xs, ys = [], []
            for mkk in _mkk_values(analysis):
                stats = analysis.constraint_stats.get((r_value, mkk, pid, tag))
                if not stats or stats.evaluated == 0:
                    continue
                xs.append(mkk / 1000.0)
                ys.append(stats.vetoed / stats.evaluated)
            if xs:
                ax.plot(xs, ys, marker="o", linewidth=1.4, label=pid)
        ax.set_title(f"r={r_value:g}")
        ax.set_xscale("log")
        ax.set_ylim(-0.03, 1.03)
        ax.grid(True, alpha=0.25)
        ax.set_xlabel("M_KK [TeV]")
        ax.set_ylabel("Veto fraction")
    _hide_unused(axes, len(r_values))
    if top_pids:
        axes[0].legend(fontsize=7, ncols=2)
    fig.suptitle(f"{tag.capitalize()} per-constraint HARD veto fraction vs M_KK")
    _save(fig, path)
    return path


def _plot_constraining_power_heatmap(analysis: Analysis, *, tag: str, path: Path) -> Path:
    pids = _top_constraint_ids(analysis, tag=tag, limit=15)
    groups = _group_keys(analysis)
    fig, ax = plt.subplots(figsize=(max(8.0, 0.42 * max(len(groups), 1)), max(4.0, 0.32 * max(len(pids), 1))))
    if not pids or not groups:
        _blank_axis(ax, f"No {tag} ranking data found")
        _save(fig, path)
        return path
    matrix = np.full((len(pids), len(groups)), np.nan, dtype=float)
    for row_idx, pid in enumerate(pids):
        for col_idx, (r_value, mkk) in enumerate(groups):
            stats = analysis.constraint_stats.get((r_value, mkk, pid, tag))
            if stats and stats.evaluated:
                matrix[row_idx, col_idx] = stats.vetoed / stats.evaluated
    image = ax.imshow(matrix, aspect="auto", vmin=0.0, vmax=1.0, cmap="magma")
    ax.set_yticks(np.arange(len(pids)), labels=pids)
    ax.set_xticks(
        np.arange(len(groups)),
        labels=[f"r={r:g}\n{m / 1000:g}T" for r, m in groups],
        rotation=90,
    )
    ax.set_title(f"{tag.capitalize()} constraining power by exact (r, M_KK) group")
    fig.colorbar(image, ax=ax, label="HARD veto fraction")
    _save(fig, path)
    return path


def _plot_yukawa_singular_values(analysis: Analysis, path: Path) -> Path:
    mkk_values = _mkk_values(analysis)
    fig, axes = _facet_axes(len(mkk_values), figsize=(13.0, max(3.4, 2.7 * math.ceil(max(len(mkk_values), 1) / 2))))
    if not mkk_values:
        _blank_axis(axes[0], "No Yukawa singular-value samples found")
    for ax, mkk in zip(axes, mkk_values):
        _boxplot_by_r(
            ax,
            analysis,
            mkk,
            series=(
                ("up", lambda rec: rec.up_singular, "#386cb0", -0.16),
                ("down", lambda rec: rec.down_singular, "#fdb462", 0.16),
            ),
        )
        ax.set_title(f"M_KK={mkk / 1000:g} TeV")
        ax.set_ylabel("Singular value")
        ax.set_yscale("log")
        ax.grid(True, axis="y", alpha=0.25)
    _hide_unused(axes, len(mkk_values))
    fig.suptitle("Yukawa singular-value distributions by exact r and M_KK")
    _save(fig, path)
    return path


def _plot_bulk_c_distributions(analysis: Analysis, path: Path) -> Path:
    mkk_values = _mkk_values(analysis)
    fig, axes = _facet_axes(len(mkk_values), figsize=(13.0, max(3.4, 2.7 * math.ceil(max(len(mkk_values), 1) / 2))))
    if not mkk_values:
        _blank_axis(axes[0], "No bulk-c samples found")
    for ax, mkk in zip(axes, mkk_values):
        _boxplot_by_r(
            ax,
            analysis,
            mkk,
            series=(
                ("c_Q", lambda rec: rec.c_q, "#7fc97f", -0.22),
                ("c_u", lambda rec: rec.c_u, "#beaed4", 0.0),
                ("c_d", lambda rec: rec.c_d, "#fdc086", 0.22),
            ),
        )
        ax.set_title(f"M_KK={mkk / 1000:g} TeV")
        ax.set_ylabel("bulk c")
        ax.set_ylim(0.25, 0.95)
        ax.grid(True, axis="y", alpha=0.25)
    _hide_unused(axes, len(mkk_values))
    fig.suptitle("Bulk localization distributions by exact r and M_KK")
    _save(fig, path)
    return path


def _plot_max_abs_yukawa(analysis: Analysis, path: Path) -> Path:
    mkk_values = _mkk_values(analysis)
    fig, axes = _facet_axes(len(mkk_values), figsize=(13.0, max(3.4, 2.7 * math.ceil(max(len(mkk_values), 1) / 2))))
    if not mkk_values:
        _blank_axis(axes[0], "No perturbativity samples found")
    for ax, mkk in zip(axes, mkk_values):
        _boxplot_by_r(
            ax,
            analysis,
            mkk,
            series=(("|Yq|max", lambda rec: () if rec.max_abs_y is None else (rec.max_abs_y,), "#80b1d3", 0.0),),
        )
        ax.axhline(4.0, color="black", linestyle=":", linewidth=1.0, label="cut=4")
        ax.set_title(f"M_KK={mkk / 1000:g} TeV")
        ax.set_ylabel("max |Y_q|")
        ax.grid(True, axis="y", alpha=0.25)
    _hide_unused(axes, len(mkk_values))
    axes[0].legend(fontsize=7)
    fig.suptitle("Perturbativity distribution by exact r and M_KK")
    _save(fig, path)
    return path


def _plot_per_r_survival(analysis: Analysis, path: Path) -> Path:
    fig, ax = plt.subplots(figsize=(8.8, 5.4))
    r_values = _r_values(analysis)
    if not r_values:
        _blank_axis(ax, "No survival data found")
    for mkk in _mkk_values(analysis):
        xs, ys = [], []
        for r_value in r_values:
            stats = analysis.group_stats.get((r_value, mkk))
            if stats and stats.evaluated:
                xs.append(r_value)
                ys.append(stats.survives_strict / stats.evaluated)
        if xs:
            ax.plot(xs, ys, marker="o", label=f"{mkk / 1000:g} TeV")
    ax.set_xscale("log")
    ax.set_xlabel("r")
    ax.set_ylabel("Strict survival fraction")
    ax.set_title("Per-r survival evolution, each curve fixed M_KK")
    ax.set_ylim(-0.03, 1.03)
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=7, ncols=2)
    _save(fig, path)
    return path


def _plot_per_r_cq(analysis: Analysis, path: Path) -> Path:
    mkk_values = _mkk_values(analysis)
    fig, axes = _facet_axes(len(mkk_values), figsize=(13.0, max(3.4, 2.7 * math.ceil(max(len(mkk_values), 1) / 2))))
    if not mkk_values:
        _blank_axis(axes[0], "No c_Q samples found")
    colors = ("#1b9e77", "#d95f02", "#7570b3")
    for ax, mkk in zip(axes, mkk_values):
        for gen in range(3):
            xs, ys = [], []
            for r_value in _r_values(analysis):
                values = [
                    rec.c_q[gen]
                    for rec in _records(analysis, r_value, mkk)
                    if len(rec.c_q) > gen
                ]
                if values:
                    xs.append(r_value)
                    ys.append(float(np.median(values)))
            if xs:
                ax.plot(xs, ys, marker="o", color=colors[gen], label=f"c_Q{gen + 1}")
        ax.set_xscale("log")
        ax.set_title(f"M_KK={mkk / 1000:g} TeV")
        ax.set_xlabel("r")
        ax.set_ylabel("median c_Q")
        ax.set_ylim(0.25, 0.95)
        ax.grid(True, alpha=0.25)
    _hide_unused(axes, len(mkk_values))
    axes[0].legend(fontsize=8)
    fig.suptitle("c_Q localization evolution with r, faceted by M_KK")
    _save(fig, path)
    return path


def _plot_per_r_yukawa(analysis: Analysis, path: Path) -> Path:
    mkk_values = _mkk_values(analysis)
    fig, axes = _facet_axes(len(mkk_values), figsize=(13.0, max(3.4, 2.7 * math.ceil(max(len(mkk_values), 1) / 2))))
    if not mkk_values:
        _blank_axis(axes[0], "No Yukawa samples found")
    colors = ("#386cb0", "#386cb0", "#386cb0", "#fdb462", "#fdb462", "#fdb462")
    linestyles = ("-", "--", ":", "-", "--", ":")
    labels = ("u s1", "u s2", "u s3", "d s1", "d s2", "d s3")
    for ax, mkk in zip(axes, mkk_values):
        for index, label in enumerate(labels):
            xs, ys = [], []
            is_up = index < 3
            component = index if is_up else index - 3
            for r_value in _r_values(analysis):
                values = []
                for rec in _records(analysis, r_value, mkk):
                    singulars = rec.up_singular if is_up else rec.down_singular
                    if len(singulars) > component:
                        values.append(singulars[component])
                if values:
                    xs.append(r_value)
                    ys.append(float(np.median(values)))
            if xs:
                ax.plot(
                    xs,
                    ys,
                    marker="o",
                    color=colors[index],
                    linestyle=linestyles[index],
                    label=label,
                )
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title(f"M_KK={mkk / 1000:g} TeV")
        ax.set_xlabel("r")
        ax.set_ylabel("median singular value")
        ax.grid(True, alpha=0.25)
    _hide_unused(axes, len(mkk_values))
    axes[0].legend(fontsize=7, ncols=2)
    fig.suptitle("Yukawa singular-value evolution with r, faceted by M_KK")
    _save(fig, path)
    return path


def _boxplot_by_r(
    ax: plt.Axes,
    analysis: Analysis,
    mkk: float,
    *,
    series: Sequence[tuple[str, Any, str, float]],
) -> None:
    r_values = _r_values(analysis)
    width = min(0.18, 0.65 / max(len(series), 1))
    handles = []
    for label, getter, color, offset in series:
        data, positions = [], []
        for index, r_value in enumerate(r_values, start=1):
            values: list[float] = []
            for rec in _records(analysis, r_value, mkk):
                values.extend(float(x) for x in getter(rec) if _finite_or_none(x) is not None)
            if values:
                data.append(values)
                positions.append(index + offset)
        if data:
            bp = ax.boxplot(
                data,
                positions=positions,
                widths=width,
                patch_artist=True,
                showfliers=False,
            )
            for patch in bp["boxes"]:
                patch.set_facecolor(color)
                patch.set_alpha(0.65)
            for median in bp["medians"]:
                median.set_color("black")
            handles.append(plt.Line2D([0], [0], color=color, linewidth=6, alpha=0.65, label=label))
    ax.set_xticks(range(1, len(r_values) + 1), labels=[f"{r:g}" for r in r_values])
    ax.set_xlabel("r")
    if handles:
        ax.legend(handles=handles, fontsize=7, loc="best")


def _write_report(analysis: Analysis, path: Path, plot_paths: Sequence[Path]) -> None:
    total_rows = sum(stats.rows for stats in analysis.group_stats.values())
    total_eval = sum(stats.evaluated for stats in analysis.group_stats.values())
    total_skip = sum(stats.skipped for stats in analysis.group_stats.values())
    lines = [
        "# WQ Quark-Only Analysis Report",
        "",
        f"- input JSONL files: {len(analysis.input_paths)}",
        f"- rows: {total_rows}",
        f"- evaluated rows: {total_eval}",
        f"- skipped rows: {total_skip}",
        f"- malformed rows skipped: {analysis.malformed_rows}",
        f"- r values: {', '.join(f'{r:g}' for r in _r_values(analysis))}",
        f"- M_KK values [TeV]: {', '.join(f'{m / 1000:g}' for m in _mkk_values(analysis))}",
        "- grouping rule: every fraction and distribution is grouped by exact (quark_fit_r, M_KK); no rows are pooled across r or M_KK.",
        f"- singular-value source counts: {dict(sorted(analysis.singular_source_counts.items()))}",
        "",
        "## Plots",
        "",
        *[f"- plots/{plot.name}" for plot in plot_paths],
        "",
        "## Survival By Group",
        "",
        "| r | M_KK [TeV] | rows | evaluated | strict survival | inclusive survival |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for r_value, mkk in _group_keys(analysis):
        stats = analysis.group_stats[(r_value, mkk)]
        lines.append(
            "| "
            f"{r_value:g} | {mkk / 1000:g} | {stats.rows} | {stats.evaluated} | "
            f"{_fraction(stats.survives_strict, stats.evaluated)} | "
            f"{_fraction(stats.survives_inclusive, stats.evaluated)} |"
        )
    lines.extend(
        [
            "",
            "## Top HARD Vetoes By Group",
            "",
            "| r | M_KK [TeV] | rigorous top vetoes | proxy top vetoes |",
            "|---:|---:|---|---|",
        ]
    )
    for r_value, mkk in _group_keys(analysis):
        rigorous = _top_for_group(analysis, r_value, mkk, tag="rigorous")
        proxy = _top_for_group(analysis, r_value, mkk, tag="proxy")
        lines.append(
            f"| {r_value:g} | {mkk / 1000:g} | "
            f"{_format_top_list(rigorous)} | {_format_top_list(proxy)} |"
        )
    if analysis.singular_source_counts.get("seed_fallback", 0):
        lines.extend(
            [
                "",
                "Note: this input lacked fitted Yukawa singular values in some rows; "
                "the singular-value plots used the drawn seed matrices only for those rows.",
            ]
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _top_for_group(
    analysis: Analysis,
    r_value: float,
    mkk: float,
    *,
    tag: str,
    limit: int = 3,
) -> list[tuple[str, float]]:
    items = []
    for (r0, m0, pid, tag0), stats in analysis.constraint_stats.items():
        if r0 == r_value and m0 == mkk and tag0 == tag and stats.evaluated:
            fraction = stats.vetoed / stats.evaluated
            if fraction > 0.0:
                items.append((pid, fraction))
    return sorted(items, key=lambda item: (-item[1], item[0]))[:limit]


def _format_top_list(items: Sequence[tuple[str, float]]) -> str:
    if not items:
        return "none"
    return ", ".join(f"{pid} {value:.3f}" for pid, value in items)


def _top_constraint_ids(analysis: Analysis, *, tag: str, limit: int) -> list[str]:
    best: dict[str, float] = defaultdict(float)
    for (_r, _mkk, pid, tag0), stats in analysis.constraint_stats.items():
        if tag0 != tag or stats.evaluated == 0:
            continue
        best[pid] = max(best[pid], stats.vetoed / stats.evaluated)
    return [
        pid
        for pid, _value in sorted(best.items(), key=lambda item: (-item[1], item[0]))[:limit]
    ]


def _records(analysis: Analysis, r_value: float, mkk: float) -> list[SampleRecord]:
    reservoir = analysis.samples.get((r_value, mkk))
    return [] if reservoir is None else reservoir.records


def _r_values(analysis: Analysis) -> list[float]:
    return sorted({key[0] for key in analysis.group_stats})


def _mkk_values(analysis: Analysis) -> list[float]:
    return sorted({key[1] for key in analysis.group_stats})


def _group_keys(analysis: Analysis) -> list[tuple[float, float]]:
    return sorted(analysis.group_stats)


def _facet_axes(count: int, *, figsize: tuple[float, float]) -> tuple[plt.Figure, list[plt.Axes]]:
    count = max(1, count)
    cols = 2 if count > 1 else 1
    rows = int(math.ceil(count / cols))
    fig, axes_raw = plt.subplots(rows, cols, figsize=figsize, squeeze=False)
    return fig, [ax for row in axes_raw for ax in row]


def _hide_unused(axes: Sequence[plt.Axes], used: int) -> None:
    for ax in axes[used:]:
        ax.set_visible(False)


def _blank_axis(ax: plt.Axes, message: str) -> None:
    ax.text(0.5, 0.5, message, ha="center", va="center", transform=ax.transAxes)
    ax.set_xticks([])
    ax.set_yticks([])


def _save(fig: plt.Figure, path: Path) -> None:
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def _float_tuple(values: Any) -> tuple[float, ...]:
    if values is None:
        return ()
    try:
        arr = np.asarray(values, dtype=float).reshape(-1)
    except (TypeError, ValueError):
        return ()
    return tuple(float(x) for x in arr if math.isfinite(float(x)))


def _finite_or_none(value: Any) -> float | None:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if math.isfinite(out) else None


def _fraction(numerator: int, denominator: int) -> str:
    if denominator <= 0:
        return "n/a"
    return f"{numerator / denominator:.6g}"


def _default_output_root(input_path: Path) -> Path:
    if input_path.is_dir():
        return input_path
    return input_path.parent


def _build_argparser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--max-samples-per-group", type=int, default=5000)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _build_argparser().parse_args(argv)
    input_path = Path(args.input)
    output_root = Path(args.output_dir) if args.output_dir else _default_output_root(input_path)
    result = analyze(input_path, max_samples_per_group=int(args.max_samples_per_group))
    plots = write_outputs(result, output_root)
    print(f"[analyze_wq] input_jsonl={len(result.input_paths)} rows={sum(s.rows for s in result.group_stats.values())}")
    print(f"[analyze_wq] plots_dir={output_root / 'plots'}")
    print(f"[analyze_wq] report={output_root / 'analysis_report.md'}")
    for plot in plots:
        print(f"[analyze_wq] plot={plot}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
