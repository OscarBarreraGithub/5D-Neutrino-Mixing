"""CFW-style comparison of the M_KK^min CDF.

Reads a run directory, computes M_KK^min from the stored per-system
Delta-F=2 ratios, and can overlay a convention-matched CFW projection
without changing the live post-audit constants in deltaf2.py.  The CFW
projection is necessarily a plot-layer rescaling/filtering of stored rows:
the scan rows do not contain enough information to rerun the CFW inverse
c-sampler or decompose epsilon_K operator-by-operator.

Usage:
    python scripts/rs_anarchy_cfw_comparison.py \
        --run scan_outputs/rs_anarchy_runC_<TS> \
        --out-dir results/figures/quark
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter

REPO = Path(__file__).resolve().parents[1]
DEFAULT_OUT = REPO / "results/figures/quark"
MKK_TICKS = [1, 2, 3, 5, 10, 20, 30, 50, 100, 200]

G_S_PERT = 1.05
DEFAULT_GS_STAR = 3.0

CFW_GENERIC_TEV = 21.0
CFW_PGB_TEV = 33.0

EPSILON_K_EXP = 2.228e-3
EPSILON_K_SM_DEFAULT = 2.161e-3
EPSILON_K_BUDGET_DEFAULT = abs(EPSILON_K_EXP - EPSILON_K_SM_DEFAULT)

CURRENT_BAGS = {
    "B_1_K": 0.5503,
    "B_4_K": 0.903,
    "B_5_K": 0.691,
}


@dataclass(frozen=True)
class CurveConfig:
    key: str
    label: str
    epsilon_k_scale: float = 1.0
    pdg_factor: float | None = None
    pdg_j_factor: float | None = None
    budget: float = EPSILON_K_BUDGET_DEFAULT
    budget_source: str = "BGS 2020/PDG 2024"
    bag_scale: float = 1.0
    bag_source: str = "FLAG 2024 defaults"


def _coerce_yaml_scalar(value: str) -> Any:
    value = value.strip()
    if value == "":
        return {}
    lower = value.lower()
    if lower in {"true", "false"}:
        return lower == "true"
    if lower in {"null", "none"}:
        return None
    try:
        return int(value)
    except ValueError:
        pass
    try:
        return float(value)
    except ValueError:
        return value.strip("'\"")


def _load_simple_yaml(path: Path) -> dict[str, Any]:
    """Load the small key/value YAML subset needed for bag overrides.

    PyYAML is intentionally not a dependency of this repo.  This parser accepts
    top-level ``key: value`` pairs plus one-level nested maps such as
    ``ratio_scales:`` followed by indented ``epsilon_K: 0.8``.
    """
    out: dict[str, Any] = {}
    current_map: dict[str, Any] | None = None
    current_indent = 0
    for lineno, raw in enumerate(path.read_text().splitlines(), start=1):
        line = raw.split("#", 1)[0].rstrip()
        if not line.strip():
            continue
        indent = len(line) - len(line.lstrip(" "))
        if ":" not in line:
            raise ValueError(f"{path}:{lineno}: expected 'key: value'")
        key, value = line.strip().split(":", 1)
        if indent == 0:
            parsed = _coerce_yaml_scalar(value)
            out[key] = parsed
            current_map = parsed if isinstance(parsed, dict) else None
            current_indent = indent
        else:
            if current_map is None or indent <= current_indent:
                raise ValueError(f"{path}:{lineno}: nested key without a parent map")
            current_map[key] = _coerce_yaml_scalar(value)
    return out


def _bag_epsilon_k_scale(path: Path | None) -> tuple[float, str]:
    if path is None:
        return 1.0, "FLAG 2024 defaults"
    data = _load_simple_yaml(path)
    ratio_scales = data.get("ratio_scales")
    if isinstance(ratio_scales, dict) and "epsilon_K" in ratio_scales:
        return float(ratio_scales["epsilon_K"]), f"{path} ratio_scales.epsilon_K"
    if "epsilon_K_scale" in data:
        return float(data["epsilon_K_scale"]), f"{path} epsilon_K_scale"

    bag_values = data.get("bags") if isinstance(data.get("bags"), dict) else data
    lr_scales = []
    for key in ("B_4_K", "B_5_K"):
        if key in bag_values:
            lr_scales.append(float(bag_values[key]) / CURRENT_BAGS[key])
    if lr_scales:
        scale = float(sum(lr_scales) / len(lr_scales))
        return scale, f"{path} mean(B_4_K,B_5_K) scale"
    if "B_1_K" in bag_values:
        return float(bag_values["B_1_K"]) / CURRENT_BAGS["B_1_K"], (
            f"{path} B_1_K scale"
        )
    return 1.0, f"{path} (no recognized epsilon_K bag keys)"


def _budget_from_args(args: argparse.Namespace) -> tuple[float, str]:
    if args.eps_k_sm is not None and args.eps_k_budget is not None:
        raise ValueError("Use either --eps-k-sm or --eps-k-budget, not both.")
    if args.eps_k_budget is not None:
        budget = float(args.eps_k_budget)
        source = f"explicit budget {budget:.6g}"
    elif args.eps_k_sm is not None:
        budget = abs(EPSILON_K_EXP - float(args.eps_k_sm))
        source = f"epsilon_K_SM={float(args.eps_k_sm):.6g}"
    else:
        budget = EPSILON_K_BUDGET_DEFAULT
        source = "BGS 2020/PDG 2024"
    if not math.isfinite(budget) or budget <= 0:
        raise ValueError(f"Invalid epsilon_K budget: {budget}")
    return budget, source


def _pdg_factor_from_args(args: argparse.Namespace) -> tuple[float | None, float | None]:
    if args.pdg_relative_tolerance is not None and args.pdg_factor is not None:
        raise ValueError("Use either --pdg-relative-tolerance or --pdg-factor, not both.")
    if args.pdg_relative_tolerance is not None:
        tol = float(args.pdg_relative_tolerance)
        if not (0.0 < tol < 1.0):
            raise ValueError("--pdg-relative-tolerance must lie between 0 and 1.")
        factor = 1.0 / (1.0 - tol)
    elif args.pdg_factor is not None:
        factor = float(args.pdg_factor)
    else:
        return None, None
    if factor <= 1.0 or not math.isfinite(factor):
        raise ValueError("PDG factor must be finite and greater than 1.")
    j_factor = float(args.pdg_j_factor) if args.pdg_j_factor is not None else factor
    return factor, j_factor


def _passes_curve_gate(row: dict[str, Any], curve: CurveConfig, require_pdg: bool) -> bool:
    if curve.pdg_factor is None:
        return (not require_pdg) or bool(row.get("passes_pdg"))
    tol = math.log(curve.pdg_factor)
    j_tol = math.log(curve.pdg_j_factor or curve.pdg_factor)
    return (
        float(row.get("up_log_max", math.inf)) <= tol
        and float(row.get("down_log_max", math.inf)) <= tol
        and float(row.get("ckm_log_max", math.inf)) <= tol
        and float(row.get("j_log", math.inf)) <= j_tol
    )


def _mkk_min_pert(row: dict[str, Any], epsilon_k_scale: float) -> float | None:
    ratios = row.get("deltaf2_ratios")
    if isinstance(ratios, dict):
        max_ratio = -math.inf
        for system, value in ratios.items():
            ratio = float(value)
            if system == "epsilon_K":
                ratio *= epsilon_k_scale
            max_ratio = max(max_ratio, ratio)
    else:
        max_ratio = float(row.get("max_ratio", math.nan)) * epsilon_k_scale
    if not math.isfinite(max_ratio) or max_ratio <= 0:
        return None
    return float(row["M_KK_GeV"]) / 1000.0 * math.sqrt(max_ratio)


def _load_curves(
    draws_path: Path,
    curves: list[CurveConfig],
    require_pdg: bool = True,
) -> dict[str, np.ndarray]:
    out: dict[str, list[float]] = {curve.key: [] for curve in curves}
    with draws_path.open() as fh:
        for line in fh:
            try:
                row = json.loads(line)
            except json.JSONDecodeError:
                continue
            if not row.get("ok"):
                continue
            for curve in curves:
                if not _passes_curve_gate(row, curve, require_pdg=require_pdg):
                    continue
                mkk = _mkk_min_pert(row, curve.epsilon_k_scale)
                if mkk is not None:
                    out[curve.key].append(mkk)
    return {key: np.asarray(values, dtype=float) for key, values in out.items()}


def _summarize(arr_pert: np.ndarray, gs_star: float) -> dict[str, float | int]:
    scaled = arr_pert * (gs_star / G_S_PERT)
    summary: dict[str, float | int] = {"n": int(arr_pert.size)}
    for q in [5, 25, 50, 75, 90, 95, 99]:
        summary[f"p{q:02d}_pert_TeV"] = float(np.percentile(arr_pert, q))
        summary[f"p{q:02d}_gs_star_TeV"] = float(np.percentile(scaled, q))
    return summary


def _parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run", type=str, required=True,
                   help="Run directory containing draws.jsonl")
    p.add_argument("--out-dir", type=str, default=str(DEFAULT_OUT))
    p.add_argument("--label-tag", type=str, default="")
    p.add_argument("--gs-star", type=float, default=DEFAULT_GS_STAR,
                   help="Dimensionless g_s^* convention for plotted M_KK values.")
    p.add_argument("--eps-k-sm", type=float, default=None,
                   help="Override epsilon_K^SM; budget is |epsilon_K^exp - value|.")
    p.add_argument("--eps-k-budget", type=float, default=None,
                   help="Direct epsilon_K NP-budget override.")
    p.add_argument("--bag-overrides", type=str, default=None,
                   help="Small YAML file with epsilon_K bag override scales.")
    p.add_argument("--pdg-factor", type=float, default=None,
                   help="Symmetric multiplicative PDG gate for masses, CKM, and J.")
    p.add_argument("--pdg-relative-tolerance", type=float, default=None,
                   help="Relative tolerance t; implemented as symmetric factor 1/(1-t).")
    p.add_argument("--pdg-j-factor", type=float, default=None,
                   help="Optional separate multiplicative gate for the Jarlskog invariant.")
    p.add_argument("--summary-out", type=str, default=None,
                   help="Write machine-readable percentile summary JSON.")
    p.add_argument("--no-plot", action="store_true",
                   help="Compute summaries without writing PDF/PNG figures.")
    return p.parse_args()


def main():
    args = _parse_args()
    run_dir = Path(args.run)
    draws = run_dir / "draws.jsonl"
    if not draws.exists():
        sys.stderr.write(f"ERROR: draws.jsonl not found in {run_dir}\n")
        sys.exit(1)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    try:
        matched_budget, matched_budget_source = _budget_from_args(args)
        matched_pdg_factor, matched_pdg_j_factor = _pdg_factor_from_args(args)
        bag_scale, bag_source = _bag_epsilon_k_scale(
            Path(args.bag_overrides) if args.bag_overrides else None
        )
    except ValueError as exc:
        sys.stderr.write(f"ERROR: {exc}\n")
        sys.exit(2)

    matched_scale = (EPSILON_K_BUDGET_DEFAULT / matched_budget) * bag_scale
    have_matched_overrides = any(
        value is not None
        for value in (
            args.eps_k_sm,
            args.eps_k_budget,
            args.bag_overrides,
            args.pdg_factor,
            args.pdg_relative_tolerance,
        )
    )

    curves = [
        CurveConfig(
            key="post_audit_default",
            label="post-audit default (BGS/FLAG, stored PDG gate)",
        )
    ]
    if have_matched_overrides:
        curves.append(
            CurveConfig(
                key="cfw_matched",
                label="CFW-matched projection",
                epsilon_k_scale=matched_scale,
                pdg_factor=matched_pdg_factor,
                pdg_j_factor=matched_pdg_j_factor,
                budget=matched_budget,
                budget_source=matched_budget_source,
                bag_scale=bag_scale,
                bag_source=bag_source,
            )
        )

    print(f"loading {draws} ...", flush=True)
    curve_arrays = _load_curves(draws, curves, require_pdg=True)
    relaxed_keys: set[str] = set()
    fallback_curves = [
        curve for curve in curves
        if curve_arrays[curve.key].size == 0 and curve.pdg_factor is None
    ]
    if fallback_curves:
        sys.stderr.write(
            "WARNING: no stored-PDG-passing draws for one or more curves; "
            "falling back to all 'ok' draws only for those stored-gate curves.\n"
        )
        fallback_arrays = _load_curves(draws, fallback_curves, require_pdg=False)
        for curve in fallback_curves:
            curve_arrays[curve.key] = fallback_arrays[curve.key]
            relaxed_keys.add(curve.key)

    nonempty_curves = [curve for curve in curves if curve_arrays[curve.key].size > 0]
    if not nonempty_curves:
        sys.stderr.write("ERROR: no usable draws for any requested curve.\n")
        sys.exit(1)

    summary = {
        "run": str(run_dir),
        "gs_pert": G_S_PERT,
        "gs_star": float(args.gs_star),
        "epsilon_K_exp": EPSILON_K_EXP,
        "epsilon_K_budget_default": EPSILON_K_BUDGET_DEFAULT,
        "curves": {},
    }

    for curve in curves:
        arr = curve_arrays[curve.key]
        if arr.size == 0:
            sys.stderr.write(f"WARNING: curve {curve.key} has zero accepted rows.\n")
            continue
        curve_summary = _summarize(arr, args.gs_star)
        curve_summary.update({
            "label": curve.label,
            "epsilon_k_scale": float(curve.epsilon_k_scale),
            "epsilon_k_budget": float(curve.budget),
            "epsilon_k_budget_source": curve.budget_source,
            "bag_scale": float(curve.bag_scale),
            "bag_source": curve.bag_source,
            "pdg_factor": curve.pdg_factor,
            "pdg_j_factor": curve.pdg_j_factor,
            "relaxed_to_all_ok": curve.key in relaxed_keys,
        })
        summary["curves"][curve.key] = curve_summary

    if args.summary_out:
        summary_path = Path(args.summary_out)
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
        print(f"wrote {summary_path}")

    if not args.no_plot:
        fig, ax = plt.subplots(figsize=(8.5, 5.5))

        colors = {"post_audit_default": "C0", "cfw_matched": "C3"}
        styles = {"post_audit_default": "-", "cfw_matched": "--"}
        for curve in nonempty_curves:
            arr_plot = np.sort(curve_arrays[curve.key] * (args.gs_star / G_S_PERT))
            cdf = (np.arange(1, arr_plot.size + 1) / arr_plot.size) * 100.0
            ax.plot(
                arr_plot,
                cdf,
                color=colors.get(curve.key, None),
                linewidth=2.4,
                linestyle=styles.get(curve.key, "-"),
                label=f"{curve.label} (n={arr_plot.size:,})",
            )

        # CFW abstract benchmarks
        ax.axvline(CFW_GENERIC_TEV, color="grey", linewidth=1.2, linestyle=":",
                   label=fr"CFW generic anarchic ({CFW_GENERIC_TEV:.0f} TeV)")
        ax.axvline(CFW_PGB_TEV, color="black", linewidth=1.2, linestyle="--",
                   label=fr"CFW PGB Higgs ({CFW_PGB_TEV:.0f} TeV)")

        ax.axhline(50.0, color="k", linewidth=0.7, alpha=0.4, linestyle=":")
        ax.axhline(95.0, color="k", linewidth=0.7, alpha=0.4, linestyle=":")

        ax.set_xscale("log")
        ax.xaxis.set_major_locator(FixedLocator(MKK_TICKS))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:g}"))
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.set_xlabel(
            fr"$M_{{\mathrm{{KK}}}}$ [TeV], $g_s^\star={args.gs_star:g}$ convention"
        )
        ax.set_ylabel(
            r"fraction of accepted draws with "
            r"$M_{\mathrm{KK}}^{\min}\leq M_{\mathrm{KK}}$ [%]"
        )
        title = "RS-anarchy comparison to Csaki--Falkowski--Weiler"
        if have_matched_overrides:
            title += " (post-audit vs convention-matched)"
        if args.label_tag:
            title = f"{title}  [{args.label_tag}]"
        ax.set_title(title, fontsize=11)
        ax.set_xlim(0.5, 250.0)
        ax.set_ylim(0, 102)
        ax.legend(fontsize=8.5, loc="lower right")
        ax.grid(True, which="both", alpha=0.25, linestyle=":")

        fig.tight_layout()

        stem = "rs_anarchy_cfw_comparison"
        if args.label_tag:
            stem = f"{stem}_{args.label_tag}"
        for ext in ("pdf", "png"):
            path = out_dir / f"{stem}.{ext}"
            fig.savefig(path, dpi=200)
            print(f"wrote {path}")

    print()
    print("=== M_KK^min percentiles (TeV) ===")
    for curve in nonempty_curves:
        print(f"\n{curve.key}: {curve.label}")
        print(
            f"  n={curve_arrays[curve.key].size:,}, "
            f"epsilon_K scale={curve.epsilon_k_scale:.6g}, "
            f"budget={curve.budget:.6g}, bag_scale={curve.bag_scale:.6g}"
        )
        for q in [5, 25, 50, 75, 90, 95, 99]:
            v_pert = np.percentile(curve_arrays[curve.key], q)
            v_gss = v_pert * (args.gs_star / G_S_PERT)
            print(f"  p{q:>2}: pert {v_pert:7.2f}    g_s*={args.gs_star:g} {v_gss:7.2f}")


if __name__ == "__main__":
    main()
