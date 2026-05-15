"""Compare typical Yukawa-element magnitudes between the envelope-fit
accepted points and the RS-anarchy ensemble draws.

Outputs:
  results/figures/quark/yukawa_size_envelope_vs_anarchic.{pdf,png}

Reads:
  scan_outputs/dense_20260506T141321/derived/accepted_points_with_yukawas.csv
  scan_outputs/rs_anarchy_20260507T030811/draws.jsonl   (streamed)
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter

REPO = Path(__file__).resolve().parents[1]
ENV_CSV = REPO / "scan_outputs/dense_20260506T141321/derived/accepted_points_with_yukawas.csv"
ANARCH_JSONL = REPO / "scan_outputs/rs_anarchy_20260507T030811/draws.jsonl"
OUT = REPO / "results/figures/quark"
OUT.mkdir(parents=True, exist_ok=True)


def envelope_yukawa_magnitudes_split() -> dict:
    df = pd.read_csv(ENV_CSV)
    diag = []
    off = []
    for f in ("u", "d"):
        for i in (1, 2, 3):
            for j in (1, 2, 3):
                mag = np.hypot(df[f"Y_{f}_{i}{j}_re"].values,
                               df[f"Y_{f}_{i}{j}_im"].values)
                (diag if i == j else off).append(mag)
    return {"diag": np.concatenate(diag), "off": np.concatenate(off)}


def anarchic_yukawa_magnitudes(max_rows: int = 50_000) -> np.ndarray:
    """Stream JSONL, accept only rows with Y_u and Y_d present."""
    mags = []
    n = 0
    with ANARCH_JSONL.open() as fh:
        for line in fh:
            try:
                row = json.loads(line)
            except json.JSONDecodeError:
                continue
            yu = row.get("Y_u")
            yd = row.get("Y_d")
            if yu is None or yd is None:
                continue
            for Y in (yu, yd):
                arr = np.asarray(Y, dtype=complex)
                mags.append(np.abs(arr).ravel())
            n += 1
            if n >= max_rows:
                break
    if not mags:
        return _anarchic_yukawa_from_prior(int(1e6))
    return np.concatenate(mags)


def _anarchic_yukawa_from_prior(n: int) -> np.ndarray:
    """Fallback: redraw from the documented prior if Y arrays not stored."""
    rng = np.random.default_rng(20260507)
    half = 1.5
    floor = 0.1
    # rejection sample
    out = []
    target = n
    while len(out) < target:
        block = target - len(out)
        re = rng.uniform(-half, half, size=block)
        im = rng.uniform(-half, half, size=block)
        mag = np.hypot(re, im)
        keep = mag >= floor
        out.extend(mag[keep].tolist())
    return np.asarray(out[:target])


def main():
    env = envelope_yukawa_magnitudes_split()
    n_pts = (env["diag"].size + env["off"].size) // 18
    print(f"envelope diag: {env['diag'].size:,}  off: {env['off'].size:,} from {n_pts:,} pts")
    print(f"  diag median {np.median(env['diag']):.3f}  p5 {np.percentile(env['diag'],5):.3f}"
          f"  p95 {np.percentile(env['diag'],95):.3f}")
    print(f"  off  median {np.median(env['off']):.3f}  p5 {np.percentile(env['off'],5):.3f}"
          f"  p95 {np.percentile(env['off'],95):.3f}")

    anar = anarchic_yukawa_magnitudes(max_rows=50_000)
    print(f"anarchic: {anar.size:,} |Y_ij| values")
    if len(anar) and not np.isfinite(anar).all():
        anar = anar[np.isfinite(anar)]

    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    bins = np.logspace(-4, 1, 80)

    ax.hist(env["diag"], bins=bins, density=True, histtype="stepfilled",
            alpha=0.45, color="C3",
            label=fr"Envelope fit, diagonal $|Y_{{ii}}|$  ({env['diag'].size:,} entries)")
    ax.hist(env["off"], bins=bins, density=True, histtype="stepfilled",
            alpha=0.45, color="C1",
            label=fr"Envelope fit, off-diagonal $|Y_{{i\neq j}}|$  ({env['off'].size:,} entries)")
    ax.hist(anar, bins=bins, density=True, histtype="step",
            color="C0", linewidth=2.0,
            label=r"Anarchic prior ($|Y|\geq 0.1$, $\mathrm{Re,Im}\sim\mathcal{U}(-1.5,1.5)$)")

    ax.set_xscale("log")
    y_ticks = [0.001, 0.01, 0.1, 1.0, 10.0]
    ax.xaxis.set_major_locator(FixedLocator(y_ticks))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:g}"))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlim(0.001, 12.0)
    ax.set_xlabel(r"$|Y_{f,ij}|$  (5D Yukawa entry magnitude)")
    ax.set_ylabel("probability density")
    ax.set_title("Envelope-fit Yukawa entries (diagonal vs.\\ off-diagonal) vs.\\ anarchic prior")
    ax.axvline(1.0, color="grey", linestyle=":", linewidth=1, alpha=0.7)
    ax.legend(fontsize=8.8, loc="upper left")
    ax.grid(True, which="both", alpha=0.25, linestyle=":")
    fig.tight_layout()

    for ext in ("pdf", "png"):
        path = OUT / f"yukawa_size_envelope_vs_anarchic.{ext}"
        fig.savefig(path, dpi=200)
        print("wrote", path)


if __name__ == "__main__":
    main()
