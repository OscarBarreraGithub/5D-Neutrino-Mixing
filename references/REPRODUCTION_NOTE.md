# Reproduction note — constraint plots vs literature

How to regenerate the per-draw plot quantities and run the comparison notebook
`notebooks/constraint_plots_vs_literature.ipynb`.

## Environment

The repo's numerics use the `ising_bootstrap` conda env. On this cluster the
env's own `libstdc++` must be on `LD_LIBRARY_PATH` (the system one is too old
for SciPy 1.17's `GLIBCXX_3.4.29`):

```bash
export LD_LIBRARY_PATH=$HOME/.conda/envs/ising_bootstrap/lib:$LD_LIBRARY_PATH
PY=$HOME/.conda/envs/ising_bootstrap/bin/python
```

(Equivalently `conda activate ising_bootstrap`, which is what the scan sbatch
scripts do.)

## Step 1 — extract + validate the per-draw quantities

`scripts/extract_plot_quantities.py` **replays** a sample of stored scan draws
deterministically through the repo's own fixed-code pipeline and extracts the
quantities the literature figures need but the scan does not store as numerics
(`delta g_L^b, delta g_R^b, S, T`, plus `|eps_K^NP|` and the D-mixing M12
pieces).

For each sampled draw it:

1. rebuilds the exact `QuarkFitSeed` from the stored anarchic seed matrices
   (`params.quark_yukawa_seed`) via the scan's own `_svd_seed_parts`;
2. re-runs `fit_quark_sector(default_quark_targets(), r=0.25, seed=..., Lambda_IR,
   k, max_nfev=120, fit_orientation=True)` (deterministic least-squares);
3. builds the point with `point_builder.build_from_rs_ew_inputs(..., ew_model=
   "minimal_rs", include_fermion_kk_mixing=True, ...)` +
   `compute_quark_kk_gluon_couplings` + `make_point`;
4. evaluates the real T010/T011/EW001/K001/C001/C002/B001/B003 constraints and
   reads off the adapter diagnostics.

```bash
$PY scripts/extract_plot_quantities.py \
    --scan-dir scan_outputs/wq_quarkonly_20260622T090807 \
    --per-tile 250 \
    --out scan_outputs/plot_quantities.parquet
```

Runtime ~15 min (one `RSEWSpectrum` + spline cache build per M_KK tile dominates;
250 draws/tile x 10 tiles = 2500 draws). Output: `scan_outputs/plot_quantities.parquet`.

### Validation gate (hard)

The script recomputes `ratio_T010` and `ratio_EW001` from the replayed adapters
and compares to the **stored** ratios per draw. Result of the reference run:

```
VALIDATION ratio_T010 : n=2500  median_rel=0  max_rel=0  max_abs=0  n(rel>1e-6)=0
VALIDATION ratio_EW001: n=2500  median_rel=0  max_rel=0  max_abs=0  n(rel>1e-6)=0
GATE PASS (threshold 1e-6)
```

i.e. the replay reproduces the exact fixed-code scan points **bit-identically**,
so the recomputed `(g_L^b, g_R^b)` and `(S, T)` are the genuine scan points.

## Step 2 — build + execute the notebook

The notebook is generated from a script (so it is reproducible and diffable):

```bash
$PY notebooks/_build_constraint_plots_vs_literature.py     # writes the .ipynb
$PY -m jupyter nbconvert --to notebook --execute --inplace \
    --ExecutePreprocessor.timeout=600 \
    notebooks/constraint_plots_vs_literature.ipynb
```

The notebook reads:
- `scan_outputs/plot_quantities.parquet` (replayed Zbb / S-T / D-mixing points),
- `scan_outputs/wq_quarkonly_20260622T090807/constraint_matrix.parquet` (full
  100k quark-sector ratios for the density scatters),
- `scan_outputs/fix100k_minimal_20260622T080053/constraint_matrix.parquet`
  (full catalog incl. `mu -> e gamma` L001).

It produces 7 embedded figures, ordered by user priority:

1. **Z -> bb (g_L^b, g_R^b) plane** — CGHNP 0807.4937 Fig. 8 (PRIORITY)
2. **S, T oblique plane** — CGHNP 0807.4937 Fig. 4 (PRIORITY)
3. epsilon_K — Bauer et al. 0912.1625 Fig. 4
4. D0 mixing — Gedalia et al. 0906.1879 Fig. 1
5. Delta m_d / Delta m_s — own scatter (Blanke 0809.1073 inputs)
6. mu -> e gamma — Agashe et al. hep-ph/0606021 / Perez-Randall 0805.4652
7. summary: median ratio vs M_KK for all constraints

Per user priority, no collider exclusion plot is built; the CMS-B2G-25-009
`M_KK >= 5.5 TeV` cutoff is drawn as a vertical line where relevant.

## Key constants used (provenance)

- Zbb SM reference point (CGHNP Fig. 8): `g_L^b = -0.42114, g_R^b = 0.077420`.
- Zbb ellipse: `R_b^0 = 0.21629(66), A_b = 0.923(20), A_FB^{0,b} = 0.0992(16)`
  with the LEP/SLC correlation matrix (CGHNP Eq. 173), mapped to the (g_L, g_R)
  plane via CGHNP Eq. 171.
- S-T ellipse (ours / live EW001 anchor): `S = 0.026(75), T = 0.047(66),
  rho = 0.90`; paper (CGHNP Eq. 148): `S = 0.07(10), T = 0.16(10), rho = 0.85`.
- EW001 proxy: `c_S = 30`, `v = 246.21965`, `sin^2 = 0.23122`, `L = 35`,
  `x_1 = 2.4487`, `M_KK = x_1 Lambda_IR` (from `quarkConstraints/oblique_stu.py`).
- mu -> e gamma limits: MEG II 2025 `1.5e-13` (scan default) and paper-era MEGA
  `1.2e-11`.
