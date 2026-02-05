# scanParams - Parameter Space Scanner

Grid-scan driver for the RS lepton-sector parameter space.

The scanner evaluates `yukawa.compute_all_yukawas()` over configured grids,
then applies:

- perturbativity,
- naturalness,
- mu->e gamma LFV bound.

Results are returned as rows (list of dicts) and can be written to CSV.

Theory/parameter-prior note: [`THEORY_PRIORS.md`](THEORY_PRIORS.md)
Task tracking for this refactor: [`IMPLEMENTATION_TASKS.md`](IMPLEMENTATION_TASKS.md)

## Quick Start

```python
import numpy as np
from scanParams import ScanConfig, run_scan

config = ScanConfig(
    # Geometry
    k=1.2209e19,
    Lambda_IR_values=np.array([3000.0, 5000.0, 7000.0]),
    xi_KK=1.0,

    # Bulk masses
    c_L_values=np.linspace(0.52, 0.70, 10),
    c_N_values=np.linspace(0.15, 0.45, 10),
    c_E_fixed=[0.75, 0.60, 0.50],

    # UV Majorana mode
    MN_mode="fixed_ratio",
    MN_over_k=0.1,

    # Neutrino mass prior (broad mode)
    lightest_nu_mass_values=np.logspace(-5, np.log10(3e-2), 8),

    # Reproducibility metadata
    rng_seed_global=20260205,
    record_git_metadata=True,
)

rows = run_scan(config, output_csv="scan_results.csv", progress_every=200)
```

## Majorana Modes

`M_N` is controlled by `MN_mode`:

- `fixed_ratio`: use one value `MN_over_k` for all points.
- `scan_ratio`: scan `MN_over_k_values`.

`M_N` is always derived internally as:

- `M_N = MN_over_k * k`

### Example: scan `MN_over_k`

```python
config = ScanConfig(
    MN_mode="scan_ratio",
    MN_over_k_values=np.logspace(-6, 0, 7),
)
```

## LFV Convention

LFV is specified by an experimental BR bound and converted to internal `C`:

- input: `br_limit`, `prefac_br`, `lfv_reference_scale`, `xi_KK`
- derived per run: `lfv_C = sqrt(br_limit / prefac_br)`

Defaults:

- `br_limit = 1.5e-13` (MEG II 2025),
- `prefac_br = 4e-8`,
- `lfv_reference_scale = 3000 GeV`,
- `xi_KK = 1.0`, with `M_KK = xi_KK * Lambda_IR`.

To reproduce the Perez-Randall paper-era bound, set:

```python
config = ScanConfig(br_limit=1.2e-11)
```

## c_E Ordering Convention

If `sort_c_E_descending=True` (default), each sampled triplet is relabeled by sorting:

- `c_E1 >= c_E2 >= c_E3`

This is a labeling convention implemented by sorting values, not reject-until-ordered sampling.

## Anarchic Yukawa Scoring Mode

Optional anarchic-prior metadata/scoring can be added per point.

```python
from scanParams import AnarchyConfig, ScanConfig

config = ScanConfig(
    anarchy=AnarchyConfig(
        magnitude_min=1/3,
        magnitude_max=3,
        yN_overall_min=0.01,
        yN_overall_max=0.2,
        w_band=1.0,
        w_cond=1.0,
        w_fit=0.0,
    ),
    anarchy_min_score=-5.0,  # optional additional filter
    rng_seed_global=7,
)
```

Rows include:

- `anarchy_score`,
- `anarchy_band_penalty`,
- `anarchy_condition_penalty`,
- `anarchy_yN_overall`,
- score weights (`anarchy_w_band`, `anarchy_w_cond`, `anarchy_w_fit`).

## API

- `ScanConfig`: scanner configuration dataclass.
- `run_scan(config, output_csv=None, extra_filters=None, progress_every=100)`.
- `AnarchyConfig`: anarchic sampling/scoring configuration dataclass.

### `extra_filters`

Each extra filter takes a `YukawaResult` and returns `(passes, label)`.

```python
def require_tau_near_one(result):
    y_tau = abs(result.Y_E_bar[2])
    return (0.5 <= y_tau <= 2.0, "tau_yukawa")

rows = run_scan(config, extra_filters=[require_tau_near_one])
```

## CSV Output

CSV includes:

- scan point parameters: `Lambda_IR`, `M_KK`, `k`, `MN_over_k`, `M_N`, `c_*`, neutrino inputs,
- LFV metadata: `lfv_model`, `br_limit`, `prefac_br`, `lfv_C`, `xi_KK`, reference scale,
- reproducibility fields: `git_commit`, `dirty_tree`, global/sample RNG seeds,
- computed outputs: Yukawas, overlaps, filter booleans, `reject_reason`,
- optional anarchic metrics.

## Notes

- v1 scan currently enforces `ordering='normal'`.
- If `record_git_metadata=True`, commit hash and dirty-tree status are embedded in each row.
