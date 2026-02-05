# scanParams — Parameter Space Scanner

Grid-scan driver for the RS lepton-sector parameter space. Loops over bulk mass parameters `(c_L, c_N)` (and optionally `c_E`) at fixed geometry and Majorana mass, evaluates Yukawa couplings, and filters for perturbativity, naturalness, and the mu-to-e-gamma bound. Results are written to CSV.

## Quick Start

```python
import numpy as np
from scanParams import ScanConfig, run_scan

config = ScanConfig(
    c_L_values=np.linspace(0.50, 0.70, 21),
    c_N_values=np.linspace(0.20, 0.50, 21),
)
results = run_scan(config, output_csv="scan_results.csv")
```

**LFV default:** `ScanConfig.lfv_C` defaults to the MEG II 2024 bound
(\(C \approx 4.33\times10^{-3}\)). Set `lfv_C=0.02` to reproduce the
Perez–Randall (paper-era) constraint.

## Full Example

```python
import numpy as np
from scanParams import ScanConfig, run_scan

config = ScanConfig(
    # Fixed physics
    Lambda_IR=5000.0,
    M_N=1e14,
    lightest_nu_mass=0.001,
    ordering='normal',

    # Scan grids
    c_L_values=np.linspace(0.50, 0.70, 21),
    c_N_values=np.linspace(0.20, 0.50, 21),

    # Fixed charged-lepton bulk masses (default)
    c_E_fixed=[0.75, 0.60, 0.50],

    # Filter thresholds
    max_Y_bar=4.0,
    naturalness_range=(0.1, 4.0),
    # Default uses MEG II 2024 bound (C ≈ 4.33e-3); use 0.02 to reproduce Perez–Randall.
    lfv_C=0.00433,
    lfv_reference_scale=3000.0,
)

results = run_scan(config, output_csv="scan_results.csv", progress_every=50)

# Count accepted points
accepted = [r for r in results if r['passes_all']]
print(f"{len(accepted)} / {len(results)} points accepted")
```

### Scanning over c_E (5D grid)

```python
config = ScanConfig(
    c_L_values=np.linspace(0.55, 0.65, 5),
    c_N_values=np.linspace(0.25, 0.35, 5),
    c_E_grid=[
        np.linspace(0.70, 0.80, 3),   # c_E1 (electron)
        np.linspace(0.55, 0.65, 3),   # c_E2 (muon)
        np.linspace(0.45, 0.55, 3),   # c_E3 (tau)
    ],
)
# 5 x 5 x 3 x 3 x 3 = 675 points
results = run_scan(config, output_csv="scan_5d.csv")
```

### Custom filters

```python
def require_tau_yukawa_near_one(result):
    """Require |Y_E_bar_tau| between 0.5 and 2."""
    y_tau = abs(result.Y_E_bar[2])
    ok = 0.5 <= y_tau <= 2.0
    return (ok, "tau_yukawa")

results = run_scan(config, extra_filters=[require_tau_yukawa_near_one])
```

## CSV Columns

| Column | Description |
|--------|-------------|
| `c_L`, `c_N`, `c_E1`, `c_E2`, `c_E3` | Bulk mass parameters |
| `Lambda_IR`, `M_N`, `lightest_nu_mass`, `ordering` | Fixed physics inputs |
| `Y_E_bar_1`, `Y_E_bar_2`, `Y_E_bar_3` | Rescaled charged-lepton Yukawas |
| `Y_N_bar_1`, `Y_N_bar_2`, `Y_N_bar_3` | Rescaled neutrino Yukawas |
| `f_L`, `f_N`, `f_N_UV` | Overlap factors |
| `max_Y_bar` | max(\|Y_E_bar\|, \|Y_N_bar\|) |
| `perturbative` | All \|Y_bar\| < max_Y_bar? |
| `natural` | All \|Y_bar\| in naturalness_range? |
| `lfv_passes` | mu-to-e-gamma constraint satisfied? |
| `lfv_ratio` | LHS/RHS of the LFV bound (>1 violates) |
| `passes_all` | All filters passed? |
| `reject_reason` | Semicolon-separated failure labels |

## Plotting

```python
import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt("scan_results.csv", delimiter=',', names=True, dtype=None, encoding='utf-8')
mask = data['passes_all'].astype(bool)
plt.scatter(data['c_L'][mask], data['c_N'][mask], c='green', label='accepted')
plt.scatter(data['c_L'][~mask], data['c_N'][~mask], c='red', alpha=0.2, label='rejected')
plt.xlabel('$c_L$'); plt.ylabel('$c_N$'); plt.legend(); plt.show()
```
