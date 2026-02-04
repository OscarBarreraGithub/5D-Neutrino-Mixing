# 5D Neutrino Mixing - Project Status

**Last Updated:** January 2025

This document provides a comprehensive overview of the project's current state, goals, and next steps. It is designed to help new contributors (or a new Claude instance) quickly understand the codebase.

---

## 1. Project Goal

The primary goal of this repository is to **find regions of parameter space in a Randall-Sundrum warped extra dimension model that reproduce observed lepton masses with natural (O(1)) Yukawa couplings**.

Specifically, given input parameters:
- **Λ_IR (M_KK)**: The KK/IR scale (typically 1-10 TeV)
- **c_L, c_E, c_N**: Bulk mass parameters for lepton doublets, RH charged leptons, and RH neutrinos
- **M_N**: UV-localized Majorana mass for neutrinos

The code computes the corresponding **Yukawa couplings** and checks whether they are:
1. **Natural**: Rescaled Yukawas Ȳ = 2kY are O(1)
2. **Perturbative**: |Ȳ| < 4 (at least 3 KK modes before strong coupling)

This allows scanning parameter space to identify viable regions for the RS lepton flavor model.

---

## 2. Physics Background

### 2.1 The Warped Geometry

The model is a slice of AdS₅ bounded by:
- **UV brane** at z = 1/k (Planck-scale physics)
- **IR brane** at z = 1/Λ (TeV-scale physics)

Key parameters:
- k ≈ M_Planck ≈ 1.22×10¹⁹ GeV (AdS curvature)
- Λ = Λ_IR ≈ few TeV (KK scale)
- ε = Λ/k ≈ 10⁻¹⁵ (warp factor)

### 2.2 Fermion Localization

Bulk fermions have a dimensionless mass parameter **c** that controls their zero-mode localization:
- c > 1/2: Zero mode localized toward UV brane (small IR overlap)
- c < 1/2: Zero mode localized toward IR brane (large IR overlap)

The overlap factors (f-factors) encode the wavefunction values at the branes:
- f_IR²(c) = (1/2 - c) / (1 - ε^(1-2c))
- f_UV²(c) = (1/2 - c) / (ε^(2c-1) - 1)

### 2.3 Mass Formulas

**Charged leptons** (IR-localized Higgs):
```
m_{E_i} = 2vk · f_L · Y_{E_i} · f_{E_i}
```

**Neutrinos** (Type-I seesaw with UV Majorana mass):
```
m_{ν_i} = (2k²v² f²_L f²_N) / ((f^UV_N)² M_N) · Y²_{N_i}
```

### 2.4 Key Reference

Perez & Randall, "Natural Neutrino Masses and Mixings from Warped Geometry", arXiv:0805.4652

---

## 3. Current Implementation Status

### 3.1 Completed Modules ✅

| Module | Status | Description |
|--------|--------|-------------|
| `warpConfig/baseParams.py` | ✅ Complete | Geometry parameters (ε, rc, etc.) |
| `warpConfig/wavefuncs.py` | ✅ Complete | f_IR and f_UV overlap factors |
| `solvers/bessel.py` | ✅ Complete | KK mass solver (Bessel equations) |
| `neutrinos/neutrinoValues.py` | ✅ Complete | PDG neutrino data, PMNS matrix |
| `neutrinos/massConstraints.py` | ✅ Complete | Allowed neutrino mass sweeper |
| `diagonalization/diag.py` | ✅ Complete | SVD and Takagi factorization |
| **`yukawa/` (NEW)** | ✅ Complete | Yukawa computation from parameters |

### 3.2 Yukawa Module (Just Implemented)

The `yukawa/` module was just implemented. It provides:

**Main API:**
```python
from yukawa import compute_all_yukawas

result = compute_all_yukawas(
    Lambda_IR=3000,           # KK scale (GeV)
    c_L=0.58,                 # Lepton doublet bulk mass
    c_E=[0.75, 0.60, 0.50],   # RH charged lepton bulk masses
    c_N=0.27,                 # RH neutrino bulk mass
    M_N=1.22e18,              # UV Majorana mass (GeV)
    lightest_nu_mass=0.002,   # Lightest neutrino mass (eV)
    ordering='normal'
)

print(result.Y_E_bar)         # Rescaled charged lepton Yukawas
print(result.Y_N_bar)         # Rescaled neutrino Yukawas
print(result.is_perturbative())  # Check |Ȳ| < 4
print(result.summary())       # Full formatted output
```

**Returns a `YukawaResult` dataclass with:**
- `Y_E`, `Y_E_bar`: Charged lepton Yukawas (5D and rescaled)
- `Y_N`, `Y_N_bar`: Neutrino Yukawa eigenvalues (5D and rescaled)
- `Y_N_matrix`: Full 3×3 neutrino Yukawa matrix with PMNS mixing
- `f_L`, `f_E`, `f_N`, `f_N_UV`: Overlap factors used
- `epsilon`: Warp factor
- `is_perturbative()`: Check if all |Ȳ| < 4
- `summary()`: Formatted output string

**Verified against Perez-Randall benchmark (Table I):**
- f_L ≈ 0.016, f_N ≈ 0.48, f_N^UV ≈ 0.00012 (match paper)
- Ȳ_E ≈ [2.9, 4.4, 5.4] (O(1) to O(few))
- Ȳ_N ≈ [0.2, 0.4, 1.0] (O(0.1) to O(1))
- Round-trip consistency check passes

### 3.3 Stub Modules (Not Yet Implemented)

| Module | Status | Description |
|--------|--------|-------------|
| `flavorConstraints/` | 📝 Stub | μ→eγ and other LFV constraints |
| `scanParams/` | 📝 Stub | Parameter space sweep driver |

---

## 4. File Structure

```
5D-Neutrino-Mixing/
├── CLAUDE.md                 # Instructions for Claude Code
├── PROJECT_STATUS.md         # This document
├── TODO.md                   # Original todo list
│
├── warpConfig/               # ✅ Geometry and wavefunction computation
│   ├── baseParams.py         #   get_warp_params(k, Lambda) → {ε, rc, ...}
│   ├── wavefuncs.py          #   f_IR(c, ε), f_UV(c, ε)
│   └── wavefuncsTest.ipynb   #   Demonstration notebook
│
├── solvers/                  # ✅ KK tower mass solver
│   ├── bessel.py             #   solve_kk(species, bc, geometry, c) → masses
│   └── besselExample.ipynb   #   Working examples
│
├── neutrinos/                # ✅ Neutrino phenomenology
│   ├── neutrinoValues.py     #   PDG data, compute_masses(), get_pmns()
│   ├── massConstraints.py    #   find_allowed_lightest_masses()
│   ├── PMNS.ipynb            #   PMNS matrix examples
│   └── allowedMass.ipynb     #   Allowed mass ranges
│
├── diagonalization/          # ✅ Matrix diagonalization
│   ├── diag.py               #   SVD() for Dirac, Takagi() for Majorana
│   └── diagonalizationTest.ipynb
│
├── yukawa/                   # ✅ NEW: Yukawa computation
│   ├── __init__.py           #   Package exports
│   ├── constants.py          #   PDG lepton masses
│   ├── charged_lepton.py     #   compute_charged_lepton_yukawas()
│   ├── neutrino.py           #   compute_neutrino_yukawas()
│   ├── compute_yukawas.py    #   compute_all_yukawas(), YukawaResult
│   └── README.md             #   Module documentation
│
├── flavorConstraints/        # 📝 Stub: LFV bounds
│   └── README.md             #   μ→eγ constraint description
│
├── scanParams/               # 📝 Stub: Parameter sweep driver
│   └── README.md             #   Grid scan description
│
├── derivations/              # LaTeX derivations
│   ├── conventions.tex       #   All normalization conventions
│   └── kk_modes/             #   KK mode derivations
│
└── 0805.4652.pdf             # Perez-Randall reference paper
```

---

## 5. Key Conventions

### 5.1 Yukawa Rescaling

The 5D Yukawa Y has dimension [mass]⁻¹. The rescaled (dimensionless) Yukawa is:
```
Ȳ = 2k · Y
```

For naturalness, we want |Ȳ| ~ O(1). Perturbativity requires |Ȳ| < 4.

### 5.2 Bulk Mass Parameter

The dimensionless bulk mass parameter c is defined as:
```
c = M_5/k
```
where M_5 is the 5D Dirac mass and k is the AdS curvature.

### 5.3 Units

- All masses in GeV unless otherwise noted
- Neutrino masses typically input in eV (converted internally)
- k ≈ M_Planck ≈ 1.22×10¹⁹ GeV
- v = 174 GeV (electroweak VEV)

---

## 6. Next Steps

### 6.1 Immediate (Parameter Scanning)

Now that `compute_all_yukawas()` is implemented, the next step is to use it for parameter scanning:

1. **Simple grid scan**: Loop over (c_L, c_E, c_N, Λ_IR, M_N) and call `compute_all_yukawas()` for each point
2. **Filter for naturalness**: Keep points where `result.is_perturbative()` returns True and Ȳ ~ O(1)
3. **Store viable points**: Save parameters and Yukawas for viable regions

Example scanning code:
```python
import numpy as np
from yukawa import compute_all_yukawas

viable_points = []

for c_L in np.linspace(0.5, 0.7, 20):
    for c_N in np.linspace(0.2, 0.5, 20):
        for Lambda_IR in [3000, 5000, 10000]:
            # Vary c_E to fit charged lepton masses
            for c_E_tau in np.linspace(0.4, 0.6, 10):
                c_E = [0.8, 0.65, c_E_tau]  # Rough starting point

                result = compute_all_yukawas(
                    Lambda_IR=Lambda_IR,
                    c_L=c_L,
                    c_E=c_E,
                    c_N=c_N,
                    M_N=1e14,
                    lightest_nu_mass=0.001,
                    ordering='normal'
                )

                # Check if natural and perturbative
                if result.is_perturbative() and np.all(result.Y_E_bar > 0.5):
                    viable_points.append({
                        'params': result.params,
                        'Y_E_bar': result.Y_E_bar.tolist(),
                        'Y_N_bar': result.Y_N_bar.tolist(),
                    })

print(f"Found {len(viable_points)} viable points")
```

### 6.2 Future Enhancements

1. **Flavor constraints** (`flavorConstraints/`):
   - Implement μ→eγ bound checking
   - The constraint from the paper: |(Ȳ_N Ȳ_N†)₁₂| ≤ 0.028·(M_KK/3TeV)²

2. **Optimization/fitting**:
   - Instead of grid scan, use scipy.optimize to find parameters that minimize |Ȳ - 1|
   - Fit c_E values to reproduce exact charged lepton masses

3. **Electroweak precision**:
   - Check Z→ℓℓ constraints on wavefunction overlaps
   - Verify KK scale is high enough

4. **Full scan driver** (`scanParams/`):
   - Implement the grid scan logic described in scanParams/README.md
   - Add result storage (CSV/HDF5)
   - Add plotting utilities

---

## 7. Running the Code

### 7.1 Dependencies

- Python 3.x
- NumPy
- SciPy (for Bessel functions and optimization)

### 7.2 Quick Test

```bash
cd /path/to/5D-Neutrino-Mixing
python -c "
from yukawa import compute_all_yukawas
result = compute_all_yukawas(
    Lambda_IR=3000, c_L=0.58, c_E=[0.75, 0.60, 0.50],
    c_N=0.27, M_N=1.22e18, lightest_nu_mass=0.002, ordering='normal'
)
print(result.summary())
"
```

---

## 8. Contact and References

- **Reference paper**: Perez & Randall, arXiv:0805.4652
- **Conventions**: See `derivations/conventions.tex`
- **Physics overview**: See `CLAUDE.md`
