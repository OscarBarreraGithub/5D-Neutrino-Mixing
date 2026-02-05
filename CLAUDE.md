# 5D Neutrino Mixing — Project Guide

This repository implements a parameter sweep for the lepton sector in a Randall-Sundrum warped extra dimension. The goal is to find regions of parameter space that reproduce observed charged lepton masses and viable light neutrino masses using geometric localization rather than hierarchical Yukawa couplings.

---

## Physics Overview

### The Warped Geometry

The spacetime is a slice of AdS₅ bounded by two branes:
- **UV brane** at z = 1/k (Planck-scale physics)
- **IR brane** at z = 1/Λ (TeV-scale physics)

Key parameters:
| Symbol | Meaning | Typical Value |
|--------|---------|---------------|
| k | AdS curvature scale | ~M_Planck ≈ 1.22 × 10¹⁹ GeV |
| Λ | IR (KK) scale | ~few TeV |
| ε = Λ/k | Warp factor | ~10⁻¹⁵ |
| r_c | Compactification radius | Fixed by ε = e^(-πkr_c) |

### Fermion Localization

Bulk fermions have a dimensionless mass parameter **c** that controls their zero-mode localization:
- **c > 1/2**: Zero mode peaks at UV brane (exponentially small IR overlap)
- **c < 1/2**: Zero mode peaks at IR brane (O(1) IR overlap)

This is the central mechanism: mass hierarchies arise from geometry, not from hierarchical 5D couplings.

**Convention:** throughout this repo we define
```
c = M_5 / k
```
so that `M_5 = k c` and the Bessel order used in KK equations is `alpha = |c + 1/2|`.

### Overlap Factors (f-factors)

The normalized zero-mode wavefunction values at the branes:

**IR brane overlap:**
$$f_{\text{IR}}^2(c) = \frac{\tfrac{1}{2} - c}{1 - \varepsilon^{1-2c}}$$

**UV brane overlap** (for right-handed neutrinos):
$$f_{\text{UV}}^2(c) = \frac{\tfrac{1}{2} - c}{\varepsilon^{2c-1} - 1}$$

These factors control all effective 4D couplings.

### Charged Lepton Masses

With Higgs and Yukawas localized on the IR brane:
$$m_{E_i} = 2\,v\,k\, f_{L_i} \cdot Y_{E_i} \cdot f_{E_i}$$

where v ≈ 174 GeV is the electroweak vev.

**Strategy**: Fix c_L (universal for doublets), vary c_{E_i} for each generation, then solve for diagonal Y_{E_i} that reproduce m_e, m_μ, m_τ.

### Neutrino Masses (Seesaw)

Right-handed neutrinos N have a UV-localized Majorana mass M_N. In the universal limit (c_L, c_N, M_N all generation-independent):

$$m_{\nu_i} \simeq \frac{2\,k^2\,v^2\, f_L^2\, f_N^2}{(f_N^{\text{UV}})^2\, M_N}\, Y_{N_i}^2$$

The factor (f_N^UV)⁻² reflects canonical normalization of the Majorana term. The Y_{N_i} eigenvalues control mass splittings; their ratios determine the neutrino spectrum hierarchy.

### PMNS Mixing

Working in the charged-lepton mass basis, mixing is encoded in the neutrino Yukawa:
$$Y_N \to V_{\text{PMNS}} \cdot \text{diag}(Y_{N_1}, Y_{N_2}, Y_{N_3})$$

---

## Code Architecture

```
5D-Neutrino-Mixing/
├── warpConfig/        # ✅ Geometry parameters and f-factor computation
│   ├── baseParams.py  #   get_warp_params(k, Lambda) → {ε, r_c, z_h, z_v, ...}
│   └── wavefuncs.py   #   f_IR(c, ε), f_UV(c, ε) — vectorized
│
├── solvers/           # ✅ KK tower mass solver
│   └── bessel.py      #   solve_kk(species, bc, geometry, c) → masses
│
├── neutrinos/         # ✅ Neutrino phenomenology
│   ├── neutrinoValues.py    #   compute_masses(), get_pmns()
│   └── massConstraints.py   #   find_allowed_lightest_masses()
│
├── diagonalization/   # ✅ Matrix diagonalization
│   └── diag.py        #   SVD() for Dirac, Takagi() for Majorana matrices
│
├── yukawa/            # ✅ Yukawa computation from parameters
│   ├── compute_yukawas.py   #   compute_all_yukawas() → YukawaResult
│   ├── charged_lepton.py    #   Charged lepton Yukawa inversion
│   ├── neutrino.py          #   Neutrino Yukawa (seesaw inversion)
│   └── constants.py         #   PDG lepton masses
│
├── flavorConstraints/ # ✅ LFV bounds (μ→eγ)
│
├── scanParams/        # ✅ Parameter space sweep driver
│
└── derivations/       # LaTeX derivations (detailed physics)
```

### Key Dependencies
- NumPy (arrays, linear algebra)
- SciPy (Bessel functions, optimization, matrix functions)

---

## Workflow

### Quick Method (Recommended)

Use the `yukawa` module to compute everything in one call:

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

# Check results
print(result.Y_E_bar)         # Rescaled charged lepton Yukawas
print(result.Y_N_bar)         # Rescaled neutrino Yukawas
print(result.is_perturbative())  # Check |Ȳ| < 4
print(result.summary())       # Full output
```

### Step-by-Step Method

For more control, you can proceed step-by-step:

1. **Set geometry**: Choose (k, Λ) → compute warp parameters
   ```python
   from warpConfig.baseParams import get_warp_params
   params = get_warp_params(k=1.22e19, Lambda_IR=3e3)
   ```

2. **Choose bulk masses**: c_L (doublets), c_N (RH neutrinos), c_{E_i} (RH charged leptons)

3. **Compute f-factors**:
   ```python
   from warpConfig.wavefuncs import f_IR, f_UV
   fL = f_IR(c_L, params['epsilon'])
   fN = f_IR(c_N, params['epsilon'])
   fN_uv = f_UV(c_N, params['epsilon'])
   ```

4. **Fit charged lepton Yukawas**: Invert m_E = 2vk f_L Y_E f_E for each generation

5. **Check neutrino mass scale**: Evaluate seesaw formula with chosen M_N and Y_N pattern

6. **Apply constraints**:
   - Yukawa perturbativity: |Ȳ| = 2k|Y| < O(few)
   - Flavor bounds: μ→eγ dipole operator (scan default uses MEG II 2024 C ≈ 4.33e-3; set C=0.02 for Perez–Randall)
   - EW precision: KK mass scale sufficiently high

7. **Collect viable points**: Store parameters that pass all filters

---

## Key Physical Inputs

| Parameter | Physical Meaning | Scan Range |
|-----------|------------------|------------|
| k | AdS curvature | Fixed ~M_Pl |
| Λ | IR/KK scale | 1–10 TeV |
| c_L | Lepton doublet bulk mass | 0.3–0.7 |
| c_{E_i} | RH charged lepton bulk masses | 0.3–0.9 |
| c_N | RH neutrino bulk mass | 0.3–0.7 |
| M_N | UV Majorana mass | 10¹⁰–10¹⁵ GeV |

**Targets**:
- m_e, m_μ, m_τ: PDG values
- Σm_ν ≤ 0.082 eV (cosmological bound)
- Δm²₂₁, Δm²₃₂: oscillation data

---

## Current Status

**See [PROJECT_STATUS.md](PROJECT_STATUS.md) for detailed implementation status and next steps.**

### Recently Completed (January 2025)

The `yukawa/` module is now fully implemented. It computes Yukawa couplings from RS parameters:

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

print(result.Y_E_bar)         # Rescaled charged lepton Yukawas Ȳ_E
print(result.Y_N_bar)         # Rescaled neutrino Yukawas Ȳ_N
print(result.is_perturbative())  # Check |Ȳ| < 4
```

### Module Status

| Module | Status | Description |
|--------|--------|-------------|
| `warpConfig/` | ✅ Complete | Geometry and f-factor computation |
| `solvers/` | ✅ Complete | KK mass solver |
| `neutrinos/` | ✅ Complete | PDG data, PMNS, mass constraints |
| `diagonalization/` | ✅ Complete | SVD and Takagi factorization |
| `yukawa/` | ✅ Complete | Yukawa computation from parameters |
| `flavorConstraints/` | ✅ Complete | μ→eγ NDA dipole bound (paper + MEG II 2024) |
| `scanParams/` | ✅ Complete | Grid-scan driver with perturbativity/naturalness/LFV filters |

### Next Steps

1. **Run scans + summarize**: use `scanParams` to map viable regions and write a short results note.
2. **Constraint upgrades**: add/confirm EW precision and any additional LFV constraints beyond μ→eγ.
3. **Optimization**: fit parameters to minimize |Ȳ - 1| after initial grid results.

---

## References

- **Research Paper.pdf**: Full derivations and theoretical motivation
- **derivations/**: LaTeX documents with step-by-step calculations
- Agashe et al. (hep-ph/0412089): RS flavor framework
- Perez et al.: Lepton sector constraints in warped models
