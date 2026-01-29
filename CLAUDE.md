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
├── warpConfig/        # Geometry parameters and f-factor computation
│   ├── baseParams.py  #   get_warp_params(k, Lambda) → {ε, r_c, z_h, z_v, ...}
│   └── wavefuncs.py   #   f_IR(c, ε), f_UV(c, ε) — vectorized
│
├── solvers/           # KK tower mass solver
│   └── bessel.py      #   solve_kk(species, bc, geometry, c) → masses
│
├── neutrinos/         # Neutrino phenomenology
│   ├── neutrinoValues.py    #   compute_masses(), get_pmns()
│   └── massConstraints.py   #   find_allowed_lightest_masses()
│
├── diagonalization/   # Matrix diagonalization
│   └── diag.py        #   SVD() for Dirac, Takagi() for Majorana matrices
│
├── flavorConstraints/ # LFV bounds (μ→eγ)
│
├── yukawa/            # Extract Y_E from observed lepton masses
│
├── scanParams/        # Parameter space sweep driver
│
└── derivations/       # LaTeX derivations (detailed physics)
```

### Key Dependencies
- NumPy (arrays, linear algebra)
- SciPy (Bessel functions, optimization, matrix functions)

---

## Workflow

A typical parameter scan proceeds as:

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
   - Flavor bounds: μ→eγ dipole operator
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

See [TODO.md](TODO.md) for outstanding items. Key areas:
- **Stable**: warpConfig, neutrinos, diagonalization
- **In progress**: KK mass derivations, seesaw precision checks
- **Planned**: High-precision diagonalization, uncertainty propagation

---

## References

- **Research Paper.pdf**: Full derivations and theoretical motivation
- **derivations/**: LaTeX documents with step-by-step calculations
- Agashe et al. (hep-ph/0412089): RS flavor framework
- Perez et al.: Lepton sector constraints in warped models
