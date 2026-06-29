# 5D Neutrino Mixing — Agent Notes

This file is agent-facing context for working in the repo. Canonical user-facing
documentation lives in `README.md` and the module `README.md` files. Superseded
planning notes should be recovered from git history when needed.

This repository originated as a lepton-sector parameter sweep in a
Randall-Sundrum warped extra dimension (geometric localization of charged-lepton
and neutrino masses), and that machinery remains. **Active documentation and the
current paper work are quark-sector / full-catalog**, not lepton: the audited
ΔF=2 + RS-electroweak constraints and the 103-constraint flavor/collider catalog.
The lepton-sector directories are follow-up scope.

All work lives on `main` (post-2026-05-25 consolidation; tagged
`v2026q2-catalog-complete`). The authoritative current-state docs are
`docs/FLOOR_SUMMARY.md`, `docs/STATE_OF_PROJECT.md`, and the June 2026
collaborator report `reports/collaborator_2026-06/CONTENT.md`. The legacy
quark-sector methodology note `docs/quark_scan_methodology_note.tex`/`.pdf` is
ΔF=2-only and pre-audit (carries a SUPERSEDED banner); its flavor-only floors are
not the current floors. The Cloudflare-deployed catalog website builds from
`main` with root `flavor_catalog/website/`.

**Corrected minimal-RS floor (post-audit, June 2026):** typical ~30 TeV from the
tunable `epsilon_K`; existence ~18-20 TeV from irreducible oblique S,T,U; Z→bb
~5 TeV after the B1 fix (the old "25-30 TeV Z→bb floor" was a B1 bug, not real).

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

The original computational stack is lepton-oriented. Post-consolidation, the
active paper work is quark-sector, with `quarkConstraints/`, `qcd/`, and
`flavor_catalog/` now first-class top-level packages alongside the lepton
modules.

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
├── quarkConstraints/  # ✅ Active quark-sector MFV scans, benchmarks, and ΔF=2 validation
│
├── qcd/               # ✅ 4-loop QCD running and threshold decoupling utilities
│
├── flavor_catalog/    # ✅ 95 PRIMARY + 8 SECONDARY constraints (103 total) plus Astro website
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
   - Flavor bounds: μ→eγ dipole operator (scan default uses the published
     MEG II 2025 limit `br_limit = 1.5e-13`, which implies `C ≈ 1.94e-3`;
     set `C=0.02` for Perez-Randall reproduction)
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

## Working Notes

- Canonical install and validation commands live in `README.md`.
- Scanner defaults should be read from `scanParams/scan.py`; the current LFV
  default is the published MEG II 2025 bound `br_limit = 1.5e-13`.
- Historical snapshots may describe superseded APIs; use git history for
  context, not as the canonical source of truth.

---

## References

- Primary paper references are linked from `README.md` and the module
  `README.md` files.
- **derivations/**: LaTeX documents with step-by-step calculations
- Agashe et al. (hep-ph/0412089): RS flavor framework
- Perez et al.: Lepton sector constraints in warped models
