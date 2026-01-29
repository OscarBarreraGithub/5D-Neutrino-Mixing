# Derivations

This directory contains all physics derivations for the 5D Neutrino Mixing project. Each derivation must be completed and approved before the corresponding code can be implemented.

## Directory Structure

| Directory | Topic | Status |
|-----------|-------|--------|
| [geometry/](geometry/) | Warp factor, brane positions, coordinates | DONE |
| [zero_modes/](zero_modes/) | f-factors, g(z) profiles | PARTIAL |
| [kk_modes/](kk_modes/) | KK mass quantization, boundary conditions | PARTIAL |
| [masses/](masses/) | Charged lepton and neutrino mass formulas | NEEDED |
| [flavor/](flavor/) | PMNS, LFV constraints, perturbativity | PARTIAL |
| [integrals/](integrals/) | Loop integrals (L_ijk) | LATER |

## Status Tags

- `[DONE]` — Derivation complete, code implemented
- `[PARTIAL]` — Started but incomplete
- `[NEEDED]` — Required before code can be written
- `[LATER]` — Not blocking current work

## Derivation File Format

Each `.tex` file should contain:

1. **Statement**: What formula/result we need
2. **Starting point**: Which 5D action terms are relevant
3. **Derivation**: Step-by-step calculation
4. **Final result**: Boxed formula in our conventions
5. **Code mapping**: Which function(s) implement this
6. **Verification**: How to check correctness (limits, literature values)

## Conventions

All derivations use the conventions defined in [conventions.tex](conventions.tex). Key points:

### Coordinates
- Conformal coordinate z: UV brane at z_h = 1/k, IR brane at z_v = e^(πkr_c)/k
- Warp factor: ε = Λ/k = e^(-πkr_c) = z_h/z_v

### Bulk Mass Parameter
- 5D fermion mass: M_5 = k(c - 1/2)
- α parameter: α = |c + 1/2|

### Boundary Conditions
- **(++)** LH zero mode: J_{α∓1}/Y_{α∓1} matching (upper sign for c > -1/2)
- **(--)** RH zero mode: J_α/Y_α matching
- **(NN)** Gauge bosons: J_0/Y_0 matching

### Mass Formulas
- Charged leptons: m_E = 2vk f_L Y_E f_E
- Neutrinos (seesaw): m_ν = (2k²v² f_L² f_N²) / ((f_N^UV)² M_N) · Y_N²

## Derivation → Code Dependencies

| Derivation | Blocks Code In | Priority |
|------------|----------------|----------|
| `zero_modes/g_profiles.tex` | `warpConfig/wavefuncs.py` | High |
| `kk_modes/fermion_pp.tex` | `solvers/bessel.py` | High |
| `kk_modes/normalization.tex` | `solvers/bessel.py` | High |
| `masses/charged_leptons.tex` | `yukawa/chargedYukawa.py` | High |
| `masses/seesaw.tex` | `yukawa/neutrinoYukawa.py` | High |
| `flavor/mu_to_e_gamma.tex` | `flavorConstraints/` | Medium |
| `flavor/yukawa_naturalness.tex` | `yukawa/` | Medium |
| `integrals/L_ijk.tex` | Future | Low |

## Building LaTeX

To compile a derivation:
```bash
cd derivations/<subdir>
pdflatex <filename>.tex
```

Or use the main build system from the `build/` directory.
