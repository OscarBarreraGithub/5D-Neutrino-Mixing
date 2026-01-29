# Mass Formula Derivations

Derivations for charged lepton and neutrino mass formulas.

## Status: NEEDED

Core mass formulas need explicit derivation.

## Files

| File | Description | Status | Blocks Code |
|------|-------------|--------|-------------|
| `charged_leptons.tex` | m_E = 2vk f_L Y_E f_E | **NEEDED** | `yukawa/chargedYukawa.py` |
| `seesaw.tex` | Type-I seesaw derivation | **NEEDED** | `yukawa/neutrinoYukawa.py` |
| `neutrino_spectrum.tex` | Mass splittings, hierarchy | DONE | — |

---

## What's Done: neutrino_spectrum.tex

The neutrino mass spectrum from oscillation data is implemented in `neutrinos/neutrinoValues.py`:
- Mass-squared differences Δm²₂₁, Δm²₃₂
- Normal vs inverted hierarchy
- Sum constraint Σm_ν ≤ 0.082 eV

---

## What's Needed: charged_leptons.tex

### Statement
Derive the 4D charged lepton mass formula from the 5D IR-brane Yukawa coupling.

### Starting Point
5D action with Higgs and Yukawa localized on IR brane:
```
S ⊃ ∫ d⁴x √(-g_IR) δ(y - πr_c) [Y_E^(5D) H L̄ E + h.c.]
```

### Required Result
```
m_E = 2vk f_L Y_E f_E
```

### Key Questions to Address
1. **Where does the factor of 2 come from?** (brane delta-function normalization?)
2. How do f_L and f_E enter from the zero-mode overlap?
3. What is the relation between Y_E^(5D) and the dimensionless Y_E?

### Verification
- Electron mass: check that typical (c_L, c_E, Y_E) values give m_e ~ 0.511 MeV
- Hierarchy: larger c_E → smaller f_E → smaller mass

### Code Mapping
`yukawa/chargedYukawa.py::extract_yukawa()` (to be implemented)

---

## What's Needed: seesaw.tex

### Statement
Derive the light neutrino mass formula from the Type-I seesaw with UV-localized Majorana mass.

### Starting Point
1. Dirac Yukawa on IR brane: couples L to N
2. Majorana mass on UV brane: gives N a large mass M_N

### Required Result
In the universal limit (c_L, c_N, M_N all generation-independent):
```
m_ν = (2k²v² f_L² f_N²) / ((f_N^UV)² M_N) × Y_N²
```

### Key Questions to Address
1. **Why (f_N^UV)^(-2)?** How does canonical normalization enter?
2. What is the full 6×6 neutral lepton mass matrix structure?
3. How does the seesaw formula emerge from diagonalization?
4. What are the corrections from the approximation?

### Verification
- Build full mass matrix and diagonalize with Takagi
- Compare eigenvalues to seesaw formula
- Check that light masses are ~ eV scale for reasonable parameters

### Code Mapping
- `yukawa/neutrinoYukawa.py::invert_seesaw()` (to be implemented)
- `neutrinos/seesawValidation.py` (validation against full diagonalization)
