# Integral Derivations

Derivations for overlap integrals needed in loop calculations.

## Status: LATER

Not blocking current work. Will be needed for advanced LFV calculations.

## Files

| File | Description | Status | Blocks Code |
|------|-------------|--------|-------------|
| `L_ijk.tex` | Triple overlap integrals | LATER | Future |

---

## What's Needed (Later): L_ijk.tex

### Statement
Derive the overlap integral:
```
L_ijk = ∫_{z_h}^{z_v} (g_i g_j g_k f_z) / (k z) dz
```

### Context
These integrals appear in loop-level calculations for flavor-violating processes. They encode the overlap of three bulk profiles (e.g., two fermions and a gauge boson KK mode).

### Dependencies
- Requires g(z) profiles from `zero_modes/g_profiles.tex`
- Requires KK wavefunctions from `kk_modes/normalization.tex`

### Priority
Low — not needed for the initial parameter sweep. Will become relevant when computing radiative corrections and precise LFV predictions.
