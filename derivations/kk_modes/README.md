# KK Mode Derivations

Derivations for Kaluza-Klein mass quantization and wavefunctions.

## Status: PARTIAL

Basic Bessel equation setup exists; boundary condition details and normalization needed.

## Files

| File | Description | Status | Blocks Code |
|------|-------------|--------|-------------|
| `bessel_equation.tex` | 5D EOM → Bessel equation | PARTIAL | — |
| `fermion_pp.tex` | (++) BC for LH zero mode | **NEEDED** | `solvers/bessel.py` |
| `fermion_mm.tex` | (--) BC for RH zero mode | **NEEDED** | `solvers/bessel.py` |
| `gauge_nn.tex` | (NN) BC for gauge bosons | **NEEDED** | `solvers/bessel.py` |
| `normalization.tex` | N_n², b_α formulas | **NEEDED** | `solvers/bessel.py` |

---

## What's Partially Done: bessel_equation.tex

The existing `KK mode masses.tex` derives how the 5D equation of motion reduces to Bessel's equation in conformal coordinates.

---

## What's Needed: fermion_pp.tex

### Statement
Derive the KK mass quantization condition for fermions with a **left-handed zero mode** (++ boundary conditions).

### Required Result
For α = |c + 1/2|:
```
J_{α∓1}(m_n z_h) / Y_{α∓1}(m_n z_h) = J_{α∓1}(m_n z_v) / Y_{α∓1}(m_n z_v)
```
- Upper sign (α-1) for c > -1/2
- Lower sign (α+1) for c < -1/2

### Key Questions to Address
1. Why does the sign of c + 1/2 determine whether we use α-1 or α+1?
2. How does this relate to the chirality of the zero mode?

### Code Mapping
`solvers/bessel.py::solve_kk(species='fermion', bc='++')`

---

## What's Needed: fermion_mm.tex

### Statement
Derive the KK mass quantization condition for fermions with a **right-handed zero mode** (-- boundary conditions).

### Required Result
```
J_α(m_n z_h) / Y_α(m_n z_h) = J_α(m_n z_v) / Y_α(m_n z_v)
```

### Code Mapping
`solvers/bessel.py::solve_kk(species='fermion', bc='--')`

---

## What's Needed: gauge_nn.tex

### Statement
Derive the KK mass quantization for gauge bosons with Neumann BCs at both branes.

### Required Result
```
J_0(m_n z_v) / Y_0(m_n z_v) = J_0(m_n z_h) / Y_0(m_n z_h)
```

### Note
The zero mode (n=0) is massless and corresponds to the 4D gauge boson.

### Code Mapping
`solvers/bessel.py::solve_kk(species='gauge', bc='NN')`

---

## What's Needed: normalization.tex

### Statement
Derive the normalization factors N_n and Bessel mixing coefficients b_α for KK wavefunctions.

### Required Results

**Fermion KK wavefunction:**
```
g_{x_i}(c, z) = (z/z_h) × (1 / N_n √(πr_c)) × [J_α(m_n z) + b_α Y_α(m_n z)]
```

**Normalization:**
```
N_n² = (1/2πr_c) × [z_v² (J_α(m_n z_v) + b_α Y_α(m_n z_v))²
                   - z_h² (J_α(m_n z_h) + b_α Y_α(m_n z_h))²]
```

**Bessel mixing coefficient:**
```
-b_α(m_n) = J_{α∓1}(m_n z_h) / Y_{α∓1}(m_n z_h)   (for ++ modes)
-b_α(m_n) = J_α(m_n z_h) / Y_α(m_n z_h)           (for -- modes)
```

### Code Mapping
Will extend `solvers/bessel.py` to return N_n and b_α alongside masses.

---

## Current Implementation Status

The current `solvers/bessel.py`:
- ✅ Solves for KK masses via root-finding
- ✅ Returns b_n coefficients
- ⚠️ Uses simplified α-1 / α for ++/-- (needs α∓1 based on sign of c+1/2)
- ❌ Does not compute N_n normalization
- ❌ Does not return full g(z) wavefunctions
