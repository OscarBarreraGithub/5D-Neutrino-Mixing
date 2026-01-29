# Zero Mode Derivations

Derivations for fermion zero-mode wavefunctions and overlap factors.

## Status: PARTIAL

f-factors implemented; g(z) profiles needed.

## Files

| File | Description | Status | Blocks Code |
|------|-------------|--------|-------------|
| `f_factors.tex` | f_IR, f_UV derivation | DONE | — |
| `g_profiles.tex` | Full g(z) wavefunctions | **NEEDED** | `warpConfig/wavefuncs.py` |

## What's Implemented

### f_IR (IR brane overlap)
```
f_IR²(c) = (1/2 - c) / (1 - ε^(1-2c))
```
- Handles c = 1/2 singularity analytically
- Implemented in `warpConfig/wavefuncs.py::f_IR()`

### f_UV (UV brane overlap)
```
f_UV²(c) = (1/2 - c) / (ε^(2c-1) - 1)
```
- Used for right-handed neutrino Majorana mass
- Implemented in `warpConfig/wavefuncs.py::f_UV()`

---

## What's Needed: g_profiles.tex

### Statement
Derive the full z-dependent zero-mode profile g(z) for a bulk fermion with parameter c.

### Starting Point
5D fermion action in AdS₅ with bulk mass M_5 = k(c - 1/2).

### Required Result
The normalized wavefunction g(z) satisfying:
```
∫_{z_h}^{z_v} dz (k/z) |g(z)|² = 1
```

### Expected Form
```
g(z) = N × (z/z_h)^(1/2 - c)  or similar
```
where N is a normalization constant.

### Code Mapping
Will be implemented as `warpConfig/wavefuncs.py::g_profile(c, z, epsilon)`

### Verification
- g(z_v) should relate to f_IR
- g(z_h) should relate to f_UV
- Normalization integral = 1
