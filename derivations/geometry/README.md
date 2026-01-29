# Geometry Derivations

Derivations for the warped 5D geometry setup.

## Status: DONE

The geometry formulas are implemented in `warpConfig/baseParams.py`.

## Files

| File | Description | Status | Blocks Code |
|------|-------------|--------|-------------|
| `warp_factor.tex` | ε, r_c, z_h, z_v relations | DONE | — |
| `brane_positions.tex` | UV/IR brane setup | DONE | — |

## Key Results

### Warp Factor
```
ε = Λ/k = e^(-πkr_c)
```

### Brane Positions (conformal coordinate)
```
z_h = 1/k        (UV brane)
z_v = 1/Λ        (IR brane)
```

### Compactification Radius
```
r_c = ln(k/Λ) / (πk)
```

## Code Mapping

- `warpConfig/baseParams.py::get_warp_params(k, Lambda_IR)`
  - Returns: `{k, Lambda_IR, epsilon, rc, z_h, z_v, warp_log}`

## Verification

- [x] ε = z_h / z_v ✓
- [x] Λ = k × ε ✓
- [x] Large hierarchy: k ~ M_Pl, Λ ~ TeV → ε ~ 10^(-15) ✓
