# Flavor Derivations

Derivations for PMNS mixing, flavor constraints, and perturbativity bounds.

## Status: PARTIAL

PMNS implemented; LFV constraints and perturbativity bounds needed.

## Files

| File | Description | Status | Blocks Code |
|------|-------------|--------|-------------|
| `pmns_convention.tex` | PMNS parameterization | DONE | — |
| `mu_to_e_gamma.tex` | μ→eγ dipole constraint | **NEEDED** | `flavorConstraints/` |
| `yukawa_naturalness.tex` | Perturbativity bounds | **NEEDED** | `yukawa/` |

---

## What's Done: pmns_convention.tex

The PMNS matrix is implemented in `neutrinos/neutrinoValues.py::get_pmns()`:
- Standard PDG parameterization
- Three mixing angles (θ₁₂, θ₂₃, θ₁₃)
- Dirac CP phase δ
- Majorana phases (α, β)

Convention:
```
|ν_α⟩ = Σᵢ V*_αi |ν_i⟩
Y_N = V_PMNS × diag(Y_N1, Y_N2, Y_N3)
```

---

## What's Needed: mu_to_e_gamma.tex

### Statement
Derive the constraint on neutrino Yukawa couplings from the μ→eγ branching ratio.

### Physics Background
In RS models, flavor-changing neutral currents arise from KK mode exchange. The radiative decay μ→eγ provides a strong constraint via a dipole operator.

### Expected Constraint Form (from literature)
```
|(Ȳ_N Ȳ_N†)₁₂| ≤ C × (M_KK / 3 TeV)²
```
where Ȳ_N = 2k Y_N and C ≈ 0.028 (to be verified).

### Key Questions to Address
1. What is the NDA (naive dimensional analysis) estimate for the dipole coefficient?
2. How does the constraint scale with M_KK?
3. What is the precise numerical coefficient C?
4. Are there other LFV processes (μ→3e, μ-e conversion) that give stronger bounds?

### Literature Reference
Perez et al. and Agashe et al. on RS flavor constraints.

### Verification
- Cross-check numerical coefficient against published values
- Test that viable parameter points pass the constraint

### Code Mapping
`flavorConstraints/muToEGamma.py::check_mu_to_e_gamma()` (to be implemented)

---

## What's Needed: yukawa_naturalness.tex

### Statement
Establish perturbativity bounds on rescaled Yukawa couplings.

### Background
The dimensionless 5D Yukawa Y couples to curvature k. The relevant expansion parameter is:
```
Ȳ = 2k Y
```

### Required Result
Perturbativity condition:
```
|Ȳ| < Λ_cutoff   (e.g., 4π or stricter)
```

### Key Questions to Address
1. What is the appropriate perturbativity bound? (4π, √4π, or something else?)
2. Is there a different bound for charged lepton vs neutrino Yukawas?
3. How does this relate to loop corrections?

### Code Mapping
`yukawa/perturbativity.py::check_perturbativity()` (to be implemented)
