# Yukawa Computation Module

Computes Yukawa couplings from RS model parameters by inverting the mass formulas.

## Quick Start

```python
from yukawa import compute_all_yukawas

result = compute_all_yukawas(
    Lambda_IR=3000,           # 3 TeV KK scale
    c_L=0.58,                 # Lepton doublet bulk mass
    c_E=[0.75, 0.60, 0.50],   # RH charged lepton bulk masses
    c_N=0.27,                 # RH neutrino bulk mass
    M_N=1.22e18,              # M_Pl/10 Majorana mass
    lightest_nu_mass=0.002,   # 2 meV
    ordering='normal'
)

print(result.summary())
print(f"Perturbative: {result.is_perturbative()}")
```

## Physics

### Charged Lepton Yukawas

The charged lepton mass with IR-localized Higgs is:
$$
m_{E_i} = 2\,v\,k\, f_L\, Y_{E_i}\, f_{E_i}
$$

Inverting:
$$
Y_{E_i} = \frac{m_{E_i}}{2\,v\,k\, f_L\, f_{E_i}}
$$

The rescaled (dimensionless) Yukawa is:
$$
\bar{Y}_{E_i} = 2k\,Y_{E_i} = \frac{m_{E_i}}{v\, f_L\, f_{E_i}}
$$

### Neutrino Yukawas (Seesaw)

In the universal limit (c_L, c_N, M_N generation-independent):
$$
m_{\nu_i} = \frac{2\,k^2\,v^2\, f_L^2\, f_N^2}{(f_N^{\text{UV}})^2\, M_N}\, Y_{N_i}^2
$$

Inverting:
$$
Y_{N_i} = \sqrt{\frac{m_{\nu_i}\, (f_N^{\text{UV}})^2\, M_N}{2\,k^2\,v^2\, f_L^2\, f_N^2}}
$$

### Perturbativity

For naturalness, we want $|\bar{Y}| \sim \mathcal{O}(1)$.

Perturbativity requires $|\bar{Y}| < 4$ (at least 3 KK modes before strong coupling).

## API Reference

### `compute_all_yukawas()`

Main entry point. Returns a `YukawaResult` dataclass containing:

- `Y_E`, `Y_E_bar`: Charged lepton Yukawas (5D and rescaled)
- `Y_N`, `Y_N_bar`: Neutrino Yukawa eigenvalues (5D and rescaled)
- `Y_N_matrix`: Full 3×3 neutrino Yukawa matrix with PMNS mixing
- `f_L`, `f_E`, `f_N`, `f_N_UV`: Overlap factors used
- `epsilon`: Warp factor
- `is_perturbative()`: Check if all |Ȳ| < 4
- `summary()`: Formatted output string

### `compute_charged_lepton_yukawas()`

Lower-level function for charged leptons only.

### `compute_neutrino_yukawas()`

Lower-level function for neutrinos only.

## Module Structure

```
yukawa/
├── __init__.py           # Package exports
├── constants.py          # PDG lepton masses
├── charged_lepton.py     # Charged lepton Yukawa computation
├── neutrino.py           # Neutrino Yukawa computation
└── compute_yukawas.py    # Main unified API
```

## Reference

Perez & Randall, "Natural Neutrino Masses and Mixings from Warped Geometry", arXiv:0805.4652
