# Constraining a Warped 5th Dimension

This repository contains numerical tools for Randall-Sundrum lepton-sector studies:
zero-mode overlaps, KK spectra, Yukawa inversion, lepton-flavor constraints, and
parameter scans.

## Install

```bash
pip install -e .[dev]
```

## Validate

```bash
ruff check .
pytest -q
python scripts/benchmark_perez_randall.py
```

## Notes

- `scanParams.ScanConfig` defaults to the published **MEG II 2025** bound
  `br_limit = 1.5e-13`; the derived scan coefficient `lfv_C` is computed per run.
- Tracked notebooks are exploratory analysis artifacts and should remain output-free.
- Historical status and planning notes live under [`docs/archive`](docs/archive).

## Packages

[`neutrinos`](neutrinos)
- Computes neutrino mass spectra, the PMNS matrix, and allowed parameter ranges under experimental constraints.

[`diagonalization`](diagonalization)
- Implements mass matrix diagonalization routines supporting both singular value decomposition (SVD) and Takagi factorization.

[`warpConfig`](warpConfig)
- Provides utilities to compute derived 5D warp parameters and the fermion zero-mode overlap
functions used throughout the project.

[`yukawa`](yukawa)
- Computes charged lepton and neutrino Yukawa couplings by inverting the RS mass formulas.

[`solvers`](solvers)
- This solver finds the Kaluza–Klein (KK) masses in Randall–Sundrum models by solving the Bessel-function boundary conditions that quantize 5D bulk fields. The bulk equations reduce to Bessel’s equation in conformal (z) coordinates; the UV/IR boundary conditions give transcendental equations whose roots are the KK masses.

[`flavorConstraints`](flavorConstraints)
- Implements the μ→eγ NDA dipole bound, keeping the Perez-Randall paper coefficient
  available while scans default to the published MEG II 2025 limit.

[`scanParams`](scanParams)
- Grid-scan driver to sweep RS lepton-sector parameters and filter by perturbativity, naturalness, and LFV bounds.

[`qcd`](qcd)
- Computes the QCD running coupling α_s(μ) with threshold matching.
