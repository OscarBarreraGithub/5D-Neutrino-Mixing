# Constraining a Warped 5th Dimension

This repository explores the constraints on warped five-dimensional Randall-Sundrum models.

[`neutrinos`](neutrinos)
- Computes neutrino mass spectra, the PMNS matrix, and allowed parameter ranges under experimental constraints.

[`diagonalization`](diagonalization)
- Implements mass matrix diagonalization routines supporting both singular value decomposition (SVD) and Takagi factorization.

[`warpConfig`](warpConfig)
- Provides utilities to compute derived 5D warp parameters and the fermion zero-mode overlap
functions used throughout the project.

[`solvers`](solvers)
- This solver finds the Kaluza–Klein (KK) masses in Randall–Sundrum models by solving the Bessel-function boundary conditions that quantize 5D bulk fields. The bulk equations reduce to Bessel’s equation in conformal (z) coordinates; the UV/IR boundary conditions give transcendental equations whose roots are the KK masses.
