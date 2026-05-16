# Constraining a Warped 5th Dimension

This repository contains numerical tools for Randall-Sundrum flavor studies
across the lepton and quark sectors: zero-mode overlaps, KK spectra, Yukawa
inversion, lepton-flavor constraints, quark-sector MFV workflows, and
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

- `Lambda_IR` denotes the geometric IR scale `1 / z_v`. Physical first-KK
  masses are sector-dependent roots times `Lambda_IR`.
- `scanParams.ScanConfig` defaults to the published **MEG II 2025** bound
  `br_limit = 1.5e-13`; the derived scan coefficient `lfv_C` is computed per run.
- The repo's default LFV convention is `M_KK = Lambda_IR` (`xi_KK = 1.0`).
  Physical first-KK mass conventions are available only as explicit utilities.
- `scripts/benchmark_perez_randall.py` keeps a historical filename, but the
  validated point is a repo-local, paper-inspired benchmark rather than a
  literal reproduction of Perez–Randall Eq. (10) / Table I.
- `scripts/audit_perez_randall_consistency.py` documents the conclusion of the
  reproduction audit: the displayed Eq. (10) neutrino Yukawas do not
  numerically reproduce Eq. (7) when combined with Eq. (6) and Eq. (11).
- Tracked notebooks live under [`notebooks/`](notebooks/) as executed analysis
  artifacts.
- Tracked paper figure exports are limited to the figures referenced by
  [`docs/quark_scan_methodology_note.tex`](docs/quark_scan_methodology_note.tex);
  older exploratory quark figures are kept out of the submission figure set.
- Tracked paper docs under [`docs/`](docs/) include the current
  [`quark_scan_methodology_note.tex`](docs/quark_scan_methodology_note.tex)
  and the older [`quark_scan_assumptions_compact.tex`](docs/quark_scan_assumptions_compact.tex).

Paper branch note: the active methodology note is
[`docs/quark_scan_methodology_note.tex`](docs/quark_scan_methodology_note.tex)
on `paper/quark-scan-2026q2`. The current paper scope is quark-sector only;
lepton-sector packages remain repo tools for follow-up work.

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

[`quarkConstraints`](quarkConstraints)
- Active quark-sector MFV program with an exploratory `repo_v1` lane, a frozen
  `paper_0710_1869` lane, and a planned `modern` production lane. See
  [`quarkConstraints/README.md`](quarkConstraints/README.md) and
  [`quarkConstraints/PAPER_0710_1869.md`](quarkConstraints/PAPER_0710_1869.md).

[`qcd`](qcd)
- Computes the QCD running coupling α_s(μ) with threshold matching.
