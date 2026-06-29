# quarkConstraints

Quark-sector implementation notes for the Randall-Sundrum MFV model of
Fitzpatrick, Perez, and Randall:

- arXiv: [0710.1869](https://arxiv.org/abs/0710.1869)
- PDF: [0710.1869](https://arxiv.org/pdf/0710.1869)
- Journal: Phys. Rev. Lett. 100, 171604 (2008)

The paper’s core claim is that 5D minimal flavor violation can collapse the
generic RS flavor/CP problem into a controlled next-to-MFV structure. The
key low-energy dial is the parameter `r`, which suppresses down-sector flavor
violation and can make KK scales near 2 TeV viable.

> **Scope note (this 2 TeV figure is the old 0710.1869/FPR MFV lane, not the
> current minimal floor).** The audited minimal non-custodial floor is *not*
> 2 TeV: it is a typical ~30 TeV from `epsilon_K` and an irreducible existence
> ~18–20 TeV from oblique S,T,U (post-B1 Z→bb is only ~5 TeV). See
> [`docs/FLOOR_SUMMARY.md`](../docs/FLOOR_SUMMARY.md). The 2 TeV claim is the
> FPR-MFV-`r` mechanism in isolation, which this package implements as one lane.

## Canonical Docs

- [`README.md`](./README.md): collaborator-facing quark package overview
- [`PAPER_0710_1869.md`](./PAPER_0710_1869.md): concise note on the published
  paper’s claim and why it matters for this repo

## Package Layout

- `repo_v1`: the top-level modules in `quarkConstraints/` provide the
  exploratory repo-local MFV implementation.
- `paper_0710_1869`: frozen paper-facing audit and reproduction helpers with
  strict scope separation from repo-local exploratory work.
- `modern`: versioned artifact/provenance machinery for the eventual
  production-quality modern-bounds lane.

## Current Scope

The current exploratory scaffold includes:

- MFV-native spurion inputs in `model.py`
- exact quark mass-matrix and CKM fitting in `fit.py`
- explicit quark KK-scale conventions in `scales.py`
- mass-basis KK-gluon couplings in `couplings.py`
- a repo-owned `Delta F = 2` exclusion slice in `deltaf2.py`
- deterministic benchmarks and validation helpers
- lightweight proxy and quark-basis misalignment diagnostics
- a scan wrapper and benchmark script

The default target table is a repo-owned bundle (`default_quark_targets()`)
whose PDG quark masses are evolved to `m_t(m_t) = 163.5 GeV`; the `mu = 3 TeV`
scale is the Wilson-coefficient *reference* scale, not the mass-input scale.
`rough_sm_targets()` remains as a compatibility helper for earlier exploratory
naming. The `modern` lane below is the current provenance/versioned-artifact
machinery (not merely "planned").

## Caveats

- the bulk-mass eigenvalue to `c` map is still a repo-local surrogate:
  `BulkMassMap` is an MFV-native numerical device, not yet a uniquely
  faithful derivation of the paper's benchmark relations
- the default benchmark is a deterministic regression point, not a literal
  0710.1869 benchmark reconstruction
- the current fitter does not explore the full spurion parameterization:
  it optimizes singular values and left rotations, while right rotations are
  inherited from the seed/template
- the `Delta F = 2` layer is a fixed-convention v1 exclusion slice with a
  repo-owned input bundle, not a full EFT/RG evolution package
- the low-level helpers still default to the repo bookkeeping convention
  `M_KK ≡ Lambda_IR`, and the benchmark/validation path now keeps that choice
  explicit through `DEFAULT_QUARK_BENCHMARK_XI_KK = 1.0` so it can be updated
  deliberately later
- the `h_RS`-style quantity is still only a proxy; the validation summary now
  reports both the current repo bookkeeping gate (`h_RS < 1`) and the stricter
  paper-scale target (`h_RS <~ 0.3`) instead of conflating them
- the off-diagonal `q`-basis diagnostics are misalignment fractions: smaller
  means more aligned. Some legacy helper outputs still use `alignment_*`
  labels for backward compatibility, but they refer to these misalignment
  fractions rather than an inverse alignment score
