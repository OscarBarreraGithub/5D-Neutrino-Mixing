# quarkConstraints

Quark-sector implementation notes for the Randall-Sundrum MFV model of
Fitzpatrick, Perez, and Randall:

- arXiv: [0710.1869](https://arxiv.org/abs/0710.1869)
- PDF: [0710.1869](https://arxiv.org/pdf/0710.1869)
- Local PDF: [`0710.1869v1.pdf`](0710.1869v1.pdf)
- Journal: Phys. Rev. Lett. 100, 171604 (2008)

The paper’s core claim is that 5D minimal flavor violation can collapse the
generic RS flavor/CP problem into a controlled next-to-MFV structure. The
key low-energy dial is the parameter `r`, which suppresses down-sector flavor
violation and can make KK scales near 2 TeV viable.

## Recommended implementation phases

The implementation should be built in stages, with strict scope control.

### Phase 1: Model spec + benchmark reproduction

Implement only the model layer that is directly supported by the paper:

- an MFV-native quark input space built from anarchic spurions `Y_u`, `Y_d`,
  their relative flavor orientation, the overall Yukawa scale, and `r`
- `C_u,d ~ Y_u,d^\dagger Y_u,d`
- `C_Q ~ r Y_u Y_u^\dagger + Y_d Y_d^\dagger`
- derived bulk-mass eigenvalues and quark zero-mode overlaps `F_Q`, `F_u`,
  `F_d`
- quark masses and CKM fit residuals from diagonalizing the full mass matrices,
  not from overlap ratios alone
- reproduction of the paper’s benchmark structure and favored `r ~ 0.1-0.4`

This phase answers whether the repo’s warp-profile machinery can represent the
paper’s MFV/NMFV flavor geometry at all.

### Phase 2: Cheap flavor diagnostics

Add fast diagnostics that follow the paper’s logic but stop short of full
hadronic phenomenology:

- KK-gluon coupling matrices in the quark mass basis
- `h_RS`-style suppression proxies, explicitly treated as proxies
- left-left and left-right `Delta F = 2` proxy summaries
- matrix-level alignment and down-sector versus up-sector misalignment
  diagnostics

This phase answers whether the model really shifts the dangerous flavor tension
out of the down sector, as claimed.

### Phase 3: External-observable integration

Only after phases 1 and 2 are stable should the repo grow a fuller constraint
stack:

- kaon, `B`, and `D` mixing observables with external EFT inputs
- EDM calculations
- scan integration with explicit provenance for non-paper ingredients

## Recommended boundaries

Do not start with a full kaon/EDM package. The first useful version should
test the structural claim of the paper, not pretend to be a complete flavor
phenomenology engine.

The first milestone should be: "can we reproduce the paper’s qualitative claim
that a small `r` suppresses down-sector flavor violation while preserving CKM
fit quality?"

## Current status

The quark-sector MFV scaffold is now present in this repo:

- MFV-native spurion inputs in `model.py`
- exact quark mass-matrix and CKM fitting in `fit.py`
- deterministic benchmarks and validation helpers
- lightweight proxy and alignment diagnostics
- a scan wrapper and benchmark script

This is usable as a repo-local exploratory implementation, not yet as a
paper-level reproduction.

## Current caveats

- the bulk-mass eigenvalue to `c` map is still a repo-local surrogate:
  `BulkMassMap` is an MFV-native numerical device, not yet a uniquely
  faithful derivation of the paper's benchmark relations
- the default benchmark is a deterministic regression point, not a literal
  0710.1869 benchmark reconstruction
- the rough SM-like target set is intended for exploratory fitting only
- the current fitter does not explore the full spurion parameterization:
  it optimizes singular values and left rotations, while right rotations are
  inherited from the seed/template
- the `h_RS`-style quantity is still a proxy and uses the repo convention
  `M_KK ≡ Lambda_IR` by default for bookkeeping and internal comparisons
- that `M_KK ≡ Lambda_IR` choice is a convention, not a literal physical
  identification of the first KK mass; physical first-KK masses remain sector
  dependent and should be supplied explicitly for quantitative comparisons

See [`IMPLEMENTATION_PLAN.md`](IMPLEMENTATION_PLAN.md) for the converged
minimal implementation proposal, and
[`AGENT_ORCHESTRATION_PLAN.md`](AGENT_ORCHESTRATION_PLAN.md) for the multi-agent
execution plan.
