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

## Program docs

- [`QUARK_COMPLETION_PROGRAM.md`](./QUARK_COMPLETION_PROGRAM.md): umbrella
  post-paper program contract for the three-lane quark architecture
- [`AGENT_ORCHESTRATION_PLAN.md`](./AGENT_ORCHESTRATION_PLAN.md):
  exploratory `repo_v1` orchestration
- [`PAPER_READY_REPRODUCTION_PLAN.md`](./PAPER_READY_REPRODUCTION_PLAN.md):
  paper-facing reproduction critical path
- [`PAPER_AGENT_WORKFLOW.md`](./PAPER_AGENT_WORKFLOW.md): persistent paper
  execution workflow to use together with the paper status and reproduction
  docs, not as a standalone replacement for them
- [`paper_0710_1869/CURRENT_STATUS.md`](./paper_0710_1869/CURRENT_STATUS.md):
  authoritative frozen paper boundary

Current broader-program status:
- the frozen published paper-facing boundary still has `next milestone: none`
  in [`paper_0710_1869/CURRENT_STATUS.md`](./paper_0710_1869/CURRENT_STATUS.md)
- the current working tree now includes a narrowed noncanonical QS1
  structural/reference slice plus a hard-gated QS1 physical bridge in
  `paper_0710_1869`
- that is still not enough for a full allowed quark Yukawa claim: QS2-QS7
  remain before the broader program can claim an allowed region, as tracked in
  [`QUARK_COMPLETION_PROGRAM.md`](./QUARK_COMPLETION_PROGRAM.md)

## Historical exploratory roadmap

This section is retained as context for how the original `repo_v1`
exploratory lane was staged. It is not the current blocker list for the
broader quark program.

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

## Historical exploratory boundary

Do not start with a full kaon/EDM package. The first useful version should
test the structural claim of the paper, not pretend to be a complete flavor
phenomenology engine.

The first milestone should be: "can we reproduce the paper’s qualitative claim
that a small `r` suppresses down-sector flavor violation while preserving CKM
fit quality?"

## Current status

The quark-sector work is active and currently split into lanes:

- `repo_v1`: the top-level `quarkConstraints/` modules are the exploratory,
  repo-local MFV implementation
- `paper_0710_1869`: the frozen paper-facing audit/reproduction lane
- `modern`: the planned canonical production lane defined in
  [`QUARK_COMPLETION_PROGRAM.md`](./QUARK_COMPLETION_PROGRAM.md)

The current working tree also carries one narrower cross-lane bridge under the
paper package:

- a noncanonical seeded structural/reference slice in
  `quarkConstraints.paper_0710_1869`
- this slice provides seeded structural MFV eigensystems and a
  benchmark-reference mass probe with explicit noncanonical provenance
- it is useful for broader-program staging, but it is not part of the frozen
  published paper-facing claim boundary

The current exploratory scaffold includes:

- MFV-native spurion inputs in `model.py`
- exact quark mass-matrix and CKM fitting in `fit.py`
- explicit quark KK-scale conventions in `scales.py`
- mass-basis KK-gluon couplings in `couplings.py`
- a repo-owned `Delta F = 2` exclusion slice in `deltaf2.py`
- deterministic benchmarks and validation helpers
- lightweight proxy and quark-basis misalignment diagnostics
- a scan wrapper and benchmark script

This is usable as a repo-local exploratory implementation, not yet as a
paper-level reproduction.

The default target table is now a repo-owned fixed-scale bundle at
`mu = 3 TeV` (`default_quark_targets()`), while `rough_sm_targets()` remains
as a compatibility helper for the earlier exploratory naming.

## Current caveats

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

See [`IMPLEMENTATION_PLAN.md`](IMPLEMENTATION_PLAN.md) for the converged
minimal implementation proposal, and
[`AGENT_ORCHESTRATION_PLAN.md`](AGENT_ORCHESTRATION_PLAN.md) for the multi-agent
execution plan. For the broader post-paper program, see
[`QUARK_COMPLETION_PROGRAM.md`](./QUARK_COMPLETION_PROGRAM.md).
