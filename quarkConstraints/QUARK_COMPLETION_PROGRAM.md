# Quark Completion Program

This document is the umbrella execution contract for the post-paper quark
program in this repo.

It sits above the existing quark docs:

- [`AGENT_ORCHESTRATION_PLAN.md`](./AGENT_ORCHESTRATION_PLAN.md):
  exploratory `repo_v1` orchestration
- [`PAPER_READY_REPRODUCTION_PLAN.md`](./PAPER_READY_REPRODUCTION_PLAN.md):
  paper-facing `0710.1869` critical path
- [`PAPER_AGENT_WORKFLOW.md`](./PAPER_AGENT_WORKFLOW.md):
  persistent paper execution playbook
- [`paper_0710_1869/CURRENT_STATUS.md`](./paper_0710_1869/CURRENT_STATUS.md):
  authoritative frozen paper boundary and handoff state

This umbrella doc adds cross-lane program coordination only. It does not
replace, weaken, or relax the stricter lane-specific contracts that already
exist.

Authority order:

1. lane-specific boundary and handoff docs remain authoritative for their lane
2. lane-specific review, test, and acceptance gates remain mandatory
3. this umbrella doc governs only the cross-lane program plan and milestone
   ordering

In particular:

- `repo_v1` work still follows
  [`AGENT_ORCHESTRATION_PLAN.md`](./AGENT_ORCHESTRATION_PLAN.md)
- `paper_0710_1869` work still follows
  [`paper_0710_1869/CURRENT_STATUS.md`](./paper_0710_1869/CURRENT_STATUS.md),
  [`PAPER_READY_REPRODUCTION_PLAN.md`](./PAPER_READY_REPRODUCTION_PLAN.md),
  and [`PAPER_AGENT_WORKFLOW.md`](./PAPER_AGENT_WORKFLOW.md) together

The repo already has an exploratory quark stack and a frozen paper-facing
package. This umbrella doc tracks the broader repo-specific program beyond the
published paper boundary:

- preserve an honest paper reproduction path
- improve the numerics where the current paper bridge is still structural
- add a separate modern-bounds path with versioned external inputs

## Lane Architecture

The quark program has three lanes. They must remain separate.

### `repo_v1`

Location:
- `quarkConstraints/` top-level modules such as `model.py`, `fit.py`,
  `couplings.py`, `deltaf2.py`, and `scan.py`

Purpose:
- fast exploratory MFV implementation
- benchmark prototyping
- proxy diagnostics
- early scan ergonomics

Allowed claims:
- repo-local exploratory behavior
- regression-tested numerical behavior under repo-owned conventions

Not allowed:
- paper-faithful claims by default
- modern exclusion claims by default

### `paper_0710_1869`

Location:
- `quarkConstraints/paper_0710_1869/`

Purpose:
- frozen audit and reproduction lane for the current paper-facing claim
  boundary
- deterministic artifacts and independent verification

Allowed claims:
- only the boundary explicitly frozen in
  [`paper_0710_1869/CURRENT_STATUS.md`](./paper_0710_1869/CURRENT_STATUS.md)

Not allowed:
- silent widening of scope
- silent import of modern inputs
- treating custom surfaces as default exported claims

### `modern`

Location:
- planned canonical production lane under `quarkConstraints/modern/`

Purpose:
- full machinery for the 0710.1869-inspired quark program with modern external
  inputs, modern provenance, and production scans

Target contents:
- paper-faithful MFV model and fit layer
- pointwise KK-gluon couplings and EFT matching
- RG and hadronic bundles with explicit scheme/scale tags
- modern observable and likelihood evaluation
- production scan, artifact, and verifier tooling

Allowed claims:
- only when backed by dated external-input bundles, frozen interpretation
  rules, and passing verifier-backed acceptance tests

## Claim Levels

Every artifact, benchmark, and scan row must declare one claim level.

### `strict_paper`

- reproduce the `0710.1869` model logic under a frozen paper contract
- use paper-lane inputs and conventions only
- do not mix in modern bounds, modern lattice averages, or later-source
  interpretation rules

### `numerically_improved_reproduction`

- stay within the paper-model boundary
- improve the numerical implementation relative to the current structural
  bridge
- allow better fitting, rotation recovery, and pointwise coupling evaluation
- do not advertise modern phenomenology unless the artifact is explicitly
  promoted into `modern_bounds`

### `modern_bounds`

- use the paper-inspired model together with versioned modern external inputs
- require explicit source bundles, provenance ids, and interpretation policy ids
- keep paper-lane goldens untouched

## Hard Separation Rules

1. `paper_0710_1869` is an audit lane, not a staging area for modern work.
2. Paper and modern lanes must use different bundle ids, manifest ids, and
   golden artifacts.
3. Modern inputs must never overwrite paper-era numbers or paper goldens.
4. `repo_v1` proxy fields such as `h_RS` remain diagnostics unless promoted by
   a reviewed milestone into a claim-bearing lane.
5. Shared helpers are allowed only through narrow adapters with tests proving
   that paper exports remain byte-stable.
6. No notebook or plotting script may contain unique physics logic.
7. The verifier must consume exported artifacts only and must not import the
   canonical matching, RG, hadronic, or observable modules it is checking.
8. Every output must carry lane id, claim level, source bundle ids, git SHA,
   and dirty-tree state.

## Milestone Graph

The post-paper program starts here:

- `QS0`: documentation and contract freeze
  - add this umbrella doc
  - correct repo-level status/docs links
  - freeze lane names and claim levels

- `QS1`: seeded structural/reference bridge and frozen physical bridge
  - status in the current working tree: a narrowed noncanonical
    structural/reference slice is implemented under `paper_0710_1869`, and a
    hard-gated QS1 physical bridge now sits beside it under the frozen
    seed-to-profile contract
  - the implemented slices cover seeded `C_Q/C_u/C_d`, ordered eigensystems,
    rotated Yukawas, Eq. (3) diagnostics, a benchmark-reference mass probe,
    and a paper-honest seed-to-physical `c/F` bridge with explicit frozen
    policy ids
  - these slices remain paper-lane only and do not widen the frozen published
    paper boundary
  - remaining broader-program blocker: QS2-QS7 still have to land before the
    repo can make a full allowed-quark-Yukawa claim

- `QS2`: pointwise KK-gluon coupling and matching bridge
  - derive mass-basis flavor couplings from fitted rotations
  - connect fitted points to EFT matching with explicit scale metadata

- `QS3`: honest `strict_paper` reproduction
  - reproduce benchmark structures and paper-style scans under frozen paper
    inputs
  - export deterministic paper artifacts and pass the independent verifier

- `QS4`: modern input registry
  - add versioned experimental, hadronic, CKM, mass, and QCD input bundles
  - keep `strict_paper_inputs` and `modern_default_inputs` separate

- `QS5`: modern `Delta F = 2` phenomenology
  - add `K`, `B_d`, `B_s`, and conservative `D0` observables
  - include `epsilon_K` before making serious modern kaon viability claims

- `QS6`: EDM and CP-odd widening
  - add dipole and EDM machinery before making EDM, CP-odd, or
    "CP-problem-improved" claims
  - keep the strong-CP assumption explicit in docs and outputs
  - `QS6` is not required for a first non-CP `modern_bounds` release
  - a first non-CP `modern_bounds` release may stop after `QS5`, proceed to
    `QS7`, and run production scans, but it must stay explicitly limited to
    non-CP observables and must not advertise EDM or CP-odd conclusions

- `QS7`: production scans and SLURM execution
  - add warm-started fits, caching, deterministic scan schemas, artifacts, and
    verifier-backed scan acceptance
  - run production scans only after the earlier claim-bearing milestones close

- `QS8`: optional second-wave observables
  - rare decays or wider global-likelihood work
  - must not block the first honest modern-bounds release

Dependency rule:
- no milestone may claim closure until the previous claim-bearing milestone has
  clean logic review, physics review, numerical review, and acceptance results

## Review Roles

The orchestrator owns ordering, ownership boundaries, and merge gates. The
orchestrator does not self-review implementation.

Required roles:

- Planner
  - writes slice-specific file lists, contracts, non-goals, and merge order
- Implementation worker
  - edits only assigned files
- Logic reviewer
  - checks control flow, contracts, determinism, and regression risk
- Physics reviewer
  - checks paper traceability, EFT conventions, scheme/scale semantics, and
    claim honesty
- Numerical checker
  - checks deterministic recomputation, tolerance policy, manifests, and
    subprocess stability
- Second-pass reviewer
  - verifies that first-wave findings were actually resolved
- Integration verifier
  - runs benchmark, artifact, and pytest acceptance jobs without authoring the
    physics logic being checked

Minimum review gates:

1. every patch gets one independent logic review
2. every formula-bearing patch gets one independent physics review
3. every benchmark/export/verifier/tolerance patch gets one independent
   numerical review
4. milestone-closing slices get second-pass review before merge

## Safe Parallelization Boundaries

Parallel work is allowed only after upstream interfaces are frozen.

Safe splits:

- docs/contracts vs code implementation
- model/fit vs external-input ingestion once lane/schema ids are frozen
- coupling/matching vs scan/artifact plumbing once data contracts are frozen
- hadronic bundles vs verifier wiring once artifact schemas are frozen
- paper-lane maintenance vs modern-lane implementation when no shared schema is
  being edited

Unsafe splits:

- two workers editing the same lane-level convention or schema file
- concurrent widening of `paper_0710_1869` defaults and `modern` defaults
- verifier workers changing canonical export semantics
- test workers changing physics formulas outside their ownership

## Acceptance Rules

Each claim-bearing lane must have its own:

- benchmark script
- pytest subset
- golden artifacts
- source/provenance manifest
- verifier report

At minimum, `modern_bounds` acceptance must demonstrate:

- paper-lane goldens unchanged
- deterministic artifacts across processes
- explicit source ids and dates for external inputs
- benchmark and scan outputs that can be recomputed from exported artifacts

If this umbrella doc is less specific than a stricter lane-specific workflow,
the lane-specific workflow wins.

## Immediate Program Order

Broader-program status:
- `QS0` is in place
- `QS1` is now present in the working tree as both a narrowed noncanonical
  structural/reference slice and a hard-gated paper-honest physical bridge in
  `paper_0710_1869`
- the remaining blocker is broader-program completion: QS2-QS7 still have to
  land before the repo can claim an allowed quark Yukawa parameter space

Reason:
- the current repo already has strong paper-artifact and verifier discipline,
  and now also has working-tree QS1 structural/reference and physical bridge
  slices
- but it still lacks the QS2-QS7 stack needed to turn that bridge into a full
  allowed-region result with pointwise flavor couplings, modern inputs, and
  production scans
- without those later milestones, neither honest paper reproduction nor
  modern scan work is end-to-end for a full allowed quark Yukawa claim
