Paper-ready execution plan for a 0710.1869 quark-sector reproduction.

This document is the implementation plan for turning the current exploratory
`quarkConstraints/` stack into a defensible, auditable reproduction path for
Fitzpatrick, Perez, and Randall, arXiv:0710.1869.

It is intentionally separate from `AGENT_ORCHESTRATION_PLAN.md`.

- `AGENT_ORCHESTRATION_PLAN.md` remains the exploratory / repo-v1 orchestration
  plan.
- This file governs a new paper-facing path and assumes multiple independent
  planning and review passes have already agreed on the main architecture.
- `PAPER_AGENT_WORKFLOW.md` is the persistent execution playbook for how paper
  planners, implementation workers, logic reviewers, numerical checkers,
  physics reviewers, and second-pass reviewers should be deployed and how their
  findings are merged.
- `paper_0710_1869/CURRENT_STATUS.md` is the authoritative handoff note for
  the current completed milestone, synced repo state, and next honest slice.

The central conclusion of those planning/review passes is:

1. The existing repo-v1 quark pipeline must remain intact for continuity.
2. A paper-ready path must be introduced as a separate mode and separate set of
   entrypoints.
3. The current surrogate bulk-mass map and proxy `Delta F = 2` stack cannot be
   the canonical paper path.
4. "Use the correct M_KK" is necessary but not sufficient; KK-gluon coupling
   normalization, EFT basis/scheme, RG running, hadronic inputs, and observable
   interpretation must all be frozen as contracts before implementation.

# Paper-Ready Reproduction Plan

## Goal

Build a paper-facing quark-sector workflow that can honestly be described as a
reproduction of the 0710.1869 model assumptions and benchmark logic, with a
deterministic numerical pipeline for:

- MFV spurions to bulk masses
- physical KK-gluon scale conventions
- KK-gluon matching at `mu = m_KK`
- `Delta F = 2` operator running
- hadronic matrix elements
- NP observables for kaon, B, and D mixing
- deterministic paper artifacts and independent verification

This plan is intentionally stricter than the current repo-v1 quark workflow.

## Scope Boundary

The first paper-facing claim should be:

- a paper-ready reproduction of the 0710.1869 quark-sector model under an
  explicit convention bundle
- with KK-gluon-driven `Delta F = 2` observables
- starting with kaon NP-only observables as the first honest milestone
- and NP-only observable interpretation unless and until a full SM+NP stack is
  added and validated

The following are out of scope for the first paper-ready claim unless added
explicitly with their own input bundles and review gates:

- KK electroweak or other non-KK-gluon mediators
- tower-summed effective matching beyond a frozen convention choice
- precision SM+NP `epsilon_K`
- long-distance SM treatment of D mixing
- global-fit interpretations
- NLO claims without validated NLO ingredients

Milestones must be named explicitly:

- Milestone 1: kaon-only NP observable reproduction with full LO RG running and
  independent verification
- Milestone 2: extend the same paper-mode contract through explicit custom
  neutral-meson slices, starting with `B_d` / `B_s` Q1 NP-only surfaces and
  then a separate `D0` Q1 NP-only amplitude surface

## Non-Negotiable Pre-Implementation Decisions

These decisions must be frozen in code, docs, and tests before paper-mode
implementation begins.

### A. Scope and Claim Definition

- Decide whether "paper-ready reproduction" means:
  - strict 0710.1869-era conventions and inputs
  - or the 0710.1869 model with modern external inputs
- Decide whether the mediator content is:
  - KK-gluon only
  - or KK-gluon plus additional RS mediators
- Decide whether the heavy-scale treatment is:
  - first KK gluon only
  - or an effective tower-summed approximation

### B. Scale Contract

Paper mode must represent these objects separately:

- `Lambda_IR_GeV`: geometric IR scale
- `m_g1_GeV`: physical first KK gluon mass
- `mu_match_GeV`: EFT matching scale
- `mu_gs_GeV`: scale used to evaluate the QCD gauge coupling in the KK-gluon
  normalization
- optional `m_KK_eff_GeV`: effective tower scale, only if explicitly adopted
- `xi_g`: mapping factor `m_g1 / Lambda_IR`

The default paper-mode contract should be frozen before any fitting or matching:

- exact meaning of `m_KK` in paper mode
- exact default `xi_g`
- exact default choice of `mu_match`
- exact default choice of `mu_gs`
- exact rule for which mass enters the heavy propagator denominator if
  `mu_match_GeV != m_g1_GeV`

### C. MFV Bulk-Mass Construction

Freeze how paper mode constructs `C_Q`, `C_u`, and `C_d`:

- which spurion terms are included
- which universal terms are included
- which coefficients are fitted vs fixed
- how Eq. (3) is interpreted operationally
- how `V_5KM` is represented in code

Paper mode must not route through the current `BulkMassMap` surrogate.

### D. EFT Contract

Freeze:

- operator basis
- operator normalization
- renormalization scheme
- perturbative order commitment
- Wilson coefficient data structure and tags

Minimum initial standard:

- one named `Delta F = 2` basis
- one named renormalization scheme
- LO QCD running with operator mixing and thresholds

### E. Observable Interpretation Contract

Freeze, per system:

- NP-only vs SM+NP interpretation
- comparison target and budget rule
- hadronic evaluation scale

Initial recommendation:

- Kaons: NP-only `Re(M12^NP)` and `Im(M12^NP)` interpretation
- B systems: NP-only `|M12^NP|` / `Delta M`
- D system: first a custom-input-only NP-only amplitude surface
  `M12_D0^NP` / `Delta m_D0^NP = 2 Re(M12_D0^NP)` under exact alignment, with
  any conservative-bound interpretation deferred to a later separately sourced
  milestone

### F. Input and Provenance Contract

Freeze the required metadata for every input bundle:

- source
- year
- citation
- table/equation reference
- scheme
- scale
- transformations applied
- bundle id
- bundle version

### G. Output and Verification Contract

Freeze:

- paper scan schema name and version
- required output keys and units
- exported intermediate artifacts for independent verification
- deterministic artifact locations and metadata rules
- verifier input format and allowed dependencies

Independent verification must be isolated by contract:

- the verifier consumes exported JSON/CSV artifacts, not live in-memory objects
- the verifier must not import the canonical paper matching, RG, or observable
  modules
- the verifier may share only frozen convention and artifact-schema
  definitions

## Parallel Work Model

Implementation is organized into one short critical path and several parallel
tracks.

### Critical Path

The following sequence is blocking and must be respected:

1. PR-1: mode split
2. PR0: freeze paper contracts
3. PR1: paper benchmark and bulk-mass fidelity
4. PR2a: coupling contract freeze
5. PR2b: coupling implementation
6. PR3: EFT matching
7. PR4a: LO RG running
8. PR5: hadronic bundles and NP observables
9. PR6: deterministic paper artifacts
10. PR7: independent verification

### Parallel Tracks

Once the relevant contracts are frozen, the following tracks can proceed in
parallel.

#### Track A: Mode, Schema, and CI

Owner: `A0 Mode/Schema Agent`

Primary files:

- `quarkConstraints/paper_0710_1869/__init__.py` (new)
- `quarkConstraints/paper_0710_1869/scan.py` (new)
- `quarkConstraints/paper_0710_1869/validation.py` (new)
- `.github/workflows/ci.yml`
- new paper-mode tests

Responsibilities:

- split `repo_v1` from `paper_0710_1869`
- define paper schema ids and versions
- prevent proxy-field leakage into paper mode
- add deterministic CI commands
- enforce separate `repo_v1_regression` and `paper_0710_1869_acceptance` lanes

#### Track B: Paper Benchmark and Model Fidelity

Owner: `A1 Paper Model Agent`

Primary files:

- `quarkConstraints/paper_0710_1869/conventions.py` (new)
- `quarkConstraints/paper_0710_1869/inputs.py` (new)
- `quarkConstraints/paper_0710_1869/model.py` (new)
- `quarkConstraints/paper_0710_1869/fit.py` (new)
- `quarkConstraints/paper_0710_1869/benchmarks.py` (new)
- benchmark tests

Responsibilities:

- remove `BulkMassMap` from paper mode
- implement paper-style `C_x` construction
- encode Table I benchmark data
- validate Eq. (3)
- freeze paper benchmark targets and scale definitions
- keep repo-v1 compatibility through explicit adapters, not shared silent
  branching in the canonical paper path

#### Track C: KK-Scale and Coupling Normalization

Owner: `A2 KK Coupling Agent`

Primary files:

- `quarkConstraints/paper_0710_1869/scales.py` (new)
- `quarkConstraints/paper_0710_1869/couplings.py` (new)
- `quarkConstraints/paper_0710_1869/kkgluon.py` (new)

Responsibilities:

- define `Lambda_IR`, `m_g1`, `mu_match`, optional `m_KK_eff`
- define `mu_gs` and propagator-mass semantics explicitly
- define coupling contract and normalization knobs
- implement KK-gluon overlap normalization
- decide and encode universal subtraction policy
- decide and encode single-mode vs effective-scale policy

#### Track D: EFT Matching and RG

Owner: `A3 EFT/RG Agent`

Primary files:

- `quarkConstraints/paper_0710_1869/eft_deltaf2/operators.py` (new)
- `quarkConstraints/paper_0710_1869/eft_deltaf2/matching_kkgluon.py` (new)
- `quarkConstraints/paper_0710_1869/eft_deltaf2/rg.py` (new)
- `qcd/` integration points

Responsibilities:

- define operator basis and Wilson-vector type
- implement cited KK-gluon matching
- implement LO anomalous dimensions and running
- add literature regression checks
- stage NLO as optional later work

#### Track E: Inputs and Observables

Owner: `A4 Inputs/Observable Agent`

Primary files:

- `quarkConstraints/paper_0710_1869/eft_deltaf2/hadronic.py` (new)
- `quarkConstraints/paper_0710_1869/eft_deltaf2/observables.py` (new)
- `quarkConstraints/paper_0710_1869/inputs/` (new)

Responsibilities:

- define input bundle format and provenance
- add hadronic / experimental bundles
- implement NP-only kaon, B, and D observables
- keep SM+NP explicitly deferred unless separately added

#### Track F: Artifacts and Verification

Owner: `A5 Artifact/Verification Agent`

Primary files:

- `scripts/benchmark_quark_0710_1869.py` (new)
- `scripts/export_quark_0710_1869_figures.py` (new)
- paper artifact tests

Responsibilities:

- produce deterministic benchmark tables
- produce deterministic paper figures
- export intermediate Wilson and observable artifacts
- define numeric golden baselines
- keep the independent verifier on an artifact-only interface

## Phase Breakdown

### PR-1: Mode Split

Goal:

- introduce `repo_v1` and `paper_0710_1869` as distinct paths

Required outcomes:

- no paper-mode code path routes through proxy-only outputs by accident
- existing repo-v1 quark tests remain stable

Exit criteria:

- paper mode has separate entrypoints
- schema version fields are defined
- contract tests fail if paper mode emits repo-v1 proxy outputs
- paper mode imports resolve through `quarkConstraints.paper_0710_1869.*`
  modules in the canonical path

### PR0: Freeze Paper Contracts

Goal:

- freeze every ambiguous physics convention before heavy implementation

Required outcomes:

- paper conventions doc/module exists
- scale objects are defined explicitly
- basis/scheme/observable interpretation are frozen
- provenance bundle format is frozen

Exit criteria:

- no unowned convention remains ambiguous
- review sign-off from physics and execution reviewers

### PR1: Benchmark and Bulk-Mass Fidelity

Goal:

- make the paper benchmark reconstructible

Required outcomes:

- Table I encoded as sourced inputs
- Eq. (3) encoded as a validator
- paper-mode `C_x` construction implemented
- paper-mode overlap checks implemented

Exit criteria:

- benchmark can reproduce paper-style `c` and `f` values within tolerance
- paper-mode tests prove negative `c` values are allowed where needed

### PR2a: Coupling Contract Freeze

Goal:

- freeze KK-gluon coupling conventions before EFT matching

Required outcomes:

- exact meaning of `m_KK` frozen
- coupling contract object defined
- universal subtraction policy frozen
- single-mode vs effective/tower policy frozen
- exact `mu_gs` choice and propagator-mass rule frozen

Exit criteria:

- matching layer can be written against a stable coupling interface

### PR2b: Coupling Implementation

Goal:

- implement paper-mode KK-gluon coupling calculation

Required outcomes:

- dimensionless flavor matrix and normalization layer separated
- scale-convention metadata carried through outputs
- no-FCNC and basis-consistency tests pass

Exit criteria:

- couplings are ready for EFT matching

### PR3: EFT Matching

Goal:

- produce physical Wilson coefficients at `mu_match`

Required outcomes:

- operator basis implementation
- cited matching formulas
- Wilson coefficient objects with tags
- scaling and symmetry tests

Exit criteria:

- `C_i(mu_match)` exported for a benchmark point
- canonical Wilson export includes basis, scheme, scale, and convention tags

### PR4a: LO RG Running

Goal:

- add the minimum honest "full RG running"

Required outcomes:

- LO evolution matrices
- threshold segmentation
- VLL and VRR LO running
- LR-MAP-1, freezing the exact paper O4/O5 scalar LR basis to BMU `Q1^LR` /
  `Q2^LR`, is a separate prerequisite milestone before any honest LR running or
  LR observables claim
- LR-RG-1, activating only the exact paper -> BMU -> LO BMU RG -> paper
  running path for LR Wilsons, is the next milestone after LR-MAP-1
- LR-HAD-1, adding only custom-input LR hadronic contracts and matrix elements
  with no default lattice dataset and no LR observables, is the next milestone
  after LR-RG-1
- LR-OBS-1, adding only a separate custom LR-only kaon observable surface
  under exact LR RG/hadronic alignment while leaving the default/exported
  Q1-only observable path, default numerics, artifacts, and verifier
  unchanged, is the next milestone after LR-HAD-1
- LR-TOTAL-1, adding only a separate custom combined kaon observable surface
  for total NP-only `Q1 + LR` under exact Wilson/custom-Q1/custom-LR
  alignment while leaving the default/exported Q1-only path, artifacts, and
  verifier unchanged, is the next milestone after LR-OBS-1
- BS-Q1-CUSTOM-1, adding only custom-input-only Q1 NP-only `B_d` and `B_s`
  observable surfaces plus custom `B_d` / `B_s` hadronic bundles under exact
  Wilson/hadronic alignment while leaving the kaon default/exported Q1 path,
  artifacts, and verifier unchanged, is the next milestone after LR-TOTAL-1
- D0-Q1-CUSTOM-1, adding only a custom-input-only Q1 NP-only `D0`
  hadronic/observable surface plus the minimal `D0` system/matching helper
  under exact Wilson/hadronic alignment while leaving the kaon
  default/exported Q1 path, artifacts, and verifier unchanged, is the next
  milestone after BS-Q1-CUSTOM-1
- LR-RCHI-FREEZE-1, freezing only the kaon LR `R_chi(mu_had)` semantics and
  provenance at `mu_had = 2.0 GeV` under the exact BV 2004 definition and a
  PDG 2024 `N_L = 4` mass-source chain while leaving default `B4/B5` LR
  freezing, LR observables, artifacts, and verifier unchanged, is the
  completed milestone after D0-Q1-CUSTOM-1
- LR-DEFAULT-HAD-1, freezing sourced default LR `B4/B5` inputs plus a full
  default LR hadronic bundle, is the next honest milestone after
  LR-RCHI-FREEZE-1
- literature-style regression checks

Exit criteria:

- benchmark `C_i(mu_had)` values reproducible in CI

### PR4b: Optional NLO Upgrade

Goal:

- improve precision only after LO is validated

Required outcomes:

- NLO inputs are sourced and scheme-consistent

Exit criteria:

- may be deferred without blocking the first paper-ready claim

### PR5: Hadronic Bundles and NP Observables

Goal:

- compute physical NP observables from the evolved Wilsons

Required outcomes:

- hadronic bundle ids and metadata
- NP-only kaon, B, and D observable functions
- observable regression baseline for at least one benchmark point

Exit criteria:

- one benchmark table with physical observables is reproducible

### PR6: Deterministic Paper Artifacts

Goal:

- generate auditable figures and tables from tested code

Required outcomes:

- paper benchmark script
- figure export script
- metadata-stamped outputs under `results/`

Exit criteria:

- artifact commands run in CI and produce stable outputs

### PR7: Independent Verification

Goal:

- verify at least one key numeric chain without reusing the exact same code path

Required outcomes:

- exported `C_i(mu_match)` and `C_i(mu_had)` intermediates
- standalone verifier computation of at least one kaon quantity
- comparison to golden numeric baselines

Exit criteria:

- independent verifier signs off on at least one end-to-end benchmark
- verifier result is produced from exported artifacts without importing the
  canonical paper matching/RG/observable modules

## Validation and Merge Gates

Every paper-mode merge must satisfy the relevant subset of these gates.

### Contract Gates

- paper mode must not use `xi_KK = 1.0` by accident
- paper mode must not emit repo-v1 proxy keys
- paper mode must not call proxy `deltaf2.py` logic in its canonical path
- basis/scheme/mu tags must match before observable evaluation
- paper mode must not import repo-v1 `model.py`, `couplings.py`, or
  `deltaf2.py` in its canonical path except through explicitly reviewed adapter
  layers

### CI Lanes

- `repo_v1_regression`: protects the current exploratory path and may continue
  to exercise the legacy `xi_KK = 1.0` bookkeeping workflow
- `paper_0710_1869_acceptance`: protects the canonical paper path and must
  never depend on repo-v1 proxy logic or legacy scale defaults

### Test Gates

- `ruff check .`
- `pytest -q`
- literature-style RG regression tests
- benchmark reproduction tests
- end-to-end paper benchmark golden-output tests

### Script Gates

- `repo_v1_regression`:
  `python scripts/benchmark_perez_randall.py`
- `repo_v1_regression`:
  `python scripts/benchmark_quark_mfv.py`
- `repo_v1_regression`:
  `python scripts/plot_quark_deltaf2_r_sweep.py --max-nfev 60 --xi-kk 1.0`
- `paper_0710_1869_acceptance`:
  `python scripts/benchmark_quark_0710_1869.py` once available
- `paper_0710_1869_acceptance`:
  `python scripts/export_quark_0710_1869_figures.py` once available

### Independent Verification Gates

- exported intermediate Wilson data exist
- independent verifier reproduces at least one kaon benchmark quantity
- artifact existence and numeric summary checks pass
- verifier imports are restricted to artifact readers, frozen convention
  definitions, and standalone numeric helpers

## Required New Tests

At minimum, paper mode needs:

- scale-contract tests
- no-proxy-leakage tests
- Table I benchmark tests
- Eq. (3) consistency tests
- coupling normalization tests
- EFT basis/scheme tagging tests
- LO RG magic-number tests
- threshold segmentation tests
- end-to-end benchmark observable regression tests
- verifier artifact-roundtrip tests

## Roles and Review Graph

Implementation roles:

- `A0`: mode/schema/CI
- `A1`: benchmark/model
- `A2`: KK-scale/coupling normalization
- `A3`: EFT/RG
- `A4`: inputs/observables
- `A5`: artifacts/verification

Independent review roles:

- `R1`: code review for `A0`/`A1`
- `R2`: code review for `A2`/`A3`
- `R3`: code review for `A4`/`A5`
- `P1`: MFV/bulk-mass physics review
- `P2`: KK-scale/coupling-normalization physics review
- `P3`: EFT/RG/hadronic/observable physics review
- `V1`: independent numeric verifier

No implementation owner reviews their own patch.

## Shortest Honest Critical Path

If schedule pressure is high, the shortest path to an honest first claim is:

1. mode split
2. freeze paper contracts
3. reproduce benchmark / Table I / Eq. (3)
4. freeze and implement first-KK-gluon coupling normalization
5. implement KK-gluon matching
6. implement LO RG with literature regressions
7. add NP-only kaon observable bundle
8. produce one paper benchmark script
9. independent verifier cross-check of one kaon benchmark quantity from
   exported artifacts only

That is the minimum paper-facing path. Everything else is an extension.

## Explicit Risks

- ambiguous `m_KK` meaning
- wrong KK-gluon normalization
- wrong or under-specified `V_5KM` basis interpretation
- unsourced TeV-scale benchmark targets
- basis/scheme mismatch between RG and hadronic inputs
- accidental reuse of proxy-era fields or logic
- fake "independent verification" that reruns the same code path

These are merge-blocking risks, not cleanup items.

## Implementation Policy

Paper mode must favor:

- explicit named conventions over inferred defaults
- separate modules over silent branching inside proxy-era code
- a dedicated `quarkConstraints.paper_0710_1869` package as the canonical
  paper path, with adapters outward rather than shared-core branching inward
- sourced bundles over repo-owned placeholders
- deterministic artifacts over notebook-only logic
- exported intermediates over opaque end-to-end scripts

If a paper-facing number cannot be traced to:

- a frozen convention object
- a cited input bundle
- a tested code path
- and an independently verifiable artifact

then it is not ready to appear in a reproduction claim.
