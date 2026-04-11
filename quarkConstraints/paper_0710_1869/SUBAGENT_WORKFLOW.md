# paper_0710_1869 Subagent Workflow

This document is the persistent orchestration contract for the paper-facing
`quarkConstraints.paper_0710_1869` path. Use it together with:

- `quarkConstraints/PAPER_READY_REPRODUCTION_PLAN.md`
- `quarkConstraints/PAPER_AGENT_WORKFLOW.md`

Its job is narrower than the higher-level plans: it freezes how the repo should
deploy subagents for PR6, PR7, and later paper-mode slices after context
compaction.

## Scope

This workflow governs only:

- `quarkConstraints/paper_0710_1869/`
- `scripts/benchmark_quark_0710_1869.py`
- `tests/test_paper_*.py`

It does not govern repo-v1 implementation except where paper acceptance checks
must prove isolation from repo-v1.

## Frozen Honest Scope

The current honest paper-mode claim remains narrower than the long-term paper
plan. The workflow below assumes:

- canonical path: `quarkConstraints.paper_0710_1869`
- benchmark scope: default kaon benchmark only
- observable scope: kaon NP-only
- supported end-to-end operator subset: `Q1_VLL`, `Q1_VRR`
- guarded subset: `Q4_LR`, `Q5_LR`
- LR-MAP-1 remains the frozen basis-map prerequisite for all later LR work
- LR-RG-1 remains the LR-running prerequisite for all later LR hadronic work
- LR-HAD-1 remains the custom-input-only LR hadronic prerequisite for LR
  observables
- landed LR-OBS-1 scope: a separate custom LR-only kaon observable surface is
  available under exact LR RG/hadronic alignment, while the default/exported
  kaon observable path, benchmark numerics, artifacts, and standalone verifier
  remain Q1-only
- landed LR-TOTAL-1 scope: a separate custom combined kaon observable surface
  for total NP-only `Q1 + LR` is available under exact
  Wilson/custom-Q1/custom-LR alignment, while the default/exported kaon Q1
  path, artifacts, and standalone verifier remain unchanged
- next narrow physics milestone after LR-TOTAL-1: BS-Q1-CUSTOM-1, adding only
  custom-input-only Q1 NP-only `B_d` and `B_s` observable surfaces plus custom
  Q1 hadronic bundles under exact Wilson/hadronic alignment
- `D0`, sourced default LR inputs, artifact/verifier widening, and
  `epsilon_K` remain later slices
- PR6 boundary: deterministic Wilson, hadronic, observable, and provenance
  artifacts for the default kaon benchmark
- PR7 boundary: JSON-only verifier input and independent recomputation of
  `M12_K^NP` without importing canonical matching, RG, hadronic, or observable
  modules

If a slice cannot satisfy that honest scope, the correct action is to keep the
feature guarded and say so explicitly. Do not widen the claim boundary by
inference.

## Roles And Ownership Boundaries

| Role | Owns | Must not own |
| --- | --- | --- |
| Orchestrator | slice decomposition, file ownership, review assignment, merge order, final verification | first-draft production code, first-pass review sign-off |
| Planner | smallest honest slice, file split, acceptance targets, open risks | production edits |
| Theory checker | source freeze before coding, allowed claim boundary, guard decisions | implementation approval of its own code |
| Implementation worker | one bounded package-code slice under `quarkConstraints/paper_0710_1869/` | review of its own patch |
| PR6 export worker | `artifacts.py` and export-only helper code | `verifier.py`, tests, CI |
| PR7 verifier worker | `verifier.py` and verifier-only helpers | `artifacts.py`, canonical physics modules, tests, CI |
| Tests/acceptance worker | `scripts/benchmark_quark_0710_1869.py`, `tests/test_paper_*.py`, CI only if explicitly assigned | canonical paper physics formulas unless ownership is reassigned |
| Logic checker | diff-level invariant, schema, import, guard, and provenance review | authoring the patch under review |
| Numerical checker | determinism, subprocess parity, scaling checks, independent recomputation | schema authorship while serving as checker |
| Physics checker | source-to-code consistency, scale semantics, operator interpretation, honest-scope review | authoring the patch under review |
| Second-pass reviewer | re-review after fixes for PR6+ or any slice with blocking findings | first-draft implementation |

Hard ownership rules:

- No worker reviews its own patch.
- Do not let two implementation workers edit the same file concurrently.
- The export worker is the schema authority for PR6 until PR6 review closes.
- The verifier worker is never allowed to redefine PR6 schemas from inside
  PR7.

## Deployment Order And Parallelism

### Phase 0: planning and theory freeze

Run in parallel:

- one planner
- one theory checker

Exit only when:

- the slice boundary is explicit
- file ownership is explicit
- the honest claim boundary is explicit
- unresolved theory questions are either frozen or converted into guards

### Phase 1: PR6 artifact export

After Phase 0 closes, run in parallel:

- one PR6 export worker
- one tests/acceptance worker

Safe split for PR6:

- export worker: `quarkConstraints/paper_0710_1869/artifacts.py` and any
  export-only helper under `quarkConstraints/paper_0710_1869/`
- tests/acceptance worker: `tests/test_paper_artifact_contracts.py`,
  related `tests/test_paper_*.py`, and
  `scripts/benchmark_quark_0710_1869.py`

Unsafe split for PR6:

- verifier edits in parallel with unresolved artifact-schema edits
- tests worker rewriting export schema details while the export worker is still
  defining them
- any worker changing canonical matching, RG, hadronic, or observable code just
  to satisfy export tests

### Phase 2: PR6 review

Once the PR6 patch is ready, run in parallel:

- one logic checker
- one numerical checker
- one physics checker

PR6 does not advance to PR7 until:

- schema shape is accepted
- deterministic bundle content is accepted
- honest-scope and provenance claims are accepted

### Phase 3: PR7 standalone verifier

PR7 starts only after PR6 artifact schemas and semantics are frozen by review.
After that freeze, run in parallel:

- one PR7 verifier worker
- one tests/acceptance worker

Safe split for PR7:

- verifier worker: `quarkConstraints/paper_0710_1869/verifier.py` and
  verifier-only helpers
- tests/acceptance worker: verifier contract tests, benchmark-script verifier
  checks, and any paper acceptance updates

Unsafe split for PR7:

- verifier worker editing `artifacts.py`
- export worker and verifier worker changing verifier inputs at the same time
- verifier worker importing canonical paper matching, RG, hadronic, or
  observable modules to "reuse" physics logic

### Phase 4: PR7 review

Once the PR7 patch is ready, run in parallel:

- one logic checker
- one numerical checker
- one physics checker

### Phase 5: fix loop and merge order

If any checker finds a blocker:

1. the orchestrator consolidates the findings into one actionable list
2. the list goes back to the implementation worker or a dedicated fix worker
3. the fix worker reports how each finding was addressed
4. the updated diff goes back to at least one original checker
5. PR6+ slices also receive one additional non-author re-review

Merge order is frozen:

1. merge PR6 export boundary first
2. rebase or restack PR7 on the accepted PR6 contract
3. merge PR7 verifier second

Future paper-mode slices use the same sequence:

1. plan and theory freeze
2. disjoint implementation and tests
3. logic/numerical/physics review
4. fix loop
5. merge only after the gates below pass

### LR-MAP-1 milestone

LR-MAP-1 freezes only the exact paper O4/O5 scalar LR basis
(`Q4_LR` / `Q5_LR` in code) to BMU `Q1^LR` / `Q2^LR`. It does not include LR
running, LR hadronic inputs, or LR observables. `lr_running_activated` must
remain `false` in this slice.

Deployment order:

1. run one planner and one theory checker in parallel
2. if the basis map is source-frozen, run one implementation worker and one
   tests/acceptance worker in parallel
3. run one logic checker, one numerical checker, and one physics checker in
   parallel
4. route findings back through one orchestrator handoff
5. require one second-pass non-author re-review before closure

Safe ownership split:

- implementation worker: `eft_deltaf2/operators.py`,
  `eft_deltaf2/matching_kkgluon.py`, `eft_deltaf2/rg_inputs.py`, and
  `eft_deltaf2/rg.py` only if needed for frozen basis metadata rather than LR
  evolution
- tests/acceptance worker: `tests/test_paper_*.py` and
  `scripts/benchmark_quark_0710_1869.py` only

Unsafe ownership split:

- any worker implementing LR running or LR observables in this slice
- any worker weakening existing LR guards to satisfy tests
- concurrent edits to the same EFT file by two implementation workers

Frozen source gate for LR-MAP-1:

- BMU eq. (2.1): `Q1^LR = (\bar s^alpha gamma_mu P_L d^alpha)
  (\bar s^beta gamma^mu P_R d^beta)`,
  `Q2^LR = (\bar s^alpha P_L d^alpha)(\bar s^beta P_R d^beta)`
- BMU eq. (3.12): `(gamma_mu P_L)_{ij} (gamma_mu P_R)_{kl} =
  2 (P_R)_{il} (P_L)_{kj}`
- BMU section 4 statement: the Fierz identity (3.12) remains valid at two
  loops in NDR-MS for VLR/SLR
- current paper operator definitions with explicit projectors:
  `O4 = (\bar q_hi^alpha P_L q_lo^alpha)(\bar q_hi^beta P_R q_lo^beta)`,
  `O5 = (\bar q_hi^alpha P_L q_lo^beta)(\bar q_hi^beta P_R q_lo^alpha)`,
  with `P_L = (1-gamma5)/2`, `P_R = (1+gamma5)/2`

Frozen map matrices for LR-MAP-1 in ordered operator bases
`[O4, O5]^T` and `[Q1^LR_BMU, Q2^LR_BMU]^T`:

- operator map, BMU to paper:
  `[O4, O5]^T = [[0, 1], [1/2, 0]] [Q1^LR_BMU, Q2^LR_BMU]^T`
- operator map, paper to BMU:
  `[Q1^LR_BMU, Q2^LR_BMU]^T = [[0, 2], [1, 0]] [O4, O5]^T`
- coefficient map, paper to BMU:
  `[C1^LR_BMU, C2^LR_BMU]^T = [[0, 1/2], [1, 0]] [C4, C5]^T`

Wording rule:

- in human-facing text, use `paper O4/O5 scalar LR basis` or
  `O4/O5 scalar LR basis`
- `SUSY` is allowed only as a historical literature label or code-compatibility
  note, not as the default human-facing name for this slice

### LR-RG-1 milestone

LR-RG-1 is the next milestone after LR-MAP-1. It activates LR running on the
paper path using the already-frozen O4/O5 scalar LR to BMU `Q1^LR` / `Q2^LR`
map. It still does not include LR hadronic inputs or LR observables.

Exact algorithm for LR-RG-1:

1. start from paper-basis LR Wilsons in ordered paper basis
   `[C4, C5]^T`
2. map paper LR Wilsons into the frozen BMU LR Wilson ordering using the
   frozen LR-MAP-1 Wilson map
3. run the BMU-basis LR pair with the BMU LO LR ADM in the BMU NDR-MS scheme
4. map the evolved BMU LR Wilsons back into the paper O4/O5 scalar LR basis
5. re-embed the evolved LR block into the public paper Wilson object

Operational rule:

- the LR map contract stays frozen in LR-RG-1
- `mapping_matrix_frozen` remains `true`
- `lr_running_activated` becomes `true` only for the LR-RG-1 slice
- LR hadronic inputs remain blocked
- LR observables remain blocked

Deployment order:

1. run one planner and one theory checker in parallel
2. after confirming LR-MAP-1 is already frozen, run one implementation worker
   and one tests/acceptance worker in parallel
3. run one logic checker, one numerical checker, and one physics checker in
   parallel
4. route findings back through one orchestrator handoff
5. require one second-pass non-author re-review before closure

Safe ownership split:

- implementation worker: `eft_deltaf2/rg.py` and `eft_deltaf2/rg_inputs.py`,
  plus any closely related paper-mode metadata file explicitly assigned for the
  LR-running contract
- tests/acceptance worker: `tests/test_paper_*.py` and
  `scripts/benchmark_quark_0710_1869.py` only

Unsafe ownership split:

- any worker rewriting the frozen LR-MAP-1 map matrices during LR-RG-1
- any worker adding LR hadronic or LR observable claims in the same slice
- concurrent edits to the same RG file by two implementation workers

### LR-HAD-1 milestone

LR-HAD-1 is the next milestone after LR-RG-1. It adds LR hadronic contracts and
LR matrix-element builders for the paper O4/O5 scalar LR basis, but only in a
custom-input mode. It does not add a default lattice dataset, LR observables,
or changes to artifacts, the standalone verifier, or default exported bundles.

Exact boundary for LR-HAD-1:

- custom-input-only LR hadronic contracts and matrix elements
- self-consistent only within the declared LR hadronic bundle scheme and
  `mu_had`
- no default lattice dataset for LR bag parameters or chiral inputs
- no LR observables
- no changes to canonical artifact schemas, the standalone verifier, or default
  benchmark bundles

Formula and source gate for LR-HAD-1:

- BV 2004, eq. (5), for the O4/O5 scalar LR matrix elements:
  `<O4(mu)> = 2 R_chi(mu) m_K^2 f_K^2 B4(mu)`
  and
  `<O5(mu)> = (2/3) R_chi(mu) m_K^2 f_K^2 B5(mu)`
- custom chiral ratio choice:
  `R_chi(mu) = [m_K / (m_s(mu) + m_d(mu))]^2`
- LR-HAD-1 must treat `R_chi(mu)` as a custom positive input supplied by the
  caller
- `B4(mu)`, `B5(mu)`, and `R_chi(mu)` must each carry source metadata with the
  same declared scheme id and the same declared `mu_had` as the LR hadronic
  bundle
- LR-HAD-1 does not require the hadronic-input layer to use BMU NDR-MS
- Wilson/hadronic scheme alignment is deferred to the later LR-OBS slice

Operational rule:

- LR-MAP-1 and LR-RG-1 contracts remain frozen
- LR hadronic support is custom-input-only in this slice
- LR hadronic inputs need only be self-consistent within the declared bundle
  scheme/scale
- LR observables remain blocked
- `artifacts.py`, `verifier.py`, and default bundle builders remain unchanged

Deployment order:

1. run one planner and one theory checker in parallel
2. after confirming LR-RG-1 is already frozen, run one implementation worker
   and one tests/acceptance worker in parallel
3. run one logic checker, one numerical checker, and one physics checker in
   parallel
4. route findings back through one orchestrator handoff
5. require one second-pass non-author re-review before closure

Safe ownership split:

- implementation worker: `eft_deltaf2/hadronic.py` only, plus any narrowly
  scoped paper-mode metadata file explicitly assigned for LR hadronic contracts
- tests/acceptance worker: `tests/test_paper_*.py` and
  `scripts/benchmark_quark_0710_1869.py` only

Unsafe ownership split:

- any worker adding default LR lattice numbers or default LR bundles in this
  slice
- any worker changing `artifacts.py`, `verifier.py`, or default benchmark
  writers to expose LR hadronic content
- any worker adding LR observables in the same slice
- concurrent edits to the same hadronic file by two implementation workers

### LR-OBS-1 milestone

LR-OBS-1 is the next milestone after LR-HAD-1. It adds only a separate custom
LR-only kaon observable surface for the paper O4/O5 scalar LR basis. It does
not widen the default/exported Q1-only observable API, default benchmark
numerics, artifacts, or the standalone verifier. Combined `Q1+LR` totals and
`epsilon_K` remain out of scope.

Exact boundary for LR-OBS-1:

- separate custom LR-only kaon observable surface only
- default/exported kaon observable path remains Q1-only and unchanged
- default benchmark numerics remain unchanged
- artifacts, standalone verifier behavior, and default exported bundles remain
  unchanged
- combined `Q1+LR` total observables remain out of scope
- `epsilon_K` remains blocked

Exact alignment required for LR-OBS-1:

- LR RG wilsons and LR hadronic inputs must agree on `system_id`
- LR RG wilsons and LR hadronic inputs must agree on operator basis id
- LR RG wilsons and LR hadronic inputs must agree on operator normalization id
- LR RG wilsons and LR hadronic inputs must agree on renormalization scheme id
- LR RG wilsons and LR hadronic inputs must agree on `mu_had` and evaluation
  scale

Active LR status rule:

- the active LR status id now means LR running is active, custom LR hadronic
  inputs are active, and the custom LR-only observable surface is active
- the same status id must still encode that the default/exported observable
  path remains Q1-only

Deployment order:

1. run one planner and one theory checker in parallel
2. after confirming LR-HAD-1 is already frozen, run one implementation worker
   and one tests/acceptance worker in parallel
3. run one logic checker, one numerical checker, and one physics checker in
   parallel
4. route findings back through one orchestrator handoff
5. require one second-pass non-author re-review before closure

Safe ownership split:

- implementation worker: `eft_deltaf2/observables.py` only, plus any narrowly
  scoped LR-observable metadata file explicitly assigned for status or contract
  ids
- tests/acceptance worker: `tests/test_paper_*.py` and
  `scripts/benchmark_quark_0710_1869.py` only

Unsafe ownership split:

- any worker widening `artifacts.py`, `verifier.py`, or default benchmark
  writers in this slice
- any worker adding combined `Q1+LR` total observables in this slice
- any worker unblocking `epsilon_K` in this slice
- concurrent edits to the same observable file by two implementation workers

### LR-TOTAL-1 milestone

LR-TOTAL-1 is the next milestone after LR-OBS-1. It adds only a separate
custom combined kaon observable surface for total NP-only `Q1 + LR`. It
requires an explicit custom Q1 hadronic bundle plus an explicit custom LR
hadronic bundle. It does not widen the default/exported Q1-only observable
path, artifacts, standalone verifier, or default result bundles. The separate
custom LR-only surface remains available. `epsilon_K` remains out of scope.

Exact boundary for LR-TOTAL-1:

- separate custom combined kaon observable surface for total NP-only
  `Q1 + LR` only
- explicit custom Q1 hadronic bundle required
- explicit custom LR hadronic bundle required
- separate custom LR-only surface remains available and unchanged
- default/exported kaon observable path remains Q1-only and unchanged
- default/exported result bundles remain unchanged
- artifacts and standalone verifier remain unchanged
- `epsilon_K` remains blocked

Exact alignment required for LR-TOTAL-1:

- Wilson snapshot, custom Q1 hadronic bundle, and custom LR hadronic bundle
  must agree on `system_id`
- Wilson snapshot, custom Q1 hadronic bundle, and custom LR hadronic bundle
  must agree on operator basis id
- Wilson snapshot, custom Q1 hadronic bundle, and custom LR hadronic bundle
  must agree on operator normalization id
- Wilson snapshot, custom Q1 hadronic bundle, and custom LR hadronic bundle
  must agree on renormalization scheme id
- Wilson snapshot, custom Q1 hadronic bundle, and custom LR hadronic bundle
  must agree on `mu_had` and evaluation scale
- Wilson snapshot, custom Q1 hadronic bundle, and custom LR hadronic bundle
  must agree on Hamiltonian convention id

Operational rule:

- LR-TOTAL-1 is custom-only and does not redefine the default/exported Q1-only
  public path
- the separate custom LR-only surface remains available
- sourced default LR inputs remain out of scope in this slice
- artifact/verifier widening remains out of scope in this slice

Deployment order:

1. run one planner and one theory checker in parallel
2. after confirming LR-OBS-1 is already frozen, run one implementation worker
   and one tests/acceptance worker in parallel
3. run one logic checker, one numerical checker, and one physics checker in
   parallel
4. route findings back through one orchestrator handoff
5. require one second-pass non-author re-review before closure

Safe ownership split:

- implementation worker: `eft_deltaf2/observables.py` only, plus any narrowly
  scoped observable metadata file explicitly assigned for custom combined-path
  contract ids
- tests/acceptance worker: `tests/test_paper_*.py` and
  `scripts/benchmark_quark_0710_1869.py` only

Unsafe ownership split:

- any worker widening `artifacts.py`, `verifier.py`, or default benchmark
  writers in this slice
- any worker introducing sourced default LR hadronic inputs in this slice
- any worker removing the separate custom LR-only surface in this slice
- any worker unblocking `epsilon_K` in this slice
- concurrent edits to the same observable file by two implementation workers

### BS-Q1-CUSTOM-1 milestone

BS-Q1-CUSTOM-1 is the next milestone after LR-TOTAL-1. It adds only
custom-input-only Q1 NP-only `B_d` and `B_s` observable surfaces plus
custom-input-only `B_d` and `B_s` Q1 hadronic bundles. It does not add `D0`.
It does not widen the kaon default/exported Q1 path, artifacts, standalone
verifier, or default result bundles. The kaon LR-only and custom combined
surfaces remain available. Sourced default LR inputs, artifact/verifier
widening, and `epsilon_K` remain out of scope.

Exact boundary for BS-Q1-CUSTOM-1:

- custom-input-only Q1 NP-only `B_d` observable surface only
- custom-input-only Q1 NP-only `B_s` observable surface only
- custom-input-only `B_d` Q1 hadronic bundle only
- custom-input-only `B_s` Q1 hadronic bundle only
- no `D0` support in this slice
- kaon default/exported Q1 path remains unchanged
- kaon LR-only and custom combined surfaces remain available and unchanged
- artifacts, standalone verifier behavior, and default result bundles remain
  unchanged
- `epsilon_K` remains blocked

Exact alignment required for BS-Q1-CUSTOM-1:

- Wilson snapshot and each custom `B_d` / `B_s` hadronic bundle must agree on
  `system_id`
- Wilson snapshot and each custom `B_d` / `B_s` hadronic bundle must agree on
  operator basis id
- Wilson snapshot and each custom `B_d` / `B_s` hadronic bundle must agree on
  operator normalization id
- Wilson snapshot and each custom `B_d` / `B_s` hadronic bundle must agree on
  renormalization scheme id
- Wilson snapshot and each custom `B_d` / `B_s` hadronic bundle must agree on
  `mu_had` and evaluation scale
- Wilson snapshot and each custom `B_d` / `B_s` hadronic bundle must agree on
  Hamiltonian convention id

Operational rule:

- BS-Q1-CUSTOM-1 is custom-only and does not redefine the kaon default/exported
  Q1 public path
- kaon LR-only and custom combined surfaces remain available
- `D0` remains out of scope in this slice
- sourced default LR inputs, artifact/verifier widening, and `epsilon_K`
  remain out of scope in this slice

Deployment order:

1. run one planner and one theory checker in parallel
2. after confirming LR-TOTAL-1 is already frozen, run one implementation
   worker and one tests/acceptance worker in parallel
3. run one logic checker, one numerical checker, and one physics checker in
   parallel
4. route findings back through one orchestrator handoff
5. require one second-pass non-author re-review before closure

Safe ownership split:

- implementation worker: `eft_deltaf2/hadronic.py` and
  `eft_deltaf2/observables.py` only, plus any narrowly scoped metadata file
  explicitly assigned for custom `B_d` / `B_s` Q1 contract ids
- tests/acceptance worker: `tests/test_paper_*.py` and
  `scripts/benchmark_quark_0710_1869.py` only

Unsafe ownership split:

- any worker widening `artifacts.py`, `verifier.py`, or default benchmark
  writers in this slice
- any worker introducing `D0` observables or hadronic bundles in this slice
- any worker changing kaon LR-only or custom combined surfaces in this slice
- any worker unblocking `epsilon_K` in this slice
- concurrent edits to the same hadronic or observable file by two
  implementation workers

## Required Review Matrix

Every PR6/PR7 slice, and every future physics-bearing paper slice, must satisfy
the following matrix.

| Lane | Required when | Minimum output |
| --- | --- | --- |
| Implementation | every slice | changed files, assumptions, checks run, guarded items that remain |
| Tests/acceptance | every slice | targeted tests, benchmark/acceptance coverage, negative tests for contract drift when applicable |
| Logic checker | every code-bearing slice | file/function or line-anchored findings, or explicit no-findings with scope inspected |
| Numerical checker | every output-bearing, export-bearing, verifier-bearing, or formula-bearing slice | determinism check, subprocess or cross-run parity, and at least one independent recomputation or scaling check |
| Physics checker | every formula-bearing, interpretation-bearing, export-bearing, or verifier-bearing slice | source anchor plus statement that the slice stays inside the honest claim boundary |

Extra rule for PR6 and later:

- one second-pass non-author re-review is mandatory before closure

## Review Handoff And Re-Review

Reviewer findings are handed back to implementation workers exactly once per
fix round through the orchestrator.

Required handoff format:

1. copy each finding exactly, including file/function/severity
2. remove duplicates without weakening the substance
3. separate blockers from hardening suggestions
4. assign the fix to the implementation worker or a dedicated fix worker
5. require a per-finding resolution note in the worker's return

Re-review triggers:

- any blocker from logic, numerical, or physics review
- any change to exported bundle fields, schema ids, provenance fields, or
  scale semantics
- any change to verifier inputs, verifier formulas, or verifier import surface
- any fix that changes claimed numbers

Closure rule:

- a finding is closed only when the original checker, or a second independent
  checker, confirms the patched diff resolves it

Passing tests alone do not close a review finding.

## Import-Isolation Expectations

### Paper mode

Paper-mode modules must remain isolated from repo-v1 and proxy-era code.

Required behavior:

- import within paper mode through relative imports or
  `quarkConstraints.paper_0710_1869.*`
- do not import or re-export repo-v1 modules such as
  `quarkConstraints.benchmarks`, `quarkConstraints.couplings`,
  `quarkConstraints.deltaf2`, `quarkConstraints.fit`,
  `quarkConstraints.model`, `quarkConstraints.proxies`,
  `quarkConstraints.scan`, `quarkConstraints.scales`, or
  `quarkConstraints.validation`
- do not reference forbidden root exports such as `BulkMassMap`

Existing enforcement surfaces:

- `tests/test_paper_mode_contracts.py`
- `scripts/benchmark_quark_0710_1869.py`

### Independent verifier

The standalone verifier is stricter than normal paper mode.

Required behavior:

- consume exported JSON artifacts only
- import only artifact schema helpers from
  `quarkConstraints.paper_0710_1869.artifacts`
- do not import canonical matching, RG, hadronic, or observable modules
- do not depend on live in-memory canonical objects
- keep verifier logic independently reviewable from the canonical pipeline

Existing enforcement surfaces:

- `tests/test_paper_artifact_contracts.py`
- logic review of `quarkConstraints/paper_0710_1869/verifier.py`

PR7 is not complete until the verifier independently recomputes `M12_K^NP`
from exported JSON inputs while preserving the isolation rules above.

## Theory Checker Read Set

The theory checker should read these files first for PR6/PR7 and any later
slice that touches the same contracts:

- `quarkConstraints/0710.1869v1.pdf`
- `quarkConstraints/PAPER_READY_REPRODUCTION_PLAN.md`
- `quarkConstraints/paper_0710_1869/artifacts.py`
- `quarkConstraints/paper_0710_1869/verifier.py`
- `quarkConstraints/paper_0710_1869/eft_deltaf2/operators.py`
- `quarkConstraints/paper_0710_1869/eft_deltaf2/matching_kkgluon.py`
- `quarkConstraints/paper_0710_1869/eft_deltaf2/rg_inputs.py`
- `quarkConstraints/paper_0710_1869/eft_deltaf2/hadronic.py`
- `quarkConstraints/paper_0710_1869/eft_deltaf2/observables.py`
- `scripts/benchmark_quark_0710_1869.py`
- `tests/test_paper_artifact_contracts.py`
- related paper tests under `tests/test_paper_*.py`

Read-order expectation:

1. source paper and reproduction plan
2. artifact and verifier contracts
3. EFT/operator/hadronic/observable modules that define what the verifier is
   allowed to reproduce independently
4. benchmark script and tests that enforce the contract

For LR-MAP-1, add these files to the read set before any coding starts:

- `quarkConstraints/paper_0710_1869/eft_deltaf2/operators.py`
- `quarkConstraints/paper_0710_1869/eft_deltaf2/matching_kkgluon.py`
- `quarkConstraints/paper_0710_1869/eft_deltaf2/rg.py`
- `quarkConstraints/paper_0710_1869/eft_deltaf2/rg_inputs.py`
- `quarkConstraints/paper_0710_1869/eft_deltaf2/hadronic.py`
- `quarkConstraints/paper_0710_1869/eft_deltaf2/observables.py`

Theory checkpoints for LR-MAP-1:

- confirm the source gate above explicitly before implementation
- freeze the exact O4/O5 scalar LR to BMU `Q1^LR` / `Q2^LR` map, including
  ordered-basis convention and coefficient-map transpose rule
- confirm that the map relies only on BMU eq. (2.1), BMU eq. (3.12), BMU's
  NDR-MS survival statement for VLR/SLR, and the current paper projector-based
  operator definitions
- record the source anchor for that map in code-facing ids or provenance fields
- if the source gate is incomplete, keep LR guarded and stop the slice
- do not claim LR running, LR hadronic inputs, or LR observables in this
  milestone

Theory checkpoints for LR-RG-1:

- confirm LR-MAP-1 is already frozen and unchanged before LR-RG-1 starts
- confirm the algorithm is exactly paper -> BMU -> LO BMU RG -> paper
- confirm the BMU LR ADM is the one already frozen under the BMU NDR-MS scheme
- confirm the LR map contract remains frozen while only the LR-running status
  changes
- keep LR hadronic inputs and LR observables blocked in this milestone

Theory checkpoints for LR-HAD-1:

- confirm LR-MAP-1 and LR-RG-1 are already frozen and unchanged before
  LR-HAD-1 starts
- confirm the LR hadronic formula gate is exactly BV 2004 eq. (5) for the
  paper O4/O5 scalar LR basis
- confirm `R_chi(mu)` is treated as a custom input choice rather than a hidden
  derived default
- confirm LR hadronic inputs are custom-input-only and carry explicit source
  metadata for `B4(mu)`, `B5(mu)`, and `R_chi(mu)`
- confirm the source metadata for `B4(mu)`, `B5(mu)`, and `R_chi(mu)` carries
  the same declared scheme id and the same declared `mu_had` as the LR
  hadronic bundle
- do not require BMU NDR-MS at the LR hadronic-input layer
- defer Wilson/hadronic scheme alignment to the later LR-OBS slice
- keep LR observables blocked in this milestone
- keep artifacts, verifier behavior, and default bundles unchanged in this
  milestone

Theory checkpoints for LR-OBS-1:

- confirm LR-MAP-1, LR-RG-1, and LR-HAD-1 are already frozen and unchanged
  before LR-OBS-1 starts
- confirm LR-OBS-1 adds only the separate custom LR-only kaon observable
  surface and does not widen the default/exported Q1-only path
- confirm the LR observable evaluation uses exact alignment between LR RG
  wilsons and LR hadronic inputs on `system_id`, operator basis id, operator
  normalization id, renormalization scheme id, and `mu_had`/evaluation scale
- confirm the active LR status id explicitly encodes custom LR hadronic active
  plus custom LR observable active while default/export remains Q1-only
- keep combined `Q1+LR` totals and `epsilon_K` blocked in this milestone
- keep artifacts, verifier behavior, and default bundles unchanged in this
  milestone

Theory checkpoints for LR-TOTAL-1:

- confirm LR-MAP-1, LR-RG-1, LR-HAD-1, and LR-OBS-1 are already frozen and
  unchanged before LR-TOTAL-1 starts
- confirm LR-TOTAL-1 adds only the separate custom combined kaon observable
  surface for total NP-only `Q1 + LR`
- confirm LR-TOTAL-1 requires an explicit custom Q1 hadronic bundle plus an
  explicit custom LR hadronic bundle
- confirm the combined observable evaluation uses exact alignment across the
  Wilson snapshot, custom Q1 hadronic bundle, and custom LR hadronic bundle on
  `system_id`, operator basis id, operator normalization id,
  renormalization scheme id, `mu_had`/evaluation scale, and Hamiltonian
  convention id
- confirm the separate custom LR-only surface remains available and unchanged
- keep the default/exported Q1-only observable path, artifacts, verifier
  behavior, and default result bundles unchanged in this milestone
- keep sourced default LR inputs, artifact/verifier widening, and `epsilon_K`
  out of scope in this milestone

Theory checkpoints for BS-Q1-CUSTOM-1:

- confirm LR-MAP-1, LR-RG-1, LR-HAD-1, LR-OBS-1, and LR-TOTAL-1 are already
  frozen and unchanged before BS-Q1-CUSTOM-1 starts
- confirm BS-Q1-CUSTOM-1 adds only custom-input-only Q1 NP-only `B_d` and
  `B_s` observable surfaces plus custom `B_d` / `B_s` Q1 hadronic bundles
- confirm the observable evaluation uses exact alignment between the Wilson
  snapshot and each custom `B_d` / `B_s` hadronic bundle on `system_id`,
  operator basis id, operator normalization id, renormalization scheme id,
  `mu_had`/evaluation scale, and Hamiltonian convention id
- confirm no `D0` support is added in this milestone
- confirm the kaon default/exported Q1 path stays unchanged
- confirm the kaon LR-only and custom combined surfaces remain available and
  unchanged
- keep sourced default LR inputs, artifact/verifier widening, and `epsilon_K`
  out of scope in this milestone

## Merge And Verification Gates

### PR6 artifact export gate

PR6 cannot merge until all of the following are true:

- the export boundary is still limited to the default kaon benchmark
- deterministic Wilson, hadronic, observable, and provenance bundles exist
- bundle schema ids, basis ids, scheme ids, and named scales are consistent
- provenance records are sufficient for the exported bundle references
- logic, numerical, and physics review all pass
- second-pass non-author re-review passes

Minimum verification commands:

- targeted `ruff check` on changed paper files
- targeted `pytest -q` for `tests/test_paper_artifact_contracts.py` and any
  touched `tests/test_paper_*.py`
- `python scripts/benchmark_quark_0710_1869.py --require-package`
- `python scripts/benchmark_quark_0710_1869.py`

### PR7 standalone verifier gate

PR7 cannot merge until all of the following are true:

- verifier inputs are JSON-only artifacts
- verifier imports only `artifacts.py` helpers plus standard-library code
- verifier independently recomputes `M12_K^NP`
- the recomputed value agrees with the exported observable contract to the
  declared precision
- no canonical matching, RG, hadronic, or observable module is imported
- logic, numerical, and physics review all pass
- second-pass non-author re-review passes

Minimum verification commands:

- targeted `ruff check` on changed paper files
- targeted `pytest -q` for `tests/test_paper_artifact_contracts.py` and any
  touched `tests/test_paper_*.py`
- `python scripts/benchmark_quark_0710_1869.py --require-package`
- `python scripts/benchmark_quark_0710_1869.py`

For any formula-bearing or export-bearing slice, run full `pytest -q` before
closure if the targeted tests do not already cover the affected paper surface.

### LR-MAP-1 gate

LR-MAP-1 cannot merge until all of the following are true:

- the slice still excludes LR running and LR observables
- the exact paper O4/O5 scalar LR basis to BMU `Q1^LR` / `Q2^LR` map is frozen
  with explicit source proof
- the exact map matrices above are recorded in code-facing docs or metadata
- `lr_running_activated` still remains `false`
- logic, numerical, and physics review all pass
- second-pass non-author re-review passes

Minimum acceptance criteria:

- targeted `ruff check` on changed paper files
- targeted `pytest -q` for touched `tests/test_paper_*.py`
- at least one negative test proving LR remains guarded beyond the frozen
  mapping boundary
- benchmark or contract summary still reports LR as unsupported beyond the
  frozen mapping slice
- reviewer-facing text uses neutral O4/O5 scalar LR wording unless citing the
  historical literature label directly

Rule for future LR unguarding:

- do not unguard LR running, LR hadronic inputs, or LR observables until the
  repo contains this frozen basis map plus explicit scheme and source proof for
  every downstream conversion

### LR-RG-1 gate

LR-RG-1 cannot merge until all of the following are true:

- LR-MAP-1 remains frozen and unchanged
- the exact running algorithm is recorded as paper -> BMU -> LO BMU RG ->
  paper
- the BMU LR ADM is the running kernel for the LR block
- the LR map contract remains frozen while `lr_running_activated` becomes
  active
- LR hadronic inputs remain blocked
- LR observables remain blocked
- logic, numerical, and physics review all pass
- second-pass non-author re-review passes

Minimum acceptance criteria:

- targeted `ruff check` on changed paper files
- targeted `pytest -q` for touched `tests/test_paper_*.py`
- at least one positive test showing nonzero LR Wilsons evolve through the
  paper -> BMU -> LO BMU RG -> paper path
- at least one negative test showing LR hadronic and LR observable paths remain
  blocked
- benchmark or contract summary reports active LR running but blocked LR
  hadronic and observable scope

### LR-HAD-1 gate

LR-HAD-1 cannot merge until all of the following are true:

- LR-MAP-1 and LR-RG-1 remain frozen and unchanged
- LR hadronic contracts are custom-input-only
- the LR hadronic formula gate is recorded as BV 2004 eq. (5) with custom
  `R_chi(mu)` input choice
- LR hadronic inputs are self-consistent within the declared bundle scheme id
  and `mu_had`
- BMU NDR-MS is not required at the LR hadronic-input layer
- Wilson/hadronic scheme alignment is deferred to LR-OBS
- no default LR lattice dataset is introduced
- LR observables remain blocked
- artifacts, verifier behavior, and default bundles remain unchanged
- logic, numerical, and physics review all pass
- second-pass non-author re-review passes

Minimum acceptance criteria:

- targeted `ruff check` on changed paper files
- targeted `pytest -q` for touched `tests/test_paper_*.py`
- at least one positive test showing custom LR hadronic inputs build
  successfully with explicit source metadata
- at least one negative test showing missing LR source metadata or missing
  custom `R_chi(mu)` input is rejected
- at least one negative test showing mismatched source-metadata scheme id or
  mismatched `mu_had` across `B4(mu)`, `B5(mu)`, and `R_chi(mu)` is rejected
- at least one negative test showing LR observables remain blocked
- benchmark or contract summary still reports unchanged artifacts/verifier and
  unchanged default bundles

### LR-OBS-1 gate

LR-OBS-1 cannot merge until all of the following are true:

- LR-MAP-1, LR-RG-1, and LR-HAD-1 remain frozen and unchanged
- only the separate custom LR-only kaon observable surface is added
- exact LR RG/hadronic alignment is enforced on `system_id`, operator basis
  id, operator normalization id, renormalization scheme id, and
  `mu_had`/evaluation scale
- the active LR status id records custom LR hadronic plus custom LR observable
  active while default/export remains Q1-only
- the default/exported Q1-only observable API and default benchmark numerics
  remain unchanged
- artifacts, verifier behavior, and default bundles remain unchanged
- combined `Q1+LR` totals remain blocked
- `epsilon_K` remains blocked
- logic, numerical, and physics review all pass
- second-pass non-author re-review passes

Minimum acceptance criteria:

- targeted `ruff check` on changed paper files
- targeted `pytest -q` for touched `tests/test_paper_*.py`
- at least one positive test showing the separate custom LR-only observable
  path succeeds under exact LR RG/hadronic alignment
- at least one negative test showing mismatched `system_id`, basis,
  normalization, renormalization scheme id, or `mu_had`/evaluation scale is
  rejected on the LR-only path
- at least one negative test showing the default/exported Q1-only path remains
  unchanged
- at least one negative test showing combined `Q1+LR` totals remain blocked
- at least one negative test showing `epsilon_K` remains blocked
- at least one negative test showing artifacts and the standalone verifier
  remain isolated from the custom LR-only observable path
- benchmark or contract summary still reports active custom LR-only
  observables with unchanged default/exported Q1-only scope

### LR-TOTAL-1 gate

LR-TOTAL-1 cannot merge until all of the following are true:

- LR-MAP-1, LR-RG-1, LR-HAD-1, and LR-OBS-1 remain frozen and unchanged
- only the separate custom combined kaon observable surface for total NP-only
  `Q1 + LR` is added
- the combined path requires both a custom Q1 hadronic bundle and a custom LR
  hadronic bundle
- exact alignment is enforced across the Wilson snapshot, custom Q1 hadronic
  bundle, and custom LR hadronic bundle on `system_id`, operator basis id,
  operator normalization id, renormalization scheme id, `mu_had`/evaluation
  scale, and Hamiltonian convention id
- the separate custom LR-only surface remains available and unchanged
- the default/exported Q1-only observable API and default/exported result
  bundles remain unchanged
- artifacts, verifier behavior, and default bundles remain unchanged
- sourced default LR inputs remain out of scope
- artifact/verifier widening remains out of scope
- `epsilon_K` remains blocked
- logic, numerical, and physics review all pass
- second-pass non-author re-review passes

Minimum acceptance criteria:

- targeted `ruff check` on changed paper files
- targeted `pytest -q` for touched `tests/test_paper_*.py`
- at least one positive test showing the separate custom combined kaon
  observable path succeeds only when the Wilson snapshot, custom Q1 hadronic
  bundle, and custom LR hadronic bundle are exactly aligned
- at least one negative test showing missing custom Q1 or missing custom LR
  hadronic input is rejected
- at least one negative test showing mismatched `system_id`, basis,
  normalization, renormalization scheme id, `mu_had`/evaluation scale, or
  Hamiltonian convention id is rejected on the custom combined path
- at least one negative test showing the default/exported Q1-only path remains
  unchanged
- at least one negative test showing the separate custom LR-only path remains
  available and unchanged
- at least one negative test showing artifacts and the standalone verifier
  remain isolated from the custom combined path
- at least one negative test showing `epsilon_K` remains blocked
- benchmark or contract summary still reports unchanged default/exported
  Q1-only scope and a separate custom combined path only

### BS-Q1-CUSTOM-1 gate

BS-Q1-CUSTOM-1 cannot merge until all of the following are true:

- LR-MAP-1, LR-RG-1, LR-HAD-1, LR-OBS-1, and LR-TOTAL-1 remain frozen and
  unchanged
- only custom-input-only Q1 NP-only `B_d` and `B_s` observable surfaces plus
  custom `B_d` / `B_s` hadronic bundles are added
- exact alignment is enforced between the Wilson snapshot and each custom
  `B_d` / `B_s` hadronic bundle on `system_id`, operator basis id, operator
  normalization id, renormalization scheme id, `mu_had`/evaluation scale, and
  Hamiltonian convention id
- no `D0` support is added
- the kaon default/exported Q1 path remains unchanged
- the kaon LR-only and custom combined surfaces remain available and unchanged
- artifacts, verifier behavior, and default result bundles remain unchanged
- sourced default LR inputs remain out of scope
- artifact/verifier widening remains out of scope
- `epsilon_K` remains blocked
- logic, numerical, and physics review all pass
- second-pass non-author re-review passes

Minimum acceptance criteria:

- targeted `ruff check` on changed paper files
- targeted `pytest -q` for touched `tests/test_paper_*.py`
- at least one positive test showing custom-input-only `B_d` and `B_s` Q1
  observable paths succeed only when the Wilson snapshot and hadronic bundle
  are exactly aligned
- at least one negative test showing missing or mismatched `system_id`, basis,
  normalization, renormalization scheme id, `mu_had`/evaluation scale, or
  Hamiltonian convention id is rejected for `B_d` / `B_s`
- at least one negative test showing `D0` remains absent or blocked
- at least one negative test showing the kaon default/exported Q1 path remains
  unchanged
- at least one negative test showing the kaon LR-only and custom combined
  surfaces remain available and unchanged
- at least one negative test showing artifacts and the standalone verifier
  remain isolated from the new `B_d` / `B_s` custom surfaces
- at least one negative test showing `epsilon_K` remains blocked
- benchmark or contract summary still reports unchanged kaon/exported scope and
  separate custom `B_d` / `B_s` Q1-only surfaces only

## Compact Checklist For Later Milestones

Reuse this checklist for every later paper-mode slice.

1. Freeze the smallest honest slice and the exact claim boundary.
2. Assign disjoint file ownership for implementation and tests.
3. Run planner and theory checker before implementation starts.
4. Keep paper-mode imports isolated from repo-v1 and keep the verifier
   artifact-only.
5. Require implementation, tests/acceptance, logic, numerical, and physics
   lanes.
6. Route reviewer findings back through one consolidated orchestrator handoff.
7. Trigger re-review whenever formulas, exported numbers, schemas, or verifier
   behavior change.
8. Merge only after the relevant artifact/verifier gates pass and one
   second-pass non-author review closes the slice.
9. For any LR-facing slice, land LR-MAP-1 first; no LR unguarding is allowed
   without the frozen O4/O5 to BMU map and downstream scheme/source proof.
10. LR-RG-1 may activate LR running only; LR hadronic and LR observables stay
    blocked until a later separately sourced slice lands.
11. LR-HAD-1 may add only custom-input LR hadronic contracts and matrix
    elements; default bundles, artifacts, verifier behavior, and LR observables
    remain out of scope.
12. LR-OBS-1 may add only the separate custom LR-only kaon observable surface;
    default/exported Q1-only numerics, artifacts, verifier behavior, combined
    `Q1+LR` totals, and `epsilon_K` remain out of scope until later slices.
13. LR-TOTAL-1 may add only the separate custom combined kaon observable
    surface for total NP-only `Q1 + LR`; sourced default LR inputs,
    default/exported Q1-only numerics, artifacts, verifier behavior, and
    `epsilon_K` remain out of scope until later slices.
14. BS-Q1-CUSTOM-1 may add only custom-input-only Q1 NP-only `B_d` and `B_s`
    observable surfaces plus custom `B_d` / `B_s` Q1 hadronic bundles; `D0`,
    sourced default LR inputs, artifact/verifier widening, and `epsilon_K`
    remain out of scope until later slices.
