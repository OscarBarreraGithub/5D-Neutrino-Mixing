# Paper Agent Workflow

This document is the persistent execution playbook for the dedicated
`paper_0710_1869` path.

It is narrower than
[`AGENT_ORCHESTRATION_PLAN.md`](/Users/oscar/Documents/Research_Code/Randall/5D-Neutrino-Mixing/quarkConstraints/AGENT_ORCHESTRATION_PLAN.md)
and should be used together with
[`PAPER_READY_REPRODUCTION_PLAN.md`](/Users/oscar/Documents/Research_Code/Randall/5D-Neutrino-Mixing/quarkConstraints/PAPER_READY_REPRODUCTION_PLAN.md).

It exists for one reason: after context compaction, a new orchestrator should be
able to resume the paper program without rediscovering the role split, review
graph, theory anchors, or verification gates.

## Scope

This workflow applies only to the paper-facing path under:

- `/Users/oscar/Documents/Research_Code/Randall/5D-Neutrino-Mixing/quarkConstraints/paper_0710_1869`
- `/Users/oscar/Documents/Research_Code/Randall/5D-Neutrino-Mixing/scripts/benchmark_quark_0710_1869.py`
- `/Users/oscar/Documents/Research_Code/Randall/5D-Neutrino-Mixing/tests/test_paper_*`

It does not govern `repo_v1` implementation work except where paper CI must
prove isolation from `repo_v1`.

## Core Rule

The orchestrator does not own the technical decisions inside each slice.
Planning, implementation, and review must be separated.

Hard rules:

1. No implementation worker reviews its own patch.
2. Every code hunk must receive one independent logic review.
3. Every formula-bearing or interpretation-bearing patch must receive one
   independent physics review.
4. Every export or verifier patch must receive one independent numerical or
   consistency review.
5. For PR6 and later, every accepted patch must receive a second-pass review by
   an agent that did not author the code and did not perform the first review.
6. No canonical paper export path may contain numerical verification logic.
7. The verifier must consume exported artifacts only and must not import the
   canonical paper matching, RG, hadronic, or observable modules.

## File Map

The paper-facing modules that reviewers must understand first are:

- `/Users/oscar/Documents/Research_Code/Randall/5D-Neutrino-Mixing/quarkConstraints/paper_0710_1869/artifacts.py`
- `/Users/oscar/Documents/Research_Code/Randall/5D-Neutrino-Mixing/quarkConstraints/paper_0710_1869/verifier.py`
- `/Users/oscar/Documents/Research_Code/Randall/5D-Neutrino-Mixing/quarkConstraints/paper_0710_1869/eft_deltaf2/operators.py`
- `/Users/oscar/Documents/Research_Code/Randall/5D-Neutrino-Mixing/quarkConstraints/paper_0710_1869/eft_deltaf2/matching_kkgluon.py`
- `/Users/oscar/Documents/Research_Code/Randall/5D-Neutrino-Mixing/quarkConstraints/paper_0710_1869/eft_deltaf2/rg.py`
- `/Users/oscar/Documents/Research_Code/Randall/5D-Neutrino-Mixing/quarkConstraints/paper_0710_1869/eft_deltaf2/rg_inputs.py`
- `/Users/oscar/Documents/Research_Code/Randall/5D-Neutrino-Mixing/quarkConstraints/paper_0710_1869/eft_deltaf2/hadronic.py`
- `/Users/oscar/Documents/Research_Code/Randall/5D-Neutrino-Mixing/quarkConstraints/paper_0710_1869/eft_deltaf2/observables.py`
- `/Users/oscar/Documents/Research_Code/Randall/5D-Neutrino-Mixing/scripts/benchmark_quark_0710_1869.py`
- `/Users/oscar/Documents/Research_Code/Randall/5D-Neutrino-Mixing/tests/test_paper_artifact_contracts.py`

## Role Taxonomy

### Orchestrator

Owns:

- task decomposition
- file ownership boundaries
- review assignment
- merge order
- final local verification

Does not own:

- first-draft implementation logic
- first-pass review conclusions
- theory sign-off

### Planner

Owns:

- slice-specific implementation plan
- file list
- contract additions
- merge order
- acceptance targets

Output requirements:

- exact files to change
- interfaces to add
- non-goals
- review dependencies
- open risks, if any

### Implementation Worker

Owns:

- code only inside its assigned files
- no review comments
- no CI gating decisions

Output requirements:

- changed files
- assumptions
- local checks run
- unresolved concerns

### Test/Acceptance Worker

Owns:

- pytest coverage
- benchmark script coverage
- CI lane updates if explicitly assigned

Must not change canonical physics formulas unless that ownership is explicitly
granted.

### Logic Reviewer

Checks:

- control flow
- contract coherence
- determinism
- data duplication vs references
- regression risk

Required output:

- line-anchored findings or explicit no-findings verdict
- each finding labeled as `change requested`, `test gap`, or `accept`

### Numerical Checker

Checks:

- deterministic serialization
- same-process and cross-process stability
- exact recomputation identities
- no silent tolerance inflation
- bundle cross-consistency

Required output:

- concrete invariants checked
- exact commands or scripts used
- explicit statement whether any value drifted

### Physics Reviewer

Checks:

- frozen scope honesty
- operator and Hamiltonian conventions
- formula correctness
- scheme and scale semantics
- provenance sufficiency
- whether claimed exports are enough for independent recomputation

Required output:

- line-anchored findings or explicit no-findings verdict
- source anchor for every nontrivial physics statement

### Second-Pass Reviewer

Checks only after fixes from the first review wave land.

Checks:

- that every earlier finding was actually resolved
- that no new inconsistency was introduced by the fix

### Integration Verifier

Owns only end-to-end validation:

- artifact generation
- independent verifier run
- benchmark acceptance
- full relevant pytest subset

## Deployment Rules

Parallel work is allowed only when file ownership is clean.

### Safe Parallel Split for PR6/PR7

- Artifact Worker A:
  - `quarkConstraints/paper_0710_1869/artifacts.py`
  - new export helper module if needed
  - export script under `scripts/`
  - no edits to `verifier.py`, tests, or CI

- Verifier Worker B:
  - `quarkConstraints/paper_0710_1869/verifier.py`
  - optional verifier-only helper module
  - no edits to canonical export helpers

- Test/Acceptance Worker C:
  - `tests/test_paper_artifact_contracts.py`
  - new artifact/verifier tests
  - `scripts/benchmark_quark_0710_1869.py`
  - `.github/workflows/ci.yml` only if explicitly assigned

Unsafe parallel split:

- two workers editing the same paper EFT module
- verifier worker editing export schemas while artifact worker is still writing
  them
- test worker changing canonical formulas to satisfy tests

## Handoff Protocol

Every worker prompt must specify:

1. exact ownership boundary
2. files the worker may edit
3. files the worker must read first
4. frozen scope and theory constraints
5. what is explicitly out of scope
6. required final report format

Every worker response must include:

- changed files
- assumptions
- checks run
- blockers or risks

The orchestrator must feed review findings back verbatim to the implementation
worker before any new edits are made.

## Review Quality Gates

An acceptable review must:

- reference the actual changed file and line
- describe the failure mode, not just a preference
- explain why the issue matters for paper honesty, correctness, or determinism
- state the expected fix direction

A reviewer must request changes when:

- a formula is ambiguous or mislabeled
- provenance is insufficient for independent recomputation
- determinism is not demonstrated across processes
- export and verifier scopes disagree
- a bundle silently encodes unsupported LR content
- the verifier imports canonical paper physics modules

A second reviewer is mandatory when:

- the patch changes artifact schemas
- the patch changes verifier inputs or acceptance logic
- the patch changes hadronic or observable semantics
- the patch adds or alters physics-facing provenance fields

## LR-MAP-1 Milestone

The next honest LR-facing slice is LR-MAP-1:

- freeze the exact paper O4/O5 scalar LR basis (`Q4_LR` / `Q5_LR` in code)
  to BMU `Q1^LR` / `Q2^LR`
- record the exact operator and coefficient map matrices
- keep `lr_running_activated = false`
- do not add LR running, LR hadronic inputs, or LR observables in this slice

Required deployment order:

1. planner and theory checker in parallel
2. implementation worker and tests/acceptance worker in parallel only after the
   basis map is source-frozen
3. logic checker, numerical checker, and physics checker in parallel
4. one orchestrator handoff of consolidated findings back to the worker
5. one second-pass non-author re-review before closure

## LR-RG-1 Milestone

The next honest LR-facing slice after LR-MAP-1 is LR-RG-1:

- activate LR running only
- keep the LR-MAP-1 map contract frozen
- run the exact algorithm paper -> BMU -> LO BMU RG -> paper
- do not add LR hadronic inputs or LR observables in this slice

Required deployment order:

1. planner and theory checker in parallel
2. implementation worker and tests/acceptance worker in parallel only after
   LR-MAP-1 is confirmed frozen and unchanged
3. logic checker, numerical checker, and physics checker in parallel
4. one orchestrator handoff of consolidated findings back to the worker
5. one second-pass non-author re-review before closure

## LR-HAD-1 Milestone

The next honest LR-facing slice after LR-RG-1 is LR-HAD-1:

- add only custom-input LR hadronic contracts and matrix elements
- require only self-consistency within the declared LR hadronic bundle scheme
  and `mu_had`
- do not add a default lattice dataset
- do not add LR observables
- keep artifacts, the standalone verifier, and default bundles unchanged

Required deployment order:

1. planner and theory checker in parallel
2. implementation worker and tests/acceptance worker in parallel only after
   LR-RG-1 is confirmed frozen and unchanged
3. logic checker, numerical checker, and physics checker in parallel
4. one orchestrator handoff of consolidated findings back to the worker
5. one second-pass non-author re-review before closure

## LR-OBS-1 Milestone

The next honest LR-facing slice after LR-HAD-1 is LR-OBS-1:

- add only a separate custom LR-only kaon observable surface
- keep the default/exported kaon observable API and numerics unchanged
- keep artifacts, the standalone verifier, and default bundles unchanged
- keep combined `Q1+LR` total observables out of scope
- keep `epsilon_K` blocked
- require exact LR RG/hadronic alignment on `system_id`, operator basis id,
  operator normalization id, renormalization scheme id, and
  `mu_had`/evaluation scale
- record that the active LR status id now means custom LR hadronic plus custom
  LR observable active while default/export remains Q1-only

Required deployment order:

1. planner and theory checker in parallel
2. implementation worker and tests/acceptance worker in parallel only after
   LR-HAD-1 is confirmed frozen and unchanged
3. logic checker, numerical checker, and physics checker in parallel
4. one orchestrator handoff of consolidated findings back to the worker
5. one second-pass non-author re-review before closure

## LR-TOTAL-1 Milestone

The next honest LR-facing slice after LR-OBS-1 is LR-TOTAL-1:

- add only a separate custom combined kaon observable surface for total
  NP-only `Q1 + LR`
- require an explicit custom Q1 hadronic bundle plus an explicit custom LR
  hadronic bundle
- require exact alignment across the Wilson snapshot, custom Q1 hadronic
  bundle, and custom LR hadronic bundle on `system_id`, operator basis id,
  operator normalization id, renormalization scheme id, `mu_had`/evaluation
  scale, and Hamiltonian convention id
- keep the separate custom LR-only surface available
- keep the default/exported Q1-only observable API and numerics unchanged
- keep artifacts, the standalone verifier, and default bundles unchanged
- keep sourced default LR inputs and artifact/verifier widening out of scope
- keep `epsilon_K` blocked

Required deployment order:

1. planner and theory checker in parallel
2. implementation worker and tests/acceptance worker in parallel only after
   LR-OBS-1 is confirmed frozen and unchanged
3. logic checker, numerical checker, and physics checker in parallel
4. one orchestrator handoff of consolidated findings back to the worker
5. one second-pass non-author re-review before closure

## BS-Q1-CUSTOM-1 Milestone

The next honest paper-mode slice after LR-TOTAL-1 is BS-Q1-CUSTOM-1:

- add only custom-input-only Q1 NP-only `B_d` and `B_s` observable surfaces
- add only custom-input-only `B_d` and `B_s` Q1 hadronic bundles
- require exact alignment between the Wilson snapshot and each custom
  hadronic bundle on `system_id`, operator basis id, operator normalization id,
  renormalization scheme id, `mu_had`/evaluation scale, and Hamiltonian
  convention id
- do not add `D0` in this slice
- keep the kaon default/exported Q1 path unchanged
- keep the kaon LR-only and custom combined surfaces available
- keep artifacts, the standalone verifier, and default result bundles unchanged
- keep sourced default LR inputs, artifact/verifier widening, and `epsilon_K`
  out of scope

Required deployment order:

1. planner and theory checker in parallel
2. implementation worker and tests/acceptance worker in parallel only after
   LR-TOTAL-1 is confirmed frozen and unchanged
3. logic checker, numerical checker, and physics checker in parallel
4. one orchestrator handoff of consolidated findings back to the worker
5. one second-pass non-author re-review before closure

## Theory Guidance

The physics reviewer for the paper path should use these anchors:

- `0710.1869` for the model-facing RS claim boundary
- `BMU, hep-ph/0005183` for LO ADM and RG conventions already frozen in the
  repo
- `Ciuchini, hep-ph/9711402` for basis and scheme caveats
- PDG/FLAG-derived hadronic references already encoded in
  `/Users/oscar/Documents/Research_Code/Randall/5D-Neutrino-Mixing/quarkConstraints/paper_0710_1869/eft_deltaf2/hadronic.py`

For LR-MAP-1, the reviewer must read first:

- `/Users/oscar/Documents/Research_Code/Randall/5D-Neutrino-Mixing/quarkConstraints/0710.1869v1.pdf`
- `/Users/oscar/Documents/Research_Code/Randall/5D-Neutrino-Mixing/quarkConstraints/PAPER_READY_REPRODUCTION_PLAN.md`
- `/Users/oscar/Documents/Research_Code/Randall/5D-Neutrino-Mixing/quarkConstraints/paper_0710_1869/eft_deltaf2/operators.py`
- `/Users/oscar/Documents/Research_Code/Randall/5D-Neutrino-Mixing/quarkConstraints/paper_0710_1869/eft_deltaf2/matching_kkgluon.py`
- `/Users/oscar/Documents/Research_Code/Randall/5D-Neutrino-Mixing/quarkConstraints/paper_0710_1869/eft_deltaf2/rg.py`
- `/Users/oscar/Documents/Research_Code/Randall/5D-Neutrino-Mixing/quarkConstraints/paper_0710_1869/eft_deltaf2/rg_inputs.py`
- `/Users/oscar/Documents/Research_Code/Randall/5D-Neutrino-Mixing/quarkConstraints/paper_0710_1869/eft_deltaf2/hadronic.py`
- `/Users/oscar/Documents/Research_Code/Randall/5D-Neutrino-Mixing/quarkConstraints/paper_0710_1869/eft_deltaf2/observables.py`

Exact source gate for LR-MAP-1:

- BMU eq. (2.1): `Q1^LR = (\bar s^alpha gamma_mu P_L d^alpha)
  (\bar s^beta gamma^mu P_R d^beta)`,
  `Q2^LR = (\bar s^alpha P_L d^alpha)(\bar s^beta P_R d^beta)`
- BMU eq. (3.12): `(gamma_mu P_L)_{ij} (gamma_mu P_R)_{kl} =
  2 (P_R)_{il} (P_L)_{kj}`
- BMU statement in section 4: the Fierz identity (3.12) remains valid at two
  loops in NDR-MS for VLR/SLR
- current paper operator definitions with explicit projectors:
  `O4 = (\bar q_hi^alpha P_L q_lo^alpha)(\bar q_hi^beta P_R q_lo^beta)`,
  `O5 = (\bar q_hi^alpha P_L q_lo^beta)(\bar q_hi^beta P_R q_lo^alpha)`,
  with `P_L = (1-gamma5)/2`, `P_R = (1+gamma5)/2`

For PR6/PR7, the reviewer must specifically check:

- no hidden `+ h.c.` factor was introduced
- `M12 = (C1_VLL + C1_VRR) * <Q1> / (2 m_K)`
- `Delta_m = 2 Re(M12)`
- `<Q1> = (8/3) f_K^2 m_K^2 B_K(mu_had)`
- kaon-only, NP-only, Q1-supported scope remains explicit
- LR remains guarded and cannot leak into the claimed artifact surface

For LR-MAP-1, the reviewer must specifically check:

- the exact source gate above is sufficient and explicitly cited before any
  coding
- the exact paper O4/O5 scalar LR to BMU `Q1^LR` / `Q2^LR` map is frozen with
  sign, color, normalization, and ordered-basis conventions stated explicitly
- the map matrices are recorded exactly:
  `[O4, O5]^T = [[0, 1], [1/2, 0]] [Q1^LR_BMU, Q2^LR_BMU]^T`,
  `[Q1^LR_BMU, Q2^LR_BMU]^T = [[0, 2], [1, 0]] [O4, O5]^T`,
  `[C1^LR_BMU, C2^LR_BMU]^T = [[0, 1/2], [1, 0]] [C4, C5]^T`
- the code and tests do not claim LR running or LR observables
- `lr_running_activated` remains `false`
- human-facing text uses neutral O4/O5 scalar LR wording; `SUSY` appears only
  as a historical label or compatibility note
- if the source proof is incomplete, LR remains guarded

Primary source proof is mandatory before future LR unguarding. A frozen basis
map alone is not sufficient unless the downstream scheme choice is also cited
and pinned.

For LR-RG-1, the reviewer must specifically check:

- LR-MAP-1 remains frozen and unchanged
- the exact algorithm is recorded and implemented as
  `paper -> BMU -> LO BMU RG -> paper`
- the running step uses the frozen BMU LO LR ADM in the BMU NDR-MS scheme
- the LR map contract remains frozen while `lr_running_activated` becomes
  active
- LR hadronic inputs remain blocked
- LR observables remain blocked
- human-facing text does not imply end-to-end LR observables support

For LR-HAD-1, the reviewer must specifically check:

- LR-MAP-1 and LR-RG-1 remain frozen and unchanged
- the LR hadronic formula gate is exactly BV 2004 eq. (5) for the O4/O5 scalar
  LR basis
- the documented formulas are
  `<O4(mu)> = 2 R_chi(mu) m_K^2 f_K^2 B4(mu)` and
  `<O5(mu)> = (2/3) R_chi(mu) m_K^2 f_K^2 B5(mu)`
- `R_chi(mu) = [m_K / (m_s(mu) + m_d(mu))]^2` is treated as a custom caller
  input with explicit source metadata, not as a hidden default dataset value
- `B4(mu)`, `B5(mu)`, and `R_chi(mu)` each carry source metadata with the same
  declared scheme id and the same declared `mu_had` as the LR hadronic bundle
- LR hadronic support is custom-input-only
- BMU NDR-MS is not required at the LR hadronic-input layer
- Wilson/hadronic scheme alignment is deferred to the later LR-OBS slice
- LR observables remain blocked
- artifacts, verifier behavior, and default bundles remain unchanged

For LR-OBS-1, the reviewer must specifically check:

- LR-MAP-1, LR-RG-1, and LR-HAD-1 remain frozen and unchanged
- the new surface is separate and LR-only rather than a widened default
  observable API
- the default/exported Q1-only kaon observable path and default benchmark
  numerics remain unchanged
- combined `Q1+LR` totals remain blocked
- `epsilon_K` remains blocked
- exact LR RG/hadronic alignment is enforced on `system_id`, operator basis id,
  operator normalization id, renormalization scheme id, and
  `mu_had`/evaluation scale
- the active LR status id explicitly means custom LR hadronic plus custom LR
  observable active while default/export remains Q1-only
- artifacts, verifier behavior, and default bundles remain unchanged

For LR-TOTAL-1, the reviewer must specifically check:

- LR-MAP-1, LR-RG-1, LR-HAD-1, and LR-OBS-1 remain frozen and unchanged
- the new surface is separate and custom combined rather than a widened
  default/exported observable API
- the new surface requires both a custom Q1 hadronic bundle and a custom LR
  hadronic bundle
- exact alignment is enforced across the Wilson snapshot, custom Q1 hadronic
  bundle, and custom LR hadronic bundle on `system_id`, operator basis id,
  operator normalization id, renormalization scheme id, `mu_had`/evaluation
  scale, and Hamiltonian convention id
- the separate custom LR-only surface remains available and unchanged
- the default/exported Q1-only kaon observable path and default benchmark
  numerics remain unchanged
- artifacts, verifier behavior, and default bundles remain unchanged
- sourced default LR inputs remain out of scope
- artifact/verifier widening remains out of scope
- `epsilon_K` remains blocked

For BS-Q1-CUSTOM-1, the reviewer must specifically check:

- LR-MAP-1, LR-RG-1, LR-HAD-1, LR-OBS-1, and LR-TOTAL-1 remain frozen and
  unchanged
- the new surfaces are custom-input-only Q1 NP-only `B_d` / `B_s` rather than
  a widened kaon/exported API
- the new slice adds only custom `B_d` / `B_s` Q1 hadronic bundles and no
  `D0`
- exact alignment is enforced between the Wilson snapshot and each custom
  `B_d` / `B_s` hadronic bundle on `system_id`, operator basis id, operator
  normalization id, renormalization scheme id, `mu_had`/evaluation scale, and
  Hamiltonian convention id
- the kaon default/exported Q1 path remains unchanged
- the kaon LR-only and custom combined surfaces remain available and unchanged
- artifacts, verifier behavior, and default result bundles remain unchanged
- sourced default LR inputs, artifact/verifier widening, and `epsilon_K`
  remain out of scope

## PR6/PR7 Verification Rules

Deliverables for the default kaon benchmark:

- `results/paper_0710_1869/default_kaon/wilsons.json`
- `results/paper_0710_1869/default_kaon/hadronic.json`
- `results/paper_0710_1869/default_kaon/observables.json`
- `results/paper_0710_1869/default_kaon/provenance.json`

Minimum exported content:

- `wilsons.json`:
  - `C_i(mu_match)` and `C_i(mu_had)` or explicit low-scale Wilson set
  - basis, scheme, normalization, scales, benchmark id, source matching tags

- `hadronic.json`:
  - `m_K0`
  - `f_K`
  - `B_K(mu_had)` or explicit conversion inputs
  - `matrix_element_formula_id`
  - `hamiltonian_convention_id`
  - `mu_had_GeV`
  - evaluation scale
  - renormalization scheme
  - provenance citations

- `observables.json`:
  - `M12_K_NP`
  - `Delta_m_K_NP`
  - explicit interpretation and scope ids

- `provenance.json`:
  - all record ids referenced by the other three bundles
  - citations and version strings sufficient to trace every frozen input

Required checks:

1. same-process deterministic export
2. cross-process deterministic export
3. typed artifact round-trip
4. cross-bundle point-id, bundle-id, and source-id consistency
5. hadronic semantics sufficient for independent recomputation
6. verifier import isolation from canonical paper physics modules
7. acceptance benchmark includes artifact export and verifier pass/fail summary

## LR-MAP-1 Verification Rules

Required checks for LR-MAP-1:

1. targeted `ruff check` on changed paper files
2. targeted `pytest -q` for touched `tests/test_paper_*.py`
3. at least one negative test showing LR remains guarded beyond the frozen
   mapping boundary
4. acceptance output still reports LR as unsupported beyond the frozen
   mapping slice
5. explicit physics review citing the source proof for the basis map
6. explicit numerical review confirming no claimed LR numbers were silently
   introduced
7. explicit logic review confirming the exact map matrices and wording rule are
   propagated consistently through docs and metadata

Future LR unguarding is blocked until the repo carries both:

- a frozen basis map into the chosen RG basis
- explicit scheme and source proof for every conversion used downstream

## LR-RG-1 Verification Rules

Required checks for LR-RG-1:

1. targeted `ruff check` on changed paper files
2. targeted `pytest -q` for touched `tests/test_paper_*.py`
3. at least one positive test showing LR Wilsons evolve through the exact
   paper -> BMU -> LO BMU RG -> paper path
4. at least one negative test showing LR hadronic inputs remain blocked
5. at least one negative test showing LR observables remain blocked
6. explicit physics review confirming the BMU LR ADM is the running kernel and
   the map contract stayed frozen
7. explicit numerical review confirming the BMU-basis and paper-basis LR
   evolution agree under the frozen maps
8. explicit logic review confirming only the LR-running status changed while
   the map contract remained frozen

## LR-HAD-1 Verification Rules

Required checks for LR-HAD-1:

1. targeted `ruff check` on changed paper files
2. targeted `pytest -q` for touched `tests/test_paper_*.py`
3. at least one positive test showing custom LR hadronic inputs build with
   explicit source metadata for `B4(mu)`, `B5(mu)`, and `R_chi(mu)`
4. at least one negative test showing missing custom `R_chi(mu)` input is
   rejected
5. at least one negative test showing missing LR source metadata is rejected
6. at least one negative test showing mismatched source-metadata scheme id or
   mismatched `mu_had` across `B4(mu)`, `B5(mu)`, and `R_chi(mu)` is rejected
7. at least one negative test showing LR observables remain blocked
8. explicit physics review confirming the BV 2004 formula gate, the custom
   `R_chi(mu)` input choice, and that BMU NDR-MS is not required at the
   hadronic-input layer
9. explicit numerical review confirming LR matrix elements follow the recorded
   O4/O5 formulas from the supplied custom inputs within the declared bundle
   scheme/scale
10. explicit logic review confirming artifacts, verifier behavior, and default
    bundles remained unchanged and Wilson/hadronic alignment is still deferred
    to LR-OBS

## LR-OBS-1 Verification Rules

Required checks for LR-OBS-1:

1. targeted `ruff check` on changed paper files
2. targeted `pytest -q` for touched `tests/test_paper_*.py`
3. at least one positive test showing the custom LR-only kaon observable path
   succeeds under exact LR RG/hadronic alignment
4. at least one negative test showing the default/exported Q1-only kaon
   observable path remains unchanged
5. at least one negative test showing artifacts and the standalone verifier
   remain unchanged and isolated from the custom LR-only observable path
6. at least one negative test showing combined `Q1+LR` totals remain blocked
7. at least one negative test showing `epsilon_K` remains blocked
8. explicit physics review confirming the exact LR RG/hadronic alignment rule,
   the separate custom LR-only scope boundary, and the unchanged default/export
   boundary
9. explicit numerical review confirming the LR-only observable is computed only
   from the LR block under the aligned LR hadronic inputs and does not drift
   into the default Q1-only observable numerics
10. explicit logic review confirming the active LR status id, API boundaries,
    and isolation rules are propagated consistently through docs, summaries,
    and acceptance outputs

## LR-TOTAL-1 Verification Rules

Required checks for LR-TOTAL-1:

1. targeted `ruff check` on changed paper files
2. targeted `pytest -q` for touched `tests/test_paper_*.py`
3. at least one positive test showing the custom combined kaon observable path
   succeeds only when the Wilson snapshot, custom Q1 hadronic bundle, and
   custom LR hadronic bundle are exactly aligned
4. at least one negative test showing missing custom Q1 or missing custom LR
   hadronic input is rejected
5. at least one negative test showing mismatched `system_id`, basis,
   normalization, renormalization scheme id, `mu_had`/evaluation scale, or
   Hamiltonian convention id is rejected on the custom combined path
6. at least one negative test showing the default/exported Q1-only kaon
   observable path remains unchanged
7. at least one negative test showing the separate custom LR-only surface
   remains available and unchanged
8. at least one negative test showing artifacts and the standalone verifier
   remain unchanged and isolated from the custom combined path
9. at least one negative test showing `epsilon_K` remains blocked
10. explicit physics review confirming the combined scope boundary, the exact
    three-way alignment rule, and that sourced default LR inputs remain out of
    scope
11. explicit numerical review confirming the combined observable is the sum of
    the aligned Q1 and LR NP pieces only and does not drift into the frozen
    default/exported numerics
12. explicit logic review confirming API boundaries, artifact/verifier
    isolation, and continued availability of the separate LR-only surface

## BS-Q1-CUSTOM-1 Verification Rules

Required checks for BS-Q1-CUSTOM-1:

1. targeted `ruff check` on changed paper files
2. targeted `pytest -q` for touched `tests/test_paper_*.py`
3. at least one positive test showing custom-input-only `B_d` and `B_s` Q1
   observable paths succeed only when the Wilson snapshot and hadronic bundle
   are exactly aligned
4. at least one negative test showing missing or mismatched `system_id`,
   basis, normalization, renormalization scheme id, `mu_had`/evaluation scale,
   or Hamiltonian convention id is rejected for `B_d` / `B_s`
5. at least one negative test showing `D0` remains absent or blocked
6. at least one negative test showing the kaon default/exported Q1 path
   remains unchanged
7. at least one negative test showing the kaon LR-only and custom combined
   surfaces remain available and unchanged
8. at least one negative test showing artifacts and the standalone verifier
   remain unchanged and isolated from the new `B_d` / `B_s` custom surfaces
9. at least one negative test showing `epsilon_K` remains blocked
10. explicit physics review confirming the custom-input-only `B_d` / `B_s`
    Q1 NP-only scope, the exact Wilson/hadronic alignment rule, and that
    `D0` plus sourced default LR inputs remain out of scope
11. explicit numerical review confirming each `B_d` / `B_s` observable path
    obeys the frozen NP-only Q1 semantics without drifting into the kaon
    default/exported numerics
12. explicit logic review confirming API boundaries, kaon/export isolation,
    and continued availability of the kaon LR-only and custom combined
    surfaces

## Environment Note

If worker execution is attempted through nested `codex exec` inside this app
environment, first run a non-mutating probe. In some desktop sandbox setups,
the nested worker may fail to use local shell/file tools with
`sandbox_apply: Operation not permitted`.

If that happens:

1. prefer any native app subagent mechanism that shares the current workspace
2. otherwise provide workers the required local context explicitly in their
   prompt
3. if neither is possible, record the limitation in the turn summary before any
   direct fallback work is done

This limitation is environmental, not a paper-path design decision.

## Compact Checklist

For a future orchestrator after compaction:

1. Read
   [`PAPER_READY_REPRODUCTION_PLAN.md`](/Users/oscar/Documents/Research_Code/Randall/5D-Neutrino-Mixing/quarkConstraints/PAPER_READY_REPRODUCTION_PLAN.md)
   and this file first.
2. Confirm the current milestone and frozen scope.
3. Rebuild clean file ownership boundaries before launching workers.
4. Launch planner, implementation, logic review, numerical review, physics
   review, then second-pass review.
5. Feed findings back to the implementation worker verbatim.
6. Run local verification only after all review findings are resolved.
7. Record which agents were used, what they owned, and what they verified.
8. For any LR-facing slice, land LR-MAP-1 first and keep
   `lr_running_activated = false` until the downstream scheme/source proof is
   frozen.
9. LR-RG-1 may change LR running status only; it must not be used to smuggle
   in LR hadronic or LR observable support.
10. LR-HAD-1 may add only custom-input LR hadronic contracts and matrix
    elements; it must not change artifacts, verifier behavior, default bundles,
    or LR observables.
11. LR-OBS-1 may add only the separate custom LR-only kaon observable surface;
    default/exported Q1-only numerics, artifacts, verifier behavior, combined
    `Q1+LR` totals, and `epsilon_K` remain for later slices.
12. LR-TOTAL-1 may add only the separate custom combined kaon observable
    surface for total NP-only `Q1 + LR`; sourced default LR inputs,
    artifact/verifier widening, default/exported Q1-only numerics, and
    `epsilon_K` remain for later slices.
13. BS-Q1-CUSTOM-1 may add only custom-input-only Q1 NP-only `B_d` and `B_s`
    observable surfaces plus custom `B_d` / `B_s` hadronic bundles; `D0`,
    sourced default LR inputs, artifact/verifier widening, and `epsilon_K`
    remain for later slices.
