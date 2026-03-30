# Quark-Sector Agent Orchestration Plan

This document is the execution plan for implementing the full quark-sector
extension in this repo with multiple Codex instances under a single
orchestrator.

It is intentionally repo-specific. The goal is not "parallelize everything,"
but to parallelize only where ownership boundaries are clean while enforcing
independent code-logic review, independent physics-logic review, and explicit
test/plot validation.

Related docs:

- [`IMPLEMENTATION_PLAN.md`](IMPLEMENTATION_PLAN.md): physics and feature scope
- [`README.md`](README.md): high-level quark-sector roadmap
- [`scanParams/README.md`](../scanParams/README.md): scan-driver conventions
- [`README.md`](../README.md): repo-wide validation commands and artifact policy

## Orchestrator rules

The orchestrator owns:

- phase ordering and task decomposition
- interface contracts between agents
- merge gates and conflict resolution
- benchmark definitions and acceptance criteria

The orchestrator does not author all code directly. Each code patch must be
checked by agents that did not write it.

Hard rules:

1. No implementation agent reviews its own patch.
2. Every code hunk must receive one independent code-logic review.
3. Every formula-bearing or physics-claim-bearing patch must receive one
   independent physics-logic review.
4. No notebook may contain unique physics logic. Notebooks may only call
   package helpers that are already covered by tests.
5. No phase advances until its required tests, benchmark script, and plot
   checks pass.

## Agent roster

Use separate Codex instances with the following roles.

### Core implementation agents

- `A1 Spurion/Model Agent`
  - Owns MFV-native parameterization and construction of `C_Q`, `C_u`, `C_d`.
  - Primary files:
    - `quarkConstraints/model.py`
    - `quarkConstraints/__init__.py`
    - `tests/test_quark_model.py`

- `A2 Profile/Diagonalization Agent`
  - Owns diagonalization of `C` matrices, derived bulk-mass eigenvalues, and
    overlap extraction using `warpConfig.wavefuncs` and `diagonalization`.
  - Primary files:
    - `quarkConstraints/model.py`
    - `quarkConstraints/benchmarks.py`
    - `tests/test_quark_model.py`
    - `tests/test_quark_benchmarks.py`

- `A3 Fit Agent`
  - Owns full quark mass-matrix construction, exact mass/CKM extraction, and
    fit residual machinery.
  - Primary files:
    - `quarkConstraints/fit.py`
    - `quarkConstraints/benchmarks.py`
    - `tests/test_quark_fit.py`
    - `tests/test_quark_benchmarks.py`

- `A4 Diagnostics Agent`
  - Owns `h_RS`-style proxy summaries, matrix-level alignment diagnostics, and
    mass-basis coupling summaries.
  - Primary files:
    - `quarkConstraints/proxies.py`
    - `tests/test_quark_proxies.py`

- `A5 Scan/CLI Agent`
  - Owns scan wrappers, CSV schema, benchmark scripts, and later integration
    with `scanParams` conventions.
  - Primary files:
    - `quarkConstraints/scan.py`
    - `scripts/benchmark_quark_mfv.py`
    - `scripts/export_quark_validation_figures.py`
    - `tests/test_quark_scan.py`

### Independent reviewer agents

- `R1 Code Reviewer for A1/A2`
- `R2 Code Reviewer for A3`
- `R3 Code Reviewer for A4/A5`

Each reviewer performs line-by-line logic review of diffs and cannot review a
patch they coauthored.

### Physics-logic reviewer agents

- `P1 Spurion/MFV Physics Reviewer`
  - Checks that flavor breaking enters only through the Yukawa spurions and
    that bulk masses are derived rather than treated as fundamental fit
    variables.

- `P2 Flavor-Suppression Physics Reviewer`
  - Checks the `r -> 0` down-sector protection story, basis/misalignment logic,
    and whether diagnostics really probe alignment rather than just overlap
    ratios.

- `P3 Observable/Phenomenology Physics Reviewer`
  - Activated in the later external-observable phase.
  - Checks operator-basis provenance, QCD running assumptions, hadronic-input
    sourcing, and whether exclusion claims are justified.

### Testing and artifact agents

- `T1 Unit/Regression Test Agent`
  - Owns deterministic pytest coverage and benchmark fixtures.

- `T2 Notebook/Plot Validation Agent`
  - Owns executed Jupyter notebooks and figure-generation scripts.
  - Ensures plots are generated from tested helpers and saved under
    `results/figures/`.

- `T3 Integration Agent`
  - Runs repo-wide validation commands, checks artifact cleanliness, and
    verifies that notebooks/scripts still work after rebases.

## File map and ownership boundaries

The quark-sector implementation should stay modular and mirror the existing
repo structure.

Recommended file set:

- `quarkConstraints/model.py`
  - dataclasses for spurion inputs and derived bulk-mass state
  - `C_Q`, `C_u`, `C_d` construction
  - diagonalization helpers
  - overlap extraction using `warpConfig.wavefuncs.f_IR`

- `quarkConstraints/fit.py`
  - exact mass-matrix construction
  - SVD-based diagonalization via `diagonalization.diag.SVD`
  - CKM extraction
  - residual/objective computation

- `quarkConstraints/proxies.py`
  - `h_RS`-style proxy calculations
  - matrix-level alignment diagnostics
  - optional mass-basis coupling summaries

- `quarkConstraints/benchmarks.py`
  - benchmark parameter points
  - benchmark targets for masses, CKM, preferred `r` window, and expected
    suppression patterns

- `quarkConstraints/scan.py`
  - quark-sector scan config dataclass
  - CSV-row construction and reproducibility metadata mirroring `scanParams`
  - optional bridge wrappers if combined scans are added later

- `scripts/benchmark_quark_mfv.py`
  - human-readable benchmark reproduction script, analogous to
    `scripts/benchmark_perez_randall.py`

- `scripts/export_quark_validation_figures.py`
  - deterministic figure-export script for notebook-equivalent plot generation

- `tests/test_quark_model.py`
- `tests/test_quark_fit.py`
- `tests/test_quark_proxies.py`
- `tests/test_quark_benchmarks.py`
- `tests/test_quark_scan.py`

- `quarkConstraints/quark_benchmark_validation.ipynb`
- `quarkConstraints/quark_alignment_validation.ipynb`
- `quarkConstraints/quark_scan_validation.ipynb`

Ownership rules:

- package logic lives in `.py` modules under `quarkConstraints/`
- scripts are thin wrappers around tested package code
- notebooks import package helpers and never define the canonical equations
- later EFT/hadronic machinery should live in new modules, not overload the
  narrow first-pass files

## Review graph

Every task has an author, a code reviewer, and a physics reviewer.

Required review pattern:

- `A1` authored patch -> `R1` code review -> `P1` physics review
- `A2` authored patch -> `R1` code review -> `P1` or `P2` physics review
- `A3` authored patch -> `R2` code review -> `P1` physics review
- `A4` authored patch -> `R3` code review -> `P2` physics review
- `A5` authored patch -> `R3` code review -> `P2` or `P3` physics review,
  depending on scope

If a patch touches both implementation and interpretation, it needs both review
types before merge.

Each reviewer must comment on every changed hunk with one of:

- `accept`
- `change requested`
- `physics concern`
- `test gap`

The orchestrator merges only when all hunks are resolved.

## Phase plan

### Phase 0: Scaffold and contracts

Goal:
- create the interface contracts before heavy implementation begins

Deliverables:
- module skeletons for `model.py`, `fit.py`, `proxies.py`, `benchmarks.py`
- exported API list in `quarkConstraints/__init__.py`
- benchmark-target schema
- initial test-file skeletons

Assigned agents:
- author: `A1`, `A3`, `A4`
- code review: `R1`, `R2`, `R3`
- physics review: `P1`

Exit criteria:
- agreed dataclasses and function signatures
- no ambiguity about fundamental inputs vs derived bulk masses
- benchmark targets defined in one canonical location

### Phase 1: MFV spurions to derived bulk masses

Goal:
- implement the MFV-native model layer

Deliverables:
- spurion input dataclass
- `C_u,d ~ Y_u,d^\dagger Y_u,d`
- `C_Q ~ r Y_u Y_u^\dagger + Y_d Y_d^\dagger`
- diagonalization to derived `c_Q`, `c_u`, `c_d`
- overlap extraction `F_Q`, `F_u`, `F_d`

Assigned agents:
- author: `A1`, `A2`
- code review: `R1`
- physics review: `P1`
- tests: `T1`

Required tests:
- Hermiticity / positive-semidefinite structure checks for `C_u`, `C_d`
- derived eigenvalues are real and reproducibly ordered
- overlap extraction matches `warpConfig.wavefuncs.f_IR`
- basis-rotation invariance tests where appropriate

Physics sign-off questions:
- are bulk masses derived strictly from spurions?
- is the basis/misalignment structure explicit rather than implicit?
- are `c_Q`, `c_u`, `c_d` outputs rather than the fit coordinates?

### Phase 2: Exact quark masses and CKM

Goal:
- reproduce masses and CKM from full mass-matrix diagonalization

Deliverables:
- exact `M_u,d ~= 2 v F_Q Y_u,d F_u,d`
- SVD-based extraction of quark masses and CKM
- residual machinery and deterministic fit driver
- benchmark reproduction script

Assigned agents:
- author: `A3`
- code review: `R2`
- physics review: `P1`
- tests: `T1`, `T3`

Required tests:
- exact diagonalization reproduces benchmark masses and CKM within stated
  tolerances
- CKM matrix is unitary to numerical tolerance
- approximate CKM scaling logic is used only as initialization or sanity check,
  not as the source of truth
- benchmark script prints stable, readable diagnostics

Physics sign-off questions:
- are observables extracted from full matrices rather than overlap-ratio fits?
- does the implementation keep the spurion logic primary?
- is any use of `V_CKM ~ f_Qi / f_Qj` clearly labeled approximate?

### Phase 3: Alignment diagnostics and fast proxies

Goal:
- test the paper's structural suppression claim before adding full hadronic
  phenomenology

Deliverables:
- down-sector `h_RS`-style proxy, explicitly labeled as a proxy
- at least one matrix-level alignment diagnostic
- optional up-sector `D - Dbar` proxy if the fit layer exposes it cheaply
- human-readable diagnostic summary object

Assigned agents:
- author: `A4`
- code review: `R3`
- physics review: `P2`
- tests: `T1`, `T2`

Required tests:
- `r -> 0` should reduce the down-sector alignment/off-diagonal measure
- benchmark window `|r| ~ 0.1-0.4` should show qualitatively suppressed
  down-sector proxy relative to generic larger-`r` points
- diagnostics must remain finite and deterministic across edge cases

Physics sign-off questions:
- do the diagnostics actually probe alignment?
- are proxy claims kept separate from observable claims?
- does the code avoid overstating phenomenological conclusions?

### Phase 4: Scan integration and batch workflows

Goal:
- make the quark-sector machinery scan-ready without breaking existing lepton
  workflows

Deliverables:
- `quarkConstraints/scan.py`
- CSV schema modeled on `scanParams/scan.py`
- reproducibility metadata: seed, git hash, benchmark label, derived fit score
- benchmark and figure-export scripts

Assigned agents:
- author: `A5`
- code review: `R3`
- physics review: `P2`
- tests: `T1`, `T3`

Required tests:
- one-point scan benchmark
- small-grid scan benchmark
- CSV schema regression tests
- script smoke tests

Exit criteria:
- quark scans are deterministic
- output rows carry provenance
- no existing scan tests regress

### Phase 5: Notebook and plot validation

Goal:
- create visual validation artifacts that show the implemented physics behaves
  as expected

Deliverables:
- `quark_benchmark_validation.ipynb`
- `quark_alignment_validation.ipynb`
- `quark_scan_validation.ipynb`
- saved figures in `results/figures/`

Assigned agents:
- author: `T2`
- code review: `R3`
- physics review: `P1`, `P2`
- integration: `T3`

Required plots and expected qualitative behavior:

- overlap/eigenvalue plot
  - derived overlaps should move monotonically with the corresponding bulk-mass
    eigenvalues in the physically relevant ranges

- exact-vs-scaling CKM plot
  - scaling estimates should track hierarchy intuition, but visible deviations
    from the exact diagonalization are acceptable and expected

- `r`-scan suppression plot
  - down-sector proxy and matrix-level alignment measure should decrease as
    `r` approaches zero

- benchmark residual plot
  - quark masses and CKM residuals should cluster within the target tolerance at
    accepted benchmark points

- optional mass-basis coupling heatmaps
  - off-diagonal down-sector structure should visibly collapse in the small-`r`
    regime relative to comparison points

Notebook rules:

- notebook cells may call package functions only
- any plotted array must be reproducible from a tested `.py` helper
- saved figures go under `results/figures/`
- notebook markdown must state whether a panel is exact, approximate, or proxy

### Phase 6: External observables and CP later-wave

Goal:
- add kaon, `B`, `D`, and EDM machinery only after the earlier structural
  phases are stable

Deliverables:
- new quark phenomenology modules for Wilson coefficients, operator basis, RG
  evolution, and observable mapping
- explicit provenance for every external input
- extended physics-review notes

Assigned agents:
- author: `A5` or new dedicated phenomenology agent
- code review: `R3`
- physics review: `P3`
- tests: `T1`, `T3`

Hard gate:
- no exclusion or viability claims tied to current experiments without sourced
  external EFT and hadronic inputs

## Test strategy

The testing hierarchy should be:

1. pure numerical invariants in `tests/test_quark_*.py`
2. benchmark scripts for human-readable reproduction
3. notebooks and figure scripts for visual checks

Do not invert this order.

Numerical tests should cover:

- matrix dimensions and dtype stability
- Hermiticity / symmetry checks
- unitary reconstruction checks for SVD-based diagonalization
- benchmark residual thresholds
- determinism under fixed seeds
- scan-row schema stability
- proxy/alignment monotonic trends in controlled `r` scans

Notebook/plot validation should not replace unit tests. It should expose:

- trends that look physically wrong even when scalar tests pass
- basis-ordering mistakes visible in heatmaps
- broken suppression patterns across the `r` scan
- benchmark points with numerically small but visually suspicious residuals

## Merge gates

No phase is complete until all applicable gates pass.

### Code gate

- independent reviewer approves each changed hunk
- public API matches the phase contract
- no duplicate implementations across scripts/notebooks/modules

### Physics gate

- physics reviewer signs off on assumptions and interpretation
- proxy language is not used as if it were a full observable prediction
- basis and alignment conventions are documented

### Test gate

- `ruff check .`
- `pytest -q`
- benchmark script for the current phase
- notebook/figure generation for the current phase

### Documentation gate

- `quarkConstraints/README.md` stays aligned with actual delivered scope
- new external inputs are called out explicitly
- benchmark conventions are written once and reused

## Execution order

Recommended order for parallel Codex work:

1. `A1` + `A2` implement the model layer while `T1` prepares test scaffolding.
2. `R1` and `P1` review Phase 1 before `A3` starts relying on the API.
3. `A3` implements exact fit logic while `A4` drafts diagnostic interfaces
   against stable fit outputs.
4. `R2`, `P1`, and `T1` clear Phase 2.
5. `A4` completes diagnostics; `T2` starts notebook prototypes from tested
   helpers.
6. `A5` adds scan/scripts only after Phases 1-3 interfaces stop moving.
7. `P3` is activated only when external observable work starts.

This preserves the narrow-first strategy while still allowing meaningful
parallelism.

## Failure modes the orchestrator must watch for

- code drifting back to a generic RS fit in terms of free `c_Q`, `c_u`, `c_d`
- notebooks becoming the only place where formulas are implemented
- proxy metrics being reported as if they were full flavor predictions
- benchmark conventions splitting across multiple files
- scan metadata not recording enough provenance to reproduce a point
- reviewer agents doing style review while missing physics or numerical logic
- physics reviewers approving intuition that is not actually what the code does

## Definition of done

The quark-sector plan is complete only when:

- the MFV-native spurion model is implemented in package code
- derived bulk masses and overlaps are reproducible and tested
- exact quark masses and CKM are extracted from full matrix diagonalization
- fast alignment/proxy diagnostics exist and are clearly labeled
- scan and benchmark tooling are deterministic
- notebook plots visibly support the expected qualitative behavior
- every merged code path has been independently reviewed for code logic and
  physics logic
- later phenomenology claims are held back until the external-input phase is
  actually implemented
