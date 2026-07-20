# Independent repository audit, status, and research roadmap

**Repository:** `OscarBarreraGithub/5D-Neutrino-Mixing`

**Snapshot audited:** `main` at `972b8a9`, synchronized with `origin/main` on 2026-07-15

**Audit date:** 2026-07-15

**Purpose:** establish what the repository currently computes, what the July Fable/Opus
audit changed, what is still unsafe to quote, what happened in the Yukawa-perturbation
study, and what would turn the project into a defensible Lisa-facing physics paper.

This is an independent read of the code, artifacts, commit history, pull request, CI,
audit ledgers, and research scripts. It does not modify any physics implementation or
relabel any existing result.

---

## 0. Bottom line

The repository contains a serious and unusually broad RS-flavor calculation stack. Its
best-developed pieces are the warped geometry and overlap machinery, quark mass/CKM
fitting, mass-basis KK-gluon matching, the core \(\Delta F=2\) evaluators, the scan
harness, and the explicit distinction between rigorous, proxy, partial, and unevaluated
constraints. The July audit also found and corrected real physics mistakes rather than
only cleaning style: Lane-C RG direction and basis maps, matrix-element normalization,
LFV and top-decay factors, a rare-decay kernel, the \(\epsilon_K\) budget, per-meson RG
stopping scales, QCD threshold constants, the Lane-A Bauer bridge, and several scan
fail-open paths.

It is nevertheless **not ready for a paper-level numerical claim**. The reason is not
one small residual. Five hard gates remain:

1. The committed million-row minimal and custodial production scans are pre-audit
   artifacts. The minimal analysis is explicitly marked `SUPERSEDED (B1 retraction)`;
   the custodial analysis is the same June vintage but lacks an equivalent top-level
   banner. No comparable post-audit production scan exists in the tree.
2. The production quark lane uses a physical \(M_{KK}\) axis but the legacy
   `perturbative_4d_legacy` KK-gluon coupling. The physical RS volume-enhanced coupling
   is implemented but has not been selected, convention-matched to a paper, or rerun.
3. The claimed Lane-A median \(\epsilon_K\) wall is contradictory: commit `2c9989b`
   and the audit ledger say the corrected value moved from about 30 TeV to about
   5.7 TeV, while every headline document still says about 30 TeV.
4. The Lane-C paper surface now accepts Q4/Q5 Wilson coefficients in a default artifact
   that contains only Q1 hadronic inputs, and `compute_kaon_np_observables()` silently
   computes Q1 only. Two unresolved post-merge review threads identify exactly this
   mismatch.
5. GitHub CI is red. Run 502 failed both jobs at lint, so its tests and benchmarks were
   skipped. The reported local “1,806 passed” result is not reproduced by the required
   workflow.

There is also a sixth correctness issue with existing scan classifications: B002 and
B004 use a convention-dependent real SM mixing amplitude even though they were included
in the production allowlist and visibly vetoed scan points. `docs/KNOWN_ISSUES.md` says
they were not scanned; that statement is false.

The right immediate action is therefore **not another large scan**. First freeze the
physics convention, repair the two silent-operator/CP-phase defects, make CI genuinely
green, reproduce exact FPR benchmark inputs, and run a small cross-convention validation
matrix. Only then should the project spend compute on new production scans.

For a Lisa-facing paper, the strongest scientific direction is not “RS still fits if we
assume alignment.” That is already the content of the original FPR construction and
several related warped flavor-protection models. A sharper project is:

> **How large is the physically distinct neighborhood of the FPR flavor-safe solution,
> which combinations of 5D Yukawa deformations destroy it, and can the local response
> geometry recover the protecting flavor symmetry and predict correlated CP signals?**

That turns the existing perturbation idea into a quantitative naturalness and symmetry
diagnostic, directly connected to the 2007 Fitzpatrick–Perez–Randall paper, while making
new predictions rather than only updating a mass floor.

---

## 1. What the repository is actually for

The package name and some top-level metadata still describe “5D neutrino mixing,” but the
active project is now substantially broader and predominantly quark-sector RS flavor.
The repository has three layers:

1. **5D geometry and flavor construction**
   - RS warp geometry, fermion zero-mode wavefunctions, Bessel KK spectrum, and overlap
     integrals;
   - lepton Yukawa inversion, a bulk-neutrino seesaw, PMNS/Takagi utilities, and LFV;
   - quark mass/CKM fitting and mappings from Yukawa spurions to bulk masses.
2. **Phenomenology**
   - KK-gluon mass-basis couplings and \(\Delta F=2\) Wilson matching;
   - LO QCD evolution and hadronic contractions for kaon, \(B_d\), \(B_s\), and charm;
   - minimal and tree-level custodial RS electroweak modules;
   - rare kaon, beauty, charm, top, Higgs, collider, LFV, EDM, and neutrino adapters.
3. **Scanning and presentation**
   - a 103-entry catalog (95 primary and 8 secondary constraints);
   - a production scan harness and strict/inclusive survival bookkeeping;
   - large historical scan artifacts, analysis notebooks, reports, and a web explorer.

The central current physics question is whether warped models can retain TeV-scale KK
states while reproducing the SM flavor data and avoiding the RS flavor/CP problem. The
dominant quark observable is usually \(\epsilon_K\), while the minimal electroweak model
also faces the familiar oblique-parameter problem. Custodial symmetry helps the latter
but does not by itself protect flavor-changing KK-gluon couplings.

### 1.1 The three lanes must remain separate

| Lane | What it is | Current implementation status | What may be quoted now |
|---|---|---|---|
| A — anarchic reproduction | Bauer/Casagrande-style anarchic Yukawa baseline and its distribution of required \(M_{KK}\) | Code was materially corrected in July; historical large runs predate the correction | No current median until a corrected rerun; commit-level exploratory result is about 5.7 TeV under its stated legacy convention |
| B — production “as run” | Simplified fitted MFV surrogate, `C_Q = r C_u + C_d`, `BulkMassMap`, no full \(V_{5KM}\) alignment | Operational, but not the exact FPR model; production scan hardcodes legacy 4D coupling | Historical/diagnostic results only, always labeled with the coupling policy and scan commit |
| C — FPR ideal / paper mode | Intended implementation of arXiv:0710.1869 with Yukawa-controlled bulk masses and down-sector alignment | Quarantined; major July fixes landed, exact Table-I calibration and default alignment remain unresolved; current Q1/LR contract is unsafe | Literature target (~2 TeV) only, not a result of this repository |

The original FPR paper explicitly makes the 5D Yukawas the only flavor-breaking sources,
uses them to control bulk masses, and identifies a limit with no down-sector flavor
violation. It reports that this mechanism can permit a KK scale around 2 TeV. That is a
model statement to reproduce, not a number to attach to Lane B by analogy. See
[Fitzpatrick, Perez, Randall, arXiv:0710.1869](https://arxiv.org/abs/0710.1869).

---

## 2. What happened in the latest commits, pushes, and merge

### 2.1 Repository and branch state

- `main` and `origin/main` both point at `972b8a9`.
- The working tree was clean before this audit document was added.
- The latest two post-audit commits add a collaborator report and VS Code exclusions;
  they do not alter physics.
- The repository has 4,041 tracked files. Of these, 1,284 are under `.orchestration`,
  181 under `flavor_catalog_constraints`, and 192 under `tests`.
- Git currently carries roughly 237 MiB of loose objects plus 54 MiB of packed data.
  This is manageable but signals that generated/research workflow state needs a clearer
  archival boundary.

### 2.2 Pull request 4

PR 4 imported the Fable audit and then applied the Opus/Codex correction campaign. Its
audited range contains 50 commits; the change from the Fable import baseline `0bd133b` to
the PR head `a68dc63` touches 230 files with 6,874 insertions and 2,219 deletions.

The local history contains those commits directly rather than a visible merge commit for
PR 4. This is consistent with a rebase/fast-forward-style integration and explains why
`git log --merges` does not show the PR even though GitHub records it as merged.

The campaign process was:

1. Fable produced a compendium and 19 module reports.
2. The reports were imported without code changes in `0bd133b`.
3. Findings were triaged in `FIX_LEDGER.md`.
4. Corrections were applied in serial groups and locally reviewed by Claude/Opus.
5. About 50 lower-severity findings were documented rather than implemented.

This was a valuable audit. The main concern is that “verified” in the ledger usually means
the correction passed a targeted local review and self-contained test set. It does not
mean the current main branch has a green end-to-end CI run or a post-fix production scan.

### 2.3 Important corrections that landed

The following are substantive and should be preserved:

- **Lane C:** corrected RG direction/transpose, LR basis sign, Q1/LR matrix-element
  normalization, localization orientation, KK volume factor, and C4 capability.
- **LFV/top/rare decays:** corrected \(\mu\to3e\) interference, tau branching-fraction
  conversion, top dipole widths, the \(B\to K^*\ell\ell\) kernel, and an
  \(R(D^{(*)})\) dimensional error.
- **\(\epsilon_K\) and \(\Delta F=2\):** unified and sign-aware budget policy,
  per-system RG stopping scales, B-mixing budgets, and explicit LO caveats.
- **QCD:** corrected decoupling and threshold behavior, \(m_t(m_t)\), top uncertainty,
  and shared \(\alpha_s(M_Z)\).
- **Fitting/scanning:** fixed a non-gauge-equivalent reported seed, double application of
  `overall_scale`, fail-open evaluated partials, and malformed HARD extras.
- **Lane A:** fixed the Bauer/repository mass bridge, restored the Wolfenstein \(A\)
  factor, and replaced salted hash seeds in the main reproduction path.
- **RS electroweak:** fixed the Higgs counterpart of the old Casagrande \(Zb\bar b\)
  convention error and reconciled the code-derived 15.96 TeV oblique proxy crossing.
- **Bessel spectrum:** repaired skipped/misordered higher roots.

### 2.4 What the merge did not establish

- It did not choose the physical KK-gluon coupling convention for production.
- It did not rerun the large scans after the audit.
- It did not make Lane B identical to the FPR model.
- It did not complete the Lane-C Table-I/source calibration.
- It did not upgrade \(\Delta F=2\) running to a scheme-consistent NLO calculation.
- It did not implement a full collider recast, full custodial fermion sector, or live
  lepton production scan.
- It did not make GitHub CI green.

---

## 3. Current component and process status

| Component/process | Status | What is in place | What is still needed before a paper claim |
|---|---|---|---|
| RS geometry and zero modes | Strong core | Stable overlap and wavefunction APIs; Bessel spectrum audited | Resolve the \(c=0.5\) numerical plateau convention and derive every normalization used by mass bridges |
| Lepton/neutrino core | Research-grade, incomplete | Seesaw, Yukawa inversion, PMNS/Takagi, LFV plumbing | Replace nonperturbative flagship point, pin public NuFIT values, document the seesaw prefactor, decide massless-neutrino naturalness policy, run a live lepton scan |
| Quark masses and CKM | Numerically strong surrogate | Fits converge and serialize diagnostics; gauge-equivalent seed reporting fixed | Update target covariance/vintage; separate “exact FPR” from `BulkMassMap`; use a likelihood rather than broad factor gates |
| Lane A anarchy | Invalidated pending rerun | Mature scan/reproduction scripts and rich historical artifacts | Rerun after M-11/M-12 and the coupling decision; regenerate all percentiles and figures |
| Lane B production | Operational but model-ambiguous | Production harness, fitted locus, catalog adapters | Freeze physical coupling and exact model definition; rerun; stop calling it FPR/MFV without the surrogate qualifier |
| Lane C FPR | Quarantined/blocking defect | Large paper-mode API, matching/RG/hadronic artifact system, many targeted tests | Fix the Q1/LR contract, vendor/extract the source paper, reproduce Table I and Eq. (3), choose the RH-down alignment, add independent benchmark oracles |
| Core \(\Delta F=2\) | Good LO research code | Mass-basis matching, sign-aware \(\epsilon_K\), per-meson RG endpoints | NLO scheme-consistent running and lattice inputs, kaon chiral scale consistency, investigate VIA additive terms, include \(\Delta m_K\) deliberately |
| Rare B/K/charm | Mixed | Many adapters, modern anchors, several independent hand checks | Several remain proxy/partial; complete sign audit, charm likelihood, nonleptonic and LFV kaon projections |
| Minimal RS electroweak | Proxy-level | \(S,T,U\), \(Zb\bar b\), Higgs/top plumbing | Replace EW001 analytic proxy with a validated global fit or clearly call it a projection; make geometry determine \(L\) |
| Custodial branch | Tree-level proxy | Diagonal \(P_{LR}\) protection and strict/proxy tagging | Full representations, top partners and loops, custodial FCNC derivation, physical partner spectrum |
| Collider constraints | Mostly proxy | Catalog entries and mass-cut logic | Model-dependent production rates, widths, branching ratios, efficiencies, and actual recasts |
| Catalog | Broad but uneven | 103 registered constraints and typed results | Only 46 were in the quark scan and none were charged-lepton constraints; classify release-quality subset and freeze source snapshots |
| Scan harness | Strong engineering, stale physics artifacts | Sharding, serialization, strict/inclusive semantics, finite-statistics helpers | Atomic output handling, fail-closed missing fit diagnostics, manifest every convention, rerun after fixes |
| Existing 1M+1M scans | Superseded/pre-audit | 1,000,000 rows each; 878,707/878,705 evaluated; common \(r\times M_{KK}\) grid | The minimal scan is explicitly pre-B1; both predate the July audit. Do not use their floors. Add a custodial stale-artifact banner and replace both with post-audit versioned scans |
| Scan Explorer/website | Useful diagnostic | Committed page and generated summary | Its inputs are stale; do not deploy as current physics until rebuilt from release artifacts |
| Tests | Large but not release-grade | About 1,800 tests claimed locally; many targeted physics checks | CI is red; reduce circular/golden tests, add external oracles and end-to-end scan checks, pin test environment |
| Reproducibility | Partial | Some manifests, checksums, source notes, and committed reduced data | Runtime extras are incomplete; raw results are mostly ignored/local; no immutable release bundle reproduces headline tables |
| Yukawa response study | Promising proof of concept | Noise scan, Im\(C_4\) gradient, ensemble PCA | It does not yet measure the fully fitted allowed set; fix metrics, gauge redundancies, refitting, full residual Jacobian, curvature and uncertainty |

### 3.1 What was actually scanned

The paired quark-only scans used:

- \(r = 0.05, 0.1, 0.25, 0.5, 1\);
- \(M_{KK} = 1,2,3,5,7,10,15,20,30,50\) TeV;
- 20,000 rows per grid cell, or 1,000,000 rows per model;
- a 46-constraint allowlist out of the 103-entry registry;
- no charged-lepton constraints and a dropped lepton sector.

Their raw runs date to June. The July modification time on the minimal analysis is the
addition of the superseded banner, not a July physics rerun. The current tree contains no
post-PR-4 production output of comparable scope.

---

## 4. Findings and required changes

### P0-1 — There is no current publication-valid production scan

**Evidence**

- Both million-row analysis reports describe June runs.
- The minimal report begins with `SUPERSEDED (B1 retraction)` and explains that its old
  25–30 TeV \(Zb\bar b\) floor came from a retracted bug.
- The audit subsequently changed Lane A normalization, \(\epsilon_K\) budgets,
  per-system RG, QCD constants, fail-closed scan behavior, LFV conventions, and catalog
  predictions.

**Impact**

Historical plots can validate pipeline mechanics and illustrate why results changed,
but they cannot establish the present floor or survival fractions. Any document that
calls these the post-audit production result is overstated.

**Required change**

1. Add a machine-readable release policy containing the model lane, code SHA, input SHA,
   coupling policy, mass convention, RG policy, hadronic scheme, budget policy, target
   snapshot, and random seed.
2. Run a 1,000-point smoke matrix after all P0 fixes.
3. Run independent 10,000-point validation samples for each planned lane and convention.
4. Only after agreement, launch new large scans and generate all figures from their
   immutable manifest.

### P0-2 — KK-gluon coupling normalization is still a human physics decision

**Evidence**

- `scripts/run_full_catalog_scan.py` repeatedly passes `g_s_star=None`.
- This resolves to `perturbative_4d_legacy` in `quarkConstraints/couplings.py`.
- The code also implements `rs_volume_sqrt2L_physical`, which uses
  \(g_{\rm eff}=g_s^{4D}\sqrt{2L}\) in the repository's profile convention.
- The audit ledger records the legacy Lane-B floor as bit-identical and separately flags
  an approximate ~59 TeV physical-coupling requote as an open decision.

**Impact**

At fixed profiles, tree-level KK-gluon Wilson coefficients scale as
\(g_{\rm eff}^2/M_{KK}^2\), so an \(M_{KK}\) floor approximately scales with
\(g_{\rm eff}\). The possible change is therefore headline-sized, not a nuisance.

**Required change**

Do not simply adopt 59 TeV. First reproduce one benchmark from a primary RS paper using
the repository's exact \(f=F/\sqrt{2}\), \(M_{KK}\), \(\Lambda_{IR}\), and operator
normalizations. Demonstrate algebraically and numerically that the volume factor is
neither missing nor double-counted. Then make the selected policy mandatory in the scan
CLI and prohibit an implicit default in release mode. Preserve the legacy result only as
a named comparison.

### P0-3 — The Lane-A headline is internally contradictory

**Evidence**

Commit `2c9989b` states that fixing the factor-two Bauer mass bridge and Wolfenstein
\(A\) changes the anarchic median \(\epsilon_K\) wall from about 30 TeV to about
5.7 TeV. `FIX_LEDGER.md` repeats that result. In contrast:

- `README.md` says ~30 TeV;
- `CLAUDE.md` says ~30 TeV;
- `docs/STATE_OF_PROJECT.md` says ~30 TeV;
- `docs/FLOOR_SUMMARY.md` says ~30 TeV;
- `docs/MODEL_CONVENTIONS.md` says ~30 TeV.

Some of those same files also retain 18–20 TeV for the EW floor after the code/audit
reconciled it to 15.96 TeV.

**Impact**

This makes the repository's most visible physics statement unreliable. It also destroys
the intended contrast between Lane A and Lane B: under the commit-level exploratory
numbers, they may be similar rather than 30 TeV versus 7 TeV.

**Required change**

Rerun Lane A under the frozen coupling and budget policies. Generate a single signed
`headline_results.json` from raw outputs. Make all human-facing documents consume or cite
that file. Until then the honest statement is: “the corrected code reported an
exploratory ~5.7 TeV value, while existing ~30 TeV documents and artifacts are stale.”

### P0-4 — Lane C can silently drop Q4/Q5 from the default kaon observable

**Evidence**

- `HadronicArtifactBundleV1` declares Q1 VLL/VRR and Q4/Q5 LR as supported, but contains
  only \(B_K\) and a Q1 matrix element.
- `_require_hadronic_compatibility()` no longer rejects nonzero Q4/Q5.
- `compute_kaon_np_observables()` still calls `_compute_q1_only_m12_value()`.
- The returned object retains the nonzero LR Wilson coefficients, making the output look
  complete.
- `verifier.py` reads `q4_lr` and `q5_lr` but never uses them; these are two of the four
  paper-slice CI lint failures.
- PR 4 has two unresolved review threads describing this exact defect.

**Impact**

The LR contribution is usually the dominant RS \(\epsilon_K\) effect. Silent omission is
publication-blocking even though Lane C is currently quarantined.

**Required change**

Split the contracts:

1. Wilson artifacts may support Q1+Q4+Q5.
2. The default Q1 hadronic artifact must advertise Q1 only.
3. The default Q1 observable must fail closed on nonzero Q4/Q5.
4. The existing explicit combined evaluator should require the separate LR hadronic
   input object and return Q1, LR, and total pieces.
5. Add a regression test with nonzero Q4/Q5 proving that the default path raises and the
   combined path matches an independent hand contraction.

### P0-5 — B002/B004 CP phases are wrong and did affect production scans

**Evidence**

`docs/KNOWN_ISSUES.md` correctly identifies the formula problem: B002/B004 form the NP
phase against a real \(M_{12}^{SM}=\Delta m/2\), which is convention-dependent. The
physical comparison uses the complex SM box amplitude proportional to
\((V_{tq}^*V_{tb})^2\).

The same document then incorrectly says B002/B004 were not in the production scan.
They are in `QUARK_ONLY_CONSTRAINT_IDS`, their required extras are wired, and the analysis
reports show nonzero B004 veto fractions, including 100% in several 1 TeV cells and 57.1%
in one 2 TeV cell.

**Impact**

Existing strict exclusion labels contain a rephasing-dependent constraint. It may or may
not set the final floor, but the classifications are contaminated.

**Required change**

Implement the complex SM box amplitude with an explicit CKM convention and source
policy, form the invariant ratio \(M_{12}^{NP}/M_{12}^{SM}\), and add random quark-field
rephasing tests. Then rerun the affected classifications.

### P0-6 — Main is red and the advertised suite did not run in CI

**Evidence**

[GitHub Actions run 502](https://github.com/OscarBarreraGithub/5D-Neutrino-Mixing/actions/runs/29449321215)
at PR head `a68dc63` failed:

- `repo_v1_regression`: 988 Ruff errors; test suite and three benchmarks skipped;
- `paper_0710_1869_acceptance`: four Ruff errors; paper tests and benchmark skipped.

The paper failures include two unused imports and the unused Q4/Q5 variables noted above.
The general lint failures include E402, E501, I001, F401, F841, and B023 across research,
scripts, and tests. A local `python -m pytest -q` audit attempt remained in collection and
was interrupted after 561.84 seconds with zero tests executed. Collection with plugin
autoload disabled also exceeded 90 seconds. This may be NFS/environment related, but it
shows that the current local environment is not a reliable independent release runner.

**Required change**

1. Pin Python, NumPy, SciPy, Pytest, and Ruff versions in a lockfile/container.
2. Decide whether the entire historical tree is lint scope. If not, make explicit
   production/paper scopes rather than accumulating broad ignores.
3. Fix all paper-slice lint errors without masking Q4/Q5 logic.
4. Add `testpaths = ["tests"]` and diagnose collection performance.
5. Make lint, unit tests, paper acceptance, and benchmark checks separately required.
6. Protect `main`: no merge on red checks, unresolved review threads, or zero reviewers.

### P0-7 — The production model is not yet the FPR model

**Evidence**

- Production routes through `BulkMassMap`.
- It has no full \(V_{5KM}\) down-sector alignment.
- `docs/MODEL_CONVENTIONS.md` explicitly admits this, while also using “our FPR-MFV”
  language elsewhere.
- Lane C still lacks exact Table-I affine coefficients and an in-repository source
  extraction from arXiv:0710.1869.

**Impact**

The paper cannot claim a modern validation of the Lisa/FPR model based on Lane B. It can
claim results for a simplified, fitted, MFV-inspired surrogate if that model is defined
as its own theory.

**Required change**

Choose one of two honest paper scopes:

- **FPR paper:** complete and benchmark Lane C, then use Lane B only as a computational
  surrogate comparison; or
- **new surrogate paper:** write a self-contained Lagrangian/parameter definition for
  Lane B and stop identifying its floor with FPR.

For a Lisa-facing project, the first option is much stronger.

---

## 5. P1 physics and reproducibility work

These items do not all block a first corrected smoke scan, but they block a durable
precision/publication claim.

### P1-1 — Upgrade \(\Delta F=2\) theory coherently

- Replace LO-only evolution with a scheme-consistent NLO implementation or interface to
  a trusted external EFT package.
- Evolve Wilson coefficients and lattice matrix elements in the same basis, scheme, and
  scale.
- Resolve the remaining kaon use of 2 GeV quark masses with 3 GeV bag parameters.
- Investigate whether VIA additive terms double-count contributions.
- Synchronize the stale core \(\Delta m_D\) input with the catalog.
- Make \(\Delta m_K\) inclusion an explicit policy rather than an accidental omission.
- Remove ambiguous `bound` fields whose semantics differ from `ratio_to_bound`.

The audit estimates the remaining LO bias in the \(\epsilon_K\) floor at roughly
10–15%; this should be a theory uncertainty until upgraded.

### P1-2 — Replace proxy global constraints with controlled likelihoods

- Use covariance-aware quark mass/CKM residuals rather than factors of 3 and 10.
- Treat experimental and SM theory uncertainties with named likelihood/nuisance models.
- Replace the C002 one-dimensional charm CP proxy with the appropriate correlated
  mixing likelihood if it is to be a HARD constraint.
- Define direction and confidence-level semantics for T010/T011.
- Ensure the `z=1.92` finite-statistics convention is never called standard one-sided
  95% coverage.

### P1-3 — Rebuild the electroweak/custodial result

- Make the warp volume \(L\) come from the active geometry instead of hardwired 35.
- Re-derive the oblique coefficient and quantify the flagged ~20% uncertainty.
- Pin a modern electroweak-fit covariance source; do not label a 2018 anchor as 2026.
- Implement the full custodial representation/top-partner spectrum and loop effects.
- Derive rather than assume custodial off-diagonal FCNC behavior.
- Keep the current tree-level \(P_{LR}\) branch labeled as a proxy until those pieces
  exist.

### P1-4 — Recast collider constraints for this model

Most CR entries reduce to mass cuts because the sigma-times-branching-ratio path is not
populated. The paper needs model-specific production, widths, branching fractions,
acceptances, and experimental likelihoods, or it must explicitly omit collider limits
from the rigorous floor. The current VLQ proxy `m_VLQ=M_KK` and contact-interaction
normalizations are not sufficient.

### P1-5 — Complete the lepton scope or narrow the paper

No charged-lepton constraint was in the paired production scans. If the paper claims a
joint quark/lepton RS model, it must:

- replace the nonperturbative lepton benchmark;
- pin current oscillation inputs;
- implement live LMFV point construction across the scan;
- resolve \(\mu N\to eN\) target-material inputs and phase policies;
- include neutrino-mass and \(0\nu\beta\beta\) constraints where applicable;
- document or derive the seesaw convention.

Otherwise the title and abstract should say “quark-sector RS flavor.”

### P1-6 — Make inputs and artifacts immutable

- Add `analysis`/`paper` dependency extras for pandas, matplotlib, pyarrow, and any other
  imported analysis package. The current declared runtime dependencies are only NumPy
  and SciPy even though tests import PyArrow and analysis scripts use pandas/matplotlib.
- Put every release input under a checksummed manifest with source URL, access date,
  scheme, scale, and covariance.
- Write output atomically so stale `.tmp` shards cannot be double counted.
- Never default missing `fit_success` to true in a release scan.
- Store a reduced but sufficient release dataset and exact regeneration command. Raw
  multi-million-row files may live in an archive, but the paper must not depend on an
  unversioned lab path.

### P1-7 — Reclassify the test suite

Tests should be labeled as:

1. algebraic identity/unit tests;
2. independent analytic or numerical oracles;
3. literature benchmark reproductions;
4. golden regression snapshots;
5. end-to-end reproducibility/statistical tests.

The audit found about 105 long-float pins and several tests that reproduce the same
formula as the implementation. Those are useful regression guards, but they must not be
presented as independent physics validation.

### P1-8 — Replace hand-maintained status prose with one generated truth source

The current documentation contains independent contradictions beyond the Lane-A floor:

- `README.md` retains an 18–20 TeV EW floor while the audited code result is 15.96 TeV.
- `docs/STATE_OF_PROJECT.md` says the Explorer is uncommitted, then says it is committed,
  and elsewhere says it is committed locally but not pushed. The audited branch is
  synchronized with `origin/main`.
- The same state file says the \(\epsilon_K\) budget differs between core and catalog,
  even though the July M-1/M-2/M-4 cycle deliberately unified the policy.
- `docs/MODEL_CONVENTIONS.md` announces “TWO lanes” immediately before a three-lane
  table.
- The project metadata still describes a lepton-only package while the active paper is
  quark flavor.

Use a small generated status block backed by `headline_results.json`, scan manifests,
and Git metadata. Keep narrative documents, but make numeric/status claims reference the
generated block rather than duplicating them.

---

## 6. Residual lower-severity inventory

The July campaign deliberately marked roughly 50 findings `DOC` rather than fixing them.
They should be converted into real GitHub issues with owner, severity, affected claims,
and a validation condition. The important residuals are:

### 6.1 Lepton/numerical core

- **LC-01:** flagship lepton benchmark violates its own perturbativity criterion.
- **LC-02:** tiny-argument `_F_exact` branch has a fake-positive sign.
- **LC-03:** fallback root scan searches only upward from the last seed.
- **LC-05:** NuFIT 6.1 values require public-source verification.
- **LC-06:** the naturalness cut excludes exactly massless-neutrino choices by design.
- **LC-08:** the seesaw prefactor is internally consistent but not derived in the docs.

### 6.2 Quark model and QCD

- **QM-01:** CKM target central values are not the claimed current set.
- **QM-03:** `np.isclose(c,0.5)` creates a finite plateau at the flat-profile point.
- **QM-04:** a legacy no-hadronic fallback applies kaon LR weights to all mesons.
- **QM-05:** diagnostic scale layout is hardwired to \(m_t(m_t)\).
- **QCD-02:** mass running ignores the requested reference-flavor scheme in one path.
- **QCD-04:** light-quark uncertainty symmetrization does not match its documentation.

### 6.3 \(\Delta F=2\), rare decays, kaon, and charm

- **DF2-01:** possible VIA additive double count.
- **DF2-02/KCH-03:** stale core charm-mixing constant despite a catalog override.
- **DF2-03:** unpublished precision in the displayed SM \(\epsilon_K\) value.
- **DF2-04:** \(\Delta m_K\) not in the default observable bundle.
- **DF2-05:** ambiguous kaon bound API.
- **RD-02:** opposite lepton-axial signs between kaon and B/charm paths require a
  convention audit.
- **RD-03:** \(B\to\tau\nu\) has silent default inputs outside the YAML path.
- **KCH-01:** K008 coefficients conflict with its YAML SM limit.
- **KCH-02:** Grossman–Nir policy is not applied to the KOTO/NA62 combination.
- **KCH-04:** direct charm CP is non-vetoing.
- **NMF-02:** K009 compresses vector and axial direct-CP weights.
- **NMF-06:** \(K_L\) LFV modes omit CP projection and an isospin Clebsch.
- **NMF-07:** exact \(\xi_{KK}\) rescaling neglects the shift in the RG start scale.
- **NMF-14:** a field named for the charged-kaon lifetime carries a \(K_L\) lifetime.
- **NMF-15:** B012 uses a bare one-sigma budget with SM set equal to measurement.

### 6.4 RS electroweak, collider, and catalog

- **RSEW-01:** a custodial SU(2)R name aliases the \(P_{LR}\) model.
- **RSEW-02:** oblique proxy hardwires \(L=35\).
- **RSEW-03:** the geometric \(S\)-coefficient may be ~20% low.
- **ECL-01:** bookkeeping `M_KK` is exported under a KK-gluon mass name.
- **ECL-03:** T010/T011 direction/CL semantics are incomplete.
- **ECL-04:** collider sigma-times-BR path is never populated.
- **ECL-05:** cosmological neutrino mass, \(0\nu\beta\beta\), \(g-2\), and robust
  hadronic EDM constraints are missing or informational.
- **ECL-07:** aluminum \(\mu\to e\) conversion borrows a gold-nucleus budget.
- **NMF-04:** vectorlike-quark mass is proxied by \(M_{KK}\).
- **NMF-05:** CR009 combines destructive endpoints and omits a contact normalization.

### 6.5 Pipeline, figures, and derivations

- **SP-02:** stale temporary files can be double counted after a crash.
- **SP-03:** the code's 5.5 TeV collider cut conflicts with 4 TeV prose.
- **SP-05:** EW003 is in a quark allowlist but needs a forbidden/missing extra.
- **SP-06:** missing fit diagnostics default to success.
- **NMF-01:** LFV unknown-phase envelopes use inconsistent constructive/destructive
  policies.
- **NMF-08:** a CFW figure fixes \(g_s\) while the run varies it.
- **NMF-09:** figure scripts disagree on fit filtering and can map nonfinite values to
  passing.
- **NMF-10:** phase calibration tolerances are based on achieved residuals, not data
  uncertainties.
- **NMF-11:** Lane-A/Lane-B mass prefactors remain insufficiently derived/documented.
- **NMF-12:** derivation READMEs mark never-written files as done.
- **NMF-13:** long-float snapshots dominate parts of the test suite.
- **NMF-17:** MEG-II comparison prose says 100x where the script implies about 80x.

These items are not all equal. DF2-01, RD-02, RSEW-02/03, ECL-03/04/07, SP-06,
NMF-07/09/11/12 should be treated as paper-relevant; spelling/field-name items may wait.

---

## 7. Short explanation of \(\epsilon_K\)

Neutral kaons mix through an off-diagonal amplitude \(M_{12}^K\). The mass difference
mostly probes its real part, while indirect CP violation is controlled by its imaginary
part. In the approximation used by this repository,

\[
  \epsilon_K^{\rm NP}\simeq
  \frac{\kappa_\epsilon}{\sqrt{2}\,\Delta m_K}
  \operatorname{Im} M_{12}^{K,\rm NP}.
\]

An RS KK gluon couples non-universally to zero-mode quarks because their 5D profiles are
different. After rotating to the mass basis, tree-level exchange generates left-left,
right-right, and mixed-chirality \(\Delta S=2\) operators. The dangerous term is usually
the LR coefficient

\[
  C_4^{sd}\propto -\frac{g_{s*}^2}{M_{KK}^2}
  (G_L^d)_{12}(G_R^d)_{12},
\]

because its hadronic matrix element is chirally enhanced and QCD running enhances the LR
sector. Generic complex Yukawas make \(\operatorname{Im}C_4\) large. One can survive by:

1. increasing \(M_{KK}\);
2. reducing the off-diagonal coupling magnitudes with degeneracy/alignment; or
3. tuning the relevant phase close to zero or \(\pi\).

The third route is fragile; the second can be symmetry-protected. Custodial symmetry
solves electroweak problems such as \(T\) and \(Zb_L\bar b_L\), but does not automatically
align color/flavor couplings. That is why \(\epsilon_K\) remains central even in a
custodial model.

The repository now has a sign-aware \(\epsilon_K\) budget, but no single current floor
should be quoted: Lane A has contradictory stale numbers, Lane B uses a legacy coupling,
Lane C is quarantined, and the large scans predate the audit.

---

## 8. What happened with the Yukawa perturbation question

### 8.1 The question that was asked

The requested experiment had four parts:

1. perturb each Yukawa entry and find how far it can move before the fit fails, including
   how badly and which observable fails;
2. add random sub-percent Gaussian noise to the complete Yukawa matrices;
3. treat the local response as a point-dependent linear map on Yukawa space;
4. study the geometry of those maps and decide whether the allowed data lie on a
   lower-dimensional submanifold.

### 8.2 What Claude/Opus actually completed

The July F7 study selected four point classes at \(M_{KK}=3\) TeV and performed:

- real multiplicative noise, by default on \(Y_d\) only;
- a central finite-difference gradient of \(\operatorname{Im}C_4\) with respect to the
  18 real components of \(Y_d\);
- PCA/SVD of normalized gradient directions over small point ensembles.

It reported:

| Class | Noise result | Local reported tuning radius | Interpretation in the study |
|---|---|---:|---|
| phase-tuned anarchy | below 50% survival by roughly \(\sigma\sim5\times10^{-3}\); median failure reaches ~36x at large noise | \(7.9\times10^{-4}\) | very fragile phase cancellation |
| typical magnitude-suppressed anarchy | robust at percent scale; ~55% at 10% noise | 0.15 | smaller local gradient |
| rank-one/U(2) toy | 100% through 3%, 85% at 10%, 63% at 30% | 0.024 | magnitude suppression plus partial directional protection |
| real-\(Y_d\) toy called “Nelson–Barr” | 100% under all tested real multiplicative noise | infinite in real directions | real noise preserves real \(Y_d\), so the result is exact by construction |

The ensemble gradients for the real-\(Y_d\) toy occupy the nine imaginary directions,
with a sharp rank drop after component nine. That is a clean demonstration that this
particular constraint and perturbation preserve a nine-real-dimensional flat. The U(2)
toy showed a more diffuse, mostly imaginary-sensitive response rather than an exact rank
cliff.

### 8.3 The clear answer: useful proof of concept, not the requested fit boundary

The result does support the original intuition:

> A point surviving \(\epsilon_K\) by phase cancellation can be destroyed by
> sub-percent Yukawa noise, while a point protected by an exact reality condition is
> insensitive to reality-preserving perturbations.

It does **not** yet answer “what is the maximum perturbation each Yukawa entry can receive
before the model no longer fits.” The reasons are concrete:

1. `multiplicative_noise()` computes `passes_pdg` but does not use it in the pass count;
   it counts only `ratio_eps_K <= 1`.
2. The profiles \(f_Q,f_u,f_d\) are held fixed after changing the Yukawas. There is no
   refit or profiling of bulk parameters, so the study measures a frozen-background
   response, not distance to the fitted solution set.
3. The CLI never exposes the internally supported `perturb="both"`; the saved F7 run
   perturbs \(Y_d\) by default, not the full \(Y_u,Y_d\) system.
4. The noise is real and multiplicative. It cannot explore generic CP phases and barely
   moves entries near zero. The “Nelson–Barr” immunity is therefore tautological for the
   chosen perturbation class and is not a demonstration of a complete Nelson–Barr model.
5. The phase-tuned point is the largest-\(|C_4|\) survivor found in a finite search,
   introducing best-of-budget selection bias.
6. Seeds use Python's salted `hash()` in this script, so the advertised seed is not fully
   reproducible across processes unless `PYTHONHASHSEED` is fixed.
7. The finite-difference step and Euclidean norm are coordinate- and normalization-
   dependent. Quark-field rephasings and weak-basis rotations have not been quotiented.
8. The “PCA effective dimension” is the span of gradient directions across an ensemble,
   not the dimension/codimension of the pointwise allowed manifold. A single scalar
   constraint generically has a 17-dimensional local kernel in an 18-dimensional
   \(Y_d\) space.
9. Only \(\operatorname{Im}C_4\) is differentiated. Masses, CKM, \(\Delta m_K\), B/D
   mixing, EW, perturbativity, and the full \(\epsilon_K\) contraction are absent.
10. No Hessian, active-set/tangent-cone treatment, uncertainty estimate, or held-out
    validation is present.
11. The calculation uses \(M_{KK}=3\) TeV and the legacy coupling/profile conventions,
    so its numerical radii are not release results.
12. The research script has no targeted tests and is part of the current lint failure.

Calling the tuned point “isolated” is therefore too strong. It is fragile along the
measured normal direction, but the pointwise kernel of one observable is large. The
interesting question is how that kernel intersects the mass/CKM-fit tangent space and
the normals of all other active constraints.

### 8.4 What the older Yukawa sheet does suggest

An older utility, `scripts/yukawa_per_element_anatomy.py`, was written to ask whether the
accepted envelope points use one collective off-diagonal suppression or tune every entry
independently. It was not incorporated into the committed results. A read-only evaluation
of its 83,961-point input during this audit gives:

- \(Y_u\) off-diagonal within-point log spread: median 0.534 dex, 95th percentile
  0.630 dex;
- \(Y_d\) off-diagonal within-point log spread: median 0.148 dex, 95th percentile
  0.339 dex;
- pooled off-diagonal medians: about 0.133 for \(Y_u\), 0.182 for \(Y_d\).

This is a genuinely interesting clue: in that old optimizer sample, the down-sector
off-diagonals move much more collectively than the up-sector ones. It is compatible with
an approximate low-dimensional down-sector suppression spurion rather than six unrelated
tunings. However, the dataset is a pre-audit existence/envelope fit, the statement is
basis-dependent, and the result has not been reproduced on corrected, independently
accepted points. It is a hypothesis generator, not yet a result.

---

## 9. The perturbation study that should be done

### 9.1 Define the physical space and residual map

Use the 36-real-dimensional pair

\[
 x=(\operatorname{Re}Y_u,\operatorname{Im}Y_u,
    \operatorname{Re}Y_d,\operatorname{Im}Y_d)
\]

plus nuisance/profile parameters \(n\). Fix a canonical weak basis or explicitly quotient
the quark-field unitary/rephasing redundancy. Define a whitened residual vector

\[
 r(x,n)=\bigl(r_{m_q},r_{CKM},r_{J/\delta},r_{\epsilon_K},r_{\Delta m_K},
 r_{B_d},r_{B_s},r_D,r_{EW},r_{\rm pert},\ldots\bigr)
\]

with each component normalized by a named covariance or likelihood scale.

Report two distinct notions of sensitivity:

- **frozen response:** perturb Yukawas and hold profiles/nuisances fixed;
- **profiled response:** perturb Yukawas, then re-minimize over allowed nuisances and
  bulk parameters while preserving the model ansatz.

The first diagnoses the local mechanism. The second answers the PI's “does it still fit?”
question.

### 9.2 Do the naive experiment completely

For every real and imaginary component of \(Y_u\) and \(Y_d\):

1. move in the positive and negative additive directions;
2. move in positive and negative fractional directions where the entry is nonzero;
3. use a bracket plus bisection to find the first fit/constraint failure;
4. record the signed distance, the failing observable, its pull at failure, and its pull
   after fixed 0.1%, 0.3%, 1%, 3%, 10%, and 30% moves;
5. repeat after profiling the nuisance/profile parameters;
6. bootstrap across many independently accepted points in every model class.

For random noise, perturb both matrices with real and complex Gaussian ensembles. Quote
the distribution of critical noise amplitudes and the 50%/90%/95% survival thresholds,
not one best point.

### 9.3 Construct the point-dependent linear map correctly

Compute the Jacobian

\[
 J=\frac{\partial r}{\partial(x,n)}.
\]

Use stable central differences with a step-convergence study, complex-step derivatives
where analytic operations allow them, or autodifferentiation after isolating
nondifferentiable fit logic. Whiten the rows by the data covariance and scale the columns
with a physically motivated metric.

For profiled sensitivity, eliminate nuisance directions with the local least-squares
projector/Schur complement. SVD of the resulting \(J_{\rm eff}\) provides:

- stiff combinations that immediately spoil the fit;
- soft/tangent combinations that keep all active observables stable;
- an effective codimension and condition spectrum;
- the observable composition of each normal direction.

Because pass/fail constraints are inequalities, the boundary is generally a stratified
set with corners, not one smooth manifold. At a boundary point, use the gradients of all
active constraints to build the tangent cone. Call a region a submanifold only where the
Jacobian has locally constant rank and the active set is unchanged.

### 9.4 Go beyond linear response

- Compute directional Hessians/second fundamental forms along the soft directions.
- Continue geodesics/tangent predictor-corrector paths and measure when the active set
  changes.
- Track local stiff/soft subspaces across points using principal angles on a Grassmann
  manifold, rather than PCA of uncentered covectors.
- Estimate the reach/thickness of each surviving sheet and the viable tube volume under a
  specified prior. This gives an interpretable naturalness measure.
- Validate the local model by predicting finite perturbation failures on held-out points.

The general information-geometry language is established outside HEP: multiparameter
models often form “hyperribbon” prediction manifolds with stiff and sloppy combinations.
See [Quinn et al., arXiv:2111.07176](https://arxiv.org/abs/2111.07176). The potential
novelty here is the gauge-invariant application to warped Yukawa flavor and the recovery
of protecting spurions, not the use of an SVD by itself.

---

## 10. Research directions worth discussing with Lisa

These are candidate directions, not claims of literature novelty. A systematic INSPIRE
review is required before choosing the paper pitch. Simple right-handed-down flavor
protection is already known (for example
[Santiago, arXiv:0806.1230](https://arxiv.org/abs/0806.1230)), and alternative warped
solutions to the mixed-chirality \(\epsilon_K\) problem already exist (for example an
extended bulk color group in
[Bauer, Malm, Neubert, arXiv:1110.0471](https://arxiv.org/abs/1110.0471)). Merely adding
“a symmetry” is not enough.

### Direction A — Response geometry as a symmetry detector

Build a coordinate-invariant atlas of the fitted Yukawa set and ask whether its soft
directions align with generators of approximate CP, U(2), or FPR/MFV spurion symmetries.
Instead of assuming the symmetry and showing it works, infer the approximate generators
from the Jacobian kernels, then test whether they reproduce the known analytic flavor
structure.

**Possible publishable result:** accidental \(\epsilon_K\) survivors, U(2)-like points,
and the exact FPR alignment have quantitatively different normal spectra, curvature, and
viable tube volumes even when they give the same mass floor.

### Direction B — Stability of the original FPR limit to MFV-breaking spurions

Complete exact FPR first. Add the most general small symmetry-breaking bulk-mass/Yukawa
spurions, organized by representations, and determine which combinations first regenerate
the LR kaon operator. Measure the scaling

\[
 \epsilon_K^{NP}\sim \delta^p\,g_{s*}^2/M_{KK}^2
\]

for each breaking direction and identify linear versus quadratic protection. This is
directly tied to the Lisa paper and answers a more important question than whether the
exactly aligned point still fits modern data: **how technically stable is the mechanism?**

### Direction C — A naturalness phase diagram, not one floor

Map \((M_{KK},g_{s*},r,\delta_{\rm break})\) to:

- viable prior volume;
- tangent-space dimension;
- smallest profiled perturbation radius;
- active constraint;
- required phase alignment.

This separates three mechanisms that a one-dimensional floor conflates: decoupling,
magnitude suppression, and phase cancellation. A “how much flavor protection is enough?”
phase diagram would be much more informative than another bound table.

### Direction D — Correlated CP signals along the \(\epsilon_K\)-safe tangent space

Once EDM and B/D CP calculations are made reliable, move along directions that remain
tangent to the \(\epsilon_K\) and mass/CKM constraints. Determine which observables grow
first: neutron/electron EDMs, \(B_{d,s}\) phases, charm CP, rare kaons, or top/Higgs FCNCs.
The old RS literature already identifies an EDM/CP problem, so the new content would be
the response-geometric correlation and symmetry classification, not the existence of an
EDM bound.

### Direction E — The collective down-sector suppression clue

Reproduce the old per-entry anatomy on a corrected ensemble and decompose the accepted
\(Y_d\) matrices into a small number of collective spurions. Test whether the 0.15-dex
within-point spread survives basis fixing, profiling, and out-of-sample prediction. If one
collective mode explains the accepted down-sector off-diagonals while \(Y_u\) requires
several, that could provide a data-driven bridge from the scan to an analytic flavor
ansatz.

### Recommended paper concept

The most coherent combination is A+B+C:

> **The geometry and robustness of flavor protection in warped space:** exact modern FPR
> reproduction; controlled symmetry-breaking deformations; a gauge-invariant
> naturalness atlas distinguishing accidental phase tuning from symmetry-protected
> Yukawa manifolds; and correlated observables along the surviving directions.

That is recognizably connected to Lisa's original construction but asks a new and useful
question. It also forces the repository to become scientifically cleaner because every
result depends on an exact model definition, physical convention, likelihood, and
reproducible local response.

---

## 11. Sequenced execution plan

### Gate 0 — Freeze claims immediately (1–2 days)

- Mark all current floor tables as historical/pending rerun.
- Add GitHub issues for every P0 item and the paper-relevant P1/minor items.
- Decide that no number without lane + coupling + mass + budget + scan SHA is quotable.
- Protect `main` and require green checks/review.

**Exit condition:** no document presents a superseded scan as current.

### Gate 1 — Convention and model note (2–5 days, human physics review)

Write a short derivation note fixing:

- \(f\) versus \(F\) normalization;
- \(v=174\) versus 246 GeV and every Lane-A/B mass prefactor;
- physical first gauge KK mass versus \(\Lambda_{IR}\);
- 4D \(g_s\), \(g_{s*}\), and volume factors;
- \(\Delta F=2\) Hamiltonian/operator signs;
- exact Lane A/B/C definitions;
- the selected \(\epsilon_K\) and statistical policies.

Reproduce at least one published KK-gluon Wilson benchmark before selecting the
production coupling.

**Exit condition:** one machine-readable policy and one derivation are signed off by a
physicist.

### Gate 2 — Repair P0 code and CI (3–7 days)

- Fix Lane-C Q1/LR contracts and default fail-closed behavior.
- Fix B002/B004 complex SM box phases and add rephasing tests.
- Remove salted-hash seeds from all research/release scripts.
- Make paper lint and tests green, then clean/scoped general CI.
- Pin a reproducible environment and repair collection performance.

**Exit condition:** green CI at the exact commit, zero unresolved review threads, all
paper benchmarks executed rather than skipped.

### Gate 3 — Exact FPR benchmark (1–2 weeks)

- Vendor a checksummed source extraction from arXiv:0710.1869.
- Implement/verify Table-I coefficients and Eq. (3) conventions.
- Reproduce a published parameter point and the down-aligned limit.
- Cross-check the calculation with an independent minimal implementation, not the
  production classes.

**Exit condition:** Lane C passes a literature benchmark with an explicit tolerance and
no circular oracle.

### Gate 4 — Small validation matrix (about 1 week)

Run common random seeds for:

- Lane A corrected anarchy;
- Lane B surrogate;
- Lane C exact FPR;
- minimal and custodial EW branches;
- legacy and selected physical coupling policies.

Use 1k smoke and 10k validation samples. Compare analytic scaling, fit rates, constraint
composition, and finite-statistics intervals.

**Exit condition:** all changes from old headlines are explained quantitatively before a
large job is submitted.

### Gate 5 — Production reruns and release artifact (1–3 weeks plus compute)

- Run new quark production ensembles with immutable manifests.
- Do not overwrite old artifacts; retain them under an explicitly historical namespace.
- Generate every table/figure from a single release command.
- Build the website from the release artifact only.
- If the paper claims leptons, run the separate live lepton lane after its benchmarks are
  fixed.

**Exit condition:** a fresh clone/container can reproduce the headline table and plots
from committed/released inputs.

### Gate 6 — Full perturbation geometry (2–4 weeks)

- Implement per-entry additive/fractional line searches and random real/complex noise.
- Add frozen and profiled modes.
- Build the full whitened residual Jacobian and active tangent cone.
- Bootstrap radii/subspaces across model classes.
- Add Hessian/geodesic continuation for promising sheets.
- Test held-out finite perturbations.

**Exit condition:** the reported dimension/radius is basis-aware, includes the full fit,
and predicts held-out failures.

### Gate 7 — Precision upgrades and research observables (parallel after the core is stable)

- NLO \(\Delta F=2\) and scheme-consistent lattice inputs;
- full custodial/top-partner physics;
- model-specific collider recasts;
- EDM and CP-correlated observables;
- current global input/likelihood updates.

These determine whether the first paper is a robust modern phenomenology paper or a
method/naturalness paper with clearly labeled leading-order constraints.

---

## 12. Definition of “solid physics repo for Lisa”

The repository is ready for Lisa when all of the following are true:

- one precise model, not a mixture of Lane B and Lane C language, is the paper subject;
- every headline number names its lane, \(M_{KK}\) definition, coupling policy, data
  snapshot, theory order, likelihood/budget, and scan SHA;
- exact FPR benchmark reproduction exists if FPR is in the title;
- the dominant \(\epsilon_K\) LR operator is never silently omitted;
- CP observables are rephasing invariant;
- current production scans postdate all physics fixes;
- GitHub CI is green and required, with tests and benchmarks actually executed;
- at least one independent oracle covers each headline calculation;
- release dependencies and data are immutable and reproducible outside the lab account;
- proxy/partial constraints never appear in a rigorous exclusion without explicit
  qualification;
- the paper makes a new claim about robustness, geometry, symmetry breaking, or
  correlated predictions—not only that a known protected RS model can still fit.

---

## 13. Key evidence locations

- Audit disposition: `docs/audits/full_repo_audit_2026-07/FIX_LEDGER.md`
- Residual findings: `docs/audits/full_repo_audit_2026-07/MINOR_FINDINGS_LEDGER.md`
- Current but inconsistent state: `docs/STATE_OF_PROJECT.md`
- Lane definitions: `docs/MODEL_CONVENTIONS.md`
- Current floor claims: `docs/FLOOR_SUMMARY.md`, `README.md`, `CLAUDE.md`
- False B002/B004 scan statement: `docs/KNOWN_ISSUES.md`
- Production allowlist: `scripts/run_full_catalog_scan.py`
- Coupling policies: `quarkConstraints/couplings.py`
- Production legacy coupling calls: `scripts/run_full_catalog_scan.py`
- Lane-C observable defect: `quarkConstraints/paper_0710_1869/eft_deltaf2/observables.py`
- Lane-C artifact defect: `quarkConstraints/paper_0710_1869/artifacts.py`
- Lane-C verifier symptom: `quarkConstraints/paper_0710_1869/verifier.py`
- Superseded scans: `scan_outputs/wq_quarkonly_1M_20128400/analysis/analysis_report.md`,
  `scan_outputs/wq_quarkonly_1M_custodial_20675555/analysis/analysis_report.md`
- Perturbation writeup: `.orchestration/runs/RS-FLAVOR-ALIGNMENT-2026-07/PERTURBATION_STUDY.md`
- Perturbation implementation: `scripts/yukawa_perturbation_study.py`
- Older per-entry anatomy: `scripts/yukawa_per_element_anatomy.py`
