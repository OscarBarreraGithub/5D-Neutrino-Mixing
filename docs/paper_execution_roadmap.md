# Roadmap for the Quark-Scan Paper Branch, v2

## A. State of the repo

The repo has a substantial quark-sector framework, but it is still a working-tree state rather than a paper-ready state. At review time the branch was `main`, with 14 modified tracked files and 36 visible untracked entries, plus ignored heavy artifacts. The visible `.claude` issue is precise: `.claude/scheduled_tasks.lock` is visible as untracked, while `.claude/settings.local.json` is ignored through the user's global gitignore, not through the repo `.gitignore`. Treat agent-local files as local state and do not let them enter the paper branch accidentally.

The new code appears to cover three linked pieces: PDG 2024 quark target construction through `qcd/mass_running.py::run_msbar_mass` and `quarkConstraints/pdg_quark_masses.py::pdg_quark_masses_at_scale`; RS-anarchy scanning through `scripts/run_rs_anarchy.py`, including `_draw_anarchic_matrix`, `_check_pdg_match`, and `_evaluate_one_draw`; and Delta F=2 constraint evaluation through `quarkConstraints/deltaf2.py`, including `compute_delta_f2_wilsons`, `_evolve_wilsons`, and `evaluate_delta_f2_constraints`.

The methodology note is a real draft: `docs/quark_scan_methodology_note.tex` is 904 lines. It is also currently ignored by repo policy because `.gitignore` ignores `docs/*` and only unignores the compact assumptions note. That means ordinary status output can hide a key paper source unless the ignore rules are amended or the file is force-added deliberately.

Several headline claims are backed by code or local run outputs, but not all are paper-grade. `pytest -q tests/test_quark_fit.py` was verified at `12 passed, 2 xfailed` in about 5.90 s, not 85 s. The scan-result claims are supported by local `scan_outputs/` JSON, but no blessed manifest yet ties a code SHA, command line, environment, output file, and figure hash to each quoted paper number.

## B. Holes to patch, prioritized

1. **Unsafe dirty-tree branch procedure**

**What's wrong:** The current paper work sits on `main` with modified files, visible untracked files, and ignored artifacts. The previous plan started logical commit splitting before proving that local `HEAD` matched `origin/main`, before saving a durable snapshot, and before checking branch-name collisions.

**Why it matters:** A mistaken pathspec or interrupted split could obscure the exact state that generated the draft claims.

**Effort:** M/L. The git commands are simple; the coordination and path inventory are not.

**Severity:** BLOCKER.

2. **No canonical run snapshot or final artifact blessing**

**What's wrong:** The methodology note cites baseline, RUNA, Run B, Run C, gate sensitivity, and follow-up crossing statistics from ignored local outputs. There is no tracked manifest identifying which artifacts are canonical.

**Why it matters:** Reproducibility requires more than hashes. Someone must decide which local outputs are final, which are scratch, and which figure was generated from which summary.

**Effort:** L unless a complete artifact inventory already exists.

**Severity:** BLOCKER.

3. **Final paper claims are not gated on bag/RG audits**

**What's wrong:** The previous plan allowed a manifest and bound language to solidify before the hadronic bag-parameter and Wilson-RG audits were complete.

**Why it matters:** If constants, schemes, thresholds, or Wilson normalizations change, every quoted `M_KK^min` p50/p95 and affected figure must be rerun or revalidated.

**Effort:** M if both audits certify no numerical change; L if rescans or deterministic re-evaluation are required.

**Severity:** BLOCKER.

4. **Two xfailed quark-fit tests under the new PDG targets**

**What's wrong:** `tests/test_quark_fit.py` still has two xfails. Removing decorators is not a validation criterion by itself, and non-strict xfail handling can hide XPASS behavior.

**Why it matters:** The failures touch benchmark logic used to interpret whether the new PDG targets preserve intended physics behavior.

**Effort:** M if only fixtures and tolerances need repair; L if benchmark construction changes.

**Severity:** BLOCKER.

5. **Hadronic bag-parameter provenance is not documented**

**What's wrong:** `quarkConstraints/deltaf2.py` encodes hadronic inputs without a complete source, scale, scheme, and uncertainty table in the methodology note.

**Why it matters:** Bag parameters are direct systematics in Delta F=2 constraints.

**Effort:** M for provenance; L if values change and outputs must be regenerated.

**Severity:** BLOCKER.

6. **Wilson-coefficient RG running has not been audited**

**What's wrong:** `_evolve_wilsons` and related evaluators need a basis, anomalous-dimension, threshold, direction, and matching-convention audit.

**Why it matters:** Epsilon_K-dominated bounds can move under convention errors.

**Effort:** L.

**Severity:** BLOCKER.

7. **CFW reproduction needs convention reconciliation, not access troubleshooting**

**What's wrong:** arXiv/PDF access is not the main risk. The real risk is reconciling Yukawa normalization, priors, observable definitions, historical inputs, and coupling conventions.

**Why it matters:** The paper should not imply a quantitative reproduction unless the old-paper conventions are actually matched.

**Effort:** L.

**Severity:** IMPORTANT.

8. **Zero-pass scan claims need finite-statistics treatment**

**What's wrong:** moreUV/moreIR c-shift scans and Run C factor-1.5 gate results report zero PDG passes, but zero observed successes in a finite ensemble is not an absolute impossibility statement.

**Why it matters:** Negative scan claims need sample size, gate, seed/tile provenance, and a one-sided upper limit.

**Effort:** M for manifest and wording; L if additional scan statistics are needed.

**Severity:** IMPORTANT.

9. **KK tower, EWPO floor, `g_s*`, and lepton scope are not closed**

**What's wrong:** The repo computes a first-mode flavor proxy, cites an EWPO floor externally, and focuses on quarks while living in a broader 5D flavor repository.

**Why it matters:** The document must distinguish computed results from cited assumptions and scope boundaries.

**Effort:** M.

**Severity:** IMPORTANT.

10. **Figure and document hygiene is not paper-ready**

**What's wrong:** Many local figures exist, but only a small named set should be tracked. README and contributor-facing notes lag the new framework.

**Why it matters:** Collaborators need a small, regenerable paper surface, not a scratch directory.

**Effort:** M.

**Severity:** IMPORTANT.

## C. Execution plan in 3 phases

### Phase 1: Repo hygiene & branching

This phase must finish before physics patching. The first goal is to preserve the exact pre-paperwork state, then split it into topic commits on a new branch. Do not commit the draft methodology note directly to `main`.

```bash
git switch main
git fetch origin && git status
# Confirm HEAD == origin/main before proceeding. If this fails, stop and record divergence.
test "$(git rev-parse HEAD)" = "$(git rev-parse origin/main)"

mkdir -p snapshots
git status --porcelain=v1 --untracked-files=all > snapshots/pre-paperwork-status.txt
git diff > snapshots/pre-paperwork-tracked.diff
git diff --staged > snapshots/pre-paperwork-staged.diff
git ls-files --others --exclude-standard -z | \
  tar --null -T - -czf snapshots/pre-paperwork-untracked.tgz
git bundle create snapshots/pre-paperwork-$(date +%Y%m%dT%H%M%S).bundle main

git branch --list paper/quark-scan-2026q2
# If the branch exists, stop and choose an explicit non-colliding name.
git switch -c paper/quark-scan-2026q2
```

After that snapshot and branch creation, do the topic-grouped commits on `paper/quark-scan-2026q2`. Recommended groups are: repo tracking rules; PDG quark mass targets; quark scan benchmarks and diagnostics; RS-anarchy scan driver and dispatch scripts; follow-up analysis utilities; plotting scripts; executed analysis notebooks; and draft paper docs. The `.gitignore` change should deliberately handle `docs/quark_scan_methodology_note.{tex,pdf}` and local agent state. Prefer ignoring `.claude/` as local state unless there is a repo reason to preserve a specific file.

Use per-hole branches off the paper branch: `fix/pdg-benchmarks`, `audit/bag-inputs`, `audit/wilson-rg`, `audit/cfw-comparison`, `scan/zero-pass-statistics`, and `docs/paper-readiness`. Merge back by reviewed PR or reviewed fast-forward. Keep `main` untouched until a release-candidate tag passes validation.

### Phase 2: Patch the holes

Start with audits that can invalidate numbers, then tests, scan statistics, and documentation. The artifact work has two stages: an early pre-audit provenance snapshot, and a final blessed artifact manifest only after the invalidation gate below.

**Xfail repair protocol:** First capture the current failing behavior:

```bash
mkdir -p tests/baselines
pytest -ra --runxfail tests/test_quark_fit.py | tee tests/baselines/pre-fix-xfail-output.txt
```

Then preserve old benchmark outputs as a regression fixture so repaired PDG-target logic can be compared against the former benchmark surface. Next, make xfail handling strict either globally through pytest config or locally with strict markers, so XPASS does not hide regressions. Final acceptance is zero xfails on `tests/test_quark_fit.py`, followed by the targeted quark regression suite. The deliverable is not merely removing decorators; it is explaining or testing why the new benchmark is equivalent, intentionally changed, or superseded.

**Bag-parameter audit:** Build a table mapping each bag parameter, decay constant, meson mass, scale, scheme, and uncertainty used by `default_delta_f2_inputs`, `_kaon_matrix_elements`, and `_meson_matrix_elements` to a cited source. If any value changes, mark all affected bound summaries and plots stale.

**Wilson-RG audit:** Check `_evolve_wilsons` and the `evaluate_*_with_running` paths against a cited formula set. Record basis, normalization, RG direction, thresholds, matching scale, and whether the implementation is LO or NLO. Add one or two tests with analytically simple Wilson vectors.

**Final-claims invalidation gate:** Insert this gate between the bag/RG audits and Phase 2 finalization: no p50/p95/figure may be quoted in the final doc until bag-parameter audit + Wilson-RG audit have either certified "no change" or rerun the affected scans. A pre-audit manifest may exist, but the final blessed manifest must not be signed until this gate passes.

**CFW reproduction:** Start from the accessible source text, not from an access-risk assumption:

```bash
curl https://arxiv.org/pdf/0804.1954 -o /tmp/cfw.pdf && pdftotext /tmp/cfw.pdf - | head -200
```

The task is to reconcile conventions: Yukawa normalization, `g_s*`, c-parameter priors, observable definitions, historical numerical inputs, and any rescaling needed to compare to this repo. Deliver either a faithful reproduction script/table or narrowed language that calls the comparison qualitative.

**Zero-pass finite-statistics task:** Add a manifest section for every zero-pass claim, including moreUV, moreIR, and Run C factor-1.5. For each entry record N draws, accepted gate definition, seed or tile provenance, output hash, binomial 95% upper limit for the pass probability, and the final wording. Use language like "no passes observed in this finite ensemble" rather than "impossible" unless supported by an analytic proof.

**Prior and gate robustness:** Bootstrap the Run B p50/p95 values, document uncertainty, and run one log-uniform `|Y|` prior scan only if the existing robustness statement still needs it after uncertainty estimates. Add CLI tests for representative `scripts/run_rs_anarchy.py` paths.

**Scope and approximation wording:** Close EWPO, `g_s*`, KK tower, and lepton-sector language. State what is computed in this repo, what is cited externally, and what is a first-mode or quark-only approximation.

### Phase 3: Paper-readiness pass

Run the final targeted suite: all `tests/test_quark*.py`, `tests/test_mass_running.py`, `tests/test_pdg_quark_masses.py`, `tests/test_diagnostics.py`, and `tests/test_rs_anarchy_priors.py`; then run full `pytest -q` if wall time is acceptable. Rebuild `docs/quark_scan_methodology_note.pdf` and `docs/quark_scan_assumptions_compact.pdf` from a clean tree.

Prune figures to a named submission set. Keep exploratory plots ignored or store them in the external artifact snapshot. Do not commit the 25 GB `scan_outputs/` tree. Create a compressed artifact bundle or institutional-storage directory with checksums, and cite the location in the tracked final manifest. Tag only after tests, PDF rebuild, final manifest checks, and the bag/RG invalidation gate pass.

## D. Risks and mitigations

**Risk 1: Branch capture misses untracked provenance.** Probability: medium. Impact: high. Mitigation: status file, tracked diff, staged diff, untracked tarball, and git bundle before topic commit splitting.

**Risk 2: Bag or Wilson-RG audit changes numerical bounds.** Probability: medium. Impact: high. Mitigation: perform audits before final claims freeze; rerun or revalidate every affected `M_KK^min` p50/p95 and figure.

**Risk 3: Benchmark repair changes physics interpretation.** Probability: medium. Impact: high. Mitigation: capture `--runxfail` baseline, preserve old outputs as fixtures, use strict xfail handling, and require human sign-off for intentional benchmark changes.

**Risk 4: CFW comparison does not survive convention reconciliation.** Probability: medium. Impact: medium to high. Mitigation: separate the new PDG result from the historical comparison, and narrow the claim if exact reproduction fails.

**Risk 5: Zero-pass language overstates finite scans.** Probability: high. Impact: medium. Mitigation: record N, gate, seed, and binomial 95% upper limits, then soften language accordingly.

**Risk 6: Artifact blessing becomes ambiguous.** Probability: medium. Impact: high. Mitigation: make the final manifest a reviewed deliverable, not a mechanical hash dump.

## E. Final deliverable

The paper should live on `paper/quark-scan-2026q2`, with short-lived topic branches merged into it. `main` should receive only a reviewed merge from the paper branch after release-candidate validation. Use tags like `quarkscan-paper-v0.1-snapshot` for the first committed capture, `quarkscan-paper-rc1` for the first fully audited draft, and `quarkscan-paper-v1.0` for submission.

Ready-to-submit checklist:

1. `paper/quark-scan-2026q2` contains all code, docs, notebooks, and whitelisted figures needed for the quark scan.
2. No xfailed tests remain in `tests/test_quark_fit.py`, and the old benchmark surface is preserved in regression fixtures.
3. Bag parameters and Delta F=2 RG conventions are source-traceable in the methodology note.
4. Bag/RG audits certify no change, or every affected p50/p95/figure has been rerun or revalidated.
5. CFW comparison is either reproduced quantitatively or explicitly narrowed.
6. `g_s*`, EWPO floor, KK tower truncation, and lepton-sector scope are unambiguous.
7. Prior/gate robustness and zero-pass statements include uncertainty or conservative finite-ensemble wording.
8. A tracked final run manifest maps every paper claim to code SHA, command, output, checksum, and figure hash.
9. Submission figures are pruned, named, and regenerable from tracked scripts.
10. Clean-tree tests, PDF builds, release tag, and external scan-output snapshot are pushed and recorded.

===PLAN_V2_END===
