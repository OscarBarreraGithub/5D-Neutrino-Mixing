# ISSUES — Catalogue of problems found during merge+review

Format per issue:
  [<UNIT>-I<N>] severity: <CRITICAL|HIGH|MEDIUM|LOW|INFO> tag: <physics|numerics|code|docs|infra>
  Title: <one-line>
  Description: <what's wrong>
  Evidence: <file:line, commit, command output>
  Recommended fix: <action>

## Open

(none — all open issues have been routed to Closed / Accepted-risk subsections
below; Tier-6 deferred work tracked in `CLEANUP_QUEUE.md` rows C02b
(sensitivity-band figure) and C02c (SLURM RUNA reruns at 3 budget edges),
both deferred to paper-finalization phase.)

## Closed / Accepted-risk

### Closed by C01

- [R03-I1] severity:HIGH tag:code  **CLOSED 2026-05-25** by commit `3ab1f8f` (C01a).
  Title: Duplicate (pre-audit) kaon hadronic constants in modern backend not updated by Phase-2 hole #5 audit
  Description: `quarkConstraints/modern/phenomenology.py:421-431` carried a private copy of the kaon hadronic constants (`_KAON_B_1 = 0.717`, `_KAON_B_4 = 0.78`, `_KAON_B_5 = 0.57`, `_KAON_EPSILON_K_SM = 1.81e-3`) that were NOT updated when `dc9c498` updated the canonical values in `quarkConstraints/deltaf2.py:618-623` (`B_1_K=0.5503`, `B_4_K=0.903`, `B_5_K=0.691`, `EPSILON_K_SM=2.161e-3`). Closed via recommended-fix option (b): the four stale literals were updated in place (option (a) `from ..deltaf2 import` is structurally impossible — the modern lane is firewalled from `deltaf2` at 3 enforcement layers: `tests/test_modern_phenomenology.py:41-52,198-249` + `quarkConstraints/modern/verifier.py:50-60`). Collaborator artifacts re-exported under `7899205` (C01b); ratio_epsilon_K columns tightened ~7x as forecast. Full diagnosis in `.orchestration/cleanup_reports/C01.md`.

- [R03-I2] severity:MEDIUM tag:numerics  **CLOSED 2026-05-25** by commit `3ab1f8f` (C01a).
  Title: Audited-constant pin test does not guard the modern backend duplicates
  Description: `tests/test_quark_deltaf2.py::test_audited_deltaf2_hadronic_constants_match_selected_sources` (lines 111-117) only asserted the four post-audit values in `quarkConstraints.deltaf2` and did not cross-pin against `quarkConstraints.modern.phenomenology._KAON_*`. Closed by adding `test_modern_phenomenology_kaon_constants_match_deltaf2_canonical` in the same file. The new test imports both modules at test-function scope (legal — the firewall guards the modern source/sys.modules surface, not tests) and asserts `np.isclose` for all 11 kaon constant pairs with an explicit "Update both literals together" failure message.

### Closed by C02a-code

- [R03-I3] severity:MEDIUM tag:docs  **CLOSED 2026-05-25** by C02a-code (commit `3e1e07c`).
  Title: Five hole #5 follow-up tasks gating the final-claim invalidation gate remain open
  Description: `docs/phase_logs/phase2_h5_signoff.md:104-136` enumerates five follow-up tasks. Closure of this issue is **bookkeeping only**: it folds together one in-this-commit closure of tasks 1, 3, and 4 plus the deferred-but-tracked follow-ups for tasks 2 and 5.
    - Tasks 1 (methodology-note band quote), 3 (band-edge $\Mkk^{\min}$ quotation `47.26^{+75.10}_{-25.05}` TeV), and 4 (methodology-note paragraph describing the band construction) are **resolved in this commit** by the new `\paragraph{Symbolic NP-budget band on $\Mkk^{\min}$}` paragraph in `docs/quark_scan_methodology_note.tex` (after line 729, immediately following the existing BGS-budget band quote). The paragraph cites the symbolic $1/\sqrt{|\Delta\epsilon_K|}$ scaling, includes the explicit ratios $2.59$ (tight) and $0.47$ (loose) at the C02a plan's documented edges $\Delta\epsilon_K\in\{1\times10^{-5},\,6.7\times10^{-5},\,3\times10^{-4}\}$, and cross-links the `--epsilon-k-budget` CLI seam.
    - Task 2 (three RUNA reruns) is **deferred to cleanup unit C02c** (SLURM-budget reason — see CLEANUP_PLAN.md C02c). The CLI seam (`scripts/run_rs_anarchy.py --epsilon-k-budget {central,low,high}`) plumbed end-to-end through `EnsembleConfig` → `_worker_init` → `evaluate_delta_f2_constraints` → `evaluate_epsilon_k` IS landed in this commit; smoke runs confirm `central` is bit-for-bit identical to the pre-flag path and `low`/`high` rescale ε_K ratios by exactly $6.7$ / $0.2233$ as predicted.
    - Task 5 (sensitivity-band figure) is **deferred to cleanup unit C02b** (paper finalization).
  Evidence: `docs/quark_scan_methodology_note.tex` (new `\paragraph{Symbolic NP-budget band on $\Mkk^{\min}$}` after line 729); `scripts/run_rs_anarchy.py` (CLI flag + `EnsembleConfig` + `_worker_init` round-trip); `quarkConstraints/deltaf2.py` (`epsilon_k_np_budget_override` keyword); `tests/test_run_rs_anarchy.py` (12 new tests, all pass); 125/125 broader regression tests pass; PDF rebuilt cleanly (20 pp). Full diagnosis in `.orchestration/cleanup_reports/C02a-code.md`.

### Closed by C03

- [R04-I1] severity:LOW tag:numerics  **CLOSED 2026-05-25** by C03.
  Title: Audit-script closed-form shortcut silently assumes upper-triangular LR ADM
  Description: `scripts/audit_wilson_rg.py::scalar_lr_segment_matrix` set the lower-left ADM entry to `0.0` literal-by-literal but did not assert the upper-triangular precondition. Closed by hoisting `gamma54 = 0.0` into a named local, building `gamma_lr_local = np.array([[gamma44, gamma45], [gamma54, gamma55]])`, and raising `ValueError("audit_wilson_rg.scalar_lr_segment_matrix requires an upper-triangular scalar LR ADM (gamma_54 == 0); update the closed-form propagator if the anomalous-dimension block is changed.")` when `gamma_lr_local[1, 0] != 0.0` or `np.allclose(gamma_lr_local - np.triu(gamma_lr_local), 0.0)` fails. Numeric output unchanged (`C4_LR=1 -> 3.53816397486`, max relative discrepancy `1.300e-16`). Full diagnosis in `.orchestration/cleanup_reports/C03.md`.

- [R04-I2] severity:LOW tag:code  **CLOSED 2026-05-25** by C03.
  Title: `_evolve_wilsons` performs a function-local import
  Description: The function-local `from .qcd_running import evolve_deltaf2_wilsons` in `quarkConstraints/deltaf2.py::_evolve_wilsons` (line 606 pre-cleanup) has been lifted to the module-top import block. Reviewer M-7 confirmed `quarkConstraints/qcd_running.py` does not import `deltaf2` at module scope (the module only imports `math` and `numpy`); the lift is therefore unconditional and non-circular. The orchestrator prompt's R04-I2 wording flipped the file (claimed the offender was in `qcd_running.py`); the plan (CLEANUP_PLAN.md line 46) and source agree the offender is in `deltaf2.py`. Closed by editing the module-top import block and replacing the function-local import line with a provenance docstring note. Focused suite `pytest -k 'wilson or rg or qcd or audit or bmu'` (103 passed) and broader sweep `pytest -k 'deltaf2 or epsilon_k or kaon or matching or modern'` (196 passed, 1 skipped) both green.

- [R04-I3] severity:INFO tag:physics  **CLOSED 2026-05-25** by C03 (verification-only).
  Title: `paper_0710_1869` LR running module carries an independent BMU map — confirm sign-chain alignment in R06
  Description: C03 cross-checked the LR-sign chain in `quarkConstraints/paper_0710_1869/eft_deltaf2/{matching_kkgluon,hadronic,observables,rg}.py` against the post-C01 audited `quarkConstraints/deltaf2.py` path and found **no mismatch**. Both modules use the conventional positive scalar O4/O5 sign with `C4_LR < 0`, `C5_LR > 0` at matching, positive bag inputs, positive matrix-element prefactors, and no extra contraction sign. The BMU LR LO ADM is stored as the Wilson-coefficient transpose `((2,12),(0,-16))` in `paper_0710_1869/eft_deltaf2/rg.py:98-101` and as the operator form `[[2,0],[12,-16]]` in `qcd_running.py:50-57` / `scripts/audit_wilson_rg.py:85-117`, which is the standard convention difference between Wilson-coefficient ADM and operator ADM (transposes). The `paper_0710_1869` default LR bag values (`B4=0.78`, `B5=0.57` from ETM 2013) are intentionally distinct from the FLAG 2024 values in the repo-v1 path (`B_4_K=0.903`, `B_5_K=0.691`); the bundle identifier `hadronic.kaon.lr.b4.etm2013.table1.ms_2gev.v1` separates them by design. No code change required. Full file:line evidence table in `.orchestration/cleanup_reports/C03.md`.

- [R04-I4] severity:INFO tag:docs  **CLOSED 2026-05-25** by C03 (docs-only tracking note).
  Title: Four R04 deferred follow-ups remain open (endpoint migration, NLO BMU, VLL/VRR NLO, methodology-note tables)
  Description: A new `## Status as of 2026-05-25 post-cleanup (C03)` section has been appended to `docs/audits/wilson_rg_inventory.md` (immediately before the existing `## Invalidation-gate note`) explicitly enumerating the four deferred follow-ups from the inventory's 10-item checklist: (item 2) Hamiltonian/CFW comparison; (item 4) per-system RG endpoint migration; (item 7) full NLO WET evolution module; (item 10) endpoint caveat for matrix-element contraction. The section also records the R04-I3 verification, the R04-I1/I2 code changes, and pins the audit-script numeric output. Per CLEANUP_PLAN.md §C C19, the substantive resolution of these methodology-note items is reassigned to **C19**; the C03 doc addendum is the bridge. (CLEANUP_QUEUE.md row for C20 also lists R04-I4 as an ACCEPTED-RISK candidate; the C19/C20 split is not litigated here.)

### Closed by C04

- [R07-I2] severity:LOW tag:numerics  **CLOSED 2026-05-25** by C04.
  Title: `wilson_upper_limit` unit tests only cover k=0 cases
  Description: `tests/test_finite_stats.py` previously exercised only `k=0` for two values of n (1000 and 1_600_000) plus three ValueError cases. The general-k formula at `quarkConstraints/finite_stats.py:24-26` was verified algebraically but no regression guarded against a future k>0 bug. Closed by adding three new test groups: (a) `test_wilson_upper_limit_k_gt_zero_matches_closed_form` parametrized over `(k, n) ∈ {(1,100), (1,1000), (5,1000), (10,1000), (50,1000)}`, cross-checking the helper against an in-file closed-form reference at z=1.92 (rel 1e-12); (b) `test_wilson_upper_limit_k_gt_zero_matches_scipy_at_z_eq_1p96` cross-checking against `scipy.stats.binomtest(k, n).proportion_ci(method='wilson', confidence_level=0.95).high` for `(1,100)`, `(5,1000)`, `(50,1000)` — the scipy oracle uses `z = norm.ppf(0.975) = 1.959964...` which is passed to the repo helper for the comparison; (c) `test_wilson_upper_limit_k_gt_zero_spot_check_k5_n1000` pinning the explicit numerical value at k=5,n=1000,z=1.92 to ~0.01146 (the dispatch-prompt hand-computation of 0.0089 was off by the `+z*sqrt(...)` margin term — see CLEANUP_PLAN.md §C C04 for the closed-form). 12/12 finite_stats tests pass, 48/48 in the focused selection. scipy>=1.10 is satisfied (in-repo scipy is 1.17.0); the closed-form group does not depend on scipy.

- [R22-I1] severity:LOW tag:tests  **CLOSED 2026-05-25** by C04.
  Title: `notation.ts` / `prose.ts` ship without automated regression tests
  Description: Closed by adding a vitest snapshot suite under `flavor_catalog/website/src/lib/__tests__/{notation,prose}.test.ts` covering the 5-entry corpus from the C04 plan (T010, CR002, K018, B015, K020 raw `process_name` fields). `notation.test.ts` snapshots all four public helpers (`normalizeNotation`, `notationToMath`, `processNameToHtml`, `processNameToPlain`) over the 5 entries (20 snapshots). `prose.test.ts` snapshots `texProseToHtml` over the same 5 entries plus 4 prose fragments that exercise the LaTeX text-mode rewriter, em-dash normalization, jargon-scrub, and shorthand-wrapping branches (9 snapshots). Snapshot files committed under `flavor_catalog/website/src/lib/__tests__/__snapshots__/{notation,prose}.test.ts.snap`. `vitest@^4.1.7` added to `devDependencies`; `npm test` and `npm run test:update` scripts added to `package.json`. `npm test` reports `Test Files 2 passed (2), Tests 29 passed (29)`.

### Closed by C05

- [R17-I1] severity:LOW tag:provenance  **CLOSED 2026-05-25** by C05.
  Title: External-research PDFs lack an explicit sha256 / MANIFEST / README in `flavor_catalog/external_research/`
  Description: Closed by adding `flavor_catalog/external_research/MANIFEST.md` (provenance log distinguishing imported artifacts vs in-house review memos, with per-file source/subject/date/sha256/description tables) and `flavor_catalog/external_research/sha256sums.txt` (six rows, validated with `sha256sum -c sha256sums.txt` → 6/6 OK). Import dates were recovered from `git log` (`022a20c` PDFs/txts 2026-05-16 18:05 EDT; `8fb5f91` may16 review 18:13; `2c00d84` may15 review 18:15). Sub-item (a) PDF binary digests — fully closed by `sha256sums.txt`. Sub-items (b) session URLs / prompts and (c) text-extraction recipe — partially closed: not retroactively recoverable; MANIFEST documents the gap and records a forward-looking convention for future external-research imports. No catalog numerical claim depends on these artifacts, so no code change or test re-run was required. Evidence in `.orchestration/cleanup_reports/C05.md`.

### Closed by C06

- [R07-I3] severity:INFO tag:infra  **CLOSED 2026-05-25** by C06.
  Title: Manifest `code_sha_at_run_time` is inferred from directory mtime, not embedded in scan output
  Description: Closed by embedding `git_sha` in the top-level `tile_summary.json` payload written by `scripts/run_rs_anarchy.py` — the sole producer of `tile_summary.json` in the repo. A new module-level helper `_resolve_git_sha()` calls `subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=_REPO_ROOT, stderr=DEVNULL)` and returns the documented `"unknown"` sentinel on `(CalledProcessError, FileNotFoundError, OSError)` or empty output (matching the manifest's existing inference-fallback convention). Fully backward compatible: the four current `tile_summary.json` readers (`scripts/{plot_rs_anarchy_summary,rs_anarchy_cfw_comparison,rs_anarchy_gate_sensitivity,rs_anarchy_mkk_min_hist_by_cvals}.py`) use plain `json.load` and pull specific keys; adding a new top-level field is invisible to them, and consumers that want to verify against a manifest can use `payload.get("git_sha", "unknown")`. Four new tests in `tests/test_run_rs_anarchy.py` cover: in-repo SHA shape, out-of-tree `"unknown"` fallback, end-to-end `main()` smoke writes the field, and legacy payloads without the field still load. Targeted: 4/4 pass in 5.42s. Broader sweep `pytest -k 'tile or summary or scan or dispatch'` 55 passed, 515 deselected, 122.46s. Evidence in `.orchestration/cleanup_reports/C06.md`.

### Closed by C07

- [R17-I2] severity:LOW tag:docs  **CLOSED 2026-05-25** by C07.
  Title: `master_compile_v02_report.md` claims "79 VERIFIED of 80" but `factcheck_status.md` at R17 endpoint only documents 75 rows
  Description: Closed by appending a "Consolidation status" block to `flavor_catalog/master_compile_v02_report.md` that reconciles the headline `79 VERIFIED + 1 PARTIAL out of 80` against the per-family audits at the v0.2 boundary. At the v0.2 tag, the consolidated `audits/factcheck_status.md` table has 75 rows (`74 VERIFIED + 1 PARTIAL`); the remaining 5 Wave-7 PRIMARY verdicts live in `audits/factcheck_top_higgs_ew.md:294-352` (T003/T004/T008/T012, all VERIFIED) and `audits/factcheck_beauty.md:92` (B012, VERIFIED). Sum: 75 + 4 + 1 = 80 = 79 VERIFIED + 1 PARTIAL. The annotation cites both addendum paths and records the forward-looking convention (regenerate the consolidated table at each master-compile bump, optionally via `tools/aggregate_factchecks.py`). Evidence in `.orchestration/cleanup_reports/C07.md`.

- [R18-I3] severity:LOW tag:docs  **CLOSED 2026-05-25** by C07.
  Title: Consolidated `factcheck_status.md` was NOT regenerated for Wave-8 — SECONDARY rows live only in per-family addenda (recurrence of R17-I2 at the SECONDARY tier)
  Description: Closed by appending a "Consolidation status" block to `flavor_catalog/master_compile_v03_report.md` that reconciles the headline `86 VERIFIED + 2 PARTIAL out of 88` against the per-family audits at the v0.3 boundary. At the v0.3 tag, the consolidated table still has 75 rows (unchanged from v0.1). The remaining 13 verdicts are: 5 Wave-7 PRIMARY (T003/T004/T008/T012/B012, all VERIFIED) and 8 Wave-8 SECONDARY in `audits/factcheck_kaon.md:254-313` (K019/K020/K021), `audits/factcheck_beauty.md:243-310` (B007/B008/B013/B014), `audits/factcheck_top_higgs_ew.md:391-427` (T014) — 7 VERIFIED + 1 PARTIAL at the v0.3 tag (K020 PARTIAL, cleared post-tag by `b5c2375`). Sum: 75 + 5 + 8 = 88 = 86 VERIFIED + 2 PARTIAL. The block also notes the orthogonal R18-I1 v0.3 tag-annotation discrepancy. Evidence in `.orchestration/cleanup_reports/C07.md`.

- [R19-I4] severity:LOW tag:docs  **CLOSED 2026-05-25** by C07.
  Title: Consolidated `factcheck_status.md` was NOT regenerated at the R19 endpoint to include the 14 Wave-9 collider_rs rows (recurrence of R17-I2 / R18-I3 at the new-family tier)
  Description: Closed by appending a "Consolidation status" block to `flavor_catalog/master_compile_v04_report.md` that reconciles the headline `101 VERIFIED + 1 PARTIAL out of 102` against the per-family audits at the v0.4 boundary. At the v0.4 tag the consolidated table has 89 rows: the original 75 DA-4 base + 14 Wave-9 collider_rs (`CR001`-`CR014`, all VERIFIED, folded in by commit `0c5dacc` at master-compile time). The remaining 13 verdicts are still the Wave-7 PRIMARY (5 VERIFIED) + Wave-8 SECONDARY (8 VERIFIED, K020 cleared by `b5c2375` before the v0.4 tag) addenda. Sum: 89 + 5 + 8 = 102 = 101 VERIFIED + 1 PARTIAL. The original ISSUES.md description framed this as "CR rows missing from consolidated table at v0.4"; in fact the CR rows ARE in the consolidated table at the v0.4 tag, and the still-uncounted entries are Wave-7/8. The annotation also notes the orthogonal R19-I3 v0.4 tag-annotation discrepancy. Closure also lands a forward-looking helper `tools/aggregate_factchecks.py` that scans `factcheck_status.md` + per-family addenda (both table-row and section-header verdict forms) and emits per-family + overall V/P/M/F counts; at HEAD it reports `102 processes / 101 VERIFIED / 1 PARTIAL`, matching `master_compile_v04_report.md`. Evidence in `.orchestration/cleanup_reports/C07.md`.

### Closed by C08

- [R10a-I2] severity:INFO tag:docs  **CLOSED 2026-05-25** by C08.
  Title: K001 `open_issues` flags the same legacy/modern epsilon_K constant mismatch tracked as R03-I1
  Description: Closed by prepending resolution markers to both K001 `open_issues` items in `flavor_catalog/processes/kaon/K001.yaml`. Item 1 (PDG review vs pdgLive choice) is marked `[RESOLVED-by-C08 2026-05-25]` with rationale citing the PDG 2026 review extract at `flavor_catalog/references/K001/pdg2026_epsilon_k.txt` as the canonical snapshot path (encoded in `pdg_or_equivalent.snapshot_path`). Item 2 (legacy/modern epsilon_K constant mismatch) is marked `[RESOLVED-by-C08 2026-05-25, tracked-in R03-I1, closed-by cleanup-C01 commit 3ab1f8f]` so the cross-link is now explicit. Original wording is preserved verbatim after `Original note: ` so the audit trail is intact. Evidence in `.orchestration/cleanup_reports/C08.md`.

- [R10a-I3] severity:INFO tag:docs  **CLOSED 2026-05-25** by C08.
  Title: K003 `open_issues` records an unresolved canonical choice for eps'/eps display value
  Description: Closed by prepending `[RESOLVED-by-C08 2026-05-25]` to the K003 `open_issues` entry in `flavor_catalog/processes/kaon/K003.yaml`. The canonical headline is the PDG listing value `Re(eps'/eps) = 0.00166 +/- 0.00023` (now in `pdg_or_equivalent.value/display_value` at K003.yaml:66,69); the S013EPS DataBlock OUR AVERAGE 1.68 +/- 0.20 x 10^-3 is retained as `supporting_value` (K003.yaml:87,90). Rationale: the listing-fit row is what downstream consumers will quote; CA Wave-1 cycle accepted this assignment implicitly when CHK-1..CHK-8 passed. Evidence in `.orchestration/cleanup_reports/C08.md`.

- [R10b-I3] severity:INFO tag:docs  **CLOSED 2026-05-25** by C08.
  Title: All 5 Wave-1 beauty yamls retain unresolved `open_issues` entries flagged by PKA/WA for CA disposition that were never formally closed
  Description: Closed by prepending resolution markers to the `open_issues` entries in all 5 Wave-1 beauty yamls. B002 (`[RESOLVED-by-C08]` HFLAV all-charmonium 0.710 canonical); B005 (`[RESOLVED-by-C08]` PDG live 3.34e-9 canonical, HFLAV Apr-2023 demoted to `hflav_historical_average`); B009 item 1 (`[RESOLVED-by-C08]` HFLAV Dec-2025 1.12e-4 canonical) + item 2 (`[DEFERRED-to-post-paper-finalization]` direct f_B|V_ub| construction not needed for operational_scan_only claim level); B011 (`[RESOLVED-by-C08]` Misiak 2020 3.40e-4 canonical SM, 2015 PDG-quote retained as `sm_prediction_pdg_review_quote`); B015 item 1 (`[CLOSED-as-N/A]` rejected-anchor provenance note, correct BaBar/Belle 1312.5364 + hep-ex/0503044 already cited) + item 2 (`[CLOSED-as-by-design]` exclusive B->K(*)ll and LFU details deferred to B016-B019 by PKA design). Original wording preserved verbatim after `Original note: ` in each. Evidence in `.orchestration/cleanup_reports/C08.md`.

### Closed by C09

- [R10b-I1] severity:INFO tag:docs  **CLOSED 2026-05-25** by C09.
  Title: B011 yaml `sm_prediction_current.observable` field carries a typo "SM prediction for B_s gamma" (should be `B -> X_s gamma`)
  Description: Closed by editing `flavor_catalog/processes/beauty/B011.yaml:94` from `observable: "SM prediction for B_s gamma"` to `observable: "SM prediction for B -> X_s gamma"`. Matches the ASCII-arrow notation used elsewhere in the same file (e.g. line 79 `"CP- and isospin-averaged branching fraction B(B -> X_s gamma)"`). Value `0.000340`, uncertainty `0.000017`, units, conditions (`E_gamma > 1.6 GeV in the decaying meson rest frame`), source (`Misiak, Rehman, Steinhauser 2020`), source_url, snapshot_path, and sha256 are all unchanged — this is a human-readable label fix only. Evidence in `.orchestration/cleanup_reports/C09.md`.

- [R10b-I2] severity:INFO tag:docs  **CLOSED 2026-05-25** by C09.
  Title: B009 yaml status_history entry for `WRITER-INITIATED -> PKA-DONE` omits the `state:` field
  Description: Closed by prepending `state: "PKA-DONE"` at `flavor_catalog/processes/beauty/B009.yaml:19`. Now matches the pattern of all six other entries in the same `status_history` list (each carries `state:` as the leading key, with `to:` echoing the same value). The canonical state was already in the `to:` field, so this is presentational alignment only; no audit information changes. Evidence in `.orchestration/cleanup_reports/C09.md`.

- [R11-I3] severity:INFO tag:docs  **CLOSED 2026-05-25** by C09.
  Title: B025 records a no-PDF-policy gap for the Belle II CKM-2025 hadronic-tag input
  Description: Closed by (a) prepending a `[DOCUMENTED-by-C09 2026-05-25]` resolution marker to the B025 `open_issues` entry at `flavor_catalog/processes/beauty/B025.yaml:233` explaining that the input is captured indirectly via the HFLAV joint-fit text snapshots and (b) adding a new citation-anchor entry at `flavor_catalog/website/_data/citation_anchors/B025.yaml` (block_key `canonical_average.belleII_hadronic_tag_input_no_pdf_policy_gap`) that uses a new `unsnapshotted_reason:` multi-line scalar and `status: UNSNAPSHOTTED-BY-POLICY` value to formally document the no-PDF-policy gap. The original wording is preserved verbatim after `Original note: ` in the B025 yaml so the audit trail is intact. Original recommended fix (snapshot the conference-talk PDF) is unchanged in priority — backfill the eventual arXiv/journal version of the hadronic-tag measurement when it appears; no urgency at the current `operational_scan_only` claim level. The new `unsnapshotted_reason:` field and `UNSNAPSHOTTED-BY-POLICY` status set a forward-looking convention for documenting future no-PDF-policy gaps. Evidence in `.orchestration/cleanup_reports/C09.md`.

### Closed by C10

- [R10c-I1] severity:INFO tag:docs  **CLOSED 2026-05-25** by C10.
  Title: T010 PKA assignment combines plan rows T010 (R_b) and T011 (A_FB^b, A_b) into a single Z bbar b pole package; T011 has no separate file
  Description: Closed via brief option (b) — created a non-intrusive stub `flavor_catalog/processes/top_higgs_ew/T011.yaml` with `schema: flavor_catalog.process.v1`, `merged_into: T010`, `canonical_home` block pointing at T010.{yaml,tex}, a one-entry `status_history` recording `CLOSED-AS-MERGED` by cleanup-C10, and a `[CLOSED-as-merged-by-C10 2026-05-25]` `open_issues` bullet. No `T011.tex` is added — `flavor_catalog/processes/top_higgs_ew/index.tex` uses `\IfFileExists`, and a hollow .tex stub would render an empty section in the master compile; the merge is already documented in `T010.tex:14-17` and in `worklogs/discovery/round_004_addendum_deferred_scope.md:25`. The plan-row -> file mapping is now explicit. Evidence in `.orchestration/cleanup_reports/C10.md`.

- [R10c-I2] severity:INFO tag:docs  **CLOSED 2026-05-25** by C10.
  Title: Reference-snapshot directories are inconsistent across the 6 R10c PKAs — T002, T010, E001 lack a top-level `sha256sums.txt` listing
  Description: Closed by generating `sha256sums.txt` in each of the three referenced dirs using `(cd flavor_catalog/references/<PID> && sha256sum *.txt source_manifest.yaml > sha256sums.txt)`: T002/ (7 entries: 6 .txt + manifest), T010/ (7 entries), E001/ (6 entries). All three validate cleanly under `sha256sum -c sha256sums.txt`. Format matches the existing files in T001/C001/L001/B011/B012/..., so the six R10c PKAs now carry uniform top-level checksum coverage. Evidence in `.orchestration/cleanup_reports/C10.md`.

- [R10c-I3] severity:INFO tag:docs  **CLOSED 2026-05-25** by C10.
  Title: Four R10c PKAs (T002, E001, C001, L001) record their WA/CA cycles in `status_history` but not in any wave-1-named worklog
  Description: Verified C08 did not address this (C08 closed R10a-I2/I3 and R10b-I3 open_issues threads). Closed by adding a `wave_assignment: "wave1"` scalar plus a `worklog_alias:` block (with `family_batch`, `writer_worklogs[]`, `checker_worklogs[]`) to the four affected yamls: T002 (w23_top_higgs_ew, 2 cycles), E001 (w23_kaon_charm_edm, 2 cycles), C001 (w23_kaon_charm_edm, 2 cycles), L001 (w23_charged_lepton, 3 cycles + `opus_arbitration_id`). A brief comment block above the new fields explains that the PKA is Wave-1 per plan v1 but was polished under the sibling Wave-2/3 batch because the wave-1 batch was already in flight. Worklog filenames themselves are left as-is (they correctly reflect the batch that ran). All five touched yamls parse cleanly via `yaml.safe_load` with schema string `flavor_catalog.process.v1` preserved. Evidence in `.orchestration/cleanup_reports/C10.md`.

### Closed by C11

- [R11-I2] severity:INFO tag:docs  **CLOSED 2026-05-26** by C11.
  Title: DA-1 worklog records four PI escalations without a resolution thread
  Description: Closed by appending a "Closure Addendum" section to `flavor_catalog/worklogs/discovery/round_001_full_scope.md` that records a per-escalation resolution table for the four PI items at lines 55-60: (a) K013 alone vs K014 alternate → RESOLVED via DA-4 Wave-7 addendum (K014 = DEFERRED-SCOPE pending PI override, `round_004_addendum_deferred_scope.md:43`); (b) EW002 standalone CKM-unitarity → RESOLVED via Wave-4 PKA + WA-v2/CA-v2 chain (`wa_w4_ew_v2.md`, `ca_w4_ew_v2.md`), now present at `flavor_catalog/processes/top_higgs_ew/EW002.{yaml,tex}`; (c) EW003 single overview vs B029/B030/B031 split → RESOLVED via Wave-4 (single overview retained at `EW003.{yaml,tex}`; B029/B030/B031 folded into EW003 per DA-4 deferred-scope addendum line 38); (d) E004/E006/E008 catalog-only vs future hard-cuts → PARTIALLY RESOLVED (drafted as catalog rows in Wave-4 via `wa_w4_kaon_edm{,_v2}.md` and `ca_w4_kaon_edm{,_v2}.md`; hard-cut intent remains a synthesis-stage decision, not a DA-side item). DOC CLOSED marker added with cross-link to DA-2 and DA-4 closure addenda. Evidence in `.orchestration/cleanup_reports/C11.md`.

- [R12-I4] severity:INFO tag:docs  **CLOSED 2026-05-26** by C11.
  Title: DA-1 PI escalations from R11-I2 partially resolved by Wave-4 but not back-annotated
  Description: Closed jointly with R11-I2 by the same DA-1 closure addendum. The addendum explicitly records EW002 as RESOLVED-by-Wave-4 with PKA + WA-v2/CA-v2 commit references (incl. `52acd5e` WA-v2, `d6e78b6` CA-v2), EW003 as RESOLVED-by-Wave-4 (single overview retained), and E004/E006/E008 as drafted in Wave-4 with the hard-cut intent left to synthesis. A symmetric Closure Addendum was also appended to `round_002_followup.md` (DA-2) capturing the K013/K014 carry-forward → RESOLVED via DA-4 Wave-7 addendum, and to `round_003_final_sweep.md` (DA-3) for symmetry. The four DA-N worklogs now all carry explicit DOC CLOSED status. Evidence in `.orchestration/cleanup_reports/C11.md`.

- [R15-I4] severity:INFO tag:docs  **CLOSED 2026-05-26** by C11.
  Title: DA-4 worklog has no downstream addendum back-linking the round-2 closure of the 25 then-pending OPUS approvals
  Description: Closed by appending a "Round-2 Closure Addendum" to `flavor_catalog/worklogs/discovery/round_004_convergence.md` recording: 21 of 25 pending sidecars approved via `flavor_catalog/signoff/round_002_index.md`; B001+B003 approved via `flavor_catalog/signoff/by_process/B001_B003.md`; B021+B023 round-2 arbitration via `flavor_catalog/signoff/by_process/B021_B023.md`; Wave-7 promotions raising catalog from 75 → 80; `catalog_master.tex` compile completed with consolidated v0.4 fact-check report at `flavor_catalog/master_compile_v04_report.md` (102 processes / 101 VERIFIED / 1 PARTIAL). DOC CLOSED. Evidence in `.orchestration/cleanup_reports/C11.md`.

### Closed by C12

- [R18-I1] severity:LOW tag:docs  **CLOSED 2026-05-26** by C12.
  Title: `flavor-catalog-v0.3` tag annotation says "87 VERIFIED + 1 PARTIAL" but compile report and runbook §8 say "86V + 2P"
  Description: Closed by appending a "Tag annotation provenance" section to `flavor_catalog/audits/factcheck_status.md` that records, in a 3-row table, the V/P count written into each `flavor-catalog-v0.{2,3,4}` tag annotation alongside the canonical post-consolidation V/P from `master_compile_v0X_report.md` and a one-sentence drift-origin note. The v0.3 row pins the canonical 86V + 2P (88 entries; tag annotation `87V+1P` projects the post-tag K020 fix at commit `b5c2375` — at the v0.3-tag boundary K020 was still PARTIAL alongside E009). Tag annotations are not rewritten (immutable history). The new provenance section also lays down a forward-looking convention binding for v0.5+ tags: the V/P count in a master-compile tag annotation MUST be sourced from the at-tag-time compile report only; post-tag corrections never flow backward into the tag message. Evidence in `.orchestration/cleanup_reports/C12.md`.

- [R19-I3] severity:LOW tag:docs  **CLOSED 2026-05-26** by C12.
  Title: `flavor-catalog-v0.4` tag annotation says "100V + 2P" but the v0.4 compile report says "101V + 1P"
  Description: Closed jointly with R18-I1 by the same `factcheck_status.md` §"Tag annotation provenance" addendum. The v0.4 row pins the canonical 101V + 1P (102 entries: 89 consolidated-table rows after Wave-9 fold-in at `0c5dacc` + 5 Wave-7 PRIMARY addenda + 8 Wave-8 SECONDARY addenda). The tag annotation `100V+2P` projects the v0.3-tagged historical K020 PARTIAL state; at the v0.4-tag boundary K020 was already VERIFIED (cleared by `b5c2375` before the v0.4 tag was placed) and only E009 remained PARTIAL. Tag is not rewritten. Cross-link to canonical sources (the C07 "Consolidation status" block in `master_compile_v04_report.md` + the `tools/aggregate_factchecks.py` HEAD output) is encoded in the new provenance section. Evidence in `.orchestration/cleanup_reports/C12.md`.

- [R21-I1] severity:LOW tag:docs  **CLOSED 2026-05-26** by C12.
  Title: `WEBSITE_RUNBOOK.md:25` "100 VERIFIED / 2 PARTIAL" tally is the v0.4 tag-annotation state; the website itself ships post-tag state (101V / 1P)
  Description: Closed by editing `flavor_catalog/website/WEBSITE_RUNBOOK.md:21` to read `101 VERIFIED / 1 PARTIAL / 0 MISMATCH (post-v0.4-tag canonical; K020 cleared by b5c2375 before the v0.4 tag — matches master_compile_v04_report.md §"Consolidation status" and the website's catalog_index.json verdict_counts)`. The runbook headline now matches what the website actually renders (`verdict_counts: {VERIFIED: 101, PARTIAL: 1}`) and aligns with the canonical row in the new `factcheck_status.md` §"Tag annotation provenance" section. No website rebuild required (the runbook is operator documentation, not ingested by the build). Evidence in `.orchestration/cleanup_reports/C12.md`.

- [R22-I3] severity:INFO tag:ux  **CLOSED 2026-05-26** by C12.
  Title: `EntryTable` still emits `<th>Tier</th>` and `data-tier` after phase-10 removes the tier badge UI
  Description: Closed by removing the dead `data-tier={r.tier}` attribute from each `<tr class="row-link">` in `flavor_catalog/website/src/components/EntryTable.astro:70`. Prior cleanup (between `7053fb7` and the C12-state HEAD) had already dropped the `<th data-sort="tier">Tier</th>` column from the thead and the dead `tier` branch in `browse.astro` (no surviving `getAttribute('data-tier')` call remains); the only residue was the unused `data-tier` row attribute, which is now gone. The TypeScript `tier: string` field in the `IndexRow` interface is preserved because `ingest_catalog.py` still emits the field in the JSON payload; dropping it from the interface would force a wider refactor with no user-visible benefit. No JavaScript reads `data-tier` after this edit (verified by `grep -n 'data-tier' flavor_catalog/website/` → empty); no test re-run required. Evidence in `.orchestration/cleanup_reports/C12.md`.

### Closed by C13

- [R02-I2] severity:INFO tag:code  **CLOSED 2026-05-26** by C13.
  Title: Re-derived spurion seed values would benefit from an in-source provenance pointer
  Description: Closed by adding two in-source provenance pointers inside `quarkConstraints/benchmarks.py::default_spurion_seed` (line 189): (a) an extended multi-line docstring stating the singular-value and rotation literals were re-derived against PDG-2024 MS-bar quark masses with a pointer to `docs/phase_logs/phase2_h4_impl.md` for re-derivation details + audit log, and a note that `tests/test_quark_fit.py` pins the values so drift surfaces as a test failure; (b) a 3-line `#`-comment immediately above the `return SpurionSeed(` literal block carrying the short form requested by the dispatch brief (`# Spurion seed values re-derived against PDG-2024 MS-bar quark masses / (see docs/phase_logs/phase2_h4_impl.md for re-derivation details and / the audit log). Tests at tests/test_quark_fit.py pin these values.`). Both pointers live inside the function body so grep/AST-jump lands on them. No literal touched: `python -c "from quarkConstraints.benchmarks import default_spurion_seed; s = default_spurion_seed(); print(s.up_singular_values, s.down_singular_values, s.overall_scale)"` reproduces the pre-edit values to machine precision (`[0.1434281 0.33368792 1.23945796]`, `[0.22903978 0.16350921 0.28800013]`, `2.8`). Docs file `docs/phase_logs/phase2_h4_impl.md` confirmed to exist; test file `tests/test_quark_fit.py` confirmed to import and reference `default_spurion_seed` (lines 38, 55, 137, 150). No test re-run required (pure-comment change). Evidence in `.orchestration/cleanup_reports/C13.md`.

### Closed by C14

- [R14-I3] severity:LOW tag:docs  **CLOSED 2026-05-26** by C14.
  Title: Codex quota pause note cites gpt-5.5 but global config and CLAUDE.md use gpt-5.4
  Description: Closed by replacing `gpt-5.5` → `gpt-5.4` in `docs/phase_logs/flavor_catalog_codex_quota_pause.md:46` (the explicit R14-I3 target) and aligning eight other prescriptive surfaces with the canonical model name set in `~/.codex/config.toml` (per `~/.claude/CLAUDE.md` Agent Orchestration section): `flavor_catalog/AGENTIC_WORKFLOW.md:42-46`, `flavor_catalog/HANDOFF_PROMPT.md:103`, `flavor_catalog/WEBSITE_BUILD_PROMPT.md:25,263`, `flavor_catalog/SESSION_NOTES.md:201,227`, `flavor_catalog/website/README.md:55,211`, `flavor_catalog/audits/factcheck_status.md:7`, `docs/paper_execution_decisions.md:45`, `docs/phase_logs/POST_COMPACTION_BRIEFING.md:28,59,125`. Historical/provenance contexts preserved verbatim: `flavor_catalog/website/WEBSITE_RUNBOOK.md` (completed-event log of Phase-2 anchor batches), `flavor_catalog/website/_data/priority/*.yaml` (`ranked_by` provenance stamps, ~99 entries), `docs/phase_logs/flavor_catalog_plan_v0.md` + `_plan_v1.md` (superseded plan drafts), and the R14/R20 review records that document the drift by definition. Pure-docs cosmetic alignment; no code change; no test re-run required. Evidence in `.orchestration/cleanup_reports/C14.md`.

- [R20-I2] severity:INFO tag:docs  **CLOSED 2026-05-26** by C14.
  Title: AGENTIC_WORKFLOW.md role table cites "Codex GPT-5.5 xhigh" model; user's global CLAUDE.md cites "gpt-5.4"
  Description: Closed by the same nine-file edit batch as R14-I3 above. The two R20-I2-named targets (`AGENTIC_WORKFLOW.md:42-46` role table — 5 rows for PKA/WA/CA/DA/Fact-Check — and `HANDOFF_PROMPT.md:103` non-negotiable rule) now read "Codex GPT-5.4 xhigh", aligning with the canonical `gpt-5.4` set in `~/.codex/config.toml`. Evidence in `.orchestration/cleanup_reports/C14.md`.

### Closed by C15

- [R05-I2] severity:LOW tag:infra  **CLOSED 2026-05-26** by C15.
  Title: `.orchestration/pytest_selection/R05.txt` was not generated by the orchestrator
  Description: Closed by backfilling `.orchestration/pytest_selection/<unit>.txt` for all 24 Phase-1 review units (R01–R09, R10a/b/c, R11–R22), generated deterministically by running `git show --name-only` against the commit list transcribed from `MERGE_PLAN.md §B.1` (lines 292–315) and collecting paths under `tests/`. Each file lists the sorted-unique `tests/*` paths touched by its unit's commits; units that touched zero test files carry a single line `# no test changes in this unit`. 24/24 files present (`ls .orchestration/pytest_selection/ | wc -l` → 24). Six units carry test paths (R01: 9, R02: 3, R03: 1, R04: 6, R06: 1, R07: 1); the other 18 carry the no-test-changes marker — consistent with the unit topology (R05/R08 = scan-rerun + figure-prune only; R09–R20 = catalog/docs; R21/R22 = website). One typo fixed inside the generator (MERGE_PLAN.md `ebdo66c` → `ebd066c`; MERGE_PLAN.md itself untouched). Five other SHAs cited in MERGE_PLAN.md §B.1 could not be resolved (R08 `e0b4e2`, R15 `7500919`, R18 `e9f3cf3`, R19 `1cf8b57`, `cd3a8fe`); all five sit in catalog/figure-only units whose resolved siblings touch zero tests, so the selection content is unaffected — flagged in the C15 report as a follow-up MERGE_PLAN bookkeeping nit, not a C15 deliverable. Pure `.orchestration/` artifact backfill; no source, test, doc, or notebook touched; no test re-run required. Evidence in `.orchestration/cleanup_reports/C15.md`.

### Closed by C16

- [R20-I1] severity:LOW tag:docs  **CLOSED 2026-05-26** by C16.
  Title: SESSION_NOTES.md §1 fact-check tallies disagree with §1 table headline at v0.4 boundary
  Description: Closed by replacing the v0.3-tagged "Fact-check status across all 88 processes: 86 VERIFIED / 2 PARTIAL" paragraph at `flavor_catalog/SESSION_NOTES.md:38-39` with the v0.4 canonical block "Fact-check status across all 102 processes at v0.4 (post-Wave-9 canonical, per C07 reconciliation): **101 VERIFIED / 1 PARTIAL / 0 MISMATCH / 0 FAILED**". E009 retained as the only standing PARTIAL; K020 documented as v0.3 PARTIAL → cleared by `b5c2375` before the v0.4 tag → VERIFIED at v0.4, with cross-link to `flavor_catalog/audits/factcheck_status.md` §"Tag annotation provenance" (the C12 closure surface for R18-I1/R19-I3 that owns at-tag vs post-tag bookkeeping). New §1 numbers match C07 §"Aggregated fact-check totals" (`VERIFIED: 101 / PARTIAL: 1 / TOTAL: 102`), `master_compile_v04_report.md`, and the website `catalog_index.json verdict_counts: {VERIFIED: 101, PARTIAL: 1}` (per C12 R21-I1 closure). §1a (Wave-8) and §1b (Wave-9) historical summaries left intact (at-the-time records, not running tallies). Pure-docs cosmetic alignment; no code change; no test re-run required. Evidence in `.orchestration/cleanup_reports/C16.md`.

- [R20-I3] severity:INFO tag:docs  **CLOSED 2026-05-26** by C16.
  Title: SESSION_NOTES.md §8 cites three Opus sign-off rounds; catalog at v0.4 has five
  Description: Closed by appending `round_004_index.md (Wave-8 SECONDARY)` and `round_005_index.md (Wave-9 collider_rs PRIMARY)` to the per-round verdict list at `flavor_catalog/SESSION_NOTES.md:373-374` and rephrasing the trailing line as "Opus's per-round verdicts (5 rounds at v0.4)". `ls flavor_catalog/signoff/round_*.md` returns 5 files, matching the new text and HANDOFF_PROMPT.md `:81-86`. The C16 §1 edit (above) also adds a "Most recent Opus round-5 sign-off: 14/14 Wave-9 collider_rs PRIMARY APPROVE (CR001-CR014)" preamble with a back-pointer to round-4. Pure-docs cosmetic alignment; no code change; no test re-run required. Evidence in `.orchestration/cleanup_reports/C16.md`.

### Closed by C17

- [R09-I2] severity:INFO tag:docs  **CLOSED 2026-05-26** by C17.
  Title: `.gitkeep` placeholders left in family subdirs alongside `index.tex`
  Description: Removed 12 redundant `.gitkeep` placeholders in directories that now contain populated content (`flavor_catalog/processes/{beauty,charged_lepton,charm,edm_neutrino,kaon,top_higgs_ew}/`, `flavor_catalog/references/`, `flavor_catalog/signoff/by_process/`, and the four `flavor_catalog/worklogs/{checker,discovery,pka,writer}/`). Retained `flavor_catalog/signoff/round_index/.gitkeep` (sole tracked file in its directory). No functional impact.

### Closed by C18

- [R12-I1] severity:INFO tag:docs  **CLOSED 2026-05-26** by C18.
  Title: R12 dispatch prompt mis-identifies SINDRUM-II Au target as L008, undercounts W4 scope
  Description: Closed via the new "Retroactive corrections (C18, 2026-05-26)" annotation block in `.orchestration/MERGE_PLAN.md` (between §B.1 granularity check and the uniqueness invariant). The block records the canonical scope (W4 = 20 PKAs, not "12+") and pins L004 as the canonical home of the SINDRUM-II Au limit. Dispatch-prompt drift only; catalog remains consistent. Evidence in `.orchestration/cleanup_reports/C18.md`.

- [R12-I2] severity:LOW tag:code  **CLOSED 2026-05-26** by C18.
  Title: PKA commit `7d3da08` bundles two independent PKA drafts (EW001 + E004)
  Description: Closed via the C18 retroactive-corrections block in `MERGE_PLAN.md` which annotates commit `7d3da08` as bundling EW001 + E004 PKAs (workflow-seam irregularity preserved as-is, no git history change). Evidence in `.orchestration/cleanup_reports/C18.md`.

- [R12-I3] severity:INFO tag:code  **CLOSED 2026-05-26** by C18.
  Title: Commit `bcd8907` triple-bundles unrelated workstreams (L001 arbitration + kaon_edm CHECKER-DONE transitions + ca_w4_kaon_edm worklog)
  Description: Closed via the C18 retroactive-corrections block in `MERGE_PLAN.md` annotating `bcd8907` as bundling L001 Opus signoff + W4 kaon_edm CA cycle-1 verdict (E004/E006/E008/K017 `CHECKER-DONE`) + `ca_w4_kaon_edm.md` worklog. Evidence in `.orchestration/cleanup_reports/C18.md`.

- [R13-I1] severity:INFO tag:docs  **CLOSED 2026-05-26** by C18.
  Title: R13 dispatch prompt mis-identifies all three suggested spot-check process IDs (K010, B023, E007)
  Description: Closed via the C18 retroactive-corrections block in `MERGE_PLAN.md` recording the canonical referents: K010 = `K_S → π⁰ e⁺e⁻`, B023 = `B → K*(892) ν ν̄`, E007 = `²²⁵Ra / ¹²⁹Xe atomic EDMs`. R13 reviewer cross-verified the actual catalog content; this is a one-time dispatch-prompt error. Evidence in `.orchestration/cleanup_reports/C18.md`.

- [R14-I1] severity:INFO tag:docs  **CLOSED 2026-05-26** by C18.
  Title: MERGE_PLAN.md R14 row mis-attributes K012/K018 across the R14/R15 seam
  Description: Closed by editing `MERGE_PLAN.md:307` R14 description from "T006 v2, B001, B016, K012, K018; codex quota pause + resumption plan." to "T006 v2 + 8 other Wave-5 CA closures + B001 + B016 PKAs; codex quota pause + resumption plan. (K012, K018 belong to R15 per R14-I1.)". K012 PKA = `6a16bf4` and K018 PKA = `98b203a`, both correctly listed in the R15 row. Evidence in `.orchestration/cleanup_reports/C18.md`.

- [R15-I1] severity:LOW tag:docs  **CLOSED 2026-05-26** by C18.
  Title: MERGE_PLAN.md R15 row lists charged_lepton_top CA-v2 commit as `7500919`, actual SHA is `7500794`
  Description: Closed by editing `MERGE_PLAN.md:308` R15 commit list to replace `7500919` (no such object) with `7500794` (canonical `flavor-catalog(charged_lepton_top): CA batch ca_w6_charged_lepton_top_v2 — verify WA polish for L023 T020`). Verified via `git log --all --format='%h %s' | grep ^7500794` and `git log --grep='ca_w6_charged_lepton_top_v2' --oneline`. Evidence in `.orchestration/cleanup_reports/C18.md`.

- [R16-I3] severity:INFO tag:docs  **CLOSED 2026-05-26** by C18.
  Title: B012 promotion from DEFERRED-SCOPE to ACTIVE could be misread as a quiet DA-rule amendment
  Description: Closed via the C18 retroactive-corrections block in `MERGE_PLAN.md` adopting the preferred framing: "DA-4 converged at 75; Wave-7 amended active coverage to 80 (4 unallocated + 1 previously-deferred B012 promoted on external-review grounds)." Documentation seam only. Evidence in `.orchestration/cleanup_reports/C18.md`.

- [R19-I1] severity:LOW tag:provenance  **CLOSED 2026-05-26** by C18.
  Title: MERGE_PLAN.md:312 R19 commit list (and the user prompt) list `1cf8b57` as the CR005 PKA SHA; the actual SHA is `1cd8b57`
  Description: Closed by editing `MERGE_PLAN.md:312` R19 commit list to replace `1cf8b57` (no such object) with `1cd8b57` (canonical `flavor-catalog(wave9): PKA-CR005 initial draft (collider_rs PRIMARY)`). Verified via `git log --all --format='%h %s' | grep ^1cd8b57`. Evidence in `.orchestration/cleanup_reports/C18.md`.

- [R19-I2] severity:LOW tag:provenance  **CLOSED 2026-05-26** by C18.
  Title: Wave-9 CA-v2 ew_tail cycle-2 verdict commit `82daa9b` is in the Wave-9 chain but is NOT listed in MERGE_PLAN.md:312 R19 commit list (or the user prompt)
  Description: Closed by inserting `82daa9b` between `950ca36` and `0c5dacc` in `MERGE_PLAN.md:312` R19 commit list and bumping the count "(33 commits ...)" → "(34 commits ...)". Granularity-check sentence at MERGE_PLAN.md:317 also updated to read "R19 is now correctly 34 commits (R03 commit `82a96f0` is excluded by design; missing `82daa9b` restored per C18)". The "82a96f0 belongs to R03 ONLY" caveat preserved. Range verification: `git log --oneline 7ed9117^..864cd6d | wc -l = 34`. Evidence in `.orchestration/cleanup_reports/C18.md`.

- [R19-I5] severity:INFO tag:docs  **CLOSED 2026-05-26** by C18.
  Title: MERGE_PLAN.md:413 R19 grep-target list mentions `claim_level` field, but the collider_rs schema uses `limit_type` (under `pdg_or_equivalent.values[]`)
  Description: Closed by editing the R19 grep-target entry at MERGE_PLAN.md (line ~419 after C18 inserts) to replace `claim_level` with `limit_type` and add a parenthetical noting the planner-era field-name guess. `grep -l "claim_level" flavor_catalog/processes/collider_rs/*.yaml` returns nothing; `grep -c "limit_type:" flavor_catalog/processes/collider_rs/*.yaml` shows the field is present in all 14 CR0XX yamls. Evidence in `.orchestration/cleanup_reports/C18.md`.

### Closed by C19

- [R08-I3] severity:INFO tag:infra  **CLOSED 2026-05-26** by C19 (accept-as-paper-archive-companion).
  Title: Historical-snapshot figure dirs (`quark_pre_audit_constants/`, `quark_baseline_800k/`) deliberately left in tree
  Description: Closed via accept-as-paper-archive-companion disposition recorded in `docs/audits/figure_prune_inventory.md` (new "Status as of 2026-05-26 post-cleanup (C19)" block). `quark_baseline_800k/` is already absent from the tree; `quark_pre_audit_constants/` is retained (58 files / ~5 MB) on `\includegraphics`-irrelevance and provenance grounds. No tracked code, no test code, no current TeX source references these dirs. Methodology PDF rebuilt (20 pp, 627505 bytes) with the C03→C19 substantive forwarding cross-reference sentence added (one line at `docs/quark_scan_methodology_note.tex`).

### Closed by C20 (auto-resolved / accepted-risk)

Sweep of all remaining `## Open` INFO/LOW issues that are non-actionable
(documentation cosmetic seams, workflow-seam observations, or items whose
original `Recommended fix` explicitly says "None", "Optional", "Accept as
adequate", or "Accept-as-is otherwise"). Two entries are AUTO-RESOLVED (the
upstream artifact those issues tracked has since been corrected). The rest
are ACCEPTED-RISK. Full C20 disposition table at
`.orchestration/cleanup_reports/C20.md`; per-issue original text preserved
below verbatim for audit trail.

- [R01-I1] severity:LOW tag:code  **CLOSED 2026-05-26** by C20 (ACCEPTED-RISK).
  Title: Dual source-of-truth for heavy-quark mass constants
  Description: Per-mille drift between `qcd/constants.py` (`M_CHARM=1.27, M_BOTTOM=4.18, M_TOP_MS=163.5`, used as threshold positions in `qcd/mass_running.py:106-108`) and `quarkConstraints/pdg_quark_masses.py:96-117` (PDG-2024 `m_c=1.273, m_b=4.183, m_t=162.5`) is acknowledged in `docs/quark_scan_methodology_note.tex:53-54`. Drift is per-mille so not load-bearing; dual-SSOT pattern accepted-as-is for the quark-scan paper. Future refactor option (re-export pattern) recorded but not gated on paper finalization.
  Evidence: `qcd/constants.py:17-19`; `quarkConstraints/pdg_quark_masses.py:96-117`; `docs/quark_scan_methodology_note.tex:53-54`; review `.orchestration/reviews/R01.md` Check 3.

- [R01-I2] severity:INFO tag:code  **CLOSED 2026-05-26** by C20 (ACCEPTED-RISK).
  Title: `MASS_TOLERANCE_FLOOR` is a placeholder pending Phase-0 calibration
  Description: `quarkConstraints/scan.py:40` sets `MASS_TOLERANCE_FLOOR = 0.003` with a TODO at lines 37-39 pointing at `scripts/calibrate_phase0.py`. The calibration script exists in-tree but has not been re-run since the kaon-constants audit (C01); Phase-0 floor calibration is deferred to the paper-finalization phase (Tier-6 deferred work). ACCEPTED-RISK rationale: the current 0.003 floor is conservative; any future re-calibration would tighten not loosen the floor, and the headline numerics in the methodology note are tested against the per-flavor log-residual structure independently.
  Evidence: `quarkConstraints/scan.py:30-40`; `scripts/calibrate_phase0.py` (in-tree); review `.orchestration/reviews/R01.md` Check 3.

- [R02-I1] severity:INFO tag:docs  **CLOSED 2026-05-26** by C20 (ACCEPTED-RISK).
  Title: MERGE_PLAN row for R02 names `pytest.ini` but the strict-xfail block landed in `pyproject.toml`
  Description: `.orchestration/MERGE_PLAN.md:293` lists `pytest.ini` as an R02 file scope, but the `xfail_strict = true` setting was added to `pyproject.toml` under `[tool.pytest.ini_options]`. Functionally equivalent for pytest; MERGE_PLAN inventory nit only. ACCEPTED-RISK.
  Evidence: `pyproject.toml:70-71`; `.orchestration/MERGE_PLAN.md:293`; commit `559b851`.

- [R03-I4] severity:INFO tag:docs  **CLOSED 2026-05-26** by C20 (ACCEPTED-RISK).
  Title: No top-level REFERENCES.md to cross-reference audit citations against
  Description: Provenance is encoded as inline arXiv citations in `docs/audits/bag_param_inventory.md` and the methodology-note appendix. All five primary eprints (2411.04268, 1911.06822, 1907.01025, 1505.06639, 1002.3612) and PDG 2024 review (Phys. Rev. D 110, 030001) are real and externally resolvable. ACCEPTED-RISK: inline-citation convention accepted as canonical for the quark-scan paper.
  Evidence: `docs/audits/bag_param_inventory.md:16-53`; `docs/quark_scan_methodology_note.tex` appendix `app:hadronic-input-provenance`.

- [R05-I1] severity:INFO tag:docs  **CLOSED 2026-05-26** by C20 (ACCEPTED-RISK).
  Title: Rerun commit subject says "8 RS-anarchy scans" but 9 dated directories exist
  Description: `db02223` subject says "regenerate 8 RS-anarchy scans"; actual rerun set is 9 SLURM jobs (RUNA, Run-3 baseline/qtop/moreUV/moreIR, Run-B narrow/wide/gaussian, Run-C). Two Run-3 variants (moreUV/moreIR) had zero PDG-passing draws, justifying the "8" count. Purely cosmetic. ACCEPTED-RISK.
  Evidence: commit `db02223` subject; `docs/phase_logs/invalidation_gate_rerun.md:14-22`; `.orchestration/MERGE_PLAN.md:296`.

- [R05-I3] severity:INFO tag:docs  **CLOSED 2026-05-26** by C20 (ACCEPTED-RISK).
  Title: Possible reader confusion between benchmark 22.49x cumulative factor and the empirical ~4.47x headline scale shift
  Description: Methodology-note 22.49 = 6.2388 * 3.6052 cumulative invalidation-factor and the empirical 4.467x headline rescaling are already distinguished in `docs/phase_logs/invalidation_gate_rerun.md:46-51`. ACCEPTED-RISK: distinction is documented at the canonical phase-log source; methodology-note prose reads adequately for the paper audience.
  Evidence: `scan_outputs/followup_crossings_summary.json`; `docs/phase_logs/invalidation_gate_rerun.md:46-51`.

- [R06-I1] severity:LOW tag:code  **CLOSED 2026-05-26** by C20 (ACCEPTED-RISK).
  Title: CFW convention constants live as private module constants of the audit driver
  Description: `CFW_RS_GS3_TEV=10.5`, `CFW_RS_GS6_TEV=21.0`, `CFW_PGB_GS3_TEV=17.0`, `CFW_PGB_GS6_TEV=33.0`, `EPSILON_K_BUDGET_DEFAULT`, `G_S_PERT` live in `scripts/rs_anarchy_cfw_comparison.py:41-78`; regression test reaches into the script via `from scripts import rs_anarchy_cfw_comparison as cfw_script`. ACCEPTED-RISK per the issue's own recommended fix ("Accept-as-is otherwise"): acceptable for a one-off audit driver. Refactor option reserved for the day a second-paper comparison driver is added.
  Evidence: `scripts/rs_anarchy_cfw_comparison.py:41-78`; `tests/test_cfw_comparison.py:15,127-129`.

- [R06-I2] severity:INFO tag:docs  **CLOSED 2026-05-26** by C20 (ACCEPTED-RISK).
  Title: CFW comparison PNG no longer ships at the publication path
  Description: PNG relocated to `results/figures/quark/exploratory/`; PDF stays at the publication path. Methodology note `\includegraphics` references the PDF only, so LaTeX build is unaffected. ACCEPTED-RISK: cosmetic discrepancy between phase-log narrative and final-tree state; no consumer is broken.
  Evidence: `results/figures/quark/rs_anarchy_cfw_comparison.pdf`; relocated PNG at `exploratory/` subdir; `docs/quark_scan_methodology_note.tex:882`.

- [R06-I3] severity:INFO tag:docs  **CLOSED 2026-05-26** by C20 (ACCEPTED-RISK).
  Title: Headline number rounding mismatch between `cfw_vs_ours.md` and methodology note
  Description: `docs/audits/cfw_vs_ours.md:12` quotes precise percentile spread `47.26^{+69.37}_{-24.98} TeV`; `docs/quark_scan_methodology_note.tex:898-903` rounds to `~47 TeV`. Numbers are consistent (precise-source vs rounded-summary). ACCEPTED-RISK per issue's own recommended fix ("Accept as adequate").
  Evidence: `docs/audits/cfw_vs_ours.md:12`; `docs/quark_scan_methodology_note.tex:898-903`.

- [R07-I1] severity:LOW tag:physics  **CLOSED 2026-05-26** by C20 (ACCEPTED-RISK).
  Title: Wilson-score z=1.92 default is non-standard; matches CP two-sided 95% at large n
  Description: `quarkConstraints/finite_stats.py:8` default `z=1.92` is documented in the helper docstring (lines 9-13) and in `docs/audits/zero_pass_inventory.md:13-14` as the "audit convention". Numerical values are correct; only the convention choice (z=1.92 vs standard Wilson z=1.96) is presentational. ACCEPTED-RISK: the docstring + audit-doc disclosure is judged sufficient for the paper audience; an explanatory methodology-note footnote remains an Optional follow-up (not gating).
  Evidence: `quarkConstraints/finite_stats.py:8-13`; `docs/quark_scan_methodology_note.tex:706-709, 748-750, 766-771, 866`; `docs/audits/zero_pass_inventory.md:13-14`.

- [R07-I4] severity:INFO tag:docs  **CLOSED 2026-05-26** by C20 (AUTO-RESOLVED).
  Title: `CLAUDE.md:12-15` paper-branch pointer to `paper/quark-scan-2026q2` will be invalidated by Phase 10
  Description: AUTO-RESOLVED. The `CLAUDE.md:12-15` paper-branch pointer was updated during Phase-1 (post-consolidation) editing; the current `CLAUDE.md:12-16` block reads "All work lives on `main` (post-2026-05-25 consolidation; tagged `v2026q2-catalog-complete`). The canonical quark-sector methodology note is `docs/quark_scan_methodology_note.tex`/`.pdf`. This paper is quark-sector only; lepton-sector directories are follow-up scope. The Cloudflare-deployed catalog website builds from `main` with root `flavor_catalog/website/`." The planned MERGE_PLAN Phase-10 update landed; the original R07-I4 tracking entry is now self-resolved.
  Evidence: `CLAUDE.md:12-16` (current state, post-Phase-1); commit `1a107c8` (original pointer); review `.orchestration/reviews/R07.md` Check 3.

- [R08-I1] severity:INFO tag:docs  **CLOSED 2026-05-26** by C20 (ACCEPTED-RISK).
  Title: `bf6186c` "final PDF rebuild after figure prune" committed a byte-identical PDF (no-op)
  Description: `bf6186c` records `docs/quark_scan_methodology_note.pdf` as `596242 -> 596242 bytes` (no change). Correct behavior — the figure prune (`e7d824d`) only relocated unreferenced files. ACCEPTED-RISK: harmless commit-log noise; meaningful post-prune rebuild ships at `e0b24e8` (620832 bytes).
  Evidence: `git show --stat bf6186c`; `git show --stat e0b24e8`; `artifacts/checksums.sha256:21`.

- [R08-I2] severity:INFO tag:numerics  **CLOSED 2026-05-26** by C20 (ACCEPTED-RISK).
  Title: rc1.1 WARNING-3 fix shifts `f_{u,2}` from `0.117` to `0.16` (37% upward)
  Description: `e0b24e8` shifts `f_{u,2}` in `docs/quark_scan_methodology_note.tex:522` from 0.117 to 0.16, simultaneously relaxing surrounding `(m_u, m_c, m_t)` estimates from explicit values to `\mathcal{O}(...)` order-of-magnitude form. The two changes are consistent. ACCEPTED-RISK per issue's own recommended fix ("accept since the surrounding numerical content was relaxed to `\mathcal{O}(...)` to match"). Spot-check at canonical c-pattern remains an Optional follow-up.
  Evidence: `docs/quark_scan_methodology_note.tex:522, 527-532`; commit `e0b24e8`.

- [R09-I1] severity:INFO tag:docs  **CLOSED 2026-05-26** by C20 (ACCEPTED-RISK).
  Title: Scaffold collapsed plan-v1 separate `edm/` + `neutrino_universality/` families into single `edm_neutrino/` dir
  Description: Plan-v1 listed 7 PRIMARY family directories; scaffold (`83c0178`) created 6, merging `edm/` + `neutrino_universality/` into `edm_neutrino/`. Collapse is recorded in `docs/phase_logs/flavor_catalog_scaffold_impl.md`; consumed without contradiction by Waves 1-9. ACCEPTED-RISK: stale 7-family listing in plan-v1 is a superseded planning artifact; current catalog tree is the canonical source.
  Evidence: `docs/phase_logs/flavor_catalog_plan_v1.md:30-55`; `flavor_catalog/catalog_master.tex:48-49`; `flavor_catalog/README.md:38`; `docs/phase_logs/flavor_catalog_scaffold_impl.md`.

- [R10a-I1] severity:INFO tag:docs  **CLOSED 2026-05-26** by C20 (ACCEPTED-RISK).
  Title: K001 sidecar `checker_agent_id: "CA"` points to a stale wave label
  Description: `flavor_catalog/processes/kaon/K001.yaml:9` declares `checker_agent_id: "CA"`; the actual CA cycle batch_id is in `status_history` at line 34. ACCEPTED-RISK: the `status_history.batch_id` field is the actual provenance source; top-level `checker_agent_id` is just unscoped; future polish-pass remains Optional.
  Evidence: `flavor_catalog/processes/kaon/K001.yaml:9, 34`; `flavor_catalog/worklogs/checker/ca_w23_kaon_charm_edm.md`.

- [R11-I1] severity:LOW tag:code  **CLOSED 2026-05-26** by C20 (ACCEPTED-RISK).
  Title: B017 introduced as new yaml+tex pair inside a WA-polish commit
  Description: B017 first appears in WA-polish commit `5031e06`, not in a dedicated PKA-draft commit (all other Wave-2/3 new processes land as dedicated PKA-draft first). File pair is complete; yaml parses against `flavor_catalog.process.v1`; DA-1 inventory counts B017 in the beauty=10 tally. ACCEPTED-RISK: workflow-seam irregularity, not a missing-process defect; no code or physics impact.
  Evidence: `git log --diff-filter=A --name-only -- flavor_catalog/processes/beauty/B017.yaml` = `5031e06`; `flavor_catalog/worklogs/discovery/round_001_full_scope.md:11`.

- [R13-I2] severity:LOW tag:docs  **CLOSED 2026-05-26** by C20 (ACCEPTED-RISK).
  Title: `signoff/round_002_index.md` description column has wrong process names for K009 and K010
  Description: K009 row reads "`K^+ → π^+ ℓ⁺ℓ⁻` form-factor" (matches K017, not K009) and K010 row reads "`K^+ → π^+ π⁰ γ / π⁰ e+e- γ` radiative" (does not match K010). APPROVE verdicts trace through CA worklog `ca_w5b_kaon_v2.md` which correctly identifies the actual process IDs. ACCEPTED-RISK: human-readable description column in signoff-index is wrong but verdict integrity is preserved via the CA worklog cross-trace; corrective edit remains Optional per issue's own recommended fix.
  Evidence: `flavor_catalog/signoff/round_002_index.md` vs `flavor_catalog/processes/kaon/{K009,K010}.yaml:4-5`; `flavor_catalog/worklogs/checker/ca_w5b_kaon_v2.md`.

- [R13-I3] severity:INFO tag:docs  **CLOSED 2026-05-26** by C20 (ACCEPTED-RISK).
  Title: DA-3 worklog does not separate Wave-4 vs Wave-5a pipeline lineage for charm family
  Description: `round_003_final_sweep.md:8` reports the 8-drafted charm-family count correctly; family-by-family Wave-4 vs Wave-5a pipeline distinction is not surfaced in the wording. ACCEPTED-RISK: cosmetic only; family inventory total is consistent with the catalog state; Wave-5a (`C005`) pipeline lineage is correctly captured in `signoff/round_001_index.md:48-53` via the `ca_w5a_kaon_charm.md` cross-link.
  Evidence: `flavor_catalog/worklogs/discovery/round_003_final_sweep.md:8-13`; `flavor_catalog/signoff/round_001_index.md:48-53`.

- [R14-I2] severity:INFO tag:docs  **CLOSED 2026-05-26** by C20 (AUTO-RESOLVED).
  Title: B001 yaml carries an unresolved `open_issue` about x_d / chi_d prose treatment
  Description: AUTO-RESOLVED per the issue's own recommended fix ("Close automatically once R15 is reviewed; no R14-side change required"). The B001 open_issue ("WA/CA should verify whether to retain only Delta m_d as the headline B001 observable or also list x_d and chi_d as companion values in prose") was resolved downstream in R15 via the Wave-6 WA → CA → Opus arbitration chain at `flavor_catalog/signoff/by_process/B001_B003.md`: the arbitration ruled `Delta m_d` is the headline prose observable; `x_d` and `chi_d` stay in `pdg_or_equivalent.companion_*` blocks. R15 was reviewed and merged.
  Evidence: `flavor_catalog/processes/beauty/B001.yaml:226` (original open_issue); `flavor_catalog/signoff/by_process/B001_B003.md` (resolution); MERGE_PLAN.md:308.

- [R15-I2] severity:INFO tag:docs  **CLOSED 2026-05-26** by C20 (ACCEPTED-RISK).
  Title: E009 carries `fact_check_verdict: PARTIAL` due to JS-only INSPIRE Weinberg-1989 URL render
  Description: E009 PARTIAL verdict is an artifact of the INSPIRE web app rendering metadata via JavaScript; local snapshot at `flavor_catalog/references/E009/weinberg1989_inspire_metadata.txt` (sha256 `d50b6648...`) is the canonical source of record; all numerics fully VERIFIED. ACCEPTED-RISK: cosmetic verdict-label artifact; updating to `VERIFIED-WITH-NOTE` remains Optional per issue's own recommended fix.
  Evidence: `flavor_catalog/processes/edm_neutrino/E009.yaml:81-83`; `flavor_catalog/references/E009/weinberg1989_inspire_metadata.txt`.

- [R15-I3] severity:INFO tag:docs  **CLOSED 2026-05-26** by C20 (ACCEPTED-RISK).
  Title: L023 (neutrino trident) filed under `family: charged_lepton` though it is technically a neutrino process
  Description: L023 classification under `family: charged_lepton` is defensible (muon current and the contact interaction are the load-bearing degrees of freedom); `L023.tex:14-18` explicitly disclaims this. ACCEPTED-RISK per issue's own recommended fix ("None — current classification is documented and consistent"). Reconsider only if >2-3 neutrino-process entries are ever added.
  Evidence: `flavor_catalog/processes/charged_lepton/L023.yaml:3`; `flavor_catalog/processes/charged_lepton/L023.tex:14-18`.

- [R16-I1] severity:INFO tag:docs  **CLOSED 2026-05-26** by C20 (ACCEPTED-RISK).
  Title: T012 single-entry consolidates three Z-pole charm observables and absorbs T013 (charm asymmetry)
  Description: T012 single-entry for `R_c^0` + `A_FB^{0,c}` + `A_c` is defensible (jointly determined from the same LEP-SLC Z-pole fit). DA-4 addendum `round_004_addendum_deferred_scope.md:26` explicitly folds T013 into T012. ACCEPTED-RISK per issue's own recommended fix ("None at present scale"). Future PI scope decision to commission an independent charm-asymmetry likelihood would re-open via a fresh DA wave.
  Evidence: `flavor_catalog/processes/top_higgs_ew/T012.yaml:3-11`; `flavor_catalog/worklogs/discovery/round_004_addendum_deferred_scope.md:26`.

- [R16-I2] severity:INFO tag:docs  **CLOSED 2026-05-26** by C20 (ACCEPTED-RISK).
  Title: R16 task-prompt anchor "BR(t->c gamma) < ~1.8e-4 ATLAS" is stale; catalog correctly carries 1.51e-5 CMS Run-2
  Description: Catalog headline at `T003.yaml:119` is the current CMS-TOP-21-013 Run-2 result `< 1.51e-5 @ 95% CL`. The R16 dispatch prompt's expected anchor was a stale Run-1 / pre-2020 limit. ACCEPTED-RISK per issue's own recommended fix ("None — catalog is current"). Future automated-prompt generation should source anchors from `pdg_or_equivalent.value_summary` rather than external recall.
  Evidence: `flavor_catalog/processes/top_higgs_ew/T003.yaml:112-185`.

- [R16-I4] severity:INFO tag:docs  **CLOSED 2026-05-26** by C20 (ACCEPTED-RISK).
  Title: B001/B003 arbitration commit `7e1b80b` chronologically predates the Wave-7 PKA cycle
  Description: Arbitration commit `7e1b80b` timestamps `2026-05-16 16:35:01`; earliest Wave-7 PKA commit is `2026-05-16T18:26:49`. MERGE_PLAN topical grouping under R16 is correct; strict-temporal reading would attach to R15. ACCEPTED-RISK per issue's own recommended fix ("None — the topical grouping is correct"). R15 already cross-references the arbitration via `signoff/by_process/B001_B003.md`; no information is lost.
  Evidence: `.orchestration/MERGE_PLAN.md:309`; `git log -1 --format="%ci" 7e1b80b`.

- [R18-I2] severity:LOW tag:provenance  **CLOSED 2026-05-26** by C20 (ACCEPTED-RISK).
  Title: Commit `37beabb` message reads "PKA-K020 initial draft" but the diff contains the B013 family files
  Description: `37beabb` is mis-labeled (commit subject says K020 but files are B013). K020 content actually landed in `bab5bd0` and a second B013 commit `ebd066c` added more references. File tree contains all expected K020 and B013 files. ACCEPTED-RISK per issue's own recommended fix ("None at present (history is immutable)"). MERGE_PLAN R18 commit list correctly attributes both commits.
  Evidence: `git show --stat 37beabb`, `git show --stat bab5bd0`; `flavor_catalog/processes/secondary/{beauty/B013,kaon/K020}.{yaml,tex}` all present.

- [R18-I4] severity:INFO tag:docs  **CLOSED 2026-05-26** by C20 (ACCEPTED-RISK).
  Title: SECONDARY tier rationale is documented in TWO places (per-YAML `priority_rationale` AND `PRIORITY_TIERS.md` §4 table)
  Description: Per-process `priority_rationale` field + wave-level `PRIORITY_TIERS.md:90-99` table both encode the SECONDARY tier reasoning; the two sources agree pairwise for all 8 SECONDARY entries. ACCEPTED-RISK per issue's own recommended fix ("None at present"): redundancy is desirable (per-process SoT + wave-level overview); drift mitigated by the wave-runbook §6 reproducibility-note pattern.
  Evidence: `flavor_catalog/processes/secondary/kaon/K020.yaml:4-6`, `K019.yaml:4-5`; `flavor_catalog/PRIORITY_TIERS.md:90-99`.

- [R21-I2] severity:INFO tag:infra  **CLOSED 2026-05-26** by C20 (ACCEPTED-RISK).
  Title: `cloudflare-pages.config.md` is a markdown doc, not a `wrangler.toml`
  Description: Cloudflare Pages git-driven deployments do not require `wrangler.toml` (unlike Cloudflare Workers); build command + output dir are entered in the dashboard. The website pairs the doc with `.node-version`, `public/_redirects`, `public/_headers` (the machine-readable Pages config files). ACCEPTED-RISK per issue's own recommended fix ("Optional ... Out of scope for R21"). Pattern is correct for Pages.
  Evidence: `flavor_catalog/website/cloudflare-pages.config.md:14-25`; `.node-version`, `public/_redirects:11-12`, `public/_headers:1-20`.

- [R21-I3] severity:LOW tag:docs  **CLOSED 2026-05-26** by C20 (ACCEPTED-RISK).
  Title: `WEBSITE_RUNBOOK.md` STOP-and-ping trigger was tripped (charm 23.7% UNRESOLVED) but handled as soft-trigger
  Description: `WEBSITE_RUNBOOK.md:35` STOP-and-ping policy vs `:53-55` soft-trigger handling under PI's "do not stop" directive is a wording seam. UNRESOLVED anchors are explicitly labelled in the citation modal so no false anchors shipped. ACCEPTED-RISK per issue's own recommended fix ("Not blocking"); the PI directive is captured in the phase-2 ledger entry and commit `6ffbb33` body.
  Evidence: `flavor_catalog/website/WEBSITE_RUNBOOK.md:35, 53-55`; commit `6ffbb33`.

- [R22-I2] severity:LOW tag:docs  **CLOSED 2026-05-26** by C20 (ACCEPTED-RISK).
  Title: methodology page remains at `/methodology/` after `5f31f2d` removes the nav link
  Description: Methodology page (`src/pages/methodology.astro`, 233 lines) is unlinked-from-nav but reachable via direct URL + family-page internal links. Intentional per commit `5f31f2d` body ("Remove Methodology from nav and home CTA"). ACCEPTED-RISK per issue's own recommended fix ("Not blocking").
  Evidence: `flavor_catalog/website/src/pages/methodology.astro`; `flavor_catalog/website/src/layouts/BaseLayout.astro` (nav after `5f31f2d`); commit `5f31f2d` body.

### Auto-resolved / no-op

- [R17-I3] severity:INFO tag:docs  **CLOSED 2026-05-26** by C16 (PRE-EXISTING-TAG-VERIFIED, no action); confirmed by C20 as the canonical location.
  Title: Catalog "v0.2" snapshot — produce or annotate
  Description: Per CR-1 in the cleanup-plan r2 revision, the original C16 v0.2 tag-creation task is DROPPED: the `flavor-catalog-v0.2` annotated tag already exists at commit `835cf48` (verified by `git show-ref --tags flavor-catalog-v0.2`). C16 took no tag action. C20 confirms this `### Auto-resolved / no-op` subsection as the canonical home for R17-I3 (the Open-section duplicate entry was removed by C20); the C16 closure marker is now final, not provisional. Evidence: `.orchestration/cleanup_reports/C16.md` §"Note on R17-I3"; `.orchestration/cleanup_reports/C20.md`.

## Infra follow-ups
- INFRA-1 severity:LOW tag:infra — Reconfigure Cloudflare to deploy from `main/flavor_catalog/website/` so the second branch can eventually be retired. **CLOSED 2026-05-25**: website branch merged into main (commit `cb58a36`); Cloudflare Pages production branch set to `main` with root `flavor_catalog/website/`; verified green deploy on commit `a809cc3`; `flavor-catalog-website/2026q2` deleted local + origin.
