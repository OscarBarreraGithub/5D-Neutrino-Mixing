# MERGE + CODE-REVIEW PLAN — Round 2 (Planner)

Repo: `/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing` (origin: `OscarBarreraGithub/5D-Neutrino-Mixing`)
Author: Claude Opus 4.7 (planner role, round 2 revision)
Status: REVISED after reviewer round 1 critique (see `PLAN_REVIEW_R1.md`) — now execution-ready

Round-2 changes are summarized in `PLAN_CHANGES_R2.md`. The substantive deltas vs. round 1 are:
- New **Phase 0 (pre-merge user gate)** for untracked-file disposition (resolves CR-3).
- New **Phase 0.5 (worktree prune + state snapshot)** before safety tagging (resolves RC-5, M-1).
- Phase 0 safety-tag push is now **verified** before any mutation proceeds (resolves RC-4).
- R19 commit list no longer includes `82a96f0` (resolves CR-1).
- All references to the non-existent SHA `82b35c0` replaced with the verified canonical `82daa9b` (resolves CR-2).
- R10 split into **R10a/R10b/R10c** by process family (resolves RC-1, OQ-3).
- Catalog-physics spot-check now requires one PRIMARY process per touched family per wave (resolves OQ-4).
- Added per-unit grep checklists for high-risk code units (RC-6).
- Added literature-manifest contract (RC-2), per-unit pytest selection (RC-3), test-timeout policy (M-5), GitHub-issue mirroring step (M-2), orchestrator commit policy (M-3), and rollback runbook reference (M-4).
- **Total unit count: 24** (R01-R09 + R10a/R10b/R10c + R11-R22).

---

## 0. Scope reminder

The user wants:
1. Branch consolidation into one (or two — website may stay) trunk.
2. Code review broken into atomic "submodules of hand-written work" with three checks per unit: (i) physics intuition vs literature, (ii) numerical correctness, (iii) repo-wide code self-consistency.
3. Reviews keyed to git commits.
4. Issue tracking.
5. GitHub mirrors local at the end.

This document is the orchestrator's blueprint. It does NOT execute anything; it has been iterated through one planner -> reviewer -> planner cycle and is now ready for execution.

---

## A. Branch consolidation strategy

### A.1 Subsumption matrix (verified via `git log <trunk>..<branch>` and `git diff --stat`)

Trunk candidate: `flavor-catalog/2026q2` (local + origin in sync, 369 commits ahead of `main`).
Website branch: `origin/flavor-catalog-website/2026q2` (382 commits ahead of `main`; 13 commits ahead of trunk, all touching ONLY `flavor_catalog/website/`).

| # | Branch | Ahead-of-trunk | Ahead-of-website | Ahead-of-main | Status | Recommended action |
|---|---|---|---|---|---|---|
| 1 | `audit/bag-inputs` | 0 | 0 | 16 | Fully subsumed by trunk | Delete local + origin |
| 2 | `audit/cfw-comparison` | 0 | 0 | 32 | Subsumed | Delete local + origin |
| 3 | `audit/wilson-rg` | 0 | 0 | 20 | Subsumed | Delete local + origin |
| 4 | `backup/ca-v2-ew-tail-399f116` | 1 | 1 | 358 | Safety snapshot. Tip commit `399f116` is an obsolete pre-CHECKER-DONE state of CR009/CR011; the canonical version lives at `82daa9b` on trunk (verified: `git log -1 --format='%h %s' 82daa9b` -> `82daa9b flavor-catalog(wave9): CA-v2 ew_tail cycle-2 verdict`, same subject as backup tip). | Convert to **annotated git tag** `safety/ca-v2-ew-tail-399f116`, delete branch local + origin |
| 5 | `backup/ca-v2-ew-tail-combined-776eff9` | 1 | 1 | 358 | Same content as #4 (sibling snapshot); canonical replacement also `82daa9b`. | Convert to tag `safety/ca-v2-ew-tail-combined-776eff9`, delete branch local + origin |
| 6 | `ca_wave1_top_higgs_ew_v2_push` | 0 | 0 | 120 | Subsumed | Delete local + origin |
| 7 | `docs/canonical-manifest` | 0 | 0 | 52 | Subsumed (no remote counterpart) | Delete local |
| 8 | `docs/paper-readiness` | 0 | 0 | 60 | Subsumed | Delete local + origin |
| 9 | `docs/scope-wording` | 0 | 0 | 48 | Subsumed | Delete local + origin |
| 10 | `fix/pdg-benchmarks` | 0 | 0 | 12 | Subsumed | Delete local + origin |
| 11 | `invalidation/gate-rerun` | 0 | 0 | 26 | Subsumed | Delete local + origin |
| 12 | `paper/quark-scan-2026q2` | 0 | 0 | 82 | Subsumed. (CLAUDE.md references this as active paper branch — outdated.) | Delete local + origin, then update `CLAUDE.md` in a follow-up commit on trunk |
| 13 | `pka-T017` | 0 | 0 | 201 | Subsumed | Delete local (no remote) |
| 14 | `pka-l008-push` | 0 | 0 | 143 | Subsumed | Delete local |
| 15 | `scan/zero-pass-statistics` | 0 | 0 | 43 | Subsumed | Delete local + origin |
| 16 | `wa-w5a-top-higgs-ew-v2-152230` | 0 | 0 | 220 | Subsumed | Delete local |
| 17 | `wa_w5b_charm_edm_work` | 0 | 0 | 211 | Subsumed (currently checked-out as a `prunable` worktree under `/tmp/wa_w5b_charm_edm_worktree`) | **Run `git worktree prune` in Phase 0.5 BEFORE deletion**, then delete local |
| 18 | `main` | n/a | n/a | n/a | Old historical pointer (last commit 2026-04-20) | **Fast-forward to trunk** at end of process |
| 19 | `flavor-catalog-website/2026q2` (origin only) | n/a | n/a | n/a | Cloudflare deploys from this branch | Keep — merge trunk into it after consolidation |

**Verification snippet for the canonical-vs-backup claim** (orchestrator should re-run):
```
git log -1 --format='%h %s' 82daa9b
git log -1 --format='%h %s' 399f116
git log -1 --format='%h %s' 776eff9
git diff 399f116 82daa9b -- flavor_catalog/processes/collider_rs/ | head
# expect: same "CA-v2 ew_tail cycle-2 verdict" subject and a small diff that the canonical commit on trunk supersedes
```

**Key finding**: 16 of 17 non-trunk feature branches contribute ZERO unique commits relative to `flavor-catalog/2026q2`. The two backup branches each carry one obsolete intermediate state, which is informational, not load-bearing. The website branch is the only branch that is genuinely ahead of trunk, and it touches a strictly disjoint subtree (`flavor_catalog/website/`).

### A.2 Target trunk recommendation

**Recommendation: fast-forward `main` to the final reviewed tip and use `main` as the single trunk.**

Rationale:
- The user said "one (or two) branches"; the cleanest realization is `main` (trunk) + `flavor-catalog-website/2026q2` (deploy).
- `main` is the GitHub default; PRs, issue search, and external links all point there.
- Naming `flavor-catalog/2026q2` as the trunk perpetuates a "quarterly branch" naming convention that the user implicitly rejects ("merge everything into one").
- The fast-forward is lossless: every commit on the historical `main` is already in `flavor-catalog/2026q2`.

End state: `main` contains all 369 commits + any new merge commit + the merged website if we choose to fold website in (we will NOT fold it; see A.3). `flavor-catalog/2026q2` itself becomes a stale label and should be deleted after the FF, OR kept as a release-snapshot tag (`v2026q2`).

Tag policy after FF:
- `safety/pre-merge-2026-05-25` — annotated tag on current trunk tip before any operation. **THIS IS THE SAFETY NET**. Must be pushed to origin AND verified before any deletion or FF (see Phase 0).
- `safety/ca-v2-ew-tail-399f116`, `safety/ca-v2-ew-tail-combined-776eff9` — tag the two backup branches before deletion.
- `v2026q2-catalog-complete` — annotated tag on the final reviewed commit, so the "2026q2 catalog v0.4 + 14 collider_rs" milestone has a permanent ref even after the branch is gone.

### A.3 Website branch decision

**Recommendation: keep `flavor-catalog-website/2026q2` as the deploy branch (option i).**

Rationale:
- Cloudflare deploys from it; reconfiguring requires (a) a working Cloudflare login, (b) changing the production target which is risky during a merge sprint, and (c) potentially re-running website CI under new branch protections.
- The website branch's diff vs trunk is **strictly additive in `flavor_catalog/website/` only** (zero files modified outside that subtree). This means merging trunk -> website is mechanically conflict-free.
- Post-consolidation workflow: after each significant catalog/data change on `main`, run `git checkout flavor-catalog-website/2026q2 && git merge main && git push`.

Future work (out of scope for this round): reconfigure Cloudflare to deploy from `main/flavor_catalog/website/` so the second branch can be retired. Note in `ISSUES.md` as `INFRA-1`.

### A.4 Consolidation command sequence

**Phase 0 — Pre-merge user gate for the 8 untracked files** (BLOCKING — must precede everything else):

The 8 untracked files are NOT auto-ignored by the current `.gitignore`. Root `.gitignore` line 71 reads `!artifacts/**` (an un-ignore), so `git add artifacts/` would sweep them all in. The orchestrator MUST resolve each file's disposition with the user before any other phase runs, and record the decision in `.orchestration/PRE_MERGE_DECISIONS.md`.

Decision matrix (each row = one user question, with planner's recommended default):

| # | Path | Recommended action | Rationale |
|---|------|--------------------|-----------|
| 1 | `scripts/export_collaborator_5tev_points.py` | **commit** | Reproducible export script useful to peers; pairs with artifact #3/#4. |
| 2 | `scripts/export_collaborator_direct_affine_points.py` | **commit** | Same. |
| 3 | `artifacts/collaborator_5tev_points.csv` | **commit** | Structured collaborator export with provenance JSON; matches `claim_level: operational_scan_only` schema. |
| 4 | `artifacts/collaborator_5tev_points.provenance.json` | **commit** | Provenance for #3. |
| 5 | `artifacts/collaborator_direct_affine_5_10tev_points.csv` | **commit** | Same family as #3. |
| 6 | `artifacts/collaborator_direct_affine_5_10tev_points.provenance.json` | **commit** | Provenance for #5. |
| 7 | `artifacts/sample_points_may18.csv` | **rm OR move to `snapshots/2026q2/`** | Different column schema from #3/#5, no provenance JSON, pattern-matches earlier scratch output. Reviewer-confirmed not peer-relevant. |
| 8 | `artifacts/sample_points_may18.txt` | **rm OR move to `snapshots/2026q2/`** | Free-form Markdown-flavoured description, not structured. Same reasoning as #7. |

For each file the user must choose one of `commit`, `gitignore-add`, or `rm` (with `rm` allowed only after user types the file name back to confirm). If the choice is `gitignore-add`, the orchestrator MUST add an explicit deny-line for that exact path **after** the `!artifacts/**` re-include line (because `!` rules are order-sensitive in `.gitignore`):
```
# Pre-merge 2026-05-25: explicit deny under !artifacts/**
artifacts/sample_points_may18.csv
artifacts/sample_points_may18.txt
```

Output of this phase: `.orchestration/PRE_MERGE_DECISIONS.md` with one resolved row per file. No subsequent phase runs until that file exists and all 8 rows are resolved.

**Phase 0.5 — Worktree prune + repo-state snapshot** (before safety tagging so the snapshot reflects the post-prune, post-decision state):

```
git worktree prune
git worktree list   # expect ONLY the main worktree at the repo root
# write .orchestration/PRE_MERGE_STATE.md containing:
#   git status
#   git branch -avv
#   git tag -l
#   git worktree list
#   git log -1 --format='%H %s' HEAD
```

**Phase 1 — Apply user decisions from Phase 0** (commits if any):
```
# (only run the relevant commands for files marked 'commit')
git add scripts/export_collaborator_5tev_points.py scripts/export_collaborator_direct_affine_points.py
git add artifacts/collaborator_5tev_points.csv artifacts/collaborator_5tev_points.provenance.json
git add artifacts/collaborator_direct_affine_5_10tev_points.csv artifacts/collaborator_direct_affine_5_10tev_points.provenance.json
git commit -m "chore(artifacts): commit collaborator export scripts and 5/10 TeV provenance bundles"

# (only if any files marked 'gitignore-add')
# Edit .gitignore to append explicit deny lines under the !artifacts/** block, then:
git add .gitignore
git commit -m "chore(gitignore): explicitly exclude scratch export iterations under artifacts/"

# (only if any files marked 'rm')
rm artifacts/sample_points_may18.csv artifacts/sample_points_may18.txt   # only the files marked rm
```

**Phase 2 — Safety net** (run after Phase 1, BEFORE any review):
```
git tag -a safety/pre-merge-2026-05-25 -m "snapshot before merge+review consolidation" flavor-catalog/2026q2
git tag -a safety/ca-v2-ew-tail-399f116 -m "preserved backup snapshot" backup/ca-v2-ew-tail-399f116
git tag -a safety/ca-v2-ew-tail-combined-776eff9 -m "preserved backup snapshot" backup/ca-v2-ew-tail-combined-776eff9
git push origin safety/pre-merge-2026-05-25 safety/ca-v2-ew-tail-399f116 safety/ca-v2-ew-tail-combined-776eff9

# Verify push succeeded — REQUIRED before proceeding. Abort if empty.
git ls-remote --tags origin | grep safety/pre-merge-2026-05-25 || { echo "ABORT: safety tag not on origin"; exit 1; }
git ls-remote --tags origin | grep safety/ca-v2-ew-tail-399f116 || { echo "ABORT: backup tag #1 not on origin"; exit 1; }
git ls-remote --tags origin | grep safety/ca-v2-ew-tail-combined-776eff9 || { echo "ABORT: backup tag #2 not on origin"; exit 1; }
```

The orchestrator must abort if any of these greps return empty. The `progress.json` field `safety_tags_pushed` flips to `true` only after the three greps pass.

**Phase 3 — Run review queue** (see Section B); only after every review unit is sealed do we proceed.

**Phase 4 — Fast-forward main**:
```
git checkout main
git pull --ff-only origin main           # sync any drift
git merge --ff-only flavor-catalog/2026q2
git push origin main
```

**Phase 5 — Update website**:
```
git fetch origin flavor-catalog-website/2026q2:flavor-catalog-website/2026q2
git checkout flavor-catalog-website/2026q2
git merge --no-ff main -m "merge: trunk consolidation 2026-05-25"
# expect zero conflicts; verify diff is empty outside flavor_catalog/website/
git push origin flavor-catalog-website/2026q2
```

**Phase 6 — Delete branches**. Order: remote first (so a failed local delete cannot leave origin polluted), then local. Backups already preserved as tags in Phase 2. Switch back to `main` before any local deletion so we are not on the branch we are about to delete.

```
git checkout main   # ensure we are not on a branch about to be deleted

# Remote deletion (only branches present on origin). Wrap each in a guard so one rejection
# (branch protection, race) does not abort the loop.
for b in audit/bag-inputs audit/cfw-comparison audit/wilson-rg \
         docs/paper-readiness docs/scope-wording fix/pdg-benchmarks \
         invalidation/gate-rerun paper/quark-scan-2026q2 scan/zero-pass-statistics \
         flavor-catalog/2026q2 ; do
  git push origin --delete "$b" 2>&1 | tee -a .orchestration/logs/remote_deletion.log || \
    echo "WARN: remote delete failed for $b — continuing" | tee -a .orchestration/logs/remote_deletion.log
done

# Local deletion (all 17 feature/backup/docs branches + trunk label)
for b in audit/bag-inputs audit/cfw-comparison audit/wilson-rg \
         backup/ca-v2-ew-tail-399f116 backup/ca-v2-ew-tail-combined-776eff9 \
         ca_wave1_top_higgs_ew_v2_push docs/canonical-manifest \
         docs/paper-readiness docs/scope-wording fix/pdg-benchmarks \
         invalidation/gate-rerun paper/quark-scan-2026q2 pka-T017 \
         pka-l008-push scan/zero-pass-statistics wa-w5a-top-higgs-ew-v2-152230 \
         wa_w5b_charm_edm_work flavor-catalog/2026q2 ; do
  git branch -D "$b" 2>/dev/null || true
done
```

**Phase 7 — Tag milestone + push tags** (single `--tags` push covers v2026q2-catalog-complete; no redundant individual push):
```
git tag -a v2026q2-catalog-complete -m "PRIMARY 94 + SECONDARY 8; post-review consolidation"
git push --tags origin
```

**Phase 8 — Verify GitHub mirrors local**:
```
git remote prune origin
git fetch --all
git branch -a            # expect: main, flavor-catalog-website/2026q2 (+remote mirrors)
git tag -l 'safety/*' 'v2026q2*'
```

**Phase 9 — Mirror local `ISSUES.md` to GitHub Issues** (user request item 5):
```
# For each entry in .orchestration/ISSUES.md that is NOT marked ACCEPTED-RISK or INFO-only,
# create a GitHub issue. Severity becomes a label.
gh issue create --title "[R##-I#] <one-line desc>" --label "severity:<crit|high|med|low>" --body "<copy from ISSUES.md + link to .orchestration/reviews/R##.md>"
# A simple awk/jq pass over ISSUES.md is sufficient. Resulting issue URLs append back into ISSUES.md
# as "Mirrored: <url>".
```

**Phase 10 — Update CLAUDE.md**:
```
# Replace lines 12 and 15 (paper/quark-scan-2026q2 references) with "main".
# Keep the lepton-sector caveat.
git add CLAUDE.md
git commit -m "docs(claude): point active-paper-branch reference at main after consolidation"
git push origin main
```

### A.5 Conflict expectations

- Trunk -> main: pure fast-forward, no conflicts.
- Trunk -> website: conflict-free by structural argument (disjoint trees). If a conflict appears, investigation required (probably means someone hand-edited a catalog `.tex` on the website branch).
- Branch deletions: none expected since all unique work is captured.

### A.6 Rollback runbook

See `.orchestration/ROLLBACK.md` (orchestrator writes this in Phase 0.5). Pre-canned commands:
```
# Full revert to the safety tag (destroys all post-tag work):
git checkout main
git reset --hard safety/pre-merge-2026-05-25
git push --force-with-lease origin main

# Restore a deleted feature branch (requires the safety tag and one commit SHA from it):
git push origin <sha>:refs/heads/<branch>

# Restore the website branch:
git push origin safety/pre-merge-2026-05-25:refs/heads/flavor-catalog-website/2026q2
```

---

## B. Atomic-review breakdown

Target unit size: each unit reviewable in one focused Opus run (~30-60 minutes wall time).
**Total: 24 units (R01-R09 + R10a/R10b/R10c + R11-R22).** Spans the 369 commits ahead of `main` plus 13 website-only commits.

The atomic seams used:
- **Code seams** (R01-R09): each `feat/`, `fix/`, `physics/`, `audit/`, `figures/`, `chore/`, `test/`, `scripts/` commit cluster on the paper rc1 work (Phase 1, Phase 2 holes #4-#10, Phase 3) is one unit. These are the only commits with executable code changes.
- **Catalog seams** (R10a/R10b/R10c, R11-R19): waves of the flavor catalog. Each wave is one unit (waves 1 — split by family, 2/3, 4, 5a, 5b, 6, 7, 8, 9 collider_rs, audits+master-compile). Catalog units have NO code; review focuses on physics-vs-literature and consistency across `processes/<family>/<id>.yaml|tex` files.
- **Docs/orchestration seams** (R20): catalog meta-docs, scaffolding, handoff notes.
- **Website seams** (R21-R22): two units for the Astro site (R21 = scaffold + ingest + content; R22 = UI polish + LaTeX/methodology + Cloudflare config).

### B.1 Unit table

| ID | Title | Commits (representative) | Files / area | What was added | Checks |
|----|-------|-------------------------|--------------|----------------|--------|
| R01 | Quark-scan foundation (feat: qcd, scan, quark, analysis) | `109ec02`, `d4f873c`, `d0a9103`, `bcca1df`, `d3bcac9`, `6ccf6d8`, `b0675cd` | `qcd/`, `quarkConstraints/`, `scanParams/`, `scripts/`, new `docs/quark_scan_methodology_note.tex` | PDG-2024 MS-bar quark running; RS-anarchy scan driver + dispatch; PDG-target gate integration; analysis utilities; methodology paper drafted. | physics, numerics, code |
| R02 | Phase-2 hole #4 — spurion seed xfail repair | `16e6b36`, `2871fce`, `559b851`, `fd83d96`, `7e679ea`, `1fec2c7` | `tests/test_*spurion*`, `quarkConstraints/`, `pytest.ini` | Captured pre-fix xfail; re-derived spurion seed against PDG; enabled strict xfail; locked passing fixtures. | numerics, code |
| R03 | Phase-2 hole #5 — bag-parameter audit (FLAG 2024 + BGS 2020) | `dc9c498`, `82a96f0`, `695f35e`, `2521f6b` | `quarkConstraints/` (kaon inputs), `docs/audits/` | Updated kaon bag parameters; documented hadronic-input provenance (`82a96f0` is the provenance-doc commit; assigned here ONLY). | physics, numerics, code |
| R04 | Phase-2 hole #6 — Wilson RG running (BMU) + sign-convention fix | `7f71908`, `c80a8e0`, `561d848`, `c83c16d` | `quarkConstraints/deltaf2.py`, `quarkConstraints/qcd_running.py`, `tests/test_qcd_running.py`, `tests/test_quark_deltaf2.py`, `tests/test_wilson_rg_audit.py`, `scripts/audit_wilson_rg.py`, `docs/audits/wilson_rg_*` | Fixed BMU LO Wilson RG; corrected BMU<->scalar-LR sign convention in ADM and matrix elements. **HIGH-RISK PHYSICS unit.** | physics, numerics, code |
| R05 | Phase-2 invalidation-gate rerun + post-audit figures | `29803ff`, `db02223`, `87c728f`, `38a4586`, `217af80`, `f19f17c` | `scripts/`, `scan_outputs/`, publication figs | Regenerated 8 RS-anarchy scans + plots under corrected DeltaF=2 constants; updated headline numbers. | numerics, code |
| R06 | Phase-2 hole #7 — CFW (Csaki-Falkowski-Weiler) comparison | `ffec986`, `06d5d85`, `da8f647`, `508ae69`, `b47aa76`, `330ffa9`, `11e0d58`, `321fe0f`, `7c284a8`, `2a3a2c0`, `f2c29c8`, `ad2e3ee` | `scripts/`, `tests/test_cfw*`, `figures/`, `docs/audits/` | Extracted CFW 0804.1954 conventions; added convention-override flags; regenerated CFW comparison; demonstrated factor-2.2 result (not 11%); regression tests. | physics, numerics, code |
| R07 | Phase-2 holes #8-#9 — zero-pass binomial bounds, scope appendix, canonical manifest | `f34f9df`, `a9d1058`, `212320f`, `1d9e17e`, `715b679`, `f911a13`, `bd8ded1`, `e359faa`, `1a107c8`, `760f23d`, `80641cf`, `cb27631`, `b86e66d`, `2729659`, `b946abf`, `5e60bc0`, `c9d9cb5` | `physics/stats` (new), `figures/zero-pass`, `docs/paper/`, `artifacts/` | Wilson-score UL helper for zero-pass; legend annotations; scope/approximations appendix; canonical-artifact manifest. | physics (UL bound only), numerics, code |
| R08 | Phase-3 hole #10 figure prune + final PDF rebuild | `e7d824d`, `3951f41`, `bf6186c`, `00b222d`, `e3c0c79`, `1faba23`, `dce4d18`, `4be26fd`, `e3a0f1e`, `e0b4e2` (rc1.1 fixes), `8541a2b`, `82ff25e`, `e2b8349`, `5471ccd`, `a4ee3be`, `2121e4e`, `8deb291`, `73d347d`, `fb37473`, `9a13d16` | `figures/`, `notebooks/`, `docs/paper/` | Pruned unreferenced figures; rebuilt PDF; rc1.1 text fixes from end-to-end Opus review; re-executed 4 notebooks under post-audit DeltaF=2 constants. | numerics (notebook re-exec), code |
| R09 | Discovery-mode catalog scaffold + planner v0/v1 | `ebe8a3c`, `744e70c`, `8981648`, `33367d9`, `95f9f63`, `c78fd73`, `83c0178`, `a516dc5` | `flavor_catalog/` initial scaffold | Defined plan-v1, scaffolded `flavor_catalog/processes/<family>/` directory structure. | code (structure only) |
| **R10a** | **Catalog Wave-1 kaon family** (K001-K006, K013) | PKAs `0155d7c`, `0a1ce41`, `6fae247`, `6fb9d78`, `e45eadc`, `21661a4`, `6027db0`; WA `50311f6`; CA `4e6b7a9`; WA-v2 `e5ba712`; CA-v2 `69c33e4` | `flavor_catalog/processes/kaon/{K001,K002,K003,K004,K005,K006,K013}.{yaml,tex}` + references + worklogs | 7 kaon PKA drafts (epsilon_K, Delta m_K, eps'/eps, K+->pi+nunu, KL->pi0nunu, KL->mumu, KL->pi0gamma gamma) + WA polish + CA verify + v2 rework. | physics, code (yaml/tex consistency) |
| **R10b** | **Catalog Wave-1 beauty family** (B002, B005, B009, B011, B015) | PKAs `c4cea8a`, `47e4339`, `a36f67e`, `fcbe870`, `9a92781`; WA `73b138a`; CA `9b49fae`; WA-v2 `25afeb0`; CA-v2 `6ffcdb2` | `flavor_catalog/processes/beauty/{B002,B005,B009,B011,B015}.{yaml,tex}` + references + worklogs | 5 beauty PKA drafts (S_psiKS, Bs->mumu, B+->tau nu, B->Xs gamma, inclusive b->s ll) + WA + CA + v2. | physics, code |
| **R10c** | **Catalog Wave-1 top/EW + EDM + charm + lepton** (T001, T002, T010, E001, C001, L001) | PKAs `18fe48f`, `be62498`, `9bb58af`, `2f443a0`, `6e60d02`, `9068956`; WA `1cf18c6`; CA `333042f`; WA-v2 `fb9a170`; CA-v2 `00393cd` | `flavor_catalog/processes/{top_higgs_ew,edm_neutrino,charm,charged_lepton}/...` + references + worklogs | 6 PKA drafts across 4 families: t->cZ, t->uZ, Zbb pole, electron EDM, D0 mixing, mu->e gamma. WA/CA cycles for top_higgs_ew family only (the other 3 families have a single representative each — verified at PKA level only). | physics, code |
| R11 | Catalog Wave-2/3 — broader rollout | `f902c33`, `56a0a7f`, `282fb9c`, `7b2c2c0`, `93c5a85`, `494309d`, `27de127`, `f93ed80`, `a2c7802`, `19a8a23`, `0fe9e31`, `5031e06`, `8795419`, `76b1fd6`, `bcb34b3`, `8b6388c`, `df697bc` + v2/v3 reworks | `flavor_catalog/processes/*` (T007, T018, B025, B032, B033, K001, K002, C001, E001, L001, L002, L007 + sweeps) | Added ~12 more processes; DA-1 discovery worklog. | physics, code |
| R12 | Catalog Wave-4 — 12+ new processes incl. EW oblique, R_K*, R_D, R_D*, mu->e conv | `f048312`, `7d3da08`, `de78fa1`, `741f7c0`, `f734e4e`, `fd49402`, `cf21e2b`, `519aec5`, `8eec2a0`, `d64d457`, `66ea613`, `fbbc3ff`, `a91f599`, `0c395f2`, `49c57c2`, `65e3ca5` + W4 WA/CA batches `1268def`, `7741e1c`, `eae8d05`, `470d498`, `0f5972b`, `bcd8907`, `a5e94cb`, `aab1862`, `c369747`, `4907544`, `a23b929`, `ff366bc`, `52acd5e`, `f863013`, `fabad04`, `049fa43`, `d6e78b6`, `8a0d6fa` | broad sweep across families | EW001 (oblique S,T), B019, B026, EW002 (CKM first-row), C002 (D mixing CP), E006, E008, K017, L003, L004, L008, L009, L010; DA-2 worklog. | physics, code |
| R13 | Catalog Wave-5a — second-round PKA + v2 (top, beauty, charged-lepton EDM, kaon-charm) | `8aa7b96`, `74ff5bc`, `121559e`, `92f3286`, `8d02741`, `a1d969e`, `82780c6`, `ca886e0`, `fd3c97e`, `140ab97`, `de608a6`, `a4519a6`, `d5dc2a6`, `74c6aa6`, `a873420`, `6a978ea` + v2 reworks `8ce0306`, `45f051a`, `a496085`, `16ef205`, `b55bed6`, `90b60ea`, `058c7f1`, `bb942a1`, `7618cd0` | broad sweep | T005, T015, T019, T006, T016, T017, K008, K009, K010, C005, C006, C008, E007, B021, B022, B023, B034, L005, L006, E002; DA-3 worklog; Opus round-1 sign-off on 50 processes. | physics, code |
| R14 | Catalog Wave-5b — completion of wave-5 followups | `7849dcc`, `d288119`, `4eca813`, `a9dd56c`, `0fcb364`, `51fe52f`, `3b76b20`, `4f07abf`, `923b8cd`, `64b5043`, `d596b3e`, `17ab228`, `bc4f253`, `f6ab87a`, `f115228`, `7f87327`, `b709912` | wave-5 leftovers, codex-quota pause notes | T006 v2, B001, B016, K012, K018; codex quota pause + resumption plan. | physics, code |
| R15 | Catalog Wave-6 — new physics drafts (E009, L023, T020, K018, K012, B001, B003, B016) | `398b939`, `50793c3`, `a28aa73`, `297b820`, `aaa22d6`, `a15edea`, `bdb8371`, `c6ec949`, `e821611` + v2 batches `56ee949`, `e9d13b7`, `4666517`, `ee5c95b`, `67912c0`, `cd8816d`, `7500919`, `6d97db8` | broad sweep | E009 (Weinberg 3-gluon), L023 (neutrino trident), T020 (h->emu), K012, K018, B001, B003, B016; DA-4 convergence; lock at 75 processes + deferred-scope tail. | physics, code |
| R16 | Catalog Wave-7 — 5 new processes + cross-cutting subtleties + arbitration | `e3ecf8b`, `55e2299`, `64d3f6a`, `1f9e16a`, `60e85b9`, `dbf2b96`, `837e709`, `99003f3`, `40248b8`, `5f21f3e`, `52f6660`, `3c21e30`, `3f46ca0`, `91358f9`, `0cfc3ed`, `a6ffc9f`, `ca2beba`, `7a6bd34`, `1aa6b37`, `7e1b80b` | new processes + DA-4 addendum | T003 (t->c gamma), T004 (t->u gamma), T008 (t->Hu), T012 (Z->ccbar), B012 (B->K* gamma); 6 cross-cutting subtleties threaded through ~22 .tex files; deferred-scope addendum; arbitration on B001+B003. | physics, code |
| R17 | Catalog master compile v0.2 + audit/factcheck consolidation (75 -> 80 processes) | `835cf48`, `9f65578`, `aa8e3e8`, `a4245e2`, `9fb5763`, `a04a101`, `7f79bab`, `0b07505`, `04b1ab1`, `994e346`, `6498fad`, `42ac647`, `022a20c`, `8fb5f91`, `2c00d84` | `flavor_catalog/catalog_master.tex`, `flavor_catalog/audits/factcheck_status.md`, `flavor_catalog/external_research/` | First master compile (75 + Wave-7 = 80); family-by-family factcheck (73-75 VERIFIED); imported external GPT Deep Research PDFs; v0.1->v0.2 milestone. | physics (audit cross-check), code |
| R18 | Catalog Wave-8 — 8 SECONDARY tier processes (T014, K019, K020, K021, B007, B008, B013, B014) | `82a2c01`, `2630168`, `bab5bd0`, `37beabb`, `ebdo66c`, `e6e1cc3`, `5b68b43`, `6b64a18`, `e348a35`, `17599d5`, `0ad4a07`, `4a6df17`, `103833f`, `d198787`, `2b23464`, `f9c7ff4`, `88f2cde`, `cd2f81b`, `b3c35d7`, `e9f3cf3`, `41069e7`, `749d789`, `5c14d2f`, `a20d75a`, `980b45b`, `118de1c`, `37514e9`, `7e121c3`, `df7af95`, `ef5a417`, `76c87ae`, `fff1898`, `6a2506c`, `b5c2375` | `flavor_catalog/processes/secondary/`, master compile v0.3 | 8 SECONDARY processes drafted/polished/verified (PRIMARY=80, SECONDARY=8); K020 cycle-3 metadata cleanup; Wave-8 close-out, v0.3 tag. | physics, code |
| R19 | Catalog Wave-9 — 14 collider_rs PRIMARY processes (CR001-CR014) + master compile v0.4 | `7ed9117`, `108cae4`, `5911918`, `f195c45`, `144cc61`, `4054e52`, `7614f6d`, `6176694`, `0a8e20a`, `28879ea`, `d09c6b3`, `835cc83`, `92b78c7`, `c22fb07`, `42428d7`, `1cf8b57`, `a13edbb`, `b96036b`, `1eac13e`, `c6d55b0`, `2286a39`, `cd3a8fe`, `e8daa00`, `7a6fa3d`, `bc35a06`, `e846f14`, `a8758ac`, `e1aec33`, `950ca36`, `0c5dacc`, `cdd0238`, `2ad34b1`, `864cd6d` (33 commits; `82a96f0` belongs to R03 ONLY and is intentionally NOT listed here) | `flavor_catalog/processes/collider_rs/` (CR001-CR014), wiring into catalog_master | NEW family `collider_rs/` with 14 collider PRIMARY entries (RS-direct collider constraints); 14 PKA + 4 WA + 4 CA cycles + Opus round-5 sign-off (14/14 APPROVE); master compile v0.4 = PRIMARY 94 + SECONDARY 8. | physics, code |
| R20 | Catalog meta-orchestration docs | `5086c99`, `5a6fae0`, `3be3a7b`, `24284d7`, `6f072d5`, `d1bbaf4`, `781bc20`, `2bda5f1`, `e3b37c7`, `6d99b17` | `flavor_catalog/AGENTIC_WORKFLOW.md`, `SESSION_NOTES.md`, `HANDOFF_PROMPT.md`, `PRIORITY_TIERS.md`, `CATALOG_METHODOLOGY.tex`/`pdf`, `WEBSITE_BUILD_PROMPT.md` | Workflow playbook + handoff prompts + tiering policy + one-page methodology PDF + website cold-boot prompt. | code (consistency only — no physics) |
| R21 | Website — scaffold, ingest, citation anchors, phases 1-5 (origin only) | `460ece9`, `6ffbb33`, `d62db33`, `5f087fd`, `344e147` | `flavor_catalog/website/` (Astro app + content/entries + citation anchors + family pages + Cloudflare config) | New Astro static site; 102-entry citation-anchor resolution; family pages; search/filters; Cloudflare config; README; methodology page. | code (consistency only — physics already reviewed at the YAML/tex source level) |
| R22 | Website — phases 6-12 polish (LaTeX normalization, priority UI, UI cleanup) | `5180988`, `eb6b8b7`, `0833d0e`, `7053fb7`, `4eba297`, `b6864f7`, `f485725`, `5f31f2d` | `flavor_catalog/website/src/lib/notation.ts`, `prose.ts`, `pages/`, runbook | Shorthand->LaTeX normalizer; constraint-priority ranking; site-wide LaTeX/em-dash/jargon sweep; methodology rewrite; Secondary nav removal; PDG column fixes. | code |

**Granularity check (round 2)**: Code units (R01-R09) average ~5 commits, R10a/b/c are 9-11 commits each, R11-R19 are ~20-35 commits with structured Writer-Checker-Opus internal gating, website units 5-8 commits. R19 is now correctly 33 commits (down from 34 after the `82a96f0` removal). Each unit fits within one Opus run.

**Uniqueness invariant**: when `COMMIT_INDEX.csv` is generated, the orchestrator MUST run:
```
awk -F, 'NR>1 {print $1}' COMMIT_INDEX.csv | sort | uniq -d
```
and abort if non-empty. Any SHA must appear in exactly one unit.

### B.2 Units that are explicitly out-of-scope (no review needed)

- Pure `docs(paper): seal Phase ... signoff` commits — these are review-of-review artifacts; the underlying work was already covered in an R-unit, so the seal commits are bookkeeping.
- Pure `chore(gitignore)` commits — captured in R07/R08 context only.

### B.3 Mapping commits -> units (machine-readable index)

The orchestrator will generate `.orchestration/COMMIT_INDEX.csv` with columns `commit_sha, unit_id, role` after this plan is approved, so each reviewer agent receives the exact `git log --oneline RANGE` for its unit. Alongside it the orchestrator generates per-unit `.orchestration/pytest_selection/<unit>.txt` with the list of `tests/*` files touched in that unit's commits (one `git log --name-only RANGE -- tests/` pass per unit), per RC-3.

---

## C. Per-unit review protocol

### C.1 Reviewer agent contract

Each Opus reviewer is dispatched with:

**Inputs provided by orchestrator** (passed in the prompt):
1. Unit ID (e.g. `R04`), title, and the row from §B.1.
2. The list of commit SHAs for the unit, with one-line subjects.
3. The unified diff: `git log --reverse -p <range>` truncated to a sensible cap, OR the file list with the agent free to `git show <sha>` per commit.
4. A pointer to `flavor_catalog/CATALOG_METHODOLOGY.tex` and the relevant section of `docs/quark_scan_methodology_note.tex` if the unit touches paper figures/numerics.
5. **Literature manifest** — `.orchestration/REFERENCES.md` (orchestrator generates this BEFORE first reviewer dispatch; see §C.2.1) with one section per unit listing exact arXiv IDs, section/eq numbers, and local PDF paths under `flavor_catalog/external_research/` where available.
6. **Per-unit pytest selection** — `.orchestration/pytest_selection/<unit>.txt`.
7. The path where the report goes: `.orchestration/reviews/<unit_id>.md`.
8. A reminder that this is read-only: the reviewer may run code (`pytest`, `python -c '...'`) to verify numerics but may NOT modify files.
9. The reviewer MAY use `gh issue list --search "..."` and `gh pr list` to cross-reference open GitHub issues touching the unit (per M-6).

**Reviewer mandate clarification (per M-6)**: The reviewer's job is NOT to re-run the project's internal Writer-Checker-Opus loop that already gated each catalog process. The reviewer's job is to (a) cross-check selected high-impact outputs against external literature, and (b) verify cross-cutting code/schema consistency that single-process reviewers cannot see.

**Outputs the reviewer must produce**:
- One report file `.orchestration/reviews/R##.md` following the template in §C.3.
- One append-only block in `.orchestration/ISSUES.md` for each issue found (severity tag + unit ID + 1-line description + pointer to report section).
- An update to `.orchestration/REVIEW_QUEUE.md` flipping the unit from `IN-PROGRESS` to `DONE` (or `FAIL` if a blocker is found).

### C.2 The three checks

#### C.2.1 Physics intuition grounded in literature

The orchestrator pre-writes `.orchestration/REFERENCES.md` with one section per unit. Each section lists exact arXiv IDs, exact eq/section numbers to consult, and any local PDF path the project already imported under `flavor_catalog/external_research/`. The reviewer is NOT expected to scavenge papers from the open internet.

**For code units (R01-R08)** — quark-scan paper machinery:
- Check formulas in `quarkConstraints/deltaf2.py`, `quarkConstraints/qcd_running.py`, `quarkConstraints/bsm_amplitudes.py` against:
  - Agashe et al. hep-ph/0412089 (RS flavor framework — main reference).
  - Csaki-Falkowski-Weiler 0804.1954 (RS flavor protection — explicitly cited in R06).
  - FLAG 2024 averages (the source `dc9c498` updated to).
  - Buras-Misiak-Urban hep-ph/0005183 (BMU operator basis — referenced in R04).
  - PDG 2024 (quark masses, CKM elements, DeltaF=2 observables).
- For R04 specifically: independently re-derive the BMU<->scalar-LR Wilson coefficient sign convention from the operator definitions and confirm against the literature; this was the subject of a discovered sign-flip bug, so the agent must NOT just trust the post-fix code.
- For R07 (Wilson-score UL): derive the 95% Clopper-Pearson and Wilson upper bound for 0/N pass observations and check agreement.

**For catalog units (R10a/b/c, R11-R19)** — each process YAML/tex:
- Sample 2-3 processes per unit (do NOT do all 14 in R19; pick high-impact: the rare-decay BR, the EDM bound, the mixing observable).
- **Tightened rule (per OQ-4 reviewer answer)**: the spot-check MUST include at least one PRIMARY-tier process from EACH family touched in the wave. Example: if a wave touches 3 kaon + 1 beauty + 1 charm, the sample is {1 kaon, 1 beauty, 1 charm} minimum, not {3 kaon, 0, 0}. R10 split mitigates this for Wave-1; later waves still need the family-coverage rule. The reviewer must document which processes were spot-checked in the report under "Sample selection".
- For each sampled process verify:
  - `pdg_or_equivalent` value matches the cited reference (PDG 2024, HFLAV 2024, LHCb/CMS/ATLAS papers, NA62, Belle II).
  - The `standard_notation` LaTeX expression matches a published convention.
  - The `theory_uncertainty` and `experimental_uncertainty` ranges are consistent with the literature (no obvious order-of-magnitude errors).
- For collider_rs (R19): cross-check the CMS/ATLAS contact-interaction limits against the cited ATLAS/CMS arXiv papers (the wave9 worklog has the exact SHA-256 of the source TXTs — the agent should `cat flavor_catalog/processes/collider_rs/CR0XX.yaml` and confirm the citation-to-value mapping).

**For docs/scaffold units (R09, R20, R21, R22)**: physics check is `N/A`. State this explicitly in the report.

#### C.2.2 Numerical correctness

**For code units**:
- Run pytest using the per-unit selection file `.orchestration/pytest_selection/<unit>.txt` (auto-generated from `git log --name-only RANGE -- tests/`, plus a default suite of `tests/test_qcd_running.py tests/test_quark_deltaf2.py tests/test_wilson_rg_audit.py tests/test_cfw_comparison.py` for high-risk units). All must pass.
- **Test-timeout policy (per M-5)**: if pytest for the unit's selection exceeds 10 min wall time, the reviewer marks numerics as PARTIAL with a note and continues; this is NOT a FAIL. (`pytest --collect-only` takes ~100s in this repo, so import overhead alone eats budget.)
- For R04 specifically, re-execute `scripts/audit_wilson_rg.py` and compare against `docs/audits/wilson_rg_reference_values.md`.
- For R05, spot-check the post-audit scan output checksums against the manifest in `artifacts/`.
- For R08, recompute one scalar in each re-executed notebook (e.g., the headline yield number) by hand against the printed cell output.

**For catalog units**: numerical check is mostly cross-referential — pick the 2-3 sampled processes from §C.2.1 and verify any arithmetic in the .tex narrative (e.g., a 4-sigma tension calculation, a 90% CL conversion to 2-sigma).

**For website units (R21, R22)**: numerical check is "do the displayed numerics match the source YAMLs?" — pick 3 random entries, compare browser-rendered (via `npm run build` + grep of generated HTML) numbers to source YAML. R22 specifically must also run `npm run build && grep -r '\\\\(' dist/ | head` to confirm LaTeX delimiters are rewritten (not raw-rendered).

#### C.2.3 Repo-wide code self-consistency

For every unit:
- `git grep` for any symbol or constant added/changed by the unit and confirm all call sites have been updated (e.g., for R04 the BMU sign fix must propagate to every consumer of those Wilson coefficients).
- `git grep` for the magic strings introduced (e.g., the `Lambda_IR` default, the `MEG II 2025 br_limit = 1.5e-13`) and confirm they appear consistently.
- For catalog units, run a YAML schema check: every new `processes/<family>/<id>.yaml` must parse and conform to whatever schema the catalog uses (the agent should look for a schema doc or sample from existing entries).
- For all code-touching units: `python -m py_compile <touched .py>` to catch syntax regressions.
- Check that `catalog_master.tex` wires include all new processes and compiles (it doesn't need to actually build the PDF, but `\input{...}` references must point at real files).

**Per-unit grep checklists** (the orchestrator extends `.orchestration/REFERENCES.md` with a `grep_targets:` block for each unit; the high-risk pre-canned lists are):
- **R04 (BMU sign fix)**: `evolve_deltaf2`, `C4_LR`, `C5_LR`, `bmu_to_scalar_lr`, `wilson_rg`. Confirm every consumer of these symbols uses the post-fix sign.
- **R06 (CFW comparison)**: `cfw_matched`, the literal `0.11`, the literal `2.2`, `factor_2p2`, the convention-override flag name (look for `cfw_convention` or `use_cfw_*` in `quarkConstraints/`).
- **R07 (zero-pass UL)**: `wilson_score`, `clopper_pearson`, `zero_pass_ul`, `0_pass_ci`.
- **R19 (collider_rs)**: each `CR0XX` ID, the `claim_level` field, the `pdg_or_equivalent` field, the SHA-256 listed in each `references/CR0XX/sha256sums.txt`.
- **MEG II constant (per OQ-7 reviewer answer)**: every unit that touches a charged-lepton process must grep for `1.5e-13`, `4.2e-13`, `2.0e-13` and flag (HIGH severity) any inconsistency with the frozen `br_limit = 1.5e-13` value. The reviewer does NOT propose code changes; the user has explicitly frozen this number.

### C.3 Report template `.orchestration/reviews/R##.md`

```
# Review R## — <Title>

**Reviewer**: Claude Opus 4.7 (single-run reviewer)
**Date**: YYYY-MM-DD
**Commits**: <range or list of SHAs>
**Files**: <high-level path list>
**Sample selection** (catalog units only): <which processes were spot-checked and why>

## Summary

<2-4 sentence summary of what this unit added>

## Check 1 — Physics intuition grounded in literature

**Verdict**: PASS | PARTIAL | FAIL | N/A
**Literature anchors consulted**: <bibliography list from REFERENCES.md>
**Evidence**:
- <bullet>
- <bullet>

## Check 2 — Numerical correctness

**Verdict**: PASS | PARTIAL | FAIL | N/A
**Commands run**: <list — include pytest selection used>
**Wall time**: <minutes> (mark PARTIAL if >10 min per M-5)
**Evidence**:
- <bullet>

## Check 3 — Repo-wide code self-consistency

**Verdict**: PASS | PARTIAL | FAIL | N/A
**Greps performed**: <list — include all targets from grep checklist>
**Evidence**:
- <bullet>

## Issues found

| ID | Severity | Description | File:line |
|----|----------|-------------|-----------|
| R##-I1 | CRITICAL/HIGH/MEDIUM/LOW/INFO | <desc> | <path:line> |

(severity: CRITICAL = wrong physics or wrong numerics; HIGH = correctness suspect; MEDIUM = consistency gap or missing test; LOW = doc/typo; INFO = noted for follow-up)

## Verdict

**Overall**: APPROVE | APPROVE-WITH-NOTES | BLOCK
```

### C.4 Issue tagging into `.orchestration/ISSUES.md`

`ISSUES.md` is grouped by severity:
```
## CRITICAL
- [R04-I1] BMU sign fix tested but not propagated to notebook X — flavor-catalog/notebooks/...
## HIGH
- ...
## MEDIUM
- ...
## LOW
- ...
## INFO
- [INFRA-1] Cloudflare deploys from flavor-catalog-website/2026q2; reconfigure to main later.
```

Each issue line has `[unit-Inum]` (or `[INFRA-num]` for orchestration-level) followed by a 1-line description and pointer. After Phase 9 each non-ACCEPTED-RISK non-INFO issue is mirrored to GitHub Issues; the URL is appended back as `Mirrored: <url>`.

A `BLOCK` verdict on any unit pauses the consolidation: the orchestrator must surface the issue to the user, await direction, then either patch on a `fix/<issue-id>` micro-branch (re-reviewed before merge) or accept-with-note (recorded in `ISSUES.md` as `ACCEPTED-RISK`).

---

## D. Tracking and resumability

`.orchestration/` layout (all files are append-only or replace-in-place; safe to inspect mid-run):

```
.orchestration/
├── MERGE_PLAN.md                # this file; planner outputs
├── PLAN_REVIEW_R1.md            # reviewer round-1 critique
├── PLAN_CHANGES_R2.md           # planner round-2 changelog
├── PRE_MERGE_DECISIONS.md       # Phase 0 user decisions on the 8 untracked files
├── PRE_MERGE_STATE.md           # Phase 0.5 snapshot (git status/branch/tag/worktree)
├── REFERENCES.md                # per-unit literature manifest + grep checklists
├── ROLLBACK.md                  # pre-canned revert commands
├── REVIEW_QUEUE.md              # status board (one row per unit)
├── COMMIT_INDEX.csv             # commit_sha -> unit_id mapping (auto-generated; uniqueness enforced)
├── pytest_selection/
│   ├── R01.txt
│   └── ...                      # one file per unit
├── ISSUES.md                    # running issue log, grouped by severity
├── reviews/
│   ├── R01.md
│   ├── R02.md
│   └── ...
├── progress.json                # machine-readable resumability state
└── logs/
    ├── remote_deletion.log
    └── <unit>.<timestamp>.log   # raw reviewer stdout (optional)
```

### D.1 `REVIEW_QUEUE.md` schema

```
| ID   | Title                              | Status      | Started      | Finished     | Verdict           | Issues |
|------|------------------------------------|-------------|--------------|--------------|-------------------|--------|
| R01  | Quark-scan foundation              | DONE        | 2026-05-25 09:10 | 2026-05-25 09:48 | APPROVE           | 0      |
| R02  | Phase-2 hole #4 — spurion seed     | IN-PROGRESS | 2026-05-25 09:50 | -            | -                 | -      |
| R10a | Wave-1 kaon family                  | PENDING     | -            | -            | -                 | -      |
| R10b | Wave-1 beauty family                | PENDING     | -            | -            | -                 | -      |
| R10c | Wave-1 top/EW + EDM + charm + lepton| PENDING     | -            | -            | -                 | -      |
```

Status values: `PENDING`, `IN-PROGRESS`, `DONE`, `BLOCKED`.

### D.2 `progress.json`

Minimal JSON the orchestrator updates after each unit:
```
{
  "schema": 2,
  "current_phase": "review",      // pre_merge_gate|worktree_prune|phase1_artifacts|safety|review|consolidate|delete|push|github_mirror|done
  "pre_merge_decisions_recorded": false,
  "worktree_pruned": false,
  "safety_tags_created": false,
  "safety_tags_pushed": false,
  "last_completed_unit": "R02",
  "next_unit": "R03",
  "blocked_units": [],
  "merge_done": false,
  "branches_deleted": false,
  "push_done": false,
  "github_issues_mirrored": false
}
```

### D.3 Resumability contract

If the orchestrator is killed mid-run, the next invocation:
1. Reads `progress.json` to determine where to resume.
2. Cross-checks `PRE_MERGE_DECISIONS.md` exists and all 8 rows resolved; if not, restart from Phase 0.
3. Cross-checks `REVIEW_QUEUE.md` for `IN-PROGRESS` rows (those need re-dispatch).
4. Cross-checks `git tag -l 'safety/*'` AND `git ls-remote --tags origin | grep safety/` to confirm safety nets exist BOTH locally and on origin (creates them if missing).
5. Proceeds from `next_unit`.

A reviewer agent that finishes its unit MUST commit nothing (the harness does that); it only writes report + queue update + issue log. The orchestrator is the sole entity that runs git mutations.

### D.4 Orchestrator commit policy (per M-3)

The orchestrator creates the following commits, all on `main` (or `flavor-catalog/2026q2` while it still exists), each authored as the user with `Co-Authored-By: Claude (orchestrator) <noreply@anthropic.com>`:

| Phase | Commit subject | Branch | Notes |
|-------|---------------|--------|-------|
| 1 | `chore(artifacts): commit collaborator export scripts and 5/10 TeV provenance bundles` | trunk | Only if any Phase-0 row chose `commit` |
| 1 | `chore(gitignore): explicitly exclude scratch export iterations under artifacts/` | trunk | Only if any Phase-0 row chose `gitignore-add` |
| 5 | `merge: trunk consolidation 2026-05-25` | website | Phase 5 trunk->website merge commit |
| 10 | `docs(claude): point active-paper-branch reference at main after consolidation` | main | Step 33 / Phase 10 |

Reviewer agents create zero commits.

---

## E. Order of operations

Default mode: **serial** (one Opus reviewer at a time) per the user's 5-hour budget constraint. Parallel speculation noted only where strictly safe.

| Step | Action | Tool / agent | Serial? |
|------|--------|--------------|---------|
| 1 | **Phase 0 pre-merge gate**: ask user about 8 untracked files (commit/gitignore-add/rm). Write `.orchestration/PRE_MERGE_DECISIONS.md`. | orchestrator + user | yes (BLOCKING) |
| 2 | **Phase 0.5**: `git worktree prune`; write `.orchestration/PRE_MERGE_STATE.md` and `.orchestration/ROLLBACK.md`. | orchestrator bash | yes |
| 3 | **Phase 1**: apply user decisions from step 1 (any commits, .gitignore edits, rm). | orchestrator bash | yes |
| 4 | **Phase 2**: create safety tags + push to origin + verify push with `git ls-remote --tags`. Set `safety_tags_pushed=true` only after verification. | orchestrator bash | yes |
| 5 | Write `.orchestration/REFERENCES.md` (literature manifest + grep checklists for all 24 units). | orchestrator | yes |
| 6 | Generate `COMMIT_INDEX.csv` from §B.1 table; assert uniqueness invariant. Generate per-unit `pytest_selection/<unit>.txt`. | orchestrator bash | yes |
| 7 | Initialize `REVIEW_QUEUE.md`, `ISSUES.md`, `progress.json`. | orchestrator bash | yes |
| 8 | Dispatch R01 reviewer. | Opus agent | yes |
| 9 | Append R01 issues to ISSUES.md, mark R01 DONE. | orchestrator | yes |
| 10-31 | Dispatch R02 through R22 (24 units total: R01-R09, R10a, R10b, R10c, R11-R22) one at a time, repeating step 9 each time. | Opus agent | yes |
| 32 | Triage `ISSUES.md`: any CRITICAL/BLOCK? If yes, surface to user; if all clear or only LOW/INFO, proceed. | orchestrator + user | yes (gating) |
| 33 | **Phase 4** — fast-forward `main` to trunk tip, push. | orchestrator bash | yes |
| 34 | **Phase 5** — merge `main` into `flavor-catalog-website/2026q2`, push. | orchestrator bash | yes |
| 35 | **Phase 6** — `git checkout main`; delete remote then local branches (guarded loop). | orchestrator bash | yes |
| 36 | **Phase 7** — tag `v2026q2-catalog-complete`, push tags via single `--tags` push. | orchestrator bash | yes |
| 37 | **Phase 8** — `git remote prune origin && git fetch --all`, final verification. | orchestrator bash | yes |
| 38 | **Phase 9** — mirror non-ACCEPTED-RISK non-INFO issues to GitHub Issues via `gh issue create`. | orchestrator bash | yes |
| 39 | **Phase 10** — update CLAUDE.md to remove stale `paper/quark-scan-2026q2` reference (one small commit on main). | orchestrator bash | yes |

**Parallelization opportunities (NOT recommended now)**:
- R10a/R10b/R10c are independent (disjoint families) and COULD run in parallel as a 3-way batch — would cut Wave-1 wall time ~3x. Default to serial.
- R11-R19 are also pairwise independent but together cover 200+ commits; parallelizing would burn ~10x token budget.
- R21 and R22 (website units) are independent — same trade-off.
- The 9 code units (R01-R09) have implicit dependencies (R04 informs R05, R06 informs R07/R08) and should stay strictly serial.

Estimated wall time at ~45 min/unit × 24 units = ~18 h serial. **This exceeds 5 h**. The orchestrator should plan to dispatch units in batches across multiple user sessions, snapshotting via `progress.json`. The first session should aim to cover R01-R08 (the code units with executable consequences) so that all the high-physics-risk review is done first.

---

## F. Open questions for the reviewer agent — STATUS

| OQ | Question | Reviewer verdict (R1) | Incorporated in R2 plan |
|----|----------|------------------------|-------------------------|
| 1 | Untracked-file disposition? | Commit the 4 collaborator artifacts + 2 provenance JSONs + 2 export scripts; do NOT commit `sample_points_may18.*` (move or rm). User-gated. | YES — new Phase 0 gate (§A.4) |
| 2 | Backup branches -> tags? | Correct; canonical replacement is `82daa9b` (NOT the planner's typo `82b35c0`). | YES — §A.1 row 4/5 corrected |
| 3 | Split R10? | Split. | YES — R10a/R10b/R10c (§B.1) |
| 4 | Catalog physics depth? | Spot-check default OK, but REQUIRE one PRIMARY per touched family per wave; document which processes checked. | YES — §C.2.1 tightened |
| 5 | Collapse R21+R22? | Keep separate; add `grep -r '\\\\(' dist/` check to R22. | YES — §C.2.2 |
| 6 | CLAUDE.md update timing? | Post-merge is fine (Step 33 / Phase 10). | YES — §E step 39 |
| 7 | MEG II `1.5e-13` scope? | Applies to code AND docs; flag inconsistencies (HIGH) but propose no code change. | YES — §C.2.3 grep checklist |

---

## Appendix A — Verification snippets the orchestrator can re-run

```
# Subsumption matrix
for b in $(git for-each-ref --format='%(refname:short)' refs/heads refs/remotes/origin | grep -v -E 'HEAD|main$|flavor-catalog/2026q2$|flavor-catalog-website/2026q2$'); do
  ahead=$(git log --oneline flavor-catalog/2026q2..$b 2>/dev/null | wc -l)
  printf "%-50s ahead-of-trunk:%d\n" "$b" "$ahead"
done

# Website branch isolation check
git diff --name-only flavor-catalog/2026q2..origin/flavor-catalog-website/2026q2 -- ':!flavor_catalog/website/'
# expect empty output

# Worktree audit before deletion (Phase 0.5)
git worktree list
git worktree prune
git worktree list      # expect only main worktree

# SHA verification (Phase 0/A.1)
git log -1 --format='%h %s' 82daa9b   # must succeed
git log -1 --format='%h %s' 82a96f0   # must succeed; subject = "docs(paper): document hadronic input provenance"
git log -1 --format='%h %s' 82b35c0   # MUST FAIL (does not exist; was the planner's R1 typo)

# Safety-tag push verification (Phase 2)
git ls-remote --tags origin | grep safety/pre-merge-2026-05-25

# COMMIT_INDEX uniqueness
awk -F, 'NR>1 {print $1}' .orchestration/COMMIT_INDEX.csv | sort | uniq -d
# expect empty
```

## Appendix B — Risk register

| Risk | Likelihood | Severity | Mitigation |
|------|-----------|----------|------------|
| Loss of unique work during branch deletion | Low | High | Phase 2 safety tags pushed AND verified before any deletion; subsumption matrix proves all unique work is preserved. |
| Cloudflare deploy breaks after website merge | Low | Medium | Trunk-into-website merge is structurally conflict-free; if Cloudflare fails, revert via the safety tag is one command (see `ROLLBACK.md`). |
| Reviewer agent finds CRITICAL physics bug | Medium | High | Orchestrator pauses and surfaces to user; consolidation does not proceed until resolved. |
| Orchestrator killed mid-review-queue | High | Low | `progress.json` resumability; reviewer state is idempotent. |
| Sign-convention regression slips through R04 review | Low | Critical | Reviewer re-derives sign from operator basis manually, does NOT trust the post-fix code; runs full pytest; grep checklist enforced. |
| Branch deletion order races (local before remote) | Low | Medium | Always remote first, then local; remote loop guarded so one failure doesn't abort. |
| Untracked scratch files silently committed | Was Medium (R1) | Now Low | Resolved by Phase 0 user gate. |
| Backup branch deleted before tag exists | Low | High | Resolved by Phase 2 push-verification gate. |
| pytest exceeds reviewer budget | Medium | Low | M-5 policy: PARTIAL not FAIL if >10 min. |
| GitHub issues never mirror local | Was High (M-2) | Now Low | Phase 9 step added; orchestrator commit policy explicit. |
