# Session Notes — Handoff to Future Orchestrator (or post-compaction self)

This file captures context-only knowledge from the orchestration session
that built the flavor-constraint catalog, that is **not** already written
into a more durable repo file. A new orchestrator (or this Claude after a
context compaction) should read this in addition to:

- `flavor_catalog/AGENTIC_WORKFLOW.md` (the reproducible pattern)
- `docs/phase_logs/flavor_catalog_plan_v1.md` (the plan)
- `docs/phase_logs/flavor_catalog_plan_v1_opus_signoff.md` (the plan
  approval and the recommended defaults on Section I)
- `docs/phase_logs/flavor_catalog_orchestrator_decisions.md` (the
  binding decisions adopted on the PI's behalf)
- The `signoff/round_001_index.md`, `round_002_index.md`,
  `round_003_index.md`, and `signoff/by_process/*.md` arbitration trail.

The structure below is **what to do if you have to pick this up tomorrow**.

---

## 1. Final state at handoff (updated 2026-05-17 after Wave-9)

| Artifact | State |
|---|---|
| Branch | `flavor-catalog/2026q2` cut from `paper/quark-scan-2026q2` at rc1.1 |
| Active processes | **102 total** — 94 PRIMARY (80 Waves 1–7 flavor + 14 Wave-9 collider_rs) + 8 SECONDARY (Wave-8), all OPUS-APPROVED |
| Deferred-scope | 43 plan-v1 Section C rows still DEFERRED-SCOPE (51 original − 8 promoted in Wave-8) |
| Out-of-scope (was) | Direct collider tails — **now IN-SCOPE via Wave-9** (custodial top partners, KK gauge resonances, VLQ, KK-graviton, EW-precision-tail) |
| Master PDF | `flavor_catalog/catalog_master.pdf`, **179 pages** (v0.4) |
| Tags | `flavor-catalog-v0` (75) → `v0.1` (fact-checked) → `v0.2` (80, +Wave-7) → `v0.3` (88, +Wave-8 SECONDARY) → **`v0.4` (102, +Wave-9 collider_rs PRIMARY)** |
| Companion paper | `quarkscan-paper-rc1.1` on `paper/quark-scan-2026q2` (frozen) |
| Priority policy | `flavor_catalog/PRIORITY_TIERS.md`. PRIMARY = `processes/<family>/`; SECONDARY = `processes/secondary/<family>/`. Wave-9 collider_rs is PRIMARY but in a new family — a new scope class opened by PI directive, not DA-deferred. |

Most recent Opus round-4 sign-off: 8/8 SECONDARY APPROVE (K019, K020, K021,
B007, B008, B013, B014, T014). See
`flavor_catalog/signoff/round_004_index.md`.

Fact-check status across all 88 processes:
**86 VERIFIED / 2 PARTIAL / 0 MISMATCH / 0 FAILED**.
- v0.2 PARTIAL: E009 (Weinberg three-gluon operator) — INSPIRE JS-only;
  APS journal + local snapshot confirm content.
- v0.3 PARTIAL: K020 (K+→π+ e μ) — NA62 2021 first-author convention
  drift; values/year/DOI/URL all verified; same class as E009; accepted
  per E009 precedent.

## 1a. Wave-8 summary (2026-05-17)

| Stage | Result |
|---|---|
| Stage 0 | `flavor_catalog/PRIORITY_TIERS.md` (109 lines) committed `24284d7`. |
| Stage 1 | 8 PKAs landed (commits `bab5bd0`, `2630168`, `ebd066c`, `e6e1cc3`, `5b68b43`, `6b64a18`, `e348a35`, `17599d5`). All YAMLs carry `priority_tier: SECONDARY`, `promoted_in_wave: 8`. |
| Stage 2 | `catalog_master.tex` wired with new Secondary section (commit `103833f`); 3 `processes/secondary/<family>/index.tex` stubs. |
| Stage 3 | 4 WA batches polished + WRITER-DONE (commits `f9c7ff4`, `88f2cde`, `2b23464`, `d198787`). |
| Stage 4 | 4 CA batches: 7/8 PASS cycle 1; K020 returned WRITER-REWORK on CHK-1 (Sher/E865 limit not promoted). WA-v2 fix (`5c14d2f`) + CA-v2 PASS (`a20d75a`). End state: 8/8 CHECKER-DONE. |
| Stage 5 | Fact-check addenda for 3 families: 7 VERIFIED / 1 PARTIAL (K020 metadata) / 0 MISMATCH / 0 FAILED (commits `118de1c`, `37514e9`, `7e121c3`). |
| Stage 6 | Opus round-4 sign-off: 8/8 APPROVE (commit `df7af95`). L001 / B001_B003 precedent uniformly applied. |
| Stage 7 | Master compile v0.3: 148 pages, 796,017 bytes PDF (commit `76c87ae`). Tag `flavor-catalog-v0.3` pushed to origin. |

Live wave state and dispatch ledger:
`flavor_catalog/worklogs/orchestration/wave_008_runbook.md`.

Wave-8 takes the catalog from 80 PRIMARY → 80 PRIMARY + 8 SECONDARY.
SECONDARY entries are repo-structurally separated under
`processes/secondary/<family>/` (PI request: easy to see + easy to
separate later). The `priority_tier` YAML field plus the dedicated master
PDF section give three independent ways to filter PRIMARY vs SECONDARY.

## 1b. Wave-9 summary (2026-05-17)

PI directive: "I think we should definitely have the custodial
constraints. What are the other constraints which are just as valid as
custodial? … lets get as many constraints as we can." Opened a new
scope class: direct-collider RS-resonance / heavy-partner searches.

| Stage | Result |
|---|---|
| Stage 0 | New family `processes/collider_rs/` + wave_009_runbook (commits `108cae4`, etc.). |
| Stage 1 | 14 PKAs landed (CR001–CR014) covering KK-gluon→tt, custodial T_{5/3}/T'/B, KK-EW dilepton + KK-W + KK-graviton spin-2, VLQ singlet + doublet, DY EFT, longitudinal VBS, diboson + diphoton resonance, 4-top. |
| Stage 2 | `catalog_master.tex` wired with new `\section{Collider RS Resonances}` (commit `b96036b`). |
| Stage 3 | 4 WA batches polished cycle-1 (commits `cd8a3fe`, `2286a39`, `c6d55b0`, `1eac13e`). |
| Stage 4 | 4 CA batches: 10/14 PASS cycle-1; 4 returned WRITER-REWORK (CR008/CR009 CHK-1 historical-numeral context; CR010 CHK-1 + CHK-2; CR011 CHK-2 citation-key mismatch). Two cycle-2 WA-v2 + CA-v2 batches (commits `a8758ac`+`950ca36`, `e1aec33`+`82daa9b`) cleared per T003 (remove historical numerals) + B023 (rename keys to manifest) precedents. End state: 14/14 CHECKER-DONE. |
| Stage 5 | Fact-check: **14/14 VERIFIED, 0 mismatches, 0 fetch exceptions** (commit `0c5dacc`, `audits/factcheck_collider_rs.md`). |
| Stage 6 | Opus round-5 sign-off: **14/14 APPROVE** (commit `cdd0238`, `signoff/round_005_index.md`). |
| Stage 7 | Master compile v0.4: 179 pages, 913,614 bytes PDF (commit `2ad34b1`). Tag `flavor-catalog-v0.4` pushed. |

Live wave state and decisions D-1..D-8:
`flavor_catalog/worklogs/orchestration/wave_009_runbook.md`.

Wave-9 sub-classes added:
- **Custodial fermion partners** (CR002 T_{5/3}, CR003 T_{2/3}=T', CR004 B): iconic RS-distinguishing custodial signatures.
- **KK gauge resonances** (CR001 KK-gluon→tt, CR005 KK-Z/γ dilepton, CR006 KK-W, CR007 KK-graviton spin-2): canonical LHC RS searches.
- **VLQ + 4-top** (CR008 singlet T, CR010 doublet (T,B), CR014 4-top): non-custodial collider partners.
- **EW-precision-tail** (CR009 DY high-mass EFT, CR011 W_LW_L scattering, CR012 diboson resonance, CR013 diphoton resonance): model-independent collider constraints with RS interpretation.

Cumulative fact-check across all 102 entries:
**100 VERIFIED / 2 PARTIAL / 0 MISMATCH / 0 FAILED.**
The two PARTIALs are both metadata-only (E009 INSPIRE JS-only at v0.2,
K020 NA62 author convention — cleared at v0.3 cleanup; v0.3 tag still
pins the original PARTIAL state for historical accuracy).

---

## 2. Branches and sync recommendation

Three branches matter:

- `main` — at session start showed many modified files (audit work,
  methodology PDF, scan code, tests). Those changes are very likely
  already on `paper/quark-scan-2026q2`, which was cut from main with that
  work brought along, but this should be verified explicitly with
  `git diff main..paper/quark-scan-2026q2` before any merge.

- `paper/quark-scan-2026q2` — the rc1.1 paper (sealed). All
  Phase 2 / Phase 3 audits, BMU sign fix, methodology note, Opus reviews
  landed here. **80+ commits ahead of main when cut from it.**

- `flavor-catalog/2026q2` — the catalog work. Cut from
  `paper/quark-scan-2026q2` at rc1.1. **Disjoint from existing repo
  content** — adds only the `flavor_catalog/` subdirectory and a few
  audit / phase-log files under `docs/`.

The PI's recommended sync sequence (the user explicitly affirmed this):

1. Verify what's on `main` that isn't on `paper/quark-scan-2026q2`.
   Most likely answer: nothing essential; main's local modifications
   from session start were the audit work that's now on the paper branch.
2. Merge `paper/quark-scan-2026q2` into `main` with `--no-ff` (history
   is the artifact — preserve the sign-off chain).
3. Run the full pytest suite on `main` post-merge. Expected: **543
   passed / 1 skipped / 0 xfails** (matches rc1.1 baseline).
4. Merge `flavor-catalog/2026q2` into `main`. Options: `--no-ff` to
   preserve worklog history, or `--squash` since `flavor_catalog/` is
   self-contained. The user did not pin a choice; both are defensible.
   Recommendation: `--no-ff` for symmetry with the paper merge.
5. Tag the merge point on `main`.
6. Keep feature branches alive for follow-up work.

The PI plans to do this manually ("I'll do it"). Do not auto-merge.

---

## 3. Open questions for the PI (worth surfacing on next touch)

These are flagged in `docs/phase_logs/flavor_catalog_plan_v1.md` Section I
and in the DA-1/DA-4 worklogs but were resolved by orchestrator on PI's
behalf via `docs/phase_logs/flavor_catalog_orchestrator_decisions.md`. The
PI may want to reaffirm or revisit:

1. **`K1 → πγ` interpretation**: orchestrator chose `K_L → π⁰γγ` (K013).
   Alternative `K_S → π⁰γγ` (K014) was deferred.
2. **Scope**: full default (quark + LFV + neutrino + EW + Higgs/top FCNC +
   EDM). PI can prune by marking specific entries DEFERRED-SCOPE.
3. **Branch policy**: `flavor-catalog/2026q2` (separate from paper). PI
   can re-merge or rebase later.
4. **License**: minimal text snapshots, no publisher PDFs. Cheap to
   upgrade; expensive to retract.
5. **rc1.1 / rc2 relationship**: companion artifact. PI may eventually
   integrate select catalog entries into the next paper as a "future
   constraints" section.
6. **Code integration target**: `quarkConstraints/modern/` for any
   future graduation from catalog to live constraint.
7. **B012 promotion** (decided by Wave-7 follow-up): B012 (B→K*γ) was
   promoted from DA-1 deferred to active because external review showed
   exclusive radiative helicity observables are RS-distinguishing.

---

## 4. Key operational quirks of this cluster

### 4.1 The codex_worker.sh wrapper

Path: `/n/home09/obarrera/bin/codex_worker.sh`.

**Critical**:
- Always invoke with `< /dev/null` to close stdin. Without this the codex
  CLI hangs silently for hours on "Reading additional input from
  stdin...". The wrapper itself does not redirect stdin.
- Do **not** use the wrapper's `-o` flag. It overwrites the output file
  with codex's summary message, clobbering any longer content. Instead,
  instruct each agent to write its deliverable file via shell heredoc.
- The wrapper's `--dangerously-bypass-approvals-and-sandbox` flag is
  mandatory for programmatic use; codex otherwise hangs on human
  approval prompts.

### 4.2 The codex usage probe

Path: `/n/home09/obarrera/bin/codex_check_usage.sh`. Built during this
session after a 4 PM (local) quota exhaustion incident. Run **before
every dispatch wave**. Returns:
- exit 0 if usable
- exit 1 if rate-limited (with stderr including the reset-time message)

### 4.3 Codex CLI is logged in via ChatGPT (not API key)

`codex login status` shows "Logged in using ChatGPT". This means usage
is subject to the ChatGPT subscription quota. There is no `codex
status` command exposed to query remaining quota; the only way to
check is `codex exec --dangerously-bypass-approvals-and-sandbox
"/status" < /dev/null` and parse stderr (which is what the probe does).

### 4.4 Model: codex GPT-5.5 xhigh (not fast)

Configured in `~/.codex/config.toml`. Do not invoke `codex --fast` —
fast mode is a lighter model and is not appropriate for the work in this
project.

### 4.5 Opus agents are dispatched via the Agent tool

Not the codex wrapper. Specifically `Agent` with
`subagent_type: general-purpose` and `model: opus`. Opus uses the
Anthropic API and is **not affected** by codex quota.

### 4.6 SLURM partition note (unrelated to catalog; for the quark scan)

The cluster's `randall_lab` partition is not available; use
`serial_requeue --requeue` for scan-output dispatch.

---

## 5. Things a new orchestrator should NOT do

- Do not push to `paper/quark-scan-2026q2` — it is **frozen at rc1.1**.
  All catalog work goes on `flavor-catalog/2026q2`.
- Do not modify any per-process `.tex` content yourself. Always delegate
  to a WA. The orchestrator only writes orchestration meta-docs.
- Do not invent agent dispatches with model names other than codex
  GPT-5.5 xhigh (for content) or Opus (for sign-off / arbitration).
  Mixing models breaks reproducibility.
- Do not override the writer/checker separation rule. Two real BLOCKER
  physics errors in the companion paper were caught only because the
  CA was genuinely independent from the WA.
- Do not skip the cycle-3 cap. If a unit fails WA/CA three times,
  escalate to Opus arbitration immediately. The L001 / B021 / B023 /
  B001 / B003 precedents are now codified.
- Do not auto-merge any feature branch into `main` without the PI's
  explicit confirmation.

---

## 6. Specific things a future Claude should know

### 6.1 The "factcheck status PARTIAL on E009" non-issue

If the fact-check status table shows E009 as PARTIAL, it is **not a
real content gap**. The INSPIRE-HEP URL renders via JavaScript so
WebFetch cannot read it; the APS journal page (Weinberg 1989) and the
local snapshot independently confirm the content. Do not "fix" it; it
is documented in `flavor_catalog/audits/factcheck_status.md`.

### 6.2 The T020 ATLAS digit-slip story

`T020` (h → eμ) had a `6.1 × 10⁻⁵` quote that should have been
`6.2 × 10⁻⁵`. Caught by the fact-check WebFetch round. Fixed in commit
`6498fad`. Mentioned here so a re-fact-check doesn't get confused
seeing the corrected value.

### 6.3 Wave naming and provenance

- **Wave-1**: 10 PI-seed processes (K003-006, K013, B009/B011/B015,
  T001/T010)
- **Wave-2**: 8 processes (K001/K002, C001, B002/B005, T002, L001, E001)
- **Wave-3**: 8 processes (B032/B033, B018, B025, L007, L002, T007, T018)
- **Wave-4a/4b**: 20 from DA-1 (8a + 12b)
- **Wave-5a/5b**: 20 from DA-2 (12a + 8b)
- **Wave-6**: 8 from DA-3 (B001/B003/B016, K012/K018, L023/T020, E009)
- **Wave-7**: 5 from external Deep Research review (T003/T004/T008/T012,
  B012)
- DA-4 converged (no Wave-7 of its own; the Wave-7 PKAs are external-review
  follow-ups, not DA-4 proposals)

### 6.4 The 5 individual Opus arbitrations (all APPROVE-OVERRIDE under
L001 precedent)

| Unit | Issue | Doc |
|---|---|---|
| L001 | "3 TeV reference scale" not in pdg_or_equivalent | `signoff/by_process/L001.md` |
| B021, B023 | Dataset metadata + citation key typo | `signoff/by_process/B021_B023.md` |
| B001, B003 | FLAG bag parameters + lattice context not in pdg_or_equivalent | `signoff/by_process/B001_B003.md` |

These set the policy precedent that downstream CAs and Opus rounds
applied uniformly.

### 6.5 The `quarkscan-paper-rc1.1` rebuild context

If you ever need to rebuild the paper PDF, the headline values to expect:
- M_KK^min p50 = **47.26 +69.4 −25.0 TeV** at g_s* = 3 (BGS 2020 + LO BMU
  + factor-3 PDG gate)
- M_KK^min p95 = **127.13 TeV**
- CFW reconciliation: factor 2.2 stronger at matched conventions
- ε_K binding fraction: 99.5%
- pytest: **543 passed / 1 skipped / 0 xfails / 0 failures**
- Methodology note: **19 pages**
- PDF sha256: `5f544e5d1654f52add06fd8adfe97fed77b968e5e1bea79e119882b1ac898883`

### 6.6 The rc1 post-compaction briefing is OBSOLETE

`docs/phase_logs/POST_COMPACTION_BRIEFING.md` covers the rc1 → rc1.1
transition. It is **no longer the right entry point** because rc1.1 is
sealed, Task A and B are done, and the catalog work is the active surface.
Read it for historical context only.

---

## 7. What I would do next if I came back tomorrow

In priority order:

1. **Have the PI confirm sync of `flavor-catalog/2026q2` and
   `paper/quark-scan-2026q2` into `main`.** Once merged, all the audit
   trail history lives on main and external collaborators can clone and
   review.

2. **Optionally**: write a short "Catalog overview" companion document
   (`flavor_catalog/CATALOG_OVERVIEW.md`) describing what the catalog
   contains in physicist-readable prose. The current
   `catalog_master.pdf` (133 pages) is too long for a quick read; a
   1–2 page summary with the top-tier observables grouped by RS
   sensitivity would be valuable for collaborators.

3. **Optionally**: dispatch a Wave-8. The full deferred-scope list with
   rationale is in
   `flavor_catalog/worklogs/discovery/round_004_addendum_deferred_scope.md`
   (51 plan-v1 Section C rows). Orchestrator-side leverage tiering for
   selecting Wave-8 candidates:

   **Top tier (high-leverage; ~8–10 candidates)** — these are the
   defensible adds if Wave-8 happens at all:
   - **K019 / K020 / K021** — LFV kaon trio (`K_L→μe`, `K⁺→π⁺μe`,
     `K_L→π⁰eμ`). Promotion trigger: scope extends to lepton-bulk RS.
   - **B007 / B008** — `B_{s,d} → e⁺e⁻, τ⁺τ⁻` rare leptonic tails;
     companions to drafted B005 / B006.
   - **B013 / B014** — exclusive `b → sγ` (`B_s → φγ`) and `b → dγ`
     (`B → ργ, ωγ`); helicity companions to B011 / B012.
   - **T014** — flavor-changing Z decays (`Z → bs, Z → bd, Z → sd`);
     distinct from drafted Z-pole / Z-LFV entries.

   **Middle tier (defensible; ~5–8 candidates)** — Wave-9 material
   unless explicitly asked:
   - B020 (`B_s → φℓ⁺ℓ⁻`), B024 (`B → πνν̄`), B027 (`R_{J/ψ}`),
     B028 (`B_c → τν`).
   - L011–L013 hadronic tau LFV anchors.
   - C009 (`D⁰ → γγ`), C010 (`D → ργ / φγ`).
   - E003 tau EDM, E005 proton / deuteron EDM prospects.

   **Tail tier (~30+ rows)** — DA-4 judged "diminishing returns";
   long-distance dominated, redundant with existing entries, or
   future-projection-only. Do not promote unless the PI specifically
   wants literal plan-v1 Section C completion (which DA-4 explicitly
   recommended against).

   Sensible Wave-8 size: **8–10 PKAs** (top tier). Pattern unchanged
   from Wave-7: PKAs in parallel → WA batches → CA → fact-check via
   WebFetch → Opus round-4 sign-off → master rebuild → tag v0.3.
   Half-day of orchestration.

4. **Optionally**: build the first MEDIUM-difficulty constraint into
   the code (e.g. b → sγ via `quarkConstraints/modern/`) as a
   proof-of-concept that the catalog can graduate to live constraints.
   This is the "top-10 horizon" mentioned in the orchestration chat
   summary: ~6–9 months with 1–2 developers.

5. **Do not** dispatch new agents without first running
   `~/bin/codex_check_usage.sh`.

---

## 8. Specific files the PI may want to read first

In rough order of importance for a quick review:

1. `flavor_catalog/AGENTIC_WORKFLOW.md` — the reproducible playbook.
2. `flavor_catalog/audits/factcheck_status.md` — the 80-row checklist.
3. `flavor_catalog/signoff/round_001_index.md`,
   `round_002_index.md`, `round_003_index.md` — Opus's per-round verdicts.
4. `flavor_catalog/signoff/by_process/L001.md`, `B021_B023.md`,
   `B001_B003.md` — the policy precedents.
5. `flavor_catalog/external_research/deepresearch_may{15,16}_review.md`
   — the external comparison reviews.
6. `flavor_catalog/worklogs/discovery/round_001_full_scope.md`,
   `round_002_followup.md`, `round_003_final_sweep.md`,
   `round_004_convergence.md`, `round_004_addendum_deferred_scope.md`
   — the discovery trail.
7. `docs/phase_logs/flavor_catalog_plan_v1.md` — the original plan.
8. `docs/phase_logs/flavor_catalog_orchestrator_decisions.md` — the
   orchestrator's adjudications on the PI's behalf.

Each is self-contained; nothing in this `SESSION_NOTES.md` is uniquely
critical for understanding the catalog's content — it's only critical
for understanding the **orchestration pattern** that produced it.

---

End of session notes. Good luck.
