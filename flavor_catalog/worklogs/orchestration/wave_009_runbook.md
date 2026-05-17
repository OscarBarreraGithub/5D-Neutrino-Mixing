# Wave-9 Orchestration Runbook (live)

**Status**: ACTIVE
**Started**: 2026-05-17 (after Wave-8 close-out and v0.3 tag)
**Orchestrator**: Claude Opus 4.7 (1M ctx), continuing from Wave-8.
**Goal**: Add 14 collider RS resonance / heavy-partner constraints to the
catalog (custodial top/bottom partners, KK gauge resonances, VLQs, EW-tail
EFT bounds). Tag `flavor-catalog-v0.4`.

This file is the canonical state file for Wave-9. Updated as agents land;
survives orchestrator replacement / context compaction (read this + the
Wave-8 runbook + `flavor_catalog/AGENTIC_WORKFLOW.md` to resume).

---

## 1. PI directive (2026-05-17)

> "I think we should definitely have the custodial constraints… now that
> this done, we should move on to do the custodial. What are the other
> constraints which are just as valid as custodial?" → after orchestrator
> proposed 14 same-tier entries: "I think that is perfect. Can you go
> ahead and act as my orchestrator with the entire review pipeline we
> have set up? lets get as many constraints as we can."

This opens a **new scope class**: direct-collider RS-resonance / heavy-
partner searches. These are complementary to (not replacements for) the
low-energy flavor entries in Waves 1-7. LHC reach ($M_{KK}\sim\text{TeV}$)
is weaker than the flavor-implied bound ($M_{KK}\sim 47$ TeV at $g_{*}=3$),
but the entries are RS-distinguishing (especially custodial top partners
and spin-2 KK-graviton) and provide cross-checks under non-anarchic RS
scenarios.

---

## 2. Wave-9 scope (decided)

14 PRIMARY entries, new family `processes/collider_rs/`, IDs CR001–CR014.

| ID | Process | Sub-class | WA batch |
|----|---------|-----------|----------|
| CR001 | KK-gluon ($g^{(1)}_{KK}$) → $t\bar{t}$ resonance | KK gauge | B |
| CR002 | $T_{5/3}$ ($X_{5/3}$) pair production → same-sign dileptons + jets | Custodial fermion partner | A |
| CR003 | $T_{2/3}$ ($T'$) → $tZ$, $tH$, $bW$ | Custodial fermion partner | A |
| CR004 | Custodial $B$ → $tW$, $bZ$, $bH$ | Custodial fermion partner | A |
| CR005 | KK-EW ($\gamma^{(1)}, Z^{(1)}$) → $\ell^+\ell^-$ high-mass tail | KK gauge | B |
| CR006 | KK-W ($W^{(1)}$) → $\ell\nu$, $tb$ | KK gauge | B |
| CR007 | KK-graviton ($G^{(1)}$) → diboson / diphoton / dilepton (spin-2) | KK gauge | B |
| CR008 | VLQ singlet $T$ → $tZ$, $tH$, $bW$ | VLQ | C |
| CR009 | Drell-Yan high-mass tail (EFT contact-op / $Z'$-like) | EW-precision-tail | D |
| CR010 | VLQ doublet $(T,B)$ → various | VLQ | C |
| CR011 | $W_L W_L \to W_L W_L$ longitudinal vector-boson scattering | EW-precision-tail | D |
| CR012 | Diboson high-mass resonance (WW/WZ/ZZ generic spin-1) | EW-precision-tail | D |
| CR013 | Diphoton high-mass resonance (spin-0/2) | EW-precision-tail | D |
| CR014 | 4-top quark ($t\bar{t}t\bar{t}$) production | VLQ adjacent | C |

---

## 3. Orchestrator decisions (binding for this wave)

### D-1: Tier = PRIMARY (implicit)
**Decision.** No `priority_tier` field in YAML sidecars. Implicit-PRIMARY
default per `flavor_catalog/PRIORITY_TIERS.md`.

**Why.** These are a new scope class opened by explicit PI directive, not
DA-deferred. SECONDARY tag is reserved for DA-deferred re-promotions.

### D-2: Family = `collider_rs`
**Decision.** New family directory `flavor_catalog/processes/collider_rs/`.
References, worklogs, signoff stay flat (per-process-ID).

**Why.** "Collider RS" is precise: these are LHC searches for RS-specific
resonances and heavy partners. Distinct from low-energy `kaon/`, `beauty/`,
etc. and from `top_higgs_ew/` (which holds top-FCNC decays, not direct
resonance searches).

### D-3: catalog_master.tex placement
**Decision.** New top-level `\section{Collider RS Resonances}` placed
after the EDM/Neutrino family section and before the existing
`\section{Secondary Entries (Wave-8+, …)}`.

**Why.** Keeps PRIMARY entries in the main flow; SECONDARY stays the
trailing section.

### D-4: WA batching (4 batches)
- **WA-A custodial-partners**: CR002, CR003, CR004 (3 entries)
- **WA-B KK-gauge**: CR001, CR005, CR006, CR007 (4 entries)
- **WA-C VLQ-and-4-top**: CR008, CR010, CR014 (3 entries)
- **WA-D EW-tail**: CR009, CR011, CR012, CR013 (4 entries)

### D-5: Fact-check = single family agent
**Decision.** One `factcheck-codex-collider_rs` agent covering all 14
entries. Output: new `flavor_catalog/audits/factcheck_collider_rs.md`.

**Why.** New family = new fact-check audit file (parallel to existing
`factcheck_kaon.md`, etc.).

### D-6: CHK-1 carve-out applies uniformly
**Decision.** Apply L001 / B001_B003 precedent: measured cross-section /
mass-exclusion limits go in `pdg_or_equivalent.values`; theoretical
cross-section predictions, integrated-luminosity dataset metadata, and
EFT-coefficient translations stay in `paper_era_reference` /
`supporting_measurements` / `auxiliary_theory_inputs`.

### D-7: Tag = `flavor-catalog-v0.4`
**Decision.** On successful Opus round-5 + master compile, annotated
tag `flavor-catalog-v0.4` at the compile commit.

### D-8: Pipeline = standard
**Decision.** PKA → WA → CA → fact-check → Opus, identical to Waves 7-8.
Cycle caps unchanged. L001 arbitration precedent available if needed.

---

## 4. Dispatch ledger

### Stage 0: scaffold (Claude inline)

| Stage | Output | Status |
|-------|--------|--------|
| 0.a | `flavor_catalog/processes/collider_rs/` dir | DONE |
| 0.b | This runbook | DONE |

### Stage 1: PKA (14 in parallel)

| Stage | Agent ID | Background ID | Output path | Status |
|-------|----------|---------------|-------------|--------|
All 14 dispatched 2026-05-17 16:23 EDT. Prompts in
`/tmp/wave9_prompts/pka_<ID>.prompt` (224 lines each, ~10KB, all clean).
Logs in `/tmp/wave9_logs/pka_<ID>.log`.

| Stage | Agent ID | Background ID | Output path | Status |
|-------|----------|---------------|-------------|--------|
| 1.a | PKA-CR001 (KK-gluon→tt) | `bp0jq3v1d` | `processes/collider_rs/CR001.{tex,yaml}` | DISPATCHED |
| 1.b | PKA-CR002 (T_{5/3}) | `b1ritjrsi` | `processes/collider_rs/CR002.{tex,yaml}` | DISPATCHED |
| 1.c | PKA-CR003 (T') | `bv9nlgn4u` | `processes/collider_rs/CR003.{tex,yaml}` | DISPATCHED |
| 1.d | PKA-CR004 (custodial B) | `brlh8bfrs` | `processes/collider_rs/CR004.{tex,yaml}` | DISPATCHED |
| 1.e | PKA-CR005 (KK-Z/γ dilepton) | `buociv7yd` | `processes/collider_rs/CR005.{tex,yaml}` | DISPATCHED |
| 1.f | PKA-CR006 (KK-W) | `byr3qw7c9` | `processes/collider_rs/CR006.{tex,yaml}` | DISPATCHED |
| 1.g | PKA-CR007 (KK-graviton) | `bd4ci8dli` | `processes/collider_rs/CR007.{tex,yaml}` | DISPATCHED |
| 1.h | PKA-CR008 (VLQ singlet) | `b6j346syj` | `processes/collider_rs/CR008.{tex,yaml}` | DISPATCHED |
| 1.i | PKA-CR009 (DY EFT) | `bs7u1j9kj` | `processes/collider_rs/CR009.{tex,yaml}` | DISPATCHED |
| 1.j | PKA-CR010 (VLQ doublet) | `ba4svf6b6` | `processes/collider_rs/CR010.{tex,yaml}` | DISPATCHED |
| 1.k | PKA-CR011 (W_L W_L scattering) | `baw22rett` | `processes/collider_rs/CR011.{tex,yaml}` | DISPATCHED |
| 1.l | PKA-CR012 (Diboson resonance) | `br6g1g507` | `processes/collider_rs/CR012.{tex,yaml}` | DISPATCHED |
| 1.m | PKA-CR013 (Diphoton resonance) | `b2s3e8ll0` | `processes/collider_rs/CR013.{tex,yaml}` | DISPATCHED |
| 1.n | PKA-CR014 (4-top production) | `b8xk0ehvo` | `processes/collider_rs/CR014.{tex,yaml}` | DISPATCHED |

### Stage 2-7: ALL DONE (see commit chain below)

| Stage | Result |
|---|---|
| 2.a | scaffold2-w9 — `catalog_master.tex` + `collider_rs/index.tex` wired (commit `b96036b`) |
| 3.a..d | 4 WA batches → 14/14 WRITER-DONE (commits `cd8a3fe`, `2286a39`, `c6d55b0`, `1eac13e`) |
| 4.cycle-1 | 4 CA batches: 10 PASS (CR001-CR007, CR012-CR014); 4 WRITER-REWORK (CR008 CHK-1, CR009 CHK-1, CR010 CHK-1+CHK-2, CR011 CHK-2) — commits `e8daa00`, `7a6fa3d`, `bc35a06`, `e846f14` |
| 4.cycle-2 | WA-v2 vlq_4top (CR008+CR010) commit `a8758ac`; WA-v2 ew_tail (CR009+CR011) commit `e1aec33`. CA-v2 vlq_4top commit `950ca36`; CA-v2 ew_tail commit `82daa9b`. T003 + B023 precedents applied — 14/14 CHECKER-DONE. |
| 5.a | Fact-check collider_rs: 14/14 VERIFIED, 0 mismatches, 0 fetch exceptions (commit `0c5dacc`) |
| 6.a | Opus round-5 (via Agent tool, `model: opus`): 14/14 APPROVE (commit `cdd0238`) |
| 7.a | Master compile v0.4: 179 pages, 913 KB PDF (commit `2ad34b1`) |
| 7.b | Tag `flavor-catalog-v0.4` annotated + pushed to origin |

---

## 6. Wave-9 close-out

**Wave status: COMPLETE.** 14 new PRIMARY collider_rs entries drafted,
polished, verified, fact-checked, signed off, compiled, tagged. New
family `collider_rs/` is a new scope class opened by PI directive,
complementary to (not replacing) the low-energy flavor entries.

Catalog state at v0.4:
- 94 PRIMARY (80 flavor Waves 1-7 + 14 collider_rs Wave-9) + 8 SECONDARY (Wave-8) = 102 total OPUS-APPROVED
- 100 VERIFIED + 2 PARTIAL (E009 v0.2 + K020 v0.3 — both metadata-only, both accepted per E009 precedent; K020 was actually cleared in the v0.3 cleanup commits but v0.3 tag pins the original PARTIAL state)
- Master PDF: 179 pages
- Tag: `flavor-catalog-v0.4` on commit `2ad34b1`

This runbook is now closed. A future Wave-10 should start a new runbook
at `flavor_catalog/worklogs/orchestration/wave_010_runbook.md` using
this one + the Wave-8 runbook as templates.

Useful Wave-9 takeaways for future orchestrators:
- 4-of-14 cycle-1 CA rework rate (~29%) is in line with Wave-8 (1-of-8 = ~13%) and Wave-1..6 historical rates. Collider entries have a slightly higher base rate of CHK-1 historical-numeral and CHK-2 citation-key issues than low-energy entries; both are cycle-2-fixable via T003 + B023 precedents.
- The new family directory + new master TeX section pattern (Wave-9) is identical to the Wave-8 `processes/secondary/` pattern in structural mechanics, differing only in tier label and rationale. A future Wave-10 opening another scope class would use the same recipe.
- Fact-check pass rate at 14/14 VERIFIED, 0 fetch exceptions is the cleanest of any wave to date. Likely because collider entries cite ATLAS / CMS arXiv abstracts and PDG live pages with WebFetch-friendly URLs (no JS-only INSPIRE renders or paywalled-only journals).

---

## 5. Recovery instructions

Same protocol as Wave-8 runbook §5. Read this file + Wave-8 runbook +
AGENTIC_WORKFLOW.md, run `git status -sb` + `~/bin/codex_check_usage.sh`,
resume per the dispatch ledger.
