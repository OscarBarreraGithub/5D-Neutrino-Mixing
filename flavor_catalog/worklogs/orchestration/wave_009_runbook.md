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
| 1.a | PKA-CR001 | (pending) | `processes/collider_rs/CR001.{tex,yaml}` | PENDING |
| 1.b | PKA-CR002 | (pending) | `processes/collider_rs/CR002.{tex,yaml}` | PENDING |
| 1.c | PKA-CR003 | (pending) | `processes/collider_rs/CR003.{tex,yaml}` | PENDING |
| 1.d | PKA-CR004 | (pending) | `processes/collider_rs/CR004.{tex,yaml}` | PENDING |
| 1.e | PKA-CR005 | (pending) | `processes/collider_rs/CR005.{tex,yaml}` | PENDING |
| 1.f | PKA-CR006 | (pending) | `processes/collider_rs/CR006.{tex,yaml}` | PENDING |
| 1.g | PKA-CR007 | (pending) | `processes/collider_rs/CR007.{tex,yaml}` | PENDING |
| 1.h | PKA-CR008 | (pending) | `processes/collider_rs/CR008.{tex,yaml}` | PENDING |
| 1.i | PKA-CR009 | (pending) | `processes/collider_rs/CR009.{tex,yaml}` | PENDING |
| 1.j | PKA-CR010 | (pending) | `processes/collider_rs/CR010.{tex,yaml}` | PENDING |
| 1.k | PKA-CR011 | (pending) | `processes/collider_rs/CR011.{tex,yaml}` | PENDING |
| 1.l | PKA-CR012 | (pending) | `processes/collider_rs/CR012.{tex,yaml}` | PENDING |
| 1.m | PKA-CR013 | (pending) | `processes/collider_rs/CR013.{tex,yaml}` | PENDING |
| 1.n | PKA-CR014 | (pending) | `processes/collider_rs/CR014.{tex,yaml}` | PENDING |

### Stage 2-7: TBD (will update as we go)

---

## 5. Recovery instructions

Same protocol as Wave-8 runbook §5. Read this file + Wave-8 runbook +
AGENTIC_WORKFLOW.md, run `git status -sb` + `~/bin/codex_check_usage.sh`,
resume per the dispatch ledger.
