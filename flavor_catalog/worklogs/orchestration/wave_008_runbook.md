# Wave-8 Orchestration Runbook (live)

**Status**: ACTIVE
**Started**: 2026-05-17
**Orchestrator**: Claude Opus 4.7 (1M ctx), succeeding the v0.2 orchestrator
**Companion handoff doc**: `HANDOFF_PROMPT.md` (committed 3be3a7b)
**Goal**: Add the 8 top-tier Wave-8 candidates from `SESSION_NOTES.md §7` to the
catalog as **SECONDARY** entries, run the full PKA→WA→CA→fact-check→Opus
pipeline, tag `flavor-catalog-v0.3`.

This file is the canonical state file for the wave. It is updated as agents
land. If the conversation is compacted, a new orchestrator can resume by
reading this file plus the standard handoff stack.

---

## 1. PI directives binding this wave

From 2026-05-17 PI message:

1. **Keep all Wave-8 additions as SECONDARY.** Non-negotiable.
2. **Make the SECONDARY distinction visible from the repo structure** — not
   just a YAML tag — so a future read can cleanly separate.
3. **Act as orchestrator until done; do not stop and ask questions; decide on
   PI's behalf and document each decision.**
4. **Notify PI when all 8 are in and reviewed + verified.**

---

## 2. Wave-8 scope (decided)

8 PKAs, drawn from `SESSION_NOTES.md §7` top tier:

| ID | Process | Family | WA batch |
|----|---------|--------|----------|
| K019 | `K_L → e^± μ^∓` | kaon | WA1 (kaon-LFV) |
| K020 | `K^+ → π^+ e^± μ^∓` | kaon | WA1 (kaon-LFV) |
| K021 | `K_L → π^0 e^± μ^∓` | kaon | WA1 (kaon-LFV) |
| B007 | `B_{s,d} → e^+ e^-` | beauty | WA2 (B-rare-leptonic) |
| B008 | `B_{s,d} → τ^+ τ^-` | beauty | WA2 (B-rare-leptonic) |
| B013 | `B_s → φ γ` | beauty | WA3 (B-radiative) |
| B014 | `B → ρ γ`, `B → ω γ` | beauty | WA3 (B-radiative) |
| T014 | `Z → bs`, `Z → bd`, `Z → sd` | top_higgs_ew | WA4 (T014 solo) |

All 8 are deferred-by-DA-4 rows from `round_004_addendum_deferred_scope.md`.

---

## 3. Orchestrator decisions on PI's behalf (binding for this wave)

### D-1: Repo-level SECONDARY separation
**Decision.** New entries land under `flavor_catalog/processes/secondary/<family>/<id>.tex` and `.yaml`, in addition to the YAML `priority_tier: SECONDARY` field. The existing 80 PRIMARY entries stay where they are (`processes/<family>/<id>.tex`).

**Why.** PI explicitly asked for repo-structure visibility. Subdirectory split satisfies "easy to see" (just `ls processes/secondary/`) and "easy to separate later" (one `git mv processes/secondary/*/<id>.* processes/<family>/` would promote any entry; one tree-deletion would purge the whole tier).

**How to apply.** Every Wave-8+ PKA writes to `processes/secondary/<family>/<id>.{tex,yaml}`. References (`references/<id>/`), worklogs (`worklogs/{pka,writer,checker}/<id>.md`), and signoff (`signoff/by_process/<id>.md`) stay flat — process IDs are globally unique so nesting them adds nothing.

### D-2: catalog_master.tex layout for SECONDARY
**Decision.** A new top-level `\section{Secondary Entries (Wave-8+, lower implementation priority)}` is appended to `catalog_master.tex`, with `\subsection*{<Family>}` and `\input{processes/secondary/<family>/index}` per affected family.

**Why.** Keeps the master PDF a single artifact (so reviewers see everything) while making the PRIMARY/SECONDARY boundary visually unmissable in both the source tree and the rendered PDF.

**How to apply.** Wave-8 affects 3 families: kaon, beauty, top_higgs_ew. A scaffold agent updates `catalog_master.tex` and creates 3 `processes/secondary/<family>/index.tex` stubs after PKAs land.

### D-3: YAML sidecar contract for SECONDARY
**Decision.** Every Wave-8 YAML sidecar carries:
```yaml
priority_tier: SECONDARY
priority_rationale: <one-line: why deferred + why now>
promoted_in_wave: 8
```
The 80 PRIMARY entries are NOT retroactively tagged this round (handoff prompt: optional; default-implicit is PRIMARY). Retroactive tag pass is a clean follow-up if PI wants.

### D-4: WA batching
4 WA batches matching the table above. Constraint: a CA never checks the batch it wrote, per AGENTIC_WORKFLOW §2 (writer/checker separation).

### D-5: Cycle-cap policy
Unchanged from plan v1 §G: max 3 WA/CA cycles per batch. On cap, escalate to Opus arbitration per L001 / B021_B023 / B001_B003 precedents (CHK-1 carve-out for theory normalization, dataset metadata, lattice averages).

### D-6: Fact-check grouping
One fact-check agent per family: `factcheck-codex-kaon-w8` (K019/K020/K021), `factcheck-codex-beauty-w8` (B007/B008/B013/B014), `factcheck-codex-tophiggsew-w8` (T014).

### D-7: Tag
On successful Opus round-4 sign-off + master compile, tag `flavor-catalog-v0.3` and push. Do NOT push to `paper/quark-scan-2026q2` (frozen).

### D-8: PRIORITY_TIERS.md ownership
Dispatched to one codex agent (not Claude), per handoff prompt instruction. Codex writes the canonical policy doc; this runbook records the dispatch metadata.

---

## 4. Dispatch ledger

Each row records one agent dispatch. Updated by the orchestrator as agents
land. Format:

| Stage | Agent ID | Started | Background ID | Output path | Status |
|-------|----------|---------|---------------|-------------|--------|

### Stage 0: scaffold

| Stage | Agent ID | Started | Background ID | Output path | Status |
|-------|----------|---------|---------------|-------------|--------|
| 0.a | priority-tiers-bootstrap | 2026-05-17 02:04 EDT | `btehv33ko` (log `/tmp/wave8_logs/priority_tiers.log`) | `flavor_catalog/PRIORITY_TIERS.md` | DONE (commit `24284d7`, 109 lines) |

### Stage 1: PKA (8 in parallel)

Prompts in `/tmp/wave8_prompts/pka_<ID>.prompt` (189 lines each, zero
unsubstituted tokens). Logs in `/tmp/wave8_logs/pka_<ID>.log`. All
dispatched 2026-05-17 02:04 EDT.

| Stage | Agent ID | Started | Background ID | Output path | Status |
|-------|----------|---------|---------------|-------------|--------|
| 1.a | PKA-K019 | 02:04 EDT | `bi755ba90` | `processes/secondary/kaon/K019.{tex,yaml}` | DONE (origin commit `e348a35`; tex 108L / yaml 169L / 12 refs) |
| 1.b | PKA-K020 | 02:04 EDT | `bzu9mclku` | `processes/secondary/kaon/K020.{tex,yaml}` | DONE (origin commit `bab5bd0`; tex 113L / yaml 156L / 8 refs) |
| 1.c | PKA-K021 | 02:04 EDT | `bjh9hox5j` | `processes/secondary/kaon/K021.{tex,yaml}` | DONE (origin commit `5b68b43`; tex 114L / yaml 180L / 11 refs) |
| 1.d | PKA-B007 | 02:04 EDT | `b5fffrzza` | `processes/secondary/beauty/B007.{tex,yaml}` | DONE (origin commit `e6e1cc3`; tex 110L / yaml 220L / 9 refs) |
| 1.e | PKA-B008 | 02:04 EDT | `bgmmjm93i` | `processes/secondary/beauty/B008.{tex,yaml}` | DONE (origin commit `2630168`; B_s 6.8e-3 / B_d 2.1e-3; tex 115L / yaml 200L / 11 refs) |
| 1.f | PKA-B013 | 02:04 EDT | `bjqgpgtf9` | `processes/secondary/beauty/B013.{tex,yaml}` | DONE (origin commit `ebd066c`; tex 113L / yaml 253L / 10 refs) |
| 1.g | PKA-B014 | 02:04 EDT | `bomkm2ot3` | `processes/secondary/beauty/B014.{tex,yaml}` | DONE (origin commit `17599d5`; tex 139L / yaml 383L / 14 refs) |
| 1.h | PKA-T014 | 02:04 EDT | `bs3v54b0w` | `processes/secondary/top_higgs_ew/T014.{tex,yaml}` | DONE (origin commit `6b64a18`; tex 111L / yaml 218L / 9 refs) |

**Stage-1 audit**: all 8 SECONDARY entries verified on origin. Every YAML
carries `priority_tier: SECONDARY` and `promoted_in_wave: 8`. Cosmetic
artefact: commit `37beabb` carries the K020 commit message but B013
content (concurrent-commit race; content-correct, just mislabelled).
Harmless — content is anchored by both `bab5bd0` (real K020) and `ebd066c`
(real B013); the noisy commit is a no-op duplicate of B013 with a wrong
message. Not worth a history rewrite.

### Stage 2: scaffold-2 (after PKAs)

| Stage | Agent ID | Started | Background ID | Output path | Status |
|-------|----------|---------|---------------|-------------|--------|
| 2.a | scaffold2-master-wire | 02:30 EDT | `btvy6zb3i` | `catalog_master.tex` + 3 `processes/secondary/<family>/index.tex` | DONE (commit `103833f`) |

### Stage 3: WA (4 batches) — DONE

All 4 batches committed WRITER-DONE on first cycle (8/8 entries: 2x WRITER-DONE blocks each = 1 transition).

| Stage | Agent ID | Started | Background ID | Output path | Status |
|-------|----------|---------|---------------|-------------|--------|
| 3.a | WA-w8-kaon-LFV | 02:30 EDT | `buwec8cf5` | K019/K020/K021 | DONE (commit `f9c7ff4`, WRITER-DONE) |
| 3.b | WA-w8-B-rare-leptonic | 02:30 EDT | `bxnxm22g9` | B007/B008 | DONE (commit `88f2cde`, WRITER-DONE) |
| 3.c | WA-w8-B-radiative | 02:30 EDT | `bp11p6e5r` | B013/B014 | DONE (commit `2b23464`, WRITER-DONE) |
| 3.d | WA-w8-T014 | 02:30 EDT | `bj3bejnp7` | T014 | DONE (commit `d198787`, WRITER-DONE) |

### Stage 4: CA — DONE

Each CA is a fresh codex session (writer/checker separation per
AGENTIC_WORKFLOW §2).

| Stage | Agent ID | Background ID | Cycle-1 result |
|-------|----------|---------------|----------------|
| 4.a | CA-w8-kaon-LFV | `b64c0s7gq` | K019 PASS, K020 **WRITER-REWORK** (CHK-1: Sher/E865 limit not in pdg_or_equivalent), K021 PASS — commit `41069e7` |
| 4.b | CA-w8-B-rare-leptonic | `by14g940i` | B007/B008 PASS, all 9 checks — commit `e9fcdf3` |
| 4.c | CA-w8-B-radiative | `bq9ndu8bo` | B013/B014 PASS, all 9 checks — commit `b3c35d7` |
| 4.d | CA-w8-T014 | `by475c2ef` | T014 PASS, all 9 checks — commit `749d789` |

### Stage 4b: Cycle-2 rework for K020 (only) — DONE

Triggered by 4.a CHK-1 fail (legitimate, not L001 carve-out). Surgical
fix per L001 / T003 precedent: promoted measured observable to
`pdg_or_equivalent.values`.

| Stage | Agent ID | Background ID | Result |
|-------|----------|---------------|--------|
| 4.e | WA-v2-K020 | `bwkk2ggxn` | DONE; added `Sher2005:K020:Kplus_piplus_mup_em_e865_only_limit` to pdg_or_equivalent.values, WRITER-DONE cycle 2 — commit `5c14d2f` |
| 4.f | CA-v2-K020 | `bkdvi62ng` | PASS on all 9 checks; CHECKER-DONE cycle 2; `checker_passed_at: 2026-05-17T12:04:08-04:00` — commit `a20d75a` |

**End-of-stage-4 state**: all 8 SECONDARY entries CHECKER-DONE.

### Stage 5: fact-check (3 family agents)

| Stage | Agent ID | Started | Background ID | Output path | Status |
|-------|----------|---------|---------------|-------------|--------|
| 5.a | factcheck-codex-kaon-w8 | 12:05 EDT | `bvw124s2n` (log `/tmp/wave8_logs/factcheck_kaon.log`) | `audits/factcheck_kaon.md` w8 addendum | DISPATCHED |
| 5.b | factcheck-codex-beauty-w8 | 12:05 EDT | `bj00yv7cx` (log `/tmp/wave8_logs/factcheck_beauty.log`) | `audits/factcheck_beauty.md` w8 addendum | DISPATCHED |
| 5.c | factcheck-codex-tophiggsew-w8 | 12:05 EDT | `b23y427q9` (log `/tmp/wave8_logs/factcheck_top.log`) | `audits/factcheck_top_higgs_ew.md` w8 addendum | DISPATCHED |

### Stage 6: Opus sign-off

| Stage | Agent ID | Started | Background ID | Output path | Status |
|-------|----------|---------|---------------|-------------|--------|
| 6.a | Opus-round-4 | (pending) | — | `signoff/round_004_index.md` | PENDING |

### Stage 7: compile + tag

| Stage | Agent ID | Started | Background ID | Output path | Status |
|-------|----------|---------|---------------|-------------|--------|
| 7.a | master-compile-v03 | (pending) | — | `catalog_master.pdf` + `master_compile_v03_report.md` | PENDING |
| 7.b | tag-v0.3 | (pending) | — | git tag `flavor-catalog-v0.3` | PENDING |

---

## 5. Recovery instructions (if orchestrator is replaced or compacted)

1. Read this file end-to-end. The dispatch ledger is authoritative for what's
   running and what's pending.
2. Read `flavor_catalog/HANDOFF_PROMPT.md`, `SESSION_NOTES.md`,
   `AGENTIC_WORKFLOW.md`, and the 5 arbitration docs.
3. Run pre-flight: `git status -sb`, `git log --oneline -5`,
   `~/bin/codex_check_usage.sh`.
4. For any PENDING agent: dispatch per its row description, using the prompt
   in `/tmp/wave8_prompts/<agent_id>.prompt` if present, or rebuild from the
   PKA/WA/CA/fact-check templates already used in Waves 1–7
   (`worklogs/pka/T003.md` is a good Wave-7 reference). All Wave-8 entries
   are SECONDARY, so reinforce the `priority_tier: SECONDARY` requirement
   and the `processes/secondary/<family>/` output paths in every prompt.
5. For any agent IN-FLIGHT (background ID present, no "completed" mark): use
   Monitor on that ID or just wait — re-dispatching duplicates work.
6. On COMPLETED stage, mark this file and proceed to the next stage.
7. NEVER push to `paper/quark-scan-2026q2`. ALWAYS invoke codex with
   `< /dev/null`. ALWAYS run `~/bin/codex_check_usage.sh` before a new
   dispatch wave.

---

## 6. Reproducibility note

This wave's protocol is the Wave-7 protocol with two adjustments:
- All deliverable paths use `processes/secondary/<family>/` not `processes/<family>/`.
- All YAML sidecars carry `priority_tier: SECONDARY`, `priority_rationale: …`,
  `promoted_in_wave: 8`.

Everything else (PKA → WA → CA → fact-check → Opus, robust git push protocol,
CHK-1..8 checklist, L001 carve-out, cycle caps) is unchanged from
`AGENTIC_WORKFLOW.md`. A future Wave-9 SECONDARY pass is the same recipe with
`promoted_in_wave: 9` and a sibling runbook.
