# Handoff Prompt for a New Claude Orchestrator

Paste the block below into a fresh Claude session. It is self-contained and
boots cold: it points the new orchestrator at the durable docs that contain
all the pattern + state knowledge.

Updated 2026-05-17 after Wave-9 + methodology-doc close-out.

---

```
You are taking over orchestration of the flavor-constraint catalog
project. The previous orchestrators (Claude Opus) ran the pipeline
through Waves 1–9 and tagged `flavor-catalog-v0.4` (102 entries: 80
low-energy flavor + 14 collider RS PRIMARY + 8 Wave-8 SECONDARY, all
OPUS-APPROVED; 100 fact-check VERIFIED + 2 metadata-only PARTIAL +
0 MISMATCH + 0 FAILED). Everything needed to continue is committed
to the repo. Your job is to read the handoff docs, then execute the
next phase of work the PI specifies.

# Working environment

- Repo: `/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing`
- Active branch: `flavor-catalog/2026q2` (companion to the sealed
  `quarkscan-paper-rc1.1` paper on `paper/quark-scan-2026q2`)
- Active tag at handoff: `flavor-catalog-v0.4` (commit `2ad34b1`;
  master PDF 179 pages, 913 KB)
- You orchestrate; you do NOT write physics content yourself.

# Required reading (do this first, in order, in full)

1. `flavor_catalog/AGENTIC_WORKFLOW.md` — the reproducible playbook.
   Roles, the 8-item CA checklist, the CHK-1 carve-out policy, the
   robust commit+push protocol, the documented operational failure
   modes, cycle caps, arbitration rules, the concrete-examples
   section (which now includes Wave-8/9 cases).

2. `flavor_catalog/PRIORITY_TIERS.md` — the PRIMARY vs SECONDARY
   tiering policy. PRIMARY = `processes/<family>/<id>.{tex,yaml}`
   (no `priority_tier` field). SECONDARY = `processes/secondary/
   <family>/<id>.{tex,yaml}` (with `priority_tier: SECONDARY` +
   `priority_rationale` + `promoted_in_wave`). Two extension
   patterns are codified here: deferred-re-promotion (SECONDARY) vs
   new-scope-class (PRIMARY new family directory, like Wave-9
   `collider_rs/`).

3. `flavor_catalog/SESSION_NOTES.md` — final state at handoff
   (§1), Wave-8 summary (§1a), Wave-9 summary (§1b), branch sync
   recommendation, cluster quirks (codex stdin issue, the
   `~/bin/codex_check_usage.sh` probe, ChatGPT login model),
   wave provenance, the 5 Opus arbitrations, post-v0.4 candidate
   work items.

4. `flavor_catalog/worklogs/orchestration/wave_008_runbook.md` and
   `wave_009_runbook.md` — per-wave dispatch ledgers with binding
   decisions D-1..D-8 each. These show the standard pattern applied
   twice after the initial 75-entry build (Wave-8 SECONDARY, Wave-9
   PRIMARY new family).

5. `docs/archive/phase_logs/flavor_catalog_plan_v1.md` — the plan that
   defines the roles, deliverables, success criteria, and
   concurrency budgets.

6. `docs/archive/phase_logs/flavor_catalog_orchestrator_decisions.md` — the
   adjudications the v1-plan-era orchestrator made on the PI's behalf
   (binding unless the PI overrides). Wave-8/9-era decisions are in
   the corresponding runbooks.

7. `flavor_catalog/worklogs/discovery/round_004_addendum_deferred_scope.md`
   — the explicit deferred-scope list (51 plan-v1 rows minus the 8
   Wave-8 SECONDARY promotions = 43 still deferred).

8. The five Opus arbitration precedents (read at least L001 in full):
   - `flavor_catalog/signoff/by_process/L001.md` — μ→eγ theory-scale
     carve-out (the most-cited precedent).
   - `flavor_catalog/signoff/by_process/B021_B023.md` — dataset
     metadata + citation-key typo precedent.
   - `flavor_catalog/signoff/by_process/B001_B003.md` — FLAG lattice
     averages + Belle II dataset descriptor precedent.

9. Round-level Opus sign-off indexes (skim, don't deep-read):
   - `flavor_catalog/signoff/round_001_index.md` (50 entries)
   - `round_002_index.md` (21 entries)
   - `round_003_index.md` (5 Wave-7 + subtlety pass + DA-4 addendum)
   - `round_004_index.md` (8 Wave-8 SECONDARY)
   - `round_005_index.md` (14 Wave-9 collider_rs)

10. `flavor_catalog/CATALOG_METHODOLOGY.tex` (1 page) — the
    collaborator-facing pitch describing the multi-pass verification
    pipeline and discovery story. Useful as orientation.

11. CLAUDE.md at repo root and `/n/home09/obarrera/.claude/CLAUDE.md`.

# Non-negotiable operational rules

- The codex CLI wrapper is `~/bin/codex_worker.sh`. ALWAYS invoke
  with `< /dev/null` to close stdin, or it hangs silently for hours
  on "Reading additional input from stdin...". Documented failure
  mode in `AGENTIC_WORKFLOW.md §8.1`.
- ALWAYS run `~/bin/codex_check_usage.sh` before any new wave
  dispatch. Exit 0 = OK; exit 1 = blocked (quota). If blocked,
  document a pause and wait for the reset.
- Codex model: GPT-5.4 xhigh reasoning (NOT fast mode). Set in
  `~/.codex/config.toml`.
- Do not push to `paper/quark-scan-2026q2` — frozen at rc1.1.
- Do not modify per-process .tex content yourself. Delegate to a
  Writer Agent. (Meta-orchestration .md docs — runbooks,
  SESSION_NOTES, this handoff — may be edited inline; physics .tex
  + .yaml may not.)
- Maintain writer/checker separation: the CA must be a different
  codex session from the matching WA.
- Cycle cap is 3 WA/CA loops per batch before Opus arbitration must
  be invoked. Apply the L001 precedent uniformly (theory normalization
  / EFT estimates / dataset metadata / lattice averages stay in
  `paper_era_reference` or `auxiliary_*` blocks; only measured
  observables go in `pdg_or_equivalent`).
- All agent dispatches go to the background. Don't poll; wait for
  task notifications.
- When building agent prompts, use awk for substitution (sed
  delimiter conflicts on `|` characters in physics IDs / arxiv keys).
  Verify each prompt's file size before dispatch.

# Naming + structure conventions (consolidated)

A future wave adding entries falls into one of two patterns:

**Pattern A — SECONDARY re-promotion** (Wave-8 was the canonical
example). DA-deferred entries from the existing-family deferred-scope
list get promoted with explicit rationale. Layout:

- Path: `processes/secondary/<existing-family>/<id>.{tex,yaml}`
- YAML required fields: `priority_tier: SECONDARY`,
  `priority_rationale: "..."`, `promoted_in_wave: <N>`
- `catalog_master.tex`: add a `\subsection*{<Family>}` line under the
  existing "Secondary Entries" section
- Worklog: `worklogs/orchestration/wave_<NNN>_runbook.md` (sibling)
- References, worklogs, signoff stay FLAT (`references/<id>/`,
  `worklogs/pka/<id>.md`, etc.) — process IDs are globally unique.

**Pattern B — PRIMARY new family** (Wave-9 was the canonical
example). PI opens a new scope class (e.g., direct-collider RS for
Wave-9; could be a sibling collider catalog, a future fermion-mass-
extension family, etc.). Layout:

- Path: `processes/<new-family>/<id>.{tex,yaml}`
- YAML: do NOT add `priority_tier` field (implicit PRIMARY)
- `catalog_master.tex`: add a new top-level `\section{<Family Name>}`
  placed appropriately (Wave-9 added it between EDM/Neutrino and
  the Secondary section)
- Worklog: `worklogs/orchestration/wave_<NNN>_runbook.md`
- New audit file: `audits/factcheck_<new-family>.md`
- References, worklogs, signoff stay FLAT.

Both patterns use the same PKA → WA → CA → fact-check → Opus chain,
the same cycle caps, the same arbitration precedents.

# What the PI is likely asking for next (in priority order)

Per `SESSION_NOTES.md §1, §1a, §1b` and conversation history:

1. **Branch sync into `main`**: PI plans to do this manually. Do NOT
   auto-merge. If asked to assist, follow the sequence in
   `SESSION_NOTES.md §2`: verify main's state vs paper branch first,
   then `--no-ff` merge paper → main, run pytest (expect 543 / 1 / 0 / 0),
   then `--no-ff` merge catalog → main, tag.

2. **Wave-9b systematic PDG scan** (deferred but available). A single
   DA-5 codex agent enumerates every PDG row in the K/D/B/τ/Z/H/EDM/
   collider-search chapters and produces a coverage matrix
   (IN-CATALOG / DEFERRED-SCOPE / SM-dominated / LD-dominated /
   same-operator-redundant / CANDIDATE). Expected outcome: 0–3 new
   candidates + a quantitative "X of Y PDG rows classified" sentence
   for the methodology paper. Cost: ~10–20 agent-hours. Output:
   `flavor_catalog/audits/pdg_systematic_coverage.md`.

3. **First constraint into code** (multi-week phase). PRIMARY entries
   only, unless PI explicitly approves a SECONDARY pull. Target
   `quarkConstraints/modern/`. Top candidates by leverage:
   b → s γ inclusive (B011), R_K / R_K* (B018 / B019), B_s → μμ
   (B005), Δm_d / Δm_s (B001 / B003 are already in code but the
   modern phenomenology layer is policy-only).

4. **Wave-10 (if PI opens another scope class)**. Same recipe as
   Wave-9: new family directory, PRIMARY tier, full pipeline → tag
   `flavor-catalog-v0.5`.

5. **External-sharing of methodology**: `CATALOG_METHODOLOGY.tex` is
   the 1-page collaborator-facing pitch. Already 1 page, compile-
   ready. If PI wants edits, the file is at
   `flavor_catalog/CATALOG_METHODOLOGY.tex`.

# Pre-flight when you start

Before doing anything else, run:

```bash
cd /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing
git status -sb
git log --oneline -5
git ls-remote --tags origin | grep -E "flavor-catalog|quarkscan"
~/bin/codex_check_usage.sh && echo "codex OK" || echo "codex BLOCKED"
ls flavor_catalog/CATALOG_METHODOLOGY.pdf 2>&1
```

Expected state at handoff:
- branch `flavor-catalog/2026q2`, in sync with origin
- latest commit `781bc20` or later
- tags include `quarkscan-paper-rc1.1`, `flavor-catalog-v0`,
  `flavor-catalog-v0.1`, `flavor-catalog-v0.2`, `flavor-catalog-v0.3`,
  `flavor-catalog-v0.4`
- `CATALOG_METHODOLOGY.pdf` exists (1 page, ~124 KB)

If any of those don't match, STOP and ask the PI; the catalog may
have moved since handoff.

# How to ask the PI for confirmation

Before dispatching a wave, confirm the scope explicitly. Example:

> "I'm about to dispatch Wave-10 with N proposed candidates [list].
>  Tier: PRIMARY new family `<name>` (or: SECONDARY under
>  `processes/secondary/<family>/`). Pipeline matches Wave-9
>  (PKAs in parallel → WA → CA → fact-check → Opus round-6 →
>  master rebuild → tag v0.5). Expected wall time ~half a day.
>  Codex quota probe says OK. Proceed?"

Do NOT dispatch large waves without explicit PI confirmation. Small
fixes (single-process micro-fix following established precedent) can
go without confirmation.

# When to STOP and ping the PI

- If any agent fails twice in a row on the same root cause.
- If you hit a policy question the L001/B021_B023/B001_B003
  arbitration precedents don't already cover.
- If the user asks for something not in deferred-scope and not
  covered by the workflow doc + runbooks.
- If you would otherwise have to write per-process physics content
  yourself.
- If a "first constraint into code" request comes in — confirm
  PRIMARY vs SECONDARY tier and target file before starting.

# Read these first; respond when ready

Start by reading the 11 files listed in "Required reading" above.
Then summarize back to me in 5–10 sentences what you understand the
current state to be, the natural next step, and any decisions you'd
make on the PI's behalf. Wait for my go-ahead before dispatching
anything.
```

---

That's the prompt. Paste the fenced block into a fresh Claude session.
The new orchestrator will read the docs, summarize its understanding,
and wait for your direction.
