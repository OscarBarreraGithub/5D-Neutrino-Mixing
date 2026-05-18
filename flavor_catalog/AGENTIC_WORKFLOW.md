# Agentic Workflow — Reproducible Playbook

This document is a generic playbook for building a research catalog (or any
similar large, verification-heavy artifact) using an orchestrator + delegated
agents. It is written so a new orchestrator (a fresh Claude session, a new
collaborator, or a new project entirely) can pick up the pattern and run it.
The concrete examples are drawn from `flavor_catalog/`, but the structure
generalizes.

The TL;DR: one orchestrator never reads physics or writes content; it
dispatches specialized agents in parallel, runs cross-review chains, and
arbitrates only when forced to. The orchestrator's tokens go to management,
not content.

---

## 1. When to use this pattern

Use it when:

- The artifact is many small independent units (e.g. one entry per process,
  one section per topic).
- Each unit must be sourced, fact-checked, and externally defensible.
- A single agent doing all the work would either run out of context, drift in
  quality, or be impossible to verify.
- You can afford 3–10× more compute by spawning parallel agents in exchange
  for higher reliability and reproducibility.

Don't use it when the artifact is small enough for one careful agent to
handle, or when the verification step is cheap relative to the writing step.

---

## 2. Roles

Each agent's responsibilities, model choice, and budget are fixed. Mixing
them defeats the cross-review safety.

| Role | Model | Purpose | Output |
|---|---|---|---|
| **Orchestrator** | Claude Opus 4.x | Plan, dispatch, arbitrate, record decisions. **Does not write content.** | Worklogs, commits, sign-off docs |
| **Process Knowledge Agent (PKA)** | Codex GPT-5.5 xhigh (not fast) | Per-unit research: download sources, extract values, write the first draft. | Per-unit `.tex` + `.yaml` sidecar + reference snapshots + worklog |
| **Writer Agent (WA)** | Codex GPT-5.5 xhigh | Polish a batch of PKA drafts: tighten language, complete metadata, fix LaTeX. Does not invent content. | Updated `.tex`/`.yaml` + batch worklog |
| **Checker Agent (CA)** | Codex GPT-5.5 xhigh, **must be a separate session from the matching WA** | Independent verification of a WA-polished batch. Runs an 8-item checklist (see below). Verdict: PASS / WRITER-REWORK. | CA worklog with per-unit pass/fail matrix + verified value table |
| **Discovery Agent (DA)** | Codex GPT-5.5 xhigh | Survey the catalog vs. a target list; propose additions or recommend convergence. Multiple rounds. | Discovery worklog |
| **Fact-Check Agent** | Codex GPT-5.5 xhigh with WebFetch | Independently re-fetch every cited source and verify that the claimed values actually appear in the source content. | Per-family fact-check report; aggregated checklist |
| **Opus sign-off agent** | Claude Opus 4.x | Final per-unit judgment: APPROVE / RETURN-TO-WA / ESCALATE-TO-PI. Reads the chain, applies the policy. | Round-level signoff index; per-unit signoff doc only for non-APPROVE |
| **Opus arbitration agent** | Claude Opus 4.x | Called when WA/CA cycle cap is hit. Applies policy carve-outs and either overrides or returns. One corrective loop before PI escalation. | Arbitration doc + sidecar status update |

The writer/checker separation rule is **non-negotiable**. Two real BLOCKER-class
physics errors in the quark-scan paper were caught only because the CA was
genuinely independent from the WA.

---

## 3. Standard wave pipeline

A "wave" is the atomic unit of orchestration. One wave processes N units
through the full chain in parallel.

```
                       (per unit, in parallel)
PKA × N    ────►  produces per-unit drafts on a feature branch
                       │
                       ▼
WA × ⌈N/batch⌉  ────►  polishes drafts in family batches of 3–5
                       │
                       ▼
CA × ⌈N/batch⌉  ────►  independently verifies; produces per-unit verdict
                       │
              ┌────────┴────────┐
              ▼                 ▼
        PASS (CHECKER-DONE)  WRITER-REWORK
              │                 │
              │                 ▼
              │           WA-v2 (cycle 2) ──► CA-v2 ──┐
              │                                       │
              │                                       ▼
              │                                 (cycles capped at 3)
              │                                       │
              │                                       ▼
              │                                 Opus arbitration
              ▼
        Fact-check pass (WebFetch every cited source)
              │
              ▼
        Opus sign-off round (per-unit verdict matrix)
              │
              ▼
        Master compile + tag vN
```

### After every wave converges
- Run a **Discovery Agent** to find gaps relative to the target list.
- If DA proposes a meaningful batch, dispatch the next wave.
- Convergence rule: **two DA rounds in a row find ≤ 1 new unit AND every
  existing unit is OPUS-APPROVED.** Then stop discovery and mark the
  remaining target-list tail explicitly DEFERRED-SCOPE.

### Status states (the YAML status_history vocabulary)
```
DRAFT
   │
   ▼
WRITER-INITIATED ── PKA wrote the initial draft
   │
   ▼
WRITER-DONE      ── WA polish (may be on cycle 1, 2, or 3)
   │
   ├──► WRITER-REWORK  ── CA returned the unit; resubmit
   │         │
   │         ▼
   │     (loop up to plan-v1 :620 max-3 cycles)
   │
   ▼
CHECKER-DONE     ── CA passed all 8 checks
   │
   ▼
OPUS-APPROVED    ── Opus signed off
   │
   ▼
FACT-CHECKED     ── independent WebFetch verification of cited sources
```

Cap-hit transitions force `BLOCKED-PI` (waiting on a human decision) or
`DEFERRED-SCOPE` (explicitly removed from active set) so no unit can be
silently open.

---

## 4. The 8-item Checker checklist (canonical)

Every CA runs these against every unit in its batch. Adapt the names but
keep the structure.

| Check | What it verifies |
|---|---|
| **CHK-1** | Every numerical claim in the unit traces to a structured `pdg_or_equivalent` (or analogous) value block with year / value / uncertainty / units / source URL / access date / sha256. |
| **CHK-2** | Every cited reference key resolves to a `source_manifest.yaml` entry; every manifest entry has a non-empty `snapshot_path` pointing to a tracked text file. |
| **CHK-3** | No publisher PDFs are tracked under the unit's `references/` directory (text snapshots and arXiv preprints only). |
| **CHK-4** | The `status_history` shows the legal transition `WRITER-INITIATED → WRITER-DONE` (or further) with ISO 8601 timestamps. |
| **CHK-5** | "In-code" coverage claims (`YES`/`PARTIAL`/`NO`) actually grep-match the cited `file:line` evidence. |
| **CHK-6** | Implementation-difficulty rating is consistent with the rubric (LOW = existing operator basis; MEDIUM = new operator, standard inputs; HIGH = new RG/lattice/mode calculation). |
| **CHK-7** | No load-bearing number from upstream artifacts (the companion paper, the methodology note) is contradicted. |
| **CHK-8** | No accidental LaTeX issues: math mode, missing braces, unresolved `\cite` / `\ref`, balanced delimiters. |

A unit must pass all 8 to advance to OPUS-APPROVED.

---

## 5. The CHK-1 carve-out (the most important policy lesson)

The single biggest source of false-positive CA failures in this project was
strict CHK-1 reading applied to numbers that aren't measured observables.
Three categories of numbers should **never** live in `pdg_or_equivalent`:

1. **Theory normalization scales / reference values.** Example: the
   Csaki-Falkowski-Weiler "3 TeV / 21 TeV KK reference scale" — it's a
   convention input, not a measurement.
2. **Order-of-magnitude SM estimates.** Example: "SM prediction is
   $\mathcal{O}(10^{-14})$" for a top FCNC mode — it's an EFT estimate.
3. **Dataset/luminosity descriptors.** Example: "LHCb 9 fb⁻¹" — that's
   metadata about a measurement, not the measurement itself.

These belong in `paper_era_reference`, `theory_inputs`, or
`supporting_measurements` with full provenance (URL, sha256, etc.), but the
top-level `pdg_or_equivalent` block is reserved for measured experimental
observables only.

The first unit to hit this was **L001** (μ → eγ): a strict CA insisted on
promoting the "3 TeV reference scale" into `pdg_or_equivalent`. Opus
arbitration ruled APPROVE-OVERRIDE and codified the policy. Subsequent
arbitrations (B021, B023, B001, B003) applied the same precedent. By round 2,
CAs had internalized it and stopped flagging false positives.

**Lesson:** the first time a recurring CA false-positive appears, dispatch an
Opus arbitration to codify the policy; don't keep iterating WA/CA forever.
Future CAs read the precedent doc and apply uniformly.

The recurring example: a wave-7 unit (T003) tried to keep
"SM expectation remains at the $10^{-14}$ level" in its prose. CHK-1 flagged
it because the number wasn't in `pdg_or_equivalent`. The fix per L001
precedent: **remove the number from prose**, keep the qualitative "many
orders of magnitude below current bounds", and leave the source under
`paper_era_reference` for traceability. That kind of cycle-3 fix is much
faster and more honest than fighting the strict rule.

---

## 6. Concurrency budgets

Budget is per orchestrator turn:

| Role | Max parallel | Notes |
|---|---|---|
| PKA | 8–12 | One per unit. Hit the upper edge only if all units are independent and the model provider is fresh. |
| WA | 4 | Batches of 3–5 units each. |
| CA | 4 | Must be different agents from WAs (writer/checker separation). |
| DA | 1 at a time | Discovery is sequential because each round depends on the previous. |
| Opus sign-off | 1 (rounds 1–2), 2 (rounds 3+) | Default to 1 to lock global consistency. |
| Opus arbitration | 1 at a time per unit | Apply policy carefully. |
| Fact-check | 6 (one per family) | All WebFetch-heavy; expect rate limits. |

**Total concurrent agents:** keep under 20. The orchestrator must track this
mentally.

---

## 7. Infrastructure

### 7.1 Wrapper for the model CLI

Wrap the codex CLI (or equivalent) in a script that:
- Sets PATH so the model binary is found.
- Passes `--dangerously-bypass-approvals-and-sandbox` (or equivalent) so the
  CLI doesn't hang waiting for human approval.
- Closes stdin: invoke via `bash wrapper.sh ... < /dev/null`. Without this,
  the CLI prints "Reading additional input from stdin..." and hangs silently
  for hours. **This is the #1 failure mode.**

Example: `~/bin/codex_worker.sh` in this repo.

### 7.2 Pre-flight quota probe

Build a probe that runs `<cli> /status` and exits non-zero on a usage-limit
error. Run it before every dispatch wave. Example:
`~/bin/codex_check_usage.sh`.

Behavior: if probe returns BLOCKED, do not dispatch; wait, log a pause
record, schedule a retry when the quota window resets. Example pause record:
`docs/phase_logs/flavor_catalog_codex_quota_pause.md`.

### 7.3 Branch policy

- The active development happens on a **dedicated feature branch**
  (`flavor-catalog/2026q2` in this project), cut from a stable release tag of
  the parent artifact.
- The parent's `main` (or release-candidate tag) is **frozen**; the feature
  branch never pushes back to it without explicit human review.
- Tag releases (e.g. `flavor-catalog-v0`, `v0.1`, `v0.2`) at each major
  convergence point.
- Multiple agents can push to the same feature branch in parallel if they
  write to **disjoint file paths**. The "robust push protocol" handles the
  rebase fan-out.

### 7.4 Robust commit + push (the per-agent boilerplate)

Every WA/CA/PKA/fact-check agent uses this loop. It tolerates concurrent
push fan-out without overwriting other agents' work:

```bash
git add <my files only>
git commit -m "<descriptive>"

for attempt in 1 2 3 4 5; do
  git fetch origin <branch>
  git rebase origin/<branch>
  if [ $? -ne 0 ]; then
    UNRESOLVED=$(git diff --diff-filter=U --name-only)
    REAL=""
    for f in $UNRESOLVED; do
      case "$f" in
        <my-write-paths>) REAL="$f"; break ;;
      esac
    done
    if [ -n "$REAL" ]; then
      git rebase --abort
      echo "FAILED: real conflict on $REAL"
      exit 1
    fi
    git rebase --abort 2>/dev/null
    sleep $((RANDOM % 5 + 1))
    continue
  fi
  if git push origin <branch>; then
    exit 0
  fi
  sleep $((RANDOM % 5 + 1))
done
exit 1
```

The critical distinction: a rebase merging in a file YOU did not touch is
**not** a conflict (the working tree just gets that file added). Only an
unresolved-file marker matching one of YOUR write paths counts as a real
collision. Earlier agents over-aborted on benign fetches; the protocol above
fixes it.

---

## 8. Operational failure modes seen in this project

These are the concrete failures we hit. Including them here is the whole
point of this document — future orchestrators avoid the same traps.

### 8.1 Silent stdin hang

The model CLI was invoked without `< /dev/null` and silently waited on stdin
for **8 hours** before being killed. Always close stdin. Always spot-check
stdout file size within 60–90 seconds of dispatch; if it's still the size of
the wrapper banner (~40 bytes), the agent has hung.

### 8.2 Shell substitution / delimiter conflicts when building prompts

When generating per-agent prompts via `sed s|placeholder|content|`, the
content sometimes contained `|` (e.g. an arXiv reference key with
`|q/p|_D`). Sed silently fails or produces empty prompts. **Use `awk -v`
substitution instead** — it doesn't care about content characters. Always
verify prompt file size before dispatching.

### 8.3 Race conditions on the shared master-status file

Six fact-check agents in parallel cannot all update one
`factcheck_status.md`. They will conflict on the rebase. Pattern:
- Each agent writes only to its own `<family>` files (process YAMLs +
  per-family report).
- A **separate aggregator agent** runs after all per-family agents land,
  reads the per-family reports, and writes the consolidated master file.

### 8.4 Quota-window collisions during long runs

The model provider's usage budget resets on a clock you can't predict
perfectly. If you hit the limit mid-batch, in-flight sessions usually
complete but new dispatches fail. Pre-flight probe before each batch.
Document the reset time when you see the error message so you can plan the
restart.

### 8.5 The CHK-1 metadata pattern

~80% of first-pass PKA drafts will fail CHK-1 in CA cycle 1 because PKAs
naturally write numbers into auxiliary blocks rather than the strict-format
`pdg_or_equivalent`. Two responses:
- **Cycle 2 WA-v2**: standard fix; promote-or-remove. About 70% pass on
  cycle 2.
- **Cycle 3 → Opus arbitration**: if the failure is on theory normalization
  or dataset metadata, apply the L001 precedent.

Budget for ~1.5 WA/CA cycles per batch on average.

### 8.6 Codex CLI rollout-recorder error

You will see lines like:
```
ERROR codex_core::session: failed to record rollout items: thread X not found
```
These are internal logging errors and do **not** affect the agent's output.
The deliverable file is still written. Ignore them.

### 8.7 "Reading additional input from stdin..." silent hang

This is the canonical symptom of the stdin issue. If you ever see this in
stdout, the agent is hung. Kill the process group, fix the invocation, and
retry.

---

## 9. The arbitration / cycle-cap rules

From plan v1 :620–622 in this project, generalized:

| Cycle | What it is | Cap |
|---|---|---|
| PKA micro-iterations | Source fixes within a single PKA call | 2 |
| WA / CA full cycles per batch | Full polish + verify loop | 3 |
| DA discovery rounds | Full catalog gap-finding | 4 total |
| Opus corrective loops | Per arbitration call | 1 before PI escalation |

When a cap is hit, the orchestrator **must** force a sidecar status
transition (`BLOCKED-PI` or `DEFERRED-SCOPE`) plus a sign-off log entry, so
no unit is silently open.

When in doubt, escalate to Opus arbitration. Opus is cheap relative to
WA/CA cycles that keep failing for the same reason.

---

## 10. External-research integration round

This pattern saved real bookkeeping gaps in this project.

1. **Get an independent review** of the catalog (e.g. ask GPT Deep Research,
   ask a collaborator, run a different model).
2. **Import the review artifact** into a dedicated subdirectory
   (`external_research/`).
3. Convert PDFs to text up front so model agents can read them.
4. **Dispatch read-only review agents** (one per source document) that
   produce a structured comparison report: Missing / Disagreements /
   Subtleties / Overall assessment.
5. Synthesize the reports for the PI; flag actionable gaps.
6. Dispatch a follow-up **Wave-N** with new PKAs / a subtlety-writer /
   a deferred-scope addendum agent — all going through the standard
   PKA→WA→CA→Opus chain.

In this project, the external review caught 5 plan-v1 process IDs that were
ambiguously closed by DA-4 (T003 / T004 / T008 / T012 / B012). They would
have been silent gaps in v0 forever.

---

## 11. Fact-check round (the final defense against hallucination)

Even after Opus sign-off, you do **not** know that the cited references
actually contain the values they're cited for. The PKA could have
hallucinated an arXiv ID, or quoted a slightly-off number.

The fact-check round dispatches one agent per family. Each:
1. Reads the family's catalog entries.
2. For every cited numerical claim, WebFetches the cited URL (arXiv abs
   page, PDG listing, HFLAV averages, journal page).
3. Verifies authors / title / year match the manifest entry.
4. Verifies the specific numerical value appears in the fetched content or
   the local text snapshot.
5. Records VERIFIED / PARTIAL / MISMATCH / FAILED per claim per unit.

Carve-out the policy items from §5 (theory normalization, EFT estimates,
dataset descriptors) — those aren't measured values; don't mark them
MISMATCH.

A single MISMATCH usually means a small digit slip (in this project: a
"6.1 × 10⁻⁵" should have been "6.2 × 10⁻⁵"). Trivial WA fix; re-fact-check;
retag.

---

## 12. Reproducibility checklist

For the artifact to be re-buildable from the source tree, you need:

- [ ] Per-unit `.tex` and `.yaml` sidecar tracked in git.
- [ ] Per-unit `references/<id>/source_manifest.yaml` + minimal text
      snapshots tracked, with sha256 recorded.
- [ ] No publisher PDFs tracked (licensing policy).
- [ ] Per-unit worklog trail: PKA / Writer / Checker / Discovery / Opus.
- [ ] Per-round sign-off index: `signoff/round_NNN_index.md`.
- [ ] Per-unit arbitration docs only where needed:
      `signoff/by_process/<id>.md`.
- [ ] Master TeX file that compiles standalone via `pdflatex` twice (no
      external dependencies).
- [ ] Fact-check status file + per-family reports.
- [ ] DA round worklogs + DA-N convergence verdict + any addendum closing
      bookkeeping.
- [ ] Plan document (this project: `flavor_catalog_plan_v1.md` in
      `docs/phase_logs/`) that captures the agent role definitions and
      target-list rationale.
- [ ] Orchestrator decisions document recording every adjudication made on
      the human's behalf.
- [ ] Tier policy document (this project: `PRIORITY_TIERS.md`) defining
      PRIMARY vs SECONDARY layout + the two extension patterns
      (SECONDARY re-promotion, PRIMARY new family).
- [ ] Per-wave runbook (`worklogs/orchestration/wave_NNN_runbook.md`)
      recording the wave's binding decisions, dispatch ledger, and
      recovery instructions. Critical for compaction-survival.
- [ ] Collaborator-facing one-page methodology pitch (this project:
      `CATALOG_METHODOLOGY.tex`) so external readers can understand the
      verification chain without reading the full playbook.
- [ ] This workflow document (`AGENTIC_WORKFLOW.md`) or equivalent so a
      new orchestrator can resume.
- [ ] Tag at every major convergence point so external readers can pin
      to a specific state.

If all of those exist, anyone can `git checkout <tag>` and re-build the
catalog from scratch with full provenance.

---

## 13. What the orchestrator never does

- Does not read full process content end-to-end (only short verdict lines).
- Does not write per-unit physics or values.
- Does not run pytest or compile LaTeX directly.
- Does not WebFetch sources directly.
- Does not arbitrate policy on its own — delegate to Opus arbitration if a
  cycle cap is hit.

The orchestrator's signature is: `dispatch + read short verdict + decide
next dispatch`. Anything heavier should be a delegated agent.

---

## 14. Concrete examples from this project

To illustrate, here are real chains:

- **B001 / B003 (already-in-code ΔF=2 mixing)**: cycle 1 CHK-1 FAIL → cycle 2
  CHK-1 FAIL on the same theory-normalization line → Opus arbitration
  APPROVE-OVERRIDE per L001 precedent. Total elapsed: ~3 hours.
  Final state: `OPUS-APPROVED via arbitration`. Sign-off doc:
  `flavor_catalog/signoff/by_process/B001_B003.md`.

- **T003 (t→cγ)**: PKA → WA → CA FAIL on "$10^{-14}$" prose → WA-v2 → CA-v2
  FAIL on same line → WA-v3 simply removed the explicit "$10^{-14}$" from
  the prose → CA-v3 PASS → fact-check VERIFIED → Opus round-3 APPROVE.
  Lesson: when the policy precedent says "don't promote theory
  normalization to pdg_or_equivalent", the right cycle-3 fix is to
  **delete the number from the .tex**, not to keep iterating against the
  CA's strict reading. Half the WA/CA cycle in this project would have
  been saved by recognizing this earlier.

- **T020 (h→eμ ATLAS limit)**: PKA / WA / CA / Opus all PASS. Fact-check
  WebFetch found a literal digit slip: "6.1 × 10⁻⁵" should have been
  "6.2 × 10⁻⁵". One-line WA-v3 fix. This is the **only** numerical
  MISMATCH the fact-check round caught across 80 processes. Lesson: the
  fact-check round is worth the cost; it catches things the
  internal-consistency chain can't.

- **External review (Wave-7)**: GPT Deep Research surfaced 5 plan-v1
  process IDs that were ambiguously closed by DA-4 — T003, T004, T008,
  T012, B012. Adding them was ~half a day's work and meaningfully
  improved coverage. Lesson: budget for at least one external-review
  round before declaring the catalog done.

- **K020 cycle-2 promotion (Wave-8)**: PKA → WA → CA cycle-1 FAIL on
  CHK-1 — the entry's prose quoted the Sher / BNL E865 measured upper
  limit `B(K^+ \to \pi^+ \mu^+ e^-) < 2.1 \times 10^{-11}` at 90% CL,
  but the structured `pdg_or_equivalent.values` block was missing that
  specific row. Unlike the L001 / T003 cases, this *is* a measured
  observable, not a theory normalization scale; the right cycle-2 fix
  was to **promote**, not delete. WA-v2 added the row; CA-v2 PASS.
  Lesson: the L001 / T003 carve-out is asymmetric — measured
  observables should be promoted into the structured block, theoretical
  normalizations should not be. Decide which by reading the source.
  Don't blanket-apply T003 to every CHK-1 fail.

- **Wave-9 collider_rs cycle-2 batch (4-of-14)**: CR008 and CR010
  failed CHK-1 because their prose quoted historical mass-exclusion
  numerals from earlier LHC searches that were superseded by the
  canonical current limit already in `pdg_or_equivalent.values`. CR010
  and CR011 also failed CHK-2 because their TeX `Key references`
  sections cited snapshot-filename stems (e.g.
  `atlas_2018_arxiv1808_02343`) instead of the manifest's canonical
  keys (e.g. `ATLAS2018_TBCombination`). Two surgical WA-v2 batches
  applied T003 (remove historical numerals from prose; keep qualitative
  language) and B023 (rename TeX keys to manifest canonical) precedents
  in parallel. CA-v2 PASS for all 4. Lesson: cycle-2 fixes can mix
  precedents — apply the appropriate one per finding, not the same
  one to all.

- **Two extension patterns (Wave-8 vs Wave-9)**: the catalog now has
  two codified ways to add entries. SECONDARY re-promotion (Wave-8)
  pulls from the DA-deferred list and lands under
  `processes/secondary/<family>/`. PRIMARY new family (Wave-9) opens
  a new scope class at PI directive and lands under
  `processes/<new-family>/` with no `priority_tier` field. Same
  pipeline, same arbitration precedents. Convention is documented in
  `flavor_catalog/PRIORITY_TIERS.md §7`.

---

End of playbook. The pattern is portable.
