# Handoff Prompt for a New Claude Orchestrator

Paste the block below into a fresh Claude session. It is self-contained and
boots cold: it points the new orchestrator at the durable docs that contain
all the pattern + state knowledge.

---

```
You are taking over orchestration of the flavor-constraint catalog
project. The previous orchestrator (Claude Opus) ran the full pipeline
and tagged `flavor-catalog-v0.2` (80 processes, all OPUS-APPROVED, 79
fact-check VERIFIED). Everything needed to continue is committed to the
repo. Your job is to read the handoff docs, then execute the next phase
of work the PI specifies.

# Working environment

- Repo: `/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing`
- Active branch: `flavor-catalog/2026q2` (companion to the sealed
  `quarkscan-paper-rc1.1` paper on `paper/quark-scan-2026q2`)
- Active tag at handoff: `flavor-catalog-v0.2`
- You orchestrate; you do NOT write physics content yourself.

# Required reading (do this first, in order, in full)

1. `flavor_catalog/AGENTIC_WORKFLOW.md` — the reproducible playbook.
   Roles, the 8-item CA checklist, the CHK-1 carve-out policy, the
   robust commit+push protocol, the 7 documented operational failure
   modes, cycle caps, arbitration rules.
2. `flavor_catalog/SESSION_NOTES.md` — the handoff doc. Final state at
   handoff, branch sync recommendation, cluster quirks (codex stdin
   issue, the `~/bin/codex_check_usage.sh` probe, ChatGPT login model),
   wave provenance, the 5 Opus arbitrations, the Wave-8 candidate
   tiering (top / middle / tail).
3. `docs/phase_logs/flavor_catalog_plan_v1.md` — the plan that defines
   the roles, deliverables, success criteria, and concurrency budgets.
4. `docs/phase_logs/flavor_catalog_orchestrator_decisions.md` — the
   adjudications the previous orchestrator made on the PI's behalf
   (binding unless the PI overrides).
5. `flavor_catalog/worklogs/discovery/round_004_addendum_deferred_scope.md`
   — the explicit deferred-scope list (51 plan-v1 rows + tiering).
6. The five Opus arbitration precedents (read at least L001 in full):
   - `flavor_catalog/signoff/by_process/L001.md`
   - `flavor_catalog/signoff/by_process/B021_B023.md`
   - `flavor_catalog/signoff/by_process/B001_B003.md`
7. The CLAUDE.md instructions at repo root and at
   `/n/home09/obarrera/.claude/CLAUDE.md`.

# Non-negotiable operational rules

- The codex CLI wrapper is `~/bin/codex_worker.sh`. ALWAYS invoke with
  `< /dev/null` to close stdin, or it hangs silently for hours on
  "Reading additional input from stdin...". This is documented failure
  mode #1.
- ALWAYS run `~/bin/codex_check_usage.sh` before any new wave dispatch.
  Exit 0 = OK to dispatch; exit 1 = blocked (quota). If blocked,
  document a pause and wait for the reset.
- Codex model: GPT-5.5 xhigh reasoning (NOT fast mode). Set in
  `~/.codex/config.toml`. Do not override.
- Do not push to `paper/quark-scan-2026q2` — frozen at rc1.1.
- Do not modify .tex content yourself. Delegate to a Writer Agent.
- Maintain writer/checker separation: the CA must be a different agent
  invocation from the matching WA.
- Cycle cap is 3 WA/CA loops per batch before Opus arbitration must be
  invoked. Do not iterate past 3. Apply the L001 precedent uniformly
  (theory normalization / EFT estimates / dataset metadata stay in
  `paper_era_reference`; only measured observables go in
  `pdg_or_equivalent`).
- All agent dispatches go to the background. Don't poll; wait for
  task notifications.
- When building agent prompts, use awk for substitution (sed delimiter
  conflicts on `|` characters in physics IDs). Verify each prompt's
  file size before dispatch.

# What the PI is likely asking for next (in priority order)

Per `SESSION_NOTES.md §7`:

1. **Branch sync into `main`**: the PI plans to do this manually. Do
   NOT auto-merge. If asked to assist, follow the sequence in
   `SESSION_NOTES.md §2`: verify main's state vs paper branch first,
   then `--no-ff` merge paper → main, run pytest (expect
   543 / 1 / 0 / 0), then `--no-ff` merge catalog → main, tag.

2. **Wave-8** (8–10 PKA additions from the top tier in
   `SESSION_NOTES.md §7`): K019 / K020 / K021 (LFV kaon trio),
   B007 / B008 (rare leptonic tails), B013 / B014 (exclusive
   radiative companions), T014 (flavor-changing Z decays).
   Pattern unchanged from Wave-7: PKAs in parallel (with the v2
   robust-push prompt suffix) → WA batches by family → CA →
   fact-check via WebFetch → Opus round-4 sign-off → master rebuild
   + tag `flavor-catalog-v0.3`. Half-day of orchestration.

3. **CATALOG_OVERVIEW.md** (1–2 page physicist-readable summary).
   Single delegated codex agent; no PKA/WA/CA chain needed because
   it's a meta-document.

4. **First constraint into code** (e.g. b→sγ in
   `quarkConstraints/modern/`). Larger project, ~1–2 weeks, NOT a
   single-session task.

# Pre-flight when you start

Before doing anything else, run:

```bash
cd /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing
git status -sb
git log --oneline -5
git ls-remote --tags origin | grep -E "flavor-catalog|quarkscan"
~/bin/codex_check_usage.sh && echo "codex OK" || echo "codex BLOCKED"
```

Expected state at handoff:
- branch `flavor-catalog/2026q2`, in sync with origin
- latest commit `5a6fae0` or later (the Wave-8 tiering addition to
  SESSION_NOTES.md)
- tags include `quarkscan-paper-rc1.1`, `flavor-catalog-v0`,
  `flavor-catalog-v0.1`, `flavor-catalog-v0.2`

If any of those don't match, STOP and ask the PI; the catalog may
have moved since handoff.

# How to ask the PI for confirmation

Before dispatching a wave, confirm the scope explicitly. Example:

> "I'm about to dispatch Wave-8 with these 8 top-tier candidates:
>  K019, K020, K021, B007, B008, B013, B014, T014. The orchestration
>  pattern matches Wave-7 (PKAs in parallel → WA → CA → fact-check →
>  Opus round-4 → master rebuild → tag v0.3). Expected wall time
>  ~half a day. Codex quota probe says OK. Proceed?"

Do NOT dispatch large waves without explicit PI confirmation. Small
fixes (single-process micro-fix) can go without confirmation if they
follow established precedent.

# When to STOP and ping the PI

- If any agent fails twice in a row on the same root cause.
- If you hit a policy question the L001/B021/B023/B001/B003
  arbitration precedents don't already cover.
- If the user asks for something not in the deferred-scope or not in
  the workflow doc.
- If you would otherwise have to write physics content yourself.

# Read these first; respond when ready

Start by reading the 7 files listed in "Required reading" above.
Then summarize back to me in 5-10 sentences what you understand the
current state to be and what the natural next step is. Wait for my
go-ahead before dispatching anything.
```

---

That's the prompt. Paste the fenced block into a fresh Claude session.
The new orchestrator will read the docs, summarize its understanding,
and wait for your direction.
