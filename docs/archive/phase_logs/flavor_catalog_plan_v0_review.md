# Flavor Catalog Plan v0 — Peer Review

**Reviewer**: independent Codex peer reviewer
**Date**: 2026-05-16
**Plan**: `docs/phase_logs/flavor_catalog_plan_v0.md` (commit `ebe8a3c`)

## Verdict
REJECT-WITH-REVISIONS

## Item-by-item (1-11)
1. PASS - Section A gives a concrete root path, `flavor_catalog/`, with a sensible process/reference/worklog/signoff layout, branch policy, figure/source caching policy, license rule, and quantified 200 MB total / 5 MB per-file guidance at lines 17-94.

2. PASS - Section B includes process name, PDG/equivalent source and value metadata, RS relevance, post-2008 developments, validity/model dependence, code coverage, implementation difficulty, and references in the YAML and TeX skeletons at lines 100-165.

3. PASS - Section C lists 128 entries, above the >=80 target, spanning kaon, charm, beauty, top/Higgs/electroweak, charged-lepton, neutrino-adjacent, and EDM-adjacent families, with one-line source hints and initial coverage annotations at lines 185-343.

4. FAIL - PKA is fully specified with model, timing, prompt, success criteria, and deliverables at lines 347-393, but WA, CA, DA, and Opus have tasks/checklists/responsibilities without explicit deliverable paths and role-specific success criteria at lines 395-448; writer/checker separation itself is present at lines 397, 410, and 456-458.

5. PASS - Section E gives an executable dispatch order from branch creation through PKA, WA, CA, discovery, Opus signoff, compilation, and merge policy, with a quantitative convergence rule at lines 452-470.

6. PASS - Section F defines the YES/PARTIAL/NO coverage approach with concrete grep commands, modules to inspect, and an implementation difficulty rubric, while Section C supplies initial coverage annotations for the full process list at lines 175-183, 187-343, and 476-509.

7. FAIL - Section G gives expected rounds and DA minimums at lines 512-517, but it does not state a hard max-iteration cap, and Section E leaves DA repetition open-ended when more than one new process appears at lines 459-460; arbitration, versioning, and sign-off log conventions are otherwise present at lines 519-538.

8. PASS - Section H provides agent-hour estimates, a risk register with mitigations, and PI-facing acceptance criteria at lines 542-573.

9. PASS - Section I explicitly asks for the `K1 -> pi gamma` clarification and for scope confirmation across quark, LFV, neutrino, electroweak, Higgs/top FCNC, and EDM-adjacent processes at lines 577-578.

10. PASS - The major cross-references and policies are internally consistent: planning does not create the catalog now, execution moves to `flavor-catalog/2026q2`, and agent writes are confined to `flavor_catalog/` at lines 5, 85-89, 357-358, and 452-463.

11. PASS - The plan substantially matches the PI's intent by keeping discovery separate from the existing repo surface, assigning PKAs only a few processes, separating WA and CA, marking code coverage, and requiring implementation difficulty at lines 78-83, 349, 397, 410, 157-161, and 494-509.

## Findings

### BLOCKER 1 - Section D - Role specifications are incomplete for non-PKA roles

Section D is strong for PKAs: it gives model choice, wall time, concurrency, an exact prompt, success criteria, and deliverables. The other roles are not specified to the same operational standard. WA has model, batch size, timing, and tasks, but no explicit deliverable list or success criteria. CA has a checklist and one worklog path, but no precise status-output contract, no failure artifact shape, and no success criteria beyond the checklist. DA has success criteria but no output path or version/update responsibility for proposed additions. Opus has responsibilities but no explicit signoff artifact, status mutation rule, or success/failure criteria beyond approving `OPUS-APPROVED`.

Fix: add short deliverables and success criteria blocks under WA, CA, DA, and Opus. For WA, specify updated `.tex`, updated `.yaml`, proposed `catalog.bib` changes, and `worklogs/writer/<batch_id>.md`. For CA, specify `worklogs/checker/<batch_id>.md`, pass/fail issue list, verified source/value table, and required sidecar status transition. For DA, specify `worklogs/discovery/round_<NNN>_<scope>.md` plus an additions/merge proposal. For Opus, specify `signoff/round_<NNN>_index.md` and `signoff/by_process/<process_id>.md` updates, with clear APPROVE / RETURN-TO-WA / ESCALATE-TO-PI outcomes. Keep the existing rule that writer and checker are separate.

### BLOCKER 2 - Section G - No hard max-iteration cap

The checklist explicitly requires a max-iteration cap. Section G states expected rounds, and Section E says to repeat discovery if more than one genuine new process appears, but neither section gives a hard ceiling or a mandatory escalation point. That leaves the orchestrator without a deterministic stop condition if WA/CA keeps cycling on the same disputed source or discovery keeps finding two marginal additions.

Fix: add a bounded cap such as: PKA source-fix micro-iterations max 2 before WA triage; WA/CA max 3 full cycles per batch before Opus arbitration; DA max 4 total rounds before PI scope decision, unless PI explicitly expands the search; Opus max 1 corrective loop before PI escalation. Tie every cap to a signoff log entry and a status like `BLOCKED-PI` or `DEFERRED-SCOPE` so no process remains silently open.

### WARNING - Sections B/E/G - Scalar status wording can become ambiguous

The YAML skeleton defines one scalar `status` with mutually exclusive values ending in `OPUS-APPROVED`, but the convergence rule says every process has `CHECKER-DONE` and `OPUS-APPROVED`. This is probably intended as a historical transition, not simultaneous scalar values, but execution agents could interpret it differently.

Fix: either clarify that `OPUS-APPROVED` implies prior checker pass, or add `status_history`, `checker_passed_at`, and `opus_approved_at` fields to the sidecar/index schema.

## Recommendation
Send back to planner with the BLOCKERs above.

===FLAVOR_CATALOG_PLAN_V0_REVIEW_END===
