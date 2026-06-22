# QUARK-FIX-2026-06-21 — orchestration status

Gate: dual-signoff. Codex (gpt-5.x) is rate-limited until 2026-06-24, so per user
decision (2026-06-21) the plan + impl reviews run with INDEPENDENT OPUS agents now,
with a **Codex cross-check reserved as a final gate BEFORE the commit is pushed**
(after 2026-06-24). Do not push until that Codex pass is done.

- Plan author: Opus agent acf3a8d85e9497bbc -> PLAN.md (v1)
- Plan review v1: Opus (independent) -> review_plan_opus_v1.md  [in progress]
- Codex plan review: BLOCKED (usage limit until 2026-06-24) -> rerun then.

## PLAN GATE CLOSED (2026-06-22)
PLAN v3 dual-approved at Opus level:
- author: acf3a8d85e9497bbc (v1) / a8da73604e68012d6 (v2) / a81a5d38d67a4dea8 (v3)
- reviewer (independent): acc8fb3814bc78c23 -> APPROVE v1; a50a86b4aa117677d -> REVISE v2
  (caught §8.2 floor-extraction rescale BLOCKER); afe97e23903e91b8b -> APPROVE v3.
- M1 gate decision: T010 R_b veto -> T011-style NP-shift budget (source-defensible). APPROVED.
- Scan recipe: ALL headline floors from full §8.1 paired-1M grid (~30min wall / ~115 core-hr).
  §8.2 cheap shortcut DROPPED as floor source (invalid for CR0xx/EW001 scaling).
- Codex plan review STILL OWED (usage limit) -> run as backup before push.

## PHASE 2 (impl) sequence: code+tests green -> Opus impl review -> run scan -> restate -> final.

## IMPL GATE CLOSED + SCAN DONE (2026-06-22)
- Impl dual-approved (Opus): reviewer a74b385559c7d46cb APPROVE; suite 1767 pass.
- Robust 100k MINIMAL quark-only scan (job 24027419) with fixes: floors EW001 20 / epsK 7 /
  Zbb 5 TeV; combined 20 TeV (proxy S,T,U); Zbb NOT leading. See FINAL_REVIEW.md.
- mu->e gamma from full-catalog lepton-viable subset (job 24025068, 4977 draws).
- Interactive notebook delivered. NOTHING committed/pushed. Codex backup still owed pre-push.
