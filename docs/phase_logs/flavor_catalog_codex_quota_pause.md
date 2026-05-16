# Flavor Catalog — Codex Quota Pause (resumption plan)

**Date**: 2026-05-16, ~14:51 EDT.
**Status**: codex CLI returned "You've hit your usage limit. Visit
https://chatgpt.com/codex/settings/usage to purchase more credits or try again
at 3:08 PM." Quota resets at **3:08 PM local (≈ 18 minutes from pause)**.

## Catalog state at pause

- **Branch**: `flavor-catalog/2026q2` (current).
- **Process drafts**: 67 across 6 families.
  - kaon: 9 (K001-K006, K008, K009, K010, K013, K017)
  - charm: 7 (C001-C008)
  - beauty: 14 (B002, B004, B005, B006, B009, B011, B015, B017-B019, B021-B026, B032-B034)
  - top_higgs_ew: 11 (T001-T019, EW001-EW003)
  - charged_lepton: 8 (L001-L010)
  - edm_neutrino: 5 (E001, E002, E004, E006-E008)
- **CHECKER-DONE**: 50 (47 Wave-1-4 + T015 + C005 + B034 from Wave-5a).
- **WRITER-REWORK queued** (Wave-5a, 9 processes): T005, T019, L005, L006,
  E002, K008, B021, B022, B023.
- **Wave-5b WA progress**: kaon done, top_higgs_ew done, charm_edm still
  running on existing codex session (may finish, may fail mid-task at quota
  limit). Need CA pass after WA lands.
- **Opus round-1 sign-off** on the 50 CHECKER-DONE: running (uses Anthropic
  API, not affected by codex quota).

## Queued work waiting on codex (do NOT dispatch until probe says OK)

1. **WA-v2 cycle-2 rework for Wave-5a (4 batches, 9 processes)**:
   - wa_w5a_top_higgs_ew_v2: T005, T019 (CHK-1 metadata pattern)
   - wa_w5a_charged_lepton_edm_v2: L005, L006, E002 (CHK-1 metadata pattern)
   - wa_w5a_kaon_charm_v2: K008 (CHK-1; C005 already PASS)
   - wa_w5a_beauty_v2: B021, B022, B023 (CHK-1, plus B021 CHK-2)
2. **CA-v2 on the above** (4 batches) once WA-v2 lands.
3. **Wave-5b CA pass** (3 batches) once Wave-5b WAs all land:
   - ca_w5b_top_higgs_ew: T006, T016, T017
   - ca_w5b_kaon: K009, K010
   - ca_w5b_charm_edm: C006, C008, E007
4. **Any Wave-5b WA-v2 rework** (likely needed; metadata pattern is common).
5. **DA-3 discovery round** (optional, ~3rd of plan's 4-round budget).
6. **Wave-6 PKAs** if DA-3 surfaces more (e.g. L020-L025, K012/K018, B029-B031).
7. **Second Opus round** if Wave-5/6 substantially grows the CHECKER-DONE set.

## Estimated codex budget for resumption

- 4 WA-v2 (~12 min each, gpt-5.5 xhigh) ≈ ~50 min serial → dispatch in
  parallel batches of 3–4. ≈ 1 wall-hour.
- 4 CA-v2 ≈ same. ≈ 1 wall-hour.
- 3 CA-w5b ≈ ~30 min.
- DA-3 ≈ ~15 min.
- Wave-6 PKAs (if dispatched, 8–12) ≈ ~1 wall-hour.
- **Total**: ~3–4 wall-hours of intermittent codex calls. Watch quota carefully.

## Pre-flight usage probe

Created `~/bin/codex_check_usage.sh`. Returns exit code 1 if quota limited.

```bash
~/bin/codex_check_usage.sh && echo OK || echo BLOCKED
```

The orchestrator should run this probe BEFORE dispatching any new codex batch:
- If OK: dispatch a wave of up to 4 codex jobs in parallel, then sleep
  until completion notifications.
- If BLOCKED: do not dispatch. Re-probe on next user prompt OR after the
  documented reset window.

## Master compile + handoff readiness

Independent of codex quota:
- Opus round-1 sign-off in flight (uses Anthropic API).
- After it lands, `flavor_catalog/catalog_master.tex` can be compiled with
  `pdflatex` to produce a draft `master_catalog.pdf` covering the 50
  CHECKER-DONE + OPUS-APPROVED entries.
- The 17 Wave-5 non-yet-Opus-approved drafts can be marked DRAFT in the
  master TOC and listed in an appendix for the PI to review.

## Resumption checklist (when codex is back at 3:08 PM)

1. Run `~/bin/codex_check_usage.sh` → confirm exit 0.
2. Dispatch the 4 WA-v2 Wave-5a rework batches in parallel.
3. After WA-v2 lands, dispatch CA-v2 on the same batches.
4. Dispatch CA-w5b once Wave-5b WAs are all in.
5. Handle any cycle-2 reworks per the established pattern.
6. Dispatch DA-3.
7. Dispatch Wave-6 if DA-3 proposes additions.
8. Dispatch Opus round-2 sign-off on the now-larger CHECKER-DONE set.
9. Rebuild master_catalog.pdf.
10. Tag `flavor-catalog-v0` and push.
