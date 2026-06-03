# W2 PHASE 3 — SUB-STEP 3d-C: rare-charm charged-lepton VECTOR rewire (codex, gpt-5.x xhigh). Repo: /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing. Implement the rare-charm part of sub-step 3d (plan `.orchestration/runs/W2-P3/plan.md` step 12), using the 3a builder `quarkConstraints/rs_semileptonic_wilsons.c_to_u_ll` (C9/C10/C9'/C10') and the 3d-B/3d-K rewire+degradation pattern (committed f04ae1c). First a SHORT plan, then implement. Codex + Opus dual-review (both must APPROVE). SCOPE = rare-CHARM ONLY; this is the LAST sub-step of Phase 3.

REWIRE these c→uℓℓ same-flavor charged-lepton consumers to consume `rs_semileptonic_wilsons.c_to_u_ll` C9/C10/C9'/C10' additively at the rare-charm core's Wilson-value consumption point (NOT via the old proxy prefactor, NO second 1/M_KK^2):
- C004, C005, C007 (verify each against the catalog; rewire only the c→uℓℓ VECTOR/axial part).
- **C006, C008 are LFV (e±μ∓)** — DO NOT touch them; they wait for Phase 4 (lepton NC). Confirm they are NOT in scope.
- Charm is theoretically dirty (large long-distance): KEEP any existing long-distance / full-q² PARTIAL / NEEDS-HUMAN diagnostics where the constraint is already partial (e.g. C007 full-q² proxy + form-factor uncertainty). Only the SHORT-DISTANCE c→u C9/C10 NP becomes rigorous; LD stays as-is.
- SM rates / form factors / anchors / budgets UNCHANGED.

GRACEFUL DEGRADATION (as 3d-B/3d-K): rigorous `rs_semileptonic_wilsons` when present; ABSENT (old-style point) ⇒ non-vetoing `evaluated=False` + `missing_extra` (no crash, no fake pass).

TESTS: rigorous-path (point via `build_from_rs_ew_inputs`, NP c→u Wilsons enter; cross-check independently of the adapter), SM-limit (universal-c ⇒ NP=0 ⇒ recover committed SM rate), absent-path (non-vetoing). Preserve C007 LD/full-q² PARTIAL diagnostics. Replace proxy-only tests with rigorous/SM/absent equivalents — enumerate, no silent coverage loss. `python -m pytest tests/ -q` stays green.

CONSTRAINTS: physics ONLY via adapters; ConstraintResult numeric fields real finite floats (complex Wilsons in diagnostics); touch ONLY C004/C005/C007 + their adapter(s)/tests; do NOT touch C006/C008.

OUTPUT (<=14 lines): short plan; rewired files + IDs; one rigorous c→u C9/C10→rate example + SM-limit recovery + absent-path; confirm C006/C008 untouched + C007 LD-partial preserved; test-count change (enumerated); pytest counts. End with: P3DC-DONE.
