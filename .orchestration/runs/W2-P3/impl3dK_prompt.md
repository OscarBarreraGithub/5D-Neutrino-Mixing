# W2 PHASE 3 — SUB-STEP 3d-K: rare-kaon charged-lepton VECTOR rewire (codex, gpt-5.x xhigh). Repo: /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing. Implement the rare-kaon part of sub-step 3d (plan `.orchestration/runs/W2-P3/plan.md` step 11), using the 3a builder `quarkConstraints/rs_semileptonic_wilsons.s_to_d_ll` (C9/C10/C9'/C10') and the 3d-B rewire+degradation pattern (committed 4834d50). First a SHORT plan, then implement. Codex + Opus dual-review (both must APPROVE). SCOPE = rare-KAON ONLY (rare-charm = 3d-C separate).

REWIRE these s→dℓℓ same-flavor charged-lepton consumers to consume `rs_semileptonic_wilsons.s_to_d_ll` C9/C10/C9'/C10', MAPPED into the existing rare-kaon core's effective inputs (the Y / y7V / y7A short-distance coefficients), additively — NOT via any old proxy prefactor, NO second 1/M_KK^2:
- K006, K008, K009, K010, K012 (verify each against the catalog; rewire only the s→dℓℓ VECTOR/axial short-distance part).
- K012 (K_S→μμ) is Im-sensitive — preserve its Im[...] structure (the NP enters via the same s→d Wilsons; keep the existing CP/Im treatment).
- K009/K010 (K_L→π0ℓℓ) are semileptonic with direct-CP y7V/y7A — map the s→d C9/C10 into y7V/y7A consistently; KEEP the long-distance / interference diagnostics where the constraint is already PARTIAL (do not claim full rigor on LD pieces).
- SM rates / form factors / anchors / budgets UNCHANGED.

GRACEFUL DEGRADATION (as 3d-B): rigorous `rs_semileptonic_wilsons` when present; ABSENT (old-style point) ⇒ non-vetoing `evaluated=False` + `missing_extra` (no crash, no fake pass).

TESTS: rigorous-path (point via `build_from_rs_ew_inputs`, NP s→d Wilsons enter the Y/y7 inputs; cross-check independently of the adapter), SM-limit (universal-c ⇒ NP=0 ⇒ recover committed SM rate), absent-path (non-vetoing). Preserve K012 Im-sensitivity test + K009/K010 LD-partial diagnostics. Replace proxy-only tests with rigorous/SM/absent equivalents — enumerate, no silent coverage loss. `python -m pytest tests/ -q` stays green.

CONSTRAINTS: physics ONLY via adapters; ConstraintResult numeric fields real finite floats (complex Wilsons in diagnostics); touch ONLY K006/K008/K009/K010/K012 + their adapters/tests.

OUTPUT (<=14 lines): short plan; rewired files + IDs; one rigorous s→d C9/C10→rate example + SM-limit recovery + absent-path; confirm K012 Im-structure + K009/K010 LD-partial preserved; test-count change (enumerated); pytest counts. End with: P3DK-DONE.
