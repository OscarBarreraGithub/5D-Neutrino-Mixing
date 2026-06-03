# W2 PHASE 3 — SUB-STEP 3d-B: rare-B charged-lepton VECTOR rewire (codex, gpt-5.x xhigh). Repo: /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing. Implement the rare-B part of sub-step 3d of the DUAL-APPROVED plan `.orchestration/runs/W2-P3/plan.md` (step 10), using the 3a builder `quarkConstraints/rs_semileptonic_wilsons` (b_to_s_ll, b_to_d_ll → C9/C10/C9'/C10') and the 3b/3c rewire+degradation pattern (committed 651389e). First a SHORT plan, then implement. Codex + Opus dual-review (both must APPROVE). SCOPE = rare-B ONLY (rare-K = 3d-K, rare-charm = 3d-C are separate).

REWIRE these rare-B same-flavor charged-lepton VECTOR (b→sℓℓ / b→dℓℓ) consumers to consume `rs_semileptonic_wilsons.b_to_s_ll`/`b_to_d_ll` C9/C10/C9'/C10' additively at the rare-B core's Wilson-value consumption point (NOT via `_wilson_prefactor`, NO second 1/M_KK²):
- Across the 5 split adapter files: `rare_b_meson.py`, `rare_b_electronic.py`, `rare_b_tauonic.py`, `rare_b_kstar_dilepton.py`, `rare_b_baryon.py`.
- Constraints: B005, B006, B007, B008, B015, B016, B017, B018, B019, B021 (verify each against the catalog; rewire only the b→sℓℓ/b→dℓℓ VECTOR C9/C10 part).
- **B017/B018/B019 are LFU RATIO observables (R_K-type)**: preserve numerator/denominator semantics; a lepton-universal Wilson largely cancels in the ratio — keep the ratio structure intact.
- **B015 also carries a C7 dipole proxy**: replace ONLY C9/C10/C9'/C10' (vector/axial); leave C7 unchanged (dipole = Phase 7).
- SM rates / form factors / anchors / budgets UNCHANGED.

GRACEFUL DEGRADATION (as 3b/3c): rigorous `rs_semileptonic_wilsons` when present; ABSENT (old-style point) ⇒ non-vetoing `evaluated=False` + `missing_extra` (no crash, no fake pass). The 100M scan always provides it.

TESTS: for the rewired constraints — rigorous-path (point via `build_from_rs_ew_inputs`, NP C9/C10 enters; cross-check independently of the adapter), SM-limit (universal-c ⇒ NP Wilsons=0 ⇒ recover committed SM rate), absent-path (non-vetoing). For B017/B018/B019 assert the ratio still behaves (lepton-universal NP ~cancels). For B015 assert C7 path untouched. Replace any proxy-only tests with rigorous/SM/absent equivalents — enumerate, no silent coverage loss. `python -m pytest tests/ -q` stays green.

CONSTRAINTS: reach physics ONLY via adapters; numeric ConstraintResult fields real finite floats (complex Wilsons in diagnostics); touch ONLY the rare-B constraints/adapters/tests for the listed IDs.

OUTPUT (<=16 lines): short plan; rewired files + IDs; one rigorous-path C9/C10→BR example + SM-limit recovery + absent-path; LFU-ratio + B015-C7 confirmations; the test-count change (enumerated); pytest counts. End with: P3DB-DONE.
