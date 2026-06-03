# W2 PHASE 4 — SUB-STEP 4c-KC: LFV rare K + charm rewire (codex, gpt-5.x xhigh). Repo: /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing. Implement the LFV-rare-K + LFV-rare-charm part of sub-step 4c (plan `.orchestration/runs/W2-P4/plan.md` step 13), using the 4a lepton builder (committed a585265: `rs_semileptonic_wilsons.lfv_llqq` off-diagonal-lepton Wilsons) and the established rewire+degradation pattern. First a SHORT plan, then implement. Codex + Opus dual-review (both must APPROVE). SCOPE = LFV rare K + charm ONLY (the LFV leptonic L00x = 4c-L is separate).

REWIRE (LFV semileptonic, off-diagonal e-μ lepton):
- **K019** (K_L→eμ), **K020** (K⁺→π⁺eμ), **K021** (K_L→π⁰eμ) — under secondary/kaon. **C006, C008** (charm c→u eμ LFV).
- These currently supply the off-diagonal e-μ lepton coupling as a DOCUMENTED PROXY (flagged NEEDS-HUMAN "lepton coupling not on ParameterPoint"). Replace that proxy with the RIGOROUS `rs_semileptonic_wilsons.lfv_llqq` (the s→d / c→u quark × off-diagonal e-μ lepton contact), feeding the existing rare cores' Y/y7 or C9/C10 inputs additively (no _wilson_prefactor, no second 1/M_KK^2).
- **v1 PHYSICS**: with the diagonal charged-lepton fit (U_e=I, universal c_L), the off-diagonal e-μ lepton coupling = 0 ⇒ tree-level LFV rate = 0 (rigorous-tree, non-vetoing). This RESOLVES the old "lepton coupling not on ParameterPoint" NEEDS-HUMAN (it IS now on the point, rigorously, and =0 for this fit). Replace that flag with an honest note: "tree-level LFV rigorous (=0 for the diagonal charged-lepton fit); nonzero only with non-diagonal lepton structure / loop-induced — deferred." Do NOT fake a nonzero bound.
- SM (LFV ⇒ SM≈0) / anchors / experimental limits UNCHANGED.

GRACEFUL DEGRADATION: rigorous `rs_semileptonic_wilsons` when present; ABSENT ⇒ non-vetoing evaluated=False + missing_extra (no crash/fake pass).

TESTS: v1-path (diagonal fit ⇒ LFV rate=0, rigorous-tree non-vetoing); LFV-LIVE (non-diagonal-U_e/non-universal-c toy ⇒ nonzero LFV rate, cross-check independent of adapter, confirm it bites + scales); absent-path. Replace proxy-only tests with these — enumerate, no silent coverage loss. `python -m pytest tests/ -q` stays green.

CONSTRAINTS: physics ONLY via adapters; ConstraintResult numeric fields real finite floats (complex in diagnostics); touch ONLY K019/K020/K021 + C006/C008 + their adapters/tests.

OUTPUT (<=14 lines): short plan; rewired files + IDs; v1-zero (rigorous-tree non-vetoing) + LFV-live-toy nonzero rate + absent; confirm old lepton-coupling NEEDS-HUMAN resolved (replaced by honest note); test-count change; pytest counts. End with: P4CKC-DONE.
