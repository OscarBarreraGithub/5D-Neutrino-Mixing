# W7 plan REVIEW — independent critique

Independently critique the plan at `.orchestration/runs/W7-MUEGAMMA/plan_codex.md`. Read the
requirements `.orchestration/runs/W7-MUEGAMMA/plan_prompt.md` and verify against the REAL repo.
Do NOT rubber-stamp. The μ→eγ model is LOCKED = LMFV (Perez-Randall, C=0.02); do not relitigate it.

Check specifically:
1. **Dipole reuse:** does the plan correctly reuse `flavorConstraints/muToEGamma.py` (LMFV spurion
   (Y_N Y_N†)_12, C=0.02, br_limit=1.5e-13) rather than re-deriving it?
2. **Lepton ParameterPoint extra:** are the fields (Y_N, PMNS, M_KK/M_N, c_L,c_E,c_N, v,k,ε)
   complete and correctly produced from a draw? Does it compose with the existing quark flow via
   the same get_extra pattern EW001 uses, WITHOUT breaking quark-only scans (byte-identical)?
3. **Constraint wiring:** does the catalog μ→eγ constraint (locate the real L00x file) read the
   extra and call the dipole correctly, with rigorous|proxy|partial tag + graceful degradation
   when the lepton extra is absent?
4. **Tests:** numeric BR oracle, quark-only byte-identity regression, absent-extra degradation,
   finite/perturbativity guards — all specified?
5. **Isolation/honesty:** files-touched list accurate; overlap with W8/W9 flagged; deferred
   loop-level / off-diagonal LFV honestly stated.

End with EXACTLY one line: `VERDICT: APPROVE` or `VERDICT: NEEDS-FIXES`, then numbered fixes.
