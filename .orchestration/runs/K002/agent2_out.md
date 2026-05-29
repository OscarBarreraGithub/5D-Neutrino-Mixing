1. NIT: Correct amplitude. K002 uses `abs_m12_np` for the verdict, not `Im(M12)`: `K002.py:172-188`, core `quarkConstraints/deltaf2.py:817-825`. Correct for a mass splitting.

2. NIT: QCD running is actually used. Verdict path is `K002.py:172-175` -> adapter `deltaf2.py:210` -> core `evaluate_delta_mk_with_running` at `quarkConstraints/deltaf2.py:863-886`. Test benchmark: unrun ratio `0.00875337`, run ratio `0.0934376`, enhancement `x10.674`.

3. NIT: Budget is right. YAML `K002.yaml:85-92` gives `0.5293 +/- 0.0009` in `10^10 hbar s^-1`; conversion gives `Delta m_K = 3.483915888e-15 GeV`, budget `Delta m_K/2 = 1.741957944e-15 GeV`. This matches core `3.484e-15 GeV` at `quarkConstraints/deltaf2.py:638` within `0.014 sigma`; not too tight/loose.

4. NIT: Long-distance treatment is correct. No SM short-distance subtraction is applied: `K002.py:192-207`; this matches the convention in `docs/quark_scan_assumptions_compact.tex:464-466`.

5. NIT: Anchor/snapshot and units check out. YAML `pdg_fit_assuming_cpt` matches snapshot `pdg_live_S013D_20260516.txt:10-13`; code converts with `hbar` in GeV seconds at `K002.py:43-44,72-113`. `Severity.HARD` at `K002.py:149` is appropriate for this veto.

PHYSICS-OK