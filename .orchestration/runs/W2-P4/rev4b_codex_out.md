1. NIT flavor_catalog_constraints/primary/top_higgs_ew/T015.py:441, T016.py:453, T017.py:463: notes say `z_delta_g_e[...]`; runtime/diagnostics correctly use `z_delta_g_L/R_e`.
2. BLOCKER/SHOULD-FIX: none found. Rewire uses adapter `physics_adapters/zpole_lfv.py:93` and indices T015 `[0,1]`, T016 `[0,2]`, T017 `[1,2]`.
3. v1-zero verdict: CORRECT. Diagonal fit gives offdiag `(0j,0j)`, BR_NP `0.0`, ratio `0.0`, `evaluated=True`, passes non-vetoing; PARTIAL/NEEDS-HUMAN and Phase-7 loop deferral retained.
4. Live recompute: T015 `2.523682343e-9`, T016 `4.328828719e-9`, T017 `8.724492866e-9`; manual BR formula matched exactly.
5. Scaling: norm high/low `0.0625`; BR high/low `0.0625000001-0.0625000005`, consistent with `(m_Z/M_KK)^4`.
6. Degradation: old-style lepton-only point returns `evaluated=False`, `missing_extra=rs_ew_couplings`, `predicted=None`, `ratio=None`, `passes=True`.
7. Coverage: targeted collect `38->30`; proxy-only tests replaced by v1-zero, live/scaling, absent old-style, invalid extra, amplified fail, finite fields, deterministic.
8. Isolation/anchors: source diff limited to adapter + T015/T016/T017 + their tests + new helper; no sidecar/anchor/limit changes.
9. Pytest: `pytest tests/ -q` => `1645 passed, 1 skipped in 784.24s`; collect-only full count was `1646` because the skipped test is collected.

P4B-OK