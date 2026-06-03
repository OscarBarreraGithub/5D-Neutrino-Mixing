1. BLOCKER: none.
2. SHOULD-FIX: none.
3. NIT: isolation has 2 extra untracked `.orchestration/runs/W2-P5/impl5b_*` notes; code/test touch set is otherwise exactly EW002/K018/K017/B009/B025 + adapter/helper/tests.
4. Rewire OK: formulas consume `rs_charged_current.epsilon`; no old `m^2/M_KK^2` proxy or second scaling found.
5. `delta_G_F/G_F`: subtracted once in Phase-5a epsilon; Phase-5b only reads epsilon, no double subtract.
6. Recompute K017: `R_K^SM=2.477e-05`, eps_e=0.1/eps_mu=0 gives `2.99717e-05 = 1.21 * R_K^SM`.
7. Recompute B025: `R_D^SM=0.296`, eps_tau=0.5/light=0 gives `0.666 = 2.25 * R_D^SM`.
8. SM-limit universal-c: EW002 shift `~0`, K018 `0.22330377397401527`, K017/B009/B025 equal SM.
9. Degradation: all five absent `rs_charged_current` return `passes=True`, `predicted=None`, `evaluated=False`, `missing_extra`.
10. v1 near-SM: EW002 `5.25e-05`, K018 `5.90e-06`, K017/B025 LFU cancel, B009 third-family `-3.06e-07` BR shift from eps_ub_tau `~-1.77e-3`.
11. Labels preserved: EW002 remains SOFT/non-veto; B025 diagnostics keep `matching_coverage=PARTIAL`.
12. Test count: touched tests old 43, new 45; rigorous/SM/absent coverage present.
13. `pytest tests/ -q`: `1673 passed, 1 skipped in 797.22s`.
14. P5B-OK