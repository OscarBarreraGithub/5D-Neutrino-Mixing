1. BLOCKER: none.
2. SHOULD-FIX: none.
3. NIT: none.
4. X_NP verdict: physical/correct, not a normalization bug. Builder maps `X_NP=C/g_SM^2`; K core adds it to `lambda_t * X_t`, not to `X_t`.
5. Independent K004 recompute: `C_LL=1.461e-10-2.017e-11j`, `C_RL=7.071e-12`, `g_SM^2=1.781e-7`, so `X_NP=8.6009e-4-1.1327e-4j`.
6. Trace: `lambda_t X_t=-4.812e-4+2.102e-4j`; adding `X_NP` gives `3.789e-4+9.697e-5j`, BR `8.4726e-11 -> 4.9286e-12`. Naive `lambda_t*(X_t+X_NP)` gives `8.4794e-11`.
7. B022 recompute: `X_NP_total=-7.4740e-2-4.11e-6j`, `R_K=5.0167`, BR `5.58e-6 -> 2.5547e-5`; B023 `R_K*=1.7293`.
8. SM-limit: universal-c gives `X_NP=0`; K004 recovers `8.472598449611133e-11`, K005 `2.95375989343059e-11`, B022/B023 recover their SM predictions.
9. SM side unchanged: `P_c=0.404`, `delta_EM=0.997`, B+ LD term retained, NA62/KOTO/HPQCD/Buras anchors still sidecar-driven.
10. Degradation: missing/invalid `rs_semileptonic_wilsons` is non-vetoing with `evaluated=False`; Majorana/Dirac light-nunu rates match.
11. Test counts: focused B022/B023/K004/K005 `47 passed`; full `pytest tests/ -q` `1664 passed, 1 skipped`.
12. P4D-OK