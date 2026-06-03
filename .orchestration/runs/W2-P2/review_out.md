1. `BLOCKER quarkConstraints/rs_ew_spectrum.py:207`: `x1` is OK (`2.450509663813735`, `m=7351.528991441206 GeV`), but the multi-root tower is wrong: expected n6 `18.120466878760`, module n6 `21.261294197360`; n10 regresses/duplicates `27.543551802036` after n9 `30.684853611711`.
2. This contaminates `a(c)`: independent correct-root Q=1024 gives `a(0.2)=21.8605473380645`, `a(0.65)=-1.45713993302811`; module/tests give `21.9165452359042`, `-1.47609745602508`.
3. `SHOULD-FIX tests/test_rs_ew_spectrum.py:69`: overlap “manual” test calls through `spectrum.omega()`/`chi_n()` and lacks sorted/unique root assertions, so it locks in the bad tower.
4. Truncation: correct-root partials pass `N=16→32` (`rel=3.651731e-05` for c=0.2, `2.518757e-04` for c=0.65); module doubling is implemented and `rel_tol=0` raises at `N=512`.
5. Other probes pass: endpoint max residual vs `wavefuncs.py` `2.220e-15`; universal `a-a_ref=0`; `a(0.2)>a(0.65)`; `(mZ²/MKK²)Δa≈3.59e-3`; cached/uncached identical; `__init__.py:42` lazy export safe.
6. Blast radius/tests: no existing constraint/scaffold/core references touched beyond the lazy export; `python -m pytest tests/ -q -p no:cacheprovider` -> `1644 passed, 1 skipped`.

PHASE2-NEEDS-FIXES