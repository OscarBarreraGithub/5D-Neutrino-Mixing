1. BLOCKER: none. Mapping verified via `rare_charm_dilepton.py:159`, `rare_charm_semileptonic.py:127`, consumed in `C004.py:395`, `C005.py:441`, `C007.py:648`.
2. SHOULD-FIX: none. Core consumption is additive: D0 uses `C10-C10'` at `rare_charm_dilepton.py:546`; D+ uses `C9+C9'`, `C10+C10'` at `rare_charm_semileptonic.py:334`.
3. NIT: none. No old proxy prefactor or second `1/M_KK^2`; diagnostics explicitly report false.
4. Recomputed C007: `C9+C9'=27.378192871893+68.049753521990j`, `C10+C10'=-364.459436526798-905.880637939158j`, `BR=3.232711560514e-9`, ratio `0.04825` vs `6.7e-8`.
5. Charm-Wilson sanity: OK; `lambda_cu=5.589e-05-1.421e-04j`, so large Wilsons are expected, and 3a `C10=-20.4408-31.3859j` is reproduced, not a bug.
6. SM-limit: OK; universal-c point gives `C004=0`, `C005=0`, `C007=1e-12`.
7. Degradation: OK; old-style quark-coupling point returns `passes=True`, `evaluated=False`, `missing_extra=rs_semileptonic_wilsons`.
8. C006/C008: untouched; no diff in LFV constraints/adapters/tests.
9. C007 partial diagnostics preserved: full-q2 proxy true, LD not subtracted, no resonance/window acceptance, form-factor diagnostics retained.
10. Isolation/numerics: only the 2 adapters, C004/C005/C007, and their tests/helper changed; physics enters through adapters; `ConstraintResult` numeric fields are real finite floats.
11. Test counts unchanged: C004 `10`, C005 `12`, C007 `11`; proxy tests replaced by rigorous core-consumption/SM/absent/no-second-MKK coverage.
12. `pytest tests/ -q`: `1646 passed, 1 skipped in 783.60s`.

P3DC-OK