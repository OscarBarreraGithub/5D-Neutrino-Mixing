1. SHOULD-FIX: Low-q2 “validation” is not anchored to B016.yaml: q2=1.1-6.0 is hardcoded and the snapshot has no bin value; test only rechecks the same formula. Correct physics needs a sourced partial/differential BR or SM benchmark. `B016.py:66-68,351-366`; `test_B016.py:262-295`; value used `1.85195e-7`.

2. SHOULD-FIX: Budget is too tight for a C9/C10-only, no-covariance total-BR proxy. Exp is `(5.76±0.40)e-7`, formula SM is `5.75019e-7`, so budget is only `4.098e-8` (7.1%). Omitted C7 alone is plausibly O(10%) in rate; nonlocal charm/form factors add more. `B016.py:373-394`; limitations admitted at `rare_b_dilepton.py:691-694`.

3. SHOULD-FIX: Neutral mode is catalogued but not constrained; the same SM machinery predicts `5.3239e-7` vs HFLAV `3.28±0.32e-7`, a `6.39σ` tension. Either split charged/neutral scope or include a neutral/covariance treatment. `B016.py:15-17,431,439-442,581-589`.

4. NIT: Mandatory Δm/CP amplitude check is N/A for B016. For this rate, the implemented P→P combinations are physically right for the stated C9/C10-only approximation: `C9+C9'`, `C10+C10'`, mod-squared, no Re/Im misuse. `rare_b_dilepton.py:1025-1034`.

5. NIT: ΔF=2 `*_with_running(mu_had=2 GeV)` check is N/A. B016 does not use deltaf2; C9/C10 semileptonic-current QCD running is not the same requirement, though C7/charm evolution is omitted and marked. `rare_b_dilepton.py:1014-1028,1091`.

6. NIT: HFLAV anchors and units load correctly: charged `5.76±0.40e-7`, neutral `3.28±0.32e-7`; BR is dimensionless and GeV units in the rate are consistent. `B016.yaml:91-107`; snapshots lines `hflav_*:11-12`.

PHYSICS-NEEDS-FIXES