1. SHOULD-FIX `tests/constraints/primary/charged_lepton/test_L007.py:195` - numeric test is independent but formula-level; it does not directly use the underlying `flavorConstraints.muToEGamma.check_mu_to_e_gamma_raw` core as requested. Fix: compute expected BR/ratio from `check_mu_to_e_gamma_raw(y, pmns[[1,2,0], :], ...)`.

2. Isolation/contract otherwise OK: `L007.py` reaches physics only via `physics_adapters.lepton_tau_mu`; no L001/lepton/core diffs; missing input returns `passes=True`, `predicted=None`, `ratio=None`, `sm_prediction=None`, `evaluated=False`; severity is HARD.

3. Anchor OK: YAML limit is `4.2e-8`; missing candidate raises `AnchorError`; mismatch probe raised `AnchorError: L007: primary_current_limit and primary_experiment differ`.

4. Direct core cross-check I ran: y=(0.10,0.20,0.30), core_lhs=0.021036774620197408, expected BR=1.7701835456839274e-11, constraint BR=same, ratio=0.00042147227278188746, pass=True.

5. Safe/excluded: safe BR=1.7701835456839277e-15 ratio=4.2147227278188757e-08 pass=True; excluded BR=4.5316698769508572e-07 ratio=10.789690183216328 pass=False.

6. Tests: `python -m pytest tests/constraints/ -q` -> 862 passed in 16.62s. Registry smoke: registry_count=84, has_L007=True, import_failures=0.

CODE-NEEDS-FIXES