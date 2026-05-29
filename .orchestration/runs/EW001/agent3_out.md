1. BLOCKER: none.
2. SHOULD-FIX: none.
3. NIT: none.
4. Isolation OK: EW001 adds only new files; no pre-existing tracked function bodies modified. Adapter boundary OK at `EW001.py:49` and `physics_adapters/oblique_stu.py:5`; no other constraint imports.
5. Contract/anchor/determinism OK: real float result fields, HARD severity, missing/invalid mass non-crashing, load_anchor route verified, mismatch/missing probes raise `AnchorError`.
6. Cross-check 6 TeV: S=0.0505200967051021, T=0.120428357556969, chi2=3.63125879555267, manual=3.63125879555267, ratio=0.606071982401272, pass=True.
7. Cross-check 3 TeV: S=0.202080386820408, T=0.481713430227877, chi2=110.844082202322, manual=110.844082202322, ratio=18.5003318188414, pass=False. SM chi2=0.962793485090136.
8. Tests: EW001 targeted `10 passed`; `python -m pytest tests/constraints/ -q` -> `862 passed in 17.83s`; registry smoke `fcc.get("EW001")` -> `top_higgs_ew PRIMARY HARD`.

CODE-OK