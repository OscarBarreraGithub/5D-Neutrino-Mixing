1. BLOCKER: none.
2. SHOULD-FIX: none.
3. NIT: none.
4. Isolation: `T015.py:43-52` imports only scaffold/base/registry and `physics_adapters.zpole_lfv`; no constraint-to-constraint import. No diff in pre-existing T010/T012/zpole/deltaf2 files.
5. Contract: `T015.py:318-341` missing-coupling path is non-crashing; evaluated numeric fields are real `float`; complex `delta_g_*` stay in diagnostics; severity is HARD.
6. Anchor: YAML CMS limit loaded as `1.9000000000000001e-07`; mismatch and missing probes raise `AnchorError` loudly.
7. Cross-check: left=`0.25+0.10j`, right=`0.05j` -> predicted=`3.508394829491628e-08`, manual=`3.508394829491628e-08`, ratio=`0.18465235944692779`.
8. Couplings: deltaL=`0.00023097717760444444+9.2390871041777783e-05j`, deltaR=`0+4.6195435520888892e-05j`, norm=`6.4020547888938088e-08`, SM total=`3.6495633333333339`.
9. Safe/fail: left=`0.10` -> BR=`4.6778599148909267e-09`, ratio=`0.024620315341531192`, pass; left=`10.0` -> BR=`4.6776411232732079e-05`, ratio=`246.19163806701093`, fail.
10. Determinism: repeated `evaluate()` results equal; point extra unchanged.
11. Pytest: targeted T015 `12 passed`; full `python -m pytest tests/constraints/ -q` -> `352 passed in 5.20s`; registry smoke ok.
CODE-OK