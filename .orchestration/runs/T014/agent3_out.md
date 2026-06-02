1. BLOCKER `flavor_catalog_constraints/physics_adapters/zpole.py:242` / `:355`: FCNC width/proxy physics is implemented inside the adapter; `quarkConstraints/zpole.py` has no append-only off-diagonal down-sector helper. Fix: move/add the core helper/dataclasses to `quarkConstraints/zpole.py` and leave the adapter as thin wrappers.

2. Isolation: T014 imports no peer constraint and reaches physics only via `flavor_catalog_constraints.physics_adapters.zpole`; AST check found no modified/removed pre-existing adapter function bodies, only added functions.

3. Contract/anchor: severity HARD; missing/invalid couplings return non-crashing unevaluated results; numeric result fields are real floats, complex only in diagnostics; anchors route through `load_anchor` and mismatch/missing probes raise.

4. Independent cross-check: limits `{bs,bd,sd}=0.0029`; SM hadronic width weight `2.5225098333333333`. For `bs_left=0.25+0.10j`, `bs_right=0.05j`: predicted/manual BR `1.5227818887525844e-07`; ratio `5.250972030181326e-05`; `delta_L=(0.00023097717760444444+9.239087104177778e-05j)`, `delta_R=4.619543552088889e-05j`; norm `6.402054788893809e-08`; FCNC width `3.8412328733362853e-07`.

5. Safe/excluded: safe point passed with ratio `0.0007001282891148734`; excluded point failed with ratio `6.861973237606202`; deterministic equality check passed.

6. Tests: targeted T014 `11 passed in 5.86s`; full `python -m pytest tests/constraints/ -q` -> `1002 passed in 26.01s`; registry smoke `T014 SECONDARY top_higgs_ew`.

CODE-NEEDS-FIXES