1. NIT `flavor_catalog_constraints/physics_adapters/bsgamma.py:26`: new `exclusive_bsphigamma_*` wrappers are not listed in `__all__`, unlike B012 wrappers; fix by adding the three names. Non-blocking for B013’s direct imports.

2. No BLOCKER/SHOULD-FIX findings. Isolation OK: B013 imports physics only via `physics_adapters.bsgamma`; adapter diff is append-only; `quarkConstraints/bsgamma.py` diff is empty.

3. Contract/anchors OK: result numeric fields are real floats, complex Wilsons stay in diagnostics, severity is HARD, missing couplings degrade without crash. Missing/mismatched anchor probes raise `AnchorError`.

4. Independent numerical cross-check: manual BR `3.359433628084557e-05`, result BR `3.359433628084557e-05`, abs diff `0.000e+00`; ratio-to-SM `0.988068714142517`.

5. Wilson check: `C7_NP=(0.0018934910536189335+0.001262327369079289j)`, `C7p_NP=(0.0063116368453964445-0.0018934910536189335j)`; LL `u77=0.511997570352028`, `u78=0.119166114187617`, `u88=0.556684863172384`.

6. Pass/fail points: right `0.1` gives BR `3.546560100925992e-05`, budget ratio `0.286110528804824`, pass; right `0.2` gives BR `3.986240403703966e-05`, budget ratio `1.144442115219293`, fail.

7. Determinism OK: same point returned equal `ConstraintResult`; registry import smoke OK: `fcc.get("B013").process_id == "B013"`.

8. Pytest: focused B013 `12 passed`; full `python -m pytest tests/constraints/ -q` -> `1033 passed in 31.62s`.

CODE-OK