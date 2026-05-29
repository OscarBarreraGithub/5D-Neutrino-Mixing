1. No BLOCKER / SHOULD-FIX / NIT findings. Isolation OK: `B004.py` imports no other constraint; `deltaf2.py` diff is append-only imports/`__all__`/new `bs_mixing_m12_np_from_wilsons_with_running`, no existing function bodies changed.

2. Contract OK: numeric `ConstraintResult` fields are real floats; complex amplitudes stay in diagnostics; missing `quark_mass_basis_couplings` returns a non-crashing pass result with `predicted=None`.

3. Anchor OK: values load via `load_anchor`; probes fail loudly:
missing key -> `AnchorError`; bad units `degrees` -> `AnchorError`.

4. Independent numerical check from `quarkConstraints.deltaf2._evolve_wilsons + compute_m12_np`:
anchors `exp=-0.041`, `sm=-0.0368`, budgets `+0.016025292509/-0.0160112460477`, `M12_SM=5.85e-12 GeV`.

5. Safe point `(1e-3+1e-3j, right=0)`: direct `M12_NP=0+4.458438064963e-15j`, `phi_np=7.621260173942e-04`, `pred=-3.603787398261e-02`, `ratio=0.3096433974355`, pass=True; constraint matches.

6. Excluded point `(1e-2+1e-2j, right=0)`: direct `M12_NP=0+4.458438064963e-13j`, `phi_np=7.606557176438e-02`, `pred=3.926557176438e-02`, `ratio=5.008680604060`, pass=False; constraint matches.

7. Tests: `python -m pytest tests/constraints/primary/beauty/test_B004.py -q` -> `10 passed`; `python -m pytest tests/constraints/ -q` -> `125 passed`.

8. Registry smoke: `registered=13`, `has_B004=True`, `family=beauty`, `import_failures=0`.

CODE-OK