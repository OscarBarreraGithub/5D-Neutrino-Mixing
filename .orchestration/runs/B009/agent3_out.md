1. No BLOCKER/SHOULD-FIX/NIT findings for B009 code/tests.
2. Isolation: B009 imports no other constraint and no direct `quarkConstraints`; it reaches physics via `physics_adapters.leptonic_tree`.
3. Contract/anchor: numeric result fields are real floats; complex amplitude stays in diagnostics; missing `kk_ew_mass_gev` returns non-crashing pass; missing-anchor probe raises `AnchorError`.
4. YAML numbers: exp=`1.120e-4`, sigma=`1.900e-5`, SM anchor=`8.650e-5`, budget=`4.493733520830467e-5`.
5. Cross-check: `f_B=0.190`, `|Vub|=0.00368`; manual/core/constraint SM all `8.630796690372665e-5`.
6. Proxy check: `MKK=3000` pred=`8.630850146600551e-5`, ratio=`1.1895727158213082e-5`, pass=True.
7. Excluded check: `MKK=10` pred=`1.4112304413683584e-4`, ratio=`1.2198114770049617`, pass=False.
8. Determinism: repeated `evaluate(500 GeV)` equal; point unchanged; pred=`8.632721218869506e-5`, ratio=`4.28269386228468e-4`.
9. Tests: B009 targeted `9 passed`; full `python -m pytest tests/constraints/ -q` => `509 passed in 15.14s`.
10. Registry smoke: `fcc.get("B009")` => `B009 beauty HARD`, registry count `49`.

CODE-OK