1. BLOCKER `flavor_catalog_constraints/physics_adapters/lfv_three_body_tau.py:141`: matrix inputs without explicit flavors use the shared core’s `mu->e` matrix default while diagnostics say `tau,mu`. Evidence: used overlap `(0.09+0.01j)` instead of expected `(0.003+0.004j)`, predicted `5.5422131956696402e-09` vs tau-mu expected `1.6896991450212318e-11` (328x). Fix by normalizing contact inputs to a pinned `LFVThreeBodyContactProxyInput` or injecting `initial_flavor="tau", final_flavor="mu"` before calling the core; add a matrix pinning test.

2. SHOULD-FIX `flavor_catalog_constraints/physics_adapters/lfv_three_body_tau.py:149`: dipole-only explicit BR still requires `m_kk_gev` via the zero-contact placeholder. Evidence: `{"dipole_parent_branching_fraction": 1e-8}` returns unevaluated `KeyError`. Fix by allowing pure dipole evaluation without MKK, using a harmless zero-contact placeholder/reference scale.

Cross-check: scalar proxy independent core recomputation matched exactly: predicted `4.9924402074646848e-11`, expected `4.9924402074646848e-11`, diff `0`, ratio `0.0026276001091919392`.

Safe/excluded: safe `box_ll=1e-5` predicted `2.0000000000000003e-10`, ratio `0.010526315789473686`, pass `True`; excluded `box_ll=5e-4` predicted `4.9999999999999998e-07`, ratio `26.315789473684209`, pass `False`.

Contract/anchors/determinism: numeric fields are real floats, complex amplitudes stay in diagnostics, missing input is non-crashing, missing/mismatched anchors raise `AnchorError`, repeated evaluation is equal and did not mutate the point.

Isolation: L009 adds only `L009.py`, `lfv_three_body_tau.py`, and `test_L009.py`; tracked shared `quarkConstraints/lfv_three_body.py` and L002 adapter are unchanged; constraint imports physics only through the tau adapter.

Registry smoke: `SMOKE L009 charged_lepton HARD BR(tau -> 3mu)`.

Pytest: `568 passed in 14.12s` for `python -m pytest tests/constraints/ -q`.

CODE-NEEDS-FIXES