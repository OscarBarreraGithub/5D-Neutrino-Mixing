1. BLOCKER: none.
2. SHOULD-FIX: none.
3. NIT: none.

Isolation: EW003 imports only `physics_adapters.semileptonic_ckm`, no other constraint; `deltaf2.py` AST check showed `modified_existing=[]`, `removed_existing=[]`, `added=['bs_mixing_m12_np_from_wilsons_with_running']`.

Contract/probes: numeric fields are real finite floats; severity is `SOFT`; empty/missing-couplings point does not crash; `evaluate()` is deterministic and does not mutate the point.

Anchor probes: 12 `load_anchor` calls to YAML list entries; missing observable raises `AnchorError`; mismatched units raises `AnchorError`.

Independent YAML recompute:
PDG `|V_cb|`: diff `2.4`, combined sigma `0.781024967591`, pull `3.07288511838951`, ratio `1.02429503946317`, fails.
PDG `|V_ub|`: diff `0.43`, combined sigma `0.301330383466`, pull `1.42700511993994`, ratio `0.475668373313313`, passes.
EW003 result: predicted `3.07288511838951`, budget `3.0`, ratio `1.02429503946317`, passes `False`.

Registry smoke: `registry_size=13`, `import_failures=0`, `EW003_present=True`.

`python -m pytest tests/constraints/ -q`: `125 passed in 3.08s`.

CODE-OK