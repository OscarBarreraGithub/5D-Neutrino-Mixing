1. SHOULD-FIX `flavor_catalog_constraints/primary/beauty/B003.py:360`/`:387`: top-level `sm_prediction` and `experimental` are half-M12 values, not `Delta m_s`; use `sm_delta_m_gev` / `experimental_delta_m_gev`, keep M12 in diagnostics.
2. SHOULD-FIX `flavor_catalog_constraints/primary/beauty/B003.py:315`: anchor probe with `units: GeV` + wrong `observable` still loads and is converted as ps^-1; validate `units == "ps^-1"` and `observable == "Delta m_s"` before budget construction.
3. NIT `flavor_catalog_constraints/primary/beauty/B003.py:206`: missing auxiliary anchor raises loudly, but message says `pdg_or_equivalent`; pass block context or raise a local auxiliary-specific `AnchorError`.

Isolation: no other constraint imports in B003; `git diff quarkConstraints/` empty; B_s adapter additions are additive and existing function bodies/signatures are unchanged.
Contract/determinism: numeric fields are real floats; complex couplings/Wilsons stay in diagnostics; missing couplings returns non-crashing pass; same input gives equal result and no coupling mutation.
Anchor probes: missing pdg key, missing aux key, bad exp numeric, bad aux numeric all raise `AnchorError`.

Cross-check safe: predicted/direct `4.76038552508622049e-14`, budget `2.63516764862967608e-13`, ratio `0.180648298697871745`, pass `True`.
Cross-check excluded: predicted/direct `4.76038552508622097e-12`, budget `2.63516764862967608e-13`, ratio `18.0648298697871752`, pass `False`.
Direct core legacy budget was `5.84399999999999974e-12`; B003 correctly uses tighter catalog budget for ratio/pass.
Tests: `python -m pytest tests/constraints/ -q` -> `51 passed in 2.42s`; B003 alone -> `8 passed in 2.24s`.
Registry smoke: `fcc.get("B003")` imports; `all_constraints()` count `5`; `import_failures {}`.
CODE-NEEDS-FIXES