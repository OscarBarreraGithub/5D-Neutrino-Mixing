1. BLOCKER: none found.
2. SHOULD-FIX: none found.
3. NIT: none found.

Evidence: CR003 imports no other constraint; physics access is via `flavor_catalog_constraints.physics_adapters.collider_resonance`; adapter AST check showed `changed_preexisting=[]`, `removed_preexisting=[]`, added only `resolve_charge_two_thirds_vlq_mkk_gev` and `vlq_t_pair_prediction_from_m_kk_gev`.

Numerical cross-check against `quarkConstraints.collider_resonance`: `M_KK=2000 GeV -> predicted=2.0 TeV, limit=1.48 TeV, ratio=0.74, pass=True`; `1800 GeV -> ratio=0.8222222222222222, pass=True`; `1200 GeV -> ratio=1.2333333333333334, pass=False`.

Contract/anchor probes: missing couplings returned non-crashing `passes=True`, `predicted=None`, `ratio=None`; mismatched anchor raised `AnchorError`; deterministic probe returned identical results and no point mutation.

Pytest: `python -m pytest tests/constraints/ -q` -> `626 passed in 14.44s`. Registry smoke: `CR003 collider_rs HARD`, registered `True`.

CODE-OK