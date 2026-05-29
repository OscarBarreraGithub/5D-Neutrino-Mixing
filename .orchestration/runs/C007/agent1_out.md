Plan executed: reused C004 `rare_charm_dilepton` Wilson proxy, added only D→πμμ short-distance rate, then C007 + tests.

Physics: vector/axial `C9+C9'`, `C10+C10'` D→π form-factor proxy; LD resonances documented, not modeled. Sources: `C007.yaml`, [LHCb 2021](https://arxiv.org/abs/2011.00217), [Sahoo et al. EPJC 2017](https://link.springer.com/article/10.1140/epjc/s10052-017-4888-4).

Budget: HARD `6.7e-8` from `flavor_catalog/processes/charm/C007.yaml` `canonical_current_limit`.

Validation: YAML SM SD scale `1.0e-12`; zero-Wilson smooth formula `0.0`; pinned NP proxy check `9.132644202844895e-12`.

NEEDS-HUMAN-PHYSICS: full RS EW/lepton/scalar/tensor/resonance/window recast still missing; diagnostics flag it.

Files changed: `quarkConstraints/rare_charm_semileptonic.py`, `physics_adapters/rare_charm_semileptonic.py`, `primary/charm/C007.py`, `tests/constraints/primary/charm/test_C007.py`.

Tests: `python -m pytest tests/constraints/primary/charm/test_C007.py -q` → 11 passed; `python -m pytest tests/constraints/ -q` → 341 passed.