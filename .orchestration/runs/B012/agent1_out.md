Plan: read B011/K001/K004/scaffold, append exclusive wrapper, add B012 constraint, add tests, run suite.

Physics: `BR(B0 -> K*0 gamma) = BR_norm * |C7_eff|^2` ratio using shared `C7/C7p` power; `BR_norm` comes from `B012.yaml` HFLAV neutral anchor.

Reused vs built: reused `quarkConstraints/bsgamma.py` C7/C8 matching + LL running; appended only exclusive B→K*γ adapter wrappers. B011 inclusive functions unchanged.

Budget: HARD `9.2e-7`, from `B012.yaml:branching_fraction_b0_to_kstar0_gamma`.

SM validation: zero-NP prediction is `4.163e-05`, within the requested `~4.3e-05` scale.

NEEDS-HUMAN-PHYSICS: flagged in diagnostics for the reused RS dipole proxy; full exclusive RS matching/form-factor/spectator inputs are not on `ParameterPoint`.

Files: `flavor_catalog_constraints/physics_adapters/bsgamma.py`, `flavor_catalog_constraints/primary/beauty/B012.py`, `tests/constraints/primary/beauty/test_B012.py`.

Tests: `python -m pytest tests/constraints/ -q` → `193 passed`. Note: git status also shows unrelated B006/B023/rare_b_nunu changes; I left them alone.