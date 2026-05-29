Plan: read scaffold/K001/K004 plus E004/C003 stubs, add E006 stub, add adapter/test, run full constraints.

Physics: `|d_Hg|` is recorded from `E006.yaml`; no Schiff-moment, atomic-structure, or RS CP-odd matching calculation is faked. Severity is `INFO`, non-vetoing.

Reused vs built: reused scaffold `load_anchor`/registry contracts; built new small adapter `physics_adapters/mercury_edm.py`.

Budget/source: `|d_Hg| < 7.4e-30 e cm` at 95% CL from Graner et al. 2016, loaded from `canonical_direct_limit`.

Validation: YAML central value `2.20e-30 e cm`, total uncertainty `3.122963e-30 e cm`, central/limit ratio `0.297297297297`.

NEEDS-HUMAN-PHYSICS: flagged both `Schiff moment + atomic structure` and `RS CP-odd matching`.

Files changed: `flavor_catalog_constraints/primary/edm_neutrino/E006.py`, `flavor_catalog_constraints/physics_adapters/mercury_edm.py`, `tests/constraints/primary/edm_neutrino/test_E006.py`.

Tests: `8 passed` for E006; `778 passed` for `python -m pytest tests/constraints/ -q`.