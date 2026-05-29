Plan implemented.

Physics: added shared `c -> u l+l-` short-distance C9/C10 machinery for `D0 -> l+l-`, with D0→μμ using the SD branching fraction against the C004 PDG upper limit. C004 YAML has LD SM context, not a nonzero SD Wilson anchor, so SM SD validates to `0.0`; YAML LD VMD context is `9.45e-13`.

Budget: HARD, `2.1e-9` from `C004.yaml` `canonical_current_limit`.

Built vs reused: built `quarkConstraints/rare_charm_dilepton.py`; wrapped via new adapter `flavor_catalog_constraints/physics_adapters/rare_charm_dilepton.py`. Constraint imports only the adapter.

NEEDS-HUMAN-PHYSICS: RS NP uses documented Z/KK-penguin C9/C10 proxy because full EW KK/Z/Z′, lepton, scalar/pseudoscalar matching is not on `ParameterPoint`.

Files changed:
`quarkConstraints/rare_charm_dilepton.py`
`flavor_catalog_constraints/physics_adapters/rare_charm_dilepton.py`
`flavor_catalog_constraints/primary/charm/C004.py`
`tests/constraints/primary/charm/test_C004.py`

Validation: safe proxy point ratio `6.46e-05`; excluded proxy point ratio `16.16`.

Tests: `python -m pytest tests/constraints/ -q` → `213 passed`.