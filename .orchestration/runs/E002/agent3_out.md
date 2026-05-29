1. NIT: No E002 code/test findings requiring changes.

Evidence: E002 imports physics only through `flavor_catalog_constraints.physics_adapters.edm`; relevant tracked EDM adapter/core/E001 files are clean, with only `E002.py` and `test_E002.py` added.

Numerical cross-check against `quarkConstraints.edm`:
`c=1.0e-6`: predicted/core `1.973269804000e-20`, ratio `0.109626100222`, pass.
`c=2.0e-6`: predicted/core `3.946539608000e-20`, ratio `0.219252200444`, pass.
`c=2.0e-5`: predicted/core `3.946539608000e-19`, ratio `2.19252200444`, fail.

Anchor probes: missing anchor raises `AnchorError`; noncanonical `pdg_combined_measurement` raises `AnchorError` requiring `canonical_direct_limit`.

Registry import smoke: OK, `E002 edm_neutrino HARD |d_mu|`.

Pytest: `tests/constraints/primary/edm_neutrino/test_E002.py -q`: 12 passed. Full `python -m pytest tests/constraints/ -q`: 568 passed.

CODE-OK