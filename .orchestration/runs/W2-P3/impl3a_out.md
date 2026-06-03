Plan executed: schema keys, 3a builder, quark Z matrices, contacts/Wilsons, gate tests; no constraint rewiring.

New API:
`flavor_catalog_constraints.point_builder.build_from_rs_ew_inputs`
New modules:
`flavor_catalog_constraints/rs_ew_builder.py`, `quarkConstraints/rs_ew_couplings.py`, `quarkConstraints/rs_semileptonic_wilsons.py`.

Asserted numbers:
`z_delta_g_R^b(3 TeV) = -1.6604227573826689e-3`; `6 TeV = -4.151056893456672e-4`; scaling ratio `0.25`.
SM-limit max checks: `z=0`, contacts `0 GeV^-2`, Wilsons `0`.
Sample manual-recompute Wilsons matched builder at `rel=1e-12`: `b->s mumu C9=0.3046272275+0.0058031605j`, `C10=-4.0552080339-0.0772518704j`; `s->d ee C9=-0.7330216826-0.2065809673j`; `c->u mumu C10=-20.4407745025-31.3859401666j`.

Tests:
`python -m pytest tests/test_rs_ew_phase3a.py -q` -> `7 passed`.
`python -m pytest tests/ -q` -> `1652 passed, 1 skipped`.

P3A-DONE.