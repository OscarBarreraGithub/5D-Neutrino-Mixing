1. No BLOCKER / SHOULD-FIX / NIT findings for L006.

Isolation: L006 adds only new files `flavor_catalog_constraints/primary/charged_lepton/L006.py`, `flavor_catalog_constraints/physics_adapters/muonium_conversion.py`, and `tests/constraints/primary/charged_lepton/test_L006.py`; no pre-existing `deltaf2.py`, `base.py`, or `point_builder.py` diff. No sibling-constraint or direct `quarkConstraints` imports.

Independent numeric check: YAML `P_limit=8.3e-11`, `G_limit/GF=0.003`. Safe `g=1e-4`: predicted `9.222222222222222e-14`, ratio `0.0011111111111111111`, pass. Excluded `g=1e-2`: predicted `9.222222222222222e-10`, ratio `11.11111111111111`, fail. Complex check `1.5e-4+2.0e-4j`: predicted `5.763888888888889e-13`, ratio `0.006944444444444444`, matches constraint.

Anchor probes: missing anchor raises `AnchorError`; mismatched MACS/PDG probability raises `AnchorError` with `8.3e-11 vs 7.7e-11`.

Focused tests: `10 passed`. Full command `python -m pytest tests/constraints/ -q`: `959 passed in 19.53s`. Registry smoke: `L006 charged_lepton PRIMARY HARD P(M -> Mbar)`, 94 registered constraints.

CODE-OK