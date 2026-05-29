Promoted T011 from the CLOSED-AS-MERGED/T010 stub to standalone YAML anchors.
Anchors: `A_FB^{0,b}=0.0992 +/- 0.0016`, `A_b=0.923 +/- 0.020` from PDG 2025 plus LEP/SLC final-combination provenance.
Added T011-local reference snapshots and source manifest under `flavor_catalog/references/T011/`.
Runtime now requires T011 `pdg_or_equivalent` anchors via `load_anchor`; legacy `canonical_home`/`merged_into` T010 fallback raises `AnchorError`.
Physics still reuses `quarkConstraints/zpole.py`: `A_b=(|g_L|^2-|g_R|^2)/(|g_L|^2+|g_R|^2)`, `A_FB^{0,b}=3/4 A_e A_b`.
SM validation: `A_b=0.935535274431`, `A_FB^{0,b}=0.103278860422`.
RS Zbb-shift proxy remains `NEEDS-HUMAN-PHYSICS`.
Tests: `python -m pytest tests/constraints/primary/top_higgs_ew/test_T011.py -q` -> 12 passed.
Tests: `python -m pytest tests/constraints/ -q` -> 863 passed.
