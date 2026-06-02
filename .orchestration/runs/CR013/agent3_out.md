1. NIT [kk_graviton_resonance.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/physics_adapters/kk_graviton_resonance.py:4): module docstring still says it adds only the CR007-specific mass proxy; fix by mentioning the appended CR013 diphoton helper too.

Isolation OK: CR013 imports no other constraint, reaches physics via adapter; AST check on pre-existing `kk_graviton_resonance.py` defs showed `changed=[]`, `missing=[]`, added only `kk_graviton_diphoton_prediction_from_lambda_ir_gev`.

Anchor/contract OK: active YAML row is `PDG2025:CR013:CMS2024_RSG_diphoton_kMPl_0p1`, value `4.8 TeV`, `ktilde=0.1`; missing and mismatched-active probes both raised `AnchorError`. Missing-input path returns non-crashing `passes=True`, `predicted=None`, `ratio=None`, real-float `experimental/budget`.

Numerical cross-check: for `kk_ew_mass_gev=3100`, direct core gives `m_G=4.85088043162 TeV`, ratio `0.989511093431`; CR013 returns the same and passes. Boundary probe: safe `m_G=6.71585298510 TeV`, ratio `0.714726783127`, pass; excluded `m_G=4.41682940298 TeV`, ratio `1.08675241040`, fail.

Tests: `python -m pytest tests/constraints/primary/collider_rs/test_CR013.py -q` -> `8 passed`; `python -m pytest tests/constraints/ -q` -> `1033 passed in 34.08s`. Registry smoke: `CR013 collider_rs PRIMARY HARD`, registered count `101`.

CODE-OK