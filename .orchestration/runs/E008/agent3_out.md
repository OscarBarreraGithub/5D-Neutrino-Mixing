1. SHOULD-FIX `flavor_catalog_constraints/primary/edm_neutrino/E008.py:44`: E008 bypasses `load_anchor` and composes `load_pdg_block`/`find_block`/`build_anchor` directly at `:249`/`:270`; fix by routing the selected qCEDM anchors through the scaffold `load_anchor` path, or an explicit dotted-anchor wrapper, and add a mismatched-anchor probe like E004.
2. NIT `tests/constraints/primary/edm_neutrino/test_E008.py:63`: missing-anchor loud fail is probed, but mismatched selected-anchor behavior is not; fix by monkeypatching the anchor loader to return the wrong block and asserting `AnchorError`.

Cross-check numbers: `result.experimental=budget=1.1e-27`; neutron recompute `1.8e-26/1.1=1.6363636363636362e-26` vs YAML `1.6e-26`; Hg recompute `7.4e-30/7.0e-3=1.0571428571428571e-27` vs YAML `1.1e-27`; safe adapter point `5.5e-28/1.1e-27=0.5 pass`; excluded `2.2e-27/1.1e-27=2.0 fail`.

Pytest: `python -m pytest tests/constraints/ -q` -> `810 passed in 17.42s`. Registry smoke: `fcc.get("E008")` OK (`edm_neutrino`, `INFO`).  
CODE-NEEDS-FIXES