1. NIT: ΔF=2 amplitude/QCD-running checks are N/A for EW002: `evaluate()` calls only CKM-sum arithmetic, no `M12`, no Wilsons, no `deltaf2` path, no GeV units (`flavor_catalog_constraints/primary/top_higgs_ew/EW002.py:225`; `physics_adapters/ckm_unitarity.py:87`).

2. NIT: Anchor arithmetic is correct: PDG sum `0.9983` vs SM target `1.0` gives `Delta_CKM=-0.0017`; budget `0.0007` gives pull `2.4286` (`EW002.py:227`; `EW002.yaml:131`).

3. NIT: Budget is defensible and uncertainty-aware: uses quoted first-row unitarity uncertainty `0.0007`, with PDG components `0.0006` and `0.0004`; not a bare residual (`EW002.py:181`; `EW002.yaml:135`).

4. NIT: Loaded inputs match YAML/snapshots: `Vud=0.97367(32)`, `Vus=0.22431(85)`, first-row sum `0.9983(6)(4)`; no standalone `Vub` YAML block, and code correctly relies on the PDG sum including tiny `|Vub|^2` (`EW002.py:232`; `EW002.yaml:86`).

5. NIT: Severity/NP handling is appropriate: `SOFT`, `np_shift_delta_ckm=0.0`, and `NEEDS-HUMAN-PHYSICS` is explicit rather than inventing an RS charged-current shift (`EW002.py:218`; `EW002.py:279`).

PHYSICS-OK