1. NIT: ΔF=2 amplitude check is N/A for E008: no `M12`, Δm, or CP-mixing observable is computed, so no wrong Re/Im component is used; verdict path records `predicted=None` only. `flavor_catalog_constraints/primary/edm_neutrino/E008.py:494`
2. NIT: QCD-running check is N/A by hard-stub design: E008 does not call a non-running ΔF=2 path, and does not compute cEDM Wilsons to run; future live cEDM matching would need explicit RG/mixing. `E008.py:497`, `physics_adapters/quark_cedm.py:93`
3. NIT: Budget is defensible only as non-vetoing reference: Hg gives `7.4e-30 / 7.0e-3 = 1.057e-27 cm`, rounded to `1.1e-27 cm`; code keeps `Severity.INFO`, `passes=True`, `predicted=None`. `E008.py:503`, `E008.py:511`
4. NIT: Anchor spot-checks match YAML/audit: neutron `1.8e-26 e cm -> 1.6e-26 cm`, Hg `7.4e-30 e cm -> 1.1e-27 cm`. `flavor_catalog/processes/edm_neutrino/E008.yaml:87`, `E008.yaml:145`
5. NIT: Units are consistent for this EDM stub: EDM anchors are `e cm`, qCEDM bounds are `cm`; GeV M12 budgets are not applicable. `E008.py:233`, `E008.py:352`
6. NIT: Wording caveat: “strongest” compares different qCEDM linear combinations, so it is only the smallest catalogued scalar reference, not a model-independent dominance statement; non-vetoing diagnostics name the chosen observable. `E008.py:171`

PHYSICS-OK