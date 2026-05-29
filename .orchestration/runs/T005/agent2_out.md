1. NIT: ΔF=2 amplitude/Re-Im checks are N/A for T005; this is a two-body BR. The code correctly uses `|ζ_L|^2+|ζ_R|^2`, giving `Γ=(2/3) αs(mt) mt (|ζ_L|^2+|ζ_R|^2)` in `quarkConstraints/top_fcnc.py:238`-`245`.

2. NIT: The `mu_had=2 GeV *_with_running` requirement is N/A here. A real RS chromodipole result would run `C_uG` from `M_KK` to `m_t`, not to 2 GeV; the proxy explicitly says no QCD running and is flagged NEEDS-HUMAN-PHYSICS in `flavor_catalog_constraints/physics_adapters/top_fcnc.py:373`-`379`.

3. NIT: Budget is defensible from YAML: active limit is ATLAS `BR(t->cg)<3.7e-4` vs CMS `4.1e-4`, with SM `4.6e-12` negligible by `8.0e7`; see `T005.py:347` and `T005.yaml:140`-`180`.

4. NIT: Anchors match snapshots: ATLAS `4.7 pb`, `0.14 TeV^-2`, `3.7e-4`; CMS `1.8e-2 TeV^-1`, `4.1e-4`; SM `4.6e-12`. One older PDG snapshot value `2.0e-4` is intentionally not in YAML; see `pdg_2025_top_tug_tcg.txt:36`-`40` and `T005.yaml:262`-`265`.

PHYSICS-OK