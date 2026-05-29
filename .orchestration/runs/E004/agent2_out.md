1. NIT: Amplitude check N/A for E004: this is not a ΔF=2 `M12` observable, and the code correctly computes no RS EDM amplitude (`predicted=None`, `sm_prediction=None`); it only records `|d_n|` bookkeeping. `flavor_catalog_constraints/primary/edm_neutrino/E004.py:398` numbers: limit `1.8e-26 e cm`, ratio `0.0`.

2. NIT: QCD-running check N/A for this honest EDM stub: no Wilsons are matched or evolved, and it does not use the non-running DeltaF=2 path. Correct physics requires CP-odd quark EDM/qCEDM/Weinberg matching plus EDM RG before any running verdict. `E004.py:51`, `E004.py:391`, `physics_adapters/neutron_edm.py:85`.

3. NIT: Budget is defensible for a non-vetoing stub: it uses the experimental upper limit itself, not a fake SM residual or hadronic translation. `E004.py:404`, `E004.py:406`; YAML `flavor_catalog/processes/edm_neutrino/E004.yaml:84` value `1.8e-26 e cm`.

4. NIT: Anchor numbers match snapshots: PDG table `<0.18` in `10^-25 e cm` converts to `1.8e-26 e cm`; Abel gives `d_n=(0.0 +/- 1.1_stat +/- 0.2_sys)e-26 e cm`. YAML `E004.yaml:78`; snapshots `pdg2026_neutron_edm_datablock.txt:18`, `abel2020_arxiv2001_11966.txt:23`.

5. NIT: Severity/units/diagnostics are physics-appropriate: `Severity.INFO`, e cm throughout, `passes=True` is non-vetoing, and both missing ingredients are flagged: hadronic matrix elements and RS CP-odd quark-dipole matching. `E004.py:382`, `E004.py:414`, `E004.py:424`, `E004.py:427`.

PHYSICS-OK