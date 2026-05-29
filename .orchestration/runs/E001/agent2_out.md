1. NIT: Correct observable projection OK: E001 is not an \(M_{12}\)/DeltaF2 amplitude; direct proxy is \(c_{\rm CPodd}=d_e/e\), complex chiral proxy uses `Im(C)`, and result compares `abs(c*hbarc)`. `quarkConstraints/edm.py:94`, `quarkConstraints/edm.py:111`, `quarkConstraints/edm.py:141`, `flavor_catalog_constraints/primary/edm_neutrino/E001.py:239`.

2. NIT: Running verdict N/A/OK: no QCD `*_with_running` DeltaF2 path is relevant for a low-energy electron EDM proxy; implemented running effect is factor 1 / 0%. Missing RS loop matching plus possible SMEFT/LEFT running is explicitly flagged `NEEDS-HUMAN-PHYSICS`. `quarkConstraints/edm.py:12`, `flavor_catalog_constraints/physics_adapters/edm.py:42`, `E001.py:224`.

3. NIT: Budget/severity OK: HARD pure-NP bound uses \(|d_e^{NP}| <= 4.1e-30 e cm\), not a central residual; SM is set to 0 at catalog precision. Equivalent proxy limit is \(4.1e-30 / 1.973269804e-14 = 2.078e-16 GeV^-1\). `E001.yaml:119`, `E001.py:145`, `E001.py:197`.

4. NIT: Anchors OK: YAML/snapshots match PDG/Roussy canonical limit `<0.041 x 10^-28 e cm = 4.1e-30 e cm` at 90% CL; ACME/Andreev 2018 `1.1e-29 e cm` is correctly superseded, not used as the budget. `E001.yaml:125`, `E001.yaml:131`, `flavor_catalog/references/E001/pdg2026_electron_edm_datablock.txt:18`.

5. NIT: Units/docs/diagnostics OK: conversion `GeV^-1 -> cm` via `hbar*c` is dimensionally consistent, notes state pure proxy and non-vetoing unevaluated semantics, and diagnostics clearly flag missing RS CP-violating dipole machinery. `quarkConstraints/edm.py:40`, `E001.py:153`, `E001.py:248`.

PHYSICS-OK