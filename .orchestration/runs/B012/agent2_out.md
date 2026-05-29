1. NIT: Rate amplitude is correct for this BR observable: unpolarized exclusive rate uses `|C7|^2 + |C7p|^2`, not Re/Im pieces; implemented via `exclusive_btokstargamma_from_couplings` at `flavor_catalog_constraints/primary/beauty/B012.py:221` and core dipole power at `quarkConstraints/bsgamma.py:301`. No deltaF2 amplitude issue applies.

2. NIT: QCD running is actually used for coupling-based verdicts: `compute_bsgamma_wilsons` runs `C7/C8` to `mu_b=4.8 GeV` at `quarkConstraints/bsgamma.py:346`; sample 3 TeV coefficients give `u77=0.512`, `u78=0.119`, `u88=0.557`, net `C7=C8` factor `0.631`. Correct for b -> s gamma; the prompt’s `2 GeV` deltaF2 check is not applicable.

3. NIT: Budget is defensible for the requested proxy: `BR_norm = 4.163e-5`, `sigma = 9.2e-7` from `B012.yaml:101`; code uses `|BR_total - BR_norm| <= sigma` at `B012.py:166`, giving a 2.21% HARD band `[4.071, 4.255]e-5`. This validates the requested SM-scale `~4.3e-5` within about 3.2%.

4. NIT: Anchor loading matches snapshots: YAML `4.163e-5 +/- 9.2e-7` and charged cross-check `3.95e-5 +/- 1.0e-6` at `B012.yaml:101` and `B012.yaml:112` match HFLAV snapshot lines `27` and `15`; audit also verifies them at `flavor_catalog/audits/factcheck_beauty.md:96`.

5. NIT: Diagnostics inherit inclusive metadata `photon_energy_cut_gev=1.6` from `quarkConstraints/bsgamma.py:409`; exclusive `B -> K* gamma` has no inclusive photon-energy cut. The exclusive diagnostics at `physics_adapters/bsgamma.py:132` correct the formula, so verdict is unaffected, but that field is potentially misleading.

PHYSICS-OK