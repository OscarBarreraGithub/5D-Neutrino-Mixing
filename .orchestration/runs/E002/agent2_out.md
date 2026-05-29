1. NIT: Correct-amplitude check OK. E002 evaluates the EDM magnitude `|d_mu^NP|`; complex chiral proxy uses `Im(C_chiral)` as the CP-odd part, not `Re` or full `|C|` (`quarkConstraints/edm.py:94`, `quarkConstraints/edm.py:183`, `flavor_catalog_constraints/primary/edm_neutrino/E002.py:248`).

2. NIT: QCD-running check is N/A for E002. This is a charged-lepton low-energy EDM proxy, not a Delta-F=2 `M_12` Wilson path; no `mu_had=2 GeV` QCD running is physically required or available here (`flavor_catalog_constraints/primary/edm_neutrino/E002.py:39`, `quarkConstraints/edm.py:3`).

3. NIT: Budget OK. Code uses the YAML canonical direct limit as a pure-NP HARD bound: `|d_mu| < 1.8e-19 e cm`; equivalent proxy-coefficient limit is `9.1219e-6 GeV^-1` using `hbar*c=1.973269804e-14 GeV cm` (`E002.py:49`, `E002.py:203`, `quarkConstraints/edm.py:40`).

4. NIT: Anchor numbers OK. YAML has `canonical_direct_limit` value `1.8e-19 e cm`, `95%`, table value `1.8` in `10^-19 e cm` (`flavor_catalog/processes/edm_neutrino/E002.yaml:100`); PDG snapshot matches `<1.8`, Bennett 2009 (`flavor_catalog/references/E002/pdg2026_muon_edm_datablock.txt:17`).

5. NIT: Severity/SM/units OK. `Severity.HARD`, SM set to `0.0`, evaluated prediction is `abs_edm_e_cm`, and ratio is `|d_mu|/1.8e-19` (`E002.py:151`, `E002.py:52`, `E002.py:244`). Units are consistently `GeV^-1 -> e cm`.

6. NIT: RS proxy caveat OK. The code does not claim a full RS one-loop prediction; it flags the explicit low-energy proxy as `NEEDS-HUMAN-PHYSICS` (`E002.py:16`, `E002.py:253`, `flavor_catalog_constraints/physics_adapters/edm.py:42`).

PHYSICS-OK