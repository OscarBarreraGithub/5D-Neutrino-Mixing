1. BLOCKER — Chirality/C7′ assignment is reversed for the stated vector-overlap proxy: with `O7=(sbar_L sigma b_R)` and `O7p=L<->R`, left-handed `b->s` overlap should feed `C7`, right-handed overlap should feed `C7p`; code maps right→`C7` and left→`C7p`, changing SM interference. `quarkConstraints/bsgamma.py:24`, `quarkConstraints/bsgamma.py:41`, `quarkConstraints/bsgamma.py:203`.

2. BLOCKER — No `b->s gamma` RG evolution is used in the HARD verdict: `C7_NP` built at `matching_scale=M_KK` is inserted directly into a `C7_SM(mu_b)=-0.304` rate. Correct physics must run `C7,C8` from `M_KK` to `mu_b`, or explicitly define the proxy as low-scale `C7_eff`. LL 3 TeV→4.8 GeV gives `C7(mu_b)≈0.50 C7(MKK)+0.121 C8(MKK)`. `quarkConstraints/bsgamma.py:13`, `quarkConstraints/bsgamma.py:48`, `quarkConstraints/bsgamma.py:203`, `flavor_catalog_constraints/primary/beauty/B011.py:240`.

3. NIT — Mandatory ΔF=2 `M12` real/imag check is N/A for B011; for the BR observable, the leading dipole rate uses the right magnitude structure `|C7+C7_NP|^2+|C7p+C7p_NP|^2`. `quarkConstraints/bsgamma.py:162`, `flavor_catalog_constraints/primary/beauty/B011.py:8`.

4. NIT — Budget is uncertainty-aware and matches YAML/snapshots: exp `(3.49±0.19)e-4`, SM `(3.40±0.17)e-4`, central gap `0.09e-4`, combined sigma `0.255e-4`, HARD budget `0.345e-4`. `flavor_catalog/processes/beauty/B011.yaml:75`, `flavor_catalog_constraints/primary/beauty/B011.py:128`.

5. NIT — Anchors/units/severity are otherwise consistent: SM limit is `3.40e-4` at `E_gamma>1.6 GeV`, validation SM is `3.36±0.23e-4`, units are branching fraction, HARD is defensible, and `NEEDS-HUMAN-PHYSICS` is surfaced. `flavor_catalog/processes/beauty/B011.yaml:91`, `flavor_catalog_constraints/primary/beauty/B011.py:317`.

PHYSICS-NEEDS-FIXES