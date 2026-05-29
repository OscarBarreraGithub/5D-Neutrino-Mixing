1. NIT: Check 1 PASS/N-A: E006 computes no ΔF=2 `M12` amplitude and no RS EDM amplitude; it records scalar `|d_Hg|` via `abs(measured)/limit`, with `predicted=None`. Correct for an EDM stub. `flavor_catalog_constraints/physics_adapters/mercury_edm.py:100`, `flavor_catalog_constraints/primary/edm_neutrino/E006.py:456`.

2. NIT: Check 2 PASS/N-A: no QCD-running verdict is expected because no CP-odd Wilson matching is implemented; the code avoids Delta-F=2 non-running paths and flags RS CP-odd matching absent. Running effect is not quantifiable in this stub. `flavor_catalog_constraints/primary/edm_neutrino/E006.py:448`, `flavor_catalog_constraints/physics_adapters/mercury_edm.py:45`.

3. NIT: Check 3 PASS: budget is the direct Graner 95% CL atom-EDM limit `7.4e-30 e cm`, not a central residual; central `2.20e-30`, total uncertainty `3.123e-30`, ratio `0.297` is bookkeeping only. `flavor_catalog/processes/edm_neutrino/E006.yaml:103`, `flavor_catalog_constraints/primary/edm_neutrino/E006.py:450`.

4. NIT: Check 4 PASS: anchor numbers match local snapshots: Graner `2.20 +/- 2.75_stat +/- 1.48_syst` and `|d_Hg| < 7.4e-30 e cm`; PDG Hg-theory rows `<0.16`, `<0.22` in `10^-25 e cm`; Sahoo context values match. `flavor_catalog/references/E006/graner2016_arxiv1601_04339.txt:20`.

5. NIT: Check 5 PASS: severity/units/diagnostics are appropriate: INFO, non-vetoing, EDM units `e cm`, dual `NEEDS-HUMAN-PHYSICS` for Schiff+atomic and RS CP-odd matching, no fake Schiff/atomic calculation. `flavor_catalog_constraints/primary/edm_neutrino/E006.py:440`, `flavor_catalog_constraints/primary/edm_neutrino/E006.py:469`.

PHYSICS-OK