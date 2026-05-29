Plan completed: inspected K001/K004, T010/T012, then added an LFV zpole extension, T015 constraint, and tests.

Physics: `BR(Z -> e mu)` is pure NP with SM set to `0.0`; rate uses zpole effective couplings with charge-summed factor `2*(|dL|^2+|dR|^2)`.  
Source/budget: CMS 2025 limit from `T015.yaml`, `BR < 1.9e-7` at 95% CL.  
Validation: independent SM total Z-width weight check gives `3.649563333333334`; zero LFV coupling gives `BR=0.0`.  
Gap: flagged `NEEDS-HUMAN-PHYSICS`; off-diagonal lepton Z coupling is proxied from supplied lepton overlap spurions because it is not on `ParameterPoint`.

Reused: existing zpole coupling/width conventions. Built: append-only `quarkConstraints/zpole_lfv.py` plus adapter `flavor_catalog_constraints/physics_adapters/zpole_lfv.py`.  
Files added: `T015.py`, `test_T015.py`, `zpole_lfv.py` core + adapter. No T010/T012 diagonal functions modified.  
Tests: `python -m pytest tests/constraints/ -q` -> `319 passed`.