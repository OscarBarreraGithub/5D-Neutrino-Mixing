1. NIT: Correct-amplitude verdict OK/N.A. C003 forms no `M12`; a live result would need CP-odd imaginary interference in `Delta C=1` decay amplitudes, not `|M12|` or `Im M12` from mixing. See `flavor_catalog_constraints/primary/charm/C003.py:418` and `physics_adapters/charm_direct_cp.py:71`.

2. NIT: QCD-running verdict OK/N.A. No non-running `deltaf2` path feeds the verdict; the stub does no Wilson calculation. A production implementation would require `Delta C=1` penguin/dipole RG, not the neutral-meson `*_with_running` evaluators. See `C003.py:410` and `charm_direct_cp.py:77`.

3. NIT: Budget verdict OK for INFO bookkeeping. Budget is full observed `|Delta A_CP| = 1.54e-3`, with `sigma = 2.9e-4` recorded, ratio `1.0`, `passes=True` advisory only. This is not an SM-subtracted exclusion; it would be too tight only if promoted to SOFT/HARD. See `C003.py:167`, `C003.py:421`, `C003.py:461`.

4. NIT: Anchor numbers match. LHCb `(-15.4 +/- 2.9)e-4 = (-0.154 +/- 0.029)%` matches YAML and snapshot; HFLAV `Delta a_CP^dir=(-0.159 +/- 0.029)%`, `a_CP^ind=(-0.010 +/- 0.012)%`, no-CPV `5.3 sigma` also match. See `flavor_catalog/processes/charm/C003.yaml:90`, `C003.yaml:126`, `references/C003/lhcb2019_arxiv1903_08726.txt:25`.

5. NIT: Severity/diagnostics OK. `Severity.INFO`, `non_vetoing=True`, `predicted=None`, `sm_prediction=None`, both SM and RS `NEEDS-HUMAN-PHYSICS` flags present; units are dimensionless asymmetry, so no GeV issue. See `C003.py:403`, `C003.py:422`, `C003.py:427`, `charm_direct_cp.py:32`.

PHYSICS-OK