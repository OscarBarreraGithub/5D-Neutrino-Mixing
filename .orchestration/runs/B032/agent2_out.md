1. NIT: No issue on amplitude-part usage. B032 is not a DeltaF=2 `M12` observable, so neither `|M12^NP|` nor `Im M12^NP` should appear; code returns `predicted=None`, `sm_prediction=None` at `flavor_catalog_constraints/primary/beauty/B032.py:589`.

2. NIT: QCD-running check is N/A and correctly not forced. This is DeltaB=1 charmless nonleptonic, not the `deltaf2.py` path; adapter explicitly avoids QCDF/SCET, penguin matching, RG, and matrix elements at `flavor_catalog_constraints/physics_adapters/charmless_b_pik.py:72`.

3. NIT: Budget is defensible for a non-vetoing stub: `Delta A_CP = 2.7% - (-8.31%) = 11.01%`, `sigma = sqrt(1.2^2+0.31^2)=1.239%`; code uses full `|Delta A_CP|=0.1101` as INFO-only room at `B032.py:511` and `B032.py:577`.

4. NIT: Anchors match YAML/snapshots: BRs `2.392, 1.322, 2.007, 1.012 x10^-5`; ACPs `-2.67, +2.7, -8.31, -1%`; PDG `C=0.00+-0.08`, `S=0.64+-0.13` in `flavor_catalog/processes/beauty/B032.yaml:107`.

5. NIT: Severity/units/diagnostics are appropriate: `Severity.INFO` at `B032.py:568`; percent ACPs converted to dimensionless at `B032.py:377`; dual `NEEDS-HUMAN-PHYSICS` flags at `B032.py:613`.

PHYSICS-OK