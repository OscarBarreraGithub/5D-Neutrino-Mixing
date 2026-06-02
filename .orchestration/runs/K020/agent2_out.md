1. NIT: Correct-amplitude verdict PASS. K020 is a branching-fraction bound, not ΔF=2; the rate uses full complex magnitudes `|YV|^2`, `|YA|^2`, not Re/Im projections, which is correct for charged `K+ -> pi+ mu+ e-` BR. `flavor_catalog_constraints/physics_adapters/rare_kaon_lfv_semileptonic.py:596-598`.

2. NIT: QCD-running verdict N/A/PASS. No `*_with_running` ΔF=2 path is appropriate here; this is a color-singlet ΔS=1 semileptonic vector/axial proxy with LO QCD factor 1.0, so running effect is 0% in this restricted model. Code does not misuse the non-running ΔF=2 verdict path. `K020.py:445-451`.

3. NIT: Budget verdict PASS. Pure-NP SM≈0 bound uses `budget = 1.3e-11`, exactly the PDG/BNL 90% CL row for `BR(K+ -> pi+ mu+ e-)`; no central residual subtraction is needed. `K020.py:124-138`, `flavor_catalog/processes/secondary/kaon/K020.yaml:121-140`.

4. NIT: Anchor-number verdict PASS. YAML/snapshot numbers match: current `mu+ e-` limit `1.3e-11`, E865-only `2.1e-11`, opposite charge `6.6e-11`, all 90% CL. `K020.yaml:125,143,164`; snapshot `flavor_catalog/references/K020/pdg2025_kplus_lfv_semileptonic_api.txt:7-8`.

5. NIT: Severity/units verdict PASS. HARD is appropriate for a null LFV upper limit; branching fractions are dimensionless, masses are GeV, lifetime is converted with `hbar_gev_s`. `K020.py:362-365,483-491`; adapter `rare_kaon_lfv_semileptonic.py:118-123,378-383`.

PHYSICS-OK