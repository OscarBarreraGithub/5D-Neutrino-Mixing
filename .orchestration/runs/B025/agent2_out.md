1. NIT: Correct-amplitude check is N/A to Δm/CP; B025 uses the declared LFU amplitude proxy `R_D = R_D^SM |1 + C_tau|^2`, with no Re/Im mixup. `quarkConstraints/semileptonic_lfu.py:191`, `quarkConstraints/semileptonic_lfu.py:216`, `flavor_catalog_constraints/primary/beauty/B025.py:373`.

2. NIT: QCD-running check is N/A; this is not a ΔF=2 `M_12` constraint and no 2 GeV hadronic Wilson evolution is expected. Verdict uses the semileptonic LFU adapter/core, not a non-running `deltaf2` path; running effect is 0 by construction of this proxy. `flavor_catalog_constraints/primary/beauty/B025.py:373`, `quarkConstraints/semileptonic_lfu.py:253`.

3. NIT: Budget is defensible and uncertainty-aware: HFLAV `R_D=0.358±0.024`, SM `0.296±0.004` gives residual `0.062`, combined σ `0.024331`, HARD budget `0.086331`; SM pull ratio is `0.718`, so this is not a bare central residual. `flavor_catalog_constraints/primary/beauty/B025.py:244`, `flavor_catalog_constraints/primary/beauty/B025.py:298`.

4. NIT: Anchor numbers match the YAML/snapshots: average `R(D)=0.358±0.024`, `rho=-0.374`, `chi2/dof=16.683/14`, SM arithmetic row `0.296±0.004`; units are dimensionless. `flavor_catalog/processes/beauty/B025.yaml:76`, `flavor_catalog/references/B025/hflav_ckm2025_rdrds.txt:171`, `flavor_catalog/references/B025/hflav_ckm2025_rdrds.txt:472`.

5. NIT: Severity/notes are appropriate: HARD for an observed LFU ratio, and the RS term is clearly marked as a charged-current stress proxy needing human physics, not a complete W'/charged-Higgs/leptoquark calculation. `flavor_catalog_constraints/primary/beauty/B025.py:16`, `flavor_catalog_constraints/primary/beauty/B025.py:29`, `quarkConstraints/semileptonic_lfu.py:21`.

PHYSICS-OK