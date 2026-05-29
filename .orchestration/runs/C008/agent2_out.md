1. NIT: Δm/CP checklist is N/A for C008: this is a ΔC=1 LFV branching-ratio bound, not an M12 observable. Code uses a positive rate from `|C9|^2 + |C10|^2`, not Re/Im M12; see [quarkConstraints/rare_charm_lfv_semileptonic.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/quarkConstraints/rare_charm_lfv_semileptonic.py:373).

2. NIT: QCD-running check is N/A here. C008 imports/calls the rare-charm LFV semileptonic adapter, not `deltaf2`; no non-running ΔF=2 path is used for the verdict. See [C008.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/primary/charm/C008.py:50) and [C008.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/primary/charm/C008.py:553). Running effect: not defined for this v1 proxy.

3. NIT: Amplitude structure is defensible: P→P semileptonic uses `C9_LFV+C9p_LFV` and `C10_LFV+C10p_LFV`, matching the C007 convention; equal-mass cross-check reproduces C007 to numerical precision. See [rare_charm_lfv_semileptonic.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/quarkConstraints/rare_charm_lfv_semileptonic.py:232).

4. NIT: Budget is correct for a pure-NP LFV upper limit: stricter per-charge 90% CL limit `min(2.1e-7, 2.2e-7)=2.1e-7`, not a central residual. Code compares the one-mode prediction to both limits via max ratio. See [C008.yaml](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog/processes/charm/C008.yaml:70) and [C008.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/primary/charm/C008.py:572).

5. NIT: Anchor numbers match snapshots/audit: PDG/LHCb `D+→π+e+μ- <2.1e-7`, `D+→π+e-μ+ <2.2e-7`; 95% companions `2.3e-7`, `2.2e-7`; BaBar predecessor `2.9e-6`, `3.6e-6`. See [factcheck_charm.md](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog/audits/factcheck_charm.md:121).

6. NIT: Severity/units/diagnostics look physically consistent: HARD is appropriate for observed LFV upper limits; masses/q2 are in GeV/GeV² and output is dimensionless BR; `NEEDS-HUMAN-PHYSICS` caveats correctly flag missing LFV lepton matching, scalar/tensor, resonance, and acceptance effects. See [C008.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/primary/charm/C008.py:603).

PHYSICS-OK