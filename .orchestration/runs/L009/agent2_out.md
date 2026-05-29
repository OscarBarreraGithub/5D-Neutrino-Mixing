1. SHOULD-FIX — L009 is flavor-pinned in labels but not in accepted spurion keys: `left_emu_overlap`/`right_emu_overlap` from the shared L002 core are accepted and interpreted as τμ. Probe: `left_emu_overlap=1e-2, m_kk=3000 GeV` gives `BR(τ→3μ)=6.76e-11`. Correct physics: τ→3μ needs μτ/τμ LFV input or a flavor matrix element, not eμ aliases. `flavor_catalog_constraints/physics_adapters/lfv_three_body_tau.py:141`, `quarkConstraints/lfv_three_body.py:593`.

2. NIT — Amplitude check: Δm/M12 real-vs-imag logic is not applicable to this charged-LFV branching ratio. L009 correctly compares the scalar pure-NP BR to the limit. `quarkConstraints/lfv_three_body.py:388`, `flavor_catalog_constraints/primary/charged_lepton/L009.py:232`.

3. NIT — QCD running check: not applicable for a purely leptonic τ→3μ BR; no hadronic Wilson running to `mu_had=2 GeV` should enter. L009 uses the LFV three-body module, not `deltaf2`. `flavor_catalog_constraints/primary/charged_lepton/L009.py:37`.

4. NIT — Budget is physically defensible: SM set to 0 and budget is the 90% CL BR upper limit, `1.9e-8`, not a central residual. `flavor_catalog/processes/charged_lepton/L009.yaml:73`, `flavor_catalog_constraints/primary/charged_lepton/L009.py:69`.

5. NIT — Anchor/severity/units check passes: PDG/Belle II `1.9e-8` at 90% CL and post-PDG LHCb `1.9e-8` at 90% CL, `2.3e-8` at 95% CL match snapshots; severity is HARD; masses are in GeV and the τ/μ dipole factor is `0.00224135`. `L009.yaml:78`, `L009.yaml:105`, `quarkConstraints/lfv_three_body.py:69`.

PHYSICS-NEEDS-FIXES