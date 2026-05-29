1. NIT: ΔF=2 amplitude check is N/A for B021: this is a ΔB=1 BR, not \(M_{12}\)/CP; verdict uses \(|BR_pred-BR_exp|/budget\), not Re/Im \(M_{12}\): [B021.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/primary/beauty/B021.py:561).

2. NIT: QCD-running-to-2-GeV check is N/A here; B021 intentionally uses the rare \(b\to s\mu\mu\) C9/C10 proxy at matching scale via `compute_rare_b_dilepton_wilsons`, not DeltaF2 `*_with_running`: [rare_b_baryon_dilepton.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/quarkConstraints/rare_b_baryon_dilepton.py:444).

3. NIT: Anchors match YAML/snapshots: LHCb \(1.18^{+0.286}_{-0.283}\times10^{-7}\,\mathrm{GeV}^{-2}\) over \(15<q^2<20\) gives \(5.90^{+1.43}_{-1.42}\times10^{-7}\); PDG total \(1.08\pm0.28\times10^{-6}\): [B021.yaml](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog/processes/beauty/B021.yaml:88), [B021.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/primary/beauty/B021.py:436).

4. NIT: SM bin validation is good: implemented SM \(BR(15,20)=5.5324\times10^{-7}\), average \(1.1065\times10^{-7}\,\mathrm{GeV}^{-2}\), within LHCb/PDG high-\(q^2\) bin; tests confirm: [test_B021.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/tests/constraints/primary/beauty/test_B021.py:322).

5. NIT: Budget is uncertainty-aware, not bare residual: \(3.6757\times10^{-8}+1.4309\times10^{-7}+0.30\times5.5324\times10^{-7}=3.4582\times10^{-7}\); defensible for HARD proxy: [B021.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/primary/beauty/B021.py:494).

6. NIT: Form-factor inputs match Detmold-Lin-Meinel-Wingate \(N_\pm,X_\pm\) values, and limitations/C7/charm/covariance are documented; tree-level `c_gamma=1,c_v=0` should be read as proxy, not full Detmold SM. Source: https://arxiv.org/abs/1212.4827; [rare_b_baryon_dilepton.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/quarkConstraints/rare_b_baryon_dilepton.py:52).

PHYSICS-OK