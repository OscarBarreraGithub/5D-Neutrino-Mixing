1. NIT: Amplitude verdict PASS. B004 uses complex `M12_Bs^NP` and `phi_s_np = phase(1 + M12_NP/M12_SM)`, not `|M12|`; `bs_mixing_from_wilsons_with_running` is diagnostics only. `flavor_catalog_constraints/primary/beauty/B004.py:379-394`.

2. NIT: QCD-running verdict PASS. The verdict path calls `bs_mixing_m12_np_from_wilsons_with_running`, whose helper evolves Wilsons before `compute_m12_np`. `B004.py:380-383`, `physics_adapters/deltaf2.py:418-443`. Representative `L_sb=1e-3(1+i)`: running changes `|M12|` from `6.115e-15` to `4.458e-15 GeV` and `phi_s^NP` from `1.045e-3` to `7.621e-4 rad` (`0.729x`).

3. NIT: Budget verdict PASS. Code compares total `phi_s = phi_s^SM + phi_s^NP` to HFLAV with combined exp+SM phase uncertainty: `sqrt(0.016^2+0.0009^2)=0.016025 rad`, lower side `0.016011 rad`; central SM-exp residual is `0.0042 rad`. This is defensible for the requested phase observable. `B004.py:250-263`, `B004.yaml:60-67`, `B004.yaml:109-115`.

4. NIT: Anchor verdict PASS. YAML values match snapshots: HFLAV all-combined `phi_s^ccs=-0.041 +/- 0.016 rad`, mode-specific `-0.050 +/- 0.017 rad`, LHCb `-0.039 +/- 0.022 +/- 0.006 rad`, SM `-2 beta_s=-0.0368^{+0.0009}_{-0.0006} rad`. `flavor_catalog/references/B004/hflav_pdg2025_phis_inputs.txt:23-27`, `flavor_catalog/references/B004/lhcb_2024_jpsikk_phis_arxiv.txt:22-33`.

5. NIT: Asymmetric-SM uncertainty side is likely reversed in the “direction-aware” selector: for positive residual `predicted-exp > 0`, the theory uncertainty that reduces tension is the SM lower error `0.0006`, not upper `0.0009`. Effect is tiny (`0.016025` vs `0.016011 rad`) and will not matter practically. `B004.py:332-344`.

PHYSICS-OK