Wrote the consensus design here: [.orchestration/rs_ew_sector_design_CONSENSUS.md](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/.orchestration/rs_ew_sector_design_CONSENSUS.md)

Summary:
- Agreed extras: `kk_ew_mass_gev`, `lepton_mass_basis_couplings`, `rs_ew_spectrum`, `rs_ew_couplings`, `rs_semileptonic_wilsons`, `rs_charged_current`, optional `rs_dipole_wilsons`, optional `rs_higgs_yukawas`.
- Folded A’s `lepton_overlaps` into `lepton_mass_basis_couplings` to avoid duplicate extras.
- Agreed count: `26` fully rigorous after fermion-KK Zbb promotion, `17` partial, `7` still human/EDM. Before the Zbb fermion-mixing subphase, full is `24`.
- Phasing: derivation pins -> EW spectrum/overlaps -> quark neutral currents -> lepton neutral currents -> charged current/G_F -> fermion-KK + Higgs -> loop dipoles/EDMs.
- Normalization pin: stored `z_delta_g_*` is dimensionless additive in existing `zpole` convention; rare-decay Wilsons use physical contacts/Wilson bundles to avoid `g_Z` double counting.
- Required tests: universal-c/a_ref subtraction gives all `delta g = 0` and existing SM-limit adapter outputs; IR `b_R` gives negative `delta g_R^b ~ 1e-3` at `M_KK ~ 3 TeV`.
- Unresolved for reviewers: exact closed-form `a(c)` prefactor, complete classic Zbb treatment, custodial/BKT defaults, finite dipole loop normalization, and observable-side long-distance/covariance choices.