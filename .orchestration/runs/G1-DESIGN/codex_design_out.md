Wrote the independent design to [codex_design.md](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/.orchestration/runs/G1-DESIGN/codex_design.md). No production code was modified.

Headline: centralize RS-EW matching in one builder that computes physical EW KK masses, mass-basis Z/W couplings, heavy neutral contact terms, and semileptonic Wilsons; leave loop dipoles/EDMs honestly partial.

Proposed `ParameterPoint` additions: fill existing `kk_ew_mass_gev` and `lepton_mass_basis_couplings`; add `rs_ew_spectrum`, `rs_ew_couplings`, `rs_semileptonic_wilsons`, `rs_charged_current`, optional `rs_dipole_wilsons`, and optional `rs_higgs_yukawas`.

Named-scope count: 26 constraints get fully rigorous tree-level RS matching, 17 remain partial, 7 EDM constraints stay NEEDS-HUMAN.

Recommended order: derivation notes/profile normalization -> EW spectrum/overlap tables -> quark neutral currents -> lepton neutral currents -> charged-current/G_F shifts -> fermion-KK mixing for Zbb/Higgs LFV -> loop dipoles/EDMs.

Top risks: coupling normalization/signs, double-counting light-Z and heavy-vector exchange, and quark/lepton basis rotations with PMNS/CKM conventions.