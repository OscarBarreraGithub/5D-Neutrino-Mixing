1. **BLOCKER**: Neutral-contact formula is incomplete for LFV/quark-diagonal channels. It needs the full NP product `(g_q^SM δij + δg_q)(g_l^SM δab + δg_l) - SM`; otherwise μ-e conversion, μ→3e, τ→3ℓ miss light-Z `g_q^SM * δg_l^LFV`.

2. **SHOULD-FIX**: The Zbb sign pin is right, but implementation must bypass old proxy helpers. Current proxy wiring gives the opposite sign for positive b-light overlap; the new vector-diagonalized path must enforce IR `b_R` ⇒ `δg_R^b < 0`, `~1e-3` at 3 TeV.

3. **SHOULD-FIX**: Universal-c test wording is too broad. All `δg=0` only if all relevant `a(c)=a_ref` or truly universal across chiral reps; per-species family-universal c removes FCNCs but can leave diagonal Z-pole shifts.

4. **SHOULD-FIX**: Lepton builder contract needs exact shapes: current `compute_all_yukawas()` returns `Y_N_bar` vector and unbarred `Y_N_matrix`; builder must define/store `Y_N_bar_matrix = 2k*Y_N_matrix`, PMNS, and identity `U_e_L/R`.

5. **SHOULD-FIX**: Wilson/contact bundles are sufficient but need typed key schemas for transition, lepton flavor, neutrino flavor, and chirality so full family coverage can be tested.

6. **SOUND**: KK masses, rotations, `zpole` additive `δg`, and rare `C9/C10`/`X_NP` normalizations are consistent if `kk_ew_mass_gev = x1*Lambda_IR` and rare adapters consume Wilson/contact bundles directly.

7. **ARBITRATION**: `a(c)` prefactor RESOLVED as numerical overlap first; no closed form approved until zero-mode/gauge-profile normalization derivation is added.

8. **ARBITRATION**: Classic Zbb RESOLVED as gauge-only PARTIAL; FULL requires fermion-KK EWSB mixing. Buras-Duling-Gori 0905.2318 supports keeping this separate.

9. **ARBITRATION**: Custodial/BKT defaults NEEDS-HUMAN because embeddings, gauge group, and brane terms are model choices.

10. **ARBITRATION**: Dipole-loop normalization NEEDS-HUMAN/LOOP-MATCHER; repo lacks finite KK-loop derivation.

11. **ARBITRATION**: Long-distance/covariance choices NEEDS-HUMAN; they are observable-side choices, not EW matching.

12. **PHASING**: Sound after fixing the contact-formula blocker: derivations → spectrum/overlaps → quark NC → lepton NC → charged current → fermion-KK/Higgs → loops.

Refs: repo `derivations/`; Casagrande et al. 0807.4937; Buras-Duling-Gori 0905.2318.

DESIGN-NEEDS-FIXES.