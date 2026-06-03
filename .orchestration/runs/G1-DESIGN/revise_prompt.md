# G1 DESIGN REVISION (codex, gpt-5.x xhigh). Repo: /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing. Do NOT modify production code. REVISE `.orchestration/rs_ew_sector_design_CONSENSUS.md` in place to incorporate the dual-review findings below, then print a tight changelog. The revision will be re-reviewed by codex + Opus (both must APPROVE).

The dual review returned codex DESIGN-NEEDS-FIXES + Opus DESIGN-OK-with-2-SHOULD-FIX. Apply ALL of:

**BLOCKER (codex #1) — full neutral-contact NP product.** The neutral-current contact/coupling formula must use the FULL product `(g_q^SM·δ_ij + δg_q)(g_l^SM·δ_ab + δg_l) − SM`, NOT just the δg×δg or δg×g_SM piece. Otherwise LFV channels (μ-e conversion, μ→3e, τ→3ℓ) miss the light-Z `g_q^SM · δg_l^LFV` cross term. Rewrite the relevant formula + matching-map entries accordingly.

**SHOULD-FIX (Opus a) — KNOWN_EXTRA_KEYS.** State explicitly that the new extras (`rs_ew_spectrum`, `rs_ew_couplings`, `rs_semileptonic_wilsons`, `rs_charged_current`, optional `rs_dipole_wilsons`/`rs_higgs_yukawas`) MUST be added to `KNOWN_EXTRA_KEYS` (fail-loud `make_point` requires it).

**SHOULD-FIX (Opus b + codex #6) — Wilson prefactor double-count.** State explicitly that the new rare-decay path does NOT reuse the existing `_wilson_prefactor` (which already bakes in 1/M_KK²); the new contact coefficients are GeV^-2 and carry their own 1/M_V². Spell out how C9_NP/C10_NP/X_NP plug into the existing `rare_*` cores WITHOUT a second 1/M_KK².

**SHOULD-FIX (codex #2) — Zbb sign via new path.** Implementation must bypass the old proxy helpers; the vector-diagonalized path must enforce IR b_R ⇒ δg_R^b < 0, ~1e-3 at M_KK~3 TeV. Encode this as a required test.

**SHOULD-FIX (codex #3) — universal-c test wording.** Tighten: all δg=0 only if all relevant a(c)=a_ref / truly universal across chiral reps; per-species family-universal c removes FCNCs but can leave diagonal Z-pole shifts. Fix the test statement.

**SHOULD-FIX (codex #4) — lepton builder shapes.** Specify exact shapes/contract: `compute_all_yukawas()` returns `Y_N_bar` (vector) and unbarred `Y_N_matrix`; the builder must define/store `Y_N_bar_matrix = 2k·Y_N_matrix`, PMNS, and identity `U_e_L/R` (charged-lepton mass basis).

**SHOULD-FIX (codex #5) — typed bundle schemas.** Define typed key schemas for the Wilson/contact bundles (transition, lepton flavor, neutrino flavor, chirality) so full family coverage is testable.

**ARBITRATIONS — bake these decisions into the doc as firm scope:**
- `a(c)` prefactor: NUMERICAL overlap (Tier-B) FIRST; NO closed form until a zero-mode/gauge-profile normalization derivation is added to `derivations/`. Mark the closed form as deferred.
- Classic Zbb: gauge-only is PARTIAL now; FULL needs fermion-KK EWSB mixing (phase 6). Keep separate.
- **Custodial/BKT defaults: NEEDS-HUMAN (model choice — embeddings/gauge group/brane terms).** Use a documented "minimal-RS" default + diagnostic; flag the custodial-protection variant as a human decision. Do NOT silently pick one.
- **Dipole-loop normalization: NEEDS-HUMAN / loop-matcher.** Repo lacks a finite KK-loop derivation → dipoles stay PARTIAL/absent, NOT promoted.
- Long-distance/covariance: observable-side, out-of-sector → diagnostics only.

Also add a short "HUMAN-INPUT ITEMS surfaced by design review" section listing the custodial/BKT choice + dipole-loop + EDM items so the orchestrator can relay them.

Final message (<=18 lines): changelog of what you changed (by finding), the corrected contact formula, the firm scope after arbitration (counts: fully-rigorous/partial/human), and the explicit list of HUMAN-INPUT items. Confirm the doc path.
