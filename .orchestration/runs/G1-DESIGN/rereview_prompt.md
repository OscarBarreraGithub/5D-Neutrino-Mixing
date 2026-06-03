# RE-REVIEW of the REVISED G1 consensus design (codex, gpt-5.x xhigh). Repo: /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing. Do NOT modify code. Read `.orchestration/rs_ew_sector_design_CONSENSUS.md` (now revised). The prior dual review (codex DESIGN-NEEDS-FIXES + Opus DESIGN-OK) raised the items below; verify EACH is now resolved and the design is sound, complete, and implementable. This design will be implemented next.

VERIFY RESOLVED:
1. **BLOCKER (was: incomplete neutral-contact formula)** — confirm the contact coefficient now uses the FULL product-minus-SM including the light-Z LFV cross-term `g_q^SM·δg_l^LFV`. The doc states: `C_AB^Z = (g_Z²/m_Z²)[(g_qA^SM δ_ij + δg_qA_ij)(g_lB^SM δ_ab + δg_lB_ab) − g_qA^SM g_lB^SM δ_ij δ_ab] + Σ_V g_Vq g_Vl/M_V²`. Confirm this is correct (it should reproduce: pure SM when δg=0; LFV (a≠b) driven by g_q^SM·δg_l; quark-FCNC (i≠j) by δg_q·g_l^SM) and that the matching-map entries for μ-e conversion / μ→3e / τ→3ℓ now consume it.
2. **No double-counting** — confirm the doc now states rare Wilsons do NOT reuse `_wilson_prefactor` (contacts are GeV^-2 already carrying 1/m_Z² or 1/M_V²).
3. **KNOWN_EXTRA_KEYS** — confirm the doc requires adding all new extras there (fail-loud make_point).
4. **Zbb sign test** — confirm the required test: vector-diagonalized path ⇒ IR b_R ⇒ δg_R^b < 0, ~1e-3 at M_KK~3 TeV (not via old proxy helpers).
5. **universal-c wording; lepton-builder shapes (Y_N_bar_matrix=2k·Y_N_matrix, PMNS, identity U_e); typed bundle schemas** — confirm all corrected.
6. **Arbitrations baked as firm scope** — a(c) numerical-first (closed form deferred); classic Zbb gauge-only PARTIAL (phase 6); custodial/BKT = minimal-RS default + diagnostic + HUMAN-flagged; dipole-loop = HUMAN, dipoles stay PARTIAL; long-distance = out-of-sector diagnostics. Confirm the doc has a "HUMAN-INPUT items" section.
7. **Implementability** — is the design now concrete enough to implement phase-by-phase (field set, point_builder logic, sweep inputs, matching map, per-phase tests)?

OUTPUT (<=16 lines): for each item RESOLVED / STILL-OPEN (with evidence); any NEW findings; final implementability call. End with: DESIGN-OK or DESIGN-NEEDS-FIXES.
