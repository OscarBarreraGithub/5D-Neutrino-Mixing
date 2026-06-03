# G1 CONSENSUS DESIGN SYNTHESIS (codex, gpt-5.x xhigh). Repo: /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing. Do NOT modify production code. WRITE a single consensus design to `.orchestration/rs_ew_sector_design_CONSENSUS.md` and print a tight summary.

Two INDEPENDENT designs for the "RS electroweak-coupling sector" now exist and largely agree:
- `.orchestration/rs_ew_sector_design.md` (design A, Opus): central insight = nearly every EW proxy is the same object, the RS Z-coupling shift `δg_f/g_Z = (m_Z²/M_KK²)[a_f(c_f) − ⟨a⟩]` from f_IR/f_UV; tiers A(closed-form)/B(numeric oracle); ~24–26 fully rigorous / ~13–16 partial / ~4–6 EDM human.
- `.orchestration/runs/G1-DESIGN/codex_design.md` (design B, codex): central builder computing EW KK spectrum, mass-basis Z/W couplings, heavy neutral contact terms, semileptonic Wilsons; ParameterPoint additions `rs_ew_spectrum`, `rs_ew_couplings`, `rs_semileptonic_wilsons`, `rs_charged_current`, optional `rs_dipole_wilsons`/`rs_higgs_yukawas`; 26 fully rigorous / 17 partial / 7 EDM human.

PRODUCE the consensus design by reading BOTH and reconciling:
1. **Physics inventory** — the unified list of RS-EW quantities + the agreed formula for each (in terms of f_IR/f_UV, ε, c-values, KK masses), with derivations/ references. Where A and B differ in form, pick the one better grounded in `derivations/` and SAY WHY (cite the derivation).
2. **ParameterPoint schema (single agreed set)** — reconcile A's and B's field proposals into ONE concrete set of new fields/extras (names, types, units, real-vs-complex→diagnostics), reusing existing KNOWN_EXTRA_KEYS (`kk_ew_mass_gev`, `lepton_mass_basis_couplings`). Keep frozen + fail-loud.
3. **point_builder logic + new sweep inputs** — the agreed build path; list the NEW sweep inputs (c_L, c_E[3], c_N, M_N, …) the current quark-only builder lacks.
4. **Matching map (single table)** — per affected adapter/constraint family: current proxy → new rigorous input → formula; mark FULLY-rigorous vs PARTIAL vs STILL-HUMAN. Reconcile the 26-vs-(24–26) and 17-vs-(13–16) counts into one agreed count.
5. **Phasing** — one agreed implementation order (sub-phases), each independently testable.
6. **Normalization & convention pin** — the single most important risk both flagged: the g_Z/cosθ_W normalization + sign conventions must match the EXISTING validated `zpole`/`rare_*`/`bsgamma` cores. Specify the exact SM-limit test (universal-c ⇒ all δg=0) and Zbb-shift sign/magnitude cross-check that implementation must pass.
7. **Open disagreements** — explicitly list any point where A and B genuinely differ and you could NOT reconcile from `derivations/`; flag for the human/reviewers to arbitrate (do NOT paper over).

Final message (<=22 lines): the agreed ParameterPoint field set, the agreed fully-rigorous/partial/human counts, the agreed phasing order, the normalization-pin test, and any unresolved disagreements. Confirm the consensus doc path.
