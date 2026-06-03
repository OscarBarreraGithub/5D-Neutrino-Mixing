# W2 PHASE 5 — SUB-STEP 5b: charged-current rewire (codex, gpt-5.x xhigh). Repo: /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing. Implement sub-step 5b of the DUAL-APPROVED plan `.orchestration/runs/W2-P5/plan.md`, using the 5a builder (committed 2b89ae2: `rs_charged_current` extra with `epsilon` shifts, `delta_g_W_ud_L`, `delta_G_F_over_G_F`) and the established rewire+degradation pattern. First a SHORT plan, then implement. Codex + Opus dual-review (both must APPROVE). SCOPE = EW002, K018, K017, B009, B025 (EW003 = 5c separate).

REWIRE (charged-current; consume the rigorous `epsilon`/`delta_G_F` shifts, NOT proxies):
- **EW002** (CKM first-row unitarity): `Delta_CKM_NP ≈ 2 Σ|V_ij|² Re(epsilon_ij)` vs the YAML first-row sum/budget. Becomes RIGOROUS but stays SOFT (non-veto).
- **K018** (|Vus| from K_l3): `|Vus|_app = |Vus|·|1+epsilon_us^l|` (or linear `δ|Vus|/|Vus|=Re epsilon_us^l`). Full for minimal LH W/W'; radiative/isospin/mode-weight stays PARTIAL if not available.
- **K017** (R_K = K→eν/μν LFU ratio — NOT |Vus|): `R_K = R_K^SM |1+epsilon_us^e|²/|1+epsilon_us^μ|²`. Full for minimal LH LFU; charged-Higgs/heavy-ν/radiative stay NEEDS-HUMAN.
- **B009** (B→τν): `BR = BR_SM |1+epsilon_ub^τ|²`, YAML f_B/|Vub| anchors unchanged. Full for minimal LH; charged-Higgs/RH/scalar PARTIAL.
- **B025** (R(D)/R(D*) LFU): `R(D) = R_SM |1+epsilon_cb^τ|²/|1+epsilon_cb^light|²`. STAYS PARTIAL — vector ratio rigorous, but scalar/RH WET + form-factor integration NOT built (do NOT claim FULL).
- Replace the old proxies (m_K²/M_KK², m_B²/M_KK², m_b m_τ/M_KK²). SM-vs-data PULL + anchors UNCHANGED.

GRACEFUL DEGRADATION: rigorous `rs_charged_current` when present; ABSENT ⇒ non-vetoing evaluated=False + missing_extra (no crash/fake pass).

v1 EXPECTATION (verify in tests): universal c_L ⇒ lepton/G_F pieces absorbed ⇒ EW002/K017/K018 first-row & LFU shifts tiny/near-SM; B025 LFU vector cancels universal pieces; third-family quark-doublet residuals ~1e-4–1e-3. So these mostly stay near-SM in v1 (not big vetoes) — that's correct, not a failure.

TESTS: rigorous-path (point via `build_from_rs_ew_inputs(...,include_charged_current=True)` ⇒ the epsilon shift enters; cross-check independent of adapter); SM-limit (universal-c ⇒ epsilon=0 ⇒ recover committed SM); absent-path; B025 PARTIAL flags preserved; EW002 stays SOFT. Replace proxy-only tests — enumerate, no silent coverage loss. `python -m pytest tests/ -q` stays green.

CONSTRAINTS: physics ONLY via adapters; ConstraintResult numeric fields real finite floats; touch ONLY EW002/K018/K017/B009/B025 + their adapters/tests.

OUTPUT (<=16 lines): short plan; rewired files + IDs; one rigorous epsilon→observable example + SM-limit recovery + absent-path; v1 near-SM confirmation; B025-PARTIAL + EW002-SOFT preserved; old proxies removed; test-count change; pytest counts. End with: P5B-DONE.
