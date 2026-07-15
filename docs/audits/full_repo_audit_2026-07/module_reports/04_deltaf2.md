# Report 04 — ΔF=2 machinery (quarkConstraints/deltaf2.py + review_local + docs/audits trail)

**Structural summary.** The core chain — KK-gluon tree matching (CFW 0804.1954 color factors {1/6, 1/6, −1, +1/3}), LO BMU RG running with 6→5→4 flavor thresholds, GGMS Eq. (8) M12-ready matrix elements, ε_K master formula with κ_ε/√2/Δm_K — is **correct in the current code**. Verified numerically: the BMU→scalar ADM similarity transform ([[2,0],[12,−16]] → [[−16,−6],[0,2]]) is exact; C4(3 TeV→2 GeV) = 3.538, C_VLL = 0.729 (right directions); the SM-box anchor with the code's ⟨O1⟩ = (1/3)f²mB̂ reproduces ε_K ≈ 1.7×10⁻³ (~2.2×10⁻³ target within rough CKM inputs — the legacy 2/3 would overshoot 2×); Δm = 2|M12| and the 1/(2m_M) each appear exactly once. The serious problems: **the entire human-review documentation layer records the pre-fix (swapped, ×2) matrix elements as verified physics**, a B-meson RG endpoint below m_b, and statistically inconsistent NP-budget conventions across systems.

### [MAJOR] Review documents record the old swapped/×2 matrix elements as confirmed physics
- **File:** review_local/constraint_formulas.tex:111–113,129; review_local/epsilon_k_review.tex:693–697; review_local/deltaF2_framework_review.tex:648–652; review_local/constraint_formulas_review.md (finding 3)
- **Category:** doc-code-mismatch
- **Claim:** Every review/tex document quotes ⟨O1⟩=(2/3)f²mB₁, ⟨O4⟩=(r/6+1/4)f²mB₄, ⟨O5⟩=(r/2+1/12)f²mB₅ *while claiming M12 normalization* — these are the pre-B3 buggy values (Q4/Q5 chiral coefficients swapped AND ×2), which the code has since fixed to GGMS hep-ph/9604387 Eq. (8): (1/3), (r/4+1/24), (r/12+1/8).
- **Evidence:** Code (deltaf2.py:726-728): `o4_lr = (m_ratio_sq*(1/4) + 1/24)…`; comment: "The previous code had these swapped AND each x2 too large". Tex: `\langle Q_4\rangle = (1/6 r_χ + 1/4) f² m B₄` with line 129 "already normalised for M12". Note (r/6+1/4) = 2×(r/12+1/8) — the tex O4 is exactly twice the code's O5. The human review md **"CONFIRMED as a standard parametrization"** — endorsing the wrong pairing; the dominant C4·⟨O4⟩ term differs by ≈0.70 (r=25.7). Tests (test_epsilon_k_physics.py:205-254) pin the correct GGMS rationals with an independent oracle.
- **Confidence:** high
- **Fix:** Rewrite the ME equations in all four review documents to the GGMS values and retract the review md's "CONFIRMED" on item 3.

### [MAJOR] B_d/B_s Wilson coefficients run below m_b while bags are quoted at μ = m_b
- **File:** quarkConstraints/deltaf2.py:464 (`mu_had: float = 2.0` applied to all systems), 663–684
- **Category:** wrong-formula (RG endpoint)
- **Claim:** B-meson Wilsons are evolved to 2 GeV (integrating out the b, a valence quark of the meson) and contracted with HPQCD/FLAG bags defined at μ = m_b, spuriously inflating the B-system LR amplitudes.
- **Evidence:** Measured with the repo's own evolver: C4(2 GeV)/C4(m_b) = **1.2496**, C_VLL(2 GeV)/C_VLL(m_b) = 0.946. So B_d/B_s C4 contributions are ~25% too large. docs/audits/wilson_rg_inventory.md item 4 acknowledges and explicitly defers this ("would invalidate all scan outputs") — documented, but it is a genuine physics error in production numbers. Same class: kaon B4/B5 are FLAG MS-bar(3 GeV) contracted at 2 GeV; r_χ uses 2 GeV masses with 3 GeV bags.
- **Confidence:** high
- **Fix:** Per-system μ_had (m_b for B, 3 GeV for D and kaon B4/B5, with r_χ masses at matching scale).

### [MAJOR] NP-budget statistics are inconsistent across systems: ε_K maximally tight, Δm_B/D maximally loose
- **File:** quarkConstraints/deltaf2.py:793–802, 956–968
- **Category:** statistics
- **Claim:** ε_K is bounded by the bare central-value gap |exp−SM| = 6.7×10⁻⁵ (3% of ε_K^exp; smaller than BGS's own 1σ = 1.8×10⁻⁴, no uncertainty propagated), while B_d/B_s use `max(Δm_exp/2, |Δm_exp−Δm_SM|/2)` = Δm_exp/2, i.e. NP allowed to be 100% of the measured mass difference (Δm_d agrees with SM to ~8%). A UTfit-style consistent treatment (~30–40% NP fraction) would loosen ε_K ~5–10× and tighten Δm_B ~5–10×.
- **Evidence:** `budget = abs(EPSILON_K_EXP - EPSILON_K_SM)` vs `_bd_budget: max(DELTA_M_BD_EXP/2, abs(DELTA_M_BD_EXP - DELTA_M_BD_SM)/2)` — the `max()` always selects the loosest option. docs/audits/epsilon_k_sm_decision.md mandates the ε_K floor be quoted as a **band** (budget 1×10⁻⁵–3×10⁻⁴, ~0.39×–2.1× on M_KK); the code default is central-only, so any single-number "lane B ~7 TeV" floor silently drops the mandated band.
- **Confidence:** high
- **Fix:** Propagate σ(SM)⊕σ(exp) into all budgets with one declared convention; always quote the ε_K band per the decision file.

### [MINOR] LO-only running understates the LR enhancement (η₄ = 3.54 vs NLO ~4.5–5)
- **File:** quarkConstraints/qcd_running.py (whole module)
- **Category:** numerics
- **Claim:** LO magic number C4(2 GeV)/C4(3 TeV) = 3.538; NLO (Ciuchini/BJW hep-ph/0102316) gives ~4.5–5, so the ε_K floor is systematically **under**-estimated by ~10–15% in M_KK.
- **Confidence:** high (LO value), medium (NLO gap size)
- **Fix:** Adopt NLO magic numbers or state the one-sided bias in the methodology note.

### [MINOR] VIA additive terms (1/24, 1/8) likely double-count against FLAG/HPQCD chiral-only bag definitions
- **File:** quarkConstraints/deltaf2.py:727–728, 918–919
- **Category:** convention-inconsistency
- **Claim:** FLAG 2024 kaon BSM bags (and lattice B-meson bags) are defined with the pure chirally-enhanced VIA denominator, so multiplying B₄/B₅ by the *full* GGMS bracket (r/4 **+1/24**, r/12 **+1/8**) overcounts the non-chiral piece: ~0.7% (O4_K), ~5.8% (O5_K), but ~48% of the subdominant ⟨O5⟩ for B mesons, a few % of |M12^B|.
- **Confidence:** medium
- **Fix:** Match each bag to its collaboration's exact denominator; drop the additive term where the lattice definition omits it.

### [MINOR] Stale hadronic inputs kept against the audit's own canonical values
- **File:** quarkConstraints/deltaf2.py:655, 683, 694–697
- **Category:** stale-data
- **Claim:** `DELTA_M_D_EXP = 6.25e-15` vs the audit's canonical 6.56×10⁻¹⁵ GeV (4.8% tightening kept in code, anti-conservative); `DELTA_M_BS_EXP = 1.1688e-11` vs PDG 17.765 ps⁻¹ = 1.1693×10⁻¹¹ (negligible); `EPSILON_K_SM = 2.161e-3` carries an unpublished digit beyond BGS's 2.16 — review FIX 4 never applied; `B_4_D = B_5_D = 1.0` vs ETM 2015's 0.91/0.97 (flagged Yellow, unfixed).
- **Confidence:** high
- **Fix:** Adopt the audit's canonical column values.

### [NOTE] Δm_K (Re M12) is not in the default constraint bundle
- **File:** quarkConstraints/deltaf2.py:208–265, 349–453
- **Category:** logic-bug (coverage gap)
- **Claim:** `evaluate_delta_mk` exists but `DEFAULT_DELTA_F2_INPUTS_V1` gates kaons only through ε_K, so CP-conserving s–d NP is unconstrained in production scans.
- **Confidence:** high
- **Fix:** Add a `delta_m_k` input row using the existing evaluator.

### [NOTE] `bound` field misleading for ε_K; legacy fallback mixes evolved Wilsons with unevolved bounds
- **File:** quarkConstraints/deltaf2.py:216 (`bound=2.0e-8`), 139–141, 543–561
- **Category:** logic-bug / convention-inconsistency
- **Claim:** For the kaon row, `DeltaF2ObservableSummary.bound` returns the legacy surrogate 2×10⁻⁸ while `ratio_to_bound` is computed against the 6.7×10⁻⁵ budget; and the `use_hadronic=False` fallback applies pre-running-calibrated `LEGACY_OPERATOR_WEIGHT_BOUNDS` to RG-*evolved* coefficients when `apply_qcd_running=True` (default).
- **Confidence:** high
- **Fix:** Populate `bound` with the hadronic budget per system; force `apply_qcd_running=False` in the legacy path.

### [NOTE] Audit-trail documents no longer match the code they certify
- **File:** docs/audits/bag_param_inventory.md:57–59; docs/audits/wilson_rg_inventory.md:46–54
- **Category:** doc-code-mismatch
- **Claim:** The audit trail records the pre-B3 ME prefactors as inspected/OK, and cites line numbers that have shifted; nothing marks these sections superseded by the B2/B3 fixes.
- **Confidence:** high
- **Fix:** Add a post-B3 addendum to both audit files.

**Verified-correct (no action):** CFW tree color/Fierz factors; γ_VLL=4 and scalar-LR ADM (re-derived, exact); running direction; threshold sequence incl. m_t (real pre-audit bug, correctly fixed); κ_ε=0.94, √2, Δm_K; f_K/f_Bd/f_Bs/f_D; single 1/(2m_M); ε_K phase convention-stable via PDG SVD rephasing; Δm=2|M12| exactly once.

**Cross-module inconsistencies:** (1) ε_K budget (central-gap, no σ) vs B/D budgets (full Δm/2) — opposite statistical philosophies in one bundle; (2) review_local/tex + docs/audits describe pre-B2/B3 code while tests pin post-fix values — the paper-facing documentation layer is the stale one; (3) default `g_s_star=None` gives g_s(M_KK)≈1.0, not the RS composite g_s*≈3–6 the lane-A/C literature assumes; (4) kaon r_χ at 2 GeV vs B4/B5 at 3 GeV vs B-meson bags at m_b, all contracted at one global μ=2 GeV.

**Not reviewed:** paper_0710_1869/ internals; test_qcd_running/test_wilson_rg_audit expected values; delta_m_s_review.tex, d0_mixing_review.tex beyond grep-level; fit.py rephasing algebra line-by-line; D0 long-distance tex discussion.
