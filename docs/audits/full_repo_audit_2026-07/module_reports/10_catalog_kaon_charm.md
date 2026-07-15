# Report 10 — Catalog kaon + charm (flavor_catalog_constraints/primary/{kaon,charm}, secondary/kaon)

**Structural summary.** Coverage map: K001=ε_K, K002=Δm_K, K003=ε′/ε (INFO stub), K004=K⁺→π⁺νν̄, K005=K_L→π⁰νν̄, K006=K_L→μμ (SD), K008=K_L→π⁰ee, K009=K_L→π⁰μμ, K010=K_S→π⁰ee, K012=K_S→μμ(SD), K013=K_L→π⁰γγ (INFO), K017=R_K(e/μ), K018=K_l3 V_us (SOFT); K019–K021=LFV kaon; C001=D mixing |M12|, C002=D indirect CPV, C003=ΔA_CP (INFO stub), C004/C005=D⁰→μμ/ee, C006=D⁰→eμ, C007/C008=D⁺→π⁺ℓℓ(′). Deep review: K001/K002/K003/K004/K005/K006/K008 chains, C001–C004, rare_charm_dilepton core.

**Verified correct (numerically spot-checked):** KK-gluon ΔF=2 matching color/Fierz; GGMS M12-ready kaon MEs (chiral 25.7, O4/O1=31.9); ε_K^NP=κ_ε·ImM12/(√2Δm_K); K002 ħ-unit conversion; g_SM² identity; κ₊=5.173e-11, P_c=0.404, X_t=1.481 (SM BR ~8.5e-11); κ_L=2.231e-10; KOTO <2.2e-9 current; κ_μ, Y_eff with correct axial minus sign; D system uses proper charm constants — **no kaon→charm copy-paste bug found**; D⁰→ℓℓ standard structure; C004 uses current 2.1e-9 @90%.

### [MAJOR] K001 ε_K HARD budget is 4.5× looser than the production-core budget that sets the ~7 TeV floor
- **File:** flavor_catalog_constraints/primary/kaon/K001.py:206-231 vs quarkConstraints/deltaf2.py:794-795
- **Category:** statistics / convention-inconsistency
- **Claim:** Catalog K001 veto uses the "loose band edge" budget |ε_exp − (ε_SM − σ_comb)| = **3.04e-4** (σ_SM=1.83e-4 BGS grouped, σ_exp=1.1e-5, SM-choice 1.5e-4 ⇒ σ_comb=2.37e-4), while core `evaluate_epsilon_k` defaults to the bare central residual = **6.7e-5**. Same observable, two HARD verdicts differing ×4.53 ⇒ M_KK floors differing ×√4.5≈2.1 (a 7 TeV core floor becomes ≈3.3 TeV under the K001 budget).
- **Evidence:** `hard_veto_budget=loose_budget` (K001.py:231); core `budget = abs(EPSILON_K_EXP - EPSILON_K_SM)`. Numeric: ratio loose/central = 4.533.
- **Confidence:** high
- **Fix:** One documented ε_K NP budget (and CL) pinned via shared constant + cross-pin test.

### [MAJOR] K001 loose-edge budget is sign-blind: grants the positive-shift room to negative NP shifts too
- **File:** K001.py:211-218 with deltaf2.py:793 (`epsilon_k_np = abs(...)`)
- **Category:** statistics / sign
- **Claim:** Band shifts SM only **down** by σ_comb and compares |ε_NP| to central+σ = 3.04e-4 both directions; for NP that *lowers* ε_K the 1σ-consistent room is ≈1.7e-4, so downward shifts up to ~1.8× too large pass. (K004 implements direction-aware budget — the logic exists but K001 doesn't use it.)
- **Confidence:** high
- **Fix:** Make K001 budget direction-aware on sign(Im M12^NP).

### [MAJOR] Mixed confidence-level conventions across HARD vetoes (68% vs 90% vs 95%)
- **File:** K004.py:327-338 (1σ ≈ 68% CL veto); K001 (1σ loose edge); K005/K006/K008/K009/C004 (90% CL); C002 (95% CL)
- **Category:** statistics
- **Claim:** K004 HARD-vetoes any point whose predicted BR differs from the NA62 central by >1σ_combined (~2.0e-11) — 68% CL — while siblings veto at 90–95%; exclusion floors not statistically commensurable.
- **Confidence:** high
- **Fix:** Standardize HARD vetoes to one CL or record and harmonize.

### [MINOR] K008 v2 hardcoded y7V/y7A coefficients do not reduce to the catalog ISU coefficients in the SM limit
- **File:** quarkConstraints/rare_kaon_dilepton.py:700-704 vs flavor_catalog/processes/kaon/K008.yaml (C_int=6.2, C_dir=2.4)
- **Category:** numerics
- **Claim:** Hardcoded C_DIR=2.67 gives SM-limit effective ISU C_dir 2.66 (+11% vs 2.4); C_INT=8.91 gives 6.50 (+5% vs 6.2); catalog values loaded but only echoed as diagnostics. Init-time KTeV-range check too wide to catch 10%.
- **Confidence:** medium
- **Fix:** Derive C_INT/C_DIR from the YAML ISU coefficients or cite the exact BMS equation.

### [MINOR] K004 HARD anchor is a *preliminary* conference value, demoting the published NA62 observation
- **File:** flavor_catalog/processes/kaon/K004.yaml (latest = Moriond-2026 preliminary 9.6e-11; published JHEP 02(2025)191 13.0e-11 relegated to diagnostic anchor)
- **Category:** stale-data / statistics
- **Claim:** HARD veto centers on unpublished preliminary; with the 1σ budget, the 9.6-vs-13.0 choice (1.2σ) directly moves the pass/fail boundary. Sidecar's own open issue flags this as unresolved.
- **Confidence:** high (policy risk)
- **Fix:** Anchor HARD veto on the published measurement.

### [MINOR] Kaon LR bag parameters at 3 GeV mixed with 2 GeV Wilsons and 2 GeV quark masses (dominant ε_K operator)
- **File:** quarkConstraints/deltaf2.py:641-652,718
- **Category:** convention-inconsistency / units (scale)
- **Claim:** O4_LR combines B4(3 GeV) with C4(2 GeV) and chiral factor at 2 GeV — O(10–20%) internal scheme mismatch on the operator driving the ε_K floor (~5–10% on M_KK). Acknowledged code TODO, but feeds the headline number.
- **Confidence:** high (mismatch), medium (size)
- **Fix:** Common scale for B4/B5, Wilsons, and masses.

### [NOTE] K005 relies on the KOTO direct limit; the Grossman-Nir bound implied by NA62 is ~5× stronger (GN with BR(K⁺)=9.6e-11 gives ≈4.1e-10 vs KOTO 2.2e-9). Flagged only in YAML open issues. Fix: add GN diagnostic/budget.
### [NOTE] C001 Δm_D catalog anchor (6.562e-15 GeV) vs stale core default 6.25e-15 — C001 overrides correctly; core paths (evaluate_d0_mixing without override, e.g. ΔF=2 bundle in run_full_catalog_scan) are ~5% tighter and inconsistent.
### [NOTE] C003 ΔA_CP (LHCb 5.3σ observation) is a non-vetoing INFO stub — loaded correctly, not faked as null, but direct-CP RS contribution effectively unconstrained; honest and documented coverage gap.

### Cross-module inconsistencies
1. ε_K NP budget: catalog K001 3.04e-4 vs core 6.7e-5 (×4.5) — single largest driver of "which floor is real".
2. Δm_D: core 6.25e-15 vs catalog 6.562e-15.
3. Δm_K: K002 consistent with core (no issue).
4. B4/B5 scale caveat shared by kaon, B, D paths via byte-identical _meson_matrix_elements — fix must be uniform.

### Not reviewed
K009/K010/K012/K013/K017/K018 beyond anchor greps (anchors look current); secondary kaon K019–K021; charm C005–C008 internals; Phase-3a/4a rs_semileptonic_wilsons C10→Y/X mapping consumed by K004–K012; evolve_deltaf2_wilsons magnitudes (sibling B-sector over-running finding likely propagates to C4^K/C4^D); tests.
