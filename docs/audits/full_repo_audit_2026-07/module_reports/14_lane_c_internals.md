# Report 14 — Gap-closing: Lane C unread internals (paper_0710_1869: kkgluon, couplings, model, fit, hadronic, observables, artifacts, verifier, tests)

**Structural summary.** ~80% frozen-ID/schema boilerplate; physics in model.py (spurions, V5KM, QS1 seed→profile map), kkgluon.py+couplings.py (couplings), fit.py (mass/CKM probes), eft_deltaf2/{hadronic,observables}.py. KK-gluon coupling chain is exactly g_s(μ_gs)=√(4πα_s) × U†diag(f²)U — suspected missing RS volume factor **CONFIRMED**. Two new independent problems: frozen QS1 affine map c=+1·λ+0 inverts the RS localization hierarchy (Table I provably unreachable), and matching-vs-hadronic operator normalizations clash by uniform ×4. ε_K is never computed anywhere in Lane C (observables stop at M12/Δm NP-only). Validation fully circular.

### [CRITICAL] CONFIRMED: no √(2πkr_c) volume enhancement anywhere in the Lane-C coupling chain — and the frozen contract forbids adding it
- **File:** paper_0710_1869/couplings.py:241-243, 271-278; kkgluon.py:809-859
- **Category:** wrong-formula (factor ~6–8.5 per coupling)
- **Claim:** Flavor matrices are g_s(μ_gs)·U†diag(f²)U with g_s=√(4πα_s(3 TeV))≈1.0005, no volume factor; physical RS coupling ≈ g_s√(2πkr_c)·f_if_j (√(2·35.94)=8.48), so |C1|,|C4| low by 2πkr_c ≈ 36–72, inferred ε_K floor low ~6–8×.
- **Evidence:** Executed: α_s(μ_gs)=0.0797, g_s=1.0005; couplings.py:241-243 **raises** unless g_s_mu_gs == √(4πα_s) — contract id `explicit_mu_gs.g_s_sqrt_4pi_alpha_s.v1` structurally excludes an enhanced g*. No kr_c/volume symbol in the package. Lane B couplings.py:115-124 documents exactly the required g_s_star = g_s·√(2kπr_c) ≈ 6 (CFW). Universal trace-subtraction is an identity shift, not the −g_s/√(2πkr_c) UV piece.
- **Confidence:** high (executed; contract-level exclusion verified)
- **Fix:** Multiply the f²-kernel by √(2πkr_c) and add the −g_s/√(2πkr_c) universal term; unfreeze the g_s contract.
- **NOTE:** the ~35× deficit in |C4| partially cancels against errors the other way (ME ×4 too big, m_g1=Λ makes 1/M² 6× too big, inverted RG enhances C1 1.37× instead of 0.73×) — defaults are accidentally within O(1–10) of right, for four wrong reasons.

### [CRITICAL] Frozen QS1 seed→profile map (c = +1.0·λ + 0.0, ascending) inverts the RS hierarchy; Table I provably unreachable
- **File:** model.py:1592-1609, 280-312; frozen at inputs.py:434-437 and model.py:315-331
- **Category:** sign / logic-bug
- **Claim:** Physical branch sets c_i equal to ascending spurion eigenvalues, so largest eigenvalue (top) gets largest c → smallest f_IR — opposite of FPR (larger Yukawa ⇒ more IR-localized); since c_u = eig(Y†Y) ≥ 0, Table I's c_u3 = −0.06 is unreachable for any seed.
- **Evidence:** Executed with Table-I-like seeds: c_Q=(0.003,0.162,0.642), F_Q=(0.705,0.581,0.0023) — first generation maximally IR-coupled, top UV-localized; probe masses m_up=(0.002, 8.65, 47.2) GeV, |V_cb|≈0.28. Affine structure needs a>1/2, b<0 to map ascending λ to descending c; frozen policy deletes both parameters. Contradicts the package's own fit.py:269-288 branch, which attaches c-descending onto λ-ascending slots (correct orientation).
- **Confidence:** high (executed; positivity argument exact)
- **Fix:** Unfreeze the affine policy to per-sector (a_x, b_x) with b_x < 0 fit to Table I.

### [MAJOR] Uniform ×4 normalization clash: P_L/P_R Wilsons contracted with (1∓γ5)-convention (Becirevic-Villadoro) matrix elements
- **File:** eft_deltaf2/hadronic.py:762-769 (Q1 = (8/3)f²m²B_K), 2233-2249 (Q4 = 2Rm²f²B4, Q5 = (2/3)Rm²f²B5); operators.py:46-55, 477-513; combined in observables.py:553-562, 878-921
- **Category:** factor-of-2 (×4) / convention-inconsistency
- **Claim:** Matching produces SUSY/CFW-basis coefficients for P_L/P_R operators (C1/(g_L²/6M²) = 1.000 exactly, executed), but hadronic layer transcribes BV Eq. (5) MEs normalized for (1−γ5)⊗(1±γ5) operators = 4× the P_L/P_R operators — every M12_NP (K, Bd, Bs, D0, LR path) overestimated by exactly 4.
- **Evidence:** Executed: Q1_ME_code/[(2/3)m_K²f_K²B_K] = 4.0; Q4_ME_code/[(R/2)B4·m_K²f_K²] = 4.0 (R_χ = 26.085 ✓). Lane B's verified pairing (deltaf2.py:722-728) corresponds to (2/3)m²f²B for the same C1 — Lane C is 4× Lane B for identical C's. operators.py's own projector_normalization_note contradicts the hadronic formulas.
- **Confidence:** high (executed; cross-checked vs Lane B)
- **Fix:** Divide the three ME formulas by 4, or re-declare/re-match in BV normalization consistently.

### [MAJOR] "strict_paper" Eq. (3) example fails to reproduce Table I (residual 0.16 in c) — failure frozen into tests as expected
- **File:** model.py:1845-1929; inputs.py:248-251 (a=0.8, r=0.3, θ=(115°,65°,70°), δ=0.6); pinned at tests/test_paper_model.py:398-410
- **Category:** paper-transcription / circular-validation
- **Claim:** With quoted example parameters, a·diag(r·V†C_uV + C_d) = (0.516, 0.604, 0.620) vs c_Q = (0.64, 0.59, 0.46) — max residual 0.160 (f_IR exponential in c ⇒ orders of magnitude in profile space), rhs ordering inverted; test suite asserts this exact residual vector rather than resolving it.
- **Evidence:** Both dagger conventions executed (residuals 0.124/−0.014/−0.160 and 0.117/−0.005/−0.162); test pins [0.12379469, −0.01376433, −0.16003036]. arXiv fetch 403-blocked — transcription checks internal-consistency only.
- **Confidence:** high (mismatch), medium (root cause)
- **Fix:** Re-extract Eq. (3) and its example from the paper.

### [MAJOR] Fully circular validation: frozen constants, golden files, verifier, and the RG test all pin the pipeline's own buggy outputs
- **File:** artifacts.py:39-52; tests/golden/paper_0710_1869/default_kaon_np_only/wilsons.json; tests/test_paper_eft_rg_lo.py:777-781, 1049-1068; verifier.py:893-916, :92
- **Category:** circular-validation / test-bug
- **Claim:** Frozen RG constant = matching × 1.372 (the inverted LO factor; correct ≈ 0.746); frozen M12 embeds the ×4 ME; the RG regression test independently re-derives the same wrong factor (α_low/α_high)^{+6/(33−2nf)} as its oracle at rel=1e-9; the "independent" verifier recomputes the identical (8/3) ME formula and hard-rejects any nonzero C4/C5 (|C4|>1e-30 fails) — wiring in the dominant LR physics breaks verification; benchmark-fit tests set targets = the probe's own outputs.
- **Confidence:** high (executed golden ratio; read test helper)
- **Fix:** After fixing rg.py/hadronic.py, regenerate golden/frozen constants; replace the test oracle with an independent reference (BBL η factors or Lane B's verified running).

### [MINOR] No mass/CKM targets or scale semantics exist in Lane C — fit quality and quark-mass scale undefined
- **File:** fit.py:361-495; scan.py:90-94
- **Claim:** No PDG targets, no scale tag; masses from M = 2v·F·Ȳ·F (~3 TeV tree relation) compared to caller-supplied numbers with no QCD running.
- **Fix:** Frozen target bundle (PDG 2024 MS-bar run to μ_match) with explicit scale id.

### [MINOR] f_IR silently returns 0 (then hard-errors) for c ≳ 10.8 — ε^(1−2c) float64 overflow; clamp masks inf→0. Fix: log-space or asymptotic branch.

### [NOTE] Scoped questions answered
- right_down exactly diagonal (default): benchmark policy `diagonal_only_until_paper_owned_rotations_exist.v1` — the *physical* builder does produce rotated right_down, but default ε_K matching consumes the benchmark object, and the verifier rejects the nonzero C4 the physical path would give.
- Suspected C_u = Y†Y misuse: NOT a bug (model.py:1622-1626 uses correct LH/RH conventions).
- scales.py m_g1 = Λ_IR: corroborated.
- ε_K/κ_ε/Im extraction: ε_K not computed anywhere in Lane C; observables end at M12_NP and Δm_NP = 2·Re M12; sign of Im M12 undeclared — moot until wired.

## Verified correct
f_IR for c > 1/2 (signs flip together; c=1/2 limit correct); V5KM/CKM parametrization (PDG form, unitarity 1e-10); basis wiring (V5KM† diag(f_Q²) V5KM transport; U†diag(F²)U internally consistent); matching prefactor 1/M², C1 = g_L²/(6M²) = 1.000 executed; B̂_K→B_K(2 GeV) conversion 0.7625→0.5557 correct direction (ironically opposite the rg.py Wilson evolution); R_χ = 26.085 with PDG-2024 masses; B4=0.78(3), B5=0.57(4) ETM 2013 plausible; Table I (c,f) pairs mutually coherent to rounding; mass-probe CKM/SVD checks correct.

## Not reviewed
artifacts.py beyond frozen constants (~2300 lines serialization); verifier.py provenance sections; conventions.py, validation.py, rg_inputs.py; B_d/B_s/D0 hadronic bundle internals; observables.py custom-combined-total path; most tests beyond the three named; direct arXiv 0710.1869 comparison (proxy 403).
