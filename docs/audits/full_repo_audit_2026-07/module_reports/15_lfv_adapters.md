# Report 15 — Gap-closing: LFV physics modules (mu_e_conversion, lfv_three_body, higgs_lfv, adapters, L002-L023, E002/E006-E009)

**Structural summary.** Three layers: reusable cores in quarkConstraints/, thin physics_adapters/ wrappers splicing in the L001 dipole and Phase-4a rs_ew_couplings tree contacts, and per-process constraint bodies loading YAML anchors + pass/fail wiring. μN→eN core is a faithful KKO implementation (overlaps, capture rates, dipole normalization verified numerically). μ→3e core reproduces Kuno–Okada for contact and dipole terms (verified by an independent spin-summed amplitude Monte Carlo) but has a wrong dipole–contact interference coefficient. The shared three-body core is silently wrong for τ→3ℓ normalization (missing BR(τ→ℓνν̄)) → L009/L010; the same missing factor afflicts the L007/L008 τ→ℓγ prefactor reuse. Stub adapters honest; tests self-referential (encode the interference bug).

### [MAJOR] μ→3e dipole–contact interference coefficient is 2√2·e; Kuno–Okada value is 8e
- **File:** quarkConstraints/lfv_three_body.py:710 (convention strings :51-59, :435-442)
- **Category:** wrong-formula (factor 2√2 ≈ 2.83)
- **Claim:** I_dc coefficient should be **8e** (KO RMP 73 (2001) 151, eq. 2.14: −8e·Re[A_R(2g₄+g₆)* + A_L(2g₃+g₅)*]); both the explicit-amplitude interference and unknown-phase envelope are ×2.83 too small.
- **Evidence:** Independent spin-summed Dirac-amplitude MC (flat Dalitz, physical m_e, identical-particle antisymmetrization) in the module's own normalization: contact-LL coefficient 2.005 (KO: 2), contact-LR 1.005 (KO: 1), dipole² 24.6 vs KO leading-log 23.2 (normalization validated), interference K = −7.999/−8.03 (|K| = 8; A_R↔(2G_LL+G_LR) pairing confirmed; sign opposite the code's `+`).
- **Confidence:** high (magnitude); sign convention-dependent but differs from KO.
- **Fix:** 2√2·e → 8e (and pin the sign against KO).

### [MAJOR] τ→3μ / τ→3e contact + interference terms missing BR(τ→ℓνν̄) ≈ 0.177 — BR overestimated ×5.65
- **File:** quarkConstraints/lfv_three_body.py:689-690, 358-390 (consumed by lfv_three_body_tau.py:163-170 and _taue.py → L009/L010, both HARD)
- **Category:** wrong-formula / units (normalization)
- **Claim:** KO-style contact rate is a width in units of G_F²m_i⁵/192π³; equals the parent total width only for the muon. For τ, BR needs × BR(τ→eνν̄) ≈ 0.177, omitted (no τ lifetime/leptonic-BR input exists).
- **Evidence:** G_F²m_τ⁵/192π³ = 4.01e-13 GeV vs Γ_τ = 2.27e-12 GeV → 0.177; overestimate ×5.65 (over-excluding). Dipole term (BR-ratio × (α/3π)(ln−11/4)) unaffected; dipole-amplitude normalization also muon-specific in interference.
- **Confidence:** high
- **Fix:** Multiply contact (and fixed interference) by BR(ℓ_i→ℓνν̄) as an input field.

### [MAJOR] L007/L008 τ→ℓγ reuse the muon NDA prefactor 4.0e-8 without the τ leptonic-BR factor
- **File:** L007.py/L008.py (load L001.yaml prefac_br), adapter lepton_tau_mu.py:248
- **Category:** wrong-formula (missing ~0.174, inside a documented proxy)
- **Claim:** BR_NP(τ→ℓγ) = prefac_br·|(YY†)|²·(3 TeV/M_KK)⁴ reuses the μ→eγ BR prefactor; τ needs an extra BR(τ→ℓνν̄) ≈ 0.17-0.18 → L007/L008 (HARD) overestimate ×~5.6.
- **Confidence:** medium-high (documented proxy, but the factor is unambiguous within the proxy's own logic)
- **Fix:** Multiply the τ radiative prefactor by BR(τ→ℓνν̄).

### [MINOR] Opposite unknown-phase envelope policies between μN→eN and ℓ→3ℓ′ (+ "conservative" mislabel)
- **File:** mu_e_conversion.py:345-369,472-476 vs lfv_three_body.py:740-754; L004.py:19-21, L005.py:20-22
- **Category:** convention-inconsistency / statistics
- **Claim:** For the same unknown dipole–contact phase, μN→eN grades on the destructive **lower** envelope (lenient) while μ→3e/τ→3ℓ grade on the constructive **upper** envelope (most excluding); "conservative" wording in L004/L005 is the wrong direction for an upper-limit veto.
- **Confidence:** high (behavior); the right policy is a judgment call.
- **Fix:** One declared envelope policy across L002-L005, L009, L010.

### [MINOR] Tests re-implement the buggy interference coefficient instead of pinning literature values
- **File:** tests/constraints/primary/charged_lepton/test_L002.py:320-330; test_L009:274
- **Category:** test-bug — the 2√2 coefficient and τ normalization errors are structurally invisible to the suite.
- **Fix:** Literature-anchored assertions (BR(μ→3e)/BR(μ→eγ)=6.1e-3 dipole-dominance, CR(Al)/BR=2.7e-3, KO interference 8e).

### [NOTE] `_is_live_complex` treats non-numeric values as "live" (adapter :412-421) → garbage extras degrade constraints to NOT-EVALUATED silently. Fix: fail loudly.
### [NOTE] μN→eN unknown-phase envelope is Cauchy–Schwarz-loose (superset interval, benign).

## Verified correct
- KKO nuclear inputs EXACT vs KKO Table (Al/Ti/Au D, V^p/n, S^p/n, ω_capt); ω_conv structure and A_R*↔g_L pairing; dipole-only CR/BR = 0.00267 (Al)/0.00414 (Ti)/0.00393 (Au) matches KKO dipole-dominance ratios (run live).
- 384π² dipole normalization; μ→3e contact coefficients and dipole-dominance factor reproduced by MC to <0.5%.
- Adapter KKO matching g_LV(q)=(C_LL+C_RL)/(√2G_F) incl. the ½ chiral→vector projection; nucleon counting g^p=2g^u+g^d.
- Z-penguin contact G_AB = 2·δg_A·g_B; SM Z couplings.
- zpole_lfv BR formula incl. charge-state factor 2; inverse algebra.
- higgs_lfv Γ = m_h(|Y_ij|²+|Y_ji|²)/(8π) (Harnik-Kopp-Zupan), m_h=125.25, Γ_h=4.07 MeV.
- L006 muonium calibration consistent with Willmann (V∓A) map.
- L023 trident: C_V^SM = 1+4s²_W, ΔC = 2v²g_νg/M² verified by SM-Z self-consistency; direction-aware two-sided Gaussian budget.
- L002-L005 wiring (anchors, Au/Ti assertions, non-vetoing unevaluated paths, Mu2e projection correctly SOFT).
- E002 unit conversion; E006/E007 INFO; E008/E009 honest stubs; vbs_longitudinal clean.

## Not reviewed
lfv_three_body_taue.py full text (findings 1-2 apply via shared core), lepton_tau_e.py beyond prefactor, zpole_lfv_etau/mutau twins, anchor-loader bodies of L006-L010/L023/E006-E009 (anchors pre-verified), physics_adapters/edm.py internals, most test files, quarkConstraints/zpole.py weights (prior audit).
