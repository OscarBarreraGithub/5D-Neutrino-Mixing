# Report 17 — Gap-closing: remaining B/K/charm constraint bodies + rare-decay adapter internals

**Structural summary.** The constraint bodies (B011–B021, K009–K021, C005–C008, B005–B008, B013/B014) are thin anchor+budget wrappers; nearly all physics lives in adapters, and bodies correctly document proxy status. Two new quantitative adapter bugs: the B→K* dilepton kernel is dimensionally wrong (missing q², missing transverse factor 2 → absolute BR ~6× low, cancels in its only consumer B019's LFU ratio), and the R(D) proxy Wilson is dimensionally inconsistent (GeV², ~7× inflated stress). Λb→Λμμ, b→sνν̄, D→πℓℓ (exact Bobeth match), and all LFV charm/kaon phase-space kernels verified correct, several numerically. test_B019 pins the buggy K* kernel by duplicating it verbatim. Radiative B013/B014 correctly use |C7|²+|C7′|² without interference.

### [MAJOR] B→K*μμ kernel missing q² factor and transverse factor 2 — absolute BR ~6× too small
- **File:** quarkConstraints/rare_b_kstar_dilepton.py:305, 353-367
- **Category:** wrong-formula | factor-of-2
- **Claim:** Differential rate is dimensionally GeV⁻⁴ instead of GeV⁻² (missing the standard q² in the transversity normalization N² ∝ q²√λβ/(3·2¹⁰π⁵m_B³)), and `h_perp²+h_parallel²+h_long²` omits the factor 2 on transverse amplitudes (A_⊥,∥ = √2·N·h).
- **Evidence:** As coded, SM BR(B0→K*0μμ) full-range = 1.16e-7, avg dBR/dq²[1.1,6] = 7.5e-9 (PDG ≈ 9.4e-7; LHCb ≈ 3.4e-8). With q²·(2h⊥²+2h∥²+h_long²): BR = 6.87e-7, avg = 2.7e-8 — correct order given C7/charm omitted. Cancels in B019 μ/e ratio; exported absolute BRs wrong.
- **Confidence:** high
- **Fix:** Multiply integrand by q²; weight transverse terms by 2.

### [MAJOR] Semileptonic-LFU proxy Wilson is dimensionful (GeV²) — double normalization
- **File:** quarkConstraints/semileptonic_lfu.py:191-193
- **Category:** units
- **Claim:** C_τ^proxy = ξ·m_b·m_τ/(2√2 G_F M_KK²) has dimension GeV² (m_b·m_τ/M² is already dimensionless; the extra 1/G_F is a second incompatible normalization), yet added to 1 in |1+C|².
- **Evidence:** At M_KK=3 TeV, ξ=1: 0.0250 "GeV²" — ×7.43 (= m_b·m_τ) larger than the four-fermion-normalized 3.4e-3, ×3e4 larger than the charged-Higgs-like 8.3e-7. Feeds B025/B026 R(D)/R(D*); over-excluding direction; documented stress proxy but matches neither documented interpretation.
- **Confidence:** high
- **Fix:** Drop either the m_b·m_τ or the 1/(2√2G_F).

### [MAJOR] B→K* uses (C9+C9′)/(C10+C10′) for ALL transversity amplitudes — primed-Wilson response wrong sign for dominant amplitudes
- **File:** quarkConstraints/rare_b_kstar_dilepton.py:453-456
- **Category:** sign | convention-inconsistency
- **Claim:** A_∥, A_0 (dominating ~80% of the SM rate) carry (C−C′); code applies the pseudoscalar (C+C′) combination to the whole rate, so sizable C9′/C10′ moves R_K* the wrong direction. Does NOT cancel in the B019 ratio (enters the NP response). Documented NEEDS-HUMAN-PHYSICS compression, but the anti-aligned choice.
- **Confidence:** high physics; medium impact (primed Wilsons RS-suppressed but nonzero)
- **Fix:** Split ⊥ vs (∥,0) blocks with (C±C′), mirroring rare_b_nunu's κ_η.

### [MAJOR] test_B019 "manual" oracle is a verbatim copy of the buggy K* kernel
- **File:** tests/constraints/primary/beauty/test_B019.py:100-170
- **Category:** test-bug
- **Claim:** `_manual_b_to_kstar_mumu` duplicates the adapter's missing-q²/missing-2 formula line for line (same 3072π⁵ prefactor), so the suite pins the wrong normalization.
- **Confidence:** high
- **Fix:** Independent oracle (literature SM BR bin) or assert the LFU-ratio response only.

### [MINOR] K009 muon direct-CP collapses ISU vector/axial phase-space weights to a single coefficient
- **File:** physics_adapters/rare_kaon_dilepton_muon.py:98-157
- **Category:** wrong-formula (proxy-adjacent)
- **Claim:** direct-CP weights |y7V|² and |y7A|² equally; ISU K_L→π0μμ vector/axial coefficients differ O(2) (axial helicity-enhanced by m_μ) — axial-only NP mis-weighted; SM limit preserved by construction. Same family as the known K008 mismatch.
- **Confidence:** medium-high
- **Fix:** Separate C_dir^{V,μ}, C_dir^{A,μ} in K009.yaml + adapter.

### [NOTE] K_L LFV modes: K_L CP projection (Y±Y′*)/√2 and K0→π0 isospin 1/√2 not applied (O(1-2) factors; currently zero-impact — tree LFV rigorously zero for diagonal lepton fit; documented proxy).
### [NOTE] `kplus_lifetime_s` field holds the K_L lifetime in K021 config (value correct, field misnamed — latent copy-paste hazard). Fix: rename parent_lifetime_s.
### [NOTE] B012 budget (bare 1σ_exp, SM≡measurement) stricter than and inconsistent with B011/B013 siblings — tightest radiative veto rests on the crudest normalization.
### [NOTE] Anchor loaders monkeypatch the global scaffold in B015/B016/B021/K019 (not thread-safe; serial today).

## Verified correct (several numerically)
- rare_b_baryon_dilepton (Λb→Λμμ): kernels dimensionally consistent; SM avg dBR/dq²[15,20] = 1.106e-7 vs Detmold-Meinel/LHCb (1.0–1.2)e-7 — numerical match; Detmold-2013 dipole parameters correct; (C+C′) compression documented.
- rare_b_nunu: g_SM² identity verified numerically; B→K (1−2η)ε², B→K* (1+1.31η)ε² with correct η; neutrino-flavor factor 3 exactly once; HPQCD LD split correct.
- rare_charm_semileptonic (D+→π+μμ): matches Bobeth-Hiller-Piranishvili term-by-term (2× bracket compensates Γ0 convention); (C+C′) correct for P→P; τ_D+ ✓; full-q² policy matches C007 semantics.
- rare_charm_lfv_dilepton (D0→eμ): unequal-mass phase space correct; charge factor 2 correct.
- rare_charm_lfv_semileptonic (D+→π+eμ): direct spin-summed lepton-tensor contraction verified; reduces to same-flavor kernel at m_e→m_μ to 1.0000000 (pure C9, pure C10, mixed).
- B011-B014 bodies: untagged |C7|²+|C7′|² (correct); B014 correct b→d (0,2) entry.
- K017 R_K construction; K018 V_us error propagation; K012 K_S→μμ SD reproduces D'Ambrosio–Kitahara ≈1.7e-13.
- K019/K020/K021: charge-state factors correct per mode; lifetimes correct; f₊(0)=0.9698; equal-mass reduction consistent with K006.
- C005-C008 bodies: upper-limit semantics, max-saturation logic, documented caveats.
- B005/B006: A_ΔΓ treatment, |C10−C10′|, direction-aware veto for B006 — sound.
- B015/B016/B021 budgets internally consistent; B015 partonic inclusive shape standard.

## Not reviewed
rare_kaon_dilepton_ks.py and radiative_kaon.py adapter internals (K010/K013 kernels); B007/B008 body internals beyond mass-swap (masses verified; inherit known rare_b_dilepton core issues); K009 body YAML-parsing regex; tests other than test_B019/test_B021; the rare_kaon_lfv √2G_Fα/(πg_SM²) mapping constant (structure consistent with audited K006 convention).
